# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# ##### END GPL LICENSE BLOCK #####

bl_info = {
    "name": "Export Quake Map (.map)",
    "author": "chedap",
    "version": (2019, 11, 26),
    "blender": (2, 81, 0),
    "location": "File > Import-Export",
    "description": "Export geometry as brushes",
    "category": "Import-Export",
}

import bpy, bmesh, math
import uuid
import os, os.path
import datetime
from shutil import copyfile
from mathutils import Vector, Matrix
from numpy.linalg import solve
from bpy_extras.io_utils import ExportHelper
from bpy.props import StringProperty, BoolProperty, FloatProperty, EnumProperty

class ExportQuakeMap(bpy.types.Operator, ExportHelper):
    bl_idname = 'export.map'
    bl_label = bl_info['name']
    bl_description = bl_info['description']
    bl_options = {'UNDO'}
    filename_ext = ".map"
    filter_glob: StringProperty(default="*.map", options={'HIDDEN'})

    option_sel: BoolProperty(name="Selection only", default=True)
    #option_tm: BoolProperty(name="Apply transform", default=True)
    option_triangulate: BoolProperty(name="Triangulate faces", default=True)
    option_backface_cull: BoolProperty(name="Apply skip to backfaces", default=True)
    option_geo: EnumProperty(name="Geo", default='Faces',
        items=( ('Brushes', "Brushes", "Export each object as a convex brush"),
                ('Faces', "Faces", "Export each face as a pyramid brush") ) )
    option_grid: FloatProperty(name="Grid", default=4.0,
        description="Snap to grid (0 for off-grid)", min=0.0, max=256.0)
    option_depth: FloatProperty(name="Depth", default=8.0,
        description="Pyramid poke offset", min=0.0, max=256.0)
    option_scale: FloatProperty(name="Scale", default=32.0,
        description="Scale to Quake units", min=0.01, max=4096.0)
    option_format: EnumProperty(name="Format", default='Valve',
        items=( ('Quake', "Standard", "Axis-aligned texture projection"),
                ('Valve', "Valve220", "Face-bound texture projection") ) )
    option_dest: EnumProperty(name="Save to", default='File',
        items=( ('File', "File", "Write data to a .map file"),
                ('Append', "Append", "Append data to an existing .map file"),
                ('Clip', "Clipboard", "Store data in system buffer") ) )
    option_skip: StringProperty(name="Fallback", default='__TB_empty',
        description="Texture to use on new and unassigned faces")

    # https://blender.stackexchange.com/questions/75332/how-to-find-the-number-of-loose-parts-with-blenders-python-api
    def walk_island(self,vert):
        ''' walk all un-tagged linked verts '''    
        vert.tag = True
        yield(vert)
        linked_verts = [e.other_vert(vert) for e in vert.link_edges
                        if not e.other_vert(vert).tag]

        for v in linked_verts:
            if v.tag:
                continue
            yield from self.walk_island(v)

    def get_islands(self,bm, verts=[]):
        def tag(verts, switch):
            for v in verts:
                v.tag = switch
        tag(bm.verts, True)
        tag(verts, False)
        ret = {"islands" : []}
        verts = set(verts)
        while verts:
            v = verts.pop()
            verts.add(v)
            island = set(self.walk_island(v))
            ret["islands"].append(list(island))
            tag(island, False) # remove tag = True
            verts -= island
        return ret
    
    def gridsnap(self, vector):
        grid = self.option_grid
        if grid:
            return [round(co/grid)*grid for co in vector]
        else:
            return vector
                
    def printvec (self, vector):
        return ' '.join([f'{co:.5g}' for co in vector])


    def get_texture_name(self, ob, face, skip=False):
        tex_name = self.option_skip
        if skip:
            return tex_name
        if ob.material_slots != []:
            mat = ob.material_slots[face.material_index].material
        for n in mat.node_tree.nodes:
            if n.type == "TEX_IMAGE":
                filepath=n.image.filepath
                tex_name= os.path.splitext(os.path.basename(filepath))[0]
                img_size_x = n.image.size[0] if n.image.size[0] > 0 else img_size_x
                img_size_y = n.image.size[1] if n.image.size[1] > 0 else img_size_y
                break
        return tex_name

    def texdata(self, face, mesh, obj, skip=False):
        mat = None
        width = height = 64

        texstring = self.get_texture_name(obj,face,skip)

        V = [loop.vert.co for loop in face.loops]
        uv_layer = mesh.loops.layers.uv.active
        if uv_layer is None:
            uv_layer = mesh.loops.layers.uv.new("dummy")
        T = [loop[uv_layer].uv for loop in face.loops]

        # UV handling ported from: https://bitbucket.org/khreathor/obj-2-map
        if self.option_format == 'Valve':
            # [ Ux Uy Uz Uoffs ] [ Vx Vy Vz Voffs ] rotation scaleU scaleV
            dummy = ' [ 1 0 0 0 ] [ 0 -1 0 0 ] 0 1 1\n'

            height = -height # workaround for flipped v

            # Set up "2d world" coordinate system with the 01 edge along X
            world01 = V[1] - V[0]
            world02 = V[2] - V[0]
            world01_02Angle = world01.angle(world02)
            if face.normal.dot(world01.cross(world02)) < 0:
                world01_02Angle = -world01_02Angle
            world01_2d = Vector((world01.length, 0.0))
            world02_2d = Vector((math.cos(world01_02Angle),
                                math.sin(world01_02Angle))) * world02.length

            # Get 01 and 02 vectors in UV space and scale them
            tex01 = T[1] - T[0]
            tex02 = T[2] - T[0]
            tex01.x *= width
            tex02.x *= width
            tex01.y *= height
            tex02.y *= height

            '''
            a = world01_2d
            b = world02_2d
            p = tex01
            q = tex02

            [ px ]   [ m11 m12 0 ] [ ax ]
            [ py ] = [ m21 m22 0 ] [ ay ]
            [ 1  ]   [ 0   0   1 ] [ 1  ]

            [ qx ]   [ m11 m12 0 ] [ bx ]
            [ qy ] = [ m21 m22 0 ] [ by ]
            [ 1  ]   [ 0   0   1 ] [ 1  ]

            px = ax * m11 + ay * m12
            py = ax * m21 + ay * m22
            qx = bx * m11 + by * m12
            qy = bx * m21 + by * m22

            [ px ]   [ ax ay 0  0  ] [ m11 ]
            [ py ] = [ 0  0  ax ay ] [ m12 ]
            [ qx ]   [ bx by 0  0  ] [ m21 ]
            [ qy ]   [ 0  0  bx by ] [ m22 ]
            '''

            # Find an affine transformation to convert 
            # world01_2d and world02_2d to their respective UV coords
            texCoordsVec = Vector((tex01.x, tex01.y, tex02.x, tex02.y))
            world2DMatrix = Matrix(((world01_2d.x, world01_2d.y, 0, 0),
                                    (0, 0, world01_2d.x, world01_2d.y),
                                    (world02_2d.x, world02_2d.y, 0, 0),
                                    (0, 0, world02_2d.x, world02_2d.y)))
            try:
                mCoeffs = solve(world2DMatrix, texCoordsVec)
            except:
                return texstring + dummy
            right_2dworld = Vector(mCoeffs[0:2])
            up_2dworld = Vector(mCoeffs[2:4])

            # These are the final scale values
            # (avoid division by 0 for degenerate or missing UVs)
            scalex = 1 / max(0.00001, right_2dworld.length)
            scaley = 1 / max(0.00001, up_2dworld.length)
            scale = Vector((scalex, scaley))

            # Get the angles of the texture axes. These are in the 2d world
            # coordinate system, so they're relative to the 01 vector
            right_2dworld_angle = math.atan2(right_2dworld.y, right_2dworld.x)
            up_2dworld_angle = math.atan2(up_2dworld.y, up_2dworld.x)

            # Recreate the texture axes in 3d world coordinates,
            # using the angles from the 01 edge
            rt = world01.normalized()
            up = rt.copy()
            rt.rotate(Matrix.Rotation(right_2dworld_angle, 3, face.normal))
            up.rotate(Matrix.Rotation(up_2dworld_angle, 3, face.normal))

            # Now we just need the offsets
            rt_full = rt.to_4d()
            up_full = up.to_4d()
            test_s = V[0].dot(rt) / (width * scale.x)
            test_t = V[0].dot(up) / (height * scale.y)
            rt_full[3] = (T[0].x - test_s) * width
            up_full[3] = (T[0].y - test_t) * height

            texstring += f" [ {self.printvec(rt_full)} ] "\
                        f"[ {self.printvec(up_full)} ]"\
                        f" 0 {self.printvec(scale)}\n"

        elif self.option_format == 'Quake':
            # offsetU offsetV rotation scaleU scaleV
            dummy = ' 0 0 0 1 1\n'

            # 01 and 02 in 3D space
            world01 = V[1] - V[0]
            world02 = V[2] - V[0]

            # 01 and 02 projected along the closest axis
            maxn = max(abs(round(crd,5)) for crd in face.normal)
            for i in [2,0,1]: # axis priority for 45 degree angles
                if round(abs(face.normal[i]),5) == maxn:
                    axis = i
                    break
            world01_2d = Vector((world01[:axis] + world01[(axis+1):]))
            world02_2d = Vector((world02[:axis] + world02[(axis+1):]))

            # 01 and 02 in UV space (scaled to texture size)
            tex01 = T[1] - T[0]
            tex02 = T[2] - T[0]
            tex01.x *= width
            tex02.x *= width
            tex01.y *= height
            tex02.y *= height
            
            # Find affine transformation between 2D and UV
            texCoordsVec = Vector((tex01.x, tex01.y, tex02.x, tex02.y))
            world2DMatrix = Matrix(((world01_2d.x, world01_2d.y, 0, 0),
                                    (0, 0, world01_2d.x, world01_2d.y),
                                    (world02_2d.x, world02_2d.y, 0, 0),
                                    (0, 0, world02_2d.x, world02_2d.y)))
            try:
                mCoeffs = solve(world2DMatrix, texCoordsVec)
            except:
                return texstring + dummy

            # Build the transformation matrix and decompose it
            tformMtx = Matrix(( (mCoeffs[0], mCoeffs[1], 0),
                                (mCoeffs[2], mCoeffs[3], 0),
                                (0,          0,          1) ))
            t0 = Vector((T[0].x * width, T[0].y * height)).to_3d()
            v0 = Vector((V[0][:axis] + V[0][(axis+1):])).to_3d()

            offset = t0 - ( tformMtx @ v0 )
            rotation = math.degrees(tformMtx.inverted_safe().to_euler().z)
            scale = tformMtx.inverted_safe().to_scale() # always positive

            # Compare normals between UV and projection to get the scale sign
            tn = tex01.to_3d().cross(tex02.to_3d())
            vn = world01_2d.to_3d().cross(world02_2d.to_3d())
            if tn.dot(vn) < 0: scale.x *= -1

            # fudge
            offset.x += width
            offset.y *= -1

            finvals = [offset.x, offset.y, rotation, scale.x, scale.y]
            texstring += f" {self.printvec(finvals)}\n"

        return texstring

    def prevent_overwrite(self, filepath):
        t = datetime.datetime.now()
        path = os.path.dirname(os.path.abspath(filepath))
        basename = os.path.basename(filepath)
        filename = ''.join(basename.split(".")[0:-1])
        ext = basename.split(".")[-1]
        backup_filepath = path + "/" + filename + t.strftime("_bak_%y%m%d_%H%M") + "." +ext
        if os.path.exists(filepath):
            copyfile(filepath, backup_filepath)
    
    def execute(self, context):
        if self.option_sel:
            objects = context.selected_objects
        else:
            objects = context.scene.objects
        objects = [obj for obj in objects if obj.type == 'MESH']

        geo = []
        fw = geo.append
        if not self.option_dest == 'Append':
            fw('// Game: Quake\n')
            if self.option_format == 'Valve':
                fw('// Format: Valve\n')
                fw('"mapversion" "220"\n')
            else:
                fw('// Format: Quake\n')
            fw('{\n"classname" "worldspawn"\n')
            fw('}\n')

        group_id=uuid.uuid4().int
        if self.option_geo == 'Faces' and objects != []:
            for obj in objects:
                ob = obj.evaluated_get(bpy.context.evaluated_depsgraph_get())
                bm = bmesh.new()
                bm.from_mesh(ob.data)
                bmesh.ops.connect_verts_concave(bm, faces=bm.faces)
                if self.option_triangulate:
                    bmesh.ops.triangulate(bm, faces=bm.faces)
                else:
                    bmesh.ops.connect_verts_concave(bm, faces=bm.faces)      
                bmesh.ops.transform(bm, matrix=obj.matrix_world*self.option_scale, verts=bm.verts)
                for vert in bm.verts:
                    vert.co = self.gridsnap(vert.co)
                bm.faces.ensure_lookup_table()

                islands = [island for island in self.get_islands(bm, verts=bm.verts)["islands"]]
                print(ob.name, "Islands:", len(islands))

                facegroups=[]
                for i in islands:
                    faces=[]
                    for f in bm.faces:
                        is_in=True
                        for v in f.verts:
                            if v not in i:
                                is_in=False
                        if is_in:    
                            faces.append(f)    
                    if len(faces)>0:
                        facegroups.append(faces)

                if len(facegroups) < 1:
                    facegroups=[bm.faces]
                for num,facegroup in enumerate(facegroups):
                    fw("{\n")
                    fw('"classname" "func_group"\n')
                    fw('"_phong" "1"\n')
                    fw('"_tb_type" "_tb_group"\n')
                    fw('"_tb_name" "' + ob.name + '_'+ str(num) +'"\n')
                    fw('"_tb_id" "' + str(group_id) +'"\n')
                    for face in facegroup[:]:
                        if face.calc_area() < 0.001:
                            continue
                        fw('//brush from face from object: ' + obj.name + ' sub: '+ str(num)+'\n')
                        fw('{\n')
                        for vert in reversed(face.verts[0:3]):
                            fw(f'( {self.printvec(vert.co)} ) ')
                        fw(self.texdata(face, bm, obj))
                        pyr = bmesh.ops.poke(bm, faces=[face],
                                     offset=-self.option_depth)
                        apex = pyr['verts'][0].co
                        pyr['verts'][0].co = self.gridsnap(apex)
                        for pyrface in pyr['faces']:
                            for vert in pyrface.verts[0:]: # backfacing
                                fw(f'( {self.printvec(vert.co)} ) ')
                            pyrface.material_index = len(obj.data.materials) - 1
                            fw(self.texdata(pyrface, bm, obj, skip=True))
                        fw('}\n') # end face
                        group_id+=1
                    fw('}\n') # end group
                    group_id+=1
                ob.to_mesh_clear()

        elif self.option_geo == 'Brushes':
            fw("{\n")
            fw('"classname" "func_group"\n')
            fw('"_phong" "1"\n')
            fw('"_tb_type" "_tb_group"\n')
            fw('"_tb_name" "' + ob.name + '_'+ str(num) +'"\n')
            fw('"_tb_id" "' + str(group_id) +'"\n')

            for obj in objects:
                bm.from_mesh(obj.data)
                #if self.option_tm:
                #bmesh.ops.transform(bm, matrix=obj.matrix_world,
                #                                    verts=bm.verts)
                bm.transform(ob.matrix_world * self.option_scale)
 
                for vert in bm.verts:
                    vert.co = self.gridsnap(vert.co)
                hull = bmesh.ops.convex_hull(bm, input=bm.verts,
                                        use_existing_faces=True)
                geom = hull['geom'] + hull['geom_holes']
                oldfaces = [face for face in bm.faces if face not in geom]
                bmesh.ops.delete(bm, geom=oldfaces, context='FACES')
                bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
                bmesh.ops.join_triangles(bm, faces=bm.faces,
                    angle_face_threshold=0.01, angle_shape_threshold=0.7)
                bmesh.ops.connect_verts_nonplanar(bm, faces=bm.faces,
                                                    angle_limit=0.0)
                fw('{\n')
                for face in bm.faces:
                    for vert in reversed(face.verts[0:3]):
                        fw(f'( {self.printvec(vert.co)} ) ')
                    fw(self.texdata(face, bm, obj))
                fw('}\n')
                bm.clear()
            fw('}\n')
        bm.free()
        #fw('}')

        if self.option_dest == 'File':
            self.prevent_overwrite(self.filepath)
            with open(self.filepath, 'w') as file:
                    file.write(''.join(geo))
        elif self.option_dest == 'Append':
            self.prevent_overwrite(self.filepath)
            with open(self.filepath, 'a') as file:
                file.write(''.join(geo))
        else:
            bpy.context.window_manager.clipboard = ''.join(geo)

        return {'FINISHED'}
    
def menu_func_export(self, context):
    self.layout.operator(ExportQuakeMap.bl_idname, text="Quake Map (.map)")

def register():
    bpy.utils.register_class(ExportQuakeMap)
    bpy.types.TOPBAR_MT_file_export.append(menu_func_export)

def unregister():
    bpy.utils.unregister_class(ExportQuakeMap)
    bpy.types.TOPBAR_MT_file_export.remove(menu_func_export)
