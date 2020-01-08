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
    "name": "Import textures as materials",
    "author": "sanderboer",
    "version": (2019, 11, 26),
    "blender": (2, 81, 0),
    "location": "File > Import-Export",
    "description": "Import a dir full of textures as materials",
    "category": "Import-Export",
}

import bpy, bmesh, math
import os
from os import listdir
from os.path import isfile, join
from bpy_extras.io_utils import ImportHelper, path_reference_mode
from bpy.props import StringProperty, BoolProperty, FloatProperty, EnumProperty,CollectionProperty


class ImportTexturesAsMaterials(bpy.types.Operator, ImportHelper):
    bl_idname = "import.textures"
    bl_label = "Import Textures"
    bl_options = {'REGISTER', 'PRESET', 'UNDO'}
    bl_description = bl_info['description']

#    filename_ext = ".png"
#    filter_glob: StringProperty(
#            default="*.jpg;*.tga;*.png",
#            options={'HIDDEN'},
#            )

#    
    # ----------------------
    # File dialog properties
    files: CollectionProperty(type=bpy.types.OperatorFileListElement, options={'HIDDEN', 'SKIP_SAVE'})

    directory: StringProperty(maxlen=1024, subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})

    filter_image: BoolProperty(default=True, options={'HIDDEN', 'SKIP_SAVE'})
    filter_movie: BoolProperty(default=True, options={'HIDDEN', 'SKIP_SAVE'})
    filter_folder: BoolProperty(default=True, options={'HIDDEN', 'SKIP_SAVE'})


    option_add_fake: BoolProperty(name="Persistant materials (add fake user)", default=False)
   
    def create_material_from_img_path(self,filepath):
        img_name= os.path.splitext(os.path.basename(filepath))[0]
        mat = bpy.data.materials.new(name=img_name)
        mat.use_fake_user = self.option_add_fake
        mat.use_nodes = True
        bsdf = mat.node_tree.nodes["Principled BSDF"]
        texImage = mat.node_tree.nodes.new('ShaderNodeTexImage')
        texImage.image = bpy.data.images.load(filepath)
        mat.node_tree.links.new(bsdf.inputs['Base Color'], texImage.outputs['Color'])
 
    def execute(self, context):
        #imgdir="/opt/quakemapping/textures/prototype_1_2"
        #files= listdir(imgdir)
       
#        print('DIR: '+str(self.directory))
#        print('FILES: '.join([ n['name'] for n in self.files]))
#        for n['name'] in self.files:
#            filepa
#        #print("dir: "+self.directory)
#        imgfiles=[]
        image_extensions=[".png",".PNG",".tga",".TGA",".jpg",".JPG"]

        for n in self.files:
            f = n['name']
            ext=os.path.splitext(f)[-1]
            filepath=join(self.directory, f)
            if ext in image_extensions :
                #imgfiles.append(filepath)
                self.create_material_from_img_path(filepath)

        return {'FINISHED'}

def menu_func_export(self, context):
    self.layout.operator(ImportTexturesAsMaterials.bl_idname, text="Batch import textures as materials")

def register():
    bpy.utils.register_class(ImportTexturesAsMaterials)
    bpy.types.TOPBAR_MT_file_import.append(menu_func_export)

def unregister():
    bpy.utils.unregister_class(ImportTexturesAsMaterials)
    bpy.types.TOPBAR_MT_file_import.remove(menu_func_export)



if __name__ == "__main__":
    #unregister()
    register()
