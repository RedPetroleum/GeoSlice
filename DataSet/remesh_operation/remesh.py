import pymeshlab
import glob, os

rootDir = os.getcwd()
print(rootDir)

files = glob.glob(rootDir + "/layers_remeshed/*")
for f in files:
    os.remove(f)

os.chdir(rootDir + "/layers_unremeshed")
for file in glob.glob("*.obj"):
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(file)
    ms.load_filter_script(rootDir + "/remesh_operation.mlx")
    ms.apply_filter_script()
    ms.save_current_mesh(rootDir + "/layers_remeshed/" + file)
    print("remesh " , file, " finished.")