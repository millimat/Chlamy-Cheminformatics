from sys import argv
from os import listdir, mkdir
from os.path import isfile, join
from shutil import copyfile

sourcedir = "latca_all/"

def main():
    if(len(argv) != 3):
        print("usage: python {0} [LATCA listing] [target dir]".format(argv[0]))
        return

    usables_path = argv[1]
    target = argv[2]
    
    file_of_usables = open(usables_path)
    usable_latids = [x.strip() for x in file_of_usables.readlines()]
    file_of_usables.close()

    all_sdfs = [obj for obj in listdir(sourcedir) if isfile(join(sourcedir, obj))]
    
    try:
        mkdir(target)
    except OSError:
        pass

    found = []
        
    for sdf_name in all_sdfs:
        if sdf_name[:-4] in usable_latids: # [:-4] strips the ".sdf"
            copyfile(join(sourcedir, sdf_name), join(target, sdf_name))
            found.append(sdf_name[:-4])
            print("Transferred {0}".format(sdf_name))

    for i in usable_latids:
        if not(i in found):
            print(i)

if __name__ == "__main__":
    main()
