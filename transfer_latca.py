from os import listdir, mkdir
from os.path import isfile, join
from shutil import copyfile

def main():
    file_of_usables = open("data/strong_purchasable_thru_45.txt")
    usable_latids = [x.strip() for x in file_of_usables.readlines()]
    file_of_usables.close()

    all_sdfs = [obj for obj in listdir("latca/") if isfile(join("latca/", obj))]
    
    try:
        mkdir("latca_strong_purchasable_thru_45/")
    except OSError:
        pass

    found = []
        
    for sdf_name in all_sdfs:
        if sdf_name[:-4] in usable_latids: # [:-4] strips the ".sdf"
            copyfile(join("latca/", sdf_name), join("latca_strong_purchasable_thru_45/", sdf_name))
            found.append(sdf_name[:-4])
            print("Transferred {0}".format(sdf_name))

    for i in usable_latids:
        if not(i in found):
            print(i)

if __name__ == "__main__":
    main()
