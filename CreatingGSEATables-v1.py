import numpy as np
import pandas as pd
import re
import xlsxwriter 
import os
import sys

def main (argv):
    if len(argv)!=2:
        usage(argv)
        sys.exit(0)

    for file in os.listdir(argv[1]):
        if file.endswith("AvgXprsn.txt"):
            if file.startswith("."):
                pass
            else:
                clust = pd.read_table(os.path.join(argv[1],file))
                outname = re.sub(".txt$", "", file)
                print ("working on ", outname)
                clust.index = clust.index.str.upper()
                total_clust_count = len(clust.columns.values)
                clust_names = clust.columns.values
                clust.reset_index(inplace=True)
                temp = clust.iloc[:,0]
                clust2 = pd.concat([clust.iloc[:,0], clust], axis=1)
                clust2.columns.values[0]="NAME"
                clust2.columns.values[1]="Description"
                line2 = "\t".join([str(len(clust2.index)), str(len(clust2.columns.values))])
                my_file = open("testfile.gct", "w")
                my_file.write("#1.2\n")
                my_file.write(line2)
                my_file.write("\n")
                my_file.close()
                with open("".join([argv[1], outname, ".gct"]), "a") as f:
                    clust2.to_csv(f, header=True, sep="\t", index=False)
                new_list ='\t'.join(str(e) for e in clust_names)
                line1 = "\t".join([str(total_clust_count), str(total_clust_count), "1"])
                line2 = "\t".join(["#", new_list])
                line3 = new_list
                clsLine = "\n".join([line1, line2, line3])
                my_clsFile = open("".join([argv[1],outname, ".cls"]), "w")
                my_clsFile.write(clsLine)
                my_clsFile.close()
            
            


def usage(argv):
    print("Enter input dir")
    

    
    
if __name__=="__main__":
    main(sys.argv)


