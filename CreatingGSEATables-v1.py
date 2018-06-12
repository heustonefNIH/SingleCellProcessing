import numpy as np
import pandas as pd
import re
import os
import sys

def usage(argv):
	print("Enter input dir")
    
    
def output_prefix (inname): #Define output prefix based on input file name
	outname = re.sub(".txt$", "", inname)
	return outname

def create_system_command(file_directory, gct_file, cls_file, output_label, number_of_clusters):
	
	for current_cluster in range(0, number_of_clusters):
		current_cluster = str(current_cluster)
		system_command = "".join(["java -cp ./gsea-3.0.jar -Xms512m -Xmx8164m xtools.gsea.Gsea -res ./",gct_file," -cls ./",cls_file,"#",current_cluster,"_versus_REST -gmx msigdb.v6.1.symbols.gmt -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute phenotype -rnd_type no_balance -scoring_scheme weighted -rpt_label C",output_label,"_",current_cluster,"_vs_REST -metric log2_Ratio_of_Classes -sort real -order descending -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out ./ -gui false"])
		with open("".join([file_directory, output_label, "_GSEA_submission.swarm"]), "a") as f:
			f.write(system_command)
			f.write("\n")
			f.close()

def main (argv):
	if len(argv)!=2:
		usage(argv)
		sys.exit(0)

	for file in os.listdir(argv[1]):
		if file.endswith("AvgXprsn.txt"):
			if file.startswith("."):
				pass
			else:
				outname = output_prefix(file)
				print ("working on ", outname)
				
				clust = pd.read_table(os.path.join(argv[1],file))
				clust.index = clust.index.str.upper()
				total_clust_count = len(clust.columns.values)
				clust_names = clust.columns.values
				
				clust.reset_index(inplace=True)
				clust2 = pd.concat([clust.iloc[:,0], clust], axis=1)
				clust2.columns.values[0]="NAME"
				clust2.columns.values[1]="Description"
				gct_file = "".join([argv[1], outname, ".gct"])
				gct_line = "\t".join([str(len(clust2.index)), str(total_clust_count)])
				with open(gct_file, "a") as f:
					f.write("#1.2\n")
					f.write(gct_line)
					f.write("\n")
					clust2.to_csv(f, header=True, sep="\t", index=False)
					f.close()

				cls_line1 = "\t".join([str(total_clust_count), str(total_clust_count), "1"])
				cls_line3 = '\t'.join(str(e) for e in clust_names)
				cls_line2 = "\t".join(["#", cls_line3])
				cls_file ="".join([argv[1],outname, ".cls"]) 
				with open(cls_file, "a") as f:
					f.write(cls_line1)
					f.write("\n")
					f.write(cls_line2)
					f.write("\n")
					f.write(cls_line3)
					f.close()
				create_system_command(argv[1], gct_file, cls_file, outname, total_clust_count)
            	
            

if __name__=="__main__":
	main(sys.argv)




