# superseeker/data_processing.py

import pandas as pd
import matplotlib.pyplot as plt

def get_CN_info(samples, facets_dir, patient_sex): # For vcf_to_pyclone_input()
    output_list = []
    x_cn_list = dict()
    for samp in samples:
        ## open facets cncf file for that sample
        cn_file = open(facets_dir+"/"+samp+".cncf.txt", "r")
        entry = []
        for line in cn_file:
            info = line.strip().split("\t")
            chrom = info[0]
            if chrom == "chrom":
                continue
            start = info[9]
            end = info[10]
            cell_fraction = info[11]
            total_cn = info[12]
            minor_cn = info[13]
            # If it's chromosome 23, then get the distance (in bp) of the entry, and then add it to a dictionary. 
            # Once all of the samples have been processed, this will find the total_copy_number that spans the 
            # most distance, and assume that is the normal copy number for the X chromosome. 
            if chrom == "23" and patient_sex != "M" and patient_sex != "Y": 
                dist = int(end) - int(start)
                if total_cn in x_cn_list.keys():
                    x_cn_list[total_cn] = x_cn_list[total_cn] + dist
                else:
                    x_cn_list[total_cn] = dist
            if minor_cn == "NA":
                continue
            major_cn = str(int(total_cn)-int(minor_cn))
            if major_cn == "0":
                major_cn = "1"
                minor_cn = "0"
            entry.append([chrom,start,end,cell_fraction,total_cn,major_cn,minor_cn])
        output_list.append(entry)
        cn_file.close()
    if patient_sex == "M":
        X_cn = "1"
    elif patient_sex == "F" or len(x_cn_list.keys()) == 0:  
        X_cn = "2"
    else:
        X_cn = max(x_cn_list, key=x_cn_list.get)
    return output_list, X_cn

def get_copy_numbers(sample_cn_data, chr, position, X_normal_cn = 2): # For vcf_to_pyclone_input()
    major = "1" 
    minor = "1"
    
    if len(chr) > 3 and chr[:3] == "chr":
        chr = chr[3:]
    for val in sample_cn_data:
        if val[0] != chr:
            continue
        if int(val[1]) < int(position) and int(val[2]) > int(position): #If start index is less than position, and end index is greater than:
            major = val[5]
            minor = val[6]
            return major, minor
    if (chr == "X"):
        minor = str(int(X_normal_cn)-1)
    return major, minor

def get_variant_output_lines(variant_line, samples, output_file, HIGH_IMPACT=False, germfilter=False, cn_override=False, cn_neutral=False, patient_sex="F", X_normal_cn=2, samples_cn_lists=[]): # For vcf_to_pyclone_input()
    fields = variant_line.strip().split("\t")
    mutation_id = fields[0]+":"+fields[1]+":"+fields[4]
    info = fields[7].split(";")
    if HIGH_IMPACT:
        for i in reversed(info):
            if i[:4] == "ANN=":
                ann = i.split("|")
                for a in ann:
                    if a == "HIGH" or a == "MODERATE":
                        break
                    if a == "MODIFIER" or a == "LOW":
                        return
                break
    Format = fields[8].split(':')
    AO_index = Format.index("AO")
    RO_index = Format.index("RO")

    # Calculate the germline AF of the variant, and exlude the variant if the 
    # AF is too high in the germline variant. This is to remove likely artifacts
    if germfilter:
        germline_alt = fields[9].split(":")[AO_index]
        germline_ref = fields[9].split(":")[RO_index]
        Germline_AF = int(germline_alt)/(int(germline_alt)+int(germline_ref))
        if Germline_AF >= 0.01:# and Germline_AF <=0.2: 
            pass
            #return
    output = ""
    i = 0
    while i < len(samples):
        if germfilter:
            sample_info = fields[10+i].split(":")
        else:
            sample_info = fields[9+i].split(":")
        sample_id = samples[i]
        ref_counts = sample_info[RO_index]
        alt_counts = sample_info[AO_index]
        
        ## You could figure out how to filter germline AFs that are high. You need to figure out how to make sure it's the germline tho.
        if cn_override:
            major_cn, minor_cn= "1","1"
            if fields[0] == "chrX" and patient_sex == "M":
                normal_cn = "1"
            else:
                normal_cn = "2"
        else:
            major_cn, minor_cn = get_copy_numbers(samples_cn_lists[i], fields[0], fields[1]) #passing in the list of lists, chr, and position
            if fields[0] == "chrX":
                normal_cn = X_normal_cn # This is a rough way to estimate. There may be a better option.
            else:
                normal_cn = "2" 
            if cn_neutral: ## If it has been specified that no variants in CNV regions should be included, this will add the minor and major CNs 
                    ## together. If it doesn't equal the expected normal copy number, skip it.
                if int(major_cn)+int(minor_cn) != int(normal_cn):
                    print("skipping: "+str(",".join([mutation_id, sample_id, ref_counts, alt_counts, major_cn, minor_cn, normal_cn])))
                    return
        output = output + str("\t".join([mutation_id, sample_id, ref_counts, alt_counts, major_cn, minor_cn, normal_cn])+"\n")
        i = i+1
    output_file.write(output)
    return

def vcf_to_pyclone_input(vcf_file_name, facets_dir, output_file_name, patient_sex, cn_neutral, cn_override, germfilter):
    ## Input: Somatic VCF file (merged?) and facets cncf file for each sample, patient sex (F or M) is optional.
    ## Important: Right now the vcf file must be decompressed!
    ## Output: Tab-delimited file with the following columns:
    # 1. mutation_id - Unique identifier for the mutation. This is free form but should match across all samples.
    # 2. sample_id - Unique identifier for the sample.
    # 3. ref_counts - Number of reads matching the reference allele.
    # 4. alt_counts - Number of reads matching the alternate allele.
    # 5. major_cn - Major copy number of segment overlapping mutation.
    # 6. minor_cn - Minor copy number of segment overlapping mutation.
    # 7. normal_cn - Total copy number of segment in healthy tissue. For autosome this will be two and male sex chromosomes one.
    # Optional:
    # 8. tumour_content - The tumour content (cellularity) of the sample. Default value is 1.0 if column is not present.
    # 9. error_rate - Sequencing error rate. Default value is 0.001 if column is not present.

    ## Steps ##
    # 1. Read in VCF file, and get the sample names that are included in VCF
    # 2. Read in the facets cncf files for each sample one at at time
    # 3. For each sample, store the CNV info for each region. Include the start and end, and the major/minor info.
    # 4. For each variant in the VCF file, iterate through each of the sample read count values. Add one line
    #    to the output per sample, indicating which sample it came from and inputting the needed info
    # 5. Once each sample value for that variant has a line in the output, move on to the next variant
    # 6. Make sure that Y chromosomes are getting a Normal CN of 1. Everything else is 2.
  
    HIGH_IMPACT = False #I think all variants are needed for clustering.
    ## Step 1. ##
    if vcf_file_name[-7:] == ".vcf.gz": ## This is not working yet. For some reason it adds b' ' to every line.
        #vcf_file = gzip.open(vcf_file_name, "rb")
        print("VCF file is compressed. Please decompress and try again.")
    elif vcf_file_name[-4:] == ".vcf":
        vcf_file = open(vcf_file_name, "r")
    else:
        print ("Improper VCF file format. Exiting...")
        exit()
    
    output_file = open(output_file_name, "w")

    ## Step 2. ##
    variant_lines = [] # This will hold the line of each variant in the VCF file for the patient.
    for line in vcf_file:
        if line[:2] == "##":
            continue
        if line[0] == "#":
            header_line = line.strip()
        else: 
            variant_lines.append(line.strip())

    columns = header_line.split("\t")
    if germfilter:
        germline_sample = columns[9]
        samples = columns[10:] # I think this way of getting sample IDs will only work for the CLL workflow. It will need to be more robust in future. Starts at 10 to skip germline sample.
    else: 
        samples = columns[9:]
    print("Samples being processed: ")
    print(samples)

    samples_cn_lists = []
    X_normal_cn = 2
    ## Step 3 ##
    if not cn_override:
        samples_cn_lists, X_normal_cn = get_CN_info(samples, facets_dir, patient_sex)

    ## print to ouput file.
    output_file.write("mutation_id\tsample_id\tref_counts\talt_counts\tmajor_cn\tminor_cn\tnormal_cn\n")
    for line in variant_lines:
        get_variant_output_lines(line, samples, output_file, HIGH_IMPACT, germfilter, cn_override, cn_neutral, patient_sex, X_normal_cn, samples_cn_lists)

    output_file.close()
    vcf_file.close()

######################################

def pyclone_to_vcf(vcf_file, clustered_file, output_file_name, HIGH_IMPACT=False):
    # The pyclone cluster assignment is added to the INFO section for each variant in the original somatic VCF file.
    # The cluster is added as ";AFCLU=0".

    # Because superseeker currently calculates the allele frequency of the cluster by itself (averaging over all cluster members, 
    # using the RO and AO values in the VCF). I would try not to use copy number variation region variants. 
    # Although pyclone can "correct" for it, it is still a bit ambiguous because it's hard to say, out of e.g. 5 total copies, 
    # how many had the reference and how many had the variant

    ## Steps ##
    # 1. Read in the VCF and pyclone output for the patient
    # 2. Go through each variant in the pyclone output, and add it to the vcf 
    # Will it be better to go through the VCF and find the pyclone output, or go through pyclone and add to VCF?

    vcf_file_name = vcf_file
    pyclone_file_name = clustered_file

    pyclone_file = open(pyclone_file_name, 'r')
    vcf_file = open(vcf_file_name, 'r')
    output_file = open(output_file_name, "w")

    # For each line in pyclone file, add the {mutation_id : cluster} to a dictionary. 
    cluster_assignments = dict()
    for line in pyclone_file:
        fields = line.strip().split('\t')
        if fields[0] == "mutation_id":
            continue
        elif fields[0] in cluster_assignments.keys():
            continue
        else:
            cluster_assignments[fields[0]] = [fields[2]]

    # for line in vcf file, 
    # if it starts with a "#", just write it directly to the output file
    # Otherwise, find the cluster the variant belongs to, add it to the info field, and then write the line.
    info_line = "##INFO=<ID=AFCLU,Number=1,Type=String,Description=\"The allele frequency cluster this variant belongs to\">"
    for line in vcf_file:
        if line.strip()[:2] == "##":
            output_file.write(line.strip()+'\n')
        elif line.strip()[0] == "#":
            output_file.write(info_line+'\n')
            output_file.write(line.strip()+'\n')
        else:
            fields = line.split('\t')
            mut_id = fields[0]+":"+fields[1]+":"+fields[4]
            if mut_id in cluster_assignments.keys():
                cluster = cluster_assignments[mut_id]
            else:
                print("Skipping " + mut_id + " from original VCF")
                continue 
            skip = False
            info = fields[7].split(";")
            if HIGH_IMPACT:
                for i in reversed(info):
                    if i[:4] == "ANN=":
                        ann = i.split("|")
                        for a in ann:
                            if a == "HIGH" or a == "MODERATE":
                                break
                            if a == "MODIFIER" or a == "LOW":
                                skip = True
                        break
            if not skip:
                info_fields = fields[7].split(";")
                if info_fields[-1].split("=")[0] == "AFCLU":
                    info_fields[-1] = "AFCLU="+cluster[0]
                else:
                    info_fields.append("AFCLU="+cluster[0])
                fields[7] = ";".join(info_fields)
                output_file.write("\t".join(fields))
    output_file.close()
    vcf_file.close()
    pyclone_file.close()


############### Below are functions for identifying patterns of evolution ##################
class subclone:
    def __init__(self, ID):
        self.ID = ID
        self.vafs = []
        self.selection = False
        self.emergence = False
        self.replacement = False

    def add_vaf(self, vaf):
        self.vafs.append(vaf)

def find_replacement(subclones):
    for ID in subclones:
        highest_at_end = True
        not_highest_at_start = False
        if subclones[ID].emergence:
            for check in subclones:
                if check == ID:
                    continue
                if subclones[check].vafs[-1] >= subclones[ID].vafs[-1]:
                    highest_at_end = False
                if subclones[check].vafs[1] > subclones[ID].vafs[1]:
                    not_highest_at_start = True
        if not_highest_at_start and highest_at_end:
            subclones[ID].replacement = True

def find_evolution(subclones): # Making this more invovled could help
    result = []
    for ID in subclones:
        i = 1 # This will start iterating the clusters at the index after germline.
        change = 0.0
        while i < len(subclones[ID].vafs)-1:
            change = change + (float(subclones[ID].vafs[i+1]) - float(subclones[ID].vafs[i]))
            if change >= 0.1:
                subclones[ID].emergence = True
            if change <= -0.1:
                subclones[ID].selection = True
            i = i+1

def identify_evolution(stats_file, output_file):
    input_file_name = stats_file
    input_file = open(input_file_name, 'r')
    subclones = dict()
    clusters = input_file.readline().strip().split("\t")
    too_small = False
    for cluster in clusters:
        subclones[cluster] = subclone(cluster)
    for line in input_file:
        if line[:4] == "Move":
            break
        sample = line.strip().split("\t")
        if len(sample) < 4: # This should exit with no evolution if there is only 1 tumor sample.
            too_small = True
            break
        i = 0
        for val in sample[1:]:
            subclones[str(i)].add_vaf(val)
            i = i+1

    if not too_small:
        find_evolution(subclones)
        find_replacement(subclones)

    selection = False
    emergence = False
    replacement = False

    selection_list = []
    emergence_list = []
    replacement_list = []

    for ID in subclones:
        if subclones[ID].selection:
            selection = True
            selection_list.append(ID)
        if subclones[ID].emergence:
            emergence = True
            emergence_list.append(ID)
        if subclones[ID].replacement:
            replacement = True
            replacement_list.append(ID)
        
    if replacement:
        print("Replacement")
        output_file.write("Replacement\n")
    elif emergence:
        print("Positive Selection")
        output_file.write("Positive Selection\n")
    elif selection:
        print("Negative Selection")
        output_file.write("Negative Selection\n")
    else:
        print("No Evolution")
        output_file.write("No Evolution\n")

    output_file.write("Subclones with Negative Selection: "+",".join(selection_list)+"\n")
    output_file.write("Subclones with Postive Selection: "+",".join(emergence_list)+"\n")
    output_file.write("Subclones with Replacement: "+",".join(replacement_list)+"\n")

    input_file.close()
    output_file.close()

def make_dot_files(subclones_vcf, tmp_graph_files):
    vcf_file_name = subclones_vcf
    tmp_dir = tmp_graph_files

    infile = open(vcf_file_name, "r")
    i = 0
    for line in infile:
        if line[0] != "#":
            break
        if line[0:10] == "##subclone":
            i=i+1
            edges = line.split("\"")[1]
            outfile = open(tmp_dir+"/solution"+str(i)+".gv", "w")
            outfile.write("digraph D{\n")
            for edge in edges.split(", "):
                outfile.write(edge+"\n")
            outfile.write("label=\"Solution "+str(i)+"\"\nlabelloc=\"t\"\n}"+"\n")
            outfile.close()

    infile.close()

def make_line_plot(vcf_file_name, title, show, save, plot_file_name="vaf.png"):

    ## Make a DataFrame of the variant allele frequencies ##    
    variant_df = pd.DataFrame(columns=["Sample", "Position", "BP_Change", "Cluster", "Allele_Frequency"])

    infile = open(vcf_file_name, "r")
    i = 0
    ## Extract info from VCF ##
    for line in infile:
        if line[0] == "#":
            header = line.strip().split()
        else:
            fields = line.strip().split()
            clu = fields[7].split(";")[-1].split("=")[1]
            Format = fields[8].split(":")
            AO_index = Format.index("AO")
            DP_index = Format.index("DP")
            for samp_col_num in range(9, len(fields)):
                sample_name = header[samp_col_num]
                dp = fields[samp_col_num].split(":")[DP_index]
                ao = fields[samp_col_num].split(":")[AO_index]
                af = float(ao)/float(dp)
                variant_df.loc[len(variant_df.index)] = [sample_name, fields[0]+":"+fields[1], fields[3]+"->"+fields[4], clu, af]
    infile.close()

    ## Make the line plot of the VAFs for each mutation in each label group (cell type) ##
    # Get the unique sample order from the original DataFrame
    sample_order = variant_df['Sample'].unique()

    # Pivot the DataFrame to make each column represent a mutation_id and each row a sample_id, with values being cellular_prevalence
    allele_freq_matrix = variant_df.pivot(index='Sample', columns='Position', values='Allele_Frequency')

    # Reindex to preserve the original sample order
    allele_freq_matrix = allele_freq_matrix.reindex(sample_order)

    # Determine unique clusters
    unique_clusters = variant_df['Cluster'].unique()
    cluster_colors = {cluster: plt.get_cmap('gist_rainbow')(i / len(unique_clusters)) for i, cluster in enumerate(unique_clusters)}
    
    # Set font sizes
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 20,
        'axes.labelsize': 18,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14
    })
    plt.figure(figsize=(16, 10))

    # Plotting each mutation_id with a color based on its cluster
    for position in allele_freq_matrix.columns:
        mutation_cluster = variant_df[variant_df['Position'] == position]['Cluster'].iloc[0]
        color = cluster_colors[mutation_cluster]
        plt.plot(allele_freq_matrix.index, allele_freq_matrix[position], label=position, 
                 linestyle='-', linewidth=2, alpha=0.7, color=color)
    
    plt.xticks(rotation=90)
    plt.title(title)
    plt.ylabel('Cellular Prevalence')
    plt.xlabel('Sample ID')
    plt.legend(title='Mutation ID', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    plt.tight_layout()
    
    if show:
        plt.show()
    if save:
        plt.savefig(plot_file_name, bbox_inches='tight')