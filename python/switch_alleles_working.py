# cd Z:\HOME\Waples\Stacks_mapping\Python
from __future__ import division
import switch_allele_functions as switchAlleles

import numpy
import itertools
import subprocess
import shlex

import cPickle as pickle
pickle_file = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/pkl/mappable_01.pkl"
INFILE = open(pickle_file, 'rb')
mappable_01 = pickle.load(INFILE)
INFILE.close()

pickle_file = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/pkl/mappable_08.pkl"
INFILE = open(pickle_file, 'rb')
mappable_08 = pickle.load(INFILE)
INFILE.close()

pickle_file = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/pkl/mappable_09.pkl"
INFILE = open(pickle_file, 'rb')
mappable_09 = pickle.load(INFILE)
INFILE.close()

pickle_file = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/pkl/chum_fams.pkl"
INFILE = open(pickle_file, 'rb')
chum_fams = pickle.load(INFILE)
INFILE.close()
##### DONE PICKLE #######

MST_exe = "/home/ipseg/Programs/MSTMap/MSTMap.exe"

# run MST map once

def run_MSTmap(MST_exe, INFILE, OUTFILE, LOGFILE):
    MSTmap_call_string = "{} {} {}".format(MST_exe, INFILE, OUTFILE)
    print(MSTmap_call_string)
    with open(LOGFILE, 'w') as LOG:
        subprocess.call(shlex.split(MSTmap_call_string),shell=False, universal_newlines=True, stdout=LOG)
        
        
# Intiall mappings   

#Chum_01_initial_MST = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_01_mstmap.txt"
#Chum_08_initial_MST = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_08_mstmap.txt"
#Chum_09_initial_MST = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_09_mstmap.txt"
#
#Chum_01_initial_MAP = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_01_mstmap.map"
#Chum_08_initial_MAP = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_08_mstmap.map"
#Chum_09_initial_MAP = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_09_mstmap.map"
#
#Chum_01_initial_LOG = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_01_mstmap.log"
#Chum_08_initial_LOG = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_08_mstmap.log"
#Chum_09_initial_LOG = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_09_mstmap.log"
#
#run_MSTmap(MST_exe = MST_exe, INFILE = Chum_01_initial_MST, OUTFILE = Chum_01_initial_MAP, LOGFILE = Chum_01_initial_LOG)
#run_MSTmap(MST_exe = MST_exe, INFILE = Chum_08_initial_MST, OUTFILE = Chum_08_initial_MAP, LOGFILE = Chum_08_initial_LOG)
#run_MSTmap(MST_exe = MST_exe, INFILE = Chum_09_initial_MST, OUTFILE = Chum_09_initial_MAP, LOGFILE = Chum_09_initial_LOG)


######## Chum_08 #########
#MSTmap_infile = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_08_mstmap.txt"
#initial_linkage_map_file= '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_08_mstmap.map'
#logfile = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_08_mstmap.map.log"
#MSTmap_call_string = "{} {} {}".format(MST_exe, MSTmap_infile, initial_linkage_map_file)
#print(MSTmap_call_string)
#with open(logfile, 'w') as log:
#    subprocess.call(shlex.split(MSTmap_call_string),shell=False, universal_newlines=True, stdout=log)


# Now set it up to run multiple times for each family
#mappable = mappable_08
#reps = 2
#NR_matrix_threshold = 840
#sum_lod_threshold = 15000
#in_linkage_map_file = Chum_08_initial_MAP

def switch_alleles_loop(mappable, family, parent, starting_map, NR_matrix_threshold, reps, MST_grouping_threshold = .0000001, map_type = 'mst', MST_exe = MST_exe):    
    for xx in range(reps):
        ini_map, sq_sum_lod, sq_NR_matrix, R, NR, lgs_longer_than_1 = switchAlleles.calculate_switch_stats(mappable, starting_map, 'mst', MST_grouping_threshold)
        where_pairs = numpy.where(sq_NR_matrix > NR_matrix_threshold)
        LGs_to_switch = switchAlleles.find_LGs_to_switch(lgs_longer_than_1, where_pairs)
        loci_to_switch = switchAlleles.find_loci_to_switch(ini_map, LGs_to_switch)
        switched_alleles_file = '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/switched_alleles_{}_iter{}.txt'.format(family, xx)
        switchAlleles.update_switched_alleles_file(switched_alleles_file, loci_to_switch)
        mappable = switchAlleles.switch_alleles(loci_to_switch, mappable)
        if len(loci_to_switch) == 0:
            print("found zero loci to switch, ending loop")
            break
        print("Starting MSTmap on family: {}, rep: {}".format(family, xx))
        #run MSTmap
        MSTmap_infile = '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_{}_iter{}_mstmap.txt'.format(family, xx)
        MSTmap_logfile = '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_{}_iter{}_mstmap.log'.format(family, xx)
        out_linkage_map_file = '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_{}_iter{}_mstmap.map'.format(family, xx)
        chum_fams.write_mstmap(parent, MSTmap_infile, mappable)
        run_MSTmap(MST_exe = MST_exe, INFILE = MSTmap_infile, OUTFILE = out_linkage_map_file, LOGFILE = MSTmap_logfile)
        print("Finished MSTmap on family: {}, rep: {}".format(family, xx))
        #renew name
        starting_map = out_linkage_map_file


switch_alleles_loop(mappable_08, family = '08', parent = "CMUW10X_0008", starting_map= "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_08_mstmap.map",
    NR_matrix_threshold = 840, reps = 2, MST_grouping_threshold = .000001)
    
switch_alleles_loop(mappable_01, family = '01', parent = "CMUW10X_0001", starting_map= "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_01_mstmap.map",
    NR_matrix_threshold = 800, reps = 3, MST_grouping_threshold = .000001)
    
switch_alleles_loop(mappable_09, family = '09', parent = "CMUW10X_0009", starting_map= "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_09_mstmap.map",
    NR_matrix_threshold = 800, reps = 3, MST_grouping_threshold = .000001)
    
    
#TODO
# Re-pickle "mappables_0X" here so as to maintain the allele switching status
to_switch_01 = list()
with open("/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/switched_alleles_01_iter0.txt") as INFILE:
    for line in INFILE:
        to_switch_01.append(line.strip())
to_switch_08 = list()
with open("/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/switched_alleles_08_iter0.txt") as INFILE:
    for line in INFILE:
        to_switch_08.append(line.strip())
to_switch_09 = list()
with open("/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/switched_alleles_09_iter0.txt") as INFILE:
    for line in INFILE:
        to_switch_09.append(line.strip())


mappable_01 = switchAlleles.switch_alleles(to_switch_01, mappable_01)
mappable_08 = switchAlleles.switch_alleles(to_switch_08, mappable_08)
mappable_09 = switchAlleles.switch_alleles(to_switch_09, mappable_09)



################################


# Calculate some summary stats from the existing recent linkage map
ini_map, sq_sum_lod, sq_NR_matrix, R, NR, lgs_longer_than_1 = switchAlleles.calculate_switch_stats(mappable, in_linkage_map_file, 'mst')
# Determine which pairs of LGs are potential switches
# Perhaps come up with a way to estiamte NR_threshold here.
where_pairs = numpy.where(sq_NR_matrix > NR_matrix_threshold)
# Determines which LGs and loci to switch 
LGs_to_switch = switchAlleles.find_LGs_to_switch(lgs_longer_than_1, where_pairs)
loci_to_switch = switchAlleles.find_loci_to_switch(ini_map, LGs_to_switch)
# updates switched alleles file
switched_alleles_file = '/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/switched_alleles_08_iter{}.txt'.format(xx)
switchAlleles.update_switched_alleles_file(switched_alleles_file, loci_to_switch)
# updates mappapble
mappable = switchAlleles.switch_alleles(loci_to_switch, mappable)
#run MSTmap
MSTmap_infile = '/home/ipseg/Data/chum/data/chum_08_iter{}_mstmap.txt'.format(xx)
MSTmap_logfile = '/home/ipseg/Data/chum/data/chum_08_iter{}_mstmap.log'.format(xx)
out_linkage_map_file = '/home/ipseg/Data/chum/data/chum_08_iter{}_mstmap.map'.format(xx)
chum_fams.write_mstmap('CMUW10X_0008', MSTmap_infile, mappable)
MSTmap_call_string = "{} {} {}".format(MST_exe, MSTmap_infile, out_linkage_map_file)
with open(MSTmap_logfile, 'w') as log:
    subprocess.call(shlex.split(MSTmap_call_string),shell=False, universal_newlines=True, stdout=log)
#renew name
in_linkage_map_file = out_linkage_map_file
    
    
    
    
    
    
    

import switch_allele_functions as switchAlleles
mappable = mappable_08
chum_fams.set_map(parent = 'CMUW10X_0008', mappable_1 = mappable[0], mappable_2 = mappable[1], source_file = '/home/ipseg/Data/chum/data/chum_08_iter0_mstmap.map', source_format = 'mst')
with open('/home/ipseg/Data/chum/data/switched_alleles_iter0.txt') as INFILE:
    loci_to_switch = list()
    for line in INFILE:
        loci_to_switch.append(line.strip())
mappable = switchAlleles.switch_alleles(loci_to_switch, mappable)
chum_fams.write_Rqtl('CMUW10X_0008', '/home/ipseg/Data/chum/data/chum_08_Rqtl_switched.txt', mappable[0], mappable[1])



#NR_matrix = get_NR_matrix(NR)
n = len(mappable[0].items()[0][1]) #number of individuals
P = .0000001
#NR_threshold = get_threshold_recombinants_for_same_LGs(n, P)
numpy.where(sq_NR_matrix > 840)

# form linkage groups (MST)



# write MSTmap file
#from psv_working import chum_fams, mappable_08, mappable_01, mappable_09
# read in mapped loci
# calculate LG pairwise LOD scores
# evaluate if we continue
# pick LGs to switch
# switch alleles in mappable
# write/edit switch alleles file

mappable = mappable_08
linkage_map_file = "Z:/HOME/Waples/Stacks_mapping/Chum_data/rqtl/chum_08_ini_map.tsv"
sum_lod_threshold = 15000
ini_map, sq_sum_lod = switchAlleles.prepare_to_switch(mappable, linkage_map_file)
LGs_to_switch, loci_to_switch = switchAlleles.switch_with_threshold(ini_map, sq_sum_lod, sum_lod_threshold)
switched_mappable_1 = switchAlleles.switch_alleles(loci_to_switch, mappable)

chum_fams.write_Rqtl('CMUW10X_0008', 'Z:/HOME/Waples/Stacks_mapping/Chum_data/psv/chum_08_switched_1_Rqtl.txt', switched_mappable_1[0], switched_mappable_1[1])
chum_fams.write_mstmap('CMUW10X_0008', 'Z:/HOME/Waples/Stacks_mapping/Chum_data/psv/chum_08_switched_1_mstmap.txt', switched_mappable_1[0], switched_mappable_1[1])
# Go to Rqtl, rebuild linkage groups

mappable = switched_mappable_1
linkage_map_file = "Z:/HOME/Waples/Stacks_mapping/Chum_data/rqtl/chum_08_switched_1_map.tsv"
ini_map, sq_sum_lod = switchAlleles.prepare_to_switch(mappable, linkage_map_file)
LGs_to_switch, loci_to_switch = switchAlleles.switch_with_threshold(ini_map, sq_sum_lod, sum_lod_threshold)
switched_mappable_2 = switchAlleles.switch_alleles(loci_to_switch, mappable)

chum_fams.write_Rqtl('CMUW10X_0008', 'Z:/HOME/Waples/Stacks_mapping/Chum_data/psv/chum_08_switched_2_Rqtl.txt', switched_mappable_2[0], switched_mappable_2[1])
chum_fams.write_mstmap('CMUW10X_0008', 'Z:/HOME/Waples/Stacks_mapping/Chum_data/psv/chum_08_switched_2_mstmap.txt', switched_mappable_2[0], switched_mappable_2[1])
# Go to Rqtl, rebuild linkage groups

mappable = switched_mappable_2
linkage_map_file = "Z:/HOME/Waples/Stacks_mapping/Chum_data/rqtl/chum_08_switched_2_map.tsv"
ini_map, sq_sum_lod = switchAlleles.prepare_to_switch(mappable, linkage_map_file)
LGs_to_switch, loci_to_switch = switchAlleles.switch_with_threshold(ini_map, sq_sum_lod, sum_lod_threshold)
switched_mappable_3 = switchAlleles.switch_alleles(loci_to_switch, mappable)


chum_fams.write_Rqtl('CMUW10X_0008', 'Z:/HOME/Waples/Stacks_mapping/Chum_data/psv/chum_08_switched_3_Rqtl.txt', switched_mappable_3[0], switched_mappable_3[1])
chum_fams.write_mstmap('CMUW10X_0008', 'Z:/HOME/Waples/Stacks_mapping/Chum_data/psv/chum_08_switched_3_mstmap.txt', switched_mappable_3[0], switched_mappable_3[1])
# Go to Rqtl, rebuild linkage groups







#Parameters


genotypes_of_locus = switchAlleles.combine_mappable_loci(mappable)

ini_map, loci_on_lg = switchAlleles.parse_map_file(linkage_map_file)

int_arr = switchAlleles.convert_genotypes_to_int_array(genotypes_of_locus, ini_map)

num_loci = int_arr.shape[0]
num_pairs =  int((num_loci * (num_loci-1))/2)
pairs = itertools.combinations(int_arr, 2)

#import timeit
#timeit.timeit('next(switchAlleles.getR(pairs))', setup = "import switch_allele_functions as switchAlleles; from __main__ import int_arr; from __main__ import pairs", number = 100000)

R = numpy.fromiter(switchAlleles.getR(pairs), dtype = numpy.float64, count = num_pairs)

pairs = itertools.combinations(int_arr, 2)
NR = numpy.fromiter(switchAlleles.getNR(pairs), dtype = numpy.float64, count = num_pairs)

ml_R_frac = switchAlleles.get_ml_R_frac(R = R, NR = NR)
Z = switchAlleles.get_LOD(R = R, NR = NR, R_frac = ml_R_frac)

rf = switchAlleles.get_rf_matrix(ml_R_frac)
lod = switchAlleles.get_lod_matrix(Z)

index_of_lg = switchAlleles.get_index_of_LG(loci_on_lg)

lgs_longer_than_1 = switchAlleles.find_LGs_with_multiple_loci(index_of_lg, loci_on_lg)

mean_rf = switchAlleles.get_LG_pairwise_mean_rf(lgs_longer_than_1, rf, index_of_lg)
mean_lod = switchAlleles.get_LG_pairwise_mean_lod(lgs_longer_than_1,lod, index_of_lg)
sum_lod = switchAlleles.get_LG_pairwise_sum_lod(lgs_longer_than_1,lod, index_of_lg)

sq_sum_lod = switchAlleles.get_square_sum_of_lod(sum_lod, lgs_longer_than_1)

where_pairs = switchAlleles.find_LG_pairs_with_sum_lod_above(sq_sum_lod, sum_lod_threshold)

LGs_to_switch = switchAlleles.find_LGs_to_switch(where_pairs)

loci_to_switch = switchAlleles.find_loci_to_switch(ini_map, LGs_to_switch)

switched_mappable = switchAlleles.switch_alleles(loci_to_switch, mappable)





