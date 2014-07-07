from parseMST import Map, write_generic_map
from switch_allele_functions import parse_map_file_MST, convert_genotypes_to_int_array
import numpy


# To summarise and plot a map in R:
  # Create two files, 
    # a map file, giving the map positions
    # a stats file, that give stats for each locus, for quality control purposes

map_chum_08 = Map('chum_08_mst', 'mst', file_path = "/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/mst/chum_08_iter0_mstmap.map")
write_generic_map("/olympus/WORK/WAPLES/Stacks_mapping/Chum_data/psv/chum_08_finalmap.txt", map_chum_08)


#Examine number of crossovers 
# Perhaps read genotype directly from the mstmap input file.
# read mst map output file as above or with parse_map_file.
# kick out LG and loci sitting alone
# How to best count the number of crossovers along a LG???
# need to switch alleles
ini_map, num_loci_on_lg = parse_map_file_MST("Y:\WORK\WAPLES\Stacks_mapping\Chum_data\mst\chum_08_iter0_mstmap.map")
int_arr = convert_genotypes_to_int_array(mappable_08, ini_map)

# This should be move into the psv_working pipeline, so that stats_08 is not needed here (it cannot be pickled).
def write_stats(filename, stats):
    with open(filename, 'w') as OUTFILE:
        for catID, values in stats.items():
            OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                catID, 
                values.epsilon,
                values.best_model,
                values.best_likelihood,
                values.second_best_model,
                values.second_best_likelihood,
                values.x1_seg_pval,
                values.x2_seg_pval
                )
            )
#from psv_working import stats_08
write_stats('/home/ipseg/Data/chum/data/chum_08.stats', stats_08)






# To get count of crossovers, per LG
def parse_mst_map_infile(infile, num_ind):
    genotypes_of_locus = dict()
    with open(infile, 'r') as INFILE:
        found_name_line = False
        for line in INFILE:
            split_line =  line.strip().split()
            if found_name_line:
                if len(split_line) != num_ind + 1: raise Exception
                locus = split_line[0]
                genotypes = split_line[1:num_ind + 1]
                genotypes_of_locus[locus] = genotypes
            else:
                if len(line.strip().split()) == num_ind + 1:
                    individuals = line.strip().split()[1:num_ind + 1]
                    found_name_line = True
    return(genotypes_of_locus)
    
    

  
ini_map, num_loci_on_lg = parse_map_file("Y:\WORK\WAPLES\Stacks_mapping\Chum_data\mst\chum_08_iter0_mstmap.map", 'mst')
genotypes_of_locus = parse_mst_map_infile("Y:\WORK\WAPLES\Stacks_mapping\Chum_data\mst\chum_08_iter0_mstmap.txt", 140)
int_arr = convert_genotypes_to_int_array(genotypes_of_locus, ini_map)
to_be_removed = 0

start_index = 0
end_index = 0

# What do we do if missing genotypes sit next to potential double crossovers?

with open ('Y:\WORK\WAPLES\Stacks_mapping\pest_out.tsv', 'w') as OUTFILE:
    start_index = 0
    end_index = 0
    OUTFILE.write(("{}\t"*13+"\n").format("LG", "length", "end_index", "start_index", "current_LG.shape", "ind_count", "ind_genotypes", "nonzero_genotypes", "count_non_zero", "crossovers", "dbl_crossovers", "flag", "corrected_crossovers"))
    for LG, length in num_loci_on_lg.items():
        end_index += length
        current_LG = int_arr[start_index:end_index]
        t_current_LG = numpy.transpose(current_LG)
        for ind_count, ind_genotypes in enumerate(t_current_LG):
            if length >1:
                other = (ind_genotypes != to_be_removed)
                nonzero_genotypes = ind_genotypes[other]
                crossovers = numpy.sum(numpy.abs(numpy.diff(nonzero_genotypes)))
                dbl_crossovers = numpy.sum(numpy.diff(numpy.where(numpy.abs(numpy.diff(nonzero_genotypes))==1)) == 1)
                trouble_1 = numpy.all(rolling_window(nonzero_genotypes, 5) == [2, 1, 2, 1, 2], axis = 1)
                trouble_2 = numpy.all(rolling_window(nonzero_genotypes, 5) == [1, 2, 1, 2, 1], axis = 1)
                flag = numpy.sum(trouble_1) + numpy.sum(trouble_2) 
                #meaningfull_crossovers = numpy.sum(numpy.abs(numpy.diff(numpy.abs(numpy.diff(nonzero_genotypes)))))/2.
                corrected_crossovers = crossovers - 2 * dbl_crossovers
                if flag:
                    corrected_crossovers = -flag
                OUTFILE.write(("{}\t"*13+"\n").format(LG, length, end_index, start_index, int(current_LG.shape[0]), ind_count, int(ind_genotypes.shape[0]), list(nonzero_genotypes), int(nonzero_genotypes.shape[0]), crossovers, dbl_crossovers, flag, corrected_crossovers))
            start_index = int(end_index)

# subset the part of int_arr that represent a LG
# transpose it
# diff it
# record results

import collections
test_seqs = collections.OrderedDict([
('single'       , numpy.array([1,1,1,1,1,2,2,2,2,2,2])), # one meaningful crossover
('double'       , numpy.array([1,1,1,1,2,1,1,1,1,1,1])), # zero meaningful crossovers
('two_single'   , numpy.array([1,1,1,2,2,2,2,2,1,1,1])), # two meaningful crossovers
('two_double'   , numpy.array([1,1,1,1,2,1,1,2,1,1,1])), # zero meaningful crossovers
('three_single' , numpy.array([1,1,2,2,1,1,1,1,2,2,2])), # three meaningful crossovers
('triple'       , numpy.array([1,1,1,1,1,2,1,2,2,2,2])), # one meaningful crossover
('two_triple'   , numpy.array([2,2,1,1,1,2,1,2,2,2,2])), # two meaningful crossover
('quad'         , numpy.array([1,1,1,1,2,1,2,1,1,1,1])), # zero (two) meaningful crossovers
('penta'        , numpy.array([1,1,1,1,2,1,2,1,2,2,2])), # one meaningful crossover(?)
('end'          , numpy.array([1,1,1,1,1,1,1,1,1,1,2])), # one meaningful crossover
('end_single'   , numpy.array([1,1,1,1,1,1,2,2,2,2,1])), # two meaningful crossovers 
('end_double'   , numpy.array([1,1,1,1,1,2,1,1,1,1,2])),  # one meaningful crossover
('near_double'  , numpy.array([1,1,1,1,1,2,2,1,1,1,1]))  # one meaningful crossover
])



# NOTE The LG numbers given in log file are +1 with repsect to the .map file!!!
from switch_allele_functions import parse_map_file, convert_genotypes_to_int_array
import numpy
import collections
def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return numpy.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def parse_MSTmap_log(infile, outfile_totals, outfile_locations, number_inds, ini_map, num_loci_on_lg, debug):
    """
    Reads a log file from MSTmap. For each LG reports, the crossovers in each indivudal.
    """
    with open(outfile_totals, 'w') as TOTALS:
        with open(outfile_locations, 'w') as LOCATIONS:
            TOTALS.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('LG', 'loci', 'ind_index', 'crossovers', 'dbl_crossovers', 'flag', 'corrected_crossovers', 'location_of_crossovers', 'location_of_dbl_crossovers', 'ind_array'))
            LOCATIONS.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format('LG', 'loci', 'ind_index', 'crossovers', 'dbl_crossovers', 'flag', 'corrected_crossovers', 'crossover_location', 'cM'))
            with open(infile, 'r') as INFILE:
                line_count = 0
                current_LG_genotypes = list()
                current_number_of_bins = None
                for line in INFILE:
                    line = line.strip()
                    line_count +=1
                    if debug: print("Looking at line: {}".format(line_count))
                    if debug: print(line)
                    if line.split()[0:2] == ['finished', 'the']:
                        log_LG =  int(line.split()[2]) 
                        map_LG = log_LG - 1  # Corrected to match the map output file, numbers are off by one.
                        # this marks the end of data from this LG, process then clear the preceeding group.
                        if debug: print('Working on LG: {}'.format(log_LG))
                        ind_genotypes = zip(*current_LG_genotypes) #transpose
                        if len(current_LG_genotypes) > 1:
                            LG_results = find_crossovers(ind_genotypes, debug = debug)
                            for ind_results in LG_results:
                                loc_cross = [int(ii) for ii in ind_results.location_of_crossovers]
                                loc_dbl_cross = [int(ii) for ii in ind_results.location_of_dbl_crossovers]
                                true_crossovers = loc_cross[:]
                                for dbl_cross in loc_dbl_cross:
                                    print(true_crossovers, dbl_cross, loc_dbl_cross)
                                    try:
                                        true_crossovers.remove(dbl_cross)
                                    except ValueError:
                                        if dbl_cross in loc_cross:
                                            pass
                                        else:
                                            raise
                                    try:    
                                        true_crossovers.remove(dbl_cross+1)
                                    except ValueError:
                                        if dbl_cross+1 in loc_cross:
                                            pass
                                        else:
                                            raise                                        
                                TOTALS.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(map_LG, current_number_of_bins, ind_results.ind_index, ind_results.crossovers, ind_results.dbl_crossovers, ind_results.flag, ind_results.corrected_crossovers, '_'.join([str(ii) for ii in loc_cross]), '_'.join([str(ii) for ii in loc_dbl_cross]), numpy.array_str(ind_results.ind_array, max_line_width = 1000)))
                                for crossover in true_crossovers:
                                    cM_location_of_crossover = cM_of_crossover(position = crossover, LG = map_LG, ini_map = ini_map, num_loci_on_lg = num_loci_on_lg)
                                    LOCATIONS.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(map_LG, current_number_of_bins, ind_results.ind_index, ind_results.crossovers, ind_results.dbl_crossovers, ind_results.flag, ind_results.corrected_crossovers, crossover, cM_location_of_crossover))
                        current_LG_genotypes = list()
                    elif (len(line) == number_inds) and len(line.translate(None, '.#')) == 0: # If length is number of inds and one "." and "#" removed its len zero.
                        locus_genotype = [jj for jj in line] # split each chr
                        current_LG_genotypes.append(locus_genotype)
                    elif line[0:15] == 'number of bins:':
                        current_number_of_bins = int(line[15:])

                
            
def find_crossovers(ind_genotypes, debug):
    to_return = list()
    Crossover_Results = collections.namedtuple(typename = 'Crossover_Results', field_names = ['ind_index', 'crossovers', 'dbl_crossovers', 'flag', 'corrected_crossovers', 'location_of_crossovers', 'location_of_dbl_crossovers', 'ind_array'])   
    geno_arr = numpy.array(ind_genotypes)
    my_a = numpy.where(geno_arr == '.')
    my_b = numpy.where(geno_arr == '#')
    int_arr = numpy.empty(shape = (geno_arr.shape), dtype = numpy.int)
    int_arr[my_a] = 1 # allele 1
    int_arr[my_b] = 2 # allele 2
    for ind_index, ind_array in enumerate(int_arr):
        # Count of raw crossovers
        crossovers = numpy.sum(numpy.abs(numpy.diff(ind_array)))
        location_of_crossovers = numpy.where(numpy.abs(numpy.diff(ind_array)) == 1)[0]
        # Count of pairs of crossovers next to each other. fyi, 3 in a row = 2 pairs
        dbl_crossovers = numpy.sum(numpy.diff(numpy.where(numpy.abs(numpy.diff(ind_array))==1)) == 1)
        location_of_dbl_crossovers = numpy.where(numpy.abs(numpy.diff(ind_array, 2))==2)[0]
        # look for three crossovers in a row, flag if seen
        flag = 0
        trouble_1 = numpy.all(rolling_window(ind_array, 4) == [2, 1, 2, 1], axis = 1)
        trouble_2 = numpy.all(rolling_window(ind_array, 4) == [1, 2, 1, 2], axis = 1)
        if numpy.sum(trouble_1) + numpy.sum(trouble_2) > 0: flag += 1 
        # look for crossovers separated by 2 postitions
        trouble_3 = numpy.all(rolling_window(ind_array, 6) == [1, 1, 2, 2, 1, 1], axis = 1)
        trouble_4 = numpy.all(rolling_window(ind_array, 6) == [2, 2, 1, 1, 2, 2], axis = 1)
        if numpy.sum(trouble_3) + numpy.sum(trouble_4) > 0: flag += 2
        corrected_crossovers = crossovers - 2 * dbl_crossovers
        if flag: # If flagged, return a negative corrected crossovers
            corrected_crossovers *= -1
        if debug: print(ind_index, crossovers, dbl_crossovers, flag, corrected_crossovers)
        results = Crossover_Results(ind_index, crossovers, dbl_crossovers, flag, corrected_crossovers, location_of_crossovers, location_of_dbl_crossovers, ind_array)
        to_return.append(results)
    return(to_return)

def cM_of_crossover(position, LG, ini_map, num_loci_on_lg):
    if num_loci_on_lg[str(LG)] < 2:
        raise ValueError('LG only len 1') 
    else:
        LG_map = [(x,y,z) for x,y,z in ini_map if LG == int(y)]
        cM_1 = float(LG_map[position][2])
        cM_2 = float(LG_map[position+1][2])
        average_cM = (cM_1 + cM_2)/2.
        return(average_cM)
                          
ini_map, num_loci_on_lg = parse_map_file("Y:\WORK\WAPLES\Stacks_mapping\Chum_data\mst\chum_08_iter0_mstmap.map", 'mst')            

mymy= parse_MSTmap_log(infile = "Y:\WORK\WAPLES\Stacks_mapping\Chum_data\mst\chum_08_iter0_mstmap.log", 
                outfile_totals =  "Y:/WORK/WAPLES/Stacks_mapping/Chum_data/crossover_totals.tsv",
                outfile_locations = "Y:/WORK/WAPLES/Stacks_mapping/Chum_data/crossover_locations.tsv",
                number_inds = 140,
                ini_map = ini_map,
                num_loci_on_lg = num_loci_on_lg,
                debug = True)
 
        
# I want to get the specific cM location for each crossover.  
# It will be the average of the closest two locus locations
        

    






# END 
with open("Y:\WORK\WAPLES\Stacks_mapping\crossover_log.data", 'r') as INFILE:
    locus_genotypes = list()
    for line in INFILE:
        locus_genotype = [jj for jj in line.strip()]    
        locus_genotypes.append(locus_genotype)
    ind_genotypes = zip(*locus_genotypes)
    geno_arr = numpy.array(ind_genotypes)
    my_a = numpy.where(geno_arr == '.')
    my_b = numpy.where(geno_arr == '#')
    int_arr = numpy.empty(shape = (geno_arr.shape), dtype = numpy.int)
    int_arr[my_a] = 1 # allele 1
    int_arr[my_b] = 2 # allele 2
    for ind_array in int_arr:
        crossovers = numpy.sum(numpy.abs(numpy.diff(ind_array)))
        dbl_crossovers = numpy.sum(numpy.diff(numpy.where(numpy.abs(numpy.diff(ind_array))==1)) == 1)
        trouble_1 = numpy.all(rolling_window(ind_array, 5) == [2, 1, 2, 1, 2], axis = 1)
        trouble_2 = numpy.all(rolling_window(ind_array, 5) == [1, 2, 1, 2, 1], axis = 1)
        flag = numpy.sum(trouble_1) + numpy.sum(trouble_2) 
        corrected_crossovers = crossovers - 2 * dbl_crossovers
        if flag:
            corrected_crossovers = -flag
        print(crossovers, dbl_crossovers, flag, corrected_crossovers)

    

    







# getting raw number of crossovers
# how to exclude adjacent double crossovers
tt = numpy.transpose(int_arr)
to_be_removed = 0
other = (tt != to_be_removed)
tt = tt[other]
aa = numpy.sum(numpy.abs(numpy.diff(tt)))
bb = numpy.sum(numpy.abs(numpy.diff(tt, n=2)) == 2)
aa - 2*bb
