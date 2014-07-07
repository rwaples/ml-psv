from __future__ import division
import collections
import scipy.stats

from CatID_Genotypes_class import CatID_Cross_Genotypes
from CatID_class import CatID
from switch_allele_functions import parse_map_file_rqtl, parse_map_file_MST

class Stacks_SQL:
	"""
	Class to represent a stacks haplotypes file as created by export_sql.pl
	Expects doubled haploid data, with one female parent per cross.
           	path = path to file
           	name = a name for this dataset
           	num_parents = count of parents, they must be the 'leftmost' individuals
           	num_offspring = how many total offspring
           	parent_offspring = a list of ints; 
                    len = num parents with the number of offspring per parent.
	"""
	def __init__(self, path, name, num_parents, num_offspring, parent_offspring):	
		self.name = name
		self.path = path
		self.index_of_parent 	 = dict()
		self.index_of_offspring  = dict()
		self.parent_of_offspring = dict()
		self.offspring_of_parent = dict()
		self.all_parents         = list()
		self.all_offspring       = list()
		self.catIDs = dict()
		self.cross_genotypes = collections.defaultdict(dict)

		with open(path, 'r') as INFILE:
			line_count = 0
			for line in INFILE:
				line_count += 1
				line_split = line.rstrip("\r\n").split('\t')
				if line_count == 1:
					column_headings = line_split
					#check column headings are as expected
					expected_heading = '# Catalog ID\tAnnotation\tChr\tBP\tConsensus Sequence\tNum Parents\tNum Progeny\tNum SNPs\tSNPs\tNum Alleles\tAlleles\tDeleveraged'
					if column_headings[0:12] != expected_heading.split("\t"):
						 print (column_headings[0:12])
						 print (expected_heading.split("\t"))
						 raise StandardError("Unexpected Column Headings, please check input file")   
					self.column_headings = column_headings
					self.all_parents = column_headings[12:(12 + num_parents)]
					self.all_offspring = column_headings[(12 + num_parents):(12 + num_parents + num_offspring)]
									
					for xx in range(len(parent_offspring)):
					   self.offspring_of_parent[self.all_parents[xx]] = list()
					   for ii in range(parent_offspring[xx]):
					       self.offspring_of_parent[self.all_parents[xx]].append(self.all_offspring.pop(0))
					for parent in self.all_parents:
					   self.index_of_parent[parent] = column_headings.index(parent)
					   for offspring in self.offspring_of_parent[parent]:
						self.index_of_offspring[offspring] = column_headings.index(offspring)
						self.parent_of_offspring[offspring] = parent
				
				elif line_count > 1:
					if len(line_split) <= 2: # don't read in the final summary portion
					   break
					else:
					   catID_str = line_split[0]
					   consensus = line_split[4]
					   num_snps, snps, num_alleles, alleles, deleveraged = line_split[7:12]				
					   current_catID = CatID(ID = catID_str, 
                              						consensus = consensus, 
                              						num_snps = num_snps, 
                              						snps = snps,
                              						num_alleles = num_alleles, 
                              						alleles = alleles, 
                              						deleveraged = deleveraged
						                      )
                                        self.add_catID(catID_str, current_catID)
				        for parent_str in self.all_parents:
				            #zip parents genotypes
				            parent_genotpye = line_split[self.index_of_parent[parent_str]]
				            zipped_parent_genotypes = [(parent_str, parent_genotpye)]
				            
				            offspring_genotypes = list()
				            for offspring_name in self.offspring_of_parent[parent_str]:
				                # zip offspring genotype
				                offspring_genotype = line_split[self.index_of_offspring[offspring_name]]
				                offspring_genotypes.append(offspring_genotype)
				            zipped_offspring_genotypes = zip(self.offspring_of_parent[parent_str], offspring_genotypes)
				            cross_genotypes = CatID_Cross_Genotypes(current_catID, zipped_offspring_genotypes, zipped_parent_genotypes)
                                            self.add_catID_cross_genotypes(current_catID, parent_str, cross_genotypes)
                                            
	def add_catID(self, catID_str, CatID):
	   """Adds info about the catID.  Info is from the entire file, rather than a single family"""
	   self.catIDs[catID_str] = CatID
	   return(True)
	
	def add_catID_cross_genotypes(self, CatID, parent, cross_genotypes):
	    catID_str = CatID.ID
	    self.cross_genotypes[catID_str][parent] = cross_genotypes
	    return(True)
	
	def purge_catID(self, catID_str):
	   """Removes the catId from all families"""
	   if catID_str in self.catIDs:
	       del self.catIDs[catID_str]
	   if catID_str in self.cross_genotypes:
	       del self.cross_genotypes[catID_str]
	   return(True)
	   
	def purge_absent_catIDs(self):
            purged = list()
            for catID in self.cross_genotypes.keys():
                if self.cross_genotypes[catID] == {}:
                    self.purge_catID(catID)
                    purged.append(catID)
            print("Purged {} CatIDs removed from all families".format(len(purged)))
            return(purged)
	
	def remove_catID_in_cross(self, catID_str, parent):
	    if catID_str in self.cross_genotypes:
	        if parent in self.cross_genotypes[catID_str]:
	            del self.cross_genotypes[catID_str][parent]
            return(True)
	   
	def add_external_genotypes(self, filepath, parent, missing_codes = ["", "-", "No Call", "Invalid"], taqMan_translate = dict({'XX': 'X/X', 'YY': 'Y/Y', 'XY': 'X/Y', 'YX': 'X/Y'})):
		"""Adds additional genotypes from taqMan assays.  Assumes a .tsv format.
		A CatId object is created with some resonable default values for this locus.
	        parent = a single parent  
	        If encountered, the genotypes on that line will be assocaited with that parent.
	        Genotype on other lines will treated as offspring.         
		Will overwrite any existing genotypes."""
		# Construct translation dict
		for jj in missing_codes:
		  taqMan_translate[jj] = "-"
		parent_genotypes_dict = dict()
		with open(filepath, 'r') as INFILE:
		  table = [line.rstrip("\r\n").split('\t') for line in INFILE]
		  #kick out parents here, to process separately
		  for ii in range(len(table)):
                    if table[ii][0] == parent:
		      parent_line = table.pop(ii)
		      locus_names_line = table[0]
		      parent_genotypes = zip(locus_names_line, parent_line)[1:]
		      for l, g in parent_genotypes:
		          parent_genotypes_dict[l] = g
		t_table = zip(*table)
		offspring = t_table[0][1:]
		for locus_line in t_table[1:]:
		  locus = locus_line[0]
		  current_catID = CatID(locus, 
				consensus = 'NA', 
				num_snps = 1, 
				snps = 'NA',
				num_alleles = 2, 
				alleles = "X;Y", 
				deleveraged = 0
  				)
  		  raw_genotypes = locus_line[1:]
		  genotypes = [taqMan_translate[jj] for jj in raw_genotypes]
		  zipped_offspring_genotypes = zip(offspring, genotypes)
		  zipped_parent_genotypes = parent_genotypes_dict.get(locus, [(parent, ("-"))])
		  cross_genotypes = CatID_Cross_Genotypes(current_catID, zipped_offspring_genotypes = zipped_offspring_genotypes, zipped_parent_genotypes = zipped_parent_genotypes , missing_codes = missing_codes)
		  self.add_catID(current_catID.ID, current_catID)
		  self.add_catID_cross_genotypes(current_catID, parent = parent, cross_genotypes = cross_genotypes)
		return(True)
		
        def remove_missing_within_cross(self, parent, max_miss):
            catIDs_removed = list()
            for catID in self.cross_genotypes.keys():
                if parent in self.cross_genotypes[catID]:
                    current_genotypes = self.cross_genotypes[catID][parent]
                    # look for missing data
                    missing_count = current_genotypes.genotype_counts[frozenset('-')]
                    overall_count = sum(current_genotypes.genotype_counts.values())
                    if float(missing_count)/overall_count > max_miss:
                        del self.cross_genotypes[catID][parent]
                        catIDs_removed.append((current_genotypes.catID.ID, float(missing_count)/overall_count))
            print('Removed {} catIDs from parent: {} for missingness'.format(len(catIDs_removed), parent))
            return(catIDs_removed) 
                
        def remove_max_alleles_within_cross(self, parent, max_alleles = 4):
            catIDs_removed = list()
            for catID in self.cross_genotypes.keys():
                if parent in self.cross_genotypes[catID]:
                    current_genotypes = self.cross_genotypes[catID][parent]
                    # Count number of alleles in this cross
                    num_alleles = len(current_genotypes.allele_counts)
                    if num_alleles > max_alleles:
                        del self.cross_genotypes[catID][parent]
                        catIDs_removed.append((current_genotypes.catID.ID, num_alleles))
            print('Removed {} catIDs from parent: {} for too many alleles'.format(len(catIDs_removed), parent))
            return(catIDs_removed)
                  
        def eval_catID_at_model(self, catID, parent, model, epsilon_bound_high, epsilon_bound_low):
            catID_cross_genotypes = self.cross_genotypes[catID][parent]
            mappings = list()
            epsilons = list()
            likelihoods = list()
            duplication_levels = list()
            for each_trans in catID_cross_genotypes.all_translations:
                translation_dict = each_trans[0]
                mappings.append(translation_dict)
                translated_genotypes = each_trans[1]
                ml_epsilon = model.estimate_epsilon(translated_genotypes)
                epsilons.append(ml_epsilon)
                duplication_levels.append(model.duplication_level)
                if ml_epsilon > epsilon_bound_high:
                    bounded_epsilon = epsilon_bound_high
                elif ml_epsilon < epsilon_bound_low:
                    bounded_epsilon = epsilon_bound_low
                else:
                    bounded_epsilon = ml_epsilon
                adj_prob = model.adjusted_probabilites(bounded_epsilon, max(len(model.alleles), len(catID_cross_genotypes.allele_counts)))
                likelihood = model.mult_likelihood(adj_prob, translated_genotypes)
                likelihoods.append(likelihood)
            results = zip(likelihoods, epsilons, mappings, duplication_levels)
            # sort by most likely
            sorted_results = sorted(results, key = lambda x: x[0], reverse = True)
            best_result = sorted_results[0]
            return(best_result)
            
        def eval_catID(self, catID, parent, models, epsilon_bound_high, epsilon_bound_low):
            results_dict = dict()
            for model_name in models.keys():
                model_to_test = models[model_name]
                result = self.eval_catID_at_model(catID = catID, parent = parent, model = model_to_test, epsilon_bound_high = epsilon_bound_high, epsilon_bound_low = epsilon_bound_low)
                results_dict[model_name] = result
            return(results_dict)
            
        def set_model_results(self, to_test, epsilon_bound_high = 0.1, epsilon_bound_low = 0.01):
            print("Found {} models to test".format(len(to_test)))
            tested_CatID_within_cross = 0
            for catID in self.cross_genotypes.keys():
                for parent in self.cross_genotypes[catID]:
                    num_alleles = len(self.cross_genotypes[catID][parent].alleles)
                    if num_alleles <= 4:
                        model_results = self.eval_catID(catID = catID, parent = parent, models = to_test, epsilon_bound_high = epsilon_bound_high, epsilon_bound_low = epsilon_bound_low)
                        self.cross_genotypes[catID][parent].set_model_results(model_results)
                        tested_CatID_within_cross += 1
                        if tested_CatID_within_cross % 1000 == 0:
                            print("{} CatIDs tested within crosses".format(tested_CatID_within_cross))
                    else:
                        self.cross_genotypes[catID][parent].set_model_results(None)

        def write_mstmap(self, parent, outfile, mappable_loci, 
            distance_function = 'kosambi', cut_off_p_value = .0000001, 
            no_map_dist = 15.0, no_map_size = 2, missing_threshold =.25, 
            estimation_before_clustering = 'yes', detect_bad_data = 'yes', 
            objective_function = 'COUNT', inc_header = 'yes'
            ):
            with open (outfile, 'w') as OUTFILE:
		header = ["population_type DH\n",
			 "population_name {}\n".format(parent), 
			 "distance_function {}\n".format(distance_function),
			 "cut_off_p_value {}\n".format(cut_off_p_value), 
			 "no_map_dist {}\n".format(no_map_dist), 
			 "no_map_size {}\n".format(no_map_size), 
			 "missing_threshold {}\n".format(missing_threshold),
			 "estimation_before_clustering {}\n".format(estimation_before_clustering),
			 "detect_bad_data {}\n".format(detect_bad_data), 
			 "objective_function {}\n".format(objective_function),
			 "number_of_loci {}\n".format(len(mappable_loci)),
			 "number_of_individual {}\n".format(len(self.offspring_of_parent[parent]))]
		if inc_header.lower() == 'yes':
	           OUTFILE.writelines(header)
	           OUTFILE.write("\n")
	        #offspring line   
	        OUTFILE.write("\t".join(['locus_name'] + self.offspring_of_parent[parent] + ["\n"]))
	        # mappped_1 
	        for catID, genos in mappable_loci.items():
	           OUTFILE.write("\t".join([catID] + genos + ["\n"]))
	
	def write_Rqtl(self, parent, outfile, mappable_loci):
	    with open (outfile, 'w') as OUTFILE:
	       # header line
	       header = ['ID'] + [catID for catID in mappable_loci.keys()]
	       OUTFILE.write(",".join(header) + "\n")
	       # LG line
	       LGs = list()
	       for locus in mappable_loci.keys():
	           LGs.append(self.map_of_parent[parent].get(locus, ('999', '1'))[0])
	       OUTFILE.write("," + ",".join(LGs) + "\n")
	       # Ind lines
	       for off_index in range(len(self.offspring_of_parent[parent])):
	           off = self.offspring_of_parent[parent][off_index]
	           OUTFILE.write(off)
	           for jj in mappable_loci.values():
	               OUTFILE.write("," + jj[off_index])
	           OUTFILE.write("\n")

        def set_map(self, parent, mappable_loci, source_file, source_format):
            try:
                self.map_of_parent[parent] = None
            except AttributeError:
                self.map_of_parent = collections.OrderedDict()
            map_position_of_locus = collections.OrderedDict()
            if source_file is None and source_format is None:
                #set all LGs to one, with arbitrary order
                position_counter = 1
                for locus in mappable_loci.keys():
                    map_position_of_locus[locus] = (str(1), str(position_counter))
                    position_counter+=1
            elif source_format.lower() == 'mst':
                ini_map, loci_on_lg = parse_map_file_rqtl(source_file, 'mst')
                for entry in ini_map:
                    locus, lg, pos = entry
                    map_position_of_locus[locus] = (str(lg), str(pos))
            else:
                Exception("Unknown source format")
            self.map_of_parent[parent] = map_position_of_locus

        def write_full_model_likelihood_file(self, parent, outfile, to_test, epsilon_bound_high = 0.1, epsilon_bound_low = 0.01):
            with open(outfile, 'w') as OUTFILE:
                for catID in self.cross_genotypes.keys():
                    if parent in self.cross_genotypes[catID]:
                        raw = self.cross_genotypes[catID][parent].genotype_counts
                        to_write = [(list(k), v) for k, v in raw.items()]
                        OUTFILE.write("{}\t{}".format(catID, to_write))
                        num_alleles = len(self.cross_genotypes[catID][parent].alleles)
                        if num_alleles <= 4:
                            model_results = self.eval_catID(catID = catID, parent = parent, models = to_test, epsilon_bound_high = epsilon_bound_high, epsilon_bound_low = epsilon_bound_low)
                            sorted_model_results = sorted(model_results.items(), key = lambda x: x[1][0], reverse = True)
                            for test in sorted_model_results:
                                model = test[0]
                                likelihood = test[1][0]
                                epsilon = test[1][1]
                                mapping = test[1][2]
                                dup_level = test[1][3]
                                OUTFILE.write("\t{}\t{}\t{}\t{}\t{}".format(model, likelihood, epsilon, mapping, dup_level))
                        OUTFILE.write("\n")

                                
        def find_mappable_markers(self, parent, models_to_test):
            mappable_1 = dict()
            mappable_2 = dict()
            catID_stats = collections.namedtuple('catID_stats',['epsilon', 'best_model',  'best_likelihood', 'second_best_model', 'second_best_likelihood', 'x1_seg_pval', 'x2_seg_pval'])
            stats_of_catID = collections.OrderedDict()
            for catID in self.cross_genotypes.keys():
                if parent in self.cross_genotypes[catID]:
                    cross_genotypes = self.cross_genotypes[catID][parent]
                    if cross_genotypes.model_results is not None:
                        sorted_model_results = sorted(cross_genotypes.model_results.items(), key = lambda x: x[1][0], reverse = True) # sorted by likelihood
                        best_results = sorted_model_results[0]
                        second_best_results = sorted_model_results[1]
                        best_model = best_results[0]
                        second_best_model = second_best_results[0]
                        best_likelihood = best_results[1][0]
                        second_best_likelihood = second_best_results[1][0]
                        epsilon = best_results[1][1]
                        mapping = best_results[1][2]
                        model_class = models_to_test[best_model]
                        translated_genotypes = cross_genotypes.get_translated_offspring_genotypes(mapping)
                        #  map_1
                        calls = list()
                        for off in self.offspring_of_parent[parent]: # fix ordering
                            if model_class.map_1_allele_1 is not None:
                                if translated_genotypes.get(off, '-') in model_class.map_1_allele_1:
                                    calls.append('a')
                                elif translated_genotypes.get(off, '-') in model_class.map_1_allele_2:
                                    calls.append('b')
                                else:
                                    calls.append('-')
                            else:
                                calls.append('-')
                        seg_test_x1 = scipy.stats.binom_test([calls.count('a'),calls.count('b')])
                        if (calls.count('a') + calls.count('b')) == 0:
                            seg_test_x1 = 'NA'
                        mappable_1[catID + "_x1"] = calls
                        
                        #map_2 line
                        calls = list()
                        for off in self.offspring_of_parent[parent]: # fix ordering
                            if model_class.map_2_allele_1 is not None:
                                if translated_genotypes.get(off, '-') in model_class.map_2_allele_1:
                                    calls.append('a')
                                elif translated_genotypes.get(off, '-') in model_class.map_2_allele_2:
                                    calls.append('b')
                                else:
                                    calls.append('-')
                            else:
                                calls.append('-')
                        seg_test_x2 = scipy.stats.binom_test([calls.count('a'), calls.count('b')])
                        if (calls.count('a') + calls.count('b')) == 0:
                            seg_test_x2 = 'NA'
                        mappable_2[catID + "_x2"] = calls
                        
                        #Record Stats
                        stats_of_catID[catID] = catID_stats(epsilon, best_model, best_likelihood, second_best_model, second_best_likelihood, seg_test_x1, seg_test_x2)

            # filter mapped_1, mapped_2
            for k, v in mappable_1.items():
                if v.count('-') == len(v):
                    del mappable_1[k]
            for k, v in mappable_2.items():
                if v.count('-') == len(v):
                    del mappable_2[k]
            for k, v in mappable_2.items():
                mappable_1[k] = v
            mappable = mappable_1        
            return(mappable, stats_of_catID)                     
