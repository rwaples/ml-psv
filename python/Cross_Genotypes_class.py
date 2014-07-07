from __future__ import division
from Locus_Genotypes_class import Locus_Genotypes
from models_to_test import models_to_test

class Cross_Genotypes():
    """genotypes_file should be an encoded output file from Stacks_Haplotypes_sql class."""
    def __init__(self, genotypes_file, delimiter = ','):
        """Constructs two dictionaries from the input file.
        One for offspring genotypes, another for parent genotypes."""
        self.offspring_genotypes_of_locus = dict()
        self.parent_genotype_of_locus = dict()
        with open(genotypes_file, 'r') as INFILE:
            self.rows = [line.split(delimiter) for line in INFILE]
            self.cols = zip(*self.rows)
            self.parent = self.cols[0][0]
            self.offspring = self.cols[0][2:]
            for locus_row in self.cols[1:]:
                locus = locus_row[0]
                self.parent_genotype_of_locus[locus] = locus_row[1]
                self.offspring_genotypes_of_locus[locus] = locus_row[2:]
    
    def load_genotypes(self, locus):
        """Creates a Locus_Genotypes object for the offspring genotypes at the target locus"""
        offspring_genotypes = self.offspring_genotypes_of_locus[locus]
        genos = Locus_Genotypes(offspring_genotypes)
        return(genos)
        
    def eval_locus(self, locus, model, epsilon_bound_high = 0.1, epsilon_bound_low = 0.01):
        """Evaluates the all permutations of the target locus with the model specified."""
        genotypes = self.load_genotypes(locus)
        dicts = list()
        genos = list()
        epsilons = list()
        likelihoods = list()
        for trans in genotypes.all_translations:
            #print(model.estimate_epsilon(trans))
            translation_dict = trans[0]
            dicts.append(translation_dict)
            translated_genotypes = trans[1]
            genos.append(translated_genotypes)
            ml_epsilon = model.estimate_epsilon(translated_genotypes)
            epsilons.append(ml_epsilon)
            if ml_epsilon > epsilon_bound_high:
                bounded_epsilon = epsilon_bound_high
            elif ml_epsilon < epsilon_bound_low:
                bounded_epsilon = epsilon_bound_low
            else:
                bounded_epsilon = ml_epsilon
            adj_prob = model.adjusted_probabilites(bounded_epsilon, max(len(model.alleles), len(genotypes.raw_allele_counts)))
            likelihood = model.mult_likelihood(adj_prob, translated_genotypes)
            likelihoods.append(likelihood)
        return(zip(dicts, genos, epsilons, likelihoods))
        
    def get_translation_with_max_likelihood(self, locus, model, epsilon_bound_high, epsilon_bound_low):
        locus_evals = self.eval_locus(locus = locus, model = model, epsilon_bound_high= epsilon_bound_high, epsilon_bound_low = epsilon_bound_low)
        # highest log likelihood
        most_likely = sorted(locus_evals, key = lambda x: x[3], reverse = True)[0] 
        return(most_likely)
        
    def eval_models_for_locus(self, locus, models, epsilon_bound_high, epsilon_bound_low):
        all_results = dict()
        for model_name in models_to_test.keys():
            result = self.get_translation_with_max_likelihood(locus = locus, model = models_to_test[model_name], epsilon_bound_high= epsilon_bound_high, epsilon_bound_low = epsilon_bound_low)
            all_results[model_name] = result
        return(all_results)
        
    def sort_results(self, all_results):
        most_likely = sorted(all_results.items(), key = lambda x: x[1][3], reverse=True)
        return(most_likely)
                                    
    def find_psvs(self, outfile, to_test = models_to_test, epsilon_bound_high = 0.1 , epsilon_bound_low = 0.01):
        """Loops through all loci and calls .eval_locus() on each, the results are written to [outfile] """
        with open(outfile, 'w') as OUTFILE:
            for locus in self.parent_genotype_of_locus.keys():
                counts = Locus_Genotypes(self.offspring_genotypes_of_locus[locus]).genotype_counts
                OUTFILE.write("{}\t{}".format(locus, counts))
                num_alleles = len(Locus_Genotypes(self.offspring_genotypes_of_locus[locus]).raw_allele_counts)
                if num_alleles <= 4:
                    results = self.eval_all_models_for_locus(locus, models = to_test, epsilon_bound_high = epsilon_bound_high, epsilon_bound_low = epsilon_bound_low)
                    sorted_results = self.sort_results(results)
                    for test in sorted_results:
                        model = test[0]
                        liklihood = test[1][3]
                        epsilon = test[1][2]
                        mapping = test[1][0]
                        OUTFILE.write("\t{}\t{}\t{}\t{}".format(model, liklihood, epsilon, mapping)) 
                OUTFILE.write("\n")