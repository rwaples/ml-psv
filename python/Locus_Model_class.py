#Likelihoods based on segregation patterns
from __future__ import division
import collections
import itertools

from Numpy_Multinomial_class import Numpy_Multinomial, log_fac_lookup
from epsilon import ml_epsilon

                     
class Locus_Model(object):
    """
    self.alleles = set of alleles
    self.duplication_level = number of duplicated positions
    self.prob_of_genotype[frozenset(["A", "B"])] = the probablity of the [A, B] genotype
    map_x_alle_x:
    Holds information on how to encode one or two mappable markers matching this model.
    Specified for each initialized model.  
    If there is one mappable position, assign a list of genotypes to each loc_1 allele, loc_2 alleles, etc. should be None
    if there are two mappable positions, assign a list of genotypes to each loc_1 & loc_2 allele
    If there are zero mappable positions, assign None to each loc_1 and loc_2 allele.
    map_1_allele_1 is a list of frozenset genotypes that all code for the first allele in map position one
    """
    def __init__(self, loc_1, loc_2 = None, loc_3 = None, loc_4 = None, 
        map_1_allele_1 = None, map_1_allele_2 = None, map_2_allele_1 = None, 
        map_2_allele_2 = None
        ):
        
        self.loc_1 = loc_1
        self.loc_2 = loc_2
        self.loc_3 = loc_3
        self.loc_4 = loc_4
        self.map_1_allele_1 = map_1_allele_1
        self.map_1_allele_2 = map_1_allele_2
        self.map_2_allele_1 = map_2_allele_1
        self.map_2_allele_2 = map_2_allele_2
        

        if self.loc_4 is None:
            if self.loc_3 is None:
                if self.loc_2 is None:
                    #only loc_1 is present
                    self.alleles = set(self.loc_1)
                    self.duplication_level = 0
                else:
                    #loc_1 & loc_2
                    self.alleles = set(self.loc_1 + self.loc_2)
                    self.duplication_level = 1
            else:
                #loc_1 & loc_2 & loc_3
                self.alleles = set(self.loc_1 + self.loc_2 + self.loc_3)
                self.duplication_level = 2
        else:
                self.alleles = set(self.loc_1 + self.loc_2 + self.loc_3 + self.loc_4)
                self.duplication_level = 3
        self.prob_of_genotype = self.genotype_probabilites(loc_1 = self.loc_1, loc_2 = self.loc_2, loc_3 = self.loc_3, loc_4 = self.loc_4)


    def find_error_genotypes_for_x_alleles(self, x):
        """Combines x different alleles in all possible ways. Reports the combinations 
        that do not appear as possibibilities for the base model"""
        error_genotypes = list()
        for jj in self.combinations_of_x_alleles(x):
            if self.prob_of_genotype[jj] == 0:
                error_genotypes.append(jj)
        return (error_genotypes)
       
    def combinations_of_x_alleles(self, x):
        """All possible genotypes with x number of alleles.  
        For x = 0 and x = 1, results for x+1 are reported."""
        if x == 0:
            allele_combinations = list(self.genotype_probabilites(loc_1 = ['A'], loc_2 = None, loc_3 = None, loc_4 = None).keys())
        elif x <= 2:
            allele_combinations = list(self.genotype_probabilites(loc_1 = ['A', 'B'], loc_2 = ['A', 'B'], loc_3 = None, loc_4 = None).keys())
        elif x == 3:
            allele_combinations = list(self.genotype_probabilites(loc_1 = ['A', 'B', 'C'], loc_2 = ['A', 'B', 'C'], loc_3 = ['A', 'B', 'C'], loc_4 = None).keys())
        elif x == 4:
            allele_combinations = list(self.genotype_probabilites(loc_1 = ['A', 'B', 'C', 'D'], loc_2 = ['A', 'B', 'C', 'D'], loc_3 = ['A', 'B', 'C', 'D'], loc_4 = ['A', 'B', 'C', 'D']).keys())
        else:
            raise StandardError("Too many alleles : {}".format(x))
        return(allele_combinations)
    
    def how_to_map(self, map_1_allele_1, map_1_allele_2, map_2_allele_1, map_2_allele_2):
        """Unused"""
        pass

            
    def genotype_probabilites(self, loc_1, loc_2, loc_3, loc_4):
        """
        Used within __int__
        Returns a collections.defaultdict with frozenset genotypes as keys 
        and the expected segregation probability as values.
        At the genotyping stage, we cannot discern dosage or phase so that:
            A/A = A
            A/A/B = A/B
            A/B = B/A
        """
        if loc_4 is None:
            if loc_3 is None:
                if loc_2 is None:
                    #only loc_1 is present
                    inherit = list(itertools.product(loc_1))
                else:
                    #loc_1 & loc_2
                    inherit = list(itertools.product(loc_1, loc_2))
            else:
                #loc_1 & loc_2 & loc_3
                inherit = list(itertools.product(loc_1, loc_2, loc_3))
        else:
                inherit = list(itertools.product(loc_1, loc_2, loc_3, loc_4))
        prob_unit = 1./len(inherit)
        visible = [frozenset(xx) for xx in inherit]
        visible_of_inherit = dict()
        prob_of_visible = collections.defaultdict(float)
        for ii in range(len(inherit)):
            visible_of_inherit[inherit[ii]] = visible[ii]
            prob_of_visible[visible[ii]] += prob_unit
        return(prob_of_visible)
        
    def estimate_epsilon(self, translated_genotypes, DEBUG = False):
        """Returns an estimate of epsilon, the error rate 
        (Epsilon = the rate at which a genotype is replaced by a random genotype).
        Epsilon is estimated by:
            Using the model to determine which combinations of alleles are possible as genotypes.
            It will be possible to classify all observed genotypes that match these as errors.
            The fraction of observed genotypes matching errors is the observable error rate.
            This rate is multiplied by (unique_genotypes/unique_error_genotypes) to account for 
            unobservable errors, ie those resulting in possible genotypes.
        """
        alleles = set()
        alleles.update(*translated_genotypes)
        err_genotypes = self.find_error_genotypes_for_x_alleles(len(alleles))
        count_of_error_genotypes = sum([1 for x in translated_genotypes if x in err_genotypes ])
        ml_epi = ml_epsilon(num_error_catagories = len(err_genotypes),
                    num_catagories = len(self.combinations_of_x_alleles(len(alleles))),
                    count_error_genotypes = count_of_error_genotypes,
                    count_total_genotypes = len(translated_genotypes)
                    )
        if DEBUG:
            print("error genotypes: {}".format(err_genotypes))
            print("alleles = {}".format(alleles))
            print("num_error_catagories = {}".format(len(err_genotypes)))
            print("num_catagories = {}".format(len(self.combinations_of_x_alleles(len(alleles)))))       
            print("count_error_genotypes = {}".format( count_of_error_genotypes))        
            print("count_total_genotypes = {}".format(len(translated_genotypes)))        
        return(ml_epi)
        
    def adjusted_probabilites(self, epsilon, num_alles):
        """num alleles should be the greater of alleles in model and alleles in genotype"""
        allele_combinations = self.combinations_of_x_alleles(num_alles)
        adjusted_probabilites = self.prob_of_genotype.copy()
        for xx in allele_combinations:
            p = self.prob_of_genotype[xx] # base probability
            adj_p = p + (epsilon/len(allele_combinations)) - (p * epsilon)
            adjusted_probabilites[xx] = adj_p
        return(adjusted_probabilites)
        
    def mult_likelihood(self, adj_prob, translated_genotypes, DEBUG = False):
        """Returns the log probability of seeing the genotypes given 
        the adjusted probabilities of observing them.
        """
        probs = list()
        counts = list()
        for ii in adj_prob.keys():
            prob = adj_prob[ii]
            count = translated_genotypes.count(ii)
            if prob == 0 and count > 1:
                raise StandardError("Genotypes Seen For impossible case: {}".format(ii))
            elif prob > 0:    
                probs.append(prob)
                counts.append(count)
        if DEBUG:
            print("Expected Ratios: {}".format(probs))
            print("Observed: {}".format(counts))    
        mult = Numpy_Multinomial(probs)
        log_pmf = mult.log_pmf(counts)
        return(log_pmf)
        