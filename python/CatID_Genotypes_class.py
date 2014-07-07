from __future__ import division
import collections
import itertools

from CatID_class import CatID

class CatID_Cross_Genotypes():
    """Class to repsent genotype data at a single catalog ID for a familily of haploid offspring """
    def __init__(self, catID, zipped_offspring_genotypes, zipped_parent_genotypes, alleles = None, missing_codes = ["", "-", "No Call", 'Invalid', 'consensus']):
        if isinstance(catID, CatID):
            self.catID = catID
        else:
            raise StandardError("catID: {} \nShould be a catID".format(catID))
        
        self.missing_codes = missing_codes
        
        if alleles is not None:
            self.alleles = set(alleles)
            self.allowed_alleles = self.alleles.union(set("-"))
        else:
            self.set_alleles(zipped_offspring_genotypes, zipped_parent_genotypes)
             
        self.offspring_genotypes = dict()    
        self.parent_genotypes = dict()            
        for offspring, genotype in zipped_offspring_genotypes:
            if genotype in self.missing_codes:
                genotype = "-"
            genotype_set = frozenset(genotype.split("/"))
            for allele in genotype_set:
                if allele not in self.allowed_alleles:
                    raise StandardError("Allele: [{}] is not expected at catID: {}".format(allele, self.catID))
            self.offspring_genotypes[offspring] = genotype_set
            
        for parent, genotype in zipped_parent_genotypes:
            if genotype in self.missing_codes:
                genotype = "-"
            genotype_set = frozenset(genotype.split("/"))
            for allele in genotype_set:
                if allele not in self.allowed_alleles:
                    raise StandardError("Allele: [{}] is not expected at catID: {}".format(allele, self.catID))
            self.parent_genotypes[parent] = genotype_set
        
        self.set_allele_counts()
        self.set_genotype_counts()
        self.set_translations()
    
    def set_allele_counts(self): 
        # Raw count of alleles within offspring
        self.allele_counts = collections.defaultdict(int)
        for offspring, genotype in self.offspring_genotypes.items():
            for allele in genotype:
                self.allele_counts[allele]+=1
        if "-" in self.allele_counts:
            del self.allele_counts["-"]
        
    def set_genotype_counts(self):
        # Count of genotypes seen       
        self.genotype_counts = collections.defaultdict(int)
        for offspring, genotype in self.offspring_genotypes.items():
            self.genotype_counts[genotype]+=1
        
    # Translations
    def set_translations(self):
        if len(self.alleles) <= 4:
                #raise StandardError("catID: {} has *more* than XX alleles: {}".format(self.catID.ID, self.alleles))                                                                                                    
            num_alleles = len(self.allele_counts)
            old_alleles = self.allele_counts.keys()
            new_alleles = 'ABCD'[:num_alleles]
            translations =  [zip(x, new_alleles) for x in itertools.permutations(old_alleles)]
            translation_dicts = [dict(ii) for ii in translations]
            all_translations = list()
            for trans_dict in translation_dicts:
                trans_dict["-"] = "-"
                single_translation = list()
                for genotype in self.offspring_genotypes.values():
                    genotype_translation = set()
                    for allele in genotype:
                        genotype_translation.add(trans_dict[allele])
                    genotype_translation.difference_update(set(['-']))
                    if len(genotype_translation) > 0:
                        single_translation.append(frozenset(genotype_translation))
                all_translations.append(single_translation)
            #list of frozensets
            self.all_translations = zip(translation_dicts, all_translations)      
            
    def get_translated_offspring_genotypes(self, mapping):
        translated_genotypes = dict()
        for off, geno in self.offspring_genotypes.items():
            trans_geno = set()
            for allele in geno:
                trans_allele = mapping[allele]
                trans_geno.update(set(trans_allele))
            translated_genotypes[off] = frozenset(trans_geno)
        return(translated_genotypes)
        

    def set_alleles(self, zipped_offspring_genotypes, zipped_parent_genotypes):
        zipped_genotypes = zipped_offspring_genotypes + zipped_parent_genotypes
        seen = set()
        #print('set_alleles')
        #print(zipped_parent_genotypes)
        #print(zipped_offspring_genotypes)
        for ind, genotype in zipped_genotypes:
            if genotype in self.missing_codes:
                genotype = "-"
            genotype_set = frozenset(genotype.split("/"))
            for allele in genotype_set:
                if allele not in self.missing_codes:
                    seen.add(allele)
        self.alleles = seen
        self.allowed_alleles = self.alleles.union(set("-"))
        
    def set_model_results(self, model_results):
        self.model_results = model_results
        
                   