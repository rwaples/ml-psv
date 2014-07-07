from Locus_Model_class import Locus_Model

models_to_test = dict()
# One Allele
models_to_test['AA_xx'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = None, loc_3 = None, loc_4 = None,
     map_1_allele_1 = None,
     map_1_allele_2 = None, 
     map_2_allele_1 = None, 
     map_2_allele_2 = None
    )


# Basic Two Allele
models_to_test['AB'] = Locus_Model(loc_1 = ["A", "B"], loc_2 = None, loc_3 = None, loc_4 = None,
     map_1_allele_1 = [frozenset(['A'])],
     map_1_allele_2 = [frozenset(['B'])],
     map_2_allele_1 = None, 
     map_2_allele_2 = None
    )


# Duplicated Two Allele
models_to_test['AA_AB'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["A", "B"], loc_3 = None, loc_4 = None,
     map_1_allele_1 = [frozenset(['A', 'A'])],
     map_1_allele_2 = [frozenset(['A', 'B'])],
     map_2_allele_1 = None, 
     map_2_allele_2 = None
    )
    
models_to_test['AA_BB'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["B", "B"], loc_3 = None, loc_4 = None)
models_to_test['AB_AB'] = Locus_Model(loc_1 = ["A", "B"], loc_2 = ["A", "B"], loc_3 = None, loc_4 = None)

models_to_test['AA_AB_AB'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["A", "B"], loc_3 = ["A", "B"], loc_4 = None)
models_to_test['AA_BB_AB'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["B", "B"], loc_3 = ["A", "B"], loc_4 = None)
models_to_test['AB_AB_AB'] = Locus_Model(loc_1 = ["A", "B"], loc_2 = ["A", "B"], loc_3 = ["A", "B"], loc_4 = None)

#models_to_test['AA_BB_AB_AB'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["B", "B"], loc_3 = ["A", "B"], loc_4 = ["A", "B"])
#models_to_test['AB_AB_AB_AB'] = Locus_Model(loc_1 = ["A", "B"], loc_2 = ["A", "B"], loc_3 = ["A", "B"], loc_4 = ["A", "B"])

# Three Allele
models_to_test['AA_BC'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["B", "C"], loc_3 = None, loc_4 = None,
     map_1_allele_1 = [frozenset(['A', 'B'])],
     map_1_allele_2 = [frozenset(['A', 'C'])],
     map_2_allele_1 = None, 
     map_2_allele_2 = None
    )

models_to_test['AB_AC'] = Locus_Model(loc_1 = ["A", "B"], loc_2 = ["A", "C"], loc_3 = None, loc_4 = None,
     map_1_allele_1 = [frozenset(['A', 'A']), frozenset(['A', 'C'])],
     map_1_allele_2 = [frozenset(['B', 'A']), frozenset(['B', 'C'])],
     map_2_allele_1 = [frozenset(['A', 'A']), frozenset(['B', 'A'])],
     map_2_allele_2 = [frozenset(['A', 'C']), frozenset(['B', 'C'])]
    )
    
#Four Allele
models_to_test['AB_CD'] = Locus_Model(loc_1 = ["A", "B"], loc_2 = ["C", "D"], loc_3 = None, loc_4 = None,
     map_1_allele_1 = [frozenset(['A', 'C']), frozenset(['A', 'D'])],
     map_1_allele_2 = [frozenset(['B', 'C']), frozenset(['B', 'D'])],
     map_2_allele_1 = [frozenset(['A', 'C']), frozenset(['B', 'C'])],
     map_2_allele_2 = [frozenset(['A', 'D']), frozenset(['B', 'D'])]
    )

 

models_to_test['AA_BB_AC'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["B", "B"], loc_3 = ["A", "C"], loc_4 = None)
models_to_test['AA_AB_AC'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["A", "B"], loc_3 = ["A", "C"], loc_4 = None)
models_to_test['AA_BB_CC'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["B", "B"], loc_3 = ["C", "C"], loc_4 = None)

#models_to_test['AA_AB_AB_AC'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["A", "B"], loc_3 = ["A", "B"], loc_4 = ["A", "C"])
#models_to_test['AB_AB_AB_AC'] = Locus_Model(loc_1 = ["A", "B"], loc_2 = ["A", "B"], loc_3 = ["A", "B"], loc_4 = ["A", "C"])
#models_to_test['AA_AB_AB_CC'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["A", "B"], loc_3 = ["A", "B"], loc_4 = ["C", "C"])
#models_to_test['AB_AB_AC_AC'] = Locus_Model(loc_1 = ["A", "B"], loc_2 = ["A", "B"], loc_3 = ["A", "C"], loc_4 = ["A", "C"])
#models_to_test['AA_BB_AC_BC'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["B", "B"], loc_3 = ["A", "C"], loc_4 = ["B", "C"])
#models_to_test['AA_BB_CC_AB'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["B", "B"], loc_3 = ["C", "C"], loc_4 = ["A", "B"])

#models_to_test['AA_BC_CD'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["B", "C"], loc_3 = ["C", "D"], loc_4 = None)
#models_to_test['AA_BB_CD'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["B", "B"], loc_3 = ["C", "D"], loc_4 = None)
#models_to_test['AB_AC_AD'] = Locus_Model(loc_1 = ["A", "B"], loc_2 = ["A", "C"], loc_3 = ["A", "D"], loc_4 = None)
#models_to_test['AA_AB_BC_CD'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["A", "B"], loc_3 = ["B", "C"], loc_4 = ["C", "D"])
#models_to_test['AA_AB_AC_AD'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["A", "B"], loc_3 = ["A", "C"], loc_4 = ["A", "D"])
#models_to_test['AA_BB_BC_AD'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["B", "B"], loc_3 = ["B", "C"], loc_4 = ["A", "D"])
#models_to_test['AA_BB_CC_AB'] = Locus_Model(loc_1 = ["A", "A"], loc_2 = ["B", "B"], loc_3 = ["C", "C"], loc_4 = ["A", "B"])

# Traps
models_to_test['singles'] = Locus_Model(loc_1 = ["A", "B", "C", "D"], loc_2 = None, loc_3 = None, loc_4 = None)
models_to_test['doubles'] = Locus_Model(loc_1 = ["A", "B", "C", "D"], loc_2 = ["A", "B", "C", "D"], loc_3 = None, loc_4 = None)
models_to_test['triples'] = Locus_Model(loc_1 = ["A", "B", "C", "D"], loc_2 = ["A", "B", "C", "D"], loc_3 = ["A", "B", "C", "D"], loc_4 = None)
models_to_test['quads'] = Locus_Model(loc_1 = ["A", "B", "C", "D"], loc_2 = ["A", "B", "C", "D"], loc_3 = ["A", "B", "C", "D"], loc_4 = ["A", "B", "C", "D"])

def write_models(outfile):
    with open(outfile, 'w') as OUTFILE:
        for name, model in models_to_test.items():
            OUTFILE.write( "{}\t{}\t{}\t{}\t{}\t{}\n".format(name, model.loc_1, model.loc_2, model.loc_3, model.loc_4, model.prob_of_genotype))