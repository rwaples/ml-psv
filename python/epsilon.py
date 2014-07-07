from __future__ import division


def probability_with_error(p, epsilon, num_catagories):
    """Adjust a genotype probability for a given 
    epsilon value and number of alternate possible categories."""
    if p > 1 or p < 0:
          raise ValueError("{} is not in range [0,1]".format(p))
    prob = p + (epsilon/num_catagories) - (p * epsilon)
    return(prob)

def ml_epsilon(num_error_catagories, num_catagories, count_error_genotypes, count_total_genotypes):
    """
    Returns an estimate of epsilon, the error rate 
    (Epsilon = the rate at which a genotype is replaced by a random genotype).
    Epsilon is estimated by:
        Using the model to determine which combinations of alleles are possible as genotypes.
        It will be possible to classify all observed genotypes that match these as errors.
        The fraction of observed genotypes matching errors is the observable error rate.
        This rate is multiplied by (unique_genotypes/unique_error_genotypes) to account for 
        unobservable errors, ie those resulting in possible genotypes.
    """
    if count_total_genotypes == 0:
        return(None)
    if num_error_catagories == 0:
        return(0)    

    else:
        inverse_fraction_errors_observed = num_catagories / num_error_catagories
        fraction_of_errors = count_error_genotypes/ count_total_genotypes
        ml_of_epsilon = inverse_fraction_errors_observed * fraction_of_errors
        return(ml_of_epsilon)      
