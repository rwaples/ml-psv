from __future__ import division
#import math
import numpy

log_fac_lookup = dict()
log_fac_lookup[0] = 0
for xx in range(1, 200):
    n = numpy.arange(1,xx+1)
    log_fac_lookup[xx] = numpy.sum(numpy.log(n))

class Numpy_Multinomial(object):
  def __init__(self, params):
    self._params = numpy.array(params, dtype = numpy.float)

  def pmf(self, counts):
    if not(len(counts)==len(self._params)):
      raise ValueError("Dimensionality of count vector is incorrect")
    prob = 1.
    for i,c in enumerate(counts):
      prob *= self._params[i]**counts[i]
    return (prob * numpy.exp(self._log_multinomial_coeff(counts)))


  def log_pmf(self,counts):
    if not(len(counts)==len(self._params)):
      raise ValueError("Dimensionality of count vector is incorrect")
    counts = numpy.array(counts, dtype = numpy.int)
    prob = 0 + numpy.sum(counts*numpy.log(self._params))
    return (prob + self._log_multinomial_coeff(counts))


  def _log_multinomial_coeff(self, counts):
    return self._log_factorial(numpy.sum(counts)) - numpy.sum(self._log_factorial(c) for c in counts)


  def _log_factorial(self, num):
    #if not num > 0:
    #  raise ValueError("Can only compute the factorial of positive ints")
    #n = numpy.arange(1,num+1)
    return log_fac_lookup[num]