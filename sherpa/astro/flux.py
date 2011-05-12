#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2009)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

#!/usr/bin/env python
import numpy
import numpy.random
from sherpa.estmethods import *
from sherpa.astro.utils import calc_energy_flux
from itertools import izip
from sherpa.utils import parallel_map
from sherpa.utils.err import EstErr

import logging
warning = logging.getLogger("sherpa").warning


__all__ = ['get_sample_uncorr', 'get_sample_corr',
           'sample_param_uncorr', 'sample_param_corr',
           'calc_flux', 'sample_flux']

def get_sample_uncorr(val, sigma, num=1, dist=numpy.random.normal):
    return dist(val, sigma, num)


def get_sample_corr(vals, cov, num=1, dist=numpy.random.multivariate_normal):
    return dist(vals, cov, num)


def sample_param_uncorr(fit, num=1, dist=numpy.random.normal):
    samples=[]
    oldestmethod = fit.estmethod

    fit.estmethod = Covariance()
    try:
        r = fit.est_errors()
    finally:
        fit.estmethod = oldestmethod

    thawedpars = [par for par in fit.model.pars if not par.frozen]
    for par, val, lo, hi in izip(thawedpars, r.parvals, r.parmins, r.parmaxes):
        sigma = None
        if lo is not None and hi is not None:
            sigma = numpy.abs(lo)
        else:
            warning("Covariance failed for '%s', trying Confidence..." %
                    par.fullname)
            fit.estmethod = Confidence()
            try:
                t = fit.est_errors(parlist = (par,))
                if t.parmins[0] is not None and t.parmaxes[0] is not None:
                    sigma = numpy.abs(t.parmins[0])
                else:
                    warning('1 sigma bounds for parameter ' +
                            par.fullname +
                            ' could not be found, using soft limit minimum')
                    sigma = numpy.abs(par.min)
            finally:
                fit.estmethod = oldestmethod
        samples.append(get_sample_uncorr(val, sigma, num, dist))
    samples = numpy.asarray(samples).transpose()
    return samples


def sample_param_corr(fit, num=1, dist=numpy.random.multivariate_normal):
    oldestmethod = fit.estmethod
    fit.estmethod = Covariance()
    try:
        r = fit.est_errors()
    finally:
        fit.estmethod = oldestmethod

    cov = r.extra_output
    if cov is None:
        raise EstErr('nocov')

    cov = numpy.asarray(cov)
    if numpy.min(numpy.linalg.eigvalsh(cov)) <= 0:
        raise TypeError("The covariance matrix is not positive definite")

    vals = fit.model.thawedpars
    samples = get_sample_corr(vals, cov, num, dist)

    return samples


def calc_flux(fit, data, src, samples, method=calc_energy_flux,
              lo=None, hi=None, numcores=None):

    def evaluate(sample):
        fit.model.thawedpars = sample
        flux = method(data, src, lo, hi)
        return [flux] + list(sample)

    old_model_vals  = fit.model.thawedpars
    try:
        fluxes = parallel_map(evaluate, samples, numcores)
    finally:
        fit.model.thawedpars = old_model_vals

    return numpy.asarray(fluxes)


def sample_flux(fit, data, src, method=calc_energy_flux, correlated=False,
                num=1, lo=None, hi=None, numcores=None):
    samples=None
    if correlated:
        samples = sample_param_corr(fit, num, numpy.random.multivariate_normal)
    else:
        samples = sample_param_uncorr(fit, num, numpy.random.normal)
    return calc_flux(fit, data, src, samples, method, lo, hi, numcores)
