#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2011)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

import numpy
from sherpa.models.parameter import Parameter, tinyval
from sherpa.models.model import ArithmeticModel, modelCacher1d
from sherpa.utils.err import ModelErr
from sherpa.utils import SherpaFloat, sao_fcmp

_tol = numpy.finfo(numpy.float32).eps

# Optical Models for SED Analysis
#
# This Sherpa Python module contains optical models for fitting to SEDs.
# These models are Python versions of models found in the Specview
# application for analyzing spectra and SEDs.  These models are meant
# to be used in conjunction with Specview to serve the VAO SED project.
#
# These models work in wavelength space (Angstroms).
#

__all__ = ('AbsorptionEdge', 'AccretionDisk', 'AbsorptionGaussian', 'AbsorptionLorentz', 'EmissionLorentz', 'OpticalGaussian', 'EmissionGaussian', 'AbsorptionVoigt', 'BlackBody', 'Bremsstrahlung', 'BrokenPowerlaw', 'CCM', 'LogAbsorption', 'LogEmission', 'Polynomial', 'Powerlaw', 'Recombination', 'EmissionVoigt', 'XGal')

# The speed of light in km/s
c_km = 2.99792458e+5


# This model sets in edge (in Angstroms) beyond which absorption
# is a significant feature to the spectrum or SED.
class AbsorptionEdge(ArithmeticModel):

    def __init__(self, name='absorptionedge'):
        self.edgew = Parameter(name, 'edgew', 5000., units='angstroms')
        self.tau = Parameter(name, 'tau', 0.5)
        self.index = Parameter(name, 'index', 3.0, alwaysfrozen=True,
                               hidden=True)

        ArithmeticModel.__init__(self, name,
                                 (self.edgew, self.tau, self.index))

    # We can turn on model caching with this commented-out feature,
    # if we find we need it.
    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        y = numpy.zeros_like(x)

        if sao_fcmp(p[0], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s tau cannot be zero' % self.name)

        idx = (x <= p[0])
        y[idx] = numpy.exp(-(p[1]*numpy.power(x[idx]/p[0], p[2])))
        return y

# This model is an accretion disk continuum function.
class AccretionDisk(ArithmeticModel):

    def __init__(self, name='accretiondisk'):

        self.ref = Parameter(name, 'ref', 5000., units='angstroms')
        self.beta = Parameter(name, 'beta', 0.5, -10, 10)
        self.ampl = Parameter(name, 'ampl', 1.)
        self.norm = Parameter(name, 'norm', 20000.0, 0, alwaysfrozen=True,
                              hidden=True)

        ArithmeticModel.__init__(self, name,
                                 (self.ref, self.beta, self.ampl, self.norm))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if 0.0 in x:
            raise ValueError('model evaluation failed, ' +
                             'x cannot be zero')

        if p[3] == 0.0:
            raise ValueError('model evaluation failed, ' +
                             'norm cannot be zero')

        return p[2]*numpy.power(x/p[3], -p[1])*numpy.exp(-p[0]/x)


# This model calculates a Gaussian function expressed in
# equivalent width, and models absorption due to this Gaussian.
class AbsorptionGaussian(ArithmeticModel):
    """Absorption Gaussian function expressed in equivalent width."""

    def __init__(self, name='absorptiongaussian'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, 
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., 0, units='angstroms')
        self.ewidth = Parameter(name, 'ewidth', 1.)
        self.limit = Parameter(name, 'limit', 4., alwaysfrozen=True,
                               hidden=True )

        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos,
                                              self.ewidth, self.limit))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if sao_fcmp(p[0], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if sao_fcmp(p[1], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        y = numpy.ones_like(x)
        sigma = p[1] * p[0] / 705951.5     # = 2.9979e5 / 2.354820044 ?
        delta = numpy.abs((x - p[1]) / sigma)
        ampl  = p[2] / sigma / 2.50662828  # document this constant

        idx = (delta < p[3])
        y[idx] = 1.0 - ampl * numpy.exp(- delta * delta / 2.0)

        return y

# This model calculates a Lorentzian function expressed in
# equivalent width, and models absorption due to this Lorentzian.
class AbsorptionLorentz(ArithmeticModel):
    """Absorption Lorentz function expressed in equivalent width."""

    def __init__(self, name='absorptionlorentz'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, 
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., 0, units='angstroms')
        self.ewidth = Parameter(name, 'ewidth', 1.)

        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos, self.ewidth))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if sao_fcmp(p[0], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if sao_fcmp(p[1], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        y = (1.0 / x - 1.0 / p[1]) * p[1] * c_km / p[0]
        y = 1.0 + 4.0 * y * y
        y *= 1.571 * p[0] * p[1] / c_km
        y = 1.0 - p[2] / y
        return y

# This model computes a Lorentzian profile for emission features.
class EmissionLorentz(ArithmeticModel):
    """Emission Lorentz function expressed in equivalent width."""

    def __init__(self, name='emissionlorentz'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval,
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., 0, units='angstroms')
        self.flux = Parameter(name, 'flux', 1.)
        self.kurt = Parameter(name, 'kurt', 1.)

        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos,
                                              self.flux, self.kurt))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if sao_fcmp(p[0], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if sao_fcmp(p[1], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        sigma = p[0] * p[1] / c_km
        arg = numpy.power(numpy.abs(x - p[1]), p[3]) + sigma/2.0 * sigma/2.0

        arg[arg<1.0e-15] = 1.0e-15
        
        return p[2] * sigma / arg / (numpy.pi*2)

# This model computes an absorption Gaussian feature expressed in
# optical depth.
class OpticalGaussian(ArithmeticModel):
    """Absorption Gaussian function expressed in optical depth."""

    def __init__(self, name='opticalgaussian'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, 
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., 0, units='angstroms')
        self.tau = Parameter(name, 'tau', 0.5)
        self.limit = Parameter(name, 'limit', 4., alwaysfrozen=True,
                               hidden=True )

        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos,
                                              self.tau, self.limit))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if sao_fcmp(p[0], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if sao_fcmp(p[1], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        y = numpy.ones_like(x)
        sigma = p[1] * p[0] / 705951.5     # = 2.9979e5 / 2.354820044 ?
        delta = numpy.abs((x - p[1]) / sigma)

        idx = (delta < p[3])
        y[idx] = numpy.exp(-p[2] * numpy.exp(- delta * delta / 2.0))

        return y

# This model computes a Gaussian profile for emission features.
class EmissionGaussian(ArithmeticModel):
    """Emission Gaussian function."""

    def __init__(self, name='emissiongaussian'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval,
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., 0, units='angstroms')
        self.flux = Parameter(name, 'flux', 1.)
        self.skew = Parameter(name, 'skew', 1.)
        self.limit = Parameter(name, 'limit', 4., alwaysfrozen=True,
                               hidden=True )

        ArithmeticModel.__init__(self, name,
                                 (self.fwhm, self.pos, self.flux,
                                  self.skew, self.limit))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if sao_fcmp(p[0], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if sao_fcmp(p[1], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        if sao_fcmp(p[3], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s skew cannot be zero' % self.name)

        y = numpy.zeros_like(x)
        sigma = p[1] * p[0] / 705951.5     # = 2.9979e5 / 2.354820044
        delta = numpy.abs((x - p[1]) / sigma)
        ampl  = p[2] / sigma / 2.50662828  # document this constant
        idx = (delta < p[3])

        arg = - delta * delta / 2.0
        if sao_fcmp(p[3], 1.0, _tol) == 0:
            y[idx] = p[0] * numpy.exp(arg) / sigma / 2.50662828

        else:
            left = (arg <= p[1])
            arg[left] = numpy.exp(arg[left])
            right = ~left
            arg[right] = numpy.exp(arg[right] / p[3] / p[3])
            y[idx] = 2.0 * p[0] * arg / sigma / 2.50662828 / (1.0 + p[3])

        return y

# This model computes absorption as a Voigt function -- i.e., with
# a Gaussian core and Lorentzian wings.
class AbsorptionVoigt(ArithmeticModel):
    """Absorption Voigt function expressed in equivalent width."""

    def __init__(self, name='absorptionvoigt'):
        self.center = Parameter(name, 'center', 5000., tinyval, hard_min=tinyval, units="angstroms")
        self.ew = Parameter(name, 'ew', 1., tinyval, hard_min=tinyval, units="angstroms")
        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, units="km/s")
        self.lg = Parameter(name, 'lg', 1., tinyval, hard_min=tinyval)

        # Create core and wings from Gaussian and Lorentz
        self._core = AbsorptionGaussian()
        self._wings = AbsorptionLorentz()
        
        ArithmeticModel.__init__(self, name, (self.center, self.ew, self.fwhm, self.lg))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        # Combining two absorption components means
        # multiplying them, it appears (at least according
        # to addAbsorption in NarrowBandFunction.java)

        core_pars = numpy.array([ p[2], p[0], p[1] / 2.0, 4.0 ])
        wing_pars = numpy.array([ p[2], p[0] * p[3], p[1] / 2.0 ])
        return self._core.calc(core_pars, x) * self._wings.calc(wing_pars, x)

# This model computes continuum emission as a blackbody function.
class BlackBody(ArithmeticModel):
    """Blackbody model."""

    def __init__(self, name='blackbody'):
        self.refer = Parameter(name, 'refer', 5000., tinyval, hard_min=tinyval, units="angstroms")
        self.ampl = Parameter(name, 'ampl', 1., tinyval, hard_min=tinyval, units="angstroms")
        self.temperature = Parameter(name, 'temperature', 3000., tinyval, hard_min=tinyval, units="Kelvin")

        self._argmin = 1.0e-3
        self._argmax = 88.0

        ArithmeticModel.__init__(self, name, (self.refer, self.ampl, self.temperature))

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        c1 = 1.438786e8
        efactor = c1 / p[2]
        #if ((efactor / p[0]) > self._argmax):
            # raise error exp too big
        #    raise ValueError('model evaluation failed, either temperature or reference wavelength too small')

        #numer = p[1] * numpy.power(p[0], 5.0) * (numpy.exp(efactor / p[0]) - 1.0)
        numer = p[1] * numpy.power(p[0], 5.0) * numpy.expm1(efactor / p[0])

        y = numpy.zeros_like(x)
        x0 = numpy.where(x > 0.0)
        if (len(x0[0]) > 0):
            arg = numpy.zeros_like(x)
            arg[x0] = efactor / x[x0]
            denon = numpy.zeros_like(x)
            denon[x0] = numpy.power(x[x0], 5)
            argmin_slice = numpy.where(arg < self._argmin)
            if (len(argmin_slice[0]) > 0): 
                denon[argmin_slice] *= arg[argmin_slice] * (1.0 + 0.5 * arg[argmin_slice])
            numpy.where(arg > self._argmax, self._argmax, arg)

            arg_slice = numpy.where(arg >= self._argmin)
            if (len(arg_slice[0]) > 0):
                denon[arg_slice] *= numpy.exp(arg[arg_slice]) - 1.0

            y[x0] = numer / denon[x0]

        return y

# This model computes continuum emission with the bremsstrahlung function.
class Bremsstrahlung(ArithmeticModel):
    """Bremsstrahlung model."""

    def __init__(self, name='bremsstrahlung'):
        self.refer = Parameter(name, 'refer', 5000., tinyval, hard_min=tinyval, units="angstroms")
        self.ampl = Parameter(name, 'ampl', 1., tinyval, hard_min=tinyval, units="angstroms")
        self.temperature = Parameter(name, 'temperature', 3000., tinyval, hard_min=tinyval, units="Kelvin")

        ArithmeticModel.__init__(self, name, (self.refer, self.ampl, self.temperature))

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        if sao_fcmp(p[0], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s refer cannot be zero' % self.name)
        
        return p[1] * numpy.power((p[0] / x), 2) * numpy.exp(-1.438779e8 / x / p[2])

# This model computes continuum emission with a broken power-law;
# that is, the power-law index changes after a break at a particular
# wavelength.
class BrokenPowerlaw(ArithmeticModel):
    """Broken power-law model."""
    
    def __init__(self, name='brokenpowerlaw'):
        self.refer = Parameter(name, 'refer', 5000., tinyval, hard_min=tinyval, units="angstroms")
        self.ampl = Parameter(name, 'ampl', 1., tinyval, hard_min=tinyval, units="angstroms")
        self.index1 = Parameter(name, 'index1', 0.1, -10.0, 10.0)
        self.index2 = Parameter(name, 'index2', -0.1, -10.0, 10.0)

        ArithmeticModel.__init__(self, name, (self.refer, self.ampl, self.index1, self.index2))

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        if sao_fcmp(p[0], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s refer cannot be zero' % self.name)
        
        x = numpy.asarray(x, dtype=SherpaFloat)
        arg = x / p[0]
        arg = p[1] * (numpy.where(arg > 1.0, numpy.power(arg, p[3]), numpy.power(arg, p[2])))
        return arg

# This model computes extinction using the function published by
# Cardelli, Clayton, and Mathis
# (ApJ, 1989, vol 345, pp 245)
class CCM(ArithmeticModel):
    """Cardelli, Clayton, and Mathis extinction curve."""
    
    def __init__(self, name='ccm'):
        self.ebv = Parameter(name, 'ebv', 0.5)
        self.r = Parameter(name, 'r', 3.2)

        ArithmeticModel.__init__(self, name, (self.ebv, self.r))

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        y = numpy.zeros_like(x)
        y2 = numpy.zeros_like(x)
        y3 = numpy.zeros_like(x)

        a = numpy.zeros_like(x)
        b = numpy.zeros_like(x)

        x = 1000.0 / (x / 10.0)

        # Infrared wavelengths
        xp = numpy.zeros_like(x)
        ir_slice = numpy.concatenate((numpy.where(x >= 0.3)[0], numpy.where(x <= 1.1)[0]), 0)
        if (len(ir_slice) > 0):
            xp[ir_slice] = numpy.power(x[ir_slice], 1.61)
            a[ir_slice] = 0.574 * xp[ir_slice]
            b[ir_slice] = -0.527 * xp[ir_slice]

        # Optical
        opt_slice = numpy.concatenate((numpy.where(x > 1.1)[0], numpy.where(x <= 3.3)[0]), 0)
        if (len(opt_slice) > 0):
            y[opt_slice] = x[opt_slice] - 1.82

            a[opt_slice] = 1.0 + 0.17699 * y[opt_slice] - 0.50477 * y[opt_slice] * y[opt_slice] - 0.02427 * numpy.power(y[opt_slice],3) + 0.72085 * numpy.power(y[opt_slice],4) + 0.01979 * numpy.power(y[opt_slice],5) - 0.77530 * numpy.power(y[opt_slice],6) + 0.32999 * numpy.power(y[opt_slice],7)

            b[opt_slice] = 0.0 + 1.41338 * y[opt_slice] + 2.28305 * y[opt_slice] * y[opt_slice] + 1.07233 * numpy.power(y[opt_slice],3) - 5.38434 * numpy.power(y[opt_slice],4) - 0.62551 * numpy.power(y[opt_slice],5) + 5.30260 * numpy.power(y[opt_slice],6) - 2.09002 * numpy.power(y[opt_slice],7)

        # Near-UV
        nuv_slice = numpy.concatenate((numpy.where(x > 3.3)[0], numpy.where(x <= 8.0)[0]), 0)
        if (len(nuv_slice) > 0):
            a[nuv_slice] = 0.0
            b[nuv_slice] = 0.0

            nuv_slice2 = numpy.concatenate((numpy.where(x >= 5.9)[0], numpy.where(x <= 8.0)[0]), 0)
            if (len(nuv_slice2) > 0):
                y[nuv_slice2] = x[nuv_slice2] - 5.9
                y2[nuv_slice2] = y[nuv_slice2] * y[nuv_slice2]
                y3[nuv_slice2] = y2[nuv_slice2] * y[nuv_slice2]

                a[nuv_slice2] = -0.04473 * y2[nuv_slice2] - 0.009779 * y3[nuv_slice2]
                b[nuv_slice2] = 0.21300 * y2[nuv_slice2] + .120700 * y3[nuv_slice2]

            a[nuv_slice] = a[nuv_slice] + 1.752 - 0.316 * x[nuv_slice] - 0.104 / (0.341 + numpy.power((x[nuv_slice]-4.67),2))

            b[nuv_slice] = b[nuv_slice] - 3.090 + 1.825 * x[nuv_slice] + 1.206 / (0.263 + numpy.power((x[nuv_slice]-4.62),2))

        # Far-UV
        fuv_slice = numpy.concatenate((numpy.where(x > 8.0)[0], numpy.where(x <= 20.0)[0]), 0)
        if (len(fuv_slice) > 0):
            y[fuv_slice] = x[fuv_slice] - 8.0
            y2[fuv_slice] = y[fuv_slice] * y[fuv_slice]
            y3[fuv_slice] = y2[fuv_slice] * y[fuv_slice]

            a[fuv_slice] = -1.073 - 0.628 * y[fuv_slice] + 0.137 * y2[fuv_slice] - 0.070 * y3[fuv_slice]
            b = 13.670 + 4.257 * y[fuv_slice] - 0.420 * y2[fuv_slice] + 0.374 * y3[fuv_slice]

        # Final extinction curve
        aext = p[1] * a + b
        return numpy.power(10.0,(-0.4 * p[0] * aext))

# This model computes absorption using a Gaussian function expressed
# in optical depth, and using the log of the FWHM.
class LogAbsorption(ArithmeticModel):
    """Log of absorption Gaussian function expressed in optical depth."""

    def __init__(self, name='logabsorption'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, 
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., 0, units='angstroms')
        self.tau = Parameter(name, 'tau', 0.5)

        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos,
                                              self.tau))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if sao_fcmp(p[0], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if sao_fcmp(p[1], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        y = numpy.ones_like(x)

        alpha = 0.69314718 / numpy.log(1.0 + p[0] / 2.9979e5 / 2.0)
        if (alpha <= 1.0):
            alpha = 1.0001

        y = numpy.where(x >= p[1], p[2] * numpy.power((x / p[1]), -alpha), p[2] * numpy.power((x / p[1]), alpha))
        
        return numpy.exp(-y)

# This model computes emission using a Gaussian function expressed
# in optical depth, and using the log of the FWHM.
class LogEmission(ArithmeticModel):
    """Log of the emission Gaussian function."""

    def __init__(self, name='logemission'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval,
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., 0, units='angstroms')
        self.flux = Parameter(name, 'flux', 1.)
        self.skew = Parameter(name, 'skew', 1.)
        self.limit = Parameter(name, 'limit', 4., alwaysfrozen=True,
                               hidden=True )

        ArithmeticModel.__init__(self, name,
                                 (self.fwhm, self.pos, self.flux,
                                  self.skew, self.limit))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if sao_fcmp(p[0], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if sao_fcmp(p[1], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        if sao_fcmp(p[3], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s skew cannot be zero' % self.name)

        arg = 0.69314718 / numpy.log (1.0 + p[0] / 2.9979e5 / 2.0);
        if (arg <= 1.0):
            arg = 1.0001

        fmax = (arg - 1.0) * p[2] / p[1] / 2.0

        if (p[3] == 1.0):
            return numpy.where(x >= p[1], fmax * numpy.power((x / p[1]), -arg), fmax * numpy.power((x / p[1]), arg)) 

        arg1 = 0.69314718 / numpy.log (1.0 + p[3] * p[0] / 2.9979e5 / 2.0)
        fmax = (arg - 1.0) * p[2] / p[1] / (1.0 + (arg - 1.0) / (arg1 - 1.0))

        return numpy.where(x <= p[1], fmax * numpy.power((x / p[1]), arg), fmax * numpy.power((x / p[1]), -arg1))

# This model computes continuum emission as a polynomial,
# y = c0 + c1 * (x - offset) + c2 * (x - offset)^2 + c3 * (x - offset)^3 + c4 * (x - offset)^4 + c5 * (x - offset)^5
class Polynomial(ArithmeticModel):
    """Polynomial model."""
    
    def __init__(self, name='polynomial'):
        pars = []
        
        for i in xrange(6):
            pars.append(Parameter(name, 'c%d' % i, 0, frozen=True))
        pars[0].val = 1
        pars[0].frozen = False
        for p in pars:
            setattr(self, p.name, p)

        self.offset = Parameter(name, 'offset', 0, frozen=True)
        pars.append(self.offset)
        ArithmeticModel.__init__(self, name, pars)

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        y = numpy.zeros_like(x)
        xtemp = x - p[6]
        y += p[5]
        for i in [4,3,2,1,0]:
            y = y*xtemp + p[i]

        return y

# This model computes continuum emission using a power-law.
class Powerlaw(ArithmeticModel):
    """Power-law model."""
    
    def __init__(self, name='powerlaw'):
        self.refer = Parameter(name, 'refer', 5000., tinyval, hard_min=tinyval, units="angstroms")
        self.ampl = Parameter(name, 'ampl', 1., tinyval, hard_min=tinyval, units="angstroms")
        self.index = Parameter(name, 'index', -0.5, -10.0, 10.0)

        ArithmeticModel.__init__(self, name, (self.refer, self.ampl, self.index))

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        if sao_fcmp(p[0], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s refer cannot be zero' % self.name)
        
        x = numpy.asarray(x, dtype=SherpaFloat)
        arg = x / p[0]
        arg = p[1] * numpy.power(arg, p[2])

        return arg

# This model computes the continuum with an optically thin
# recombination function.
class Recombination(ArithmeticModel):
    """Optically thin recombination model."""
    
    def __init__(self, name='recombination'):
        self.refer = Parameter(name, 'refer', 5000., tinyval, hard_min=tinyval, units="angstroms")
        self.ampl = Parameter(name, 'ampl', 1., tinyval, hard_min=tinyval, units="angstroms")
        self.temperature = Parameter(name, 'temperature', 3000., tinyval, hard_min=tinyval, units="Kelvin")
        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, units="km/s")

        ArithmeticModel.__init__(self, name, (self.refer, self.ampl, self.temperature, self.fwhm))

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        if sao_fcmp(p[0], 0.0, _tol) == 0:
            raise ValueError('model evaluation failed, ' +
                             '%s refer cannot be zero' % self.name)
        
        x = numpy.asarray(x, dtype=SherpaFloat)
        sigma = p[0] * p[3] / 705951.5; # = 2.9979e5 / 2.354820044
        delta = 1.440e8 * (1.0 / x - 1.0 / p[0]) / p[2]

        return numpy.where(delta < 0.0, p[1] * numpy.exp(-numpy.power((x - p[0]), 2.0) / numpy.power(sigma, 2.0) / 2.0), p[1] * numpy.power((p[0] / x), 2.0) * numpy.exp(-delta))

# This model computes emission as a Voigt function -- i.e., with
# a Gaussian core and Lorentzian wings.
class EmissionVoigt(ArithmeticModel):
    """Emission Voigt function."""

    def __init__(self, name='emissionvoigt'):
        self.center = Parameter(name, 'center', 5000., tinyval, hard_min=tinyval, units="angstroms")
        self.flux = Parameter(name, 'flux', 1.)
        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, units="km/s")
        self.lg = Parameter(name, 'lg', 1., tinyval, hard_min=tinyval)

        # Create core and wings from Gaussian and Lorentz
        self._core = EmissionGaussian()
        self._wings = EmissionLorentz()
        
        ArithmeticModel.__init__(self, name, (self.center, self.flux, self.fwhm, self.lg))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        # Combining two emission components means
        # adding them, it appears (at least according
        # to addEmission in NarrowBandFunction.java)

        core_pars = numpy.array([ p[2], p[0], p[1], 1.0])
        wing_pars = numpy.array([ p[3]*p[2], p[0], p[1], 2.0])
        return self._core.calc(core_pars, x) + self._wings.calc(wing_pars, x)

# This model computes the extragalactic extinction function of
# Calzetti, Kinney and Storchi-Bergmann, 1994, ApJ, 429, 582
class XGal(ArithmeticModel):
    """Extragalactic extinction function of Calzetti, Kinney and Storchi-Bergmann"""
    
    def __init__(self, name='xgal'):
        self.ebv = Parameter(name, 'ebv', 0.5)

        ArithmeticModel.__init__(self, name, (self.ebv,))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        x = 1000.0 / x

        # Formula from paper with zero point moved to (x = 0)
        ext = ((0.011 * x - 0.198) * x + 1.509) * x

        # Normalize the result according to Kailash Sahu's calculations
        ext *= 2.43

        return numpy.power(10.0,(-0.4 * p[0] * ext))
