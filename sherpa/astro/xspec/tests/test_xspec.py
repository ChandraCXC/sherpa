#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
import numpy
import sherpa.astro.xspec as xs
from sherpa.utils import SherpaTestCase

import logging
error = logging.getLogger(__name__).error

def is_proper_subclass(obj, cls):
    if type(cls) is not tuple:
        cls = (cls,)
    if obj in cls:
        return False
    return issubclass(obj, cls)


class test_xspec(SherpaTestCase):

    def test_create_model_instances(self):
        count = 0

        for cls in dir(xs):
            if not cls.startswith('XS'):
                continue

            cls = getattr(xs, cls)

            if is_proper_subclass(cls, (xs.XSAdditiveModel,
                                        xs.XSMultiplicativeModel)):
                m = cls()
                count += 1

        self.assertEqual(count, 133)

    def test_evaluate_model(self):
        m = xs.XSbbody()
        out = m([1,2,3,4])
        if m.calc.__name__.startswith('C_'):
            otype = numpy.float64
        else:
            otype = numpy.float32
        self.assert_(out.dtype.type is otype)
        self.assertEqual(int(numpy.flatnonzero(out == 0.0)), 3)


    def test_xspec_models(self):
        models = [model for model in dir(xs) if model[:2] == 'XS']
        models.remove('XSModel')
        models.remove('XSMultiplicativeModel')
        models.remove('XSAdditiveModel')
        models.remove('XSTableModel')

        xx = numpy.arange(0.1, 11.01, 0.01, dtype=float)
        xlo = numpy.array(xx[:-1])
        xhi = numpy.array(xx[1:])
        for model in models:
            cls = getattr(xs, model)
            foo = cls('foo')
            vals = foo(xlo,xhi)
            try:
                self.assert_(not numpy.isnan(vals).any() and
                             not numpy.isinf(vals).any() )
            except AssertionError:
                error('XS%s model evaluation failed' % model)
                raise


if __name__ == '__main__':

    from sherpa.utils import SherpaTest
    SherpaTest(xs).test()
