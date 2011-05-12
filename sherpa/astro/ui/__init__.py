#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
import sherpa.all
import sherpa.astro.all
import sherpa.astro.ui.utils
from sherpa.utils import calc_mlr, calc_ftest, rebin, histogram1d, \
    histogram2d, gamma, lgam, erf, igamc, igam, incbet, multinormal_pdf, \
    multit_pdf
from sherpa.data import Data1D, Data1DInt, Data2D, Data2DInt
from sherpa.astro.data import DataARF, DataRMF, DataPHA, DataIMG, DataIMGInt
from sherpa.logposterior import Prior


# We build up __all__ as we go along
__all__ = ['DataARF', 'DataRMF','DataPHA', 'DataIMG', 'DataIMGInt', 'Data1D',
           'Data1DInt','Data2D', 'Data2DInt', 'calc_mlr', 'calc_ftest', 'rebin',
           'histogram1d', 'histogram2d', 'gamma', 'lgam', 'erf', 'igamc',
           'igam', 'incbet', 'Prior', 'multinormal_pdf', 'multit_pdf']

_session = utils.Session()
_session._add_model_types(sherpa.models.basic)
_session._add_model_types(sherpa.astro.models)
# To add PSFModel to list -- doesn't inherit from ArithmeticModel
_session._add_model_types(sherpa.instrument,baselist=(sherpa.models.Model,))
# Get RMFModel, ARFModel in list of models
_session._add_model_types(sherpa.astro.instrument)

if hasattr(sherpa.astro, 'xspec'):
    _session._add_model_types(sherpa.astro.xspec,
                              (sherpa.astro.xspec.XSAdditiveModel,
                               sherpa.astro.xspec.XSMultiplicativeModel))

    from sherpa.astro.xspec import get_xsabund, get_xscosmo, get_xsxsect, \
         set_xsabund, set_xscosmo, set_xsxsect, set_xsxset, get_xsxset, \
         get_xschatter, set_xschatter
    __all__.extend(('get_xsabund', 'get_xschatter', 'get_xscosmo',
                    'get_xsxsect', 'set_xsabund', 'set_xschatter',
                    'set_xscosmo', 'set_xsxsect', 'set_xsxset',
                    'get_xsxset'))

__all__.extend(_session._export_names(globals()))


__all__ = tuple(__all__)
