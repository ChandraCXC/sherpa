#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
import os.path
import sherpa
from sherpa.utils import SherpaTestCase


class test_sherpa(SherpaTestCase):

    def test_include_dir(self):
        incdir = os.path.join(sherpa.get_include(), 'sherpa')
        self.assert_(os.path.isdir(incdir))
