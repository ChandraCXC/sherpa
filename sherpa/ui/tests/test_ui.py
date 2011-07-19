#/usr/bin/env python

from sherpa.utils import SherpaTest, SherpaTestCase, needs_data
from sherpa.models import ArithmeticModel, Parameter
import sherpa.ui as ui
import numpy

class UserModel(ArithmeticModel):

    def __init__(self, name='usermodel'):
        self.param1 = Parameter(name, 'param1', 1, min=0, max=100)
        self.param2 = Parameter(name, 'param2', 1, min=-100, max=100)

        ArithmeticModel.__init__(self, name, (self.param1,
                                              self.param2))

    def calc(self, p, x, *args, **kwargs):
        return p[0]*x+p[1]



class test_ui(SherpaTestCase):

    @needs_data
    def setUp(self):
        self.ascii = self.datadir + '/threads/ascii_table/sim.poisson.1.dat'
        self.single = self.datadir + '/single.dat'
        self.double = self.datadir + '/double.dat'
        self.filter = self.datadir + '/filter_single_integer.dat'
        self.func = lambda x: x
        
        ui.dataspace1d(1,1000,dstype=ui.Data1D)

    @needs_data
    def test_ascii(self):
        ui.load_data(1, self.ascii)
        ui.load_data(1, self.ascii, 2)
        ui.load_data(1, self.ascii, 2, ("col2", "col1"))


    # Test table model
    @needs_data
    def test_table_model_ascii_table(self):
        ui.load_table_model('tbl', self.single)
        ui.load_table_model('tbl', self.double)


    # Test user model
    @needs_data
    def test_user_model_ascii_table(self):
        ui.load_user_model(self.func, 'mdl', self.single)
        ui.load_user_model(self.func, 'mdl', self.double)


    @needs_data
    def test_filter_ascii(self):
        ui.load_filter(self.filter)
        ui.load_filter(self.filter, ignore=True)

    @needs_data
    def test_add_model(self):
        ui.add_model(UserModel)
        ui.set_model('usermodel.user1')

    @needs_data
    def test_set_full_model(self):
        ui.load_psf('psf1', 'gauss2d.g1')
        ui.set_full_model('psf1(gauss2d.g2)+const2d.c1')
        ui.get_model()
        ui.get_source()



class test_psf_ui(SherpaTestCase):

    models1d = ['gauss1d', 'delta1d', 'normgauss1d']
    models2d = ['gauss2d', 'delta2d', 'normgauss2d']

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_psf_model2d(self):
        ui.dataspace1d(1, 10)
        for model in self.models1d:
            try:
                ui.load_psf('psf1d', model+'.mdl')
                ui.set_psf('psf1d')
                mdl = ui.get_model_component('mdl')
                self.assert_( (numpy.array(mdl.get_center()) ==
                               numpy.array([4])).all() )
            except:
                print model
                raise


    def test_psf_model2d(self):
        ui.dataspace2d([216,261])
        for model in self.models2d:
            try:
                ui.load_psf('psf2d', model+'.mdl')
                ui.set_psf('psf2d')
                mdl = ui.get_model_component('mdl')
                self.assert_( (numpy.array(mdl.get_center()) ==
                               numpy.array([108,130])).all() )
            except:
                print model
                raise


if __name__ == '__main__':

    import sys
    if len(sys.argv) > 1:
        SherpaTest(ui).test(datadir=sys.argv[1])
