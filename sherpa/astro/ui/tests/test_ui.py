#/usr/bin/env python

from sherpa.utils import SherpaTest, SherpaTestCase, needs_data
import sherpa.astro.ui as ui

class test_ui(SherpaTestCase):

    @needs_data
    def setUp(self):
        self.ascii = self.datadir + '/threads/ascii_table/sim.poisson.1.dat'
        self.fits = self.datadir + '/1838_rprofile_rmid.fits'
        self.singledat = self.datadir + '/single.dat'
        self.singletbl = self.datadir + '/single.fits'
        self.doubledat = self.datadir + '/double.dat'
        self.doubletbl = self.datadir + '/double.fits'
        self.img = self.datadir + '/img.fits'
        self.filter_single_int_ascii = self.datadir + '/filter_single_integer.dat'
        self.filter_single_int_table = self.datadir + '/filter_single_integer.fits'
        self.filter_single_log_table = self.datadir + '/filter_single_logical.fits'

        self.func = lambda x: x
        ui.dataspace1d(1,1000,dstype=ui.Data1D)


    @needs_data
    def test_ascii(self):
        ui.load_ascii(1, self.ascii)
        ui.load_ascii(1, self.ascii, 2)
        ui.load_ascii(1, self.ascii, 2, ("col2", "col1"))


    @needs_data
    def test_table(self):
        ui.load_table(1, self.fits)
        ui.load_table(1, self.fits, 3)
        ui.load_table(1, self.fits, 3, ["RMID","SUR_BRI","SUR_BRI_ERR"])
        ui.load_table(1, self.fits, 4, ('R',"SUR_BRI",'SUR_BRI_ERR'),
                      ui.Data1DInt)


    # Test table model
    @needs_data
    def test_table_model_ascii_table(self):
        ui.load_table_model('tbl', self.singledat)
        ui.load_table_model('tbl', self.doubledat)


    @needs_data
    def test_table_model_fits_table(self):
        ui.load_table_model('tbl', self.singletbl)
        ui.load_table_model('tbl', self.doubletbl)


    @needs_data
    def test_table_model_fits_image(self):
        ui.load_table_model('tbl', self.img)


    # Test user model
    @needs_data
    def test_user_model_ascii_table(self):
        ui.load_user_model(self.func, 'mdl', self.singledat)
        ui.load_user_model(self.func, 'mdl', self.doubledat)


    @needs_data
    def test_user_model_fits_table(self):
        ui.load_user_model(self.func, 'mdl', self.singletbl)
        ui.load_user_model(self.func, 'mdl', self.doubletbl)


    @needs_data
    def test_filter_ascii(self):
        ui.load_filter(self.filter_single_int_ascii)
        ui.load_filter(self.filter_single_int_ascii, ignore=True)


    # Test load_filter
    @needs_data
    def test_filter_table(self):
        ui.load_filter(self.filter_single_int_table)
        ui.load_filter(self.filter_single_int_table, ignore=True)

        ui.load_filter(self.filter_single_log_table)
        ui.load_filter(self.filter_single_log_table, ignore=True)


if __name__ == '__main__':

    import sys
    if len(sys.argv) > 1:
        SherpaTest(ui).test(datadir=sys.argv[1])
