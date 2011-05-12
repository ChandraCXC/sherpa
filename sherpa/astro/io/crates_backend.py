#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2008)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
from itertools import izip
import os.path
import copy
import numpy
import pycrates
from sherpa.utils.err import IOErr
from sherpa.utils import SherpaInt, SherpaUInt, SherpaFloat, is_binary_file
from sherpa.astro.utils import resp_init
from sherpa.astro.io.meta import *

import logging
warning = logging.getLogger(__name__).warning
error = logging.getLogger(__name__).error
info    = logging.getLogger(__name__).info

transformstatus = False
try:
    from sherpa.astro.io.wcs import WCS
    transformstatus = True
except:
    warning('failed to import WCS module; WCS routines will not be ' +
            'available')

__all__ = ('get_table_data', 'get_image_data', 'get_arf_data', 'get_rmf_data',
           'get_pha_data','set_table_data', 'set_image_data', 'set_pha_data',
           'get_column_data', 'get_ascii_data')

def _open_crate(type, args):
    crate = type(*args)
    if crate.get_status() == pycrates.dmFAILURE:
        raise IOErr('openfailed',crate.get_status_message().strip('ERROR').strip(' - ').strip())
    return crate

def _try_hdr_key(crate, name, dtype=str):
    name = name.upper().strip()
    if not crate.key_exists(name):
        return None

    return _try_key(crate.get_key(name), dtype)
    
def _try_key(key, dtype=str):
    if key is None:
        return None
    key = key.get_value()

    if str(key).find('none') != -1:
        return None

    return dtype( key )


def _get_meta_data(crate):
    meta = Meta()
    if crate.get_keynames is not None and crate.get_nkeys() > 0:
        for key in crate.get_keynames():
            val = crate.get_key_value(key)

            # empty numpy strings are not recognized by load pickle!
            if type(val) is numpy.str_ and val == '':
                val = ''

            meta[key] = val
    return meta


def _set_key(crate, key, name, val, fix_type=False, dtype=str):
    key.name = str(name).strip().upper()
    
    if fix_type:
        val = dtype(val)

    key.set_value( val )
    pycrates.add_key(crate, key)

def _set_column(crate, col, name, val):
    col.name = name.upper()
    col.set_nsets( len(val) )
    col.set_values( numpy.asarray(val) )    
    pycrates.add_col(crate, col)

def _set_pha_column(crate, col, name, val):
    col.name = name.upper()

    if col.load( numpy.asarray(val), True) != pycrates.dmSUCCESS:
        raise IOErr('setcolfailed', col.name, col.get_status_message())
    
    exec('crate.' + name + ' = col')

def _require_key(crate, name, key, dtype=str):
    key = _try_key(key, dtype)
    if key is None:
        raise IOErr('nokeyword', crate.get_filename(), name)
    return key

def _require_hdr_key(crate, name, dtype=str):
    key = _try_hdr_key(crate, name, dtype)
    if key is None:
        raise IOErr('nokeyword', crate.get_filename(), name)
    return key
   
def _try_unit(col, dtype=str):
    if col is None:
        return None

    unit = col.get_unit()

    if str(unit) == '':
        return None

    return str(unit)

def _try_col(col, make_copy=False, fix_type=False, dtype=SherpaFloat):
    if col is None:
        return None

    if make_copy:
        # Make a copy if a filename passed in
        data = numpy.array(col.get_values()).ravel()

    else:
        # Use a reference if a crate passed in
        data = numpy.asarray(col.get_values()).ravel()

    if fix_type:
        data = data.astype(dtype)
 
    return data

def _require_tbl_col(crate, colname, cnames, make_copy=False,
                     fix_type=False, dtype=SherpaFloat):

    name = str(colname).strip()
    
    if pycrates.col_exists(crate, name) == pycrates.dmFAILURE:
        raise IOErr('reqcol', name, cnames)

    data = pycrates.get_colvals(crate, name)

    if make_copy:
        # Make a copy if a filename passed in
        data = numpy.array(data)
    else:
        # Use a reference if a crate passed in
        data = numpy.asarray(data)
    
    if fix_type:
        data = data.astype(dtype)
 
    return numpy.column_stack(data)

def _try_key_list(key, num, dtype=SherpaFloat, fix_type=False):
    if key is None:
        return numpy.array([None]*num)
    
    # Make a copy of the data, since we don't know that pycrates will
    # do something sensible wrt reference counting
    key = numpy.array([key.get_value()]*num)
    
    if fix_type:
        key = key.astype(dtype)
        
    return key

def _try_col_list(col, num, make_copy=False, fix_type=False,
                  dtype=SherpaFloat):
    if col is None:
        return numpy.array([None]*num)
    
    if make_copy:
        # Make a copy if a filename passed in
        col = numpy.array(col.get_values())
    else:
        # Use a reference if a crate passed in
        col = numpy.asarray(col.get_values())
    
    if fix_type:
        col = col.astype(dtype)
        
    return col

def _require_col_list(col, num, make_copy=False, fix_type=False,
                      dtype=SherpaFloat):
    if col is None:
        raise IOErr('badcol', col.name)
    return _try_col_list(col, num, make_copy, fix_type, dtype)

def _require_col(col, make_copy=False, fix_type=False, dtype=SherpaFloat):
    data = _try_col(col, make_copy, fix_type, dtype)
    if data is None:
        raise IOErr('badcol', col.name)
    return data

def _try_image(crate, make_copy=False, fix_type=False, dtype=SherpaFloat):

    if make_copy:
        # Make a copy if a filename passed in
        dat = pycrates.copy_piximgvals(crate).squeeze()
    else:
        # Use a reference if a crate passed in
        dat = pycrates.get_piximgvals(crate).squeeze()

    if fix_type:
        dat = dat.astype(dtype)
        
    # FITS standard is FORTRAN filling, Sherpa is c-style
#    return dat.reshape(dat.shape[::-1]) 

    # Crates now returns image in c-style filling
    return dat


def _require_image(crate, make_copy=False, fix_type=False, dtype=SherpaFloat):
    dat = _try_image(crate, make_copy, fix_type, dtype)
    if dat is None:
        raise IOErr('badimg', crate.get_filename())
    return dat


## Read Functions ##

def get_header_data( arg, blockname=None, hdrkeys=None ):

    filename = ''
    if type(arg) == str and pycrates.Crate(arg).is_image()==0:

        filename = arg

        isbinary = True
        colnames = True
        dmsyn = ''

        if '[' in filename and ']' in filename:
            parts = filename.split('[')
            filename = parts.pop(0)
            if parts:
                dmsyn = parts.pop(0).lower()

        if not is_binary_file(filename):
            isbinary = False
            fd = open(filename, 'r')
            try:
                last=None
                line = fd.readline().strip()
                while len(line) > 0 and line[0] in '#%':
                    last = line
                    line = fd.readline().strip()
                if (last is not None and
                    (len(last.split(' ')) != len(line.split(' ')) )):
                    colnames = False
            finally:
                fd.close()

        if blockname is not None:
            arg += "[%s]" % str(blockname).upper()

        if (not isbinary) and (not colnames) and (not 'cols' in dmsyn):
            arg += "[opt colnames=none]"

        try:
            tbl = _open_crate(pycrates.TABLECrate, [arg])
        except Exception, e:
            try:
                tbl =  _open_crate(pycrates.IMAGECrate, [arg])
            except:
                raise e

        filename = tbl.get_filename()

        # Make a copy of the data, since we don't know that pycrates will
        # do something sensible wrt reference counting
    elif isinstance(arg, pycrates.TABLECrate):
        tbl = arg
        filename = arg.get_filename()
        make_copy=False
    else:
        raise IOErr('badfile', arg, 'TABLECrate obj')

    hdr = {}
    if hdrkeys is not None:
        for key in hdrkeys:
            hdr[key] = _require_hdr_key(tbl, key)
    else:
        for key in tbl.get_keynames():
           hdr[key] = _require_hdr_key(tbl, key) 

    return hdr


def get_column_data( *args ):
    """
    get_column_data( *NumPy_args )

    get_column_data( *CrateData_args )
    
    get_column_data( *NumPy_and_CrateData_args ) 
    """
    # args is passed as type list
    if len(args) == 0:
        raise IOErr('noarrays')

    cols = []
    for arg in args:
        if isinstance(arg, pycrates.CrateData):
            vals = arg.get_values()
        elif arg is None or type(arg) in (numpy.ndarray, list, tuple):
            vals = arg
        else:
            raise IOErr('badarray', arg)

        if arg is not None:
            vals = numpy.asarray( vals )
            for col in numpy.column_stack(vals):
                cols.append( col )
        else:
            cols.append( vals )

    return cols

def get_ascii_data(filename, ncols=2, colkeys=None, **kwargs):
    """
    get_table_data( filename [, ncols=2 [, colkeys=None [, **kwargs ]]] ) 
    """
    return get_table_data( filename, ncols, colkeys )[:3]


def get_table_data( arg, ncols=1, colkeys=None, make_copy=True, fix_type=True,
                    blockname=None, hdrkeys=None):
    """
    get_table_data( filename , ncols=1 [, colkeys=None [, make_copy=True [,
                    fix_type=True [, blockname=None [, hdrkeys=None ]]]]])

    get_table_data( TABLECrate , ncols=1 [, colkeys=None [, make_copy=True [,
                    fix_type=True [, blockname=None [, hdrkeys=None ] ]]]])
    """
    filename = ''
    if type(arg) == str and pycrates.Crate(arg).is_image()==0:

        filename = arg

        isbinary = True
        colnames = True
        dmsyn = ''

        if '[' in filename and ']' in filename:
            parts = filename.split('[')
            filename = parts.pop(0)
            if parts:
                dmsyn = parts.pop(0).lower()

        if not is_binary_file(filename):
            isbinary = False
            fd = open(filename, 'r')
            try:
                last=None
                line = fd.readline().strip()
                while len(line) > 0 and line[0] in '#%':
                    last = line
                    line = fd.readline().strip()
                if (last is not None and
                    (len(last.split(' ')) != len(line.split(' ')) )):
                    colnames = False
            finally:
                fd.close()

        if blockname is not None:
            arg += "[%s]" % str(blockname).upper()

        if (not isbinary) and (not colnames) and (not 'cols' in dmsyn):
            arg += "[opt colnames=none]"

        tbl = _open_crate(pycrates.TABLECrate, [arg])

        filename = tbl.get_filename()

        # Make a copy of the data, since we don't know that pycrates will
        # do something sensible wrt reference counting
    elif isinstance(arg, pycrates.TABLECrate):
        tbl = arg
        filename = arg.get_filename()
        make_copy=False
    else:
        raise IOErr('badfile', arg, 'TABLECrate obj')

    cnames = list(pycrates.get_col_names(tbl, vectors=False, rawonly=True))

    if colkeys is not None:
        colkeys = [str(name).strip() for name in list(colkeys)]

    elif (type(arg) == str and (not os.path.isfile(arg))
          and '[' in arg and ']' in arg):
        colkeys = cnames

    # Try Channel, Counts or X,Y before defaulting to first two table cols
    elif 'CHANNEL' in cnames and 'COUNTS' in cnames:
        colkeys = ['CHANNEL','COUNTS']

    elif 'X' in cnames and 'Y' in cnames:
        colkeys = ['X','Y']

    else:
        colkeys = cnames[:ncols]

    cols = []
    for name in colkeys:
        for col in _require_tbl_col(tbl, name, cnames, make_copy, fix_type):
            cols.append(col)

    hdr={}
    if hdrkeys is not None:
        for key in hdrkeys:
            hdr[key] = _require_hdr_key(tbl, key)

    return colkeys, cols, filename, hdr


def get_image_data(arg, make_copy=True, fix_type=True):
    """
    get_image_data ( filename [, make_copy=True, fix_type=True ])

    get_image_data ( IMAGECrate [, make_copy=True, fix_type=True ])
    """
    filename = ''
    if type(arg) == str and pycrates.Crate(arg).is_image()==1:
        #img = pycrates.read_file(arg)
        img = _open_crate(pycrates.IMAGECrate, [arg])
        filename = arg

    elif isinstance(arg, pycrates.IMAGECrate):
        img = arg
        filename = arg.get_filename()
        make_copy=False
    else:
        raise IOErr('badfile', arg, "IMAGECrate obj")

    data = {}

    data['y'] = _require_image(img, make_copy, fix_type)

    sky = None
    if img.get_axis('sky') is not None:
        sky = pycrates.get_transform(img, 'sky')

    elif img.get_axis('SKY') is not None:
        sky = pycrates.get_transform(img, 'SKY')

    elif img.get_axis('pos') is not None:
        sky = pycrates.get_transform(img, 'pos')
        
    elif img.get_axis('POS') is not None:
        sky = pycrates.get_transform(img, 'POS')
    
    wcs = None
    if img.get_axis('EQPOS') is not None:
        wcs = pycrates.get_transform(img, 'EQPOS')

    if sky is not None and transformstatus:
        linear = pycrates.WCSTANTransform()
        linear.set_flavor("LINEAR")
        linear.set_transform_matrix(sky.get_transform_matrix())
        cdelt = numpy.array(linear.get_parameter_value('CDELT'))
        crpix = numpy.array(linear.get_parameter_value('CRPIX'))
        crval = numpy.array(linear.get_parameter_value('CRVAL'))
        data['sky'] = WCS('physical', 'LINEAR', crval, crpix, cdelt)

    if wcs is not None and transformstatus:
        cdelt = numpy.array(wcs.get_parameter_value('CDELT'))
        crpix = numpy.array(wcs.get_parameter_value('CRPIX'))
        crval = numpy.array(wcs.get_parameter_value('CRVAL'))
        crota = SherpaFloat(wcs.get_parameter_value('CROTA'))
        equin = SherpaFloat(wcs.get_parameter_value('EQUINOX'))
        epoch = SherpaFloat(wcs.get_parameter_value('EPOCH'))
        data['eqpos'] = WCS('world', 'WCS', crval, crpix, cdelt,
                            crota, epoch, equin)

    data['header'] = _get_meta_data(img)

    keys = ['MTYPE1','MFORM1','CTYPE1P','CTYPE2P','WCSNAMEP','CDELT1P',
            'CDELT2P','CRPIX1P','CRPIX2P','CRVAL1P','CRVAL2P',
            'MTYPE2','MFORM2','CTYPE1','CTYPE2','CDELT1','CDELT2','CRPIX1',
            'CRPIX2','CRVAL1','CRVAL2','CUNIT1','CUNIT2','EQUINOX']
#            'WCSTY1P', 'WCSTY2P']
    
    for key in keys:
        try:
            data['header'].pop(key)
        except KeyError:
            pass

    return data, filename


def get_arf_data(arg, make_copy=True):
    """
    get_arf_data( filename [, make_copy=True ])

    get_arf_data( ARFCrate [, make_copy=True ])
    """
    filename = ''
    if type(arg) == str:
        #arf = pycrates.read_arf(arg)
        arf = _open_crate(pycrates.ARFCrate, [arg])
        filename = arg

        # Make a copy of the data, since we don't know that pycrates will
        # do something sensible wrt reference counting
    elif pycrates.is_arf(arg) == pycrates.dmSUCCESS:
        arf = arg
        filename = arg.get_filename()
        make_copy=False
    else:
        raise IOErr('badfile', arg, "ARFCrate obj")

    if arf is None or arf.get_colnames() is None:
        raise IOErr('filenotfound', arg)

    data = {}

    if arf.energ_lo is None:
        raise IOErr('reqcol', 'ENERG_LO', filename)
    if arf.energ_hi is None:
        raise IOErr('reqcol', 'ENERG_HI', filename)
    if arf.specresp is None:
        raise IOErr('reqcol', 'SPECRESP', filename)

    data['energ_lo'] = _require_col(arf.energ_lo, make_copy, fix_type=True)
    data['energ_hi'] = _require_col(arf.energ_hi, make_copy, fix_type=True)
    data['specresp'] = _require_col(arf.specresp, make_copy,fix_type=True)
    data['exposure'] = _try_key(arf.exposure, dtype=SherpaFloat)
    data['bin_lo']   = _try_col(arf.bin_lo, make_copy, fix_type=True)
    data['bin_hi']   = _try_col(arf.bin_hi, make_copy, fix_type=True)
    data['header']     = _get_meta_data(arf)
    data['header'].pop('EXPOSURE')

    return data, filename


def get_rmf_data(arg, make_copy=True):
    """
    get_rmf_data( filename [, make_copy=True ])

    get_rmf_data( RMFCrate [, make_copy=True ])
    """
    filename = ''
    if type(arg) == str:
        #rmf = pycrates.read_rmf(arg)
        rmf = _open_crate(pycrates.RMFCrate, [arg])
        filename = arg

        # Make a copy of the data, since we don't know that pycrates will
        # do something sensible wrt reference counting
    elif pycrates.is_rmf(arg) == pycrates.dmSUCCESS:
        rmf = arg
        filename = arg.get_filename() 
        make_copy = False
    else:
        raise IOErr('badfile', arg, "RMFCrate obj")

    if rmf is None or rmf.get_colnames() is None:
        raise IOErr('filenotfound', arg)

    data = {}

    if rmf.energ_lo is None:
        raise IOErr('reqcol', 'ENERG_LO', filename)
    if rmf.energ_hi is None:
        raise IOErr('reqcol', 'ENERG_HI', filename)
    if rmf.matrix is None:
        raise IOErr('reqcol', 'MATRIX', filename)
    if rmf.n_grp is None:
        raise IOErr('reqcol', 'N_GRP', filename)
    if rmf.f_chan is None:
        raise IOErr('reqcol', 'F_CHAN', filename)
    if rmf.n_chan is None:
        raise IOErr('reqcol', 'N_CHAN', filename)

    data['detchans'] = _require_hdr_key(rmf, 'DETCHANS', SherpaInt)
    data['energ_lo'] = _require_col(rmf.energ_lo, make_copy, fix_type=True)
    data['energ_hi'] = _require_col(rmf.energ_hi, make_copy, fix_type=True)
    data['n_grp'] = _require_col(rmf.n_grp, make_copy,
                                 dtype=SherpaUInt, fix_type=True)
    fcbuf = _require_col(rmf.f_chan, make_copy)
    ncbuf = _require_col(rmf.n_chan, make_copy)
    respbuf = _require_col_list(rmf.matrix, 1, make_copy)
    data['e_min']    = _try_col(rmf.e_min, make_copy, fix_type=True)
    data['e_max']    = _try_col(rmf.e_max, make_copy, fix_type=True)
    data['header']     = _get_meta_data(rmf)
    data['header'].pop('DETCHANS')

    offset = rmf.f_chan.get_tlmin()

    if offset < 0:
        error("Failed to locate TLMIN keyword for F_CHAN" +
              " column in RMF file '%s'; "  % filename +
              'Update the offset value in the RMF data set to' +
              ' appropriate TLMIN value prior to fitting')

    if offset < 0 and rmf.channel is not None:
        offset = rmf.channel.get_tlmin()
        
    # If response is non-OGIP, tlmin is -(max of type), so resort to default
    if not (offset < 0):
        data['offset'] = offset

    #
    # FIXME:
    #
    # Currently, CRATES does something screwy:  If n_grp is zero in a bin,
    # it appends a zero to f_chan, n_chan, and matrix.  I have no idea what
    # the logic behind this is -- why would you add data that you know you
    # don't need?  Although it's easy enough to filter the zeros out of
    # f_chan and n_chan, it's harder for matrix, since zero is a legitimate
    # value there.
    #
    # I think this crazy behavior of CRATES should be changed, but for the
    # moment we'll just punt in this case.  (If we don't, the calculation
    # in rmf_fold() will be trashed.)


    # CRATES does not support variable length arrays, so here we condense
    # the array of tuples into the proper length array

    chan_width = data['n_grp'].max()
    resp_width = 0
    if len(respbuf.shape) > 1:
        resp_width = respbuf.shape[1]
    
    (data['f_chan'], data['n_chan'],
     data['matrix'] ) = resp_init( data['n_grp'], fcbuf, ncbuf,
                                   chan_width, respbuf.ravel(), resp_width )

    return data, filename


def get_pha_data(arg, make_copy=True, use_background=False):
    """
    get_pha_data( filename [, make_copy=True [, use_background=False]])

    get_pha_data( PHACrate [, make_copy=True [, use_background=False]])
    """
    filename = ''
    if type(arg) == str:
        #pha = pycrates.read_pha(arg, use_background)
        pha = _open_crate(pycrates.PHACrate, [arg, use_background])

        filename = arg

        # Make a copy of the data, since we don't know that pycrates will
        # do something sensible wrt reference counting
    elif pycrates.is_pha(arg) == pycrates.dmSUCCESS:
        pha = arg
        filename = arg.get_filename()
        make_copy=False
    else:
        raise IOErr('badfile', arg, "PHACrate obj")

    if pha is None or pha.get_colnames() is None:
        raise IOErr('filenotfound', arg)

    keys = ['BACKFILE','ANCRFILE','RESPFILE',
            'BACKSCAL','AREASCAL','EXPOSURE']

    datasets = []
    if (pha.pha1_type_flag or
        (not pha.pha1_type_flag and len(pha.spec_num.get_values())==1)):
        data = {}

        # Keywords
        data['exposure'] = _try_key(pha.exposure, SherpaFloat)
        #data['poisserr'] = _try_key(pha.poisserr, bool)
        data['backfile'] = _try_key(pha.backfile_key)
        data['arffile']  = _try_key(pha.ancrfile_key)
        data['rmffile']  = _try_key(pha.respfile_key)

        # Keywords or columns
        data['backscal'] = _try_key(pha.backscal_key, SherpaFloat)
        if data['backscal'] is None:
            data['backscal'] = _try_col(pha.backscal, make_copy)

        data['backscup'] = _try_key(pha.backscup_key, SherpaFloat)
        if data['backscup'] is None:
            data['backscup'] = _try_col(pha.backscup, make_copy)

        data['backscdn'] = _try_key(pha.backscdn_key, SherpaFloat)
        if data['backscdn'] is None:
            data['backscdn'] = _try_col(pha.backscdn, make_copy)

        data['areascal'] = _try_key(pha.areascal_key, SherpaFloat)
        if data['areascal'] is None:
            data['areascal'] = _try_col(pha.areascal, make_copy)

        data['header']            = _get_meta_data(pha)

        for key in keys:
            try:
                data['header'].pop(key)
            except KeyError:
                pass

        # Columns

        if pha.channel is None:
            raise IOErr('reqcol', 'CHANNEL', filename)

        data['channel']         = _require_col(pha.channel, make_copy, fix_type=True)
        # Make sure channel numbers, not indices
        if pha.channel.get_tlmin() == 0:
            data['channel'] = data['channel']+1
        data['counts'] = None
        if pha.counts is not None:
            data['counts'] = _require_col(pha.counts, make_copy, fix_type=True)
        else:
            if pha.rate is None:
                raise IOErr('reqcol', 'COUNTS or RATE', filename)
            data['counts'] = _require_col(pha.rate, make_copy, fix_type=True)*data['exposure']
        data['staterror']       = _try_col(pha.stat_err, make_copy)
        data['syserror']        = _try_col(pha.sys_err, make_copy)
        data['background_up']   = _try_col(pha.background_up, make_copy, fix_type=True)
        data['background_down'] = _try_col(pha.background_down, make_copy, fix_type=True)
        data['bin_lo']          = _try_col(pha.bin_lo, make_copy,fix_type=True)
        data['bin_hi']          = _try_col(pha.bin_hi, make_copy,fix_type=True)
        data['grouping']        = _try_col(pha.grouping, make_copy)
        data['quality']         = _try_col(pha.quality, make_copy)

        datasets.append(data)

    else:
        # Type 2 PHA file support
        data = {}
        num = pha.spec_num.get_nsets()

        # Keywords
        exposure = _try_key(pha.exposure, SherpaFloat)
        #poisserr = _try_key(pha.poisserr, bool)
        backfile = _try_key(pha.backfile_key)
        arffile  = _try_key(pha.ancrfile_key)
        rmffile  = _try_key(pha.respfile_key)

        # Keywords or columns
        backscal = _try_key_list(pha.backscal_key, num)
        if pha.backscal_key is None:
            backscal = _try_col_list(pha.backscal, num, make_copy)

        backscup = _try_key_list(pha.backscup_key, num)
        if pha.backscup_key is None:
            backscup = _try_col_list(pha.backscup, num, make_copy)

        backscdn = _try_key_list(pha.backscdn_key, num)
        if pha.backscdn_key is None:
            backscdn = _try_col_list(pha.backscdn, num, make_copy)

        areascal = _try_key_list(pha.areascal_key, num)
        if pha.areascal_key is None:
            areascal = _try_col_list(pha.areascal, num, make_copy)

        # Columns

        if pha.channel is None:
            raise IOErr('reqcol', 'CHANNEL', filename)

        channel         = _require_col_list(pha.channel, num, make_copy, fix_type=True)
        counts = None
        if pha.counts is not None:
            counts      = _require_col_list(pha.counts, num, make_copy, fix_type=True)
        else:
            if pha.rate is None:
                raise IOErr('reqcol', 'COUNTS or RATE', filename)
            counts      = _require_col_list(pha.rate, num, make_copy, fix_type=True) * exposure
        staterror       = _try_col_list(pha.stat_err, num, make_copy)
        syserror        = _try_col_list(pha.sys_err, num, make_copy)
        background_up   = _try_col_list(pha.background_up, num, make_copy, fix_type=True)
        background_down = _try_col_list(pha.background_down, num, make_copy, fix_type=True)
        bin_lo          = _try_col_list(pha.bin_lo, num, make_copy,
                                        fix_type=True)
        bin_hi          = _try_col_list(pha.bin_hi, num, make_copy,
                                        fix_type=True)
        grouping        = _try_col_list(pha.grouping, num, make_copy)
        quality         = _try_col_list(pha.quality, num, make_copy) 

        orders          = _try_key_list(pha.tg_m_key,num)
        if pha.tg_m_key is None:
            orders      = _try_col_list(pha.tg_m, num, make_copy)

        parts           = _try_key_list(pha.tg_part_key,num)
        if pha.tg_part_key is None:
            parts       = _try_col_list(pha.tg_part, num, make_copy)

        specnums        = _try_col_list(pha.spec_num, num, make_copy)
        srcids          = _try_col_list(pha.tg_srcid, num, make_copy)

        # Iterate over all rows of channels, counts, errors, etc
        # Populate a list of dictionaries containing individual dataset info
        for (bscal, bscup, bscdn, arsc, chan, cnt, staterr, syserr,
             backup, backdown, binlo, binhi, grp, qual, ordr, prt,
             specnum, srcid
             ) in izip(backscal, backscup, backscdn, areascal, channel,
                       counts, staterror, syserror, background_up,
                       background_down, bin_lo, bin_hi, grouping, quality,
                       orders, parts, specnums, srcids):

            data = {}

            data['exposure'] = exposure
            #data['poisserr'] = poisserr
            data['backfile'] = backfile
            data['arffile']  = arffile
            data['rmffile']  = rmffile

            data['backscal'] = bscal
            data['backscup'] = bscup
            data['backscdn'] = bscdn
            data['areascal'] = arsc

            data['channel']         = chan
            data['counts']          = cnt
            data['staterror']       = staterr
            data['syserror']        = syserr
            data['background_up']   = backup
            data['background_down'] = backdown
            data['bin_lo']          = binlo
            data['bin_hi']          = binhi
            data['grouping']        = grp
            data['quality']         = qual
            data['header']            = _get_meta_data(pha)
            data['header']['TG_M']     = ordr
            data['header']['TG_PART']  = prt
            data['header']['SPEC_NUM'] = specnum
            data['header']['TG_SRCID'] = srcid

            for key in keys:
                try:
                    data['header'].pop(key)
                except KeyError:
                    pass

            datasets.append(data)

    return datasets, filename


#
### Write/Pack Functions ####
#


def set_image_data(filename, data, header, ascii=False, clobber=False,
                   packup=False):

    if not packup and os.path.isfile(filename) and not clobber:
        raise IOErr('filefound', filename)

    img = pycrates.IMAGECrate()

    # Write Image Header Keys
    for key in header.keys():
        if header[key] is None:
            continue
        _set_key(img, pycrates.CrateKey(), key, header[key])

    # Write Image WCS Header Keys
    if data['eqpos'] is not None:
        cdeltw = data['eqpos'].cdelt
        crvalw = data['eqpos'].crval
        crpixw = data['eqpos'].crpix
        equin  = data['eqpos'].equinox

    if data['sky'] is not None:
        cdeltp = data['sky'].cdelt
        crvalp = data['sky'].crval
        crpixp = data['sky'].crpix

        _set_key(img, pycrates.CrateKey(), 'MTYPE1', 'sky     ')
        _set_key(img, pycrates.CrateKey(), 'MFORM1', 'x,y     ')
        _set_key(img, pycrates.CrateKey(), 'CTYPE1P','x       ')
        _set_key(img, pycrates.CrateKey(), 'CTYPE2P','y       ')
        _set_key(img, pycrates.CrateKey(), 'WCSNAMEP','PHYSICAL')
        _set_key(img, pycrates.CrateKey(), 'CDELT1P', cdeltp[0])
        _set_key(img, pycrates.CrateKey(), 'CDELT2P', cdeltp[1])
        _set_key(img, pycrates.CrateKey(), 'CRPIX1P', crpixp[0])
        _set_key(img, pycrates.CrateKey(), 'CRPIX2P', crpixp[1])
        _set_key(img, pycrates.CrateKey(), 'CRVAL1P', crvalp[0])
        _set_key(img, pycrates.CrateKey(), 'CRVAL2P', crvalp[1])

        if data['eqpos'] is not None:
            # Simply the inverse of read transformations in get_image_data
            cdeltw = cdeltw * cdeltp
            crpixw = ((crpixw - crvalp) /  cdeltp + crpixp )

    if data['eqpos'] is not None:
        _set_key(img, pycrates.CrateKey(), 'MTYPE2', 'EQPOS   ')
        _set_key(img, pycrates.CrateKey(), 'MFORM2', 'RA,DEC  ')
        _set_key(img, pycrates.CrateKey(), 'CTYPE1', 'RA---TAN')
        _set_key(img, pycrates.CrateKey(), 'CTYPE2', 'DEC--TAN')
        _set_key(img, pycrates.CrateKey(), 'CDELT1', cdeltw[0])
        _set_key(img, pycrates.CrateKey(), 'CDELT2', cdeltw[1])
        _set_key(img, pycrates.CrateKey(), 'CRPIX1', crpixw[0])
        _set_key(img, pycrates.CrateKey(), 'CRPIX2', crpixw[1])
        _set_key(img, pycrates.CrateKey(), 'CRVAL1', crvalw[0])
        _set_key(img, pycrates.CrateKey(), 'CRVAL2', crvalw[1])
        _set_key(img, pycrates.CrateKey(), 'EQUINOX', equin)

    # Write Image pixel values
    shape = data['pixels'].shape
    pix_col = pycrates.CrateData()
    pix_col.set_dimarr( shape, len(shape) )
    pix_col.set_nsets(1)
    pix_col.set_values( data['pixels'].ravel() )

    pycrates.add_piximg(img, pix_col)

    if packup:
        return img

    if ascii and '[' not in filename and ']' not in filename:
        #filename += "[opt kernel=text/simple]"
        raise IOErr('writenoimg')

    pycrates.write_file(img, filename)


def set_table_data(filename, data, col_names, hdr=None, hdrnames=None,
                   ascii=False, clobber=False, packup=False):

    if not packup and os.path.isfile(filename) and not clobber:
        raise IOErr("filefound", filename)

    tbl = pycrates.TABLECrate()

    col_names = [name for name in col_names if data[name] != None]
    col_names.remove('name')    
    try:
        cols = [pycrates.CrateData() for i in range(len(col_names))]
        for name, col in izip(col_names, cols):
            if data[name] is None:
                continue
            _set_column( tbl, col, name, data[name] )

    finally:
        if packup:
            return tbl

        if ascii and '[' not in filename and ']' not in filename:
            filename += "[opt kernel=text/simple]"

        pycrates.write_file(tbl, filename)


def set_pha_data(filename, data, col_names, header=None,
                 ascii=False, clobber=False, packup=False):

    if not packup and os.path.isfile(filename) and not clobber:
        raise IOErr("filefound", filename)

    pha = pycrates.PHACrate()

    cols = [pycrates.CrateData() for i in range(len(col_names))]

    try:

        # Write header values using CrateKey objects
        for key in header.keys():
            if header[key] is None:
                continue
            _set_key( pha, pycrates.CrateKey(), key, header[key] )

        # Write column values using CrateData objects
        for name, col in izip(col_names, cols):
            if data[name] is None:
                continue
            _set_pha_column( pha, col, name, data[name] )

    finally:
        if packup:
            return pha

        if ascii and '[' not in filename and ']' not in filename:
            filename += "[opt kernel=text/simple]"

        pycrates.write_pha(pha, filename)


def set_arrays(filename, args, fields=None, ascii=True, clobber=False):

    if os.path.isfile(filename) and not clobber:
        raise IOErr("filefound", filename)

    if not numpy.iterable(args) or len(args) == 0:
        raise IOErr('noarrayswrite')

    if not numpy.iterable(args[0]):
        raise IOErr('noarrayswrite')

    size = len(args[0])
    for arg in args:
        if not numpy.iterable(arg):
            raise IOErr('noarrayswrite')
        elif len(arg) != size:
            raise IOErr('arraysnoteq')

    if ascii and '[' not in filename and ']' not in filename:
        filename += "[opt kernel=text/simple]"

    tbl = pycrates.TABLECrate()

    cols = [pycrates.CrateData() for ii in range(len(args))]

    if fields is None:
        fields = ['col%i' % (ii+1) for ii in range(len(args))]

    if len(args) != len(fields):
        raise IOErr('toomanycols', str(len(fields)), str(len(args)))

    for col, val, name in izip(cols, args, fields):
        _set_column(tbl, col, name, val)

    pycrates.write_file(tbl, filename)
