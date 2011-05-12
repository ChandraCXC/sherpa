#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
import numpy
from sherpa.utils import SherpaFloat, get_num_args, is_binary_file
from sherpa.utils.err import IOErr
from sherpa.data import Data1D, BaseData
from itertools import izip
from exceptions import ValueError
import os


_all__ = ('read_data', 'write_data', 'get_ascii_data', 'read_arrays',
          'write_arrays')


def _is_subclass(t1, t2):
    return isinstance(t1, type) and issubclass(t1, t2) and (t1 is not t2)

def _check_args( size, dstype ):
    # Find the number of required args minus self, filename
    req_args = (get_num_args( dstype.__init__ )[1] - 2)

    if size < req_args:
        #raise IOErr('badargs', dstype.__name__, req_args)
        raise TypeError("data set '%s' takes at least %s args" % 
                        (dstype.__name__, req_args))

def read_file_data(filename, sep=' ', comment='#'):

    bad_chars = '\t\n\r,;: |'
    fp = open(filename, 'r')
    raw_names = []
    rows = []
    try:
        for line in fp:
            for char in bad_chars:
                if char in line:
                    # replace any bad chars in line with sep for tokenize
                    line = line.replace(char, sep)

            line = line.strip()

            #look for last commented line before data 
            if len(line) > 0 and line[0] == comment:
                line = line.replace(comment, ' ') #slice off comment
                raw_names = line.strip().split(sep)
                
            elif line == '':
                continue
            else:
                # split line at sep
                elems = line.strip().split(sep)
                row = [elem for elem in elems if elem != '']
                
                # make list of row elements
                rows.append(row)

    finally:
        fp.close()

    # rotate rows into list of columns
    cols = numpy.column_stack(rows)

    # cast columns to appropriate type
    args = []
    for col in cols:
        try:
            args.append(col.astype(SherpaFloat))
        except ValueError:
            args.append(col)

    names = [name.strip(bad_chars) for name in raw_names if name != '']
  
    if len(names) == 0:
        names = ['col%i' % (i+1) for i in xrange(len(args))]

    return names, args

def get_column_data( *args ):
   """
   get_column_data( *NumPy_args )
   """
   if len(args) == 0:
       raise IOErr('noarrays')
   
   cols = []
   for arg in args:
       if arg is None or type(arg) in (numpy.ndarray, list, tuple):
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


def get_ascii_data(filename, ncols=1, colkeys=None, sep=' ', dstype=Data1D,
                   comment='#'):

    if is_binary_file(filename):
        raise IOErr('notascii', filename)

    names, args = read_file_data(filename, sep, comment)

    if colkeys is None:
        kwargs = []
        if ncols != 1:
            _check_args(ncols, dstype)
        kwargs.extend(args[:ncols])
        return (names, kwargs, filename)

    kwargs = []
    colkeys = list(colkeys)

    if len(names) > len(args):
        raise IOErr('toomanycols')

    assert( len(names) <= len(args) )

    for key in colkeys:
        if key not in names:
            raise IOErr('reqcol', key, numpy.asarray(names, numpy.string_))
        kwargs.append(args[names.index(key)])

    _check_args(len(kwargs), dstype)
    return (colkeys, kwargs, filename)


def read_data( filename, ncols=2, colkeys=None, sep=' ', dstype=Data1D,
               comment='#'):

    colnames, args, name = get_ascii_data(filename, ncols, colkeys,
                                          sep, dstype, comment)
    return dstype(name, *args)


def read_arrays(*args): 
    """
    read_table( *NumPy_args [, dstype = Data1D ] )
    """
    args = list(args)
    if len(args) == 0:
        raise IOErr('noarrays')

    dstype=Data1D
    if _is_subclass(args[-1], BaseData):
        dstype = args.pop()

    args = get_column_data(*args)

    # Determine max number of args for dataset constructor
    _check_args( len(args), dstype )

    return dstype('', *args)


def write_arrays(filename, args, fields=None, sep=' ', comment='#',
                 clobber=False, linebreak='\n', format='%g'):

    if os.path.isfile(filename) and not clobber:
        raise IOErr("filefound",filename)

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

    args = numpy.column_stack(numpy.asarray(args))

    f = file(filename, 'w')

    if fields is not None:
        f.write(comment + sep.join(fields) + linebreak)

    lines = []
    for arg in args:
        line = [format % elem for elem in arg]
        lines.append(sep.join(line))

    f.write(linebreak.join(lines))

    # add a newline at end
    f.write(linebreak)
    f.close()


def write_data(filename, dataset, fields=None, sep=' ', comment='#',
               clobber=False, linebreak='\n', format='%g'):

    if fields is None:
        fields = dataset._fields

    cols = []
    col_names = []

    for name in fields:
        field = getattr(dataset, name, None)
        if field is not None and name != 'name':
            col_names.append(name.upper())
            cols.append(field)

    write_arrays(filename, cols, col_names, sep, comment, clobber,
                 linebreak, format)
