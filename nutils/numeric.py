# -*- coding: utf8 -*-
#
# Module NUMERIC
#
# Part of Nutils: open source numerical utilities for Python. Jointly developed
# by HvZ Computational Engineering, TU/e Multiscale Engineering Fluid Dynamics,
# and others. More info at http://nutils.org <info@nutils.org>. (c) 2014

"""
The numeric module provides methods that are lacking from the numpy module.
"""

from __future__ import print_function, division
import numpy

_abc = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ' # indices for einsum

def grid( shape ):
  shape = tuple(shape)
  grid = numpy.empty( (len(shape),)+shape, dtype=int )
  for i, sh in enumerate( shape ):
    grid[i] = numpy.arange(sh)[(slice(None),)+(numpy.newaxis,)*(len(shape)-i-1)]
  return grid

def round( arr ):
  arr = numpy.asarray( arr )
  return arr if arr.dtype == int \
    else ( arr - numpy.less(arr,0) + .5 ).astype( int )

def sign( arr ):
  arr = numpy.asarray( arr )
  return (arr>=0).astype(int) - (arr<=0).astype(int)

def floor( arr ):
  ass = numpy.asarray( arr )
  return arr if arr.dtype == int \
    else ( arr - numpy.less(arr,0) ).astype( int )

def overlapping( arr, axis=-1, n=2 ):
  'reinterpret data with overlaps'

  arr = numpy.asarray( arr )
  if axis < 0:
    axis += arr.ndim
  assert 0 <= axis < arr.ndim
  shape = arr.shape[:axis] + (arr.shape[axis]-n+1,n) + arr.shape[axis+1:]
  strides = arr.strides[:axis] + (arr.strides[axis],arr.strides[axis]) + arr.strides[axis+1:]
  overlapping = numpy.lib.stride_tricks.as_strided( arr, shape, strides )
  overlapping.flags.writeable = False
  return overlapping

def normdim( ndim, n ):
  'check bounds and make positive'

  assert isint(ndim) and ndim >= 0, 'ndim must be positive integer, got %s' % ndim
  if n < 0:
    n += ndim
  assert 0 <= n < ndim, 'argument out of bounds: %s not in [0,%s)' % (n,ndim)
  return n

def align( arr, trans, ndim ):
  '''create new array of ndim from arr with axes moved accordin
  to trans'''

  arr = numpy.asarray( arr )
  assert arr.ndim == len(trans)
  if not len(trans):
    return arr[(numpy.newaxis,)*ndim]
  transpose = numpy.empty( ndim, dtype=int )
  trans = numpy.asarray( trans )
  nnew = ndim - len(trans)
  if nnew > 0:
    remaining = numpy.ones( ndim, dtype=bool )
    remaining[trans] = False
    inew, = remaining.nonzero()
    trans = numpy.hstack([ inew, trans ])
    arr = arr[(numpy.newaxis,)*nnew]
  transpose[trans] = numpy.arange(ndim)
  return arr.transpose( transpose )

def get( arr, axis, item ):
  'take single item from array axis'

  arr = numpy.asarray( arr )
  axis = normdim( arr.ndim, axis )
  return arr[ (slice(None),) * axis + (item,) ]

def expand( arr, *shape ):
  'expand'

  newshape = list( arr.shape )
  for i, sh in enumerate( shape ):
    if sh != None:
      if newshape[i-len(shape)] == 1:
        newshape[i-len(shape)] = sh
      else:
        assert newshape[i-len(shape)] == sh
  expanded = numpy.empty( newshape )
  expanded[:] = arr
  return expanded

def linspace2d( start, stop, steps ):
  'linspace & meshgrid combined'

  start = tuple(start) if isinstance(start,(list,tuple)) else (start,start)
  stop  = tuple(stop ) if isinstance(stop, (list,tuple)) else (stop, stop )
  steps = tuple(steps) if isinstance(steps,(list,tuple)) else (steps,steps)
  assert len(start) == len(stop) == len(steps) == 2
  values = numpy.empty( (2,)+steps )
  values[0] = numpy.linspace( start[0], stop[0], steps[0] )[:,numpy.newaxis]
  values[1] = numpy.linspace( start[1], stop[1], steps[1] )[numpy.newaxis,:]
  return values

def contract( A, B, axis=-1 ):
  'contract'

  A = numpy.asarray( A )
  B = numpy.asarray( B )

  maxdim = max( A.ndim, B.ndim )
  m = _abc[maxdim-A.ndim:maxdim]
  n = _abc[maxdim-B.ndim:maxdim]

  axes = sorted( [ normdim(maxdim,axis) ] if isinstance(axis,int) else [ normdim(maxdim,ax) for ax in axis ] )
  o = _abc[:maxdim-len(axes)] if axes == range( maxdim-len(axes), maxdim ) \
    else ''.join( _abc[a+1:b] for a, b in zip( [-1]+axes, axes+[maxdim] ) if a+1 != b )

  return numpy.einsum( '%s,%s->%s' % (m,n,o), A, B )

def contract_fast( A, B, naxes ):
  'contract last n axes'

  assert naxes >= 0
  A = numpy.asarray( A )
  B = numpy.asarray( B )

  maxdim = max( A.ndim, B.ndim )
  m = _abc[maxdim-A.ndim:maxdim]
  n = _abc[maxdim-B.ndim:maxdim]
  o = _abc[:maxdim-naxes]

  return numpy.einsum( '%s,%s->%s' % (m,n,o), A, B )

def dot( A, B, axis=-1 ):
  '''Transform axis of A by contraction with first axis of B and inserting
     remaining axes. Note: with default axis=-1 this leads to multiplication of
     vectors and matrices following linear algebra conventions.'''

  A = numpy.asarray( A )
  B = numpy.asarray( B )

  m = _abc[:A.ndim]
  x = _abc[A.ndim:A.ndim+B.ndim-1]
  n = m[axis] + x
  o = m[:axis] + x
  if axis != -1:
    o += m[axis+1:]

  return numpy.einsum( '%s,%s->%s' % (m,n,o), A, B )

def fastrepeat( A, nrepeat, axis=-1 ):
  'repeat axis by 0stride'

  A = numpy.asarray( A )
  assert A.shape[axis] == 1
  shape = list( A.shape )
  shape[axis] = nrepeat
  strides = list( A.strides )
  strides[axis] = 0
  return numpy.lib.stride_tricks.as_strided( A, shape, strides )

def fastmeshgrid( X, Y ):
  'mesh grid based on fastrepeat'

  return fastrepeat(X[numpy.newaxis,:],len(Y),axis=0), fastrepeat(Y[:,numpy.newaxis],len(X),axis=1)

def meshgrid( *args ):
  'multi-dimensional meshgrid generalisation'

  args = [ numpy.asarray(arg) for arg in args ]
  shape = [ len(args) ] + [ arg.size for arg in args if arg.ndim ]
  dtype = int if all( isintarray(a) for a in args ) else float
  grid = numpy.empty( shape, dtype=dtype )
  n = len(shape)-1
  for i, arg in enumerate( args ):
    if arg.ndim:
      n -= 1
      grid[i] = arg[(slice(None),)+(numpy.newaxis,)*n]
    else:
      grid[i] = arg
  assert n == 0
  return grid

def appendaxes( A, shape ):
  'append axes by 0stride'

  shape = (shape,) if isinstance(shape,int) else tuple(shape)
  A = numpy.asarray( A )
  return numpy.lib.stride_tricks.as_strided( A, A.shape + shape, A.strides + (0,)*len(shape) )

def takediag( A ):
  diag = A[...,0] if A.shape[-1] == 1 \
    else A[...,0,:] if A.shape[-2] == 1 \
    else numpy.einsum( '...ii->...i', A )
  diag.flags.writeable = False
  return diag

def reshape( A, *shape ):
  'more useful reshape'

  newshape = []
  i = 0
  for s in shape:
    if isinstance( s, (tuple,list) ):
      assert numpy.product( s ) == A.shape[i]
      newshape.extend( s )
      i += 1
    elif s == 1:
      newshape.append( A.shape[i] )
      i += 1
    else:
      assert s > 1
      newshape.append( numpy.product( A.shape[i:i+s] ) )
      i += s
  assert i <= A.ndim
  newshape.extend( A.shape[i:] )
  return A.reshape( newshape )

def mean( A, weights=None, axis=-1 ):
  'generalized mean'

  return A.mean( axis ) if weights is None else dot( A, weights / weights.sum(), axis )

def normalize( A, axis=-1 ):
  'devide by normal'

  s = [ slice(None) ] * A.ndim
  s[axis] = numpy.newaxis
  return A / numpy.linlg.norm( A, axis )[ tuple(s) ]

def cross( v1, v2, axis ):
  'cross product'

  if v1.ndim < v2.ndim:
    v1 = v1[ (numpy.newaxis,)*(v2.ndim-v1.ndim) ]
  elif v2.ndim < v1.ndim:
    v2 = v2[ (numpy.newaxis,)*(v1.ndim-v2.ndim) ]
  return numpy.cross( v1, v2, axis=axis )

def times( A, B ):
  """Times
  Multiply such that shapes are concatenated."""
  A = numpy.asarray( A )
  B = numpy.asarray( B )

  o = _abs[:A.ndim+B.ndim]
  m = o[:A.ndim]
  n = o[A.ndim:]

  return numpy.einsum( '%s,%s->%s' % (m,n,o), A, B )

def stack( arrays, axis=0 ):
  'powerful array stacker with singleton expansion'

  arrays = [ numpy.asarray(array,dtype=float) for array in arrays ]
  shape = [1] * max( array.ndim for array in arrays )
  axis = normdim( len(shape)+1, axis )
  for array in arrays:
    for i in range(-array.ndim,0):
      if shape[i] == 1:
        shape[i] = array.shape[i]
      else:
        assert array.shape[i] in ( shape[i], 1 )
  stacked = numpy.empty( shape[:axis]+[len(arrays)]+shape[axis:], dtype=float )
  for i, arr in enumerate( arrays ):
    stacked[(slice(None),)*axis+(i,)] = arr
  return stacked

def bringforward( arg, axis ):
  'bring axis forward'

  arg = numpy.asarray(arg)
  axis = normdim(arg.ndim,axis)
  if axis == 0:
    return arg
  return arg.transpose( [axis] + range(axis) + range(axis+1,arg.ndim) )

def diagonalize( arg ):
  'append axis, place last axis on diagonal of self and new'

  diagonalized = numpy.zeros( arg.shape + (arg.shape[-1],) )
  diag = numpy.einsum( '...ii->...i', diagonalized )
  assert diag.base is diagonalized
  diag.flags.writeable = True
  diag[:] = arg
  return diagonalized

def isbool( a ):
  return isboolarray( a ) and a.ndim == 0 or numpy.issubdtype( type(a), numpy.bool )

def isboolarray( a ):
  return isinstance( a, numpy.ndarray ) and numpy.issubdtype( a.dtype, numpy.bool )

def isint( a ):
  return isintarray( a ) and a.ndim == 0 or numpy.issubdtype( type(a), numpy.integer )

def isintarray( a ):
  return isinstance( a, numpy.ndarray ) and numpy.issubdtype( a.dtype, numpy.integer )

# vim:shiftwidth=2:foldmethod=indent:foldnestmax=2
