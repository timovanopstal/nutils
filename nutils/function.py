# -*- coding: utf8 -*-
#
# Module FUNCTION
#
# Part of Nutils: open source numerical utilities for Python. Jointly developed
# by HvZ Computational Engineering, TU/e Multiscale Engineering Fluid Dynamics,
# and others. More info at http://nutils.org <info@nutils.org>. (c) 2014

"""
The function module defines the :class:`Evaluable` class and derived objects,
commonly referred to as nutils functions. They represent mappings from a
:mod:`nutils.topology` onto Python space. The notabe class of :class:`Array`
objects map onto the space of Numpy arrays of predefined dimension and shape.
Most functions used in nutils applicatons are of this latter type, including the
geometry and function bases for analysis.

Nutils functions are essentially postponed python functions, stored in a tree
structure of input/output dependencies. Many :class:`Array` objects have
directly recognizable numpy equivalents, such as :class:`Sin` or
:class:`Inverse`. By not evaluating directly but merely stacking operations,
complex operations can be defined prior to entering a quadrature loop, allowing
for a higher level style programming. It also allows for automatic
differentiation and code optimization.

It is important to realize that nutils functions do not map for a physical
xy-domain but from a topology, where a point is characterized by the combination
of an element and its local coordinate. This is a natural fit for typical finite
element operations such as quadrature. Evaluation from physical coordinates is
possible only via inverting of the geometry function, which is a fundamentally
expensive and currently unsupported operation.
"""

from . import util, numpy, numeric, log, core, cache, transform, expression, _
import sys, warnings, itertools, functools, operator, inspect, numbers, builtins, re, types, collections.abc

isevaluable = lambda arg: isinstance(arg, Evaluable)

class Evaluable( cache.Immutable ):
  'Base class'

  def __init__(self, args:tuple):
    assert all(isevaluable(arg) for arg in args)
    self.__args = args

  def evalf( self, *args ):
    raise NotImplementedError( 'Evaluable derivatives should implement the evalf method' )

  @cache.property
  def dependencies(self):
    '''collection of all function arguments'''
    args = set()
    for func in self.__args:
      if func not in args:
        args |= func.dependencies
        args.add(func)
    return args

  @property
  def isconstant(self):
    return EVALARGS not in self.dependencies

  @cache.property
  def ordereddeps(self):
    '''collection of all function arguments such that the arguments to
    dependencies[i] can be found in dependencies[:i]'''
    return tuple([EVALARGS] + sorted(self.dependencies - {EVALARGS}, key=lambda f: len(f.dependencies)))

  @cache.property
  def dependencytree(self):
    '''lookup table of function arguments into ordereddeps, such that
    ordereddeps[i].__args[j] == ordereddeps[dependencytree[i][j]], and
    self.__args[j] == ordereddeps[dependencytree[-1][j]]'''
    args = self.ordereddeps
    return tuple(tuple(map(args.index, func.__args)) for func in args+(self,))

  @property
  def serialized(self):
    return zip(self.ordereddeps[1:]+(self,), self.dependencytree[1:])

  def asciitree( self, seen=None ):
    'string representation'

    if seen is None:
      seen = []
    try:
      index = seen.index( self )
    except ValueError:
      pass
    else:
      return '%{}'.format( index )
    asciitree = self._asciitree_str()
    if core.getprop( 'richoutput', False ):
      select = '├ ', '└ '
      bridge = '│ ', '  '
    else:
      select = ': ', ': '
      bridge = '| ', '  '
    for iarg, arg in enumerate( self.__args ):
      n = iarg >= len(self.__args) - 1
      asciitree += '\n' + select[n] + ( ('\n' + bridge[n]).join( arg.asciitree( seen ).splitlines() ) if isevaluable(arg) else '<{}>'.format(arg) )
    index = len(seen)
    seen.append( self )
    return '%{} = {}'.format( index, asciitree )

  def _asciitree_str(self):
    return str(self)

  def __str__( self ):
    return self.__class__.__name__

  def eval(self, **evalargs):
    values = [evalargs]
    for op, indices in self.serialized:
      try:
        args = [values[i] for i in indices]
        retval = op.evalf(*args)
      except KeyboardInterrupt:
        raise
      except:
        etype, evalue, traceback = sys.exc_info()
        excargs = etype, evalue, self, values
        raise EvaluationError(*excargs).with_traceback(traceback)
      values.append(retval)
    return values[-1]

  @log.title
  def graphviz( self ):
    'create function graph'

    import os, subprocess, hashlib

    dotpath = core.getprop( 'dot', True )
    if not isinstance( dotpath, str ):
      dotpath = 'dot'

    lines = []
    lines.append( 'digraph {' )
    lines.append( 'graph [ dpi=72 ];' )
    lines.extend( '%d [label="%d. %s"];' % (i, i, name._asciitree_str()) for i, name in enumerate(self.ordereddeps+(self,)) )
    lines.extend( '%d -> %d;' % (j,i) for i, indices in enumerate(self.dependencytree) for j in indices )
    lines.append( '}' )
    imgdata = '\n'.join(lines).encode()

    imgtype = core.getprop( 'imagetype', 'png' )
    imgpath = 'dot_{}.{}'.format(hashlib.sha1(imgdata).hexdigest(), imgtype)
    if not os.path.exists( imgpath ):
      with core.open_in_outdir( imgpath, 'w' ) as img:
        with subprocess.Popen( [dotpath,'-T'+imgtype], stdin=subprocess.PIPE, stdout=img ) as dot:
          dot.communicate( imgdata )

    log.info( imgpath )

  def stackstr( self, nlines=-1 ):
    'print stack'

    lines = ['  %0 = EVALARGS']
    for op, indices in self.serialized:
      args = [ '%%%d' % idx for idx in indices ]
      try:
        code = op.evalf.__code__
        offset = 1 if getattr( op.evalf, '__self__', None ) is not None else 0
        names = code.co_varnames[ offset:code.co_argcount ]
        names += tuple( '%s[%d]' % ( code.co_varnames[ code.co_argcount ], n ) for n in range( len(indices) - len(names) ) )
        args = [ '%s=%s' % item for item in zip( names, args ) ]
      except:
        pass
      lines.append( '  %%%d = %s( %s )' % (len(lines), op._asciitree_str(), ', '.join(args)) )
      if len(lines) == nlines+1:
        break
    return '\n'.join( lines )

  @cache.property
  def simplified(self):
    return self.edit(lambda arg: arg.simplified if isevaluable(arg) else arg)

class EvaluationError( Exception ):
  'evaluation error'

  def __init__( self, etype, evalue, evaluable, values ):
    'constructor'

    self.etype = etype
    self.evalue = evalue
    self.evaluable = evaluable
    self.values = values

  def __repr__( self ):
    return 'EvaluationError%s' % self

  def __str__( self ):
    'string representation'

    return '\n%s --> %s: %s' % ( self.evaluable.stackstr( nlines=len(self.values) ), self.etype.__name__, self.evalue )

EVALARGS = Evaluable(args=())

class Cache(Evaluable):
  def __init__(self):
    super().__init__(args=[EVALARGS])
  def evalf(self, evalargs):
    try:
      return evalargs['_cache']
    except:
      return cache.WrapperDummyCache()

CACHE = Cache()

class Trans(Evaluable):
  def __init__(self, n):
    self.n = n
    super().__init__(args=[EVALARGS])
  def evalf(self, evalargs):
    trans = evalargs['_transforms'][self.n]
    assert isinstance(trans, transform.TransformChain)
    return trans

TRANS = Trans(0)
OPPTRANS = Trans(1)

class Points(Evaluable):
  def __init__(self, opposite=False):
    super().__init__(args=[EVALARGS])
  def evalf(self, evalargs):
    points = evalargs['_points']
    assert numeric.isarray(points) and points.ndim == 2
    return numeric.const(points)

POINTS = Points()

class Tuple( Evaluable ):

  def __init__(self, items:tuple):
    self.items = items
    args = []
    indices = []
    for i, item in enumerate(self.items):
      if isevaluable(item):
        args.append(item)
        indices.append(i)
    self.indices = tuple(indices)
    super().__init__(args)

  @cache.property
  def simplified(self):
    return Tuple([item.simplified if isevaluable(item) else item for item in self.items])

  def edit(self, op):
    return Tuple([op(item) for item in self.items])

  def evalf( self, *items ):
    'evaluate'

    T = list(self.items)
    for index, item in zip( self.indices, items ):
      T[index] = item
    return tuple( T )

  def __iter__( self ):
    'iterate'

    return iter(self.items)

  def __len__( self ):
    'length'

    return len(self.items)

  def __getitem__( self, item ):
    'get item'

    return self.items[item]

  def __add__( self, other ):
    'add'

    return Tuple( self.items + tuple(other) )

  def __radd__( self, other ):
    'add'

    return Tuple( tuple(other) + self.items )

class SelectChain( Evaluable ):
  def __init__(self, trans, first:bool):
    self.trans = trans
    self.first = first
    super().__init__(args=[trans])
  def evalf( self, trans ):
    assert isinstance( trans, transform.TransformChain )
    bf = trans[0]
    assert isinstance( bf, transform.Bifurcate )
    ftrans = bf.trans1 if self.first else bf.trans2
    return transform.TransformChain( ftrans + trans[1:] )

class Promote(Evaluable):
  def __init__(self, ndims:int, trans):
    self.ndims = ndims
    super().__init__(args=[trans])
  def evalf(self, trans):
    head, tail = trans.canonical.promote(self.ndims)
    return transform.TransformChain(head + tail)

# ARRAYFUNC
#
# The main evaluable. Closely mimics a numpy array.

def add(a, b):
  a, b = _numpy_align(a, b)
  return Add([a, b])

def multiply(a, b):
  a, b = _numpy_align(a, b)
  return Multiply([a, b])

def sum(arg, axis=None):
  arg = asarray(arg)
  if axis is None:
    axis = numpy.arange(arg.ndim)
  elif not util.isiterable(axis):
    axis = numeric.normdim(arg.ndim, axis),
  else:
    axis = _norm_and_sort(arg.ndim, axis)
    assert numpy.greater(numpy.diff(axis), 0).all(), 'duplicate axes in sum'
  summed = arg
  for ax in reversed(axis):
    summed = Sum(summed, ax)
  return summed

def product(arg, axis):
  arg = asarray(arg)
  axis = numeric.normdim(arg.ndim, axis)
  shape = arg.shape[:axis] + arg.shape[axis+1:]
  trans = [i for i in range(arg.ndim) if i != axis] + [axis]
  return Product(transpose(arg, trans))

def power(arg, n):
  arg, n = _numpy_align(arg, n)
  return Power(arg, n)

def dot(a, b, axes=None):
  if axes is None:
    a = asarray(a)
    b = asarray(b)
    assert b.ndim == 1 and b.shape[0] == a.shape[0]
    for idim in range(1, a.ndim):
      b = insertaxis(b, idim, a.shape[idim])
    axes = 0,
  else:
    a, b = _numpy_align(a, b)
  if not util.isiterable(axes):
    axes = axes,
  axes = _norm_and_sort(a.ndim, axes)
  return Dot([a, b], axes)

def transpose(arg, trans=None):
  arg = asarray(arg)
  if trans is None:
    normtrans = range(arg.ndim-1, -1, -1)
  else:
    normtrans = _normdims(arg.ndim, trans)
    assert sorted(normtrans) == list(range(arg.ndim))
  return Transpose(arg, normtrans)

def swapaxes(arg, axis1, axis2):
  arg = asarray(arg)
  trans = numpy.arange(arg.ndim)
  trans[axis1], trans[axis2] = trans[axis2], trans[axis1]
  return transpose(arg, trans)

asdtype = lambda arg: arg if arg in (bool, int, float) else {'f': float, 'i': int, 'b': bool}[numpy.dtype(arg).kind]
asarray = lambda arg: arg if isarray(arg) else Constant(arg) if numeric.isarray(arg) or numpy.asarray(arg).dtype != object else stack(arg, axis=0)
asarrays = lambda args: tuple(asarray(arg) for arg in args)

class Array( Evaluable ):
  'array function'

  __array_priority__ = 1. # http://stackoverflow.com/questions/7042496/numpy-coercion-problem-for-left-sided-binary-operator/7057530#7057530

  def __init__(self, args:tuple, shape:tuple, dtype:asdtype):
    assert all(numeric.isint(sh) or isarray(sh) and sh.ndim == 0 and sh.dtype == int for sh in shape)
    self.shape = tuple(sh if not isarray(sh) else sh.eval()[0] if sh.isconstant else sh.simplified for sh in shape)
    self.ndim = len(shape)
    self.dtype = dtype
    super().__init__(args=args)

  def __getitem__(self, item):
    if not isinstance(item, tuple):
      item = item,
    iell = None
    nx = self.ndim - len(item)
    for i, it in enumerate(item):
      if it is ...:
        assert iell is None, 'at most one ellipsis allowed'
        iell = i
      elif it is _:
        nx += 1
    array = self
    axis = 0
    for it in item + (slice(None),)*nx if iell is None else item[:iell] + (slice(None),)*(nx+1) + item[iell+1:]:
      if numeric.isint(it):
        array = get(array, axis, item=it)
      elif it is _:
        array = expand_dims(array, axis)
        axis += 1
      elif it is slice(None):
        axis += 1
      elif isinstance(it, slice):
        assert it.step == None or it.step == 1
        start = 0 if it.start is None else it.start if it.start >= 0 else it.start + array.shape[axis]
        stop = array.shape[axis] if it.stop is None else it.stop if it.stop >= 0 else it.stop + array.shape[axis]
        array = take(array, index=Range(stop-start, start), axis=axis)
        axis += 1
      else:
        array = take(array, index=it, axis=axis)
        axis += 1
    assert axis == array.ndim
    return array

  def __len__(self):
    if self.ndim == 0:
      raise TypeError('len() of unsized object')
    return self.shape[0]

  def __iter__(self):
    if not self.shape:
      raise TypeError('iteration over a 0-d array')
    return ( self[i,...] for i in range(self.shape[0]) )

  size = property(lambda self: util.product(self.shape) if self.ndim else 1)
  T = property(lambda self: transpose(self))

  __add__ = __radd__ = add
  __sub__ = lambda self, other: subtract(self, other)
  __rsub__ = lambda self, other: subtract(other, self)
  __mul__ = __rmul__ = multiply
  __truediv__ = lambda self, other: divide(self, other)
  __rtruediv__ = lambda self, other: divide(other, self)
  __neg__ = lambda self: negative(self)
  __pow__ = power
  __abs__ = lambda self: abs(self)
  __mod__  = lambda self, other: mod(self, other)
  __str__ = __repr__ = lambda self: 'Array<{}>'.format(','.join(map(str, self.shape)) if hasattr(self, 'shape') else '?')

  sum = sum
  prod = product
  vector = lambda self, ndims: vectorize([self] * ndims)
  dot = dot
  normalized = lambda self, axis=-1: normalized(self, axis)
  normal = lambda self, exterior=False: normal(self, exterior)
  curvature = lambda self, ndims=-1: curvature(self, ndims)
  swapaxes = swapaxes
  transpose = transpose
  grad = lambda self, geom, ndims=0: grad(self, geom, ndims)
  laplace = lambda self, geom, ndims=0: grad(self, geom, ndims).div(geom, ndims)
  add_T = lambda self, axes=(-2,-1): add_T(self, axes)
  symgrad = lambda self, geom, ndims=0: symgrad(self, geom, ndims)
  div = lambda self, geom, ndims=0: div(self, geom, ndims)
  dotnorm = lambda self, geom, axis=-1: dotnorm(self, geom, axis)
  tangent = lambda self, vec: tangent(self, vec)
  ngrad = lambda self, geom, ndims=0: ngrad(self, geom, ndims)
  nsymgrad = lambda self, geom, ndims=0: nsymgrad(self, geom, ndims)

  @property
  def blocks(self):
    return [(tuple(Range(n) for n in self.shape), self)]

  def _asciitree_str(self):
    return '{}({})'.format(type(self).__name__, ','.join(['?' if isarray(sh) else str(sh) for sh in self.shape]))

  # simplifications
  _multiply = lambda self, other: None
  _transpose = lambda self, axes: None
  _insertaxis = lambda self, axis, length: None
  _dot = lambda self, other, axes: None
  _get = lambda self, i, item: None
  _power = lambda self, n: None
  _add = lambda self, other: None
  _sum = lambda self, axis: None
  _take = lambda self, index, axis: None
  _determinant = lambda self: None
  _inverse = lambda self: None
  _takediag = lambda self, axis, rmaxis: None
  _kronecker = lambda self, axis, length, pos: None
  _diagonalize = lambda self, axis, newaxis: None
  _product = lambda self: None
  _cross = lambda self, other, axis: None
  _sign = lambda self: None
  _eig = lambda self, symmetric: None
  _inflate = lambda self, dofmap, length, axis: None
  _mask = lambda self, maskvec, axis: None
  _unravel = lambda self, axis, shape: None
  _ravel = lambda self, axis: None

class Normal( Array ):
  'normal'

  def __init__(self, lgrad:asarray):
    assert lgrad.ndim == 2 and lgrad.shape[0] == lgrad.shape[1]
    self.lgrad = lgrad
    super().__init__(args=[lgrad], shape=(len(lgrad),), dtype=float)

  def evalf( self, lgrad ):
    n = lgrad[...,-1]
    if n.shape[-1] == 1: # geom is 1D
      return numpy.sign(n)
    # orthonormalize n to G
    G = lgrad[...,:-1]
    GG = numeric.contract( G[:,:,_,:], G[:,:,:,_], axis=1 )
    v1 = numeric.contract( G, n[:,:,_], axis=1 )
    v2 = numpy.linalg.solve( GG, v1 )
    v3 = numeric.contract( G, v2[:,_,:], axis=2 )
    return numeric.normalize( n - v3 )

  def _derivative(self, var, seen):
    if len(self) == 1:
      return zeros(self.shape + var.shape)
    G = self.lgrad[...,:-1]
    GG = matmat(G.T, G)
    Gder = derivative(G, var, seen)
    nGder = matmat(self, Gder)
    return -matmat(G, inverse(GG), nGder)

class ArrayFunc( Array ):
  'deprecated ArrayFunc alias'

  def __init__(self, args:tuple, shape:tuple):
    warnings.warn( 'function.ArrayFunc is deprecated; use function.Array instead', DeprecationWarning )
    super().__init__(args=args, shape=shape, dtype=float)

class Constant( Array ):

  def __init__(self, value:numeric.const):
    self.value = value
    super().__init__(args=[], shape=value.shape, dtype=value.dtype)

  @cache.property
  def simplified(self):
    if not self.value.any():
      return zeros_like(self)
    return self

  def evalf( self ):
    return self.value[_]

  @cache.property
  def _isunit( self ):
    return numpy.equal(self.value, 1).all()

  def _derivative(self, var, seen):
    return zeros(self.shape + var.shape)

  def _transpose(self, axes):
    return Constant(self.value.transpose(axes))

  def _sum(self, axis):
    return Constant(numpy.sum(self.value, axis))

  def _get(self, i, item):
    if item.isconstant:
      item, = item.eval()
      return Constant(numeric.get(self.value, i, item))

  def _add(self, other):
    if isinstance(other, Constant):
      return Constant(numpy.add(self.value, other.value))

  def _inverse(self):
    return Constant(numpy.linalg.inv(self.value))

  def _product(self):
    return Constant(self.value.prod(-1))

  def _multiply(self, other):
    if self._isunit:
      return other
    if isinstance(other, Constant):
      return Constant(numpy.multiply(self.value, other.value))

  def _takediag(self, axis, rmaxis):
    return Constant(numeric.takediag(self.value, axis, rmaxis))

  def _take( self, index, axis ):
    if isinstance(index, Constant):
      return Constant(self.value.take(index.value, axis))

  def _power(self, n):
    if isinstance(n, Constant):
      return Constant(numeric.power(self.value, n.value))

  def _dot( self, other, axes ):
    if isinstance(other, Constant):
      return Constant(numeric.contract(self.value, other.value, axes))
    if self._isunit:
      summed = other
      for axis in reversed(sorted(axes)):
        summed = Sum(summed, axis)
      return summed

  def _cross(self, other, axis):
    if isinstance(other, Constant):
      return Constant(numeric.cross(self.value, other.value, axis))

  def _eig(self, symmetric):
    eigval, eigvec = (numpy.linalg.eigh if symmetric else numpy.linalg.eig)(self.value)
    return Tuple((Constant(eigval), Constant(eigvec)))

  def _sign(self):
    return Constant(numeric.sign(self.value))

  def _unravel(self, axis, shape):
    shape = self.value.shape[:axis] + shape + self.value.shape[axis+1:]
    return Constant(self.value.reshape(shape))

  def _mask(self, maskvec, axis):
    return Constant(self.value[(slice(None),)*axis+(maskvec,)])

class DofMap(Array):

  def __init__(self, dofs:tuple, index:asarray):
    assert index.ndim == 0 and index.dtype == int
    self.dofs = dofs
    self.index = index
    length = get([len(d) for d in dofs] + [0], iax=0, item=index)
    super().__init__(args=[index], shape=(length,), dtype=int)

  @property
  def dofmap(self):
    return self.index.asdict(self.dofs)

  def evalf(self, index):
    index, = index
    return (self.dofs[index] if index < len(self.dofs) else numpy.empty([0], dtype=int))[_]

class ElementSize( Array):
  'dimension of hypercube with same volume as element'

  def __init__(self, geometry:asarray, ndims:int=0):
    assert geometry.ndim == 1
    self.ndims = len(geometry)+ndims if ndims <= 0 else ndims
    iwscale = jacobian( geometry, self.ndims )
    super().__init__(args=[iwscale], shape=(), dtype=float)

  def evalf( self, iwscale ):
    volume = iwscale.sum()
    return numeric.power( volume, 1/self.ndims )[_]

class InsertAxis(Array):

  def __init__(self, func:asarray, axis:int, length:asarray):
    assert length.ndim == 0 and length.dtype == int
    assert 0 <= axis <= func.ndim
    self.func = func
    self.axis = axis
    self.length = length
    super().__init__(args=[func, length], shape=func.shape[:axis]+(length,)+func.shape[axis:], dtype=func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    retval = func._insertaxis(self.axis, self.length)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return InsertAxis(func, self.axis, self.length)

  def evalf(self, func, length):
    length, = length
    return numeric.const(func).insertaxis(self.axis+1, length)

  def _derivative(self, var, seen):
    return insertaxis(derivative(self.func, var, seen), self.axis, self.length)

  def _get(self, i, item):
    if i == self.axis:
      if item.isconstant and self.length.isconstant:
        assert item.eval()[0] < self.length.eval()[0]
      return self.func
    return InsertAxis(Get(self.func, i-(i>self.axis), item), self.axis-(i<self.axis), self.length)

  def _sum(self, i):
    if i == self.axis:
      return Multiply([self.func, _inflate_scalar(self.length, self.func.shape)])
    return InsertAxis(Sum(self.func, i-(i>self.axis)), self.axis-(i<self.axis), self.length)

  def _product(self):
    if self.axis == self.ndim-1:
      return Power(self.func, _inflate_scalar(self.length, self.func.shape))
    return InsertAxis(Product(self.func), self.axis, self.length)

  def _power(self, n):
    if isinstance(n, InsertAxis) and self.axis == n.axis:
      assert n.length == self.length
      return InsertAxis(Power(self.func, n.func), self.axis, self.length)

  def _add(self, other):
    if isinstance(other, InsertAxis) and self.axis == other.axis:
      assert self.length == other.length
      return InsertAxis(Add([self.func, other.func]), self.axis, self.length)

  def _multiply(self, other):
    if isinstance(other, InsertAxis) and self.axis == other.axis:
      assert self.length == other.length
      return InsertAxis(Multiply([self.func, other.func]), self.axis, self.length)

  def _insertaxis(self, axis, length):
    if (not length.isconstant, axis) < (not self.length.isconstant, self.axis):
      return InsertAxis(InsertAxis(self.func, axis-(axis>self.axis), length), self.axis+(axis<=self.axis), self.length)

  def _take(self, index, axis):
    if axis == self.axis:
      return InsertAxis(self.func, self.axis, index.shape[0])
    return InsertAxis(Take(self.func, index, axis-(axis>self.axis)), self.axis, self.length)

  def _takediag(self, axis, rmaxis):
    if self.axis == rmaxis:
      return self.func
    elif self.axis == axis:
      return Transpose(self.func, list(range(axis))+[rmaxis-1]+list(range(axis, rmaxis-1))+list(range(rmaxis, self.func.ndim)))
    else:
      return InsertAxis(TakeDiag(self.func, axis-(self.axis<axis), rmaxis-(self.axis<rmaxis)), self.axis-(self.axis>rmaxis), self.length)

  def _dot(self, other, axes):
    if self.axis in axes:
      assert other.shape[self.axis] == self.shape[self.axis]
      return Dot([self.func, Sum(other, self.axis)], [ax-(ax>self.axis) for ax in axes if ax != self.axis])

  def _mask(self, maskvec, axis):
    if axis == self.axis:
      assert len(maskvec) == self.shape[self.axis]
      return InsertAxis(self.func, self.axis, maskvec.sum())
    return InsertAxis(Mask(self.func, maskvec, axis-(self.axis<axis)), self.axis, self.length)

  def _transpose(self, axes):
    i = axes.index(self.axis)
    return InsertAxis(Transpose(self.func, [ax-(ax>self.axis) for ax in axes[:i]+axes[i+1:]]), i, self.length)

  def _unravel(self, axis, shape):
    if axis == self.axis:
      return InsertAxis(InsertAxis(self.func, self.axis, shape[1]), self.axis, shape[0])
    else:
      return InsertAxis(Unravel(self.func, axis-(axis>self.axis), shape), self.axis+(axis<self.axis), self.length)

class Transpose(Array):

  def __init__(self, func:asarray, axes:tuple):
    assert sorted(axes) == list(range(func.ndim))
    self.func = func
    self.axes = axes
    super().__init__(args=[func], shape=[func.shape[n] for n in axes], dtype=func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    if self.axes == tuple(range(self.ndim)):
      return func
    retval = func._transpose(self.axes)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Transpose(func, self.axes)

  def evalf(self, arr):
    return arr.transpose([0] + [n+1 for n in self.axes])

  def _transpose(self, axes):
    newaxes = [self.axes[i] for i in axes]
    return Transpose(self.func, newaxes)

  def _takediag(self, axis, rmaxis):
    if self.axes[axis] < self.axes[rmaxis]:
      axes = self.axes
    else:
      axes = list(self.axes)
      axes[axis], axes[rmaxis] = axes[rmaxis], axes[axis]
    assert axes[axis] < axes[rmaxis]
    return Transpose(TakeDiag(self.func, axes[axis], axes[rmaxis]), [ax-(ax>axes[rmaxis]) for ax in axes[:rmaxis]+axes[rmaxis+1:]])

  def _get(self, i, item):
    axis = self.axes[i]
    axes = [ax-(ax>axis) for ax in self.axes if ax != axis]
    return Transpose(Get(self.func, axis, item), axes)

  def _sum(self, i):
    axis = self.axes[i]
    axes = [ax-(ax>axis) for ax in self.axes if ax != axis]
    return Transpose(Sum(self.func, axis), axes)

  def _derivative(self, var, seen):
    return transpose(derivative(self.func, var, seen), self.axes+tuple(range(self.ndim, self.ndim+var.ndim)))

  def _multiply(self, other):
    if isinstance(other, Transpose) and self.axes == other.axes:
      return Transpose(Multiply([self.func, other.func]), self.axes)
    other_trans = other._transpose(_invtrans(self.axes))
    if other_trans is not None:
      return Transpose(Multiply([self.func, other_trans]), self.axes)

  def _add(self, other):
    if isinstance(other, Transpose) and self.axes == other.axes:
      return Transpose(Add([self.func, other.func]), self.axes)
    other_trans = other._transpose(_invtrans(self.axes))
    if other_trans is not None:
      return Transpose(Add([self.func, other_trans]), self.axes)

  def _take(self, indices, axis):
    return Transpose(Take(self.func, indices, self.axes[axis]), self.axes)

  def _dot(self, other, axes):
    sumaxes = [self.axes[axis] for axis in axes]
    trydot = self.func._dot(transpose(other, _invtrans(self.axes)), sumaxes)
    if trydot is not None:
      trans = [axis - builtins.sum(ax<axis for ax in sumaxes) for axis in self.axes if axis not in sumaxes]
      return Transpose(trydot, trans)

  def _mask(self, maskvec, axis):
    return Transpose(Mask(self.func, maskvec, self.axes[axis]), self.axes)

class Get(Array):

  def __init__(self, func:asarray, axis:int, item:asarray):
    assert item.ndim == 0 and item.dtype == int
    self.func = func
    self.axis = axis
    self.item = item
    assert 0 <= axis < func.ndim, 'axis is out of bounds'
    if item.isconstant and numeric.isint(func.shape[axis]):
      assert 0 <= item.eval()[0] < func.shape[axis], 'item is out of bounds'
    super().__init__(args=[func, item], shape=func.shape[:axis]+func.shape[axis+1:], dtype=func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    item = self.item.simplified
    retval = func._get(self.axis, item)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Get(func, self.axis, item)

  def evalf(self, arr, item):
    if len(item) == 1:
      item, = item
      p = slice(None)
    else:
      p = numpy.arange(len(item))
    return arr[(p,)+(slice(None),)*self.axis+(item,)]

  def _derivative(self, var, seen):
    f = derivative(self.func, var, seen)
    return get(f, self.axis, self.item)

  def _get(self, i, item):
    tryget = self.func._get(i+(i>=self.axis), item)
    if tryget is not None:
      return Get(tryget, self.axis, self.item)

  def _take(self, indices, axis):
    return Get(Take(self.func, indices, axis+(axis>=self.axis) ), self.axis, self.item)

class Product( Array ):

  def __init__(self, func:asarray):
    self.func = func
    super().__init__(args=[func], shape=func.shape[:-1], dtype=func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    retval = func._product()
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Product(func)

  def evalf( self, arr ):
    assert arr.ndim == self.ndim+2
    return numpy.product( arr, axis=-1 )

  def _derivative(self, var, seen):
    grad = derivative(self.func, var, seen)
    funcs = Stack([util.product(self.func[...,j] for j in range(self.func.shape[-1]) if i != j) for i in range(self.func.shape[-1])], axis=self.ndim)
    return (grad * funcs[(...,)+(_,)*var.ndim]).sum(self.ndim)

    ## this is a cleaner form, but is invalid if self.func contains zero values:
    #ext = (...,)+(_,)*len(shape)
    #return self[ext] * ( derivative(self.func,var,shape,seen) / self.func[ext] ).sum( self.ndim )

  def _get(self, i, item):
    func = Get(self.func, i, item)
    return Product(func)

class RootCoords( Array ):

  def __init__(self, ndims:int, trans=TRANS):
    self.trans = trans
    super().__init__(args=[POINTS,trans], shape=[ndims], dtype=float)

  def evalf( self, points, chain ):
    'evaluate'

    ndims = len(self)
    head, tail = chain.promote( ndims )
    while head and head[0].todims != ndims:
      head = head[1:]
    return transform.apply( head + tail, points )

  def _derivative(self, var, seen):
    if isinstance(var, LocalCoords) and len(var) > 0:
      return RootTransform(len(self), len(var), self.trans)
    return zeros(self.shape+var.shape)

class RootTransform( Array ):

  def __init__(self, ndims:int, nvars:int, trans):
    super().__init__(args=[Promote(ndims, trans)], shape=(ndims,nvars), dtype=float)

  def evalf(self, chain):
    todims, fromdims = self.shape
    while chain and chain[0].todims != todims:
      chain = chain[1:]
    return transform.linearfrom(chain, fromdims)[_]

  def _derivative(self, var, seen):
    return zeros(self.shape+var.shape)

class Function( Array ):

  def __init__(self, stds:tuple, depth:int, trans, index:asarray, derivs:tuple=()):
    assert index.ndim == 0 and index.dtype == int
    self.stds = stds
    self.depth = depth
    self.trans = trans
    self.index = index
    nshapes = get([std.nshapes for std in stds] + [0], iax=0, item=index)
    super().__init__(args=(CACHE,POINTS,trans,index), shape=(nshapes,)+derivs, dtype=float)

  @property
  def stdmap(self):
    return self.index.asdict(self.stds)

  def evalf(self, cache, points, trans, index):
    index, = index
    if index == len(self.stds):
      return numpy.empty((1,0)+self.shape[1:])
    tail = trans[self.depth:]
    if tail:
      points = cache[transform.apply](tail, points)
    fvals = cache[self.stds[index].eval](points, self.ndim-1)
    assert fvals.ndim == self.ndim+1
    if tail:
      for i, ndims in enumerate(self.shape[1:]):
        linear = cache[transform.linearfrom](tail, ndims)
        fvals = numeric.dot(fvals, linear, axis=i+2)
    return fvals

  def _derivative(self, var, seen):
    if isinstance(var, LocalCoords):
      return Function(self.stds, self.depth, self.trans, self.index, self.shape[1:]+var.shape)
    return zeros(self.shape+var.shape, dtype=self.dtype)

class Inverse( Array ):

  def __init__(self, func:asarray):
    assert func.ndim >= 2 and func.shape[-1] == func.shape[-2]
    self.func = func
    super().__init__(args=[func], shape=func.shape, dtype=float)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    retval = func._inverse()
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Inverse(func)

  def evalf( self, arr ):
    assert arr.ndim == self.ndim+1
    try:
      inv = numpy.linalg.inv( arr )
    except numpy.linalg.LinAlgError:
      inv = numpy.empty_like( arr )
      flat = (-1,) + arr.shape[-2:]
      for arri, invi in zip( arr.reshape(flat), inv.reshape(flat) ):
        try:
          invi[...] = numpy.linalg.inv(arri)
        except numpy.linalg.LinAlgError:
          invi[...] = numpy.nan
    return inv

  def _derivative(self, var, seen):
    G = derivative(self.func, var, seen)
    n = var.ndim
    a = slice(None)
    return -sum(self[(...,a,a,_,_)+(_,)*n] * G[(...,_,a,a,_)+(a,)*n] * self[(...,_,_,a,a)+(_,)*n], [-2-n, -3-n])

  def _eig(self, symmetric):
    eigval, eigvec = Eig(self.func, symmetric)
    return Tuple((reciprocal(eigval), eigvec))

class Concatenate(Array):

  def __init__(self, funcs:tuple, axis:int=0):
    ndim = funcs[0].ndim
    assert all(isarray(func) and func.ndim == ndim for func in funcs)
    assert 0 <= axis < ndim
    assert all(func.shape[:axis] == funcs[0].shape[:axis] and func.shape[axis+1:] == funcs[0].shape[axis+1:] for func in funcs[1:])
    length = util.sum(func.shape[axis] for func in funcs)
    shape = funcs[0].shape[:axis] + (length,) + funcs[0].shape[axis+1:]
    dtype = _jointdtype(*[func.dtype for func in funcs])
    self.funcs = funcs
    self.axis = axis
    super().__init__(args=funcs, shape=shape, dtype=dtype)

  def edit(self, op):
    return Concatenate([op(func) for func in self.funcs], self.axis)

  @cache.property
  def _withslices(self):
    return tuple((Range(func.shape[self.axis], n), func) for n, func in zip(util.cumsum(func.shape[self.axis] for func in self.funcs), self.funcs))

  @cache.property
  def simplified(self):
    funcs = tuple(func.simplified for func in self.funcs if func.shape[self.axis] != 0)
    if all(iszero(func) for func in funcs):
      return zeros_like(self)
    if len(funcs) == 1:
      return funcs[0]
    return Concatenate(funcs, self.axis)

  def evalf(self, *arrays):
    shape = list(builtins.max(arrays, key=len).shape)
    shape[self.axis+1] = builtins.sum(array.shape[self.axis+1] for array in arrays)
    retval = numpy.empty(shape, dtype=self.dtype)
    n0 = 0
    for array in arrays:
      n1 = n0 + array.shape[self.axis+1]
      retval[(slice(None),)*(self.axis+1)+(slice(n0,n1),)] = array
      n0 = n1
    assert n0 == retval.shape[self.axis+1]
    return retval

  @cache.property
  def blocks(self):
    return _concatblocks(((ind[:self.axis], ind[self.axis+1:]), (ind[self.axis]+n, f))
      for n, func in zip(util.cumsum(func.shape[self.axis] for func in self.funcs), self.funcs)
        for ind, f in func.blocks)

  def _get(self, i, item):
    if i != self.axis:
      axis = self.axis - (self.axis > i)
      return Concatenate([Get(f, i, item) for f in self.funcs], axis=axis)
    if item.isconstant:
      item, = item.eval()
      for f in self.funcs:
        if item < f.shape[i]:
          return Get(f, i, item)
        item -= f.shape[i]
      raise Exception

  def _derivative(self, var, seen):
    funcs = [derivative(func, var, seen) for func in self.funcs]
    return concatenate(funcs, axis=self.axis)

  def _multiply(self, other):
    funcs = [Multiply([func, Take(other, s, self.axis)]) for s, func in self._withslices]
    return Concatenate(funcs, self.axis)

  def _cross(self, other, axis):
    if axis != self.axis:
      funcs = [Cross(func, Take(other, s, self.axis), axis) for s, func in self._withslices]
      return Concatenate(funcs, self.axis)

  def _add(self, other):
    if isinstance(other, Concatenate) and self.axis == other.axis:
      if [f1.shape[self.axis] for f1 in self.funcs] == [f2.shape[self.axis] for f2 in other.funcs]:
        funcs = [add(f1, f2) for f1, f2 in zip(self.funcs, other.funcs)]
      else:
        if isarray(self.shape[self.axis]):
          raise NotImplementedError
        funcs = []
        beg1 = 0
        for func1 in self.funcs:
          end1 = beg1 + func1.shape[self.axis]
          beg2 = 0
          for func2 in other.funcs:
            end2 = beg2 + func2.shape[self.axis]
            if end1 > beg2 and end2 > beg1:
              mask = numpy.zeros(self.shape[self.axis], dtype=bool)
              mask[builtins.max(beg1, beg2):builtins.min(end1, end2)] = True
              funcs.append(Add([Mask(func1, mask[beg1:end1], self.axis), Mask(func2, mask[beg2:end2], self.axis)]))
            beg2 = end2
          beg1 = end1
    else:
      funcs = [Add([func, Take(other, s, self.axis)]) for s, func in self._withslices]
    return Concatenate(funcs, self.axis)

  def _sum(self, axis):
    funcs = [Sum(func, axis) for func in self.funcs]
    if axis == self.axis:
      while len(funcs) > 1:
        funcs[-2:] = Add(funcs[-2:]),
      return funcs[0]
    return Concatenate(funcs, self.axis - (axis<self.axis))

  def _transpose(self, axes):
    funcs = [Transpose(func, axes) for func in self.funcs]
    axis = axes.index(self.axis)
    return Concatenate(funcs, axis)

  def _insertaxis(self, axis, length):
    funcs = [InsertAxis(func, axis, length) for func in self.funcs]
    return Concatenate(funcs, self.axis+(axis<=self.axis))

  def _takediag(self, axis, rmaxis):
    if self.axis == axis:
      funcs = [TakeDiag(Take(func, s, rmaxis), axis, rmaxis) for s, func in self._withslices]
      return Concatenate(funcs, axis=axis)
    elif self.axis == rmaxis:
      funcs = [TakeDiag(Take(func, s, axis), axis, rmaxis) for s, func in self._withslices]
      return Concatenate(funcs, axis=axis)
    else:
      return Concatenate([TakeDiag(f, axis, rmaxis) for f in self.funcs], axis=self.axis-(self.axis>rmaxis))

  def _take(self, indices, axis):
    if axis != self.axis:
      return Concatenate([Take(func, indices, axis) for func in self.funcs], self.axis)
    if not indices.isconstant:
      return
    indices, = indices.eval()
    assert numpy.logical_and(numpy.greater_equal(indices, 0), numpy.less(indices, self.shape[axis])).all()
    ifuncs = numpy.hstack([ numpy.repeat(ifunc,func.shape[axis]) for ifunc, func in enumerate(self.funcs) ])[indices]
    splits, = numpy.nonzero( numpy.diff(ifuncs) != 0 )
    funcs = []
    for i, j in zip( numpy.hstack([ 0, splits+1 ]), numpy.hstack([ splits+1, len(indices) ]) ):
      ifunc = ifuncs[i]
      assert numpy.equal(ifuncs[i:j], ifunc).all()
      offset = builtins.sum(func.shape[axis] for func in self.funcs[:ifunc])
      funcs.append(Take(self.funcs[ifunc], indices[i:j] - offset, axis))
    if len(funcs) == 1:
      return funcs[0]
    return Concatenate(funcs, axis=axis)

  def _dot(self, other, axes):
    funcs = [Dot([func, Take(other, s, self.axis)], axes) for s, func in self._withslices]
    if self.axis in axes:
      while len(funcs) > 1:
        funcs[-2:] = Add(funcs[-2:]),
      return funcs[0]
    return Concatenate(funcs, self.axis - builtins.sum(axis < self.axis for axis in axes))

  def _power(self, n):
    return Concatenate([Power(func, Take(n, s, self.axis)) for s, func in self._withslices], self.axis)

  def _diagonalize(self, axis, newaxis):
    if self.axis != axis:
      return Concatenate([Diagonalize(func, axis, newaxis) for func in self.funcs], self.axis+(newaxis<=self.axis))

  def _kronecker(self, axis, length, pos):
    return Concatenate([kronecker(func,axis,length,pos) for func in self.funcs], self.axis+(axis<=self.axis))

  def _mask( self, maskvec, axis ):
    if axis != self.axis:
      return Concatenate([Mask(func,maskvec,axis) for func in self.funcs], self.axis)
    if all(s.isconstant for s, func in self._withslices):
      return Concatenate([Mask(func, maskvec[s.eval()[0]], axis) for s, func in self._withslices], axis)

  def _unravel(self, axis, shape):
    if axis != self.axis:
      return Concatenate([Unravel(func, axis, shape) for func in self.funcs], self.axis+(self.axis>axis))

class Interpolate( Array ):
  'interpolate uniformly spaced data; stepwise for now'

  def __init__(self, x:asarray, xp:numeric.const, fp:numeric.const, left=None, right=None):
    assert xp.ndim == fp.ndim == 1
    if not numpy.greater(numpy.diff(xp), 0).all():
      warnings.warn( 'supplied x-values are non-increasing' )
    assert x.ndim == 0
    self.xp = xp
    self.fp = fp
    self.left = left
    self.right = right
    super.__init__(args=[x], shape=(), dtype=float)

  def evalf( self, x ):
    return numpy.interp( x, self.xp, self.fp, self.left, self.right )

class Cross( Array ):

  def __init__(self, func1:asarray, func2:asarray, axis:int):
    assert func1.shape == func2.shape
    assert 0 <= axis < func1.ndim and func2.shape[axis] == 3
    self.func1 = func1
    self.func2 = func2
    self.axis = axis
    super().__init__(args=(func1,func2), shape=func1.shape, dtype=_jointdtype(func1.dtype, func2.dtype))

  @cache.property
  def simplified(self):
    func1 = self.func1.simplified
    func2 = self.func2.simplified
    retval = func1._cross(func2, self.axis)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    retval = func2._cross(func1, self.axis)
    if retval is not None:
      assert retval.shape == self.shape
      return negative(retval).simplified
    return Cross(func1, func2, self.axis)

  def evalf( self, a, b ):
    assert a.ndim == b.ndim == self.ndim+1
    return numeric.cross( a, b, self.axis+1 )

  def _derivative(self, var, seen):
    ext = (...,)+(_,)*var.ndim
    return cross(self.func1[ext], derivative(self.func2, var, seen), axis=self.axis) \
         - cross(self.func2[ext], derivative(self.func1, var, seen), axis=self.axis)

  def _take(self, index, axis):
    if axis != self.axis:
      return Cross(Take(self.func1, index, axis), Take(self.func2, index, axis), self.axis)

class Determinant( Array ):

  def __init__(self, func:asarray):
    assert isarray(func) and func.ndim >= 2 and func.shape[-1] == func.shape[-2]
    self.func = func
    super().__init__(args=[func], shape=func.shape[:-2], dtype=func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    retval = func._determinant()
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Determinant(func)

  def evalf( self, arr ):
    assert arr.ndim == self.ndim+3
    return numpy.linalg.det( arr )

  def _derivative(self, var, seen):
    Finv = swapaxes(inverse(self.func), -2, -1)
    G = derivative(self.func, var, seen)
    ext = (...,)+(_,)*var.ndim
    return self[ext] * sum(Finv[ext] * G, axis=[-2-var.ndim,-1-var.ndim])

class Multiply(Array):

  def __init__(self, funcs:util.frozenmultiset):
    self.funcs = funcs
    func1, func2 = funcs
    assert isarray(func1) and isarray(func2) and func1.shape == func2.shape
    super().__init__(args=self.funcs, shape=func1.shape, dtype=_jointdtype(func1.dtype,func2.dtype))

  def edit(self, op):
    return Multiply([op(func) for func in self.funcs])

  @cache.property
  def simplified(self):
    func1, func2 = [func.simplified for func in self.funcs]
    if func1 == func2:
      return power(func1, 2).simplified
    retval = func1._multiply(func2)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    retval = func2._multiply(func1)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Multiply([func1, func2])

  def evalf( self, arr1, arr2 ):
    return arr1 * arr2

  def _sum(self, axis):
    func1, func2 = self.funcs
    return Dot([func1, func2], [axis])

  def _get(self, axis, item):
    func1, func2 = self.funcs
    return Multiply([Get(func1, axis, item), Get(func2, axis, item)])

  def _add(self, other):
    func1, func2 = self.funcs
    if other == func1:
      return Multiply([func1, Add([func2, ones_like(func2)])])
    if other == func2:
      return Multiply([func2, Add([func1, ones_like(func1)])])
    if isinstance(other, Multiply) and not self.funcs.isdisjoint(other.funcs):
      f = next(iter(self.funcs & other.funcs))
      return Multiply([f, Add(self.funcs + other.funcs - [f,f])])

  def _determinant( self ):
    func1, func2 = self.funcs
    if self.shape[-2:] == (1,1):
      return Multiply([Determinant(func1), Determinant(func2)])

  def _product(self):
    func1, func2 = self.funcs
    return Multiply([Product(func1), Product(func2)])

  def _multiply(self, other):
    func1, func2 = self.funcs
    func1_other = func1._multiply(other)
    if func1_other is not None:
      return Multiply([func1_other, func2])
    func2_other = func2._multiply(other)
    if func2_other is not None:
      return Multiply([func1, func2_other])

  def _derivative(self, var, seen):
    func1, func2 = self.funcs
    ext = (...,)+(_,)*var.ndim
    return func1[ext] * derivative(func2, var, seen) \
         + func2[ext] * derivative(func1, var, seen)

  def _takediag(self, axis, rmaxis):
    func1, func2 = self.funcs
    return Multiply([TakeDiag(func1, axis, rmaxis), TakeDiag(func2, axis, rmaxis)])

  def _take(self, index, axis):
    func1, func2 = self.funcs
    return Multiply([Take(func1, index, axis), Take(func2, index, axis)])

  def _power(self, n):
    func1, func2 = self.funcs
    func1pow = func1._power(n)
    func2pow = func2._power(n)
    if func1pow is not None and func2pow is not None:
      return Multiply([func1pow, func2pow])

  def _dot(self, other, axes):
    func1, func2 = self.funcs
    trydot1 = func1._dot(Multiply([func2, other]), axes)
    if trydot1 is not None:
      return trydot1
    trydot2 = func2._dot(Multiply([func1, other]), axes)
    if trydot2 is not None:
      return trydot2

class Add(Array):

  def __init__(self, funcs:util.frozenmultiset):
    self.funcs = funcs
    func1, func2 = funcs
    assert isarray(func1) and isarray(func2) and func1.shape == func2.shape
    super().__init__(args=self.funcs, shape=func1.shape, dtype=_jointdtype(func1.dtype,func2.dtype))

  def edit(self, op):
    return Add([op(func) for func in self.funcs])

  @cache.property
  def simplified(self):
    func1, func2 = [func.simplified for func in self.funcs]
    if iszero(func1):
      return func2
    if iszero(func2):
      return func1
    if func1 == func2:
      return multiply(func1, 2).simplified
    retval = func1._add(func2)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    retval = func2._add(func1)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Add([func1, func2])

  def evalf( self, arr1, arr2=None ):
    return arr1 + arr2

  def _sum(self, axis):
    return Add([Sum(func, axis) for func in self.funcs])

  def _derivative(self, var, seen):
    func1, func2 = self.funcs
    return derivative(func1, var, seen) + derivative(func2, var, seen)

  def _get(self, axis, item):
    func1, func2 = self.funcs
    return Add([Get(func1, axis, item), Get(func2, axis, item)])

  def _takediag(self, axis, rmaxis):
    func1, func2 = self.funcs
    return Add([TakeDiag(func1, axis, rmaxis), TakeDiag(func2, axis, rmaxis)])

  def _take(self, index, axis):
    func1, func2 = self.funcs
    return Add([Take(func1, index, axis), Take(func2, index, axis)])

  def _add(self, other):
    func1, func2 = self.funcs
    func1_other = func1._add(other)
    if func1_other is not None:
      return Add([func1_other, func2])
    func2_other = func2._add(other)
    if func2_other is not None:
      return Add([func1, func2_other])

  def _mask(self, maskvec, axis):
    func1, func2 = self.funcs
    return Add([Mask(func1, maskvec, axis), Mask(func2, maskvec, axis)])

class BlockAdd( Array ):
  'block addition (used for DG)'

  def __init__(self, funcs:util.frozenmultiset):
    self.funcs = funcs
    shapes = set(func.shape for func in funcs)
    assert len(shapes) == 1, 'multiple shapes in BlockAdd'
    shape, = shapes
    super().__init__(args=funcs, shape=shape, dtype=_jointdtype(*[func.dtype for func in self.funcs]))

  def edit(self, op):
    return BlockAdd([op(func) for func in self.funcs])

  @cache.property
  def simplified(self):
    funcs = []
    for func in self.funcs:
      func = func.simplified
      if isinstance(func, BlockAdd):
        funcs.extend(func.funcs)
      elif not iszero(func):
        funcs.append(func)
    return BlockAdd(funcs) if len(funcs) > 1 else funcs[0] if funcs else zeros_like(self)

  def evalf(self, *args):
    return util.sum(args)

  def _add(self, other):
    return BlockAdd(tuple(self.funcs) + tuple(other.funcs if isinstance(other, BlockAdd) else [other]))

  def _dot( self, other, axes ):
    return BlockAdd([Dot([func, other], axes) for func in self.funcs])

  def _sum(self, axis):
    return BlockAdd([sum(func, axis) for func in self.funcs])

  def _derivative(self, var, seen):
    return BlockAdd([derivative(func, var, seen) for func in self.funcs])

  def _get(self, i, item):
    return BlockAdd([Get(func, i, item) for func in self.funcs])

  def _takediag(self, axis, rmaxis):
    return BlockAdd([TakeDiag(func, axis, rmaxis) for func in self.funcs])

  def _take(self, indices, axis):
    return BlockAdd([take(func, indices, axis) for func in self.funcs])

  def _transpose(self, axes):
    return BlockAdd([Transpose(func, axes) for func in self.funcs])

  def _insertaxis(self, axis, length):
    return BlockAdd([InsertAxis(func, axis, length) for func in self.funcs])

  def _multiply(self, other):
    return BlockAdd([multiply(func, other) for func in self.funcs])

  def _kronecker(self, axis, length, pos):
    return BlockAdd([kronecker(func, axis, length, pos) for func in self.funcs])

  def _mask(self, maskvec, axis):
    return BlockAdd([Mask(func, maskvec, axis) for func in self.funcs])

  def _unravel(self, axis, shape):
    return BlockAdd([unravel(func, axis, shape) for func in self.funcs])

  @cache.property
  def blocks(self):
    gathered = tuple((ind, util.sum(f)) for ind, f in util.gather(block for func in self.funcs for block in func.blocks))
    if len(gathered) > 1:
      for idim in range(self.ndim):
        gathered = _concatblocks(((ind[:idim], ind[idim+1:]), (ind[idim], f)) for ind, f in gathered)
    return gathered

class Dot(Array):

  def __init__(self, funcs:util.frozenmultiset, axes:tuple):
    self.funcs = funcs
    func1, func2 = funcs
    assert isarray(func1) and isarray(func2) and func1.shape == func2.shape
    self.axes = axes
    assert all(0 <= ax < func1.ndim for ax in axes)
    assert all(ax1 < ax2 for ax1, ax2 in zip(axes[:-1], axes[1:]))
    shape = func1.shape
    self.axes_complement = list(range(func1.ndim))
    for ax in reversed(self.axes):
      shape = shape[:ax] + shape[ax+1:]
      del self.axes_complement[ax]
    _abc = numeric._abc[:func1.ndim+1]
    self._einsumfmt = '{0},{0}->{1}'.format(_abc, ''.join(a for i, a in enumerate(_abc) if i-1 not in axes))
    super().__init__(args=funcs, shape=shape, dtype=_jointdtype(func1.dtype,func2.dtype))

  def edit(self, op):
    return Dot([op(func) for func in self.funcs], self.axes)

  @cache.property
  def simplified(self):
    func1, func2 = [func.simplified for func in self.funcs]
    if len(self.axes) == 0:
      return multiply(func1, func2).simplified
    if iszero(func1) or iszero(func2):
      return zeros(self.shape)
    for i, axis in enumerate(self.axes):
      if func1.shape[axis] == 1:
        return dot(sum(func1,axis), sum(func2,axis), self.axes[:i] + tuple(axis-1 for axis in self.axes[i+1:])).simplified
    retval = func1._dot(func2, self.axes)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    retval = func2._dot(func1, self.axes)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Dot([func1, func2], self.axes)

  def evalf( self, arr1, arr2 ):
    return numpy.einsum(self._einsumfmt, arr1, arr2)

  def _get(self, axis, item):
    func1, func2 = self.funcs
    funcaxis = self.axes_complement[axis]
    return Dot([Get(func1, funcaxis, item), Get(func2, funcaxis, item)], [ax-(ax>=funcaxis) for ax in self.axes])

  def _derivative(self, var, seen):
    func1, func2 = self.funcs
    ext = (...,)+(_,)*var.ndim
    return dot(derivative(func1, var, seen), func2[ext], self.axes) \
         + dot(func1[ext], derivative(func2, var, seen), self.axes)

  def _add(self, other):
    if isinstance(other, Dot) and self.axes == other.axes and not self.funcs.isdisjoint(other.funcs):
      f = next(iter(self.funcs & other.funcs))
      return Dot([f, Add(self.funcs + other.funcs - [f,f])], self.axes)

  def _takediag(self, axis, rmaxis):
    func1, func2 = self.funcs
    faxis = self.axes_complement[axis]
    frmaxis = self.axes_complement[rmaxis]
    return Dot([TakeDiag(func1, faxis, frmaxis), TakeDiag(func2, faxis, frmaxis)], [ax-(ax>frmaxis) for ax in self.axes])

  def _sum(self, axis):
    funcaxis = self.axes_complement[axis]
    func1, func2 = self.funcs
    return Dot([func1, func2], sorted(self.axes + (funcaxis,)))

  def _take(self, index, axis):
    func1, func2 = self.funcs
    funcaxis = self.axes_complement[axis]
    return Dot([Take(func1, index, funcaxis), Take(func2, index, funcaxis)], self.axes)

class Sum( Array ):

  def __init__(self, func:asarray, axis:int):
    self.axis = axis
    self.func = func
    assert 0 <= axis < func.ndim, 'axis out of bounds'
    shape = func.shape[:axis] + func.shape[axis+1:]
    super().__init__(args=[func], shape=shape, dtype=int if func.dtype == bool else func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    retval = func._sum(self.axis)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Sum(func, self.axis)

  def evalf( self, arr ):
    assert arr.ndim == self.ndim+2
    return numpy.sum(arr, self.axis+1)

  def _sum(self, axis):
    trysum = self.func._sum(axis+(axis>=self.axis))
    if trysum is not None:
      return Sum(trysum, self.axis-(axis<self.axis))

  def _derivative(self, var, seen):
    return sum(derivative(self.func, var, seen), self.axis)

class Debug( Array ):
  'debug'

  def __init__(self, func:asarray):
    self.func = func
    super().__init__(args=[func], shape=func.shape, dtype=func.dtype)

  def evalf( self, arr ):
    'debug'

    assert arr.ndim == self.ndim+1
    log.debug( 'debug output:\n%s' % arr )
    return arr

  def __str__( self ):
    'string representation'

    return '{DEBUG}'

  def _derivative(self, var, seen):
    return Debug(derivative(self.func, var, seen))

class TakeDiag( Array ):

  def __init__(self, func:asarray, axis:int, rmaxis:int):
    assert func.shape[axis] == func.shape[rmaxis]
    assert 0 <= axis < rmaxis < func.ndim
    self.func = func
    self.axis = axis
    self.rmaxis = rmaxis
    super().__init__(args=[func], shape=func.shape[:rmaxis]+func.shape[rmaxis+1:], dtype=func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    if self.shape[self.axis] == 1:
      return get(func, self.rmaxis, 0).simplified
    retval = func._takediag(self.axis, self.rmaxis)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return TakeDiag(func, self.axis, self.rmaxis)

  def evalf(self, arr):
    assert arr.ndim == self.ndim+2
    return numeric.takediag(arr, self.axis+1, self.rmaxis+1)

  def _derivative(self, var, seen):
    return TakeDiag(derivative(self.func, var, seen), self.axis, self.rmaxis)

  def _sum(self, axis):
    if axis != self.axis:
      return TakeDiag(Sum(self.func, axis+(axis>=self.rmaxis)), self.axis-(axis<self.axis), self.rmaxis-(axis<self.rmaxis))

class Take( Array ):

  def __init__(self, func:asarray, indices:asarray, axis:int):
    assert indices.ndim == 1 and indices.dtype == int
    assert 0 <= axis < func.ndim
    self.func = func
    self.axis = axis
    self.indices = indices
    shape = func.shape[:axis] + indices.shape + func.shape[axis+1:]
    super().__init__(args=[func,indices], shape=shape, dtype=func.dtype)

  @cache.property
  def simplified(self):
    if self.shape[self.axis] == 0:
      return zeros(self.shape, dtype=self.dtype)
    func = self.func.simplified
    indices = self.indices.simplified
    length = self.func.shape[self.axis]
    if indices == Range(length):
      return func
    if indices.isconstant and numeric.isint(length):
      indices_, = indices.eval()
      if numpy.greater(numpy.diff(numpy.mod(indices_, length)), 0).all():
        mask = numpy.zeros(length, dtype=bool)
        mask[indices_] = True # note: includes proper bounds check
        return Mask(func, mask, self.axis).simplified
    retval = func._take(indices, self.axis)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Take(func, indices, self.axis)

  def evalf( self, arr, indices ):
    if indices.shape[0] != 1:
      raise NotImplementedError( 'non element-constant indexing not supported yet' )
    return numpy.take( arr, indices[0], self.axis+1 )

  def _derivative(self, var, seen):
    return take(derivative(self.func, var, seen), self.indices, self.axis)

  def _take(self, index, axis):
    if axis == self.axis:
      return Take(self.func, self.indices[index], axis)
    trytake = self.func._take(index, axis)
    if trytake is not None:
      return Take(trytake, self.indices, self.axis)

class Power(Array):

  def __init__(self, func:asarray, power:asarray):
    assert func.shape == power.shape
    self.func = func
    self.power = power
    super().__init__(args=[func,power], shape=func.shape, dtype=float)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    power = self.power.simplified
    if iszero(power):
      return ones_like(self).simplified
    retval = func._power(power)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Power(func, power)

  def evalf( self, base, exp ):
    return numeric.power( base, exp )

  def _derivative(self, var, seen):
    ext = (...,)+(_,)*var.ndim
    if self.power.isconstant:
      p, = self.power.eval()
      return zeros(self.shape + var.shape) if p == 0 \
        else multiply(p, power(self.func, p-1))[ext] * derivative(self.func, var, seen)
    # self = func**power
    # ln self = power * ln func
    # self` / self = power` * ln func + power * func` / func
    # self` = power` * ln func * self + power * func` * func**(power-1)
    return (self.power * power(self.func, self.power - 1))[ext] * derivative(self.func, var, seen) \
         + (ln(self.func) * self)[ext] * derivative(self.power, var, seen)

  def _power(self, n):
    func = self.func
    newpower = Multiply([self.power, n])
    if iszero(self.power % 2) and not iszero(newpower % 2):
      func = abs(func)
    return Power(func, newpower)

  def _get(self, axis, item):
    return Power(Get(self.func, axis, item), Get(self.power, axis, item))

  def _sum(self, axis):
    if self == (self.func**2):
      return Dot([self.func, self.func], [axis])

  def _takediag(self, axis, rmaxis):
    return Power(TakeDiag(self.func, axis, rmaxis), TakeDiag(self.power, axis, rmaxis))

  def _take(self, index, axis):
    return Power(Take(self.func, index, axis), Take(self.power, index, axis))

  def _multiply(self, other):
    if isinstance(other, Power) and self.func == other.func:
      return Power(self.func, Add([self.power, other.power]))
    if other == self.func:
      return Power(self.func, Add([self.power, ones_like(self.power)]))

  def _sign( self ):
    if iszero(self.power % 2):
      return ones_like(self)

class Pointwise( Array ):

  deriv = None

  def __init__(self, *args:asarrays):
    retval = self.evalf(*[numpy.ones((), dtype=arg.dtype) for arg in args])
    shapes = set(arg.shape for arg in args)
    assert len(shapes) == 1, 'pointwise arguments have inconsistent shapes'
    shape, = shapes
    self.args = args
    super().__init__(args=args, shape=shape, dtype=retval.dtype)

  @cache.property
  def simplified(self):
    args = [arg.simplified for arg in self.args]
    if all(arg.isconstant for arg in args):
      retval, = self.evalf(*[arg.eval() for arg in args])
      return Constant(retval).simplified
    return self.__class__(*args)

  def _derivative(self, var, seen):
    if self.deriv is None:
      raise NotImplementedError('derivative is not defined for this operator')
    return util.sum(deriv(*self.args)[(...,)+(_,)*var.ndim] * derivative(arg, var, seen) for arg, deriv in zip(self.args, self.deriv))

  def _takediag(self, axis, rmaxis):
    return self.__class__(*[TakeDiag(arg, axis, rmaxis) for arg in self.args])

  def _get(self, axis, item):
    return self.__class__(*[Get(arg, axis, item) for arg in self.args])

  def _take(self, index, axis):
    return self.__class__(*[Take(arg, index, axis) for arg in self.args])

class Cos(Pointwise):
  evalf = numpy.cos
  deriv = lambda x: -Sin(x),

class Sin(Pointwise):
  evalf = numpy.sin
  deriv = Cos,

class Tan(Pointwise):
  evalf = numpy.tan
  deriv = lambda x: Cos(x)**-2,

class ArcSin(Pointwise):
  evalf = numpy.arcsin
  deriv = lambda x: reciprocal(sqrt(1-x**2)),

class ArcCos(Pointwise):
  evalf = numpy.arccos
  deriv = lambda x: -reciprocal(sqrt(1-x**2)),

class Exp(Pointwise):
  evalf = numpy.exp
  deriv = lambda x: Exp(x),

class Log(Pointwise):
  evalf = numpy.log
  deriv = lambda x: reciprocal(x),

class Mod(Pointwise):
  evalf = numpy.mod

class ArcTan2(Pointwise):
  evalf = numpy.arctan2
  deriv = lambda x, y: y / (x**2 + y**2), lambda x, y: -x / (x**2 + y**2)

class Greater(Pointwise):
  evalf = numpy.greater
  deriv = (lambda a, b: Zeros(a.shape, dtype=int),) * 2

class Equal(Pointwise):
  evalf = numpy.equal
  deriv = (lambda a, b: Zeros(a.shape, dtype=int),) * 2

class Less(Pointwise):
  evalf = numpy.less
  deriv = (lambda a, b: Zeros(a.shape, dtype=int),) * 2

class Minimum(Pointwise):
  evalf = numpy.minimum
  deriv = Less, lambda x, y: 1 - Less(x, y)

class Maximum(Pointwise):
  evalf = numpy.maximum
  deriv = lambda x, y: 1 - Less(x, y), Less

class Int(Pointwise):
  evalf = staticmethod(lambda a: a.astype(int))
  deriv = lambda a: Zeros(a.shape, int),

class Sign( Array ):

  def __init__(self, func:asarray):
    self.func = func
    super().__init__(args=[func], shape=func.shape, dtype=func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    retval = func._sign()
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Sign(func)

  def evalf( self, arr ):
    assert arr.ndim == self.ndim+1
    return numpy.sign( arr )

  def _derivative(self, var, seen):
    return zeros(self.shape + var.shape)

  def _takediag(self, axis, rmaxis):
    return Sign(TakeDiag(self.func, axis, rmaxis))

  def _get(self, axis, item):
    return Sign(Get(self.func, axis, item))

  def _take(self, index, axis):
    return Sign(Take(self.func, index, axis))

  def _sign( self ):
    return self

  def _power(self, n):
    if iszero(n % 2):
      return ones_like(self)

class Sampled( Array ):
  'sampled'

  def __init__(self, data:util.frozendict, trans=TRANS):
    self.data = data.copy()
    self.trans = trans
    items = iter(self.data.items())
    trans0, (values0,points0) = next(items)
    shape = values0.shape[1:]
    assert all( transi.fromdims == trans0.fromdims and valuesi.shape == pointsi.shape[:1]+shape for transi, (valuesi,pointsi) in items )
    super().__init__(args=[trans,POINTS], shape=shape, dtype=float)

  def evalf( self, trans, points ):
    (myvals,mypoints), tail = trans.lookup_item( self.data )
    evalpoints = tail.apply( points )
    assert mypoints.shape == evalpoints.shape and numpy.equal(mypoints, evalpoints).all(), 'Illegal point set'
    return myvals

class Elemwise( Array ):
  'elementwise constant data'

  def __init__(self, fmap:util.frozendict, shape:tuple, default=None, trans=TRANS):
    self.fmap = fmap
    self.default = default
    self.trans = trans
    super().__init__(args=[trans], shape=shape, dtype=float)

  def evalf( self, trans ):
    try:
      value, tail = trans.lookup_item( self.fmap )
    except KeyError:
      value = self.default
      if value is None:
        raise
    value = numpy.asarray( value )
    assert value.shape == self.shape, 'wrong shape: {} != {}'.format( value.shape, self.shape )
    return value[_]

  def _derivative(self, var, seen):
    return zeros(self.shape+var.shape)

class Eig( Evaluable ):

  def __init__(self, func:asarray, symmetric:bool=False):
    assert func.ndim >= 2 and func.shape[-1] == func.shape[-2]
    self.symmetric = symmetric
    self.func = func
    super().__init__(args=[func])

  def __len__(self):
    return 2

  def __iter__(self):
    yield ArrayFromTuple(self, index=0, shape=self.func.shape[:-1], dtype=float)
    yield ArrayFromTuple(self, index=1, shape=self.func.shape, dtype=float)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    retval = func._eig(self.symmetric)
    if retval is not None:
      assert len(retval) == 2
      return retval.simplified
    return Eig(func, self.symmetric)

  def evalf(self, arr):
    return (numpy.linalg.eigh if self.symmetric else numeric.eig)(arr)

class ArrayFromTuple(Array):

  def __init__(self, arrays, index:int, shape:tuple, dtype:asdtype):
    assert isevaluable(arrays)
    assert 0 <= index < len(arrays)
    self.arrays = arrays
    self.index = index
    super().__init__(args=[arrays], shape=shape, dtype=dtype)

  def evalf(self, arrays):
    assert isinstance(arrays, tuple)
    return arrays[self.index]

class Zeros( Array ):
  'zero'

  def __init__(self, shape:tuple, dtype:asdtype):
    super().__init__(args=[asarray(sh) for sh in shape], shape=shape, dtype=dtype)

  def evalf(self, *shape):
    if shape:
      shape, = zip(*shape)
    return numpy.zeros((1,)+shape, dtype=self.dtype)

  @property
  def blocks(self):
    return ()

  def _derivative(self, var, seen):
    return zeros(self.shape+var.shape, dtype=self.dtype)

  def _add(self, other):
    return other

  def _multiply(self, other):
    return self

  def _dot(self, other, axes):
    shape = [sh for axis, sh in enumerate(self.shape) if axis not in axes]
    return Zeros(shape, dtype=_jointdtype(self.dtype,other.dtype))

  def _cross(self, other, axis):
    return self

  def _diagonalize(self, axis, newaxis):
    return Zeros(self.shape[:newaxis]+(self.shape[axis],)+self.shape[newaxis:], dtype=self.dtype)

  def _sum(self, axis):
    return Zeros(self.shape[:axis] + self.shape[axis+1:], dtype=self.dtype)

  def _transpose(self, axes):
    shape = [self.shape[n] for n in axes]
    return Zeros(shape, dtype=self.dtype)

  def _insertaxis(self, axis, length):
    return Zeros(self.shape[:axis]+(length,)+self.shape[axis:], self.dtype)

  def _get(self, i, item):
    return Zeros(self.shape[:i] + self.shape[i+1:], dtype=self.dtype)

  def _takediag(self, axis, rmaxis):
    return Zeros(self.shape[:rmaxis]+self.shape[rmaxis+1:], dtype=self.dtype)

  def _take(self, index, axis):
    return Zeros(self.shape[:axis] + index.shape + self.shape[axis+1:], dtype=self.dtype)

  def _inflate(self, dofmap, length, axis):
    assert not isinstance( self.shape[axis], int )
    return Zeros(self.shape[:axis] + (length,) + self.shape[axis+1:], dtype=self.dtype)

  def _power(self, n):
    return self

  def _kronecker(self, axis, length, pos):
    return Zeros(self.shape[:axis]+(length,)+self.shape[axis:], dtype=self.dtype)

  def _mask(self, maskvec, axis):
    return Zeros(self.shape[:axis] + (maskvec.sum(),) + self.shape[axis+1:], dtype=self.dtype)

  def _unravel( self, axis, shape ):
    shape = self.shape[:axis] + shape + self.shape[axis+1:]
    return Zeros(shape, dtype=self.dtype)

  def _ravel(self, axis):
    return Zeros(self.shape[:axis] + (self.shape[axis]*self.shape[axis+1],) + self.shape[axis+2:], self.dtype)

class Inflate( Array ):

  def __init__(self, func:asarray, dofmap:asarray, length:int, axis:int):
    assert not dofmap.isconstant
    self.func = func
    self.dofmap = dofmap
    self.length = length
    self.axis = axis
    assert 0 <= axis < func.ndim
    assert func.shape[axis] == dofmap.shape[0]
    shape = func.shape[:axis] + (length,) + func.shape[axis+1:]
    super().__init__(args=[func,dofmap], shape=shape, dtype=func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    dofmap = self.dofmap.simplified
    retval = func._inflate(dofmap, self.length, self.axis)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Inflate(func, dofmap, self.length, self.axis)

  def evalf(self, array, indices):
    assert indices.shape[0] == 1
    indices, = indices
    assert array.ndim == self.ndim+1
    warnings.warn( 'using explicit inflation; this is usually a bug.' )
    shape = list(array.shape)
    shape[self.axis+1] = self.length
    inflated = numpy.zeros(shape, dtype=self.dtype)
    inflated[(slice(None),)*(self.axis+1)+(indices,)] = array
    return inflated

  @property
  def blocks(self):
    for ind, f in self.func.blocks:
      assert ind[self.axis] == Range(self.func.shape[self.axis])
      yield (ind[:self.axis] + (self.dofmap,) + ind[self.axis+1:]), f

  def _mask(self, maskvec, axis):
    if axis != self.axis:
      return Inflate(Mask(self.func, maskvec, axis), self.dofmap, self.length, self.axis)
    newlength = maskvec.sum()
    selection = Take(maskvec, self.dofmap, axis=0)
    renumber = numpy.empty( len(maskvec), dtype=int )
    renumber[:] = newlength # out of bounds
    renumber[maskvec] = numpy.arange(newlength)
    newdofmap = Take(renumber, Take(self.dofmap, Find(selection), axis=0 ), axis=0)
    newfunc = Take(self.func, Find(selection), axis=self.axis)
    return Inflate(newfunc, newdofmap, newlength, self.axis)

  def _inflate( self, dofmap, length, axis ):
    assert axis != self.axis
    if axis > self.axis:
      return
    return Inflate(Inflate(self.func, dofmap, length, axis), self.dofmap, self.length, self.axis)

  def _derivative(self, var, seen):
    return inflate(derivative(self.func, var, seen), self.dofmap, self.length, self.axis)

  def _transpose(self, axes):
    axis = axes.index(self.axis)
    return Inflate(Transpose(self.func, axes), self.dofmap, self.length, axis)

  def _insertaxis(self, axis, length):
    return Inflate(InsertAxis(self.func, axis, length), self.dofmap, self.length, self.axis+(axis<=self.axis))

  def _get(self, axis, item):
    assert axis != self.axis
    return Inflate(Get(self.func,axis,item), self.dofmap, self.length, self.axis-(axis<self.axis))

  def _dot(self, other, axes):
    if isinstance(other, Inflate) and other.axis == self.axis:
      assert self.dofmap == other.dofmap
      other = other.func
    else:
      other = Take(other, self.dofmap, self.axis)
    arr = Dot([self.func, other], axes )
    if self.axis in axes:
      return arr
    return Inflate(arr, self.dofmap, self.length, self.axis - builtins.sum(axis < self.axis for axis in axes))

  def _multiply(self, other):
    if isinstance(other, Inflate) and self.axis == other.axis:
      assert self.dofmap == other.dofmap and self.length == other.length
      take_other = other.func
    else:
      take_other = Take(other, self.dofmap, self.axis)
    return Inflate(Multiply([self.func, take_other]), self.dofmap, self.length, self.axis)

  def _add(self, other):
    if isinstance(other, Inflate) and self.axis == other.axis and self.dofmap == other.dofmap:
      return Inflate(Add([self.func, other.func]), self.dofmap, self.length, self.axis)
    return BlockAdd([self, other])

  def _cross(self, other, axis):
    if isinstance(other, Inflate) and self.axis == other.axis:
      assert self.dofmap == other.dofmap
      other = other.func
    else:
      other = Take(other, self.dofmap, self.axis)
    return Inflate(Cross(self.func,other,axis), self.dofmap, self.length, self.axis)

  def _power(self, n):
    return Inflate(Power(self.func, n), self.dofmap, self.length, self.axis)

  def _takediag(self, axis, rmaxis):
    if self.axis == axis:
      return Inflate(TakeDiag(take(self.func, self.dofmap, rmaxis), axis, rmaxis), self.dofmap, self.length, axis)
    elif self.axis == rmaxis:
      return Inflate(TakeDiag(take(self.func, self.dofmap, axis), axis, rmaxis), self.dofmap, self.length, axis)
    else:
      return Inflate(TakeDiag(self.func, axis, rmaxis), self.dofmap, self.length, self.axis-(self.axis>rmaxis))

  def _take(self, index, axis):
    if axis != self.axis:
      return Inflate(Take(self.func, index, axis), self.dofmap, self.length, self.axis)
    if index == self.dofmap:
      return self.func

  def _diagonalize(self, axis, newaxis):
    if self.axis != axis:
      return Inflate(Diagonalize(self.func, axis, newaxis), self.dofmap, self.length, self.axis+(newaxis<=self.axis))

  def _sum(self, axis):
    arr = Sum(self.func, axis)
    if axis == self.axis:
      return arr
    return Inflate(arr, self.dofmap, self.length, self.axis-(axis<self.axis))

  def _kronecker(self, axis, length, pos):
    return Inflate(kronecker(self.func,axis,length,pos), self.dofmap, self.length, self.axis+(axis<=self.axis))

  def _unravel(self, axis, shape):
    if axis != self.axis:
      return Inflate(Unravel(self.func, axis, shape), self.dofmap, self.length, self.axis+(self.axis>axis))

class Diagonalize( Array ):

  def __init__(self, func:asarray, axis=int, newaxis=int):
    assert 0 <= axis < newaxis <= func.ndim
    self.func = func
    self.axis = axis
    self.newaxis = newaxis
    super().__init__(args=[func], shape=func.shape[:newaxis]+(func.shape[axis],)+func.shape[newaxis:], dtype=func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    if func.shape[self.axis] == 1:
      return insertaxis(func, self.newaxis, 1).simplified
    retval = func._diagonalize(self.axis, self.newaxis)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Diagonalize(func, self.axis, self.newaxis)

  def evalf( self, arr):
    assert arr.ndim == self.ndim
    return numeric.diagonalize(arr, self.axis+1, self.newaxis+1)

  def _derivative(self, var, seen):
    return diagonalize(derivative(self.func, var, seen), self.axis, self.newaxis)

  def _get(self, i, item):
    if i != self.axis and i != self.newaxis:
      return Diagonalize(Get(self.func, i-(i>self.newaxis), item), self.axis-(i<self.axis), self.newaxis-(i<self.newaxis))
    if item.isconstant:
      pos, = item.eval()
      funcs = [Zeros(self.func.shape[:self.axis]+self.func.shape[self.axis+1:], dtype=self.func.dtype)] * self.func.shape[self.axis]
      funcs[pos] = Get(self.func, self.axis, item)
      return Stack(funcs, axis=self.axis if i == self.newaxis else self.newaxis-1)

  def _inverse(self):
    if self.axis == self.func.ndim-1 and self.newaxis == self.ndim-1:
      return Diagonalize(reciprocal(self.func), self.axis, self.newaxis)

  def _determinant(self):
    if self.axis == self.func.ndim-1 and self.newaxis == self.ndim-1:
      return Product(Transpose(self.func, list(range(self.axis))+list(range(self.axis+1,self.func.ndim))+[self.axis]))

  def _multiply(self, other):
    return Diagonalize(Multiply([self.func, TakeDiag(other, self.axis, self.newaxis)]), self.axis, self.newaxis)

  def _dot(self, other, axes):
    faxes = [axis-(axis>self.newaxis) for axis in axes if axis not in (self.axis, self.newaxis)]
    assert self.axis not in faxes
    if len(faxes) < len(axes): # one of or both diagonalized axes are summed
      if len(faxes) == len(axes) - 2:
        faxes.append(self.axis)
      retval = Dot([self.func, TakeDiag(other, self.axis, self.newaxis)], faxes)
      if len(faxes) == len(axes) - 1 and self.axis in axes:
        axis = self.axis-builtins.sum(ax<self.axis for ax in faxes)
        newaxis = self.newaxis-builtins.sum(ax<self.newaxis for ax in faxes)
        retval = Transpose(retval, list(range(axis))+list(range(axis+1,newaxis))+[axis]+list(range(newaxis,retval.ndim)))
      return retval
    return Diagonalize(Dot([self.func, TakeDiag(other, self.axis, self.newaxis)], faxes), self.axis-builtins.sum(ax<self.axis for ax in faxes), self.newaxis-builtins.sum(ax<self.newaxis for ax in faxes))

  def _add(self, other):
    if isinstance(other, Diagonalize) and other.axis == self.axis and other.newaxis == self.newaxis:
      return Diagonalize(Add([self.func, other.func]), self.axis, self.newaxis)

  def _sum(self, axis):
    if axis == self.newaxis:
      return self.func
    if axis == self.axis:
      return Transpose(self.func, list(range(self.axis))+list(range(self.axis+1,self.newaxis))+[self.axis]+list(range(self.newaxis,self.func.ndim)))
    return Diagonalize(Sum(self.func, axis-(axis>self.newaxis)), self.axis-(axis<self.axis), self.newaxis-(axis<self.newaxis))

  def _transpose(self, axes):
    axis = axes.index(self.axis)
    newaxis = axes.index(self.newaxis)
    if newaxis < axis:
      axes = list(axes)
      axes[axis] = self.newaxis
      axes[newaxis] = self.axis
      axis, newaxis = newaxis, axis
    newaxes = [ax-(ax>self.newaxis) for ax in axes[:newaxis]+axes[newaxis+1:]]
    return Diagonalize(Transpose(self.func, newaxes), axis, newaxis)

  def _insertaxis(self, axis, length):
    return Diagonalize(InsertAxis(self.func, axis-(axis>self.newaxis), length), self.axis+(axis<=self.axis), self.newaxis+(axis<=self.newaxis))

  def _takediag(self, axis, rmaxis):
    if self.axis == axis and self.newaxis == rmaxis:
      return self.func

  def _take(self, index, axis):
    if axis not in (self.axis, self.newaxis):
      return Diagonalize(Take(self.func, index, axis-(axis>self.newaxis)), self.axis, self.newaxis)
    if numeric.isint(self.func.shape[self.axis]):
      diag = Diagonalize(Take(self.func, index, self.axis), self.axis, self.newaxis)
      return Inflate(diag, index, self.func.shape[self.axis], self.newaxis if axis == self.axis else self.axis)

  def _mask(self, maskvec, axis):
    if axis not in (self.axis, self.newaxis):
      return Diagonalize(Mask(self.func, maskvec, axis-(axis>self.newaxis)), self.axis, self.newaxis)
    indices, = numpy.where(maskvec)
    if not numpy.equal(numpy.diff(indices), 1).all():
      return
    # consecutive sub-block
    ax = self.axis if axis == self.newaxis else self.newaxis
    masked = Diagonalize(Mask(self.func, maskvec, self.axis), self.axis, self.newaxis)
    return Concatenate([Zeros(masked.shape[:ax] + (indices[0],) + masked.shape[ax+1:], dtype=self.dtype), masked, Zeros(masked.shape[:ax] + (self.shape[ax]-(indices[-1]+1),) + masked.shape[ax+1:], dtype=self.dtype)], axis=ax)

  def _unravel(self, axis, shape):
    if axis == self.axis:
      return Ravel(Diagonalize(Diagonalize(Unravel(self.func, self.axis, shape), self.axis, self.newaxis+1), self.axis+1, self.newaxis+2), self.newaxis+1)

class Guard( Array ):
  'bar all simplifications'

  def __init__(self, fun:asarray):
    self.fun = fun
    super().__init__(args=[fun], shape=fun.shape, dtype=fun.dtype)

  @staticmethod
  def evalf( dat ):
    return dat

  def _derivative(self, var, seen):
    return Guard(derivative(self.fun, var, seen))

class TrigNormal( Array ):
  'cos, sin'

  def __init__(self, angle:asarray):
    assert angle.ndim == 0
    self.angle = angle
    super().__init__(args=[angle], shape=(2,), dtype=float)

  def _derivative(self, var, seen):
    return trigtangent(self.angle)[(...,)+(_,)*var.ndim] * derivative(self.angle, var, seen)

  def evalf( self, angle ):
    return numpy.array([ numpy.cos(angle), numpy.sin(angle) ]).T

  def _dot( self, other, axes ):
    assert axes == (0,)
    if isinstance( other, (TrigTangent,TrigNormal) ) and self.angle == other.angle:
      return Constant(1. if isinstance(other,TrigNormal) else 0.)

class TrigTangent( Array ):
  '-sin, cos'

  def __init__(self, angle:asarray):
    assert angle.ndim == 0
    self.angle = angle
    super().__init__(args=[angle], shape=(2,), dtype=float)

  def _derivative(self, var, seen):
    return -trignormal(self.angle)[(...,)+(_,)*var.ndim] * derivative(self.angle, var, seen)

  def evalf( self, angle ):
    return numpy.array([ -numpy.sin(angle), numpy.cos(angle) ]).T

  def _dot(self, other, axes):
    assert axes == (0,)
    if isinstance( other, (TrigTangent,TrigNormal) ) and self.angle == other.angle:
      return Constant(1. if isinstance(other,TrigTangent) else 0.)

class Find( Array ):
  'indices of boolean index vector'

  def __init__(self, where:asarray):
    assert isarray(where) and where.ndim == 1 and where.dtype == bool
    self.where = where
    super().__init__(args=[where], shape=[where.sum()], dtype=int)

  def evalf( self, where ):
    assert where.shape[0] == 1
    where, = where
    index, = where.nonzero()
    return index[_]

class Stack( Array ):

  def __init__(self, funcs:tuple, axis:int):
    shapes = set(func.shape for func in funcs if func is not None)
    assert shapes, 'cannot determine shape of stack'
    assert len(shapes) == 1, 'multiple shapes in stack'
    shape, = shapes
    dtype = _jointdtype(*[func.dtype for func in funcs if func is not None])
    assert 0 <= axis <= len(shape)
    self.funcs = funcs
    self.axis = axis
    self.nz = tuple(ifunc for ifunc, func in enumerate(funcs) if func is not None)
    super().__init__(args=[funcs[i] for i in self.nz], shape=shape[:axis]+(len(funcs),)+shape[axis:], dtype=dtype)

  def edit(self, op):
    return Stack([op(func) if func is not None else None for func in self.funcs], self.axis)

  @cache.property
  def simplified(self):
    if len(self.funcs) == 1:
      return InsertAxis(self.funcs[0], axis=self.axis, length=1).simplified
    krons = Zeros(self.shape, self.dtype)
    funcs = [None] * len(self.funcs)
    for ifunc in self.nz:
      func = self.funcs[ifunc].simplified
      kron = func._kronecker(self.axis, len(self.funcs), ifunc)
      if kron is not None:
        assert kron.shape == self.shape
        krons = Add([krons, kron]).simplified
      elif not iszero(func):
        funcs[ifunc] = func
    if tuple(funcs) == self.funcs: # avoid recursion
      assert iszero(krons)
      return self
    if all(func is None for func in funcs):
      return krons
    return Add([Stack(funcs, self.axis), krons]).simplified

  def evalf(self, *funcs):
    shape = builtins.max(funcs, key=len).shape
    array = numpy.zeros(shape[:self.axis+1] + (len(self.funcs),) + shape[self.axis+1:], dtype=self.dtype)
    for i, func in zip(self.nz, funcs):
      array[(slice(None),)*(self.axis+1)+(i,)] = func
    return array

  def _derivative(self, var, seen):
    return Stack([derivative(func, var, seen) if func is not None else None for func in self.funcs], self.axis)

  def _get(self, i, item):
    if i != self.axis:
      return Stack([Get(func,i-(i>self.axis),item) if func is not None else None for func in self.funcs], self.axis-(i<self.axis))
    if item.isconstant:
      item, = item.eval()
      func = self.funcs[item]
      if func is None:
        return Zeros(self.shape[:self.axis]+self.shape[self.axis+1:], dtype=self.dtype)
      return func

  def _add(self, other):
    if isinstance(other, Stack) and other.axis == self.axis:
      return Stack([func1 if func2 is None else func2 if func1 is None else Add([func1, func2]) for func1, func2 in zip(self.funcs, other.funcs)], self.axis)

  def _multiply(self, other):
    return Stack([Multiply([func, Get(other, self.axis, ifunc)]) if func is not None else None for ifunc, func in enumerate(self.funcs)], self.axis)

  def _dot(self, other, axes):
    newaxis = self.axis
    newaxes = []
    for ax in axes:
      if ax < self.axis:
        newaxis -= 1
        newaxes.append( ax )
      elif ax > self.axis:
        newaxes.append( ax-1 )
    return util.sum(Dot([self.funcs[ifunc], Get(other, self.axis, ifunc)], newaxes) for ifunc in self.nz) if len(newaxes) < len(axes) \
      else Stack([Dot([func, Get(other, self.axis, ifunc)], newaxes) if func is not None else None for ifunc, func in enumerate(self.funcs)], newaxis)

  def _sum(self, axis):
    if axis == self.axis:
      if any(func is not None and func.dtype == bool for func in self.funcs):
        raise NotImplementedError
      return util.sum(func for func in self.funcs if func is not None)
    return Stack([Sum(func, axis-(axis>self.axis)) if func is not None else None for func in self.funcs], self.axis-(axis<self.axis))

  def _transpose(self, axes):
    newaxis = axes.index(self.axis)
    newaxes = [ax-(ax>self.axis) for ax in axes if ax != self.axis]
    return Stack([Transpose(func, newaxes) if func is not None else None for func in self.funcs], newaxis)

  def _takediag(self, axis, rmaxis):
    if self.axis == rmaxis:
      return Stack([Get(func, axis, ifunc) if func is not None else None for ifunc, func in enumerate(self.funcs)], axis)
    elif self.axis == axis:
      return Stack([Get(func, rmaxis-1, ifunc) if func is not None else None for ifunc, func in enumerate(self.funcs)], axis)
    else:
      return Stack([TakeDiag(func, axis-(self.axis<axis), rmaxis-(self.axis<rmaxis)) if func is not None else None for func in self.funcs], self.axis-(rmaxis<self.axis))

  def _take(self, index, axis):
    if axis != self.axis:
      return Stack([Take(func, index, axis-(axis>self.axis)) if func is not None else None for func in self.funcs], self.axis)
    # TODO select axis in index

  def _power(self, n):
    return Stack([Power(func, Get(n, self.axis, ifunc)) if func is not None else None for ifunc, func in enumerate(self.funcs)], self.axis)

  def _mask(self, maskvec, axis):
    if axis != self.axis:
      return Stack([Mask(func, maskvec, axis-(axis>self.axis)) if func is not None else None for func in self.funcs], self.axis)
    newlength = maskvec.sum()
    funcs = [func for ifunc, func in enumerate(self.funcs) if maskvec[ifunc]]
    if all(func is None for func in funcs):
      return Zeros(self.shape[:axis]+(len(funcs),)+self.shape[axis+1:], self.dtype)
    return Stack(funcs, self.axis)

  def _insertaxis(self, axis, length):
    return Stack([InsertAxis(func, axis-(axis>self.axis), length) if func is not None else None for func in self.funcs], self.axis+(self.axis>=axis))

  def _product(self):
    if self.axis == self.ndim-1:
      if len(self.nz) < len(self.funcs):
        return Zeros(self.shape[:-1], self.dtype)
      return util.product(self.funcs)

class DerivativeTargetBase( Array ):
  'base class for derivative targets'

  @property
  def isconstant( self ):
    return False

class Argument(DerivativeTargetBase):
  '''Array argument, to be substituted before evaluation.

  The :class:`Argument` is an :class:`Array` with a known shape, but whose
  values are to be defined later, before evaluation, e.g. using
  :func:`replace_arguments`.

  It is possible to take the derivative of an :class:`Array` to an
  :class:`Argument`:

  >>> from nutils import function
  >>> a = function.Argument('x', [])
  >>> b = function.Argument('y', [])
  >>> f = a**3 + b**2
  >>> function.derivative(f, a).simplified == (3.*a**2).simplified
  True

  Furthermore, derivatives to the local cooardinates are remembered and applied
  to the replacement when using :func:`replace_arguments`:

  >>> from nutils import mesh
  >>> domain, x = mesh.rectilinear([2,2])
  >>> basis = domain.basis('spline', degree=2)
  >>> c = function.Argument('c', basis.shape)
  >>> replace_arguments(c.grad(x), dict(c=basis)) == basis.grad(x)
  True

  Args
  ----
  name : :class:`str`
      The Identifier of this argument.
  shape : :class:`tuple` of :class:`int`\\s
      The shape of this argument.
  nderiv : :class:`int`, non-negative
      Number of times a derivative to the local coordinates is taken.  Default:
      ``0``.
  '''

  def __init__(self, name, shape:tuple, nderiv:int=0):
    self._name = name
    self._nderiv = nderiv
    super().__init__(args=[EVALARGS], shape=shape, dtype=float)

  def evalf(self, evalargs):
    assert self._nderiv == 0
    try:
      value = evalargs[self._name]
    except KeyError:
      raise ValueError('argument {!r} missing'.format(self._name))
    else:
      assert numeric.isarray(value) and value.shape == self.shape
      return value[_]

  def _derivative(self, var, seen):
    if isinstance(var, Argument) and var._name == self._name:
      assert var._nderiv == 0 and self.shape[:self.ndim-self._nderiv] == var.shape
      if self._nderiv:
        return zeros(self.shape+var.shape)
      result = _inflate_scalar(1., self.shape)
      for i, sh in enumerate(self.shape):
        result = diagonalize(result, i, i+self.ndim)
      return result
    elif isinstance(var, LocalCoords):
      return Argument(self._name, self.shape+var.shape, self._nderiv+1)
    else:
      return zeros(self.shape+var.shape)

  def __str__(self):
    return '{} {!r} <{}>'.format(self.__class__.__name__, self._name, ','.join(map(str, self.shape)))

class LocalCoords( DerivativeTargetBase ):
  'local coords derivative target'

  def __init__(self, ndims:int):
    super().__init__(args=[], shape=[ndims], dtype=float)

  def evalf( self ):
    raise Exception( 'LocalCoords should not be evaluated' )

class Ravel( Array ):

  def __init__(self, func:asarray, axis:int):
    assert 0 <= axis < func.ndim-1
    self.func = func
    self.axis = axis
    super().__init__(args=[func], shape=func.shape[:axis]+(func.shape[axis]*func.shape[axis+1],)+func.shape[axis+2:], dtype=func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    if func.shape[self.axis] == 1:
      return get(func, self.axis, 0 ).simplified
    if func.shape[self.axis+1] == 1:
      return get(func, self.axis+1, 0 ).simplified
    retval = func._ravel(self.axis)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Ravel(func, self.axis)

  def evalf( self, f ):
    return f.reshape( f.shape[:self.axis+1] + (f.shape[self.axis+1]*f.shape[self.axis+2],) + f.shape[self.axis+3:] )

  def _multiply(self, other):
    if isinstance(other, Ravel) and other.axis == self.axis and other.func.shape[self.axis:self.axis+2] == self.func.shape[self.axis:self.axis+2]:
      return Ravel(Multiply([self.func, other.func]), self.axis)
    return Ravel(Multiply([self.func, Unravel(other, self.axis, self.func.shape[self.axis:self.axis+2])]), self.axis)

  def _add(self, other):
    if isinstance(other, Ravel) and other.axis == self.axis and other.func.shape[self.axis:self.axis+2] == self.func.shape[self.axis:self.axis+2]:
      return Ravel(Add([self.func, other.func]), self.axis)
    return Ravel(Add([self.func, Unravel(other, self.axis, self.func.shape[self.axis:self.axis+2])]), self.axis)

  def _get(self, i, item):
    if i != self.axis:
      return Ravel(Get(self.func, i+(i>self.axis), item), self.axis-(i<self.axis))
    if item.isconstant and numeric.isint(self.func.shape[self.axis+1]):
      item, = item.eval()
      i, j = divmod(item, self.func.shape[self.axis+1])
      return Get(Get(self.func, self.axis, i), self.axis, j)

  def _dot( self, other, axes ):
    newaxes = [ax+(ax>=self.axis) for ax in axes]
    if self.axis in axes:
      return Dot([self.func, Unravel(other, self.axis, self.func.shape[self.axis:self.axis+2])], sorted(newaxes+[self.axis]))
    else:
      return Ravel(Dot([self.func, Unravel(other, self.axis, self.func.shape[self.axis:self.axis+2])], newaxes), self.axis-builtins.sum(ax<self.axis for ax in axes))

  def _sum(self, axis):
    if axis == self.axis:
      return Sum(Sum(self.func, axis), axis)
    return Ravel(Sum(self.func, axis+(axis>self.axis) ), self.axis-(axis<self.axis))

  def _derivative(self, var, seen):
    return ravel(derivative(self.func, var, seen), axis=self.axis)

  def _transpose(self, axes):
    ravelaxis = axes.index(self.axis)
    funcaxes = [ax+(ax>self.axis) for ax in axes]
    funcaxes = funcaxes[:ravelaxis+1] + [self.axis+1] + funcaxes[ravelaxis+1:]
    return Ravel(Transpose(self.func, funcaxes), ravelaxis)

  def _kronecker(self, axis, length, pos):
    return Ravel(kronecker(self.func, axis+(axis>self.axis), length, pos), self.axis+(axis<=self.axis))

  def _takediag(self, axis, rmaxis):
    if not {self.axis, self.axis+1} & {axis, rmaxis}:
      return Ravel(TakeDiag(self.func, axis+(axis>self.axis), rmaxis+(rmaxis>self.axis)), self.axis-(self.axis>rmaxis))

  def _diagonalize(self, axis, newaxis):
    if axis != self.axis:
      return Ravel(Diagonalize(self.func, axis+(axis>self.axis), newaxis+(newaxis>self.axis)), self.axis+(self.axis>=newaxis))

  def _take(self, index, axis):
    if axis not in (self.axis, self.axis+1):
      return Ravel(Take(self.func, index, axis+(axis>self.axis)), self.axis)

  def _unravel( self, axis, shape ):
    if axis != self.axis:
      return Ravel(Unravel(self.func, axis+(axis>self.axis), shape), self.axis+(self.axis>axis))
    elif shape == self.func.shape[axis:axis+2]:
      return self.func

  def _insertaxis(self, axis, length):
    return Ravel(InsertAxis(self.func, axis+(axis>self.axis), length), self.axis+(axis<=self.axis))

  def _mask(self, maskvec, axis):
    if axis != self.axis:
      return Ravel(Mask(self.func, maskvec, axis+(axis>self.axis)), self.axis)

  @property
  def blocks(self):
    for ind, f in self.func.blocks:
      newind = ravel(ind[self.axis][:,_] * self.func.shape[self.axis+1] + ind[self.axis+1][_,:], axis=0)
      yield (ind[:self.axis] + (newind,) + ind[self.axis+2:]), ravel(f, axis=self.axis)

class Unravel( Array ):

  def __init__(self, func:asarray, axis:int, shape:tuple):
    assert 0 <= axis < func.ndim
    assert func.shape[axis] == numpy.product(shape)
    assert len(shape) == 2
    self.func = func
    self.axis = axis
    self.unravelshape = shape
    super().__init__(args=[func]+[asarray(sh) for sh in shape], shape=func.shape[:axis]+shape+func.shape[axis+1:], dtype=func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    if self.shape[self.axis] == 1:
      return InsertAxis(func, self.axis, 1).simplified
    if self.shape[self.axis+1] == 1:
      return InsertAxis(func, self.axis+1, 1).simplified
    retval = func._unravel(self.axis, self.unravelshape)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Unravel(func, self.axis, self.unravelshape)

  def _derivative(self, var, seen):
    return unravel(derivative(self.func, var, seen), axis=self.axis, shape=self.unravelshape)

  def evalf(self, f, sh1, sh2):
    sh1, = sh1
    sh2, = sh2
    return f.reshape(f.shape[:self.axis+1]+(sh1, sh2)+f.shape[self.axis+2:])

  def _ravel(self, axis):
    if axis == self.axis:
      return self.func

class Mask( Array ):

  def __init__(self, func:asarray, mask:numeric.const, axis:int):
    assert len(mask) == func.shape[axis]
    self.func = func
    self.axis = axis
    self.mask = mask
    super().__init__(args=[func], shape=func.shape[:axis]+(mask.sum(),)+func.shape[axis+1:], dtype=func.dtype)

  @cache.property
  def simplified(self):
    func = self.func.simplified
    if self.mask.all():
      return func
    if not self.mask.any():
      return zeros_like(self)
    retval = func._mask(self.mask, self.axis)
    if retval is not None:
      assert retval.shape == self.shape
      return retval.simplified
    return Mask(func, self.mask, self.axis)

  def evalf( self, func ):
    return func[(slice(None),)*(self.axis+1)+(self.mask,)]

  def _derivative(self, var, seen):
    return mask(derivative(self.func, var, seen), self.mask, self.axis)

  def _get(self, i, item):
    if i != self.axis:
      return Mask(Get(self.func, i, item), self.mask, self.axis-(i<self.axis))
    if item.isconstant:
      item, = item.eval()
      return Get(self.func, i, numpy.arange(len(self.mask))[self.mask][item])

  def _sum(self, axis):
    if axis != self.axis:
      return Mask(sum(self.func, axis), self.mask, self.axis-(axis<self.axis))
    if self.shape[axis] == 1:
      (item,), = self.mask.nonzero()
      return Get(self.func, axis, item)

  def _take(self, index, axis):
    if axis != self.axis:
      return Mask(Take(self.func, index, axis), self.mask, self.axis)

  def _product(self):
    if self.axis != self.ndim-1:
      return Mask(Product(self.func), self.mask, self.axis)

  def _mask(self, maskvec, axis):
    if axis == self.axis:
      newmask = numpy.zeros(len(self.mask), dtype=bool)
      newmask[self.mask] = maskvec
      return Mask(self.func, newmask, self.axis)

  def _takediag(self, axis, rmaxis):
    if self.axis not in (axis, rmaxis):
      return Mask(TakeDiag(self.func, axis, rmaxis), self.mask, self.axis-(rmaxis<self.axis))

class FindTransform(Array):

  def __init__(self, transforms:tuple, trans):
    self.transforms = transforms
    bits = []
    bit = 1
    while bit <= len(transforms):
      bits.append(bit)
      bit <<= 1
    self.bits = numpy.array(bits[::-1])
    super().__init__(args=[trans], shape=(), dtype=int)

  def asdict(self, values):
    assert len(self.transforms) == len(values)
    return dict(zip(self.transforms, values))

  def evalf(self, trans):
    n = len(self.transforms)
    index = 0
    for bit in self.bits:
      i = index|bit
      if i <= n and trans >= self.transforms[i-1]:
        index = i
    index -= 1
    if index < 0 or trans[:len(self.transforms[index])] != self.transforms[index]:
      index = len(self.transforms)
    return numpy.array(index)[_]

class Range(Array):

  def __init__(self, length:asarray, offset:asarray=Zeros((), int)):
    assert length.ndim == 0 and length.dtype == int
    assert offset.ndim == 0 and offset.dtype == int
    self.length = length
    self.offset = offset
    super().__init__(args=[length, offset], shape=[length], dtype=int)

  def evalf(self, length, offset):
    length, = length
    offset, = offset
    return numpy.arange(offset, offset+length)[_]

# AUXILIARY FUNCTIONS (FOR INTERNAL USE)

_ascending = lambda arg: numpy.greater(numpy.diff(arg), 0).all()
_normdims = lambda ndim, shapes: tuple( numeric.normdim(ndim,sh) for sh in shapes )

def _jointdtype( *dtypes ):
  'determine joint dtype'

  type_order = bool, int, int, float
  kind_order = 'biuf'
  itype = builtins.max( kind_order.index(dtype.kind) if isinstance(dtype,numpy.dtype)
           else type_order.index(dtype) for dtype in dtypes )
  return type_order[itype]

def _matchndim( *arrays ):
  'introduce singleton dimensions to match ndims'

  arrays = [ asarray(array) for array in arrays ]
  ndim = builtins.max( array.ndim for array in arrays )
  return tuple(array[(_,)*(ndim-array.ndim)] for array in arrays)

def _obj2str( obj ):
  'convert object to string'

  if numeric.isarray(obj):
    if obj.size < 6:
      return _obj2str(obj.tolist())
    return 'array<%s>' % 'x'.join( str(n) for n in obj.shape )
  if isinstance( obj, list ):
    if len(obj) < 6:
      return '[%s]' % ','.join( _obj2str(o) for o in obj )
    return '[#%d]' % len(obj)
  if isinstance( obj, (tuple,set) ):
    if len(obj) < 6:
      return '(%s)' % ','.join( _obj2str(o) for o in obj )
    return '(#%d)' % len(obj)
  if isinstance(obj, collections.abc.Mapping):
    return '{#%d}' % len(obj)
  if isinstance( obj, slice ):
    I = ''
    if obj.start is not None:
      I += str(obj.start)
    if obj.step is not None:
      I += ':' + str(obj.step)
    I += ':'
    if obj.stop is not None:
      I += str(obj.stop)
    return I
  if obj is Ellipsis:
    return '...'
  return str(obj)

def _invtrans(trans):
  trans = numpy.asarray(trans)
  assert trans.dtype == int
  invtrans = numpy.empty(len(trans), dtype=int)
  invtrans[trans] = numpy.arange(len(trans))
  return tuple(invtrans)

def _norm_and_sort( ndim, args ):
  'norm axes, sort, and assert unique'

  normargs = tuple( sorted( numeric.normdim( ndim, arg ) for arg in args ) )
  assert _ascending( normargs ) # strict
  return normargs

def _concatblocks(items):
  gathered = util.gather(items)
  order = [ind for ind12, ind_f in gathered for ind, f in ind_f]
  blocks = []
  for (ind1, ind2), ind_f in gathered:
    if len(ind_f) == 1:
      ind, f = ind_f[0]
    else:
      inds, fs = zip(*sorted(ind_f, key=lambda item: order.index(item[0])))
      ind = Concatenate(inds, axis=0)
      f = Concatenate(fs, axis=len(ind1))
    blocks.append(((ind1+(ind,)+ind2), f))
  return tuple(blocks)

def _numpy_align(*arrays):
  '''reshape arrays according to Numpy's broadcast conventions'''
  arrays = [asarray(array) for array in arrays]
  if len(arrays) > 1:
    ndim = builtins.max([array.ndim for array in arrays])
    for idim in range(ndim):
      lengths = [array.shape[idim] for array in arrays if array.ndim == ndim and array.shape[idim] != 1]
      length = lengths[0] if lengths else 1
      assert all(l == length for l in lengths), 'incompatible shapes: {}'.format(' != '.join(str(l) for l in lengths))
      for i, a in enumerate(arrays):
        if a.ndim < ndim:
          arrays[i] = insertaxis(a, idim, length)
        elif a.shape[idim] != length:
          arrays[i] = repeat(a, length, idim)
  return arrays

def _inflate_scalar(arg, shape):
  arg = asarray(arg)
  assert arg.ndim == 0
  for idim, length in enumerate(shape):
    arg = insertaxis(arg, idim, length)
  return arg

# FUNCTIONS

isarray = lambda arg: isinstance( arg, Array )
iszero = lambda arg: isinstance( arg, Zeros )
zeros = lambda shape, dtype=float: Zeros( shape, dtype )
zeros_like = lambda arr: zeros(arr.shape, arr.dtype)
ones = lambda shape, dtype=float: _inflate_scalar(numpy.ones((), dtype=dtype), shape)
ones_like = lambda arr: ones(arr.shape, arr.dtype)
reciprocal = lambda arg: power( arg, -1 )
grad = lambda arg, coords, ndims=0: asarray( arg ).grad( coords, ndims )
symgrad = lambda arg, coords, ndims=0: asarray( arg ).symgrad( coords, ndims )
div = lambda arg, coords, ndims=0: asarray( arg ).div( coords, ndims )
negative = lambda arg: multiply( arg, -1 )
nsymgrad = lambda arg, coords: ( symgrad(arg,coords) * coords.normal() ).sum(-1)
ngrad = lambda arg, coords: ( grad(arg,coords) * coords.normal() ).sum(-1)
sin = lambda x: Sin(x)
cos = lambda x: Cos(x)
rotmat = lambda arg: Stack([trignormal(arg), trigtangent(arg)], 0)
tan = lambda x: Tan(x)
arcsin = lambda x: ArcSin(x)
arccos = lambda x: ArcCos(x)
exp = lambda x: Exp(x)
ln = lambda x: Log(x)
mod = lambda arg1, arg2: Mod(*_numpy_align(arg1, arg2))
log2 = lambda arg: ln(arg) / ln(2)
log10 = lambda arg: ln(arg) / ln(10)
sqrt = lambda arg: power( arg, .5 )
arctan2 = lambda arg1, arg2: ArcTan2(*_numpy_align(arg1, arg2))
greater = lambda arg1, arg2: Greater(*_numpy_align(arg1, arg2))
equal = lambda arg1, arg2: Equal(*_numpy_align(arg1, arg2))
less = lambda arg1, arg2: Less(*_numpy_align(arg1, arg2))
min = lambda a, b: Minimum(*_numpy_align(a, b))
max = lambda a, b: Maximum(*_numpy_align(a, b))
abs = lambda arg: arg * sign(arg)
sinh = lambda arg: .5 * ( exp(arg) - exp(-arg) )
cosh = lambda arg: .5 * ( exp(arg) + exp(-arg) )
tanh = lambda arg: 1 - 2. / ( exp(2*arg) + 1 )
arctanh = lambda arg: .5 * ( ln(1+arg) - ln(1-arg) )
piecewise = lambda level, intervals, *funcs: Get(Stack(asarrays(funcs), axis=0), axis=0, item=util.sum(Int(greater(level, interval)) for interval in intervals))
trace = lambda arg, n1=-2, n2=-1: sum(takediag(arg, n1, n2), numeric.normdim(arg.ndim, n1))
normalized = lambda arg, axis=-1: divide(arg, expand_dims(norm2(arg, axis=axis), axis))
norm2 = lambda arg, axis=-1: sqrt( sum( multiply( arg, arg ), axis ) )
heaviside = lambda arg: Int(greater(arg, 0))
divide = lambda arg1, arg2: multiply( arg1, reciprocal(arg2) )
subtract = lambda arg1, arg2: add( arg1, negative(arg2) )
mean = lambda arg: .5 * ( arg + opposite(arg) )
jump = lambda arg: opposite(arg) - arg
add_T = lambda arg, axes=(-2,-1): swapaxes( arg, axes ) + arg
blocks = lambda arg: asarray(arg).simplified.blocks
rootcoords = lambda ndims: RootCoords( ndims )
sampled = lambda data, ndims: Sampled( data )
opposite = cache.replace(initcache={TRANS: OPPTRANS, OPPTRANS: TRANS})
bifurcate1 = cache.replace(initcache={TRANS: SelectChain(TRANS, True), OPPTRANS: SelectChain(OPPTRANS, True)})
bifurcate2 = cache.replace(initcache={TRANS: SelectChain(TRANS, False), OPPTRANS: SelectChain(OPPTRANS, False)})
bifurcate = lambda arg1, arg2: ( bifurcate1(arg1), bifurcate2(arg2) )
curvature = lambda geom, ndims=-1: geom.normal().div(geom, ndims=ndims)
laplace = lambda arg, geom, ndims=0: arg.grad(geom, ndims).div(geom, ndims)
symgrad = lambda arg, geom, ndims=0: multiply(.5, add_T(arg.grad(geom, ndims)))
div = lambda arg, geom, ndims=0: trace(arg.grad(geom, ndims))
tangent = lambda geom, vec: subtract(vec, multiply(dot(vec, normal(geom), -1)[...,_], normal(geom)))
ngrad = lambda arg, geom, ndims=0: dotnorm(grad(arg, geom, ndims), geom)
nsymgrad = lambda arg, geom, ndims=0: dotnorm(symgrad(arg, geom, ndims), geom)
expand_dims = lambda arg, n: InsertAxis(arg, n, 1)

def trignormal( angle ):
  angle = asarray( angle )
  assert angle.ndim == 0
  if iszero( angle ):
    return kronecker( 1, axis=0, length=2, pos=0 )
  return TrigNormal( angle )

def trigtangent( angle ):
  angle = asarray( angle )
  assert angle.ndim == 0
  if iszero( angle ):
    return kronecker( 1, axis=0, length=2, pos=1 )
  return TrigTangent( angle )

eye = lambda n, dtype=float: diagonalize(ones([n], dtype=dtype))

def insertaxis(arg, n, length):
  arg = asarray(arg)
  n = numeric.normdim(arg.ndim+1, n)
  return InsertAxis(arg, n, length)

stack = lambda args, axis=0: Stack(_numpy_align(*args), axis)

def chain( funcs ):
  'chain'

  funcs = [ asarray(func) for func in funcs ]
  shapes = [ func.shape[0] for func in funcs ]
  return [ concatenate( [ func if i==j else zeros( (sh,) + func.shape[1:] )
             for j, sh in enumerate(shapes) ], axis=0 )
               for i, func in enumerate(funcs) ]

vectorize = lambda args: concatenate([kronecker(arg, axis=-1, length=len(args), pos=iarg) for iarg, arg in enumerate(args)])

def repeat(arg, length, axis):
  arg = asarray(arg)
  assert arg.shape[axis] == 1
  return insertaxis(get(arg, axis, 0), axis, length)

def get(arg, iax, item):
  arg = asarray(arg)
  item = asarray(item)
  iax = numeric.normdim(arg.ndim, iax)
  sh = arg.shape[iax]
  if numeric.isint(sh) and item.isconstant:
    item = numeric.normdim(sh, item.eval()[0])
  return Get(arg, iax, item)

def align(arg, axes, ndim):
  keep = numpy.zeros(ndim, dtype=bool)
  keep[list(axes)] = True
  renumber = keep.cumsum()-1
  transaxes = _invtrans(renumber[numpy.asarray(axes)])
  retval = transpose(arg, transaxes)
  for axis in numpy.where(~keep)[0]:
    retval = expand_dims(retval, axis)
  for i, j in enumerate(axes):
    assert arg.shape[i] == retval.shape[j]
  return retval

def bringforward( arg, axis ):
  'bring axis forward'

  arg = asarray(arg)
  axis = numeric.normdim(arg.ndim,axis)
  if axis == 0:
    return arg
  return transpose( args, [axis] + range(axis) + range(axis+1,args.ndim) )

def jacobian( geom, ndims ):
  assert geom.ndim == 1
  J = localgradient( geom, ndims )
  cndims, = geom.shape
  assert J.shape == (cndims,ndims), 'wrong jacobian shape: got %s, expected %s' % ( J.shape, (cndims, ndims) )
  assert cndims >= ndims, 'geometry dimension < topology dimension'
  detJ = abs( determinant( J ) ) if cndims == ndims \
    else 1. if ndims == 0 \
    else abs( determinant( ( J[:,:,_] * J[:,_,:] ).sum(0) ) )**.5
  return detJ

def matmat( arg0, *args ):
  'helper function, contracts last axis of arg0 with first axis of arg1, etc'
  retval = asarray( arg0 )
  for arg in args:
    arg = asarray( arg )
    assert retval.shape[-1] == arg.shape[0], 'incompatible shapes'
    retval = dot( retval[(...,)+(_,)*(arg.ndim-1)], arg[(_,)*(retval.ndim-1)], retval.ndim-1 )
  return retval

def determinant(arg, axes=(-2,-1)):
  arg = asarray(arg)
  ax1, ax2 = _norm_and_sort(arg.ndim, axes)
  assert ax2 > ax1 # strict
  trans = [i for i in range(arg.ndim) if i not in (ax1, ax2)] + [ax1, ax2]
  arg = transpose(arg, trans)
  return Determinant(arg)

def inverse(arg, axes=(-2,-1)):
  arg = asarray( arg )
  ax1, ax2 = _norm_and_sort(arg.ndim, axes)
  assert ax2 > ax1 # strict
  trans = [i for i in range(arg.ndim) if i not in (ax1, ax2)] + [ax1, ax2]
  arg = transpose(arg, trans)
  return transpose(Inverse(arg), _invtrans(trans))

def takediag(arg, axis=-2, rmaxis=-1):
  arg = asarray(arg)
  axis = numeric.normdim(arg.ndim, axis)
  rmaxis = numeric.normdim(arg.ndim, rmaxis)
  assert axis < rmaxis
  return TakeDiag(arg, axis, rmaxis)

def derivative(func, var, seen=None):
  'derivative'

  assert isinstance(var, DerivativeTargetBase), 'invalid derivative target {!r}'.format(var)
  if seen is None:
    seen = {}
  func = asarray(func)
  if func in seen:
    result = seen[func]
  else:
    result = func._derivative(var, seen)
    seen[func] = result
  assert result.shape == func.shape+var.shape, 'bug in {}._derivative'.format(func)
  return result

def localgradient(arg, ndims):
  'local derivative'

  return derivative(arg, LocalCoords(ndims))

def dotnorm( arg, coords ):
  'normal component'

  return sum( arg * coords.normal(), -1 )

normal = lambda geom: geom.normal()

def kronecker(arg, axis, length, pos):
  arg = asarray(arg)
  axis = numeric.normdim(arg.ndim+1, axis)
  funcs = [None] * length
  funcs[pos] = arg
  return Stack(funcs, axis)

def diagonalize(arg, axis=-1, newaxis=-1):
  arg = asarray(arg)
  axis = numeric.normdim(arg.ndim, axis)
  newaxis = numeric.normdim(arg.ndim+1, newaxis)
  assert axis < newaxis
  return Diagonalize(arg, axis, newaxis)

def concatenate(args, axis=0):
  args = _matchndim(*args)
  axis = numeric.normdim(args[0].ndim, axis)
  return Concatenate(args, axis)

def cross(arg1, arg2, axis):
  arg1, arg2 = _numpy_align(arg1, arg2)
  axis = numeric.normdim(arg1.ndim, axis)
  assert arg1.shape[axis] == 3
  return Cross(arg1, arg2, axis)

def outer( arg1, arg2=None, axis=0 ):
  'outer product'

  if arg2 is not None and arg1.ndim != arg2.ndim:
    warnings.warn( 'varying ndims in function.outer; this will be forbidden in future', DeprecationWarning )
  arg1, arg2 = _matchndim( arg1, arg2 if arg2 is not None else arg1 )
  axis = numeric.normdim( arg1.ndim, axis )
  return expand_dims(arg1,axis+1) * expand_dims(arg2,axis)

def sign(arg):
  arg = asarray(arg)
  return Sign(arg)

def eig(arg, axes=(-2,-1), symmetric=False):
  arg = asarray(arg)
  ax1, ax2 = _norm_and_sort( arg.ndim, axes )
  assert ax2 > ax1 # strict
  trans = [i for i in range(arg.ndim) if i not in (ax1, ax2)] + [ax1, ax2]
  transposed = transpose(arg, trans)
  eigval, eigvec = Eig(transposed, symmetric)
  return Tuple([transpose(diagonalize(eigval), _invtrans(trans)), transpose(eigvec, _invtrans(trans))])

def function(fmap, nmap, ndofs):
  transforms = sorted(fmap)
  depth, = set(len(transform) for transform in transforms)
  fromdims, = set(transform.fromdims for transform in transforms)
  promote = Promote(fromdims, trans=TRANS)
  index = FindTransform(transforms, promote)
  dofmap = DofMap([nmap[trans] for trans in transforms], index=index)
  func = Function(stds=[fmap[trans] for trans in transforms], depth=depth, trans=promote, index=index)
  return Inflate(func, dofmap, ndofs, axis=0)

def elemwise( fmap, shape, default=None ):
  return Elemwise( fmap=fmap, shape=shape, default=default )

def take(arg, index, axis):
  arg = asarray(arg)
  axis = numeric.normdim(arg.ndim, axis)
  index = asarray(index)
  assert index.ndim == 1
  if index.dtype == bool:
    assert index.shape[0] == arg.shape[axis]
    if index.isconstant:
      mask, = index.eval()
      return Mask(arg, mask, axis)
    index = find(index)
  return Take(arg, index, axis)

def find( arg ):
  'find'

  arg = asarray( arg )
  assert arg.ndim == 1 and arg.dtype == bool

  if arg.isconstant:
    arg, = arg.eval()
    index, = arg.nonzero()
    return asarray( index )

  return Find( arg )

def inflate(arg, dofmap, length, axis):
  arg = asarray(arg)
  dofmap = asarray(dofmap)
  axis = numeric.normdim(arg.ndim, axis)
  shape = arg.shape[:axis] + (length,) + arg.shape[axis+1:]
  if dofmap.isconstant:
    n = arg.shape[axis]
    assert numeric.isint(n), 'constant inflation only allowed over fixed-length axis'
    index, = dofmap.eval()
    assert len(index) == n
    assert numpy.greater_equal(index, 0).all() and numpy.less(index, length).all()
    assert numpy.equal(numpy.diff(index), 1).all(), 'constant inflation must be contiguous'
    if n == length:
      retval = arg
    else:
      parts = []
      if index[0] > 0:
        parts.append( zeros( arg.shape[:axis] + (index[0],) + arg.shape[axis+1:], dtype=arg.dtype ) )
      parts.append( arg )
      if index[0] + n < length:
        parts.append( zeros( arg.shape[:axis] + (length-index[0]-n,) + arg.shape[axis+1:], dtype=arg.dtype ) )
      retval = concatenate( parts, axis=axis )
    assert retval.shape == tuple(shape)
    return retval
  return Inflate(arg, dofmap, length, axis)

def mask(arg, mask, axis=0):
  arg = asarray(arg)
  axis = numeric.normdim(arg.ndim, axis)
  assert numeric.isarray(mask) and mask.ndim == 1 and mask.dtype == bool
  assert arg.shape[axis] == len(mask)
  return Mask(arg, mask, axis)

def J( geometry, ndims=None ):
  if ndims is None:
    ndims = len(geometry)
  elif ndims < 0:
    ndims += len(geometry)
  return jacobian( geometry, ndims )

def unravel(func, axis, shape):
  func = asarray(func)
  axis = numeric.normdim(func.ndim, axis)
  shape = tuple(shape)
  assert func.shape[axis] == numpy.product(shape)
  return Unravel(func, axis, tuple(shape))

def ravel(func, axis):
  func = asarray(func)
  axis = numeric.normdim( func.ndim-1, axis )
  return Ravel(func, axis)

@cache.replace
def replace_arguments(value, arguments):
  '''Replace :class:`Argument` objects in ``value``.

  Replace :class:`Argument` objects in ``value`` according to the ``arguments``
  map, taking into account derivatives to the local coordinates.

  Args
  ----
  value : :class:`Array`
      Array to be edited.
  arguments : :class:`collections.abc.Mapping` with :class:`Array`\\s as values
      :class:`Argument`\\s replacements.  The key correspond to the ``name``
      passed to an :class:`Argument` and the value is the replacement.

  Returns
  -------
  :class:`Array`
      The edited ``value``.
  '''
  if isinstance(value, Argument) and value._name in arguments:
    v = asarray(arguments[value._name])
    assert value.shape[:value.ndim-value._nderiv] == v.shape
    for ndims in value.shape[value.ndim-value._nderiv:]:
      v = localgradient(v, ndims)
    return v

@cache.replace
def zero_argument_derivatives(arg):
  if isinstance(arg, Argument) and arg._nderiv > 0:
    return zeros_like(arg)

def _eval_ast(ast, functions):
  '''evaluate ``ast`` generated by :func:`nutils.expression.parse`'''

  op, *args = ast
  if op is None:
    value, = args
    return value

  args = (_eval_ast(arg, functions) for arg in args)
  if op == 'group':
    array, = args
    return array
  elif op == 'arg':
    name, *shape = args
    return Argument(name, shape)
  elif op == 'substitute':
    array, arg, value = args
    assert isinstance(arg, Argument) and arg._nderiv == 0
    return replace_arguments(array, {arg._name: value})
  elif op == 'call':
    func, arg = args
    return functions[func](arg)
  elif op == 'eye':
    length, = args
    return eye(length)
  elif op == 'normal':
    geom, = args
    return normal(geom)
  elif op == 'getitem':
    array, dim, index = args
    return get(array, dim, index)
  elif op == 'trace':
    array, n1, n2 = args
    return trace(array, n1, n2)
  elif op == 'sum':
    array, axis = args
    return sum(array, axis)
  elif op == 'concatenate':
    return concatenate(args, axis=0)
  elif op == 'grad':
    array, geom = args
    return grad(array, geom)
  elif op == 'surfgrad':
    array, geom = args
    return grad(array, geom, len(geom)-1)
  elif op == 'derivative':
    func, target = args
    return derivative(func, target)
  elif op == 'append_axis':
    array, length = args
    return repeat(asarray(array)[..., None], length, -1)
  elif op == 'transpose':
    array, trans = args
    return transpose(array, trans)
  elif op == 'jump':
    array, = args
    return jump(array)
  elif op == 'mean':
    array, = args
    return mean(array)
  elif op == 'neg':
    array, = args
    return -asarray(array)
  elif op in ('add', 'sub', 'mul', 'truediv', 'pow'):
    left, right = args
    return getattr(operator, '__{}__'.format(op))(asarray(left), asarray(right))
  else:
    raise ValueError('unknown opcode: {!r}'.format(op))

class Namespace:
  '''Namespace for :class:`Array` objects supporting assignments with tensor expressions.

  The :class:`Namespace` object is used to store :class:`Array` objects.

  >>> from nutils import function
  >>> ns = function.Namespace()
  >>> ns.A = function.zeros([3, 3])
  >>> ns.x = function.zeros([3])
  >>> ns.c = 2

  In addition to the assignment of :class:`Array` objects, it is also possible
  to specify an array using a tensor expression string — see
  :func:`nutils.expression.parse` for the syntax.  All attributes defined in
  this namespace are available as variables in the expression.  If the array
  defined by the expression has one or more dimensions the indices of the axes
  should be appended to the attribute name.  Examples:

  >>> ns.cAx_i = 'c A_ij x_j'
  >>> ns.xAx = 'x_i A_ij x_j'

  It is also possible to simply evaluate an expression without storing its
  value in the namespace by passing the expression to the method ``eval_``
  suffixed with appropriate indices:

  >>> ns.eval_('2 c')
  Array<>
  >>> ns.eval_i('c A_ij x_j')
  Array<3>
  >>> ns.eval_ij('A_ij + A_ji')
  Array<3,3>

  For zero and one dimensional expressions the following shorthand can be used:

  >>> '2 c' @ ns
  Array<>
  >>> 'A_ij x_j' @ ns
  Array<3>

  When evaluating an expression through this namespace the following functions
  are available: ``opposite``, ``sin``, ``cos``, ``tan``, ``sinh``, ``cosh``,
  ``tanh``, ``arcsin``, ``arccos``, ``arctan2``, ``arctanh``, ``exp``, ``abs``,
  ``ln``, ``log``, ``log2``, ``log10``, ``sqrt`` and ``sign``.

  Args
  ----
  default_geometry_name : :class:`str`
      The name of the default geometry.  This argument is passed to
      :func:`nutils.expression.parse`.  Default: ``'x'``.

  Attributes
  ----------
  arg_shapes : :class:`types.MappingProxyType`
      A readonly map of argument names and shapes.
  default_geometry_name : :class:`str`
      The name of the default geometry.  See argument with the same name.
  '''

  __slots__ = '_attributes', '_arg_shapes', 'arg_shapes', 'default_geometry_name'

  _re_assign = re.compile('^([a-zA-Zα-ωΑ-Ω][a-zA-Zα-ωΑ-Ω0-9]*)(_[a-z]+)?$')

  _functions = dict(
    opposite=opposite, sin=sin, cos=cos, tan=tan, sinh=sinh, cosh=cosh,
    tanh=tanh, arcsin=arcsin, arccos=arccos, arctan2=arctan2, arctanh=arctanh,
    exp=exp, abs=abs, ln=ln, log=ln, log2=log2, log10=log10, sqrt=sqrt,
    sign=sign,
  )
  _functions_nargs = {k: len(inspect.signature(v).parameters) for k, v in _functions.items()}

  def __init__(self, *, default_geometry_name='x'):
    if not isinstance(default_geometry_name, str):
      raise ValueError('default_geometry_name: Expected a str, got {!r}.'.format(default_geometry_name))
    if '_' in default_geometry_name or not self._re_assign.match(default_geometry_name):
      raise ValueError('default_geometry_name: Invalid variable name: {!r}.'.format(default_geometry_name))
    super().__setattr__('_attributes', {})
    super().__setattr__('_arg_shapes', {})
    super().__setattr__('arg_shapes', types.MappingProxyType(self._arg_shapes))
    super().__setattr__('default_geometry_name', default_geometry_name)
    super().__init__()

  @property
  def default_geometry(self):
    ''':class:`nutils.function.Array`: The default geometry, shorthand for ``getattr(ns, ns.default_geometry_name)``.'''
    return getattr(self, self.default_geometry_name)

  def __or__(self, subs):
    '''Return a copy with arguments replaced by ``subs``.

    Return a copy of this namespace with :class:`Argument` objects replaced
    according to ``subs``.

    Args
    ----
    subs : :class:`dict` of :class:`str` and :class:`nutils.function.Array` objects
        Replacements of the :class:`Argument` objects, identified by their names.

    Returns
    -------
    ns : :class:`Namespace`
        The copy of this namespace with replaced :class:`Argument` objects.
    '''

    if not isinstance(subs, collections.abc.Mapping):
      return NotImplemented
    ns = Namespace(default_geometry_name=self.default_geometry_name)
    for k, v in self._attributes.items():
      setattr(ns, k, replace_arguments(v, subs))
    return ns

  def copy_(self, *, default_geometry_name=None):
    '''Return a copy of this namespace.'''

    if default_geometry_name is None:
      default_geometry_name = self.default_geometry_name
    ns = Namespace(default_geometry_name=default_geometry_name)
    for k, v in self._attributes.items():
      setattr(ns, k, v)
    return ns

  def __getattr__(self, name):
    '''Get attribute ``name``.'''

    if name.startswith('eval_'):
      return lambda expr: _eval_ast(expression.parse(expr, variables=self._attributes, functions=self._functions_nargs, indices=name[5:], arg_shapes=self._arg_shapes, default_geometry_name=self.default_geometry_name)[0], self._functions)
    try:
      return self._attributes[name]
    except KeyError:
      pass
    raise AttributeError(name)

  def __setattr__(self, name, value):
    '''Set attribute ``name`` to ``value``.'''

    if name in self.__slots__:
      raise AttributeError('readonly')
    m = self._re_assign.match(name)
    if not m or m.group(2) and len(set(m.group(2))) != len(m.group(2)):
      raise AttributeError('{!r} object has no attribute {!r}'.format(type(self), name))
    else:
      name, indices = m.groups()
      indices = indices[1:] if indices else ''
      if isinstance(value, str):
        ast, arg_shapes = expression.parse(value, variables=self._attributes, functions=self._functions_nargs, indices=indices, arg_shapes=self._arg_shapes, default_geometry_name=self.default_geometry_name)
        value = _eval_ast(ast, self._functions)
        self._arg_shapes.update(arg_shapes)
      else:
        assert not indices
      self._attributes[name] = asarray(value)

  def __delattr__(self, name):
    '''Delete attribute ``name``.'''

    if name in self.__slots__:
      raise AttributeError('readonly')
    elif name in self._attributes:
      del self._attributes[name]
    else:
      raise AttributeError('{!r} object has no attribute {!r}'.format(type(self), name))

  def __rmatmul__(self, expr):
    '''Evaluate zero or one dimensional ``expr``.'''

    if not isinstance(expr, str):
      return NotImplemented
    try:
      ast = expression.parse(expr, variables=self._attributes, functions=self._functions_nargs, indices=None, arg_shapes=self._arg_shapes, default_geometry_name=self.default_geometry_name)[0]
    except expression.AmbiguousAlignmentError:
      raise ValueError('`expression @ Namespace` cannot be used because the expression has more than one dimension.  Use `Namespace.eval_...(expression)` instead')
    return _eval_ast(ast, self._functions)

def normal(arg, exterior=False):
  assert arg.ndim == 1
  if not exterior:
    lgrad = localgradient(arg, len(arg))
    return Normal(lgrad)
  lgrad = localgradient(arg, len(arg)-1)
  if len(arg) == 2:
    return asarray([lgrad[1,0], -lgrad[0,0]]).normalized()
  if len(arg) == 3:
    return cross(lgrad[:,0], lgrad[:,1], axis=0).normalized()
  raise NotImplementedError

def grad(self, geom, ndims=0):
  assert geom.ndim == 1
  if ndims <= 0:
    ndims += geom.shape[0]
  J = localgradient(geom, ndims)
  if J.shape[0] == J.shape[1]:
    Jinv = inverse(J)
  elif J.shape[0] == J.shape[1] + 1: # gamma gradient
    G = dot(J[:,:,_], J[:,_,:], 0)
    Ginv = inverse(G)
    Jinv = dot(J[_,:,:], Ginv[:,_,:], -1)
  else:
    raise Exception( 'cannot invert %sx%s jacobian' % J.shape )
  return dot(localgradient(self, ndims)[...,_], Jinv, -2)

def dotnorm(arg, geom, axis=-1):
  axis = numeric.normdim(arg.ndim, axis)
  assert geom.ndim == 1 and geom.shape[0] == arg.shape[axis]
  return dot(arg, normal(geom)[(slice(None),)+(_,)*(arg.ndim-axis-1)], axis)

# vim:shiftwidth=2:softtabstop=2:expandtab:foldmethod=indent:foldnestmax=2
