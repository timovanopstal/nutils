from . import util, element, numpy, function, _

class AddPart( function.ArrayFunc ):
  'add partial'

  needxi = True

  def __init__( self, func1, func2 ):
    'constructor'

    self.args = func1, func2, func2.shape[0]
    self.shape = func1.shape

  @staticmethod
  def eval( xi, func1, func2, I ):
    'evaluate'

    retval = func1.copy()
    retval[ I ] += func2
    return retval

class CacheFuncColoc( object ):
  'cached function evaluation for colocation'

  def __init__( self, coords, *payload ):
    'compare'

    self.data = coords, payload

  def iterdata( self, myelem ):
    'iterate'

    yield self.data
  
class CacheFuncND( object ):
  'cached function evaluation'

  def __init__( self, topo, coords, *payload ):
    'compare'

    self.F = function.Tuple([ coords, function.Tuple( payload ) ])
    self.data = [ (elem,) + self.F(elem.eval('gauss7')) for elem in topo ]

  def iterdata( self, myelem ):
    'iterate'

    for elem, y, w_payload in self.data:
      if elem is myelem:
        y, w_payload = self.F( elem.eval('uniform345') )
      yield y, w_payload

  def __eq__( self, other ):
    'compare'

    return isinstance( other, CacheFuncND ) and self.F == other.F

class Convolution( function.Evaluable ):
  'integrate iterator'

  needxi = True

  def __init__( self, coords, cache ):
    'constructor'

    self.args = coords, cache

  @staticmethod
  def eval( xi, x, cachefunc ):
    'convolute shapes'

    iterdata = []
    for y, w_payload in cachefunc.iterdata( xi.elem ):
      d = x[:,:,_] - y[:,_,:] # FIX count number of axes in x
      r2 = util.contract( d, d, 0 )
      r = numpy.sqrt( r2 )
      logr = .5 * numpy.log( r2 )
      iterdata.append( (d/r,r,r2,logr) + w_payload )
    return xi.points.coords.shape[1:], iterdata

class Convolution3D( function.Evaluable ):
  'integrate iterator'

  needxi = True

  def __init__( self, mycoords, topo, coords, func, *other ):
    'constructor'

    self.args = mycoords, CacheFunc3D( topo, coords, func, other )

  @staticmethod
  def eval( xi, x, cachefunc ):
    'convolute shapes'

    iterdata = []
    for elem, y, funcs in cachefunc.iterdata( xi.elem ):
      d = x[:,:,_] - y[:,_,:] # FIX count number of axes in x
      r2 = util.contract( d, d, 0 )
      r1 = numpy.sqrt( r2 )
      d_r = d / r1
      iterdata.append( (d_r,r1,r2) + funcs )
    return xi.points.coords.shape[1:], iterdata

# LAPLACE

class Laplacelet( function.ArrayFunc ):
  'Laplacelet'

  def __init__( self, mycoords, topo, coords, funcsp ):
    'constructor'

    self.topo = topo
    self.mycoords = mycoords
    self.coords = coords
    self.funcsp = funcsp
    self.shape = int(funcsp.shape[0]),
    self.args = int(funcsp.shape[0]), Convolution( mycoords, topo, coords, funcsp, funcsp.shape[0] )

  @staticmethod
  def eval( ndofs, (shape,iterdata) ):
    'evaluate'

    retval = numpy.zeros( (ndofs,)+shape )
    for D, R2, logR, wf, I in iterdata:
      kernel = logR / (-2*numpy.pi)
      retval[ I ] += numpy.tensordot( wf, kernel, (-1,-1) )
    return retval

  def grad( self ):
    'gradient'

    return LaplaceletGrad( self.mycoords, self.topo, self.coords, self.funcsp )

  def flux( self ):
    'flux'

    return AddPart( ( self.grad() * self.coords.normal() ).sum(), .5 * self.funcsp )

  def reconstruct( self, bval, flux ):
    'reconstruct'

    return LaplaceletReconstruct( self.mycoords, self.topo, self.coords, bval, flux )

class LaplaceletGrad( function.ArrayFunc ):
  'laplacelet gradient'

  def __init__( self, mycoords, topo, coords, funcsp ):
    'constructor'

    self.shape = int(funcsp.shape[0]), 2
    self.args = int(funcsp.shape[0]), Convolution( mycoords, topo, coords, funcsp, funcsp.shape[0] )

  @staticmethod
  def eval( ndofs, (shape,iterdata) ):
    'evaluate'

    retval = numpy.zeros( (ndofs,2) + shape )
    for D, R2, logR, wf, I in iterdata:
      kernel = D / ( R2 * (-2*numpy.pi) )
      retval[ I ] += numpy.tensordot( wf, kernel, (-1,-1) )
    return retval

class LaplaceletReconstruct( function.ArrayFunc ):
  'laplacelet reconstruction'

  def __init__( self, mycoords, topo, coords, bval, flux ):
    'constructor'

    self.shape = ()
    self.args = Convolution( mycoords, topo, coords, 1, bval, flux, coords.normal() ),

  @staticmethod
  def eval( (shape,iterdata) ):
    'evaluate'

    retval = 0
    for D, R2, logR, w, bval, flux, normal in iterdata:
      retval += numpy.dot( logR * flux + util.contract( D, normal[:,_,:], 0 ) * bval / R2, w ) # TODO fix for arbitrary axes
    return retval / (-2*numpy.pi)

# STOKES 2D

class Stokeslet( function.ArrayFunc ):
  'stokeslet'

  def __init__( self, mycoords, topo, coords, funcsp, mu ):
    'constructor'

    self.topo = topo
    self.coords = coords
    self.mycoords = mycoords
    self.funcsp = funcsp
    self.mu = mu
    self.shape = int( funcsp.shape[0] ) * 2, 2

    if isinstance( funcsp, function.Evaluable ):
      iweights = coords.iweights( topo.ndims ) * funcsp
      cache = CacheFuncND( topo, coords, iweights, funcsp.shape[0] )
    else:
      assert isinstance( funcsp, numpy.ndarray )
      cache = CacheFuncColoc( funcsp.T, numpy.eye(funcsp.shape[0]), slice(None) )

    self.args = int(funcsp.shape[0]), funcsp.shape[0], Convolution( mycoords, cache ), mu

  @staticmethod
  def eval( ndofs, N, (shape,iterdata), mu ):
    'evaluate'

    retval = numpy.zeros( (ndofs*2,2) + shape )
    retval_swap = retval.reshape( 2, ndofs, 2, *shape ).swapaxes(0,1) # follows ordering Vectorize
    for D_R, R, R2, logR, w, M in iterdata:
      kernel = D_R[:,_] * D_R[_,:]
      util.ndiag( kernel, [0,1] )[:] -= logR
      retval_swap[M] += numpy.tensordot( w, kernel, (-1,-1) )
    retval /= 4 * numpy.pi * mu
    return retval

  def traction( self ):
    'flux'

    return StokesletTrac( self.mycoords, self.topo, self.coords, self.funcsp )

  def curvature( self ):
    'curvature'

    grad = StokesletGrad( self.mycoords, self.topo, self.coords, self.funcsp, self.mu )
    normal = self.coords.normal()
    return grad.trace(1,2) - ( ( grad * normal ).sum() * normal ).sum()

  def reconstruct( self, velo, trac, surftens=0 ):
    'reconstruct'

    return StokesletReconstruct( self.mycoords, self.topo, self.coords, velo, trac, surftens, self.mu )

class StokesletTrac( function.ArrayFunc ):
  'stokeslet stress'

  def __init__( self, mycoords, topo, coords, funcsp ):
    'constructor'

    self.shape = int( funcsp.shape[0] ) * 2, 2

    assert isinstance( funcsp, function.Evaluable )
    iweights = coords.iweights( topo.ndims )
    cache = CacheFuncND( topo, coords, iweights, funcsp, coords.normal(), funcsp.shape[0] )
    self.args = int(funcsp.shape[0]), funcsp.shape[0], funcsp, mycoords.normal(), Convolution( mycoords, cache )

  @staticmethod
  def eval( ndofs, N, func, norm, (shape,iterdata) ):
    'evaluate'

    retval = numpy.zeros( (ndofs*2,2) + shape )
    retval_swap = retval.reshape( 2, ndofs, 2, *shape ).swapaxes(0,1) # follows ordering Vectorize
    tot = numpy.zeros( (len(N),2,2) + shape )
    for D_R, R, R2, logR, w, f, n, M in iterdata:
      kernel = D_R[:,_,_] * D_R[_,:,_] * D_R[_,_,:]
      kernel /= R
      trac    =    f[:,_,_,_,:] * util.contract( kernel, norm[:,:,_], -3 )
      kernsub = func[:,_,_,:,_] * util.contract( kernel,    n[:,_,:], -3 )
      if N is M:
        retval_swap[M] += numpy.dot( trac + kernsub, w )
      else:
        retval_swap[M] += numpy.dot( trac, w )
        retval_swap[N] += numpy.dot( kernsub, w )
    retval /= -numpy.pi
    return retval

class StokesletGrad( function.ArrayFunc ):
  'stokeslet gradient'

  def __init__( self, mycoords, topo, coords, funcsp, mu ):
    'constructor'

    self.shape = int( funcsp.shape[0] ) * 2, 2, 2

    iweights = coords.iweights( topo.ndims ) * funcsp
    self.args = int(funcsp.shape[0]), Convolution( mycoords, topo, coords, iweights, funcsp.shape[0] ), mu

  @staticmethod
  def eval( ndofs, (shape,iterdata), mu ):
    'evaluate'

    retval = numpy.zeros( (ndofs*2,2,2) + shape )
    retval_swap = retval.reshape( 2, ndofs, 2, 2, *shape ).swapaxes(0,1) # follows ordering Vectorize
    for D_R, R, R2, logR, w, I in iterdata:
      kernel = D_R[:,_,_] * D_R[_,:,_] * D_R[_,_,:]
      kernel *= -2
      util.ndiag( kernel, [1,2] )[:] += D_R
      util.ndiag( kernel, [0,1] )[:] -= D_R
      util.ndiag( kernel, [2,0] )[:] += D_R
      kernel /= R
      retval_swap[I] += numpy.tensordot( w, kernel, (-1,-1) )
    retval /= 4 * numpy.pi * mu
    return retval

class StokesletReconstruct( function.ArrayFunc ):
  'stokeslet reconstruction'

  def __init__( self, mycoords, topo, coords, velo, trac, surftens, mu ):
    'constructor'

    self.shape = 2,
    assert velo is 0 or velo.shape == (2,)
    assert trac is 0 or trac.shape == (2,)

    iweights = coords.iweights( topo.ndims )
    cache = CacheFuncND( topo, coords, iweights, velo, trac, coords.normal() )
    self.args = Convolution( mycoords, cache ), surftens, mu

  @staticmethod
  def eval( (shape,iterdata), surftens, mu ):
    'evaluate'

    retval = 0
    for D_R, R, R2, logR, w, velo, trac, norm in iterdata:
      Dnorm = util.contract( D_R, norm[:,_,:], axis=0 )
      d = 0
      if velo is not 0:
        Dvelo = util.contract( D_R, velo[:,_,:], axis=0 )
        d -= (D_R/R) * Dnorm * Dvelo
      if trac is not 0:
        Dtrac = util.contract( D_R, trac[:,_,:], axis=0 )
        d += ( D_R * Dtrac - trac[:,_,:] * logR ) / (4*mu)
      if surftens is not 0:
        d += (D_R/R) * ( 1 - Dnorm**2 ) * (surftens/(2*mu))
      retval += numpy.dot( d, w )
    return retval / numpy.pi

# STOKES 3D

class Stokeslet3D( function.ArrayFunc ):
  'stokeslet'

  def __init__( self, mycoords, topo, coords, funcsp, mu ):
    'constructor'

    self.topo = topo
    self.coords = coords
    self.mycoords = mycoords
    self.funcsp = funcsp
    self.mu = mu
    self.shape = int( funcsp.shape[0] ) * 3, 3
    self.args = int(funcsp.shape[0]), funcsp.shape[0], funcsp, coords.normal(), Convolution3D( mycoords, topo, coords, funcsp, funcsp.shape[0] ), mu

  @staticmethod
  def eval( ndofs, N, func, norm, (shape,iterdata), mu ):
    'evaluate'

    retval = numpy.zeros( (ndofs*3,3) + shape )
    retval_swap = retval.reshape( 3, ndofs, 3, *shape ).swapaxes(0,1) # follows ordering Vectorize
    for D_R, R, R2, f, w, n, M in iterdata:
      kernel = D_R[:,_] * D_R[_,:]
      kernel[0,0] += 1
      kernel[1,1] += 1
      kernel[2,2] += 1
      kernel /= R
      retval_swap[M] += numpy.tensordot( w * f, kernel, (-1,-1) )

#     flux = util.contract( kernel, n[:,_,:], axis=-3 )
#     shift = util.contract( norm[:,:,_], flux, axis=0 )

#     if M is N:
#       #tmp = f[:,_,_,_,:] * kernel
#       #tmp[:,0,0] -= shift
#       #tmp[:,1,1] -= shift
#       #tmp[:,2,2] -= shift
#       if n[0,0] > .3:
#         from matplotlib import pylab

#         tmp = ( kernel - flux / norm[0,0] ).reshape( 3, 3, 13, 111, 111 )
#         pylab.imshow( tmp[0,0,0] )
#         pylab.colorbar()
#         pylab.show()

#     retval_swap[N,0,0] -= numpy.dot( shift, w )
#     retval_swap[N,1,1] -= numpy.dot( shift, w )
#     retval_swap[N,2,2] -= numpy.dot( shift, w )

    retval /= 8 * numpy.pi * mu
    return retval

  def traction( self ):
    'flux'

    return StokesletTrac3D( self.mycoords, self.topo, self.coords, self.funcsp )

  def curvature( self ):
    'curvature'

    grad = StokesletGrad3D( self.mycoords, self.topo, self.coords, self.funcsp, self.mu )
    normal = self.coords.normal()
    return grad.trace(1,2) - ( ( grad * normal ).sum() * normal ).sum()

  def reconstruct( self, bval, flux ):
    'reconstruct'

    return StokesletReconstruct( self.mycoords, self.topo, self.coords, bval, flux, self.mu )

class StokesletGrad3D( function.ArrayFunc ):
  'stokeslet gradient'

  def __init__( self, mycoords, topo, coords, funcsp, mu ):
    'constructor'

    self.shape = int( funcsp.shape[0] ) * 3, 3, 3
    self.args = int(funcsp.shape[0]), Convolution3D( mycoords, topo, coords, funcsp, funcsp.shape[0] ), mu

  @staticmethod
  def eval( ndofs, (shape,iterdata), mu ):
    'evaluate'

    retval = numpy.zeros( (ndofs*3,3,3) + shape )
    retval_swap = retval.reshape( 3, ndofs, 3, 3, *shape ).swapaxes(0,1) # follows ordering Vectorize
    for D_R, R, R2, f, w, n, I in iterdata:
      kernel = D_R[:,_,_] * D_R[_,:,_] * D_R[_,_,:]
      kernel *= -3
      kernel[:,0,0] += D_R
      kernel[:,1,1] += D_R
      kernel[:,2,2] += D_R
      kernel[0,:,0] += D_R
      kernel[1,:,1] += D_R
      kernel[2,:,2] += D_R
      kernel[0,0,:] -= D_R
      kernel[1,1,:] -= D_R
      kernel[2,2,:] -= D_R
      kernel /= R2
      retval_swap[I] += numpy.tensordot( w * f, kernel, (-1,-1) )
    retval /= 8 * numpy.pi * mu
    return retval

class StokesletPres3D( function.ArrayFunc ):
  'stokeslet gradient'

  def __init__( self, mycoords, topo, coords, funcsp ):
    'constructor'

    self.shape = int( funcsp.shape[0] ) * 3
    self.args = int(funcsp.shape[0]), Convolution3D( mycoords, topo, coords, funcsp, funcsp.shape[0] )

  @staticmethod
  def eval( ndofs, (shape,iterdata) ):
    'evaluate'

    retval = numpy.zeros( (ndofs*3,) + shape )
    retval_swap = retval.reshape( 3, ndofs, *shape ).swapaxes(0,1) # follows ordering Vectorize
    for D_R, R, R2, f, w, n, I in iterdata:
      retval_swap[I] += numpy.tensordot( w * f, D_R / R2, (-1,-1) )
    retval /= 4 * numpy.pi
    return retval

class StokesletTrac3D( function.ArrayFunc ):
  'stokeslet stress'

  def __init__( self, mycoords, topo, coords, funcsp ):
    'constructor'

    self.shape = int( funcsp.shape[0] ) * 3, 3
    self.args = int(funcsp.shape[0]), funcsp.shape[0], funcsp, mycoords.normal(), Convolution3D( mycoords, topo, coords, funcsp, funcsp.shape[0] )

  @staticmethod
  def eval( ndofs, N, func, norm, (shape,iterdata) ):
    'evaluate'

    retval = numpy.zeros( (ndofs*3,3) + shape )
    retval_swap = retval.reshape( 3, ndofs, 3, *shape ).swapaxes(0,1) # follows ordering Vectorize
    for D_R, R, R2, f, w, n, M in iterdata:
      kernel = D_R[:,_,_] * D_R[_,:,_] * D_R[_,_,:]
      kernel /= R2
      trac    =    f[:,_,_,_,:] * util.contract( kernel, norm[:,:,_], -3 )
      kernsub = func[:,_,_,:,_] * util.contract( kernel,    n[:,_,:], -3 )
      if N is M:
        retval_swap[M] += numpy.dot( trac + kernsub, w )
      else:
        retval_swap[M] += numpy.dot( trac, w )
        retval_swap[N] += numpy.dot( kernsub, w )
    retval *= -.75 / numpy.pi
    ## instead of kernsub:
    #retval_swap[N,0,0] += .5 * func
    #retval_swap[N,1,1] += .5 * func
    #retval_swap[N,2,2] += .5 * func
    return retval

# tests

def testgrad( domain, coords, funcsp ):
  'numeric differentiation test'

  import mesh
  testdomain, testcoords = mesh.rectilinear( [0,1], [0,1] )
  elem, = testdomain
  eps = 1e-10
  p = elem.eval( numpy.array([[.5],[.5]]) )
  px = elem.eval( numpy.array([[.5-.5*eps,.5+.5*eps],[.5,.5]]) )
  py = elem.eval( numpy.array([[.5,.5],[.5-.5*eps,.5+.5*eps]]) )

  fval = Stokeslet( testcoords, domain, coords, funcsp, 1. )
  grad = StokesletGrad( testcoords, domain, coords, funcsp, 1. )

  dx = fval(px)
  print numpy.hstack([ ( dx[...,1] - dx[...,0] ) / eps, grad( p )[...,0,0] ])
  print ( dx[...,1] - dx[...,0] ) / eps - grad( p )[...,0,0]
  print

  dy = fval(py)
  print numpy.hstack([ ( dy[...,1] - dy[...,0] ) / eps, grad( p )[...,1,0] ])
  print ( dy[...,1] - dy[...,0] ) / eps - grad( p )[...,1,0]

  raise SystemExit

def testgrad3D( domain, coords, funcsp ):
  'numeric differentiation test'

  import mesh
  testdomain, testcoords = mesh.rectilinear( [0,1], [0,1], [0,1] )
  elem, = testdomain
  eps = 1e-10
  MU = 1.

  fval = Stokeslet3D( testcoords, domain, coords, funcsp, MU )
  grad = StokesletGrad3D( testcoords, domain, coords, funcsp, MU )

  p = elem.eval( numpy.array([[.5],[.5],[.5]]) )
  gradp = grad( p )[...,0]
  diffp = gradp.copy()

  print gradp[:,0,0] + gradp[:,1,1] + gradp[:,2,2]

  px = elem.eval( numpy.array([[.5-.5*eps,.5+.5*eps],[.5,.5],[.5,.5]]) )
  dx = fval(px)
  #print numpy.hstack([ ( dx[...,1] - dx[...,0] ) / eps, gradp[...,0] ])

  py = elem.eval( numpy.array([[.5,.5],[.5-.5*eps,.5+.5*eps],[.5,.5]]) )
  dy = fval(py)
  #print numpy.hstack([ ( dy[...,1] - dy[...,0] ) / eps, gradp[...,1] ])

  pz = elem.eval( numpy.array([[.5,.5],[.5,.5],[.5-.5*eps,.5+.5*eps]]) )
  dz = fval(pz)
  #print numpy.hstack([ ( dz[...,1] - dz[...,0] ) / eps, gradp[...,2] ])

  diffp[...,0] -= ( dx[...,1] - dx[...,0] ) / eps
  diffp[...,1] -= ( dy[...,1] - dy[...,0] ) / eps
  diffp[...,2] -= ( dz[...,1] - dz[...,0] ) / eps

  print diffp

  pres = StokesletPres3D( testcoords, domain, coords, funcsp )
  presp = pres( p )[...,0]

  stress = MU * ( gradp + gradp.swapaxes(-2,-1) )
  stress[...,0,0] -= presp
  stress[...,1,1] -= presp
  stress[...,2,2] -= presp
  print stress[0]

  stress = StokesletStress3D( testcoords, domain, coords, funcsp )( p )[...,0]
  print stress[0]

  raise SystemExit

#---- COLOCATION ----

class Stokeslet3DColoc( function.ArrayFunc ):
  'stokeslet'

  def __init__( self, coords, points, mu ):
    'constructor'

    self.coords = coords
    self.points = points
    self.mu = mu
    npoints, ndims = points.shape
    self.shape = npoints * ndims, ndims
    self.args = coords, points, mu

  @staticmethod
  def eval( x, points, mu ):
    'evaluate'

    npoints, ndims = points.shape

    D = x[_,:,:] - points[:,:,_]
    R2 = util.contract( D, D, 1 )
    R = numpy.sqrt( R2 )

    U = D[:,:,_,:] * D[:,_,:,:]
    U *= 1. / R2[:,_,_,:]
    U[:,0,0,:] += 1
    U[:,1,1,:] += 1
    U[:,2,2,:] += 1
    U *= ( .125 / ( numpy.pi * mu ) ) / R[:,_,_,:]

    return U.swapaxes(0,1).reshape( npoints * ndims, ndims, -1 )

  def grad( self ):
    'gradient'

    return StokesletGrad3DColoc( self.coords, self.points, self.mu )

  def pres( self ):
    'pressure'

    return StokesletPres3DColoc( self.coords, self.points )

  def stress( self ):
    'gradient'

    return StokesletStress3DColoc( self.coords, self.points )

  def traction( self ):
    'flux'

    return ( self.stress() * self.coords.normal() ).sum()

class StokesletPres3DColoc( function.ArrayFunc ):
  'stokeslet'

  def __init__( self, coords, points ):
    'constructor'

    npoints, ndims = points.shape
    self.shape = npoints * ndims
    self.args = coords, points

  @staticmethod
  def eval( x, points ):
    'evaluate'

    npoints, ndims = points.shape

    D = x[_,:,:] - points[:,:,_]
    R2 = util.contract( D, D, 1 )
    R = numpy.sqrt( R2 )
    R3 = R * R2

    U = D / ( ( 4 * numpy.pi ) * R3[:,_,:] )

    return U.swapaxes(0,1).reshape( npoints * ndims, -1 )

class StokesletGrad3DColoc( function.ArrayFunc ):
  'stokeslet'

  def __init__( self, coords, points, mu ):
    'constructor'

    npoints, ndims = points.shape
    self.shape = npoints * ndims, ndims, ndims
    self.args = coords, points, mu

  @staticmethod
  def eval( x, points, mu ):
    'evaluate'

    npoints, ndims = points.shape

    D = x[_,:,:] - points[:,:,_]
    R2 = util.contract( D, D, 1 )
    R = numpy.sqrt( R2 )
    R3 = R * R2

    U = D[:,:,_,_,:] * D[:,_,:,_,:] * D[:,_,_,:,:]
    U *= -3. / R2[:,_,_,_,:]
    U[:,0,0,:,:] -= D
    U[:,1,1,:,:] -= D
    U[:,2,2,:,:] -= D
    U[:,0,:,0,:] += D
    U[:,1,:,1,:] += D
    U[:,2,:,2,:] += D
    U[:,:,0,0,:] += D
    U[:,:,1,1,:] += D
    U[:,:,2,2,:] += D
    U *= ( .125 / ( numpy.pi * mu ) ) / R3[:,_,_,_,:]

    return U.swapaxes(0,1).reshape( npoints * ndims, ndims, ndims, -1 )

class StokesletStress3DColoc( function.ArrayFunc ):
  'stokeslet'

  def __init__( self, coords, points ):
    'constructor'

    npoints, ndims = points.shape
    self.shape = npoints * ndims, ndims, ndims
    self.args = coords, points

  @staticmethod
  def eval( x, points ):
    'evaluate'

    npoints, ndims = points.shape

    D = x[_,:,:] - points[:,:,_]
    R2 = util.contract( D, D, 1 )
    R = numpy.sqrt( R2 )
    R4 = R2**2
    R5 = R * R4

    U = D[:,:,_,_,:] * D[:,_,:,_,:] * D[:,_,_,:,:]
    U *= ( -.75 / numpy.pi ) / R5[:,_,_,_,:]

    return U.swapaxes(0,1).reshape( npoints * ndims, ndims, ndims, -1 )

def testgrad3Dcoloc( points ):
  'numeric differentiation test'

  import mesh
  domain, coords = mesh.rectilinear( [0,1], [0,1], [0,1] )
  elem, = domain
  eps = 1e-10
  MU = 1.

  fval = Stokeslet3DColoc( coords, points, MU )
  grad = fval.grad()

  p = elem.eval( numpy.array([[.5],[.5],[.5]]) )
  gradp = grad( p )[...,0]
  diffp = gradp.copy()

  print 'DIVERGENCE'
  print gradp[:,0,0] + gradp[:,1,1] + gradp[:,2,2]
  print

  px = elem.eval( numpy.array([[.5-.5*eps,.5+.5*eps],[.5,.5],[.5,.5]]) )
  dx = fval(px)
  #print numpy.hstack([ ( dx[...,1] - dx[...,0] ) / eps, gradp[...,0] ])

  py = elem.eval( numpy.array([[.5,.5],[.5-.5*eps,.5+.5*eps],[.5,.5]]) )
  dy = fval(py)
  #print numpy.hstack([ ( dy[...,1] - dy[...,0] ) / eps, gradp[...,1] ])

  pz = elem.eval( numpy.array([[.5,.5],[.5,.5],[.5-.5*eps,.5+.5*eps]]) )
  dz = fval(pz)
  #print numpy.hstack([ ( dz[...,1] - dz[...,0] ) / eps, gradp[...,2] ])

  diffp[...,0] -= ( dx[...,1] - dx[...,0] ) / eps
  diffp[...,1] -= ( dy[...,1] - dy[...,0] ) / eps
  diffp[...,2] -= ( dz[...,1] - dz[...,0] ) / eps

  print 'FINITE DIFFERENCE ERROR'
  print util.norm2( diffp.reshape( -1, 9 ), axis=1 )
  print

  pres = fval.pres()
  presp = pres( p )[...,0]

  stress = MU * ( gradp + gradp.swapaxes(-2,-1) )
  stress[...,0,0] -= presp
  stress[...,1,1] -= presp
  stress[...,2,2] -= presp

  print 'STRESS ERROR'
  print util.norm2( ( stress - fval.stress()( p )[...,0] ).reshape( -1, 9 ), axis=1 )
  print

  raise SystemExit

# vim:shiftwidth=2:foldmethod=indent:foldnestmax=2
