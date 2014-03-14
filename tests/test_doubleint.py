#!/usr/bin/env python

from nutils import *
import numpy, warnings
almostEquals = lambda val, places=7: numeric.abs( val ) < 10.**(-places)
infnorm = lambda f: numpy.linalg.norm( f, numeric.inf )
grid = lambda n: numeric.linspace( -n/2., n/2., n+1 )
pi = numeric.pi
def V( x, y ):
  rInv = function.norm2( x-y )**-1.
  return 0.125*pi**-1. * (function.eye(3)*rInv + (x-y)[:,_]*(x-y)[_,:]*rInv**3)
def K( x, y ):
  rInv = function.norm2( x-y )**-1.
  return 0.75*pi**-1. * (x-y)[:,_]*(x-y)[_,:] * ((x-y)*y.normal()).sum() * rInv**5
def inputfile( case='sphere' ):
  f = open( '%s.sdv'%case, 'w' )
  if case=='sphere':
    f.write( 'Verts 26\n0. 0. -1.\n0. 0. 1.\n1. 0. 0.\n0. -1. 0.\n-1. 0. 0.\n0. 1. 0.\n0.75 0. -0.75\n0. 0.75 -0.75\n0.75 0.75 0.\n0. -0.75 -0.75\n0.75 -0.75 0.\n-0.75 0. -0.75\n-0.75 -0.75 0.\n-0.75 0.75 0.\n0.75 0. 0.75\n0. 0.75 0.75\n0. -0.75 0.75\n-0.75 0. 0.75\n0.555556 0.555556 -0.555556\n0.555556 -0.555556 -0.555556\n-0.555556 -0.555556 -0.555556\n-0.555556 0.555556 -0.555556\n0.555556 0.555556 0.555556\n0.555556 -0.555556 0.555556\n-0.555556 -0.555556 0.555556\n-0.555556 0.555556 0.555556\nRings 24 0\n4 3 18 8 5 7 0 6 2 20 9 19 10 11 21 13\n4 3 19 10 2 6 0 9 3 21 11 20 12 7 18 8\n4 3 20 12 3 9 0 11 4 18 7 21 13 6 19 10\n4 3 21 13 4 11 0 7 5 19 6 18 8 9 20 12\n4 3 22 8 2 14 1 15 5 24 17 25 13 16 23 10\n4 3 25 13 5 15 1 17 4 23 16 24 12 14 22 8\n4 3 24 12 4 17 1 16 3 22 14 23 10 15 25 13\n4 3 23 10 3 16 1 14 2 25 15 22 8 17 24 12\n4 3 18 7 0 6 2 8 5 23 14 22 15 10 19 9\n4 3 22 15 5 8 2 14 1 19 10 23 16 6 18 7\n4 3 23 16 1 14 2 10 3 18 6 19 9 8 22 15\n4 3 19 9 3 10 2 6 0 22 8 18 7 14 23 16\n4 3 19 6 0 9 3 10 2 24 16 23 14 12 20 11\n4 3 23 14 2 10 3 16 1 20 12 24 17 9 19 6\n4 3 24 17 1 16 3 12 4 19 9 20 11 10 23 14\n4 3 20 11 4 12 3 9 0 23 10 19 6 16 24 17\n4 3 20 9 0 11 4 12 3 25 17 24 16 13 21 7\n4 3 24 16 3 12 4 17 1 21 13 25 15 11 20 9\n4 3 25 15 1 17 4 13 5 20 11 21 7 12 24 16\n4 3 21 7 5 13 4 11 0 24 12 20 9 17 25 15\n4 3 22 14 1 15 5 8 2 21 7 18 6 13 25 17\n4 3 18 6 2 8 5 7 0 25 13 21 11 15 22 14\n4 3 21 11 0 7 5 13 4 22 15 25 17 8 18 6\n4 3 25 17 4 13 5 15 1 18 8 22 14 7 21 11\nFaceGroups 1\nGroupI 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23' )
  if case=='tripod':
    f.write( 'Verts 74\n 0. 0. -1.\n 1. 0. 2.\n 2. 0. 1.\n 0. -1. 0.\n -1. 0. 0.\n 1. 2. 0.\n 3. 0. 0.\n 2. 1. 0.\n 2. 0. -1.\n 2. -1. 0.\n 0. 0. 3.\n 0. 1. 2.\n 0. -1. 2.\n -1. 0. 2.\n 0. 3. 0.\n -1. 2. 0.\n 0. 2. -1.\n 0. 2. 1.\n 1. 0. -1.\n 0. 1. -1.\n 1.25 1.25 0.\n 0. -0.75 -0.75\n 1. -1. 0.\n -0.75 0. -0.75\n -0.75 -0.75 0.\n -1. 1. 0.\n 1.25 0. 1.25\n 0. 1.25 1.25\n 0. -1. 1.\n -1. 0. 1.\n 2.75 0. -0.75\n 2.75 0.75 0.\n 2.75 -0.75 0.\n 2.75 0. 0.75\n 2. -0.75 0.75\n 2. 0.75 0.75\n 2. 0.75 -0.75\n 2. -0.75 -0.75\n 0.75 0. 2.75\n 0. 0.75 2.75\n 0. -0.75 2.75\n -0.75 0. 2.75\n 0.75 0.75 2.\n 0.75 -0.75 2.\n -0.75 0.75 2.\n -0.75 -0.75 2.\n 0. 2.75 -0.75\n 0.75 2.75 0.\n -0.75 2.75 0.\n 0. 2.75 0.75\n 0.75 2. -0.75\n 0.75 2. 0.75\n -0.75 2. 0.75\n -0.75 2. -0.75\n 1. 1. -0.84\n 1. -0.75 -0.75\n -0.5555556 -0.5555556 -0.5555556\n -0.75 1. -0.75\n 1. 1. 1.\n 1. -0.84 1.\n -0.75 -0.75 1.\n -0.84 1. 1.\n 2.555556 0.555556 -0.5555556\n 2.555556 -0.555556 -0.5555556\n 2.555556 0.555556 0.5555556\n 2.555556 -0.555556 0.5555556\n 0.555556 0.555556 2.5555556\n 0.555556 -0.555556 2.5555556\n -0.555556 -0.555556 2.5555556\n -0.555556 0.555556 2.5555556\n 0.555556 2.555556 -0.5555556\n -0.555556 2.555556 -0.5555556\n 0.555556 2.555556 0.5555556\n -0.555556 2.555556 0.5555556\n' )
    f.write( 'Rings 72 0\n 4 5 54 50 16 19 0 18 8 36 7 20 5 56 21 55 37 23 57 53\n 5 4 0 23 56 21 55 18 54 19 57 9 37 8 36 22 3 24\n 4 3 56 24 3 21 0 23 4 54 19 57 25 18 55 22\n 5 4 0 18 54 19 57 23 56 21 55 15 25 4 24 53 16 50\n 4 5 59 28 12 43 1 26 2 34 9 22 3 66 42 58 35 38 67 40\n 4 6 58 35 2 26 1 42 11 27 17 51 5 20 7 67 38 66 39 43 59 34\n 4 3 66 39 11 42 1 38 10 59 43 67 40 26 58 27\n 4 3 67 40 10 38 1 43 12 58 26 59 28 42 66 39\n 4 6 58 20 7 35 2 26 1 42 11 27 17 51 5 65 34 59 43 33 64 31\n 4 5 59 43 1 26 2 34 9 22 3 28 12 64 33 65 32 35 58 42\n 4 3 65 32 9 34 2 33 6 58 35 64 31 26 59 22\n 4 3 64 31 6 33 2 35 7 59 26 58 20 34 65 32\n 5 4 3 28 59 22 55 21 56 24 60 8 18 0 23 37 9 34\n 4 5 59 34 9 22 3 28 12 43 1 26 2 56 24 60 45 21 55 37\n 5 4 3 21 56 24 60 28 59 22 55 13 45 12 43 29 4 23\n 4 3 56 23 4 24 3 21 0 59 22 55 18 28 60 29\n 4 3 56 21 0 23 4 24 3 61 29 60 28 25 57 19\n 5 4 4 25 61 29 60 24 56 23 57 12 28 3 21 45 13 44\n 4 5 61 44 13 29 4 25 15 52 17 27 11 56 23 57 53 24 60 45\n 5 4 4 24 56 23 57 25 61 29 60 16 53 15 52 19 0 21\n 4 6 58 27 17 51 5 20 7 35 2 26 1 42 11 70 50 54 36 47 72 49\n 4 5 54 36 7 20 5 50 16 19 0 18 8 72 47 70 46 51 58 35\n 4 3 70 46 16 50 5 47 14 58 51 72 49 20 54 19\n 4 3 72 49 14 47 5 51 17 54 20 58 27 50 70 46\n 4 3 62 36 8 30 6 31 7 65 33 64 35 32 63 37\n 4 3 64 35 7 31 6 33 2 63 32 65 34 30 62 36\n 4 3 65 34 2 33 6 32 9 62 30 63 37 31 64 35\n 4 3 63 37 9 32 6 30 8 64 31 62 36 33 65 34\n 4 5 54 18 8 36 7 20 5 50 16 19 0 64 35 58 51 31 62 30\n 4 6 58 51 5 20 7 35 2 26 1 42 11 27 17 62 31 64 33 36 54 50\n 4 3 64 33 2 35 7 31 6 54 36 62 30 20 58 26\n 4 3 62 30 6 31 7 36 8 58 20 54 18 35 64 33\n 5 4 8 36 54 18 55 37 63 30 62 3 22 9 32 21 0 19\n 4 5 54 19 0 18 8 36 7 20 5 50 16 63 30 62 31 37 55 21\n 4 3 62 31 7 36 8 30 6 55 37 63 32 18 54 20\n 4 3 63 32 6 30 8 37 9 54 18 55 22 36 62 31\n 4 5 59 26 2 34 9 22 3 28 12 43 1 63 37 55 21 32 65 33\n 5 4 9 32 63 37 55 22 59 34 65 0 21 3 28 18 8 30\n 4 3 63 30 8 37 9 32 6 59 34 65 33 22 55 18\n 4 3 65 33 6 32 9 34 2 55 22 59 26 37 63 30\n 4 3 66 42 1 38 10 39 11 68 41 69 44 40 67 43\n 4 3 69 44 11 39 10 41 13 67 40 68 45 38 66 42\n 4 3 68 45 13 41 10 40 12 66 38 67 43 39 69 44\n 4 3 67 43 12 40 10 38 1 69 39 66 42 41 68 45\n 4 6 58 26 1 42 11 27 17 51 5 20 7 35 2 69 44 61 52 39 66 38\n 4 5 61 52 17 27 11 44 13 29 4 25 15 66 39 69 41 42 58 51\n 4 3 69 41 13 44 11 39 10 58 42 66 38 27 61 29\n 4 3 66 38 10 39 11 42 1 61 27 58 26 44 69 41\n 5 4 12 43 59 28 60 45 68 40 67 4 29 13 41 24 3 22\n 4 5 59 22 3 28 12 43 1 26 2 34 9 68 40 67 38 45 60 24\n 4 3 67 38 1 43 12 40 10 60 45 68 41 28 59 26\n 4 3 68 41 10 40 12 45 13 59 28 60 29 43 67 38\n 4 5 61 27 11 44 13 29 4 25 15 52 17 68 45 60 24 41 69 39\n 5 4 13 41 68 45 60 29 61 44 69 3 24 4 25 28 12 40\n 4 3 68 40 12 45 13 41 10 61 44 69 39 29 60 28\n 4 3 69 39 10 41 13 44 11 60 29 61 27 45 68 40\n 4 3 72 51 17 49 14 47 5 71 46 70 50 48 73 52\n 4 3 70 50 5 47 14 46 16 73 48 71 53 49 72 51\n' )
    f.write( '4 3 71 53 16 46 14 48 15 72 49 73 52 47 70 50\n 4 3 73 52 15 48 14 49 17 70 47 72 51 46 71 53\n 5 4 15 52 61 25 57 53 71 48 73 0 19 16 46 23 4 29\n 4 5 61 29 4 25 15 52 17 27 11 44 13 71 48 73 49 53 57 23\n 4 3 73 49 17 52 15 48 14 57 53 71 46 25 61 27\n 4 3 71 46 14 48 15 53 16 61 25 57 19 52 73 49\n 4 5 54 20 5 50 16 19 0 18 8 36 7 71 53 57 23 46 70 47\n 5 4 16 46 71 53 57 19 54 50 70 4 23 0 18 25 15 48\n 4 3 71 48 15 53 16 46 14 54 50 70 47 19 57 25\n 4 3 70 47 14 46 16 50 5 57 19 54 20 53 71 48\n 4 5 61 25 15 52 17 27 11 44 13 29 4 72 51 58 42 49 73 48\n 4 6 58 42 11 27 17 51 5 20 7 35 2 26 1 73 49 72 47 52 61 44\n 4 3 72 47 5 51 17 49 14 61 52 73 48 27 58 20\n 4 3 73 48 14 49 17 52 15 58 27 61 25 51 72 47\n' )
    f.write( 'FaceGroups 1\n GroupI 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71')
  f.close()

class TestGaussDoubleInt( object ):
  'Gauss quadrature on product domain.'
  def distance( self, val ):
    assert almostEquals( val - 0. )

  def distancesquared( self, val ):
    assert almostEquals( val - 1./6 )

  def test_polynomials( self ):
    domain, geom = mesh.rectilinear( [[0,.5,1]] )
    ddomain = domain * domain

    iw = function.iwscale( geom, domain.ndims )
    iweights = iw * function.opposite( iw ) * function.IWeights()

    x = geom[0]
    y = function.opposite( geom[0] )

    self.distance( ddomain.integrate( x-y, iweights=iweights, ischeme='gauss2' ) )
    self.distancesquared( ddomain.integrate( (x-y)**2, iweights=iweights, ischeme='gauss2' ) )

class TestSingularDoubleInt( object ):
  'Regularized quadrature on product domain.'
  def patch( self, val ):
    assert almostEquals( val - 1. )

  def distance( self, val ):
    assert almostEquals( val - 1./3 )

  def test_Integration( self ):
    grid = numeric.linspace( 0., 1., 4 )
    domain, geom = mesh.rectilinear( 2*(grid,) )
    ddomain = domain * domain

    iw = function.iwscale( geom, domain.ndims )
    iweights = iw * function.opposite( iw ) * function.IWeights()

    x = geom
    y = function.opposite( geom[0] )
    r = function.norm2( x-y )
    
    self.patch( ddomain.integrate( 1, iweights=iweights, ischeme='singular2' ) )
    self.distance( ddomain.integrate( r**2, iweights=iweights, ischeme='singular3' ) )

class TestNormalInKernelOfV( object ):
  'Convolute normal with single-layer to verify it is in the kernel, note that this works with all gauss schemes!'
  def __init__( self ):
    'Geometry definitions.'
    domain, geom = mesh.rectilinear( (grid(4),grid(2)), periodic=(0,) )
    self.geom, self.domain, self.ddomain = geom, domain, domain * domain

  def template( self, degree, geometry, dump=False ):
    'Template for Vn = 0 tests on different geometries.'
    trac = self.domain.splinefunc( degree=2*(2,) ).vector(3)
    if dump:
      geo = domain.projection( geometry, onto=trac, geometry=self.geom )
      refine = 3 if geometry is sphere else 0
      plot.writevtu( 'geometry.vtu', domain.refine(3), geometry )
      if refine: warnings.warn( 'The plotted geometry is a projection, and a bad approximation.' )

    iw = function.iwscale( geometry, self.domain.ndims )
    iweights = iw * function.opposite( iw ) * function.IWeights()

    x = geometry
    y = function.opposite( geometry )

    return self.ddomain.integrate( (V(x,y)*x.normal()).sum(), iweights=iweights, ischeme='singular{0}'.format(degree) )

  def test_SphericalGeometry( self ):
    'n in ker(V), case: sphere.'
    cos, sin, pi = function.cos, function.sin, numeric.pi
    phi, theta = .5*pi*self.geom # |phi| < pi, |theta| < pi/2
    self.sphere = function.stack( [cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)] )

    err = self.template( 4, self.sphere )
    assert almostEquals( infnorm(err) - 0., places=3 )

  def test_PolyhedronGeometry( self ):
    'n in ker(V), case: polyhedron'
    # raise NotImplemented( 'Choose has no .opposite yet' )

    abs = function.abs
    xi, eta = self.geom
    self.octahedron = function.Concatenate( [[(1.-abs( eta ))*function.piecewise( xi, (-1., 0., 1.), 1, -2*xi-1, -1, 2*xi-3 )], [(1.-abs( eta ))*function.piecewise( xi, (-1., 0., 1.), 2*xi+3, 1, 1-2*xi, -1 )], [numeric.sqrt(2)*eta]] )

    err = self.template( 4, self.octahedron )
    assert almostEquals( infnorm( err ) - 0., places=3 )

class TestKroneckerKernelGivesSurface( object ):
  'Convolute a uniform velocity field with the identity to verify it gives the surface.'

  def test_SphericalGeometry( self ):
    'Integrating Id gives surf, case: sphere.'
    domain, geom = mesh.rectilinear( (grid(4),grid(2)), periodic=(0,) )
    cos, sin, pi = function.cos, function.sin, numeric.pi
    phi, theta = .5*pi*geom # |phi| < pi, |theta| < pi/2
    sphere = function.stack( [cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)] )
    velo = domain.splinefunc( degree=2*(2,) ).vector(3)
    vinf = function.stack( (0.,0.,1.) )
    val = domain.integrate( (velo*vinf).sum(-1), geometry=sphere, ischeme='gauss8' ).sum()

    surf = 4.*pi
    assert almostEquals( val-surf )

class TestOneInKernelOfK( object ):
  'Convolute a uniform velocity field with the dual-layer to verify it is in the kernel.'
  def __init__( self ):
    'Geometry definitions.'
    domain, geom = mesh.rectilinear( (grid(4),grid(2)), periodic=(0,) )
    self.geom, self.domain, self.ddomain = geom, domain, domain * domain

  def template( self, geometry, degree ):
    trac = self.domain.splinefunc( degree=2*(2,) ).vector(3)
    vinf = function.stack( (0.,0.,1.) )

    iw = function.iwscale( geometry, self.domain.ndims )
    iweights = iw * function.opposite( iw ) * function.IWeights()

    x = geometry
    y = function.opposite( geometry )

    Kvinf = (K(x, y)[:,:]*vinf[_,:]).sum(-1)
    doublelayer = self.ddomain.integrate( (trac[:,:]*Kvinf[_,:]).sum(), iweights=iweights, ischeme='singular{0}'.format(degree) )
    identity = self.domain.integrate( (trac*vinf).sum(), geometry=geometry, ischeme='gauss4' )

    return .5*identity + doublelayer

  def test_SphericalGeometry( self ):
    '1 in ker(K), case: sphere.'
    cos, sin, pi = function.cos, function.sin, numeric.pi
    phi, theta = .5*pi*self.geom # |phi| < pi, |theta| < pi/2
    self.sphere = function.stack( [cos(phi)*cos(theta), sin(phi)*cos(theta), sin(theta)] )

    err = infnorm( self.template( self.sphere, 4 ) )
    assert almostEquals( err - 0., places=2 )

class TestShearFlow( object ):
  'Torus cut of shear flow.'
  def _bemtest( self, domain, coords, ischeme, tol=1e-2, visual=False ):
    'Interior Stokes-Dirichlet BEM test case: shear flow.'
    # domain, space
    funcsp = domain.splinefunc( degree=3 ).vector(3)
    l2norm = lambda func, title='l2norm': numpy.sqrt( sum( domain.integrate( func**2., 'gauss6', coords, title=title ) ) )
    
    # boundary data
    velo = function.stack( [coords[2], 0., 0.] )
    trac = function.stack( [coords.normal()[2], 0., coords.normal()[0]] )
    # boundary data compatibility
    influx = (velo*coords.normal()).sum(-1)
    assert numpy.abs( domain.integrate( influx, geometry=coords, ischeme='gauss9', title='compatibility' ) ) < 1.e-5, 'int v.n = 0 condition violated.'
  
    # Matric/vector assembly
    ddomain = domain*domain
    iw = function.iwscale( coords, 2 )
    iweights = iw * function.opposite( iw ) * function.IWeights()
    x, y = coords, function.opposite( coords )
    kwargs = {'iweights':iweights, 'force_dense':True}
    rhs = 0.5*domain.integrate( (funcsp*velo).sum(-1), geometry=x, ischeme=ischeme['F'], title='bem[F]' ) \
        + ddomain.integrate( (funcsp*(K(x,y)*function.opposite(velo)).sum()).sum(),
          ischeme=ischeme['K'], title='bem[K]', **kwargs )
    mat = ddomain.integrate_symm( (funcsp*(V(x,y)*function.opposite(funcsp)[:,_,_,:]).sum()).sum(),
          ischeme=ischeme['V'], title='bem[V]', **kwargs )

    # Solve
    lhs = mat.solve( rhs, tol=1.e-8 )
    trach = funcsp.dot(lhs)
    tracp = domain.projection( trac, onto=funcsp, geometry=coords, ischeme='gauss6', title='project exact sol' )
    err = l2norm(trach-tracp, 'err')/l2norm(numpy.ones(1), 'surf')
    if visual:
      plot.writevtu( visual, domain.refine(2), coords, sdv=1.e-4 ) # , pointdata={'trac0':trac, 'trach':trach} )
      log.info( 'L2 err per unit surface area: %.2e' % err )
    assert err < tol, 'Traction computed in BEM example exceeds tolerance.'

  def test_InteriorProblem( self, visual=False ):
    'Interior Dirichlet on torus with b-splines.'
    # Topology and torus geometry
    domain, coords = mesh.rectilinear( 2*(range(-2,3),), periodic=(0,1) )
    cos, sin = function.cos, function.sin
    R, r = 3, 1
    phi, theta = .5*pi*coords
    torus = function.stack( [(r*cos(theta) + R)*cos(phi), (r*cos(theta) + R)*sin(phi), r*sin(theta)] )

    self._bemtest( domain, torus, {'F':'gauss3', 'K':'singular5', 'V':'singular3'}, tol=2e-2, visual='./torus.vtu' if visual else False )

  def test_SubdivisionGeometry( self, visual=False, case='sphere' ):
    'Interior Dirichlet on approximate sphere with subdivision surfaces.'
    tol = 1e-2 if case=='sphere' else 7e-2
    adaptive = 'singular3' # {'adaptive':'singular%i', 'qmin':2, 'qmax':17, 'step':2, 'TOL':1.e-5}
    ischemes = {'F':'gauss4', 'K':adaptive, 'V':adaptive}
    inputfile( case=case )
    domain, coords = mesh.blender( '%s.sdv'%case )
    self._bemtest( domain, coords, ischemes, tol=tol, visual='./%s.vtu'%case if visual else False )

def main():
  a = TestShearFlow()
  # a.test_InteriorProblem( visual=True )
  a.test_SubdivisionGeometry( visual=True, case='sphere' )

if __name__ == '__main__':
  util.run( main )

# vim:shiftwidth=2:foldmethod=indent:foldnestmax=2
