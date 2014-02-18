#!/usr/bin/env python

from nutils import *
import numpy, copy

grid = numeric.linspace( 0., 1., 4 )
def inputfile():
  'Create approximate sphere: smallest closed subdivision surface with extraordinary points and only quad elements.'
  f = open( 'sphere.sdv', 'w' )
  f.write( 'Verts 26\n0. 0. -1.\n0. 0. 1.\n1. 0. 0.\n0. -1. 0.\n-1. 0. 0.\n0. 1. 0.\n0.75 0. -0.75\n0. 0.75 -0.75\n0.75 0.75 0.\n0. -0.75 -0.75\n0.75 -0.75 0.\n-0.75 0. -0.75\n-0.75 -0.75 0.\n-0.75 0.75 0.\n0.75 0. 0.75\n0. 0.75 0.75\n0. -0.75 0.75\n-0.75 0. 0.75\n0.555556 0.555556 -0.555556\n0.555556 -0.555556 -0.555556\n-0.555556 -0.555556 -0.555556\n-0.555556 0.555556 -0.555556\n0.555556 0.555556 0.555556\n0.555556 -0.555556 0.555556\n-0.555556 -0.555556 0.555556\n-0.555556 0.555556 0.555556\nRings 24 0\n4 3 18 8 5 7 0 6 2 20 9 19 10 11 21 13\n4 3 19 10 2 6 0 9 3 21 11 20 12 7 18 8\n4 3 20 12 3 9 0 11 4 18 7 21 13 6 19 10\n4 3 21 13 4 11 0 7 5 19 6 18 8 9 20 12\n4 3 22 8 2 14 1 15 5 24 17 25 13 16 23 10\n4 3 25 13 5 15 1 17 4 23 16 24 12 14 22 8\n4 3 24 12 4 17 1 16 3 22 14 23 10 15 25 13\n4 3 23 10 3 16 1 14 2 25 15 22 8 17 24 12\n4 3 18 7 0 6 2 8 5 23 14 22 15 10 19 9\n4 3 22 15 5 8 2 14 1 19 10 23 16 6 18 7\n4 3 23 16 1 14 2 10 3 18 6 19 9 8 22 15\n4 3 19 9 3 10 2 6 0 22 8 18 7 14 23 16\n4 3 19 6 0 9 3 10 2 24 16 23 14 12 20 11\n4 3 23 14 2 10 3 16 1 20 12 24 17 9 19 6\n4 3 24 17 1 16 3 12 4 19 9 20 11 10 23 14\n4 3 20 11 4 12 3 9 0 23 10 19 6 16 24 17\n4 3 20 9 0 11 4 12 3 25 17 24 16 13 21 7\n4 3 24 16 3 12 4 17 1 21 13 25 15 11 20 9\n4 3 25 15 1 17 4 13 5 20 11 21 7 12 24 16\n4 3 21 7 5 13 4 11 0 24 12 20 9 17 25 15\n4 3 22 14 1 15 5 8 2 21 7 18 6 13 25 17\n4 3 18 6 2 8 5 7 0 25 13 21 11 15 22 14\n4 3 21 11 0 7 5 13 4 22 15 25 17 8 18 6\n4 3 25 17 4 13 5 15 1 18 8 22 14 7 21 11\nFaceGroups 1\nGroupI 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23' )
  f.close()

class ConnectivityStructuredBase( object ):
  'Tests StructuredTopology.neighbor(), also handles periodicity.'

  def test_1DConnectivity( self ):
    domain = mesh.rectilinear( 1*(grid,), periodic=[0] if self.periodic else [] )[0]
    elem = domain.structure
    assert domain.neighbor( elem[0], elem[0] ) ==  0, 'Failed to identify codim 0 neighbors'
    assert domain.neighbor( elem[1], elem[2] ) ==  1, 'Failed to identify codim 1 neighbors'
    if self.periodic:
      assert domain.neighbor( elem[0], elem[2] ) ==  1, 'Failed to identify periodicity neighbors'
    else:
      assert domain.neighbor( elem[0], elem[2] ) == -1, 'Failed to identify non-neighbors'

  def test_2DConnectivity( self ):
    domain = mesh.rectilinear( 2*(grid,), periodic=[0] if self.periodic else [] )[0]
    elem = domain.structure
    assert domain.neighbor( elem[0,0], elem[0,0] ) ==  0, 'Failed to identify codim 0 neighbors'
    assert domain.neighbor( elem[1,1], elem[1,2] ) ==  1, 'Failed to identify codim 1 neighbors'
    assert domain.neighbor( elem[0,0], elem[1,1] ) ==  2, 'Failed to identify codim 2 neighbors'
    assert domain.neighbor( elem[1,1], elem[0,0] ) ==  2, 'Failed to identify codim 2 neighbors'
    if self.periodic:
      assert domain.neighbor( elem[2,1], elem[0,1] ) ==  1, 'Failed to identify periodicity neighbors'
      assert domain.neighbor( elem[2,1], elem[0,0] ) ==  2, 'Failed to identify periodicity neighbors'
    else:
      assert domain.neighbor( elem[2,1], elem[0,1] ) == -1, 'Failed to identify non-neighbors'

  def test_3DConnectivity( self ):
    domain = mesh.rectilinear( 3*(grid,), periodic=[0] if self.periodic else [] )[0]
    elem = domain.structure
    assert domain.neighbor( elem[1,1,1], elem[1,1,1] ) ==  0, 'Failed to identify codim 0 neighbors'
    assert domain.neighbor( elem[1,1,1], elem[1,1,2] ) ==  1, 'Failed to identify codim 1 neighbors'
    assert domain.neighbor( elem[1,1,1], elem[1,2,2] ) ==  2, 'Failed to identify codim 2 neighbors'
    assert domain.neighbor( elem[1,1,1], elem[2,2,2] ) ==  3, 'Failed to identify codim 3 neighbors'
    if self.periodic:
      assert domain.neighbor( elem[0,2,2], elem[2,2,2] ) ==  1, 'Failed to identify periodicity neighbors'
      assert domain.neighbor( elem[0,2,2], elem[2,1,2] ) ==  2, 'Failed to identify periodicity neighbors'
    else:
      assert domain.neighbor( elem[0,2,2], elem[2,2,2] ) == -1, 'Failed to identify non-neighbors'

class TestConnectivityStructured( ConnectivityStructuredBase ):
  periodic = False

class TestConnectivityStructuredPeriodic( ConnectivityStructuredBase ):
  periodic = True

class TestConnectivitySubdivision( object ):
  'Tests neighbor() and orientation for product domain of blender input.'
  def __init__( self ):
    inputfile()
    self.domain, self.coords = mesh.blender( 'sphere.sdv' )
    self.ddomain = self.domain * self.domain

  def orientations( self, visual=False ):
    'Test the orientation of product elements.'
    if not visual: raise NotImplementedError( 'Have only the visual inspection.' )

    idx = lambda el: self.domain.elements.index(el)
    def plotelem( fig, elem, style='.k:', verts=False ):
      'style is ordered: t, c, l'
      x, y, z = self.coords( elem, 'bezier2' ).T
      fig.plot( x, y, style[:2] )
      if verts:
        ha, va = ('left', 'bottom') if verts==1 else ('right', 'top')
        for ip, (px, py) in enumerate( zip( x, y ) ):
          fig.text( px, py, '%i'%ip, color=style[1], ha=ha, va=va )
        fig.text( numpy.mean(x), numpy.mean(y), 'elem%i'%verts, color=style[1],
            ha='center', va='center')
      for edge in elem.edges:
        x, y, z = self.coords( edge, 'bezier9' ).T
        fig.plot( x, y, style[1:] )
    def plotgrid( fig ):
      for elem in self.domain:
        if self.coords( elem, 'gauss1' )[0,2] < 0:
          continue # Skip bottom half of sphere
        plotelem( fig, elem )

    for elem in self.ddomain:
      if self.coords( elem.elem1, 'gauss1' )[0,2] < 0 or \
         self.coords( elem.elem2, 'gauss1' )[0,2] < 0:
        continue # Skip bottom half of sphere
      with plot.PyPlot( 'fig' ) as fig:
        inp = (idx(elem.elem1),idx(elem.elem2)) + elem.orientation
        fig.title( 'elem1:%i elem2:%i: %i/%i/%i'%inp )
        plotgrid( fig )
        plotelem( fig, elem.elem1, 'xr-', 1 )
        plotelem( fig, elem.elem2, '+g-', 2 )

class TestStructure2D( object ):
  'Test coordinate evaluation for StructuredTopology.'

  def verify_connectivity( self, structure, geom ):
    (e00,e01), (e10,e11) = structure

    geom = geom.compiled()

    a0 = geom.eval( e00, numeric.array([0,1]) )
    a1 = geom.eval( e01, numeric.array([0,0]) )
    numpy.testing.assert_array_almost_equal( a0, a1 )

    b0 = geom.eval( e10, numeric.array([1,1]) )
    b1 = geom.eval( e11, numeric.array([1,0]) )
    numpy.testing.assert_array_almost_equal( b0, b1 )

    c0 = geom.eval( e00, numeric.array([1,0]) )
    c1 = geom.eval( e10, numeric.array([0,0]) )
    numpy.testing.assert_array_almost_equal( c0, c1 )

    d0 = geom.eval( e01, numeric.array([1,1]) )
    d1 = geom.eval( e11, numeric.array([0,1]) )
    numpy.testing.assert_array_almost_equal( d0, d1 )

    x00 = geom.eval( e00, numeric.array([1,1]) )
    x01 = geom.eval( e01, numeric.array([1,0]) )
    x10 = geom.eval( e10, numeric.array([0,1]) )
    x11 = geom.eval( e11, numeric.array([0,0]) )
    numpy.testing.assert_array_almost_equal( x00, x01 )
    numpy.testing.assert_array_almost_equal( x10, x11 )
    numpy.testing.assert_array_almost_equal( x00, x11 )

  def testMesh( self ):
    domain, geom = mesh.rectilinear( [[-1,0,1]]*2 )
    self.verify_connectivity( domain.structure, geom )

  def testBoundaries( self ):
    domain, geom = mesh.rectilinear( [[-1,0,1]]*3 )
    for grp in 'left', 'right', 'top', 'bottom', 'front', 'back':
      bnd = domain.boundary[grp]
      self.verify_connectivity( bnd.structure, geom )
      xn = bnd.elem_eval( geom.dotnorm(geom), ischeme='gauss1', separate=False )
      numpy.testing.assert_array_less( 0, xn, 'inward pointing normals' )

class TestTopologyGlueing( object ):
  'Test glueing of compatible topologies along prespecified boundary.'

  def __init__( self ):
    'Create half dome geometry for glueing.'
    # Aliases
    pi, sqrt, sin, cos, abs = numeric.pi, function.sqrt, function.sin, function.cos, function.abs
    grid = numeric.linspace( -.25*pi, .25*pi, 5 )

    # Half dome
    self.topo0, (xi, eta) = mesh.rectilinear( 2*(grid,) )
    x, y = sqrt(2)*sin(xi)*cos(eta), sqrt(2)*sin(eta)*cos(xi)
    self.geom0 = function.stack( [x, y, abs(1-x**2-y**2)] ) # Don't take sqrt, upsets BEM conv.

    # Plane, rotated to ensure singular-integration-compatible connectivity
    self.topo1, (xi, eta) = mesh.rectilinear( 2*(grid,) )
    for elem in self.topo1: # relabel vertices
      elem.vertices = tuple( vertex+"/" for vertex in elem.vertices )
    x, y = sin(xi)*cos(eta)-sin(eta)*cos(xi), sin(xi)*cos(eta)+sin(eta)*cos(xi)
    self.geom1 = function.stack( [x, -y, 0] ) # minus to get normal downwards

    # Merged function object and coordinate function
    # For one single merged coordinate system we need the cascades to match up, so we project, this
    # turns out to preserve the matching of element edges up to errors of 10^-4 at the vertices.
    splines_on = lambda topo: topo.splinefunc(degree=4)
    self.funcsp = function.merge( [splines_on(self.topo0), splines_on(self.topo1)] ).vector(3)
    dofs = self.topo0.project( self.geom0, self.funcsp, self.geom0, exact_boundaries=True, ischeme='gauss8' ) \
         | self.topo1.project( self.geom1, self.funcsp, self.geom1, exact_boundaries=True, ischeme='gauss8' )
    self.geom = self.funcsp.dot( dofs )

    # Glue boundary definition
    self.topo0.boundary.groups['__glue__'] = self.topo0.boundary
    self.topo1.boundary.groups['__glue__'] = self.topo1.boundary
    self.topo = topology.glue( self.topo0, self.topo1, self.geom, tol=1e-4 )

  def plot_gauss_on_circle( self, elem, ischeme='singular2', title='' ):
    'Given a product element on our 4x4 circular domain (see __init__), plot gauss points'
    dom, coo = mesh.rectilinear( 2*([0,1],) )
    circumf = dom.boundary.elem_eval( coo, 'bezier5', separate=True )
    quadpoints = elem.eval( ischeme )[0]
    with plot.PyPlot( 'quad', figsize=(6,5) ) as fig:
      # Glued grids
      for partition, style in (('__master__', 'r-'), ('__slave__', 'g-')):
        for cell in self.topo.groups[partition]:
          pts = self.geom( cell, circumf )
          fig.plot( pts[:,0], pts[:,1], style )
      # Quad points on elem pair
      for element, points in zip( (elem.elem1, elem.elem2), (quadpoints[:,:2], quadpoints[:,2:]) ):
        style = 'rx' if element in self.topo.groups['__master__'] else 'g+'
        pts = self.geom( element, points )
        fig.plot( pts[:,0], pts[:,1], style )

      fig.title( title + ' | n:%i, t:%i, %i'%elem.orientation )

  def _integrate( self, func, ecoll, qset=range(1,9), qmax=16, slopes=None, plot_quad_points=False ):
    '''Test convergence of approximation on all product element types.
    I: func,    integrand,
       ddomain, product domain over which to perform integration test,
       dcoords, tuple of corresponding coordinate functions,
       ecoll,   product elements over which to perform integration test,
       qset,    set of quadrature orders, length (1,2, >2) determines type of test,
       qmax,    reference quadrature level,
       slopes,  expected rate of convergence.'''
    devel = len(qset) > 2

    # This could be handled underwater by Topology.integrate only if geom can be glued.
    iw = function.iwscale( self.geom, 2 )
    iweights = iw * function.opposite( iw ) * function.IWeights()

    # integrands and primitives
    for neighbor, elems in enumerate( ecoll ):
      if devel: errs = {}
      for key, elem in elems.iteritems():
        topo = topology.UnstructuredTopology( [elem], ndims=2 )
        integrate = lambda q: topo.integrate( func, iweights=iweights, ischeme='singular%i'%q )
        F = integrate( qmax )

        if devel:
          # Devel mode (default), visual check of convergence
          errs[key] = []
          for q in qset:
            Fq = integrate( q )
            errs[key].append( numeric.abs(F/Fq-1) )

        elif len(qset) == 1:
          # Test assertions on exact quadrature
          Fq = integrate( qset[0] )
          err = numeric.abs(F/Fq-1)
          assert err < 1.e-12, 'Nonexact quadrature, err = %.1e' % err

        elif len(qset) == 2:
          # Test assertions on convergence rate of quadrature
          q0, q1 = tuple( qset )
          F0 = integrate( q0 )
          F1 = integrate( q1 )
          err0 = numeric.abs(F/F0-1)
          err1 = numeric.abs(F/F1-1)
          slope = numeric.log10(err1/err0)/(q1-q0)
          assert slope <= (-2. if slopes is None else slopes[neighbor]) or err1 < 1.e-12, \
              'Insufficient quadrature convergence (is func analytic?), slope = %.2f' % slope

        else:
          raise ValueError( 'Range of quadrature orders should contain >=1 value.' )

      if devel and len(elems):
        with plot.PyPlot( 'conv' ) as fig:
          for val in errs.itervalues():
            fig.semilogy( qset, val )
          i = len(qset)//2
          slope = fig.slope_triangle( qset[i::i-1][::-1], val[i::i-1][::-1], slopefmt='{0:.2f}' )
          fig.title( 'n-type: %i'%(-1 if neighbor is 3 else neighbor) )

  def test_2DNodeRelabelingCorrect( self, visual=False ):
    'Topology glueing should not raise any errors.'
    # 0. Test if glue passes without errors: done in __init__(), is the resulting glued topology up to specs?
    assert len(self.topo) == 32
    keys_provided = set( self.topo.groups )
    keys_required = set(['master', 'slave'])
    assert keys_provided == keys_required, 'Something went awry with copying groups into union topology.'

    bkeys_provided = set( self.topo.boundary.groups )
    bkeys_required = set(['master', 'master_bottom', 'master_left', 'master_right', 'master_top',
                          'slave', 'slave_bottom', 'slave_left', 'slave_right', 'slave_top' ])
    assert bkeys_provided == bkeys_required, 'Something went awry with copying boundary groups into union topology.'

    # 1. The connectivity should still be correct, cf. test_quadrature.TestSingularQuadrature.test_connectivity
    elem = self.topo.elements
    neighbor_tests = [[(0,1),1], [(0,2),1], [(0,0),2], [(0,3),2], [(1,1),-1]]
    index = lambda alpha, offset=16: offset+numpy.ravel_multi_index( alpha, (4,4) ) # Return index of multiindex
    for alpha, n in neighbor_tests:
      assert elem[3].neighbor( elem[index(alpha)] ) == n, \
          'Failed to identify codim %i neighbor over boundary' % n

    # 2. Orientation information should also be preserved, cf. test_quadrature.TestSingularQuadrature.test_orientations
    ddom = topology.UnstructuredTopology( elem[3:4], ndims=2 ) * \
           topology.UnstructuredTopology( elem[16:], ndims=2 )
    orientations = {str( elem[index((0,3))] ): [2, (0,0), (0,7), (7,0), (0,7)],
                    str( elem[index((0,2))] ): [1, (3,7), (7,3)],
                    str( elem[index((0,1))] ): [1, (2,7), (6,3)],
                    str( elem[index((0,0))] ): [2, (2,3), (2,6), (5,3), (5,6)]}
    ecoll = [{}, {}, {}, {}] # Collect product elements to be used in integration test below.
    for alpha, pelem in enumerate( ddom ):
      orientation = orientations.get( str( pelem.elem2 ), [-1, (0,0)] )
      ecoll[orientation[0]][pelem.__repr__()] = pelem
      if visual: self.plot_gauss_on_circle( pelem )
      assert pelem.orientation[0] == orientation[0], 'Incorrect neighbor type.'
      assert pelem.orientation[1:] in orientation[1:], 'Incorrect transformation.'

    # 3. Integration should converge, cf. test_quadrature.TestSingularQuadrature.test_stronglysingularfunc
    kwargs = {'qset':(2,4), 'qmax':8, 'slopes':(-0.0, -1.0, -1.0, -1.0)}
    ddom = topology.UnstructuredTopology( elem[:16], ndims=2 ) * \
           topology.UnstructuredTopology( elem[16:], ndims=2 )
    func = function.norm2( self.geom-function.opposite(self.geom) )**-2
    if visual:
      kwargs.update( {'qset':range(1,10), 'qmax':16, 'plot_quad_points':True} )
      # Fill inspection set ecoll and count number of product elements for each neighbor type
      ecoll = [{}, {}, {}, {}]
      for pelem in ddom: ecoll[pelem.orientation[0]][pelem.__repr__()] = pelem
      for i, coll in enumerate(ecoll): log.warning( 'n: %i, #el: %i' % (i, len(coll)) )
    self._integrate( func, ecoll, **kwargs )

  def test_2DNodeRelabelingBigMaster( self ):
    'This should raise an AssertionError, as there are too many master elements.'
    self.topo0.boundary.groups['__glue__'] = self.topo0.boundary + self.topo0[:1,:1].boundary['right'] # For some strange reason deepcopy skips boundary['__glue__']
    numpy.testing.assert_raises( AssertionError, topology.glue, self.topo0, self.topo1, self.geom )

  def test_2DNodeRelabelingBigSlave( self ):
    'This should raise an AssertionError, as there are too many slave elements.'
    self.topo1.boundary.groups['__glue__'] = self.topo1.boundary + self.topo1[:1,:1].boundary['right']
    numpy.testing.assert_raises( AssertionError, topology.glue, self.topo0, self.topo1, self.geom )

  def StokesBEM( self, visual=False ):
    'The singular integration scheme depends on the correct functioning of Topology.glue().'
    # Aliases and definitions
    pi, sqrt, sin, cos, abs = numeric.pi, function.sqrt, function.sin, function.cos, function.abs
    def V( x, y ):
      rInv = function.norm2( x-y )**-1.
      return 0.125*pi**-1. * (function.eye(3)*rInv + (x-y)[:,_]*(x-y)[_,:]*rInv**3)
    def K( x, y ):
      rInv = function.norm2( x-y )**-1.
      return 0.75*pi**-1. * (x-y)[:,_]*(x-y)[_,:] * ((x-y)*y.normal()).sum() * rInv**5
    l2norm = lambda func: sqrt( sum( self.topo.integrate( func**2., 'gauss10', self.geom ) ) )

    # Boundary data
    velo = function.stack( [self.geom[2], 0., 0.] )
    trac = function.stack( [self.geom.normal()[2], 0., self.geom.normal()[0]] )

    # Matrix/vector assembly (integration schemes optimized)
    prod_topo = self.topo*self.topo
    iw = function.iwscale(self.geom,2)
    iweights = iw * function.opposite(iw) * function.IWeights()
    x, y = self.geom, function.opposite(self.geom)
    kwargs = {'iweights':iweights, 'ischeme':'singular6', 'force_dense':True}
    integrand = (self.funcsp*(V(x,y)*function.opposite(self.funcsp)[:,_,_,:]).sum()).sum()
    mat = prod_topo.integrate_symm( integrand, title='int[mat]', **kwargs )
    integrand = (self.funcsp*(K(x,y)*function.opposite(velo)).sum()).sum()
    vec = 0.5 * self.topo.integrate( (self.funcsp*velo).sum(), geometry=x, ischeme='gauss6' ) \
        + prod_topo.integrate( integrand, title='int[vec]', **kwargs )

    # Solve
    sol = mat.solve( vec, tol=0, symmetric=True )
    trach = self.funcsp.dot( sol )
    if visual:
      plot.writevtu( './dome.vtu', self.topo.refine(2), self.geom, sdv=1.e-3,
          pointdata={'trac0':trac, 'trach':trach} )
    relerr = l2norm(trach-trac)/l2norm(trac)
    log.info( 'rel. L2 err: %.2e' % relerr )
    assert relerr < 1.e-2, 'Traction computed in BEM example exceeds tolerance.'

def visualinspect():
  'Visual inspection of StokesBEM test case.'
  # visual = TestConnectivitySubdivision()
  # visual.orientations( visual=True )
  visual = TestTopologyGlueing()
  visual.StokesBEM( visual=True )

if __name__ == '__main__':
  util.run( visualinspect )

# vim:shiftwidth=2:foldmethod=indent:foldnestmax=2
