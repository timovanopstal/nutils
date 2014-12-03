#!/usr/bin/env python

from nutils import *
import numpy

class FiniteCellTestBase( object ):

  def __init__ ( self, ndims, nelems, maxrefine, voldec, surfdec ):
    domain, self.geom = mesh.rectilinear( (numpy.linspace(0,1,nelems+1),)*ndims )
    self.radius = numpy.sqrt( .5 )
    levelset = self.radius**2 - ( self.geom**2 ).sum()
    self.trimdomain, complement = domain.trim( levelset=levelset, maxrefine=maxrefine )
    V = 1.
    Vprev = 1. / (numpy.pi*self.radius)
    for idim in range( ndims ):
      S = Vprev * (2*numpy.pi*self.radius)
      Vprev = V
      V = S * (self.radius/(idim+1))
    self.volume = V / 2**ndims
    self.trimsurface = S / 2**ndims
    self.totalsurface = self.trimsurface + Vprev / (2**(ndims-1)) * ndims
    self.voldec = voldec
    self.surfdec = surfdec

  def all( self ):
    self.test_volume()
    self.test_surface()

  def test_volume( self ):
    vol = self.trimdomain.volume( self.geom )
    numpy.testing.assert_almost_equal( vol, self.volume, decimal=self.voldec )

  def test_divergence( self ):
    self.trimdomain.volume_check( self.geom, decimal=15 )
 
  def test_surface( self ):
    trimsurface = self.trimdomain.boundary['trimmed'].volume( self.geom )
    numpy.testing.assert_almost_equal( trimsurface, self.trimsurface, decimal=self.surfdec )
    totalsurface = self.trimdomain.boundary.volume( self.geom )
    numpy.testing.assert_almost_equal( totalsurface, self.totalsurface, decimal=self.surfdec )

class TestCircle( FiniteCellTestBase ):

  def __init__( self ):
    FiniteCellTestBase.__init__( self, ndims=2, nelems=2, maxrefine=5, voldec=3, surfdec=3 )
  
class TestSphere( FiniteCellTestBase ):

  def __init__( self ):
    FiniteCellTestBase.__init__( self, ndims=3, nelems=2, maxrefine=4, voldec=2, surfdec=2 )
  
class TestHierarchical():

  def test_hierarchical( self, makeplots=False ):

    # Topologies:
    # ref0    [  .  .  .  |  .  .  .  ]
    # ref1    [  .  .  .  |  .  |  .  ]
    # ref2    [  .  .  .  |  |  |  .  ]
    # trimmed [  .  .  .  |]

    ref0, geom = mesh.rectilinear( [[0,1,2]] )
    e1, e2 = ref0
    ref1 = ref0.refined_by( [e2] )
    e1, e2, e3 = ref1
    ref2 = ref1.refined_by( [e2] )

    basis = ref2.basis( 'std', degree=1 )
    assert basis.shape == (5,)
    x, y = ref2.elem_eval( [ geom[0], basis ], ischeme='bezier2', separate=False )
    assert numpy.all( y == .25 * numpy.array(
      [[4,0,0,0,0],
       [0,4,0,0,0],
       [0,3,2,0,4],
       [0,2,4,0,0],
       [0,0,0,4,0]] )[[0,1,1,2,2,3,3,4]] )

    if makeplots:
      with plot.PyPlot( 'basis' ) as plt:
        plt.plot( x, y )

    levelset = 1.125 - geom[0]
    trimmed, complement = ref2.trim( levelset, maxrefine=3 )
    trimbasis = trimmed.basis( 'std', degree=1 )
    x, y = trimmed.simplex.elem_eval( [ geom[0], trimbasis ], ischeme='bezier2', separate=False )
    assert numpy.all( y == .125 * numpy.array(
      [[8,0,0],
       [0,8,0],
       [0,7,4]] )[[0,1,1,2]] )

    if makeplots:
      with plot.PyPlot( 'basis' ) as plt:
        plt.plot( x, y )


if __name__ == '__main__':
  def hierarchical(): return TestHierarchical().test_hierarchical()
  def two_D(): return TestCircle().all()
  def three_D(): return TestSphere().all()
  util.run( hierarchical, two_D, three_D )
