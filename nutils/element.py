from . import log, util, cache, numeric, transform, function, _
import warnings, numpy

PrimaryVertex = str
HalfVertex = lambda vertex1, vertex2, xi=.5: '%s<%.3f>%s' % ( (vertex1,xi,vertex2) if vertex1 < vertex2 else (vertex2,1-xi,vertex1) )
ProductVertex = lambda *vertices: ','.join( vertices )

class Element( object ):
  '''Element base class.

  Represents the topological shape.'''

  __slots__ = 'vertices', 'ndims', 'parent', 'context', 'interface', 'root_transform'

  def __init__( self, ndims, vertices, parent=None, context=None, interface=None ):
    'constructor'

    #assert all( isinstance(vertex,Vertex) for vertex in vertices )
    self.vertices = tuple(vertices)
    self.ndims = ndims
    self.parent = parent
    self.context = context
    self.interface = interface

    if parent:
      pelem, trans = parent
      self.root_transform = pelem.root_transform * trans
    else:
      self.root_transform = transform.Identity( self.ndims )

  def __mul__( self, other ):
    'multiply elements'

    return ProductElement( self, other )

  def neighbor( self, other ):
    'level of neighborhood; 0=self'

    if self == other:
      return 0
    ncommon = len( set(self.vertices) & set(other.vertices) )
    return self.neighbormap[ ncommon ]

  def eval( self, where, cache=lambda f, *args: f(*args) ):
    'get points'

    if isinstance( where, str ):
      points, weights = cache( self.getischeme, self.ndims, where )
    else:
      points = where
      weights = None
    return points, weights

  def __str__( self ):
    'string representation'

    return '%s(%s)' % ( self.__class__.__name__, self.vertices )

  __repr__ = __str__

  def __hash__( self ):
    'hash'

    return hash(str(self))

  def __eq__( self, other ):
    'equal'

    return self is other or (
          self.__class__ == other.__class__
      and self.vertices == other.vertices
      and self.parent == other.parent
      and self.context == other.context
      and self.interface == other.interface )

  def intersected( self, levelset, lscheme, evalrefine=0 ):
    '''check levelset intersection:
      +1 for levelset > 0 everywhere
      -1 for levelset < 0 everywhere
       0 for intersected element'''

    levelset = function.ascompiled( levelset )
    elems = iter( [self] )
    for irefine in range(evalrefine):
      elems = ( child for elem in elems for child in elem.children )

    inside = numeric.greater( levelset.eval( elems.next(), lscheme ), 0 )
    if inside.all():
      for elem in elems:
        inside = numeric.greater( levelset.eval( elem, lscheme ), 0 )
        if not inside.all():
          return 0
      return 1
    elif not inside.any():
      for elem in elems:
        inside = numeric.greater( levelset.eval( elem, lscheme ), 0 )
        if inside.any():
          return 0
      return -1
    return 0

  def trim( self, levelset, maxrefine, lscheme, finestscheme, evalrefine ):
    'trim element along levelset'

    levelset = function.ascompiled( levelset )
    intersected = self.intersected( levelset, lscheme, evalrefine )

    if intersected > 0:
      return self

    if intersected < 0:
      return None

    parent = self, transform.Identity(self.ndims)
    return TrimmedElement( elem=self, vertices=self.vertices, levelset=levelset, maxrefine=maxrefine, lscheme=lscheme, finestscheme=finestscheme, evalrefine=evalrefine, parent=parent )

  def get_simplices ( self, maxrefine ):
    'divide in simple elements'

    return self,

  def get_trimmededges ( self, maxrefine ):
    return []

class ProductElement( Element ):
  'element product'

  __slots__ = 'elem1', 'elem2'

  @staticmethod
  def getslicetransforms( ndims1, ndims2 ):
    ndims = ndims1 + ndims2
    slice1 = transform.Slice( ndims, 0, ndims1 )
    slice2 = transform.Slice( ndims, ndims1, ndims )
    return slice1, slice2

  def __init__( self, elem1, elem2 ):
    'constructor'

    self.elem1 = elem1
    self.elem2 = elem2
    slice1, slice2 = self.getslicetransforms( elem1.ndims, elem2.ndims )
    iface1 = elem1, slice1
    iface2 = elem2, slice2
    vertices = [] # TODO [ ProductVertex(vertex1,vertex2) for vertex1 in elem1.vertices for vertex2 in elem2.vertices ]
    Element.__init__( self, ndims=elem1.ndims+elem2.ndims, vertices=vertices, interface=(iface1,iface2) )
    self.root_transform = transform.Scale( numeric.array([elem1.root_transform.det * elem2.root_transform.det]) ) # HACK

  @staticmethod
  def get_tri_bem_ischeme( ischeme, neighborhood ):
    'Some cached quantities for the singularity quadrature scheme.'
    points, weights = QuadElement.getischeme( ndims=4, where=ischeme )
    eta1, eta2, eta3, xi = points.T
    if neighborhood == 0:
      temp = xi*eta1*eta2*eta3
      pts0 = xi*eta1*(1 - eta2)
      pts1 = xi - pts0
      pts2 = xi - temp
      pts3 = xi*(1 - eta1)
      pts4 = pts0 + temp
      pts5 = xi*(1 - eta1*eta2)
      pts6 = xi*eta1 - temp
      points = numeric.asarray(
        [[1-xi,   1-pts2, 1-xi,   1-pts5, 1-pts2, 1-xi  ],
         [pts1, pts3, pts4, pts0, pts6, pts0],
         [1-pts2, 1-xi,   1-pts5, 1-xi,   1-xi,   1-pts2],
         [pts3, pts1, pts0, pts4, pts0, pts6]]).reshape( 4, -1 ).T
      points = numeric.asarray( points * [-1,1,-1,1] + [1,0,1,0] ) # flipping in x -GJ
      weights = numeric.concatenate( 6*[xi**3*eta1**2*eta2*weights] )
    elif neighborhood == 1:
      A = xi*eta1
      B = A*eta2
      C = A*eta3
      D = B*eta3
      E = xi - B
      F = A - B
      G = xi - D
      H = B - D
      I = A - D
      points = numeric.asarray(
        [[1-xi, 1-xi, 1-E,  1-G,  1-G ],
         [C,  G,  F,  H,  I ],
         [1-E,  1-G,  1-xi, 1-xi, 1-xi],
         [F,  H,  D,  A,  B ]] ).reshape( 4, -1 ).T
      temp = xi*A
      weights = numeric.concatenate( [A*temp*weights] + 4*[B*temp*weights] )
    elif neighborhood == 2:
      A = xi*eta2
      B = A*eta3
      C = xi*eta1
      points = numeric.asarray(
        [[1-xi, 1-A ],
         [C,  B ],
         [1-A,  1-xi],
         [B,  C ]] ).reshape( 4, -1 ).T
      weights = numeric.concatenate( 2*[xi**2*A*weights] )
    else:
      assert neighborhood == -1, 'invalid neighborhood %r' % neighborhood
      points = numeric.asarray([ eta1*eta2, 1-eta2, eta3*xi, 1-xi ]).T
      weights = eta2*xi*weights
    return points, weights
  
  @staticmethod
  def get_quad_bem_ischeme( ischeme, neighborhood ):
    'Some cached quantities for the singularity quadrature scheme.'
    points, weights = QuadElement.getischeme( ndims=4, where=ischeme )
    eta1, eta2, eta3, xi = points.T
    if neighborhood == 0:
      xe = xi*eta1
      A = (1 - xi)*eta3
      B = (1 - xe)*eta2
      C = xi + A
      D = xe + B
      points = numeric.asarray(
        [[A, B, A, D, B, C, C, D],
         [B, A, D, A, C, B, D, C],
         [C, D, C, B, D, A, A, B],
         [D, C, B, C, A, D, B, A]]).reshape( 4, -1 ).T
      weights = numeric.concatenate( 8*[xi*(1-xi)*(1-xe)*weights] )
    elif neighborhood == 1:
      ox = 1 - xi
      A = xi*eta1
      B = xi*eta2
      C = ox*eta3
      D = C + xi
      E = 1 - A
      F = E*eta3
      G = A + F
      points = numeric.asarray(
        [[D,  C,  G,  G,  F,  F ],
         [B,  B,  B,  xi, B,  xi],
         [C,  D,  F,  F,  G,  G ],
         [A,  A,  xi, B,  xi, B ]]).reshape( 4, -1 ).T
      weights = numeric.concatenate( 2*[xi**2*ox*weights] + 4*[xi**2*E*weights] )
    elif neighborhood == 2:
      A = xi*eta1
      B = xi*eta2
      C = xi*eta3
      points = numeric.asarray(
        [[xi, A,  A,  A ], 
         [A,  xi, B,  B ],
         [B,  B,  xi, C ], 
         [C,  C,  C,  xi]]).reshape( 4, -1 ).T
      weights = numeric.concatenate( 4*[xi**3*weights] )
    else:
      assert neighborhood == -1, 'invalid neighborhood %r' % neighborhood
    return points, weights

  @staticmethod
  def concat( ischeme1, ischeme2 ):
    coords1, weights1 = ischeme1
    coords2, weights2 = ischeme2
    if weights1 is not None:
      assert weights2 is not None
      weights = numeric.asarray( ( weights1[:,_] * weights2[_,:] ).ravel() )
    else:
      assert weights2 is None
      weights = None
    npoints1,ndims1 = coords1.shape  
    npoints2,ndims2 = coords2.shape 
    coords = numeric.empty( [ coords1.shape[0], coords2.shape[0], ndims1+ndims2 ] )
    coords[:,:,:ndims1] = coords1[:,_,:]
    coords[:,:,ndims1:] = coords2[_,:,:]
    coords = numeric.asarray( coords.reshape(-1,ndims1+ndims2) )
    return coords, weights
  
  @property
  def orientation( self ):
    '''Neighborhood of elem1 and elem2 and transformations to get mutual overlap in right location
    O: neighborhood,  as given by Element.neighbor(),
       transf1,       required rotation of elem1 map: {0:0, 1:pi/2, 2:pi, 3:3*pi/2},
       transf2,       required rotation of elem2 map (is indep of transf1 in UnstructuredTopology.'''
    neighborhood = self.elem1.neighbor( self.elem2 )
    common_vertices = list( set(self.elem1.vertices) & set(self.elem2.vertices) )
    vertices1 = [self.elem1.vertices.index( ni ) for ni in common_vertices]
    vertices2 = [self.elem2.vertices.index( ni ) for ni in common_vertices]
    if neighborhood == 0:
      # test for strange topological features
      assert self.elem1 == self.elem2, 'Topological feature not supported: try refining here, possibly periodicity causes elems to touch on both sides.'
      transf1 = transf2 = 0
    elif neighborhood == -1:
      transf1 = transf2 = 0
    elif isinstance( self.elem1, QuadElement ):
      # define local map rotations
      if neighborhood==1:
        trans = [0,2], [2,3], [3,1], [1,0], [2,0], [3,2], [1,3], [0,1]
      elif neighborhood==2:
        trans = [0], [2], [3], [1]
      else:
        raise ValueError( 'Unknown neighbor type %i' % neighborhood )
      transf1 = trans.index( vertices1 )
      transf2 = trans.index( vertices2 )
    elif isinstance( self.elem1, TriangularElement ):
      raise NotImplementedError( 'Pending completed implementation and verification.' )
      # define local map rotations
      if neighborhood==1:
        trans = [0,1], [1,2], [0,2]
      elif neighborhood==2:
        trans = [0], [1], [2]
      else:
        raise ValueError( 'Unknown neighbor type %i' % neighborhood )
      transf1 = trans.index( vertices1 )
      transf2 = trans.index( vertices2 )
    else:
      raise NotImplementedError( 'Reorientation not implemented for element of class %s' % type(self.elem1) )
    return neighborhood, transf1, transf2

  @staticmethod
  def singular_ischeme_tri( orientation, ischeme ):
    neighborhood, transf1, transf2 = orientation
    points, weights = ProductElement.get_tri_bem_ischeme( ischeme, neighborhood )
    transfpoints = points#numeric.empty( points.shape )
    #   transfpoints[:,0] = points[:,0] if transf1 == 0 else \
    #                       points[:,1] if transf1 == 1 else \
    #                     1-points[:,0] if transf1 == 2 else \
    #                     1-points[:,1]
    #   transfpoints[:,1] = points[:,1] if transf1 == 0 else \
    #                     1-points[:,0] if transf1 == 1 else \
    #                     1-points[:,1] if transf1 == 2 else \
    #                       points[:,0]
    #   transfpoints[:,2] = points[:,2] if transf2 == 0 else \
    #                       points[:,3] if transf2 == 1 else \
    #                     1-points[:,2] if transf2 == 2 else \
    #                     1-points[:,3]
    #   transfpoints[:,3] = points[:,3] if transf2 == 0 else \
    #                     1-points[:,2] if transf2 == 1 else \
    #                     1-points[:,3] if transf2 == 2 else \
    #                       points[:,2]
    return numeric.asarray( transfpoints ), numeric.asarray( weights )
    
  @staticmethod
  def singular_ischeme_quad( orientation, ischeme ):
    neighborhood, transf1, transf2 = orientation
    points, weights = ProductElement.get_quad_bem_ischeme( ischeme, neighborhood )
    transfpoints = numeric.empty( points.shape )
    def flipxy( points, orientation ):
      x, y = points[:,0], points[:,1]
      tx = x if orientation in (0,1,6,7) else 1-x
      ty = y if orientation in (0,3,4,7) else 1-y
      return function.stack( (ty, tx) if orientation%2 else (tx, ty), axis=1 )
    transfpoints[:,:2] = flipxy( points[:,:2], transf1 )
    transfpoints[:,2:] = flipxy( points[:,2:], transf2 )
    return numeric.asarray( transfpoints ), numeric.asarray( weights )
    
  def eval( self, where, cache=lambda f, *args: f(*args) ):
    'get integration scheme'
    
    if where.startswith( 'singular' ):
      assert type(self.elem1) == type(self.elem2), 'mixed element-types case not implemented'
      assert self.elem1.ndims == 2 and self.elem2.ndims == 2, 'singular quadrature only for bivariate surfaces'
      gauss = 'gauss%d'% (int(where[8:])*2-2)
      if isinstance( self.elem1, QuadElement ):
        xw = cache( self.singular_ischeme_quad, self.orientation, gauss )
      elif isinstance( self.elem1, TriangularElement ):
        if self.elem1 == self.elem2:
          xw = cache( self.get_tri_bem_ischeme, gauss, neighborhood=0 )
        else:
          xw = self.concat( self.elem1.eval(gauss,cache), self.elem2.eval(gauss,cache) )
      else:
        raise Exception, 'invalid element type %r' % type(self.elem1)
    else:
      where1, where2 = where.split( '*' ) if '*' in where else ( where, where )
      xw = self.concat( self.elem1.eval(where1,cache), self.elem2.eval(where2,cache) )
    return xw

class TrimmedElement( Element ):
  'trimmed element'

  __slots__ = 'elem', 'levelset', 'maxrefine', 'lscheme', 'finestscheme', 'evalrefine'

  def __init__( self, elem, levelset, maxrefine, lscheme, finestscheme, evalrefine, parent, vertices ):
    'constructor'

    assert not isinstance( elem, TrimmedElement )
    self.elem = elem
    self.levelset = function.ascompiled( levelset )
    self.maxrefine = maxrefine
    self.lscheme = lscheme
    self.finestscheme = finestscheme if finestscheme != None else 'simplex1'
    self.evalrefine = evalrefine

    Element.__init__( self, ndims=elem.ndims, vertices=vertices, parent=parent )

  def eval( self, ischeme, cache=lambda f, *args: f(*args) ):
    'get integration scheme'

    assert isinstance( ischeme, str )

    if ischeme[:7] == 'contour':
      n = int(ischeme[7:] or 0)
      points, weights = self.elem.eval( 'contour{}'.format(n), cache )
      inside = numeric.greater_equal( self.levelset.eval( self.elem, points ), 0 )
      points = points[inside]
      return points, None

    if self.maxrefine <= 0:
      if self.finestscheme.startswith('simplex'):

        points  = []
        weights = []

        for simplex in self.get_simplices( 0 ):

          spoints, sweights = simplex.eval( ischeme, cache )
          pelem, trans = simplex.parent

          assert pelem is self 

          points.append( trans.apply( spoints ) )
          weights.append( sweights * trans.det )

        if len(points) == 0:
          points = numeric.zeros((0,self.ndims))
          weights = numeric.zeros((0,))
        else:
          points = numeric.concatenate( points, axis=0 )
          weights = numeric.concatenate( weights )

        return points, weights

      else:
        
        if self.finestscheme.endswith( '.all' ):
          points, weights = self.elem.eval( self.finestscheme[:-4], cache )
        elif self.finestscheme.endswith( '.none' ):
          points, weights = self.elem.eval( self.finestscheme[:-5], cache )
          inside = numeric.zeros_like(weights,dtype=bool) # what is .none supposed to do? -GJ
        else:  
          points, weights = self.elem.eval( self.finestscheme, cache )
          inside = numeric.greater( self.levelset.eval( self.elem, points ), 0 )
        points = points[inside]
        if weights is not None:
          weights = weights[inside]
        return points, weights

    allcoords = []
    allweights = []
    for child in self.children:
      if child is None:
        continue
      points, weights = child.eval( ischeme, cache )
      pelem, trans = child.parent
      assert pelem == self
      allcoords.append( trans.apply(points) )
      allweights.append( weights * trans.det )

    coords = numeric.concatenate( allcoords, axis=0 )
    weights = numeric.concatenate( allweights, axis=0 )
    return coords, weights

  @property
  def children( self ):
    'all 1x refined elements'

    children = []
    for ielem, child in enumerate( self.elem.children ):
      isect = child.intersected( self.levelset, self.lscheme, self.evalrefine-1 )
      pelem, trans = child.parent
      parent = self, trans
      if isect < 0:
        child = None
      elif isect > 0:
        child = QuadElement( vertices=child.vertices, ndims=self.ndims, parent=parent )
      else:
        child = TrimmedElement( vertices=child.vertices, elem=child, levelset=self.levelset, maxrefine=self.maxrefine-1, lscheme=self.lscheme, finestscheme=self.finestscheme, evalrefine=self.evalrefine-1, parent=parent )
      children.append( child )
    return tuple( children )

  def edge( self, iedge ):
    'edge'

    # TODO fix trimming of edges once refine/edge operations commute
    edge = self.elem.edge( iedge )
    pelem, trans = edge.context

    # transform = self.elem.edgetransform( self.ndims )[ iedge ]
    return QuadElement( vertices=edge.vertices, ndims=self.ndims-1, context=(self,trans) )

  def get_simplices ( self, maxrefine ):
    'divide in simple elements'

    if maxrefine > 0 or self.evalrefine > 0:
      return [ simplex for child in filter(None,self.children) for simplex in child.get_simplices( maxrefine=maxrefine-1 ) ]

    simplices, trimmededges = self.triangulate()

    return simplices

  def get_trimmededges ( self, maxrefine ):

    if maxrefine > 0 or self.evalrefine > 0:
      return [ trimmededge for child in filter(None,self.children) for trimmededge in child.get_trimmededges( maxrefine=maxrefine-1 ) ]

    simplices, trimmededges = self.triangulate()

    return trimmededges

  def triangulate ( self ):

    assert self.finestscheme.startswith('simplex'), 'Expected simplex scheme'
    order = int(self.finestscheme[7:])

    lvltol  = numeric.spacing(100)
    xistol  = numeric.spacing(100)
    ischeme = self.elem.getischeme( self.elem.ndims, 'bezier2' )
    where   = numeric.greater( self.levelset.eval( self.elem, ischeme ), -lvltol )
    points  = ischeme[0][where]
    vertices   = numeric.array(self.vertices)[where].tolist()
    norig   = sum(where)

    if not where.any():
      return []

    if where.all():
	    lines = []
    else:		
    	lines = self.elem.ribbons

    for line in lines:
      
      ischeme = line.getischeme( line.ndims, 'bezier'+str(order+1) )
      vals    = self.levelset.eval( line, ischeme )
      pts     = ischeme[0]
      where   = numeric.greater( vals, -lvltol )
      zeros   = numeric.less(vals, lvltol) & where

      if order == 1:

        if where[0] == where[1] or zeros.any():
          continue
        xi = vals[0]/(vals[0]-vals[1])

      elif order == 2:

        #Check whether the levelset is linear
        if abs(2*vals[1]-vals[0]-vals[2]) < lvltol:

          if where[0] == where[2] or zeros[0] or zeros[2]:
            continue
          xi = vals[0]/(vals[0]-vals[2])
          
        else:

          disc = vals[0]**2+(-4*vals[1]+vals[2])**2-2*vals[0]*(4*vals[1]+vals[2])

          if disc < -lvltol or abs(disc) < lvltol:
            #No intersections or minimum at zero
            continue
          else:
            #Two intersections
            num2 = numeric.sqrt( disc )
            num1 = 3*vals[0]-4*vals[1]+vals[2]
            denr = 1./(4*(vals[0]-2*vals[1]+vals[2]))

            xis = [(num1-num2)*denr,\
                   (num1+num2)*denr ]

            intersects = [(xi > xistol and xi < 1-xistol) for xi in xis]

            if sum(intersects) == 0:
              continue
            elif sum(intersects) == 1:
              xi = xis[intersects[0] == False]
            else:
              raise Exception('Found multiple ribbon intersections. MAXREFINE should be increased.')

      else:
        #TODO General order scheme based on bisection
        raise NotImplementedError('Simplex generation only implemented for order 1 and 2')

      assert ( xi > xistol and xi < 1.-xistol ), 'Illegal local coordinate'
 
      elem, trans = line.context

      pts = trans.apply( pts )

      newpoint = pts[0] + xi * ( pts[-1] - pts[0] )

      points = numeric.concatenate( [ points, newpoint[_] ], axis=0 ) 
      v1, v2 = line.vertices
      vertices.append( HalfVertex( v1, v2, xi=xi ) )

    try:
      submesh = util.delaunay( points )
    except RuntimeError:
      return [], []

    Simplex = TriangularElement if self.ndims == 2 else TetrahedronElement

    convex_hull = [ [vertices[iv] for iv in tri] for tri in submesh.convex_hull if numeric.greater_equal(tri,norig).all() ]

    ##########################################
    # Extract the simplices from the submesh #
    ##########################################

    simplices = []
    degensim  = []
    for tri in submesh.vertices:

      for j in range(2): #Flip two points in case of negative determinant
        if self.ndims == 3:
          # TODO renumber tetrahedron vertices to match triangle, removes the following:
          offset = points[ tri[0] ]
          affine = numeric.array( [ points[ tri[ii+1] ] - offset for ii in range(self.ndims) ] ).T
        else:
          offset = points[ tri[-1] ]
          affine = numeric.array( [ points[ tri[ii] ] - offset for ii in range(self.ndims) ] ).T

        trans = transform.Linear( affine ) + offset

        if trans.det > numeric.spacing(100):
          break

        tri[-2:] = tri[:-3:-1]

      else:
        if abs(trans.det) < numeric.spacing(100):
          degensim.append( [ vertices[ii] for ii in tri ] )
          continue

        raise Exception('Negative determinant with value %12.10e could not be resolved by flipping two vertices' % trans.det )

      simplices.append( Simplex( vertices=[ vertices[ii] for ii in tri ], parent=(self,trans) ) )

    assert len(simplices)+len(degensim)==submesh.vertices.shape[0], 'Simplices should be stored in either of the two containers'

    #############################################################
    # Loop over the edges of the simplex and check whether they #
    # reside in the part of the convex hull on the levelset     #
    #############################################################
      
    trimmededges = []  
              
    import itertools          
    for simplex in simplices:
      for iedge in range(self.ndims+1):

        #The edge potentially to be added to the trimmededges
        sedge = simplex.edge(iedge) 

        #Create lists to store edges which are to be checked on residence in the
        #convex hull, or which have been checked
        checkedges = [ sedge.vertices ]
        visitedges = []

        while checkedges:
          #Edge to be check on residence in the convex hull
          checkedge = checkedges.pop(0)
          visitedges.append( checkedge )

          #Check whether this edge is in the convex hull
          for hull_edge in convex_hull:
            #The checkedge is found in the convex hull. Append trimmededge and
            #terminate loop
            if all(checkvertex in hull_edge for checkvertex in checkedge):
              trimmededges.append( sedge )
              checkedges = []
              break
          else:
            #Check whether the checkedge is in a degenerate simplex
            for sim in degensim:
              if all(checkvertex in sim for checkvertex in checkedge):
                #Append all the edges to the checkedges pool
                for jedge in itertools.combinations(sim,self.ndims):
                  dedge = list(jedge)
                  for cedge in visitedges:
                    #The dedge is already in visitedges
                    if all(dvertex in cedge for dvertex in dedge):
                      break
                  else:
                    #The dedge is appended to to pool
                    checkedges.append( dedge )


    return simplices, trimmededges
  
class QuadElement( Element ):
  'quadrilateral element'

  __slots__ = ()

  def __init__( self, ndims, vertices, parent=None, context=None, interface=None ):
    'constructor'

    assert len(vertices) == 2**ndims
    Element.__init__( self, ndims, vertices, parent=parent, context=context, interface=interface )

  @property
  def neighbormap( self ):
    'maps # matching vertices --> codim of interface: {0: -1, 1: 2, 2: 1, 4: 0}'
    return dict( [ (0,-1) ] + [ (2**(self.ndims-i),i) for i in range(self.ndims+1) ] )

  def children_by( self, N ):
    'divide element by n'

    assert len(N) == self.ndims
    vertices = numeric.empty( [ ni+1 for ni in N ], dtype=object )
    vertices[ tuple( slice(None,None,ni) for ni in N ) ] = numeric.asarray( self.vertices ).reshape( [2]*self.ndims )
    for idim in range(self.ndims):
      s1 = tuple( slice(None) for ni in N[:idim] )
      s2 = tuple( slice(None,None,ni) for ni in N[idim+1:] )
      for i in range( 1, N[idim] ):
        vertices[s1+(i,)+s2] = numeric.objmap( HalfVertex, vertices[s1+(0,)+s2], vertices[s1+(2,)+s2], float(i)/N[idim] )

    elemvertices = [ vertices[ tuple( slice(i,i+2) for i in index ) ].ravel() for index in numeric.ndindex(*N) ]
    return tuple( QuadElement( vertices=elemvertices[ielem], ndims=self.ndims, parent=(self,trans) )
      for ielem, trans in enumerate( self.refinedtransform(N) ) )

  @property
  def children( self ):
    'all 1x refined elements'

    return self.children_by( (2,)*self.ndims )

  @staticmethod
  def edgetransform( ndims ):
    'edge transforms'

    transforms = []
    for idim in range( ndims ):
      offset = numeric.zeros( ndims )
      offset[:idim] = 1
      matrix = numeric.zeros(( ndims, ndims-1 ))
      matrix.flat[ :(ndims-1)*idim :ndims] = -1
      matrix.flat[ndims*(idim+1)-1::ndims] = 1
      # flip normal:
      offset[idim-1] = 1 - offset[idim-1]
      matrix[idim-1] *= -1
      transforms.append( transform.Linear(matrix.copy()) + offset )
      # other side:
      offset[idim] = 1
      # flip normal back:
      offset[idim-1] = 1 - offset[idim-1]
      matrix[idim-1] *= -1
      transforms.append( transform.Linear(matrix) + offset )
    return transforms

  @property
  def ribbons( self ):
    'ribbons'

    if self.ndims == 2:
      return [ self.edge(iedge) for iedge in range(4) ]

    if self.ndims != 3:
      raise NotImplementedError('Ribbons not implemented for ndims=%d'%self.ndims)

    ndvertices = numeric.asarray( self.vertices ).reshape( [2]*self.ndims )
    ribbons = []
    for i1, i2 in numeric.array([[[0,0,0],[1,0,0]],
                                 [[0,0,0],[0,1,0]],
                                 [[0,0,0],[0,0,1]],
                                 [[1,1,1],[0,1,1]],
                                 [[1,1,1],[1,0,1]],
                                 [[1,1,1],[1,1,0]],
                                 [[1,0,0],[1,1,0]],
                                 [[1,0,0],[1,0,1]],
                                 [[0,1,0],[1,1,0]],
                                 [[0,1,0],[0,1,1]],
                                 [[0,0,1],[1,0,1]],
                                 [[0,0,1],[0,1,1]]] ):
      trans = transform.Linear((i2-i1)[:,_]) + i1
      vertices = ndvertices[tuple(i1)], ndvertices[tuple(i2)]
      ribbons.append( QuadElement( vertices=vertices, ndims=1, context=(self,trans) ) )

    return ribbons

  @property
  def edges( self ):
    return [ self.edge(iedge) for iedge in range(2**self.ndims) ]

  def edge( self, iedge ):
    'edge'
    trans = self.edgetransform( self.ndims )[ iedge ]
    idim = iedge // 2
    iside = iedge % 2
    s = (slice(None,None, 1 if iside else -1),) * idim + (iside,) \
      + (slice(None,None,-1 if iside else  1),) * (self.ndims-idim-1)
    vertices = numeric.asarray( numeric.asarray( self.vertices ).reshape( (2,)*self.ndims )[s] ).ravel() # TODO check
    return QuadElement( vertices=vertices, ndims=self.ndims-1, context=(self,trans) )

  @staticmethod
  def refinedtransform( N ):
    'refined transform'

    Nrcp = numeric.reciprocal( N, dtype=float )
    scale = transform.Scale(Nrcp)
    return [ scale + i*Nrcp for i in numeric.ndindex(*N) ]

  def refine( self, n ):
    'refine n times'

    elems = [ self ]
    for i in range(n):
      elems = [ child for elem in elems for child in elem.children ]
    return elems

  @staticmethod
  def getgauss( degree ):
    'compute gauss points and weights'

    assert isinstance( degree, int ) and degree >= 0
    k = numeric.arange( 1, degree // 2 + 1 )
    d = k / numeric.sqrt( 4*k**2-1 )
    x, w = numeric.eigh( numeric.diagflat(d,-1) ) # eigh operates (by default) on lower triangle
    return (x+1) * .5, w[0]**2

  @classmethod
  def getischeme( cls, ndims, where ):
    'get integration scheme'

    x = w = None
    if ndims == 0:
      coords = numeric.zeros([1,0])
      weights = numeric.array([1.])
    elif where.startswith( 'gauss' ):
      N = eval( where[5:] )
      if isinstance( N, tuple ):
        assert len(N) == ndims
      else:
        N = [N]*ndims
      x, w = zip( *map( cls.getgauss, N ) )
    elif where.startswith( 'uniform' ):
      N = eval( where[7:] )
      if isinstance( N, tuple ):
        assert len(N) == ndims
      else:
        N = [N]*ndims
      x = [ numeric.arange( .5, n ) / n for n in N ]
      w = [ numeric.appendaxes( 1./n, n ) for n in N ]
    elif where.startswith( 'bezier' ):
      N = int( where[6:] )
      x = [ numeric.linspace( 0, 1, N ) ] * ndims
      w = [ numeric.appendaxes( 1./N, N ) ] * ndims
    elif where.startswith( 'subdivision' ):
      N = int( where[11:] ) + 1
      x = [ numeric.linspace( 0, 1, N ) ] * ndims
      w = None
    elif where.startswith( 'vtk' ):
      # TODO: if-block and subdiv fix unnecessary: coords = numpy.array( numpy.ndindex( ndims*(2,) ) )
      if ndims == 1:
        coords = numeric.array([[0,1]]).T
      elif ndims == 2:
        eps = 0 if not len(where[3:]) else float(where[3:]) # subdivision fix (avoid extraordinary point)
        coords = numeric.array([[eps,eps],[1-eps,eps],[1-eps,1-eps],[eps,1-eps]])
      elif ndims == 3:
        coords = numeric.array([ [0,0,0], [1,0,0], [0,1,0], [1,1,0], [0,0,1], [1,0,1], [0,1,1], [1,1,1] ])
      else:
        raise Exception, 'contour not supported for ndims=%d' % ndims
    elif where.startswith( 'contour' ):
      N = int( where[7:] )
      p = numeric.linspace( 0, 1, N )
      if ndims == 1:
        coords = p[_].T
      elif ndims == 2:
        coords = numeric.array([ p[ range(N) + [N-1]*(N-2) + range(N)[::-1] + [0]*(N-1) ],
                                 p[ [0]*(N-1) + range(N) + [N-1]*(N-2) + range(0,N)[::-1] ] ]).T
      elif ndims == 3:
        assert N == 0
        coords = numeric.array([ [0,0,0], [1,0,0], [0,1,0], [1,1,0], [0,0,1], [1,0,1], [0,1,1], [1,1,1] ])
      else:
        raise Exception, 'contour not supported for ndims=%d' % ndims
    else:
      raise Exception, 'invalid element evaluation %r' % where
    if x is not None:
      coords = numeric.empty( map( len, x ) + [ ndims ] )
      for i, xi in enumerate( x ):
        coords[...,i] = xi[ (slice(None),) + (_,)*(ndims-i-1) ]
      coords = coords.reshape( -1, ndims )
    if w is not None:
      weights = reduce( lambda weights, wi: ( weights * wi[:,_] ).ravel(), w )
    else:
      weights = None
    return coords, weights

  def select_contained( self, points, eps=0 ):
    'select points contained in element'

    selection = numeric.ones( points.shape[0], dtype=bool )
    for idim in range( self.ndims ):
      newsel = numeric.greater_equal( points[:,idim], -eps ) & numeric.less_equal( points[:,idim], 1+eps )
      selection[selection] &= newsel
      points = points[newsel]
      if not points.size:
        return None, None
    return points, selection

class TriangularElement( Element ):
  '''triangular element
     conventions: reference elem:   unit simplex {(x,y) | x>0, y>0, x+y<1}
                  vertex numbering: {(1,0):0, (0,1):1, (0,0):2}
                  edge numbering:   {bottom:0, slanted:1, left:2}
                  edge local coords run counter-clockwise.'''

  __slots__ = ()

  neighbormap = -1, 2, 1, 0
  edgetransform = (
    transform.Linear( numeric.array([[ 1],[ 0]]) ),
    transform.Linear( numeric.array([[-1],[ 1]]) ) + [1,0],
    transform.Linear( numeric.array([[ 0],[-1]]) ) + [0,1] )

  def __init__( self, vertices, parent=None, context=None ):
    'constructor'

    assert len(vertices) == 3
    Element.__init__( self, ndims=2, vertices=vertices, parent=parent, context=context )

  @property
  def children( self ):
    'all 1x refined elements'

    t1, t2, t3, t4 = self.refinedtransform( 2 )
    v1, v2, v3 = self.vertices
    h1, h2, h3 = HalfVertex(v1,v2), HalfVertex(v2,v3), HalfVertex(v3,v1)
    return tuple([ # TODO check!
      TriangularElement( vertices=[v1,h1,h3], parent=(self,t1) ),
      TriangularElement( vertices=[h1,v2,h2], parent=(self,t2) ),
      TriangularElement( vertices=[h3,h2,v3], parent=(self,t3) ),
      TriangularElement( vertices=[h2,h3,h1], parent=(self,t4) ) ])
      
  @property
  def edges( self ):
    return [ self.edge(iedge) for iedge in range(3) ]

  def edge( self, iedge ):
    'edge'

    trans = self.edgetransform[ iedge ]
    vertices = [ self.vertices[::-2], self.vertices[:2], self.vertices[1:] ][iedge]
    return QuadElement( vertices=vertices, ndims=1, context=(self,trans) )

  @staticmethod
  def refinedtransform( n ):
    'refined transform'

    transforms = []
    scale = transform.Scale( numeric.asarray([1./n,1./n]) )
    negscale = transform.Scale( numeric.asarray([-1./n,-1./n]) )
    for i in range( n ):
      transforms.extend( scale + numeric.array( [i,j], dtype=float ) / n for j in range(0,n-i) )
      transforms.extend( negscale + numeric.array( [n-j,n-i], dtype=float ) / n for j in range(n-i,n) )
    return transforms

  def refined( self, n ):
    'refine'

    assert n == 2
    if n == 1:
      return self
    return [ TriangularElement( id=self.id+'.child({})'.format(ichild), parent=(self,trans) ) for ichild, trans in enumerate( self.refinedtransform( n ) ) ]

  @staticmethod
  def getischeme( ndims, where ):
    '''get integration scheme
    gaussian quadrature: http://www.cs.rpi.edu/~flaherje/pdf/fea6.pdf
    '''

    assert ndims == 2
    if where.startswith( 'contour' ):
      n = int( where[7:] or 0 )
      p = numeric.arange( n+1, dtype=float ) / (n+1)
      z = numeric.zeros_like( p )
      coords = numeric.hstack(( [1-p,p], [z,1-p], [p,z] ))
      weights = None
    elif where.startswith( 'vtk' ):
      coords = numeric.array([[0,0],[1,0],[0,1]]).T
      weights = None
    elif where == 'gauss1':
      coords = numeric.array( [[1],[1]] ) / 3.
      weights = numeric.array( [1] ) / 2.
    elif where in 'gauss2':
      coords = numeric.array( [[4,1,1],[1,4,1]] ) / 6.
      weights = numeric.array( [1,1,1] ) / 6.
    elif where == 'gauss3':
      coords = numeric.array( [[5,9,3,3],[5,3,9,3]] ) / 15.
      weights = numeric.array( [-27,25,25,25] ) / 96.
    elif where == 'gauss4':
      A = 0.091576213509771; B = 0.445948490915965; W = 0.109951743655322
      coords = numeric.array( [[1-2*A,A,A,1-2*B,B,B],[A,1-2*A,A,B,1-2*B,B]] )
      weights = numeric.array( [W,W,W,1/3.-W,1/3.-W,1/3.-W] ) / 2.
    elif where == 'gauss5':
      A = 0.101286507323456; B = 0.470142064105115; V = 0.125939180544827; W = 0.132394152788506
      coords = numeric.array( [[1./3,1-2*A,A,A,1-2*B,B,B],[1./3,A,1-2*A,A,B,1-2*B,B]] )
      weights = numeric.array( [1-3*V-3*W,V,V,V,W,W,W] ) / 2.
    elif where == 'gauss6':
      A = 0.063089014491502; B = 0.249286745170910; C = 0.310352451033785; D = 0.053145049844816; V = 0.050844906370207; W = 0.116786275726379
      VW = 1/6. - (V+W) / 2.
      coords = numeric.array( [[1-2*A,A,A,1-2*B,B,B,1-C-D,1-C-D,C,C,D,D],[A,1-2*A,A,B,1-2*B,B,C,D,1-C-D,D,1-C-D,C]] )
      weights = numeric.array( [V,V,V,W,W,W,VW,VW,VW,VW,VW,VW] ) / 2.
    elif where == 'gauss7':
      A = 0.260345966079038; B = 0.065130102902216; C = 0.312865496004875; D = 0.048690315425316; U = 0.175615257433204; V = 0.053347235608839; W = 0.077113760890257
      coords = numeric.array( [[1./3,1-2*A,A,A,1-2*B,B,B,1-C-D,1-C-D,C,C,D,D],[1./3,A,1-2*A,A,B,1-2*B,B,C,D,1-C-D,D,1-C-D,C]] )
      weights = numeric.array( [1-3*U-3*V-6*W,U,U,U,V,V,V,W,W,W,W,W,W] ) / 2.
    elif where[:7] == 'uniform':
      N = int( where[7:] )
      points = ( numeric.arange( N ) + 1./3 ) / N
      NN = N**2
      C = numeric.empty( [2,N,N] )
      C[0] = points[:,_]
      C[1] = points[_,:]
      coords = C.reshape( 2, NN )
      flip = coords[0] + numeric.greater( coords[1], 1 )
      coords[:,flip] = 1 - coords[::-1,flip]
      weights = numeric.appendaxes( .5/NN, NN )
    elif where[:6] == 'bezier':
      N = int( where[6:] )
      points = numeric.linspace( 0, 1, N )
      coords = numeric.array([ [x,y] for i, y in enumerate(points) for x in points[:N-i] ]).T
      weights = None
    else:
      raise Exception, 'invalid element evaluation: %r' % where
    return coords.T, weights

  def select_contained( self, points, eps=0 ):
    'select points contained in element'

    selection = numeric.ones( points.shape[0], dtype=bool )
    for idim in 0, 1, 2:
      points_i = points[:,idim] if idim < 2 else 1-points.sum(1)
      newsel = numeric.greater_equal( points_i, -eps )
      selection[selection] &= newsel
      points = points[newsel]
      if not points.size:
        return None, None

    return points, selection

class TetrahedronElement( Element ):
  'tetrahedron element'

  __slots__ = ()

  neighbormap = -1, 3, 2, 1, 0
  #Defined to create outward pointing normal vectors for all edges (i.c. triangular faces)
  edgetransform = (
    transform.Linear( numeric.array([[ 0, 1],[1,0],[0,0]]) ),
    transform.Linear( numeric.array([[ 1, 0],[0,0],[0,1]]) ),
    transform.Linear( numeric.array([[ 0, 0],[0,1],[1,0]]) ),
    transform.Linear( numeric.array([[-1,-1],[1,0],[0,1]]) ) + [1,0,0] )

  def __init__( self, vertices, parent=None, context=None ):
    'constructor'

    assert len(vertices) == 4
    Element.__init__( self, ndims=3, vertices=vertices, parent=parent, context=context )

  @property
  def children( self ):
    'all 1x refined elements'
    raise NotImplementedError( 'Children of tetrahedron' )  
      
  @property
  def edges( self ):
    return [ self.edge(iedge) for iedge in range(4) ]

  def edge( self, iedge ):
    'edge'

    trans = self.edgetransform[ iedge ]
    v1, v2, v3, v4 = self.vertices
    vertices = [ [v1,v3,v2], [v1,v2,v4], [v1,v4,v3], [v2,v3,v4] ][ iedge ] # TODO check!
    return TriangularElement( vertices=vertices, context=(self,trans) )

  @staticmethod
  def refinedtransform( n ):
    'refined transform'
    raise NotImplementedError( 'Transformations for refined tetrahedrons' )  

  def refined( self, n ):
    'refine'
    raise NotImplementedError( 'Refinement tetrahedrons' )  

  @staticmethod
  def getischeme( ndims, where ):
    '''get integration scheme
       http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html'''

    assert ndims == 3
    if where.startswith( 'vtk' ):
      coords = numeric.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]).T
      weights = None
    elif where == 'gauss1':
      coords = numeric.array( [[1],[1],[1]] ) / 4.
      weights = numeric.array( [1] ) / 6.
    elif where == 'gauss2':
      coords = numeric.array([[0.5854101966249685,0.1381966011250105,0.1381966011250105],
                              [0.1381966011250105,0.1381966011250105,0.1381966011250105],
                              [0.1381966011250105,0.1381966011250105,0.5854101966249685],
                              [0.1381966011250105,0.5854101966249685,0.1381966011250105]]).T
      weights = numeric.array([1,1,1,1]) / 24.
    elif where == 'gauss3':
      coords = numeric.array([[0.2500000000000000,0.2500000000000000,0.2500000000000000],
                              [0.5000000000000000,0.1666666666666667,0.1666666666666667],
                              [0.1666666666666667,0.1666666666666667,0.1666666666666667],
                              [0.1666666666666667,0.1666666666666667,0.5000000000000000],
                              [0.1666666666666667,0.5000000000000000,0.1666666666666667]]).T
      weights = numeric.array([-0.8000000000000000,0.4500000000000000,0.4500000000000000,0.4500000000000000,0.4500000000000000]) / 6.
    elif where == 'gauss4':
      coords = numeric.array([[0.2500000000000000,0.2500000000000000,0.2500000000000000],
                              [0.7857142857142857,0.0714285714285714,0.0714285714285714],
                              [0.0714285714285714,0.0714285714285714,0.0714285714285714],
                              [0.0714285714285714,0.0714285714285714,0.7857142857142857],
                              [0.0714285714285714,0.7857142857142857,0.0714285714285714],
                              [0.1005964238332008,0.3994035761667992,0.3994035761667992],
                              [0.3994035761667992,0.1005964238332008,0.3994035761667992],
                              [0.3994035761667992,0.3994035761667992,0.1005964238332008],
                              [0.3994035761667992,0.1005964238332008,0.1005964238332008],
                              [0.1005964238332008,0.3994035761667992,0.1005964238332008],
                              [0.1005964238332008,0.1005964238332008,0.3994035761667992]]).T
      weights = numeric.array([-0.0789333333333333,0.0457333333333333,0.0457333333333333,0.0457333333333333,0.0457333333333333,0.1493333333333333,0.1493333333333333,0.1493333333333333,0.1493333333333333,0.1493333333333333,0.1493333333333333]) / 6.
    elif where == 'gauss5':
      coords = numeric.array([[0.2500000000000000,0.2500000000000000,0.2500000000000000],
                            [0.0000000000000000,0.3333333333333333,0.3333333333333333],
                            [0.3333333333333333,0.3333333333333333,0.3333333333333333],
                            [0.3333333333333333,0.3333333333333333,0.0000000000000000],
                            [0.3333333333333333,0.0000000000000000,0.3333333333333333],
                            [0.7272727272727273,0.0909090909090909,0.0909090909090909],
                            [0.0909090909090909,0.0909090909090909,0.0909090909090909],
                            [0.0909090909090909,0.0909090909090909,0.7272727272727273],
                            [0.0909090909090909,0.7272727272727273,0.0909090909090909],
                            [0.4334498464263357,0.0665501535736643,0.0665501535736643],
                            [0.0665501535736643,0.4334498464263357,0.0665501535736643],
                            [0.0665501535736643,0.0665501535736643,0.4334498464263357],
                            [0.0665501535736643,0.4334498464263357,0.4334498464263357],
                            [0.4334498464263357,0.0665501535736643,0.4334498464263357],
                            [0.4334498464263357,0.4334498464263357,0.0665501535736643]]).T
      weights = numeric.array([0.1817020685825351,0.0361607142857143,0.0361607142857143,0.0361607142857143,0.0361607142857143,0.0698714945161738,0.0698714945161738,0.0698714945161738,0.0698714945161738,0.0656948493683187,0.0656948493683187,0.0656948493683187,0.0656948493683187,0.0656948493683187,0.0656948493683187]) / 6.
    elif where == 'gauss6':
      coords = numeric.array([[0.3561913862225449,0.2146028712591517,0.2146028712591517],
                            [0.2146028712591517,0.2146028712591517,0.2146028712591517],
                            [0.2146028712591517,0.2146028712591517,0.3561913862225449],
                            [0.2146028712591517,0.3561913862225449,0.2146028712591517],
                            [0.8779781243961660,0.0406739585346113,0.0406739585346113],
                            [0.0406739585346113,0.0406739585346113,0.0406739585346113],
                            [0.0406739585346113,0.0406739585346113,0.8779781243961660],
                            [0.0406739585346113,0.8779781243961660,0.0406739585346113],
                            [0.0329863295731731,0.3223378901422757,0.3223378901422757],
                            [0.3223378901422757,0.3223378901422757,0.3223378901422757],
                            [0.3223378901422757,0.3223378901422757,0.0329863295731731],
                            [0.3223378901422757,0.0329863295731731,0.3223378901422757],
                            [0.2696723314583159,0.0636610018750175,0.0636610018750175],
                            [0.0636610018750175,0.2696723314583159,0.0636610018750175],
                            [0.0636610018750175,0.0636610018750175,0.2696723314583159],
                            [0.6030056647916491,0.0636610018750175,0.0636610018750175],
                            [0.0636610018750175,0.6030056647916491,0.0636610018750175],
                            [0.0636610018750175,0.0636610018750175,0.6030056647916491],
                            [0.0636610018750175,0.2696723314583159,0.6030056647916491],
                            [0.2696723314583159,0.6030056647916491,0.0636610018750175],
                            [0.6030056647916491,0.0636610018750175,0.2696723314583159],
                            [0.0636610018750175,0.6030056647916491,0.2696723314583159],
                            [0.2696723314583159,0.0636610018750175,0.6030056647916491],
                            [0.6030056647916491,0.2696723314583159,0.0636610018750175]]).T
      weights = numeric.array([0.0399227502581679,0.0399227502581679,0.0399227502581679,0.0399227502581679,0.0100772110553207,0.0100772110553207,0.0100772110553207,0.0100772110553207,0.0553571815436544,0.0553571815436544,0.0553571815436544,0.0553571815436544,0.0482142857142857,0.0482142857142857,0.0482142857142857,0.0482142857142857,0.0482142857142857,0.0482142857142857,0.0482142857142857,0.0482142857142857,0.0482142857142857,0.0482142857142857,0.0482142857142857,0.0482142857142857]) / 6.
    elif where == 'gauss7':
      coords = numeric.array([[0.2500000000000000,0.2500000000000000,0.2500000000000000],
                            [0.7653604230090441,0.0782131923303186,0.0782131923303186],
                            [0.0782131923303186,0.0782131923303186,0.0782131923303186],
                            [0.0782131923303186,0.0782131923303186,0.7653604230090441],
                            [0.0782131923303186,0.7653604230090441,0.0782131923303186],
                            [0.6344703500082868,0.1218432166639044,0.1218432166639044],
                            [0.1218432166639044,0.1218432166639044,0.1218432166639044],
                            [0.1218432166639044,0.1218432166639044,0.6344703500082868],
                            [0.1218432166639044,0.6344703500082868,0.1218432166639044],
                            [0.0023825066607383,0.3325391644464206,0.3325391644464206],
                            [0.3325391644464206,0.3325391644464206,0.3325391644464206],
                            [0.3325391644464206,0.3325391644464206,0.0023825066607383],
                            [0.3325391644464206,0.0023825066607383,0.3325391644464206],
                            [0.0000000000000000,0.5000000000000000,0.5000000000000000],
                            [0.5000000000000000,0.0000000000000000,0.5000000000000000],
                            [0.5000000000000000,0.5000000000000000,0.0000000000000000],
                            [0.5000000000000000,0.0000000000000000,0.0000000000000000],
                            [0.0000000000000000,0.5000000000000000,0.0000000000000000],
                            [0.0000000000000000,0.0000000000000000,0.5000000000000000],
                            [0.2000000000000000,0.1000000000000000,0.1000000000000000],
                            [0.1000000000000000,0.2000000000000000,0.1000000000000000],
                            [0.1000000000000000,0.1000000000000000,0.2000000000000000],
                            [0.6000000000000000,0.1000000000000000,0.1000000000000000],
                            [0.1000000000000000,0.6000000000000000,0.1000000000000000],
                            [0.1000000000000000,0.1000000000000000,0.6000000000000000],
                            [0.1000000000000000,0.2000000000000000,0.6000000000000000],
                            [0.2000000000000000,0.6000000000000000,0.1000000000000000],
                            [0.6000000000000000,0.1000000000000000,0.2000000000000000],
                            [0.1000000000000000,0.6000000000000000,0.2000000000000000],
                            [0.2000000000000000,0.1000000000000000,0.6000000000000000],
                            [0.6000000000000000,0.2000000000000000,0.1000000000000000]]).T
      weights = numeric.array([0.1095853407966528,0.0635996491464850,0.0635996491464850,0.0635996491464850,0.0635996491464850,-0.3751064406859797,-0.3751064406859797,-0.3751064406859797,-0.3751064406859797,0.0293485515784412,0.0293485515784412,0.0293485515784412,0.0293485515784412,0.0058201058201058,0.0058201058201058,0.0058201058201058,0.0058201058201058,0.0058201058201058,0.0058201058201058,0.1653439153439105,0.1653439153439105,0.1653439153439105,0.1653439153439105,0.1653439153439105,0.1653439153439105,0.1653439153439105,0.1653439153439105,0.1653439153439105,0.1653439153439105,0.1653439153439105,0.1653439153439105]) / 6.
    elif where == 'gauss8':
      coords = numeric.array([[0.2500000000000000,0.2500000000000000,0.2500000000000000],
                            [0.6175871903000830,0.1274709365666390,0.1274709365666390],
                            [0.1274709365666390,0.1274709365666390,0.1274709365666390],
                            [0.1274709365666390,0.1274709365666390,0.6175871903000830],
                            [0.1274709365666390,0.6175871903000830,0.1274709365666390],
                            [0.9037635088221031,0.0320788303926323,0.0320788303926323],
                            [0.0320788303926323,0.0320788303926323,0.0320788303926323],
                            [0.0320788303926323,0.0320788303926323,0.9037635088221031],
                            [0.0320788303926323,0.9037635088221031,0.0320788303926323],
                            [0.4502229043567190,0.0497770956432810,0.0497770956432810],
                            [0.0497770956432810,0.4502229043567190,0.0497770956432810],
                            [0.0497770956432810,0.0497770956432810,0.4502229043567190],
                            [0.0497770956432810,0.4502229043567190,0.4502229043567190],
                            [0.4502229043567190,0.0497770956432810,0.4502229043567190],
                            [0.4502229043567190,0.4502229043567190,0.0497770956432810],
                            [0.3162695526014501,0.1837304473985499,0.1837304473985499],
                            [0.1837304473985499,0.3162695526014501,0.1837304473985499],
                            [0.1837304473985499,0.1837304473985499,0.3162695526014501],
                            [0.1837304473985499,0.3162695526014501,0.3162695526014501],
                            [0.3162695526014501,0.1837304473985499,0.3162695526014501],
                            [0.3162695526014501,0.3162695526014501,0.1837304473985499],
                            [0.0229177878448171,0.2319010893971509,0.2319010893971509],
                            [0.2319010893971509,0.0229177878448171,0.2319010893971509],
                            [0.2319010893971509,0.2319010893971509,0.0229177878448171],
                            [0.5132800333608811,0.2319010893971509,0.2319010893971509],
                            [0.2319010893971509,0.5132800333608811,0.2319010893971509],
                            [0.2319010893971509,0.2319010893971509,0.5132800333608811],
                            [0.2319010893971509,0.0229177878448171,0.5132800333608811],
                            [0.0229177878448171,0.5132800333608811,0.2319010893971509],
                            [0.5132800333608811,0.2319010893971509,0.0229177878448171],
                            [0.2319010893971509,0.5132800333608811,0.0229177878448171],
                            [0.0229177878448171,0.2319010893971509,0.5132800333608811],
                            [0.5132800333608811,0.0229177878448171,0.2319010893971509],
                            [0.7303134278075384,0.0379700484718286,0.0379700484718286],
                            [0.0379700484718286,0.7303134278075384,0.0379700484718286],
                            [0.0379700484718286,0.0379700484718286,0.7303134278075384],
                            [0.1937464752488044,0.0379700484718286,0.0379700484718286],
                            [0.0379700484718286,0.1937464752488044,0.0379700484718286],
                            [0.0379700484718286,0.0379700484718286,0.1937464752488044],
                            [0.0379700484718286,0.7303134278075384,0.1937464752488044],
                            [0.7303134278075384,0.1937464752488044,0.0379700484718286],
                            [0.1937464752488044,0.0379700484718286,0.7303134278075384],
                            [0.0379700484718286,0.1937464752488044,0.7303134278075384],
                            [0.7303134278075384,0.0379700484718286,0.1937464752488044],
                            [0.1937464752488044,0.7303134278075384,0.0379700484718286]]).T
      weights = numeric.array([-0.2359620398477557,0.0244878963560562,0.0244878963560562,0.0244878963560562,0.0244878963560562,0.0039485206398261,0.0039485206398261,0.0039485206398261,0.0039485206398261,0.0263055529507371,0.0263055529507371,0.0263055529507371,0.0263055529507371,0.0263055529507371,0.0263055529507371,0.0829803830550589,0.0829803830550589,0.0829803830550589,0.0829803830550589,0.0829803830550589,0.0829803830550589,0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0134324384376852,0.0134324384376852,0.0134324384376852,0.0134324384376852,0.0134324384376852,0.0134324384376852,0.0134324384376852,0.0134324384376852,0.0134324384376852,0.0134324384376852,0.0134324384376852,0.0134324384376852]) / 6.
    else:
      raise Exception, 'invalid element evaluation: %r' % where
    return coords.T, weights

  def select_contained( self, points, eps=0 ):
    'select points contained in element'
    raise NotImplementedError( 'Determine whether a point resides in the tetrahedron' )  


## STDELEMS


class StdElem( cache.WeakCacheObject ):
  'stdelem base class'

  __slots__ = 'ndims', 'nshapes'

  def __init__( self, ndims, nshapes ):
    self.ndims = ndims
    self.nshapes = nshapes

  def __mul__( self, other ):
    'multiply elements'

    return PolyProduct( self, other )

  def __pow__( self, n ):
    'repeated multiplication'

    assert n >= 1
    return self if n == 1 else self * self**(n-1)

  def extract( self, extraction ):
    'apply extraction matrix'

    return ExtractionWrapper( self, extraction )

class PolyProduct( StdElem ):
  'multiply standard elements'

  __slots__ = 'std',

  def __init__( self, *std ):
    'constructor'

    std1, std2 = self.std = std
    StdElem.__init__( self, ndims=std1.ndims+std2.ndims, nshapes=std1.nshapes*std2.nshapes )

  def eval( self, points, grad=0 ):
    'evaluate'

    assert isinstance( grad, int ) and grad >= 0

    assert points.shape[-1] == self.ndims
    std1, std2 = self.std

    s1 = slice(0,std1.ndims)
    p1 = points[...,s1]
    s2 = slice(std1.ndims,None)
    p2 = points[...,s2]

    E = Ellipsis,
    S = slice(None),
    N = _,

    shape = points.shape[:-1] + (std1.nshapes * std2.nshapes,)
    G12 = [ ( std1.eval( p1, grad=i )[E+S+N+S*i+N*j]
            * std2.eval( p2, grad=j )[E+N+S+N*i+S*j] ).reshape( shape + (std1.ndims,) * i + (std2.ndims,) * j )
            for i,j in zip( range(grad,-1,-1), range(grad+1) ) ]

    data = numeric.empty( shape + (self.ndims,) * grad )

    s = (s1,)*grad + (s2,)*grad
    R = numeric.arange(grad)
    for n in range(2**grad):
      index = n>>R&1
      n = index.argsort() # index[s] = [0,...,1]
      shuffle = range(points.ndim) + list( points.ndim + n )
      iprod = index.sum()
      data.transpose(shuffle)[E+s[iprod:iprod+grad]] = G12[iprod]

    return data

  def __str__( self ):
    'string representation'

    return '%s*%s' % self.std

class PolyLine( StdElem ):
  'polynomial on a line'

  __slots__ = 'degree', 'poly'

  @classmethod
  def bernstein_poly( cls, degree ):
    'bernstein polynomial coefficients'

    # magic bernstein triangle
    revpoly = numeric.zeros( [degree+1,degree+1], dtype=int )
    for k in range(degree//2+1):
      revpoly[k,k] = root = (-1)**degree if k == 0 else ( revpoly[k-1,k] * (k*2-1-degree) ) / k
      for i in range(k+1,degree+1-k):
        revpoly[i,k] = revpoly[k,i] = root = ( root * (k+i-degree-1) ) / i
    return revpoly[::-1]

  @classmethod
  def spline_poly( cls, p, n ):
    'spline polynomial coefficients'

    assert p >= 0, 'invalid polynomial degree %d' % p
    if p == 0:
      assert n == -1
      return numeric.array( [[[1.]]] )

    assert 1 <= n < 2*p
    extractions = numeric.empty(( n, p+1, p+1 ))
    extractions[0] = numeric.eye( p+1 )
    for i in range( 1, n ):
      extractions[i] = numeric.eye( p+1 )
      for j in range( 2, p+1 ):
        for k in reversed( range( j, p+1 ) ):
          alpha = 1. / min( 2+k-j, n-i+1 )
          extractions[i-1,:,k] = alpha * extractions[i-1,:,k] + (1-alpha) * extractions[i-1,:,k-1]
        extractions[i,-j-1:-1,-j-1] = extractions[i-1,-j:,-1]

    poly = cls.bernstein_poly( p )
    return numeric.contract( extractions[:,_,:,:], poly[_,:,_,:], axis=-1 )

  @classmethod
  def spline_elems( cls, p, n ):
    'spline elements, minimum amount (just for caching)'

    return map( cls, cls.spline_poly(p,n) )

  @classmethod
  def spline_elems_neumann( cls, p, n ):
    'spline elements, neumann endings (just for caching)'

    polys = cls.spline_poly(p,n)
    poly_0 = polys[0].copy()
    poly_0[:,1] += poly_0[:,0]
    poly_e = polys[-1].copy()
    poly_e[:,-2] += poly_e[:,-1]
    return cls(poly_0), cls(poly_e)

  @classmethod
  def spline_elems_curvature( cls ):
    'spline elements, curve free endings (just for caching)'

    polys = cls.spline_poly(1,1)
    poly_0 = polys[0].copy()
    poly_0[:,0] += 0.5*(polys[0][:,0]+polys[0][:,1])
    poly_0[:,1] -= 0.5*(polys[0][:,0]+polys[0][:,1])

    poly_e = polys[-1].copy()
    poly_e[:,-2] -= 0.5*(polys[-1][:,-1]+polys[-1][:,-2])
    poly_e[:,-1] += 0.5*(polys[-1][:,-1]+polys[-1][:,-2])

    return cls(poly_0), cls(poly_e)

  @classmethod
  def spline( cls, degree, nelems, periodic=False, neumann=0, curvature=False ):
    'spline elements, any amount'

    p = degree
    n = 2*p-1
    if periodic:
      assert not neumann, 'periodic domains have no boundary'
      assert not curvature, 'curvature free option not possible for periodic domains'
      if nelems == 1: # periodicity on one element can only mean a constant
        elems = cls.spline_elems( 0, n )
      else:
        elems = cls.spline_elems( p, n )[p-1:p] * nelems
    else:
      elems = cls.spline_elems( p, min(nelems,n) )
      if len(elems) < nelems:
        elems = elems[:p-1] + elems[p-1:p] * (nelems-2*(p-1)) + elems[p:]
      if neumann:
        elem_0, elem_e = cls.spline_elems_neumann( p, min(nelems,n) )
        if neumann & 1:
          elems[0] = elem_0
        if neumann & 2:
          elems[-1] = elem_e
      if curvature:
        assert neumann==0, 'Curvature free not allowed in combindation with Neumann'
        assert degree==2, 'Curvature free only allowed for quadratic splines'
        elem_0, elem_e = cls.spline_elems_curvature()
        elems[0] = elem_0
        elems[-1] = elem_e

    return numeric.array( elems )

  def __init__( self, poly ):
    '''Create polynomial from order x nfuncs array of coefficients 'poly'.
       Evaluates to sum_i poly[i,:] x**i.'''

    self.poly = numeric.asarray( poly, dtype=float )
    order, nshapes = self.poly.shape
    self.degree = order - 1
    StdElem.__init__( self, ndims=1, nshapes=nshapes )

  def eval( self, points, grad=0 ):
    'evaluate'

    assert points.shape[-1] == 1
    x = points[...,0]

    if grad > self.degree:
      return numeric.appendaxes( 0., x.shape+(self.nshapes,)+(1,)*grad )

    poly = self.poly
    for n in range(grad):
      poly = poly[1:] * numeric.arange( 1, poly.shape[0] )[:,_]

    polyval = numeric.empty( x.shape+(self.nshapes,) )
    polyval[:] = poly[-1]
    for p in poly[-2::-1]:
      polyval *= x[...,_]
      polyval += p

    return polyval[(Ellipsis,)+(_,)*grad]

  def extract( self, extraction ):
    'apply extraction'

    return PolyLine( numeric.dot( self.poly, extraction ) )

  def __repr__( self ):
    'string representation'

    return 'PolyLine#%x' % id(self)

class PolyTriangle( StdElem ):
  '''poly triangle (linear for now)
     conventions: dof numbering as vertices, see TriangularElement docstring.'''

  __slots__ = ()

  def __init__( self, order ):
    'constructor'

    assert order == 1
    StdElem.__init__( self, ndims=2, nshapes=3 )

  def eval( self, points, grad=0 ):
    'eval'

    npoints, ndim = points.shape
    if grad == 0:
      x, y = points.T
      data = numeric.array( [ x, y, 1-x-y ] ).T
    elif grad == 1:
      data = numeric.array( [[1,0],[0,1],[-1,-1]], dtype=float )
    else:
      data = numeric.array( 0 ).reshape( (1,) * (grad+ndim) )
    return data

  def __repr__( self ):
    'string representation'

    return '%s#%x' % ( self.__class__.__name__, id(self) )

class BubbleTriangle( StdElem ):
  '''linear triangle + bubble function
     conventions: dof numbering as vertices (see TriangularElement docstring), then barycenter.'''

  __slots__ = ()

  def __init__( self, order ):
    'constructor'

    assert order == 1
    StdElem.__init__( self, ndims=2, nshapes=4 )

  def eval( self, points, grad=0 ):
    'eval'

    npoints, ndim = points.shape
    if grad == 0:
      x, y = points.T
      data = numeric.array( [ x, y, 1-x-y, 27*x*y*(1-x-y) ] ).T
    elif grad == 1:
      x, y = points.T
      const_block = numeric.array( [1,0,0,1,-1,-1]*npoints, dtype=float ).reshape( npoints,3,2 )
      grad1_bubble = 27*numeric.array( [y*(1-2*x-y),x*(1-x-2*y)] ).T.reshape( npoints,1,2 )
      data = numeric.concatenate( [const_block, grad1_bubble], axis=1 )
    elif grad == 2:
      x, y = points.T
      zero_block = numeric.zeros( (npoints,3,2,2) )
      grad2_bubble = 27*numeric.array( [-2*y,1-2*x-2*y, 1-2*x-2*y,-2*x] ).T.reshape( npoints,1,2,2 )
      data = numeric.concatenate( [zero_block, grad2_bubble], axis=1 )
    elif grad == 3:
      zero_block = numeric.zeros( (3,2,2,2) )
      grad3_bubble = 27*numeric.array( [0,-2,-2,-2,-2,-2,-2,0], dtype=float ).reshape( 1,2,2,2 )
      data = numeric.concatenate( [zero_block, grad3_bubble], axis=0 )
    else:
      assert ndim==2, 'Triangle takes 2D coordinates' # otherwise tested by unpacking points.T
      data = numeric.array( 0 ).reshape( (1,) * (grad+2) )
    return data

  def __repr__( self ):
    'string representation'

    return '%s#%x' % ( self.__class__.__name__, id(self) )

class ExtractionWrapper( StdElem ):
  'extraction wrapper'

  __slots__ = 'stdelem', 'extraction'

  def __init__( self, stdelem, extraction ):
    'constructor'

    self.stdelem = stdelem
    assert extraction.shape[0] == stdelem.nshapes
    self.extraction = extraction
    StdElem.__init__( self, ndims=stdelem.ndims, nshapes=extraction.shape[1] )

  def eval( self, points, grad=0 ):
    'call'

    return numeric.dot( self.stdelem.eval( points, grad ), self.extraction, axis=1 )

  def __repr__( self ):
    'string representation'

    return '%s#%x:%s' % ( self.__class__.__name__, id(self), self.stdelem )

class CatmullClarkElem( StdElem ):
  '''subdivision surface element
     implemented by Pieter Barendrecht/Timo van Opstal 2013.'''

  basepoly = PolyLine( numeric.array(
      [[ 1./6,  2./3,  1./6, 0.  ],
       [-1./2,  0.,    1./2, 0.  ],
       [ 1./2, -1.,    1./2, 0.  ],
       [-1./6,  1./2, -1./2, 1./6]] ) )**2

  @staticmethod
  def X( n, etype ):
    if etype == 1: # Interpolating corner
      Q0 = [0,0,1,2,0,0,1,2,3,3,4,5,6,6,7,8]
      X = numeric.zeros( (16, 9) )
      for R in range(16):
        X[ R,Q0[R] ] = 1
      # Phantom points (9)
      X[0,[0,1,3,4]] = 4,-2,-2,1
      X[1,[0,3]] = 2,-1
      X[2,[1,4]] = 2,-1
      X[3,[2,5]] = 2,-1
      X[4,[0,1]] = 2,-1
      X[8,[3,4]] = 2,-1
      X[12,[6,7]] = 2,-1

    elif etype == 3: # Regular boundary
      Q0 = [4,3,0,11,4,3,0,11,5,2,1,10,6,7,8,9]
      X = numeric.zeros( (16, 12) )
      for R in range(16):
        X[ R,Q0[R] ] = 1;
      # Phantom points (4)
      X[0,[4,5]] = 2,-1
      X[1,[3,2]] = 2,-1
      X[2,[0,1]] = 2,-1
      X[3,[11,10]] = 2,-1

    elif etype >= 4: # Interior (both regular and extraordinary)
      # These rows from Abar will be used
      Q00 = 1 if n==3 else 7
      Q0 = [Q00, 0, 3,     2*n+6,  6, 5, 4,     2*n+5,  2*n+4, 2*n+3, 2*n+2, 2*n+1,  2*n+12, 2*n+11, 2*n+10, 2*n+9 ]
      Q1 = [0,   3, 2*n+6, 2*n+15, 5, 4, 2*n+5, 2*n+14, 2*n+3, 2*n+2, 2*n+1, 2*n+13, 2*n+11, 2*n+10, 2*n+9,  2*n+8 ]
      Q2 = [1,   2, 2*n+7, 2*n+16, 0, 3, 2*n+6, 2*n+15, 5,     4,     2*n+5, 2*n+14, 2*n+3,  2*n+2,  2*n+1,  2*n+13]
      X = numeric.empty( [3, 16, 2*n+17] )
      X0 = numeric.zeros( (16, 2*n+17) )
      X1 = numeric.zeros( (16, 2*n+17) )
      X2 = numeric.zeros( (16, 2*n+17) )
      for R in range(16):
        X0[ R, Q0[R] ] = 1;
        X1[ R, Q1[R] ] = 1;
        X2[ R, Q2[R] ] = 1;
      X[0,:,:] = X0
      X[1,:,:] = X1
      X[2,:,:] = X2
 
    return X

  @staticmethod
  def Abar( n, etype ):
    'extended subdivision matrix'
    assert etype >= 4, 'Extended subdivision matrix should not be required.'

    # Matrix blocks
    i1, i2, i3 = numeric.cumsum( [ 2*n+1, 7, 9 ] )
    Abar = numeric.empty( [i3,i2] )
    Abar[0:i1,i1:i2] = 0
    S   = Abar[ 0:i1, 0:i1]
    S11 = Abar[i1:i2, 0:i1]
    S12 = Abar[i1:i2,i1:i2]
    S21 = Abar[i2:i3, 0:i1]
    S22 = Abar[i2:i3,i1:i2]

    # Extraordinary coefficients
    a = 1. - 7./(4.*n)
    b = 3./(2.*n**2)
    c = 1./(4.*n**2)
    d = 3./8.
    e = 1./16.
    f = 1./4.

    S[0,0] = a
    S[0,1:] = [b, c]*n
    S[1:,0] = [d, f]*n
    S[1,1:4] = d, e, e
    S[1,4:-2] = 0
    S[1,-2:] = e, e
    S[2,1:4] = f, f, f
    S[2,4:] = 0
    for i in range(3,2*n+1):
      for j in range(1,2*n+1):
        S[i,j] = S[i-2,-2] if j==1 else S[i-2,-1] if j==2 else S[i-2,j-2] # Circulant

    # Regular coefficients
    a = 9./16.
    b = 3./32.
    c = 1./64.

    if n == 3:
      S11[0] = c,0,0,b,a,b,0
      S11[1] = e,0,0,e,d,d,0
      S11[2] = b,c,0,c,b,a,b
      S11[3] = e,e,0,0,0,d,d
      S11[4] = e,0,0,d,d,e,0
      S11[5] = b,c,b,a,b,c,0
      S11[6] = e,e,d,d,0,0,0
    else:
      S11[0,:8] = c,0,0,b,a,b,0,0
      S11[1,:8] = e,0,0,e,d,d,0,0
      S11[2,:8] = b,0,0,c,b,a,b,c
      S11[3,:8] = e,0,0,0,0,d,d,e
      S11[4,:8] = e,0,0,d,d,e,0,0
      S11[5,:8] = b,c,b,a,b,c,0,0
      S11[6,:8] = e,e,d,d,0,0,0,0
      S11[:,8:] = 0

    S12[0] = c,b,c,0,b,c,0
    S12[1] = 0,e,e,0,0,0,0
    S12[2] = 0,c,b,c,0,0,0
    S12[3] = 0,0,e,e,0,0,0
    S12[4] = 0,0,0,0,e,e,0
    S12[5] = 0,0,0,0,c,b,c
    S12[6] = 0,0,0,0,0,e,e

    S21[0,:7] = 0,0,0,0,f,0,0
    S21[1,:7] = 0,0,0,0,d,e,0
    S21[2,:7] = 0,0,0,0,f,f,0
    S21[3,:7] = 0,0,0,0,e,d,e
    S21[4,:7] = 0,0,0,0,0,f,f
    S21[5,:7] = 0,0,0,e,d,0,0
    S21[6,:7] = 0,0,0,f,f,0,0
    S21[7,:7] = 0,0,e,d,e,0,0
    S21[8,:7] = 0,0,f,f,0,0,0
    if n != 3:
      S21[:,7:] = 0

    S22[0] = f,f,0,0,f,0,0
    S22[1] = e,d,e,0,e,0,0
    S22[2] = 0,f,f,0,0,0,0
    S22[3] = 0,e,d,e,0,0,0
    S22[4] = 0,0,f,f,0,0,0
    S22[5] = e,e,0,0,d,e,0
    S22[6] = 0,0,0,0,f,f,0
    S22[7] = 0,0,0,0,e,d,e
    S22[8] = 0,0,0,0,0,f,f
    return Abar

  @staticmethod
  def LimitStencil( n, etype ):
    'limit stencil'
    assert etype >= 4, 'Limit stencil should not be required.'

    # Multiplied by n with respect to the coefficients defined earlier
    b = 3./(2.*n)
    c = 1./(4.*n)
    # cf. Lai & Cheng paper
    Denom = 5. + 14.*b + 16.*c
    Alpha = 5. / Denom
    Beta = (12.*b + 8.*c) / (Denom*n)
    Gamma = (2.*b + 8.*c) / (Denom*n)
    return numeric.array( [Alpha] + [Beta, Gamma]*n + [0.]*7 )

  def __init__( self, n, etype ):
    'constructor'
    self.valence = n
    self.etype = etype

  def eval( self, points, grad=0 ):
    'evaluate'
    size = {1:9, 3:12}.get( self.etype, 2*self.valence+8 )
    result = numeric.empty( [numpy.size(points,0), size] + [2]*grad )
    X = self.X( self.valence, self.etype )

    if self.etype > 3:
      Abar = self.Abar( self.valence, self.etype )
      A = Abar[:2*self.valence+8]
      XAbar = numeric.dot( X, Abar )
      transf = lambda x, shift=0: (2*x-shift)[:,_]

      # Compute levels
      mvec = numpy.amax( numpy.asarray(points), axis=1 )
      where = mvec != 0. # relabel inf, would be lost by int conversion
      lvec = numpy.array( mvec.shape[0]*[-1] )
      lvec[where] = numpy.floor( -numpy.log2(mvec[where]) ).astype(int)
      assert all( lvec>-2 ), 'Point outside standard elem (i.e. outside [0,1]**d)'

      # Extraordinary point (l=-1)
      l_indices, = numpy.where( lvec==-1 )
      if len(l_indices):
        assert not grad, 'grad in extraordinary point unbounded.'
        result[l_indices,:] = self.LimitStencil( self.valence, self.etype )[_,:]

      # All other levels
      for level in range( max(lvec)+1 ):
        l_indices, = numpy.where( lvec==level )
        Al = numpy.dot( Al, A ) if level>0 else numeric.eye(len(A)) # At level 0, Al = Id
        if not len(l_indices): continue # No points to evaluate at this level
        coeff = 2**(grad*(level+1))

        ubar = 2**level*numpy.asarray(points)[l_indices,][:,0]
        vbar = 2**level*numpy.asarray(points)[l_indices,][:,1]
        # Evaluate points at 'level' for quadrant 'k' in [0,1,2]
        for k in range(3):
          k_indices, = numpy.where( vbar<.5 ) if k==0 else \
                       numpy.where( ubar<.5 ) if k==2 else \
                       numpy.where( numpy.logical_and( ubar>=.5, vbar>=.5 ) )
          if not len(k_indices): continue # No points to evaluate in this quadrant
          pts = numeric.concatenate( [transf( ubar[k_indices,], not k==2 ),
                                    transf( vbar[k_indices,], not k==0 )], axis=1 )
          temp0 = numeric.dot( XAbar[k], Al ).T
          temp1 = self.basepoly.eval( pts, grad=grad ).swapaxes(0,1)
          lk_indices = l_indices[k_indices,]
          result[lk_indices,] = coeff * numeric.dot( temp0, temp1 ).swapaxes(0,1)

    else: # etype in (1,3)
      # TODO: includes two DIRTY hacks verified only for etype=3
      newpoints = numeric.roll( points, 1, 1 ) if self.etype==3 and grad==0 else points
      temp = self.basepoly.eval( newpoints, grad=grad ).swapaxes(0,1)
      if self.etype==3 and grad==1: temp = temp.reshape(4,4,-1,2).swapaxes(0,1).reshape(16,-1,2)
      result = numeric.dot( X.T, temp ).swapaxes(0,1)

    return result

  #----------------------------------------------------------------------------------------------------
  @staticmethod
  def X_tmp( n, etype ):
    if etype == 1: # Interpolating corner
      Q0 = [0,0,1,2,0,0,1,2,3,3,4,5,6,6,7,8]
      X = numeric.zeros( (16, 9) )
      for R in range(16):
        X[ R,Q0[R] ] = 1
      # Phantom points (9)
      X[0,[0,1,3,4]] = 4,-2,-2,1
      X[1,[0,3]] = 2,-1
      X[2,[1,4]] = 2,-1
      X[3,[2,5]] = 2,-1
      X[4,[0,1]] = 2,-1
      X[8,[3,4]] = 2,-1
      X[12,[6,7]] = 2,-1

    elif etype == 3: # Regular boundary
      Q0 = [4,3,0,11,4,3,0,11,5,2,1,10,6,7,8,9]
      X = numeric.zeros( (16, 12) )
      for R in range(16):
        X[ R,Q0[R] ] = 1;
      # Phantom points (4)
      X[0,[4,5]] = 2,-1
      X[1,[3,2]] = 2,-1
      X[2,[0,1]] = 2,-1
      X[3,[11,10]] = 2,-1

    elif etype >= 4: # Interior (both regular and extraordinary)
      # These rows from Abar will be used
      Q00 = 1 if n == 3 else 7
      Q0 = [Q00, 6, 2*n + 4, 2*n + 12, 0, 5, 2*n + 3, 2*n + 11, 3, 4, 2*n + 2, 2*n + 10, 2*n + 6, 2*n + 5, 2*n + 1, 2*n + 9]
      Q1 = [0, 5, 2*n + 3, 2*n + 11, 3, 4, 2*n + 2, 2*n + 10, 2*n + 6, 2*n + 5, 2*n + 1, 2*n + 9, 2*n + 15, 2*n + 14, 2*n + 13, 2*n + 8]
      Q2 = [1, 0, 5, 2*n + 3, 2, 3, 4, 2*n + 2, 2*n + 7, 2*n + 6, 2*n + 5, 2*n + 1, 2*n + 16, 2*n + 15, 2*n + 14, 2*n + 13]
      X = numeric.empty( [3, 16, 2*n+17] )
      X0 = numeric.zeros( (16, 2*n+17) )
      X1 = numeric.zeros( (16, 2*n+17) )
      X2 = numeric.zeros( (16, 2*n+17) )
      for R in range(16):
        X0[ R, Q0[R] ] = 1;
        X1[ R, Q1[R] ] = 1;
        X2[ R, Q2[R] ] = 1;
      X[0,:,:] = X0
      X[1,:,:] = X1
      X[2,:,:] = X2
    #elif etype >= 4:
    #
    #  # Picking vectors
    #  # Dit zijn de rijen die uit Abar worden gebruikt.
    #  Q0 = [7, 6, 2*n + 4, 2*n + 12, 0, 5, 2*n + 3, 2*n + 11, 3, 4, 2*n + 2, 2*n + 10, 2*n + 6, 2*n + 5, 2*n + 1, 2*n + 9]
    #  if n == 3:
    #    Q0[0] = 1        
    #  Q1 = [0, 5, 2*n + 3, 2*n + 11, 3, 4, 2*n + 2, 2*n + 10, 2*n + 6, 2*n + 5, 2*n + 1, 2*n + 9, 2*n + 15, 2*n + 14, 2*n + 13, 2*n + 8]
    #  Q2 = [1, 0, 5, 2*n + 3, 2, 3, 4, 2*n + 2, 2*n + 7, 2*n + 6, 2*n + 5, 2*n + 1, 2*n + 16, 2*n + 15, 2*n + 14, 2*n + 13]
    #  
    #  X = numpy.empty( [3, 16, 2*n+17] )
    #  X0 = numpy.zeros( (16, 2*n+17) )
    #  X1 = numpy.zeros( (16, 2*n+17) )
    #  X2 = numpy.zeros( (16, 2*n+17) )
    #  
    #  for R in range(16):
    #    X0[ R, Q0[R] ] = 1;
    #    X1[ R, Q1[R] ] = 1;
    #    X2[ R, Q2[R] ] = 1;
    #  
    #  X[0,:,:] = X0
    #  X[1,:,:] = X1
    #  X[2,:,:] = X2

    #'''
    ## Interpolating corner
    #if etype == 1:
    #  Q0 = [0,0,1,2,0,0,1,2,3,3,4,5,6,6,7,8]
    #  #Q0 = [1,1,2,3,1,1,2,3,4,4,5,6,7,7,8,9]-1
    #  X = numpy.zeros( (16, 9) ) 
    #  
    #  for R in range(16):
    #    X[ R,Q0[R] ] = 1;
    #  
    #  # Phantom points (9)
    #  X[0,[0,1,3,4]] = 4,-2,-2,1
    #  X[1,[0,3]] = 2,-1
    #  X[2,[1,4]] = 2,-1
    #  X[3,[2,5]] = 2,-1
    #  X[4,[0,1]] = 2,-1
    #  X[8,[3,4]] = 2,-1
    #  X[12,[6,7]] = 2,-1
    #  
    ## Regular boundary
    #elif etype == 3:
    #  Q0 = [4,3,0,11,4,3,0,11,5,2,1,10,6,7,8,9]
    #  X = numpy.zeros( (16, 12) )
    #  
    #  for R in range(16):
    #    X[ R,Q0[R] ] = 1;
    #    
    #  # Phantom points (4)      
    #  X[0,[4,5]] = 2,-1
    #  X[1,[3,2]] = 2,-1
    #  X[2,[0,1]] = 2,-1
    #  X[3,[11,10]] = 2,-1
    #
    ## Interior (both regular and extraordinary)
    #elif etype >= 4:
    #
    #  # Picking vectors
    #  # Dit zijn de rijen die uit Abar worden gebruikt.
    #  Q0 = [7, 6, 2*n + 4, 2*n + 12, 0, 5, 2*n + 3, 2*n + 11, 3, 4, 2*n + 2, 2*n + 10, 2*n + 6, 2*n + 5, 2*n + 1, 2*n + 9]
    #  if n == 3:
    #    Q0[0] = 1        
    #  Q1 = [0, 5, 2*n + 3, 2*n + 11, 3, 4, 2*n + 2, 2*n + 10, 2*n + 6, 2*n + 5, 2*n + 1, 2*n + 9, 2*n + 15, 2*n + 14, 2*n + 13, 2*n + 8]
    #  Q2 = [1, 0, 5, 2*n + 3, 2, 3, 4, 2*n + 2, 2*n + 7, 2*n + 6, 2*n + 5, 2*n + 1, 2*n + 16, 2*n + 15, 2*n + 14, 2*n + 13]
    #  
    #  X = numpy.empty( [3, 16, 2*n+17] )
    #  X0 = numpy.zeros( (16, 2*n+17) )
    #  X1 = numpy.zeros( (16, 2*n+17) )
    #  X2 = numpy.zeros( (16, 2*n+17) )
    #  
    #  for R in range(16):
    #    X0[ R, Q0[R] ] = 1;
    #    X1[ R, Q1[R] ] = 1;
    #    X2[ R, Q2[R] ] = 1;
    #  
    #  X[0,:,:] = X0
    #  X[1,:,:] = X1
    #  X[2,:,:] = X2
    #'''
      
    return X

  def eval_bsplines( self, s, t, grad, direction ):
    BSplineS = numpy.array([-s**3 + 3*s**2 - 3*s + 1, 3*s**3 - 6*s**2 + 4, -3*s**3 + 3*s**2 + 3*s + 1, s**3]) / 6.0
    BSplineT = numpy.array([-t**3 + 3*t**2 - 3*t + 1, 3*t**3 - 6*t**2 + 4, -3*t**3 + 3*t**2 + 3*t + 1, t**3]) / 6.0
    
    dBSplineS = numpy.array([-3*s**2 + 6*s - 3, 9*s**2 - 12*s, -9*s**2 + 6*s + 3, 3*s**2]) / 6.0
    dBSplineT = numpy.array([-3*t**2 + 6*t - 3, 9*t**2 - 12*t, -9*t**2 + 6*t + 3, 3*t**2]) / 6.0
    
    ddBSplineS = numpy.array([ -6*s + 6, 18*s - 12, -18*s + 6, 6*s ]) / 6.0
    ddBSplineT = numpy.array([ -6*t + 6, 18*t - 12, -18*t + 6, 6*t ]) / 6.0
    
    if grad == 0:
      BSplineST = ( BSplineT[:,_] * BSplineS[_,:] ).reshape(16)
      
    elif grad == 1:
      if direction == 0:
        # S
        BSplineST = ( BSplineT[:,_] * dBSplineS[_,:] ).reshape(16)
      elif direction == 1:
        # T
        BSplineST = ( dBSplineT[:,_] * BSplineS[_,:] ).reshape(16)
      
    elif grad == 2:
      if direction == 0:
        # SS
        BSplineST = ( BSplineT[:,_] * ddBSplineS[_,:] ).reshape(16)
      elif direction == 1:
        # ST of TS
        BSplineST = ( dBSplineT[:,_] * dBSplineS[_,:] ).reshape(16)
      elif direction == 2:
        # TT
        BSplineST = ( ddBSplineT[:,_] * BSplineS[_,:] ).reshape(16)
    
    return BSplineST
    
  def eval_tmp( self, points, grad=0 ):
    X = self.X(self.valence,self.etype)
    if self.etype >= 4:
      result = numpy.empty( [ numpy.size(points,0), 2*self.valence+8] + [2]*grad )
      Abar = self.Abar(self.valence,self.etype)
      A = Abar[:2*self.valence+8]
    
      for p in range(numpy.size(points,0)):
    
        u = points[p,0]
        v = points[p,1] 
        
        m = max(u,v)
        
        if m == 0.:
          # Limitpoint!
          # S = self.Vinv[0].sum() is gelijk aan 1/Self.V[0,0]
          # Niet nodig als de eerste eigenvector gelijk is aan [1,1,...,1].T
          
          if grad == 0:
            # Voor het geval dat V[0,0] niet gelijk is aan 1:
            result[p,:] = self.LimitStencil(self.valence,self.etype) #self.Vinv[0] * self.V[0,0]          
            ''' Als grad > 0, is het resultaat 0. Dat komt omdat de afgeleiden van de B-splines een 
                partition of nullity vormen (is dat altijd zo met POU functies)? '''          
          else:
            # TODO Exception
            #print("Een grad > 0 in (0,0)")
            result[:] = numpy.nan
          
        else:
          # l = level
          l = int( numpy.floor(-numpy.log2(m))+1 )
          
          ubar = 2**(l-1)*u
          vbar = 2**(l-1)*v
          
          if vbar < .5:
            k = 0
            s = 2*ubar-1
            t = 2*vbar
          elif ubar > .5:
            k = 1
            s = 2*ubar-1
            t = 2*vbar-1
          else:
            k = 2
            s = 2*ubar
            t = 2*vbar-1
            
          # Efficienter implementeren  
          if grad == 0:
            result[p,:] = numpy.dot( numpy.linalg.matrix_power(A, (l-1)).T, numpy.dot( Abar.T, numpy.dot( X[k,:,:].T, self.eval_bsplines(s,t,grad,-1) ) ) )
          elif grad == 1:
            result[p,:,0] = 2**l * numpy.dot( numpy.linalg.matrix_power(A, (l-1)).T, numpy.dot( Abar.T, numpy.dot( X[k,:,:].T, self.eval_bsplines(s,t,grad,0) ) ) )
            result[p,:,1] = 2**l * numpy.dot( numpy.linalg.matrix_power(A, (l-1)).T, numpy.dot( Abar.T, numpy.dot( X[k,:,:].T, self.eval_bsplines(s,t,grad,1) ) ) )
          elif grad == 2:
            result[p,:,0,0] = 2**(2*l) * numpy.dot( numpy.linalg.matrix_power(A, (l-1)).T, numpy.dot( Abar.T, numpy.dot( X[k,:,:].T, self.eval_bsplines(s,t,grad,0) ) ) )
            result[p,:,0,1] = 2**(2*l) * numpy.dot( numpy.linalg.matrix_power(A, (l-1)).T, numpy.dot( Abar.T, numpy.dot( X[k,:,:].T, self.eval_bsplines(s,t,grad,1) ) ) )
            result[p,:,1,0] = 2**(2*l) * numpy.dot( numpy.linalg.matrix_power(A, (l-1)).T, numpy.dot( Abar.T, numpy.dot( X[k,:,:].T, self.eval_bsplines(s,t,grad,1) ) ) )
            result[p,:,1,1] = 2**(2*l) * numpy.dot( numpy.linalg.matrix_power(A, (l-1)).T, numpy.dot( Abar.T, numpy.dot( X[k,:,:].T, self.eval_bsplines(s,t,grad,2) ) ) )
    
    # etype = 1 of etype = 3
    else:
      if self.etype == 1:
        result = numpy.empty( [ numpy.size(points,0), 9] + [2]*grad )
      elif self.etype == 3:
        result = numpy.empty( [ numpy.size(points,0), 12] + [2]*grad )
    
      for p in range(numpy.size(points,0)):
    
        u = points[p,0]
        v = points[p,1]
        
        # Geen factoren voor de afgeleiden nodig! Oftewel, l=0 voor de hele patch.
    
        if grad == 0:
          # old = self.eval_bsplines(u,v,grad,-1)
          # new = self.basepoly.eval(numpy.array([v,u]),grad) # u,v swap
          # print 'grad=0, err: ', numpy.linalg.norm(old-new)
          result[p,:] = numpy.dot( X.T, self.eval_bsplines(u,v,grad,-1) )
        elif grad == 1:
          # err = []
          # for i in range(2):
          #   old = self.eval_bsplines(u,v,grad,i)
          #   new = self.basepoly.eval(numpy.array([u,v]),grad)[...,i].reshape(4,4).T.flatten() # sort
          #   err.append( numpy.linalg.norm(old-new) )
          # print 'grad=1, err: ', err
          result[p,:,0] = numpy.dot( X.T, self.eval_bsplines(u,v,grad,0) )
          result[p,:,1] = numpy.dot( X.T, self.eval_bsplines(u,v,grad,1) )
        elif grad == 2:
          result[p,:,0,0] = numpy.dot( X.T, self.eval_bsplines(u,v,grad,0) )
          result[p,:,0,1] = numpy.dot( X.T, self.eval_bsplines(u,v,grad,1) )
          result[p,:,1,0] = numpy.dot( X.T, self.eval_bsplines(u,v,grad,1) )
          result[p,:,1,1] = numpy.dot( X.T, self.eval_bsplines(u,v,grad,2) )
    temp_result = self.eval_tmp( points, grad )
    err = numpy.linalg.norm(temp_result-result)
    if err>1.e-10: print 'err: ', err, self.valence, self.etype, grad
                  
    return temp_result # numeric.array( result )

# vim:shiftwidth=2:foldmethod=indent:foldnestmax=2
