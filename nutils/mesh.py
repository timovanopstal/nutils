from . import topology, function, util, element, transform, numeric, _, log
import os, warnings

# MESH GENERATORS

class GridFunc( function.ElemFunc ):
  __slots__ = 'structure', 'grid'
  def __init__( self, domainelem, structure, grid ):
    self.structure = structure
    self.grid = grid
    function.ElemFunc.__init__( self, domainelem )
  def finditer( self, x ):
    assert x.ndim == 2
    assert x.shape[1] == len(self.grid)
    N = numeric.array([ numeric.searchsorted(gi,xi)-1 for gi, xi in zip(self.grid,x.T) ]).T
    I = numeric.arange( x.shape[0] )
    while N.size:
      n = N[0]
      GN = zip(self.grid,n)
      assert all( 0 <= ni < len(gi)-1 for gi, ni in GN )
      w = numeric.equal( N, n ).all( axis=1 )
      x0 = numeric.array([ gi[ni] for gi, ni in GN ])
      dx = numeric.array([ gi[ni+1]-gi[ni] for gi, ni in GN ])
      yield self.structure[tuple(n)], (x[w]-x0)/dx, I[w]
      N = N[~w]
      I = I[~w]
      x = x[~w]

def rectilinear( vertices, periodic=(), name='rect' ):
  'rectilinear mesh'

  ndims = len(vertices)
  indices = numeric.grid( len(n)-1 for n in vertices )
  domainelem = element.Element( ndims=ndims, vertices=[] )

  vertexfmt = name + '(' + ','.join( '%%%dd' % len(str(len(n)-1)) for n in vertices ) + ')'
  vertexobjs = numeric.objmap( lambda *index: element.PrimaryVertex(vertexfmt%index), *numeric.grid( len(n) for n in vertices ) )
  for idim in periodic:
    tmp = numeric.bringforward( vertexobjs, idim )
    tmp[-1] = tmp[0]

  structure = numeric.objmap( lambda *index: element.QuadElement(
    ndims=ndims,
    parent=( domainelem, transform.Scale( numeric.array([ n[i+1]-n[i] for n,i in zip(vertices,index) ]) ) + [ n[i] for n,i in zip(vertices,index) ] ),
    vertices=vertexobjs[tuple(slice(i,i+2) for i in index)].ravel() ), *indices )
  topo = topology.StructuredTopology( structure )
  coords = GridFunc( domainelem, structure, vertices )
  if periodic:
    topo = topo.make_periodic( periodic )
  return topo, coords

def revolve( topo, coords, nelems, degree=3, axis=0 ):
  'revolve coordinates'

  # This is a hack. We need to be able to properly multiply topologies.
  DEGREE = (2,) # Degree of coords element

  structure = numeric.array([ [ element.QuadElement( ndims=topo.ndims+1 ) for elem in topo ] for ielem in range(nelems) ])
  revolved_topo = topology.StructuredTopology( structure.reshape( nelems, *topo.structure.shape ), periodic=0 )
  if nelems % 2 == 0:
    revolved_topo.groups[ 'top' ] = revolved_topo[:nelems//2]
    revolved_topo.groups[ 'bottom' ] = revolved_topo[nelems//2:]

  print 'topo:', revolved_topo.structure.shape
  revolved_func = revolved_topo.splinefunc( degree=(degree,)+DEGREE )

  assert isinstance( coords, function.StaticDot )
  assert coords.array.ndim == 2
  nvertices, ndims = coords.array.shape

  phi = ( .5 + numeric.arange(nelems) - .5*degree ) * ( 2 * numeric.pi / nelems )
  weights = numeric.empty(( nelems, nvertices, ndims+1 ))
  weights[...,:axis] = coords.array[:,:axis]
  weights[...,axis] = numeric.cos(phi)[:,_] * coords.array[:,axis]
  weights[...,axis+1] = numeric.sin(phi)[:,_] * coords.array[:,axis]
  weights[...,axis+2:] = coords.array[:,axis+1:]
  weights = numeric.reshape( weights, 2, 1 )

  return revolved_topo, revolved_func.dot( weights )

def gmesh( path, btags={}, name=None ):
  'gmesh'

  if name is None:
    name = os.path.basename(path)

  if isinstance( btags, str ):
    btags = { i+1: btag for i, btag in enumerate( btags.split(',') ) }

  lines = iter( open( path, 'r' ) )

  assert lines.next() == '$MeshFormat\n'
  version, filetype, datasize = lines.next().split()
  assert lines.next() == '$EndMeshFormat\n'

  assert lines.next() == '$Nodes\n'
  nvertices = int( lines.next() )
  coords = numeric.empty(( nvertices, 3 ))
  for iNode in range( nvertices ):
    items = lines.next().split()
    assert int( items[0] ) == iNode + 1
    coords[ iNode ] = map( float, items[1:] )
  assert lines.next() == '$EndNodes\n'

  assert numeric.less( abs( coords[:,2] ), 1e-5 ).all(), 'ndims=3 case not yet implemented.'
  coords = coords[:,:2]

  boundary = []
  elements = []
  connected = [ set() for i in range( nvertices ) ]
  nmap = {}
  fmap = {}

  assert lines.next() == '$Elements\n'
  domainelem = element.Element( ndims=2, vertices=[] )
  vertexobjs = numeric.array( [ element.PrimaryVertex( '%s(%d)' % (name,ivertex) ) for ivertex in range(nvertices) ], dtype=object )
  for ielem in range( int( lines.next() ) ):
    items = lines.next().split()
    assert int( items[0] ) == ielem + 1
    elemtype = int( items[1] )
    ntags = int( items[2] )
    tags = [ int(tag) for tag in set( items[3:3+ntags] ) ]
    elemvertices = numeric.asarray( items[3+ntags:], dtype=int ) - 1
    elemvertexobjs = vertexobjs[ elemvertices ]
    elemcoords = coords[ elemvertices ]
    if elemtype == 1: # boundary edge
      boundary.append(( elemvertices, tags ))
    elif elemtype in (2,4):
      if elemtype == 2: # interior element, triangle
        parent = domainelem, transform.Linear( (elemcoords[:2]-elemcoords[2]).T ) + elemcoords[2]
        elem = element.TriangularElement( vertices=elemvertexobjs, parent=parent )
        stdelem = element.PolyTriangle( 1 )
      else: # interior element, quadrilateral
        raise NotImplementedError
        elem = element.QuadElement( ndims=2 )
        stdelem = element.PolyQuad( (2,2) )
      elements.append( elem )
      fmap[ elem ] = stdelem
      nmap[ elem ] = elemvertices
      for n in elemvertices:
        connected[ n ].add( elem )
    elif elemtype == 15: # boundary vertex
      pass
    else:
      raise Exception, 'element type #%d not supported' % elemtype
  assert lines.next() == '$EndElements\n'

  belements = []
  bgroups = {}
  for vertices, tags in boundary:
    n1, n2 = vertices
    elem, = connected[n1] & connected[n2]
    loc_vert_indices = [elem.vertices.index(vertexobjs[v]) for v in vertices] # in [0,1,2]
    match = numeric.array( loc_vert_indices ).sum()-1
    iedge = [1, 0, 2][match]
    belem = elem.edge( iedge )
    belements.append( belem )
    for tag in tags:
      bgroups.setdefault( tag, [] ).append( belem )

  structured = True
  for i, el in enumerate( belements ):
    if not set(belements[i-1].vertices) & set(el.vertices):
      structured = False
      warnings.warn( 'Boundary elements are not sorted: boundary group will be an UnstructuredTopology.' )
      break

  linearfunc = function.function( fmap, nmap, nvertices, 2 )
  # Extend linearfunc by bubble functions for the P^1+bubble basis
  fmap_b, nmap_b = {}, {}
  for i, (key,val) in enumerate( nmap.iteritems() ): # enumerate bubble functions
    fmap_b[key] = element.BubbleTriangle( 1 )
    nmap_b[key] = numeric.concatenate( [val, [nvertices+i]] )
  bubblefunc = function.function( fmap_b, nmap_b, nvertices+len(nmap), 2 )

  namedfuncs = { 'spline1': linearfunc, 'bubble1': bubblefunc }
  topo = topology.UnstructuredTopology( elements, ndims=2, namedfuncs=namedfuncs )
  topo.boundary = topology.StructuredTopology( belements, periodic=(0,) ) if structured else \
                  topology.UnstructuredTopology( belements, ndims=1 )
  topo.boundary.groups = {}
  for tag, group in bgroups.items():
    try:
      tag = btags[tag]
    except:
      pass
    topo.boundary.groups[tag] = topology.UnstructuredTopology( group, ndims=1 )

  return topo, function.ElemFunc( domainelem )

def blender( path ):
  '''Read an sdv file generated by the blender plugin for nutils subdivision surfaces.
     Implemented by Pieter Barendrecht August 2013.'''
  # Some definitions
  log.context( 'sdv' )
  log.info( 'importing: %s' % path )
  lines = iter( open( path, 'r' ) )
  # Non-spline element types
  corner_elem = element.CatmullClarkElem(4,etype=1)
  boundary_elem = element.CatmullClarkElem(4,etype=3)
  interior_elem = [element.CatmullClarkElem(valence,etype=4) for valence in range(3,9)]
  regular_elem = element.CatmullClarkElem(4,etype=5)

  # Read vertices
  items = lines.next().split()
  assert items[0] == 'Verts', 'Corrupted input file'
  numverts = int(items[1])
  dofs = numeric.zeros([numverts,3])
  for v in range(numverts):
    items = lines.next().split()
    dofs[v,0] = float(items[0])
    dofs[v,1] = float(items[1])
    dofs[v,2] = float(items[2])

  # Read rings
  items = lines.next().split()
  assert items[0] == 'Rings', 'Corrupted input file'
  numrings = int(items[1])
  neighbor_flag = bool(int(items[2]))

  # Read elements
  elements = []
  nmap = {}
  fmap = {}
  label = path.rsplit('/',1)[-1].rsplit('.')[0]
  localverts = {1:[2,3,5,6], 2:[], 3:[5,4,2,3], 4:[2,5,7,6], 5:[2,5,7,6]}
  for r in range(numrings):
    items = map( int, lines.next().split() )
    etype = items[0]
    valence = items[1]
    vertices = [ element.PrimaryVertex( label+'(%i)'%items[i] ) for i in localverts[etype] ]
    elem = element.QuadElement( 2, vertices )
    elements.append( elem )
    if neighbor_flag: lines.next() # neighbor information extracted from shared vertices

    # Local and global DOF numbering
    if etype == 1: # CornerBnd
      fmap[elem] = corner_elem
      warnings.warn( 'Corners are interpolary.' )
    elif etype == 2: # ExtBnd
      pass
    elif etype == 3: # RegBnd
      fmap[elem] = boundary_elem
    elif etype == 4: # ExtInt
      fmap[elem] = interior_elem[valence-3]      # Use Stam element
    elif etype == 5: # RegInt
      fmap[elem] = regular_elem
    nmap[elem] = numeric.array( items[2:] )

  # Read groups
  items = lines.next().split()
  assert items[0] == 'FaceGroups', 'Corrupted input file'
  numgroups = int(items[1]) # TODO: consider the case that there are no groups in the input file
  groups = {}
  for g in range(numgroups):
    items = lines.next().split()
    groupname = items[0]
    groupelem = [elements[int(item)] for item in items[1:]]
    groups[groupname] = topology.UnstructuredTopology( groupelem, ndims=2 )
    log.info( 'group %s contains %i elements' % (groupname, len(groupelem)) )

  # Create nutils objects
  cubicfunc = function.function( fmap, nmap, numverts, ndims=2 )
  namedfuncs = {'spline3': cubicfunc}
  topo = topology.UnstructuredTopology( elements, ndims=2, namedfuncs=namedfuncs )
  topo.groups = groups
  coords = (cubicfunc[:,_]*dofs).sum(0)
  return topo, coords

def triangulation( vertices, nvertices ):
  'triangulation'

  bedges = {}
  nmap = {}
  I = numeric.array( [[2,0],[0,1],[1,2]] )
  for n123 in vertices:
    elem = element.TriangularElement()
    nmap[ elem ] = n123
    for iedge, (n1,n2) in enumerate( n123[I] ):
      try:
        del bedges[ (n2,n1) ]
      except KeyError:
        bedges[ (n1,n2) ] = elem, iedge

  dofaxis = function.DofAxis( nvertices, nmap )
  stdelem = element.PolyTriangle( 1 )
  linearfunc = function.Function( dofaxis=dofaxis, stdmap=dict.fromkeys(nmap,stdelem) )
  namedfuncs = { 'spline1': linearfunc }

  connectivity = dict( bedges.iterkeys() )
  N = list( connectivity.popitem() )
  while connectivity:
    N.append( connectivity.pop( N[-1] ) )
  assert N[0] == N[-1]

  structure = []
  for n12 in zip( N[:-1], N[1:] ):
    elem, iedge = bedges[ n12 ]
    structure.append( elem.edge( iedge ) )
    
  topo = topology.UnstructuredTopology( list(nmap), ndims=2, namedfuncs=namedfuncs )
  topo.boundary = topology.StructuredTopology( structure, periodic=(1,) )
  return topo

def igatool( path, name=None ):
  'igatool mesh'

  if name is None:
    name = os.path.basename(path)

  import vtk

  reader = vtk.vtkXMLUnstructuredGridReader()
  reader.SetFileName( path )
  reader.Update()

  mesh = reader.GetOutput()

  FieldData = mesh.GetFieldData()
  CellData = mesh.GetCellData()

  NumberOfPoints = int( mesh.GetNumberOfPoints() )
  NumberOfElements = mesh.GetNumberOfCells()
  NumberOfArrays = FieldData.GetNumberOfArrays()

  points = util.arraymap( mesh.GetPoint, float, range(NumberOfPoints) )
  Cij = FieldData.GetArray( 'Cij' )
  Cv = FieldData.GetArray( 'Cv' )
  Cindi = CellData.GetArray( 'Elem_extr_indi')

  elements = []
  degree = 3
  ndims = 2
  nmap = {}
  fmap = {}

  poly = element.PolyLine( element.PolyLine.bernstein_poly( degree ) )**ndims

  for ielem in range(NumberOfElements):

    cellpoints = vtk.vtkIdList()
    mesh.GetCellPoints( ielem, cellpoints )
    nids = util.arraymap( cellpoints.GetId, int, range(cellpoints.GetNumberOfIds()) )

    assert mesh.GetCellType(ielem) == vtk.VTK_HIGHER_ORDER_QUAD
    nb = (degree+1)**2
    assert len(nids) == nb

    n = range( *util.arraymap( Cindi.GetComponent, int, ielem, [0,1] ) )
    I = util.arraymap( Cij.GetComponent, int, n, 0 )
    J = util.arraymap( Cij.GetComponent, int, n, 1 )
    Ce = numeric.zeros(( nb, nb ))
    Ce[I,J] = util.arraymap( Cv.GetComponent, float, n, 0 )

    vertices = [ element.PrimaryVertex( '%s(%d:%d)' % (name,ielem,ivertex) ) for ivertex in range(2**ndims) ]
    elem = element.QuadElement( vertices=vertices, ndims=ndims )
    elements.append( elem )

    fmap[ elem ] = element.ExtractionWrapper( poly, Ce.T )
    nmap[ elem ] = nids

  splinefunc = function.function( fmap, nmap, NumberOfPoints, ndims )
  namedfuncs = { 'spline%d' % degree: splinefunc }

  boundaries = {}
  elemgroups = {}
  vertexgroups = {}
  renumber   = (0,3,1,2)
  for iarray in range( NumberOfArrays ):
    name = FieldData.GetArrayName( iarray )
    index = name.find( '_group_' )
    if index == -1:
      continue
    grouptype = name[:index]
    groupname = name[index+7:]
    A = FieldData.GetArray( iarray )
    I = util.arraymap( A.GetComponent, int, range(A.GetSize()), 0 )
    if grouptype == 'edge':
      belements = [ elements[i//4].edge( renumber[i%4] ) for i in I ]
      boundaries[ groupname ] = topology.UnstructuredTopology( belements, ndims=ndims-1 )
    elif grouptype == 'vertex':
      vertexgroups[ groupname ] = I
    elif grouptype == 'element':
      elemgroups[ groupname ] = topology.UnstructuredTopology( [ elements[i] for i in I ], namedfuncs=namedfuncs, ndims=2 )
    else:
      raise Exception, 'unknown group type: %r' % grouptype

  topo = topology.UnstructuredTopology( elements, namedfuncs=namedfuncs, ndims=ndims )
  topo.groups = elemgroups
  topo.boundary = topology.UnstructuredTopology( elements=[], ndims=ndims-1 )
  topo.boundary.groups = boundaries

  for group in elemgroups.values():
    myboundaries = {}
    for name, boundary in boundaries.iteritems():
      belems = [ belem for belem in boundary.elements if belem.parent[0] in group ]
      if belems:
        myboundaries[ name ] = topology.UnstructuredTopology( belems, ndims=ndims-1 )
    group.boundary = topology.UnstructuredTopology( elements=[], ndims=ndims-1 )
    group.boundary.groups = myboundaries

  funcsp = topo.splinefunc( degree=degree )
  coords = ( funcsp[:,_] * points ).sum( 0 )
  return topo, coords #, vertexgroups

def fromfunc( func, nelems, ndims, degree=1 ):
  'piecewise'

  if isinstance( nelems, int ):
    nelems = [ nelems ]
  assert len( nelems ) == func.func_code.co_argcount
  topo, ref = rectilinear( [ numeric.linspace(0,1,n+1) for n in nelems ] )
  funcsp = topo.splinefunc( degree=degree ).vector( ndims )
  coords = topo.projection( func, onto=funcsp, coords=ref, exact_boundaries=True )
  return topo, coords

def demo( xmin=0, xmax=1, ymin=0, ymax=1 ):
  'demo triangulation of a rectangle'

  phi = numeric.arange( 1.5, 13 ) * (2*numeric.pi) / 12
  P = numeric.array([ numeric.cos(phi), numeric.sin(phi) ])
  P /= abs(P).max(axis=0)
  phi = numeric.arange( 1, 9 ) * (2*numeric.pi) / 8
  Q = numeric.array([ numeric.cos(phi), numeric.sin(phi) ])
  Q /= 2 * numeric.sqrt( abs(Q).max(axis=0) / numeric.sqrt(2) )
  R = numeric.zeros([2,1])

  coords = [.5*(xmin+xmax),.5*(ymin+ymax)] \
         + [.5*(xmax-xmin),.5*(ymax-ymin)] * numeric.hstack( [P,Q,R] ).T
  vertices = numeric.array(
    [ ( i, (i+1)%12, 12+(i-i//3)%8 )   for i in range(12) ]
  + [ ( 12+(i+1)%8, 12+i, i+1+(i//2) ) for i in range( 8) ]
  + [ ( 12+i, 12+(i+1)%8, 20 )         for i in range( 8) ] )
  
  domainelem = element.Element( ndims=2, vertices=[] )
  elements = []
  vertexobjs = numeric.array([ element.PrimaryVertex( 'demo.%d' % ivertex ) for ivertex in range(len(vertices)) ])
  for ielem, elemvertices in enumerate( vertices ):
    elemcoords = coords[ numeric.array(elemvertices) ]
    parent = domainelem, transform.Linear( (elemcoords[:2]-elemcoords[2]).T ) + elemcoords[2]
    elem = element.TriangularElement( vertices=vertexobjs[elemvertices], parent=parent )
    elements.append( elem )

  fmap = dict.fromkeys( elements, element.PolyTriangle(1) )
  nmap = dict( zip( elements, vertices ) )
  belems = [ elem.edge(1) for elem in elements[:12] ]
  bgroups = { 'top': belems[0:3], 'left': belems[3:6], 'bottom': belems[6:9], 'right': belems[9:12] }

  linearfunc = function.function( fmap, nmap, ndofs=21, ndims=2 )
  namedfuncs = { 'spline1': linearfunc }
  topo = topology.UnstructuredTopology( elements, ndims=2, namedfuncs=namedfuncs )
  topo.boundary = topology.UnstructuredTopology( belems, ndims=1 )
  topo.boundary.groups = dict( ( tag, topology.UnstructuredTopology( group, ndims=1 ) ) for tag, group in bgroups.items() )

  return topo, function.ElemFunc( domainelem )

#----------------------------------------------------------------------------------------------------
#def sdv_reader( path ):  
#  lines = iter( open( path, 'r' ) )
#  
#  # Read Verts  
#  items = lines.next().split()
#  numverts = int(items[1])
#  
#  #print "Number of vertices/nodes:", numverts
#  Verts = numeric.zeros([numverts,3])
#  
#  for v in range(numverts):
#    items = lines.next().split()
#    Verts[v,0] = float(items[0])
#    Verts[v,1] = float(items[1])
#    Verts[v,2] = float(items[2])
#  
#  # Read valencies (van de vertices of van de faces)
#  #items = lines.next().split()
#  #valencies = int(items[1:])
#  
#  # Read Rings
#  items = lines.next().split()
#  numrings = int(items[1])
#  
#  neighbour_flag = bool(int(items[2]))
#  #neighbour_flag = False
#  
#  #print "Number of faces/elements:", numrings
#  
#  elements = []
#  
#  # Valentie is natuurlijk 2 i.p.v. 4, later aanpassen.
#  corner_elem = CatmullClarkElem(4,etype=1)
#  # Valentie is natuurlijk 3 i.p.v. 4, later aanpassen.
#  boundary_elem = CatmullClarkElem(4,etype=3)
#  interior_elem = [ CatmullClarkElem(valence,etype=4) for valence in range(3,33) ]
#  regular_elem  = CatmullClarkElem(4,etype=5)  
#
#  nmap = {}
#  fmap = {} 
#  
#  # Cubic
#  p = 4
#  # Wat is n
#  n = 2*(p-1)-1
#  #poly = element.PolyLine.spline_poly( p, n )[p-2]
#  #stdelemcub = element.PolyLine( poly )**2
#  
#  for r in range(numrings):
#    items = map( int, lines.next().split() )
#    
#    nodes = [ element.PrimaryVertex( '%s(%d:%d)' % ('{}.quad({})'.format(path.rsplit('.',1)[0],r),r,inode) ) for inode in range(2**2) ]
#    elem = element.QuadElement( 2, nodes )
#    
#    #elem = element.QuadElement( id='{}.quad({})'.format(path.rsplit('.',1)[0],r),ndims=2 )
#    elements.append( elem )
#    
#    etype = items[0]
#    valence = items[1]
#    
#    #print ":: Loading, etype=", etype, "valence=", valence
#    
#    if etype == 1:
#      # 1 = CornerBnd
#      fmap[ elem ] = corner_elem
#    elif etype == 2:
#      # 2 = ExtBnd
#      pass
#    elif etype == 3:
#      # 3 = RegBnd
#      fmap[ elem ] = boundary_elem
#    elif etype == 4:
#      # 4 = ExtInt      
#      # Gebruik overal Stam elementen voor.
#      fmap[ elem ] = interior_elem[valence-3]      
#    elif etype == 5:
#      # 5 = RegInt
#      #fmap[ elem ] = stdelemcub
#      #nmap[ elem ] = numpy.array( [items[i+1] for i in (13,12,11,10,7,6,5,14,8,1,4,15,9,2,3,16)] )      
#      fmap[ elem ] = regular_elem  
#    
#    #if etype != 5:      
#    nmap[ elem ] = numeric.array( items[2:] )
#    
#    if neighbour_flag:
#      # Neighbours!
#      Neighbours = map( int, lines.next().split() )
#      #print "Neighbours of element", r, "(n=", valence, ") are", Neighbours
#      
#  # Het is niet zeker dat er groepen in het bestand aanwezig zijn. Check hier op!
#  # TODO      
#      
#  # Read Groups
#  items = lines.next().split()
#  numgroups = int(items[1])
#  
#  #print "Number of groups: ", numgroups
#  
#  groups = {}
#  for g in range(numgroups):
#    items = lines.next().split()    
#    groupname = items[0]
#    groupelem = [elements[int(item)] for item in items[1:]]
#    groups[groupname] = topology.UnstructuredTopology( groupelem, ndims=2 )
#    #print "Group ", groupname, " contains ", len(groupelem), " elements."
#    
#  cubicfunc = function.function( fmap, nmap, numverts, ndims=2 )
#  namedfuncs = { 'spline3': cubicfunc }
#  topo = topology.UnstructuredTopology( elements, ndims=2, namedfuncs=namedfuncs )
#  
#  topo.groups = groups
#  
#  coords = (cubicfunc*Verts.T).sum()
#  
#  # Verts toegevoegd
#  topo.verts = Verts
#  return topo, coords
#
#from numpy import *
#class CatmullClarkElem( element.StdElem ):
#  def __init__( self, n, etype ):
#    self.valence = n
#    self.etype = etype
#
#    # Cumulatief!
#    i1, i2, i3 = cumsum( [ 2*n+1, 7, 9 ] )
#    j1, j2 = cumsum( [2*n+1, 7] )    
#    
#    # Extended subdivision matrix
#    Abar = empty( [i3,j2] )
#    A = Abar[:i2]
#    
#    S = Abar[0:i1,0:j1]
#    Z = Abar[0:i1,j1:j2]
#    S11 = Abar[i1:i2,0:j1]
#    S12 = Abar[i1:i2,j1:j2]
#    S21 = Abar[i2:i3,0:j1]
#    S22 = Abar[i2:i3,j1:j2]
#    
#    # Extraordinary coefficients
#    a = 1. - 7./(4.*n)
#    b = 3./(2.*n**2)
#    c = 1./(4.*n**2)
#    d = 3./8.
#    e = 1./16.
#    f = 1./4.
#    
#    S[0,0] = a
#    # First row
#    S[0,1:] = [b, c]*n
#    # First column
#    S[1:,0] = [d, f]*n
#    
#    S[1,1:4] = d, e, e
#    S[1,4:-2] = 0
#    S[1,-2:] = e, e
#    S[2,1:4] = f, f, f
#    S[2,4:] = 0
#    
#    # Circulant
#    for i in range(3,2*n+1):
#      for j in range(1,2*n+1):
#        if j == 1:
#          S[i,j] = S[i-2,-2]
#        elif j == 2:
#          S[i,j] = S[i-2,-1]
#        else:
#          S[i,j] = S[i-2,j-2]
#    
#    Z[:] = 0
#    
#    # Regular coefficients
#    a = 9./16.
#    b = 3./32.
#    c = 1./64.
#    
#    if n == 3:
#      S11[0] = c,0,0,b,a,b,0
#      S11[1] = e,0,0,e,d,d,0
#      S11[2] = b,c,0,c,b,a,b
#      S11[3] = e,e,0,0,0,d,d
#      S11[4] = e,0,0,d,d,e,0
#      S11[5] = b,c,b,a,b,c,0
#      S11[6] = e,e,d,d,0,0,0      
#     
#    else:
#      # S11 = zeros( [7,2*n+1], dtype=float )
#      S11[0,:8] = c,0,0,b,a,b,0,0
#      S11[1,:8] = e,0,0,e,d,d,0,0
#      S11[2,:8] = b,0,0,c,b,a,b,c
#      S11[3,:8] = e,0,0,0,0,d,d,e
#      S11[4,:8] = e,0,0,d,d,e,0,0
#      S11[5,:8] = b,c,b,a,b,c,0,0
#      S11[6,:8] = e,e,d,d,0,0,0,0
#      
#      # Additional zeros      
#      S11[:,8:] = 0
#          
#    S12[0] = c,b,c,0,b,c,0
#    S12[1] = 0,e,e,0,0,0,0
#    S12[2] = 0,c,b,c,0,0,0
#    S12[3] = 0,0,e,e,0,0,0
#    S12[4] = 0,0,0,0,e,e,0
#    S12[5] = 0,0,0,0,c,b,c
#    S12[6] = 0,0,0,0,0,e,e
#
#    S21[0,:7] = 0,0,0,0,f,0,0
#    S21[1,:7] = 0,0,0,0,d,e,0
#    S21[2,:7] = 0,0,0,0,f,f,0
#    S21[3,:7] = 0,0,0,0,e,d,e
#    S21[4,:7] = 0,0,0,0,0,f,f
#    S21[5,:7] = 0,0,0,e,d,0,0
#    S21[6,:7] = 0,0,0,f,f,0,0
#    S21[7,:7] = 0,0,e,d,e,0,0
#    S21[8,:7] = 0,0,f,f,0,0,0
#    
#    # zeros
#    if n != 3:
#      S21[:,7:] = 0
#   
#    S22[0] = f,f,0,0,f,0,0
#    S22[1] = e,d,e,0,e,0,0
#    S22[2] = 0,f,f,0,0,0,0
#    S22[3] = 0,e,d,e,0,0,0
#    S22[4] = 0,0,f,f,0,0,0
#    S22[5] = e,e,0,0,d,e,0
#    S22[6] = 0,0,0,0,f,f,0
#    S22[7] = 0,0,0,0,e,d,e
#    S22[8] = 0,0,0,0,0,f,f
#    
#    # Interpolating corner
#    if etype == 1:
#      Q0 = [0,0,1,2,0,0,1,2,3,3,4,5,6,6,7,8]
#      #Q0 = [1,1,2,3,1,1,2,3,4,4,5,6,7,7,8,9]-1
#      X = zeros( (16, 9) ) 
#      
#      for R in range(16):
#        X[ R,Q0[R] ] = 1;
#      
#      # Phantom points (9)
#      X[0,[0,1,3,4]] = 4,-2,-2,1
#      X[1,[0,3]] = 2,-1
#      X[2,[1,4]] = 2,-1
#      X[3,[2,5]] = 2,-1
#      X[4,[0,1]] = 2,-1
#      X[8,[3,4]] = 2,-1
#      X[12,[6,7]] = 2,-1
#      
#    # Regular boundary
#    elif etype == 3:
#      Q0 = [4,3,0,11,4,3,0,11,5,2,1,10,6,7,8,9]
#      X = zeros( (16, 12) )
#      
#      for R in range(16):
#        X[ R,Q0[R] ] = 1;
#        
#      # Phantom points (4)      
#      X[0,[4,5]] = 2,-1
#      X[1,[3,2]] = 2,-1
#      X[2,[0,1]] = 2,-1
#      X[3,[11,10]] = 2,-1
#    
#    # Interior (both regular and extraordinary)
#    elif etype >= 4:
#    
#      self.A = A
#      self.Abar = Abar
#      
#      # Picking vectors
#      # Dit zijn de rijen die uit Abar worden gebruikt.
#      Q0 = [7, 6, 2*n + 4, 2*n + 12, 0, 5, 2*n + 3, 2*n + 11, 3, 4, 2*n + 2, 2*n + 10, 2*n + 6, 2*n + 5, 2*n + 1, 2*n + 9]
#      if n == 3:
#        Q0[0] = 1        
#      Q1 = [0, 5, 2*n + 3, 2*n + 11, 3, 4, 2*n + 2, 2*n + 10, 2*n + 6, 2*n + 5, 2*n + 1, 2*n + 9, 2*n + 15, 2*n + 14, 2*n + 13, 2*n + 8]
#      Q2 = [1, 0, 5, 2*n + 3, 2, 3, 4, 2*n + 2, 2*n + 7, 2*n + 6, 2*n + 5, 2*n + 1, 2*n + 16, 2*n + 15, 2*n + 14, 2*n + 13]
#      
#      X = empty( [3, 16, 2*n+17] )
#      X0 = zeros( (16, 2*n+17) )
#      X1 = zeros( (16, 2*n+17) )
#      X2 = zeros( (16, 2*n+17) )
#      
#      for R in range(16):
#        X0[ R, Q0[R] ] = 1;
#        X1[ R, Q1[R] ] = 1;
#        X2[ R, Q2[R] ] = 1;
#      
#      X[0,:,:] = X0
#      X[1,:,:] = X1
#      X[2,:,:] = X2
#      
#      # Multiplied by n with respect to the ones defined earlier
#      a = 1. - 7./(4.)
#      b = 3./(2.*n)
#      c = 1./(4.*n)
#      
#      # Paper van Lai en Cheng
#      Denom = 5. + 14.*b + 16.*c
#      Alpha = 5. / Denom
#      Beta = (12.*b + 8.*c) / (Denom*n)
#      Gamma = (2.*b + 8.*c) / (Denom*n)
#      
#      self.LimitStencil = [Alpha] + [Beta, Gamma]*n + [0.]*7
#      
#    #print ":: Initializing @ CatmullClarkElem for etype=", etype, " and valence=", n
#    self.X = X
#    print '@pieterinit: ', self.valence, self.etype, linalg.norm(self.X)
#    if etype>3: print '             ', linalg.norm(self.Abar)
#    
## ---------------------------------------------  
#
#  def eval_bsplines( self, s, t, grad, direction ):
#  
#    # Lambda-notatie gebruiken?
#    BSplineS = array([-s**3 + 3*s**2 - 3*s + 1, 3*s**3 - 6*s**2 + 4, -3*s**3 + 3*s**2 + 3*s + 1, s**3]) / 6.0
#    BSplineT = array([-t**3 + 3*t**2 - 3*t + 1, 3*t**3 - 6*t**2 + 4, -3*t**3 + 3*t**2 + 3*t + 1, t**3]) / 6.0
#    
#    dBSplineS = array([-3*s**2 + 6*s - 3, 9*s**2 - 12*s, -9*s**2 + 6*s + 3, 3*s**2]) / 6.0
#    dBSplineT = array([-3*t**2 + 6*t - 3, 9*t**2 - 12*t, -9*t**2 + 6*t + 3, 3*t**2]) / 6.0
#    
#    ddBSplineS = array([ -6*s + 6, 18*s - 12, -18*s + 6, 6*s ]) / 6.0
#    ddBSplineT = array([ -6*t + 6, 18*t - 12, -18*t + 6, 6*t ]) / 6.0
#    
#    if grad == 0:
#      BSplineST = ( BSplineT[:,_] * BSplineS[_,:] ).reshape(16)
#      
#    elif grad == 1:
#      if direction == 0:
#        # S
#        BSplineST = ( BSplineT[:,_] * dBSplineS[_,:] ).reshape(16)
#      elif direction == 1:
#        # T
#        BSplineST = ( dBSplineT[:,_] * BSplineS[_,:] ).reshape(16)
#      
#    elif grad == 2:
#      if direction == 0:
#        # SS
#        BSplineST = ( BSplineT[:,_] * ddBSplineS[_,:] ).reshape(16)
#      elif direction == 1:
#        # ST of TS
#        BSplineST = ( dBSplineT[:,_] * dBSplineS[_,:] ).reshape(16)
#      elif direction == 2:
#        # TT
#        BSplineST = ( ddBSplineT[:,_] * BSplineS[_,:] ).reshape(16)
#    
#    return BSplineST
#    
## --------------------------------------------- 
#
#  def eval( self, points, grad=0 ):
#    # Caching! Executed only once per valence.
#    
#    # Size van result hangt af van grad
#    # 1 alleen basis
#    # 2 basis naar x en y
#    # 4 basis naar xx, xy, yx en yy
#    
#    # grad = 0: oude situatie
#    # grad = 1: .. , .. , 2
#    # grad = 2:  .., .., 2, 2
#  
#    # u1 v1
#    # u2 v2
#    # .. ..
#    # un vn
#    
#    if self.etype >= 4:
#      result = empty( [ size(points,0), 2*self.valence+8] + [2]*grad )
#    
#      for p in range(size(points,0)):
#    
#        u = points[p,0]
#        v = points[p,1] 
#        
#        m = max(u,v)
#        
#        if m == 0.:
#          # Limitpoint!
#          # S = self.Vinv[0].sum() is gelijk aan 1/Self.V[0,0]
#          # Niet nodig als de eerste eigenvector gelijk is aan [1,1,...,1].T
#          
#          if grad == 0:
#            # Voor het geval dat V[0,0] niet gelijk is aan 1:
#            result[p,:] = self.LimitStencil #self.Vinv[0] * self.V[0,0]          
#            ''' Als grad > 0, is het resultaat 0. Dat komt omdat de afgeleiden van de B-splines een 
#                partition of nullity vormen (is dat altijd zo met POU functies)? '''          
#          else:
#            # TODO Exception
#            #print("Een grad > 0 in (0,0)")
#            result[:] = numpy.nan
#          
#        else:
#          # l = level
#          l = int( floor(-log2(m))+1 )
#          
#          ubar = 2**(l-1)*u
#          vbar = 2**(l-1)*v
#          
#          if vbar < .5:
#            k = 0
#            s = 2*ubar-1
#            t = 2*vbar
#          elif ubar > .5:
#            k = 1
#            s = 2*ubar-1
#            t = 2*vbar-1
#          else:
#            k = 2
#            s = 2*ubar
#            t = 2*vbar-1
#            
#          # Efficienter implementeren  
#          if grad == 0:
#            result[p,:] = dot( linalg.matrix_power(self.A, (l-1)).T, dot( self.Abar.T, dot( self.X[k,:,:].T, self.eval_bsplines(s,t,grad,-1) ) ) )
#          elif grad == 1:
#            result[p,:,0] = 2**l * dot( linalg.matrix_power(self.A, (l-1)).T, dot( self.Abar.T, dot( self.X[k,:,:].T, self.eval_bsplines(s,t,grad,0) ) ) )
#            result[p,:,1] = 2**l * dot( linalg.matrix_power(self.A, (l-1)).T, dot( self.Abar.T, dot( self.X[k,:,:].T, self.eval_bsplines(s,t,grad,1) ) ) )
#          elif grad == 2:
#            result[p,:,0,0] = 2**(2*l) * dot( linalg.matrix_power(self.A, (l-1)).T, dot( self.Abar.T, dot( self.X[k,:,:].T, self.eval_bsplines(s,t,grad,0) ) ) )
#            result[p,:,0,1] = 2**(2*l) * dot( linalg.matrix_power(self.A, (l-1)).T, dot( self.Abar.T, dot( self.X[k,:,:].T, self.eval_bsplines(s,t,grad,1) ) ) )
#            result[p,:,1,0] = 2**(2*l) * dot( linalg.matrix_power(self.A, (l-1)).T, dot( self.Abar.T, dot( self.X[k,:,:].T, self.eval_bsplines(s,t,grad,1) ) ) )
#            result[p,:,1,1] = 2**(2*l) * dot( linalg.matrix_power(self.A, (l-1)).T, dot( self.Abar.T, dot( self.X[k,:,:].T, self.eval_bsplines(s,t,grad,2) ) ) )
#    
#    # etype = 1 of etype = 3
#    else:
#    
#      if self.etype == 1:
#        result = empty( [ size(points,0), 9] + [2]*grad )
#      elif self.etype == 3:
#        result = empty( [ size(points,0), 12] + [2]*grad )
#    
#      for p in range(size(points,0)):
#    
#        u = points[p,0]
#        v = points[p,1]
#        
#        # Geen factoren voor de afgeleiden nodig! Oftewel, l=0 voor de hele patch.
#    
#        if grad == 0:
#          result[p,:] = dot( self.X.T, self.eval_bsplines(u,v,grad,-1) )
#        elif grad == 1:
#          result[p,:,0] = dot( self.X.T, self.eval_bsplines(u,v,grad,0) )
#          result[p,:,1] = dot( self.X.T, self.eval_bsplines(u,v,grad,1) )
#        elif grad == 2:
#          result[p,:,0,0] = dot( self.X.T, self.eval_bsplines(u,v,grad,0) )
#          result[p,:,0,1] = dot( self.X.T, self.eval_bsplines(u,v,grad,1) )
#          result[p,:,1,0] = dot( self.X.T, self.eval_bsplines(u,v,grad,1) )
#          result[p,:,1,1] = dot( self.X.T, self.eval_bsplines(u,v,grad,2) )
#                  
#    return numeric.array( result )

# vim:shiftwidth=2:foldmethod=indent:foldnestmax=1
