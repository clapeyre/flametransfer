#!/usr/bin/env python
# Author : A. Ruiz G. Staffebach
# Date : 01/02/2015

'''
Creates a XDMF_ file that describes the data stored in an AVBP HDF5 mesh and solution file.
The generated XDMF file can be used for visualizing the data with different 
visualization packages, e.g. Paraview, visit or ensight.

.. note:: This module can be directly executed. To see the available options type,

    >>> python xdmf.py -h

.. _XDMF: http://www.xdmf.org/index.php/Main_Page
.. _Paraview: http://www.paraview.org/
'''


# Define which datasets should be divided by rho
# All members of the RhoSpecies and FictiveSpecies will be divided by rho 
didiveByRhoVar = [ "rhoE","rhou","rhov","rhow","N2","O2","CO2","CO","CH4","KERO_LUCHE","H2O","O2","C3H8","OH","AIR","C8H18","C7H16"  ]


# ==========================================
# Define some generic string, useful later on
_header="""<?xml version="1.0" ?>
<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
<Domain>
"""

_header_bnd="""  <Grid GridType="Collection" CollectionType="Temporal" >
"""

_time = """    <Time Value="{0}"/>
"""


_footer = """</Domain>
</Xdmf>
"""

_footer_bnd = """  </Grid>
"""

# The topology ... in avbp, connectivity starts with node index 1, in 
# xdmf, should be 0 ... remove 1 from connectivity arrays

_Start_grid_collection = """  <Grid Name="{0}" GridType="Collection" CollectionType="Spatial" >
"""

_End_grid_collection = """  </Grid>
"""

_End_grid = """    </Grid>
"""

_Start_grid = """    <Grid Name="{0}" GridType="Uniform" >
"""

_Topology = """       <Topology TopologyType="{0}" Dimensions="{1}">
         <DataItem ItemType="Function" Dimensions="{2}" Function="$0 - 1">
           <DataItem Format="HDF" Dimensions="{2}" NumberType="Int" Precision="8" >
           {3}:{4}
           </DataItem>
         </DataItem>
       </Topology>
"""

_cutTopology = """       <Topology TopologyType="{0}" Dimensions="{1}">
         <DataItem Format="HDF" Dimensions="{2}" NumberType="Int" Precision="8" >
         {3}:{4}
         </DataItem>
       </Topology>
"""

_Topology_prism = """       <Topology TopologyType="{0}" Dimensions="{1}" Order="0 5 3 1 4 2">
         <DataItem ItemType="Function" Dimensions="{2}" Function="$0 - 1">
           <DataItem Format="HDF" Dimensions="{2}" NumberType="Int" Precision="8" >
           {3}:{4}
           </DataItem>
         </DataItem>
       </Topology>
"""

_Patch = """       <Grid Collection="{0}" Name="{1}">
         <Topology TopologyType="{2}" Dimensions="{3}">
              <DataItem ItemType="Function" Dimensions="{4}" Function="$0 - 1">
            <DataItem ItemType="Hyperslab" Dimensions="{4}" Type="HyperSlab" >
              <DataItem Dimensions="3"  Format="XML">
              {5} 1 {6}
              </DataItem>
                 <DataItem Format="HDF" Dimensions="{7}" NumberType="Int" Precision="8" >
                   {8}:{9}
                 </DataItem>
              </DataItem>
            </DataItem>
         </Topology>
         <Geometry Reference="/Xdmf/Domain/Grid/Grid/Grid/Geometry" />
"""

_geometry2D = """
       <Geometry GeometryType="X_Y">
         <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
         {0}:/Coordinates/x 
         </DataItem>
         <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
         {0}:/Coordinates/y
         </DataItem>
       </Geometry>
"""

_geometry3D = """
       <Geometry GeometryType="X_Y_Z">
         <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
         {0}:/Coordinates/x 
         </DataItem>
         <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
         {0}:/Coordinates/y
         </DataItem>
         <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
         {0}:/Coordinates/z
         </DataItem>
       </Geometry>
"""

_cutgeometry3D = """
       <Geometry GeometryType="X_Y_Z">
         <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
         {0}:/mesh/x 
         </DataItem>
         <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
         {0}:/mesh/y
         </DataItem>
         <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
         {0}:/mesh/z
         </DataItem>
       </Geometry>
"""

_bndgeometry2D = """
        <Geometry GeometryType="X_Y">
          <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
            {0}:{2}/Coordinates/x 
          </DataItem>
          <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
            {0}:{2}/Coordinates/y
          </DataItem>
        </Geometry>
"""

_bndgeometry3D = """
        <Geometry GeometryType="X_Y_Z">
          <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
            {0}:{2}/Coordinates/x 
          </DataItem>
          <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
            {0}:{2}/Coordinates/y
          </DataItem>
          <DataItem NumberType="Float" ItemType="Uniform" Dimensions="{1}" Precision="8" Format="HDF"> 
            {0}:{2}/Coordinates/z
          </DataItem>
        </Geometry>
"""

_attribute = """
       <Attribute Name="{0}" AttributeType="Scalar" Center="Node">
         <DataItem Format="HDF" Dimensions="{1}" NumberType="Float" Precision="8"> {2}:{3} </DataItem>
       </Attribute>
"""

_attribute_divideByRho = """
       <Attribute Name="{0}" AttributeType="Scalar" Center="Node">
         <DataItem ItemType="Function" Dimensions="{1}" Function="$0 / $1">
           <DataItem Format="HDF" Dimensions="{1}" NumberType="Float" Precision="8"> {2}:{3} </DataItem>
           <DataItem Format="HDF" Dimensions="{1}" NumberType="Float" Precision="8"> {2}:/{4}/rho </DataItem>
         </DataItem>
       </Attribute>        
"""

_attribute_3DvectordivideByRho = """
       <Attribute Name="Velocity" AttributeType="Vector" Center="Node">
         <DataItem ItemType="Function" Dimensions="{0} 3" Function="join($0 , $1 , $2 )">
           <DataItem ItemType="Function" Dimensions="{0}" Function="$0 / $1">
             <DataItem Format="HDF" Dimensions="{0}" NumberType="Float" Precision="8"> {1}:/GaseousPhase/rhou </DataItem>
             <DataItem Format="HDF" Dimensions="{0}" NumberType="Float" Precision="8"> {1}:/GaseousPhase/rho </DataItem>
           </DataItem>
           <DataItem ItemType="Function" Dimensions="{0}" Function="$0 / $1">
             <DataItem Format="HDF" Dimensions="{0}" NumberType="Float" Precision="8"> {1}:/GaseousPhase/rhov </DataItem>
             <DataItem Format="HDF" Dimensions="{0}" NumberType="Float" Precision="8"> {1}:/GaseousPhase/rho </DataItem>
           </DataItem>
           <DataItem ItemType="Function" Dimensions="{0}" Function="$0 / $1">
             <DataItem Format="HDF" Dimensions="{0}" NumberType="Float" Precision="8"> {1}:/GaseousPhase/rhow </DataItem>
             <DataItem Format="HDF" Dimensions="{0}" NumberType="Float" Precision="8"> {1}:/GaseousPhase/rho </DataItem>
           </DataItem>
         </DataItem>
       </Attribute>        
"""

# ==========================================


def outputAttributes_bnd(xdmf,solutbound,PatchPATH,nnodes):
    import h5py as h5
    import numpy as np
    import os.path

    solb  = h5.File(solutbound, 'r')

    # Define which datasets should be divided by rho
    # All members of the RhoSpecies will be divided by rho so don't put
    # them here
    groupPath = PatchPATH.strip()

    for datasetName in solb[groupPath].keys():
    # Put the dataset only if the datasetname is present in Listvar
    # if Listvar is an empty list ... put everything
      datasetPath = groupPath+"/"+datasetName
      # Only write the dataset if it has length nnodes
      # (skip parameters like dtsum)
      try:
        if solb[datasetPath].len()==nnodes:

          label = datasetName
          xdmf.write( _attribute.format(label, nnodes, solutbound, datasetPath) )

      except:
        pass # cannot use len() for a string dataset like versionstring ... 


def outputAttributes(xdmf,solName,Listvar,nnodes,ndim):
    import h5py as h5
    import numpy as np
    import os.path

    # Check consistency between solution and mesh[0]
    if Listvar == [ ]:

      sol  = h5.File(solName, 'r')

      groupNames = sol["/"].keys() # groups in file

      if ndim == 3 and "GaseousPhase" in groupNames: xdmf.write( _attribute_3DvectordivideByRho.format(nnodes, solName) )

      if "data" in groupNames:
        rhogroupname = 'data'
      else :
        rhogroupname = 'GaseousPhase'

      if "Parameters" in groupNames:
        dsets = sol["/Parameters"].keys()
        if "nnode" in dsets:
          sol_nnodes = sol["/Parameters/nnode"][0]
          if sol_nnodes!=nnodes:
            print "Incompatiblity : solution nnodes = {0}, mesh nnodes = {1}".format(sol_nnodes,nnodes)
            sys.exit()


# Treat groups 
      for groupName in groupNames:
        groupPath = "/"+groupName

        for datasetName in sol[groupPath].keys():
        # Put the dataset only if the datasetname is present in Listvar
        # if Listvar is an empty list ... put everything

          datasetPath = groupPath+"/"+datasetName
        # Only write the dataset if it has length nnodes
        # (skip parameters like dtsum)
          try:
            if sol[datasetPath].len()==nnodes:
#debug
#              print "storing {0}".format(datasetName)

            # Check if the datasetname needs to be divided by rho
              if datasetName in didiveByRhoVar or groupName == "RhoSpecies" or groupName == "FictiveSpecies":
                # rhou is changed into u
                label = datasetName.replace('rho', '')
                xdmf.write( _attribute_divideByRho.format(label, nnodes, solName, datasetPath,rhogroupname) )
              else:
                label = datasetName
                xdmf.write( _attribute.format(label, nnodes, solName, datasetPath) )

              Listvar.append(label)
          except:
            pass # cannot use len() for a string dataset like versionstring ...

    sol.close()

def toXDMF(xdmf, meshName , solName, varList ):
    """ 
        Creates XDMF_ file that describes data in HDF5 mesh and solution file.

        PARAMETERS:
            - xdmf:     xdmf file name
            - meshName: path to hdf5 mesh
            - solName:  path to hdf5 solution
            - varList:  if present, only variables in this list are included in XDMF file.

        .. note: XMDF file is created in the directory that contains HDF5 source file(s). 
    """ 
    import h5py as h5
    import numpy as np
    import os.path

    print " "

    mesh = h5.File(meshName, 'r')
    if  solName != None:
      print " -> Processing file "+ os.path.abspath(solName)
      sol  = h5.File(solName, 'r')
      # Read time value
      status="/Parameters/dtsum" in sol
      if not status : 
        print "    Warning dtsum not found in solution. Default used"
        dtsum = 0.0
      else: 
        dtsum=sol["/Parameters/dtsum"].value[0]
    else :
      print " -> Processing file "+ os.path.abspath(meshName)
      dtsum = 0.0

    # Coordinates
    nnodes = mesh["/Coordinates/x"].len()

# debug
#    print " The Mesh has {0} nodes".format(nnodes)

    totelements = 0

    ndim = 0
    for coord  in mesh["/Coordinates"].keys():
      ndim = ndim + 1

    for elementName in mesh["/Connectivity"].keys():
      totelements = totelements + 1

    if ndim == 2:
      geo=_geometry2D
    else:
      geo=_geometry3D

# debug
#    print " Found {0} element types".format(totelements)

    # Connectivity
    # Loop over element type
    Listvar = [ ]
    nelement = 0

    if solName != None:
      xdmf.write(_Start_grid_collection.format(solName))
    else:
      xdmf.write(_Start_grid_collection.format(meshName))

    xdmf.write(_time.format(dtsum))

    for elementName in mesh["/Connectivity"].keys():
      if '->node' in elementName:
        if elementName == 'tri->node':
          nvertex = 3
          cellType = 'Triangle'
          Toponame = 'Tri'
        elif elementName == 'qua->node':
          nvertex = 4
          cellType = 'Quadrilateral'
          Toponame = 'Quad'
        elif elementName == 'tet->node':
          nvertex = 4
          cellType = 'Tetrahedron'
          Toponame = 'Tetrahedra'
        elif elementName == 'pyr->node':
          nvertex = 5
          cellType = 'Pyramid'
          Toponame = 'Pyramid'
        elif elementName == 'hex->node':
          nvertex = 8
          cellType = 'Hexahedron'
          Toponame = 'Hexahedra'
        elif elementName == 'pri->node':
          nvertex = 6
          cellType = 'Wedge'
          Toponame = 'Wedge'
        else:
          print "Unknown element type :",elementName
          sys.exit()

        connectivityPath="/Connectivity/"+elementName
        ncells = mesh[connectivityPath].len() / nvertex

        xdmf.write( _Start_grid.format(Toponame) )

        if cellType == 'Wedge':
          xdmf.write( _Topology_prism.format(cellType, ncells, ncells*nvertex, meshName, connectivityPath) )
        else:
          xdmf.write( _Topology.format(cellType, ncells, ncells*nvertex, meshName, connectivityPath) )

        xdmf.write( geo.format(meshName, nnodes) )
        if solName != None:
          Listvar = [ ]
          outputAttributes(xdmf = xdmf,solName =solName, Listvar = Listvar ,nnodes = nnodes, ndim = ndim,)

        xdmf.write(_End_grid)

    xdmf.write(_End_grid_collection)

    mesh.close()

def toXDMFcut(xdmf, meshName , solName, varList ):
    """ 
        Creates XDMF_ file that describes data in HDF5 mesh and solution file.

        PARAMETERS:
            - xdmf:     xdmf file name
            - meshName: path to hdf5 mesh
            - solName:  path to hdf5 solution
            - varList:  if present, only variables in this list are included in XDMF file.

        .. note: XMDF file is created in the directory that contains HDF5 source file(s). 
    """ 
    import h5py as h5
    import numpy as np
    import os.path

    print " "
    print " -> Processing file "+ os.path.abspath(solName)

    mesh = h5.File(meshName, 'r')
    sol  = h5.File(solName, 'r')

    # Read time value
    dtsum=sol["/data/dtsum"].value[0]
    xdmf.write(_Start_grid_collection.format(solName))
    xdmf.write(_time.format(dtsum))

    # Coordinates
    nnodes = mesh["/mesh/x"].len()

    # Cuts are only available on 3d meshes
    ndim = 3
    geo=_cutgeometry3D
 
    # Looking for grid connectvity

    connec_store = []
    totelements = 0

    #only triangles for the cut 
    for dsetname in mesh["/mesh"].keys():
        if dsetname == 'triangles':
           totelements = totelements + 1
           connec_store.append(dsetname)

    if totelements==0:
        print 'Connectivity not found.'
        sys.exit()

    # Connectivity
    # Loop over element type: In this context only triangles are valid
    Listvar = [ ]
    nelement = 0
    output = 0
    while nelement < totelements : 
      elementName = connec_store[nelement]
      nelement +=1
      if elementName == 'triangles':
        nvertex = 3
        cellType = 'Triangle'
        Toponame = 'Tri'
      else:
        print "Unknown element type :",elementName
        sys.exit()

      connectivityPath="/mesh/"+elementName
      ncells = mesh[connectivityPath].len() / nvertex

      xdmf.write( _Start_grid.format(Toponame) )

      xdmf.write( _cutTopology.format(cellType, ncells, ncells*nvertex, meshName, connectivityPath) )

      xdmf.write( geo.format(meshName, nnodes) )

      Listvar = [ ]
      outputAttributes(xdmf = xdmf,solName =solName, Listvar = Listvar ,nnodes = nnodes, ndim = ndim)

      xdmf.write(_End_grid)

    xdmf.write(_End_grid_collection)

    sol.close()


def toXDMF_bnd(xdmf, meshName , solutbound ):
    """ 
        Creates XDMF_ file that describes data in HDF5 mesh and solutbound file.

        PARAMETERS:
            - xdmf:     xdmf file name
            - meshName: path to hdf5 mesh
            - solutbound:  path to hdf5 solutbound

        .. note: XMDF file is created in the directory that contains HDF5 source file(s). 
    """ 
    import h5py as h5
    import numpy as np
    import os.path

    print " "
    print " -> Processing file "+ os.path.abspath(solutbound)

    mesh = h5.File(meshName, 'r')

# Dummy time value

    Patches = mesh["/Patch"].keys()
    labels =  mesh["/Boundary/PatchLabels"]


    xdmf.write(_Start_grid_collection.format("Collection.h5"))

    for groupName in Patches:
       groupPath = "/Patch/"+groupName
       print "  >> Processing patch ", groupPath
       nnodes = mesh[groupPath+"/Coordinates/x"].len()
       if 'z' in mesh[groupPath+"/Coordinates"]:
         ndim = 3
         geo=_bndgeometry3D
       else:
         ndim = 2
         geo=_bndgeometry2D

       connec_store = []
       totelements = 0

       for dsetname in mesh[groupPath+"/Boundary"].keys():
         if dsetname == 'bnd_tri->node':
           totelements = totelements + 1
           connec_store.append(dsetname)      
         elif dsetname == 'bnd_qua->node':
           totelements = totelements + 1
           connec_store.append(dsetname)      
         elif dsetname == 'bnd_bi->node':
           totelements = totelements + 1
           connec_store.append(dsetname)      

# Connectivity
# Loop over element type: In this context only triangles are valid
       Listvar = [ ]
       nelement = 0
       while nelement < totelements : 
         elementName = connec_store[nelement]
         nelement +=1
         if elementName == 'bnd_tri->node':
           nvertex = 3
           cellType = 'Triangle'
           Toponame = 'Tri'
         elif elementName == 'bnd_qua->node':
           nvertex = 4
           cellType = 'Quadrilateral'
           Toponame = 'Quad'
         elif elementName == 'bnd_bi->node':
           nvertex = 2
           cellType = 'Polyline'
           Toponame = 'Bi'
         else:
           print "Unknown element type :",elementName
           sys.exit()

         connectivityPath=groupPath+"/Boundary/"+elementName
         ncells = mesh[connectivityPath].len() / nvertex

         if totelements > 1 : 
           xdmf.write(_Start_grid.format(labels[int(groupName)-1]+Toponame))
         else:
           xdmf.write(_Start_grid.format(labels[int(groupName)-1]))


         xdmf.write( _Topology.format(cellType, ncells, ncells*nvertex, meshName, connectivityPath) )

         xdmf.write( geo.format(meshName, nnodes, groupPath) )

         if int(groupName) < 10: add = '00'
         elif int(groupName) < 100: add = '0'

         PatchPATH="Patch_"+add+groupName+'-'+labels[int(groupName)-1]

         outputAttributes_bnd(xdmf = xdmf,solutbound = solutbound,  PatchPATH = PatchPATH ,nnodes = nnodes)


         xdmf.write(_End_grid)

    xdmf.write(_End_grid_collection)

    mesh.close()

# ===============================================================================
def main():
    from optparse import OptionParser
    import sys
    import os
    import glob
    import os.path

    print "\n Tool to generate visualization files \n"
    parser = OptionParser()
    parser.add_option('-m', '--mesh', dest = 'mesh', help = 'Path to hdf5 mesh file.'    , default = None)
    parser.add_option('-s', '--sol' , dest = 'sol' , help = 'Path to hdf5 solution (inst or average) file.', default = None)
    parser.add_option('-c', '--cut' , dest = 'cut' , help = 'Path to hdf5 cut file.', default = None)
    parser.add_option('-b', '--bnd' , dest = 'bnd' , help = 'Path to hdf5 solutbound file.', default = None)
    parser.add_option('-v', '--vars', dest = "vars", help = 'List of variables to export. For instance "rho H2O" ', default = None)

    soltype=False
    cutype=False

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(-1)
    else:
        (options, args) = parser.parse_args(sys.argv)

    if options.mesh == None:
        print """No mesh is specified"""
        sys.exit(-1)

    if options.sol == None:
        soltype = False
    else:
        soltype = True

    if options.cut == None:
        cutype = False
    else:
        cutype = True

    if options.bnd == None:
        bndtype = False
    else:
        bndtype = True

    # The varList, if present, has the following syntax
    # python xdmf.py -m mesh.mesh.h5 -s init.h5 -v "rho H2O rhoE"
    varList = None
    if options.vars != None:    
        options.vars = options.vars.replace('"','')
        varList = options.vars.split() # Split by space

    # If no meshName is specified, try to read the ../run.dat and use
    # this mesh
    meshName = None
    if options.mesh == None:    
      print "ERROR - cannot determine mesh file name, use --mesh or fill ../run.dat"
      parser.print_help()
      sys.exit(-1)
    else:
      meshName=options.mesh

    print 100*"-"
    print " "
    print " -> GRID : "+ os.path.abspath(meshName)
    print " "

    if cutype == False and soltype == False and bndtype == False and options.mesh != None:
      xdmfName = meshName.replace(".h5",".xmf")
      xdmf = open(xdmfName, "w")
      xdmf.write(_header)
      toXDMF(xdmf = xdmf, meshName = meshName, solName = None, varList = None)
      xdmf.write(_footer)
      xdmf.close()

      print "    - Xdmf file generated : ",xdmfName
      print " "

      print 100*"-"
      sys.exit(0)

    # Print argument -s to check wildcards are correctly taken into
    # account
    # Check if a wildcard is given in argument
    # if so, use glob to expand wildcard into file list
    if soltype :
      if contains(options.sol,"*"):
        if contains(options.sol,".h5"):
          solutionList = glob.glob(options.sol)
        else: 
          print 'You are using the * Wildcard. Please add .h5 to avoid conflits ex: -s "Mysol*.h5"'
          sys.exit(-1)
      else:
        # Create a List by splitting multiple solution by space
        solutionList = options.sol.split()

    if cutype :
      if contains(options.cut,"*"):
        if contains(options.cut,".h5"):
          solutionList = glob.glob(options.cut)
        else: 
          print 'You are using the * Wildcard. Please add .h5 to avoid conflits ex: -s "Mysol*.h5"'
          sys.exit(-1)
      else:
        # Create a List by splitting multiple solution by space
        solutionList = options.cut.split()

    if bndtype : 
      solutionList = options.bnd.split()

    # Sort solution list
    solutionList = sorted(solutionList)
    print " " 
    print " -> Files to be processed : \n"
    print "    ", solutionList
    print " " 

    # Loop over solution list
    # Note : it seems that the grid information is necessary for each time step in xdmf :(
    # It might more rapid if the mesh is read only once ... I prefer to
    # open it each time, as the xdmf data structure appears more clearly
    # this way
    for solName in solutionList:
      if not(os.path.isfile(solName)):
        print solName + " does not exist"
        sys.exit(-1)

      xdmfName = solName.replace(".h5",".xmf")
      xdmf = open(xdmfName, "w")

      xdmf.write(_header)
      if soltype: 
        toXDMF(xdmf = xdmf, meshName = meshName, solName = solName, varList = varList)
      if cutype: 
        toXDMFcut(xdmf = xdmf, meshName = meshName, solName = solName, varList = varList)
      if bndtype: 
        xdmf.write(_header_bnd)
        toXDMF_bnd(xdmf = xdmf, meshName = meshName, solutbound = solName)
        xdmf.write(_footer_bnd)

      xdmf.write(_footer)
      xdmf.close()
      print "    - Xdmf file generated : ",xdmfName
      print " "

    print 100*"-"


def contains(theString, theQueryValue):
  """ check if theString is contained into theQueryValue
x = "abcdefg"
print contains(x, "abc") // True
print contains(x, "xyz") //False """
  return theString.find(theQueryValue) > -1


# ==============================================================================
if __name__== "__main__":
    main()
    print "\n Warning : \n"
    print "  If variable is not present in the following list, it is not scaled by rho!!"
    print "  Scaled Variables are: \n"
    print  didiveByRhoVar
