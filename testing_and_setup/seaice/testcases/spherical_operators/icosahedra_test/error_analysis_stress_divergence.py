from netCDF4 import Dataset
import numpy as np

#-----------------------------------------------------------------------

def error_analysis_stress_divergence():

    # taylor series terms:
    # 0: f(x0,y0)
    # 1: df/dx(x0,y0) \Delta x
    # 2: df/dy(x0,y0) \Delta y
    # 3: 1/2 d2f/dx2(x0,y0) \Delta x^2
    # 4: 1/2 d2f/dy2(x0,y0) \Delta y^2
    # 5: d2f/dxdy(x0,y0) \Delta x \Delta y

    nTerms = 6
    dxPow = [0,1,0,2,1,0]
    dyPow = [0,0,1,0,1,2]
    terms = ["f","df/dx","df/dy","d2f/dx2","d2f/dxdy","d2f/dy2"]

    gridTypes = ["hex"]
    grids = {"hex" :"icos4"}

    for gridType in gridTypes:

        # grid file
        filegrid = Dataset("icos4/base_mesh.nc","r")

        nVertices = len(filegrid.dimensions["nVertices"])
        vertexDegree = len(filegrid.dimensions["vertexDegree"])

        nEdgesOnCell = filegrid.variables["nEdgesOnCell"][:]
        cellsOnVertex = filegrid.variables["cellsOnVertex"][:]
        cellsOnVertex[:] = cellsOnVertex[:] - 1
        cellsOnEdge = filegrid.variables["cellsOnEdge"][:]
        cellsOnEdge[:] = cellsOnEdge[:] - 1
        verticesOnCell = filegrid.variables["verticesOnCell"][:]
        verticesOnCell[:] = verticesOnCell[:] - 1
        verticesOnEdge = filegrid.variables["verticesOnEdge"][:]
        verticesOnEdge[:] = verticesOnEdge[:] - 1
        edgesOnVertex = filegrid.variables["edgesOnVertex"][:]
        edgesOnVertex[:] = edgesOnVertex[:] - 1
        edgesOnCell = filegrid.variables["edgesOnCell"][:]
        edgesOnCell[:] = edgesOnCell[:] - 1

        xVertex = filegrid.variables["xVertex"][:]
        yVertex = filegrid.variables["yVertex"][:]
        xCell = filegrid.variables["xCell"][:]
        yCell = filegrid.variables["yCell"][:]
        xEdge = filegrid.variables["xEdge"][:]
        yEdge = filegrid.variables["yEdge"][:]
        dcEdges = filegrid.variables["dcEdge"][:]
        dvEdge = filegrid.variables["dvEdge"][:]
        areaCell = filegrid.variables["areaCell"][:]
        areaTriangle = filegrid.variables["areaTriangle"][:]

        filegrid.close()

        with open('./namelist.seaice.stress_divergence') as namelistFile:
            namelistTxt = namelistFile.read().split('\n')
            for line in namelistTxt:
                if 'config_use_c_grid = ' in line:
                    config_use_c_grid = line.split()[-1]
                    break

        if (config_use_c_grid == 'true') :

           #iVertexTest = 1123 #for hex
           #print("  iVertexTest: ",iVertexTest)
           #for iCellOnVertex in range(0, vertexDegree):
           #    print("  iCell: ", iCellOnVertex, cellsOnVertex[iVertexTest,iCellOnVertex])

           #print("Edges involved:")
           #print(edgesOnVertex[iVertexTest,0])
           #print(edgesOnVertex[iVertexTest,1])
           #print(edgesOnVertex[iVertexTest,2])

           iEdgeTest = 3016 #2893, 2997, 3016
           print("  iEdgeTest: ",iEdgeTest)
           for iCellOnEdge in range(0, 2):
              print("  iCell: ", iCellOnEdge, cellsOnEdge[iEdgeTest,iCellOnEdge])

           #print("xVertex[iVertexTest]=",xVertex[iVertexTest],"yVertex[iVertexTest]=",yVertex[iVertexTest])
           print("xEdge[iEdgeTest]=",xEdge[iEdgeTest],"yEdge[iEdgeTest]=",yEdge[iEdgeTest])

           basises = ["wachspress","pwl"]
           for basis in basises:

              print("   ",basis)

              # variational
              filein = Dataset("./output_%s_icos4/output.2000.nc" % (basis),"r")

              basisIntegralsU = filein.variables["basisIntegralsU"][:]
              basisIntegralsV = filein.variables["basisIntegralsV"][:]
              basisIntegralsMetric = filein.variables["basisIntegralsMetric"][:]
              basisIntegralsUTri = filein.variables["basisIntegralsUTri"][:]
              basisIntegralsVTri = filein.variables["basisIntegralsVTri"][:]
              basisIntegralsMetricTri = filein.variables["basisIntegralsMetricTri"][:]
              cellEdgesAtEdge = filein.variables["cellEdgesAtEdge"][:]
              cellEdgesAtEdge[:] = cellEdgesAtEdge[:] - 1
              triangleEdgesAtEdge = filein.variables["triangleEdgesAtEdge"][:]
              triangleEdgesAtEdge[:] = triangleEdgesAtEdge[:] - 1
              varDenCGrid = filein.variables["variationalDenominatorCGrid"][:]
              areaCellCGrid = filein.variables["areaCellCGrid"][:]
              areaTriCGrid = filein.variables["areaTriCGrid"][:]

              filein.close()

              # gradient taylor terms
              dfdx = np.zeros(nTerms) # ideal: (0,1,0,0,0,0)
              dfdy = np.zeros(nTerms) # ideal: (0,0,1,0,0,0)

              for iTerm in range(0,nTerms):

                 # loop over surrounding cells
                 for iSurroundingCell in range(0,2):
                    termCheck = 0

                    # get the cell number of this cell
                    iCell = cellsOnEdge[iEdgeTest,iSurroundingCell]

                    # get the local index of iEdge on iSurroundingCell
                    iVelocityEdge = cellEdgesAtEdge[iEdgeTest,iSurroundingCell]

                    # loop over the vertices of the surrounding cell
                    for iStressEdge in range(0,nEdgesOnCell[iCell]):

                        iEdge = edgesOnCell[iCell,iStressEdge]

                        dx = xEdge[iEdge] - xEdge[iEdgeTest]
                        dy = yEdge[iEdge] - yEdge[iEdgeTest]

                        fTerm = pow(dx,dxPow[iTerm]) * pow(dy,dyPow[iTerm])

                        dfdx[iTerm] -= (fTerm * basisIntegralsU[iCell, iVelocityEdge, iStressEdge]) #/ varDenCGrid[iEdgeTest]
                        dfdy[iTerm] -= (fTerm * basisIntegralsV[iCell, iVelocityEdge, iStressEdge]) #/ varDenCGrid[iEdgeTest]

                    #CHECK
                    areaCheck = 0.0
                    for i in range(0,nEdgesOnCell[iCell]):
                       for j in range(0,nEdgesOnCell[iCell]):
                          areaCheck += basisIntegralsMetric[iCell,i,j]
                    if iTerm == 0:
                        print("cell=", iCell, "areaCell=", areaCell[iCell], "areaCheck=", areaCheck, "areaCellCGrid=", areaCellCGrid[iCell], "diff=", abs(areaCheck-areaCellCGrid[iCell]))
          

                 # loop over surrounding triangles
                 for iSurroundingTri in range(0,2):

                    termCheck = 0
                    # get the vertex number of this triangle
                    iVertex = verticesOnEdge[iEdgeTest,iSurroundingTri]

                    # get the local index of iEdge on iSurroundingTri
                    iVelocityEdge = triangleEdgesAtEdge[iEdgeTest,iSurroundingTri]

                    # loop over the edges of the surrounding triangle
                    for iStressEdge in range(0,vertexDegree):

                       iEdge = edgesOnVertex[iVertex,iStressEdge]

                       dx = xEdge[iEdge] - xEdge[iEdgeTest]
                       dy = yEdge[iEdge] - yEdge[iEdgeTest]

                       fTerm = pow(dx,dxPow[iTerm]) * pow(dy,dyPow[iTerm])

                       dfdx[iTerm] -= (fTerm * basisIntegralsUTri[iVertex, iVelocityEdge, iStressEdge]) #/ varDenCGrid[iEdgeTest]
                       dfdy[iTerm] -= (fTerm * basisIntegralsVTri[iVertex, iVelocityEdge, iStressEdge]) #/ varDenCGrid[iEdgeTest]

                    #CHECK
                    areaCheck = 0.0
                    for i in range(0,vertexDegree):
                       for j in range(0,vertexDegree):
                         areaCheck += basisIntegralsMetricTri[iVertex,i,j]
                    if iTerm == 0:
                       print("triangle=", iVertex, "areaCheck=", areaCheck, "areaTriCGrid=", areaTriCGrid[iVertex], "diff=", abs(areaCheck-areaTriCGrid[iVertex]))

              dfdx[1] = abs(dfdx[1] - varDenCGrid[iEdgeTest])
              dfdy[2] = abs(dfdy[2] - varDenCGrid[iEdgeTest])
              dfdxAvg = np.mean(dfdx,axis=0)
              dfdyAvg = np.mean(dfdy,axis=0)
              print("          # Term       df/dx                  df/dy")
              print("          ---------------------------------------------------------")
              for iTerm in range(0,nTerms):
                 print("      %5i %-8s: % 20.15e % 20.15e" %(iTerm, terms[iTerm], dfdx[iTerm], dfdy[iTerm]))


        else :

           iVertexTest = 1123
           print("  iVertexTest: ",iVertexTest)
           for iCellOnVertex in range(0, vertexDegree):
              print("  iCell: ", iCellOnVertex, cellsOnVertex[iVertexTest,iCellOnVertex])

           iEdge = edgesOnVertex[iVertexTest,0]
           dcEdge = dcEdges[iEdge]


           basises = ["wachspress","pwl"]
           for basis in basises:

              print("   ",basis)

              # variational
              filein = Dataset("./output_%s_icos4/output.2000.nc" % (basis),"r")

              basisIntegralsU = filein.variables["basisIntegralsU"][:]
              basisIntegralsV = filein.variables["basisIntegralsV"][:]
              basisIntegralsMetric = filein.variables["basisIntegralsMetric"][:]
              cellVerticesAtVertex = filein.variables["cellVerticesAtVertex"][:]
              cellVerticesAtVertex[:] = cellVerticesAtVertex[:] - 1

              filein.close()

              # gradient taylor terms
              dfdx = np.zeros(nTerms) # ideal: (0,1,0,0,0,0)
              dfdy = np.zeros(nTerms) # ideal: (0,0,1,0,0,0)

              for iTerm in range(0,nTerms):

                 # loop over surrounding cells
                 for iSurroundingCell in range(0,vertexDegree):

                    # get the cell number of this cell
                    iCell = cellsOnVertex[iVertexTest,iSurroundingCell]

                    # get the vertexOnCell number of the iVertex velocity point from cell iCell
                    iVelocityVertex = cellVerticesAtVertex[iVertexTest,iSurroundingCell]

                    # loop over the vertices of the surrounding cell
                    for iStressVertex in range(0,nEdgesOnCell[iCell]):

                       iVertex = verticesOnCell[iCell,iStressVertex]

                       dx = xVertex[iVertex] - xVertex[iVertexTest]
                       dy = yVertex[iVertex] - yVertex[iVertexTest]

                       fTerm = pow(dx,dxPow[iTerm]) * pow(dy,dyPow[iTerm])

                       dfdx[iTerm] -= (fTerm * basisIntegralsU[iCell, iVelocityVertex, iStressVertex]) #/ areaTriangle[iVertexTest]
                       dfdy[iTerm] -= (fTerm * basisIntegralsV[iCell, iVelocityVertex, iStressVertex]) #/ areaTriangle[iVertexTest]

                    #CHECK
                    areaCheck = 0.0
                    for i in range(0,nEdgesOnCell[iCell]):
                       for j in range(0,nEdgesOnCell[iCell]):
                          areaCheck += basisIntegralsMetric[iCell,i,j]
                    if iTerm == 0:
                        print("areaCell=", areaCell[iCell], "areaCheck=", areaCheck, "diff=", abs(areaCheck-areaCell[iCell]))

              dfdx[1] = (dfdx[1] -  areaTriangle[iVertexTest])
              dfdy[2] = (dfdy[2] -  areaTriangle[iVertexTest])
              dfdxAvg = np.mean(dfdx,axis=0)
              dfdyAvg = np.mean(dfdy,axis=0)
              print("          # Term       df/dx                  df/dy")
              print("          ---------------------------------------------------------")
              for iTerm in range(0,nTerms):
                 print("      %5i %-8s: % 20.15e % 20.15e" %(iTerm, terms[iTerm], dfdx[iTerm], dfdy[iTerm]))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    error_analysis_stress_divergence()
