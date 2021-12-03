import os
import numpy as np
import xarray
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

from pyremap import MpasMeshDescriptor, Remapper

def stress_divergence_scaling_numerical():

    tests = ['wachspress', 'pwl']
    schemes = ['C-grid']
    nTests = len(tests)
    grids = ['icos4','icos5','icos6','icos7']
    nGrids = len(grids)
    dx = [2562, 10242, 40962, 163842]
    dx1 = ['4','5','6','7']
    dxRates = [2562, 10242, 40962]
    order1 = [0.00004, 0.00001, 0.0000025]
    sqrtOrder1 = [0.00008, 0.00004, 0.00002]
    Linf = np.zeros([nTests, nGrids])
    L2 = np.zeros([nTests, nGrids])
    
    for j in range(nTests):
        test = tests[j]
        for s in range(nGrids):
            grid = grids[s]

#-- setup remapping descriptors for the MPAS meshes to process  
            src_descriptor = MpasMeshDescriptor(
                #"../../" + grid + "/base_mesh.nc", meshName="oQUn")
                "icos9/base_mesh.nc", meshName="oQU9")
                #test + "/convergence/icos9/base_mesh.nc", meshName="oQU8") 

            dst_descriptor = MpasMeshDescriptor(
                #"../../icos8/base_mesh.nc", meshName="oQU8")
                #test + "/convergence/" + grid + "/base_mesh.nc", meshName="oQUn")
                grid + "/base_mesh.nc", meshName="oQUn")

#-- create the remapping object and save the weights to a file
            remapper = Remapper(
                src_descriptor, dst_descriptor, grid + "Weights.nc")

#-- build the actual remapping weights, uses weights.nc if it
#-- exists already, so be careful to delete if you want to 
#-- build new weights from scratch...
#-- METHOD="CONSERVE" is a finite-volume style remap
#-- METHOD="BILINEAR" is a finite-difference (non-conservative)
#-- approach
#-- See ESMF_regrid for more information...
            remapper.build_mapping_file(method="conserve2nd")

#-- setup some example data to remap
#-- load the mesh files just to get the cell lon-lat positions
            #msrc = xarray.open_dataset("../../" + grid + "/output.nc")
            #mdst = xarray.open_dataset("../../icos8/output.nc")
            msrc = xarray.open_dataset("output_" + test + "_icos9/output.2000.nc")
            mdst = xarray.open_dataset("output_" + test + "_" + grid + "/output.2000.nc")
   
            dsrc = xarray.Dataset()
            dsrc["stressDivergenceUCGrid"] = msrc["stressDivergenceUCGrid"] 

#-- In practice, this is where you'd read in the time series
#-- data that you'd like to remap...

#-- remap any field in DSRC to DDST using the weights computed
#-- above
            ddst = remapper.remap(dsrc) 

#-- check conservation of remapping operation
        #ssrc = float(numpy.sum(
        #    dsrc["layerThickness"] * msrc["areaCell"]))

        #sdst = float(numpy.sum(
        #    ddst["layerThickness"] * mdst["areaCell"]))

        #print("integral-src =", ssrc)
        #print("integral-dst =", sdst)
        #print("rel.-error = ", (ssrc - sdst) / max(abs(ssrc), abs(sdst)))

            areas = mdst.variables['variationalDenominatorCGrid'][:]
            var = mdst.variables['stressDivergenceUCGrid'][0,:]
            sol = ddst.variables['stressDivergenceUCGrid'][0,:]
            dif = abs(var - sol)
            multDen = (sol**2)*areas
            multNum = (dif**2)*areas 
            Linf[j, s] = np.max(dif[:])/np.max(abs(sol[:]))
            denL2 = np.sum(multDen[:])/np.sum(areas[:])
            numL2 = np.sum(multNum[:])/np.sum(areas[:])
            L2[j, s] = np.sqrt(numL2)/np.sqrt(denL2)
    
    plt1 = plt.figure(1)
    plt.loglog(dx, Linf[0, :], '-x', label='L_inf_HCt')
    plt.loglog(dx, L2[0, :], '-x', label='L_2_HCt')
    plt.loglog(dxRates, order1, 'k', label='slope=-2')
    plt.loglog(dxRates, sqrtOrder1, 'k--', label='slope=-1')
    plt.title('Convergence rate of div of stress')
    plt.xlim(1.e3, 1.e6)
    plt.ylim(1.e-6, 1.e-1)
    plt.ylabel('L_2 and L_infinity error norms')
    plt.legend()
    plt.grid()
    plt.xticks(dx, dx)
    plt.xlabel('# of cells')
   
    plt1.savefig('convergence_div_of_stress.png')

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    stress_divergence_scaling_numerical()

