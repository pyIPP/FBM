import os, sys, datetime
import numpy as np
import read_ac
from scipy.io import netcdf


def ac2cdf(f_ac, gc=True):

    print(f_ac)
    fbm = read_ac.READ_FBM(f_ac)

    n_cells = len(fbm.r2d)
    
    trdir, basename = os.path.split(f_ac)
    tmp = basename.split('.')
    runid = tmp[0]
    t_id = tmp[1][4:]
    f_cdf = '%s/%s_fitest_%s.cdf' %(trdir, runid, t_id)
    
# Create NetCDF file
    f = netcdf.netcdf_file(f_cdf, 'w', mmap=False)
    f.history = "Created " + datetime.datetime.today().strftime("%d/%m/%y")

# Grids

    f.createDimension('INDEX_2D', n_cells)
    f.createDimension('TIME', 1)

    mc = f.createVariable('INDEX_2D', np.int32, ('INDEX_2D', ))
    mc[:] = 1 + np.arange(n_cells)
    mc.units = ''
    mc.long_name = 'Energy grid'

    time = f.createVariable('TIME', np.float32, ('TIME',))
    time[:] = fbm.time
    time.units = 's'
    time.long_name = 'Time'

    r2d = f.createVariable('R2D', np.float32, ('INDEX_2D', ))
    r2d[:] = fbm.r2d
    r2d.units = 'cm'
    r2d.long_name = 'R of MC cells'

    z2d = f.createVariable('Z2D', np.float32, ('INDEX_2D', ))
    z2d[:] = fbm.z2d
    z2d.units = 'cm'
    z2d.long_name = 'Z of MC cells'

    x2d = f.createVariable('X2D', np.float32, ('INDEX_2D', ))
    x2d[:] = fbm.x2d
    x2d.units = ''
    x2d.long_name = 'rho_tor of MC cells'

    th2d = f.createVariable('TH2D', np.float32, ('INDEX_2D', ))
    th2d[:] = fbm.th2d
    th2d.units = ''
    th2d.long_name = 'Theta of MC cells'

    bmvol = f.createVariable('BMVOL', np.float32, ('INDEX_2D', ))
    bmvol[:] = fbm.bmvol
    bmvol.units = '1/cm**3'
    bmvol.long_name = 'Volume of MC cells'

    for spc_lbl in fbm.species:

        n_cells, n_pit, n_E = fbm.fdist[spc_lbl].shape

        Elbl  = 'E_%s'  %spc_lbl
        Albl  = 'A_%s'  %spc_lbl
        EBlbl = 'EB_%s' %spc_lbl
        ABlbl = 'AB_%s' %spc_lbl
        Flbl  = 'F_%s'  %spc_lbl

        ab = np.zeros(n_pit + 1)
        ab[1:-1] = 0.5*(fbm.a_d[spc_lbl][1:] + fbm.a_d[spc_lbl][:-1])
        ab[0]  = 2.*ab[1]  - ab[2]
        ab[-1] = 2.*ab[-2] - ab[-3]

        f.createDimension(Albl, n_pit)
        f.createDimension(Elbl, n_E)
        f.createDimension(ABlbl, n_pit+1)
        f.createDimension(EBlbl, n_E+1)

        E = f.createVariable(Elbl, np.float32, (Elbl, ))
        E[:] = fbm.e_d[spc_lbl]
        E.units = 'eV'
        E.long_name = 'Energy grid'

        EB = f.createVariable(EBlbl, np.float32, (EBlbl, ))
        EB[:] = fbm.eb_d[spc_lbl]
        EB.units = 'eV'
        EB.long_name = 'Energy boundary grid'

        A = f.createVariable(Albl, np.float32, (Albl, ))
        A[:] = fbm.a_d[spc_lbl]
        A.units = ''
        A.long_name = 'Pitch angle grid'

        AB = f.createVariable(ABlbl, np.float32, (ABlbl, ))
        AB[:] = ab
        AB.units = ''
        AB.long_name = 'Pitch angle boundary grid'

# Variables

        fb = f.createVariable(Flbl, np.float32, ('INDEX_2D', Albl, Elbl))
        fb[:] = 2.*fbm.fdist[spc_lbl]
        fb.units = '1/cm**3/s/sterad'
        fb.long_name = '%s fast ion distribution function' %spc_lbl

    f.close()
    print('Stored %s' %f_cdf)


if __name__ == '__main__':


    f_cdf1 = '/afs/ipp/home/g/git/tr_client/AUGD/29783/A01/29783A01_fi_1.cdf'
    cv1 = netcdf.netcdf_file(f_cdf1, 'r', mmap=False)
    for key, val in cv1.dimensions.items():
        print(key, val)
    for key, val in cv1.variables.items():
        print(key, val)
#        for attr in val.__dict__.keys():
        print(val.dimensions)
        print(val._attributes)
#        print(val.long_name)

    f_ac = '/afs/ipp/home/g/git/tr_client/AUGD/29783/A01/29783A01.DATA1'
    ac2cdf(f_ac)
    f_cdf2 = '/afs/ipp/home/g/git/tr_client/AUGD/29783/A01/29783A01_fitest_1.cdf'
    cv2 = netcdf.netcdf_file(f_cdf2, 'r', mmap=False).variables
    for key, val in cv2.items():
        print(key)
#        for attr in val.__dict__.keys():
#            print(attr)
#        print(val.units)
#        print(val.long_name)
