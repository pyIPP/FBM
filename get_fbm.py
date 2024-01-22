#!/usr/bin/env python

import os, datetime, argparse
import numpy as np
import read_ac
from scipy.io import netcdf_file

flt_type = np.float64

def ac2cdf(f_ac, gc=True):

    print(f_ac)
    fbm = read_ac.READ_FBM(f_ac)

    n_cells = len(fbm.r2d)
    nrho_surf, nthe_surf = fbm.rsurf.shape
    n_species = len(fbm.species)

    trdir, basename = os.path.split(f_ac)
    tmp = basename.split('.')
    runid = tmp[0]
    t_id = tmp[1][4:]
    f_cdf = '%s/%s_fi_%s.cdf' %(trdir, runid, t_id)

# Create NetCDF file
    f = netcdf_file(f_cdf, 'w', mmap=False)
    f.history = "Created " + datetime.datetime.today().strftime("%d/%m/%y")

# Grids

    f.createDimension('INDEX_2D', n_cells)
    f.createDimension('TIME', 1)
    f.createDimension('XSURF' , nrho_surf)
    f.createDimension('THSURF', nthe_surf)
    f.createDimension('nchar', 8)

    mc = f.createVariable('INDEX_2D', np.int32, ('INDEX_2D', ))
    mc[:] = 1 + np.arange(n_cells)
    mc.units = ''
    mc.long_name = 'Energy grid'

    time = f.createVariable('TIME', flt_type, ('TIME', ))
    time[:] = fbm.time
    time.units = 's'
    time.long_name = 'Time'

    xsurf = f.createVariable('XSURF', flt_type, ('XSURF', ))
    xsurf[:] = fbm.xsurf
    xsurf.units = ''
    xsurf.long_name = 'xsurf grid (for flux surfaces)'

    thsurf = f.createVariable('THSURF', flt_type, ('THSURF', ))
    thsurf[:] = fbm.thsurf
    thsurf.units = 'Radians'
    thsurf.long_name = 'thsurf grid (for flux surfaces)'

# Variables
 
    gcf = f.createVariable('GCFLAG', np.int32, ('TIME', ))
    gcf[:] = gc
    gcf.units = ''
    gcf.long_name = 'Distributoin at gyro-center'
   
    dt_avg = f.createVariable('DT_AVG', flt_type, ('TIME', ))
    dt_avg[:] = fbm.dt
    dt_avg.units = 's'
    dt_avg.long_name = 'Averaging time'

    nshot = f.createVariable('nshot', np.int32, ('TIME', ))
    nshot[:] = fbm.nshot
    nshot.units = ''
    nshot.long_name = 'Shot number'

    r2d = f.createVariable('R2D', flt_type, ('INDEX_2D', ))
    r2d[:] = fbm.r2d
    r2d.units = 'cm'
    r2d.long_name = 'R of MC cells'

    z2d = f.createVariable('Z2D', flt_type, ('INDEX_2D', ))
    z2d[:] = fbm.z2d
    z2d.units = 'cm'
    z2d.long_name = 'Z of MC cells'

    x2d = f.createVariable('X2D', flt_type, ('INDEX_2D', ))
    x2d[:] = fbm.x2d
    x2d.units = ''
    x2d.long_name = 'rho_tor of MC cells'

    th2d = f.createVariable('TH2D', flt_type, ('INDEX_2D', ))
    th2d[:] = fbm.th2d
    th2d.units = ''
    th2d.long_name = 'Theta of MC cells'

    bmvol = f.createVariable('BMVOL', flt_type, ('INDEX_2D', ))
    bmvol[:] = fbm.bmvol
    bmvol.units = 'cm**3'
    bmvol.long_name = 'Volume of MC cells'

    rsurf = f.createVariable('RSURF', flt_type, ('XSURF', 'THSURF'))
    rsurf[:] = fbm.rsurf
    rsurf.units = 'cm'
    rsurf.long_name = 'flux surface R locations'

    zsurf = f.createVariable('ZSURF', flt_type, ('XSURF', 'THSURF'))
    zsurf[:] = fbm.zsurf
    zsurf.units = 'cm'
    zsurf.long_name = 'flux surface Z locations'

    bdens2 = f.createVariable('bdens2', flt_type, ('INDEX_2D', ))
    bdens2[:] = fbm.bdens2
    bdens2.units = '1/cm**3'
    bdens2.long_name = 'Fast ion density in a MC cell'

    trunid = f.createVariable('TRANSP_RUNID', 'S1', ('nchar', ))
    trunid[:] = [x for x in runid]
    trunid.long_name = 'TRANSP runID'

    for jspc, spc_lbl in enumerate(fbm.species):

        n_cells, n_pit, n_E = fbm.fdist[spc_lbl].shape
        nSpcLabel = len(spc_lbl)

        Elbl  = 'E_%s'  %spc_lbl
        Albl  = 'A_%s'  %spc_lbl
        EBlbl = 'EB_%s' %spc_lbl
        ABlbl = 'AB_%s' %spc_lbl
        Flbl  = 'F_%s'  %spc_lbl
        NTOTlbl = 'NTOT_%s' %spc_lbl
        
        ab = np.zeros(n_pit + 1)
        ab[1:-1] = 0.5*(fbm.a_d[spc_lbl][1:] + fbm.a_d[spc_lbl][:-1])
        ab[0]  = 2.*ab[1]  - ab[2]
        ab[-1] = 2.*ab[-2] - ab[-3]

        f.createDimension(Albl, n_pit)
        f.createDimension(Elbl, n_E)
        f.createDimension(ABlbl, n_pit+1)
        f.createDimension(EBlbl, n_E+1)
        f.createDimension('SpcStr', nSpcLabel)

        SpcIndex = 'SPECIES_%d' %(jspc+1)
        SpcLabel = f.createVariable(SpcIndex, 'S1', ('SpcStr', ))
        SpcLabel[:] = [x for x in spc_lbl]
        SpcLabel.units = ''
        SpcLabel.long_name = 'Species label'

        E = f.createVariable(Elbl, flt_type, (Elbl, ))
        E[:] = fbm.e_d[spc_lbl]
        E.units = 'eV'
        E.long_name = 'Energy grid'

        EB = f.createVariable(EBlbl, flt_type, (EBlbl, ))
        EB[:] = fbm.eb_d[spc_lbl]
        EB.units = 'eV'
        EB.long_name = 'Energy boundary grid'

        A = f.createVariable(Albl, flt_type, (Albl, ))
        A[:] = fbm.a_d[spc_lbl]
        A.units = ''
        A.long_name = 'Pitch angle grid'

        AB = f.createVariable(ABlbl, flt_type, (ABlbl, ))
        AB[:] = ab
        AB.units = ''
        AB.long_name = 'Pitch angle boundary grid'

# Total amount of fast ions

        Ntot = f.createVariable(NTOTlbl, flt_type, ('TIME', ))
        Ntot[:] = fbm.n_tot[spc_lbl]
        Ntot.units = ''
        Ntot.long_name = 'total number of ions'

# Distribution function

        fb = f.createVariable(Flbl, flt_type, ('INDEX_2D', Albl, Elbl))
        fb[:] = 2.*fbm.fdist[spc_lbl]
        fb.units = '1/cm**3/s/sterad'
        fb.long_name = '%s fast ion distribution function' %spc_lbl

    f.close()
    print('Stored %s' %f_cdf)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert AC fbm to NetCDF')
    parser.add_argument('-r', '--runid', help='TRANSP runID', required=True)
    parser.add_argument('-t', '--t_id' , help='t_id', type=int, required=False, default=1)

    args = parser.parse_args()

    shot = args.runid[:5]
    tail = args.runid[5:]
    f_ac = '%s/AUGD/%s/%s/%s.DATA%d' %(os.getenv('TR_CLIENT'), shot, tail, args.runid, args.t_id)
    ac2cdf(f_ac)
