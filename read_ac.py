import os
import numpy as np
import matplotlib.pylab as plt
import parse_ac, mom2rz

iso_d  = {(1, 1): 'H', (2, 1) : 'D', (3, 1): 'T', (3, 2): 'HE3', (4, 2): 'HE4'}

class READ_FBM:


    def __init__(self, f_ac, gc=True):

# FBM for GC, FBM_PART for real position
        if gc:
            fbm_lbl = 'FBM'
        else:
            fbm_lbl = 'FBM_PTCL'
        list_read = ['BMVOL', 'EFBM', 'EFBMB', 'XXKSID', fbm_lbl, 'RSURF1', 'YSURF1', 'RMC', 'YMC', 'XZBEAMS', 'ABEAMS', 'NLFPROD']

        fbm_d = parse_ac.parse_ac(f_ac, list_read=list_read)

# Dimensions
        n_cells = fbm_d['NFBZNSI']
        n_pit   = fbm_d['NZNBMA']
        n_E     = fbm_d['NZNBME'] #NZNBMEA divided for species
        n_zones = fbm_d['NZONE_FB']
        n_mom   = fbm_d['NMOM']

# Time
        self.time = fbm_d['time']
        
        self.bmvol = fbm_d['BMVOL'][:n_cells]
        vol = 1.e-6*np.sum(self.bmvol)
        print('Volume is %8.4f m^-3' %vol)

        rho_lab = np.zeros(n_cells, dtype=int)      # rho index, takes values 0:n_zones-1

#-------------------------
# Space-time grids from AC
#-------------------------

        if fbm_d['NLSYM']: # up-down symmetric equilibrium
            n_sym = 1
            thbdy0 = 0.
            nthsurf = 101
        else: # asymmetric
            n_sym = 2
            thbdy0 = -np.pi
            nthsurf = 201
        thbdy1 = np.pi

        nxsurf = 2*n_zones + 1
        print('N cells', n_cells, n_zones)

        nthe_tr, nrho_tr = fbm_d['RSURF1'].shape
        nrho_step = nrho_tr//n_zones
        self.rsurf = 100*fbm_d['RSURF1'][:, ::nrho_step].T
        self.zsurf = 100*fbm_d['YSURF1'][:, ::nrho_step].T

        self.x2d  = np.zeros(n_cells)
        self.th2d = np.zeros(n_cells)
        self.r2d  = np.zeros(n_cells)
        self.z2d  = np.zeros(n_cells)

        rcos = fbm_d['RMC'][3: nrho_tr: nrho_step, :n_mom, 0] #rho, nmom, cos/sin
        rsin = fbm_d['RMC'][3: nrho_tr: nrho_step, :n_mom, 1]
        zcos = fbm_d['YMC'][3: nrho_tr: nrho_step, :n_mom, 0]
        zsin = fbm_d['YMC'][3: nrho_tr: nrho_step, :n_mom, 1]

        self.rho_grid  = (0.5 + np.arange(n_zones))/float(n_zones)
        rmaj_min = np.zeros(n_zones)         # min(R(rho))
        thsurf = np.linspace(thbdy0, thbdy1, nthe_tr)
        xsurf  = np.linspace(0, 1, nxsurf)
        rho_lab = np.zeros(n_cells, dtype=int)      # rho index, takes values 0:n_zones-1

        thb_grid2d = []
        min_ind = 0
        for jzone in range(n_zones):
            nth = 2*n_sym*(1 + jzone)
            the_step = (thbdy1 - thbdy0)/float(nth)
            th_loc = thbdy0 + the_step*(0.5 + np.arange(nth))
            r2d_loc, z2d_loc = mom2rz.mom2rz(rcos[jzone], rsin[jzone], zcos[jzone], zsin[jzone], theta=th_loc)
            rmaj_min[jzone] = min(r2d_loc)
            thb_grid2d.append(th_loc - 0.5*the_step)

            ind = range(min_ind, min_ind + nth)
            self.x2d[ind] = self.rho_grid[jzone]
            rho_lab  [ind] = jzone
            self.th2d[ind] = th_loc
            self.r2d [ind] = r2d_loc
            self.z2d [ind] = z2d_loc
            min_ind += nth

# MC cells: grid bars

        self.rbar = []
        self.zbar = []
        for jrho in range(n_zones):
            for ythe in thb_grid2d[jrho]:
                ithe = np.min(np.where(thsurf > ythe))
                ind = [ithe - 1, ithe ]
                th_ref = thsurf[ind]
                rbar = np.zeros(2)
                zbar = np.zeros(2)
                rbar[0] = np.interp(ythe, th_ref, self.rsurf[jrho  , ind])
                zbar[0] = np.interp(ythe, th_ref, self.zsurf[jrho  , ind])
                rbar[1] = np.interp(ythe, th_ref, self.rsurf[jrho+1, ind])
                zbar[1] = np.interp(ythe, th_ref, self.zsurf[jrho+1, ind])
                self.rbar.append(rbar)
                self.zbar.append(zbar)
               
#----------------------------------------------
# Read distribution, energy, pitch data from AC
#----------------------------------------------

        self.fdist  = {}
        self.a_d    = {}
        self.e_d    = {}
        self.eb_d   = {}
        self.bdens  = {}
        self.btrap  = {}
        self.trap_pit   = {}
        self.int_en_pit = {}
        self.int_en_pit_frac_trap = {}
        self.dens_zone  = {}
        self.trap_zone  = {}
        self.dens_vol   = {}
        self.trap_vol   = {}

        n_spec = len(fbm_d['ABEAMS' ])
        self.species = []
        for jspec in range(n_spec):
            if jspec < fbm_d['NSBEAM']:
                if fbm_d['NLFPROD'][jspec]:
                    prod_lbl = 'FUS'
                else:
                    prod_lbl = 'NBI'
            else:
                prod_lbl = 'RFI'

            mass   = int(fbm_d['ABEAMS' ][jspec] + 0.2)
            charge = int(fbm_d['XZBEAMS'][jspec] + 0.001)

            iso_lbl  = iso_d[(mass, charge)]

            spc_lbl = '%s_%s' %(iso_lbl, prod_lbl)
            self.species.append(spc_lbl)
            print(spc_lbl)
            self.e_d [spc_lbl] = fbm_d['EFBM'] [:, jspec]
            self.eb_d[spc_lbl] = fbm_d['EFBMB'][:, jspec]
            self.a_d[spc_lbl] = fbm_d['XXKSID']
            dE  = np.diff(self.eb_d[spc_lbl])
            dpa = np.repeat(self.a_d[spc_lbl][1] - self.a_d[spc_lbl][0], n_pit) # assuming regular p.a, grid
            dpa_dE = np.outer(dpa, dE)

            if gc:
                fbm = 0.5*fbm_d[fbm_lbl][:, :, :, jspec].transpose((2, 1, 0)) #4th: species
            else:
                fbm = 0.5*fbm_d[fbm_lbl][:, :, :, jspec, 0].transpose((2, 1, 0))
            self.fdist[spc_lbl] = fbm

            fbm_trap = np.zeros((n_cells, n_pit, n_E))
            self.trap_pit[spc_lbl] = 1. - rmaj_min[rho_lab[:]]/self.r2d[:]
            for jcell in range(n_cells):
                ind_trap = (self.a_d[spc_lbl]**2 <= self.trap_pit[spc_lbl][jcell])
                fbm_trap[jcell, ind_trap, :] = fbm[jcell, ind_trap, :]

# Integrals

            self.int_en_pit[spc_lbl] = np.tensordot(fbm, dpa_dE, axes=((1, 2), (0, 1)))
            int_en_pit_trap = np.tensordot(fbm_trap, dpa_dE, axes=((1, 2), (0, 1)))

            self.dens_zone[spc_lbl] = np.zeros((n_zones, n_pit, n_E))
            self.trap_zone[spc_lbl] = np.zeros((n_zones, n_pit, n_E))

            self.dens_vol[spc_lbl] = np.tensordot(fbm     , self.bmvol, axes=(0, 0))/np.sum(self.bmvol)
            self.trap_vol[spc_lbl] = np.tensordot(fbm_trap, self.bmvol, axes=(0, 0))/np.sum(self.bmvol)
            vol_zone = np.zeros(n_zones)
            for jrho in range(n_zones):
                indrho = (rho_lab == jrho)
                vol_zone[jrho] = np.sum(self.bmvol[indrho])
                self.dens_zone[spc_lbl][jrho, :, :] = np.tensordot(fbm[     indrho, :, :], self.bmvol[indrho], axes = (0, 0))
                self.trap_zone[spc_lbl][jrho, :, :] = np.tensordot(fbm_trap[indrho, :, :], self.bmvol[indrho], axes = (0, 0))
                self.dens_zone[spc_lbl][jrho, :, :] *= 1/vol_zone[jrho]
                self.trap_zone[spc_lbl][jrho, :, :] *= 1/vol_zone[jrho]

            self.bdens[spc_lbl] = np.tensordot(self.dens_zone[spc_lbl], dpa_dE, axes=((1, 2), (0, 1)))
            self.btrap[spc_lbl] = np.tensordot(self.trap_zone[spc_lbl], dpa_dE, axes=((1, 2), (0, 1)))

            part_tot = np.sum(self.bdens[spc_lbl]*vol_zone)
            trap_tot = np.sum(self.btrap[spc_lbl]*vol_zone)
            print('Trapped #%12.4e    Total #%12.4e    Fraction %12.4e' %(trap_tot, part_tot, trap_tot/part_tot))
            print('Volume averaged fast ion density #12.4e m^-3' %(part_tot/vol))
            self.int_en_pit_frac_trap[spc_lbl] = int_en_pit_trap/self.int_en_pit[spc_lbl]


if __name__ == '__main__':

    runid = '29783A01'

    shot = runid[:-3]
    tail = runid[-3:]

    run_dir = '%s/tr_client/AUGD/%s/%s' %(os.getenv('HOME'), shot, tail)
    f_ac  = '%s/%s.DATA1'    %(run_dir, runid)

    fbm = READ_FBM(f_ac)
