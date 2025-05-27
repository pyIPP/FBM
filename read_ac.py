import os, logging
import numpy as np
import parse_ac, mom2rz

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('readAC')
if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)

iso_d  = {(1, 1): 'H', (2, 1) : 'D', (3, 1): 'T', (3, 2): 'HE3', (4, 2): 'HE4'}


class READ_FBM:


    def __init__(self, f_ac, gc=True):

# FBM for GC, FBM_PART for real position
        if gc:
            fbm_lbl = 'FBM'
        else:
            fbm_lbl = 'FBM_PTCL'
        list_read = ['AVGTIM', 'BMVOL', 'EFBM', 'EFBMB', 'XXKSID', fbm_lbl, 'RSURF1', 'YSURF1',
            'RMC', 'YMC', 'XZBEAMS', 'ABEAMS', 'NLFPROD', 'RAXIS', 'YAXIS', 'BDENS2', 'BDENS2',
            'RLIM_PTS', 'YLIM_PTS', 'GOOSE_LOST','RMJION_LOST', 'SPECIES_LOST','WGHT_LOST',
            'XKSID_LOST', 'YNBIEN_LOST','YNBSCE_LOST','YYION_LOST','YZE_LOST',
            'GTIME_LOST','GTIME1_LOST','GTIME2_LOST']

        fbm_d = parse_ac.parse_ac(f_ac, list_read=list_read)

# Dimensions
        n_cells = fbm_d['NFBZNSI']
        n_pit   = fbm_d['NZNBMA']
        n_E     = fbm_d['NZNBME'] #NZNBMEA divided for species
        n_zones = fbm_d['NZONE_FB']
        n_mom   = fbm_d['NMOM']

# Time
        self.fileName = f_ac
        self.dt     = fbm_d['AVGTIM']
        self.time   = fbm_d['time']
        self.nshot  = fbm_d['NSHOT']
        self.bmvol  = fbm_d['BMVOL'][:n_cells]
        vol = 1.e-6*np.sum(self.bmvol)
        logger.info('Volume is %8.4f m^-3', vol)

        rho_lab = np.zeros(n_cells, dtype=int)      # rho index, takes values 0:n_zones-1

#-------------------------
# Space-time grids from AC
#-------------------------

        if 'NLSYM' in fbm_d.keys() and fbm_d['NLSYM']: # up-down symmetric equilibrium
            n_sym = 1
            thbdy0 = 0.
            nthsurf = 101
        else: # asymmetric
            n_sym = 2
            thbdy0 = -np.pi
            nthsurf = 201
        thbdy1 = np.pi

        nthe_tr, nrho_tr = fbm_d['RSURF1'].shape
        nrho_step1 = nrho_tr//n_zones
        self.r_surf = 100*fbm_d['RSURF1'][:, ::nrho_step1].T
        self.z_surf = 100*fbm_d['YSURF1'][:, ::nrho_step1].T
        th_surf = np.linspace(thbdy0, thbdy1, nthe_tr)

        self.x2d  = np.zeros(n_cells)
        self.th2d = np.zeros(n_cells)
        self.r2d  = np.zeros(n_cells)
        self.z2d  = np.zeros(n_cells)

        rcos = fbm_d['RMC'][nrho_step1-1: nrho_tr: nrho_step1, :n_mom, 0] #rho, nmom, cos/sin
        rsin = fbm_d['RMC'][nrho_step1-1: nrho_tr: nrho_step1, :n_mom, 1]
        zcos = fbm_d['YMC'][nrho_step1-1: nrho_tr: nrho_step1, :n_mom, 0]
        zsin = fbm_d['YMC'][nrho_step1-1: nrho_tr: nrho_step1, :n_mom, 1]

        self.rho_grid  = (0.5 + np.arange(n_zones))/float(n_zones)
        rmaj_min = np.zeros(n_zones)         # min(R(rho))

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

# Surface contours

        nxsurf = 2*n_zones + 1
        nrho_step2 = nrho_tr//(nxsurf-1)
        if nrho_step2 == 0 :
            nrho_step2 = 1

        self.thsurf = np.linspace(thbdy0, thbdy1, nthsurf)
        self.xsurf = np.linspace(0, 1, nxsurf)
        self.rsurf = np.zeros((nxsurf, nthsurf))
        self.zsurf = np.zeros((nxsurf, nthsurf))
        rcos = fbm_d['RMC'][(nrho_step2-1):: nrho_step2, :n_mom, 0] #rho, nmom, cos/sin
        rsin = fbm_d['RMC'][(nrho_step2-1):: nrho_step2, :n_mom, 1]
        zcos = fbm_d['YMC'][(nrho_step2-1):: nrho_step2, :n_mom, 0]
        zsin = fbm_d['YMC'][(nrho_step2-1):: nrho_step2, :n_mom, 1]
        self.rsurf, self.zsurf = mom2rz.mom2rz(rcos, rsin, zcos, zsin, theta=self.thsurf)
        self.rsurf[0, :] = fbm_d['RAXIS']
        self.zsurf[0, :] = fbm_d['YAXIS']
        self.rlim_pts= fbm_d['RLIM_PTS']
        self.ylim_pts= fbm_d['YLIM_PTS']

# MC cells: grid bars

        self.rbar = []
        self.zbar = []
        for jrho in range(n_zones):
            for ythe in thb_grid2d[jrho]:
                ithe = np.min(np.where(th_surf > ythe))
                ind = [ithe - 1, ithe ]
                th_ref = th_surf[ind]
                rbar = np.zeros(2)
                zbar = np.zeros(2)
                rbar[0] = np.interp(ythe, th_ref, self.r_surf[jrho  , ind])
                zbar[0] = np.interp(ythe, th_ref, self.z_surf[jrho  , ind])
                rbar[1] = np.interp(ythe, th_ref, self.r_surf[jrho+1, ind])
                zbar[1] = np.interp(ythe, th_ref, self.z_surf[jrho+1, ind])
                self.rbar.append(rbar)
                self.zbar.append(zbar)
               
#----------------------------------------------
# Read distribution, energy, pitch data from AC
#----------------------------------------------

        self.fdist = {}
        self.a_d   = {}
        self.e_d   = {}
        self.eb_d  = {}
        self.btrap = {}
        self.bdens = {}
        self.n_tot = {}
        self.trap_pit   = {}
        self.int_en_pit = {}
        self.int_en_pit_frac_trap = {}
        self.dens_zone  = {}
        self.trap_zone  = {}
        self.dens_vol   = {}
        self.trap_vol   = {}
        self.gtimelost = {}
        self.wghtlost = {}
        self.rmjionlost= {}
        self.yyionlost= {}
        self.xksidlost= {}
        self.yzelost= {}
        self.ynbscelost= {}
        self.gooselost= {}
        self.ta = fbm_d['GTIME1_LOST']
        self.tb= fbm_d['GTIME2_LOST']

        n_spec = len(fbm_d['ABEAMS' ])
        self.species = []
        self.bdens2  = [] # Nr. of fast ions in a cell

        for jspec in range(n_spec):

            if gc:
                fbm = 0.5*fbm_d[fbm_lbl][:, :, :, jspec].transpose((2, 1, 0)) #4th: species
            else:
                fbm = 0.5*fbm_d[fbm_lbl][:, :, :, jspec, 0].transpose((2, 1, 0))
# Skip population with zero fast ions
            if np.max(np.abs(fbm)) == 0:
                continue

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
            self.bdens2.append(fbm_d['BDENS2'][:n_cells, jspec])
            self.fdist[spc_lbl] = fbm
            self.e_d [spc_lbl] = fbm_d['EFBM' ][:, jspec]
            self.eb_d[spc_lbl] = fbm_d['EFBMB'][:, jspec]
            self.a_d[spc_lbl] = np.linspace(-1, 1, num=n_pit)
            (indx, )= np.where(jspec+1==fbm_d['SPECIES_LOST'])
            self.gtimelost[spc_lbl]  = fbm_d['GTIME_LOST'][indx]
            self.wghtlost[spc_lbl]   = fbm_d['WGHT_LOST'][indx]
            self.rmjionlost[spc_lbl] = fbm_d['RMJION_LOST'][indx]
            self.yyionlost[spc_lbl]  = fbm_d['YYION_LOST'][indx]
            self.xksidlost[spc_lbl]  = fbm_d['XKSID_LOST'][indx]
            self.yzelost[spc_lbl]    = fbm_d['YZE_LOST'][indx]
            self.ynbscelost[spc_lbl] = fbm_d['YNBSCE_LOST'][indx]
            self.gooselost[spc_lbl]  = fbm_d[ 'GOOSE_LOST'][indx]
             
# Trapped particles

            fbm_trap = np.zeros((n_cells, n_pit, n_E))
            self.trap_pit[spc_lbl] = 1. - rmaj_min[rho_lab[:]]/self.r2d[:]
            a_squared = self.a_d[spc_lbl]**2
            trap_thresholds = self.trap_pit[spc_lbl][:, np.newaxis]  # shape: (n_cells, 1)
            mask = a_squared <= trap_thresholds  # shape: (n_cells, N)
            fbm_trap[mask] = fbm[mask]

# Integrals

            dE  = np.diff(self.eb_d[spc_lbl])
            dpa = np.repeat(self.a_d[spc_lbl][1] - self.a_d[spc_lbl][0], n_pit) # assuming regular p.a, grid
            dpa_dE = np.outer(dpa, dE)

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
            self.n_tot[spc_lbl] = np.sum(self.bdens[spc_lbl]*vol_zone)
            if self.n_tot[spc_lbl] > 0:
                trap_tot = np.sum(self.btrap[spc_lbl]*vol_zone)
                logger.info('Trapped #%12.4e    Total #%12.4e    Fraction %12.4e', trap_tot, self.n_tot[spc_lbl], trap_tot/self.n_tot[spc_lbl])
                logger.info('Volume averaged fast ion density #12.4e m^-3', self.n_tot[spc_lbl]/vol)
                self.int_en_pit_frac_trap[spc_lbl] = int_en_pit_trap/self.int_en_pit[spc_lbl]
            else:
                logger.error( 'WARNING! FBM data missing in input file')
        
        self.bdens2 = np.array(self.bdens2)


if __name__ == '__main__':

    import config
    runid = '29783A01'

    shot = runid[:-3]
    tail = runid[-3:]

    f_ac  = '%s/%s/%s/%s.DATA1'    %(config.tr_clientDir, shot, tail, runid)

    fbm = READ_FBM(f_ac)
    print(fbm.dt)
