#/p/transpusers/mgorelen/work/AUGD/41385M01/41385M01.DATA1
import numpy as np
import pandas as pd
from scipy.io import netcdf_file
try:
    import Tkinter as tk
    import ttk
except:
    import tkinter as tk
    from tkinter import ttk

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
try:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as nt2tk
except:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg as nt2tk


def get_lost(runid, fbm, topframe=None):

# Vessel compoments for plot
    xext = ( fbm.rlim_pts.min() + fbm.rlim_pts.max())*0.05
    yext = (-fbm.ylim_pts.min() + fbm.ylim_pts.max())*0.05
    xtop_lim = np.array ([-fbm.rlim_pts.max() - xext, fbm.rlim_pts.max() + xext])
    xtop_lim = np.round(xtop_lim).astype(int)
    xpol_lim = np.array ([fbm.rlim_pts.min() - xext, fbm.rlim_pts.max() + xext])
    ypol_lim = np.array ([fbm.ylim_pts.min() - yext, fbm.ylim_pts.max() + yext])
    xpol_lim = np.round(xpol_lim).astype(int)
    ypol_lim = np.round(ypol_lim).astype(int)
    rlin = fbm.rlim_pts
    zlin = fbm.ylim_pts

# Get 
    nr = 141
    nz = 101
    ntheta = 101
    Rmaj = fbm.rsurf[0, 0]
    Rmin = fbm.rlim_pts.min()
    Rmax = fbm.rlim_pts.max()
    zmin = fbm.zsurf.min()
    zmax = fbm.zsurf.max()
    R_grid = np.linspace(Rmin, Rmax, nr)
    z_grid = np.linspace(zmin, zmax, nz)

    phi_tor = np.linspace(0, 2*np.pi, ntheta)
    cosp = np.cos(phi_tor)
    sinp = np.sin(phi_tor)
    Rlbl = 'R [cm]'
    zlbl = 'Z [cm]'
    fsize = 12
    colors = ('b','c', 'g', 'y', 'r')
    sizes=(7,8,9,12,14)
    nrows=1
    n_pwr = 5

    xgrid, ygrid = np.meshgrid(R_grid, z_grid, indexing='ij')

    tbm1=fbm.ta
    tbm2=fbm.tb

    print('PLOT_LOST')
    src_arr = fbm.species
    
    print('Sources: ', src_arr)

    comp_lbl = {1: 'P', 2: 'NP', 3: 'Total'}
    n_src  = len(src_arr)
# One tab for each source

    if topframe is None:
        topframe = tk.Toplevel()
        xmax = topframe.winfo_screenwidth()
        ymax = topframe.winfo_screenheight()
        width  = min(1500, int(0.95*xmax))
        height = min(960, int(0.88*ymax)) 
        topframe.title('Losses Dist.')
        topframe.geometry('%dx%d' %(width, height))
    nbsource = ttk.Notebook(topframe, name='nb source')
    nbsource.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    for jsrc in src_arr:
        weight = np.array(fbm.wghtlost[jsrc])
        Rj = np.array(fbm.rmjionlost[jsrc])
        zj = np.array(fbm.yyionlost[jsrc])
        pwr = np.array(fbm.xksidlost[jsrc])
        Elost = np.array(fbm.yzelost[jsrc])
        j_nbi = np.array(fbm.ynbscelost[jsrc]).astype(int)
        n_lost = len(Rj)
        print('# of MC particles lost: %d' %n_lost)    
        if n_lost == 0:
            continue
        j_comp = np.zeros(n_lost, dtype=np.int32)
        j_comp_tot = np.zeros(n_lost, dtype=np.int32)
    
        (indx_p , ) = np.where(fbm.gooselost[jsrc] < 1.000001)
        (indx_np, ) = np.where(fbm.gooselost[jsrc] >= 1.000001)
        j_comp[indx_p ] = 1
        j_comp[indx_np] = 2
        j_comp_tot [indx_p ] = 3
        j_comp_tot [indx_np] = 3
        comp_arr = np.unique(j_comp)
        n_comp = len(comp_arr)
        ncols = n_comp + 1

        pwr_lost_tot = 0.0
        pwr_lost = np.zeros(n_lost)
        il = 0
        while il < n_lost: 
            pwr_lost[il] = np.array(weight[il]*Elost[il]/(tbm2 - tbm1)*1.6021766*1.e-19)
            pwr_lost_tot += pwr_lost[il] 
            il += 1

        pwr_min = np.min(pwr_lost[np.nonzero(pwr_lost)])
        pwr_max = np.max(pwr_lost)
        print(jsrc + 'Total power lost, [Watts] ', pwr_lost_tot)
        pwr_edges = np.linspace(np.log10(pwr_min), np.log10(pwr_max), n_pwr+1)

        frame_source = tk.Frame(nbsource)
        lbl = jsrc
        
        nblost = ttk.Notebook(frame_source, name='nblost')
        nblost.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
        frame_lost = tk.Frame(nblost)

        fig_lost = Figure(figsize=(14., 8.45), dpi=100)
        can_lost = FigureCanvasTkAgg(fig_lost, master=frame_lost)
        can_lost._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        fig_lost.clf()
        fig_lost.subplots_adjust(left=0.05, bottom=0.08, right=0.97, top=0.97, \
                                  wspace=0.15, hspace=0)

        fig_lost.text(0.5, 0.01, jsrc+' Power losses for %s, t =%6.3f s ' %(runid, tbm1), ha='center',fontsize=16,fontweight='bold')
 
        jsplot = 1

# Overplot 3 species mix
# Poloidal section

        axpol = fig_lost.add_subplot(nrows, ncols, jsplot, aspect='equal')
        axpol.set_title(jsrc+' total losses',fontsize=fsize,fontweight='bold')
# Lost locations
            
        for jcol, jcomp in enumerate(comp_arr):
            ind = (j_comp == jcomp)
            axpol.plot(Rj[ind]  , zj[ind]  , '%so' %colors[jcol], label=comp_lbl[jcomp],markersize=sizes[jcol])

        axpol.legend(loc=2, numpoints=1, prop={'size': 10})

# Tokamak components
        try:
            import aug_sfutils as sf
            gc_d = sf.getgc()
            m2cm = 100.
        except:
            pass

        axpol.set_xlabel(Rlbl, fontsize=fsize,fontweight='bold')
        axpol.set_ylabel(zlbl, fontsize=fsize,fontweight='bold')
        axpol.set_xlim(xpol_lim)
        axpol.set_ylim(ypol_lim)
        
        if 'gc_d' in locals():
            for gc in gc_d.values():
                axpol.plot(m2cm*gc.r, m2cm*gc.z, 'b-')
        
        if 'fbm' in locals():
            for irho in range(fbm.rsurf.shape[0]):
                axpol.plot(fbm.rsurf[irho, :], fbm.zsurf[irho, :], 'r-', linewidth=0.15)
            for jbar, myr in enumerate(fbm.rbar):
                axpol.plot(myr, fbm.zbar[jbar], 'r-', linewidth=0.5)

        if 'rlin' in locals():
            axpol.plot(rlin, zlin, 'g-', linewidth=2.5)
            
# For each species overplot 5 pwr angle range

        jsplot = 2

        for jcomp in comp_arr:
            axpol = fig_lost.add_subplot(nrows, ncols, jsplot, aspect='equal')
            
            axpol.set_title('%s Pwr. losses, Watts' %comp_lbl[jcomp], \
                            fontsize=fsize,fontweight='bold')

            axpol.set_xlabel(Rlbl, fontsize=fsize,fontweight='bold')
            axpol.set_ylabel(zlbl, fontsize=fsize,fontweight='bold')
            axpol.set_xlim(xpol_lim)
            axpol.set_ylim(ypol_lim)
            ind1 = (j_comp == jcomp)
            jpwr =n_pwr-1
            
            while jpwr >=0 :
                p1 = pwr_edges[jpwr]
                p2 = pwr_edges[jpwr+1]
                ind = ind1 & (pwr_lost>  np.float_power(10,p1)) & (pwr_lost <= np.float_power(10,p2))
                axpol.plot(Rj[ind], zj[ind], '%so' %colors[jpwr], \
                           label='%s < Pwr < %s' %(round(np.float_power(10,p1),2), \
                                                   round(np.float_power(10,p2),2)),markersize=sizes[jpwr])
                jpwr =jpwr-1

            if 'gc_d' in locals():
                for gc in gc_d.values():
                    axpol.plot(m2cm*gc.r, m2cm*gc.z, 'b-')            

            if 'fbm' in locals():
                for irho in range(fbm.rsurf.shape[0]):
                    axpol.plot(fbm.rsurf[irho, :], fbm.zsurf[irho, :], 'r-', linewidth=0.15)
                for jbar, myr in enumerate(fbm.rbar):
                    axpol.plot(myr, fbm.zbar[jbar], 'r-', linewidth=0.5)

            if 'rlin' in locals():
                axpol.plot(rlin, zlin, 'g-', linewidth=2.5)
 
            axpol.legend(loc=2, numpoints=1, prop={'size': 10})
            jsplot += 1

        toolbar = nt2tk(can_lost, frame_lost)
        toolbar.update()

# Add tabs

        nblost.add(frame_lost, text='Lost location')
        nbsource.add(frame_source, text=lbl)


if __name__ == "__main__":

    flost = '/afs/ipp/home/g/git/tr_client/AUGD/29783/A01/29783A01_lost.cdf1'
    lost_tk = tk.Tk()
    lost_tk.title('Lost location')
    lost_tk.geometry('1500x940')
    lost_tk.option_add("*Font", "Helvetica")
    read_lost(flost, topframe=lost_tk)
    lost_tk.mainloop()
