import os, sys, traceback
import numpy as np
from scipy.io import netcdf_file
try:
    import Tkinter as tk
    import ttk
    import tkMessageBox as tkmb
except:
    import tkinter as tk
    from tkinter import ttk
    from tkinter import messagebox as tkmb

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
try:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as nt2tk
except:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg as nt2tk


def read_birth(birth_file, fbm, topframe=None):

    print('Reading %s' %birth_file)

    if not os.path.isfile(birth_file):
        print('Error %s not found' %birth_file)
        return

    birthdir, birthfile  = os.path.split(birth_file)
    tmp = birthfile.split('_')
    runid = tmp[0]
    t_id  = birthfile.split('cdf')[-1]

# Get 
    ntheta = 101
    Rmaj = fbm.rsurf[0, 0]
    Rtor_in  = np.min(fbm.rsurf[np.nonzero(fbm.rsurf)])
    Rtor_out = fbm.rsurf.max()
    phi_tor = np.linspace(0, 2*np.pi, ntheta)
    cosp = np.cos(phi_tor)
    sinp = np.sin(phi_tor)
 
    cv = netcdf_file(birth_file, 'r', mmap=False).variables

    print('PLOT_BIRTH')
    print(birth_file)

    mcl = cv['mclabel'].data
    mc_label = "".join([x.decode('UTF-8') for x in mcl[0]]).strip()

# Read data from cdf

# r, z : deposition location for each particle
# time
# einj = injection energy
# xksid = pitch angle
# zeta = toroidal angle
# ib = #NBI source associated to each MC marker

    Rj      = cv['bs_rgc_%s'   %mc_label].data
    zj      = cv['bs_zgc_%s'   %mc_label].data
    Einj    = cv['bs_einj_%s'  %mc_label].data
    pitch   = cv['bs_xksid_%s' %mc_label].data
    weight  = cv['bs_wght_%s'  %mc_label].data
    tor_ang = cv['bs_zeta_%s'  %mc_label].data
    t_birth = cv['bs_time_%s'  %mc_label][0]
    j_nbi   = cv['bs_ib_%s'    %mc_label].data

    phi_dep = np.radians(tor_ang)
    xtop = Rj*np.cos(phi_dep)
    ytop = Rj*np.sin(phi_dep)
    n_birth = len(Rj)
    print('# of MC particles: %d' %n_birth)
    src_arr = np.unique(j_nbi)
    
    print('Sources: ', src_arr)
    
    # Vessel compoments for plot
    try:
        import plot_aug
        import aug_sfutils as sf
        xlin, ylin, rlin, zlin = plot_aug.nbi_plot(nbis=src_arr, runid=runid)
        gc_d = sf.getgc()
        tor_d = sf.getgc_tor(rotate=False)
        m2cm = 100.
        xpol_lim = (90, 230)
        ypol_lim = (-125, 125)
        xtop_lim = (-600, 400)
    except:
        print(traceback.format_exc())
        rlin = fbm.rlim_pts
        zlin = fbm.ylim_pts
        xext = ( fbm.rlim_pts.min() + fbm.rlim_pts.max())*0.025
        yext = (-fbm.ylim_pts.min() + fbm.ylim_pts.max())*0.025
        xtop_lim = np.array ([-fbm.rlim_pts.max() - xext, fbm.rlim_pts.max() + xext])
        xtop_lim = np.round(xtop_lim).astype(int)
        xpol_lim = np.array ([fbm.rlim_pts.min() - xext, fbm.rlim_pts.max() + xext])
        ypol_lim = np.array ([fbm.ylim_pts.min() - yext, fbm.ylim_pts.max() + yext])
        xpol_lim = np.round(xpol_lim).astype(int)
        ypol_lim = np.round(ypol_lim).astype(int)

    j_comp = np.zeros(n_birth, dtype=np.int32)
    ind_nbi = {}
    for jsrc in src_arr:
        (index, ) = np.where(j_nbi == jsrc)
        E = np.max(Einj[index])
        j_comp[index] = np.int32(np.max(E)/Einj[index] + 1e-3)
        ind_nbi[jsrc] = index
 
    comp_arr = np.unique(j_comp)
    print('Energy components', comp_arr)
    comp_lbl = {1: 'Full', 2: ' Half', 3: 'Third'}
    n_src  = len(src_arr)
    n_comp = len(comp_arr)

# Determine 1/2 and 1/3 fractions

    nr = 141
    nz = 101
    Rmin = Rtor_in
    Rmax = Rtor_out
    zmin = fbm.zsurf.min()
    zmax = fbm.zsurf.max()
    R_grid = np.linspace(Rmin, Rmax, nr)
    z_grid = np.linspace(zmin, zmax, nz)
    dep_matrix = {}

    for jsrc in src_arr:
        dep_matrix[jsrc] = {}
        for jcomp in comp_arr:
            ind = (j_nbi == jsrc) & (j_comp == jcomp)
            dep_matrix[jsrc][jcomp], Redge, zedge = \
            np.histogram2d(Rj[ind], zj[ind], bins=[nr, nz], \
            range=[[Rmin, Rmax], [zmin, zmax]], weights = weight[ind])

    res_R = {}
    for jsrc in src_arr:
        res_R[jsrc] = {}
        for jcomp in comp_arr:
            dep_R = np.sum(dep_matrix[jsrc][jcomp], axis=1) # z-sum
            res_R[jsrc][jcomp] = np.cumsum(dep_R)
            print('Deposited particles for source %d, component %s: %10.3e/s' %(jsrc, jcomp, res_R[jsrc][jcomp][-1]) )

# NBI geometry for plots
    print('RUNID = %s' %runid)

#------
# Plots
#------

    Rlbl = 'R [cm]'
    zlbl = 'z [cm]'
    fsize = 12
    colors = ('r', 'b', 'g', 'm', 'y')

    nrows = 2
    ncols = n_comp + 1

    n_pitch = 5
    pitch_edges = np.linspace(0, 1, n_pitch+1)

    xgrid, ygrid = np.meshgrid(R_grid, z_grid, indexing='ij')

# One tab for each source

    if topframe is None:
        topframe = tk.Toplevel()
        xmax = topframe.winfo_screenwidth()
        ymax = topframe.winfo_screenheight()
        width  = min(1500, int(0.95*xmax))
        height = min(940, int(0.88*ymax)) 
        topframe.title('Birth location')
        topframe.geometry('%dx%d' %(width, height))
    nbsource = ttk.Notebook(topframe, name='nb source')
    nbsource.pack(side=tk.TOP, fill=tk.X)

    for jnb, jsrc in enumerate(src_arr):
        frame_source = tk.Frame(nbsource)
        lbl = ' NBI #%d ' %jsrc

        nbbirth = ttk.Notebook(frame_source, name='nbbirth')
        nbbirth.pack(side=tk.TOP, fill=tk.X)

        frame_birth = tk.Frame(nbbirth)

        fig_birth = Figure(figsize=(14., 8.45), dpi=100)
        can_birth = FigureCanvasTkAgg(fig_birth, master=frame_birth)
        can_birth._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        fig_birth.clf()
        fig_birth.subplots_adjust(left=0.05, bottom=0.08, right=0.97, top=0.97,  \
                                  wspace=0.15, hspace=0.)

        fig_birth.text(0.5, 0.95, '%s, t =%6.3f s, top view' %(runid, t_birth), ha='center')
        fig_birth.text(0.5, 0.55, 'Poloidal section'  , ha='center')

        jsplot = 1

# Overplot 3 species mix
# Above view
        axtop = fig_birth.add_subplot(nrows, ncols, jsplot, aspect='equal')
        axtop.set_title('All energy components', fontsize=fsize)

# Poloidal section

        axpol = fig_birth.add_subplot(nrows, ncols, jsplot+ncols, aspect='equal')
        axpol.set_title('All energy components', fontsize=fsize)

# Birth locations

        for jcol, jcomp in enumerate(comp_arr):
            ind = (j_comp == jcomp) & (j_nbi == jsrc)
            axtop.plot(xtop[ind], ytop[ind], '%so' %colors[jcol], label=comp_lbl[jcomp])
            axpol.plot(Rj[ind]  , zj[ind]  , '%so' %colors[jcol], label=comp_lbl[jcomp])

        axtop.legend(loc=2, numpoints=1, prop={'size': 8})
        axpol.legend(loc=2, numpoints=1, prop={'size': 8})

# Tokamak components
  
        axpol.set_xlabel(Rlbl, fontsize=fsize)
        axpol.set_ylabel(zlbl, fontsize=fsize)
        axpol.set_xlim(xpol_lim)
        axpol.set_ylim(ypol_lim)
        axtop.set_xlim(xtop_lim)

        if 'gc_d' in locals():
            for gc in gc_d.values():
                axpol.plot(m2cm*gc.r, m2cm*gc.z, 'b-')
        if 'tor_d' in locals():
            for tor in tor_d.values():
                axtop.plot(m2cm*tor.x, m2cm*tor.y, 'b-')
        
        if 'fbm' in locals():
            axtop.plot(Rtor_in *cosp, Rtor_in *sinp, 'r-')
            axtop.plot(Rtor_out*cosp, Rtor_out*sinp, 'r-')
            axtop.plot(Rmaj*cosp, Rmaj*sinp, 'r--')

            for irho in range(fbm.rsurf.shape[0]):
                axpol.plot(fbm.rsurf[irho, :], fbm.zsurf[irho, :], 'r-', linewidth=0.1)
            for jbar, myr in enumerate(fbm.rbar):
                axpol.plot(myr, fbm.zbar[jbar], 'r-', linewidth=0.1)

        if 'rlin' in locals():
            axpol.plot(rlin, zlin, 'g-', linewidth=2.5)

# For each species overplot 5 pitch angle range

        jsplot = 2

        for jcomp in comp_arr:
            axtop = fig_birth.add_subplot(nrows, ncols, jsplot      , aspect='equal')
            axpol = fig_birth.add_subplot(nrows, ncols, jsplot+ncols, aspect='equal')

            axtop.set_title('%s energy' %jcomp, fontsize=fsize)
            axpol.set_title('%s energy' %jcomp, fontsize=fsize)

            axpol.set_xlabel(Rlbl, fontsize=fsize)
            axpol.set_ylabel(zlbl, fontsize=fsize)
            axpol.set_xlim(xpol_lim)
            axpol.set_ylim(ypol_lim)
            axtop.set_xlim(xtop_lim)

            ind1 = (j_comp == jcomp) & (j_nbi == jsrc)
            for jpitch in range(n_pitch):
                p1 = pitch_edges[jpitch]
                p2 = pitch_edges[jpitch+1]
                ind = ind1 & (pitch >  p1) & (pitch <= p2)
                axtop.plot(xtop[ind], ytop[ind], '%so' %colors[jpitch], \
                           label='%3.1f < p.a. < %3.1f' %(p1, p2))
                axpol.plot(Rj[ind], zj[ind], '%so' %colors[jpitch], \
                           label='%3.1f < p.a. < %3.1f' %(p1, p2))

            if 'gc_d' in locals():
                for gc in gc_d.values():
                    axpol.plot(m2cm*gc.r, m2cm*gc.z, 'b-')
            if 'tor_d' in locals():
                for tor in tor_d.values():
                    axtop.plot(m2cm*tor.x, m2cm*tor.y, 'b-')

            if 'fbm' in locals():
                axtop.plot(Rtor_in *cosp, Rtor_in *sinp, 'r-')
                axtop.plot(Rtor_out*cosp, Rtor_out*sinp, 'r-')
                axtop.plot(Rmaj*cosp, Rmaj*sinp, 'r--')

                for irho in range(fbm.rsurf.shape[0]):
                    axpol.plot(fbm.rsurf[irho, :], fbm.zsurf[irho, :], 'r-', linewidth=0.1)
                for jbar, myr in enumerate(fbm.rbar):
                    axpol.plot(myr, fbm.zbar[jbar], 'r-', linewidth=0.1)

            if 'rlin' in locals():
                axpol.plot(rlin, zlin, 'g-', linewidth=2.5)

            axtop.legend(loc=2, numpoints=1, prop={'size': 8})
            axpol.legend(loc=2, numpoints=1, prop={'size': 8})
            jsplot += 1

        toolbar = nt2tk(can_birth, frame_birth)
        toolbar.update()

#-------------------------
# Deposition & attenuation
#-------------------------

        frame_dep = tk.Frame(nbbirth)

        fig_pol = Figure(figsize=(6., 3.9), dpi=100)
        fig_pol.clf()
        frame_pol = tk.Frame(frame_dep, height=600)
        frame_pol.pack(side=tk.TOP, fill=tk.BOTH)
        frame_pol.pack_propagate(0)
        can_pol = FigureCanvasTkAgg(fig_pol, master=frame_pol)
        can_pol._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        fig_pol.subplots_adjust(left=0.05, bottom=0.1, right=0.98, top=0.92)

        fig_pol.text(0.33, 0.95, '%s, t =%6.3f s' %(runid, t_birth), ha='center')
# 2D deposition, poloidal section
        jsplot = 1
        for jcomp in comp_arr:
            zgrid = dep_matrix[jsrc][jcomp]
            ind = np.where(zgrid == 0)
            zgrid[ind] = None
            axpol = fig_pol.add_subplot(1, n_comp, jsplot, aspect='equal')

            axpol.set_title('%s energy' %comp_lbl[jcomp], fontsize=fsize)
            axpol.set_xlim(xpol_lim)
            axpol.set_ylim(ypol_lim)
            ctr = axpol.contourf(xgrid, ygrid, zgrid)
            fig_pol.colorbar(ctr, shrink=0.9, aspect=10)
            if 'gc_d' in locals():
                for gc in gc_d.values():
                    axpol.plot(m2cm*gc.r, m2cm*gc.z, 'b-')
            if 'fbm' in locals():
                for irho in range(fbm.rsurf.shape[0]):
                    axpol.plot(fbm.rsurf[irho, :], fbm.zsurf[irho, :], 'r-', linewidth=0.1)
                for jbar, myr in enumerate(fbm.rbar):
                    axpol.plot(myr, fbm.zbar[jbar], 'r-', linewidth=0.1)
            if 'rlin' in locals():
                axpol.plot(rlin, zlin, 'g-', linewidth=2.5)

            jsplot += 1

        toolbar = nt2tk(can_pol, frame_pol)
        toolbar.update()

# Attenutation

        fig_att = Figure(figsize=(3., 2.4), dpi=100)
        fig_att.clf()
        frame_att = tk.Frame(frame_dep)
        frame_att.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        frame_att.pack_propagate(0)
        can_att = FigureCanvasTkAgg(fig_att, master=frame_att)
        can_att._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        fig_att.subplots_adjust(left=0.05, bottom=0.2, right=0.98, \
                                top=0.9)

        jsplot = 1
        for jcomp in comp_arr:
            axatt = fig_att.add_subplot(1, n_comp, jsplot)
            axatt.set_title('%s energy' %comp_lbl[jcomp], fontsize=fsize)
            axatt.set_xlabel('R [cm]', fontsize=fsize)
            axatt.set_ylabel('NBI attenuation', fontsize=fsize)
            axatt.plot(R_grid, res_R[jsrc][jcomp])
            jsplot += 1

        toolbar = nt2tk(can_att, frame_att)
        toolbar.update()

# Add tabs

        nbbirth.add(frame_birth, text='Birth location')
        nbbirth.add(frame_dep  , text='NBI penetration')

        nbsource.add(frame_source, text=lbl)


if __name__ == "__main__":

    fbirth = '/afs/ipp/home/g/git/tr_client/AUGD/29783/A01/29783A01_birth.cdf1'
    birth_tk = tk.Tk()
    birth_tk.title('Birth location')
    birth_tk.geometry('1500x940')
    birth_tk.option_add("*Font", "Helvetica")
    read_birth(fbirth, topframe=birth_tk)
    birth_tk.mainloop()
