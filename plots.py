import os, logging
from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication, QGridLayout, QMenu, QAction, QLabel, QPushButton, QLineEdit, QCheckBox, QSpinBox, QDoubleSpinBox, QFileDialog, QRadioButton, QButtonGroup, QTabWidget, QVBoxLayout, QHBoxLayout
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from scipy.interpolate import griddata
from scipy.io import netcdf_file
import numpy as np

try:
    import aug_sfutils as sf
    gc_d = sf.getgc()
    m2cm = 100.
    xpol_lim = (90, 230)
    ypol_lim = (-125, 125)
    tor_d = sf.getgc_tor(rotate=False)
except:
    pass

fsize = 12
cellLinewidth = 0.5
limLinewidth = 2.5
titSize = 14

def clear_layout(layout):

    while layout.count():
        item = layout.takeAt(0)
        widget = item.widget()
        if widget is not None:
            widget.setParent(None)
            widget.deleteLater()
        elif item.layout():  # In case it's a nested layout
            clear_layout(item.layout())


def contourPlotRZ(fbm, fig, r_grid, z_grid, f_in, title=''):

    fbmfile = fbm.fileName.split('/')[-1]
    runid   = fbmfile[:8]
    title += ', run %s, t=%6.3f s' %(runid, fbm.time)
    fig.clf()
    ax = fig.add_subplot(111, aspect='equal')
    fig.text(0.5, 0.95, title, ha='center', fontsize=titSize)
    f_grid = griddata((fbm.r2d, fbm.z2d), f_in, (r_grid, z_grid), method='cubic')
    plot2d = ax.contourf(r_grid, z_grid, f_grid, levels=50, cmap='viridis')
    fig.colorbar(plot2d, aspect=20)
    addGrids(fbm, ax)
    canvas = fig.canvas
    canvas.draw()

    return ax


def addGrids(fbm, ax):
    
    if 'gc_d' in globals():
        for gc in gc_d.values():
            ax.plot(m2cm*gc.r, m2cm*gc.z, 'b-')
    for irho in range(fbm.r_surf.shape[0]):
        ax.plot(fbm.r_surf[irho, :], fbm.z_surf[irho, :], 'r-', linewidth=cellLinewidth)
    for jbar, myr in enumerate(fbm.rbar):
        ax.plot(myr, fbm.zbar[jbar], 'r-')
# Draw limiter, set plot boundaries
    Rmax = fbm.rlim_pts.max()
    Rmin = fbm.rlim_pts.min()
    Zmax = fbm.ylim_pts.max()
    Zmin = fbm.ylim_pts.min()
    Rext = 0.025*(Rmax - Rmin)
    Zext = 0.025*(Zmax - Zmin)
    Rpol_lim = [int(Rmin - Rext), int(Rmax + Rext)]
    Zpol_lim = [int(Zmin - Zext), int(Zmax + Zext)]
    Rlin = fbm.rlim_pts
    Zlin = fbm.ylim_pts
    ax.set_xlabel('R [cm]')
    ax.set_ylabel('Z [cm]')
    ax.set_xlim(Rpol_lim)
    ax.set_ylim(Zpol_lim)
    ax.plot(Rlin, Zlin, 'g-', linewidth=limLinewidth)

    return ax


def plotTrapped(fbm, r_grid, z_grid, trapLayout):

    clear_layout(trapLayout)

    specTabs = QTabWidget()
    specTabs.setStyleSheet("QTabBar::tab { width: 120 }")
    trapLayout.addWidget(specTabs)
    for spc_lbl, f_in in fbm.int_en_pit_frac_trap.items():
        qspec = QWidget()
        specTabs.addTab(qspec, spc_lbl)
        tabLayout = QHBoxLayout()
        qspec.setLayout(tabLayout)
        trap_right_widget = QWidget()
        trap_left_widget  = QWidget()
        trap_left_layout  = QHBoxLayout()
        trap_right_layout = QVBoxLayout()
        trap_left_widget.setLayout(trap_left_layout)
        trap_right_widget.setLayout(trap_right_layout)
        figTrap = Figure()
        canvasTrap = FigureCanvas(figTrap)
        axtrap = figTrap.add_subplot(111, aspect='equal')
        trap_left_layout.addWidget(canvasTrap)
        tabLayout.addWidget(trap_left_widget)
        tabLayout.addWidget(trap_right_widget)
        contourPlotRZ(fbm, figTrap, r_grid, z_grid, f_in, title='Trapped fast ion fraction')
        canvasTrap.draw()


def plotLost(fbm, r_grid, z_grid, lossLayout):

    clear_layout(lossLayout)

    lossTabs = QTabWidget()
    lossTabs.setStyleSheet("QTabBar::tab { width: 120 }")
    lossLayout.addWidget(lossTabs)

    colors = ('b','c', 'g', 'y', 'r')
    sizes = (7, 8, 9, 12, 14)
    nrows = 1
    n_pwr = 5
    comp_lbl = ['P', 'NP', 'Total']
    tbm1 = fbm.ta
    tbm2 = fbm.tb
    fbmfile = fbm.fileName.split('/')[-1]
    runid   = fbmfile[:8]
    for spc_lbl, weight in fbm.wghtlost.items():
        qspec = QWidget()
        lossTabs.addTab(qspec, spc_lbl)
        tabLayout = QHBoxLayout()
        qspec.setLayout(tabLayout)
        figLost = Figure()
        canvasLost = FigureCanvas(figLost)
        figLost.subplots_adjust(left=0.05, bottom=0.08, right=0.99, top=0.97, wspace=0.2, hspace=0)
        figLost.text(0.5, 0.95, '%s Power losses for %s, t =%6.3f s ' %(spc_lbl, runid, tbm1), ha='center', fontsize=titSize)
        tabLayout.addWidget(canvasLost)

        Rj = np.array(fbm.rmjionlost[spc_lbl])
        zj = np.array(fbm.yyionlost[spc_lbl])
        pwr = np.array(fbm.xksidlost[spc_lbl])
        Elost = np.array(fbm.yzelost[spc_lbl])
        j_nbi = np.array(fbm.ynbscelost[spc_lbl]).astype(int)
        n_lost = len(Rj)
        print('# of MC particles lost: %d' %n_lost)    
        if n_lost == 0:
            continue
        j_comp     = np.zeros(n_lost, dtype=np.int32)
        j_comp_tot = np.zeros(n_lost, dtype=np.int32)
    
        (indx_p , ) = np.where(fbm.gooselost[spc_lbl] < 1.000001)
        (indx_np, ) = np.where(fbm.gooselost[spc_lbl] >= 1.000001)
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
        print('%s Total power lost, [Watts] ' %spc_lbl, pwr_lost_tot)
        pwr_edges = np.linspace(np.log10(pwr_min), np.log10(pwr_max), n_pwr+1)

# Left plot: lost locations
        jsplot = 1     
        axpol = figLost.add_subplot(nrows, ncols, jsplot, aspect='equal')
        axpol.set_title('%s total losses' %spc_lbl, fontsize=fsize, fontweight='bold')
        for jcol, jcomp in enumerate(comp_arr):
            ind = (j_comp == jcomp)
            axpol.plot(Rj[ind], zj[ind], '%so' %colors[jcol], label=comp_lbl[jcol], markersize=sizes[jcol])
        axpol.legend(loc=2, numpoints=1, prop={'size': 10})
        addGrids(fbm, axpol)

# For each species overplot 5 pwr angle range
        jsplot = 2
        for jcomp in comp_arr:
            axpol = figLost.add_subplot(nrows, ncols, jsplot, aspect='equal')
            axpol.set_title('%s Pwr. losses, Watts' %comp_lbl[jcomp], \
                            fontsize=fsize,fontweight='bold')
            addGrids(fbm, axpol)
            ind1 = (j_comp == jcomp)
            for jpwr in range(n_pwr-1, -1, -1):
                p1 = pwr_edges[jpwr]
                p2 = pwr_edges[jpwr+1]
                ind = ind1 & (pwr_lost>  np.float_power(10,p1)) & (pwr_lost <= np.float_power(10,p2))
                axpol.plot(Rj[ind], zj[ind], '%so' %colors[jpwr], \
                           label='%s < Pwr < %s' %(round(np.float_power(10, p1), 2), \
                                                   round(np.float_power(10, p2), 2)), markersize=sizes[jpwr])
            axpol.legend(loc=2, numpoints=1, prop={'size': 10})
            jsplot += 1


def plotBirth(fbm, r_grid, z_grid, birthLayout):

    print(r_grid.shape)
    nZ, nR = r_grid.shape
    fbmdir, fbmfile  = os.path.split(fbm.fileName)
    tmp = fbmfile.split('.')
    runid = tmp[0]
    t_id = tmp[1][4:]
    birth_file =  '%s/%s_birth.cdf%s' %(fbmdir, runid, t_id)

    if not os.path.isfile(birth_file):
        print('Error %s not found' %birth_file)
        return
    
    ntheta = 101
    Rmaj = fbm.rsurf[0, 0]
    Rtor_in  = np.min(fbm.rsurf[np.nonzero(fbm.rsurf)])
    Rtor_out = fbm.rsurf.max()
    phi_tor = np.linspace(0, 2*np.pi, ntheta)
    cosp = np.cos(phi_tor)
    sinp = np.sin(phi_tor)
    Rmin = Rtor_in
    Rmax = Rtor_out
    zmin = fbm.zsurf.min()
    zmax = fbm.zsurf.max()

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

    dep_matrix = {}
    for jsrc in src_arr:
        dep_matrix[jsrc] = {}
        for jcomp in comp_arr:
            ind = (j_nbi == jsrc) & (j_comp == jcomp)
            dep_matrix[jsrc][jcomp], Redge, zedge = \
            np.histogram2d(Rj[ind], zj[ind], bins=[nR, nZ], \
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

    clear_layout(birthLayout)
    sourceTabs = QTabWidget()
    sourceTabs.setStyleSheet("QTabBar::tab { width: 120 }")
    birthLayout.addWidget(sourceTabs)
    for jnb, jsrc in enumerate(src_arr):
        lbl = ' NBI #%d ' %jsrc
        qsource = QWidget()
        sourceTabs.addTab(qsource, lbl)
        tabLayout = QHBoxLayout()
        qsource.setLayout(tabLayout)
        figBirth = Figure()
        canvasBirth = FigureCanvas(figBirth)
        figBirth.subplots_adjust(left=0.05, bottom=0.08, right=0.97, top=0.97,  \
                                  wspace=0.15, hspace=0.)
        figBirth.text(0.5, 0.95, '%s, t =%6.3f s, top view' %(runid, t_birth), ha='center')
        figBirth.text(0.5, 0.55, 'Poloidal section'  , ha='center')

        jsplot = 1
# Overplot 3 species mix
# Above view
        axtop = figBirth.add_subplot(nrows, ncols, jsplot, aspect='equal')
        axtop.set_title('All energy components', fontsize=fsize)

# Poloidal section
        axpol = figBirth.add_subplot(nrows, ncols, jsplot+ncols, aspect='equal')
        axpol.set_title('All energy components', fontsize=fsize)

# Birth locations
        for jcol, jcomp in enumerate(comp_arr):
            ind = (j_comp == jcomp) & (j_nbi == jsrc)
            axtop.plot(xtop[ind], ytop[ind], '%so' %colors[jcol], label=comp_lbl[jcomp])
            axpol.plot(Rj[ind]  , zj[ind]  , '%so' %colors[jcol], label=comp_lbl[jcomp])

        axtop.legend(loc=2, numpoints=1, prop={'size': 8})
        axpol.legend(loc=2, numpoints=1, prop={'size': 8})

        addGrids(fbm, axpol)
        if 'tor_d' in globals():
            for tor in tor_d.values():
                axtop.plot(m2cm*tor.x, m2cm*tor.y, 'b-')
        axtop.plot(Rtor_in *cosp, Rtor_in *sinp, 'r-')
        axtop.plot(Rtor_out*cosp, Rtor_out*sinp, 'r-')
        axtop.plot(Rmaj*cosp, Rmaj*sinp, 'r--')

        for irho in range(fbm.rsurf.shape[0]):
            axpol.plot(fbm.rsurf[irho, :], fbm.zsurf[irho, :], 'r-', linewidth=0.1)
        for jbar, myr in enumerate(fbm.rbar):
            axpol.plot(myr, fbm.zbar[jbar], 'r-', linewidth=0.1)

# For each species overplot 5 pitch angle range

        jsplot = 2

        for jcomp in comp_arr:
            axtop = figBirth.add_subplot(nrows, ncols, jsplot      , aspect='equal')
            axpol = figBirth.add_subplot(nrows, ncols, jsplot+ncols, aspect='equal')

            axtop.set_title('%s energy' %jcomp, fontsize=fsize)
            axpol.set_title('%s energy' %jcomp, fontsize=fsize)

            axpol.set_xlabel(Rlbl, fontsize=fsize)
            axpol.set_ylabel(zlbl, fontsize=fsize)
            axpol.set_xlim(xpol_lim)
            axpol.set_ylim(ypol_lim)

            ind1 = (j_comp == jcomp) & (j_nbi == jsrc)
            for jpitch in range(n_pitch):
                p1 = pitch_edges[jpitch]
                p2 = pitch_edges[jpitch+1]
                ind = ind1 & (pitch >  p1) & (pitch <= p2)
                axtop.plot(xtop[ind], ytop[ind], '%so' %colors[jpitch], \
                           label='%3.1f < p.a. < %3.1f' %(p1, p2))
                axpol.plot(Rj[ind], zj[ind], '%so' %colors[jpitch], \
                           label='%3.1f < p.a. < %3.1f' %(p1, p2))

            axtop.legend(loc=2, numpoints=1, prop={'size': 8})
            axpol.legend(loc=2, numpoints=1, prop={'size': 8})
            addGrids(fbm, axpol)
            if 'tor_d' in globals():
                for tor in tor_d.values():
                    axtop.plot(m2cm*tor.x, m2cm*tor.y, 'b-')
            axtop.plot(Rtor_in *cosp, Rtor_in *sinp, 'r-')
            axtop.plot(Rtor_out*cosp, Rtor_out*sinp, 'r-')
            axtop.plot(Rmaj*cosp, Rmaj*sinp, 'r--')
            jsplot += 1

        toolbar = NavigationToolbar(canvasBirth)
        tabLayout.addWidget(canvasBirth)
        canvasBirth.draw()
