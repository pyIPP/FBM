from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication, QGridLayout, QMenu, QAction, QLabel, QPushButton, QLineEdit, QCheckBox, QSpinBox, QDoubleSpinBox, QFileDialog, QRadioButton, QButtonGroup, QTabWidget, QVBoxLayout, QHBoxLayout
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from scipy.interpolate import griddata
import numpy as np

try:
    import aug_sfutils as sf
    gc_d = sf.getgc()
    m2cm = 100.
    xpol_lim = (90, 230)
    ypol_lim = (-125, 125)
except:
    pass

fsize = 12
cellLinewidth = 0.5
limLinewidth = 2.5

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
    title += ', Run %s, t =%6.3f s' %(runid, fbm.time)
    fig.clf()
    ax = fig.add_subplot(111, aspect='equal')
    fig.text(0.5, 0.95, title, ha='center')
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
    Rpol_lim = np.array([Rmin - Rext, Rmax + Rext])
    Zpol_lim = np.array([Zmin - Zext, Zmax + Zext])
    Rpol_lim = np.round(Rpol_lim).astype(int)
    Zpol_lim = np.round(Zpol_lim).astype(int)
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
        figLost.subplots_adjust(left=0.05, bottom=0.08, right=0.97, top=0.97, wspace=0.15, hspace=0)
        figLost.text(0.5, 0.95, '%s Power losses for %s, t =%6.3f s ' %(spc_lbl, runid, tbm1), ha='center', fontsize=16)
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
