from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication, QGridLayout, QMenu, QAction, QLabel, QPushButton, QLineEdit, QCheckBox, QSpinBox, QDoubleSpinBox, QFileDialog, QRadioButton, QButtonGroup, QTabWidget, QVBoxLayout, QHBoxLayout
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from scipy.interpolate import griddata
import numpy as np


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

    fbmfile  = fbm.fileName.split('/')[-1]
    runid    = fbmfile[:8]
    title += ', Run %s, t =%6.3f s' %(runid, fbm.time)
    fig.clf()
    ax = fig.add_subplot(111, aspect='equal')
    fig.text(0.5, 0.95, title, ha='center')
    f_grid = griddata((fbm.r2d, fbm.z2d), f_in, (r_grid, z_grid), method='cubic')
    plot2d = ax.contourf(r_grid, z_grid, f_grid, levels=50, cmap='viridis')
    fig.colorbar(plot2d, aspect=20)
    ax.set_xlabel('R [cm]')
    ax.set_ylabel('Z [cm]')
    addGrids(fbm, ax)
    canvas = fig.canvas
    canvas.draw()

    return ax


def addGrids(fbm, ax):
    
    if 'gc_d' in globals():
        for gc in gc_d.values():
            ax.plot(m2cm*gc.r, m2cm*gc.z, 'b-')
    for irho in range(fbm.r_surf.shape[0]):
        ax.plot(fbm.r_surf[irho, :], fbm.z_surf[irho, :], 'r-', linewidth=0.5)
    for jbar, myr in enumerate(fbm.rbar):
        ax.plot(myr, fbm.zbar[jbar], 'r-')
# Draw limiter, set plot boundaries
    xext =  (fbm.rlim_pts.min() + fbm.rlim_pts.max())*0.025
    yext = (-fbm.ylim_pts.min() + fbm.ylim_pts.max())*0.025
    xpol_lim = np.array ([fbm.rlim_pts.min() - xext, fbm.rlim_pts.max() + xext])
    ypol_lim = np.array ([fbm.ylim_pts.min() - yext, fbm.ylim_pts.max() + yext])
    xpol_lim = np.round(xpol_lim).astype(int)
    ypol_lim = np.round(ypol_lim).astype(int)
    rlin = fbm.rlim_pts
    zlin = fbm.ylim_pts
    ax.set_xlim(xpol_lim)
    ax.set_ylim(ypol_lim)
    ax.plot(rlin, zlin, 'g-', linewidth=2.5)

    return ax


def plotTrapped(fbm, r_grid, z_grid, trapLayout):

    clear_layout(trapLayout)

    trapTabs = QTabWidget()
    trapTabs.setStyleSheet("QTabBar::tab { width: 120 }")
    trapLayout.addWidget(trapTabs)
    for spc_lbl, f_in in fbm.int_en_pit_frac_trap.items():
        qspec = QWidget()
        trapTabs.addTab(qspec, spc_lbl)
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
 
