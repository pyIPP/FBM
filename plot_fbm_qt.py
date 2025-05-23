#!/usr/bin/env python

__author__  = 'Giovanni Tardini (Tel. +49 89 3299-1898)'
__version__ = '0.0.1'
__date__    = '29.05.2025'

import os, sys, logging, webbrowser
from scipy.io import netcdf_file
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from scipy.interpolate import griddata

try:
    from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication, QGridLayout, QMenu, QAction, QLabel, QPushButton, QLineEdit, QCheckBox, QSpinBox, QDoubleSpinBox, QFileDialog, QRadioButton, QButtonGroup, QTabWidget, QVBoxLayout, QHBoxLayout
    from PyQt5.QtGui import QPixmap, QIcon
    from PyQt5.QtCore import Qt, QRect, QSize, QLocale
    qt5 = True
except:
    from PyQt4.QtCore import Qt, QRect, QSize, QLocale
    from PyQt4.QtGui import QPixmap, QIcon, QMainWindow, QWidget, QApplication, QGridLayout, QMenu, QAction, QLabel, QPushButton, QLineEdit, QCheckBox, QSpinBox, QDoubleSpinBox, QFileDialog, QRadioButton, QButtonGroup, QTabWidget, QVBoxLayout, QHBoxLayout
    qt5 = False

try:
    import aug_sfutils as sf
    gc_d = sf.getgc()
    m2cm = 100.
    xpol_lim = (90, 230)
    ypol_lim = (-125, 125)
except:
    pass
import read_ac, plot_birth, config, plot_lost

os.environ['BROWSER'] = '/usr/bin/firefox'

usLocale = QLocale('us')

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger = logging.getLogger('DPSD_GUI')
logger.addHandler(hnd)
logger.setLevel(logging.INFO)

fbm_dir = os.path.dirname(os.path.realpath(__file__))

bdens_d = {'D_NBI': 'BDENS', 'H_NBI': 'BDENS', 'HE3_FUS': 'FDENS_3', \
           'H_FUS': 'FDENS_P', 'T_FUS': 'FDENS_T', 'HE4_FUS': 'FDENS_4'}

lframe_wid = 670
rframe_wid = 800
fsize = 12

def about():
    webbrowser.open('http://www.aug.ipp.mpg.de/aug/manuals/transp/fbm/plot_fbm.html')


class FBM(QMainWindow):


    def __init__(self):

        if sys.version_info[0] == 3:
            super().__init__()
        else:
            super(QMainWindow, self).__init__()

        self.setLocale(usLocale)

        xwin  = lframe_wid + rframe_wid
        yhead = 44
        yline = 30
        ybar  = 48
        ywin  = 900 + yhead + ybar

        qhead  = QWidget(self)
        qbar   = QWidget(self)
        qtabs  = QTabWidget(self)
        qhead.setGeometry(QRect(0,    0, xwin, yhead))
        qbar.setGeometry(QRect(0, yhead, xwin, ybar))
        qtabs.setGeometry(QRect(0, yhead+ybar, xwin, ywin-yhead-ybar))
        qtabs.setStyleSheet("QTabBar::tab { width: 120 }")
        header_grid = QGridLayout(qhead) 
        tbar_grid   = QGridLayout(qbar) 

#-----
# Tabs
#-----

# Distribution TAB
        qdist = QWidget()
        dist_layout = QHBoxLayout()
        qdist.setLayout(dist_layout)
        qtabs.addTab(qdist, '2D dist')
        dist_right_widget = QWidget()
        dist_left_widget  = QWidget()
        dist_left_layout  = QHBoxLayout()
        dist_right_layout = QVBoxLayout()
        dist_left_widget.setLayout(dist_left_layout)
        dist_right_widget.setLayout(dist_right_layout)

        figdist1 = Figure()
        self.distCanvas1 = FigureCanvas(figdist1)
        self.axdist1 = figdist1.add_subplot(111, aspect='equal')
# Plot vessel even before reading FBM file
        if 'gc_d' in globals():
            for gc in gc_d.values():
                self.axdist1.plot(m2cm*gc.r, m2cm*gc.z, 'b-')
        dist_left_layout.addWidget(self.distCanvas1)
        dist_layout.addWidget(dist_left_widget)
        dist_layout.addWidget(dist_right_widget)

# Trapped particle tab
        qtrap = QWidget()
        trap_layout = QHBoxLayout()
        qtrap.setLayout(trap_layout)
        qtabs.addTab(qtrap, 'Trap. Frac')
        trap_right_widget = QWidget()
        trap_left_widget  = QWidget()
        trap_left_layout  = QHBoxLayout()
        trap_right_layout = QVBoxLayout()
        trap_left_widget.setLayout(trap_left_layout)
        trap_right_widget.setLayout(trap_right_layout)

        figtrap1 = Figure()
        self.trapCanvas1 = FigureCanvas(figtrap1)
        self.axtrap1 = figtrap1.add_subplot(111, aspect='equal')
        trap_left_layout.addWidget(self.trapCanvas1)
        trap_layout.addWidget(trap_left_widget)
        trap_layout.addWidget(trap_right_widget)

# Neutron tab
        qneut = QWidget()
        neut_layout = QHBoxLayout()
        qneut.setLayout(neut_layout)
        qtabs.addTab(qneut, 'Neut. Frac')

        figneut = Figure()
        self.neutCanvas = FigureCanvas(figneut)
        self.axneut1 = figneut.add_subplot(121, aspect='equal')
        self.axneut2 = figneut.add_subplot(122, aspect='equal')
        neut_layout.addWidget(self.neutCanvas)

# Fast ions birth location
        qbirth = QWidget()
        birth_layout = QGridLayout()
        qbirth.setLayout(birth_layout)
        qtabs.addTab(qbirth, 'Birth profile')

# NBI losses
        qloss = QWidget()
        loss_layout = QGridLayout()
        qloss.setLayout(loss_layout)
        qtabs.addTab(qloss, 'Losses')

#--------
# Menubar
#--------

        menubar = self.menuBar()
        fileMenu = QMenu('&File', self)
        helpMenu = QMenu('&Help', self)
        menubar.addMenu(fileMenu)
        menubar.addMenu(helpMenu)

        runAction  = QAction('&Read FBM', fileMenu)
        exitAction = QAction('&Exit'    , fileMenu)
        runAction.triggered.connect(self.load_fbm)
        exitAction.triggered.connect(sys.exit)
        fileMenu.addAction(runAction)
        fileMenu.addSeparator()
        fileMenu.addAction(exitAction)

        aboutAction = QAction('&Web docu', helpMenu)
        aboutAction.triggered.connect(about)
        helpMenu.addAction(aboutAction)

        header_grid.addWidget(menubar, 0, 0, 1, 10)

        self.setStyleSheet("QLabel { width: 4 }")
        self.setStyleSheet("QLineEdit { width: 4 }")
        self.setGeometry(10, 10, xwin, ywin)
        self.setWindowTitle('DPSD')
        self.show()


    def load_fbm(self):

        dir_init = config.tr_clientDir + '/29795/A05'
        ffbm =  QFileDialog.getOpenFileName(self, 'Open file', \
            '%s/' %dir_init, "AC files (*.DATA*)")
        if qt5:
            f_ac = ffbm[0]
        else:
            f_ac = str(ffbm)
        print(f_ac)
        self.read_all(f_ac)


    def read_all(self, f_ac):

        fbmdir, fbmfile  = os.path.split(f_ac)

        tmp = fbmfile.split('.')
        runid = tmp[0]
        t_id = tmp[1][4:]
        self.fbmr = read_ac.READ_FBM(f_ac)

# Read FBM 

        self.species = []
        for spc_lbl in self.fbmr.int_en_pit_frac_trap.keys():
            self.species.append(spc_lbl)

# Read CDF
        tr_file = '%s/%s.CDF' %(fbmdir, runid)
        cv_all = netcdf_file(tr_file, 'r', mmap=False).variables
        bdlist = [x for x in bdens_d.values()]
        sigs = ['BTNT2_DD', 'BBNT2_DD', 'TIME3', 'X'] + bdlist
        tdist = (cv_all['TIME3'].data - self.fbmr.time)**2
        jtclose = np.argmin(tdist)
        self.cv = {}
        for lbl in sigs:
            if lbl in cv_all.keys():
                self.cv[lbl] = cv_all[lbl][jtclose]

        btneut = self.cv['BTNT2_DD']
        bbneut = self.cv['BBNT2_DD']
        if (np.max(btneut) == 0) and (np.max(bbneut) == 0):
            print('Zero neutrons')
            return
        indbt = (btneut == 0)
        indbb = (btneut == 0)
        btneut[indbt] = np.nan
        bbneut[indbb] = np.nan

#--------
# Plots
#--------

        r_grid, z_grid = np.meshgrid(
        np.linspace(self.fbmr.r2d.min(), self.fbmr.r2d.max(), 100),
        np.linspace(self.fbmr.z2d.min(), self.fbmr.z2d.max(), 100))

# FBM integrated over E, mu
        f_in = self.fbmr.int_en_pit[spc_lbl]
        self.contourPlot(self.axdist1, r_grid, z_grid, f_in)

# Trapped particle fraction
        f_in = self.fbmr.int_en_pit_frac_trap[spc_lbl]
        self.contourPlot(self.axtrap1, r_grid, z_grid, f_in)

# Neutron
        f_in = btneut
        self.contourPlot(self.axneut1, r_grid, z_grid, f_in)
        f_in = bbneut
        self.contourPlot(self.axneut2, r_grid, z_grid, f_in)


    def contourPlot(self, ax, r_grid, z_grid, f_in):
        if 'gc_d' in globals():
            for gc in gc_d.values():
                ax.plot(m2cm*gc.r, m2cm*gc.z, 'b-')
        for irho in range(self.fbmr.r_surf.shape[0]):
            ax.plot(self.fbmr.r_surf[irho, :], self.fbmr.z_surf[irho, :], 'r-', linewidth=0.5)
        for jbar, myr in enumerate(self.fbmr.rbar):
            ax.plot(myr, self.fbmr.zbar[jbar], 'r-')
        f_grid = griddata((self.fbmr.r2d, self.fbmr.z2d), f_in, (r_grid, z_grid), method='cubic')
        ax.contourf(r_grid, z_grid, f_grid, levels=50, cmap='viridis')
        ax.set_xlabel('R [cm]')
        ax.set_ylabel('Z [cm]')
        canvas = ax.figure.canvas
        canvas.draw()


if __name__ == '__main__':


    app = QApplication(sys.argv)
    main = FBM()
    app.exec_()
