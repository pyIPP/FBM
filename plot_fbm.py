#!/usr/bin/env python

__author__  = 'Giovanni Tardini (Tel. +49 89 3299-1898)'
__version__ = '0.0.2'
__date__    = '29.05.2025'

import os, sys, logging, webbrowser
from scipy.io import netcdf_file
import numpy as np
import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from PyQt5.QtWidgets import QMainWindow, QWidget, QApplication, QGridLayout, QMenu, QAction, QLabel, QCheckBox, QDoubleSpinBox, QFileDialog, QRadioButton, QButtonGroup, QTabWidget, QVBoxLayout, QHBoxLayout
from PyQt5.QtCore import Qt, QRect, QLocale
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavTool

try:
    import aug_sfutils as sf
    gc_d = sf.getgc()
    m2cm = 100.
except:
    pass
import read_ac, config
from plots import contourPlotRZ, plotTrapped, clear_layout, plotLost, plotBirth, plotNeutrons

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
logger = logging.getLogger('FBM')
if len(logger.handlers) == 0:
    hnd = logging.StreamHandler()
    hnd.setFormatter(fmt)
    logger.addHandler(hnd)

#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

logger.info('Using version %s', __version__)

os.environ['BROWSER'] = '/usr/bin/firefox'

usLocale = QLocale('us')

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger = logging.getLogger('FBM_GUI')
logger.addHandler(hnd)
logger.setLevel(logging.INFO)

fbm_dir = os.path.dirname(os.path.realpath(__file__))

bdens_d = {'D_NBI': 'BDENS', 'H_NBI': 'BDENS', 'HE3_FUS': 'FDENS_3', \
           'H_FUS': 'FDENS_P', 'T_FUS': 'FDENS_T', 'HE4_FUS': 'FDENS_4'}
rblist = ['Single cell', 'Theta averaged', 'Volume averaged']

yhead = 42
widgetHeight = 36
textWidth = 100
fsize = 10
fac = 0.58
ywin = int(fac*1150)
xwin = int(fac*1500)
distRightCanHeight = int(fac*500)
canCellHeight = int(fac*450)


def about():
    webbrowser.open('http://www.aug.ipp.mpg.de/aug/manuals/transp/fbm/plot_fbm.html')


class FBM(QMainWindow):


    def __init__(self):

        if sys.version_info[0] == 3:
            super().__init__()
        else:
            super(QMainWindow, self).__init__()

        self.setLocale(usLocale)

        qhead  = QWidget(self)
        qtabs  = QTabWidget(self)
        qhead.setGeometry(QRect(0,     0, xwin, yhead))
        qtabs.setGeometry(QRect(0, yhead, xwin, ywin-yhead))
        qtabs.setStyleSheet("QTabBar::tab { width: 120 }")
        header_grid = QGridLayout(qhead) 

#------#
# Tabs #
#------#

#-----------------
# Distribution TAB
        qdist = QWidget()
        self.distLayout = QHBoxLayout()
        qdist.setLayout(self.distLayout)
        qtabs.addTab(qdist, '2D dist')

# Initial plot: vessel components
        self.plotDistribution(self.distLayout)

#-----------------
# Trapped particle tab
        qtrap = QWidget()
        self.trapLayout = QHBoxLayout()
        qtrap.setLayout(self.trapLayout)
        qtabs.addTab(qtrap, 'Trap. Frac')

#-----------------
# Neutron tab
        qneut = QWidget()
        self.neutLayout = QHBoxLayout()
        qneut.setLayout(self.neutLayout)
        qtabs.addTab(qneut, 'Neutrons')

#-----------------
# Fast ions birth location
        qbirth = QWidget()
        self.birthLayout = QGridLayout()
        qbirth.setLayout(self.birthLayout)
        qtabs.addTab(qbirth, 'Birth profile')

#-----------------
# NBI losses
        qloss = QWidget()
        self.lossLayout = QGridLayout()
        qloss.setLayout(self.lossLayout)
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

        self.setGeometry(10, 10, xwin, ywin)
        self.setWindowTitle('FBM viewer')
        self.show()


    def load_fbm(self):

        dir_init = config.tr_clientDir
        if os.getenv('USER') == 'git': # For Giovanni Tardini's quick testing
            dir_init += '/29795/A05'
        ffbm =  QFileDialog.getOpenFileName(self, 'Open file', dir_init, "AC files (*.DATA*)")
        f_ac = ffbm[0]
        self.read_all(f_ac)


    def read_all(self, f_ac):

        fbmdir, fbmfile  = os.path.split(f_ac)

        tmp = fbmfile.split('.')
        runid = tmp[0]
        t_id = tmp[1][4:]

# Read FBM
        self.fbmr = read_ac.READ_FBM(f_ac)

# Read CDF
        tr_file = '%s/%s.CDF' %(fbmdir, runid)
        cv_all = netcdf_file(tr_file, 'r', mmap=False).variables
        bdlist = [x for x in bdens_d.values()]
        sigs = ['BTNT2_DD', 'BBNT2_DD', 'TIME3', 'X'] + bdlist
        tdist = (cv_all['TIME3'].data - self.fbmr.time)**2
        jtclose = np.argmin(tdist)
        self.cv = {}
        for lbl in sigs:
            if lbl in cv_all:
                self.cv[lbl] = cv_all[lbl][jtclose]

#--------
# Plots
#--------

        nR = 100
        nZ = 130
        self.r_grid, self.z_grid = np.meshgrid(
        np.linspace(self.fbmr.r2d.min(), self.fbmr.r2d.max(), nR),
        np.linspace(self.fbmr.z2d.min(), self.fbmr.z2d.max(), nZ))

        self.plotDistribution(self.distLayout)

# Trapped particle fraction
        plotTrapped(self.fbmr, self.r_grid, self.z_grid, self.trapLayout)

# Neutron
        plotNeutrons(self.fbmr, self.cv, self.r_grid, self.z_grid, self.neutLayout)

# Lost particles and power
        plotLost(self.fbmr, self.r_grid, self.z_grid, self.lossLayout)

# Particles birth
        plotBirth(self.fbmr, self.r_grid, self.z_grid, self.birthLayout)


    def my_call(self, event):

        if event.button in (2, 3):
            dist2 = (self.fbmr.r2d - event.xdata)**2 + (self.fbmr.z2d - event.ydata)**2
            jcell = np.argmin(dist2)
            self.plot_fbm_cell(jcell)


    def plotDistribution(self, distLayout):

        clear_layout(distLayout)

        distTabs = QTabWidget()
        distTabs.setStyleSheet("QTabBar::tab { width: 120 }")
        distLayout.addWidget(distTabs)
        self.figDist = {}
        if hasattr(self, 'fbmr'):
            mydic = self.fbmr.int_en_pit
        else:
            mydic = {'D_NBI': None}
        for spc_lbl, f_in in mydic.items():
            qspec = QWidget()
            distTabs.addTab(qspec, spc_lbl)
            tabLayout = QHBoxLayout()
            qspec.setLayout(tabLayout)
            dist_left_widget  = QWidget()
            dist_right_widget = QWidget()
            dist_left_layout  = QVBoxLayout()
            dist_right_layout = QVBoxLayout()
            dist_right1_widget = QWidget()
            dist_right2_widget = QWidget()
            dist_right3_widget = QWidget()
            dist_left_widget.setLayout(dist_left_layout)
            dist_right_widget.setLayout(dist_right_layout)
            self.figDist[spc_lbl] = Figure()
            self.figDist[spc_lbl].subplots_adjust(left=0.05, bottom=0.1, right=0.97, top=0.9)
            canvasDist = FigureCanvas(self.figDist[spc_lbl])
            toolbar = NavTool(canvasDist)
            dist_left_layout.addWidget(canvasDist)
            dist_left_layout.addWidget(toolbar, alignment=Qt.AlignBottom)
            tabLayout.addWidget(dist_left_widget, 45)
            tabLayout.addWidget(dist_right_widget, 55)
# Right column
            dist_right_layout.addWidget(dist_right1_widget, 14)
            dist_right_layout.addWidget(dist_right2_widget, 43)
            dist_right_layout.addWidget(dist_right3_widget, 43)

# Buttons in dist_right1_widget
            dist_right1_layout = QVBoxLayout()
            dist_right2_layout = QVBoxLayout()
            dist_right3_layout = QVBoxLayout()
            dist_right1_widget.setLayout(dist_right1_layout)
            dist_right2_widget.setLayout(dist_right2_layout)
            dist_right3_widget.setLayout(dist_right3_layout)

            dist_right1_widget.setFixedHeight(4*widgetHeight + 5)
            dist_right2_widget.setFixedHeight(distRightCanHeight)
            mouseLbl = QLabel('Right-mouse click on a cell for local FBM plot')
            buttons_widget = QWidget()
            buttons_layout = QHBoxLayout()
            buttons_widget.setLayout(buttons_layout)
            self.butt_d = {}
            self.butt_d['integral'] = QButtonGroup()
            for jcol, val in enumerate(rblist):
                but = QRadioButton(val)
                if jcol == 0:
                    but.setChecked(True)
                self.butt_d['integral'].addButton(but)
                self.butt_d['integral'].setId(but, jcol)
                buttons_layout.addWidget(but)
            buttons_layout.addStretch()
            energy_widget = QWidget()
            energy_layout = QHBoxLayout()
            Elbl = QLabel('Emax [keV]')
            self.butt_d['Emax'] = QDoubleSpinBox()
            self.butt_d['Emax'].setValue(100.)
            self.butt_d['Emax'].setDecimals(1)
            self.butt_d['Emax'].setRange(0., 1000.)
            energy_widget.setLayout(energy_layout)
            energy_layout.addWidget(Elbl)
            energy_layout.addWidget(self.butt_d['Emax'], alignment=Qt.AlignLeft)
            energy_layout.addStretch()

            log_widget = QWidget()
            log_layout = QHBoxLayout()
            log_widget.setLayout(log_layout)
            self.butt_d['logScale'] = QCheckBox('Log scale')
            self.butt_d['logScale'].setChecked(True)
            logMinLbl = QLabel('    f_log_min')
            logMaxLbl = QLabel('    f_log_max')
            self.butt_d['logMin'] = QDoubleSpinBox()
            self.butt_d['logMax'] = QDoubleSpinBox()
            self.butt_d['logMin'].setValue(4.5)
            self.butt_d['logMax'].setValue(8.5)
            self.butt_d['logMin'].setDecimals(2)
            self.butt_d['logMax'].setDecimals(2)
            self.butt_d['logMin'].setRange(-1., 20.)
            self.butt_d['logMax'].setRange(-1., 20.)
            for lbl in ('Emax', 'logMin', 'logMax'):
                self.butt_d[lbl].setFixedWidth(textWidth)
            log_layout.addWidget(self.butt_d['logScale'])
            log_layout.addWidget(logMinLbl)
            log_layout.addWidget(self.butt_d['logMin'])
            log_layout.addWidget(logMaxLbl)
            log_layout.addWidget(self.butt_d['logMax'])
            log_layout.addStretch()

            for layout in dist_left_layout, dist_right_layout, dist_right1_layout, dist_right2_layout, dist_right3_layout, log_layout, energy_layout:
                layout.setContentsMargins(0, 0, 0, 0)
            for layout in dist_right1_layout, log_layout, energy_layout:
                layout.setSpacing(0)
# BDENS canvas in dist_right2_widget

            figBdens = Figure()
            canvasBdens = FigureCanvas(figBdens)

# Bdens
            axBdens = figBdens.add_subplot(111)
            figBdens.subplots_adjust(left=0.1, bottom=0.3, right=0.97, top=0.9)
            if f_in is not None:
                axBdens.plot(self.fbmr.rho_grid, self.fbmr.bdens[spc_lbl], 'r-', label='From FBM', linewidth=2.5)
                axBdens.plot(self.cv['X'], self.cv[bdens_d[spc_lbl]], 'g-', label='From CDF', linewidth=2.5)
                axBdens.set_ylabel(r'%s [1/cm$^3$]' %bdens_d[spc_lbl], fontsize=fsize)
            axBdens.set_xlabel(r'$\rho_{tor}$', fontsize=fsize)
            axBdens.set_xlim([0, 1])
            axBdens.set_ylim(ymin = 0)
            axBdens.legend()

# FBM integrated over E, mu
            if f_in is None:
                if 'gc_d' in globals():
                    ax = self.figDist[spc_lbl].add_subplot(111, aspect='equal')
                    for gc in gc_d.values():
                        ax.plot(m2cm*gc.r, m2cm*gc.z, 'b-')
            else:
                ax = contourPlotRZ(self.fbmr, self.figDist[spc_lbl], self.r_grid, self.z_grid, f_in, title=r'$\int$ f($\rho$, $\theta$, E, $\mu$) dE d$\mu$')
                self.cell_mark, = ax.plot([], [], 'ro')
                self.currentSpecies = spc_lbl
                self.figDist[spc_lbl].canvas.mpl_connect('button_press_event', self.my_call)

# Cell canvas in dist_right2_widget

            self.figCell = Figure()
            canvasCell = FigureCanvas(self.figCell)
            toolbar = NavTool(canvasCell)

            for wid in mouseLbl, buttons_widget, energy_widget, log_widget:
                dist_right1_layout.addWidget(wid)
                wid.setFixedHeight(widgetHeight)
            dist_right2_layout.addWidget(canvasBdens)
            dist_right3_layout.addWidget(canvasCell)
            dist_right3_layout.addWidget(toolbar, alignment=Qt.AlignBottom)
            canvasCell.setFixedHeight(canCellHeight)
            if f_in is not None:
                self.plot_fbm_cell(0)


    def plot_fbm_cell(self, jcell):

        self.figCell.clf()
        spc_lbl = self.currentSpecies
        ax = self.figCell.add_subplot(111)
        self.figCell.subplots_adjust(left=0.11, bottom=0.3, right=0.85, top=0.9)
        ax.set_ylim((-1, 1))
        ax.set_xlabel('Energy [keV]', fontsize=fsize)
        ax.set_ylabel('Pitch angle' , fontsize=fsize)

        bid = self.butt_d['integral'].checkedId()
        integrFlag = rblist[bid]

        if integrFlag == 'Volume averaged':
            self.cell_mark.set_data(self.fbmr.r2d, self.fbmr.z2d)
            tit_lbl = 'Volume averaged, t=%6.3f' %self.fbmr.time
            zarr_lin = self.fbmr.dens_vol[spc_lbl]
        elif integrFlag == 'Theta averaged':
            jrho = np.where(self.fbmr.rho_grid == self.fbmr.x2d[jcell])[0][0]
            ind = np.where(self.fbmr.x2d == self.fbmr.x2d[jcell])
            self.cell_mark.set_data(self.fbmr.r2d[ind], self.fbmr.z2d[ind])
            tit_lbl = r'$\rho_{tor}$' + \
                      r' = %8.4f, $\theta$ averaged, t=%6.3f' %(self.fbmr.x2d[jcell], self.fbmr.time)
            zarr_lin = self.fbmr.dens_zone[spc_lbl][jrho]
        else:
            self.cell_mark.set_data([self.fbmr.r2d[jcell]], [self.fbmr.z2d[jcell]])
            tit_lbl = r'$\rho_{tor} = $ %8.4f $\theta = $%8.4f deg, t=%6.3f s' \
                %(self.fbmr.x2d[jcell], np.degrees(self.fbmr.th2d[jcell]), self.fbmr.time)
            zarr_lin = self.fbmr.fdist[spc_lbl][jcell]
        self.figDist[spc_lbl].canvas.draw() # update canvas for the red marker

        zmax_lin = np.max(zarr_lin[~np.isnan(zarr_lin)])
        zmin_lin = 1e-8*zmax_lin
        flog_max = np.log10(zmax_lin)
        flog_min = np.log10(zmin_lin)
        flog_min = self.butt_d['logMin'].value()
        flog_max = self.butt_d['logMax'].value()
        logscale = self.butt_d['logScale'].isChecked()
        indzero = np.where(zarr_lin <= zmin_lin)
        zarr_lin[indzero] = np.nan
        zarr_log = np.log10(zarr_lin)
        n_levels = 15
        if logscale:
            zmin = flog_min
            zmax = flog_max
            zarr = zarr_log
        else:
            zmin = 0
            zmax = zmax_lin
            zarr = zarr_lin

        bounds = np.linspace(zmin, zmax, n_levels)

        Emax = self.butt_d['Emax'].value()
        ax.set_xlim((0, Emax))
        ax.set_title(tit_lbl, fontsize=fsize)

        ctr = ax.contourf(1e-3*self.fbmr.e_d[spc_lbl], self.fbmr.a_d[spc_lbl], zarr, bounds)
        norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
        cb_ax = self.figCell.add_axes([0.89, 0.15, 0.03, 0.84]) #l, b, w, h
        mpl.colorbar.ColorbarBase(cb_ax, norm=norm, boundaries=bounds, ticks=bounds)
        ax.plot([0, Emax], [0, 0], 'k-')
        if integrFlag == 'Single cell':
            ax.plot([0, Emax], [ self.fbmr.trap_pit[spc_lbl][jcell],  self.fbmr.trap_pit[spc_lbl][jcell]], 'g-')
            ax.plot([0, Emax], [-self.fbmr.trap_pit[spc_lbl][jcell], -self.fbmr.trap_pit[spc_lbl][jcell]], 'g-')

        ax.figure.canvas.draw()


if __name__ == '__main__':


    app = QApplication(sys.argv)
    main = FBM()
    app.exec_()
