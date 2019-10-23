import numpy as np
import rot_matrix


def nbi_plot(nbis=[], runid=None, raug=False):

    import nbi_geo, los

    nbis = np.atleast_1d(nbis)
    aug_geo = {}
    aug = nbi_geo.AUG_TR(0, runid=runid, raug=raug)
    xend = np.array(4*[-60] + [40, 250, 250, 40])

    xlin = []
    ylin = []
    rlin = []
    zlin = []
    for jnbi in nbis:
        jnb = jnbi - 1
        aug_geo['xend']  = xend[jnb]
        aug_geo['x0']    = aug.xsrc[jnb]
        aug_geo['y0']    = aug.ysrc[jnb]
        aug_geo['z0']    = aug.xybsca[jnb]
        aug_geo['theta'] = aug.theta_los[jnb]
        aug_geo['phi']   = aug.phi_los[jnb]
        xlin.append([los.PHILOS(aug_geo).xline[0], los.PHILOS(aug_geo).xline[-1]])
        ylin.append([los.PHILOS(aug_geo).yline[0], los.PHILOS(aug_geo).yline[-1]])
        rlin.append(los.PHILOS(aug_geo).rline)
        zlin.append(los.PHILOS(aug_geo).zline)
    return xlin, ylin, rlin, zlin
    

class STRUCT:


    coils = {'R':2.3, 'phi_beg':1.402, 'dphi':0.730, 'phi_shift':0.7854}


    def __init__(self, f_in= '/afs/ipp/aug/ads-diags/common/diaggeom.data/tor.data'):

# Toroidal

        print('Reading structure data from file %s' %f_in)
        f = open(f_in,'r')

        xtor_struct = {}
        ytor_struct = {}
        jstr = 0
        xtor_struct[jstr] = []
        ytor_struct[jstr] = []
        for line in f.readlines():
            if (line.strip() != ''):
                xval, yval = line.split()
                xtor_struct[jstr].append(float(xval))
                ytor_struct[jstr].append(float(yval))
            else:
                jstr += 1
                xtor_struct[jstr] = []
                ytor_struct[jstr] = []
        f.close()
        nstr = jstr

# Rotate

        gamma = -3*np.pi/8 # 3 sectors out of 16
        self.tor_str = {}
        self.tor_old = {}
        for jstr in range(nstr):
            x_in = np.array(xtor_struct[jstr])
            y_in = np.array(ytor_struct[jstr])
            self.tor_str[jstr] = rot_matrix.ROT_MATRIX(gamma, x_in, y_in)
            self.tor_old[jstr] = rot_matrix.ROT_MATRIX(0,     x_in, y_in)


if __name__ == "__main__":

    import sys
    sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib/')
    import matplotlib.pylab as plt
    import get_gc

    nshot = 28746
    aug = plt.figure('AUG', figsize=(18, 11))
    plt.figtext(0.5, 0.95, 'AUG #%d' %nshot, ha='center', fontsize=20)
    aug.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.9, hspace=0., wspace=0.3)

# Poloidal

    plt.subplot2grid((1, 5),(0, 0),aspect='equal', colspan=2)
    plt.title('Poloidal cross section')
    plt.xlabel('R [m]')
    plt.ylabel('z [m]')

    gc_r, gc_z = get_gc.get_gc()
    for key in gc_r.keys():
        print key
        if key.strip() in ('VESiR',):
            plt.plot(gc_r[key], gc_z[key], 'b-')
        else:
            if key.strip() in ('ICRHa',):
                plt.fill(gc_r[key], gc_z[key], fc='#ccbbbb',alpha=1, ec='b')
            else:
                plt.fill(gc_r[key], gc_z[key], fc='#dddddd',alpha=1, ec='b')

# Above

    plt.subplot2grid((1,5),(0,2),aspect='equal', colspan=3)
    plt.title('View from above')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
#    dic = STRUCT().tor_str
    dic1 = STRUCT().tor_old
    for key in dic1.keys():
        plt.plot(dic1[key].x,dic1[key].y,'g-')

# B-coils
    coil_d = STRUCT().coils
    nphi = 21
    phi_min = coil_d['phi_beg']
    Rlbl = 2
    n_coils = 8
    for jcoil in range(n_coils):
        phi_plot = np.linspace(phi_min, phi_min + coil_d['dphi'])
        plt.plot(coil_d['R']*np.cos(phi_plot), coil_d['R']*np.sin(phi_plot), 'b-', linewidth=1)
        phi_lbl = phi_min + 0.5*coil_d['dphi']
        xlbl = Rlbl*np.cos(phi_lbl)
        ylbl = Rlbl*np.sin(phi_lbl)
        plt.text(xlbl, ylbl, 'B-%d' %(jcoil + 1), fontsize=14, ha='center', va='center')
        phi_min += coil_d['phi_shift']

    plt.show()
