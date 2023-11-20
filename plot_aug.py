import numpy as np
import rot_matrix


def nbi_plot(nbis=[], runid=None):

    import nbi_geo, los

    nbis = np.atleast_1d(nbis)
    aug_geo = {}
    aug = nbi_geo.AUG_TR(0, runid=runid)
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
        geom = los.PHILOS(aug_geo)
        xlin.append([geom.xline[0], geom.xline[-1]])
        ylin.append([geom.yline[0], geom.yline[-1]])
        rlin.append(geom.rline)
        zlin.append(geom.zline)
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
