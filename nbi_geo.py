import numpy as np
import rot_matrix

rot = rot_matrix.ROT_MATRIX


def dd2xy(inj_tor_box, src_slit, op, rabsca, phi_box):

    R_src = src_slit + op
    tmp   = rot(-inj_tor_box, R_src, rabsca, x_cen=op)
    alpha = np.arctan(tmp.y/tmp.x)
    beta  = np.arctan(tmp.y/(tmp.x - op))
    rm = rot(phi_box, tmp.x, tmp.y)
    phi = alpha - beta
    return rm.x, rm.y, phi


def tr2xy(rtcena, midpl2, alpha, phi_box):

    tmp = rot(-alpha, midpl2, rtcena)
    rm  = rot(phi_box, tmp.x, tmp.y)
    phi = np.arctan(rtcena/midpl2)
    return rm.x, rm.y, phi


class AUG_TR:


    def __init__(self, nbi=None, nshot=None, runid='', raug=True):

# distances: cm; angles: rad
# op: distance box-torus axis
# phi_box: absolute angle of NBI box (does not matter in poloidal cross section)
# inj_tor_box: toroidal angle between mid source and torus axis
# inj: toroidal angle between source and box aligned
# src_slit: distance ion source - aperture, for the center of the beamlines
# src_hw: horizontal semi-displacement of sources
# xybsca: source vertical elevation w.r.t. midplane
# rabsca: source horizontal half width w.r.t. mid-point
# RTCENA: R_tang
# vert_incl: vertical inclination
# XLBAPA: length = midpl1 / COS(vert_incl)
# midpl2: Grid - P(RTCENA)   projected on horizont. midplane
# XLBTNA: length = midpl2 / COS(vert_incl)
# XYBAPA: elevation at slit


# Rectangular half aperture in point P
      self.rapedga = 16
      self.xzpedga = 19

      self.phi_box = np.array(4*[np.radians(33.75)] + 4*[np.radians(209)])

      if raug:
# Read AUG shotfile

          if (nshot is None) and (nbi is None):
              return
          if nbi is None:
              import read_nbi
              print('Test')
              nbi = read_nbi.NBI(nshot)
          for key in nbi.__dict__.keys():
              print(key)
          self.op = np.array(4*[284.20] + 4*[329.63])
          self.inj_tor_box = np.array(4*[np.radians(15)] + 4*[np.radians(18.9)])
          src_slit = 650.
          self.src_slit = np.array(8*[src_slit])
          self.src_hw  = 47.
          mag_ang = np.arctan(self.src_hw/src_slit)
          self.midpl1  = np.hypot(self.src_slit, self.src_hw)

          self.magic_angle = np.array(2*[-mag_ang, mag_ang, mag_ang, -mag_ang])

          self.xybsca    = np.array([60, 60, -60, -60, 60, 70, -70, -60])
          self.rabsca    = self.src_hw*np.array(2*[1, -1, -1, 1])
          self.vert_incl = np.radians(nbi.beta_sf_deg)
          self.theta_los = np.array(2*[-1, -1, 1, 1])*self.vert_incl
          self.alpha     = self.inj_tor_box + self.magic_angle
          self.rtcena    = self.op*np.sin(self.alpha)
          self.midpl2    = self.midpl1 + np.sqrt(self.op**2 - self.rtcena**2)
          self.xlbtna    = self.midpl2/np.cos(self.vert_incl)
          self.xlbapa    = self.midpl1/np.cos(self.vert_incl)
          self.xybapa    = self.xybsca - np.sign(self.xybsca)*self.src_slit*np.tan(self.vert_incl)

          self.xsrc, self.ysrc, self.phi_los = dd2xy(self.inj_tor_box, self.src_slit, self.op, \
                                self.rabsca, self.phi_box)
          tmp = np.degrees(np.arctan(self.ysrc/self.xsrc))
          ind = (self.phi_box > np.pi)
          tmp[ind] += 180
          self.xbzeta = tmp

# Species mix

          self.ffulla = nbi.part_mix[:, 0]
          self.fhalfa = nbi.part_mix[:, 1]
          self.ffulla[np.isnan(self.ffulla)] = 1.
          self.fhalfa[np.isnan(self.fhalfa)] = 0.
          self.einja  = 1e3*nbi.einj_kev

      else:

# Read from namelist

          import tr_path
          from parsenml import parsenml

          tr = tr_path.TR_PATH(runid)
          print('Namelist: %s' %tr.fnml)

          self.rtcena   = parsenml(tr.fnml, 'RTCENA', fmt=5)
          self.xlbapa   = parsenml(tr.fnml, 'XLBAPA', fmt=5)
          self.xybapa   = parsenml(tr.fnml, 'XYBAPA', fmt=5)
          self.xybsca   = parsenml(tr.fnml, 'XYBSCA', fmt=5)
          self.xlbtna   = parsenml(tr.fnml, 'XLBTNA', fmt=5)
          self.xbzeta   = parsenml(tr.fnml, 'XBZETA', fmt=5)
          self.src_slit = parsenml(tr.fnml, 'FOCLRA', fmt=5)
          self.ffulla   = parsenml(tr.fnml, 'FFULLA', fmt=5)
          self.fhalfa   = parsenml(tr.fnml, 'FHALFA', fmt=5)
          self.einja    = parsenml(tr.fnml, 'EINJA' , fmt=5)

# Derived TRANSP variables
    
          self.theta_los = np.arcsin((self.xybapa - self.xybsca)/self.xlbapa)
          self.vert_incl = np.abs(self.theta_los)

          self.midpl2    = self.xlbtna*np.cos(self.vert_incl)
          self.midpl1    = self.xlbapa*np.cos(self.vert_incl)
          self.op        = np.hypot((self.midpl2 - self.midpl1), self.rtcena)
          self.alpha     = np.arcsin(self.rtcena/self.op)

          
          self.xsrc, self.ysrc, self.phi_los = tr2xy(self.rtcena, self.midpl2, \
                                                     self.alpha, self.phi_box)
          self.magic_angle = np.arccos(self.src_slit/self.midpl1) # missing sign
          self.rabsca = self.midpl1*np.sin(self.magic_angle)      # missing sign
          self.inj_tor_box = self.alpha - self.magic_angle

if __name__ == "__main__":

    aug = AUG_TR(28053)
    print(aug.xybsca)
    print(aug.xybapa)
    print(aug.vert_incl)

    aug1 = AUG_TR(23076, runid='23076W01', raug=False)
    print(aug1.xybsca)
    print(aug1.xybapa)
    print(aug1.vert_incl)
