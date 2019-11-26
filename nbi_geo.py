import numpy as np
import rot_matrix

rot = rot_matrix.ROT_MATRIX


def tr2xy(rtcena, midpl2, alpha, phi_box):

    tmp = rot(-alpha, midpl2, rtcena)
    rm  = rot(phi_box, tmp.x, tmp.y)
    phi = np.arctan2(rtcena, midpl2)
    return rm.x, rm.y, phi


class AUG_TR:


    def __init__(self, nbi=None, nshot=None, runid=''):

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
