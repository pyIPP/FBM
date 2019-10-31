import sys
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
import numpy as np
import dd_20180130


sf = dd_20180130.shotfile()

kv_def = {'INJ1':60, 'INJ2':93}


class NBI:


    def __init__(self, nshot, tbeg=0., tend=10., nbi_diff=1.e3, nbi_min=3e4):


        diag = 'NIS'
        sig = 'PNIQ'

        beta_deg = {}
        einj_nis = {}
        mass_amu = {}
        pow_mix  = {}
        if sf.Open(diag, nshot):
            tim = sf.GetTimebase(sig)
            dat = sf.GetSignal(sig)
            for box in ('INJ1','INJ2'):
                beta_deg[box] = sf.GetParameter(box, 'BETA')
                mass_amu[box] = sf.GetParameter(box, 'M')
                einj_nis[box] = sf.GetParameter(box, 'UEXQ')
                pow_mix [box] = sf.GetParameter(box, 'SPEC')
                ind = (einj_nis[box] <= 0)
                einj_nis[box][ind] = kv_def[box]
            sf.Close()
        else:
            print('No NIS shotfile found for shot %d' %nshot)
            return

# Parameters

        self.n_box = np.array( [ len(einj_nis['INJ1']), len(einj_nis['INJ2']) ] )
        n_nbi = np.sum(self.n_box)

        if nshot > 14050:
            self.beta_sf_deg = np.append(beta_deg['INJ1'], beta_deg['INJ2'])
        else:
            self.beta_sf_deg = np.array(8*[4.8891] )
