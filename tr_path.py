import os
import config

class TR_PATH:

    def __init__(self, runid, nshot=None, t_id=1):

        if nshot is None:
            shot = runid[:-3]
        else:
            shot = str(nshot)
        tail = runid[-3:]
        self.udb     = '%s/%s' %(config.udbDir, shot)
        self.run_dir = '%s/%s' %(sonfig.tr_clientDir, shot)
        self.cdf_dir = '%s/%s'    %(self.run_dir, tail)
        self.fnml = '%s/%sTR.DAT' %(self.run_dir, runid)
        self.fcdf = '%s/%s.CDF'   %(self.cdf_dir, runid)
        self.ffbm = '%s/%s_fi_%d.cdf' %(self.cdf_dir, runid, int(t_id))
