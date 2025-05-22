import time
import numpy as np
import b64conv


def parse_ac(f_ac, list_read=None, list_no=None):

    f = open(f_ac, 'r')
    lines = f.readlines()
    f.close()

    lines = [line.strip() for line in lines if (line.strip() != '')]

    ac_d = {}

    nlin = len(lines)
    print('#lines = %d' %nlin)

    tmp = lines[0].split()
    ac_d['runid']  = tmp[0]
    ac_d['t_id']   = int(tmp[1])
    ac_d['time']   = float(tmp[2])
    ac_d['encode'] = int(tmp[3])

    for jl, line in enumerate(lines):
        pieces = line.split()
        if 'NUMAVG' in pieces:
            try:
                ac_d['numavg'] = int(pieces[2])
            except:
                pass
        if 'MTHDAVG' in pieces:
            ac_d['mthdavg']  = int(pieces[2])
        if 'AVGTIM' in pieces:
            ac_d['avgtim']  =  b64conv.tra2dbl(pieces[2])
        if 'AVGSAMP' in pieces:
            ac_d['avgsamp'] = b64conv.tra2dbl(pieces[2])
        else:
            if line[0] != '*': # non data
                continue
            desc = pieces[0].strip()
            dtyp = desc[1]
            ndim  = int(desc[2])
            lbl = pieces[1].strip()
            if ndim == 0: # scalars
                if dtyp not in ('C', ):
                    str64 = pieces[2]
                if dtyp == 'L':
                    ac_d[lbl] = b64conv.tra2log(str64)
                elif dtyp == 'I':
                    ac_d[lbl] = b64conv.tra2int(str64)
                elif dtyp == 'R':
                    ac_d[lbl] = b64conv.tra2flt(str64)
                elif dtyp == 'D':
                    ac_d[lbl] = b64conv.tra2dbl(str64)
                else:
                    ac_d[lbl] = None
            else: # ndim > 0
                str64 = lines[jl + 1].strip()
                tmp = str64.split()
                size = [b64conv.tra2int(sval) for sval in tmp]
# If there is a list_read, list_no is ignored
                if list_read is not None:
                    if lbl not in list_read:
                        continue
                else:
                    if list_no is not None and lbl in list_no:
                        continue
                datarr = []
                jlin = jl + 2
                lin2 = lines[jlin]
                if dtyp == 'L':
                    for lin in lines[jl + 2:]:
                        if lin[0] == '*':
                            break
                        datarr.extend(b64conv.tra2log(sval) for sval in lin)
                    ac_d[lbl] = np.array(datarr, dtype=bool)
                elif dtyp == 'I':
                    lin = lines[jlin]
                    lin = lin.replace('-',' -')
                    while lin[0] != '*':
                        tmp = lin.split()
                        datarr += [b64conv.tra2int(sval) for sval in tmp]
                        jlin += 1
                        if lin[0] == '*':
                            break
                        lin = lines[jlin].replace('-',' -')
                    ac_d[lbl] = np.array(datarr, dtype=np.int32)
                if dtyp in ('R', 'D'):
                    if dtyp == 'R':
                        strlen = 6
                        nptype = np.float32
                    else:
                        strlen = 12
                        nptype = np.float64
                    while lin2[0] != '*':
                        str_arr = [lin2[start:start+strlen] for start in range(0, len(lin2), strlen)]
                        arr = []
                        for sval in str_arr:
                            if sval[0] == '_':
                                zstr = sval[3:]
                                n_zero = b64conv.tra2int(zstr)
                                arr += n_zero*[0]
                            elif len(sval.strip()) < strlen:
                                arr.append(0)
                            else:
                                if dtyp == 'R':
                                    arr.append(b64conv.tra2flt(sval))
                                else: # 'D'
                                    arr.append(b64conv.tra2dbl(sval))
                        datarr += arr
                        jlin += 1
                        if jlin == nlin:
                            break
                        lin2 = lines[jlin]
                    ac_d[lbl] = np.array(datarr, dtype=nptype)
        if ndim > 1:
            ac_d[lbl] = ac_d[lbl].reshape(size[::-1]).T

    return ac_d


if __name__ == '__main__':

    import os, config
    runid = '29783A01'

    shot = runid[:-3]
    tail = runid[-3:]

    run_dir = '%s/%s/%s' %(config.tr_clientDir, shot, tail)
    f_ac  = '%s/%s.DATA1'    %(run_dir, runid)
    list_no = ['FBM_PTCL', 'NSTAT_TRACK_XJA', \
        'TRACK_DE_FLR' , 'TRACK_DR_FLR', 'TRACK_DVPV_FLR', 'TRACK_DZ_FLR', \
        'TRACK_EINJ'   , 'TRACK_PDEP'  , 'TRACK_PHIIN'   , 'TRACK_PHIOUT', \
        'TRACK_RIN'    , 'TRACK_ROUT'  , 'TRACK_SINJ'    , 'TRACK_TIME'  , \
        'TRACK_VPV_DEP', 'TRACK_XL'    , 'TRACK_ZIN'     , 'TRACK_ZOUT']

    t1 = time.time()
    fbm_d = parse_ac(f_ac, list_no=list_no)
    t2 = time.time()
    print(t2-t1)
#    for key in fbm_d:
    for key in ('BDENS', ):
        print(key)
        print(fbm_d[key])

