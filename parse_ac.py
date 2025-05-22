import time
import numpy as np
import b64conv

strlen_d = {'R': 6, 'D': 12}
nptype_d = {'R': np.float32, 'D': np.float64}

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

    jlin = 0
    while jlin < nlin:
        line = lines[jlin].strip()
        pieces = line.split()
        if 'NUMAVG' in pieces:
            jlin += 1
            try:
                ac_d['numavg'] = int(pieces[2])
            except:
                pass
        if 'MTHDAVG' in pieces:
            ac_d['mthdavg']  = int(pieces[2])
            jlin += 1
        elif 'AVGTIM' in pieces:
            ac_d['avgtim']  =  b64conv.tra2dbl(pieces[2])
            jlin += 1
        elif 'AVGSAMP' in pieces:
            ac_d['avgsamp'] = b64conv.tra2dbl(pieces[2])
            jlin += 1
        else:
            if line[0] != '*': # data line
                jlin += 1
                continue
            desc = pieces[0].strip()
            dtyp = desc[1]
            ndim = int(desc[2])
            lbl = pieces[1].strip()
            if ndim == 0: # scalars, values on the same line
                if dtyp not in ('C', ):
                    str64 = pieces[2]
                if dtyp == 'L':
                    ac_d[lbl] = str64.strip().upper() == 'T'
                elif dtyp == 'I':
                    ac_d[lbl] = b64conv.tra2int(str64)
                elif dtyp == 'R':
                    ac_d[lbl] = b64conv.tra2flt(str64)
                elif dtyp == 'D':
                    ac_d[lbl] = b64conv.tra2dbl(str64)
                else:
                    ac_d[lbl] = None
                jlin += 1
            else: # ndim > 0, start collecting data from the following line
                jlin += 1
                line = lines[jlin].strip()
                pieces = line.split()
                size = [b64conv.tra2int(sval) for sval in pieces]
#                print(lbl, 'ciao', jlin, '|', size)
# If there is a list_read, list_no is ignored
                if list_read is not None:
                    if lbl not in list_read:
                        jlin += 1
                        continue
                else:
                    if list_no is not None and lbl in list_no:
                        jlin += 1
                        continue
                if dtyp not in ('L', 'I', 'R', 'D'):
                    jlin += 1
                    continue
                line_arr = []
                for djlin, lin in enumerate(lines[jlin+1: ]):
                    if lin[0] == '*':
                        break
                    line_arr.append(lin)
                if dtyp == 'L':
                    strval = ''.join(line_arr).upper()
                    ac_d[lbl] = np.array(list(strval)) == 'T'
                elif dtyp == 'I':
                    words = ' '.join(line_arr).replace('-',' -').split()
                    ac_d[lbl] = b64conv.base64_to_int_vec(words)
                    print(lbl, ac_d[lbl])
                if dtyp in ('R', 'D'):
                    datarr = []
                    strlen = strlen_d[dtyp]
                    nptype = nptype_d[dtyp]
                    for lin in line_arr:
                        str_arr = [lin[start:start+strlen] for start in range(0, len(lin), strlen)]
                        arr = []
                        for sval in str_arr:
                            if sval[0] == '_':
                                zstr = sval[3:]
                                n_zero = b64conv.tra2int(zstr)
                                arr.extend(n_zero*[0])
                            elif len(sval.strip()) < strlen:
                                arr.append(0)
                            else:
                                if dtyp == 'R':
                                    arr.append(b64conv.tra2flt(sval))
                                else: # 'D'
                                    arr.append(b64conv.tra2dbl(sval))
                        datarr.extend(arr)
                    ac_d[lbl] = np.array(datarr, dtype=nptype)
                jlin += djlin + 1
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
    print(len(fbm_d))

