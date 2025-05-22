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

    jlin = 0
    while jlin < nlin:
        line = lines[jlin].strip()
        if line[0] != '*': # skip data line, look for next medata line
            jlin += 1
            continue
        pieces = line.split()
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
# If there is a list_read, list_no is ignored
            if list_read is not None:
                if lbl not in list_read:
                    jlin += 1
                    continue
            else:
                if list_no is not None and lbl in list_no:
                    jlin += 1
                    continue
            jlin += 1
            line = lines[jlin].strip()
            pieces = line.split()
            size = [b64conv.tra2int(sval) for sval in pieces]
            if dtyp not in ('L', 'I', 'R', 'D'):
                jlin += 1
                continue
            line_arr = []
            for djlin, lin in enumerate(lines[jlin+1: ]):
                if lin[0] == '*':
                    break
                line_arr.append(lin.strip())
            if dtyp == 'L':
                strval = ''.join(line_arr).upper()
                ac_d[lbl] = np.array(list(strval)) == 'T'
            elif dtyp == 'I':
                words = ' '.join(line_arr).replace('-',' -').split()
                ac_d[lbl] = b64conv.base64_to_int_vec(words)
            elif dtyp  == 'R':
                datarr = []
                strlen = 6
                loc_lines = [line_loc.ljust((len(line_loc) + strlen-1) // strlen * strlen).replace(' ', '0') for line_loc in line_arr]
                full_string = ''.join(loc_lines)
                str_arr = [full_string[start:start+strlen] for start in range(0, len(full_string), strlen)]
#                    myarr = b64conv.tra2flt_numpy(str_arr)
                for jval, sval in enumerate(str_arr):
                    if sval[0] == '_':
                        zstr = sval[3:]
                        n_zero = b64conv.tra2int(zstr)
                        datarr.extend(n_zero*[0])
                    else:
                        datarr.append(b64conv.tra2flt(sval))
                ac_d[lbl] = np.array(datarr, dtype=np.float32)
                print(lbl, ac_d[lbl])
            elif dtyp == 'D':
                datarr = []
                strlen = 12
                for lin in line_arr:
                    str_arr = [lin[start:start+strlen] for start in range(0, len(lin), strlen)]
                    for sval in str_arr:
                        if sval[0] == '_':
                            zstr = sval[3:]
                            n_zero = b64conv.tra2int(zstr)
                            datarr.extend(n_zero*[0])
                        elif len(sval.strip()) < strlen:
                            datarr.append(0)
                        else:
                            datarr.append(b64conv.tra2dbl(sval))
                ac_d[lbl] = np.array(datarr, dtype=np.float64)
            jlin += 1
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

    mylist = ['FBM', 'BDENS']
    t1 = time.time()
    fbm_d = parse_ac(f_ac, list_read=mylist)
    t2 = time.time()
    print(t2-t1)
#    for key in fbm_d:
    for key in ('BDENS', 'NUMAVG', 'AVGTIM', 'AVGSAMP'):
        print(key)
        val = fbm_d[key]
        if hasattr(val, 'dtype'):
            print(val.dtype)
        print(val)
    print(len(fbm_d))

