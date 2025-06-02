import numpy as np
import b64conv


def decodeFloatArray(line_arr, strlen=12, transpMap=b64conv.tra2dbl):
    datarr = []
    str_arr = [lin[start:start+strlen] for lin in line_arr for start in range(0, len(lin), strlen)]
    for sval in str_arr:
        if sval.startswith('_'):
            n_zero = b64conv.tra2int(sval[3:])
            datarr.extend(n_zero*[0])
        elif len(sval.strip()) < strlen:
            datarr.append(0)
        else:
            datarr.append(transpMap(sval))
    return datarr


def decoder(f_ac, list_read=None, list_no=None):

    f = open(f_ac, 'r')
    lines = f.readlines()
    f.close()

    lines = [line.strip() for line in lines if line.strip()]
    nlin = len(lines)
    indStars = [i for i, line in enumerate(lines) if line.startswith('*')]
    n_star = len(indStars)

    ac_d = {}
    piecesFirst = lines[0].split()
    ac_d['runid']  = piecesFirst[0]
    ac_d['t_id']   = int(piecesFirst[1])
    ac_d['time']   = float(piecesFirst[2])
    ac_d['encode'] = int(piecesFirst[3])

    for j, jlin in enumerate(indStars):
        pieces = lines[jlin].split()
        desc = pieces[0]
        dtyp = desc[1]
        if dtyp == 'C': # skip complex data
            continue
        ndim = int(desc[2])
        lbl = pieces[1]
        if ndim == 0: # scalars, value on the same line
            str64 = pieces[2]
            if dtyp == 'L':
                ac_d[lbl] = str64 == 'T'
            elif dtyp == 'I':
                ac_d[lbl] = b64conv.tra2int(str64)
            elif dtyp == 'R':
                ac_d[lbl] = b64conv.tra2flt(str64)
            elif dtyp == 'D':
                ac_d[lbl] = b64conv.tra2dbl(str64)
            else:
                ac_d[lbl] = None
        else: # ndim > 0, start collecting data from the following line
# If the keyword list_read is set, the list_no argument is ignored
            if list_read is not None:
                if lbl not in list_read:
                    continue
            else:
                if list_no is not None and lbl in list_no:
                    continue

            if j == n_star-1:
                line_arr = lines[jlin+2: ]
            else:
                line_arr = lines[jlin+2: indStars[j+1]]

            if dtyp == 'L':
                strval = ''.join(line_arr)
                ac_d[lbl] = np.array(list(strval)) == 'T'
            elif dtyp == 'I':
                words = ' '.join(line_arr).replace('-',' -').split()
                datarr = [b64conv.tra2int(word) for word in words]
                ac_d[lbl] = np.array(datarr, dtype=np.int64)
            elif dtyp == 'R':
                datarr = decodeFloatArray(line_arr, strlen=6, transpMap=b64conv.tra2flt)
                ac_d[lbl] = np.array(datarr, dtype=np.float32)
            elif dtyp == 'D':
                datarr = decodeFloatArray(line_arr)
                ac_d[lbl] = np.array(datarr, dtype=np.float64)
        if ndim > 1:
            pieces = lines[jlin+1].split()
            size = [b64conv.tra2int(sval) for sval in pieces]
            ac_d[lbl] = ac_d[lbl].reshape(size[::-1]).T

    return ac_d


if __name__ == '__main__':

    import os, config, time
    runid = '29783A01'

    shot = runid[:-3]
    tail = runid[-3:]

    run_dir = '%s/%s/%s' %(config.tr_clientDir, shot, tail)
    f_ac  = '%s/%s.DATA1'    %(run_dir, runid)
    list_no = ['NSTAT_TRACK_XJA', \
        'TRACK_DE_FLR' , 'TRACK_DR_FLR', 'TRACK_DVPV_FLR', 'TRACK_DZ_FLR', \
        'TRACK_EINJ'   , 'TRACK_PDEP'  , 'TRACK_PHIIN'   , 'TRACK_PHIOUT', \
        'TRACK_RIN'    , 'TRACK_ROUT'  , 'TRACK_SINJ'    , 'TRACK_TIME'  , \
        'TRACK_VPV_DEP', 'TRACK_XL'    , 'TRACK_ZIN'     , 'TRACK_ZOUT']

    mylist = ['FBM', 'BDENS']
    t1 = time.time()
#    fbm_d = decoder(f_ac)
    fbm_d = decoder(f_ac, list_no=list_no)
#    fbm_d = decoder(f_ac, list_read=mylist)
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
