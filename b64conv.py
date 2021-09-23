import base64, struct
import numpy as np

tra_b64 = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ<>='
rfc3548 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/='

tr2def = {}
for pos, char in enumerate(tra_b64):
    tr2def[char] = rfc3548[pos]

def tra2ieee(str_in):

    str_out = ''
    for char in str_in:
        str_out += tr2def[char]

    return(str_out)

def tra2log(slog):

    return (slog.strip().upper() == 'T')

def tra2int(sint):

    sint = sint.strip()
    if sint[0] == '-':
        sign = -1
        sint = sint[1:]
    else:
        sign = 1

    num = 0
    for jexp, char in enumerate(sint[::-1]):
        num += 64**jexp * tra_b64.find(char)
    return sign*num

def tra2flt(sflt, fmt='>f'):

    str_ieee = tra2ieee(sflt)

# Fix padding
    while(len(str_ieee) % 8):
        str_ieee += '='
    sval = base64.b64decode(str_ieee.encode())
    num = struct.unpack_from('>f', sval)
#    num = np.ndarray(sval, '>f')

    return -num[0] #TRANSP convention

def tra2dbl(sdbl):

    exp = tra2int(sdbl[:2])
    ic1 = tra2int(sdbl[2:7])
    ic2 = tra2int(sdbl[7:12])
    if exp > 1000:
        exp -= 1000
        sgn = 1
    else:
        sgn = -1
    exp -= 500
    if ic1 == 0:
        z = 0.
    else:
        z  = 10**exp *ic1 *1.e-9
        z += 10**exp *ic2 *1.e-18
        z *= sgn

    return z


if __name__ == '__main__':

    str64 = '7Q'
    print(tra2int(str64))
    
    print('\npi')
    str64 = 'M4AfT0'
    num = tra2flt(str64)
    print(num)

    print('\nEPS32')
    seps32 = 'J03pmg'
    num = tra2flt(seps32)
    print(num)

# Double precision

    print('\nRAXIS')
    str_tra = 'nva8vZx2VrBI'
    print(tra2dbl(str_tra))
