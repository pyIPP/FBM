import base64, struct
import numpy as np

tra_b64 = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ<>='
rfc3548 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/='

tr2def = {c: r for c, r in zip(tra_b64, rfc3548)}
char2pos = {c: i for i, c in enumerate(tra_b64)}

def tra2ieee(str_in):
    return ''.join(tr2def[char] for char in str_in)

def tra2int(sint):
    sint = sint.strip()
    sign = -1 if sint.startswith('-') else 1
    if sint[0] in '+-':
        sint = sint[1:]
    num = sum(char2pos[char] * (64 ** i) for i, char in enumerate(reversed(sint)))
    return sign * num

def base64_to_int_vec(strs):
    def decode_one(s):
        s = s.strip()
        sign = -1 if s.startswith('-') else 1
        if s[0] in '+-':
            s = s[1:]
        return sign * sum(char2pos[c] * (64 ** i) for i, c in enumerate(reversed(s)))
    
    return np.array([decode_one(s) for s in strs], dtype=np.int64)

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
