import base64, struct
import numpy as np

tra_b64 = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ<>='
rfc3548 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/='

tr2def = {c: r for c, r in zip(tra_b64, rfc3548)}
char2pos = {c: i for i, c in enumerate(tra_b64)}

transp2def = str.maketrans(tra_b64, rfc3548)

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
    str_ieee = tra2ieee(sflt) + '=='
    sval = base64.b64decode(str_ieee.encode())
    num = struct.unpack_from(fmt, sval)[0]
    return -num #TRANSP convention

def tra2flt_numpy(sflt_list):
    ieee_strs = [s.translate(tr2def) + '==' for s in sflt_list]
    decoded = b''.join(base64.b64decode(s) for s in ieee_strs)
    return np.frombuffer(decoded, dtype='>f')

def tra2flt_batch(sflt_list, fmt='>f'):
# Translate all strings and pad with '=='
    ieee_strs = [s.translate(transp2def) + '==' for s in sflt_list]    
# Decode all base64 strings
    decoded = [base64.b64decode(s) for s in ieee_strs]
# Unpack floats using struct
    return [struct.unpack_from(fmt, d)[0] for d in decoded]

def tra2dbl(sdbl):
    (exp, ic1, ic2) = base64_to_int_vec((sdbl[:2], sdbl[2:7], sdbl[7:12]))
    if exp > 1000:
        exp -= 1000
        sgn = 1
    else:
        sgn = -1
    exp -= 500
    if ic1 == 0:
        z = 0.
    else:
        z  = 10.**exp *ic1 *1.e-9
        z += 10.**exp *ic2 *1.e-18
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
