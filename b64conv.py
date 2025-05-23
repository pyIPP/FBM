import base64, struct

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
    return [decode_one(s) for s in strs]

def tra2flt(sflt, fmt='>f'):
    str_ieee = tra2ieee(sflt) + '=='
    sval = base64.b64decode(str_ieee.encode())
    num = struct.unpack_from(fmt, sval)[0]
    return -num #TRANSP convention

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
        scale = 10.0 ** exp
        z = scale * (ic1 * 1e-9 + ic2 * 1e-18) * sgn
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
