"""
This file implements base64-like encoding of DT codes.

We use the following DT code as an example here::

   >>> code = [(6, -8), (-10, 12), (-2, 4)]

(To see the link determined by this code, and the way that the crossings
are numbered, run the command Link('DT:[(6,-8),(-10,12),(-2,4)]').view()
in SnapPy and select the "Info->DT Labels" menu item from the PLink window.)

Associated to a DT code are the flips, which contain information about the signs
of the crossings.  The flips can be deduced from the DT code, but the algorithm
for laying out a link projection from a DT code runs faster if the flips are
known in advance. Here, the associated flips are::
  
   >>> from spherogram import DTcodec
   >>> flips = DTcodec(code).flips
   >>> flips
   [False, True, False, True, True, False]

To use the base64-like encoding (without or with flips)::

   >>> encode_base64_like_DT_code(code)
   '1dcgOcQmcIe'
   >>> encode_base64_like_DT_code(code, flips)
   '1dcgOcQmcIeA'

To decode::

   >>> decode_base64_like_DT_code('1dcgOcQmcIe')
   ([(6, -8), (-10, 12), (-2, 4)], None)
   >>> decode_base64_like_DT_code('1dcgOcQmcIeA')
   ([(6, -8), (-10, 12), (-2, 4)], [False, True, False, True, True, False])

The base64-like encoding of DT codes is based on base64-like encoding used in
Burton's isomorphism signatures ("base64-like" because two of the characters
differ from the base64 encoding specified in RFC 1521).

To encode an integer, we write it in base 64 and translate each base 64-digit to
an ASCII character (see _base64LikeEncoding). We use big endian (with groups of
6 bits instead of bytes) and zero's complement.

The very first character of a base64-like encoding of a DT code is always a
digit 1-9 which specifies the number of ASCII characters used to encode each of
the (even) integer crossing labels.  For the example above, the first character
is "1" because the crossing label with the highest absolute value, namely 12,
fits into 5 bits.

The next integer in the encoding specifies the number of components, here 3.

This is followed by the encodings of the components, each of which consists of
the length of the component followed by the crossing labels associated to
the component.

In this example the encoded integers following the first character "1" are::
3 2 6 -8 2 -10 12 2 -2 4.

If there are flips, they are simply stored as bit fields following the encoding
of the last component.
"""

_base64LikeEncoding = (
    'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789+-')

def _unsigned_int_to_char(i):
    """
    Convert an integer 0-63 to ASCII character.
    """
    return _base64LikeEncoding[i]

def _char_to_unsigned_int(char):
    """
    Convert an ASCII character to an integer 0-63
    """

    i = ord(char)
    if i >= 97 and i <= 122:  # a z
        return i - 97
    if i >= 65 and i <= 90: # A Z
        return i - 39
    if i >= 48 and i <= 57: # 0 9
        return i + 4
    if i == 43:
        return 62
    return 63

def _chars_to_int(chars):
    """
    Take a string of ASCII characters and convert it to integer using the
    base64-like scheme described above.
    """
    
    value = 0

    for pos, char in enumerate(chars):
        i = _char_to_unsigned_int(char)
        if pos == 0:
            # The very first bit gives the sign
            if i & (1 << 5):
                sign = -1
            else:
                sign = +1
            i = i & 31
        value = (value << 6) + i

    return sign * value

def _consume_int_and_advance(chars, pos, num_chars):
    """
    Read the num_chars characters from the string chars and interpret it as
    base64-like encoded integer. Also return the end of that integer. This
    is supposed to be in a pattern like this::

        >>> chars = "abcdef"
        >>> pos = 0 # Start at the beginning
        >>> first_integer,  pos = _consume_int_and_advance(chars, pos, 2)
        >>> second_integer, pos = _consume_int_and_advance(chars, pos, 2)
        >>> third_integer,  pos = _consume_int_and_advance(chars, pos, 2)
    """

    end = pos + num_chars
    return _chars_to_int(chars[pos:end]), end

def _int_to_chars(value, num_chars):
    """
    Encode the given integer using base64-like encoding with num_chars
    characters.
    """
    
    abs_value = abs(value)

    chars = ""
    for pos in range(num_chars):
        if pos == num_chars - 1:
            if value < 0:
                abs_value |= (1 << 5)

        chars = _unsigned_int_to_char(abs_value & 63) + chars
        abs_value >>= 6

    return chars

def _get_num_chars(code):
    """
    Given a DT code, determine the minimal number of characters to encode
    each integer so that the code fits.
    """
    max_number = max([
            max([abs(crossing) for crossing in comp])
            for comp in code])
    num_chars = 1
    while max_number > 31:
        num_chars += 1
        max_number >>= 6
    return num_chars

def _encode_DT_code(code):
    """
    Given a DT code, convert to base64-like encoding.
    """
    num_chars = _get_num_chars(code)
    
    # 1 -> "1", 2 -> "2", ...
    chars  = _unsigned_int_to_char(num_chars + 52)

    chars += _int_to_chars(len(code), num_chars)
    for comp in code:
        chars += _int_to_chars(len(comp), num_chars)
        for crossing in comp:
            chars += _int_to_chars(crossing, num_chars)
        
    return chars

def _boolean_to_int(b):
    return 1 if b else 0

def _encode_flips(flips):
    """
    Given flips, encode as bit field.
    """
    chars = ""

    padded_flips = flips + 5 * [ 0 ]

    for pos in range(len(padded_flips) // 6):
        val = 0
        for i in range(6):
            val |= _boolean_to_int(padded_flips[6 * pos + i]) << i
        chars += _unsigned_int_to_char(val)

    return chars

def encode_base64_like_DT_code(code, flips = None):
    """
    Given a DT code and optionally flips, convert to base64-like encoding.
    """
    if flips:
        return _encode_DT_code(code) + _encode_flips(flips)
    else:
        return _encode_DT_code(code)

def _decode_DT_code(chars):
    """
    Given a base64-like encoding, return the DT code and the position where
    the remaining characters not consumed yet start.
    """
    code = []

    # "1" -> 1, "2" -> 2
    num_chars = _char_to_unsigned_int(chars[0]) - 52

    num_components, pos = _consume_int_and_advance(chars, 1, num_chars)
    
    for i in range(num_components):
        component = []
        
        num_crossings, pos = _consume_int_and_advance(chars, pos, num_chars)
        for j in range(num_crossings):
            crossing, pos = _consume_int_and_advance(chars, pos, num_chars)
            component.append(crossing)

        code.append(tuple(component))

    return code, pos

def _decode_flips(chars):
    """
    Read a bit field from base64-like encoding.
    """
    flips = []
    for char in chars:
        val = _char_to_unsigned_int(char)
        for j in range(6):
            flips.append(bool((val >> j) & 1))
    return flips

def _empty_to_none(l):
    if not l:
        return None
    return l

def decode_base64_like_DT_code(chars):
    """
    Given a base64-like encoding, return the DT code and, if present in the
    encoding as well, return the flips, otherwise None.
    """
    code, pos = _decode_DT_code(chars)

    num_crossings = sum([len(component) for component in code])
    flips = _decode_flips(chars[pos:])[:num_crossings]

    return code, _empty_to_none(flips)
    
