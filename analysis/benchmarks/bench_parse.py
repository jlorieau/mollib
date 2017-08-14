import timeit


line = "ATOM   1617  C   SER A 244      21.367   6.240  37.294  1.00 19.56           C"


# Regex Implementation - 2.196 us per match, 0.138 s per 62875 matches
setup = """
import re
re_atom = re.compile((r"(?P<type>ATOM  |HETATM)"
                            "(?P<number>[\s\d]{5}) "
                            "(?P<name>[\s\w]{4})"
                            "(?P<alt_loc>[\w\s])"
                            "(?P<residue_name>[\w\s]{3}) "
                            "(?P<chain>[\s\w]{1})"
                            "(?P<residue_number>[\s\w]{4})"
                            "(?P<icode>[\w\s])   "
                            "(?P<x>[\d\s\.\-]{8})"
                            "(?P<y>[\d\s\.\-]{8})"
                            "(?P<z>[\d\s\.\-]{8})"
                            "(?P<occupancy>[\d\s\.\-]{6})"
                            "(?P<B_factor>[\d\s\.\-]{6})          "
                            "(?P<element>[\s\w]{2})"
                            "(?P<charge>[\d\s\.\-]{2})?"))
line = """
setup += '"' + line + '"\n'

number = 1000000
time = timeit.timeit("re_atom.match(line)", setup, number=number)
print("Regex: {:.3f} us per match, {:.3f} s per 62875 matches".format(time*1000000/number,
                                                                      time*62875/number))

# Split line implementation - 5.555 us per match, 0.349 s per 62875 matches
#  1 -  6        Record name   "ATOM  "
#  7 - 11        Integer       serial       Atom  serial number.
# 13 - 16        Atom          name         Atom name.
# 17             Character     altLoc       Alternate location indicator.
# 18 - 20        Residue name  resName      Residue name.
# 22             Character     chainID      Chain identifier.
# 23 - 26        Integer       resSeq       Residue sequence number.
# 27             AChar         iCode        Code for insertion of residues.
# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
# 55 - 60        Real(6.2)     occupancy    Occupancy.
# 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
# 77 - 78        LString(2)    element      Element symbol, right-justified.
# 79 - 80        LString(2)    charge       Charge  on the atom.
setup = """
fields = ((1,6), (7,11), (13, 16), (17,17), (18,20), (22,22), (23,26), (27,27),
          (31,38), (39,46), (47,54), (55,60), (61,66), (77,78), (79,80))
fields = tuple([(i-1,j) for i,j in fields])

line = """
setup += '"' + line + '"\n'

number = 1000000
time = timeit.timeit("list(map(lambda f: line[slice(*f)], fields))", setup, number=number)
print("Map/split: {:.3f} us per match, {:.3f} s per 62875 matches".format(time*1000000/number,
                                                                      time*62875/number))

# (BEST) Struct and binary implementaion - 0.701 us per match, 0.044 s per 62875 matches
setup = """
import struct

fieldwidths = (6, 5, -1, 4, 1, 3, -1, 1, 4, 1, -3, 8, 8, 8, 6, 6, -10, 2, 2 )
fmtstring = ' '.join('{}{}'.format(abs(fw), 'x' if fw < 0 else 's')
                     for fw in fieldwidths)
fieldstruct = struct.Struct(fmtstring)

line = """
setup += '"' + line + '  ".encode()\n'

number = 1000000
time = timeit.timeit("atom_type, number, name, alt_loc, residue_name, chain_id, residue_number, icode, x, y, z, occupancy, B_factor, element, charge=fieldstruct.unpack_from(line)", setup, number=number)
print("struct: {:.3f} us per match, {:.3f} s per 62875 matches".format(time*1000000/number,
                                                                      time*62875/number))

# (2nd BEST) nametuple and binary implementation - 1.042 us per match, 0.066 s per 62875 matches
setup = """
import struct
from collections import namedtuple

fieldwidths = (6, 5, -1, 4, 1, 3, -1, 1, 4, 1, -3, 8, 8, 8, 6, 6, -10, 2, 2 )
fmtstring = ' '.join('{}{}'.format(abs(fw), 'x' if fw < 0 else 's')
                     for fw in fieldwidths)
atom = namedtuple('atom', ['type', 'number', 'name', 'alt_loc', 'residue_name',
                  'chain', 'residue_number', 'icode', 'x', 'y', 'z',
                  'occupancy', 'B_factor', 'element', 'charge'])

line = """
setup += '"' + line + '  ".encode()\n'

number = 1000000
time = timeit.timeit("atom._make(struct.unpack(fmtstring, line))", setup, number=number)
print("namedtuple: {:.3f} us per match, {:.3f} s per 62875 matches".format(time*1000000/number,
                                                                      time*62875/number))