import timeit

setup="""
import re
re_int = re.compile(r'-?\d+')

value = " 12312416587 "
"""

# Direct conversion, single value - Direct (int): 0.612us per iteration
number = 1000000
time = timeit.timeit("int(value)", setup, number=number)
print("Direct (int): {:.3f}us per iteration.".format(time*1000000./number))

# Regex conversion, single value - Regex (int): 1.385us per iteration
number = 1000000
time = timeit.timeit("int(re_int.search(value).group())", setup, number=number)
print("Regex (int): {:.3f}us per iteration.".format(time*1000000./number))

# converting a struct item-by-item - 1.562 us per iteration, 0.098 s per 62875 matches
setup="""
import struct

fieldwidths = (6, 5, -1, 4, 1, 3, -1, 1, 4, 1, -3, 8, 8, 8, 6, 6, -10, 2, 2 )
fmtstring = ' '.join('{}{}'.format(abs(fw), 'x' if fw < 0 else 's')
                     for fw in fieldwidths)
fieldstruct = struct.Struct(fmtstring)

line = "ATOM   1617  C   SER A 244      21.367   6.240  37.294  1.00 19.56           C  "

atom_type, number, name, alt_loc, residue_name, chain_id, residue_number, icode, x, y, z, occupancy, B_factor, element, charge=fieldstruct.unpack_from(line)
"""

method="""
atom_type = atom_type.strip()
number = int(number)
name = name.strip()
alt_loc = alt_loc.strip()
residue_name = residue_name.strip()
chain_id = chain_id.strip()
residue_number = int(residue_number)
icode = icode.strip()
x = float(x)
y = float(y)
z = float(z)
occupancy = float(occupancy)
B_factor = float(B_factor)
element = element.strip()
#charge = float(charge) if charge.strip() else None
"""

number = 1000000
time = timeit.timeit(method, setup, number=number)
print("struct one-by-one: {:.3f} us per iteration, {:.3f} s per 62875 matches".format(time*1000000/number,
                                                                                     time*62875/number))

# converting a namedtuple with properties - 12.882 us per iteration, 0.810 s per 62875 matches
setup="""
import struct
from collections import namedtuple

fieldwidths = (6, 5, -1, 4, 1, 3, -1, 1, 4, 1, -3, 8, 8, 8, 6, 6, -10, 2, 2 )
fmtstring = ' '.join('{}{}'.format(abs(fw), 'x' if fw < 0 else 's')
                     for fw in fieldwidths)
fieldstruct = struct.Struct(fmtstring)

line = "ATOM   1617  C   SER A 244      21.367   6.240  37.294  1.00 19.56           C  "

class atom(namedtuple('atom',['atype', 'anumber', 'aname', 'aalt_loc', 'aresidue_name',
                  'achain', 'aresidue_number', 'aicode', 'ax', 'ay', 'az',
                  'aoccupancy', 'aB_factor', 'aelement', 'acharge'])):
    @property
    def type(self): return self.atype.strip()

    @property
    def number(self): return int(self.anumber)

    @property
    def name(self): return self.aname.strip()

    @property
    def alt_loc(self): return self.aalt_loc.strip()

    @property
    def residue_name(self): return self.aresidue_name.strip()

    @property
    def chain(self): return self.achain.strip()

    @property
    def residue_number(self): return int(self.aresidue_number)

    @property
    def icode(self): return self.aicode.strip()

    @property
    def x(self): return float(self.ax)

    @property
    def y(self): return float(self.ay)

    @property
    def z(self): return float(self.az)

    @property
    def occupancy(self): return float(self.aoccupancy)

    @property
    def B_factor(self): return float(self.aB_factor)

    @property
    def charge(self): return float(self.acharge) if self.acharge.strip() else ''


a = atom(*fieldstruct.unpack_from(line))
"""

method = """
[getattr(a, attr) for attr in ('type', 'number', 'name', 'alt_loc',
                               'residue_name', 'chain', 'residue_number',
                               'icode', 'x', 'y', 'z', 'occupancy',
                               'B_factor', 'charge')]
"""

number = 100000
time = timeit.timeit(method, setup, number=number)
print("namedtuple properties: {:.3f} us per iteration, {:.3f} s per 62875 matches".format(time*1000000/number,
                                                                                     time*62875/number))