"""
Test different approaches to parsing files using regexes
"""
import timeit

# Reference, one-by-one

setup1 = """
import gzip
import re
import io

f = io.BufferedReader(gzip.open('1htq.pdb.gz'))
re_atom = re.compile(r"(?P<type>ATOM  |HETATM)"
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
                           "(?P<charge>[\d\s\.\-]{2})?")
"""

stmt1 = """
count = 0
for line in f:
    if type(line) == bytes:
        line = line.decode('latin-1')
    m = re_atom.match(line)
    if m:
        count += 1
print("matches:", count)
"""

print('-'*80)
print("Reference Implementation")
time = timeit.timeit(stmt=stmt1, setup=setup1, number=1)
print("time: {:.1f} ms/iter".format(time * 1E3))

# Findall implementation

setup2 = """
import gzip
import re
import io

contents = gzip.open('1htq.pdb.gz').read().decode('latin-1')
re_atom = re.compile(r"(?P<type>ATOM  |HETATM)"
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
                           "(?P<charge>[\d\s\.\-]{2})?")
count = 0
"""

stmt2 = """
count = 0
for m in re_atom.finditer(contents):
    count += 1
print("matches:", count)
"""

print('-'*80)
print("Finditer Implementation")
time = timeit.timeit(stmt=stmt2, setup=setup2, number=1)
print("time: {:.1f} ms/iter".format(time * 1E3))

# Findall with simpler regex

setup3 = """
import gzip
import re
import io

contents = gzip.open('1htq.pdb.gz').read().decode('latin-1')
re_atom = re.compile(r"^ATOM.*")
count = 0
"""

stmt3 = """
count = 0
for m in re_atom.finditer(contents):
    count += 1
print("matches:", count)
"""

print('-'*80)
print("Finditer Implementation with simpler regex")
time = timeit.timeit(stmt=stmt3, setup=setup3, number=1)
print("time: {:.1f} ms/iter".format(time * 1E3))

# List comprehension with 'startswith' matching

setup4 = """
import gzip
import re
import io

contents = gzip.open('1htq.pdb.gz').read().decode('latin-1')
re_atom = re.compile(r"(?P<type>ATOM  |HETATM)"
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
                           "(?P<charge>[\d\s\.\-]{2})?")
count = 0
"""

stmt4 = """
count = 0
# print(dir(contents))
atom_lines = [line for line in contents.splitlines() 
              if line.startswith('ATOM') or
              line.startswith('HETATM') or
              line.startswith('MODEL') or
              line.startswith('CONECT')]
count = len(atom_lines)
print("matches:", count)
"""

print('-'*80)
print("List comprehension with 'startswith' matching")
time = timeit.timeit(stmt=stmt4, setup=setup4, number=1)
print("time: {:.1f} ms/iter".format(time * 1E3))

# List comprehension with 'startswith' matching and grouping

setup5 = """
import gzip
import re
import io
from itertools import groupby

contents = gzip.open('1htq.pdb.gz').read().decode('latin-1')
re_atom = re.compile(r"(?P<type>ATOM  |HETATM)"
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
                           "(?P<charge>[\d\s\.\-]{2})?")
count = 0
"""

stmt5 = """
count = 0
# print(dir(contents))
atom_lines = [line for line in contents.splitlines() 
              if line.startswith('ATOM') or
              line.startswith('HETATM') or
              line.startswith('MODEL') or
              line.startswith('CONECT')]
groups = groupby(atom_lines, key=lambda x: x[0:7])
count = sum(len(list(l)) for k,l in groups if k.startswith('ATOM'))
print("atoms:", count)
"""

print('-'*80)
print("List comprehension with 'startswith' matching and grouping")
time = timeit.timeit(stmt=stmt5, setup=setup5, number=1)
print("time: {:.1f} ms/iter".format(time * 1E3))

# List comprehension with 'startswith' matching and manual grouping

setup6 = """
import gzip
import re
import io
from itertools import groupby

contents = gzip.open('1htq.pdb.gz').read().decode('latin-1')
re_atom = re.compile(r"(?P<type>ATOM  |HETATM)"
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
                           "(?P<charge>[\d\s\.\-]{2})?")
count = 0
"""

stmt6 = """
count = 0
atom_lines = [line for line in contents.splitlines() 
              if line.startswith('ATOM') or
              line.startswith('HETATM') or
              line.startswith('MODEL') or
              line.startswith('CONECT')]
              
models = []
current_molecule = []
models.append(current_molecule)
connections = []
for line in contents.splitlines():
    if line.startswith('MODEL') and len(current_molecule) > 0:
        current_molecule = []
        models.append(current_molecule)
    elif line.startswith('ATOM') or line.startswith('HETATM'):
        current_molecule.append(line)
    elif line.startswith('CONECT'):
        connections.append(line)
        
print("atoms:", len(models))
"""

print('-'*80)
print("List comprehension with 'startswith' matching and manual grouping")
time = timeit.timeit(stmt=stmt6, setup=setup6, number=1)
print("time: {:.1f} ms/iter".format(time * 1E3))


print('-'*80)

# ATOM: (6, 5, -1, 4, 1, 3, -1, 1, 4, 1, -3, 8, 8, 8, 6, 6, 10, 2, 2)
#  def split(string, fieldwidths):
# ...     fmtstring = ' '.join('{}{}'.format(abs(fw), 'x' if fw < 0 else 's')
# ...                             for fw in fieldwidths)
# ...     fieldstruct = struct.Struct(fmtstring)
# ...     parse = fieldstruct.unpack_from
# ...     fields = parse(string)
# ...     print('fields: {}'.format(fields))
#'6s 5s 1x 4s 1s 3s 1x 1s 4s 1s 3x 8s 8s 8s 6s 6s 10s 2s 2s'