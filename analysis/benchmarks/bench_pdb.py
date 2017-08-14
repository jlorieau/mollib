import timeit

# Tested on Mac OS X 2.6 GHz Intel Core i7 (4 core) w 16MB RAM
# PDB file distribution: https://plot.ly/~fomightez/0/file-size-distribution-of-all-118362-protein-data-bank-entries-as-of-may-15-2016.embed
#
# Summary. Most of the PDB (1.5%) consists of files below 1MB.
#     - The chunk implementation is just as fast as the non-chunk version for
#       all files (within +/-15%).
#     - Single process_measurement, chunk is faster by 15-20% for files <1MB, and it is
#       slower by 15% for files around 9MB.
#     - Multitprocess, chunk is slower by 10% for files smaller than chunk, and
#       it is 2-3 faster for files larger than chunk size.
#     - Conclusion: The chunk approach is faster for files smaller than the
#       chunk size if mp is avoided.
#     - Conclusion: The chunk approach is faster for files larger than the
#       chunk size if mp is used.
#     - Conclusion: A chunk size of 512KB appears to be ideal for all files. A
#       chunk size of 256KB significantly slows down larger files (3H0G(.
#
# 3H0G (9MB file, 0.01% of PDB) - 1MB chunks
# (5) chunks & type conv.                 : 1199.39ms per iter. (14% slower than (4))
# (6) chunks/mp & type conv.              : 450.54ms per iter. (43% faster than (4))
# (7) chunks/mp/mm & type conv.           : 473.47ms per iter.
# ------------------------------------------------------------------------------
# (4) generator & type conv.              : 1051.98ms per iter.
# (3) readline, line skip & type conv.    : 1318.99ms per iter.
# (2) readline & type conv.               : 1428.30ms per iter.
# (1) readline/no type conv.              : 1020.67ms per iter.

# 3H0G (9MB file, 0.01% of PDB) - 512KB chunks
# (5) chunks & type conv.                 : 1513.29ms per iter. (26% slower than 1MB chunks)
# (6) chunks/mp & type conv.              : 504.27ms per iter. (50% faster than (4))
# (7) chunks/mp/mm & type conv.           : 480.31ms per iter.
# ------------------------------------------------------------------------------
# (4) generator & type conv.              : 1002.89ms per iter.
# (3) readline, line skip & type conv.    : 1205.19ms per iter.
# (2) readline & type conv.               : 1177.30ms per iter.
# (1) readline/no type conv.              : 853.88ms per iter.

# 5CJP (1MB file, 1.5% of PDB) - 1MB chunks
# (5) chunks & type conv.                 : 127.30ms per iter.
# (6) chunks/mp & type conv.              : 162.44ms per iter.
# (7) chunks/mp/mm & type conv.           : 161.19ms per iter.
# ------------------------------------------------------------------------------
# (4) generator & type conv.              : 152.97ms per iter.
# (3) readline, line skip & type conv.    : 176.38ms per iter.
# (2) readline & type conv.               : 167.50ms per iter.
# (1) readline/no type conv.              : 128.24ms per iter.

# 5CJP (1MB file, 1.5% of PDB) - 512KB chunks
# (5) chunks & type conv.                 : 126.36ms per iter.
# (6) chunks/mp & type conv.              : 91.92ms per iter. (40% faster than (4))
# (7) chunks/mp/mm & type conv.           : 97.54ms per iter.
# ------------------------------------------------------------------------------
# (4) generator & type conv.              : 149.36ms per iter.
# (3) readline, line skip & type conv.    : 178.79ms per iter.
# (2) readline & type conv.               : 169.88ms per iter.
# (1) readline/no type conv.              : 117.83ms per iter.

# 2KXA - 1MB chunks
# (5) chunks & type conv.                 : 43.12ms per iter.
# (6) chunks/mp & type conv.              : 55.06ms per iter. (28% slower than (5))
# (7) chunks/mp/mm & type conv.           : 54.66ms per iter.
# ------------------------------------------------------------------------------
# (4) generator & type conv.              : 47.14ms per iter.
# (3) readline, line skip & type conv.    : 51.09ms per iter.
# (2) readline & type conv.               : 50.88ms per iter.
# (1) readline/no type conv.              : 38.60ms per iter.

# # 2KXA - 512KB chunks
# (5) chunks & type conv.                 : 36.13ms per iter.
# (6) chunks/mp & type conv.              : 48.20ms per iter.
# (7) chunks/mp/mm & type conv.           : 65.78ms per iter.
# ------------------------------------------------------------------------------
# (4) generator & type conv.              : 44.06ms per iter.
# (3) readline, line skip & type conv.    : 62.75ms per iter.
# (2) readline & type conv.               : 55.07ms per iter.
# (1) readline/no type conv.              : 39.29ms per iter.

setup="""
import re
import gzip
import struct

fieldwidths = (6, 5, -1, 4, 1, 3, -1, 1, 4, 1, -3, 8, 8, 8, 6, 6, -10, 2, 2 )
fmtstring = ' '.join('{}{}'.format(abs(fw), 'x' if fw < 0 else 's')
                     for fw in fieldwidths)
fieldstruct = struct.Struct(fmtstring)

# regex to match atom lines
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
"""

msg = "{:40s}: {:.2f}ms per iter."

# Implementation #5
setup5 = setup + """
def getchunks(file, size=512*1024):
    f = gzip.open(file)
    while 1:
        start = f.tell()
        f.seek(size, 1)
        s = f.readline()
        yield start, f.tell() - start
        if not s:
            break

def process_measurement(file, chunk):
    f = gzip.open(file)
    f.seek(chunk[0])

    parsed = []
    for line in f.read(chunk[1]).splitlines():
        if "ATOM" in line or "HETATM" in line:
            m = re_atom.match(line)
            if m:
                d = m.groupdict()
                parsed.append((
                str(d['type']),
                int(d['number']),
                str(d['name']),
                str(d['alt_loc']),
                str(d['residue_name']),
                str(d['chain']),
                int(d['residue_number']),
                str(d['icode']),
                float(d['x']),
                float(d['y']),
                float(d['z']),
                float(d['occupancy']),
                str(d['element'])))
    return parsed

"""

implementation5 = """
parsed = []
for chunk in getchunks('2KXA.pdb.gz'):
    parsed += process_measurement('2KXA.pdb.gz', chunk)
print(len(parsed))
"""
number = 1
time = timeit.timeit(implementation5, setup5, number=number)
print(msg.format('(5) chunks & type conv.', time*1000./number, time,
                 number))

# Implementation #6

setup6 = setup5 + """
import multiprocessing as mp

# job queue
queue = mp.JoinableQueue()

# result queue
result = mp.Queue()

def process_queue(q, result):
    while True:
        rv = q.get()
        if rv is None:
          q.task_done()
          break
        result.put(process_measurement(rv[0], rv[1]))
        q.task_done()

num_workers = mp.cpu_count()

workers = [mp.Process(target=process_queue, args=(queue, result))
           for i in range(num_workers)]
for w in workers:
    w.start()
"""

# Implementation 6
implementation6 = """
count = 0
for chunk in getchunks('2KXA.pdb.gz'):
    queue.put(('2KXA.pdb.gz', chunk))
    count += 1

for i in range(num_workers):
    queue.put(None)

queue.join()

parsed = []
while count > 0:
    parsed += result.get()
    count -= 1

print(len(parsed))
"""
number = 1
time = timeit.timeit(implementation6, setup6, number=number)
print(msg.format('(6) chunks/mp & type conv.', time*1000./number, time,
                 number))


# Implementation 7
setup7 = setup6 + """
import mmap

def process_measurement(file, chunk):
    fa = gzip.open(file)
    f = mmap.mmap(
            fa.fileno(),
            os.path.getsize(file),
            access=mmap.ACCESS_READ
        )
    f.seek(chunk[0])

    parsed = []
    for line in f.read(chunk[1]).splitlines():
        if "ATOM" in line or "HETATM" in line:
            m = re_atom.match(line)
            if m:
                d = m.groupdict()
                parsed.append((
                str(d['type']),
                int(d['number']),
                str(d['name']),
                str(d['alt_loc']),
                str(d['residue_name']),
                str(d['chain']),
                int(d['residue_number']),
                str(d['icode']),
                float(d['x']),
                float(d['y']),
                float(d['z']),
                float(d['occupancy']),
                str(d['element'])))
    return parsed
"""

number = 1
time = timeit.timeit(implementation6, setup7, number=number)
print(msg.format('(7) chunks/mp/mm & type conv.', time*1000./number, time,
                 number))

print "-"*80

# Implementation #4
implementation4 = """
parsed = []
with gzip.open('2KXA.pdb.gz', 'rb') as f:
    matches = (re_atom.match(line)
               for line in f if "ATOM" in line or "HETATM" in line)
    mapp    = (match.groupdict() for match in matches if match)

    for d in mapp:
        parsed.append((
        str(d['type']),
        int(d['number']),
        str(d['name']),
        str(d['alt_loc']),
        str(d['residue_name']),
        str(d['chain']),
        int(d['residue_number']),
        str(d['icode']),
        float(d['x']),
        float(d['y']),
        float(d['z']),
        float(d['occupancy']),
        str(d['element'])))
print(len(parsed))
"""
number = 1
time = timeit.timeit(implementation4, setup, number=number)
print(msg.format('(4) generator & type conv.', time*1000./number, time,
                 number))


# Implementation #3
implementation3 = """
parsed = []
with gzip.open('2KXA.pdb.gz', 'rb') as f:
    for line in f.readlines():
        line = line.decode('latin-1')
        if "ATOM" not in line and "HETATM" not in line:
            continue
        m = re_atom.match(line)
        if m:
            d = m.groupdict()
            parsed.append((
            str(d['type']),
            int(d['number']),
            str(d['name']),
            str(d['alt_loc']),
            str(d['residue_name']),
            str(d['chain']),
            int(d['residue_number']),
            str(d['icode']),
            float(d['x']),
            float(d['y']),
            float(d['z']),
            float(d['occupancy']),
            str(d['element'])))
print(len(parsed))
"""
number = 1
time = timeit.timeit(implementation3, setup, number=number)
print(msg.format('(3) readline, line skip & type conv.', time*1000./number, time,
                 number))


# Implementation #2.

implementation2 = """
parsed = []
with gzip.open('2KXA.pdb.gz', 'rb') as f:
    for line in f.readlines():
        line = line.decode('latin-1')
        m = re_atom.match(line)
        if m:
            d = m.groupdict()
            parsed.append((
            str(d['type']),
            int(d['number']),
            str(d['name']),
            str(d['alt_loc']),
            str(d['residue_name']),
            str(d['chain']),
            int(d['residue_number']),
            str(d['icode']),
            float(d['x']),
            float(d['y']),
            float(d['z']),
            float(d['occupancy']),
            str(d['element'])))
print(len(parsed))
"""
number = 1
time = timeit.timeit(implementation2, setup, number=number)
print(msg.format('(2) readline & type conv.', time*1000./number, time,
                 number))


# Implementation #1.

implementation1 = """
number = []
with gzip.open('2KXA.pdb.gz', 'rb') as f:
    for line in f.readlines():
        line = line.decode('latin-1')
        m = re_atom.match(line)
        if m:
            number.append(int(m.groupdict()['number']))
print(len(number))
"""
number = 1
time = timeit.timeit(implementation1, setup, number=number)
print(msg.format('(1) readline/no type conv.', time*1000./number, time,
                 number))
