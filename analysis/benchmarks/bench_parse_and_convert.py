"""
Test different approaches to parsing files and converting types
"""
import timeit

run_tests = {1, 2, 3, 4, 5}

# List comprehension with 'startswith' matching, struct segmentation and zip
# conversion

setup1 = """
import gzip
import re
import io
import struct
from itertools import groupby

contents = gzip.open('1htq.pdb.gz').read().decode('latin-1')

fieldwidths = (6, 5, -1, 4, 1, 3, -1, 1, 4, 1, -3, 8, 8, 8, 6, 6, 10, 2, 2)
fmtstring = ' '.join('{}{}'.format(abs(fw), 'x' if fw < 0 else 's')
                                   for fw in fieldwidths)
fieldstruct = struct.Struct(fmtstring)


conversions = (str.strip,                                # ATOM 
               int,                                      # atom number
               str.strip,                                # atom name
               str.strip,                                # alt loc
               str.strip,                                # residue name
               str.strip,                                # chain id
               int,                                      # residue number
               lambda x: int(x) if x.isdigit() else '',  # residue seq number
               float,                                    # x
               float,                                    # y
               float,                                    # z
               float,
               float,
               str,
               str)
models = []

"""

# Without conversion

stmt1 = """
atom_lines = [line for line in contents.splitlines() 
              if line.startswith('ATOM') or
              line.startswith('HETATM')]
groups = groupby(atom_lines, key=lambda x: x[0:6])

for k, lines in groups:

    if not k.startswith('ATOM') and not k.startswith('HETATM'):
        continue

    model = list(lines)

    models.append(model)
print("No. models: {}".format(len(models)))
"""

if 1 in run_tests:
    print('-' * 80)
    print("1. List comprehension with 'startswith', no conversion.")
    time = timeit.timeit(stmt=stmt1, setup=setup1, number=1)
    print("time: {:.1f} ms/iter".format(time * 1E3))


#  Conversion with Struct and zip
#  Doesn't work in Python3

stmt2 = """
atom_lines = [line for line in contents.splitlines() 
              if line.startswith('ATOM') or
              line.startswith('HETATM')]
groups = groupby(atom_lines, key=lambda x: x[0:6])

for k, lines in groups:

    if not k.startswith('ATOM') and not k.startswith('HETATM'):
        continue

    # flat list
    #model = [x(y) for l in lines 
    #         for x,y in zip(conversions, fieldstruct.unpack_from(l)) ]

    model = [[x(y) for x,y in zip(conversions, fieldstruct.unpack_from(l))]
             for l in lines]

    models.append(model)
print("No. models: {}, No. atoms: {}".format(len(models), len(models[0])))
"""

# List comprehension with 'startswith', manual conversion

if 2 in run_tests:
    print('-' * 80)
    print("2. List comprehension with 'startswith', struct segmentation and "
          "zip conv.")
    time = timeit.timeit(stmt=stmt2, setup=setup1, number=1)
    print("time: {:.1f} ms/iter".format(time * 1E3))

stmt3 = """
atom_lines = [line for line in contents.splitlines() 
              if line.startswith('ATOM') or
              line.startswith('HETATM')]
groups = groupby(atom_lines, key=lambda x: x[0:6])

for k, lines in groups:

    if not k.startswith('ATOM') and not k.startswith('HETATM'):
        continue


    model = [[str(l[0:6]).strip(),
              int(l[7:11]),
              str(l[13:16]).strip(),
              str(l[16:17]).strip(),
              str(l[17:20]).strip(),
              str(l[21:22]),
              int(l[23:26]),
              str(l[27:28]),
              float(l[31:38]),
              float(l[39:46]),
              float(l[47:54]),
              float(l[55:60]),
              float(l[61:66]),
              str(l[76:78]),
              '',
              ]
             for l in lines]

    models.append(model)
print("No. models: {}, No. atoms: {}".format(len(models), len(models[0])))
"""

if 3 in run_tests:
    print('-' * 80)
    print("3. List comprehension with 'startswith', manual conversion.")
    time = timeit.timeit(stmt=stmt3, setup=setup1, number=1)
    print("time: {:.1f} ms/iter".format(time * 1E3))


# MP (Queue) List comprehension with 'startswith', manual conversion


setup4 = """
import gzip
import re
import io
import multiprocessing as mp
from itertools import groupby

contents = gzip.open('1htq.pdb.gz').read().decode('latin-1')

_re_model = re.compile(r'^MODEL\s+(?P<model_id>[\d]+)\s+$')

line_chunks = 1000

# job queue
queue = mp.JoinableQueue()

# result queue
result = mp.Queue()

def process_measurement(model_id, lines):

    processed = [[str(l[0:6]).strip(),
                  int(l[7:11]),
                  str(l[13:16]).strip(),
                  str(l[16:17]).strip(),
                  str(l[17:20]).strip(),
                  str(l[21:22]),
                  int(l[23:26]),
                  str(l[27:28]),
                  float(l[31:38]),
                  float(l[39:46]),
                  float(l[47:54]),
                  float(l[55:60]),
                  float(l[61:66]),
                  str(l[76:78]),
                  '',
                ]
                for l in lines]

    return model_id, processed

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

stmt4 = """
atom_lines = [line for line in contents.splitlines() 
              if line.startswith('ATOM') or
              line.startswith('HETATM') or
              line.startswith('MODEL') or
              line.startswith('CONECT')]
              
groups = groupby(atom_lines, key=lambda x: x[0:6])

###

count = 0
current_model_id = None
for k, lines in groups:

    if k.startswith('MODEL'):
        match = _re_model.match(next(lines))
        if match:
            d = match.groupdict()
            current_model_id = int(d['model_id'])
        continue
    
    if not k.startswith('ATOM') and not k.startswith('HETATM'):
        continue

    # Chunk the groups
    atom_lines = list(lines)
    for i in range(0, len(atom_lines), line_chunks):
        chunk = atom_lines[i:i + line_chunks]
        queue.put((current_model_id, chunk))
        count += 1

for i in range(num_workers):
    queue.put(None)


# collect the results
queue.join()

models = {}
while count > 0:
    rv = result.get()
    model_id, processed = rv
    count -= 1
    
    model_list = models.setdefault(model_id, [])
    model_list += processed
    

print("No. models: {}, No. atoms: {}".format(len(models), len(models[1])))
"""

if 4 in run_tests:
    print('-' * 80)
    print("4. MP(Queue) List comprehension with 'startswith', manual "
          "conversion.")
    time = timeit.timeit(stmt=stmt4, setup=setup4, number=1)
    print("time: {:.1f} ms/iter".format(time * 1E3))


# MP using a Manager


setup5 = """
import gzip
import re
import io
import multiprocessing as mp
from itertools import groupby

contents = gzip.open('1htq.pdb.gz').read().decode('latin-1')

_re_model = re.compile(r'^MODEL\s+(?P<model_id>[\d]+)\s+$')

line_chunks = 10000


def process_measurement(atom_dict, key):
    lines = atom_dict[key]
    
    processed = [[str(l[0:6]).strip(),
                  int(l[7:11]),
                  str(l[13:16]).strip(),
                  str(l[16:17]).strip(),
                  str(l[17:20]).strip(),
                  str(l[21:22]),
                  int(l[23:26]),
                  str(l[27:28]),
                  float(l[31:38]),
                  float(l[39:46]),
                  float(l[47:54]),
                  float(l[55:60]),
                  float(l[61:66]),
                  str(l[76:78]),
                  '',
                ]
                for l in lines]
    print(processed[0])
    atom_dict[key] = processed


mgr = mp.Manager()
atom_dict = mgr.dict()


"""

stmt5 = """
atom_lines = [line for line in contents.splitlines() 
              if line.startswith('ATOM') or
              line.startswith('HETATM') or
              line.startswith('MODEL') or
              line.startswith('CONECT')]

groups = groupby(atom_lines, key=lambda x: x[0:6])

jobs = []

current_model_id = None

for k, lines in groups:

    if k.startswith('MODEL'):
        match = _re_model.match(next(lines))
        if match:
            d = match.groupdict()
            current_model_id = int(d['model_id'])
        continue

    if not k.startswith('ATOM') and not k.startswith('HETATM'):
        continue

    # Chunk the groups
    atom_lines = list(lines)
    for i in range(0, len(atom_lines), line_chunks):
        chunk = atom_lines[i:i + line_chunks]
        key = (current_model_id, k, i)

        atom_dict[key] = chunk
        
        jobs.append(mp.Process(target=process_measurement, 
                    args=(atom_dict, (current_model_id, k, i))))

        

for j in jobs:
    j.start()
for j in jobs:
    j.join()

models = {}
for k,v in atom_dict.items():

    model_id, line_type, chunk_num = k
    l = models.setdefault(model_id, [])
    l += v

print("No. models: {}, No. atoms: {}".format(len(models), len(models[1])))
"""

if 5 in run_tests:
    print('-' * 80)
    print("5. MP(Manager) List comprehension with 'startswith', manual "
          "conversion.")
    time = timeit.timeit(stmt=stmt5, setup=setup5, number=1)
    print("time: {:.1f} ms/iter".format(time * 1E3))

print('-'*80)
