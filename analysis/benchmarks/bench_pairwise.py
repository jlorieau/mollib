"""
Bench pairwise searches.
"""
import timeit
from mollib import Molecule

base_setup = """
from mollib import Molecule

mol = Molecule('{pdb}')
point_list = [a.pos for a in mol.atoms]
"""

### Linear Search Implementation ###

setup_linear = """
from math import sqrt

point = mol['A'][10]['N'].pos

point_list = [a.pos for a in mol.atoms]
"""

stmt_linear = """
a = [p for p in point_list if
     sqrt((p[0]-point[0])**2 + (p[1]-point[1])**2 + (p[2]-point[2])**2) 
     < {dist}]
"""

### KD-tree ###

setup_kd = """
from scipy import spatial

point = mol['A'][10]['N'].pos
atom_nums = {{count:a for count, a in enumerate(mol.atoms)}}
tree = spatial.cKDTree(point_list)
"""

stmt_kd = """
a_nums = tree.query_ball_point(point, {dist})
a = [atom_nums[i] for i in a_nums] 
"""

### Geometry Box Implementation ###

setup_box = """
from mollib.core import within_distance
"""

stmt_box = """
# Find all atoms within 5A of an atom
a = within_distance(mol['A'][10]['N'], cutoff={dist})
"""

for pdb_id, dist, N in (('2KXA', 5.0, 5000),
                        ('2MJB', 5.0, 2500),
                        ('3H0G', 5.0, 10),):
    mol = Molecule(pdb_id)

    print('\n')
    print('-' * 80)
    print("Molecule: ('{}') {} atoms.".format(pdb_id, mol.atom_size))
    print('-' * 80)

    # Linear search
    setup = '\n'.join((base_setup, setup_linear))
    time = timeit.timeit(stmt=stmt_linear.format(dist=dist),
                         setup=setup.format(pdb=pdb_id),
                         number=N)
    print("Linear Search: {:.3f} iter/s".format(N / time))

    # KD-tree
    setup = '\n'.join((base_setup, setup_kd))
    time = timeit.timeit(stmt=stmt_kd.format(dist=dist),
                         setup=setup.format(pdb=pdb_id),
                         number=N)
    print("KD-tree Search: {:.3f} iter/s".format(N / time))

    # Geometry box
    setup = '\n'.join((base_setup, setup_box))
    time = timeit.timeit(stmt=stmt_box.format(dist=dist),
                         setup=setup.format(pdb=pdb_id),
                         number=N)
    print("Geometry box: {:.3f} iter/s".format(N / time))
