"""
Functions to collect statistics on molecules.
"""
import shelve
import csv
import logging

from mollib import Molecule
from . import settings


class Statistics(object):
    """Class to collect statistics on molecules.

    Attributes
    ----------
    molecule_files: str
        The filenames for the file containing the list of molecule ids to
        collect statistics on.
    measurement_file_base: str
        The base path and filename for the output measurement data files.
    output_ext: str
        The extension for the output data files.
    """

    measurement_file_base = settings.output_measurement_file_base

    def __init__(self, *filenames):
        filenames = set(*filenames)
        filenames.add(*settings.input_molecule_files)

        self.molecule_files = filenames
        self.functions = []

        # Read in the molecule identifiers
        self.read_identifiers()

    def process(self, molecule):
        """Process the molecule.

        Parameters
        ----------
        molecule: :obj:`mollib.Molecule` or str
            A molecule or a identifier string for a molecule.

        Returns
        -------
        :obj:`Molecule`
            (parent class) The molecule to process. Subclasses should call this
            parent class method to get the molecule.
        dict
            (subclasses) Return a dict containing values to collect statistics
            on.
        """
        if isinstance(molecule, str):
            molecule = Molecule(molecule)
        print(molecule.name)
        #logging.info('Processing measurements statistics '
        #             'for {}'.format(molecule.name))
        return molecule

    def process_measurements(self):
        """Load all of the molecules and process each one, saving the data to
        the output measurement file.

        Results are stored

        Returns
        -------

        """
        # Open output data file
        try:
            db = self.get_measurement_file()

            # Find the identifiers that need to be processed
            processed_ids = set(db.keys())
            to_process_ids = processed_ids ^ self.identifiers

            # Iterate over the identifiers and process each
            for identifier in sorted(to_process_ids):
                molecule = Molecule(identifier)
                return_dict = self.process(molecule)
                db[identifier] = return_dict
        except Exception as e:
            raise e
        finally:
            db.close()

    def read_identifiers(self):
        """Read in the list of identifiers from the molecule files."""
        self.identifiers = set()

        for filename in self.molecule_files:
            try:
                f = open(filename, 'r')
            except IOError:
                logging.warning("Statistics file '{}' not "
                                "found.".format(filename))
                continue

            with f:
                reader = csv.reader(f, skipinitialspace=True)
                for items in reader:
                    if len(items) > 0:
                        self.identifiers.update(items)

    def get_measurement_filename(self):
        "Return the filename for the measurement output file."
        base = self.measurement_file_base
        class_name = self.__class__.__name__.lower()
        return base + class_name

    def get_measurement_file(self):
        "Return the dict object for the measurement file."
        filename = self.get_measurement_filename()
        db = shelve.open(filename)
        return db



import mollib.hbonds
import mollib.hydrogens

class RamachandranStatistics(Statistics):
    """Collect statistics on Ramachandran angles.
    """

    def process(self, molecule):
        """Process the molecule and return the Ramachandran angle statistics.

        Returns
        -------
        dict
            A dict with the following keys/values:

            - key: (str) secondary structure classification
            - value: (list of tuples) listing of Ramachandran angles (phi, psi)
        """
        molecule = super(RamachandranStatistics, self).process(molecule)

        # Hydrogenate the molecule
        mollib.hydrogens.add_hydrogens(molecule)

        # Get the hydrogen bonds for the molecule
        hbonds = mollib.hbonds.find_hbond_partners(molecule)

        # Group the hbonds into (chain.id, residue.number) and
        #  major_classifications
        residue_classes = dict()

        for hbond in hbonds:
            try:
                donor_residue = hbond.donor.atom2.residue
                acceptor_residue = hbond.acceptor.atom2.residue
            except AttributeError:
                continue
            major_classification = hbond.major_classification
            minor_classification = hbond.minor_classification

            if major_classification == mollib.hbonds.settings.major_bb_bb_amide:
                key = (donor_residue.chain.id, donor_residue.number)
                residue_classes[key] = minor_classification

                key = (acceptor_residue.chain.id, acceptor_residue.number)
                residue_classes[key] = minor_classification

        # Get the Ramachandran angles and group them.
        return_dict = dict()

        for residue in molecule.residues:
            phi_psi = residue.ramachandran_angles
            key = (residue.chain.id, residue.number)
            group = residue_classes.get(key, 'Other')

            return_dict.setdefault(group, list()).append(phi_psi)

        return return_dict

