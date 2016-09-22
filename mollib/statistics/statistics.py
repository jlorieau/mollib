"""
Functions to collect statistics on molecules.
"""
import csv
import logging
import json
import os
import tarfile
import io
import sys

from setuptools import Command

from mollib import Molecule
from mollib.core import settings


class Statistics(object):
    """Class to collect statistics on molecules.

    Attributes
    ----------
    molecule_files: str
        The filenames for the file containing the list of molecule ids to
        collect statistics on.
    data_path: str
        The base path and filename for data files.
    output_ext: str
        The extension for the output data files.
    """

    data_path = settings.dataset_path

    def __init__(self, *filenames):
        filenames = set(*filenames)
        filenames.add(*settings.model_molecule_identifiers)

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
        measured_dict = self.get_measurement_dict()
        processed_ids = set(measured_dict.keys())
        to_process_ids = processed_ids ^ self.identifiers

        tarfile_object = self.open_measurement_tarfile()

        try:

            for count, identifier in enumerate(sorted(to_process_ids),
                                               len(processed_ids) + 1):
                # Skip empty entries
                if not identifier.strip():
                    continue

                print('{}'.format(count).rjust(4) + '.' +
                     ' {}'.format(identifier))
                molecule = Molecule(identifier)

                return_dict = self.process(molecule)
                self.write_measurement_dict(identifier=identifier,
                                            d=return_dict,
                                            tarfile_object=tarfile_object)
            # Iterate over the identifiers and process each
        except Exception as e:
            raise e
        finally:
            tarfile_object.close()

    def read_identifiers(self):
        """Read in the list of identifiers from the molecule files."""
        self.identifiers = set()

        for filename in self.molecule_files:
            try:
                f = open(filename, 'r')
            except IOError:
                logging.warning("Input identifiers file '{}' not "
                                "found.".format(filename))
                continue

            with f:
                reader = csv.reader(f, skipinitialspace=True)
                for items in reader:
                    if len(items) > 0:
                        self.identifiers.update(items)

    def get_measurement_filename(self):
        "Return the filename for the measurement data file"
        base = self.data_path
        class_name = self.__class__.__name__.lower()
        base = os.path.join(base, class_name)

        # Create the directory, if needed
        if not os.path.isdir(base):
            os.makedirs(base)

        filename = os.path.join(base, 'measurements.tar')
        return filename

    def open_measurement_tarfile(self, mode='a'):
        "Open the measurement tarfile."
        filename = self.get_measurement_filename()
        # Create the file if needed
        try:
            tfile = tarfile.open(name=filename, mode=mode)
        except (IOError, ValueError, tarfile.ReadError):
            tfile = tarfile.open(name=filename, mode='w')
            tfile.close()
            tfile = tarfile.open(name=filename, mode=mode)

        return tfile

    def get_measurement_dict(self):
        "Return the dict object for the measurement file."
        tfile = self.open_measurement_tarfile(mode='r')
        measurement_dict = {}

        with tfile:
            # Extract the files
            for member in tfile.getmembers():
                with tfile.extractfile(member) as f:
                    try:
                        identifier = member.name.strip('.json')
                        string = f.read().decode()
                        return_dict = json.loads(string)
                        measurement_dict[identifier] = return_dict
                    except:
                        continue
        return measurement_dict

    def write_measurement_dict(self, identifier, d, tarfile_object):
        "Write the measurement dict."
        try:
            string = json.dumps(d)
        except:
            return None

        # Write the string buffer to the tar file
        tarinfo = tarfile.TarInfo(identifier + '.json')
        tarinfo.size = len(string)
        tarfile_object.addfile(tarinfo=tarinfo,
                               fileobj=io.BytesIO(string.encode()))
        return tarfile_object

    def process_data(self, data_dict):
        pass


import mollib.hbonds
import mollib.hydrogens

class RamachandranStatistics(Statistics):
    """Collect statistics on Ramachandran angles.
    """

    data_path = settings.dataset_path

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
            group = residue_classes.get(key, 'No hydrogen bonds')

            # None values are not saved
            if all(i is not None for i in phi_psi):
                return_dict.setdefault(group, list()).append(phi_psi)

        return return_dict

class BuildData(Command):
    """Build the Statistics data"""

    description = "Build the statistics data"

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        "Process the measurements of all subclasses"
        for Subclass in Statistics.__subclasses__():
            subclass = Subclass()

            msg = "Processing data for '{}'"
            print(msg.format(subclass.__class__.__name__))

            subclass.process_measurements()
