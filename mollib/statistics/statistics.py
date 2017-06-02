"""
Functions to collect statistics on molecules.
"""

import csv
import logging
import json
import os
import tarfile
import io

from setuptools import Command

try:
    from tqdm import tqdm
except ImportError:
    pass

from mollib import Molecule, __path__ as module_path
from mollib.core import settings, load_settings

class Statistics(object):
    """Class to collect statistics on molecules.

    Statistics are compiled first by processing measurements, then by processing
    resulting datasets. Subclasses implement :meth:`process_measurement` and
    :meth:`process_data`, and optionally should set the :attr:`data_path`.

    Attributes
    ----------
    molecule_files: str
        The filenames for the file containing the list of molecule ids to
        collect statistics on.
    data_path: str
        The base path and filename for data files.
    """

    def __init__(self, *filenames):
        # Set the data paths
        self.abs_path = module_path[0]
        self.data_path = os.path.join(self.abs_path, settings.dataset_path)

        # Find all of the files that contain identifiers.
        filenames = set(filenames)
        settings_filenames = [os.path.join(self.data_path, i)
                              for i in settings.model_molecule_identifiers]
        filenames.add(*settings_filenames)

        self.molecule_files = filenames
        self.functions = []

        # Read in the molecule identifiers
        self.read_identifiers()

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

    def process_measurement(self, molecule):
        """Process the measurements for the molecule.

        Parameters
        ----------
        molecule: :obj:`mollib.Molecule` or str
            A molecule or a identifier string for a molecule.

        Returns
        -------
        :obj:`mollib.Molecule`
            (parent class) The molecule to process. Subclasses should call this
            parent class method to get the molecule.
        dict
            (subclasses) Return a dict containing values to collect statistics
            on.
        """
        if isinstance(molecule, str):
            molecule = Molecule(molecule)
        return molecule

    def process_measurements(self):
        """Load all of the molecules and process_measurement each one, saving
        the data to the output measurement file and returning the values in the
        measurement dict.

        Returns
        -------
        measurement_dict: dict
            Collected processed values from :meth:`process_measurement`.
            These are organized as follows:

            - **key**: identifier of the molecule, str
            - **value**: measurement data for that molecule, dict
        """
        msg = "Processing measurements for {}."
        print(msg.format(self.__class__.__name__))

        measurement_dict = self.get_measurement_dict()
        processed_ids = set(measurement_dict.keys())
        to_process_ids = self.identifiers - processed_ids

        tarfile_object = self.open_measurement_tarfile()

        # Create the iterator, optionally using the tqdm counter
        it = enumerate(sorted(to_process_ids), len(processed_ids) + 1)
        if 'tqdm' in globals():
            it = tqdm(it, initial=len(processed_ids),
                      total=len(self.identifiers))

        identifier = None
        try:
            for count, identifier in it:
                # Skip empty entries
                if not identifier.strip():
                    continue

                molecule = Molecule(identifier)

                return_dict = self.process_measurement(molecule)
                self.write_measurement_dict(identifier=identifier,
                                            d=return_dict,
                                            tarfile_object=tarfile_object)
                # Update the measurement_dict
                measurement_dict[identifier] = return_dict

            # Iterate over the identifiers and process_measurement each
        except Exception as e:
            print("Failed at identifier '{}'".format(identifier))
            raise e
        finally:
            tarfile_object.close()
        return measurement_dict

    def get_measurement_filename(self):
        "Return the filename for the measurement data file"
        base = self.data_path

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
                f = tfile.extractfile(member)
                try:
                    identifier = member.name.strip('.json')
                    string = f.read().decode()
                    return_dict = json.loads(string)
                    measurement_dict[identifier] = return_dict
                except:
                    continue
                finally:
                    f.close()
        return measurement_dict

    def write_measurement_dict(self, identifier, d, tarfile_object):
        "Write the measurement dict."
        # Dump the dict into a json string. Note that the dict's keys have
        # to be strings
        string = json.dumps(d)


        # Write the string buffer to the tar file
        tarinfo = tarfile.TarInfo(identifier + '.json')
        tarinfo.size = len(string)
        tarfile_object.addfile(tarinfo=tarinfo,
                               fileobj=io.BytesIO(string.encode()))
        return tarfile_object

    def process_data(self, measurement_dict):
        """Process the data dict.

        Parameters
        ----------
        measurement_dict: dict
            The data dict produced by :meth:`process_measurement`.
        """
        msg = "Processing data for {}."
        print(msg.format(self.__class__.__name__))


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

        # Load the settings to correctly retrieve paths
        load_settings()

        # Disable the logger for the following loop
        logger = logging.getLogger()
        logger.disabled = True

        # Run all of the Statistics submodules
        for Subclass in Statistics.__subclasses__():
            subclass = Subclass()

            result = subclass.process_measurements()
            subclass.process_data(result)
