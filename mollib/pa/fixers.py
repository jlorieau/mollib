"""
Tools to fix mistakes and errors in partial alignment (RDC and RACS) data.
"""
from copy import deepcopy

from mollib import Molecule
from mollib.utils.interactions import interaction_type
from .process_molecule import Process
from .svd import calc_pa_SVD
from .analysis import find_outliers
from . import settings


class Fixer(object):
    """Fix mistakes in the input data and fit parameters for a SVD fit.

    The Fixer is a Chain-of-Responsibility pattern.
    
    Attributes
    ----------
    molecules: list of :obj:`mollib.Molecule`
        A list of molecule objects.
    enabled: bool
        If True, then the Fixer will be executed.
    order: int
        The order in which to execute the fixer.
    """

    enabled = False
    order = 0

    def __init__(self, molecules):

        # Initialize the molecules attribute
        if isinstance(molecules, list):
            self.molecules = molecules
        elif isinstance(molecules, Molecule):
            self.molecules = [molecules, ]

        # Init subclasses
        self._subclass_instances = []
        for sub_cls in sorted(self.__class__.__subclasses__(),
                              key=lambda x: x.order):
            new_instance = sub_cls(molecules)
            self._subclass_instances.append(new_instance)

            # Determine whether this fixer is enabled or not from the settings
            subclass_name = sub_cls.__name__.lower()

            if hasattr(settings, 'enable_' + subclass_name):
                new_instance.enabled = getattr(settings,
                                               'enable_' + subclass_name)

    def fit(self, data):
        """Fit the RDC and RACS data with and SVD.
        
        Parameters
        ----------
        data: dict
            A dict with the data (to be fixed). 
            - **key**: interaction labels (str). ex: '14N-H'
            - **value**: RDC or RACS datum objects (:obj:`RDC` :obj:`RACS`)
        
        Returns
        -------
        RMS, Q, data_pred: float, float, dict
            - **RMS**: The root-mean-square deviation of the fit.
            - **Q**: The fit Q-factor (in percent)
            - **data_pred**: The predicted data dictionary.
        """
        # Prepare the magnetic interactions for the molecules
        labels = data.keys()
        process = Process(self.molecules)
        magnetic_interactions = process.process(labels=labels)

        # Conduct the SVD on the data
        (data_pred, Saupe_components,
         stats) = calc_pa_SVD(magnetic_interactions, data)

        return stats['Overall']['RMS'], stats['Overall']['Q (%)'], data_pred

    def fix(self, data):
        """Fix the data to improve the SVD fit.

        Parameters
        ----------
        data: dict
        A dict with the data (to be fixed). 
        - **key**: interaction labels (str). ex: '14N-H'
        - **value**: RDC or RACS datum objects (:obj:`RDC` :obj:`RACS`)

        Returns
        -------
        data_fixed: dict or None
            A dict with the data. 
            
            - **key**: interaction labels (str). ex: '14N-H'
            - **value**: RDC or RACS datum objects (:obj:`RDC` :obj:`RACS`)
            
            None is returned if none of the fixes worked
        fixes: list or str
            A list of strings of the fixes conducted to generate data_fixed.
        """
        data_fixed = None
        data_returned = None
        fixes = []
        # Process all of the subclasses and store their results
        for instance in self._subclass_instances:
            if not instance.enabled:
                continue
            data_fixed = data_fixed if data_fixed is not None else data
            data_returned, f = instance.fix(data_fixed)
            data_fixed = (data_returned if data_returned is not None
                          else data_fixed)
            fixes += f

        return data_fixed, fixes

    def copy_data(self, data):
        """Get a deep copy of the data, which can be modified without
        consequence."""
        return deepcopy(data)


class SignFixer(Fixer):
    """Fix the sign of RDCs and RACS"""

    order = 10

    def fix(self, data):
        # Prepare the fixed message
        msg = ("Inverting the sign of '{}' interactions improved the overall "
               "Q-factor from {:.1f}% to {:.1f}%.")

        # Get the reference RMS
        RMS_ref, Q_ref, data_pred = self.fit(data)

        # Setup the fixed data and fixes return values
        data_fixed = None
        fixes = []

        # Get the difference interaction types for the data
        interaction_types = {interaction_type(i) for i in data.keys()}

        # Process the N-H RDC sign first. This one is likely to be wrong
        # if the user calculated the |J+D| - |J| instead of using the J-
        # coupling
        if 'N-H' in interaction_types or 'H-N' in interaction_types:
            # Remove the 'N-H' or 'H-N' interaction type from the
            # interactions_types
            interaction_types -= {'N-H', 'H-N'}

            # Try inverting the sign of N-H couplings in a new dataset
            data_copy = self.copy_data(data)

            for k,v in data_copy.items():
                k_type = interaction_type(k)  # get the interaction type of k
                if k_type == 'N-H' or k_type == 'H-N':
                    v.value *= -1.0

            # Calculate the updated fit
            new_RMS, new_Q, data_pred = self.fit(data_copy)

            # See if it's an improvement. If it is, keep it.
            if new_Q < Q_ref:
                fixes.append(msg.format('N-H', Q_ref, new_Q))
                Q_ref = new_Q
                RMS_ref = new_RMS
                data_fixed = data_copy

        # Process the other interaction types. These must be processed
        # sequentially
        for int_type in interaction_types:
            # Copy the dataset and try inverting the signs
            data_copy = (self.copy_data(data_fixed) if data_fixed is not None
                         else self.copy_data(data))

            # Invert all of the values for the given interaction type
            for k,v in data_copy.items():
                k_type = interaction_type(k)  # get the interaction type of k
                if k_type == int_type:
                    v.value *= -1.0

            # Calculate the updated fit
            new_RMS, new_Q, data_pred = self.fit(data_copy)

            # See if it's an improvement. If it is, keep it.
            if new_Q < Q_ref:
                fixes.append(msg.format(int_type, Q_ref, new_Q))
                Q_ref = new_Q
                RMS_ref = new_RMS
                data_fixed = data_copy

        return data_fixed, fixes


class NHScaleFixer(Fixer):
    """Fix the RDCs of couplings that have been scaled to match the magnitude
    of NH couplings. 
    
    All of the RDCs for a given interaction type are scaled together.
    """

    order = 30

    def fix(self, data):
        # Prepare the fixed message
        msg = ("Re-scaling the '{}' couplings from N-H values improved the "
               "overall Q-factor from {:.1f}% to {:.1f}%.")

        # Get the reference RMS
        RMS_ref, Q_ref, data_pred = self.fit(data)

        # Setup the fixed data and fixes return values
        data_fixed = None
        fixes = []

        # Get the difference interaction types for the data. We will not
        # rescale the 'N-H' interactions, so these can be removed.
        interaction_types = {interaction_type(i) for i in data.keys()}
        interaction_types -= {'N-H', 'H-N'}

        # Calculate the A-matrix for all of the interactions. This is done to
        # get the default scaling constant for each interaction type.
        labels = data.keys()
        process = Process(self.molecules)
        magnetic_interactions = process.process(labels=labels)
        amatrix = magnetic_interactions[0]  # only need the first molecule

        # Process the interaction types. These must be processed
        # sequentially
        for int_type in interaction_types:
            # Copy the dataset and try scaling the signs altogether
            data_copy = (self.copy_data(data_fixed) if data_fixed is not None
                         else self.copy_data(data))

            # Scale all of the values for the given interaction type
            for k, v in data_copy.items():
                k_type = interaction_type(k)  # get the interaction type of k
                if k_type == int_type and k in amatrix:
                    scale, _ = amatrix[k]

                    # The scale includes the factor of 2 needed for RDCs,
                    # whereas the default_predicted_rdcs['N-H'] value does not.
                    # For this reason, the denominator is multiplied by 2.
                    # The absolute value is taken here because this should not
                    # change the sign of the RDC/RACS
                    # TODO: move factor of 2 in the scale to the SVD calculation
                    cnst = abs(scale /
                               (settings.default_predicted_rdcs['N-H'] * 2.))
                    v.value *= cnst

            # Calculate the updated fit
            new_RMS, new_Q, data_pred = self.fit(data_copy)

            # See if it's an improvement. If it is, keep it.
            if new_Q < Q_ref:
                fixes.append(msg.format(int_type, Q_ref, new_Q))
                Q_ref = new_Q
                RMS_ref = new_RMS
                data_fixed = data_copy

        return data_fixed, fixes


# Not Implemented
# class StereoFixer(Fixer):
#     order = 20


class OutlierFixer(Fixer):
    """Removes outliers from the data."""

    order = 100

    def fix(self, data):
        # Prepare the fixed message
        msg = ("Removing outlier data points {} improved the overall Q-factor "
               "from {:.1f}% to {:.1f}%.")

        # Get the reference RMS
        RMS_ref, Q_ref, data_pred = self.fit(data)

        # Setup the fixed data and fixes return values
        data_fixed = None
        fixes = []

        warning, bad = find_outliers(data, data_pred)

        # See if there are any outliers
        if len(warning) > 0 or len(bad) > 0:
            # Copy the data and remove the outliers
            data_copy = self.copy_data(data)
            data_copy = {k: v for k,v in data_copy.items()
                         if k not in warning and k not in bad}

            # Recalculate the fit
            new_RMS, new_Q, new_data_pred = self.fit(data_copy)

            # See if it's an improvement
            if new_Q < Q_ref:
                outliers = ", ".join(bad + warning)
                fixes.append(msg.format(outliers, Q_ref, new_Q))
                data_fixed = data_copy

        return data_fixed, fixes


# Not Implemented: This will not be implement until more CSA datasets are
# produced and available
# class CSAOptimizer(Fixer):
#     """Optimize the CSA tensor parameters."""


class SplitFixer(Fixer):
    """Splits the dataset into contiguous pieces and fits them individually.
    """
    pass