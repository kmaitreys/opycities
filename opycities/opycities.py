import numpy as np


class Model:
    """
    This class contains all the parameters that are used in the model.
    A model intialises the default values for all the parameters that will
    go in the calculation of opacities.

    Parameters
    ----------
    nlam : int
        Number of wavelengths

    Attributes
    ----------
    amin : float
        Minimum grain size in microns
    amax : float
        Maximum grain size in microns
    apow : float
        Power law index for grain size distribution
    amean : float
        Mean grain size in microns
    asig : float
        Standard deviation of grain size distribution in microns
    na : int
        Number of grain sizes
    blendonly : bool
        If True, only calculate the blend opacity
    split : bool
        If True, split the grain size distribution into two parts
    quiet : bool
        If True, don't print anything to the screen
    verbose : bool
        If True, print more information to the screen
    debug : bool
        If True, print even more information to the screen
    write_grd : bool
        If True, write the grain size distribution to a file
    mmfss : bool
        If True, use the Mie theory to calculate the opacities
    lam : array
        Wavelengths in microns
    """
    def __init__(self, nlam):
        self.amin = 0.05
        self.amax = 3000.0
        self.apow = 3.5
        self.amean = 0.0
        self.asig = 0.0
        self.na = 0
        self.blendonly = False
        self.split = False
        self.quiet = False
        self.verbose = False
        self.debug = False
        self.write_grd = False
        self.mmfss = False
        self.lam = np.zeros(nlam)




model = Model(14)


print(model.lam)