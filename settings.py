"""
constants.py
A fast and accurate semi-analytic method for predicting small molecule adsorption in nanoporous materials ideally suited for high-throughput screening applications. Courtesy of R.Q.Snurr Research Group, Northwestern Univerisity.

Handles the primary functions
"""


def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format)

    Replace this function and doc string for your own project

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
# import numpy as np
# import scipy.integrate as intg 

# # constants
# R=8.314 # the ideal gas constant
# theta=0 # The Temkin approximatoin ot BET constant
# epsilon = 36.7 # lJ well depth for the gust-guest correction


# # The adsorption/desorption conditions
# Tads=77
# Tdes=160

# # Linearly spaced pressure prediction points
# press_power = np.linspace(low_lim, high_lim, num=number)
# Pressures   = np.power(10, press_power)
# Pressures   = np.around(Pressures, decimals=8)