import os
import numpy as np
import astropy.units as u

from openAFR.afrfilter import AFRFilter
from openAFR.hdf import ensure_hdf5
from openAFR.afrlineout import AFRLineout

directory = os.path.join("C://", "Users", "pvheu", "Box", "pheuer",
                         "Research","Experiments and Projects",
                         "2022_4_27 FoilMHD-EP-22A","Data","s36755")
afr1_path = ensure_hdf5(os.path.join(directory, 'AFR1_CCD-s36755_afr1_4wp.hdf'))

afr3filter = AFRFilter("AFR3")

a = AFRLineout()
a.set_filter('afr1', afr3filter) # filter is AFR3

a.load_data(afr1_path, 'afr1')
a.rotate_data('afr1', 87) #Angle chosen to get the fringes to be level


a.plot_data()

# Define the starting position and size of the lineout box
start = np.array([-4500, -1000])*u.um
size = np.array([5200, 500])*u.um

a.lineout( 'afr1', start, size)

a.plot_data()
a.plot_line()

# Threshold is used to turn the AFR into a binary signal
# 'monotonic' assumes that each filter transition corresponds to an increase
# in density
a.process('afr1', threshold=0.1, filters='monotonic')
a.plot_line()

a.solve('afr1')
a.plot_solution()
