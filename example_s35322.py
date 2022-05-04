import os
import numpy as np
import astropy.units as u

from openAFR.afrfilter import AFRFilter
from openAFR.hdf import ensure_hdf5
from openAFR.afrlineout import AFRLineout

directory = os.path.join("C://", "Users", "pvheu", "Box", "pheuer",
                         "Research","Experiments and Projects",
                         "2021_8_19 MagShckHeat-EP","Data","s35322")
afr1_path = ensure_hdf5(os.path.join(directory, 'AFR1_CCD-s35322_afr1_4wp.hdf'))
shadowgraph_path = ensure_hdf5(os.path.join(directory, 'AFR2_CCD-s35322_afr2_4wp.hdf'))

afr3filter = AFRFilter("AFR3")

a = AFRLineout()
a.load_data(afr1_path, 'afr1')

a.set_filter('afr1', afr3filter)
a.rotate_data('afr1', -30)
a.lineout( 'afr1', np.array([-5000, -250])*u.um,
          np.array([6000, 500])*u.um, reverse=True)

a.plot_data()
a.plot_line(oneplot=False)

a.process('afr1', threshold=0.3, filters=[0,1,4,1,2,3,4,5,6])

a.solve('afr1')
a.plot_solution()