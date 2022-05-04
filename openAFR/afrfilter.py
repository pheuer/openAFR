# -*- coding: utf-8 -*-
"""
@author: Peter Heuer
"""

import numpy as np

class AFRFilter:
    """
    This class describes the geometry of a given AFR filter, based on a name
    e.g. "AFR3".
    
    
    name -> Name of filter, used to select pre-defined filters
    
    filter_radii -> List of edge transition radii
    
    filter_binary -> Binary mask values for each region
    
    filter_mask -> Function of mask over the radii array
    
    radii -> Radii at which to calculate the mask. By defualt, spans the filter_radii
    
    """
    
    
    def __init__(self, name, filter_radii=None, 
                 filter_binary=None, filter_mask=None, radii=None):
    
        # Go through the pre-defined names first
        if name == 'AFR3':
            filter_radii  = np.array( [0.00, 0.250, 0.950, 1.776, 2.729, 3.807, 5.012, 6.343, 7.800, 9.384,
                     11.093, 12.929, 14.891, 16.979, 19.194, 21.534, 24.000]) #AF3 filter edge radii[mm]
            filter_binary = np.array([0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1])
            filter_mask = None

        else:
            if filter_radii is None:
                raise ValueError(f"Unknown filter name: {name}")
                
                
        if radii is None:
            self.radii = np.linspace(0, np.max(filter_radii), num=200)
        else:
            self.radii = radii
        
        # If a filter mask array is not provided, calculate it from the binary array
        if filter_mask is None:
            if filter_binary is None:
                raise ValueError("Must specify either filter_mask or filter_binary, and not both.")
            
            filter_mask = np.zeros(self.radii.size)
            for i in range(len(filter_radii)-1):
                a = np.argmin(np.abs(self.radii - filter_radii[i]))
                b = np.argmin(np.abs(self.radii - filter_radii[i+1]))
                filter_mask[a:b] = filter_binary[i]
        
        # If a filter mask is provided, create a binary mask using 0.5 as the threshold 
        else:
            if filter_binary is not None:
                raise ValueError("Must specify either filter_mask or filter_binary, and not both.")
                
            filter_binary = np.zeros(filter_radii.size - 1)
            for i in range(len(filter_radii)-1):
                a = np.argmin(np.abs(self.radii - filter_radii[i]))
                b = np.argmin(np.abs(self.radii - filter_radii[i+1]))
                if np.mean(filter_mask[a:b]) > 0.5:
                    filter_binary[i] = 1
                else:
                    filter_binary[i] = 0
        
    
        # Make an array of the centers of each region
        self.filter_centers = np.zeros(filter_radii.size - 1)
        for i in range(len(filter_radii)-1):
            self.filter_centers[i] = 0.5*(filter_radii[i] + filter_radii[i+1])
            
            
        self.filter_regions = np.arange(filter_binary.size)
        self.filter_radii = filter_radii
        self.filter_binary = filter_binary
        self.filter_mask = filter_mask