# -*- coding: utf-8 -*-
"""
@author: Peter Heuer
"""

import numpy as np
import h5py

import astropy.units as u

import matplotlib.pyplot as plt
import scipy.ndimage

from scipy.interpolate import interp1d
from scipy.signal import find_peaks

import matplotlib.patches as patches

from openAFR.afrfilter import AFRFilter


class AFRLineout:
    
    def __init__(self,*args):
        
        # Interpret one argument as the AFR filter
        if len(args) >= 1:
            self.set_filter(args[0])
        
        # Interpret as second arg as an AFR file path
        if len(args) >= 2:
            self.load_afr(args[1])
        
        # Interpret as second arg as an AFR file path
        if len(args) >= 3:
            self.load_shadowgraph(args[2])
        
        
        
        # These variables are dictionaries which hold various data sources
        # under the names 
        # 'afr'
        # 'shadowgraph'
        
        self.data = {}
        self.data_hax  = {}
        self.data_vax  = {}
        self.line  = {}
        self.line_ax  = {}
        # Box that outlines the region of the lineout in coordinates
        self.lineout_box = {} 
        self.binary_signal = {}
        # The locations (in data space) of each transition edge, + the first
        # and last element
        self.edges = {}
        # AFRFilter objets for each dataset
        self.AFRFilters = {}
        self.angle = {}
        # The filters corresponding to each segment
        # where nsegments = self.edges.size - 2
        self.segment_filters = {}
        # An array of the filter chosen at each location, for plotting
        self.filters_arr = {}
        # Points at which the angle is measured
        self.point_loc = {}
        self.point_val = {}
        # Solution
        self.ne = {}

        
    def set_filter(self, name,  filter_obj):
        self.AFRFilters[name] = filter_obj
        
    def load_data(self, path, name,
                  dhax=5*u.um, dvax=5*u.um):
        """
        path -> Path to the file
        
        name -> One of 
            - 'afr'
            - 'shadowgraph'
            
        dhax, dvax -> pizel size in horizontal and vertical directions
            
        """
        with h5py.File(path, 'r') as f:
            self.data[name] = f['Streak_array'][0,:,:].T

        nx, ny = self.data[name].shape
        self.data_hax[name] = (np.arange(nx) - int(nx/2))*dhax
        self.data_vax[name] = (np.arange(ny) - int(ny/2))*dvax
        
        
    def rotate_data(self, name, rotation):
        
        self.data[name] = scipy.ndimage.rotate(self.data[name], rotation, cval=np.nan)
        
        # Make new axes, assuming the pixel size hasn't changed
        nx, ny = self.data[name].shape
        dhax = np.mean(np.gradient(self.data_hax[name]))
        dvax = np.mean(np.gradient(self.data_vax[name]))
        self.data_hax[name] = (np.arange(nx) - int(nx/2))*dhax
        self.data_vax[name] = (np.arange(ny) - int(ny/2))*dvax
        

    def filter2d(self, name, *, cutoff=None):
        
        """
        name -> name of the dataset to filter
        
        cutoff -> critical frequencies for top and bottom of bandpass
        In 1/m
         
        """
        data = self.data[name]
        dx = np.mean(np.gradient(self.data_hax[name])).to(u.m).value
        dy = np.mean(np.gradient(self.data_vax[name])).to(u.m).value
        
    
        xfreq = np.fft.fftfreq(data.shape[0], d=dx)
        xfreq = np.fft.fftshift(xfreq)
        yfreq = np.fft.fftfreq(data.shape[1], d=dy)
        yfreq = np.fft.fftshift(yfreq)

        xfreq, yfreq = np.meshgrid(xfreq, yfreq, indexing='ij')
        rfreq = np.sqrt(xfreq**2 + yfreq**2)
        nx, ny = rfreq.shape
        dfreq = np.mean(np.gradient(xfreq))

        # Create the mask
        mask = np.zeros(rfreq.shape)

        i = (rfreq>=cutoff[0]) & (rfreq<=cutoff[1])
        dC = cutoff[1] - cutoff[0]
        X = (rfreq[i]-cutoff[0])/dC
        mask[i] = np.sin(np.pi*X)
        
        
        fig, ax = plt.subplots()
        ax.set_xscale('log')
        ax.plot(xfreq, mask[:, 1000])
        plt.show()
    
    
        
        fft = np.fft.fftshift(np.fft.fft2(data))
        fft *= mask
        fft = np.fft.fftshift(fft)
        self.data[name] = np.real(np.fft.ifft2(fft))
        
        
    def filter2d_flattop(self, name, *, cutoff=None, transition=2e3):
        
            """
            name -> name of the dataset to filter
            
            cutoff -> critical frequencies for top and bottom of bandpass
            In 1/m
             
            """
            cutoff = (2e3, 5e3)
            data = self.data[name]
            dx = np.mean(np.gradient(self.data_hax[name])).to(u.m).value
            dy = np.mean(np.gradient(self.data_vax[name])).to(u.m).value
            
            
            
            xfreq = np.fft.fftfreq(data.shape[0], d=dx)
            xfreq = np.fft.fftshift(xfreq)
            yfreq = np.fft.fftfreq(data.shape[1], d=dy)
            yfreq = np.fft.fftshift(yfreq)
    
            xfreq, yfreq = np.meshgrid(xfreq, yfreq, indexing='ij')
            rfreq = np.sqrt(xfreq**2 + yfreq**2)
            nx, ny = rfreq.shape
            dfreq = np.mean(np.gradient(xfreq))
    
    
            # Create the mask
            mask = np.zeros(rfreq.shape)
            
            c0 = cutoff[0]-transition
            c1 = cutoff[1]+transition
            
            i = (rfreq>=cutoff[0]) & (rfreq<=cutoff[1])
            mask[i] = 1
            
            i = (rfreq>c0) & (rfreq<cutoff[0])
            X = (rfreq[i]-c0)/transition
            mask[i] = np.sin(0.5*np.pi*X)
            
            i = (rfreq<c1) & (rfreq>cutoff[1])
            X = (rfreq[i]-c1)/transition
            mask[i] = np.sin(0.5*np.pi*X + np.pi)
            
            fig, ax = plt.subplots()
            ax.set_xscale('log')
            ax.plot(xfreq, mask[:, 1000])
            plt.show()
            
            
            
            fft = np.fft.fftshift(np.fft.fft2(data))
            
            
            
            plt.pcolormesh(xfreq, yfreq, np.real(fft))
            plt.show()
            
            fft *= mask
            
            plt.pcolormesh(xfreq, yfreq, np.real(fft))
            plt.show()
            
            fft = np.fft.fftshift(fft)
            data = np.real(np.fft.ifft2(fft))
            
            plt.pcolormesh(data)
            plt.show()
            
            self.data[name] = data
             
        

    def plot_data(self, names=None):
        
        if names is None:
            names =list(self.data.keys())

        fig, axarr = plt.subplots(ncols=len(names))
        
        if not isinstance(axarr, np.ndarray):
            axarr = np.array([axarr,])
        
        for i, name in enumerate(names):
            
            if name == 'afr':
                title = "AFR (AFR1)"
            elif name == 'shadowgraph':
                title = "Shadowgraph (AFR2)"
            else:
                title = name
            
            ax = axarr[i]
            ax.set_title(title)
            ax.set_aspect('equal')
            ax.pcolormesh(self.data_hax[name].to(u.um).value,
                          self.data_vax[name].to(u.um).value,
                          self.data[name].T)
            
            if name in list(self.lineout_box.keys()):
                rect = patches.Rectangle(*self.lineout_box[name], linewidth=1, edgecolor='red',
                                         facecolor='none')
                ax.add_patch(rect)

        plt.show()
        
        
    def lineout(self, name, start, size, reverse=False, scale=True):
        """
        This function pulls a lineout from a dataset
        
        start -> lower left corner (x, y), np.ndarray
        size -> Size of box (width, height), np.ndarray

        """

        
        x0 = np.argmin(np.abs(start[0] - self.data_hax[name]))
        x1 = np.argmin(np.abs(start[0] + size[0] - self.data_hax[name]))
        y0 = np.argmin(np.abs(start[1] - self.data_vax[name]))
        y1 = np.argmin(np.abs(start[1] + size[1] - self.data_vax[name]))
        

        line = np.mean(self.data[name][x0:x1, y0:y1], axis=-1)
        
        if scale:
            line += - line[-1]
            line *= 1/np.max(line)

        line_ax = self.data_hax[name][x0:x1]
        
        if reverse:
            line = line[::-1]
            
        self.line[name] = line
        self.line_ax[name] = line_ax
        self.lineout_box[name] = (start.to(u.um).value, 
                                  size[0].to(u.um).value, 
                                  size[1].to(u.um).value)


    def calc_ne_shadowgraph(self, name='shadowgraph',
                            zero_region=slice(0,200,None)):
        Lz = 0.3 # cm
        Loi = 0.3 # cm
        ncrit = 1.6e22 # n_c for 4w 263nm light
        dx = 5*1e-4 # 5 um -> cm
        
        line = self.line[name] - np.mean(self.line[name][zero_region])
        
        dndx = 2*ncrit/(Lz*Loi)*np.cumsum(line)*dx
        dndx += -np.mean(dndx[zero_region]) # Assume zero gradient at the beginning

        ne = np.cumsum(dndx)*dx
        ne += -np.mean(ne[zero_region])
        
        plt.plot(self.line_ax[name], line/np.max(np.abs(line)))
        plt.plot(self.line_ax[name],dndx/np.max(np.abs(dndx)))
        plt.plot(self.line_ax[name], ne/np.max(np.abs(ne)))
        plt.show()
        
        self.shadowgraph_ne = ne
        
        
    def plot_line(self, names=None, plot_binary=True,
                  plot_filters=True, plot_segments=True,
                  plot_points=True, oneplot=False):
        """
        
        

        Parameters
        ----------
        names : TYPE, optional
            Names of the datasets to plot.
        plot_binary : TYPE, optional
            Plot the binary thresholded signal?
        plot_filters : TYPE, optional
            Plot the trace of filters at each position?
        plot_segments : TYPE, optional
            Plot the segments axis?
        plot_points : TYPE, optional
            Plot the points? 
        oneplot : bool, optional
            If true, plot all datasets on one plot. If false, make one plot
            per dataset. Default False.

        Returns
        -------
        None.

        """
        
        if names is None:
            names = list(self.line.keys())
            
        if oneplot:
            fig, ax = plt.subplots(figsize=(6,3))
            ax.set_xlabel("um")
            ax.set_ylim(0,1.3)
            
        for name in names:
            
            
            if not oneplot:
                fig, ax = plt.subplots(figsize=(6,4))
                fig.subplots_adjust(top=0.8)
                fig.suptitle(name)
                ax.set_xlabel("um")
                ax.set_ylim(0,1.3)
            
            ax.plot(self.line_ax[name], self.line[name], 
                    color='black', linestyle='dashed', label=name)
            
            if name in list(self.binary_signal.keys()) and plot_binary:
                ax.plot(self.line_ax[name], self.binary_signal[name],
                        color='lime', label=f'Binary {name}')
                
                
            if name in list(self.edges.keys()) and plot_segments:
                # Calculate the center of each segment
                segment_centers = []
                for i in range(len(self.edges[name])-1):
                    segment_centers.append((self.edges[name][i] + 
                                            self.edges[name][i+1]).to(u.um).value/2)
                
                segment_centers = np.array(segment_centers)
                segments = np.arange(segment_centers.size)
                
                ax2 = ax.twiny()
                ax2.set_xlim(ax.get_xlim())
                ax2.set_xticks(segment_centers)
                ax2.set_xticklabels(segments)
                ax2.set_xlabel("Segment")
                
                
            if name in list(self.filters_arr.keys()) and plot_filters:
                ax2 = ax.twinx()
                ax2.plot(self.line_ax[name], self.filters_arr[name], color='red')
                ax2.set_ylabel("Filter Region")
                ax2.set_ylim(0, 1.3*np.max(self.filters_arr[name]))
                
                    
            if name in list(self.point_loc.keys()) and plot_points:
                ax2.scatter(self.point_loc[name], self.point_val[name],
                            color='blue')
            
        """
        if self.shadowgraph_ne is not None:
            #ax.plot(self.shadowgraph_ax, self.shadowgraph_line, color='green', linestyle='dashed')
            ne_norm = self.shadowgraph_ne/np.max(self.shadowgraph_ne)
            ax.plot(self.line_ax['shadowgraph'], ne_norm, color='green', linestyle='dashed')
        """
        
        plt.show()
        
        
    
    def _find_edges(self, name, threshold=None):
        if threshold is None:
            threshold = input("Input threshold value (float 0-1.0):")
            threshold = float(threshold)
        
        self.binary_signal[name] = np.where(self.line[name] > threshold, 1, 0)
        
        # Find edges by looking at peaks in the absolute value of the 
        # derivative
        grad_binary_signal = np.gradient(self.binary_signal[name])
        peaks, _ = find_peaks(np.abs(grad_binary_signal))
        
        # Create a list of the edges, adding on the beginning and end of the
        # array
        self.edges[name] = list(self.line_ax[name][peaks])
        # Insert/append the beginning and ends of the array
        self.edges[name].insert(0, self.line_ax[name][0])
        self.edges[name].append(self.line_ax[name][-1])
        
        
    def _monotonic_filters(self, name):
        self.segment_filters[name] = np.arange(len(self.edges[name])-1)
        
        
        
    def _find_filters(self, name, filters=None):
        
        # Get user input (or use keyword value) for the filter band
        # for each region
        
        # Validate that the number of filter regions provided matches the number of
        # segments found
        if filters is not None:
            if len(filters) != len(self.edges[name])-1:
                raise ValueError("Wrong number of filter for segments")
        
        
        self.segment_filters[name]= np.zeros(len(self.edges[name])-1).astype(np.int32)
        for i in range(len(self.segment_filters[name])):
            
            # Find the binary value of the signal in the center of this region
            a  = np.argmin(np.abs(self.line_ax[name] - 
                                  0.5*(self.edges[name][i]+self.edges[name][i+1])))
            binary_value = self.binary_signal[name][a]
            
            if filters is None:
                # Get user input
                while True:
                    x = input(f"Specify filter # for segment {i}: ")
                    x = int(x)
                    if np.isclose(binary_value, self.AFRFilters[name].filter_binary[x]):
                        self.segment_filters[name][i] = x
                        break
                    else:
                        print(f"Bad filter match: segment {i} is binary "
                              f"{binary_value} but AFR filter region "
                              f"{x} is binary value "
                              f"{self.AFRFilters[name].filter_binary[x]}")
            else:
                x = int(filters[i])
                if np.isclose(binary_value, self.AFRFilters[name].filter_binary[x]):
                    self.segment_filters[name][i] = x
                else:
                    raise ValueError(f"Bad filter match: segment {i} is binary "
                          f"{binary_value} but AFR filter region "
                          f"{x} is binary value "
                          f"{self.AFRFilters[name].filter_binary[x]}")


    def _assign_points(self, name):
        """
        With edges found and filters chosen for each segment, chose points on
        either side of each transition

        """
        
        self.point_loc[name] = []
        self.point_val[name] = []
        for i, edge in enumerate(self.edges[name]):
            # If this is the first or last edge, just include one point
            # since this isn't actually a transition
            if i==0:
                self.point_loc[name].append(edge.to(u.um).value)
                self.point_val[name].append(self.segment_filters[name][i])
            elif i==len(self.edges[name])-1:
                self.point_loc[name].append(edge.to(u.um).value)
                self.point_val[name].append(self.segment_filters[name][i-1])
               
            # Otherwise, add points both before and after
            else:
                # TODO: Chose displacement for each transition based on the
                # lineshape: maybe 90% and 10% points of that section?
                
                # Displacement is the amount to offset each point from the edge
                displacement = 5 #um
                # Add the lower point
                self.point_loc[name].append(edge.to(u.um).value - displacement)
                self.point_val[name].append(self.segment_filters[name][i-1])
                
                # Add the higher point
                self.point_loc[name].append(edge.to(u.um).value + displacement)
                self.point_val[name].append(self.segment_filters[name][i])
        
        # Cast as np.array
        self.point_loc[name] = np.array(self.point_loc[name])
        self.point_val[name] = np.array(self.point_val[name])

    
    def process(self, name, threshold=None, filters=None):
        """
        This function thresholds the lineout, determines (with user input)
        which filter ring each region belongs to, then figures out the
        radius corresponding to each transition.
        
        
        threshold : float, < 1
            Thresholding value for determing the binary AFR signal
            
        
        filters : list of ints
            A list of AFR filter regions for each segment identified in the
            AFR signal. Must exactly match the number of segments found. Optional. 
            

        
        """
        self._find_edges(name, threshold=threshold)
        self.plot_line()

        if filters == 'monotonic':
            self._monotonic_filters(name)
        else:
            self._find_filters(name, filters=filters)
        
        # Make an array that shows which filter ring corresponds to each point
        self.filters_arr[name] = np.zeros(self.line_ax[name].size)
        for i in range(len(self.edges[name])-1):
            a = np.argmin(np.abs(self.edges[name][i]  - self.line_ax[name]))
            b = np.argmin(np.abs(self.edges[name][i+1]  - self.line_ax[name]))
            self.filters_arr[name][a:b] = self.segment_filters[name][i]
        self.filters_arr[name][a:] = self.segment_filters[name][-1] # Force end value to last filter
        
        self.plot_line()
        
        self._assign_points(name)
        self.plot_line()
        

    def plot_solution(self, names=None):
        
        if names is None:
            names = list(self.ne.keys())
        
        fig, ax = plt.subplots(figsize=(6,3))
        
        for name in names:
            
            if name in list(self.angle.keys()):
                ax.plot(self.line_ax[name], self.angle[name], 
                        color='C0', label=f'Angle {name}')
                ax.set_ylabel("Angle (deg)")
            
            ax2 = ax.twinx()
            ax2.plot(self.line_ax[name], self.ne[name], 
                     color='C1', label=f'$n_e$ {name}')
            ax2.set_ylabel("$n_e$ (cm$^{-3}$")
            ax2.set_ylim(0, 2e20)
        
            ax.plot([], [], color='C1', label=f'$n_e$ {name}')
        ax.legend()
        
        
    
    
    def solve(self, name,  Lz=0.3 * u.cm):
    
        transition_radii = self.AFRFilters[name].filter_radii[self.point_val[name].astype(np.int32)]

        
        # Constant is for 4w setup geometry: 0.368 deg/mm
        transition_angle  = 0.368 * transition_radii
        
        
        # linearly interpolate from these values over the whole axis
        interp = interp1d(self.point_loc[name], transition_angle, kind='linear', 
                          bounds_error=False, fill_value = (0, transition_angle[-1]))
        self.angle[name] = interp(self.line_ax[name])
        

        # Integrate to get the density
        ncrit = 1.6e22 # n_c for 4w 263nm light
        Lz = Lz.to(u.cm).value# cm
        dx = 5*1e-4 # 5 um -> cm
        
        self.ne[name] = np.cumsum(np.deg2rad(self.angle[name]))*dx*2*ncrit/Lz
        
        self.plot_solution()
  