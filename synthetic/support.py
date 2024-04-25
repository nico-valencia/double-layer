from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import copy
import matplotlib.gridspec as gridspec
import os.path
try:
    import matplotlib
    matplotlib.use('TkAgg')
except:
    pass
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
# import splitwavepy 
from splitwavepy.core.data import Data, WindowPicker

from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")

#################### particle motion for plot_trace ####################

def trace_ppm(self, ax, **kwargs):
    """Plot particle motion on *ax* matplotlib axis object, specifically for plot_trace code.
    """
    
    data = self.copy()
    # data.rotateto(0)
    x, y = data.x , data.y
    t = data.t()

    # middle third of uncut trace 
    x = x[int(len(x)*0.33):int(len(x)*0.66)]  
    y = y[int(len(y)*0.33):int(len(y)*0.66)]
    t = t[int(len(t)*0.33):int(len(t)*0.66)]

            
    # plot data
    # ax.plot(self.chop().y,self.chop().x)
    
    # multi-colored
    norm = plt.Normalize(t.min(), t.max())
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='plasma', norm=norm, alpha=0.7)
    lc.set_array(t)
    lc.set_linewidth(2)
    line = ax.add_collection(lc)
    # plt.colorbar(line)

    # set limit
    lim = np.abs(self.data()).max() * 1.1
    if 'lims' not in kwargs: kwargs['lims'] = [-lim, lim] 
    ax.set_aspect('equal')
    ax.set_xlim(kwargs['lims'])
    ax.set_ylim(kwargs['lims'])

    # set labels
    if 'cmplabels' not in kwargs: kwargs['cmplabels'] = data.cmplabels
    ax.set_xlabel('Radial')
    ax.set_ylabel('Transverse')
    
    # turn off tick annotation
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    return

############################################################



#################### particle motion for eigenvalue method ####################

def _pppm(self, baz,calcbaz, ax, **kwargs):
    """Plot particle motion on *ax* matplotlib axis object.
    """
    
    data = self.copy()
    # data.rotateto(0)
    x, y = data.x , data.y
    t = data.t()
            
    # plot data
    # ax.plot(self.chop().y,self.chop().x)
    
    # multi-colored
    norm = plt.Normalize(t.min(), t.max())
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='plasma', norm=norm, alpha=0.7)
    lc.set_array(t)
    lc.set_linewidth(2)
    line = ax.add_collection(lc)
    # plt.colorbar(line)

    # set limit
    lim = np.abs(self.data()).max() * 1.1
    if 'lims' not in kwargs: kwargs['lims'] = [-lim, lim] 
    ax.set_aspect('equal')
    ax.set_xlim(kwargs['lims'])
    ax.set_ylim(kwargs['lims'])

    # set labels
    if 'cmplabels' not in kwargs: kwargs['cmplabels'] = data.cmplabels
    ax.set_xlabel('Radial'+ '\n' + 
                  '(Actual baz: ' + str(round(float(baz),3))+')' + '\n' + 
                  '(Calculated baz: ' + str(round(calcbaz,3)) + ' or ' + str(round((calcbaz+180),3)) + ')'  )
    ax.set_ylabel('Transverse')
    
    # turn off tick annotation
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    return

############################################################



#################### particle motion for transverse minimization method ####################

def trans_pppm(self, baz,calcbaz, ax, **kwargs):
    """Plot particle motion on *ax* matplotlib axis object.
    """
    
    data = self.copy()
    data.rotateto(0)
    x, y = data.x , data.y
    t = data.t()
            
    # plot data
    # ax.plot(self.chop().y,self.chop().x)
    
    # multi-colored
    norm = plt.Normalize(t.min(), t.max())
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='plasma', norm=norm, alpha=0.7)
    lc.set_array(t)
    lc.set_linewidth(2)
    line = ax.add_collection(lc)
    # plt.colorbar(line)

    # set limit
    lim = np.abs(self.data()).max() * 1.1
    if 'lims' not in kwargs: kwargs['lims'] = [-lim, lim] 
    ax.set_aspect('equal')
    ax.set_xlim(kwargs['lims'])
    ax.set_ylim(kwargs['lims'])

    # set labels
    if 'cmplabels' not in kwargs: kwargs['cmplabels'] = data.cmplabels
    ax.set_xlabel('Radial'+ '\n' + 
                  '(Actual baz: ' + str(round(float(baz),3))+')' 
                  )
    ax.set_ylabel('Transverse')
    
    # turn off tick annotation
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    return

############################################################



#################### plots trace and particle motion only ####################

def plot_trace(data,baz, **kwargs):
    """
    Plot trace data and particle motion
    """

    plt.figure(figsize=(12, 3))     
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
    
    # trace
    ax0 = plt.subplot(gs[0])
    trace = copy.deepcopy(data)
    # trace.rotateto(baz)
    trace._ptr(ax0, **kwargs)
    
    # particle  motion
    ax1 = plt.subplot(gs[1])
    # data._ppm( ax1, **kwargs) 
    trace_ppm(data,ax1,**kwargs) 
    # trans_pppm(data, baz,0, ax1, **kwargs) 
    
                                
    # show
    plt.tight_layout()
    
    if not 'figname' in kwargs:
        kwargs['figname'] = "dataPlot.png"

    if not 'dpi' in kwargs:
        kwargs['dpi'] = 300

    plt.savefig(kwargs['figname'], bbox_inches='tight', dpi=kwargs['dpi'])

############################################################



#################### returns maximum lambda ratio ####################

def lambda_filter_2(self,**kwargs):             # used to filter only the most quality eigenvalue measurements (when including all data, not filtering during the measurement process)
                                                                    # uses all the data from eigenvalue method including errorbars and lambda ratio 
    if 'vals' not in kwargs: 
        kwargs['vals'] = self.lam1 / self.lam2
        
    return np.max(kwargs['vals']), self.dfast, self.dlag

############################################################



#################### returns maximum cross correlation (w/ preset constraints) ####################

def rotation_filter(self,minxc,maxfast,maxlag,**kwargs):             # used to filter only the most quality eigenvalue measurements 
                                                                    # uses all the data from eigenvalue method including errorbars and lambda ratio 
    if 'vals' not in kwargs: 
        kwargs['vals'] = self.xc
        
    if np.max(kwargs['vals']) >= minxc:
        if self.dlag <= maxlag:
            if self.dfast <= maxfast:
                        # return True
                        return np.max(kwargs['vals']), self.dfast, self.dlag
            else:
                        return False
        else:
                        return False 
    else:
                        return False
    
############################################################



#################### eigenvalue method splitting corrections and gridsearch ####################

def plot_eigen_measure(m,baz,title,**kwargs):
    # setup figure and subplots
    fig = plt.figure(figsize=(12,6)) 
    fig.suptitle("Eigenvalue Method for: " + title)
    gs = gridspec.GridSpec(2, 3,
                        width_ratios=[1,1,2]
                        )    
    ax0 = plt.subplot(gs[0,0])
    ax1 = plt.subplot(gs[0,1])
    ax2 = plt.subplot(gs[1,0])
    ax3 = plt.subplot(gs[1,1])
    ax4 = plt.subplot(gs[:,2])
    
    # data to plot
    d1 = m.data.chop()
    d1f = m.srcpoldata().chop()
    d2 = m.data_corr().chop()
    d2s = m.srcpoldata_corr().chop()

    # display predicted back azimuth
    calcbaz = m.srcpol()  
    if calcbaz < 0:
           calcbaz = calcbaz + 180


    print(type(d1))
    print(calcbaz)
    
    # flip polarity of slow wave in panel one if opposite to fast
    # d1f.y = d1f.y * np.sign(np.tan(m.srcpol()-m.fast))
    
    # get axis scaling
    lim = np.abs(d2s.data()).max() * 1.1
    ylim = [-lim,lim]

    # original
    d1f._ptr(ax0,ylim=ylim,cmplabels=('Radial','Transverse'),**kwargs)
    _pppm(d1f,baz,calcbaz,ax1,lims=ylim,**kwargs)
    # corrected
    d2s._ptr(ax2,ylim=ylim,cmplabels=('Radial','Transverse'),**kwargs)
    _pppm(d2s,baz,calcbaz,ax3,lims=ylim,**kwargs)

    # error surface
    if 'vals' not in kwargs:
        # kwargs['vals'] = (m.lam1 - m.lam2) / m.lam2
        # kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'
        kwargs['vals'] = m.lam1 / m.lam2
        kwargs['title'] = r'$\lambda_1 / \lambda_2$'
    
    # add marker and info box by default
    if 'marker' not in kwargs: kwargs['marker'] = True
    if 'info' not in kwargs: kwargs['info'] = True
    if 'conf95' not in kwargs: kwargs['conf95'] = True
    m._psurf(ax4,**kwargs)
    
    # title
    if m.name != 'Untitled':
        plt.suptitle(m.name)
    
    # neaten
    plt.tight_layout()
    # plt.show()
    if not 'figname' in kwargs:
        kwargs['figname'] = "eigenM.png"

    if not 'dpi' in kwargs:
        kwargs['dpi'] = 300
    # plt.savefig(.png')

    plt.savefig(kwargs['figname'], bbox_inches='tight', dpi=kwargs['dpi'])

############################################################



#################### transverse minimization method splitting corrections and grid search ####################

def plot_trans_measure(m,baz,dist,depth,title,**kwargs):
    # setup figure and subplots
    fig = plt.figure(figsize=(14,7)) 
    # fig.suptitle("Transversal Min. Method for: " + title + "")
    fig.suptitle(f"Transversal Min. Method for {title}, Distance of {round(dist/1000,0)} km, Depth of {depth} km.")
    gs = gridspec.GridSpec(2, 3,
                        width_ratios=[1,1,2]
                        )    
    ax0 = plt.subplot(gs[0,0])
    ax1 = plt.subplot(gs[0,1])
    ax2 = plt.subplot(gs[1,0])
    ax3 = plt.subplot(gs[1,1])
    ax4 = plt.subplot(gs[:,2])  
    


    # data to plot
    d1 = m.data.chop()
    d1f = m.srcpoldata().chop()
    d2 = m.data_corr().chop()
    d2s = m.srcpoldata_corr().chop()

    # display predicted back azimuth
    calcbaz = m.srcpol()  
    if calcbaz < 0:
            calcbaz = calcbaz + 180


    # print(type(d1))
    # print(calcbaz)

    # flip polarity of slow wave in panel one if opposite to fast
    # d1f.y = d1f.y * np.sign(np.tan(m.srcpol()-m.fast))

    # get axis scaling
    lim = np.abs(d2s.data()).max() * 1.1
    ylim = [-lim,lim]

    # original
    d1f._ptr(ax0,ylim=ylim,cmplabels=('Radial','Transverse'),**kwargs)
    trans_pppm(d1f,baz,calcbaz,ax1,lims=ylim,**kwargs)
    # corrected
    d2s._ptr(ax2,ylim=ylim,cmplabels=('Radial','Transverse'),**kwargs)
    trans_pppm(d2s,baz,calcbaz,ax3,lims=ylim,**kwargs)

    # error surface
    if 'vals' not in kwargs:
        # kwargs['vals'] = (m.lam1 - m.lam2) / m.lam2
        # kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'
        kwargs['vals'] = m.lam1 / m.lam2
        kwargs['title'] = r'$\lambda_1 / \lambda_2$'

    # add marker and info box by default
    if 'marker' not in kwargs: kwargs['marker'] = True
    if 'info' not in kwargs: kwargs['info'] = True
    if 'conf95' not in kwargs: kwargs['conf95'] = True
    m._psurf(ax4,**kwargs)

    # title
    if m.name != 'Untitled':
        plt.suptitle(m.name)

    # neaten
    plt.tight_layout()
    # plt.show()
    if not 'figname' in kwargs:
        kwargs['figname'] = "eigenM.png"

    if not 'dpi' in kwargs:
        kwargs['dpi'] = 300
    # plt.savefig(.png')


    plt.savefig(kwargs['figname'], bbox_inches='tight', dpi=kwargs['dpi'])

############################################################



#################### select smaller window for trace ####################

def plotdata(self,baz,calcbaz,**kwargs):
       
    """
    Plot trace data and particle motion
    """

    self.rotateto(baz)

    fig = plt.figure(figsize=(12, 3))     
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
    
    # trace
    ax0 = plt.subplot(gs[0])
    self._ptr( ax0,cmplabels=('Radial','Transverse'), **kwargs)
    
        # d1f._ptr(ax0,ylim=ylim,cmplabels=('Radial','Transverse'),**kwargs)
    # trans_pppm(d1f,baz,calcbaz,ax1,lims=ylim,**kwargs)

    # particle  motion
    ax1 = plt.subplot(gs[1])
    # self._ppm( ax1, **kwargs)  
    trans_pppm(self,baz,calcbaz,ax1,**kwargs) #lims = ylim???? (normalizing data)
    
    # optional pick window
    if 'pick' in kwargs and kwargs['pick'] == True:
        windowpicker = WindowPicker(self,fig,ax0)
        windowpicker.connect()
                                
    # show
    plt.tight_layout()
    plt.show()

    self.rotateto(0)
    plt.close()

############################################################



# ########## supplementary data plots ##########

# def extradata(m,baz,depth,angle,title,**kwargs):
#     # setup figure and subplots
#     fig = plt.figure(figsize=(16,8)) 
#     fig.suptitle("extra data for: " + title)
#     gs = gridspec.GridSpec(3, 3,
#                         width_ratios=[1,1,2]
#                         )    
#     ax0 = plt.subplot(gs[0,0])
#     ax1 = plt.subplot(gs[0,1])
#     ax2 = plt.subplot(gs[1,0])
#     ax3 = plt.subplot(gs[1,1])
#     ax4 = plt.subplot(gs[:,2]) 
#     ax5 = plt.subplot(gs[2,0:1]) 
#     ax6 = plt.subplot(gs[2,1:2])  

#     # data to plot
#     d1 = m.data.chop()
#     d1f = m.srcpoldata().chop()
#     d2 = m.data_corr().chop()
#     d2s = m.srcpoldata_corr().chop()

#     # display predicted back azimuth
#     calcbaz = m.srcpol()  
#     if calcbaz < 0:
#             calcbaz = calcbaz + 180


#     print(type(d1))
#     print(calcbaz)

#     # flip polarity of slow wave in panel one if opposite to fast
#     # d1f.y = d1f.y * np.sign(np.tan(m.srcpol()-m.fast))

#     # get axis scaling
#     lim = np.abs(d2s.data()).max() * 1.1
#     ylim = [-lim,lim]

#     # original
#     d1f._ptr(ax0,ylim=ylim,cmplabels=('Radial','Transverse'),**kwargs)
#     trans_pppm(d1f,baz,calcbaz,ax1,lims=ylim,**kwargs)
#     # corrected
#     d2s._ptr(ax2,ylim=ylim,cmplabels=('Radial','Transverse'),**kwargs)
#     trans_pppm(d2s,baz,calcbaz,ax3,lims=ylim,**kwargs)

#     # error surface
#     if 'vals' not in kwargs:
#         # kwargs['vals'] = (m.lam1 - m.lam2) / m.lam2
#         # kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'
#         kwargs['vals'] = m.lam1 / m.lam2
#         kwargs['title'] = r'$\lambda_1 / \lambda_2$'

#     # add marker and info box by default
#     if 'marker' not in kwargs: kwargs['marker'] = True
#     if 'info' not in kwargs: kwargs['info'] = True
#     if 'conf95' not in kwargs: kwargs['conf95'] = True
#     m._psurf(ax4,**kwargs)

#     # ray paths and arrival times
#     traveltimes = model.get_travel_times(source_depth_in_km=depth, distance_in_degree=angle) #,phase_list=["P","S","SKS", "SKKS","ScS"])
#     arrivals = model.get_ray_paths(source_depth_in_km=depth, distance_in_degree=angle) #,phase_list=["P","S","SKS", "SKKS","ScS"])
#     arrivals.plot_rays(plot_type='cartesian',legend=False,show=False,fig=fig,ax=ax5)
#     arrivals.plot_times(legend=False,show=False,fig=fig,ax=ax6)
    




    

#     # ax6
#     # plt.table(traveltimes[0],traveltimes[1],loc=ax6)
#     # plt.axis('off')
#     print(traveltimes)
#     print(traveltimes[0])
#     print(len(traveltimes))
#     print(type(traveltimes))


#     # title
#     if m.name != 'Untitled':
#         plt.suptitle(m.name)

#     # neaten
#     plt.tight_layout()
#     # plt.show()
#     if not 'figname' in kwargs:
#         kwargs['figname'] = "eigenM.png"

#     if not 'dpi' in kwargs:
#         kwargs['dpi'] = 300
#     # plt.savefig(.png')


#     plt.savefig(kwargs['figname'], bbox_inches='tight', dpi=kwargs['dpi'])

# ##############################




#################### old code ####################
    
# def plot_measure(m,**kwargs):
#     # setup figure and subplots
#     fig = plt.figure(figsize=(12,6)) 
#     gs = gridspec.GridSpec(2, 3,
#                         width_ratios=[1,1,2]
#                         )    
#     ax0 = plt.subplot(gs[0,0])
#     ax1 = plt.subplot(gs[0,1])
#     ax2 = plt.subplot(gs[1,0])
#     ax3 = plt.subplot(gs[1,1])
#     ax4 = plt.subplot(gs[:,2])
    
#     # data to plot
#     d1 = m.data.chop()
#     d1f = m.srcpoldata().chop()
#     d2 = m.data_corr().chop()
#     d2s = m.srcpoldata_corr().chop()
    
#     # flip polarity of slow wave in panel one if opposite to fast
#     # d1f.y = d1f.y * np.sign(np.tan(m.srcpol()-m.fast))
    
#     # get axis scaling
#     lim = np.abs(d2s.data()).max() * 1.1
#     ylim = [-lim,lim]

#     # original
#     d1f._ptr(ax0,ylim=ylim,**kwargs)
#     d1._ppm(ax1,lims=ylim,**kwargs)
#     # corrected
#     d2s._ptr(ax2,ylim=ylim,**kwargs)
#     d2._ppm(ax3,lims=ylim,**kwargs)

#     # error surface
#     if 'vals' not in kwargs:
#         # kwargs['vals'] = (m.lam1 - m.lam2) / m.lam2
#         # kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'
#         kwargs['vals'] = m.lam1 / m.lam2
#         kwargs['title'] = r'$\lambda_1 / \lambda_2$'
    
#     # add marker and info box by default
#     if 'marker' not in kwargs: kwargs['marker'] = True
#     if 'info' not in kwargs: kwargs['info'] = True
#     if 'conf95' not in kwargs: kwargs['conf95'] = True
#     m._psurf(ax4,**kwargs)
    
#     # title
#     if m.name != 'Untitled':
#         plt.suptitle(m.name)
    
#     # neaten
#     plt.tight_layout()
#     # plt.show()
#     if not 'figname' in kwargs:
#         kwargs['figname'] = "eigenM.png"

#     if not 'dpi' in kwargs:
#         kwargs['dpi'] = 300

#     # plt.savefig(kwargs['figname'], bbox_inches='tight', dpi=kwargs['dpi'])
#     plt.show()


# def lambda_filter(self,minlam,maxfast,maxlag,**kwargs):             # used to filter only the most quality eigenvalue measurements 
#                                                                     # uses all the data from eigenvalue method including errorbars and lambda ratio 
#     if 'vals' not in kwargs: 
#         kwargs['vals'] = self.lam1 / self.lam2
        
#     if np.max(kwargs['vals']) >= minlam:
#         if self.dlag <= maxlag:
#             if self.dfast <= maxfast:
#                         # return True
#                         return np.max(kwargs['vals']), self.dfast, self.dlag
#             else:
#                         return False
#         else:
#                         return False 
#     else:
#                         return False

############################################################