#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright 2012-2013 Thoeger Emil Rivera-Thorsen.
# Distributed under the terms of the GNU General Public License v. 3
# http://www.gnu.org/licenses/gpl-3.0.html

""" Standalone version of the profile builder class used in grism.
The module consists of a class and a convenience function to load it quicker.

Classes:
========

ProfileEditor
-------------
Args:
-----

wavlens:  numpy array
    wavelength (or velocity) array.

indata:    numpy array
    Fluxes.

inerrs:    numpy array
    Array of standard deviations.

linecen:    float
    Wavelength to consider center of the transition. Use 0 in velocity space.

fitrange:    List of tuples.
    Optional list of (min, max) tuples of wavelength/velocity ranges to
    include in fit. If none given, will fall back on linecen +/- 15.

Functions:
==========

load_profile
------------

Convenience function.
Args:
-----

data:    pandas dataframe
    Run this module as program to see the format of dataframe.

"""

# IMPORTS:
# NumPy and other functionality
import scipy as sp
# Gaussian and other functions
import scipy.stats as stats
# Collection of constants in SI units
import scipy.constants as con
# Advanced, labelled data structures
import pandas as pd
# Traits - for the model handling and GUI
from traits.api import HasTraits, Float, List, Dict, Bool, Array, DelegatesTo,\
        PrototypedFrom, Instance, Button, Str, Range, Enum, Int, Property, \
        on_trait_change, Any
from traitsui.api import View, Group, HGroup, VGroup, Item, Spring, EnumEditor
from traitsui.menu import OKButton,  CancelButton, RevertButton, LiveButtons,\
    ModalButtons, UndoButton, ApplyButton
# Enable & chaco - for plotting in the model editor GUI
from enable.component_editor import ComponentEditor
from chaco.example_support import COLOR_PALETTE
from chaco.chaco_plot_editor import ChacoPlotItem
from chaco.api import ArrayPlotData, PlotLabel, Plot, HPlotContainer, \
        GridContainer, ImageData, bone, gist_heat, gist_rainbow, DataRange1D, \
        ScatterInspectorOverlay, ColorBar, LinearMapper, Legend, \
        color_map_name_dict, KeySpec, create_line_plot, LinePlot, \
        add_default_axes, add_default_grids, OverlayPlotContainer
from chaco.tools.api import ScatterInspector, ZoomTool, PanTool, \
        BroadcasterTool, LegendTool, RangeSelection, RangeSelectionOverlay


class ProfileEditor(HasTraits):
    """ The line profile intitial guess editor class.
    This is the line profile editor module.
    It can be used to provide initial guesses to any fitting package according
    to the user's tastes.

    Usage: ProfileEditor(wave, data, errors, center)

    * wave: one-dimensional wavelength-like array, can be either
            wavelength or velocity units.
    * data: one-dimensional spectrum.
    * errors: one-dimensional noise spectrum.
    * center: Float. Central wavelength. 0 if in velocity space.
    """
    # TODO: Implement the model in a different way. Make a class for each, add
    # them as new instances, use DelegatesTo for the important parts, maybe
    # even not, maybe just use the 'object.MyInstance.attribute' notation and
    # only store the plotdata of the given model in this class...?
    # Probably needs a new way to store the parameters, though. Either the
    # Components Dict can take an extra "Kind" keyword in it, or restructure
    # the whole thing into a DataFrame object...? The latter will require a
    # major amount of work. On the other hand, it could mean much better
    # modularity.
    # TODO: Eeerrrmmm That there thing that I just though about and then forgot
    # about almost immediately again. Yikes.

    CompNum = Int(1)
    Components = Dict
    Locks = Dict
    CompoList = List()

    x = Array
    mod_x = Array
    Feedback = Str

    Sigma = Range(.1, 200.)
    Centr = Range(-100., 100., 0.)
    Heigh = Range(0., 20000.)

    # Define vars to regulate whether the above vars are locked
    #    down in the GUI:
    # TODO: Want to keep these options? Or is it just unnecessary clutter?
    # ANSWER: Keep! For when fitting.
    LockSigma = Bool()
    LockCentr = Bool()
    LockHeigh = Bool()
    LockConti = Bool()

    continuum_estimate = Range(0., 2000.)
    plots = {}
    plotrange = ()
    resplot = Instance(Plot)
    Model = Array
    Resids = Property(Array, depends_on='Model')
    y = {}
    # Define buttons for interface:
    add_profile = Button(label='Add component')
    remove_profile = Button(label='Remove latest')
    Go_Button = Button(label='Fit model')
    plwin = Instance(GridContainer)
    select = Str
    line_center = Float()

    def _line_center_changed(self):
        self.build_plot()  # Define a lightweight'er "update_plot" func later.
        # Set ranges to change automatically when plot values change.
        #plot.value_range.low_setting,\
        #    plot.value_range.high_setting = (-minval, maxval)
        #plot.index_range.low_setting,\
        #    plot.index_range.high_setting = (self.line_center - 20.,
        #                                     self.line_center + 20.)
        #resplot.value_range.low_setting,\
        #    resplot.value_range.high_setting = (-resmin, resmin)
        #resplot.index_range.low_setting,\
        #    resplot.index_range.high_setting = (plot.index_range.low_setting,
        #                                        plot.index_range.high_setting)


    def _get_Resids(self):
        intmod = sp.interp(self.x, self.mod_x, self.Model)
        resids = (self.indata - intmod) / self.errs
        return resids

    def _Components_default(self):
        return {'Contin': self.continuum_estimate,
                'Comp1': [0., .1, 0., 'a']}  # Center, Sigma, Height, Identifier

    def _Locks_default(self):
        return {'Comp1': [False, False, False, False]}

    def _CompoList_default(self):
        return ['Comp1']

    def _y_default(self):
        return {}

    def _select_default(self):
        return 'Comp1'

    def build_plot(self):
        global plotdata  # FIXME: Don't use global!
        onearray = Array
        onearray = sp.ones(self.indata.shape[0])
        minuses = onearray * (-1.)

        # Define index array for fit function:
        self.mod_x = sp.arange(self.line_center - 50.,
                               self.line_center + 50., .01)
        self.Model = sp.zeros(self.mod_x.shape[0])

        # Establish continuum array in a way that opens for other, more
        #   elaborate continua.
        self.contarray = sp.ones(self.mod_x.shape[0]) * \
                self.Components['Contin']
        self.y = {}

        for comp in self.CompoList:
            self.y[comp] = stats.norm.pdf(
                self.mod_x,
                self.Components[comp][0] + self.line_center,
                self.Components[comp][1]
            ) * self.Components[comp][1] * sp.sqrt(2. * sp.pi) * \
                    self.Components[comp][2]
        self.Model = self.contarray + self.y[self.select]

        broca = BroadcasterTool()

        # Define the part of the data to show in initial view:
        plotrange = sp.where((self.x > self.line_center - 20) &
                             (self.x < self.line_center + 20))
        # Define the y axis max value in initial view (can be panned/zoomed):
        maxval = float(self.indata[plotrange].max() * 1.1)
        minval = maxval / 15.
        maxerr = self.errs[plotrange].max() * 1.3
        resmin = max(self.Resids[plotrange].max(), 5.) * 1.2
        cenx = sp.array([self.line_center, self.line_center])
        ceny = sp.array([-minval, maxval])
        cenz = sp.array([-maxval, maxval])
        # Build plot of data and model
        plotdata = ArrayPlotData(
            wl=self.x,
            data=self.indata,
            xs=self.mod_x,
            cont=self.contarray,
            ones=onearray,
            minus=minuses,
            model=self.Model,
            errors=self.errs,
            ceny=ceny,
            cenz=cenz,
            cenx=cenx,
            Residuals=self.Resids,
        )
        for comp in self.CompoList:
            plotdata.set_data(comp, self.y[comp])
        olplot = GridContainer(shape=(2, 1), padding=10,
                               fill_padding=True,
                               bgcolor='transparent',
                               spacing=(5, 10))
        plot = Plot(plotdata)
        plot.y_axis.title = 'Flux density'
        resplot = Plot(plotdata, tick_visible=True, y_auto=True)
        resplot.x_axis.title = u'Wavelength [Å]'
        resplot.y_axis.title = u'Residuals/std. err.'

        # Create initial plot: Spectrum data, default first component,
        #   default total line profile.
        self.comprenders = []
        self.datarender = plot.plot(('wl', 'data'), color='black',
                                    name='Data',
                                    render_style='connectedhold')
        self.contrender = plot.plot(('xs', 'cont'), color='darkgray',
                                    name='Cont')
        self.modlrender = plot.plot(('xs', 'model'), color='blue',
                                    line_width=1.6, name='Model')
        self.centrender = plot.plot(('cenx', 'ceny'),
                                    color='black',
                                    type='line',
                                    line_style='dot',
                                    name='Line center',
                                    line_width=1.)

        # There may be an arbitrary number of gaussian components, so:
        for comp in self.CompoList:
            self.comprenders.append(
                plot.plot(('xs', comp),
                          type='line',
                          #color='auto',
                          color=tuple(COLOR_PALETTE[self.CompNum]),
                          line_style='dash',
                          name=comp)
            )

        # Create panel with residuals:
        resplot.plot(('wl', 'Residuals'), color='black', name='Resids')
        resplot.plot(('wl', 'ones'), color='green')
        resplot.plot(('wl', 'minus'), color='green')
        resplot.plot(('cenx', 'cenz'), color='red',
                     type='line',
                     line_style='dot',
                     line_width=.5)
        plot.x_axis.visible = False

        # Set ranges to change automatically when plot values change.
        plot.value_range.low_setting,\
            plot.value_range.high_setting = (-minval, maxval)
        plot.index_range.low_setting,\
            plot.index_range.high_setting = (self.line_center - 20.,
                                             self.line_center + 20.)
        resplot.value_range.low_setting,\
            resplot.value_range.high_setting = (-resmin, resmin)
        resplot.index_range.low_setting,\
            resplot.index_range.high_setting = (plot.index_range.low_setting,
                                                plot.index_range.high_setting)
        #resplot.index_range = plot.index_range  # Yes or no? FIXME
        plot.overlays.append(ZoomTool(plot, tool_mode='box',
                                      drag_button='left',
                                      always_on=False))

        resplot.overlays.append(ZoomTool(resplot, tool_mode='range',
                                         drag_button='left',
                                         always_on=False))

        # List of renderers to tell the legend what to write
        self.plots['Contin'] = self.contrender
        self.plots['Center'] = self.centrender
        self.plots['Model'] = self.modlrender
        for i in sp.arange(len(self.comprenders)):
            self.plots[self.CompoList[i]] = self.comprenders[i]

        # Build Legend:
        legend = Legend(component=plot, padding=10, align="ur")
        legend.tools.append(LegendTool(legend, drag_button="right"))
        legend.plots = self.plots
        plot.overlays.append(legend)
        olplot.tools.append(broca)
        pan = PanTool(plot)
        respan = PanTool(resplot, constrain=True, constrain_direction='x')
        broca.tools.append(pan)
        broca.tools.append(respan)
        plot.overlays.append(ZoomTool(plot, tool_mode='box',
                                      always_on=False))
        olplot.add(plot)
        olplot.add(resplot)
        olplot.components[0].set(resizable='hv', bounds=[500, 400])
        olplot.components[1].set(resizable='h', bounds=[500, 100])
        self.plot = plot
        self.resplot = resplot
        self.plwin = olplot
        self.legend = legend
        self.plotrange = plotrange

    # Select component 1 by default:

    def __init__(self, wavlens, indata, inerrs, linecen, fitrange=None):
        self.fitrange = fitrange
        if fitrange is None:
            self.fitrange = ()
        # Define index array for data:
        self.x = wavlens
        self.indata = indata
        self.errs = inerrs
        self.line_center = linecen
        self.build_plot()

    ### =======================================================================
    #     Reactive functions: What happens when buttons are pressed, parameters
    #     are changes etc.

    # Add component to model

    def _add_profile_fired(self):  # FIXME Don't set plotdata global!
        """ Add new component to model
        """
        global plotdata
        self.CompNum += 1
        Name = 'Comp' + str(self.CompNum)
        self.CompoList.append(Name)
        print "Added component nr. " + Name
        self.Components[Name] = [0., .1, 0., chr(self.CompNum+96)]
        self.Locks[Name] = [False, False, False, False]
        self.select = Name
        # And the plotting part:
        #    Add y array for this component.
        self.y[self.select] = stats.norm.pdf(
            self.mod_x,
            self.Centr + self.line_center,
            self.Sigma) * self.Sigma * sp.sqrt(2. * sp.pi) * self.Heigh
        plotdata[self.select] = self.y[self.select]
        render = self.plot.plot(('xs', self.select), type='line',
                                line_style='dash',
                                color=tuple(COLOR_PALETTE[self.CompNum]),
                                name=Name)
        self.plots[self.select] = render
        self.legend.plots = self.plots
        return

    def _remove_profile_fired(self):
        """ Remove the last added component
        """
        global plotdata
        if self.CompNum > 1:
            oldName = 'Comp' + str(self.CompNum)
            newName = 'Comp' + str(self.CompNum - 1)
            self.plot.delplot(oldName)
            plotdata.del_data(oldName)
            del self.y[oldName]
            del self.plots[oldName]
            del self.Components[oldName]
            del self.Locks[oldName]
            self.select = newName
            print 'Removed component nr. ' + str(self.CompNum)
            self.legend.plots = self.plots
            self.CompoList.pop()
            self.CompNum -= 1
        else:
            print 'No more components to remove'

    def _Go_Button_fired(self):

        # Transform the internal dict holding the model to a Pandas dataframe
        # that the lmfit wrapper will digest:
        tofit = pd.DataFrame.from_dict(self.Components).T
        print tofit
        print tofit.index
        tofit.columns = ['Pos', 'Sigma', 'Ampl', 'Identifier']
        tofit['Line center'] = self.line_center
        tofit.set_value('Contin', 'Lock', self.LockConti)
        for lines in self.Components.keys():
            if lines == 'Contin':
                continue
            tofit.set_value(lines, 'Lock', self.Locks[lines][:3])
        print self.Components
        print tofit

        # Make sure no dataarrays belonging to the parent class are altered.
        x = self.x.copy()
        data = self.indata.copy()
        errs = self.errs.copy()

        # If fitranges set, use them, otherwise use all data.
        fitrange = []
        if len(self.fitrange) == 0:
            self.fitrange = [(self.line_center - 15, self.line_center + 15)]
            print 'Fitrange 1: ', fitrange
        if len(self.fitrange) > 0:
            for ran in self.fitrange:
                print 'Ran', ran
                rmin, rmax = ran[0], ran[1]
                fitrange += sp.where((self.x > rmin) & (self.x < rmax))
            fitrange = sp.hstack(fitrange[:])
            fitrange.sort()
            x = x[fitrange]
            data = data[fitrange]
            errs = errs[fitrange]
        try:
            import lmfit_wrapper as lw
        except ImportError:
            print 'Could not import LMfit'
            return

        params = lw.load_params(tofit)
        result = lw.fit_it(
            #lw.build_model,
            params,
            args=(self.x[fitrange], self.indata[fitrange], self.errs[fitrange]))

        output = lw.params_to_grism(params, output_format='df')
        print output
        print tofit
        output['Identifier'] = tofit['Identifier']
        print output
        output.set_value('Contin', 'Identifier', sp.float64('nan'))
        output['Pos'] -= tofit['Line center']
        print output

        outdict = {}
        for i in output.index:
            row = output.ix[i]
            if i == 'Contin':
                outdict[i] = row['Ampl']
            else:
                outdict[i] = [row['Pos'], row['Sigma'], \
                              row['Ampl'], row['Identifier']]
        self.Components = outdict
        self.import_model()
        return result

    # Define what to do when a new component is selected.

    def _select_changed(self):
        # First, show the values of current component in sliders!
        #print "Selection changed, selected component is: " + self.select
        self.Centr = self.Components[self.select][0]
        self.Sigma = self.Components[self.select][1]
        self.Heigh = self.Components[self.select][2]
        self.LockCentr = self.Locks[self.select][0]
        self.LockSigma = self.Locks[self.select][1]
        self.LockHeigh = self.Locks[self.select][2]
        self.plot.request_redraw()
        return

    # Every time one of the parameters in the interactive window is changed,
    #   write the change to the parameters list of the selected component.
    #   Do this one-by-one, as it is otherwise going to mess up the
    #   creation and selection of components.

    def _Centr_changed(self):
        self.Components[self.select][0] = self.Centr
        self.update_plot()
        return

    def _Sigma_changed(self):
        self.Components[self.select][1] = self.Sigma
        self.update_plot()
        return

    def _Heigh_changed(self):
        self.Components[self.select][2] = self.Heigh
        self.update_plot()
        return

    def _continuum_estimate_changed(self):
        self.Components['Contin'] = self.continuum_estimate
        self.update_plot()

    def _LockCentr_changed(self):
        self.Locks[self.select][0] = self.LockCentr
        return

    def _LockSigma_changed(self):
        self.Locks[self.select][1] = self.LockSigma
        return

    def _LockHeigh_changed(self):
        self.Locks[self.select][2] = self.LockHeigh
        return

    ###========================================================================
    # Define the graphical user interface

    view = View(
        Group(
            Group(
                Item('plwin', editor=ComponentEditor(),
                     show_label=False, springy=True),
                Group(
                    Group(
                        Item('Centr', label='Center',
                             enabled_when='LockCentr==False'),
                        Item('Sigma', label='Sigma',
                             enabled_when='LockSigma==False'),
                        Item('Heigh', label=u'Strength ',
                             enabled_when='LockHeigh==False'),
                        springy=True,
                        show_border=False,
                        orientation='vertical'),
                    Group(
                        Item('LockCentr', label='Lock'),
                        Item('LockSigma', label='Lock'),
                        Item('LockHeigh', label='Lock'),
                        orientation='vertical'),
                    orientation='horizontal',
                    show_border=True,
                    label='Component parameters'),
                HGroup(
                    Item('continuum_estimate',
                         enabled_when='LockConti==False',
                         label='Contin.  ',
                         springy=True
                         ),
                    Item('LockConti', label='Lock'),
                    show_border=True,
                    springy=True
                ),
                show_border=True),
            Group(Item('add_profile'),
                  Item('remove_profile'),
                  Item('Go_Button'),
                  Item('Feedback', style='readonly'),
                  Item('Feedback', style='readonly'),
                  Item(name='select', editor=EnumEditor(name='CompoList'),
                       style='custom'),
                  orientation='vertical',
                  show_labels=False,
                  show_border=True),
            orientation='horizontal'),
        resizable=True,
        height=700, width=1000,  # ),
        buttons=[UndoButton, ApplyButton, CancelButton, OKButton],
        close_result=True,
        kind='live',  # Works but not perfect.
        title="Grism - line profile editor")

    def import_model(self):
        #global plotdata
        #print '    '
        #print 'This is the model importing method of ProfileEditor: '
        self.CompoList = sorted(self.Components.keys())[:-1]
        self.CompNum = len(self.CompoList)
        #for com in self.CompoList:
        #    self.Locks[com] = [False] * 4
        #print self.Components
        self.continuum_estimate = self.Components['Contin']
        self._select_changed()
        self.build_plot()
        self.update_plot()
        print '    '


    def update_plot(self):
        self.y[self.select] = stats.norm.pdf(
            self.mod_x,
            self.Centr + self.line_center,
            self.Sigma) * self.Sigma * sp.sqrt(2. * sp.pi) * self.Heigh
        ys = sp.asarray(self.y.values()).sum(0)
        self.contarray = sp.ones(self.mod_x.shape[0]) * self.continuum_estimate
        self.Model = self.contarray + ys
        plotdata.set_data('cont', self.contarray)
        plotdata.set_data(self.select, self.y[self.select])
        plotdata.set_data('model', self.Model)
        plotdata.set_data('Residuals', self.Resids)

    @on_trait_change('Resids')
    def update_resid_window(self):
        resmin = max(self.Resids[self.plotrange].max(), 5.) * 1.2
        self.resplot.value_range.low_setting,\
            self.resplot.value_range.high_setting = (-resmin, resmin)
        self.resplot.request_redraw()

###===========================================================================
#            Convenience- and helper functions
###===========================================================================

def load_profile(dataframe, centroid):
    wave = dataframe['wave'].values
    data = dataframe['data'].values
    errs = dataframe['errs'].values

    lb = ProfileEditor(wave, data, errs, centroid)
    lb.configure_traits()
    return lb

###============================================================================
#             What to do if script is called from command line rather than
#             imported from an external python module or shell.
###============================================================================

if __name__ == '__main__':
    df = pd.read_csv('testdata1d.dat')
    pe = load_profile(df, 6699.5)

