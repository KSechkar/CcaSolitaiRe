# PLATE_READER_ANTOOLS.PY - analysis tools for plate reader data
# By Kirill Sechkar

# PACKAGE IMPORTS
import numpy as np

from bokeh import plotting as bkplot, models as bkmodels, layouts as bklayouts, io as bkio, transform as bktransform, palettes as bkpalettes
from bokeh.colors import RGB as bkRGB
from bokeh import io as bkio
import bokeh.resources as bkresources


# AUXILIARY FUNCTIONS --------------------------------------------------------------------------------------------------

# sample path length for a given plate
def get_path_length(sample_volume,  # volume of the sample in the plate, uL=mm^3
                    plate  # plate type, currently only 'Grenier' is supported (96-well plate with 34mm^2 well area)
                    ):

    # well top and bottom diameters & total heights, mm
    well_d_top = {'Grenier': 6.96, 'Costar': 6.86}  # top diameter of the well, mm
    well_d_bot = {'Grenier': 6.39, 'Costar': 6.35}  # bottom diameter of the well, mm
    well_height = {'Grenier': 10.9, 'Costar': 10.67}  # height of a Grenier 96-well plate well, mm

    # get well top and bottom radii and height
    r_top = well_d_top[plate] / 2
    r_bot = well_d_bot[plate] / 2
    height = well_height[plate]

    # get an auxiliary coefficient
    D = (r_top - r_bot) / height

    # path length is the real positive root of a cubic polynomial
    polypath_roots = np.polynomial.polynomial.polyroots(
        np.array([
            -sample_volume,  # constant term
            np.pi * r_bot ** 2,  # linear term
            np.pi * D * r_bot,  # quadratic term
            np.pi * D ** 2 / 3  # cubic term
        ])
    )
    path_length_mm =np.real(max(polypath_roots[np.isreal(polypath_roots)]))

    # return the path length in cm
    return path_length_mm/10

# DATA STRUCTURES FOR THE MEASUREMENTS ---------------------------------------------------------------------------------
# blanks
class Blanks:
    def __init__(self, timepoints,
                 gen_config, compet_ind, # genetic configuration, competitor induction for the sample
                 probe_colour='darkgreen', compet_colour='crimson'  # probe and competitor plot colours
                ):
        # title summarising the sample considered
        self.title = gen_config + ', ' + compet_ind

        # measurement timepoints
        self.timepoints = timepoints
        num_timepoints = len(timepoints)
        # OD600
        self.od600 = np.zeros((1,num_timepoints))
        # probe fluorescence
        self.probe = np.zeros((1,num_timepoints))
        # compet fluorescence
        self.compet = np.zeros((1,num_timepoints))

        # colours for plotting
        self.probe_colour = probe_colour
        self.compet_colour = compet_colour
        return

    # figure for different measurements over time
    def mot_fig(self,
                what,   # which measurement to plot, 'od600', 'probe' or 'compet'
                y_range=None
        ):
        rep_dashes = ['solid', 'dashed', 'dotted']  # dashes for the repeats

        # determine the y-axis label
        ylabel = what


        # create the figure
        fig = bkplot.figure(title=self.title,
                            x_axis_label='Time, h', y_axis_label=ylabel,
                            tools='pan,box_zoom,reset,save',
                            width=155, height=150)
        fig.output_backend = 'svg'
        # plot individual repeats
        for i in range(0, self.od600.shape[0]):
            if (what == 'od600'):
                fig.line(self.timepoints, self.od600[0,:],
                         line_color='black',
                         line_alpha=0.3, line_width=2, line_dash=rep_dashes[i])
            elif (what == 'probe'):
                fig.line(self.timepoints, self.probe[0,:],
                         line_color=self.probe_colour,
                         line_alpha=0.5, line_width=2, line_dash=rep_dashes[i])
            elif (what == 'compet'):
                fig.line(self.timepoints, self.compet[0,:],
                         line_color=self.compet_colour,
                         line_alpha=0.5, line_width=2, line_dash=rep_dashes[i])


        # axis settings
        fig.xaxis.ticker = bkmodels.BasicTicker(desired_num_ticks=3)

        # font settings
        fig.title.text_font_size = '8pt'
        fig.yaxis.axis_label_text_font_size = '8pt'

        # set y axis range if specified
        if (y_range != None):
            fig.y_range = bkmodels.Range1d(y_range[0], y_range[1])

        # return the figure
        return fig

# sample measurements on the plate
class PlateMeas:
    # INITIALISATION
    # initialise
    def __init__(self,
                 num_reps,  # number of replicates
                 timepoints,    # measurement timepoints, h
                 gen_config, compet_ind, light,  # genetic configuration, competitor induction and light conditions for the sample
                 probe_colour='darkgreen', compet_colour='crimson' # probe and competitor plot colours
                 ):
        # title summarising the sample considered
        self.title = gen_config + ', ' + compet_ind + ', ' + light

        # measurement timepoints
        self.timepoints = timepoints
        num_timepoints = len(timepoints)
        # OD600
        self.od600_raw = np.zeros((num_reps,num_timepoints))
        self.od600 = np.zeros((num_reps,num_timepoints))
        self.od600_mean = np.zeros(num_timepoints)
        # probe fluorescence
        self.probe_raw = np.zeros((num_reps,num_timepoints))
        self.probe = np.zeros((num_reps,num_timepoints))
        self.probe_mean = np.zeros(num_timepoints)
        # compet fluorescence
        self.compet_raw = np.zeros((num_reps,num_timepoints))
        self.compet = np.zeros((num_reps,num_timepoints))
        self.compet_mean = np.zeros(num_timepoints)

        # mask for different repeats' measurements - originally all True
        self.measmask = np.ones((num_reps,num_timepoints), dtype=bool)
        # means using masked measurements only
        self.od600_mean_masked = np.zeros(num_timepoints)
        self.probe_mean_masked = np.zeros(num_timepoints)
        self.compet_mean_masked = np.zeros(num_timepoints)
        # timepoints for the measurements using masked measurements only
        self.mean_masked_timepoints = np.zeros(num_timepoints)

        # logistic growth curve fittinbg outcomes
        self.od600_0 = np.zeros(num_reps)   # initial OD600
        self.k = np.zeros(num_reps) # medium carrying capacity
        self.l = np.zeros(num_reps) # exponential growth rate

        # measurements for the timepoint closest to exponential growth
        self.egt_idcs = np.zeros(num_reps,dtype=np.int32)  # timepoint indices closest to exponential growth
        self.probe_egt = np.zeros(num_reps)  # probe fluorescence at the exponential growth timepoint
        self.compet_egt = np.zeros(num_reps)  # compet fluorescence at the exponential growth timepoint
        self.probe_egt_mean = 0.0  # mean probe fluorescence at the exponential growth timepoint
        self.compet_egt_mean = 0.0 # mean compet fluorescence at the exponential growth timepoint

        # plotting styles
        self.probe_colour = probe_colour
        self.compet_colour = compet_colour
        return

    # RAW DATA PROCESSING
    # find the culture's OD600 from raw measurements and blanks
    def find_od600(self,
                   blank, # blank measurement object
                   sample_volume,    # volume of the sample in the plate, uL=mm^3
                   plate = 'Grenier'  # plate type, currently only 'Grenier' is supported (96-well plate with 34mm^2 well area)
                   ):
        # if sample volume is a single value, use it for all timepoints
        if(np.isscalar(sample_volume)):
            # get the OD with path length correction and clipping for positivity
            self.od600 = (self.od600_raw-blank.od600)
            self.od600 = np.clip(np.float64(self.od600), 0.005*np.ones(self.od600.shape), None)
        # if sample volumes may differ, go through each measurement individually
        else:
            for i in range(0,len(sample_volume)):
                for j in range(0,self.od600_raw.shape[0]):
                    # get the OD with path length correction and clipping for positivity
                    self.od600[j,i] = (self.od600_raw[j,i]-blank.od600[0,i])
                    self.od600[j,i] = np.clip(np.float64(self.od600[j,i]), 0.005, None)

        # get the mean
        self.od600_mean = np.mean(self.od600, axis=0)
        return

    # find the probe and compet fluorescence intensities normalised to od600
    def find_normflus(self,
                      blank,  # blank measurement object
                      sample_volume  # volume of the sample in the plate, uL=mm^3
                      ):

        # get the normalised probe fluorescence with blanks subtracted, clipped for non-negativity
        self.probe = (self.probe_raw - blank.probe) / self.od600
        self.probe = np.clip(np.float64(self.probe), np.zeros(self.od600.shape), None)

        # get the mean probe fluorescence
        self.probe_mean = np.mean(self.probe, axis=0)

        # get the normalised compet fluorescence with blanks subtracted, clipped for non-negativity
        self.compet = (self.compet_raw - blank.compet) / self.od600
        self.compet = np.clip(np.float64(self.compet), np.zeros(self.od600.shape), None)

        # get the mean compet fluorescence
        self.compet_mean = np.mean(self.compet, axis=0)

        return

    # recalculate the mean values using unmasked measurements only
    def find_masked_means(self):
        # initialise masked means from a clean slate
        self.od600_mean_masked = np.zeros(self.od600_mean.shape)
        self.probe_mean_masked = np.zeros(self.od600_mean.shape)
        self.compet_mean_masked = np.zeros(self.od600_mean.shape)
        # initialise masked stdev from a clean slate
        self.od600_stdev_masked = np.zeros(self.od600_mean.shape)
        self.probe_stdev_masked = np.zeros(self.od600_mean.shape)
        self.compet_stdev_masked = np.zeros(self.od600_mean.shape)

        num_measurements_included = np.zeros(
            self.od600.shape[1])  # how many measurements included in the mean for each timepoint
        for i in range(0, self.od600.shape[1]):
            # get means
            for rep in range(1, self.od600_0.shape[0]+1):
                if(self.measmask[rep-1, i]):
                    num_measurements_included[i] += 1
                    self.od600_mean_masked[i] += self.od600[rep-1, i]
                    self.probe_mean_masked[i] += self.probe[rep-1, i]
                    self.compet_mean_masked[i] += self.compet[rep-1, i]
            if(num_measurements_included[i] > 0):
                self.od600_mean_masked[i] /= num_measurements_included[i]
                self.probe_mean_masked[i] /= num_measurements_included[i]
                self.compet_mean_masked[i] /= num_measurements_included[i]
                self.mean_masked_timepoints[i] = self.timepoints[i]
            # get stdevs
            for rep in range(1, self.od600_0.shape[0]+1):
                if(self.measmask[rep-1, i]):
                    num_measurements_included[i] += 1
                    self.od600_stdev_masked[i] += (self.od600[rep-1, i]-self.od600_mean_masked[i])**2
                    self.probe_stdev_masked[i] += (self.probe[rep-1, i]-self.probe_mean_masked[i])**2
                    self.compet_stdev_masked[i] += (self.compet[rep-1, i]-self.compet_mean_masked[i])**2
            if(num_measurements_included[i] > 1):
                self.od600_stdev_masked[i] = np.sqrt(self.od600_stdev_masked[i]/num_measurements_included[i])
                self.probe_stdev_masked[i] = np.sqrt(self.probe_stdev_masked[i]/num_measurements_included[i])
                self.compet_stdev_masked[i] = np.sqrt(self.compet_stdev_masked[i]/num_measurements_included[i])

        # delete timepoints with no measurements at all
        self.mean_masked_timepoints = self.mean_masked_timepoints[num_measurements_included > 0]
        self.od600_mean_masked = self.od600_mean_masked[num_measurements_included > 0]
        self.probe_mean_masked = self.probe_mean_masked[num_measurements_included > 0]
        self.compet_mean_masked = self.compet_mean_masked[num_measurements_included > 0]

        return

    # PLOTTING
    # figure for different measurements over time
    def mot_fig(self,
                what,   # which measurement to plot, 'od600', 'probe' or 'compet'
                y_range=None,
                masked=False, # whether to mask weird measurements
        ):
        rep_dashes = ['solid', 'dashed', 'dotted']  # dashes for the repeats

        if(masked):
            measmask = self.measmask
        else:
            measmask = np.ones(self.od600.shape, dtype=bool)

        # determine the y-axis label
        if (what == 'od600'):
            ylabel = what
        else:
            ylabel = what + ' flu./OD600'

        # create the figure
        fig = bkplot.figure(title=self.title,
                            x_axis_label='Time, h', y_axis_label=ylabel,
                            tools='pan,box_zoom,reset,save',
                            width=155, height=150)
        fig.output_backend = 'svg'
        # plot individual repeats
        for i in range(0, self.od600.shape[0]):
            if (what == 'od600'):
                fig.line(self.timepoints[measmask[i, :]], self.od600[i, measmask[i, :]],
                         line_color='black',
                         line_alpha=0.3, line_width=2, line_dash=rep_dashes[i])
            elif (what == 'probe'):
                fig.line(self.timepoints[measmask[i, :]], self.probe[i, measmask[i, :]],
                         line_color=self.probe_colour,
                         line_alpha=0.5, line_width=2, line_dash=rep_dashes[i])
            elif (what == 'compet'):
                fig.line(self.timepoints[measmask[i, :]], self.compet[i, measmask[i, :]],
                         line_color=self.compet_colour,
                         line_alpha=0.5, line_width=2, line_dash=rep_dashes[i])

        # plot the mean
        if(masked):
            if (what == 'od600'):
                fig.line(self.mean_masked_timepoints, self.od600_mean_masked,
                         line_color='black',
                         line_width=2)
            elif (what == 'probe'):
                fig.line(self.mean_masked_timepoints, self.probe_mean_masked,
                         line_color=self.probe_colour,
                         line_width=2)
            elif (what == 'compet'):
                fig.line(self.mean_masked_timepoints, self.compet_mean_masked,
                         line_color=self.compet_colour,
                         line_width=2)
        else:
            if (what == 'od600'):
                fig.line(self.timepoints, self.od600_mean,
                         line_color='black',
                         line_width=2)
            elif (what == 'probe'):
                fig.line(self.timepoints, self.probe_mean,
                         line_color=self.probe_colour,
                         line_width=2)
            elif (what == 'compet'):
                fig.line(self.timepoints, self.compet_mean,
                         line_color=self.compet_colour,
                         line_width=2)

        # axis settings
        fig.xaxis.ticker = bkmodels.BasicTicker(desired_num_ticks=3)

        # font settings
        fig.title.text_font_size = '8pt'
        fig.yaxis.axis_label_text_font_size = '8pt'

        # set y axis range if specified
        if (y_range != None):
            fig.y_range = bkmodels.Range1d(y_range[0], y_range[1])

        # return the figure
        return fig

# READING PLATE READER DATA --------------------------------------------------------------------------------------------
# based on the condition/replicate-to-well-matches, get a translator from well names to conditions/replicates, and a blank translator
def get_w2cr_w2b(cr2w):
    w2cr = {}
    w2b = {}
    for gen_config in cr2w.keys():
        if (gen_config != 'Blanks'):
            for compet_ind in cr2w[gen_config].keys():
                for light in cr2w[gen_config][compet_ind].keys():
                    for rep in cr2w[gen_config][compet_ind][light].keys():
                        well = cr2w[gen_config][compet_ind][light][rep]
                        w2cr[well] = (gen_config, compet_ind, light, rep)
        else:
            for blank_gen_config in cr2w['Blanks'].keys():
                for compet_ind in cr2w['Blanks'][blank_gen_config].keys():
                    well = cr2w['Blanks'][blank_gen_config][compet_ind]
                    w2b[well] = (blank_gen_config, compet_ind)
    return w2cr, w2b


# reading data from the plate reader upstairs - for the advanced experiments
def get_block_of_measurements_upstairs_adv(plate_log,  # plate log dataframe
                                           timepoint_cntr,  # timepoint counter as we go through the files
                                           w2b,  # dictionary of well to blanks mapping, {well: (gen_config, compet_ind)}
                                           w2cr, # dictionary of well to culture replicates mapping, {well: (gen_config, compet_ind, light, rep)}
                                           pms, blnks,  # dictionaries of PlateMeas and Blanks objects to be filled
                                           what2block, # mapping of measurement types to blocks in record
                                           what,  # which measurement to read, 'od600', 'probe' or 'compet'
                                           postfix = ''  # postfix to well indices if the same timepoint's measurements are on different plates
                                           ):
    # determine what we are looking for, see what is the corresponding block in record
    block_in_record = what2block[what]

    # locate the log row with well NUMBERS for the corresponding block
    wellnos_log_rows_found = 0 # how many well log rows already found
    wellnos_log_row_cntr = -1 # row counter as we look for the right well log row
    while(wellnos_log_rows_found <= block_in_record):
        wellnos_log_row_cntr+=1
        if(plate_log.iloc[wellnos_log_row_cntr,0]=='<>'):
            wellnos_log_rows_found+=1
    wellnos_log_row=plate_log.iloc[wellnos_log_row_cntr,:]

    # go along the plate till we run out of well numbers
    wellnos_log_col_cntr = 1
    while(wellnos_log_col_cntr < len(wellnos_log_row)):
        if(float(wellnos_log_row.iloc[wellnos_log_col_cntr]) not in [1,2,3,4,5,6,7,8,9,10,11,12]):
            break
        else:
            wellnos_log_col_cntr += 1

    # go down the plate till we run out of well letters
    welllets_log_row_cntr=wellnos_log_row_cntr+1
    while (str(plate_log.iloc[welllets_log_row_cntr,0]) in 'ABCDEFGH'):
        welllets_log_row_cntr+=1

    # go through the columns and rows to get the blanks
    for wellno_cntr in range(1,wellnos_log_col_cntr): # cycling through well numbers
        for welllet_cntr in range(wellnos_log_row_cntr+1,welllets_log_row_cntr):
            well = str(plate_log.iloc[welllet_cntr,0])+str(int(wellnos_log_row.iloc[wellno_cntr]))+postfix
            if(well in w2b.keys()):
                # get the genetic configuration and competitor induction
                gen_config, compet_ind = w2b[well]
                # record the value
                if(what=='od600'):
                    blnks[gen_config][compet_ind].od600[0,timepoint_cntr] = np.clip(float(plate_log.iloc[welllet_cntr,wellno_cntr]),0.0, None)
                elif(what=='probe'):
                    blnks[gen_config][compet_ind].probe[0,timepoint_cntr] = np.clip(float(plate_log.iloc[welllet_cntr,wellno_cntr]), 0.0, None)
                elif(what=='compet'):
                    blnks[gen_config][compet_ind].compet[0,timepoint_cntr] = np.clip(float(plate_log.iloc[welllet_cntr,wellno_cntr]), 0.0, None)

    # go through the columns and rows to get the blanks
    for wellno_cntr in range(1,wellnos_log_col_cntr): # cycling through well numbers
        for welllet_cntr in range(wellnos_log_row_cntr+1,welllets_log_row_cntr):
            well = str(plate_log.iloc[welllet_cntr,0])+str(int(wellnos_log_row.iloc[wellno_cntr]))+postfix
            if(well in w2cr):
                # get the genetic configuration, competitor induction, light and repeat
                gen_config, compet_ind, light, rep = w2cr[well]
                # record the raw measurements
                if(what=='od600'):
                    pms[gen_config][compet_ind][light].od600_raw[rep-1,timepoint_cntr] = np.clip(float(plate_log.iloc[welllet_cntr,wellno_cntr]), 0.0, None)
                elif(what=='probe'):
                    pms[gen_config][compet_ind][light].probe_raw[rep-1,timepoint_cntr] = np.clip(float(plate_log.iloc[welllet_cntr,wellno_cntr]), 0.0, None)
                elif(what=='compet'):
                    pms[gen_config][compet_ind][light].compet_raw[rep-1,timepoint_cntr] = np.clip(float(plate_log.iloc[welllet_cntr,wellno_cntr]),0.0, None)

    # return updated pms and blnks
    return pms, blnks