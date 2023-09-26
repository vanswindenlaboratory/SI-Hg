"""
Author: Federica Gugole

This script has been developed for the data processing of the measurements 
according to the 19NRM03 SI-Hg protocol “Protocol for the SI-traceable calibration 
of elemental mercury gas generators used in the field”. 
This SI-Hg protocol specifies the procedures for establishing traceability 
to the SI units for the quantitative output of elemental mercury generators 
that are employed in regulatory applications for emission monitoring or testing. 
This protocol provides methods for:
-   Calibrating the output of a mercury gas generator by comparison 
    with a reference standard
-   Calculating the uncertainty of the mercury concentration generated 
    with the candidate generator in relation to the known uncertainty 
    of the reference standard.

The mercury concentration in a gas mixture prepared with a mercury gas generator 
is determined by comparison with a metrologically traceable reference standard 
to calibrate the output of a candidate generator.

The comparison can be performed at one concentration level (single-point calibration) 
or at several concentration levels (multipoint calibration) using a bracketing sequence. 
When applying the bracketing sequence, the outputs from the generators 
are introduced alternately to an analyser so that each response from the candidate generator 
is bracketed by a pair of responses from the reference generator. 
At each concentration level, the bracketing procedure requires a minimum 
of four responses from the certified reference standard generator 
and three responses from the candidate generator

The software is used to calculate different parameters based on the responses obtained 
with the bracketing measurement sequence: 
-   Output mercury concentration(s) of the candidate generator (Y(cand(i))) 
    at the setpoint(s) (c(cand(i)))
-   In case of a multi-point calibration: the interpolation function 
    for the output mercury concentration of the candidate generator 
    as a function of the setpoint
-   Deviation between the candidate generator setpoint 
    and the calculated output mercury concentration
-   Uncertainty of the calculated candidate generator output mercury concentration (U(Y(cand)))

The single point analysis is done following 19NRM03 SI-Hg protocol 
“Protocol for the SI-traceable calibration of elemental mercury gas generators 
used in the field”.

Everything else follows the steps defined in the SI-Hg protocol, thus all the procedures 
explained in the SI-Hg protocol apply also to this procedure. 
For more details on the mathematical formulas, please check the SI-Hg protocol.

The codes have been written trying to highlight the blocks of code 
performing some specific calculations described in the SI-Hg protocol. 
Variables and custom functions have been named trying to match the naming 
used in the SI-Hg protocol or, if not possible, a (hopefully) self-explanatory name 
has been assigned to the variable.
"""

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

"""
Import and general settings
"""

import xlwings as xw
import os
import sys
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
plt.rcParams['figure.figsize'] = 9,5
import seaborn as sns

from scipy.interpolate import interp1d

logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

"""
Define auxiliary functions
"""

def output_ratio(m_cand, m_refA, m_refB, t_cand, t_refA, t_refB):
    """
    Input:
        - m_cand : int, float, list or array; candidate generator measurement
        - m_refA : int, float, list or array; reference generator measurement after candidate
        - m_refB : int, float, list or array; reference generator measurement before candidate
        - t_cand : int, float, list or array; time of candidate generator measurement
        - t_refA : int, float, list or array; time of reference generator measurement after candidate
        - t_refB : int, float, list or array; time of reference generator measurement before candidate generator
    
    Output:
        output ratio (eq. 3b in the SI-Hg protocol)
    """
    # check if input variables are in the correct format
    input_vars = [m_cand, m_refA, m_refB, t_cand, t_refA, t_refB]
    for var in input_vars:
        if not isinstance(var, (int, float, list, np.ndarray)):
            print('Received variable of type {}. Acceptable formats are: int, float, list and numpy.ndarray.'.format(type(var)))
            # terminate execution
            exit()
    
    if isinstance(m_cand, list): 
        m_cand = np.array(m_cand)
    if isinstance(m_refA, list): 
        m_refA = np.array(m_refA)
    if isinstance(m_refB, list): 
        m_refB = np.array(m_refB)
    if isinstance(t_cand, list): 
        t_cand = np.array(t_cand)
    if isinstance(t_refA, list): 
        t_refA = np.array(t_refA)
    if isinstance(t_refB, list): 
        t_refB = np.array(t_refB)

    return m_cand / ( m_refB*(t_refA-t_cand)/(t_refA-t_refB) + m_refA*(t_cand-t_refB)/(t_refA-t_refB) )

# ----------------------------------------------------------------------------

def SE_stability(m_var, t_var, mt_cov, n):
    """
    Input:
        - m_var : int or float; measurements variance
        - t_var : int or float; time variance of measurements
        - mt_cov : int or float; covariance between measurement values and measurement times
        - n : int or float; number of measurements in the bracketing procedure per MuT or reference
    
    Output:
        stability standard error (Annex 1 of the SI-Hg protocol)
    """
    # check if input variables are in the correct format
    input_vars = [m_var, t_var, mt_cov, n]
    for var in input_vars:
        if not isinstance(var, (int, float)):
            print('Received variable of type {}. Acceptable formats are: int, float.'.format(type(var)))
            # terminate execution
            exit()

    return np.sqrt((m_var - mt_cov**2/t_var)/(n-2))

# ----------------------------------------------------------------------------

def u_stability(R, SE_cand, SE_refA, SE_refB, m_cand, m_refA, m_refB, t_cand, t_refA, t_refB):
    """
    Input:
        - R : int, float or array; output ratio of candidate measurement
        - SE_cand : int, float or array; standard error of candidate generator measurement
        - SE_refA : int, float or array; standard error or reference generator measurement after candidate
        - SE_refB : int, float or array; standard error of reference generator measurement before candidate
        - m_cand : int, float or array; candidate generator measurement
        - m_refA : int, float or array; reference generator measurement after candidate
        - m_refB : int, float or array; reference generator measurement before candidate
        - t_cand : int, float or array; time of candidate generator measurement
        - t_refA : int, float or array; time of reference generator measurement after candidate
        - t_refB : int, float or array; time of reference generator measurement before candidate generator
    
    Output:
        stability uncertainty component (Annex 1 of the SI-Hg protocol)
    """
    # check if input variables are in the correct format
    input_vars = [R, SE_cand, SE_refA, SE_refB, m_cand, m_refA, m_refB, t_cand, t_refA, t_refB]
    for var in input_vars:
        if not isinstance(var, (int, float, np.ndarray)):
            print('Received variable of type {}. Acceptable formats are: int, float and numpy.ndarray.'.format(type(var)))
            # terminate execution
            exit()

    SE_B2 = (SE_refB / ( m_refB*(t_refA-t_refB)/(t_refA-t_cand)) )**2
    SE_cand2 = (SE_cand/m_cand)**2
    SE_A2 = (SE_refA / ( m_refA*(t_refA-t_refB)/(t_cand-t_refB)) )**2
    return R*np.sqrt(SE_B2 + SE_cand2 + SE_A2)

# ----------------------------------------------------------------------------

def u_repeatability(K, L, R_cand, uncer_stab):
    """
    Input:
        - K : int; number of candidate measurements per bracketing procedure
        - L : int; see Annex 1 of the SI-Hg protocol
        - R_cand : list or array; output ratio of candidate measurements
        - uncer_stab : float; stability uncertainty component as calculated by function u_stability
    
    Output:
        repeatability uncertainty component (Annex 1 of the SI-Hg protocol)
    """
    # check if input variables are in the correct format
    input_vars = [K, L]
    for var in input_vars:
        if not isinstance(var, int):
            print('Received variable of type {}. Acceptable formats are: int.'.format(type(var)))
            # terminate execution
            exit()

    if not isinstance(R_cand, (list, np.ndarray)):
        print('Received variable of type {}. Acceptable formats are: list and numpy.ndarray.'.format(type(R_cand)))
        # terminate execution
        exit()
    
    if not isinstance(uncer_stab, float):
        print('Received variable of type {}. Acceptable formats are: float.'.format(type(uncer_stab)))
        # terminate execution
        exit()
    
    R_cand_avg = np.mean(R_cand)
    s1_squared = K*uncer_stab**2
    s2_squared = np.sum((R_cand-R_cand_avg)**2)/(K-1)
    return np.sqrt(np.max([0, s2_squared-s1_squared/L]))/np.sqrt(K)

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------


def main():
    # ----------------------------------------------------------------------------

    """
    Set up log file and file with console output (i.e. what Python would normally print to console)
    """

    # get current working directory
    dir_path = os.path.dirname(os.path.realpath(__file__))

    # create log file
    logfilename = os.path.join(dir_path, 'logfile_main.log')
    logging.basicConfig(filename=logfilename,
                    format='%(levelname)s %(asctime)s :: %(message)s',
                    filemode='w',
                    level=logging.DEBUG)

    logging.debug('Inside of main()')

    # write console outputs to file (for debugging purposes)
    sys.stdout = open(os.path.join(dir_path, 'console_output_main.txt'), 'w')

    # create Excel workbook
    wb = xw.Book.caller()

    # ----------------------------------------------------------------------------
    """Read input info"""
    data_sheet = 'Data'  

    ref_concentration = wb.sheets['Input_info'].range('D2').value
    cand_concentration = wb.sheets['Input_info'].range('D3').value
    extended_uncertainty_ref_generator = wb.sheets['Input_info'].range('D4').value/100.
    u_ref_concentration = ref_concentration*extended_uncertainty_ref_generator/2.0

    datetime_column = 'DateTime' 
    name_column = 'Name' 
    valid_column = 'Valid' 
    value_column = 'Pk Area' 

    Ref_A = wb.sheets['Input_info'].range('D5').value
    Ref_B = wb.sheets['Input_info'].range('D6').value
    Cand_A = wb.sheets['Input_info'].range('D7').value
    Cand_B = wb.sheets['Input_info'].range('D8').value
    Zero_A = wb.sheets['Input_info'].range('D9').value
    Zero_B = wb.sheets['Input_info'].range('D10').value

    apply_zero_correction = bool(wb.sheets['Input_info'].range('D11').value)
    reproducibility_uncertainty = wb.sheets['Input_info'].range('D12').value/100.
    datetime_format = wb.sheets['Input_info'].range('D13').value

    # ----------------------------------------------------------------------------
    """Load data"""
    # cel = wb.sheets[data_sheet].range('A1')
    df = wb.sheets[data_sheet].range('A1').options(pd.DataFrame, header=1, index=False, expand='table').value
    print(df.info())
    if isinstance(df[datetime_column], object):
        df[datetime_column] = pd.to_datetime(df[datetime_column], format=datetime_format)
        print(df.info())

    # ----------------------------------------------------------------------------
    """Pre-process data"""

    # convert Valid from int64 to boolean 
    df[valid_column] = df[valid_column].astype('bool')

    # if present, drop invalid measurements
    df = df.drop(df[df[valid_column] == False].index)

    # if present, drop nan values in value_column
    df = df.dropna(subset=[value_column])
    df = df.reset_index(drop=True)

    # convert datetime to number of seconds elapsed since the first recorded measurement 
    time_elapsed = df[datetime_column] - df[datetime_column].iloc[0]
    time_elapsed = time_elapsed.dt.total_seconds().values
    df['Elapsed_seconds'] = time_elapsed
    df['Elapsed_minutes'] = time_elapsed/60.0

    # find indices for each of the above ref, cand and zero
    idx_ref_A  = df[name_column].str.contains(Ref_A).fillna(value=False)
    idx_ref_B  = df[name_column].str.contains(Ref_B).fillna(value=False)
    idx_cand_A = df[name_column].str.contains(Cand_A).fillna(value=False)
    idx_cand_B = df[name_column].str.contains(Cand_B).fillna(value=False)
    idx_zero_A = df[name_column].str.contains(Zero_A).fillna(value=False)
    idx_zero_B = df[name_column].str.contains(Zero_B).fillna(value=False)

    df['Channel']  = 'A'
    df['Channel'].iloc[idx_ref_B]  = 'B'
    df['Channel'].iloc[idx_cand_B] = 'B'
    df['Channel'].iloc[idx_zero_B] = 'B'
    df.Channel = df.Channel.astype('category')

    df['Zero'] = False
    df['Zero'].iloc[idx_zero_A] = True
    df['Zero'].iloc[idx_zero_B] = True

    df['Generator'] = 'Cand'
    df['Generator'].iloc[idx_ref_A] = 'Ref'
    df['Generator'].iloc[idx_ref_B] = 'Ref'
    df.Generator = df.Generator.astype('category')

    df['Bracket'] = 0
    count = 0
    for idx in range(1, len(df.index)):
        if df[value_column].iloc[idx-1]>(10*df[value_column].iloc[idx]):
            # this assumes that the setpoint being measured is larger than 10*(zero response)
            # i.e., we are measuring reasonably far away from the zero point
            count += 1
        df['Bracket'].iloc[idx] = count
    df.Bracket = df.Bracket.astype('category')

    # zero measurements subseries
    zero_A_measurements = df[idx_zero_A].copy()
    zero_A_measurements = zero_A_measurements.reset_index(drop=True)
    zero_B_measurements = df[idx_zero_B].copy()
    zero_B_measurements = zero_B_measurements.reset_index(drop=True)

    # ref measurements subseries
    ref_A_measurements = df[idx_ref_A].copy()
    ref_A_measurements = ref_A_measurements.reset_index(drop=True)
    ref_B_measurements = df[idx_ref_B].copy()
    ref_B_measurements = ref_B_measurements.reset_index(drop=True)

    # cand measurements subseries
    cand_A_measurements = df[idx_cand_A].copy()
    cand_A_measurements = cand_A_measurements.reset_index(drop=True)
    cand_B_measurements = df[idx_cand_B].copy()
    cand_B_measurements = cand_B_measurements.reset_index(drop=True)

    # ----------------------------------------------------------------------------
    """Plot data"""

    # add one sheet to display the plots of the data
    try:
        logging.debug('Add sheet to the Excel workbook')
        wb.sheets.add('Plot data')
    except ValueError as V:
        logging.debug('Sheet already present in workbook')
        print("Error:", V)

    # zero concentration measurements
    fig = plt.figure('zero measurements')
    ax = fig.add_subplot(
        111, 
        xlabel='Time [min]', 
        #ylabel=r'Peak area [mm$^2$]', 
        title='Zero concentration')
    sns.scatterplot(
        data=df[df['Zero']==True], 
        x='Elapsed_minutes', 
        y=value_column, 
        hue='Channel')
    plt.tight_layout()

    # add figure to Excel sheet
    wb.sheets['Plot data'].pictures.add(
        fig, 
        name='Zero concentration', 
        update=True, 
        left=wb.sheets['Plot data'].range('A1').left, 
        top=wb.sheets['Plot data'].range('A1').top)


    # non-zero concentration measurements
    fig = plt.figure('measurements')
    ax = fig.add_subplot(
        111, 
        xlabel='Time [min]', 
        #ylabel=r'Peak area [mm$^2$]', 
        title='Mercury concentration')
    sns.scatterplot(
        data=df[df['Zero']==False], 
        x='Elapsed_minutes', 
        y=value_column, 
        hue='Channel', 
        style='Generator')
    #ax.set_ylim([6000, 8000])
    ax.legend(loc='center left', ncol=1, bbox_to_anchor=(1., .5, 0.5, 0.5))
    plt.tight_layout()

    # add figure to Excel sheet
    wb.sheets['Plot data'].pictures.add(
        fig, 
        name='Non-zero concentration', 
        update=True, 
        left=wb.sheets['Plot data'].range('N1').left, 
        top=wb.sheets['Plot data'].range('N1').top)

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    """Zero correction"""
    if apply_zero_correction:
        # get mean and std of zero readings per channel per bracket 
        # channel A
        mean_zero_A_meas = zero_A_measurements.groupby('Bracket').median()[value_column]
        std_sero_A_meas = zero_A_measurements.groupby('Bracket').std()[value_column]
        time_mean_zero_A_meas = zero_A_measurements.groupby('Bracket').mean()['Elapsed_seconds']

        # channel B
        mean_zero_B_meas = zero_B_measurements.groupby('Bracket').median()[value_column]
        std_sero_B_meas = zero_B_measurements.groupby('Bracket').std()[value_column]
        time_mean_zero_B_meas = zero_B_measurements.groupby('Bracket').mean()['Elapsed_seconds']

        # interpolate the mean zero readings on the mean time
        # channel A
        zero_correction_A = interp1d(time_mean_zero_A_meas, mean_zero_A_meas, kind='linear')
        zero_values_ref_A = zero_correction_A(ref_A_measurements.Elapsed_seconds)
        zero_values_cand_A = zero_correction_A(cand_A_measurements.Elapsed_seconds)

        # channel B
        zero_correction_B = interp1d(time_mean_zero_B_meas, mean_zero_B_meas, kind='linear')
        zero_values_ref_B = zero_correction_B(ref_B_measurements.Elapsed_seconds)
        zero_values_cand_B = zero_correction_B(cand_B_measurements.Elapsed_seconds)

        # compute zero corrected value_column
        value_column_zero_cor = value_column + ' Zero Cor'
        # channel A
        cand_A_measurements[value_column_zero_cor] = cand_A_measurements[value_column].values - zero_values_cand_A
        ref_A_measurements[value_column_zero_cor] = ref_A_measurements[value_column].values - zero_values_ref_A

        # channel B
        cand_B_measurements[value_column_zero_cor] = cand_B_measurements[value_column].values - zero_values_cand_B
        ref_B_measurements[value_column_zero_cor] = ref_B_measurements[value_column].values - zero_values_ref_B

        # drop mislabelled non-zero measurements that have negative peak area after correction
        cand_A_measurements = cand_A_measurements.drop(cand_A_measurements[cand_A_measurements[value_column_zero_cor]<0].index)
        cand_A_measurements = cand_A_measurements.reset_index(drop=True)

        cand_B_measurements = cand_B_measurements.drop(cand_B_measurements[cand_B_measurements[value_column_zero_cor]<0].index)
        cand_B_measurements = cand_B_measurements.reset_index(drop=True)

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    """Compute output ratio"""

    cand_A_ratio = []
    cand_B_ratio = []

    for bracket in df.Bracket.cat.categories:
        # channel A
        if apply_zero_correction:
            cand_A_meas = cand_A_measurements[value_column_zero_cor][cand_A_measurements.Bracket == bracket].values
        else:
            cand_A_meas = cand_A_measurements[value_column][cand_A_measurements.Bracket == bracket].values
        
        if (len(cand_A_meas)>0) :
            cand_A_time = cand_A_measurements['Elapsed_seconds'][cand_A_measurements.Bracket == bracket].values
            if apply_zero_correction:
                ref_A_meas = ref_A_measurements[value_column_zero_cor][ref_A_measurements.Bracket == bracket].values
            else:
                ref_A_meas = ref_A_measurements[value_column][ref_A_measurements.Bracket == bracket].values
            #
            time_ref_A = ref_A_measurements['Elapsed_seconds'][ref_A_measurements.Bracket == bracket].values
            
            for idx in range(len(cand_A_meas)):
                try:
                    cand_A_out_ratio = output_ratio(cand_A_meas[idx], ref_A_meas[idx+1], ref_A_meas[idx], 
                                                    cand_A_time[idx], time_ref_A[idx+1], time_ref_A[idx])
                    cand_A_ratio.append(cand_A_out_ratio)
                except IndexError:
                    if not (len(ref_A_meas)>len(cand_A_meas)):
                        print('Channel A: candidate generator has more measurements than reference generator.')
                        cand_A_ratio.append(np.nan)

        # channel B
        if apply_zero_correction:
            cand_B_meas = cand_B_measurements[value_column_zero_cor][cand_B_measurements.Bracket == bracket].values
        else:
            cand_B_meas = cand_B_measurements[value_column][cand_B_measurements.Bracket == bracket].values
        
        if (len(cand_B_meas)>0) :
            cand_B_time = cand_B_measurements['Elapsed_seconds'][cand_B_measurements.Bracket == bracket].values
            if apply_zero_correction:
                ref_B_meas = ref_B_measurements[value_column_zero_cor][ref_B_measurements.Bracket == bracket].values
            else:
                ref_B_meas = ref_B_measurements[value_column][ref_B_measurements.Bracket == bracket].values
            #
            time_ref_B = ref_B_measurements['Elapsed_seconds'][ref_B_measurements.Bracket == bracket].values
            
            for idx in range(len(cand_B_meas)):
                try:
                    cand_B_out_ratio = output_ratio(cand_B_meas[idx], ref_B_meas[idx+1], ref_B_meas[idx], 
                                                    cand_B_time[idx], time_ref_B[idx+1], time_ref_B[idx])
                    cand_B_ratio.append(cand_B_out_ratio)
                except IndexError:
                    if not (len(ref_B_meas)>len(cand_B_meas)):
                        print('Channel B: candidate generator has more measurements than reference generator.')
                        cand_B_ratio.append(np.nan)

    cand_A_measurements['Output_ratio'] = cand_A_ratio
    cand_B_measurements['Output_ratio'] = cand_B_ratio

    # average ratio
    # channel A
    avg_ratio_A = cand_A_measurements.groupby('Bracket').mean()['Output_ratio']
    # channel B
    avg_ratio_B = cand_B_measurements.groupby('Bracket').mean()['Output_ratio']

    mean_ratio_A = np.nanmean(avg_ratio_A)
    mean_ratio_B = np.nanmean(avg_ratio_B)

    ratio_RSD_A = {}
    ratio_RSD_B = {}
    # ratio RSD
    for bracket in df.Bracket.cat.categories:
        # channel A
        out_ratio_A = cand_A_measurements.Output_ratio[cand_A_measurements.Bracket == bracket].values
        # remove nan
        out_ratio_A = out_ratio_A[~np.isnan(out_ratio_A)]
        if (len(out_ratio_A)>0):
            ratio_RSD_A[bracket] = (np.sqrt(np.sum((out_ratio_A-avg_ratio_A[bracket])**2)/
                                            (len(out_ratio_A)-1)) / avg_ratio_A[bracket])
            print('Ratio RSD for channel A at bracket {}: {:2.2f}%'.format(bracket,ratio_RSD_A[bracket]*100))
        # channel B
        out_ratio_B = cand_B_measurements.Output_ratio[cand_B_measurements.Bracket == bracket].values
        # remove nan
        out_ratio_B = out_ratio_B[~np.isnan(out_ratio_B)]
        if (len(out_ratio_B)>0) :
            ratio_RSD_B[bracket] = (np.sqrt(np.sum((out_ratio_B-avg_ratio_B[bracket])**2)/
                                            (len(out_ratio_B)-1)) / avg_ratio_B[bracket])
            print('Ratio RSD for channel B at bracket {}: {:2.2f}%'.format(bracket,ratio_RSD_B[bracket]*100))

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    """Compute measurement stability uncertainty"""
    SE_stab_cand_A = {}
    SE_stab_ref_A = {}

    SE_stab_cand_B = {}
    SE_stab_ref_B = {}

    for channel in df.Channel.cat.categories:
        for bracket in df.Bracket.cat.categories:
            for generator in df.Generator.cat.categories:
                meas = df[value_column][(df.Channel == channel) & (df.Zero == False) &
                                        (df.Bracket == bracket) & (df.Generator == generator)].values
                if (len(meas)>2): # need minimum 3 measurements
                    time = df['Elapsed_seconds'][(df.Channel == channel) & (df.Zero == False) &
                                                 (df.Bracket == bracket) & (df.Generator == generator)].values
                    # compute averages
                    meas_avg = np.mean(meas)
                    time_avg = np.mean(time)
                
                    # compute sum of variances
                    meas_var = np.sum((meas-meas_avg)**2)
                    time_var = np.sum((time-time_avg)**2)
                
                    # compute covariance squared
                    meas_time_cov = np.sum((meas-meas_avg)*(time-time_avg))
                
                    if channel == 'A':
                        if generator == 'Cand':
                            SE_stab_cand_A[bracket] = SE_stability(meas_var, time_var, meas_time_cov, len(meas))
                        else: # generator == 'Ref'
                            SE_stab_ref_A[bracket] = SE_stability(meas_var, time_var, meas_time_cov, len(meas))
                    else: # channel == 'B'
                        if generator == 'Cand':
                            SE_stab_cand_B[bracket] = SE_stability(meas_var, time_var, meas_time_cov, len(meas))
                        else: # generator == 'Ref'
                            SE_stab_ref_B[bracket] = SE_stability(meas_var, time_var, meas_time_cov, len(meas))

    # ----------------------------------------------------------------------------
    u_stability_A = {}
    u_stability_B = {}

    for bracket in cand_A_measurements.Bracket.cat.categories:
        # channel A    
        ratio_A = cand_A_measurements['Output_ratio'][cand_A_measurements.Bracket == bracket].values
        # remove nan
        ratio_A = ratio_A[~np.isnan(ratio_A)]
    
        if (len(ratio_A)>0):
            u_A = []
            if apply_zero_correction:
                cand_A_meas = cand_A_measurements[value_column_zero_cor][cand_A_measurements.Bracket == bracket].values
                ref_A_meas = ref_A_measurements[value_column_zero_cor][ref_A_measurements.Bracket == bracket].values
            else:
                cand_A_meas = cand_A_measurements[value_column][cand_A_measurements.Bracket == bracket].values
                ref_A_meas = ref_A_measurements[value_column][ref_A_measurements.Bracket == bracket].values
        
            cand_A_time = cand_A_measurements['Elapsed_seconds'][cand_A_measurements.Bracket == bracket].values
            ref_A_time = ref_A_measurements['Elapsed_seconds'][ref_A_measurements.Bracket == bracket].values
    
            for idx in range(len(ratio_A)):
                try:
                    u_A.append(u_stability(ratio_A[idx], SE_stab_cand_A[bracket], SE_stab_ref_A[bracket], 
                                           SE_stab_ref_A[bracket], cand_A_meas[idx], ref_A_meas[idx+1], 
                                           ref_A_meas[idx], cand_A_time[idx], ref_A_time[idx+1], 
                                           ref_A_time[idx]))
                except IndexError:
                    if not (len(ref_A_meas)>len(cand_A_meas)):
                        print('Channel A: candidate generator has more measurements than reference generator.')
            u_A = np.array(u_A)
            u_stability_A[bracket] = np.sqrt(np.sum(u_A**2)/(len(u_A)**2))

        # channel B    
        ratio_B = cand_B_measurements['Output_ratio'][cand_B_measurements.Bracket == bracket].values
        # remove nan
        ratio_B = ratio_B[~np.isnan(ratio_B)]
    
        if (len(ratio_B)>0):
            u_B = []
            if apply_zero_correction:
                cand_B_meas = cand_B_measurements[value_column_zero_cor][cand_B_measurements.Bracket == bracket].values
                ref_B_meas = ref_B_measurements[value_column_zero_cor][ref_B_measurements.Bracket == bracket].values
            else:
                cand_B_meas = cand_B_measurements[value_column][cand_B_measurements.Bracket == bracket].values
                ref_B_meas = ref_B_measurements[value_column][ref_B_measurements.Bracket == bracket].values
        
            cand_B_time = cand_B_measurements['Elapsed_seconds'][cand_B_measurements.Bracket == bracket].values
            ref_B_time = ref_B_measurements['Elapsed_seconds'][ref_B_measurements.Bracket == bracket].values
    
            for idx in range(len(ratio_B)):
                try:
                    u_B.append(u_stability(ratio_B[idx], SE_stab_cand_B[bracket], SE_stab_ref_B[bracket], 
                                           SE_stab_ref_B[bracket], cand_B_meas[idx], ref_B_meas[idx+1], 
                                           ref_B_meas[idx], cand_B_time[idx], ref_B_time[idx+1], 
                                           ref_B_time[idx]))
                except IndexError:
                    if not (len(ref_B_meas)>len(cand_B_meas)):
                        print('Channel B: candidate generator has more measurements than reference generator.')
            u_B = np.array(u_B)
            u_stability_B[bracket] = np.sqrt(np.sum(u_B**2)/(len(u_B)**2))

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    """Compute measurement repeatability uncertainty"""

    u_repeatability_A = {}
    u_repeatability_B = {}

    L = 3 # in this set of measurements there is 1 cand measurement bracketed between 2 ref measurements

    for bracket in df.Bracket.cat.categories:
        # channel A    
        ratio_A = cand_A_measurements['Output_ratio'][cand_A_measurements.Bracket == bracket].values
        # remove possible nan values
        ratio_A = ratio_A[~np.isnan(ratio_A)]
        K_A = len(ratio_A)
        if (K_A>0):
            u_repeatability_A[bracket] = u_repeatability(K_A, L, ratio_A, u_stability_A[bracket])
    
        # channel B
        ratio_B = cand_B_measurements['Output_ratio'][cand_B_measurements.Bracket == bracket].values
        # remove possible nan values
        ratio_B = ratio_B[~np.isnan(ratio_B)]
        K_B = len(ratio_B)
        if (K_B>0):
            u_repeatability_B[bracket] = u_repeatability(K_B, L, ratio_B, u_stability_B[bracket])

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    """
    Combine repeatability and stability uncertainty
    Compute reproducibility uncertainty contribution
    Get combined uncertainty
    """

    u_stab_repeat_A = {}
    u_stab_repeat_B = {}

    for bracket in u_stability_A.keys():
        # channel A
        u_stab_repeat_A[bracket] = np.sqrt(u_stability_A[bracket]**2 + u_repeatability_A[bracket]**2)

        # channel B
        u_stab_repeat_B[bracket] = np.sqrt(u_stability_B[bracket]**2 + u_repeatability_B[bracket]**2)

    u_comparison_A = np.sqrt( np.sum((ref_concentration*np.array(list(u_stab_repeat_A.values())))**2)/
                             (len(u_stab_repeat_A.values())**2) )
    u_comparison_B = np.sqrt( np.sum((ref_concentration*np.array(list(u_stab_repeat_B.values())))**2)/
                             (len(u_stab_repeat_B.values())**2) )

    u_reprod_A = ref_concentration*(np.nanmax(avg_ratio_A)-np.nanmin(avg_ratio_A))/np.sqrt(12)
    if not u_reprod_A > 1e-16:
        # single bracket
        u_reprod_A = reproducibility_uncertainty*ref_concentration*mean_ratio_A

    u_reprod_B = ref_concentration*(np.nanmax(avg_ratio_B)-np.nanmin(avg_ratio_B))/np.sqrt(12)
    if not u_reprod_B > 1e-16:
        # single bracket
        u_reprod_B = reproducibility_uncertainty*ref_concentration*mean_ratio_B

    u_reference_A = u_ref_concentration*mean_ratio_A
    u_reference_B = u_ref_concentration*mean_ratio_B

    u_combined_A = np.sqrt(u_comparison_A**2 + u_reprod_A**2 + u_reference_A**2)
    u_combined_B = np.sqrt(u_comparison_B**2 + u_reprod_B**2 + u_reference_B**2)

    print('Relative expanded certification uncertainty (k=2) for channel A: {:2.2f}%'
          .format(100*(2*u_combined_A)/(ref_concentration*mean_ratio_A)))
    print('Relative expanded certification uncertainty (k=2) for channel B: {:2.2f}%'
          .format(100*(2*u_combined_B)/(ref_concentration*mean_ratio_B)))

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    """Write main results to file"""

    # channel A
    # create dataframe where to store uncertainty components
    col_names = ['Symbol', 'Unit']
    for bracket in u_repeatability_A.keys():
        col_names.append('Bracket_{}'.format(bracket+1))
    col_names.append('All_brackets')

    index_names = [
        'Mean output ratio', 
        'RSD output ratio', 
        'Stability uncertainty', 
        'Repeatability uncertainty', 
        'Comparison uncertainty', 
        'Reproducibility uncertainty', 
        'Reference uncertainty', 
        'Combined std uncertainty', 
        'Expanded uncertainty (k=2)', 
        'Relative expanded uncertainty',
        'Candidate concentration',
        'Relative difference in concentration']

    results_A = pd.DataFrame(columns=col_names, index=index_names)

    results_A['Symbol'] = [
        'R', 
        '', 
        'u_stability', 
        'u_repeatability',
        'u_comparison',
        'u_reproducibility',
        'u_reference',
        'u(c)',
        'U(c)',
        'U(c)',
        'c',
        'D_i']

    results_A['Unit'] = [
        '',
        '%',
        'ng/m3',
        'ng/m3',
        'ng/m3',
        'ng/m3',
        'ng/m3',
        'ng/m3',
        'ng/m3',
        '%',
        'ng/m3',
        '%']

    for bracket in u_repeatability_A.keys():
        column_name = 'Bracket_{}'.format(bracket+1)
        results_A.loc['Mean output ratio', column_name] = avg_ratio_A[bracket]
        results_A.loc['RSD output ratio', column_name] = 100.0*ratio_RSD_A[bracket]
        results_A.loc['Stability uncertainty', column_name] = u_stability_A[bracket]
        results_A.loc['Repeatability uncertainty', column_name] = u_repeatability_A[bracket]

    results_A.loc['Mean output ratio', 'All_brackets'] = mean_ratio_A

    results_A.loc['Comparison uncertainty', 'All_brackets'] = u_comparison_A

    results_A.loc['Reproducibility uncertainty', 'All_brackets'] = u_reprod_A

    results_A.loc['Reference uncertainty', 'All_brackets'] = u_reference_A

    results_A.loc['Combined std uncertainty', 'All_brackets'] = u_combined_A

    results_A.loc['Expanded uncertainty (k=2)', 'All_brackets'] = u_combined_A*2

    results_A.loc['Relative expanded uncertainty', 'All_brackets'] = 100.0*(u_combined_A*2)/(ref_concentration*mean_ratio_A)

    results_A.loc['Candidate concentration', 'All_brackets'] = ref_concentration*mean_ratio_A

    results_A.loc['Relative difference in concentration', 'All_brackets'] = 100.0*(ref_concentration*mean_ratio_A/cand_concentration - 1.0)

    # ----------------------------------------------------------------------------
    # channel B
    # create dataframe where to store uncertainty components
    results_B = pd.DataFrame(columns=col_names, index=index_names)

    results_B['Symbol'] = [
        'R', 
        '', 
        'u_stability', 
        'u_repeatability',
        'u_comparison',
        'u_reproducibility',
        'u_reference',
        'u(c)',
        'U(c)',
        'U(c)',
        'c',
        'D_i']

    results_B['Unit'] = [
        '',
        '%',
        'ng/m3',
        'ng/m3',
        'ng/m3',
        'ng/m3',
        'ng/m3',
        'ng/m3',
        'ng/m3',
        '%',
        'ng/m3',
        '%']

    for bracket in u_repeatability_B.keys():
        column_name = 'Bracket_{}'.format(bracket+1)
        results_B.loc['Mean output ratio', column_name] = avg_ratio_B[bracket]
        results_B.loc['RSD output ratio', column_name] = 100.0*ratio_RSD_B[bracket]
        results_B.loc['Stability uncertainty', column_name] = u_stability_B[bracket]
        results_B.loc['Repeatability uncertainty', column_name] = u_repeatability_B[bracket]

    results_B.loc['Mean output ratio', 'All_brackets'] = mean_ratio_B

    results_B.loc['Comparison uncertainty', 'All_brackets'] = u_comparison_B

    results_B.loc['Reproducibility uncertainty', 'All_brackets'] = u_reprod_B

    results_B.loc['Reference uncertainty', 'All_brackets'] = u_reference_B

    results_B.loc['Combined std uncertainty', 'All_brackets'] = u_combined_B

    results_B.loc['Expanded uncertainty (k=2)', 'All_brackets'] = u_combined_B*2

    results_B.loc['Relative expanded uncertainty', 'All_brackets'] = 100.0*(u_combined_B*2)/(ref_concentration*mean_ratio_B)

    results_B.loc['Candidate concentration', 'All_brackets'] = ref_concentration*mean_ratio_B

    results_B.loc['Relative difference in concentration', 'All_brackets'] = 100.0*(ref_concentration*mean_ratio_B/cand_concentration - 1.0)

    # ----------------------------------------------------------------------------

    # write to Excel
    # add one sheet to display the plots of the data
    try:
        logging.debug('Add sheet to the Excel workbook')
        wb.sheets.add('Channel A')
        wb.sheets.add('Channel B')
    except ValueError as V:
        logging.debug('Sheet already present in workbook')
        print("Error:", V)

    # set bold font on first row
    wb.sheets['Channel A'].range('A1:AAA1').api.Font.Bold = True
    wb.sheets['Channel B'].range('A1:AAA1').api.Font.Bold = True

    # set bold font on first column
    wb.sheets['Channel A'].range('A1:A20').api.Font.Bold = True
    wb.sheets['Channel B'].range('A1:A20').api.Font.Bold = True

    # write dataset to Excel sheet
    wb.sheets['Channel A'].range('A1').options(index=True).value = results_A
    wb.sheets['Channel B'].range('A1').options(index=True).value = results_B
