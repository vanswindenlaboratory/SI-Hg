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
from scipy.special import gammainc
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
plt.rcParams['figure.figsize'] = 9,5
import seaborn as sns

from scipy.interpolate import interp1d
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std

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
        output ratio (eq. 3b in SI-Hg protocol)
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

def regress_function(deg, xdata, ydata, sigma, include_inverse=False):
    """
    Input:
        - deg : integer; degree of polynomial to be fitted
        - xdata : array or list, shape(M,); independent variable of the M sample points
        - ydata : array or list, shape(M,); dependent variable of the M sample points
        - sigma : array or list, shape(M,); weights to apply to the dependent variable of the M sample points
        - include_inverse : boolean; if True, the term 1/xdata is added to the features (default value is False)
    
    Output:
        - params : array; values of the function parameters
        - cov : array; variance and covariance of function parameter values
    """
    # check if input variables are in the correct format
    if not isinstance(deg, int):
        print('Received variable deg of type {}. Acceptable formats are: integer.'.format(type(deg)))
        # terminate execution
        exit()        
    
    input_vars = [xdata, ydata, sigma]
    for var in input_vars: 
        if not isinstance(var, (list, np.ndarray)):
            print('Received variable {} of type {}. Acceptable formats are: list and numpy.ndarray.'.format(var,type(var)))
            # terminate execution
            exit()

    # if type(sigma) is list, convert it into np.ndarray
    if isinstance(sigma, list):
        sigma = np.array(sigma)
    
    if deg == 0:
        features = np.ones(len(ydata)) # only intercept
    elif deg > 0:
        features = sm.add_constant(xdata) # intercept and slope
        if deg > 1:
            for i in range(2, deg+1):
                # higher degree terms
                features = np.column_stack((features, np.power(xdata,i)))
    else:
         # deg < 0
         print('Received deg = {}. Requirement: deg>=0'.format(deg))
         # terminate execution
         exit()
    
    if include_inverse:
        features = np.column_stack((features, np.power(xdata,-1)))
    
    return sm.WLS(ydata, features, np.power(sigma,-2)).fit()

# ----------------------------------------------------------------------------

def predict(deg, wls_res, xdata, include_inverse=False):
    """
    Input:
        - deg : integer; degree of polynomial to be fitted
        - wls_res : statsmodels WLS object; regressed function obtained using statsmodels library 
        - xdata : array or list, shape(M,); independent variable of the M sample points
        - include_inverse : boolean; if True, the term 1/xdata is added to the features (default value is False)
    
    Output:
        - pred : array, shape(M,); prediction for the M sample points
    """
    if not isinstance(deg, int):
        print('Received variable deg of type {}. Acceptable formats are: integer.'.format(type(deg)))
        # terminate execution
        exit()        

    if not isinstance(xdata, (list, np.ndarray)):
        print('Received variable xdata of type {}. Acceptable formats are: list and numpy.ndarray.'.format(type(xdata)))
        # terminate execution
        exit()

    if not isinstance(wls_res, sm.regression.linear_model.RegressionResultsWrapper):
        print('Input regression is not a statsmodels object')
        exit()

    if isinstance(xdata, list):
        xdata = np.array(xdata)

    if deg == 0:
        features = np.ones(len(xdata)) # only intercept
    elif deg > 0:
        features = sm.add_constant(xdata) # intercept and slope
        if deg > 1:
            for i in range(2, deg+1):
                # higher degree terms
                features = np.column_stack((features, np.power(xdata,i)))
    else:
         # deg < 0
         print('Received deg = {}. Requirement: deg>=0'.format(deg))
         # terminate execution
         exit()
    
    if include_inverse:
        features = np.column_stack((features, np.power(xdata,-1)))
    
    return wls_res.predict(features)

# ----------------------------------------------------------------------------

def chi2_WLS(y, x, params, sigma_y):
    """
    Input:
        - y : array, list or float, shape(M,); target variable of the M sample points
        - x : array, list or float, shape(M,); independent variable of the M sample points
        - params : array, shape(n_params,); parameters of the polynomial y = p(x)
        - sigma_y : array, list or float, shape(M,); uncertainties associated to the M sample points
    
    Output:
        - chi2 : float; chi square value 
    """
    y_pred = 0
    for deg in range(len(params)):
        y_pred += params[deg]*np.power(x, deg)

    return sum(np.power((y-y_pred)/sigma_y, 2))

# ----------------------------------------------------------------------------

def plot_poly(x, x_pred, wls_regression, wls_pred, axis, color_line, label_name, conf_lvl=0.95, transparancy=1):    
    # get WLS 95% CI
    wls_pred_err, wls_lower, wls_upper = wls_prediction_std(wls_regression, alpha=1-conf_lvl)
        
    # plot regression line
    axis.plot(x_pred, wls_pred, lw=2, color=color_line, alpha=transparancy, label='{}'.format(label_name))
    axis.plot(x, wls_lower, lw=0.5, color=color_line, alpha=transparancy, ls='--')
    axis.plot(x, wls_upper, lw=0.5, color=color_line, alpha=transparancy, ls='--')

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

    ref_setpoints = wb.sheets['Input_info'].range('D2').options(ndim=1, expand='horizontal').value
    cand_setpoints = wb.sheets['Input_info'].range('D3').options(ndim=1, expand='horizontal').value

    if len(ref_setpoints)==len(cand_setpoints):
        wb.sheets['Input_info'].range('L11').delete()
        ref_concentration = {cand_setpoints[i]:ref_setpoints[i] for i in range(len(cand_setpoints))}
    else:
        wb.sheets['Input_info'].range('L11').value = (
            'The number of setpoints given for the reference generator is not the same as the number of setpoints given for the candidate. Please fix.')
        exit()

    extended_uncertainty_ref_generator = wb.sheets['Input_info'].range('D4').value/100.
    u_ref_concentration = {}
    for setpoint in ref_concentration.keys():
        u_ref_concentration[setpoint] = extended_uncertainty_ref_generator*ref_concentration[setpoint]/2.0 # in this example, U_ref = 3% (k=2)

    data_sheet = 'Data'

    datetime_column = 'DateTime' 
    name_column = 'Name' 
    valid_column = 'Valid' 
    value_column = 'Pk Area' 
    setpoint_column = 'Setpoint_candidate' 

    Ref_A = wb.sheets['Input_info'].range('D5').value
    Ref_B = wb.sheets['Input_info'].range('D6').value
    Cand_A = wb.sheets['Input_info'].range('D7').value
    Cand_B = wb.sheets['Input_info'].range('D8').value
    Zero_A = wb.sheets['Input_info'].range('D9').value
    Zero_B = wb.sheets['Input_info'].range('D10').value

    apply_zero_correction = bool(wb.sheets['Input_info'].range('D11').value)
    reproducibility_uncertainty = wb.sheets['Input_info'].range('D12').value/100.0

    max_mins_between_meas = int(wb.sheets['Input_info'].range('D13').value)
    datetime_format = wb.sheets['Input_info'].range('D14').value

    # ----------------------------------------------------------------------------
    """Load and pre-process data."""
    df = wb.sheets[data_sheet].range('A1').options(pd.DataFrame, header=1, index=False, expand='table').value

    # convert valid_column to boolean 
    df[valid_column] = df[valid_column].astype('bool')

    # if present, drop invalid measurements
    df = df.drop(df[df[valid_column] == False].index)

    # drop completely empty lines
    df = df.dropna(how='all')
    # reset index to avoid holes in the index
    df = df.reset_index(drop=True)
    
    # convert setpoint_column from float to category 
    df[setpoint_column] = df[setpoint_column].astype('category')

    # convert Date - Time from object to datetime
    if isinstance(df[datetime_column], object):
        df[datetime_column] = pd.to_datetime(df[datetime_column], format=datetime_format)
        print(df.info())

    # if present, drop nan values in value_column
    df = df.dropna(subset=[value_column])
    df = df.reset_index(drop=True)

    # ----------------------------------------------------------------------------
    df['elapsed_mins'] = 0.0
    # convert datetime to number of seconds elapsed since the first recorded measurement 
    time_elapsed = df[datetime_column] - df[datetime_column].iloc[0]
    time_elapsed = time_elapsed.dt.total_seconds().values/60.0
    df['elapsed_mins'] = time_elapsed

    # ----------------------------------------------------------------------------
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

    df['Date'] = df[datetime_column].dt.date

    # ----------------------------------------------------------------------------
    # create sub-datasets for ref, cand and zero measurements

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
    fig = plt.figure('zero measurements', figsize=[13,5])
    ax = fig.add_subplot(
        111, 
        xlabel='Time [min]', 
        ylabel=r'Peak area [mm$^2$]', 
        title='Zero concentration')
    sns.scatterplot(
        data=df[df['Zero']==True], 
        x='elapsed_mins', 
        y=value_column, 
        hue='Channel', 
        size='Date')
    ax.legend(
        loc='center left', 
        ncol=1, 
        bbox_to_anchor=(1., .5, 0.5, 0.5))
    plt.tight_layout()

    # add figure to Excel sheet
    wb.sheets['Plot data'].pictures.add(
        fig, 
        name='Zero concentration', 
        update=True, 
        left=wb.sheets['Plot data'].range('A1').left, 
        top=wb.sheets['Plot data'].range('A1').top)

    # non-zero concentration measurements
    fig = plt.figure('measurements', figsize=[13,7])
    ax = fig.add_subplot(
        111, 
        xlabel='Time [min]', 
        ylabel=r'Peak area [mm$^2$]', 
        title='Mercury concentration')
    sns.scatterplot(
        data=df[df['Zero']==False], 
        x='elapsed_mins', 
        y=value_column, 
        hue='Channel', 
        style='Generator', 
        size='Date')
    ax.legend(
        loc='center left', 
        ncol=1, 
        bbox_to_anchor=(1., .5, 0.5, 0.5))
    plt.tight_layout()

    # add figure to Excel sheet
    wb.sheets['Plot data'].pictures.add(
        fig, 
        name='Mercury concentration', 
        update=True, 
        left=wb.sheets['Plot data'].range('T1').left, 
        top=wb.sheets['Plot data'].range('T1').top)

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------

    """Perform single point analysis for each point in setpoint_column: apply zero correction"""

    if apply_zero_correction:
        # initialize column that will contain the zero-corrected peak area values
        value_column_zero_cor = value_column + ' Zero Cor'
        df[value_column_zero_cor] = 0.0

        # channel A

        # divide zero measurements in groups
        group = 0
        zero_group = [group]
        for idx in range(len(zero_A_measurements.index)-1):
            if (zero_A_measurements.elapsed_mins.iloc[idx+1]-zero_A_measurements.elapsed_mins.iloc[idx])>max_mins_between_meas:
                group += 1
            zero_group.append(group)
        zero_A_measurements['Zero_group'] = zero_group
    
        # get the median peak area for each zero group
        medians_zero_A = zero_A_measurements.groupby('Zero_group').median()[value_column]
        # get the time corresponding to the medians
        time_medians_zero_A = zero_A_measurements.groupby('Zero_group').median()['elapsed_mins']
        # interpolate between the values of the median
        zero_correction_A = interp1d(time_medians_zero_A, medians_zero_A, kind='linear')
    
        # correct measurements of reference analyzer
        zero_values_ref_A = zero_correction_A(ref_A_measurements.elapsed_mins)
        ref_A_measurements[value_column_zero_cor] = ref_A_measurements[value_column]-zero_values_ref_A

        # correct measurements of candidate analyzer
        zero_values_cand_A = zero_correction_A(cand_A_measurements.elapsed_mins)
        cand_A_measurements[value_column_zero_cor] = cand_A_measurements[value_column]-zero_values_cand_A

        # ----------------------------------------------------------------------------
        # channel B

        # divide zero measurements in groups
        group = 0
        zero_group = [0]
        for idx in range(len(zero_B_measurements.index)-1):
            if (zero_B_measurements.elapsed_mins.iloc[idx+1]-zero_B_measurements.elapsed_mins.iloc[idx])>max_mins_between_meas:
                group += 1
            zero_group.append(group)
        zero_B_measurements['Zero_group'] = zero_group
    
        # get the median peak area for each zero group
        medians_zero_B = zero_B_measurements.groupby('Zero_group').median()[value_column]
        # get the time corresponding to the medians
        time_medians_zero_B = zero_B_measurements.groupby('Zero_group').median()['elapsed_mins']
        # interpolate between the values of the median
        zero_correction_B = interp1d(time_medians_zero_B, medians_zero_B, kind='linear')
    
        # correct measurements of reference analyzer
        zero_values_ref_B = zero_correction_B(ref_B_measurements.elapsed_mins)
        ref_B_measurements[value_column_zero_cor] = ref_B_measurements[value_column]-zero_values_ref_B

        # correct measurements of candidate analyzer
        zero_values_cand_B = zero_correction_B(cand_B_measurements.elapsed_mins)
        cand_B_measurements[value_column_zero_cor] = cand_B_measurements[value_column]-zero_values_cand_B

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------

    """Perform single point analysis for each point in setpoint_column: compute output ratio"""

    cand_A_measurements['output_ratio'] = 0.0
    cand_B_measurements['output_ratio'] = 0.0

    cand_A_ratio = []
    cand_B_ratio = []

    for setpoint in df[setpoint_column].cat.categories:
        if setpoint == 0.0:
            continue
        else:
            # ----------------------------------------------------------------------------
            # channel A
            if apply_zero_correction:
                cand_A_meas = cand_A_measurements[value_column_zero_cor][cand_A_measurements[setpoint_column] == setpoint].values
            else:
                cand_A_meas = cand_A_measurements[value_column][cand_A_measurements[setpoint_column] == setpoint].values
        
            if (len(cand_A_meas)>0) :
                cand_A_time = cand_A_measurements.elapsed_mins[cand_A_measurements[setpoint_column] == setpoint].values
                if apply_zero_correction:
                    ref_A_meas = ref_A_measurements[value_column_zero_cor][ref_A_measurements[setpoint_column] == setpoint].values
                else:
                    ref_A_meas = ref_A_measurements[value_column][ref_A_measurements[setpoint_column] == setpoint].values
                #
                time_ref_A = ref_A_measurements.elapsed_mins[ref_A_measurements[setpoint_column] == setpoint].values
                for idx in range(len(cand_A_meas)):
                    try:
                        cand_A_out_ratio = output_ratio(
                            cand_A_meas[idx], 
                            ref_A_meas[idx+1], 
                            ref_A_meas[idx], 
                            cand_A_time[idx], 
                            time_ref_A[idx+1], 
                            time_ref_A[idx])
                        cand_A_ratio.append(cand_A_out_ratio)
                    except IndexError:
                        if not (len(ref_A_meas)>len(cand_A_meas)):
                            print('Channel A: candidate generator has more measurements than reference generator.')
                            cand_A_ratio.append(np.nan)

                # add output ratio to corresponding setpoint measurements
                cand_A_measurements.output_ratio[cand_A_measurements[setpoint_column] == setpoint] = cand_A_ratio
                # reset temp storage list
                cand_A_ratio = []

            # ----------------------------------------------------------------------------
            # channel B
            if apply_zero_correction:
                cand_B_meas = cand_B_measurements[value_column_zero_cor][cand_B_measurements[setpoint_column] == setpoint].values
            else:
                cand_B_meas = cand_B_measurements[value_column][cand_B_measurements[setpoint_column] == setpoint].values

            if (len(cand_B_meas)>0) :
                cand_B_time = cand_B_measurements.elapsed_mins[cand_B_measurements[setpoint_column] == setpoint].values
                if apply_zero_correction:
                    ref_B_meas = ref_B_measurements[value_column_zero_cor][ref_B_measurements[setpoint_column] == setpoint].values
                else:
                    ref_B_meas = ref_B_measurements[value_column][ref_B_measurements[setpoint_column] == setpoint].values
                #
                time_ref_B = ref_B_measurements.elapsed_mins[ref_B_measurements[setpoint_column] == setpoint].values
                for idx in range(len(cand_B_meas)):
                    try:
                        cand_B_out_ratio = output_ratio(
                            cand_B_meas[idx], 
                            ref_B_meas[idx+1], 
                            ref_B_meas[idx], 
                            cand_B_time[idx], 
                            time_ref_B[idx+1], 
                            time_ref_B[idx])
                        cand_B_ratio.append(cand_B_out_ratio)
                    except IndexError:
                        if not (len(ref_B_meas)>len(cand_B_meas)):
                            print('Channel B: candidate generator has more measurements than reference generator.')
                            cand_B_ratio.append(np.nan)

                # add output ratio to corresponding setpoint measurements
                cand_B_measurements.output_ratio[cand_A_measurements[setpoint_column] == setpoint] = cand_B_ratio
                # reset temp storage list
                cand_B_ratio = []

    # average ratio per setpoint
    # channel A
    avg_ratio_A = cand_A_measurements.groupby(setpoint_column).mean()['output_ratio']
    avg_ratio_A = avg_ratio_A[~np.isnan(avg_ratio_A)]
    avg_ratio_A = avg_ratio_A.to_dict()
    # channel B
    avg_ratio_B = cand_B_measurements.groupby(setpoint_column).mean()['output_ratio']
    avg_ratio_B = avg_ratio_B[~np.isnan(avg_ratio_B)]
    avg_ratio_B = avg_ratio_B.to_dict()

    # ----------------------------------------------------------------------------
    # compute output ratio RSD
    ratio_RSD_A = {}
    ratio_RSD_B = {}

    for setpoint in df[setpoint_column].cat.categories:
        if setpoint == 0.0:
            continue
        else:
            # ----------------------------------------------------------------------------
            # channel A
            ratio_A = cand_A_measurements.output_ratio[cand_A_measurements[setpoint_column] == setpoint]
            if (len(ratio_A)>0) :
                ratio_RSD_A[setpoint] = (np.sqrt(np.sum((ratio_A-avg_ratio_A[setpoint])**2)/
                                                        (len(ratio_A)-1)) / avg_ratio_A[setpoint])
                print('Ratio RSD for channel A at setpoint {}: {:2.2f}%'.format(setpoint,ratio_RSD_A[setpoint]*100))
            # ----------------------------------------------------------------------------
            # channel B
            ratio_B = cand_B_measurements.output_ratio[cand_B_measurements[setpoint_column] == setpoint]
            if (len(ratio_B)>0) :
                ratio_RSD_B[setpoint] = (np.sqrt(np.sum((ratio_B-avg_ratio_B[setpoint])**2)/
                                                        (len(ratio_B)-1)) / avg_ratio_B[setpoint])
                print('Ratio RSD for channel B at setpoint {}: {:2.2f}%'.format(setpoint,ratio_RSD_B[setpoint]*100)) 

    # ----------------------------------------------------------------------------
    # compute concentration of candidate
    cand_A_concentrations = {}
    cand_B_concentrations = {}
    for setpoint in ref_concentration.keys():
        cand_A_concentrations[setpoint] = avg_ratio_A[setpoint]*ref_concentration[setpoint]
        cand_B_concentrations[setpoint] = avg_ratio_B[setpoint]*ref_concentration[setpoint]

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------

    """Perform single point analysis for each point in setpoint_column: compute measurement stability uncertainty"""

    SE_stab_cand_A = {}
    SE_stab_ref_A = {}

    SE_stab_cand_B = {}
    SE_stab_ref_B = {}

    for channel in df.Channel.cat.categories:
        for setpoint in df[setpoint_column].cat.categories:
            for generator in df.Generator.cat.categories:
                meas = df[value_column][(df.Channel == channel) & (df.Zero == False) &
                                        (df[setpoint_column] == setpoint) & (df.Generator == generator)].values
                if (len(meas)>2): # need minimum 3 measurements
                    time = df.elapsed_mins[(df.Channel == channel) & (df.Zero == False) &
                                           (df[setpoint_column] == setpoint) & (df.Generator == generator)].values
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
                            SE_stab_cand_A[setpoint] = SE_stability(meas_var, time_var, meas_time_cov, len(meas))
                        else: # generator == 'Ref'
                            SE_stab_ref_A[setpoint] = SE_stability(meas_var, time_var, meas_time_cov, len(meas))
                    else: # channel == 'B'
                        if generator == 'Cand':
                            SE_stab_cand_B[setpoint] = SE_stability(meas_var, time_var, meas_time_cov, len(meas))
                        else: # generator == 'Ref'
                            SE_stab_ref_B[setpoint] = SE_stability(meas_var, time_var, meas_time_cov, len(meas))

    # ----------------------------------------------------------------------------

    u_stability_A = {}
    u_stability_B = {}

    for setpoint in cand_A_measurements[setpoint_column].cat.categories:
        if setpoint == 0.0:
            continue
        else:
            # ----------------------------------------------------------------------------
            # channel A    
            ratio_A = cand_A_measurements.output_ratio[cand_A_measurements[setpoint_column] == setpoint].values
            # remove nan
            ratio_A = ratio_A[~np.isnan(ratio_A)]
    
            if (len(ratio_A)>0):
                u_A = []
                if apply_zero_correction:
                    cand_A_meas = cand_A_measurements[value_column_zero_cor][cand_A_measurements[setpoint_column] == setpoint].values
                    ref_A_meas = ref_A_measurements[value_column_zero_cor][ref_A_measurements[setpoint_column] == setpoint].values
                else:
                    cand_A_meas = cand_A_measurements[value_column][cand_A_measurements[setpoint_column] == setpoint].values
                    ref_A_meas = ref_A_measurements[value_column][ref_A_measurements[setpoint_column] == setpoint].values
        
                cand_A_time = cand_A_measurements.elapsed_mins[cand_A_measurements[setpoint_column] == setpoint].values
                ref_A_time = ref_A_measurements.elapsed_mins[ref_A_measurements[setpoint_column] == setpoint].values
    
                for idx in range(len(ratio_A)):
                    try:
                        u_A.append(u_stability(
                            ratio_A[idx], 
                            SE_stab_cand_A[setpoint], 
                            SE_stab_ref_A[setpoint], 
                            SE_stab_ref_A[setpoint], 
                            cand_A_meas[idx], 
                            ref_A_meas[idx+1], 
                            ref_A_meas[idx], 
                            cand_A_time[idx], 
                            ref_A_time[idx+1], 
                            ref_A_time[idx]))
                    except IndexError:
                        if not (len(ref_A_meas)>len(cand_A_meas)):
                            print('Channel A: candidate generator has more measurements than reference generator.')
                u_A = np.array(u_A)
                u_stability_A[setpoint] = np.sqrt(np.sum(u_A**2)/(len(u_A)**2))

            # ----------------------------------------------------------------------------
            # channel B    
            ratio_B = cand_B_measurements.output_ratio[cand_B_measurements[setpoint_column] == setpoint].values
            # remove nan
            ratio_B = ratio_B[~np.isnan(ratio_B)]
    
            if (len(ratio_B)>0):
                u_B = []
                if apply_zero_correction:
                    cand_B_meas = cand_B_measurements[value_column_zero_cor][cand_B_measurements[setpoint_column] == setpoint].values
                    ref_B_meas = ref_B_measurements[value_column_zero_cor][ref_B_measurements[setpoint_column] == setpoint].values
                else:
                    cand_B_meas = cand_B_measurements[value_column][cand_B_measurements[setpoint_column] == setpoint].values
                    ref_B_meas = ref_B_measurements[value_column][ref_B_measurements[setpoint_column] == setpoint].values
        
                cand_B_time = cand_B_measurements.elapsed_mins[cand_B_measurements[setpoint_column] == setpoint].values
                ref_B_time = ref_B_measurements.elapsed_mins[ref_B_measurements[setpoint_column] == setpoint].values
    
                for idx in range(len(ratio_B)):
                    try:
                        u_B.append(u_stability(
                            ratio_B[idx], 
                            SE_stab_cand_B[setpoint], 
                            SE_stab_ref_B[setpoint], 
                            SE_stab_ref_B[setpoint], 
                            cand_B_meas[idx], 
                            ref_B_meas[idx+1], 
                            ref_B_meas[idx], 
                            cand_B_time[idx], 
                            ref_B_time[idx+1], 
                            ref_B_time[idx]))
                    except IndexError:
                        if not (len(ref_B_meas)>len(cand_B_meas)):
                            print('Channel B: candidate generator has more measurements than reference generator.')
                u_B = np.array(u_B)
                u_stability_B[setpoint] = np.sqrt(np.sum(u_B**2)/(len(u_B)**2))

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------

    """Perform single point analysis for each point in setpoint_column: 
    compute measurement repeatability uncertainty"""

    u_repeatability_A = {}
    u_repeatability_B = {}

    L = 3 # in this set of measurements there is 1 cand measurement bracketed between 2 ref measurements

    for setpoint in df[setpoint_column].cat.categories:
        if setpoint == 0.0:
            continue
        else:
            # channel A    
            ratio_A = cand_A_measurements.output_ratio[cand_A_measurements[setpoint_column] == setpoint].values
            # remove possible nan values
            ratio_A = ratio_A[~np.isnan(ratio_A)]
            K_A = len(ratio_A)
            if (K_A>0):
                u_repeatability_A[setpoint] = u_repeatability(K_A, L, ratio_A, u_stability_A[setpoint])
    
            # channel B
            ratio_B = cand_B_measurements.output_ratio[cand_B_measurements[setpoint_column] == setpoint].values
            # remove possible nan values
            ratio_B = ratio_B[~np.isnan(ratio_B)]
            K_B = len(ratio_B)
            if (K_B>0):
                u_repeatability_B[setpoint] = u_repeatability(K_B, L, ratio_B, u_stability_B[setpoint])

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------

    """
    Perform single point analysis for each point in setpoint_column:
        combine repeatability and stability uncertainty
        compute reproducibility uncertainty contribution
        get combined uncertainty
    """

    # ----------------------------------------------------------------------------
    # channel A
    u_stab_repeat_A = {}
    u_comparison_A = {}
    u_reprod_A = {}
    u_reference_A = {}
    u_combined_A = {}

    # ----------------------------------------------------------------------------
    # channel B
    u_stab_repeat_B = {}
    u_comparison_B = {}
    u_reprod_B = {}
    u_reference_B = {}
    u_combined_B = {}

    for setpoint in u_stability_A.keys():
        if setpoint == 0.0:
            continue
        else:
            # ----------------------------------------------------------------------------
            # channel A
            u_stab_repeat_A[setpoint] = np.sqrt(
                u_stability_A[setpoint]**2 
                + u_repeatability_A[setpoint]**2)
            u_comparison_A[setpoint] = setpoint*u_stab_repeat_A[setpoint]
            u_reprod_A[setpoint] = reproducibility_uncertainty*cand_A_concentrations[setpoint]
            u_reference_A[setpoint] = u_ref_concentration[setpoint]*avg_ratio_A[setpoint]
            
            # total uncertainty
            u_combined_A[setpoint] = np.sqrt(
                u_comparison_A[setpoint]**2 
                + u_reprod_A[setpoint]**2 
                + u_reference_A[setpoint]**2)

            # ----------------------------------------------------------------------------
            # channel B
            u_stab_repeat_B[setpoint] = np.sqrt(
                u_stability_B[setpoint]**2 
                + u_repeatability_B[setpoint]**2)
            u_comparison_B[setpoint] = setpoint*u_stab_repeat_B[setpoint]
            u_reprod_B[setpoint] = reproducibility_uncertainty*cand_B_concentrations[setpoint]
            u_reference_B[setpoint] = u_ref_concentration[setpoint]*avg_ratio_B[setpoint]
            
            # total uncertainty
            u_combined_B[setpoint] = np.sqrt(
                u_comparison_B[setpoint]**2 
                + u_reprod_B[setpoint]**2 
                + u_reference_B[setpoint]**2)

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    """Write main results to file"""

    # channel A
    # create dataframe where to store uncertainty components
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

    results_A = pd.DataFrame(u_combined_A, index=index_names)
    cols = results_A.columns

    results_A['c_cand'] = [
        'R', 
        '', 
        'u(R)_stability', 
        'u_repeatability',
        'u_comparison',
        'u_reproducibility',
        'u_reference',
        'u(c)',
        'U(c)',
        'U(c)',
        'c',
        'D_i']

    results_A['ng/m3'] = [
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

    cols = cols.insert(0,'ng/m3')
    cols = cols.insert(0, 'c_cand') 
    results_A = results_A[cols]

    for setpoint in u_combined_A.keys():
        results_A.loc['Mean output ratio', setpoint] = avg_ratio_A[setpoint]
        results_A.loc['RSD output ratio', setpoint] = 100.0*ratio_RSD_A[setpoint]
        results_A.loc['Stability uncertainty', setpoint] = u_stability_A[setpoint]
        results_A.loc['Repeatability uncertainty', setpoint] = u_repeatability_A[setpoint]
        results_A.loc['Comparison uncertainty', setpoint] = u_comparison_A[setpoint]
        results_A.loc['Reproducibility uncertainty', setpoint] = u_reprod_A[setpoint]
        results_A.loc['Reference uncertainty', setpoint] = u_reference_A[setpoint]
        results_A.loc['Combined std uncertainty', setpoint] = u_combined_A[setpoint]
        results_A.loc['Expanded uncertainty (k=2)', setpoint] = 2*u_combined_A[setpoint]
        results_A.loc['Relative expanded uncertainty', setpoint] = 100.0*(2*u_combined_A[setpoint])/cand_A_concentrations[setpoint]
        results_A.loc['Candidate concentration', setpoint] = cand_A_concentrations[setpoint]
        results_A.loc['Relative difference in concentration', setpoint] = 100.0*(cand_A_concentrations[setpoint]/setpoint - 1.0)

    # ----------------------------------------------------------------------------
    # channel B
    results_B = pd.DataFrame(u_combined_B, index=index_names)

    results_B['c_cand'] = [
        'R', 
        '', 
        'u(R)_stability', 
        'u_repeatability',
        'u_comparison',
        'u_reproducibility',
        'u_reference',
        'u(c)',
        'U(c)',
        'U(c)',
        'c',
        'D_i']

    results_B['ng/m3'] = [
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

    results_B = results_B[cols]

    for setpoint in u_combined_B.keys():
        results_B.loc['Mean output ratio', setpoint] = avg_ratio_B[setpoint]
        results_B.loc['RSD output ratio', setpoint] = 100.0*ratio_RSD_B[setpoint]
        results_B.loc['Stability uncertainty', setpoint] = u_stability_B[setpoint]
        results_B.loc['Repeatability uncertainty', setpoint] = u_repeatability_B[setpoint]
        results_B.loc['Comparison uncertainty', setpoint] = u_comparison_B[setpoint]
        results_B.loc['Reproducibility uncertainty', setpoint] = u_reprod_B[setpoint]
        results_B.loc['Reference uncertainty', setpoint] = u_reference_B[setpoint]
        results_B.loc['Combined std uncertainty', setpoint] = u_combined_B[setpoint]
        results_B.loc['Expanded uncertainty (k=2)', setpoint] = 2*u_combined_B[setpoint]
        results_B.loc['Relative expanded uncertainty', setpoint] = 100.0*(2*u_combined_B[setpoint])/cand_B_concentrations[setpoint]
        results_B.loc['Candidate concentration', setpoint] = cand_B_concentrations[setpoint]
        results_B.loc['Relative difference in concentration', setpoint] = 100.0*(cand_A_concentrations[setpoint]/setpoint - 1.0)

    # ----------------------------------------------------------------------------

    # write to Excel
    try:
        logging.debug('Add sheet to the Excel workbook')
        wb.sheets.add('Channel A')
        wb.sheets.add('Channel B')
    except ValueError as V:
        logging.debug('Sheet already present in workbook')
        print("Error:", V)

    # # set bold font on first row
    # wb.sheets['Channel A'].range('A1:AAA1').api.Font.Bold = True
    # wb.sheets['Channel B'].range('A1:AAA1').api.Font.Bold = True

    # set bold font on first column
    wb.sheets['Channel A'].range('A1:A20').api.Font.Bold = True
    wb.sheets['Channel B'].range('A1:A20').api.Font.Bold = True

    # write dataset to Excel sheet
    wb.sheets['Channel A'].range('A1').options(index=True).value = results_A
    wb.sheets['Channel B'].range('A1').options(index=True).value = results_B

    wb.sheets['Channel A'].range('A1').value = 'Setpoint candidate generator'
    wb.sheets['Channel B'].range('A1').value = 'Setpoint candidate generator'

    # ----------------------------------------------------------------------------

    """Plot calculated concentrations vs setpoints"""

    dic_color = {'A': 'tab:blue', 'B': 'tab:orange'}
    dic_marker = {'A': 's', 'B': '^'}

    fig = plt.figure('results_singlepoint', figsize=[8,7])
    ax = fig.add_subplot(
        111, 
        xlabel=r'Setpoint [$ng\ m^{-3}$]', 
        ylabel=r'Calc. concentration [$ng\ m^{-3}$]')
    ax.errorbar(
        x=list(cand_A_concentrations.keys()), 
        y=list(cand_A_concentrations.values()), 
        yerr=2*np.array(list(u_combined_A.values())), 
        fmt=dic_marker['A'], 
        color=dic_color['A'], 
        label='A')
    ax.errorbar(
        x=list(cand_B_concentrations.keys()), 
        y=list(cand_B_concentrations.values()), 
        yerr=2*np.array(list(u_combined_B.values())), 
        fmt=dic_marker['B'], 
        color=dic_color['B'], 
        label='B')
    ax.legend(loc='best')
    plt.tight_layout()

    # ----------------------------------------------------------------------------

    # add figure to Excel sheet
    wb.sheets['Channel A'].pictures.add(
        fig, 
        name='single point calibration', 
        update=True, 
        left=wb.sheets['Plot data'].range('L1').left, 
        top=wb.sheets['Plot data'].range('L1').top)

    wb.sheets['Channel B'].pictures.add(
        fig, 
        name='single point calibration', 
        update=True, 
        left=wb.sheets['Plot data'].range('L1').left, 
        top=wb.sheets['Plot data'].range('L1').top)

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------

    """Regression calculated concentrations vs setpoints
    apply WLS to estimated measurement as a function of nominal value"""

    # add one sheet to display the regression results
    sheet_name = 'Regression results'    
    try:
        logging.debug('Add sheet to the Excel workbook')
        wb.sheets.add(sheet_name)
    except ValueError as V:
        logging.debug('Sheet already present in workbook')
        print("Error:", V)

    n_setpoints = len(cand_A_concentrations.keys())

    list_functions = []

    # ----------------------------------------------------------------------------
    # channel A
    wls_regression_A = {}
    wls_pred_A = {}
    wls_aic_A = {}
    wls_bic_A = {}
    wls_aicc_A = {}
    wls_ssr_A = {}

    # ----------------------------------------------------------------------------
    # channel B
    wls_regression_B = {}
    wls_pred_B = {}
    wls_aic_B = {}
    wls_bic_B = {}
    wls_aicc_B = {}
    wls_ssr_B = {}

    # set of points for better display of the regressed line
    setpoint_predict_A = np.arange(
        start=min(list(cand_A_concentrations.keys())), 
        stop=max(list(cand_A_concentrations.keys()))+100, 
        step=100)
    setpoint_predict_B = np.arange(
        start=min(list(cand_B_concentrations.keys())), 
        stop=max(list(cand_B_concentrations.keys()))+100, 
        step=100)

    # create plot of WLS regression lines
    col = {'poly0':'tab:blue', 'poly1':'tab:orange', 'poly2':'tab:green', 'poly3':'tab:red', 'poly4':'tab:purple', 'constant':'tab:cyan'}
    fig = plt.figure('WLS regression lines', figsize=[9,10])
    ax1 = fig.add_subplot(
        212, 
        title='Channel B', 
        xlabel=r'Setpoint $[ng\ m^{-3}]$', 
        ylabel=r'Calc. concentration $[ng\ m^{-3}]$')
    ax0 = fig.add_subplot(
        211, 
        title='Channel A', 
        ylabel=r'Calc. concentration $[ng\ m^{-3}]$', 
        sharex=ax1)

    # WLS 
    # remove reference uncertainty from combined uncertainty to be used in WLS regression
    wls_weights_A = np.sqrt(
        np.power(list(u_combined_A.values()),2)-np.power(list(u_reference_A.values()),2)
    )
    wls_weights_B = np.sqrt(
        np.power(list(u_combined_B.values()),2)-np.power(list(u_reference_B.values()),2)
    )

    # write data to file for validation of regression code
    # input_data = pd.DataFrame({
    #     'setpoint' : list(cand_A_concentrations.keys()), 
    #     'concentration_chA' : list(cand_A_concentrations.values()),
    #     'uncertainty_chA' : wls_weights_A,
    #     'concentration_chB' : list(cand_B_concentrations.values()),
    #     'uncertainty_chB' : wls_weights_B})

    # input_data.to_csv('input_data.txt', sep=' ', index=False)
    # print(os.getcwd())

    if n_setpoints > 2: # ISO6143 asks for a minimum of 3 points 
        for poly_deg in range(max(2,n_setpoints-2)):  
            # -2 is necessary to avoid division by zero in the calculation of AICc
            # but we still want to fit a polynomial order 1 with 3 data points
            function = 'poly{}'.format(poly_deg)
            list_functions.append(function)
        
            # ----------------------------------------------------------------------------
            # channel A
            # apply WLS regression
            wls_regression_A[function] = regress_function(poly_deg, 
                list(cand_A_concentrations.keys()), 
                list(cand_A_concentrations.values()), 
                wls_weights_A)
    
            # get WLS predictions
            wls_pred_A[function] = predict(
                poly_deg, 
                wls_regression_A[function], 
                setpoint_predict_A) 

            plot_poly(
                x=list(cand_A_concentrations.keys()), 
                x_pred=setpoint_predict_A,
                wls_regression=wls_regression_A[function],
                wls_pred=wls_pred_A[function],
                axis=ax0,
                color_line=col[function],
                label_name=function)
        
            # get AIC, BIC and AICc
            wls_aic_A[function] = wls_regression_A[function].aic
            wls_bic_A[function] = wls_regression_A[function].bic
            n_params = len(wls_regression_A[function].params) 
            try:
                wls_aicc_A[function] = wls_regression_A[function].aic + (2*n_params*(n_params+1))/(n_setpoints-n_params-1)
            except ZeroDivisionError as Z:
                logging.debug('Zero division in the calculation of AICc. Use AIC instead.')
                print("Error:", Z)
                wls_aicc_A[function] = wls_regression_A[function].aic
            wls_ssr_A[function] = wls_regression_A[function].ssr
    
            # ----------------------------------------------------------------------------
            # channel B
            # apply WLS regression
            wls_regression_B[function] = regress_function(poly_deg, 
                list(cand_B_concentrations.keys()), 
                list(cand_B_concentrations.values()), 
                wls_weights_B)
    
            # get WLS predictions
            wls_pred_B[function] = predict(
                poly_deg, 
                wls_regression_B[function], 
                setpoint_predict_B) 

            plot_poly(
                x=list(cand_B_concentrations.keys()), 
                x_pred=setpoint_predict_B,
                wls_regression=wls_regression_B[function],
                wls_pred=wls_pred_B[function],
                axis=ax1,
                color_line=col[function],
                label_name=function)    
    
            # get AIC, BIC and AICc
            wls_aic_B[function] = wls_regression_B[function].aic
            wls_bic_B[function] = wls_regression_B[function].bic
            n_params = len(wls_regression_A[function].params) 
            try:
                wls_aicc_B[function] = wls_regression_B[function].aic + (2*n_params*(n_params+1))/(n_setpoints-n_params-1)
            except ZeroDivisionError as Z:
                logging.debug('Zero division in the calculation of AICc. Use AIC instead.')
                print("Error:", Z)
                wls_aicc_B[function] = wls_regression_B[function].aic 
            wls_ssr_B[function] = wls_regression_B[function].ssr
    
    else:
        print('Number of setpoints is not enough to fit any function. Minimum requirement: 3 points')
    
    # add data points to the plot
    ax0.errorbar(
        x=list(cand_A_concentrations.keys()), 
        y=list(cand_A_concentrations.values()), 
        yerr=2*np.array(list(u_combined_A.values())), 
        fmt='o', 
        color='black', 
        elinewidth=2, 
        barsabove=True, 
        label='data')
    ax1.errorbar(
        x=list(cand_B_concentrations.keys()), 
        y=list(cand_B_concentrations.values()), 
        yerr=2*np.array(list(u_combined_B.values())), 
        fmt='o', 
        color='black', 
        elinewidth=2, 
        barsabove=True)

    ax0.set_ylim([0.0, 1.05*max(cand_A_concentrations.values())])
    ax1.set_ylim([0.0, 1.05*max(cand_B_concentrations.values())])

    ax0.tick_params('x', labelbottom=False)
    ax0.legend(loc='best', ncol=1, bbox_to_anchor=(1., .5, 0.5, 0.5))
    plt.tight_layout()

    # add figure to Excel sheet
    wb.sheets[sheet_name].pictures.add(
        fig, 
        name='WLS regression', 
        update=True, 
        left=wb.sheets['Plot data'].range('A1').left, 
        top=wb.sheets['Plot data'].range('A1').top)

    # ----------------------------------------------------------------------------

    # get key corresponding to optimal polynomial
    best_poly_A = min(wls_aicc_A, key=wls_aicc_A.get)
    best_poly_B = min(wls_aicc_B, key=wls_aicc_B.get)

    nparams_best_poly_A = len(wls_regression_A[best_poly_A].params)
    nparams_best_poly_B = len(wls_regression_B[best_poly_B].params)

    # create list of indices to reflect the terms involved in the calibration function
    index_df_poly_coef = []
    for idx in range(max(nparams_best_poly_A, nparams_best_poly_B)):
        index_df_poly_coef.append('b_{}'.format(idx))

    # make a string with the formula of the optimal polynomials
    formula_A = 'b_0*x^0 '
    for idx in range(1,nparams_best_poly_A):
        formula_A += '+ b_{}*x^{} '.format(idx, idx)

    formula_B = 'b_0*x^0 '
    for idx in range(1,nparams_best_poly_B):
        formula_B += '+ b_{}*x^{} '.format(idx, idx)

    # make dataframe with the parameters of the optimal polynomials
    df_poly_coef_A = pd.DataFrame(
        index=index_df_poly_coef, 
        columns=['Parameters', 'Standard error', 'P-values'])

    df_poly_coef_B = pd.DataFrame(
        index=index_df_poly_coef, 
        columns=['Parameters', 'Standard error', 'P-values'])

    if nparams_best_poly_A < nparams_best_poly_B:
        df_poly_coef_A['Parameters'] = np.append(wls_regression_A[best_poly_A].params,
            [0]*(nparams_best_poly_B-nparams_best_poly_A))
        df_poly_coef_A['Standard error'] = np.append(wls_regression_A[best_poly_A].bse,
            [0]*(nparams_best_poly_B-nparams_best_poly_A))
        df_poly_coef_A['P-values'] = np.append(wls_regression_A[best_poly_A].pvalues,
            [0]*(nparams_best_poly_B-nparams_best_poly_A))
    else:
        df_poly_coef_A['Parameters'] = wls_regression_A[best_poly_A].params
        df_poly_coef_A['Standard error'] = wls_regression_A[best_poly_A].bse
        df_poly_coef_A['P-values'] = wls_regression_A[best_poly_A].pvalues

    if nparams_best_poly_B < nparams_best_poly_A:
        df_poly_coef_B['Parameters'] = np.append(wls_regression_B[best_poly_B].params,
            [0]*(nparams_best_poly_A-nparams_best_poly_B))
        df_poly_coef_B['Standard error'] = np.append(wls_regression_B[best_poly_B].bse,
            [0]*(nparams_best_poly_A-nparams_best_poly_B))
        df_poly_coef_B['P-values'] = np.append(wls_regression_B[best_poly_B].pvalues,
            [0]*(nparams_best_poly_A-nparams_best_poly_B))
    else:
        df_poly_coef_B['Parameters'] = wls_regression_B[best_poly_B].params
        df_poly_coef_B['Standard error'] = wls_regression_B[best_poly_B].bse
        df_poly_coef_B['P-values'] = wls_regression_B[best_poly_B].pvalues

    # ----------------------------------------------------------------------------

    """Plot weighted residuals"""
    fig = plt.figure('weighted residuals', figsize=[9,10])
    ax1 = fig.add_subplot(
        212, 
        title='Channel B', 
        xlabel=r'Setpoint $[ng\ m^{-3}]$', 
        ylabel=r'Weighted residuals $[ng\ m^{-3}]$')
    ax0 = fig.add_subplot(
        211, 
        title='Channel A', 
        ylabel=r'Weighted residuals $[ng\ m^{-3}]$', 
        sharex=ax1)

    for poly_deg in range(max(2,n_setpoints-2)):
        function = list_functions[poly_deg]
        # ----------------------------------------------------------------------------
        # channel A
        ax0.scatter(
            x=list(cand_A_concentrations.keys()), 
            y=(wls_regression_A[function].resid)/wls_weights_A, 
            color=col[function], 
            label='{}'.format(function))

        # ----------------------------------------------------------------------------
        # channel B
        ax1.scatter(
            x=list(cand_B_concentrations.keys()), 
            y=(wls_regression_B[function].resid)/wls_weights_B, 
            color=col[function], 
            label='{}'.format(function))

    ax0.hlines(
        y=0, 
        xmin=min(list(cand_A_concentrations.keys())), 
        xmax=max(list(cand_A_concentrations.keys())), 
        lw=2, 
        color='black')
    ax1.hlines(
        y=0, 
        xmin=min(list(cand_B_concentrations.keys())), 
        xmax=max(list(cand_B_concentrations.keys())), 
        lw=2, 
        color='black')

    ax0.set_ylim([-1.0, +1.0])
    ax1.set_ylim([-1.0, +1.0])

    ax0.legend(loc='best', ncol=1, bbox_to_anchor=(1., .5, 0.5, 0.5))
    plt.tight_layout()

    # add figure to Excel sheet
    wb.sheets[sheet_name].pictures.add(
        fig, 
        name='WLS residuals', 
        update=True, 
        left=wb.sheets['Plot data'].range('N1').left, 
        top=wb.sheets['Plot data'].range('N1').top)

    # ----------------------------------------------------------------------------

    # write AICc and SSR values
    df_chA = pd.DataFrame.from_dict(wls_aicc_A, orient='index')
    df_chA.columns = ['AICc']
    df_chA['Sum squared residuals'] = wls_ssr_A.values()

    df_chB = pd.DataFrame.from_dict(wls_aicc_B, orient='index')
    df_chB.columns = ['AICc']
    df_chB['Sum squared residuals'] = wls_ssr_B.values()

    # set bold font on header
    wb.sheets[sheet_name].range('Y1:AO1').api.Font.Bold = True
    wb.sheets[sheet_name].range('AE15:AN15').api.Font.Bold = True
    wb.sheets[sheet_name].range('Y10:AB10').api.Font.Bold = True

    # set bold font on index
    wb.sheets[sheet_name].range('Y1:Y20').api.Font.Bold = True
    wb.sheets[sheet_name].range('AE1:AE23').api.Font.Bold = True

    # write AICc and SSR dataset to Excel sheet
    wb.sheets[sheet_name].range('Y1').options(index=True).value = df_chA
    wb.sheets[sheet_name].range('Y1').value = 'Channel A'

    wb.sheets[sheet_name].range('Y10').options(index=True).value = df_chB
    wb.sheets[sheet_name].range('Y10').value = 'Channel B'

    wb.sheets[sheet_name].range('AE1').options(index=True).value = df_poly_coef_A
    wb.sheets[sheet_name].range('AE1').value = 'Channel A'
    # change background color of cell
    wb.sheets[sheet_name].range('AE1').color = (255, 165, 0) # RGB

    wb.sheets[sheet_name].range('AN1').value = 'Covariance matrix'
    wb.sheets[sheet_name].range('AP1').options(ndim=1, expand='horizontal').value = df_poly_coef_A.index.values
    wb.sheets[sheet_name].range('AO2').options(transpose=True).value = df_poly_coef_A.index.values
    wb.sheets[sheet_name].range('AP2').options(expand='table').value = wls_regression_A[best_poly_A].cov_params()

    wb.sheets[sheet_name].range('AE8').options(index=True).value = 'Formula of optimal polynomial for channel A:'
    wb.sheets[sheet_name].range('AJ8').options(index=True).value = formula_A

    # ----------------------------------------------------------------------------

    wb.sheets[sheet_name].range('AE15').options(index=True).value = df_poly_coef_B
    wb.sheets[sheet_name].range('AE15').value = 'Channel B'
    # change background color of cell
    wb.sheets[sheet_name].range('AE15').color = (255, 165, 0) # RGB

    wb.sheets[sheet_name].range('AN15').value = 'Covariance matrix'
    wb.sheets[sheet_name].range('AP15').options(ndim=1, expand='horizontal').value = df_poly_coef_B.index.values
    wb.sheets[sheet_name].range('AO16').options(transpose=True).value = df_poly_coef_B.index.values
    wb.sheets[sheet_name].range('AP16').options(expand='table').value = wls_regression_B[best_poly_B].cov_params()
    
    wb.sheets[sheet_name].range('AE23').options(index=True).value = 'Formula of optimal polynomial for channel B:'
    wb.sheets[sheet_name].range('AJ23').options(index=True).value = formula_B

    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------

    """Check if the regressions of the two channels agree with each other
    and, if they do, compute the average regression line"""

    # add one sheet to display the regression results
    sheet_name = 'Calibration results'
    try:
        logging.debug('Add sheet to the Excel workbook')
        wb.sheets.add(sheet_name)
    except ValueError as V:
        logging.debug('Sheet already present in workbook')
        print("Error:", V)

    # show overlap of the two regressed functions
    fig = plt.figure('WLS regression lines', figsize=[10,5])
    ax0 = fig.add_subplot(
        111, 
        title='Channel A and Channel B', 
        xlabel=r'Setpoint $[ng\ m^{-3}]$',
        ylabel=r'Calc. concentration $[ng\ m^{-3}]$')
    # Channel A
    plot_poly(
        x=list(cand_A_concentrations.keys()), 
        x_pred=setpoint_predict_A,
        wls_regression=wls_regression_A[best_poly_A],
        wls_pred=wls_pred_A[best_poly_A],
        axis=ax0,
        color_line='tab:blue',
        label_name='Channel A',
        transparancy=.5)
    # Channel B
    plot_poly(
        x=list(cand_B_concentrations.keys()), 
        x_pred=setpoint_predict_B,
        wls_regression=wls_regression_B[best_poly_B],
        wls_pred=wls_pred_B[best_poly_B],
        axis=ax0,
        color_line='tab:green',
        label_name='Channel B',
        transparancy=.5)

    ax0.set_ylim([.95*min(cand_A_concentrations.values()), 1.05*max(cand_A_concentrations.values())])

    ax0.legend(loc='best', ncol=1, bbox_to_anchor=(1., .5, 0.5, 0.5))
    plt.tight_layout()

    # add figure to Excel sheet
    wb.sheets[sheet_name].pictures.add(
        fig, 
        name='WLS regression comparison', 
        update=True, 
        left=wb.sheets[sheet_name].range('A1').left, 
        top=wb.sheets[sheet_name].range('A1').top)

    #########################################################################
    # check if both channels have the same type of optimal function
    if (len(wls_regression_A[best_poly_A].params) != 
        len(wls_regression_B[best_poly_B].params)):
        # channels have different optimal function 
        wb.sheets[sheet_name].range('A25').options(index=True).value = (
            'Channels have different optimal functions.')
        wb.sheets[sheet_name].range('A26').options(index=True).value = (
            'Optimal polynomial for channel A has degree {}.'.format(len(wls_regression_A[best_poly_A].params)-1))
        wb.sheets[sheet_name].range('A27').options(index=True).value = (
            'Optimal polynomial for channel B has degree {}.'.format(len(wls_regression_B[best_poly_B].params)-1))

        # use simplest polynomial to compute calibration function
        degree_poly_calibration = min(len(wls_regression_A[best_poly_A].params)-1,
                              len(wls_regression_B[best_poly_B].params)-1)
        wb.sheets[sheet_name].range('A28').options(index=True).value = (
            'Using polynomial of degree {} to compute calibration function.'.format(degree_poly_calibration))
    else:
        # channels have same optimal function
        degree_poly_calibration = len(wls_regression_A[best_poly_A].params)-1

    # compute average regression line
    wls_regression_avg = regress_function(
        degree_poly_calibration, 
        list(cand_A_concentrations.keys())+list(cand_B_concentrations.keys()),
        list(cand_A_concentrations.values())+list(cand_B_concentrations.values()),
        np.concatenate((wls_weights_A, wls_weights_B)))

    # get WLS predictions
    wls_pred_avg = predict(
        len(wls_regression_A[best_poly_A].params)-1, 
        wls_regression_avg, 
        setpoint_predict_A) 

    # get WLS 95% CI
    wls_pred_err_avg, wls_lower_avg, wls_upper_avg = wls_prediction_std(wls_regression_avg, alpha=0.05)

    # plot regression line
    fig = plt.figure('WLS average regression line', figsize=[8,5])
    ax0 = fig.add_subplot(
        111, 
        title='Average regression line', 
        xlabel=r'Setpoint $[ng\ m^{-3}]$',
        ylabel=r'Calc. concentration $[ng\ m^{-3}]$')
    # Channel A
    ax0.plot(
        setpoint_predict_A,
        wls_pred_avg,
        lw=2,
        color='tab:orange')
    ax0.plot(
        list(cand_A_concentrations.keys())+list(cand_B_concentrations.keys()),
        wls_lower_avg,
        lw=1,
        color='tab:orange',
        ls='--')
    ax0.plot(
        list(cand_A_concentrations.keys())+list(cand_B_concentrations.keys()), 
        wls_upper_avg, 
        lw=1, 
        color='tab:orange', 
        ls='--')

    ax0.set_ylim([.95*min(cand_A_concentrations.values()), 1.05*max(cand_A_concentrations.values())])

    plt.tight_layout()

    # add figure to Excel sheet
    wb.sheets[sheet_name].pictures.add(
        fig, 
        name='WLS average regression', 
        update=True, 
        left=wb.sheets[sheet_name].range('P1').left, 
        top=wb.sheets[sheet_name].range('P1').top)

    # compute chi2
    chi2 = chi2_WLS(
        list(cand_A_concentrations.values())+list(cand_B_concentrations.values()),
        list(cand_A_concentrations.keys())+list(cand_B_concentrations.keys()),
        wls_regression_avg.params,
        list(u_combined_A.values())+list(u_combined_B.values()))

    n_params_chi2 = len(list(cand_A_concentrations.keys())+list(cand_B_concentrations.keys())) - len(wls_regression_avg.params)

    prob_chi2 = 1-gammainc(n_params_chi2/2,chi2/2)

    # make dataframe with the parameters of the optimal polynomials
    df_poly_coef_avg = pd.DataFrame(
        index=index_df_poly_coef[:(degree_poly_calibration+1)], 
        columns=['Parameters', 'Standard error', 'P-values', 'Sum squared residuals'])

    df_poly_coef_avg['Parameters'] = wls_regression_avg.params
    df_poly_coef_avg['Standard error'] = wls_regression_avg.bse
    df_poly_coef_avg['P-values'] = wls_regression_avg.pvalues
    df_poly_coef_avg['Sum squared residuals'] = wls_regression_avg.ssr

    # write to Excel
    wb.sheets[sheet_name].range('P25').options(index=True).value = df_poly_coef_avg
    wb.sheets[sheet_name].range('P25').value = 'Average regression'
    # change background color of cell
    wb.sheets[sheet_name].range('P25').color = (255, 165, 0) # RGB

    wb.sheets[sheet_name].range('X25').value = 'Covariance matrix'
    wb.sheets[sheet_name].range('Z25').options(ndim=1, expand='horizontal').value = df_poly_coef_avg.index.values
    wb.sheets[sheet_name].range('Y26').options(transpose=True).value = df_poly_coef_avg.index.values
    wb.sheets[sheet_name].range('Z26').options(expand='table').value = wls_regression_avg.cov_params()

    wb.sheets[sheet_name].range('P25:Y25').api.Font.Bold = True

    # chi2
    wb.sheets[sheet_name].range('P30').value = 'Chi squared = {}'.format(chi2)
    wb.sheets[sheet_name].range('P31').value = 'Probability that such a chi squared value should occur by chance = {}'.format(prob_chi2)

    if (prob_chi2>=0.05):
        wb.sheets[sheet_name].range('P34').value = 'The average regressed polynomial is acceptable.'
    else:
        wb.sheets[sheet_name].range('P34').value = 'The average regressed polynomial is NOT acceptable.'

