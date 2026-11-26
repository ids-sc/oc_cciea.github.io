import os
import xarray as xr
import pandas as pd
import numpy as np
import seawater as sw
from numpy import ma

# ----------------------------------------------------------------------
# --BEGIN: Change These
# ----------------------------------------------------------------------
# Creates two netcdf files
# 1) 'prelim_all_ctd.nc' -- has most the data interopolated on bottle depth
# 2) 'prelim_ctd.nc' -- has the best available columns interpolated 1 m 
#                       intervals, data to be appended to the bottle database



# --max depth of data matrix (time x depth) for each station
dpth1_int = np.array([0, 2, 4, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 125, 150, 200, 250, 300, 350, 400, 450, 500, 520, 540])

# -- max depth of data matrix (time x depth) used to append to the bottle database
max_b_dpth = 500
dpth2_int = np.arange(0, max_b_dpth+1)

# --output file type extension
file_type = 'nc'
file_s = './data_gha/CalCofiCSV/calcofi_stations66.csv'
dir_ctd_d = './data_gha/CalCofiCSV/data_from_calcofi_org/'
# file_s = './oc_cciea.github.io/data_gha/CalCofiCSV/calcofi_stations66.csv'
# dir_ctd_d = './oc_cciea.github.io/data_gha/CalCofiCSV/data_from_calcofi_org/'

study_d = ['2107SR', '2111SR',
           '2204SH', '2208BH', '2211',
           '2301RL', '2304SH', '2307SR', '2311SR',
           '2401RL', '2404SH', '2408SR', '2411SR',
           '2501RL', '2504SH', '2507', '2511',
           '2601', '2604', '2607']

yr_d = [2021, 2021,
        2022, 2022, 2022,
        2023, 2023, 2023, 2023,
        2024, 2024, 2024, 2024,
        2025, 2025, 2025, 2025,
        2026, 2026, 2026]

qrtr_d = [3, 4,
          2, 3, 4,
          1, 2, 3, 4,
          1, 2, 3, 4,
          1, 2, 3, 4,
          1, 2, 3]

in_var2_b_list = [5, 0,
                  4, 0, 0,
                  5, 4, 0, 2,
                  1, 0, 0, 0,
                  0, 6, 3, 3,
                  0, 0, 0]

# ---------------------STATION
# line_s = 'Line '
line_s = 'Line'
sttn_s = 'Station'

# ---------------------CTD
line_d = 'Line'
sttn_d = 'Sta'
dpth_d = 'Depth'
cc_d = 'Ord_Occ'
time_d = 'Date_Time_PST'
lon_d = 'Lon_Dec'
lat_d = 'Lat_Dec'
var1_d_wnt = [line_d, sttn_d, lat_d, lon_d]

# can have different columns based on the CTD CSV and the available columns,
# another one can be created if these don't account for columns with data
# NOTE: Rasmus cautioned using the O2 station avg data, suggested using Ox1 or Ox1_StaCorr
var1_d_lbl = ['line', 'sttn', 'lat', 'lon', 'quarter', 'month']
var2_d_wnt = ['Temp1', 'Temp2', 'TempAve', 'Salt1', 'Salt2', 'Salt1_Corr', 'Salt2_Corr', 'SigThetaTS1', 'SigThetaTS2', 'DynHt', 'Ox1', 'Ox1Q', 'Ox1_CruiseCorr', 'Ox1_StaCorr', 'Ox2', 'Ox2Q', 'Ox2_CruiseCorr', 'Ox2_StaCorr', 'OxAve_StaCorr', 'OxSat1', 'OxSat2', 'OxB', 'SIL', 'PO4', 'NO3', 'NO2', 'Chl-a', 'Phaeo']

# use the bottle labels for the 
var2_b_lbl = ['R_POTEMP', 'R_SALINITY', 'R_SIGMA', 'R_DYNHT', 'R_O2', 'R_O2Sat', 'R_SIO3', 'R_PO4', 'R_NO3', 'R_NO2', 'R_CHLA', 'R_PHAEO']


# ----------------------Nutrient CSV
file_n = './data_gha/CalCofiCSV/data_from_calcofi_org/compiled_finalized_nutrients_2111_2504_ALE.csv'
file_n = './oc_cciea.github.io/data_gha/CalCofiCSV/data_from_calcofi_org/compiled_finalized_nutrients_2111_2504_ALE.csv'

# ------------------Nutrient
line_n = 'Line'
sttn_n = 'Sta'
dpth_n = 'Depth'
cc_n = 'Ord_Occ'
time_n = 'Date_Time_PST'
lon_n = 'Lon_Dec'
lat_n = 'Lat_Dec'
var1_n_wnt = [line_n, sttn_n, lat_n, lon_n]

var2_n_wnt = ['NO3', 'NO2', 'PO4', 'SIL']

# -- output NetCDF file
var2_d_lbl = var2_d_wnt
var1_d_lbl = ['line', 'sttn', 'lat', 'lon', 'quarter', 'month']

# output directory, save two data files
# fn1_out: can be used for diagnostic as it has all the data saved, especially oxygen
# fn2_out: will be a curated dataset using the best columns so that it can be append to the bottle
dir_out = './data_gha/CalCofiCSV/'
# dir_out = './data_p14/CalCofiCSV/'
fn1_out = 'prelim_all_ctd.nc'
fn2_out = 'prelim_ctd.nc'

# ----------------------------------------------------------------------
# --END: Change These
# ----------------------------------------------------------------------

# --get the size of input variables
num_var1_d = len(var1_d_wnt)
num_qrtr_d = len(qrtr_d)
num_var2_n = len(var2_n_wnt)
num_int1 = np.size(dpth1_int)
num_int2 = np.size(dpth2_int)
num_var2_d = len(var2_d_wnt)

# --------------------------------------------------------------------
# --Lines and Stations
# --------------------------------------------------------------------

# --get line/sttn of 66 "core" stations
df_s = pd.read_csv(file_s)
sttn_ps = df_s[line_s].astype(str)+'_'+df_s[sttn_s].astype(str)
sttn_wnt = sttn_ps.values
num_sttn = len(sttn_wnt)

# --------------------------------------------------------------------
# --Process CTD data
# --------------------------------------------------------------------

# -- create CTD data vector (sttns X variables_1d X quarters)
#    note: +2 to num_var1_d is for the quarter and month variable have to input
#          these by hand since there is no column data with this info
var1_mtrx_d = np.zeros([num_sttn, len(var1_d_lbl), num_qrtr_d])*np.nan

# --create ctd data matrix (sttns X variables_2d X depth X quarters)
var2_mtrx_d = np.zeros([num_sttn, num_var2_d, num_int1, num_qrtr_d])*np.nan

# -- create Year and Julian date vectors to create a date, do this to
#    match year and jd columns in the pd df of the Cast file
yr_vec_d = np.zeros([num_sttn, num_qrtr_d])*np.nan
jd_vec_d = np.zeros([num_sttn, num_qrtr_d])*np.nan

for ii in range(0, num_qrtr_d):
    # --create ctd filename and open as pd Dataframe
    d_file = '{}{}_CTD_cast_{}_q{}_subset.csv'.format(dir_ctd_d, study_d[ii], yr_d[ii], qrtr_d[ii])
    df_d = pd.read_csv(d_file, low_memory=False, na_values=' ')

    # --change the "generic" Dataframe index to an datetime index
    # df_d['Datetime'] = pd.to_datetime(df_d[time_d].astype(str), infer_datetime_format=True)
    df_d['Datetime'] = pd.to_datetime(df_d[time_d].astype(str), format='mixed')
    df_d = df_d.set_index('Datetime')

    # --get year and jd from Datetime index, need these if more than one
    # --cast per station, will take average of the dates for final date
    jd = df_d.index.dayofyear
    year = df_d.index.year
    month = df_d.index.month

    # --
    for jj in range(0, num_sttn):
        # --indices of all casts for a give station
        in_d_wnt = ((df_d[line_d] == df_s[line_s][jj]) &
                    (df_d[sttn_d] == df_s[sttn_s][jj]))
        # indx_d_wnt = in_d_wnt.nonzero()[0]
        indx_d_wnt = in_d_wnt.values.nonzero()[0]

        # --unique consecative cast number, ie number of casts at the station
        # cc_d_unq = np.unique(df_d[cc_d][indx_d_wnt])
        cc_d_unq = np.unique(df_d[cc_d].values[indx_d_wnt])
        num_d_unq = cc_d_unq.size

        if num_d_unq > 0:
            # --create data matrix (variables X interp depths X unique dates)
            unq_mtrx_d = np.zeros([num_var2_d, num_int1, num_d_unq])*np.nan
            yr_unq = np.zeros([num_d_unq])
            mon_unq = np.zeros([num_d_unq])
            jd_unq = np.zeros([num_d_unq])
            lon_unq = np.zeros([num_d_unq])
            lat_unq = np.zeros([num_d_unq])

            for kk in range(0, num_d_unq):
                # --indices of an individual cast
                in_dq = (df_d[cc_d][in_d_wnt] == cc_d_unq[kk])
                indx_dq = in_dq.values.nonzero()[0]

                # --depths of profile, used to interpolate each variable
                dpth_dq = df_d[dpth_d].iloc[indx_d_wnt[indx_dq]].values
                yr_unq[kk] = year[indx_d_wnt[in_dq]].values[0]
                mon_unq[kk] = month[indx_d_wnt[in_dq]].values[0]
                jd_unq[kk] = jd[indx_d_wnt[in_dq]].values[0]
                lon_unq[kk] = df_d[lon_d].iloc[indx_d_wnt[in_dq]].values[0]
                lat_unq[kk] = df_d[lat_d].iloc[indx_d_wnt[in_dq]].values[0]

                for ll in range(0, num_var2_d):
                    lbl_var_wnt = var2_d_wnt[ll]

                    # interpolate
                    var_ll = df_d[lbl_var_wnt].iloc[indx_d_wnt[indx_dq]].values
                    # --masked array for missing values
                    vard = ma.array(var_ll, mask=np.isnan(var_ll))
                    dpthd = ma.array(dpth_dq, mask=np.isnan(var_ll))
                    if dpthd.compressed().size > 0:
                        data_int = np.interp(dpth1_int, dpthd.compressed(
                        ), vard.compressed(), left=np.nan, right=np.nan)
                    else:
                        data_int = np.zeros(num_int1)*np.nan

                    unq_mtrx_d[ll, :, kk] = data_int

                # recalculate potential temperature and replace the one that
                # was in the CSV file, want to do this because sometimes the
                # CSV file doesn't have data in the PoT1 column
                t1 = np.squeeze(unq_mtrx_d[0, :, kk])
                s1 = np.squeeze(unq_mtrx_d[1, :, kk])
                datad = sw.ptmp(s1, t1, dpth1_int, 0)
                unq_mtrx_d[0, :, kk] = datad

            # --2D variables, take mean of all casts to get a single profile
            # --for a given station
            var2_mtrx_d[jj, :, :, ii] = np.nanmean(unq_mtrx_d, axis=2)

            # --1D variables, take mean of yr and jd to get single date
            yr_vec_d[jj, ii] = np.mean(yr_unq)
            jd_vec_d[jj, ii] = np.mean(jd_unq)
            var1_mtrx_d[jj, 0, ii] = df_s[line_s][jj]
            var1_mtrx_d[jj, 1, ii] = df_s[sttn_s][jj]
            var1_mtrx_d[jj, 2, ii] = np.mean(lat_unq)
            var1_mtrx_d[jj, 3, ii] = np.mean(lon_unq)
            var1_mtrx_d[jj, 4, ii] = qrtr_d[ii]
            var1_mtrx_d[jj, 5, ii] = np.mean(mon_unq)
# --------------------------------------------------------------------
# --Process Nutrient data
# --------------------------------------------------------------------
# Use the nutrient CSV to add the nutrients (most of the CTD csv files
# do not have nutrient data, but some do -- such as Summer 2021).
df_n = pd.read_csv(file_n)
study_longname_n = np.unique(df_n['Study'].values)

study_all_n = []
for i in range(len(study_longname_n)):
    study_all_n.append(study_longname_n[i][0:4])

study_all_d = []
for i in range(num_qrtr_d):
    study_type1 = study_d[i][0:4]
    study_all_d.append(study_type1)

num_qrtr_n = len(study_all_n)

for ii in range(num_qrtr_n):
    df_n1 = df_n.loc[df_n['Study']==study_all_n[ii]]

    # index that the nutrient data matches the quarter of the ctd_prelim
    in_qrtr_mtrx = np.where(np.array(study_all_d)==study_all_n[ii])[0]

    # --change the "generic" Dataframe index to an datetime index
    df_n1['Datetime'] = pd.to_datetime(df_n1[time_d].astype(str))
    df_n1 = df_n1.set_index('Datetime')

    # --get year and jd from Datetime index, need these if more than one
    # --cast per station, will take average of the dates for final date
    jd = df_n1.index.dayofyear
    year = df_n1.index.year
    month = df_n1.index.month

    # --
    for jj in range(0, num_sttn):
        # --indices of all casts for a give station
        in_n_wnt = ((df_n1[line_d] == df_s[line_s][jj]) &
                    (df_n1[sttn_d] == df_s[sttn_s][jj]))
        indx_n_wnt = in_n_wnt.values.nonzero()[0]

        # --unique consecative cast number, ie number of casts at the station
        cc_n_unq = np.unique(df_n1[cc_n].iloc[indx_n_wnt])
        num_n_unq = cc_n_unq.size

        if num_n_unq > 0:
            # --create data matrix (variables X interp depths X unique dates)
            unq_mtrx_n = np.zeros([num_var2_n, num_int1, num_n_unq])*np.nan

            for kk in range(0, num_n_unq):
                # --indices of an individual cast
                in_nq = (df_n1[cc_n][in_n_wnt] == cc_n_unq[kk])
                indx_nq = in_nq.values.nonzero()[0]

                # --depths of profile, used to interpolate each variable
                dpth_nq = df_n1[dpth_n].iloc[indx_n_wnt[indx_nq]].values

                for ll in range(0, num_var2_n):
                    var_ll = df_n1[var2_n_wnt[ll]].iloc[indx_n_wnt[indx_nq]].values
                    # --masked array for missing values
                    vard = ma.array(var_ll, mask=np.isnan(var_ll))
                    dpthn = ma.array(dpth_nq, mask=np.isnan(var_ll))
                    if dpthn.compressed().size > 0:
                        data_int = np.interp(
                            dpth_int, dpthn.compressed(), vard.compressed(),
                            left=np.nan, right=np.nan)
                    else:
                        data_int = np.zeros(num_int1)*np.nan

                    unq_mtrx_n[ll, :, kk] = data_int

            unq_mtrx_n_mn = np.nanmean(unq_mtrx_n, axis=2)
            for ll in range(num_var2_n):
                in_var_mtrx = np.where(np.array(var2_d_wnt) == var2_n_wnt[ll])[0]
                var2_mtrx_d[jj, in_var_mtrx, :, in_qrtr_mtrx] = unq_mtrx_n_mn[ll, :]

# Create a netcdf file with the prelim data
indx_time_d = np.arange(0, num_qrtr_d)
da1d_out = xr.DataArray(var1_mtrx_d, coords=[sttn_wnt, var1_d_lbl, indx_time_d], dims=['station', 'var1', 'index_time'])
da2d_out = xr.DataArray(var2_mtrx_d, coords=[sttn_wnt, var2_d_lbl, dpth1_int, indx_time_d], dims=['station', 'var2', 'depth', 'index_time'])
da3d_out = xr.DataArray(yr_vec_d, coords=[sttn_wnt, indx_time_d], dims=['station', 'index_time'])
da4d_out = xr.DataArray(jd_vec_d, coords=[sttn_wnt, indx_time_d], dims=['station', 'index_time'])


ds1d_out = da1d_out.to_dataset(name='var1_mtrx_d')
ds1d_out['var2_mtrx_d'] = da2d_out
ds1d_out['yr_d'] = da3d_out
ds1d_out['jd_d'] = da4d_out

# -------------------------------------------------------
# --create the cdt cast data that will be append to the bottle
# -------------------------------------------------------
# -- these are used to align with the variables in the data bottle base
# -- these are self curated as some sensors might have errors or not report data
var2_b_wnt0 = ['TempAve', 'Salt2_Corr', 'SigThetaTS1', 'DynHt',
               'Ox1_StaCorr', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']

var2_b_wnt1 = ['Temp1', 'Salt1_Corr', 'SigThetaTS1', 'DynHt',
               'Ox1_StaCorr', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']
var2_b_wnt2 = ['Temp2', 'Salt2_Corr', 'SigThetaTS2', 'DynHt',
               'Ox1_StaCorr', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']

var2_b_wnt3 = ['TempAve', 'Salt2', 'SigThetaTS1', 'DynHt',
               'Ox1', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']

var2_b_wnt4 = ['TempAve', 'Salt2_Corr', 'SigThetaTS1', 'DynHt',
               'Ox1_StaCorr', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']

var2_b_wnt5 = ['TempAve', 'Salt2_Corr', 'SigThetaTS1', 'DynHt',
               'Ox1_StaCorr', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']

var2_b_wnt6 = ['TempAve', 'Salt2_Corr', 'SigThetaTS1', 'DynHt',
               'Ox1_CruiseCorr', 'OxSat1', 'SIL', 'PO4', 'NO3', 'NO2',
               'Chl-a', 'Phaeo']
var2_b_list = [var2_b_wnt0, var2_b_wnt1, var2_b_wnt2, var2_b_wnt3, var2_b_wnt4, var2_b_wnt5, var2_b_wnt6]

num_var2_b = len(var2_b_wnt0)
var2_mtrx_b = np.zeros([num_sttn, num_var2_b, num_int2, num_qrtr_d])*np.nan
for ii in range(0, num_qrtr_d):
    # get the column variable names
    var2_b_wnt = var2_b_list[in_var2_b_list[ii]]
    
    for jj in range(0, num_sttn):
        da1jj = ds1d_out['var2_mtrx_d'][jj, :, :, ii]
        for ll in range(0, num_var2_b):
            # get the data wanted
            lbl_var_wnt = var2_b_wnt[ll]

            # Special exceptions to the oxygen data, mostly for "unimportant"
            # stations. This is very ugly code, but it was constructed
            # by inspecting each ctd data file using Libre office and identifying
            # oxygen columns that would work and inspecting plots of all 66
            # stations over all depths to find stations that had extreme DO values.
            if study_d[ii] == '2107SR' and df_s['Line'].values[jj] == 93.3 and df_s['Station'].values[jj] == 28.0 and ll == 4:
                lbl_var_wnt = 'Ox2_StaCorr'
            if study_d[ii] == '2107SR' and df_s['Line'].values[jj] == 93.3 and df_s['Station'].values[jj] == 30.0 and ll == 4:
                lbl_var_wnt = 'Ox2_StaCorr'
            if study_d[ii] == '2208BH' and df_s['Line'].values[jj] == 93.3 and df_s['Station'].values[jj] == 55.0 and ll == 4:
                lbl_var_wnt = 'Ox2_StaCorr'
            if study_d[ii] == '2301RL' and df_s['Line'].values[jj] == 93.3 and df_s['Station'].values[jj] == 28.0 and ll == 4:
                lbl_var_wnt = 'Ox2'
            if study_d[ii] == '2307SR' and df_s['Line'].values[jj] == 93.3 and df_s['Station'].values[jj] == 28.0 and ll == 4:
                lbl_var_wnt = 'Ox2'
            if study_d[ii] == '2307SR' and df_s['Line'].values[jj] == 90.0 and df_s['Station'].values[jj] == 110.0 and ll == 4:
                lbl_var_wnt = 'Ox2'
            if study_d[ii] == '2311SR' and df_s['Line'].values[jj] == 83.3 and df_s['Station'].values[jj] == 55.0 and ll == 4:
                lbl_var_wnt = 'Ox2'
            if study_d[ii] == '2311SR' and df_s['Line'].values[jj] == 90.0 and df_s['Station'].values[jj] == 120.0 and ll == 4:
                lbl_var_wnt = 'Ox2'
            if study_d[ii] == '2404SH' and df_s['Line'].values[jj] == 90.0 and df_s['Station'].values[jj] == 53.0 and ll == 4:
                lbl_var_wnt = 'Ox1'
            if study_d[ii] == '2408SR' and df_s['Line'].values[jj] == 93.3 and df_s['Station'].values[jj] == 30.0 and ll == 4:
                lbl_var_wnt = 'Ox1'

            da1ll = da1jj.sel(var2=lbl_var_wnt)

            # interpolate
            dpth_ll = da1ll.depth.data
            var_ll = da1ll.data

            # --masked array for missing values
            var_ll_m = ma.array(var_ll, mask=np.isnan(var_ll))
            dpth_ll_m = ma.array(dpth_ll, mask=np.isnan(var_ll))

            if dpthd.compressed().size > 0:
                data_int = np.interp(dpth2_int, dpthd.compressed(
                ), vard.compressed(), left=np.nan, right=np.nan)
            else:
                data_int = np.zeros(num_int2)*np.nan

            var2_mtrx_b[jj, ll, :, ii] = data_int

# Create a netcdf file with the prelim data
da1b_out = xr.DataArray(var1_mtrx_d, coords=[sttn_wnt, var1_d_lbl, indx_time_d], dims=['station', 'var1', 'index_time'])
da2b_out = xr.DataArray(var2_mtrx_b, coords=[sttn_wnt, var2_b_lbl, dpth2_int, indx_time_d], dims=['station', 'var2', 'depth', 'index_time'])
da3b_out = xr.DataArray(yr_vec_d, coords=[sttn_wnt, indx_time_d], dims=['station', 'index_time'])
da4b_out = xr.DataArray(jd_vec_d, coords=[sttn_wnt, indx_time_d], dims=['station', 'index_time'])


ds1b_out = da1b_out.to_dataset(name='var1_mtrx_d')
ds1b_out['var2_mtrx_d'] = da2b_out
ds1b_out['yr_d'] = da3b_out
ds1b_out['jd_d'] = da4b_out

# -------------------------------------------------------
# --create output directory and save xarray dataset to netcdf file
# -------------------------------------------------------

# --check if directory exist, if it doesn't then create
try:
    os.makedirs(dir_out)
except OSError:
    if not os.path.isdir(dir_out):
        raise

# --Save Dataset to a netcdf file
dir1_fn = dir_out + '/' + fn1_out
dir2_fn = dir_out + '/' + fn2_out
ds1d_out.to_netcdf(dir1_fn)
ds1b_out.to_netcdf(dir2_fn)
