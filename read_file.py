import os
import xarray as xr
import numpy as np


def find_file(varname, path):
    result = []
    for file in os.listdir(path):
        dumy_name = str("." + varname + ".") # make sure to choose the specific var
        if "CESM2-WACCM-SSP245" in path:
            if dumy_name in file and "2015" in file:
                   result.append(path + file)
        else:
            if dumy_name in file:
                result.append(path + file)
    
    if len(result) > 1:
         print("ERR: multiple files containing "+varname + ":")
         print(result)
    return result[0]


def read_file(filename, varname, time_span, anom=True):
    DS = xr.open_dataset(filename)
    #change time steps 
    DS = DS.assign_coords(time = DS.time_bnds[:,0])
    v = getattr(DS, varname)
    v = v.sel(time = slice(time_span[0], time_span[1]))
    return v

def read_data(varname, PATH, CASE, time_span):
    for i in range(10):
        print("read " + varname + " in case:" + str(i+1))
        if i == 9:
            case = PATH + CASE + "010/atm/proc/tseries/month_1/"
        else:
            case = PATH + CASE + "00" + str(i+1) + "/atm/proc/tseries/month_1/"

        fname = find_file(varname, case)
        
        temp = read_file(fname, varname, time_span)
        if i == 0:
            sst = temp
        else:
            sst = xr.concat([sst, temp], "ens")
    
    return sst

#-----------------------------------------------
# read sulfate injection amount from log file
#-----------------------------------------------
def read_sai(PATH):
    for i in range(10):
        
        # find file
        if i == 9:
            path = PATH + "010/controller/"
        else:
            path = PATH + "00" + str(i+1) + "/controller/"
        
        for file in os.listdir(path):
            dumy_name = "ControlLog_"
            if dumy_name in file:
                fname = file
                break;
        print("read SAI in case" + str(i+1) + ': ' + fname)
          
        # read ascii file
        year = np.loadtxt(path + fname, skiprows = 1)[:, 0]
        
        t0_tmp = xr.DataArray(np.loadtxt(path + fname, skiprows = 1)[:, 1],
                          dims = 'year',
                          coords = {'year':year},
                          attrs=dict(units = 'K',
                          )
                         )
        
        t1_tmp = xr.DataArray(np.loadtxt(path + fname, skiprows = 1)[:, 3],
                  dims = 'year',
                  coords = {'year':year},
                  attrs=dict(units = 'K',
                  )
                 )
        
        t2_tmp = xr.DataArray(np.loadtxt(path + fname, skiprows = 1)[:, 5],
                  dims = 'year',
                  coords = {'year':year},
                  attrs=dict(units = 'K',
                  )
                 )
        
        s_tmp = xr.DataArray(np.loadtxt(path + fname, skiprows = 1)[:, -4:],
                         dims = ['year', 'loc'],
                         coords = {'year':year, 'loc':['30S','15S','15N','30N']},
                         attrs = dict(untis = 'Tg')
                        )
        if i == 0:
            t0 = t0_tmp
            t1 = t1_tmp
            t2 = t2_tmp
            s = s_tmp
        else:
            t0 = xr.concat([t0, t0_tmp], "ens")
            t1 = xr.concat([t1, t1_tmp], "ens")
            t2 = xr.concat([t2, t2_tmp], "ens")
            s = xr.concat([s, s_tmp], "ens")
    
    return t0, t1, t2, s




def read_ssp(varname, time_span):
    path = '/glade/scratch/dchenrui/CESM2-WACCM-SSP245/' + varname + '/'
    
    for i in range(10):
        print("read " + varname + " in case:" + str(i+1))
        
        if i == 9: 
            fname = path + 'b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.010.cam.h0.TREFHT.201501-206912.nc'
        else:
            fname = path + 'b.e21.BWSSP245cmip6.f09_g17.CMIP6-SSP2-4.5-WACCM.00'+str(i+1)+'.cam.h0.TREFHT.201501-206912.nc'
        
        temp = read_file(fname, varname, time_span)
        
        if i ==0:
            var = temp
        else:
            var = xr.concat([var, temp], 'ens')
    return var