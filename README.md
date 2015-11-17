# supermag_import_data_matlab
Small fuction to convert supermag csv data into a decent usable matlab format.

This code takes magnetometer csv datasets provided from supermag.jhuapl.edu/, reads them, and stores the variables as convenient matrixes in matlab.

The supermagdata2matrix function deletes all of the stations with data gaps for that day and returns all of the variables as double matrixes. This is NOT IDEAL as it obviously deletes a lot of data.

The supermagdataall2matrix function does not delete any data, replaces data gaps with NaN's, and returns double arrays of each variable as (1440 by No. of stations), making this function extremely convenient. Column index is equal to station number as given by variable 'IAGA', and specific station data is given in the stationdata.m file.
