# supermag_data_tool_matlab
Tool to effectively download, import, and use Supermag Magnetometer data.

Start with issiworkshop.m as it calls all the other functions.

The supermag_data_download.m downloads one csv file per day within the specified range (startTime,endTime) to a directory.

The supermag_combine_data.m function organises all of the data for each imported file using the supermagdataall2matrix.m function, taking into account the station and time index. This creates a time series to work with, and outputs all of the new data to a txt file where the first column is jd2000. Once imported, all gaps contain NaNs and work as regular double arrays.

The supermagdataall2matrix.m function does not delete any data, replaces all data gaps with NaN's, and returns double arrays of each variable as (1440 by No. of stations), making this function extremely convenient. Column index is equal to station number as given by variable 'IAGA', and specific station data is given in the stationdata.m file.

Output to txt files is to handle LARGE time ranges (tested with 548 days of continues data for 18 different stations at one minute resolution). Note, large datasets may result in issues with memory or will take a long time to run..
