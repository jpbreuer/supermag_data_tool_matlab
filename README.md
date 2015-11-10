# supermag_import_data_matlab
Small fuction to convert supermag data into a decent usable matlab format.

This code takes magnetometer datasets provided from supermag.jhuapl.edu/ reads them and stores the variables as convenient matrixes in matlab.

supermagdata2matrix deletes all of the stations with data gaps for that day and returns all of the variables as matrixes
supermagdataall2matrix does not delete anything but returns each station as a cell for each variable. Difference is that accessing points from a cell requires slightly different indexing.
