clear all

file = './supermag-stations.csv';
data = csvimport(file);

stncode = data(:,1);               %IAGA
stngeolon = data(:,2);             %GLON
stngeolat = data(:,3);             %GLAT
stnmaglon = data(:,4);             %MLON
stnmaglat = data(:,5);             %MLAT
stnstation_name = data(:,6);       %STATIONNAME
stnoperator_number = data(:,7);    %OPERATORNUM
stnoperators = data(:,8);          %OPERATORS