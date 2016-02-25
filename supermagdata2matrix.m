function [Date_UTC, jd2000, cellNos, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = data2matrix(filepath)
% Import Data
olddata = readtable(filepath);
olddata = sortrows(olddata,2);

% Delete Incomplete Station Data
min_day = 60*24;
IAGA = olddata(:,2);
[cellNos cellStartInd enumCells ] = unique(IAGA);
cellEndInd = cellStartInd(2:end) -1;
cellEndInd(end+1) = length(enumCells);

stnstart_diff = zeros(length(cellStartInd));
for ii = 2:length(cellStartInd)
    stnstart_diff(ii-1) = cellStartInd(ii) - cellStartInd(ii-1);
end
stnstart_diff = stnstart_diff(:,1);
stnstart_diff = abs(stnstart_diff - min_day);
stnstart_diffindex = find(stnstart_diff > 1);

indexrange = [];
for ii = 2:length(stnstart_diffindex)
    startindex(ii-1) = cellStartInd(stnstart_diffindex(ii-1));
    endindex(ii-1) = cellStartInd(stnstart_diffindex(ii-1)+1);
    indexrange = [indexrange startindex(ii-1):endindex(ii-1)-1];
end
lastindexrange = [endindex(end):cellEndInd(end)];
indexrange = [indexrange lastindexrange];
[data PS] = removerows(olddata,indexrange);
 
IAGA = data(:,2);
[cellNos cellStartInd enumCells ] = unique(IAGA);
cellEndInd = cellStartInd(2:end) -1;
cellEndInd(end+1) = length(enumCells);
[~,sortIdx] = sort(enumCells);
% IAGA = reshape(IAGA(sortIdx),min_day,length(cellNos));

% Time Data
Date_UTC = data(:,1);
Date_UTC = table2array(Date_UTC);

% idealsize = 1440*height(cellNos);
% idealStartrange = 1:1440:idealsize;
% idealStartrange = idealStartrange';
% idealtmatrix = zeros(length(idealsize));
% acdiff = idealStartrange - cellStartInd;

for d = 1:min_day
    time(d) = Date_UTC(d);
end
timevec = datetime_JP(time);
jd2000 = julian_JP(timevec);

% Import other Data
MLT = data(:,3);
MLAT = data(:,4);
IGRF_DECL = data(:,5);
SZA = data(:,6);
N = data(:,7);
E = data(:,8);
Z = data(:,9);

MLT = table2array(MLT);
MLAT = table2array(MLAT);
IGRF_DECL = table2array(IGRF_DECL);
SZA = table2array(SZA);
N = table2array(N);
E = table2array(E);
Z = table2array(Z);

% Replace 999999 with NaN
MLT(MLT == 999999) = NaN;
MLAT(MLAT == 999999) = NaN;
IGRF_DECL(IGRF_DECL == 999999) = NaN;
SZA(SZA == 999999) = NaN;
N(N == 999999) = NaN;
E(E == 999999) = NaN;
Z(Z == 999999) = NaN;

% Reshape data
MLT = reshape(MLT,min_day,height(cellNos));
MLAT = reshape(MLAT,min_day,height(cellNos));
IGRF_DECL = reshape(IGRF_DECL,min_day,height(cellNos));
SZA = reshape(SZA,min_day,height(cellNos));
N = reshape(N,min_day,height(cellNos));
E = reshape(E,min_day,height(cellNos));
Z = reshape(Z,min_day,height(cellNos));
end