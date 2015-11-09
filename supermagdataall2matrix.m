function [jd2000, IAGA, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = supermagdataall2matrix(filepath)
% Import Data
olddata = readtable(filepath);
olddata = sortrows(olddata,2);

% Time Data
min_day = 60*24;

Date_UTC = olddata(:,1);
IAGA = olddata(:,2);
MLT = olddata(:,3);
MLAT = olddata(:,4);
IGRF_DECL = olddata(:,5);
SZA = olddata(:,6);
N = olddata(:,7);
E = olddata(:,8);
Z = olddata(:,9);

Date_UTC = table2array(Date_UTC);
IAGA = table2array(IAGA);
MLT = table2array(MLT);
MLAT = table2array(MLAT);
IGRF_DECL = table2array(IGRF_DECL);
SZA = table2array(SZA);
N = table2array(N);
E = table2array(E);
Z = table2array(Z);

MLT(MLT == 999999) = NaN;
MLAT(MLAT == 999999) = NaN;
IGRF_DECL(IGRF_DECL == 999999) = NaN;
SZA(SZA == 999999) = NaN;
N(N == 999999) = NaN;
E(E == 999999) = NaN;
Z(Z == 999999) = NaN;

[IAGA,~,index] = unique(IAGA);
MLT = accumarray(index(:),MLT,[],@(x) {x});
MLAT = accumarray(index(:),MLAT,[],@(x) {x});
IGRF_DECL = accumarray(index(:),IGRF_DECL,[],@(x) {x});
SZA = accumarray(index(:),SZA,[],@(x) {x});
N = accumarray(index(:),N,[],@(x) {x});
E = accumarray(index(:),E,[],@(x) {x});
Z = accumarray(index(:),Z,[],@(x) {x});


[Date_UTC,~,indexUTC] = unique(Date_UTC);
% Date_UTC = accumarray(indexUTC(:),Date_UTC,[],@(x) {x});
for d = 1:min_day
    jd2000(d) = julian(datetime(Date_UTC(d)));
end
mapObj = containers.Map(Date_UTC,jd2000);

Date_UTC = olddata(:,1);
Date_UTC = table2array(Date_UTC);
jd2000 = values(mapObj,Date_UTC);
jd2000 = table2array(jd2000);
jd2000 = accumarray(index(:),jd2000,[],@(x) {x});
end

