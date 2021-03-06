function [jd2000, IAGA, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = supermagdataall2matrix(filepath)
% clear all
% Import Data
% filepath = '~/Work/Internships/NASA_Goddard/GCR_Neutron/SUPERMAG_data/supermag_2012-09-30.csv'%'./supermag_data/csv/2014-Jan-01.csv';

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
%plot(1:1440,N(10081:11520));hold on;plot(1:1440,E(10081:11520));plot(1:1440,Z(10081:11520));

[IAGA,~,index] = unique(IAGA);
MLT = accumarray(index(:),MLT,[],@(x) {x});
MLAT = accumarray(index(:),MLAT,[],@(x) {x});
IGRF_DECL = accumarray(index(:),IGRF_DECL,[],@(x) {x});
SZA = accumarray(index(:),SZA,[],@(x) {x});
N = accumarray(index(:),N,[],@(x) {x});
E = accumarray(index(:),E,[],@(x) {x});
Z = accumarray(index(:),Z,[],@(x) {x});


[Date_UTC,~,indexUTC] = unique(Date_UTC);

datevec = datetime_JP(Date_UTC);
decimal_min = datevec(:,5)./60;

jd2000exact = jd2000_new(datevec(:,1),datevec(:,2),datevec(:,3),(datevec(:,4)+decimal_min));

mapObj = containers.Map(Date_UTC,jd2000exact);

Date_UTC = olddata(:,1);
Date_UTC = table2array(Date_UTC);

jd2000 = values(mapObj,Date_UTC);
jd2000 = cell2mat(jd2000);
jd2000 = accumarray(index(:),jd2000,[],@(x) {x});

% fill gaps with NaN
tdiffindex = [];
for d = 1:length(jd2000);
    if length(jd2000{d,1}) == 1440
        tdiffindex(d) = 0;
    else
        tdiffindex(d) = 1;
    end
end
tdiffindex = find(tdiffindex > 0);

% diff = jd2000{1,1}(2) - jd2000{1,1}(1);

% jd2000no = [];
jd = jd2000;
% datenew = Date_UTC;
% iaganew = IAGA;
mltnew = MLT;
mlatnew = MLAT;
igrfnew = IGRF_DECL;
szanew = SZA;
nnew = N;
enew = E;
znew = Z;


for ii = 1:length(tdiffindex)
%     jd2000no = [jd2000no jd2000{tdiffindex(ii),1}'];
    [tf, loc] = ismember(jd2000{tdiffindex(ii),1}, jd2000exact);
     jd{tdiffindex(ii),1} = nan(size(jd2000exact));
    
% define empty nan arrays
%     datenew{tdiffindex(ii),1} = nan(size(jd2000exact));
%     iaganew{tdiffindex(ii),1} = nan(size(jd2000exact));
    mltnew{tdiffindex(ii),1} = nan(size(jd2000exact));
    mlatnew{tdiffindex(ii),1} = nan(size(jd2000exact));
    igrfnew{tdiffindex(ii),1} = nan(size(jd2000exact));
    szanew{tdiffindex(ii),1} = nan(size(jd2000exact));
    nnew{tdiffindex(ii),1} = nan(size(jd2000exact));
    enew{tdiffindex(ii),1} = nan(size(jd2000exact));
    znew{tdiffindex(ii),1} = nan(size(jd2000exact));
    
% replace values in correct location leaving nans
     jd{tdiffindex(ii),1}(loc) = jd2000{tdiffindex(ii),1};
%     datenew{tdiffindex(ii),1}(loc) = Date_UTC{tdiffindex(ii),1};
%     iaganew{tdiffindex(ii),1}(loc) = IAGA{tdiffindex(ii),1};
    mltnew{tdiffindex(ii),1}(loc) = MLT{tdiffindex(ii),1};
    mlatnew{tdiffindex(ii),1}(loc) = MLAT{tdiffindex(ii),1};
    igrfnew{tdiffindex(ii),1}(loc) = IGRF_DECL{tdiffindex(ii),1};
    szanew{tdiffindex(ii),1}(loc) = SZA{tdiffindex(ii),1};
    nnew{tdiffindex(ii),1}(loc) = N{tdiffindex(ii),1};
    enew{tdiffindex(ii),1}(loc) = E{tdiffindex(ii),1};
    znew{tdiffindex(ii),1}(loc) = Z{tdiffindex(ii),1};
end

jd2000 = cell2mat(jd);
MLT = cell2mat(mltnew);
MLAT = cell2mat(mlatnew);
IGRF_DECL = cell2mat(igrfnew);
SZA = cell2mat(szanew);
N = cell2mat(nnew);
E = cell2mat(enew);
Z = cell2mat(znew);
% 
MLT(MLT == 999999) = NaN;
MLAT(MLAT == 999999) = NaN;
IGRF_DECL(IGRF_DECL == 999999) = NaN;
SZA(SZA == 999999) = NaN;
N(N == 999999) = NaN;
E(E == 999999) = NaN;
Z(Z == 999999) = NaN;
% 
[IAGA,~,index] = unique(IAGA);

jd2000 = reshape(jd2000,min_day,length(IAGA));jd2000 = jd2000(:,1);
MLT = reshape(MLT,min_day,length(IAGA));
MLAT = reshape(MLAT,min_day,length(IAGA));
IGRF_DECL = reshape(IGRF_DECL,min_day,length(IAGA));
SZA = reshape(SZA,min_day,length(IAGA));
N = reshape(N,min_day,length(IAGA));
E = reshape(E,min_day,length(IAGA));
Z = reshape(Z,min_day,length(IAGA));
% 
% MLT_dat = dataset(MLT);MLT_dat.Properties.VarNames = IAGA;

end


% MLT = accumarray(index(:),MLT,[],@(x) {x});
% MLAT = accumarray(index(:),MLAT,[],@(x) {x});
% IGRF_DECL = accumarray(index(:),IGRF_DECL,[],@(x) {x});
% SZA = accumarray(index(:),SZA,[],@(x) {x});
% N = accumarray(index(:),N,[],@(x) {x});
% E = accumarray(index(:),E,[],@(x) {x});
% Z = accumarray(index(:),Z,[],@(x) {x});