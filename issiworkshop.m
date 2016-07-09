clear all;close all;clear;clc;

%% Download dataset locally
username = 'jp.breuer'; % must be SUPERMAG username
startTime = '2012-09-30T00:00:00.000Z';
endTime = '2014-04-01T00:00:00.000Z';
% stations_all = 'ALE,THL,SVS,RES,KUV,NRD,UPN,CY0,IGC,TAL,GHC,DMH,CBB,UMQ,RPB,NAL,GDH,LYR,DNB,PGC,ATU,CHC,HRN,CDC,BLC,VIZ,CNL,HOP,STF,RAN,IQA,SKT,MCR,SCO,BJN,INK,EKP,T35,KAV,DED,GHB,BRW,JAN,C07,YKC,AMK,FCC,DIK,FHB,NOR,SMI,T29,FSP,FYU,RAL,SOR,T31,GIM,TRO,EAG,AND,KEV,DAW,MAS,NAQ,TIK,PBQ,KIL,PKR,PBK,T41,AMD,CGO,CMO,LEK,ABK,CHD,LRV,MUO,KIR,FMC,ISL,HLL,LOZ,T42,T38,SOD,T33,PEL,JCK,GAK,DON,T39,T36,MEA,RVK,T40,OUL,S01,LYC,T28,C06,OUJ,PIN,FAR,T32,T22,T30,SIT,RED,T37,DOB,C12,C04,HAN,C11,SOL,T43,T03,LER,T15,M01,LET,NUR,C08,UPS,KAR,OTT,T24,KVI,LOV,CLK,T17,C10,NEW,GML,T23,YAK,TAR,BOR,BOX,MGD,VIC,GTF,STJ,C01,SHU,ESK,MSH,T21,ARS,BRZ,BFE,M04,ROE,KNZ,T25,YOR,NVS,HLP,WNG,SUW,VAL,BOU,FRD,SZC,M05,NGK,HAD,BEL,DSO,IRT,MAB,KIV,DOU,PET,ZAG,LVV,M06,T16,BDV,VYH,BFO,CLF,FUR,HRB,NCK,PAG,FRN,THY,ODE,CST,BSL,KHB,SUA,P01,GCK,TUC,M08,DLR,WMQ,RNC,AAA,CNH,MSR,MMB,AQU,RIK,ISK,DUR,IZN,BMT,EBR,SPT,MIZ,PEG,LZH,CYG,E05,KAK,SJG,A02,SFS,ONW,HTY,MID,KAG,KNY,HON,JAI,CBI,LNP,GZH,GUI,PHU,ABG,HYB,MUT,KOU,TAM,A12,GUA,A04,DLT,AAE,HUA,MBO,A13,A07,A01,A06,A03,BNG,VRE,GAN,KTB,BIK,ASC,API,PPT,SER,VSS,PIL,IPM,EWA,KDU,CKI,SHE,OSO,PAC,TAN,CTA,TRW,A05,TSU,LRM,ASP,HBK,KMH,PST,TDC,HER,GNG,GNA,CNB,ESC,ORC,LIV,OHI,AMS,PAL,EYR,AIA,B03,CZT,PAF,B11,B12,MCQ,DMC,B14,B07,B16,B18,MAW,PG4,B10,B19,B23,PG3,B21,SPA,PG2,B22,PG1,B20,SBA,MCM,DRV,CSY,VOS';
stations = 'DOU,NUR,SOD,NGK,WNG,FUR,BFO,VAL,ESK,HAD,LER,CLF,MAB,ABG,FRD,PG3,PG1,PG2,';

% path_to_download = './supermag_data/csv/';
path_to_download = '~/Work/Internships/NASA_Goddard/GCR_Neutron/SUPERMAG_data/';

supermag_data_download(username,startTime,endTime,stations,path_to_download)

%% Specify dataset
% filepath = './supermag_data/csv/2014-Jan-01.csv';
list_temp = dir(path_to_download);% '*.csv']);
list = {list_temp.name}';
list{1} = {};list{2} = {};list{3} = {};list = list(~cellfun(@isempty,list));
for ii = 1:length(list)
    filepath{ii} = sprintf([path_to_download list{ii}]);
end
%% Load data

IAGA_all = {'SPA','B23','B21','B22','B19','B20','PG1','B18','PG2','PG3','B17','B16','PG4','B14','B15','B27','B13','SBA','VOS','B12','MCM','B11','B10','DMC','B09','B08','B07','B24','NVL','VNA','B06','DVS','PRG','DRV','MAW','B04','B05','B03','CSY','MIR','AIA','B02','PAL','B01','LIV','OHI','ESC','ORC','KEP','MCQ','PNT','PST','ENP','PAF','CZT','EYR','TRW','LEM','VLD','OSO','PAC','TDC','AMS','CNB','CAN','HER','KAT','ADL','CER','GNA','SUT','PIL','GNG','SER','IPM','KMH','DAL','HBK','BSV','LMM','ASP','VSS','ANT','LRM','EWA','CTA','TAN','TSU','PUT','A05','VRE','PPT','SHE','NMP','API','ASA','DRW','WEP','KDU','HUA','CKI','A10','WTK','ASC','WEW','BIK','A11','PTN','KTB','KOU','TND','GAN','BNG','A03','A06','KOR','A08','A13','AAE','YAP','CRP','A01','A07','DLT','GUA','MBO','MUT','A04','A09','A12','PNL','HYB','ABG','SJG','TEO','HON','PHU','M11','TAM','GZH','SON','LNP','M10','CBI','M09','JAI','GUI','MID','FIT','CDP','BSL','DLR','JAX','ELT','MLT','M08','TUC','KAG','YMK','KNY','BGY','ONW','HTY','QSB','USC','M07','M06','T26','T27','LZH','FRN','SMA','KAK','TUL','DSO','SFS','E05','CYG','A02','FRD','ASH','E04','PEG','BOU','APL','ESA','MIZ','M05','T16','BMT','SPT','TOL','BJI','ISK','TKT','EBR','T20','E02','E03','M04','IZN','E01','AQU','DUR','C01','MMB','AAA','PPI','RIK','MSH','GTF','M03','OTT','SUA','MSR','CLK','GCK','SBL','T21','T23','NKK','ODE','M02','T24','CNH','WMQ','RNC','STJ','THY','M01','C08','T17','T25','T15','T18','CST','P01','VIC','NEW','CLF','FUR','HRB','NCK','YSS','BDV','VLO','BFO','C10','C11','T19','PAG','KHB','VYH','WIC','MAB','DOU','PIN','MZH','GLN','LET','T50','T51','KIV','KGD','LVV','VAL','HAD','BEL','ROT','T03','C04','C12','T30','T32','BRD','ZAG','WHS','T49','NGK','IRT','RED','C13','T43','C05','MEA','WNG','SL','LAN','YOR','EDM','SAS','MSK','C06','T28','PET','SZC','PBQ','ESK','HLP','MNK','BFE','ROE','NVS','RSV','MOS','C09','T33','T36','T37','SUW','T52','T42','T45','T48','SIT','GIM','FMC','NAN','FSJ','KNZ','SHU','T31','ARS','BRZ','BOX','BOR','CRK','GML','LOV','FCC','RAL','FVE','C02','T22','LER','SMI','KAR','TAR','LNN','YAK','MGD','KVI','HOM','C03','T29','T44','AMU','NAQ','EKP','NUR','UPS','T38','GRK','T46','T53','YKC','FSP','MEK','FHB','RAN','DOB','SOL','HAN','GAK','FAR','HLM','TLK','TRP','T39','T47','BLC','LRV','GHB','DAW','IQA','CHC','HLL','S01','T35','T40','SKT','AMK','RVK','LYC','OUJ','EAG','CMO','CGO','PKR','ARK','CDC','OUL','C07','T34','MCR','CNL','JCK','DON','ZYK','KOT','BET','FYU','CWE','ZGN','PGC','RPB','ATU','STF','SOD','PEL','T41','CBB','ARC','INK','GDH','AND','KAU','IVA','ABK','LEK','MUO','LOZ','KIR','CPS','CKA','LOP','GHC','IGC','PBC','DED','NOK','UMQ','SCO','TAL','KAV','NOR','JAN','SOR','TRO','ALT','KEV','MAS','KIL','CPY','AMD','BRW','JCO','TIK','CHD','PBK','CY0','UPN','MCE','MCW','MCG','SAH','DIK','RES','KUV','DNB','MCN','BJN','TAB','SVS','KTN','MBC','THL','DMH','HOP','HRN','CCS','HOR','NAL','LYR','BBG','VIZ','HIS','EUA','ALE','NRD'};
IAGA_all = sort(IAGA_all);

% jd2000_vec = [];MLT_vec = [];MLAT_vec = [];IGRF_DECL_vec = [];SZA_vec = [];N_vec = [];E_vec = [];Z_vec = [];B_vec = [];
% jd2000 = [];IAGA = {};MLT = [];MLAT = [];IGRF_DECL = [];SZA = [];N = [];E = [];Z = [];B = [];
MLT_matrix = NaN(1440*length(list),length(IAGA_all));jd2000_vec = [];MLAT_matrix = MLT_matrix;IGRF_DECL_matrix = MLT_matrix;SZA_matrix = MLT_matrix;N_matrix = MLT_matrix;E_matrix = MLT_matrix;Z_matrix = MLT_matrix;B_matrix = MLT_matrix;

iterIndex = 1;
for ii = 1:length(list)
    [jd2000, IAGA, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = supermagdataall2matrix(filepath{ii});
    B = sqrt(N.^2 + E.^2 + Z.^2);
    
    [~, matchIndex] = ismember(IAGA,IAGA_all);
    
%     jd2000_vec = [jd2000_vec; jd2000];
%     MLT_vec = [MLT_vec; MLT];
%     MLAT_vec = [MLAT_vec; MLAT];
%     IGRF_DECL_vec = [IGRF_DECL_vec; IGRF_DECL];
%     SZA_vec = [SZA_vec; SZA];
%     N_vec = [N_vec; N];
%     E_vec = [E_vec; E];
%     Z_vec = [Z_vec; Z];
%     B_vec = [B_vec; B];
    jd2000_vec = [jd2000_vec;jd2000];    
    MLT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLT;
    MLAT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLAT;
    IGRF_DECL_matrix(iterIndex:length(MLT)*ii,matchIndex) = IGRF_DECL;
    SZA_matrix(iterIndex:length(MLT)*ii,matchIndex) = SZA;
    N_matrix(iterIndex:length(MLT)*ii,matchIndex) = N;
    E_matrix(iterIndex:length(MLT)*ii,matchIndex) = E;
    Z_matrix(iterIndex:length(MLT)*ii,matchIndex) = Z;
    B_matrix(iterIndex:length(MLT)*ii,matchIndex) = B;
    
    iterIndex = (ii * 1440)+1;
end
% IAGA variable contains relative indexed information for all available
% stations (alphabetized). Each station index corresponds to a column in any of the data
% vectors.

%% Create Table

% for ii = 1:length(IAGA)
%     MLT_matrix(:,matchIndex(ii)) = MLT(:,ii);
%     MLAT_matrix(:,matchIndex(ii)) = MLAT(:,ii);
%     IGRF_DECL_matrix(:,matchIndex(ii)) = IGRF_DECL(:,ii);
%     SZA_matrix(:,matchIndex(ii)) = SZA(:,ii);
%     N_matrix(:,matchIndex(ii)) = N(:,ii);
%     E_matrix(:,matchIndex(ii)) = E(:,ii);
%     Z_matrix(:,matchIndex(ii)) = Z(:,ii);
%     B_matrix(:,matchIndex(ii)) = B(:,ii);
% end

% T = table(jd2000,'VariableNames','A13');%IAGA_all);