function supermag_combine_data(startTime,endTime,stations,path_to_files,path_to_output_files)

STATIONS = sort(strsplit(stations,','));

list_temp = dir(path_to_files);% '*.csv']);
list = {list_temp.name}';
list{1} = {};list{2} = {};list{3} = {};list = list(~cellfun(@isempty,list));
days_of_data = length(list);
for ii = 1:length(list)
    filepath{ii} = sprintf([path_to_files list{ii}]);
end

IAGA_all = {'SPA','B23','B21','B22','B19','B20','PG1','B18','PG2','PG3','B17','B16','PG4','B14','B15','B27','B13','SBA','VOS','B12','MCM','B11','B10','DMC','B09','B08','B07','B24','NVL','VNA','B06','DVS','PRG','DRV','MAW','B04','B05','B03','CSY','MIR','AIA','B02','PAL','B01','LIV','OHI','ESC','ORC','KEP','MCQ','PNT','PST','ENP','PAF','CZT','EYR','TRW','LEM','VLD','OSO','PAC','TDC','AMS','CNB','CAN','HER','KAT','ADL','CER','GNA','SUT','PIL','GNG','SER','IPM','KMH','DAL','HBK','BSV','LMM','ASP','VSS','ANT','LRM','EWA','CTA','TAN','TSU','PUT','A05','VRE','PPT','SHE','NMP','API','ASA','DRW','WEP','KDU','HUA','CKI','A10','WTK','ASC','WEW','BIK','A11','PTN','KTB','KOU','TND','GAN','BNG','A03','A06','KOR','A08','A13','AAE','YAP','CRP','A01','A07','DLT','GUA','MBO','MUT','A04','A09','A12','PNL','HYB','ABG','SJG','TEO','HON','PHU','M11','TAM','GZH','SON','LNP','M10','CBI','M09','JAI','GUI','MID','FIT','CDP','BSL','DLR','JAX','ELT','MLT','M08','TUC','KAG','YMK','KNY','BGY','ONW','HTY','QSB','USC','M07','M06','T26','T27','LZH','FRN','SMA','KAK','TUL','DSO','SFS','E05','CYG','A02','FRD','ASH','E04','PEG','BOU','APL','ESA','MIZ','M05','T16','BMT','SPT','TOL','BJI','ISK','TKT','EBR','T20','E02','E03','M04','IZN','E01','AQU','DUR','C01','MMB','AAA','PPI','RIK','MSH','GTF','M03','OTT','SUA','MSR','CLK','GCK','SBL','T21','T23','NKK','ODE','M02','T24','CNH','WMQ','RNC','STJ','THY','M01','C08','T17','T25','T15','T18','CST','P01','VIC','NEW','CLF','FUR','HRB','NCK','YSS','BDV','VLO','BFO','C10','C11','T19','PAG','KHB','VYH','WIC','MAB','DOU','PIN','MZH','GLN','LET','T50','T51','KIV','KGD','LVV','VAL','HAD','BEL','ROT','T03','C04','C12','T30','T32','BRD','ZAG','WHS','T49','NGK','IRT','RED','C13','T43','C05','MEA','WNG','SL','LAN','YOR','EDM','SAS','MSK','C06','T28','PET','SZC','PBQ','ESK','HLP','MNK','BFE','ROE','NVS','RSV','MOS','C09','T33','T36','T37','SUW','T52','T42','T45','T48','SIT','GIM','FMC','NAN','FSJ','KNZ','SHU','T31','ARS','BRZ','BOX','BOR','CRK','GML','LOV','FCC','RAL','FVE','C02','T22','LER','SMI','KAR','TAR','LNN','YAK','MGD','KVI','HOM','C03','T29','T44','AMU','NAQ','EKP','NUR','UPS','T38','GRK','T46','T53','YKC','FSP','MEK','FHB','RAN','DOB','SOL','HAN','GAK','FAR','HLM','TLK','TRP','T39','T47','BLC','LRV','GHB','DAW','IQA','CHC','HLL','S01','T35','T40','SKT','AMK','RVK','LYC','OUJ','EAG','CMO','CGO','PKR','ARK','CDC','OUL','C07','T34','MCR','CNL','JCK','DON','ZYK','KOT','BET','FYU','CWE','ZGN','PGC','RPB','ATU','STF','SOD','PEL','T41','CBB','ARC','INK','GDH','AND','KAU','IVA','ABK','LEK','MUO','LOZ','KIR','CPS','CKA','LOP','GHC','IGC','PBC','DED','NOK','UMQ','SCO','TAL','KAV','NOR','JAN','SOR','TRO','ALT','KEV','MAS','KIL','CPY','AMD','BRW','JCO','TIK','CHD','PBK','CY0','UPN','MCE','MCW','MCG','SAH','DIK','RES','KUV','DNB','MCN','BJN','TAB','SVS','KTN','MBC','THL','DMH','HOP','HRN','CCS','HOR','NAL','LYR','BBG','VIZ','HIS','EUA','ALE','NRD'};
IAGA_all = sort(IAGA_all);

[~, columnIndex] = ismember(STATIONS,IAGA_all);
% jd2000_vec = [];MLT_vec = [];MLAT_vec = [];IGRF_DECL_vec = [];SZA_vec = [];N_vec = [];E_vec = [];Z_vec = [];B_vec = [];

% MLT_matrix = NaN(1440*days_of_data,length(IAGA_all));

%% MLT

if ~exist(sprintf([path_to_output_files '%s_%s_supermag_MLT.txt'],startTime(1:10),endTime(1:10)), 'file')
    jd2000_vec = [];
    % MLAT_matrix = MLT_matrix;
    % IGRF_DECL_matrix = MLT_matrix;
    % SZA_matrix = MLT_matrix;
    % N_matrix = MLT_matrix;
    % E_matrix = MLT_matrix;
    % Z_matrix = MLT_matrix;
    % B_matrix = MLT_matrix;
    
    
    iterIndex = 1;
    for ii = 1:length(list)
        [jd2000, IAGA, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = supermagdataall2matrix(filepath{ii});
        B = sqrt(N.^2 + E.^2 + Z.^2);
        [~, matchIndex] = ismember(IAGA,IAGA_all);
        %     if ii == 1
        %         stationIndex = matchIndex;
        %     end
        % %     [~, missingIndex]
        %     if length(IAGA) ~= length(stationIndex)
        %         missingIndex = setdiff(stationIndex,matchIndex);
        %         NanIndexstart = [find(missingIndex == stationIndex)];
        %         if length(NanIndexstart) > 1
        %                 error('must create method');break
        %         end
        %
        %             MLT = [MLT NaN(length(MLT),length(NanIndexstart))];
        %             MLT(:,NanIndexstart+1:end) = MLT(:,NanIndexstart:end-1);
        %             MLT(:,NanIndexstart) = NaN;
        %
        %             MLAT = [MLAT NaN(length(MLAT),length(NanIndexstart))];
        %             MLAT(:,NanIndexstart+1:end) = MLAT(:,NanIndexstart:end-1);
        %             MLAT(:,NanIndexstart) = NaN;
        %
        %             IGRF_DECL = [IGRF_DECL NaN(length(IGRF_DECL),length(NanIndexstart))];
        %             IGRF_DECL(:,NanIndexstart+1:end) = IGRF_DECL(:,NanIndexstart:end-1);
        %             IGRF_DECL(:,NanIndexstart) = NaN;
        %
        %             SZA = [SZA NaN(length(SZA),length(NanIndexstart))];
        %             SZA(:,NanIndexstart+1:end) = SZA(:,NanIndexstart:end-1);
        %             SZA(:,NanIndexstart) = NaN;
        %
        %             N = [N NaN(length(N),length(NanIndexstart))];
        %             N(:,NanIndexstart+1:end) = N(:,NanIndexstart:end-1);
        %             N(:,NanIndexstart) = NaN;
        %
        %             E = [E NaN(length(E),length(NanIndexstart))];
        %             E(:,NanIndexstart+1:end) = E(:,NanIndexstart:end-1);
        %             E(:,NanIndexstart) = NaN;
        %
        %             Z = [Z NaN(length(Z),length(NanIndexstart))];
        %             Z(:,NanIndexstart+1:end) = Z(:,NanIndexstart:end-1);
        %             Z(:,NanIndexstart) = NaN;
        %
        %             B = [B NaN(length(B),length(NanIndexstart))];
        %             B(:,NanIndexstart+1:end) = B(:,NanIndexstart:end-1);
        %             B(:,NanIndexstart) = NaN;
        %     end
        
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
        %     MLAT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLAT;
        %     IGRF_DECL_matrix(iterIndex:length(MLT)*ii,matchIndex) = IGRF_DECL;
        %     SZA_matrix(iterIndex:length(MLT)*ii,matchIndex) = SZA;
        %     N_matrix(iterIndex:length(MLT)*ii,matchIndex) = N;
        %     E_matrix(iterIndex:length(MLT)*ii,matchIndex) = E;
        %     Z_matrix(iterIndex:length(MLT)*ii,matchIndex) = Z;
        %     B_matrix(iterIndex:length(MLT)*ii,matchIndex) = B;
        
        iterIndex = (ii * 1440)+1;
    end
    
    MLT_new = MLT_matrix(:,columnIndex);
    MLT_new(find(MLT_new == 0)) = NaN;
    
    MLT_final = horzcat(jd2000_vec,MLT_new);
    
    % Save data
    filename_out = sprintf([path_to_output_files '%s_%s_supermag_MLT.txt'],startTime(1:10),endTime(1:10));
    fid_out = fopen(filename_out,'w');
    
    fprintf(fid_out,'%s' ,'Time (jd2000) ');
    for mm = 1:length(columnIndex)
        fprintf(fid_out,'%s ',IAGA_all{columnIndex(mm)});%,IAGA_all{columnIndex(2)},IAGA_all{columnIndex(3)},IAGA_all{columnIndex(4)},IAGA_all{columnIndex(5)},IAGA_all{columnIndex(6)},IAGA_all{columnIndex(7)},IAGA_all{columnIndex(8)},IAGA_all{columnIndex(9)},IAGA_all{columnIndex(10)},IAGA_all{columnIndex(11)},IAGA_all{columnIndex(12)},IAGA_all{columnIndex(13)},IAGA_all{columnIndex(14)},IAGA_all{columnIndex(15)},IAGA_all{columnIndex(16)},IAGA_all{columnIndex(17)},IAGA_all{columnIndex(18)});
    end
    fprintf(fid_out,'\n');
    [rows cols] = size(MLT_final);
    x = repmat('%10.5f\t',1,(cols-1));
    fprintf(fid_out,[x,'%10.5f\n'],MLT_final');

%     for mm = 1:length(columnIndex)
%         for kk = 1:length(MLT_new)
%             fprintf(fid_out,'%s ',jd2000_vec(kk),MLT_new(kk,mm));%,MLT_new(kk,2),MLT_new(kk,3),MLT_new(kk,4),MLT_new(kk,5),MLT_new(kk,6),MLT_new(kk,7),MLT_new(kk,8),MLT_new(kk,9),MLT_new(kk,10),MLT_new(kk,11),MLT_new(kk,12),MLT_new(kk,13),MLT_new(kk,14),MLT_new(kk,15),MLT_new(kk,16),MLT_new(kk,17),MLT_new(kk,18));
%         end
%     end
    
    fclose(fid_out);
    
    clear('MLT_matrix','MLT_new');
    
end

%% MLAT

jd2000_vec = [];
if ~exist(sprintf([path_to_output_files '%s_%s_supermag_MLAT.txt'],startTime(1:10),endTime(1:10)), 'file')
    
    iterIndex = 1;
    for ii = 1:length(list)
        [jd2000, IAGA, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = supermagdataall2matrix(filepath{ii});
        B = sqrt(N.^2 + E.^2 + Z.^2);
        [~, matchIndex] = ismember(IAGA,IAGA_all);
        
        %     MLT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLT;
        MLAT_matrix(iterIndex:length(MLAT)*ii,matchIndex) = MLAT;
        %     IGRF_DECL_matrix(iterIndex:length(MLT)*ii,matchIndex) = IGRF_DECL;
        %     SZA_matrix(iterIndex:length(MLT)*ii,matchIndex) = SZA;
        %     N_matrix(iterIndex:length(MLT)*ii,matchIndex) = N;
        %     E_matrix(iterIndex:length(MLT)*ii,matchIndex) = E;
        %     Z_matrix(iterIndex:length(MLT)*ii,matchIndex) = Z;
        %     B_matrix(iterIndex:length(MLT)*ii,matchIndex) = B;
        
        jd2000_vec = [jd2000_vec; jd2000];
        
        iterIndex = (ii * 1440)+1;
    end
    
    MLAT_new = MLAT_matrix(:,columnIndex);
    MLAT_new(find(MLAT_new == 0)) = NaN;
    
    MLAT_final = horzcat(jd2000_vec,MLAT_new);
    % Save data
    
    filename_out = sprintf([path_to_output_files '%s_%s_supermag_MLAT.txt'],startTime(1:10),endTime(1:10));
    fid_out = fopen(filename_out,'w');
    
     fprintf(fid_out,'%s' ,'Time (jd2000) ');
    for mm = 1:length(columnIndex)
        fprintf(fid_out,'%s ',IAGA_all{columnIndex(mm)});%,IAGA_all{columnIndex(2)},IAGA_all{columnIndex(3)},IAGA_all{columnIndex(4)},IAGA_all{columnIndex(5)},IAGA_all{columnIndex(6)},IAGA_all{columnIndex(7)},IAGA_all{columnIndex(8)},IAGA_all{columnIndex(9)},IAGA_all{columnIndex(10)},IAGA_all{columnIndex(11)},IAGA_all{columnIndex(12)},IAGA_all{columnIndex(13)},IAGA_all{columnIndex(14)},IAGA_all{columnIndex(15)},IAGA_all{columnIndex(16)},IAGA_all{columnIndex(17)},IAGA_all{columnIndex(18)});
    end
    
    fprintf(fid_out,'\n');
    [rows cols] = size(MLAT_final);
    x = repmat('%10.5f\t',1,(cols-1));
    fprintf(fid_out,[x,'%10.5f\n'],MLAT_final');

    fclose(fid_out);
    
    clear('MLAT_matrix','MLAT_new');
    
end
%% IGRF_DECL

jd2000_vec = [];
if ~exist(sprintf([path_to_output_files '%s_%s_supermag_IGRF_DECL.txt'],startTime(1:10),endTime(1:10)), 'file')
    iterIndex = 1;
    for ii = 1:length(list)
        [jd2000, IAGA, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = supermagdataall2matrix(filepath{ii});
        B = sqrt(N.^2 + E.^2 + Z.^2);
        [~, matchIndex] = ismember(IAGA,IAGA_all);
        
        %     MLT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLT;
        %     MLAT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLAT;
        IGRF_DECL_matrix(iterIndex:length(IGRF_DECL)*ii,matchIndex) = IGRF_DECL;
        %     SZA_matrix(iterIndex:length(MLT)*ii,matchIndex) = SZA;
        %     N_matrix(iterIndex:length(MLT)*ii,matchIndex) = N;
        %     E_matrix(iterIndex:length(MLT)*ii,matchIndex) = E;
        %     Z_matrix(iterIndex:length(MLT)*ii,matchIndex) = Z;
        %     B_matrix(iterIndex:length(MLT)*ii,matchIndex) = B;
        
        jd2000_vec = [jd2000_vec; jd2000];
        
        iterIndex = (ii * 1440)+1;
    end
    
    IGRF_DECL_new = IGRF_DECL_matrix(:,columnIndex);
    IGRF_DECL_new(find(IGRF_DECL_new == 0)) = NaN;
    
    IGRF_DECL_final = horzcat(jd2000_vec,IGRF_DECL_new);
    % Save data
    
    filename_out = sprintf([path_to_output_files '%s_%s_supermag_IGRF_DECL.txt'],startTime(1:10),endTime(1:10));
    fid_out = fopen(filename_out,'w');
    
     fprintf(fid_out,'%s' ,'Time (jd2000) ');
    for mm = 1:length(columnIndex)
        fprintf(fid_out,'%s ',IAGA_all{columnIndex(mm)});%,IAGA_all{columnIndex(2)},IAGA_all{columnIndex(3)},IAGA_all{columnIndex(4)},IAGA_all{columnIndex(5)},IAGA_all{columnIndex(6)},IAGA_all{columnIndex(7)},IAGA_all{columnIndex(8)},IAGA_all{columnIndex(9)},IAGA_all{columnIndex(10)},IAGA_all{columnIndex(11)},IAGA_all{columnIndex(12)},IAGA_all{columnIndex(13)},IAGA_all{columnIndex(14)},IAGA_all{columnIndex(15)},IAGA_all{columnIndex(16)},IAGA_all{columnIndex(17)},IAGA_all{columnIndex(18)});
    end
    fprintf(fid_out,'\n');
    [rows cols] = size(IGRF_DECL_final);
    x = repmat('%10.5f\t',1,(cols-1));
    fprintf(fid_out,[x,'%10.5f\n'],IGRF_DECL_final');

    fclose(fid_out);
    
    clear('IGRF_DECL_matrix','IGRF_DECL_new');
    
end
%% SZA

jd2000_vec = [];
if ~exist(sprintf([path_to_output_files '%s_%s_supermag_SZA.txt'],startTime(1:10),endTime(1:10)), 'file')
    iterIndex = 1;
    for ii = 1:length(list)
        [jd2000, IAGA, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = supermagdataall2matrix(filepath{ii});
        B = sqrt(N.^2 + E.^2 + Z.^2);
        [~, matchIndex] = ismember(IAGA,IAGA_all);
        
        %     MLT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLT;
        %     MLAT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLAT;
        %     IGRF_DECL_matrix(iterIndex:length(MLT)*ii,matchIndex) = IGRF_DECL;
        SZA_matrix(iterIndex:length(SZA)*ii,matchIndex) = SZA;
        %     N_matrix(iterIndex:length(MLT)*ii,matchIndex) = N;
        %     E_matrix(iterIndex:length(MLT)*ii,matchIndex) = E;
        %     Z_matrix(iterIndex:length(MLT)*ii,matchIndex) = Z;
        %     B_matrix(iterIndex:length(MLT)*ii,matchIndex) = B;
        
        jd2000_vec = [jd2000_vec; jd2000];
        
        iterIndex = (ii * 1440)+1;
    end
    
    SZA_new = SZA_matrix(:,columnIndex);
    SZA_new(find(SZA_new == 0)) = NaN;
    
    SZA_final = horzcat(jd2000_vec,SZA_new);
    % Save data
    
    filename_out = sprintf([path_to_output_files '%s_%s_supermag_SZA.txt'],startTime(1:10),endTime(1:10));
    fid_out = fopen(filename_out,'w');
    
     fprintf(fid_out,'%s' ,'Time (jd2000) ');
    for mm = 1:length(columnIndex)
        fprintf(fid_out,'%s ',IAGA_all{columnIndex(mm)});%,IAGA_all{columnIndex(2)},IAGA_all{columnIndex(3)},IAGA_all{columnIndex(4)},IAGA_all{columnIndex(5)},IAGA_all{columnIndex(6)},IAGA_all{columnIndex(7)},IAGA_all{columnIndex(8)},IAGA_all{columnIndex(9)},IAGA_all{columnIndex(10)},IAGA_all{columnIndex(11)},IAGA_all{columnIndex(12)},IAGA_all{columnIndex(13)},IAGA_all{columnIndex(14)},IAGA_all{columnIndex(15)},IAGA_all{columnIndex(16)},IAGA_all{columnIndex(17)},IAGA_all{columnIndex(18)});
    end
    fprintf(fid_out,'\n');
    [rows cols] = size(SZA_final);
    x = repmat('%10.5f\t',1,(cols-1));
    fprintf(fid_out,[x,'%10.5f\n'],SZA_final');

    fclose(fid_out);
    
    clear('SZA_matrix','SZA_new');
    
end
%% N

jd2000_vec = [];
if ~exist(sprintf([path_to_output_files '%s_%s_supermag_N.txt'],startTime(1:10),endTime(1:10)), 'file')
    iterIndex = 1;
    for ii = 1:length(list)
        [jd2000, IAGA, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = supermagdataall2matrix(filepath{ii});
        B = sqrt(N.^2 + E.^2 + Z.^2);
        [~, matchIndex] = ismember(IAGA,IAGA_all);
        
        %     MLT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLT;
        %     MLAT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLAT;
        %     IGRF_DECL_matrix(iterIndex:length(MLT)*ii,matchIndex) = IGRF_DECL;
        %     SZA_matrix(iterIndex:length(MLT)*ii,matchIndex) = SZA;
        N_matrix(iterIndex:length(N)*ii,matchIndex) = N;
        %     E_matrix(iterIndex:length(MLT)*ii,matchIndex) = E;
        %     Z_matrix(iterIndex:length(MLT)*ii,matchIndex) = Z;
        %     B_matrix(iterIndex:length(MLT)*ii,matchIndex) = B;
        
        jd2000_vec = [jd2000_vec; jd2000];
        
        iterIndex = (ii * 1440)+1;
    end
    
    N_new = N_matrix(:,columnIndex);
    N_new(find(N_new == 0)) = NaN;
    
    N_final = horzcat(jd2000_vec,N_new);
    % Save data
    
    filename_out = sprintf([path_to_output_files '%s_%s_supermag_N.txt'],startTime(1:10),endTime(1:10));
    fid_out = fopen(filename_out,'w');
    
     fprintf(fid_out,'%s' ,'Time (jd2000) ');
    for mm = 1:length(columnIndex)
        fprintf(fid_out,'%s ',IAGA_all{columnIndex(mm)});%,IAGA_all{columnIndex(2)},IAGA_all{columnIndex(3)},IAGA_all{columnIndex(4)},IAGA_all{columnIndex(5)},IAGA_all{columnIndex(6)},IAGA_all{columnIndex(7)},IAGA_all{columnIndex(8)},IAGA_all{columnIndex(9)},IAGA_all{columnIndex(10)},IAGA_all{columnIndex(11)},IAGA_all{columnIndex(12)},IAGA_all{columnIndex(13)},IAGA_all{columnIndex(14)},IAGA_all{columnIndex(15)},IAGA_all{columnIndex(16)},IAGA_all{columnIndex(17)},IAGA_all{columnIndex(18)});
    end
    fprintf(fid_out,'\n');
    [rows cols] = size(N_final);
    x = repmat('%10.5f\t',1,(cols-1));
    fprintf(fid_out,[x,'%10.5f\n'],N_final');

    fclose(fid_out);
    
    clear('N_matrix','N_new');
    
end
%% E

jd2000_vec = [];
if ~exist(sprintf([path_to_output_files '%s_%s_supermag_E.txt'],startTime(1:10),endTime(1:10)), 'file')
    iterIndex = 1;
    for ii = 1:length(list)
        [jd2000, IAGA, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = supermagdataall2matrix(filepath{ii});
        B = sqrt(N.^2 + E.^2 + Z.^2);
        [~, matchIndex] = ismember(IAGA,IAGA_all);
        
        %     MLT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLT;
        %     MLAT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLAT;
        %     IGRF_DECL_matrix(iterIndex:length(MLT)*ii,matchIndex) = IGRF_DECL;
        %     SZA_matrix(iterIndex:length(MLT)*ii,matchIndex) = SZA;
        %     N_matrix(iterIndex:length(MLT)*ii,matchIndex) = N;
        E_matrix(iterIndex:length(E)*ii,matchIndex) = E;
        %     Z_matrix(iterIndex:length(MLT)*ii,matchIndex) = Z;
        %     B_matrix(iterIndex:length(MLT)*ii,matchIndex) = B;
        
        jd2000_vec = [jd2000_vec; jd2000];
        
        iterIndex = (ii * 1440)+1;
    end
    
    E_new = E_matrix(:,columnIndex);
    E_new(find(E_new == 0)) = NaN;
    
    E_final = horzcat(jd2000_vec,E_new);
    % Save data
    
    filename_out = sprintf([path_to_output_files '%s_%s_supermag_E.txt'],startTime(1:10),endTime(1:10));
    fid_out = fopen(filename_out,'w');
    
     fprintf(fid_out,'%s' ,'Time (jd2000) ');
    for mm = 1:length(columnIndex)
        fprintf(fid_out,'%s ',IAGA_all{columnIndex(mm)});%,IAGA_all{columnIndex(2)},IAGA_all{columnIndex(3)},IAGA_all{columnIndex(4)},IAGA_all{columnIndex(5)},IAGA_all{columnIndex(6)},IAGA_all{columnIndex(7)},IAGA_all{columnIndex(8)},IAGA_all{columnIndex(9)},IAGA_all{columnIndex(10)},IAGA_all{columnIndex(11)},IAGA_all{columnIndex(12)},IAGA_all{columnIndex(13)},IAGA_all{columnIndex(14)},IAGA_all{columnIndex(15)},IAGA_all{columnIndex(16)},IAGA_all{columnIndex(17)},IAGA_all{columnIndex(18)});
    end
    fprintf(fid_out,'\n');
    [rows cols] = size(E_final);
    x = repmat('%10.5f\t',1,(cols-1));
    fprintf(fid_out,[x,'%10.5f\n'],E_final');

    fclose(fid_out);
    
    clear('E_matrix','E_new');
    
end
%% Z

jd2000_vec = [];
if ~exist(sprintf([path_to_output_files '%s_%s_supermag_Z.txt'],startTime(1:10),endTime(1:10)), 'file')
    iterIndex = 1;
    for ii = 1:length(list)
        [jd2000, IAGA, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = supermagdataall2matrix(filepath{ii});
        B = sqrt(N.^2 + E.^2 + Z.^2);
        [~, matchIndex] = ismember(IAGA,IAGA_all);
        
        %     MLT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLT;
        %     MLAT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLAT;
        %     IGRF_DECL_matrix(iterIndex:length(MLT)*ii,matchIndex) = IGRF_DECL;
        %     SZA_matrix(iterIndex:length(MLT)*ii,matchIndex) = SZA;
        %     N_matrix(iterIndex:length(MLT)*ii,matchIndex) = N;
        %     E_matrix(iterIndex:length(MLT)*ii,matchIndex) = E;
        Z_matrix(iterIndex:length(Z)*ii,matchIndex) = Z;
        %     B_matrix(iterIndex:length(MLT)*ii,matchIndex) = B;
        
        jd2000_vec = [jd2000_vec; jd2000];
        
        iterIndex = (ii * 1440)+1;
    end
    
    Z_new = Z_matrix(:,columnIndex);
    Z_new(find(Z_new == 0)) = NaN;
    
    Z_final = horzcat(jd2000_vec,Z_new);
    % Save data
    
    filename_out = sprintf([path_to_output_files '%s_%s_supermag_Z.txt'],startTime(1:10),endTime(1:10));
    fid_out = fopen(filename_out,'w');
    
     fprintf(fid_out,'%s' ,'Time (jd2000) ');
    for mm = 1:length(columnIndex)
        fprintf(fid_out,'%s ',IAGA_all{columnIndex(mm)});%,IAGA_all{columnIndex(2)},IAGA_all{columnIndex(3)},IAGA_all{columnIndex(4)},IAGA_all{columnIndex(5)},IAGA_all{columnIndex(6)},IAGA_all{columnIndex(7)},IAGA_all{columnIndex(8)},IAGA_all{columnIndex(9)},IAGA_all{columnIndex(10)},IAGA_all{columnIndex(11)},IAGA_all{columnIndex(12)},IAGA_all{columnIndex(13)},IAGA_all{columnIndex(14)},IAGA_all{columnIndex(15)},IAGA_all{columnIndex(16)},IAGA_all{columnIndex(17)},IAGA_all{columnIndex(18)});
    end
    fprintf(fid_out,'\n');
    [rows cols] = size(Z_final);
    x = repmat('%10.5f\t',1,(cols-1));
    fprintf(fid_out,[x,'%10.5f\n'],Z_final');

    fclose(fid_out);
    
    clear('Z_matrix','Z_new');
    
end

%% B

jd2000_vec = [];
if ~exist(sprintf([path_to_output_files '%s_%s_supermag_B.txt'],startTime(1:10),endTime(1:10)), 'file')
    iterIndex = 1;
    for ii = 1:length(list)
        [jd2000, IAGA, MLT, MLAT, IGRF_DECL, SZA, N, E, Z] = supermagdataall2matrix(filepath{ii});
        B = sqrt(N.^2 + E.^2 + Z.^2);
        [~, matchIndex] = ismember(IAGA,IAGA_all);
        
        %     MLT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLT;
        %     MLAT_matrix(iterIndex:length(MLT)*ii,matchIndex) = MLAT;
        %     IGRF_DECL_matrix(iterIndex:length(MLT)*ii,matchIndex) = IGRF_DECL;
        %     SZA_matrix(iterIndex:length(MLT)*ii,matchIndex) = SZA;
        %     N_matrix(iterIndex:length(MLT)*ii,matchIndex) = N;
        %     E_matrix(iterIndex:length(MLT)*ii,matchIndex) = E;
        %     Z_matrix(iterIndex:length(MLT)*ii,matchIndex) = Z;
        B_matrix(iterIndex:length(B)*ii,matchIndex) = B;
        
        jd2000_vec = [jd2000_vec; jd2000];
        
        iterIndex = (ii * 1440)+1;
    end
    
    B_new = B_matrix(:,columnIndex);
    B_new(find(B_new == 0)) = NaN;
    
    B_final = horzcat(jd2000_vec,B_new);
    % Save data
    
    filename_out = sprintf([path_to_output_files '%s_%s_supermag_B.txt'],startTime(1:10),endTime(1:10));
    fid_out = fopen(filename_out,'w');
    
     fprintf(fid_out,'%s' ,'Time (jd2000) ');
    for mm = 1:length(columnIndex)
        fprintf(fid_out,'%s ',IAGA_all{columnIndex(mm)});%,IAGA_all{columnIndex(2)},IAGA_all{columnIndex(3)},IAGA_all{columnIndex(4)},IAGA_all{columnIndex(5)},IAGA_all{columnIndex(6)},IAGA_all{columnIndex(7)},IAGA_all{columnIndex(8)},IAGA_all{columnIndex(9)},IAGA_all{columnIndex(10)},IAGA_all{columnIndex(11)},IAGA_all{columnIndex(12)},IAGA_all{columnIndex(13)},IAGA_all{columnIndex(14)},IAGA_all{columnIndex(15)},IAGA_all{columnIndex(16)},IAGA_all{columnIndex(17)},IAGA_all{columnIndex(18)});
    end
    fprintf(fid_out,'\n');
    [rows cols] = size(B_final);
    x = repmat('%10.5f\t',1,(cols-1));
    fprintf(fid_out,[x,'%10.5f\n'],B_final');

    fclose(fid_out);
    
    clear('B_matrix','B_new');
end
end
%%

% IAGA variable contains relative indexed information for all available
% stations (alphabetized). Each station index corresponds to a column in any of the data
% vectors.


