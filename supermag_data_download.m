function supermag_data_download(username,startTime,endTime,stations,path_to_download)
 % 2012-09-30T00:00:00.000Z % 24:00 (URLcolon=%3A) (URLcomma=%2C)
% ALE,THL,SVS,RES,KUV,NRD,UPN,CY0,IGC,TAL,GHC,DMH,CBB,UMQ,RPB,NAL,GDH,LYR,DNB,PGC,ATU,CHC,HRN,CDC,BLC,VIZ,CNL,HOP,STF,RAN,IQA,SKT,MCR,SCO,BJN,INK,EKP,T35,KAV,DED,GHB,BRW,JAN,C07,YKC,AMK,FCC,DIK,FHB,NOR,SMI,T29,FSP,FYU,RAL,SOR,T31,GIM,TRO,EAG,AND,KEV,DAW,MAS,NAQ,TIK,PBQ,KIL,PKR,PBK,T41,AMD,CGO,CMO,LEK,ABK,CHD,LRV,MUO,KIR,FMC,ISL,HLL,LOZ,T42,T38,SOD,T33,PEL,JCK,GAK,DON,T39,T36,MEA,RVK,T40,OUL,S01,LYC,T28,C06,OUJ,PIN,FAR,T32,T22,T30,SIT,RED,T37,DOB,C12,C04,HAN,C11,SOL,T43,T03,LER,T15,M01,LET,NUR,C08,UPS,KAR,OTT,T24,KVI,LOV,CLK,T17,C10,NEW,GML,T23,YAK,TAR,BOR,BOX,MGD,VIC,GTF,STJ,C01,SHU,ESK,MSH,T21,ARS,BRZ,BFE,M04,ROE,KNZ,T25,YOR,NVS,HLP,WNG,SUW,VAL,BOU,FRD,SZC,M05,NGK,HAD,BEL,DSO,IRT,MAB,KIV,DOU,PET,ZAG,LVV,M06,T16,BDV,VYH,BFO,CLF,FUR,HRB,NCK,PAG,FRN,THY,ODE,CST,BSL,KHB,SUA,P01,GCK,TUC,M08,DLR,WMQ,RNC,AAA,CNH,MSR,MMB,AQU,RIK,ISK,DUR,IZN,BMT,EBR,SPT,MIZ,PEG,LZH,CYG,E05,KAK,SJG,A02,SFS,ONW,HTY,MID,KAG,KNY,HON,JAI,CBI,LNP,GZH,GUI,PHU,ABG,HYB,MUT,KOU,TAM,A12,GUA,A04,DLT,AAE,HUA,MBO,A13,A07,A01,A06,A03,BNG,VRE,GAN,KTB,BIK,ASC,API,PPT,SER,VSS,PIL,IPM,EWA,KDU,CKI,SHE,OSO,PAC,TAN,CTA,TRW,A05,TSU,LRM,ASP,HBK,KMH,PST,TDC,HER,GNG,GNA,CNB,ESC,ORC,LIV,OHI,AMS,PAL,EYR,AIA,B03,CZT,PAF,B11,B12,MCQ,DMC,B14,B07,B16,B18,MAW,PG4,B10,B19,B23,PG3,B21,SPA,PG2,B22,PG1,B20,SBA,MCM,DRV,CSY,VOS

%     startTime = '2012-09-30T00:00:00.000Z';
% %     endTime = '2014-04-01T00:00:00.000Z';
%     endTime = '2012-10-03T00:00:00.000Z';
%     stations = 'DOU,NUR,SOD,NGK,WNG,FUR,BFO,VAL,ESK,HAD,LER,CLF,MAB,ABG,FRD,PG3,PG1,PG2,';
    
    start_num = datenum(regexprep(startTime,'[A-Z]',' '));
    end_num = datenum(regexprep(endTime,'[A-Z]',' '));
    
    iterations = end_num-start_num;
    
    counter = 0;
    for ii = 1:iterations
        startTime = datestr(start_num + counter,'yyyy-mm-ddTHH:MM:SS.000Z');
        if exist(sprintf([path_to_download 'supermag_%s.csv'],startTime(1:10)), 'file')
            counter = counter+1;
            continue
        else
            api = 'http://supermag.jhuapl.edu/mag/lib/services/?';
            url = sprintf([api '?user=%s&start=%s&interval=23:59&service=mag&stations=%s&delta=none&baseline=all&options=+mlt+sza+decl&format=csv'],username,startTime,stations);
            path = sprintf([path_to_download 'supermag_%s.csv'],startTime(1:10));
            outfilename = websave(path,url);
            counter = counter+1;
        end
    end
end    
