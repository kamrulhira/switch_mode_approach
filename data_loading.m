function [LOD,ECG] = data_loading(ki)


if ki==1
    load DATA_01_TYPE01.mat
    load DATA_01_TYPE01_BPMtrace.mat
elseif ki==2
    load DATA_02_TYPE02.mat
    load DATA_02_TYPE02_BPMtrace.mat
elseif ki==3
    load DATA_03_TYPE02.mat
    load DATA_03_TYPE02_BPMtrace.mat
elseif ki==4
    load DATA_04_TYPE02.mat
    load DATA_04_TYPE02_BPMtrace.mat
elseif ki==5
    load DATA_05_TYPE02.mat
    load DATA_05_TYPE02_BPMtrace.mat
elseif ki==6
    load DATA_06_TYPE02.mat
    load DATA_06_TYPE02_BPMtrace.mat
elseif ki==7
    load DATA_07_TYPE02.mat
    load DATA_07_TYPE02_BPMtrace.mat
elseif ki==8
    load DATA_08_TYPE02.mat
    load DATA_08_TYPE02_BPMtrace.mat
elseif ki==9
    load DATA_09_TYPE02.mat
    load DATA_09_TYPE02_BPMtrace.mat
elseif ki==10
    load DATA_10_TYPE02.mat
    load DATA_10_TYPE02_BPMtrace.mat
elseif ki==11
    load DATA_11_TYPE02.mat
    load DATA_11_TYPE02_BPMtrace.mat
elseif ki==12
    load DATA_12_TYPE02.mat
    load DATA_12_TYPE02_BPMtrace.mat
    
elseif ki==13
    load DATA_S04_T01.mat
    load BPM_S04_T01.mat
    
elseif ki==14
    load TEST_S01_T01.mat
    load True_S01_T01.mat
elseif ki==15
    load TEST_S02_T01.mat
    load True_S02_T01.mat
elseif ki==16
    load TEST_S02_T02.mat
    load True_S02_T02.mat
elseif ki==17
    load TEST_S03_T02.mat
    load True_S03_T02.mat
elseif ki==18
    load TEST_S04_T02.mat
    load True_S04_T02.mat
elseif ki==19
    load TEST_S05_T02.mat
    load True_S05_T02.mat
elseif ki==20
    load TEST_S06_T01.mat
    load True_S06_T01.mat
elseif ki==21
    load TEST_S06_T02.mat
    load True_S06_T02.mat
elseif ki==22
    load TEST_S07_T02.mat
    load True_S07_T02.mat
elseif ki==23
    load TEST_S08_T01.mat
    load True_S08_T01.mat
end
LOD=sig;
ECG=BPM0;
end

