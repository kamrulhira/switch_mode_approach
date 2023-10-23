function [p_3,fmax] = find3peaks(N,NFFT,Fs)

L=length(N);
no_pk=3;
W=blackman(L);
 %[pxx_N,f_N]=periodogram(N,W,NFFT,Fs,'psd');
[pxx_N,f_N]=pwelch(N,round(L/2),round(.45*L),NFFT,Fs,'onesided'); 
 %figure(1),plot(f_N*60,pxx_N),xlim([60 400]),hold on;
b48=round(NFFT*1/Fs);
b200=round(NFFT*3/Fs);
pxx=zeros(1,b200);
pxx(1:b48)=pxx_N(b48);
pxx(b48:b200)=pxx_N(b48:b200); %% truncated part of periodogram 1hz to 3.3 hz
% figure,plot(f_N*60,pxx_N)

[val_pk_max,lcn]=findpeaks(pxx/max(pxx),'NPEAKS',1);

if(val_pk_max)
[pkn,lcn]=findpeaks(pxx/max(pxx),'MINPEAKDISTANCE',5,'NPEAKS',no_pk,...
    'SORTSTR','descend','MINPEAKHEIGHT',0.15*val_pk_max);
else
    lcn=ones(1,no_pk);
end

freq_n=f_N(lcn)*60;

p_3=zeros(1,no_pk); %% this is for if 3 pks are not found

p_3(1:length(freq_n))=freq_n;
fmax=p_3(1);

end

