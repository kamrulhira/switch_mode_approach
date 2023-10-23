%% $$$ Start
%%
clear all;close all;
clc;tic;format short G
a=[];
tic 
pth=[pwd '\NewData'];
addpath(pth);
warning off
[bh,ah] = butter(4,2*.5/25,'high');
[bl,al] = butter(4,2*3.25/25,'low');

%% $$$ Initialization
%%
save_on=0;
check=1;Xp=0;Yp=0;XYcount=1;cnt=1;
factorr=2;s_NS=[];s_val=[];
%% $$$ Loop Start
%%
sub=[1];%2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
for ki=sub    %% ki means ki th subject
    %% $$$ some initialization and data loading
    %%
    tic
    
    LOD=[];ECG=[];
    str=['Subject_' num2str(ki) '.mat'];
    load(str)
    LOD=[sigPPG;sigGyro;sigAcc];
    %LOD=[sigPPG;sigAcc;sigGyro];
    LOD=(downsample(LOD',factorr))';
    
    frm_tm=8;
    frm_ovr=2;
    ECG= GT_for_new_data(bpmECG,timeECG,frm_tm,frm_ovr);
    h = waitbar(0,'Please Wait');
    iter=0;              %why
    flag_cnt=0;
    lastf=0;
    hr_prev=0;
    up_noise=[];
    dn_noise=[];peak1=[];peak2=[];peak3=[];
    Noisy_show=[];TR(1)=0;
    tmm=[];tmm2=[];pnn=zeros(416,3);
    framee=1:length(ECG);
    se=[];sd=[];sp=[];sq=[];
    sx=[];sy=[];sz=[];
    sgx=[];sgy=[];sgz=[];
    P_All=[];
    for i=framee,waitbar(i/length(framee),h,sprintf('Please Wait \t Subject: %d \n Completed Frames:  %d/%d',ki,i,length(framee)))
        [ki i];
        iter=iter+1; %% counting how many frame is run
        NFFT=4096;
        Fs=50/factorr;
        
        alpha = 200;        % moderate bandwidth constraint
        tau = 0.5;            % noise-tolerance (no strict fidelity enforcement)
        
        K = 3;              % 3 modes
        DC = 0;             % no DC part imposed
        init = 1;           % initialize omegas uniformly
        tol = 1e-3;
        Wn = 1/Fs;
        
        E=LOD(1,((i-1)*Fs*frm_ovr+1:(i-1)*Fs*frm_ovr+Fs*frm_tm));
        P=LOD(2,((i-1)*Fs*frm_ovr+1:(i-1)*Fs*frm_ovr+Fs*frm_tm));
        Q=LOD(3,((i-1)*Fs*frm_ovr+1:(i-1)*Fs*frm_ovr+Fs*frm_tm));
        X=LOD(4,((i-1)*Fs*frm_ovr+1:(i-1)*Fs*frm_ovr+Fs*frm_tm));
        Y=LOD(5,((i-1)*Fs*frm_ovr+1:(i-1)*Fs*frm_ovr+Fs*frm_tm));
        Z=LOD(6,((i-1)*Fs*frm_ovr+1:(i-1)*Fs*frm_ovr+Fs*frm_tm));
        VZ(i)=var(X)+var(Y)+var(Z); %why??
        
        sigf=zscore(E)+zscore(P)+zscore(Q);% why??
        
        sigf=filter(bh,ah,sigf);
        sigf=filter(bl,al,sigf);
        %
        if VZ(i)>5e6 && iter>1
            siggn=zscore(X)+zscore(Z)+zscore(Z);
         
            [n_3,fmax2] = find3peaks(siggn,NFFT,Fs);
            
            siggn=filter(bh,ah,siggn);
            siggn=filter(bl,al,siggn);
            [imf, u_hat, omega] =MVMD(siggn, alpha, tau, K, DC, init, tol);
            
            size(imf);
            imf=imf';
            
            for ni=1:K
                [nnn_3,n_imf(ni)] = find3peaks(imf(:,ni),NFFT,Fs);
            end
            
            [ECG(i) lastf(i-1) n_imf]
            
            pp=sum(sum(imf.^2));
            
            [val,idd]=sort(abs(n_imf-lastf(i-1)));
           
            imf2=imf;
            n_imf2=n_imf;
            for ni=1:K
                imf(:,ni)=imf2(:,idd(ni));
                n_imf(ni)=n_imf2(:,idd(ni));
            end
            
            n_imf;
%             W=blackman(size(imf,1));
%              [kk,ff]=periodogram(imf,W,4096,Fs);
%               pecg=round(ECG(i)*NFFT/(Fs*60));
%                
%                          if i>50 && i<100
%                          W=blackman(length(siggn));
%                             
%                            [pxx_N,f_N]=periodogram(siggn,W,NFFT,Fs,'psd');
%                            figure,plot(f_N*60,pxx_N),xlim([0 240]),ylim([0 35]);hold on;
%                           plot(f_N(pecg)*60,pxx_N(pecg),'r^');%hold on;
%                           %[pxx_N,f_N]=periodogram(sigf,W,NFFT,Fs,'psd');
%                             %plot(f_N*60,pxx_N,'g'),xlim([0 240]);hold on;plot(f_N(pecg)*60,pxx_N(pecg),'g*');%hold on;
%                        xlabel('Bit Per-Minute (BPM)');
%                        ylabel('Magnitude');
%                        legend('Gyroscopic Signal');
%                        %legend('Noisy Signal','Actual HR','Noise Free Signal','Actual HR');
%                             [pxx_N,f_N]=periodogram(imf(:,1),W,NFFT,Fs,'psd');%hold on;
%                             figure,plot(f_N*60,pxx_N),xlim([0 240]),ylim([0 35]);hold on; plot(f_N(pecg)*60,pxx_N(pecg),'r^')
%                        xlabel('Bit Per-Minute (BPM)');
%                        ylabel('Magnitude');
%                        legend('Gyroscopic Signal');
%                              [pxx_N,f_N]=periodogram(imf(:,2),W,NFFT,Fs,'psd');
%                              figure,plot(f_N*60,pxx_N),xlim([0 240]),ylim([0 35]);hold on; plot(f_N(pecg)*60,pxx_N(pecg),'r^')
%                        xlabel('Bit Per-Minute (BPM)');
%                        ylabel('Magnitude');
%                        legend('Gyroscopic Signal');
%                              [pxx_N,f_N]=periodogram(imf(:,3),W,NFFT,Fs,'psd');
%                              figure,plot(f_N*60,pxx_N),xlim([0 240]),ylim([0 35]);hold on; plot(f_N(pecg)*60,pxx_N(pecg),'r^')
%                        xlabel('Bit Per-Minute (BPM)');
%                        ylabel('Magnitude');
%                        legend('Gyroscopic Signal');
%                           %return
%                          end 
                         
            e=sigf/max(sigf);
            
            for jj=1:K
                
                if abs(lastf(i-1)-n_imf(jj))<5
                   M=25;
                 stp=0.001;
                 
                 N11=imf(:,jj);
                N11=N11/max(N11);
                hh=dsp.LMSFilter('Length',M,'StepSizeSource','Property','StepSize',stp);
                %hh= dsp.RLSFilter('Length',M);
                e=reshape(e,length(e),1);
                N11=reshape(N11,length(N11),1);
                N11=N11/max(N11);
                [y,e]=hh(N11,e);
                
                
                 elseif abs(lastf(i-1)-n_imf(jj))<15
                   M=25;
                  stp=0.01; 
                 %stp=0.0012;
                 
                 N11=imf(:,jj);
                N11=N11/max(N11);
                hh=dsp.LMSFilter('Length',M,'StepSizeSource','Property','StepSize',stp);
                %hh= dsp.RLSFilter('Length',M);
                e=reshape(e,length(e),1);
                N11=reshape(N11,length(N11),1);
                N11=N11/max(N11);
                [y,e]=hh(N11,e);
                           
%                 
%                 elseif abs(lastf(i-1)-n_imf(jj))<35
%                    M=25;
%                  %stp=0.01; 
%                  stp=0.0012;
%                  
%                  N11=imf(:,jj);
%                 N11=N11/max(N11);
%                 hh=dsp.LMSFilter('Length',M,'StepSizeSource','Property','StepSize',stp);
%                 %hh= dsp.RLSFilter('Length',M);
%                 e=reshape(e,length(e),1);
%                 N11=reshape(N11,length(N11),1);
%                 N11=N11/max(N11);
%                 [y,e]=hh(N11,e);
%                 
                else
                    M=2*25;
                   stp=0.01;
                   
                   N11=imf(:,jj);
                N11=N11/max(N11);
                hh=dsp.LMSFilter('Length',M,'StepSizeSource','Property','StepSize',stp);
                %hh= dsp.RLSFilter('Length',M);
                e=reshape(e,length(e),1);
                N11=reshape(N11,length(N11),1);
                N11=N11/max(N11);
                [y,e]=hh(N11,e);
                end

            end
            
            e=xcorr(e,e);
            [p_3,fmax] = find3peaks(e,NFFT,Fs);
            
            fout=fmax;
            pkk(i,:)=[p_3];
            pnn(i,:)=[n_3];
            
        else
            
            siggn=zscore(X)+zscore(Y)+zscore(Z); %noise
            [n_3,fnmax] = find3peaks(siggn,NFFT,Fs);
            e=xcorr(sigf,sigf); %ppg
            [p_3,fmax] = find3peaks(e,NFFT,Fs);
            
            fout=fmax;
            pkk(i,:)=[p_3];
            pnn(i,:)=[n_3];
        end
        pkno1=pkk(:,1)';
        TR(1)=0;
%        pp_3(i,:)=p_3;
%         %[fout,TR,peak1,peak2,peak3] = Track_n_verify_new1(i,iter,p_3,hr_prev,n_3,TR,pkno1)
%         if i>=3
%          pp=[pp_3(i-1,:) pp_3(i,:) pp_3(i+1,:) ];
%         end
[fout,TR,peak1,peak2,peak3] = TrackingVerification_data04_01(iter,p_3,lastf,i,n_3,TR,pkno1);        
lastf(i)=fout;
                    pkk(i,:)=[peak1 peak2 peak3];
        
    end
    
 A(i,:)=[fout ECG(i)];
        
        Xp(XYcount)=ECG(i);
        Yp(XYcount)=fout;
        XYcount=XYcount+1;
        
    ERR(ki)=mean(abs(lastf-ECG))
        ERR2(ki)=mean(abs(lastf-ECG)./ECG)*100
 ERR3(ki)=std(lastf-ECG)
    a = [a ERR(ki)]
    figure,plot(pkk(:,1),'r*');hold on;plot(pkk(:,2),'bd');%plot(pkk(:,3),'k+');
    
    plot(ECG(framee));%plot(lastf);plot(pnn,'go')
    xlabel('Frame No.');
    ylabel('BPM');
    legend('Tracked BPM','Ground Truth');
    

PEARSON =Plot2figure(Xp,Yp,str,ERR3,sub,save_on);
   figure,plot(VZ)
    delete(h);
end
toc

%save 'All_vmd_related_pks.mat' ALL_DS ALL_X ALL_ECG