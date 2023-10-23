function [fout,TR,peak1,peak2,peak3] = TrackingVerification4(iter,Noisy_show,e_3,lastf,i,Noises,TR,pkno1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

   
if e_3(1)>60 && e_3(1)<210
    fout=e_3(1);
elseif iter>1
    fout=lastf(i-1);
else fout=e_3(2);
end
peak1=e_3(1);
peak2=e_3(2);
peak3=e_3(3);
if iter>1  %%&& sum(Noisy_show(i,:))>1 %% for PPG in NOISE
    
    hvr_tol=10;
    %% checking first pk is very close to any noise freq
    %% if very close then swap 1st and 2nd pk
    %%

    %     if e_3(1)==0 && (e_3(3)<60 || e_3(3)>180)
    %         e_3(1)=e_3(2);
    %         e_3(2)=0;
    %         e_3(3)=0;
    %     end
    
    % nos_tol=2;
    % tmpp=e_3;
    %     if (sum(abs(e_3(1)-Noises)<nos_tol) && e_3(1)>60 && e_3(1)<210)...
    %             && (sum(abs(e_3(2)-Noises)<nos_tol)  || e_3(2)<60  || e_3(2)>210) &&...
    %            (sum(abs(e_3(3)-Noises)<nos_tol)  || e_3(3)<60  || e_3(3)>210)
    %            e_3=e_3;
    %     else
    % for ii=3:-1:1
    %     if sum(abs(e_3(ii)-Noises)<nos_tol) || e_3(ii)==0
    %         e_3(ii)=0;
    %     end
    % end
    %     end
    % if ~length(e_3)
    %     e_3(2)=tmpp(1);
    % end
    %     e_3=[e_3 zeros(1,3-length(e_3)) ];
    
     nos_tol=2;
    if sum(abs(e_3(1)-Noises)<nos_tol) && e_3(1)>60 && e_3(1)<210
        if sum(abs(e_3(2)-Noises)<nos_tol) && e_3(2)>60 && e_3(2)<210
            if sum(abs(e_3(3)-Noises)<nos_tol) && e_3(3)>60 && e_3(3)<210
                e_3(3)=e_3(2);
                e_3(2)=e_3(1);
                e_3(1)=e_3(1);
            else
                tmp=e_3(3);
                e_3(3)=e_3(2);
                e_3(2)=e_3(1);
                e_3(1)=tmp;
            end
        else
            tmp=e_3(1);
            e_3(1)=e_3(2);
            e_3(2)=tmp;
            
        end
    end
    
    peak1=e_3(1);
    peak2=e_3(2);
    peak3=e_3(3);
    
    %% SS after swapping pks priority then the nearest pks of last BPM is set as HR for further investigation
    %%
    for ii=1:3
        nos_frm(ii)=sum(abs(e_3(ii)-Noises)<nos_tol);
    end
    vval=abs(e_3-lastf(i-1 ));
    [vv,indd]=sort(vval,'ascend');
    
    fout=0;
    for ii=3:-1:1
        if ~nos_frm(ii)
            if ii~=1 && abs(e_3(ii)-lastf(i-1 ))<5 && abs(abs(e_3(1)-lastf(i-1))-abs(e_3(ii)-lastf(i-1 )))>5
                fout=e_3(ii);
                vvin=indd(ii);
                break
            elseif ii==1 && abs(e_3(ii)-lastf(i-1 ))<25
                fout=e_3(ii);
                vvin=indd(ii);
                break
            end
        end
    end
    
    if fout==0
        [vval,vvin]=min(abs(e_3-lastf(i-1 )));
        fout=e_3(vvin);
    end
    %% Resolving the pks priority and their relevance
    prop_detec=0;
    if  ((fout-lastf(i-1))>20 || (lastf(i-1)-fout)>10 || fout<60 || fout>200) && vvin==1
        fout=lastf(i-1);
    elseif (vvin==2 || vvin==3) && ((fout-lastf(i-1))>5 || (lastf(i-1)-fout)>5 || fout<60 || fout>200)
        fout=lastf(i-1);
    elseif (vvin==2 || vvin==3) &&  abs(abs(e_3(1)-lastf(i-1))-abs(fout-lastf(i-1)))<5
        fout=e_3(1);
    end
    
    
    if iter>2
        TR(i)=.8*(fout-lastf(i-1))+.2*(TR(i-1)-TR(i-2));
    else
        TR(i)=0;
    end
    
    %% Using trends and other condition for accurate HR
    %%
    
    if iter>1
        if TR(i-1)==0
            if abs( fout-lastf(i-1))>hvr_tol
                fout=lastf(i-1)+5*sign(fout-lastf(i-1));
            end
        else
            
            if abs( fout-lastf(i-1))>abs(2*TR(i-1)) && abs( fout-lastf(i-1))>5
                fout=lastf(i-1)+TR(i);
            elseif abs( fout-lastf(i-1))<abs(2*TR(i-1)) || abs( fout-lastf(i-1))<5
                fout=fout;
            else abs( fout-lastf(i-1))>2*hvr_tol
                fout=lastf(i-1);
            end
        end
    end
    
    %% if three consecutive 1st pks are found very close then HR will go there
    %%
    pkno1(i)=peak1;
    no_of_ip=3;
    if fout~=e_3(1) && iter>no_of_ip && prop_detec~=1
        pkk=pkno1(i-no_of_ip+1:i);
        yyy=(pkk==0);
        pkk;
        corrr=corr([1:length(pkk)]',pkk');
        if sum(yyy)==0 && std(abs(diff(pkk)))<3 && mean(abs(diff(pkk)))<5 ...
                && e_3(1)< lastf(i-1)*1.9 && e_3(1)> lastf(i-1)*.7
            std(abs(diff(pkk)));
            fout=e_3(1);
            TR(1:i)=zeros(1,i);
        end
    end
else TR(1:i)=zeros(1,i);
end
end

