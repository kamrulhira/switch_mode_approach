function ECG= GT_for_new_data(bpmECG,timeECG,frm_tm,frm_ovr)
    frm_no=floor(timeECG(end)/frm_ovr-(frm_tm/frm_ovr-1));
    ECG=[];
    for i=1:frm_no
%         ECG(i)=median( bpmECG( find(timeECG>(i-1)*2 & timeECG<(i-1)*2+8 ) ) );
kk= sort(bpmECG( find(timeECG>(i-1)*frm_ovr & timeECG<(i-1)*2+frm_tm ) )) ;
        ECG(i)=mean(kk);

    end
end

