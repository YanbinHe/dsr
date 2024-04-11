clc;
close all;
clear;
% load the configuration
run('Configuration.m')
%%
for avg = 1 : AVG % for AVG trials

    % produce different realizations of channel
%     run('channel_realization.m')

    filename = ['./Results/compare2_', num2str(avg),'.mat'];
    load(filename)

    error = [];
    ser = [];
    time = [];
    supprecovery = [];
    noise_var_error = [];

    %% start simulation
    for jj = 1:JJ % snr

        SNR_10 = SNRl_10(jj);
        noise_var = (signal_power)/SNR_10;

        x_pilot = X(:,1:K1);
        IRS = irs_pattern(:,1:K2);
        % different irs reflection patterns, the number is equal to #overheads,
        H_p1 = IRS.'*kr(A_irs_a.',A_irs_d').'*10;
        H_p1 = H_p1(:,1:Res1);
        H_p1_ori = H_p1;
        H_p2 = x_pilot.'*conj(A1)/10;
        H_p2_ori = H_p2;
        H = kron(kron(H_p1,H_p2),A2);
        %% SNR part
        % received noisy signals
        y_tilde = vec(y_bar(1:K1*bs_ante,1:K2));
        % generate noise
        noise = sqrt(noise_var / 2)*(randn(size(y_tilde))+1i*randn(size(y_tilde)));
        y = y_tilde + noise;
        %% part 2: channel estimation with different techniques and compute the symbol error rate
        % true channel for each IRS pattern
        for i = 1:K2
            Htrue(:,:,i) = vec(H2*diag(IRS(:,i))*H1);
        end
        inPut = randi([0, modulating_scheme-1], 1, timeslots);
        Xun = qammod(inPut, modulating_scheme);

        %% different techniques
        %% classicSBL
        [time_csbl,x_rec,nvc] = classicSBL(noise_var,numItr,Res1,H,y,SNRl(jj));
        % metrics compute
        noise_var_error(4,jj) = abs(nvc-noise_var);
        [error_csbl, Hre_csbl] = ce_error(x_rec,Res1,K2,IRS,A_irs_a,A_irs_d,A1,A2,H1,H2);
        error(4,jj) = error_csbl;
        time(4,jj) = time_csbl;
        supprecovery(4,jj) = recover_rate(suppTrue,x_rec);
        ser(4,jj) = ser_compute(Hre_csbl,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2,inPut);
        %% decoKroSBL
        A{1} = H_p1;
        A{2} = H_p2;
        A{3} = A2;
        [time_dsbl,x_red] = deco_kroSBL(numItr,Res1,A,y);
        % metrics compute
        [error_dksbl, Hre_dksbl] = ce_error(x_red,Res1,K2,IRS,A_irs_a,A_irs_d,A1,A2,H1,H2);
        error(3,jj) = error_dksbl;
        time(3,jj) = time_dsbl;
        supprecovery(3,jj) = recover_rate(suppTrue,x_red);
        ser(3,jj) = ser_compute(Hre_dksbl,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2,inPut);
        %% SVD
        [time_svdsbl,x_res,nvs] = kroSBL_svd(noise_var,Res1,numItr,H,H_p1,H_p2,A2,y,SNRl(jj));
        % metrics compute
        noise_var_error(1,jj) = abs(nvs-noise_var);
        [error_svdsbl, Hre_svdsbl] = ce_error(x_res,Res1,K2,IRS,A_irs_a,A_irs_d,A1,A2,H1,H2);
        error(1,jj) = error_svdsbl;
        time(1,jj) = time_svdsbl;
        supprecovery(1,jj) = recover_rate(suppTrue,x_res);
        ser(1,jj) = ser_compute(Hre_svdsbl,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2,inPut);
        %% Alternating
        [time_am,x_rea,nva] = kroSBL_am(noise_var,Res1,numItr,H,H_p1,H_p2,A2,y,SNRl(jj));
        % metrics compute
        noise_var_error(2,jj) = abs(nva-noise_var);
        [error_am, Hre_am] = ce_error(x_rea,Res1,K2,IRS,A_irs_a,A_irs_d,A1,A2,H1,H2);
        error(2,jj) = error_am;
        time(2,jj) = time_am;
        supprecovery(2,jj) = recover_rate(suppTrue,x_rea);
        ser(2,jj) = ser_compute(Hre_am,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2,inPut);
        %% OMP
        [timeo,x_reo] = omp2(H,y,epi,40);
        % metrics compute
        [erroro, Hreo] = ce_error(x_reo,Res1,K2,IRS,A_irs_a,A_irs_d,A1,A2,H1,H2);
        error(5,jj) = erroro;
        time(5,jj) = timeo;
        supprecovery(5,jj) = recover_rate(suppTrue,x_reo);
        ser(5,jj) = ser_compute(Hreo,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2,inPut);
        %% decoOMP
        [time_do,x_redo] = deco_OMP(Res1,A,y,epi);
        % metrics compute
        [error_do, Hre_do] = ce_error(x_redo,Res1,K2,IRS,A_irs_a,A_irs_d,A1,A2,H1,H2);
        error(6,jj) = error_do;
        time(6,jj) = time_do;
        supprecovery(6,jj) = recover_rate(suppTrue,x_redo);
        ser(6,jj) = ser_compute(Hre_do,Htrue,SNRl(jj),Xun,symbolmatrix,ms_ante,bs_ante,K2,inPut);
    end
    filename = ['./Results/ce_', num2str(avg),'.mat'];
    save(filename,'noise_var_error','error','time','supprecovery','ser')
end