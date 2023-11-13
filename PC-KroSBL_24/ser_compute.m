function [ser] = ser_compute(H_est,H_true,SNR,X,symbolmatrix,ms_ante,bs_ante,K2,inPut)

% H_est: estimated channel matrix (in the form of column vector)
% H_true: the true channel matrix (in the form of column vector)
% SNR: scenario snr
% X: transmitted symbols (row vector)
% inPut: the true label of transmitting symbols
% symbolmatrix: #codes in codebook x timeslot

% We first generate a certain number of symbols x according to the
% variable sym. Then we transmit these symbols through the channel by
% multiplying them and obtain y. The detector T is designed by Linear MMSE 
% detector. The received symbols y_mmse are obtained by Ty. To reterive the
% symbols, we project the received symbols on the discrete constellation by
% computing the Euclidean distance.


ser = 0;

for i = 1:K2
    Htrue = reshape(H_true(:,:,i),bs_ante,ms_ante);
    Hest = reshape(H_est(:,:,i),bs_ante,ms_ante);
    
    % generate beamformer
    [s,~,d] = svd(Hest);
    precoder = d(:,1);
    decoder = s(:,1);
    
    SNRl = 10^(SNR/10);
 
    Yo = decoder'*Htrue*precoder*(1/sqrt(ms_ante))*X;
    
    signal_power = norm(vec(Yo))^2 / numel(Yo);
    sigmasq = signal_power/SNRl;
    noise = sqrt(sigmasq / 2)*(randn(size(Yo))+1i*randn(size(Yo)));
    Y = Yo + noise; % received signal: row vector

    matcompare = decoder'*Hest*precoder*(1/sqrt(ms_ante))*symbolmatrix;

    for j = 1:size(matcompare,1)
        matcompare(j,:) = matcompare(j,:) - Y;
    end

    [~,decode_sol] = min(abs(matcompare),[],1);

    error = length(find(decode_sol ~= (inPut+1)));
    
    ser = ser + error/numel(X);
end
ser = ser/K2;
end
