nFFT = 64; nDSC = 52; nSym = 1e4;  
EbN0dB = 0:10;  
EsN0dB = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80);  
nErr = zeros(size(EbN0dB));
for ii = 1:length(EbN0dB)
    ipBit = randi([0, 1], 1, nDSC * nSym);  
    ipMod = 2 * ipBit - 1;  
    ipMod = reshape(ipMod, nDSC, nSym).';  
    xF = [zeros(nSym, 6), ipMod(:, 1:nDSC/2), zeros(nSym, 1), ipMod(:, nDSC/2+1:end), zeros(nSym, 5)];
    xt = (nFFT/sqrt(nDSC)) * ifft(fftshift(xF.')).';  
    xt = [xt(:, 49:64), xt];  
    xt = xt(:).';  
    nt = (randn(size(xt)) + 1j * randn(size(xt))) / sqrt(2);  
    yt = sqrt(80/64) * xt + 10^(-EsN0dB(ii)/20) * nt;  
    yt = reshape(yt, 80, []).';  
    yt = yt(:, 17:80);  
    yF = (sqrt(nDSC)/nFFT) * fftshift(fft(yt.')).';  
    yMod = yF(:, [6+(1:nDSC/2), 7+(nDSC/2+1:end)]);  
    ipBitHat = reshape((sign(real(yMod)) + 1) / 2, 1, []);  
    nErr(ii) = sum(ipBitHat ~= ipBit);
end
simBer = nErr / (nSym * nDSC);  
theoryBer = 0.5 * erfc(sqrt(10.^(EbN0dB/10)));  
semilogy(EbN0dB, theoryBer, 'bs-', 'LineWidth', 2);
hold on;
semilogy(EbN0dB, simBer, 'mx-', 'LineWidth', 2);
grid on;
legend('Theory', 'Simulation');
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate');
title('BER for BPSK-OFDM');


