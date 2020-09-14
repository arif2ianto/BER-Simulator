%A script to compare the BER of BASK, BFSK, BPSK, and DPSK from simulation and theory

%Written by: ARIF RAHMADIAN ARIFIANTO - 18117016

%Reference:
%#https://www.mathworks.com/matlabcentral/fileexchange/67649-matlab-code-for-ber-performance-of-binary-ask-digital-modulation
%#https://www.mathworks.com/matlabcentral/fileexchange/28649-ber-of-bfsk-in-awgn-channel
%#https://www.mathworks.com/matlabcentral/fileexchange/66310-ber-for-bpsk-over-awgn-channel-vs-theoretical

%Variables: 
%Sig_Length - bit length of the simulated signal
%Eb - bit energy of the simulated signal

Sig_Length = 200000;
Eb = 1;

%Eb/No for theoretical calculations (high resolution)
EtoN_dB = linspace(0,20,100);
EtoN = 10.^(EtoN_dB/10);

%Eb/No for simulation (five data points)
EtoN_dB_sim = linspace(0,20,5);
EtoN_sim = 10.^(EtoN_dB_sim/10);

%No for simulation purposes
No = Eb./EtoN_sim;

%Theoretical BER calculations
BER_BASK_te = (1/2)*erfc(sqrt(EtoN/4));
BER_BFSK_te = (1/2)*erfc(sqrt(EtoN/2));
BER_BPSK_te = (1/2)*erfc(sqrt(EtoN));
BER_DPSK_te = (1/2)*exp(-EtoN);

for iter = 1:length(EtoN_sim)  
    %Initialize error as 0
    E_BASK = 0;
    E_BFSK = 0;
    E_BPSK = 0;
    E_DPSK = 0;
    
    %Generate random string of bits with the length specified above
    bit = randi([0 1],1,Sig_Length);
    
    %Modulation of BFSK signal
    x_BFSK = bit+1j*(~bit);
    
    %Modulation of BASK signal
    x_BASK = bit;
    
    %Modulation of BPSK signal
    x_BPSK = 2*bit-1;
    
    %Modulation of DPSK signal
    %Defining the first bit of DPSK as the first bit of the information
    x_DPSK(1) = bit(1);
    %Differential encoding of the information
    for iter_DPSK = 1:Sig_Length
        x_DPSK(iter_DPSK + 1) = xor(x_DPSK(iter_DPSK),bit(iter_DPSK));
    end   
    %Modulating the output of differential encoding
    x_DPSK = 2*x_DPSK-1;
    
    %Creating complex AWGN channel
    N_ril = sqrt(No(iter)/2)*randn(1,Sig_Length);
    N_imj = sqrt(No(iter)/2)*randn(1,Sig_Length);
    N = N_ril + 1j*(N_imj);
    
    %Creating complex AWGN channel for DPSK
    N_ril_DPSK = sqrt(No(iter)/2)*randn(1,Sig_Length+1);
    N_imj_DPSK = sqrt(No(iter)/2)*randn(1,Sig_Length+1);
    N_DPSK = N_ril_DPSK + 1j*(N_imj_DPSK);
    
    %Adding AWGN noise to the modulated signal
    Y_BASK = x_BASK + N;
    Y_BFSK = x_BFSK + N;
    Y_BPSK = x_BPSK + N;
    Y_DPSK = x_DPSK + N_DPSK;
    
    %Detector BPSK for DPSK signal
    for iter_DPSK = 1:Sig_Length+1
        if Y_DPSK(iter_DPSK) > 0
            Y_DPSK(iter_DPSK) = 1;
        else
            Y_DPSK(iter_DPSK) = 0;
        end
    end
    
    %Differential decoder for the first bit of DPSK signal
    if Y_DPSK(1) > 0.5
        Z_pre_DPSK(1) = 1;
        Z_DPSK(1) = 1;
    else
        Z_pre_DPSK(1) = 0;
        Z_DPSK(1) = 0;
    end
    
    for iter2 = 1:Sig_Length
        %BASK detector
        Z_BASK(iter2) = (Y_BASK(iter2));
        %Decision circuit BASK
        if (Z_BASK(iter2) > 0.5 && bit(iter2) == 0) || (Z_BASK(iter2) < 0.5 && bit(iter2) == 1);
            E_BASK = E_BASK + 1;
        end
        
        %BFSK detector
        Z_BFSK(iter2) = real(Y_BFSK(iter2))-imag(Y_BFSK(iter2));
        %Decision circuit BFSK
        if (Z_BFSK(iter2) > 0 && bit(iter2) == 0) || (Z_BFSK(iter2) < 0 && bit(iter2) == 1);
            E_BFSK = E_BFSK + 1;
        end
        
        %BPSK detector
        Z_BPSK(iter2) = Y_BPSK(iter2);
        %Decision circuit BPSK
        if (Z_BPSK(iter2) > 0 && bit(iter2) == 0) || (Z_BPSK(iter2) < 0 && bit(iter2) == 1);
            E_BPSK = E_BPSK + 1;
        end
        
        %Differential decoder for the rest of DPSK signal
        if Y_DPSK(iter2+1) > 0.5
            Z_pre_DPSK(iter2+1) = 1;
        else
            Z_pre_DPSK(iter2+1) = 0;
        end
    
        %Detektor diff. encoding
        Z_DPSK(iter2)=(xor(Z_pre_DPSK(iter2),Z_pre_DPSK(iter2+1)));
        %Decision circuit DPSK
        if (Z_DPSK(iter2) > 0.5 && bit(iter2) == 0) || (Z_DPSK(iter2) < 0.5 && bit(iter2) == 1);
            E_DPSK = E_DPSK + 1;
        end
    end
    
    %Simulated BER calculations
    BER_BASK_sim(iter) = E_BASK/Sig_Length;
    BER_BFSK_sim(iter) = E_BFSK/Sig_Length;
    BER_BPSK_sim(iter) = E_BPSK/Sig_Length;
    BER_DPSK_sim(iter) = E_DPSK/Sig_Length;
end

%Making the graph
semilogy(EtoN_dB,BER_BASK_te,'k','color','red');
hold on
semilogy(EtoN_dB_sim,BER_BASK_sim,'k*','color','red');
semilogy(EtoN_dB,BER_BFSK_te,'k','color','green');
semilogy(EtoN_dB_sim,BER_BFSK_sim,'k*','color','green');
semilogy(EtoN_dB,BER_BPSK_te,'k','color','blue');
semilogy(EtoN_dB_sim,BER_BPSK_sim,'k*','color','blue');
semilogy(EtoN_dB,BER_DPSK_te,'k','color','cyan');
semilogy(EtoN_dB_sim,BER_DPSK_sim,'k*','color','cyan');
legend('BASK teori','BASK simulasi','BFSK teori','BFSK simulasi','BPSK teori','BPSK simulasi','DPSK teori','DPSK simulasi','location','best');
axis([min(EtoN_dB) max(EtoN_dB) 10^(-6) 1]);
xlabel('Eb/No (dB)');
ylabel('BER');
title('Grafik BER terhadap Eb/No');
grid on;
hold off