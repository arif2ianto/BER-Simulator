%A script to compare the BER of BPSK through simulation and theory

%Written by: ARIF RAHMADIAN ARIFIANTO - 18117016

%Reference:
%#https://www.mathworks.com/matlabcentral/fileexchange/67649-matlab-code-for-ber-performance-of-binary-ask-digital-modulation
%#https://www.mathworks.com/matlabcentral/fileexchange/28649-ber-of-bfsk-in-awgn-channel
%#https://www.mathworks.com/matlabcentral/fileexchange/66310-ber-for-bpsk-over-awgn-channel-vs-theoretical

%Variables:
%Sig_Length - bit length of the simulated signal
%Eb - bit energy of the simulated signal

%%%%%%%%%%%%%%%%%%%%%
%                   %
%   Task 1 part 1   %
%                   %
%%%%%%%%%%%%%%%%%%%%%
clear all;

Sig_Length = 2000000;
Eb = 1;

%Eb/No for theoretical calculations (high resolution)
EtoN_dB = linspace(0,9,100);
EtoN = 10.^(EtoN_dB/10);

%Eb/No for simulation (five data points)
EtoN_dB_sim = [0, 2, 4, 6, 8];
EtoN_sim = 10.^(EtoN_dB_sim/10);

%No for simulation purposes
No = Eb./EtoN_sim;

%Theoretical BER calculations
BER_BPSK_te = (1/2)*erfc(sqrt(EtoN));

for iter1 = 1:length(EtoN_sim)
    %Initialize error as 0
    E_BPSK = 0;
    
    %Generate random string of bits with the length specified above
    bit = randi([0 1],1,Sig_Length);
    
    %Modulation of BPSK signal
    x_BPSK = 2*bit-1;
    
    %Creating complex AWGN channel
    n_ril = sqrt(No(iter1)/2)*randn(1,Sig_Length);
    n_imj = sqrt(No(iter1)/2)*randn(1,Sig_Length);
    n = n_ril + 1j*(n_imj);
    
    %Adding AWGN noise to the modulated signal
    y_BPSK_eq = x_BPSK + n;
    
    for iter3 = 1:Sig_Length
        %BPSK detector
        z_BPSK(iter3) = real(y_BPSK_eq(iter3));
        %Decision circuit BPSK
        if (z_BPSK(iter3) > 0 && bit(iter3) == 0) || (z_BPSK(iter3) < 0 && bit(iter3) == 1);
            E_BPSK = E_BPSK + 1;
        end
    end
    
    %Simulated BER calculations
    BER_BPSK_sim(iter1) = E_BPSK/Sig_Length;
end

%Making the graph
figure();
semilogy(EtoN_dB,BER_BPSK_te,'k','color','blue');
hold on
legend('BPSK teori','location','best');
axis([min(EtoN_dB) max(EtoN_dB) 10^(-6) 1]);
xlabel('Eb/No (dB)');
ylabel('BER');
title('Grafik BER Teoretis terhadap Eb/No');
grid on;
hold off

figure();
semilogy(EtoN_dB,BER_BPSK_te,'k','color','blue');
hold on
semilogy(EtoN_dB_sim,BER_BPSK_sim,'k*','color','blue');
legend('BPSK teori','BPSK simulasi','location','best');
axis([min(EtoN_dB) max(EtoN_dB) 10^(-6) 1]);
xlabel('Eb/No (dB)');
ylabel('BER');
title('Grafik BER Teoretis dan Simulasi terhadap Eb/No');
grid on;
hold off

%%%%%%%%%%%%%%%%%%%%%
%                   %
%   Task 1 part 2   %
%                   %
%%%%%%%%%%%%%%%%%%%%%
clear all;

FD = 30;
Rb = 128000;

Sig_Length = 2000000;
Eb = 1;

%Eb/No for theoretical calculations (high resolution)
EtoN_dB = linspace(0,26,100);
EtoN = 10.^(EtoN_dB/10);

%Eb/No for simulation (five data points)
EtoN_dB_sim = [0, 5, 10, 15, 20, 25];
EtoN_sim = 10.^(EtoN_dB_sim/10);

%No for simulation purposes
No = Eb./EtoN_sim;

%Theoretical BER calculations
BER_BPSK_te = (1/2)*erfc(sqrt(EtoN));
BER_BPSK_fade_te = (1/2)*(1-sqrt(EtoN./(1+EtoN)));

for iter1 = 1:length(EtoN_sim)
    %Initialize error as 0
    E_BPSK = 0;
    E_BPSK_fade = 0;
    
    %Generate random string of bits with the length specified above
    bit = randi([0 1],1,Sig_Length);
    
    %Modulation of BPSK signal
    x_BPSK = 2*bit-1;
    
    %Creating complex AWGN channel
    n_ril = sqrt(No(iter1)/2)*randn(1,Sig_Length);
    n_imj = sqrt(No(iter1)/2)*randn(1,Sig_Length);
    n = n_ril + 1j*n_imj;
    
    h_ril = 1/sqrt(2)*randn(1,Sig_Length);
    h_imj = 1/sqrt(2)*randn(1,Sig_Length);
    h = h_ril + 1j*h_imj;
    
    %Adding AWGN noise to the modulated signal
    y_BPSK_eq = x_BPSK + n;
    y_BPSK_fade = x_BPSK.*abs(h) + n;
    
    for iter3 = 1:Sig_Length
        %BPSK detector
        z_BPSK(iter3) = real(y_BPSK_eq(iter3));
        %Decision circuit BPSK
        if (z_BPSK(iter3) > 0 && bit(iter3) == 0) || (z_BPSK(iter3) < 0 && bit(iter3) == 1);
            E_BPSK = E_BPSK + 1;
        end
    end
    
    for iter4 = 1:Sig_Length
        %BPSK detector
        z_BPSK_fade(iter4) = real(y_BPSK_fade(iter4));
        %Decision circuit BPSK
        if (z_BPSK_fade(iter4) > 0 && bit(iter4) == 0) || (z_BPSK_fade(iter4) < 0 && bit(iter4) == 1);
            E_BPSK_fade = E_BPSK_fade + 1;
        end
    end
    
    %Simulated BER calculations
    BER_BPSK_sim(iter1) = E_BPSK/Sig_Length;
    BER_BPSK_fade(iter1) = E_BPSK_fade/Sig_Length;
end

%Making the graph
figure();
semilogy(EtoN_dB,BER_BPSK_te,'k','color','blue');
hold on
semilogy(EtoN_dB,BER_BPSK_fade_te,'k','color','green');
legend('AWGN Teori','Rayleigh Fading Teori','location','best');
axis([min(EtoN_dB) max(EtoN_dB) 10^(-6) 1]);
xlabel('Eb/No (dB)');
ylabel('BER');
title('Grafik BER Teoretis terhadap Eb/No');
grid on;
hold off

figure();
semilogy(EtoN_dB,BER_BPSK_te,'k','color','blue');
hold on
semilogy(EtoN_dB,BER_BPSK_fade_te,'k','color','green');
semilogy(EtoN_dB_sim,BER_BPSK_sim,'k*','color','blue');
semilogy(EtoN_dB_sim,BER_BPSK_fade, 'k*','color','green');
legend('AWGN Teori','Rayleigh Fading Teori','AWGN Simulasi','Rayleigh Fading Simulasi','location','best');
axis([min(EtoN_dB) max(EtoN_dB) 10^(-6) 1]);
xlabel('Eb/No (dB)');
ylabel('BER');
title('Grafik BER Teoretis dan Simulasi terhadap Eb/No');
grid on;
hold off

%%%%%%%%%%%%%%%%%%%%%
%                   %
%   Task 2 part 1   %
%                   %
%%%%%%%%%%%%%%%%%%%%%
clear all;

Sig_Length = 100000;

%Make the Rayleigh Fading channel
h = 1/sqrt(2)*(randn(2,Sig_Length) + 1j*randn(2,Sig_Length));

%Calculate the signal power
h_power = h.*conj(h);

%Choose the highest power between two antenna
[hMaxVal ind] = max(h_power,[],1);
hMaxValMat = kron(ones(2,1),hMaxVal);

h_detected = h(h_power==hMaxValMat);
hPowerSel = h_detected.*conj(h_detected);

%Making the graph
figure();
subplot(3,1,1),plot(abs(h_power(1,:)));
title('Daya Sinyal pada Antena 1');
subplot(3,1,2),plot(abs(h_power(2,:)));
title('Daya Sinyal pada Antena 2');
subplot(3,1,3),plot(abs(hPowerSel));
title('Daya Sinyal setelah SDC');

%%%%%%%%%%%%%%%%%%%%%
%                   %
%   Task 2 part 2   %
%                   %
%%%%%%%%%%%%%%%%%%%%%
clear all;

Sig_Length = 100000;
Eb = 1;

EtoN_dB = linspace(0,26,100);
EtoN = 10.^(EtoN_dB/10);

%Eb/No for simulation (five data points)
EtoN_dB_sim = [0, 5, 10, 15, 20, 25];
EtoN_sim = 10.^(EtoN_dB_sim/10);

%No for simulation purposes
No = Eb./EtoN_sim;

%Theoretical BER calculations
BER_BPSK_te = (1/2)*erfc(sqrt(EtoN));
BER_BPSK_fade_te = (1/2)*(1-sqrt(EtoN./(1+EtoN)));
BER_BPSK_sdc_te = 0.5.*(1-2*(1+1./EtoN).^(-0.5) +(1+2./EtoN).^(.5));

Ant = [1 2];

%Generate random bit as transmitted stream
bit = rand(1,Sig_Length)>0.5;

%BPSK Modulation
X_BPSK = 2*bit-1;

for iter1 = 1:length(Ant)
    for iter2 = 1:length(EtoN_dB_sim)
        %Generate AWGN and Rayleigh Fading channel
        n = 1/sqrt(2)*(randn(Ant(iter1),Sig_Length) +   1j*randn(Ant(iter1),Sig_Length)); %white gaussian noise, 0dB variance
        h = 1/sqrt(2)*(randn(Ant(iter1),Sig_Length) + 1j*randn(Ant(iter1),Sig_Length)); % Kanal dan Penambahan Noise
        sD = kron(ones(Ant(iter1),1),X_BPSK);
        y_BPSK = h.*sD + 10^(-EtoN_dB_sim(iter2)/20)*n;
        
        %Get the signal power
        h_power = h.*conj(h);
        
        %Get the maximum power of the channel
        [hMaxVal ind] = max(h_power,[],1);
        hMaxValMat = kron(ones(Ant(iter1),1),hMaxVal);
        
        %Choose the highest power level
        y_BPSK = y_BPSK(h_power==hMaxValMat);
        h_detected = h(h_power==hMaxValMat);
        
        %Equalization
        y_BPSK_eq = y_BPSK./h_detected;
        y_BPSK_eq = reshape(y_BPSK_eq,1,Sig_Length); %to get the matrix dimension proper
        
        %Decoding
        z_BPSK = real(y_BPSK_eq)>0;
        
        %Error counting
        E_BPSK(iter1,iter2) = size(find([bit- z_BPSK]),2);
    end
end
% Simulated BER
BER_BPSK = E_BPSK/Sig_Length;

%Making the Chart
figure
semilogy(EtoN_dB_sim,BER_BPSK(1,:),'k*-','color','red');
hold on;
semilogy(EtoN_dB_sim,BER_BPSK(2,:),'k*-','color','green');
axis([0 27 1e-6 1e0])
grid on
legend('Satu Antena','Dua Antena (SDC)');
title('Grafik BER Normal dan SDC terhadap Eb/No');
xlabel('Eb/No (dB)');
ylabel('Bit Error Rate (BER)');