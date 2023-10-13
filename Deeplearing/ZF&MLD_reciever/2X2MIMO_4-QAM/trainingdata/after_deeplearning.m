close all
clear
clc

% Parameters of Transmitter
Tx                           = 2;
Rx                           = 2;
ModulationOrder              = 4;
nu                           = log2(ModulationOrder);
PowerNormalizationFactor_QAM = sqrt(2/3 .* (ModulationOrder-1));
LengthBitSequence            = nu * Tx;
NumberIteration              = 10^7;

%SNR (Eb/N0)
EbN0_dB = 0 : 5 : 25;
EbN0 = 10 .^ (EbN0_dB / 10); %db2pow(EbNodB)
LengthEbN0 = length(EbN0_dB);
EsN0 = EbN0 * nu;
EsN0_dB = pow2db(EsN0);

ErrorCount_ZF = zeros(1, length(EbN0_dB));

count1=0;

for iTotal = 1 : NumberIteration
    tic

    BitSequence = randi([0 1], 1, LengthBitSequence);
    BitSequence_Save(iTotal,:) = BitSequence;
    
    IntegerSymbol = zeros(Tx, 1);
    QAMSymbol     = zeros(Tx, 1);
    for indx_Symbol = 1 : Tx
        IntegerSymbol(indx_Symbol) = bi2de(BitSequence((indx_Symbol - 1) * nu + 1 : indx_Symbol * nu), 'left-msb');
    end
    
    QAMSymbol = qammod(IntegerSymbol,ModulationOrder) / PowerNormalizationFactor_QAM;
    
            
    Noise = (randn(Rx, 1) + 1j * randn(Rx, 1)) ./ sqrt(2); 
    H = (randn(Tx, Rx) + 1j * randn(Tx, Rx)) ./ sqrt(2);
    condition_number = cond(H);

    for indx_EsN0 = 1 : length(EsN0)

        % Recieve
        y = sqrt(EsN0(indx_EsN0) / Tx) * H * QAMSymbol + Noise;

        % ZF
        w_ZF = 1 / sqrt(EsN0(indx_EsN0) / Tx) * pinv(H);
        z_ZF = w_ZF * y;
        IntegerSymbol_Detected_ZF = qamdemod(z_ZF * PowerNormalizationFactor_QAM, ModulationOrder);
        dis = abs(z_ZF * PowerNormalizationFactor_QAM - w_ZF * H * IntegerSymbol_Detected_ZF);

        BitSequence_Detected_ZF_tmp = de2bi(IntegerSymbol_Detected_ZF, nu, 'left-msb');
        
        for indx_Symbol = 1 : Tx
            BitSequence_Detected_ZF((indx_Symbol - 1) * nu + 1 : indx_Symbol * nu) = BitSequence_Detected_ZF_tmp(indx_Symbol, : );
        end
        
        ErrorCount_Tmp_ZF = sum(IntegerSymbol ~= IntegerSymbol_Detected_ZF);
        
        ErrorCount_Tmp_ZF_1 = sum(IntegerSymbol(1,1) ~= IntegerSymbol_Detected_ZF(1,1));
        ErrorCount_Tmp_ZF_2 = sum(IntegerSymbol(2,1) ~= IntegerSymbol_Detected_ZF(2,1));
        
        %두 심볼 다 맞았을 때
        if ErrorCount_Tmp_ZF_1 == 0 && ErrorCount_Tmp_ZF_2 == 0 
            num3 = 0;
        end
        % y1 심볼만 맞았을 때
        if ErrorCount_Tmp_ZF_1 == 0 && ErrorCount_Tmp_ZF_2 == 1 
            num3 = 1;
        end
        % y2 심볼만 맞았을 때
        if ErrorCount_Tmp_ZF_1 == 1 && ErrorCount_Tmp_ZF_2 == 0 
            num3 = 2;
        end
        %두 심볼 다 틀렸을 때
        if ErrorCount_Tmp_ZF_1 == 1 && ErrorCount_Tmp_ZF_2 == 1 
            num3 = 3;
        end
        if count1 < 1000 && EbN0_dB(indx_EsN0)==0 && num3==3
            learning = [real(H(1,1)), imag(H(1,1)),real(H(1,2)), imag(H(1,2)),real(H(2,1)), imag(H(2,1)),real(H(2,2)), imag(H(2,2)) ,real(y(1,1)), imag(y(1,1)),real(y(2,1)), imag(y(2,1)), real(w_ZF(1,1)), imag(w_ZF(1,1)),real(w_ZF(1,2)), imag(w_ZF(1,2)),real(w_ZF(2,1)), imag(w_ZF(2,1)),real(w_ZF(2,2)), imag(w_ZF(2,2)),dis(1,1),dis(2,1),condition_number];
            filename = sprintf('learning_EsN0(0dB).csv');
            writematrix(learning,filename, 'WriteMode', 'append');
            
            filename = sprintf('datalabel_EsN0(0dB).csv');
            writematrix(num3,filename, 'WriteMode', 'append');
%             after_learing = [real(H(1,1)), imag(H(1,1)),real(H(1,2)), imag(H(1,2)),real(H(2,1)), imag(H(2,1)),real(H(2,2)), imag(H(2,2)) ,real(y(1,1)), imag(y(1,1)),real(y(2,1)), imag(y(2,1))];
%             filename = sprintf('calculate_EsN0.csv');
%             writematrix(after_learing,filename, 'WriteMode', 'append');
%  
%             filename = sprintf('correct_EsN0.csv');
%             writematrix(IntegerSymbol',filename, 'WriteMode', 'append');
             count1 = count1 + 1;
             fprintf('%d\n',count1);
        end
        if count1 > 999
            exit;
        end
       
        ErrorCount_ZF(1, indx_EsN0) = ErrorCount_ZF(1, indx_EsN0) + ErrorCount_Tmp_ZF;
    end
    t = toc;
    if mod(iTotal, NumberIteration * 0.01) == 0
                   fprintf('%5.1f초 남았습니다.(%4.1f%% 수행)\n', t * (NumberIteration - iTotal),iTotal * 100 / NumberIteration);
    end
end

BER_ZF = ErrorCount_ZF / (2* NumberIteration);

% Plot
figure()
semilogy(EbN0_dB, BER_ZF, 'ro-');

axis([0 30 10^-5 0.5])
grid on
legend('Simulation ZF (Rayleigh)');
xlabel('Eb/No [dB]');
ylabel('ICS BER');
title('ICS BER for (M='+string(ModulationOrder)+')');