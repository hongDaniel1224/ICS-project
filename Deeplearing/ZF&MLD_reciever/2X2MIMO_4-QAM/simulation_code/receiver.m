

% Parameters of Transmitter
Tx                           = 2;
Rx                           = 2;
ModulationOrder              = 4;
nu                           = log2(ModulationOrder);
PowerNormalizationFactor_QAM = sqrt(2/3 .* (ModulationOrder-1));
LengthBitSequence            = nu * Tx;
NumberIteration              = 10^4;

%SNR (Eb/N0)
EbN0_dB = 0 : 5 : 25;
EbN0 = 10 .^ (EbN0_dB / 10); %db2pow(EbNodB)
LengthEbN0 = length(EbN0_dB);
EsN0 = EbN0 * nu;


CandidateSample = zeros(ModulationOrder^Tx,Rx);
Candidate = zeros(ModulationOrder^Tx,Rx);

CandidateSample = 0 : 1 :  ModulationOrder^Rx - 1;
CandidateStr = dec2base(CandidateSample.',ModulationOrder,Rx);

for indx_col = 1 : ModulationOrder^Rx
    for indx_row = 1 : Rx
        if (ModulationOrder == 4)
            Candidate(indx_col, indx_row) = str2double(CandidateStr(indx_col, indx_row));
        end
        if (ModulationOrder == 16)
            Candidate(indx_col, indx_row) = hex2dec(CandidateStr(indx_col, indx_row));
        end
    end
end
QammodArray = qammod(Candidate,ModulationOrder,'InputType','integer') .* (1 ./ sqrt(2/3 .* (ModulationOrder-1)));


ErrorCount_ZF = zeros(1, length(EbN0_dB));
ErrorCount_MMSE = zeros(1, length(EbN0_dB));
ErrorCount_MLD = zeros(1, length(EbN0_dB));
ErrorCount_ZF_MLD = zeros(1, length(EbN0_dB));

ErrorCount_ZF_SER = zeros(1, length(EbN0_dB));
ErrorCount_MLD_SER = zeros(1, length(EbN0_dB));


for iTotal = 1 : NumberIteration
    tic

    BitSequence = randi([0 1], 1, LengthBitSequence);
    
    IntegerSymbol = zeros(Tx, 1);
    QAMSymbol     = zeros(Tx, 1);
    for indx_Symbol = 1 : Tx
        IntegerSymbol(indx_Symbol) = bi2de(BitSequence((indx_Symbol - 1) * nu + 1 : indx_Symbol * nu), 'left-msb');
    end

    QAMSymbol = qammod(IntegerSymbol,ModulationOrder) / PowerNormalizationFactor_QAM;

    Noise = (randn(Rx, 1) + 1j * randn(Rx, 1)) ./ sqrt(2); 
    H = (randn(Tx, Rx) + 1j * randn(Tx, Rx)) ./ sqrt(2);
  
   
%     for indx_EsN0 = 1 : length(EsN0)
%         
%         % Recieve
%         y = sqrt(EsN0(indx_EsN0) / Tx) * H * QAMSymbol + Noise;
%         
%         % ZF
%         w_ZF = 1 / sqrt(EsN0(indx_EsN0) / Tx) * pinv(H);
%         z_ZF = w_ZF * y;
%         IntegerSymbol_Detected_ZF = qamdemod(z_ZF * PowerNormalizationFactor_QAM, ModulationOrder);
%         
%         ErrorCount1(indx_EsN0) = sum(IntegerSymbol_Detected_ZF(1,1) ~= IntegerSymbol(1,1));
%         ErrorCount2(indx_EsN0) = sum(IntegerSymbol_Detected_ZF(2,1) ~= IntegerSymbol(2,1));
%         ErrorCount3(indx_EsN0) = sum(IntegerSymbol_Detected_ZF(3,1) ~= IntegerSymbol(3,1));
%         ErrorCount4(indx_EsN0) = sum(IntegerSymbol_Detected_ZF(4,1) ~= IntegerSymbol(4,1));
% 
%         Error = [sum(ErrorCount1),sum(ErrorCount2),sum(ErrorCount3),sum(ErrorCount4)];
%         [~, index] = min(Error);
%     end
% 
%     idx = find(QammodArray(:, index) == QAMSymbol(index,1));
%     QammodArray_dis = QammodArray(idx, :);

    for indx_EsN0 = 1 : length(EsN0)

        % Recieve
        y = sqrt(EsN0(indx_EsN0) / Tx) * H * QAMSymbol + Noise;

        % ZF
        w_ZF = 1 / sqrt(EsN0(indx_EsN0) / Tx) * pinv(H);
        z_ZF = w_ZF * y;
        IntegerSymbol_Detected_ZF = qamdemod(z_ZF * PowerNormalizationFactor_QAM, ModulationOrder);

%         ZF_MLDShortDis = norm((y - sqrt(EsN0(indx_EsN0) / Tx) * H * QammodArray_dis(1,:).'),'fro');
%         ZF_MLD_Symbol = QammodArray_dis(1,:);
%         for ZF_MLD_indx = 1 : ModulationOrder^(Tx-1)
%             ZF_MLDShortDisnew = norm((y - sqrt(EsN0(indx_EsN0) / Tx) * H * QammodArray_dis(ZF_MLD_indx,:).'),'fro');
%             if (ZF_MLDShortDisnew <= ZF_MLDShortDis)
%                 ZF_MLDShortDis = ZF_MLDShortDisnew;
%                 ZF_MLD_Symbol = QammodArray_dis(ZF_MLD_indx,:).';
%             end
%         end
%         IntegerSymbol_Detected_ZF_MLD = qamdemod(ZF_MLD_Symbol * PowerNormalizationFactor_QAM, ModulationOrder);

        w_MMSE = 1 / sqrt(EsN0(indx_EsN0) / Tx) * pinv(H' * H + Tx / (EsN0(indx_EsN0)) * eye(Tx)) * H';
        z_MMSE = w_MMSE * y;
        IntegerSymbol_Detected_MMSE = qamdemod(z_MMSE * PowerNormalizationFactor_QAM, ModulationOrder);
        
        MLDShortDisFirst = norm((y - sqrt(EsN0(indx_EsN0) / Tx) * H * QammodArray(1,:).'),'fro');
        UsedMLD = QammodArray(1,:);
        for MLD_indx = 1 : ModulationOrder^Tx
            MLDShortDisNew = norm((y - sqrt(EsN0(indx_EsN0) / Tx) * H * QammodArray(MLD_indx,:).'),'fro');
            if (MLDShortDisNew <= MLDShortDisFirst)
                MLDShortDisFirst = MLDShortDisNew;
                UsedMLD = QammodArray(MLD_indx,:).';
            end
        end
        IntegerSymbol_Detected_MLD = qamdemod(UsedMLD * PowerNormalizationFactor_QAM, ModulationOrder);

        BitSequence_Detected_ZF_tmp = de2bi(IntegerSymbol_Detected_ZF, nu, 'left-msb');
        BitSequence_Detected_MMSE_tmp = de2bi(IntegerSymbol_Detected_MMSE, nu, 'left-msb');
        BitSequence_Detected_MLD_tmp = de2bi(IntegerSymbol_Detected_MLD, nu, 'left-msb');
%         BitSequence_Detected_ZF_MLD_tmp = de2bi(IntegerSymbol_Detected_ZF_MLD, nu, 'left-msb');
        
        for indx_Symbol = 1 : Tx
            BitSequence_Detected_ZF((indx_Symbol - 1) * nu + 1 : indx_Symbol * nu) = BitSequence_Detected_ZF_tmp(indx_Symbol, : );
            BitSequence_Detected_MMSE((indx_Symbol - 1) * nu + 1 : indx_Symbol * nu) = BitSequence_Detected_MMSE_tmp(indx_Symbol, : );
            BitSequence_Detected_MLD((indx_Symbol - 1) * nu + 1 : indx_Symbol * nu) = BitSequence_Detected_MLD_tmp(indx_Symbol, : );
%             BitSequence_Detected_ZF_MLD((indx_Symbol - 1) * nu + 1 : indx_Symbol * nu) = BitSequence_Detected_ZF_MLD_tmp(indx_Symbol, : );
        end
        
        ErrorCount_Tmp_ZF_SER = sum(IntegerSymbol_Detected_ZF ~= IntegerSymbol);
        ErrorCount_Tmp_MLD_SER = sum(IntegerSymbol_Detected_MLD ~= IntegerSymbol);
        
        ErrorCount_Tmp_ZF = sum(BitSequence_Detected_ZF ~= BitSequence);
        ErrorCount_Tmp_MMSE = sum(BitSequence_Detected_MMSE ~= BitSequence);
        ErrorCount_Tmp_MLD = sum(BitSequence_Detected_MLD ~= BitSequence);
%         ErrorCount_Tmp_ZF_MLD = sum(BitSequence_Detected_ZF_MLD ~= BitSequence);
        
        ErrorCount_ZF_SER(1, indx_EsN0) = ErrorCount_ZF_SER(1, indx_EsN0) + ErrorCount_Tmp_ZF_SER;
        ErrorCount_MLD_SER(1, indx_EsN0) = ErrorCount_MLD_SER(1, indx_EsN0) + ErrorCount_Tmp_MLD_SER;

        ErrorCount_ZF(1, indx_EsN0) = ErrorCount_ZF(1, indx_EsN0) + ErrorCount_Tmp_ZF;
        ErrorCount_MMSE(1, indx_EsN0) = ErrorCount_MMSE(1, indx_EsN0) + ErrorCount_Tmp_MMSE;
        ErrorCount_MLD(1, indx_EsN0) = ErrorCount_MLD(1, indx_EsN0) + ErrorCount_Tmp_MLD;
%         ErrorCount_ZF_MLD(1, indx_EsN0) = ErrorCount_ZF_MLD(1, indx_EsN0) + ErrorCount_Tmp_ZF_MLD;
    end
end

BER_ZF = ErrorCount_ZF / (LengthBitSequence * NumberIteration);

BER_MMSE = ErrorCount_MMSE / (LengthBitSequence * NumberIteration);

BER_MLD = ErrorCount_MLD /  (LengthBitSequence * NumberIteration);

% BER_ZF_MLD = ErrorCount_ZF_MLD /  (LengthBitSequence * NumberIteration);

% Plot
figure()
semilogy(EbN0_dB, BER_ZF, 'ro-');
hold on
semilogy(EbN0_dB, BER_MMSE, 'go-');
semilogy(EbN0_dB, BER_MLD, 'co-');
% semilogy(EbN0_dB, BER_ZF_MLD, 'ko-');

axis([0 30 10^-5 0.5])
grid on
legend('SimulationZF (Rayleigh)','SimulationMMSE (Rayleigh)','SimulationMLD (Rayleigh)','SimulationZFMLD (Rayleigh)');
xlabel('Eb/No [dB]');
ylabel('ICS BER');
title('ICS BER for (M='+string(ModulationOrder)+')');