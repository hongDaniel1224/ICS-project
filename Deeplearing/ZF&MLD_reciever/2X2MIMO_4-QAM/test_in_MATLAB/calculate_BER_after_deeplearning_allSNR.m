% CSV 파일 경로 지정
file_path = 'calculate_allSNR.csv';

% CSV 파일 읽기
variable_data = readmatrix(file_path);

% CSV 파일 경로 지정
file_path = 'correct_allSNR.csv';

% CSV 파일 읽기
correct_data = readmatrix(file_path);

% CSV 파일 경로 지정
file_path = 'my_list_learningat25dB.csv';

% CSV 파일 읽기
num = readmatrix(file_path);

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
EsN0_dB = pow2db(EsN0);

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
ErrorCount_onlyZF= zeros(1, length(EbN0_dB));
ErrorCount_onlyMLD= zeros(1, length(EbN0_dB));
count = 0;
H = zeros(2,2);
y = zeros(2,1);

for iTotal =1:NumberIteration
    tic
    for indx_EsN0 = 1 : length(EsN0)
        data_row = (iTotal-1)*6+indx_EsN0;
        H(1,1) = variable_data(data_row,1)+variable_data(data_row,2)*1j;
        H(1,2) = variable_data(data_row,3)+variable_data(data_row,4)*1j;
        H(2,1) = variable_data(data_row,5)+variable_data(data_row,6)*1j;
        H(2,2) = variable_data(data_row,7)+variable_data(data_row,8)*1j;

        y(1,1) = variable_data(data_row,9)+variable_data(data_row,10)*1j;
        y(2,1) = variable_data(data_row,11)+variable_data(data_row,12)*1j;

        w_ZF = 1 / sqrt(EsN0(indx_EsN0) / Tx) * pinv(H);
        z_ZF = w_ZF * y;
        IntegerSymbol_Detected_ZF = qamdemod(z_ZF * PowerNormalizationFactor_QAM, ModulationOrder);
        ErrorCount_Tmp_onlyZF = sum(correct_data(data_row,1:2) ~= IntegerSymbol_Detected_ZF');
        ErrorCount_onlyZF(1,indx_EsN0) = ErrorCount_onlyZF(1,indx_EsN0) + ErrorCount_Tmp_onlyZF;

        MLDShortDisFirst = norm((y - sqrt(EsN0(indx_EsN0) / Tx) * H * QammodArray(1,:).'),'fro');
        UsedMLD = QammodArray(1,:);
        for MLD_indx = 1 : ModulationOrder^Tx
            MLDShortDisNew = norm((y - sqrt(EsN0(indx_EsN0) / Tx) * H * QammodArray(MLD_indx,:).'),'fro');
            if (MLDShortDisNew <= MLDShortDisFirst)
                MLDShortDisFirst = MLDShortDisNew;
                UsedMLD = QammodArray(MLD_indx,:).';
            end
        end
        IntegerSymbol_Detected_onlyMLD = qamdemod(UsedMLD * PowerNormalizationFactor_QAM, ModulationOrder);
        ErrorCount_Tmp_onlyMLD = sum(correct_data(data_row,1:2) ~= IntegerSymbol_Detected_onlyMLD');
        ErrorCount_onlyMLD(1,indx_EsN0) = ErrorCount_onlyMLD(1,indx_EsN0) + ErrorCount_Tmp_onlyMLD;

        if num(data_row,1) == 0
             % ZF
             w_ZF = 1 / sqrt(EsN0(indx_EsN0) / Tx) * pinv(H);
             z_ZF = w_ZF * y;
             IntegerSymbol_Detected_onlyZF = qamdemod(z_ZF * PowerNormalizationFactor_QAM, ModulationOrder);
             ErrorCount_Tmp_ZF = sum(correct_data(data_row,1:2) ~= IntegerSymbol_Detected_onlyZF');
             ErrorCount_ZF(1,indx_EsN0) = ErrorCount_ZF(1,indx_EsN0) + ErrorCount_Tmp_ZF;
             count = count + 1;
         end

        if num(data_row,1)==1
            % ZF
            w_ZF = 1 / sqrt(EsN0(indx_EsN0) / Tx) * pinv(H);
            z_ZF = w_ZF * y;
            IntegerSymbol_Detected_ZF = qamdemod(z_ZF * PowerNormalizationFactor_QAM, ModulationOrder);
            idx = find(IntegerSymbol_Detected_ZF(1, 1)' == Candidate(:,1));
            QammodArray_dis = QammodArray(idx, :);
            ZF_MLDShortDis = norm((y - sqrt(EsN0(indx_EsN0) / Tx) * H * QammodArray_dis(1,:).'),'fro');
            ZF_MLD_Symbol = QammodArray_dis(1,:);
            for ZF_MLD_indx = 1 : ModulationOrder^(Tx-1)
                ZF_MLDShortDisnew = norm((y - sqrt(EsN0(indx_EsN0) / Tx) * H * QammodArray_dis(ZF_MLD_indx,:).'),'fro');
                if (ZF_MLDShortDisnew <= ZF_MLDShortDis)
                    ZF_MLDShortDis = ZF_MLDShortDisnew;
                    ZF_MLD_Symbol = QammodArray_dis(ZF_MLD_indx,:).';
                end
            end
            IntegerSymbol_Detected_firstZF_secondMLD = qamdemod(ZF_MLD_Symbol * PowerNormalizationFactor_QAM, ModulationOrder);
            ErrorCount_Tmp_ZF = sum(correct_data(data_row,1:2) ~= IntegerSymbol_Detected_firstZF_secondMLD');
            ErrorCount_ZF(1,indx_EsN0) = ErrorCount_ZF(1,indx_EsN0) + ErrorCount_Tmp_ZF;
            count = count + 4;
        end
        if num(data_row,1)==2
            % ZF
            w_ZF = 1 / sqrt(EsN0(indx_EsN0) / Tx) * pinv(H);
            z_ZF = w_ZF * y;
            IntegerSymbol_Detected_ZF = qamdemod(z_ZF * PowerNormalizationFactor_QAM, ModulationOrder);
            idx = find(IntegerSymbol_Detected_ZF(2, 1)' == Candidate(:,2));
            QammodArray_dis = QammodArray(idx, :);
            ZF_MLDShortDis = norm((y - sqrt(EsN0(indx_EsN0) / Tx) * H * QammodArray_dis(1,:).'),'fro');
            ZF_MLD_Symbol = QammodArray_dis(1,:);
            for ZF_MLD_indx = 1 : ModulationOrder^(Tx-1)
                ZF_MLDShortDisnew = norm((y - sqrt(EsN0(indx_EsN0) / Tx) * H * QammodArray_dis(ZF_MLD_indx,:).'),'fro');
                if (ZF_MLDShortDisnew <= ZF_MLDShortDis)
                    ZF_MLDShortDis = ZF_MLDShortDisnew;
                    ZF_MLD_Symbol = QammodArray_dis(ZF_MLD_indx,:).';
                end
            end
            IntegerSymbol_Detected_firstMLD_secondZF = qamdemod(ZF_MLD_Symbol * PowerNormalizationFactor_QAM, ModulationOrder);
            ErrorCount_Tmp_ZF = sum(correct_data(data_row,1:2) ~= IntegerSymbol_Detected_firstMLD_secondZF');
            ErrorCount_ZF(1,indx_EsN0) = ErrorCount_ZF(1,indx_EsN0) + ErrorCount_Tmp_ZF;
            count = count + 4;
        end
        if num(data_row,1)==3
            MLDShortDisFirst = norm((y - sqrt(EsN0(indx_EsN0) / Tx) * H * QammodArray(1,:).'),'fro');
            UsedMLD = QammodArray(1,:);
            for MLD_indx = 1 : ModulationOrder^Tx
                MLDShortDisNew = norm((y - sqrt(EsN0(indx_EsN0) / Tx) * H * QammodArray(MLD_indx,:).'),'fro');
                if (MLDShortDisNew <= MLDShortDisFirst)
                    MLDShortDisFirst = MLDShortDisNew;
                    UsedMLD = QammodArray(MLD_indx,:).';
                end
            end
            IntegerSymbol_Detected_onlyMLD = qamdemod(UsedMLD * PowerNormalizationFactor_QAM, ModulationOrder);
            ErrorCount_Tmp_ZF = sum(correct_data(data_row,1:2) ~= IntegerSymbol_Detected_onlyMLD');
            ErrorCount_ZF(1,indx_EsN0) = ErrorCount_ZF(1,indx_EsN0) + ErrorCount_Tmp_ZF;
            count = count + 16;
        end
    end
    t = toc;
    if mod(iTotal, NumberIteration * 0.01) == 0
                   fprintf('%5.1f초 남았습니다.(%4.1f%% 수행)\n', t * (NumberIteration - iTotal),iTotal * 100 / NumberIteration);
    end
end
ErrorcountZF = [7194,3715,1426,512,163,59];
ErrorcountMLD = [5348,1878,354,48,1,0];

semilogy(EbN0_dB, ErrorcountZF/20000, 'ro-');
hold on
semilogy(EbN0_dB, ErrorcountMLD/20000, 'co-');
semilogy(EbN0_dB, ErrorCount_ZF/20000, 'kpentagram-');

axis([0 25 10^-5 0.5])
grid on
legend('ZF','MLD','ZF-MLD');
xlabel('Eb/No [dB]');
ylabel('SER');
title('2X2 MIMO 4-QAM');