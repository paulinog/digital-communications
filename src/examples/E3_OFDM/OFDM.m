%Single frame size: 96 bits
%Total no. of frames: 100
%coding used: Convulutinal coding
%Modulation: 16-QAM
%NO. of Pilots 4
%No of carries 64
%Cylic Extension: 25% (16)
close all;
clear all;
clc;

%%Generating and coding data
t_data = randi ([0 1],1,9600);
NoPilots = 4;
BitInFrame = 1;
BERrow = 1;

%% TX

for d = 1:100
    data = t_data(BitInFrame:BitInFrame+95);
    BitInFrame = BitInFrame+96;
    trellis = poly2trellis(7, [171 133]); % used for convulationally enconde/decode
    % 7 is the contrain legth of the convulational code; [171 133] is the
    % polynomial the input/output relation in OCT
    
    [cext_data NoCarriers] = Transmitter(data, NoPilots, trellis);
    
    %% Channel % SNR
    o=1;
    for snr = 0:2:50
        ofdm_sig = awgn(cext_data,snr,'measured'); %add AWGN
        figure;
        index=1:80;
        plot(index,cext_data,'b',index,ofdm_sig, 'r'); %plot both signals
        legend('Original Signal to be transmitted','Signal with AWGN');
       
        %% RX
        
        rxed_data = Receiver(ofdm_sig, NoCarriers, NoPilots, trellis);
        
        %% Calculating BER
        
        rxed_data = rxed_data(:)';
        c = xor(data,rxed_data);
        erros = nnz(c);
        figure;
        %subplot(311); plot(1:96,data); title('original signal');
        %subplot(312); plot(1:96,rxed_data); title('received signal');
        %subplot(313); plot(1:96,data, '--',1:96, rxed_data, ':');
        %legend('original , 'received','comparision');
        
        BER(Berrow,o) = erros;length(data);
        
    end % SNR loop ends here
    BERrow = BERrow+1;
end % main data loop

%% Time averaging for optimum results

for col = 1:25;  % change if SNR loop changed
    ber(1,col) = 0;
    for row = 1:100; 
        ber(1,col) = ber(1,col)+BER(row,col);
    end
    
end

ber = ber./100;
%%
figure;
i = 0:2:48; 
semilogy(i, ber);
title('BER VS SNR');
ylabel('BER');
xlabel('SNR (dB)');
grid on

    
        