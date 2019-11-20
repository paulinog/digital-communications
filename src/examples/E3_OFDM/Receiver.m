function data_receive = Receiver(data, NoCarriers, NoPilots, trellis)
%% Receiver of the OFDM system

%%Removing Cyclic Extension
rxed_sig = CyclicPrefixRemove(data);

%% FFT 
ff_sig = fft(rexd_sig,NoCarriers);

%%Pilot Synch
synched_sig = PilotSynch(ff_sig,NoPilots);
% scatterplot(synched_sig);

%%Demodulation 16-QAM/16-PSK
dem_data = qamdemod(synched_sig,16);
% hDemod = comm.PSKDemodulator(M, 'PhasheOffset',pi/16)
% dem_data = step(hDemod, synched_sig');

%%Decimal to binary conversion
bin = de2bi(dem_data','left-msb');
bin = bin';

%%De-Interleaving 
QAMbit = 4; %16-QAM> 2^4=16
deintlvddata = matdeintrlv(bin,2,QAMbit/2); %De-Interleave
deintlvddata = deintlvddata';
deintlvddata = deintlvddata(:)';

%%Deconding data
data_receive = vitdec(deintlvddata,trellis,5,'trunc','hard'); %viterbi decoder
 

