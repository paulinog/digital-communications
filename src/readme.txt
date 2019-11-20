Disciplina IE533A, 2s2019 (FEEC/UNICAMP)
Guilherme Paulino, RA 117119
(guilherme.paulino@msn.com)

=======================
Projeto Final - Parte 1
=======================
Descri��o: Modulador Ac�stico FSK bin�rio

Op��es de execu��o:
	1. Simula��o, usando loopback na placa de som do Windows.
	2. Transmissor e Receptor, usando dois laptops ou um celular, ou outro dispositivo externo, para gravar o som.

==================
Op��o 1: Simulador
==================
Passos:
1. Menu iniciar do Windows, pesquise por "Gerenciador de dispositivos de �udio" ("Manage audio devices").
2. Na Tab "Recording", desabilite o microfone, habilite a placa Stereo Mix e a defina como padr�o.
3. Abra o MATLAB, e rode o script "simulate.m"

==================
Op��o 2: TX/RX
==================
OBS. N�o � necess�rio habilitar o loopback na placa de som. Apenas abra o MATLAB e execute os scripts.

1. Primeiro, execute o script "transmitter.m", e grave o �udio com o dispositivo externo
2. Depois, execute o script "receiver.m" para gravar o �udio e demodular e decodificar a informa��o
