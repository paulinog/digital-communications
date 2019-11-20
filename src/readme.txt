Disciplina IE533A, 2s2019 (FEEC/UNICAMP)
Guilherme Paulino, RA 117119
(guilherme.paulino@msn.com)

=======================
Projeto Final - Parte 1
=======================
Descrição: Modulador Acústico FSK binário

Opções de execução:
	1. Simulação, usando loopback na placa de som do Windows.
	2. Transmissor e Receptor, usando dois laptops ou um celular, ou outro dispositivo externo, para gravar o som.

==================
Opção 1: Simulador
==================
Passos:
1. Menu iniciar do Windows, pesquise por "Gerenciador de dispositivos de áudio" ("Manage audio devices").
2. Na Tab "Recording", desabilite o microfone, habilite a placa Stereo Mix e a defina como padrão.
3. Abra o MATLAB, e rode o script "simulate.m"

==================
Opção 2: TX/RX
==================
OBS. Não é necessário habilitar o loopback na placa de som. Apenas abra o MATLAB e execute os scripts.

1. Primeiro, execute o script "transmitter.m", e grave o áudio com o dispositivo externo
2. Depois, execute o script "receiver.m" para gravar o áudio e demodular e decodificar a informação
