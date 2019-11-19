clc;
%% Poly to Trellis
% Convert convolutional code polynomial to trellis description

trellis = poly2trellis(3, [5 7]);

%% Tx
tam=1000;
msg=randsrc(1, tam, [0 1]);

msgENC = convenc(msg, trellis);

%% Rx
msgDEC = vitdec(msgENC, trellis, 20,  'trunc', 'hard');

err = sum(msg(1:end)~=msgDEC(1:end));
if (err > 0)
    error([num2str(err) ' symbol errors'])
else
    disp('no errors')
end