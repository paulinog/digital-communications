%% Record
info = audiodevinfo;

recObj = audiorecorder
%%
record(recObj, 10);
disp('Start speaking.')
pause(5);
stop(recObj);
disp('End of Recording.');

disp('Playing sound.');
play(recObj);

y = getaudiodata(recObj);
plot(y)