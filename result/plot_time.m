ionosphere_time = load('ionosphere_time_arr.mat');
ionosphere_time = ionosphere_time.time_arr / 60;

digits1v7_time = load('digits1v7_time_arr.mat');
digits1v7_time = digits1v7_time.time_arr / 60;

ringnorm_time = load('ringnorm_time_arr.mat');
ringnorm_time = ringnorm_time.time_arr / 60;

ionosphere_x = [0.1 0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01];
digits1v7_x = [0.1 0.09 0.08 0.07 0.06 0.05 0.04 0.03];
ringnorm_x = [0.1 0.09 0.08 0.07 0.06 0.05 0.04];

plot(ionosphere_x , ionosphere_time);
hold on
plot(digits1v7_x , digits1v7_time);
plot(ringnorm_x , ringnorm_time);
legend('ionosphere' , 'digits1v7' , 'ringnorm');