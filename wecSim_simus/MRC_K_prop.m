close all; clearvars; clc;

filename = "Extracted_energy_ss7.mat";
K_init  =   1.0;
k_end   =   100.0;
step    =   0.1;

N       = K_init:step:k_end;
numN    = numel(N);
gain    = zeros(1,numN);
Ener_ext = zeros(1,numN);

for run_ind = 1:numN
    prop_gain = N(run_ind);
    wecSim
    gain(run_ind) = prop_gain;
    Ener_ext(run_ind) = Output_energy.signals.values(end);
end
clearvars -except gain Ener_ext filename
save(filename,"gain", "Ener_ext")