
gain = [];
Energy_e = [];


for i = 0:0.1:42
    prop_gain = i;
    wecSim;++
    Ener = Output_energy.signals.values(end);

    gain = [gain, prop_gain];
    Energy_e = [Energy_e, Ener];
end


T = table(gain', Energy_e', 'VariableNames', {'Prop_gain', 'Energy_extracted'});


writetable(T, 'Energy_extracted.csv')

figure;                           % Crea una nueva ventana de figura
plot(gain, Energy_e, 'b-', 'LineWidth', 2); % Grafica con línea azul y grosor 2
xlabel('Prop Gain');              % Etiqueta del eje X
ylabel('Energy Extracted(J)');       % Etiqueta del eje Y
title('Relación entre Prop Gain y Energy Extracted'); % Título de la gráfica
grid on; 
