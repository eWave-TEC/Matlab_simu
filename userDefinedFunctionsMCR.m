function userDefinedFunctionsMCR(output, simu, waves, body, ptos, constraint, imcr)
    % Cargar datos del caso actual
    energyExtracted = trapz(output.pto.time, output.pto.powerInternalMechanics);

    % Guardar resultados en un archivo único por caso
    save(sprintf('output_case_%d.mat', imcr), 'energyExtracted', 'output', 'simu', 'waves', 'body', 'ptos', 'constraint');
end
