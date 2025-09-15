function Tzs = Tz_from_FM(tabFP)
    % Get CoPs
    cops = cop_from_FM(tabFP);

    % Disregard nans
    tabFP = tabFP(4:end, 3:end);

    % Rename variables for clarity
    Fxs = [tabFP(:, 1), tabFP(:, 18), tabFP(:, 35), tabFP(:, 52)];
    Fys = [tabFP(:, 2), tabFP(:, 19), tabFP(:, 36), tabFP(:, 53)];
    Mzs = [tabFP(:, 6), tabFP(:, 23), tabFP(:, 40), tabFP(:, 57)];
    
    % Calculate Tzs
    Tzs = [Mzs(:, 1) - Fys(:, 1).*cops(:, 1) + Fxs(:, 1).*cops(:, 2), ...
        Mzs(:, 2) - Fys(:, 2).*cops(:, 4) + Fxs(:, 2).*cops(:, 5), ...
        Mzs(:, 3) - Fys(:, 3).*cops(:, 7) + Fxs(:, 3).*cops(:, 8), ...
        Mzs(:, 4) - Fys(:, 4).*cops(:, 10) + Fxs(:, 4).*cops(:, 11)];
end