function cops = cop_from_FM(tabFP)
    % Disregard nans
    tabFP = tabFP(4:end, 3:end);

    % Guess for some parameters
    az0 = 0; % -> Mx and My returned by Nexus are actually Mx' and My'
    
    % Rename variables for clarity
    Fxs = [tabFP(:, 1), tabFP(:, 18), tabFP(:, 35), tabFP(:, 52)];
    Fys = [tabFP(:, 2), tabFP(:, 19), tabFP(:, 36), tabFP(:, 53)];
    Fzs = [tabFP(:, 3), tabFP(:, 20), tabFP(:, 37), tabFP(:, 54)];
    Mxs = [tabFP(:, 4), tabFP(:, 21), tabFP(:, 38), tabFP(:, 55)];
    Mys = [tabFP(:, 5), tabFP(:, 22), tabFP(:, 39), tabFP(:, 56)];
    copzs = [tabFP(:, 9), tabFP(:, 26), tabFP(:, 43), tabFP(:, 60)];

    % Middles of force plates
    npts = size(tabFP, 1);
    middle_FP1 = repmat(tabFP(end, 7:9), npts, 1);
    middle_FP2 = repmat(tabFP(end, 24:26), npts, 1);
    middle_FP3 = repmat(tabFP(end, 41:43), npts, 1);
    middle_FP4 = repmat(tabFP(end, 58:60), npts, 1);
    
    % Calculate intermediate parameters
    Mxs_accent = Mxs + Fys*az0;
    Mys_accent = Mys - Fxs*az0;
    
    % Calculate cop coordinates
    copxs = -Mys_accent./Fzs;
    copys = Mxs_accent./Fzs;
    
    % Bundle and return
    cops = [copxs(:, 1) + middle_FP1(:, 1), copys(:, 1) + middle_FP1(:, 2), copzs(:, 1), ...
        copxs(:, 2) + middle_FP2(:, 1), copys(:, 2) + middle_FP2(:, 2), copzs(:, 2), ...
        copxs(:, 3) + middle_FP3(:, 1), copys(:, 3) + middle_FP3(:, 2), copzs(:, 3), ...
        copxs(:, 4) + middle_FP4(:, 1), copys(:, 4) + middle_FP4(:, 2), copzs(:, 4)];
end