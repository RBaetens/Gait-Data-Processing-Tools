function CoPraw = cop_from_raw(tabFP)
    % Check if force plate is being stood on
    Th = 20;
    FP1_ON = (abs(tabFP(4:end, 5)) > Th);
    FP2_ON = (abs(tabFP(4:end, 22)) > Th);
    FP3_ON = (abs(tabFP(4:end, 39)) > Th);
    FP4_ON = (abs(tabFP(4:end, 56)) > Th);
    
    % Middles of force plates
    npts = size(FP1_ON, 1);
    middle_FP1 = repmat(tabFP(end, 9:11), npts, 1);
    middle_FP2 = repmat(tabFP(end, 26:28), npts, 1);
    middle_FP3 = repmat(tabFP(end, 43:45), npts, 1);
    middle_FP4 = repmat(tabFP(end, 60:62), npts, 1);
    
    % Rename the raw measurements for clarity
    FzsFP1 = tabFP(4:end, 16:19);
    FzsFP2 = tabFP(4:end, 33:36);
    FzsFP3 = tabFP(4:end, 50:53);
    FzsFP4 = tabFP(4:end, 67:70);
    
    Fx12FP1 = tabFP(4:end, 12);
    Fx34FP1 = tabFP(4:end, 13);
    Fy14FP1 = tabFP(4:end, 14);
    Fy23FP1 = tabFP(4:end, 15);
    
    Fx12FP2 = tabFP(4:end, 29);
    Fx34FP2 = tabFP(4:end, 30);
    Fy14FP2 = tabFP(4:end, 31);
    Fy23FP2 = tabFP(4:end, 32);
    
    Fx12FP3 = tabFP(4:end, 46);
    Fx34FP3 = tabFP(4:end, 47);
    Fy14FP3 = tabFP(4:end, 48);
    Fy23FP3 = tabFP(4:end, 49);
    
    Fx12FP4 = tabFP(4:end, 63);
    Fx34FP4 = tabFP(4:end, 64);
    Fy14FP4 = tabFP(4:end, 65);
    Fy23FP4 = tabFP(4:end, 66);
    
    % Guess for some parameters
    az0 = [-40, -40, -40, -40];
    a = [420/2, 420/2, 420/2, 420/2];
    b = [520/2, 520/2, 520/2, 520/2];
    
    % Variable to store 'raw' CoPs
    CoPraw = zeros(npts, 3*4);
    
    % COP FP1
    ax_ay = cop_xy(az0(1), a(1), b(1), FzsFP1(:, 1), FzsFP1(:, 2), FzsFP1(:, 3), FzsFP1(:, 4), Fx12FP1, Fx34FP1, Fy14FP1, Fy23FP1);
    ax_ay(:, 1) = -ax_ay(:, 1); % not rotated
    ax_ay_az = [ax_ay, zeros(npts, 1)] + middle_FP1;
    CoPraw(:, 1:3) = (ax_ay_az.*FP1_ON) + (middle_FP1.*(~FP1_ON));
    
    % COP FP2
    ax_ay = cop_xy(az0(2), a(2), b(2), FzsFP2(:, 1), FzsFP2(:, 2), FzsFP2(:, 3), FzsFP2(:, 4), Fx12FP2, Fx34FP2, Fy14FP2, Fy23FP2);
    ax_ay(:, 2) = -ax_ay(:, 2); % rotated
    ax_ay_az = [ax_ay, zeros(npts, 1)] + middle_FP2;
    CoPraw(:, 4:6) = (ax_ay_az.*FP2_ON) + (middle_FP2.*(~FP2_ON));
    
    % COP FP3
    ax_ay = cop_xy(az0(3), a(3), b(3), FzsFP3(:, 1), FzsFP3(:, 2), FzsFP3(:, 3), FzsFP3(:, 4), Fx12FP3, Fx34FP3, Fy14FP3, Fy23FP3);
    ax_ay(:, 2) = -ax_ay(:, 2); % rotated
    ax_ay_az = [ax_ay, zeros(npts, 1)] + middle_FP3;
    CoPraw(:, 7:9) = (ax_ay_az.*FP3_ON) + (middle_FP3.*(~FP3_ON));
    
    % COP FP4
    ax_ay = cop_xy(az0(4), a(4), b(4), FzsFP4(:, 1), FzsFP4(:, 2), FzsFP4(:, 3), FzsFP4(:, 4), Fx12FP4, Fx34FP4, Fy14FP4, Fy23FP4);
    ax_ay(:, 2) = -ax_ay(:, 2); % rotated
    ax_ay_az = [ax_ay, zeros(npts, 1)] + middle_FP4;
    CoPraw(:, 10:12) = (ax_ay_az.*FP4_ON) + (middle_FP4.*(~FP4_ON));
end