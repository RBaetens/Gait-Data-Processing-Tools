% calculate cop for many rows of data
function ax_ay = cop_xy(az0, a, b, Fz1, Fz2, Fz3, Fz4, Fx12, Fx34, Fy14, Fy23)
    npts = size(Fz1, 1);
    ax_ay = zeros(npts, 2);
    for i = 1:npts
        ax_ay(i, 1) = cop_x(az0, a, Fz1(i), Fz2(i), Fz3(i), Fz4(i), Fx12(i), Fx34(i));
        ax_ay(i, 2) = cop_y(az0, b, Fz1(i), Fz2(i), Fz3(i), Fz4(i), Fy14(i), Fy23(i));
    end
end