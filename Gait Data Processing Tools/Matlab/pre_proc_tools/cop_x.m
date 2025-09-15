% calculate x-coordinate of cop
function ax = cop_x(az0, a, Fz1, Fz2, Fz3, Fz4, Fx12, Fx34)
    ax = ((Fx12+Fx34)*az0 + a*(Fz1-Fz2-Fz3+Fz4))/(Fz1+Fz2+Fz3+Fz4);
end