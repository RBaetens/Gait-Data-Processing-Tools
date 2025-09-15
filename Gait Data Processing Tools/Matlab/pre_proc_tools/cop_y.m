% calculate y-coordinate of cop
function ay = cop_y(az0, b, Fz1, Fz2, Fz3, Fz4, Fy14, Fy23)
    ay = ((Fy14+Fy23)*az0 + b*(Fz1+Fz2-Fz3-Fz4))/(Fz1+Fz2+Fz3+Fz4);
end