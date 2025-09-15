% inversion of transformation matrix
function inv_T = trans_mat_inv(T)
    R = T(1:3, 1:3);
    t = T(1:3, 4);
    inv_T = zeros(4, 4);
    inv_T(1:3, 1:3) = transpose(R);
    inv_T(1:3, 4) = -1*transpose(R)*t;
    inv_T(4, 4) = 1;
end