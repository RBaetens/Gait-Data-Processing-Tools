% stretching + rotation + translation matrix from defining parameters
function map = map_from_params(tx, ty, tz, Rz, ax, ay)
    stretch_mat = eye(4);
    stretch_mat(1, 1) = ax;
    stretch_mat(2, 2) = ay;

    trans_mat = eye(4);
    trans_mat(1, 4) = tx;
    trans_mat(2, 4) = ty;
    trans_mat(3, 4) = tz;
    trans_mat(1, 1) = cos(Rz);
    trans_mat(1, 2) = -sin(Rz);
    trans_mat(2, 1) = sin(Rz);
    trans_mat(2, 2) = cos(Rz);

    map = trans_mat*stretch_mat;
end