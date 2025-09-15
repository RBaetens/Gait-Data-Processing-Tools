% checks if point is in rectangle
function check = pnt_in_rect(point, corners)
    check = false;
    xs = corners(:, 1);
    ys = corners(:, 2);
    if (point(1) <= max(xs)) && (point(1) >= min(xs)) && (point(2) <= max(ys)) && (point(2) >= min(ys))
        check = true;
    end
end