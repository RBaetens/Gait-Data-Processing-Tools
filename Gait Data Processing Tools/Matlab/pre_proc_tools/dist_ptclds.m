% sum of distances between corresponding points of two point clouds
function dist = dist_ptclds(ptcld1, ptcld2)
    dist = 0;
    npts = size(ptcld1, 1);
    for i=1:npts
        dist = dist + sqrt(sum((ptcld1(i, :)-ptcld2(i, :)).^2));
    end
end