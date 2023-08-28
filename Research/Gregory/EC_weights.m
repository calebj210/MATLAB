function w = EC_weights(al, h, z)
    dm = size(z); z = z(:); sq = (0 : length(z) - 1);
    A = rot90(vander(z)); rhs = zeta(al + 1 - sq).*h.^sq;
    w = A \ rhs'; w = reshape(w, dm);
end