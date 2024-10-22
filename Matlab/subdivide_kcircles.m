function new_circles = subdivide_kcircles(z0,r0,k)

    % function to divide initial circle into smaller circles
    con_circle = [z0, r0/2];
    r0 = 5*r0/12;
    k_all = (0:1:k-1).';
    z0 = z0 + (3/4)*r0*exp(2*pi*1i.*k_all/k);
    new_circles = [con_circle; ...
                    z0, ones(5,1).*r0];
end