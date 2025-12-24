function y = V_compare(V1, V2)
% comparison metric for eigenvectors
k = 0;
for i = 1:2
    for j = 1:2
        k = k + 1;
        dV = [(-1)^i*V2(:,1), (-1)^j*V2(:,2)];
        rho(k) = norm(V1'*dV - eye(2),'fro');
    end
end

V2 = V2(:,[2 1]);
for i = 1:2
    for j = 1:2
        k = k + 1;
        dV = [(-1)^i*V2(:,1), (-1)^j*V2(:,2)];
        rho(k) = norm(V1'*dV - eye(2),'fro');
    end
end

y = min(rho);

end