function psi = mode_rot(V)
psi_p{1} = [atan2(V(2,1),V(1,1)), atan2(-V(2,1),-V(1,1))];
psi_p{2} = [atan2(V(2,2),V(1,2)), atan2(-V(2,2),-V(1,2))];
psi_q{1} = [atan2(-V(1,1),V(2,1)), atan2(V(1,1),-V(2,1))];
psi_q{2} = [atan2(-V(1,2),V(2,2)), atan2(V(1,2),-V(2,2))];

[~,ind11] = min(abs(psi_p{1}));
psi11 = psi_p{1}(ind11);
[~,ind12] = min(abs(psi_p{2}));
psi12 = psi_p{2}(ind12);
[~,ind21] = min(abs(psi_q{1}));
psi21 = psi_q{1}(ind21);
[~,ind22] = min(abs(psi_q{2}));
psi22 = psi_q{2}(ind22);
if abs(psi11)<abs(psi21)
    psi1 = psi11;
    psi2 = psi22;
else
    psi1 = psi12;
    psi2 = psi21;
end
psi = psi1;

end