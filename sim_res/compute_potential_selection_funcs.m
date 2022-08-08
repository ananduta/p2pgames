%%
load('sim_B_test5.mat')

% for i=1:3
%     [o1{i}] = com_potential(np,o{i});
% end

for k=1:3
    q = o{k};
    phi = 0;
    for i=1:np.b
        for jj = 1:length(np.B{i})
            j = np.B{i}(jj);
            phi = phi + (norm(q.p_l{i,j},2))^2;
            
        end
    end
    phis(k) = phi;
end