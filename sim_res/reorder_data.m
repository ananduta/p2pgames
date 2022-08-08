%%
%load('sim_A9_inert06-Jun-2022.mat')
ord = [5, 4, 3, 2, 1];
q1 = q;
for i = 1:length(nAg)
    ii = ord(i);
    q.erStd{i} = q1.erStd{ii};
    q.erOvr{i} = q1.erOvr{ii};
    q.erIne{i} = q1.erIne{ii};
end