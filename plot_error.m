% %%
% figure
% hold on, grid on, box on
% for l = 1:3
%     plot(o{l}.error,'LineWidth',1.5)
% end
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% %label(['standard','overrelaxed','inertial'])



%%
figure
hold on, grid on, box on
%for l = 1:3
%    plot(o{l}.error,'LineWidth',1.5)
%end
i = 1
plot(q.erStd{i},'LineWidth',1.5)
plot(q.erOvr{i},'LineWidth',1.5)
plot(q.erIne{i},'LineWidth',1.5)
set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
legend('standard','overrelaxed','inertial')