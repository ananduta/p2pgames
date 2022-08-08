%%
clear all
clc
close all
ndata = 29;


for i =1:ndata
    load(['sim_A_',num2str(i),'.mat'])
    for j = 1:5
        it{1}(i,j) = length(q.erStd{j});
        it{2}(i,j) = length(q.erOvr{j});
        it{3}(i,j) = length(q.erIne{j});
    end
end

figure
hold on
grid on
box on
t = it{1};
t_a = mean(t);

plot([40:10:80],t_a,'-o','LineWidth',2,'color',[1,0,0])

t = it{2};
t_a = mean(t);

plot([40:10:80],t_a,'-x','LineWidth',2,'color',[0,0,1])

t = it{3};
t_a = mean(t);

plot([40:10:80],t_a,'-^','LineWidth',2,'color',[0,0.6,0])

t = it{1};
t_min = min(t);
t_max = max(t);
[ph,msg]=jbfill([40:10:80],t_min,t_max,[1, 0, 0],[1, 0, 0],0,0.3);


t = it{2};
t_min = min(t);
t_max = max(t);
[ph,msg]=jbfill([40:10:80],t_min,t_max,[0, 0, 1],[0, 0, 1],0,0.3);


t = it{3};
t_min = min(t);
t_max = max(t);
[ph,msg]=jbfill([40:10:80],t_min,t_max,[0, 1, 0],[0, 1, 0],0,0.3);
t_a = mean(t);

t = it{1};
t_a = mean(t);

plot([40:10:80],t_a,'-o','LineWidth',2,'color',[1,0,0])

t = it{2};
t_a = mean(t);

plot([40:10:80],t_a,'-x','LineWidth',2,'color',[0,0,1])

t = it{3};
t_a = mean(t);

plot([40:10:80],t_a,'-^','LineWidth',2,'color',[0,0.6,0])

legend({'Standard','Over-relaxed','Inertial'},'FontSize',12,'Interpreter','Latex')
xticks([40:10:80])
xlabel('Number of prosumers','FontSize',13,'Interpreter','Latex')
ylabel('Number of iterations','FontSize',13,'Interpreter','Latex')


%%
for j = 1:3
    for k=1:ndata
        imp(k,j) = mean((it{j}(k,:)-it{1}(k,:))./it{1}(k,:));
    end
end
mean(imp)
