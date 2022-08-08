clear all
clc
close all

% load('case_sim_B_03-Jun-2022.mat')
% np = build_mat_exDSO(np);
% 
% for k = 1:3
%     load(['sim_B_',num2str(k),'b.mat'])
% 
%     phi(k) = 0;
%     c=1;
%     for y = 1:np.b
%         for zz = 1:length(np.B{y})
%             z = np.B{y}(zz);
%             pl{k}(y,z) = (norm(r.p_l{y,z},2))^2;
% 
%             r.q_l{y,z} = np.Sql{y,z}*r.u_DSO{y};
%             ql{k}(y,z) = (norm(r.q_l{y,z},2))^2;
% 
%             phi(k) = phi(k) + pl{k}(y,z) + ql{k}(y,z);
%             if z > y
%                 sl(k,c) = pl{k}(y,z) + ql{k}(y,z);
%                 c= c+1;
%             end
%         end
%     end
% end


%load('case_sim_B_06-Jun-2022.mat')
%np = build_mat_exDSO(np);
load('sim_B_test6.mat')
np = build_mat_exDSO(np);
for k = 1:3
    

    phi(k) = 0;
    c=1;
    for y = 1:np.b
        for zz = 1:length(np.B{y})
            z = np.B{y}(zz);
            pl{k}(y,z) = (norm(o{k}.p_l{y,z},2))^2;

            o{k}.q_l{y,z} = np.Sql{y,z}*o{k}.u_DSO{y};
            ql{k}(y,z) = (norm(o{k}.q_l{y,z},2))^2;

            phi(k) = phi(k) + pl{k}(y,z) + ql{k}(y,z);
            if z > y
                sl(k,c) = pl{k}(y,z) + ql{k}(y,z);
                c= c+1;
            end
        end
    end
    o{k}.pmg_all = zeros(np.h,1);
    o{k}.Ptr_t = zeros(np.h,1);
    o{k}.Ptr_b = zeros(np.h,1);
    o{k}.Ptr_e = zeros(np.h,1);
    for i = 1:np.n
        o{k}.pmg_all = o{k}.pmg_all + o{k}.p_mg{i};
        for jj = 1:length(np.N{i})
            j = np.N{i}(jj);
            %if j > i
                %if c==3
                %    Ptr_t = Ptr_t + max(0,o.p_tr{i,j}(:,end));
               % else
                    o{k}.Ptr_t = o{k}.Ptr_t + max(0,o{k}.p_tr{i,j}(:,end));
                %end
            %end
            if np.B_n(i) == np.B_n(j)
                o{k}.Ptr_b = o{k}.Ptr_b + max(0,o{k}.p_tr{i,j}(:,end));
            else
                o{k}.Ptr_e = o{k}.Ptr_e + max(0,o{k}.p_tr{i,j}(:,end));
            end
        end
    end
end


%% Selection function values (normalized)

phin = phi/phi(1);

sl = sqrt(sl');

%%
figure
grid on
box on
hold on
bar(sl)
legend({'baseline','eq. selection','modified game'},'interpreter','Latex','FontSize',12)
xlabel('Line','interpreter','Latex','FontSize',13)
ylabel('Apparent line power [kVA]','interpreter','Latex','FontSize',13)

%%
figure
grid on
box on
hold on
for k=1:3
    plot(o{k}.Ptr_t,'LineWidth',1.5)
end

%% Potential function evaluation

