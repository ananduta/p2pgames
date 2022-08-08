% Simulations
% extended P2P market, single simulation
% Test inertial and overrelaxed PPP
% W. Ananduta
% 02/05/2022


% clear all
% close all
% clc
% 
% run('pathdef.m')
% rng(240522)
% % Add path of folder 'functions'
% addpath([pwd,'/functions'])
% %addpath([pwd,'/functions/osqp'])
% 
% ty = [0]; %type of case study: (0)heterogenous  or (1)uniform  agents
% tc = [1]; %uniform trading cost
% 
% % set the number of agents
% n_agents = 40;
% n_passive = 50;

% generate case
%run('case_37bus_N.m')

% identify set of neighbors
%np.N = id_neigh(np.Adj);
%np.B = id_neigh(np.Adj_p);


% selections of line capacity constraints
%sb = [((np.n+np.pas_ag)/np.b+2)*600 ((np.n+np.pas_ag)/np.b+2)*300]; % PARAMETERS VARIED

%np.sb_set = sb;

load('case_sim_B_24-May-2022.mat')
cc = 2;
%load('sim_B_hsdm25-May-2022_1.mat')
for i = 1:np.n
    np.u_init{i} = o{cc}.u{i};
    
end
np.lambda_mg_init = o{cc}.lambda_mg;

for i=1:np.n
    for jj=1:length(np.N{i})
        j = np.N{i}(jj);
        np.mu_tr_init{i,j} = o{cc}.mu_tr{i,j} ;
    end
end
for y=1:np.b
    np.u_DSO_init{y} = o{cc}.u_DSO{y};
    np.mu_pb_init{y} = o{cc}.mu_pb{y};
end
np.mu_tg_init = o{cc}.mu_tg;
np.t_st = o{cc}.iter;
save(['case_sim_B_',date],'np')
%% 
for cc = 2%length(sb)
        
    % set line capacity constraint
    %np.s_bar = sb(cc)*ones(np.b);

    
    if cc ==1
        [s1,sl1,np1] = ppp(np);
        % compute total cost
        [s1,o{cc}] = com_cost(s1,np1);     
      %  o{cc}.comp_time = s1.comp_time;
        o{cc}.error = s1.error;
        o{cc}.error_v = s1.error_v;
        
    elseif cc == 2
        [s2,sl2,np2] = ppp_hsdm3(np);
        % compute total cost
        [s2,o{cc}] = com_cost(s2,np2);     
      %  o{cc}.comp_time = s2.comp_time;
        o{cc}.error = s2.error;
        o{cc}.error_v = s2.error_v;
    end
    
    r = o{cc};
    save(['sim_B_hsdm',date,'_',num2str(cc)],'r','o')
    %clearvars('s','sl');

end
save(['sim_B_hsdm',date],'o','np')
