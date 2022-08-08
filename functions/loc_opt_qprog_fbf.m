function s = loc_opt_qprog_fbf(np,s,t,i)
% Local optimization (using quadprog)
% W. Ananduta
% 02/12/2019

%for i = 1:np.n
    
    % Construct coefficient of the linear cost
    
    % linear term associated with primal cost function
    %cc = 1;
    c = [np.c_dg(i); np.c_st(i); np.c_st(i); 0; zeros(sum(np.Adj(i,:)),1); np.q_tr*ones(sum(np.Adj(i,:)),1)];
    
    %mu_tr=[];
    for jj = 1:length(np.N{i})
        j = np.N{i}(jj);
        c(np.no_localDecision+jj,1) = np.c_tr(i,j);
        
     %   mu_tr = [mu_tr;sl.mu_tr{i,j}(:,t)];
        
    end
    
    ch = kron(ones(np.h,1),c);
    
    pVar.u = s.u{i}(:,t);
    pVar.p_mg = s.p_mg{i}(:,t);
    pVar.sigma_mg = s.sigma_mg(:,t);
    dVar.mu_pb = s.mu_pb{np.B_n(i)}(:,t);
    dVar.lambda_mg = s.lambda_mg(:,t);
    dVar.mu_tg = s.mu_tg(:,t);
    for jj=1:length(np.N{i})
        j = np.N{i}(jj);
        dVar.mu_tr{i,j} =  s.mu_tr{i,j}(:,t);
    end
            
    nablaL = gradxL(np,pVar,dVar,i);
    
    np.c{i} = ch + nablaL;
    
    
    
%end

A = np.A_ineq{i};
b = np.b_ineq{i};
Aeq = np.Aeq{i};
beq = np.beq{i};

% %quadprog
% H = np.H{i};
% f = np.c{i};
% options = optimset('Display','off');
% %options = optimset('Display','on');
% [u_i,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq,[],[],[],options);

% if exitflag ~= 1 && exitflag ~= 0 %quadprog
%     %disp(st.problem): in case solver cannot find the solution
%    exitflag
%    disp('not solved properly')
%    st.problem
   
%lsqlin
%C = np.H_half{i};
%d = -np.H_half_inv{i}*np.c{i};
%options = optimset('Display','off');
%[u_i,~,~,exitflag] = lsqlin(C,d,A,b,Aeq,beq,[],[],[],options);

%OSQP
H = np.H{i};
f = np.c{i};
A = [A;Aeq;-Aeq];
b = [b;beq;-beq];

% Create an OSQP object
prob = osqp;

settings = prob.default_settings();
settings.eps_abs= 1e-8;
settings.eps_rel= 1e-8;
settings.max_iter = 1e5;
settings.verbose = 0;
% Setup workspace and change alpha parameter
prob.setup(H, f, A, [], b, settings);
res = prob.solve();
u_i = res.x;



if res.info.status_val ~= 1 && res.info.status_val ~= -2 %&& res.info.status_val ~= 2    
    res.info.status_val
    disp('not solved properly')
    pause
    a = s.pmg{i}(:,t);
    s.u_tilde{i}(:,t) = s.u{i}(:,t);
    s.p_mg_tilde{i}(:,t) = s.p_mg{i}(:,t);
    s.p_di_tilde{i}(:,t) = s.p_di{i}(:,t);
    s.p_ch_tilde{i}(:,t) = s.p_ch{i}(:,t);
    s.p_ds_tilde{i}(:,t) = s.p_ds{i}(:,t);
    for jj = 1:length(np.N{i})
        j = np.N{i}(jj);
        s.p_tr_tilde{i,j}(:,t) = s.p_tr{i,j}(:,t);
    end
else
    %Assigning the decisions of agent i
    %s.ph_mg(:,t+1) = 0;
    %dim_prev_ag = 0;
    %for i = 1:np.n
    %    dim_i = np.h*(3+sum(np.Adj(i,:)));
     %   u_i = u_all(dim_prev_ag+1:dim_prev_ag+dim_i);
        s.u_tilde{i}(:,t) = u_i;
        s.p_di_tilde{i}(:,t) = np.Sdi{i}*u_i;
        s.p_ch_tilde{i}(:,t) = np.Sch{i}*u_i;
        s.p_ds_tilde{i}(:,t) = np.Sds{i}*u_i;
        s.p_mg_tilde{i}(:,t) = np.Smg{i}*u_i;
        for jj = 1:length(np.N{i})
            j = np.N{i}(jj);
            s.p_tr_tilde{i,j}(:,t+1) = np.Str{i,j}*u_i;
        end
    %    s.ph_mg(:,t+1) = s.ph_mg(:,t+1) + s.pmg{i}(:,t+1);
    %    dim_prev_ag = dim_prev_ag+dim_i;
    %end
end


end