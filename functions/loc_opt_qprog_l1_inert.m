function s = loc_opt_qprog_l1_inert(np,s,t,i)
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
    
    c1 = [];
    
    
    pmg_i_t = np.Smg{i}*s.u_tilde{i}(:,t);
    for h = 1:np.h
        c_1 = c + [0; 0; 0; np.d_mg*(s.sigma_mg_tilde(h,t)-pmg_i_t(h)); zeros(2*sum(np.Adj(i,:)),1)]; % Nash eq.
        %c_1 = c + [0; 0; 0; np.d_mg*(s.sigma_mg_tilde(h,t)); zeros(2*sum(np.Adj(i,:)),1)];            % Wadrop eq.
        c1 = [c1;c_1];
    end
    
    % linear term associated with p_di
    c2 = np.Sdi{i}'*-s.mu_pb_tilde{np.B_n(i)}(:,t);
    
    % linear term associated with p_ch
    c3a = np.Sch{i}'*s.mu_pb_tilde{np.B_n(i)}(:,t);
    
    % linear term associated with p_ds
    c3b = np.Sds{i}'*-s.mu_pb_tilde{np.B_n(i)}(:,t);
    
    c3 = c3a+c3b;
    
    % linear term associated with grid constraints
    c4 = s.lambda_mg_tilde(:,t)'*[np.Smg{i};-np.Smg{i}]+s.mu_tg_tilde(:,t)'*np.Smg{i};
    c4 = c4';
    
    % linear term associated with reciprocity constraints
    c5 = zeros(np.h*(np.no_localDecision+2*sum(np.Adj(i,:))),1);
    for j=1:np.n
        if np.Adj(i,j) == 1
            c5 = c5 + np.Str{i,j}'*s.mu_tr_tilde{i,j}(:,t);
        end     
    end
    %c3 = c3;
    
    % linear term associated with proximal term
    c6 = -np.A{i}*s.u_tilde{i}(:,t);
    
    np.c{i} = c1+ c2 + c3 + c4 + c5 + c6;
    
    
    
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
    s.u{i}(:,t+1) = s.u{i}(:,t);
    a = s.pmg{i}(:,t);
    s.p_mg{i}(:,t+1) = s.p_mg{i}(:,t);
    s.p_di{i}(:,t+1) = s.p_di{i}(:,t);
    s.p_ch{i}(:,t+1) = s.p_ch{i}(:,t);
    s.p_ds{i}(:,t+1) = s.p_ds{i}(:,t);
    for jj = 1:length(np.N{i})
        j = np.N{i}(jj);
        s.p_tr{i,j}(:,t+1) = s.p_tr{i,j}(:,t);
    end
else
    %Assigning the decisions of agent i
        s.u{i}(:,t+1) = u_i;
        s.p_di{i}(:,t+1) = np.Sdi{i}*u_i;
        s.p_ch{i}(:,t+1) = np.Sch{i}*u_i;
        s.p_ds{i}(:,t+1) = np.Sds{i}*u_i;
        s.p_mg{i}(:,t+1) = np.Smg{i}*u_i;
        for jj = 1:length(np.N{i})
            j = np.N{i}(jj);
            s.p_tr{i,j}(:,t+1) = np.Str{i,j}*u_i;
        end

end


end