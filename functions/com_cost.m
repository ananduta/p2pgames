function [s,o] = com_cost(s,np)
    s.t = length(s.res1);
    
    
    for i = 1:np.n
        Qh = np.Qh{i};
        ch = np.ch{i};
        z = s.u{i}(:,s.t);
        pmg_i_t = np.Smg{i}*z;
        f_mg = np.d_mg*(s.sigma_mg(:,s.t)'-pmg_i_t')*np.Smg{i}*z;
        %f_mg = np.d_mg*(s.sigma_mg(:,s.t)')*np.Smg{i}*z;
        s.J(i,1) = z'*Qh*z + ch'*z + f_mg;
        s.J_mg(i,1) = f_mg;
        
        
        
        
        % save output (primal and dual variables)

        o.u{i}= z;
        o.J(i,1) = s.J(i,1);
        o.J_mg(i,1) = f_mg;
        
        o.p_di{i}(:,1) = s.p_di{i}(:,end) ;
        o.p_ch{i}(:,1) = s.p_ch{i}(:,end) ;
        o.p_ds{i}(:,1) = s.p_ds{i}(:,end) ;
        o.p_mg{i}(:,1) = s.p_mg{i}(:,end) ;
        
        
        
        for jj = 1:length(np.N{i})
            j = np.N{i}(jj);
            o.p_tr{i,j}(:,1) = s.p_tr{i,j}(:,end) ;
            
            o.mu_tr{i,j}(:,1) = s.mu_tr{i,j}(:,end);
        end
        
        % Compute potential function
        negQ = diag([0 0 0 np.d_mg zeros(1,length(np.N{i})) zeros(1,length(np.N{i}))]);
        negQ = kron(eye(np.h),negQ);
        Qloc = Qh - negQ;
        floc = o.u{i}'*Qloc*o.u{i} + ch'*o.u{i};
        
        D = np.Smg{i}'*np.d_mg*np.Smg{i};
        
        o.P(i) = 0.5*(o.J(i) + floc + o.u{i}'*D*o.u{i});
        
    end
    o.sigma_mg = s.sigma_mg(:,s.t);
    
    for y=1:np.b
        o.u_DSO{y} = s.u_DSO{y}(:,s.t);
        
        o.mu_pb{y} = s.mu_pb{y}(:,end);
        
        for zz=1:length(np.B{y})
            z = np.B{y}(zz);
            o.p_l{y,z}(:,1) = s.p_l{y,z}(:,end);
        end
    end
    o.mu_tg(:,1) = s.mu_tg(:,end);
    o.lambda_mg(:,1) = s.lambda_mg(:,end);
    
    % Compute cost
    
    s.Jt = sum(s.J);
    s.Jt_mg= sum(s.J_mg);
    s.J_pd = np.d_mg*(s.sigma_mg(:,s.t)')*np.sumPd;
    
    o.J = s.J;
    o.Jt = s.Jt;
    o.Jt_mg = s.Jt_mg;
    o.J_pd = s.J_pd;
    o.Jall = o.Jt + o.J_pd;
    o.Pt = sum(o.P);
    o.iter = s.t;
end