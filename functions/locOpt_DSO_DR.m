function s = locOpt_DSO_DR(s,a_DSO,np,t)
% Douglas-Rachford splitting to compute optimal decision of DSO
% W. Ananduta
% 13/05/2022
kmax =1e4;
er_max =3e-5;
tic
for y = 1:np.b
    if t==0
        x{y}(:,1) = np.init*ones(np.h*(3+2*length(np.B{y})),1);
    else
        x{y}(:,1) = s.u_DSO{y}(:,t)- inv(np.A_DSO{y})*a_DSO{y};
    end
    xi{y}(:,1) = x{y};
end
xi_all(:,1) = cat(1,xi{:});

%cost function (quadratic term)
A_DSO_all = blkdiag(np.A_DSO{:});
H = np.Spqall'*np.Spqall + 0.5*A_DSO_all;

H = H;
       
Aeq= np.A_eq_dso;
beq= zeros(np.h*np.n_eq_dso,1);
%A = [Aeq;-Aeq];
%b = [beq;-beq];

OMin = pseudoinverse([H Aeq'; Aeq zeros(size(Aeq,1))]);
for k=1:kmax
    
    % project onto C1 (step 1)
    for y=1:np.b
        x_aux = xi{y}(:,k);
        z_aux = zeros(np.h*(3+2*length(np.B{y})),1);
        % project theta (18)
        z_aux(1:np.h,1) = max(min(x_aux(1:np.h,1),np.theta_max(y)),np.theta_min(y));
        
        % project v (19)
        z_aux(np.h+1:2*np.h,1) = max(min(x_aux(np.h+1:2*np.h,1),np.v_max(y)),np.v_min(y));
        
        % project p_tg (20)
        if ~isempty(find(np.B_mg==y))
            z_aux(2*np.h+1:3*np.h,1) = x_aux(2*np.h+1:3*np.h,1);
        %else
        %    z_aux(2*np.h+1:3*np.h,1) = zeros(np.h,1);
        end
        
        % project p_l and q_l (17)
        for zz = 1:length(np.B{y})
            z = np.B{y}(zz);
            p_l = x_aux(np.h*(3+zz-1)+1:np.h*(3+zz),1);
            q_l = x_aux(np.h*(3+length(np.B{y})+zz-1)+1:np.h*(3+length(np.B{y})+zz),1);
            
            for h = 1:np.h
                c_aux = np.s_bar(y,z)/max(norm([p_l(h);q_l(h)],2),np.s_bar(y,z));
                z_aux(np.h*(3+zz-1)+h,1) = c_aux*p_l(h);
                z_aux(np.h*(3+length(np.B{y})+zz-1)+h,1) = c_aux*q_l(h);
            end
        end
        z_b{y} = z_aux;
        %xi_b{y} = xi{y}(:,k);
    end
    
    zh_all(:,k) = cat(1,z_b{:});
    xih_all(:,k) = 2*zh_all(:,k) - xi_all(:,k);
    
    % proximal operation with C2
    % Update xi (step 3) 
    
    

    f = - A_DSO_all*xih_all(:,k); 
    
    omega = OMin*[-f;-beq];
    
    z_all(:,k+1) = omega(1:size(H,1));
    
    %sol_er = norm(Aeq*z_all - beq,inf)
%     % Create an OSQP object
%     prob = osqp;
% 
%     settings = prob.default_settings();
%     settings.eps_abs= 1e-8;
%     settings.eps_rel= 1e-8;
%     settings.max_iter = 1e4;
%     settings.verbose = 0;
%     % Setup workspace and change alpha parameter
%     prob.setup(H, f, A, [], b, settings);
%     res = prob.solve();
%     z_all(:,k+1) = res.x;
    
%     options = optimset('Display','off');
%    [z_all(:,k+1),fval,exitflag] = quadprog(H,f,[],[],Aeq,beq,[],[],[],options);
    %xi_proj = np.Aaug_dso_inv*[2*z_all(:,k)-xi_all(:,k)];% zeros(size(np.A_eq_dso,1),1)];
    
    %xi_proj = (2*z_all(:,k)-xi_all(:,k)) - np.A_eq_dso_A*(2*z_all(:,k)-xi_all(:,k));
    %equation in step 3
    xi_all(:,k+1) = xi_all(:,k) + z_all(:,k+1) - zh_all(:,k);
    
    %assignment for each bus
    for y=1:np.b
        xi{y}(:,k+1) = xi_all(np.h*sum(np.le_u(1:y-1))+1 : np.h*sum(np.le_u(1:y)), k+1);
    end
    
    %stopping criteria
    er(1,k) = norm(xi_all(:,k+1)-xi_all(:,k),inf);
    er(3,k) = norm(xi_all(:,k) - z_all(:,k),inf);
    %er(3,k)
    if k>1
        er(2,k) = norm(z_all(:,k+1)-z_all(:,k),inf);
        
    end
%    if k>1 && er(2,k) < er_max && er(1,k) < er_max
     if k>1 && er(3,k) < er_max && er(2,k) < er_max
        disp('terminate at k:')
        k+1
        break
     end
    error_DRmeth = norm(er(:,k),inf)
    
end
if k == kmax
    disp('max iter reached in projDSO')
end
%locOpt_DSO_time = toc
for y=1:np.b
    s.u_DSO{y}(:,t+1) = xi{y}(:,end);
    s.p_tg{y}(:,t+1) = s.u_DSO{y}(np.h*2+1:np.h*3,t+1);
        for zz = 1:length(np.B{y})
            z = np.B{y}(zz);
            s.p_l{y,z}(:,t+1) =  s.u_DSO{y}(np.h*(3+zz-1)+1:np.h*(3+zz),t+1);
        end
end



