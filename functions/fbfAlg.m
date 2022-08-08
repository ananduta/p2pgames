function [s,sl,np] = fbfAlg(np)
%% Distributed algorithm
% extended P2P market
% W. Ananduta
% 12/11/2020

        %% INITIALIZATION
        % Assigning parameters of the algorithm
        %np = alg_param(np);
        %np = alg_param_37b_l2(np); 
        np = alg_param_37b_fbf(np);
        np.t_max = 1e4;
        np.er_max = 1e-2; 

        

        % Generate matrices for cost function and constraints       
        np = build_mat_exP2P_tr1(np);
        np = build_mat_exDSO(np);
        
        %% Initialization of the variables
        t=1;
        s.sigma_mg = np.sumPd;
        for i=1:np.n
            % decision variables
            %s = initialize_u_quadp(np,s,i);
            s.u{i}(:,1) = np.init*ones(size(np.A_ineq{i},2),1);
            s.p_di{i}(:,1) = np.Sdi{i}*s.u{i}(:,1);
            s.p_ch{i}(:,1) = np.Sch{i}*s.u{i}(:,1);
            s.p_ds{i}(:,1) = np.Sds{i}*s.u{i}(:,1);
            s.p_mg{i}(:,1) = np.Smg{i}*s.u{i}(:,1);
            for jj = 1:length(np.N{i})
                j = np.N{i}(jj);
                s.p_tr{i,j} = np.Str{i,j}*s.u{i}(:,1);
            end
            sl.b{i}(:,1) = np.Pd(i,1:np.h)' - s.p_di{i}(:,1) - s.p_ds{i}(:,1)+ s.p_ch{i}(:,1);

            s.sigma_mg(:,1) = s.sigma_mg(:,1) + s.p_mg{i}(:,1);
        end
        s.lambda_mg= zeros(2*np.h,1);

        for i=1:np.n
            for jj=1:length(np.N{i})
                j = np.N{i}(jj);
                sl.c_tr{i,j}(:,1) = s.p_tr{i,j}(:,1) + s.p_tr{j,i}(:,1);
                
                s.mu_tr{i,j}(:,1) = np.init*ones(np.h,1);
            end
        end

        % Initialization of DSO's decision variables         
        s = projDSO(s,0,np,0);                                     % 
        
        %s.sum_p_pd = zeros(np.h,1);
        
        sl.b_DSO(:,1) = kron([1;-1],s.sigma_mg(:,1))-[np.pmg_max*ones(np.h,1);-np.pmg_min*ones(np.h,1)];
        
        for y = 1:np.b
            
            sl.c_pb{y}(:,1) =  - s.p_tg{y}(:,1);
            for ii = 1:length(np.Pasag_b{y})
                    i = np.Pasag_b{y}(ii);
                    sl.c_pb{y}(:,1) = sl.c_pb{y}(:,1) + np.Pd(i+np.n,:)';
            end
            
            for ii = 1:length(np.N_b{y})
                i = np.N_b{y}(ii);
                sl.c_pb{y}(:,1) = sl.c_pb{y}(:,1) + sl.b{i}(:,1);
            end

            for zz = 1:length(np.B{y})
                z = np.B{y}(zz);
                sl.c_pb{y}(:,1) = sl.c_pb{y}(:,1) - s.p_l{y,z}(:,1);
            end

            s.mu_pb{y}(:,1) = np.init*ones(np.h,1);

            %s.sum_p_pd = s.sum_p_pd + np.Pd(np.n+y,1:np.h)';

        end
        %s.sigma_mg(:,1) = s.sigma_mg(:,1) + s.sum_p_pd;

        sl.c_tg(:,1) = s.sigma_mg(:,1);
        for yy = 1:length(np.B_mg)
            y = np.B_mg(yy);
            sl.c_tg(:,1) = sl.c_tg(:,1) - s.p_tg{y}(:,1);                
        end
        s.mu_tg(:,1) = np.init*ones(np.h,1);



        %% Iteration
        
        for t = 1:np.t_max
            %tic
            t
            % 1) Strategy update
            s.sigma_mg_tilde(:,t) = np.sumPd;
            
            
           for i=1:np.n
                
                % primal update of prosumer (quadratic prog.)
                tic
                
                s = loc_opt_qprog_fbf(np,s,t,i);
                s.comp_p(i) = toc;
                % local load imbalance of prosumer i
                sl.b_tilde{i}(:,t) = np.Pd(i,1:np.h)' - s.p_di_tilde{i}(:,t) - s.p_ds_tilde{i}(:,t) + s.p_ch_tilde{i}(:,t);

                % forward p_mg to DSO
                s.sigma_mg_tilde(:,t) = s.sigma_mg_tilde(:,t) + s.p_mg_tilde{i}(:,t);
                
                % compute error

                res1(i,t+1)= norm(s.u_tilde{i}(:,t)-s.u{i}(:,t),inf);
            end

            res2(t) = 0; % Extra: compute norm of residual
            dres2(t) = 0;

            % dual variables (reciprocity constraint) Prosumers
            for i=1:np.n

                for j=1:np.n
                    if np.Adj(i,j) == 1

                        % Dual variables (reciprocity constraints) update:
                        
                        % Aux vector
                        sl.c_tr_tilde{i,j}(:,t) = s.p_tr_tilde{i,j}(:,t) + s.p_tr_tilde{j,i}(:,t);
                        
                        % Reflected dual ascent
                        s.mu_tr_tilde{i,j}(:,t) = s.mu_tr{i,j}(:,t) + np.beta_tr(i,j)*(sl.c_tr{i,j}(:,t));

                        % Extra: compute norm of residual
                        res = sl.c_tr_tilde{i,j}(:,t);
                        res2(t) = norm([res2(t);res],inf);

                    end
                end
            end


            % DSO=========================================================================================================================================
            % Primal update 
            
            for y=1:np.b
                a_DSO{y} = [zeros(2*np.h,1); -s.mu_tg(:,t)-s.mu_pb{y}(:,t);kron(ones(length(np.B{y}),1),-s.mu_pb{y}(:,t));zeros(np.h*length(np.B{y}),1)]; 
            end
            tic
            s = projDSO_fbf(s,a_DSO,np,t);                                         
            s.comp_DNO = toc;                       
            
            % Dual variable (grid constraints)
            sl.b_DSO_tilde(:,t) = kron([1;-1],s.sigma_mg_tilde(:,t))-[np.pmg_max*ones(np.h,1);-np.pmg_min*ones(np.h,1)];
            
            s.lambda_mg_tilde(:,t) = max(0, s.lambda_mg(:,t) + np.gamma_mg*sl.b_DSO(:,t));        

            % Dual variable (local power balance of busses)
            res_dso2(t)=0;
            dres_dso2(t)=0;
            for y=1:np.b
                sl.c_pb_tilde{y}(:,t) = - s.p_tg_tilde{y}(:,t);
                for ii = 1:length(np.Pasag_b{y})
                    i = np.Pasag_b{y}(ii);
                    sl.c_pb_tilde{y}(:,t) = sl.c_pb_tilde{y}(:,t) + np.Pd(i+np.n,:)';
                end
                for ii = 1:length(np.N_b{y})
                    i = np.N_b{y}(ii);
                    sl.c_pb_tilde{y}(:,t) = sl.c_pb_tilde{y}(:,t) + sl.b_tilde{i}(:,t);
                end

                for zz = 1:length(np.B{y})
                    z = np.B{y}(zz);
                    sl.c_pb_tilde{y}(:,t) = sl.c_pb_tilde{y}(:,t) - s.p_l_tilde{y,z}(:,t);
                end
                
                
                s.mu_pb_tilde{y}(:,t) = s.mu_pb{y}(:,t) + np.beta_pb(y)*(sl.c_pb{y}(:,t) ) ;
                
                % error
                res_dso2(t) = norm([res_dso2(t);sl.c_pb{y}(:,t)],inf);
                
                
                res1(np.n+y,t+1)= norm(s.u_DSO_tilde{y}(:,t)-s.u_DSO{y}(:,t),inf);
            end

            % Dual variable (trading with the main grid)
            s.sigma_tg_tilde(:,t) = zeros(np.h,1) ;
            for yy = 1:length(np.B_mg)
                y = np.B_mg(yy);
                s.sigma_tg_tilde(:,t) = s.sigma_tg_tilde(:,t) + s.p_tg_tilde{y}(:,t);                
            end
            sl.c_tg_tilde(:, t) = s.sigma_mg_tilde(:,t) - s.sigma_tg_tilde(:,t);
            
            s.mu_tg_tilde(:,t) = s.mu_tg(:,t) + np.beta_tg*(sl.c_tg(:,t));
            
            res_gridTrading = norm(sl.c_tg(:,t),inf);
            
            
            
            % Second Forward steps
            s.sigma_mg(:,t+1) = np.sumPd;
            
            
            % prosumers
            for i = 1:np.n
                
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

                nablaLi = gradxL(np,pVar,dVar,i);
                
                
                pVar.u = s.u_tilde{i}(:,t);
                pVar.p_mg = s.p_mg_tilde{i}(:,t);
                pVar.sigma_mg = s.sigma_mg_tilde(:,t);
                dVar.mu_pb = s.mu_pb_tilde{np.B_n(i)}(:,t);
                dVar.lambda_mg = s.lambda_mg_tilde(:,t);
                dVar.mu_tg = s.mu_tg_tilde(:,t);
                for jj=1:length(np.N{i})
                    j = np.N{i}(jj);
                    dVar.mu_tr{i,j} =  s.mu_tr_tilde{i,j}(:,t);
                end
                
                nablaLi_tilde = gradxL(np,pVar,dVar,i);
                
                s.u{i}(:,t+1) = s.u_tilde{i}(:,t) - np.gamma_mg*(nablaLi_tilde - nablaLi);
                s.p_di{i}(:,t+1) = np.Sdi{i}*s.u{i}(:,t+1);
                s.p_ch{i}(:,t+1) = np.Sch{i}*s.u{i}(:,t+1);
                s.p_ds{i}(:,t+1) = np.Sds{i}*s.u{i}(:,t+1);
                s.p_mg{i}(:,t+1) = np.Smg{i}*s.u{i}(:,t+1);
                for jj = 1:length(np.N{i})
                    j = np.N{i}(jj);
                    s.p_tr{i,j}(:,t+1) = np.Str{i,j}*s.u{i}(:,t+1);
                end
                
                % local load imbalance of prosumer i
                sl.b{i}(:,t+1) = np.Pd(i,1:np.h)' - s.p_di{i}(:,t+1) - s.p_ds{i}(:,t+1) + s.p_ch{i}(:,t+1);

                % forward p_mg to DSO
                s.sigma_mg(:,t+1) = s.sigma_mg(:,t+1) + s.p_mg{i}(:,t+1);
                
                % compute error

                res1a(i,t+1)= norm(s.u{i}(:,t+1)-s.u{i}(:,t),inf);
                                
            end
            
            for i=1:np.n

                for j=1:np.n
                    if np.Adj(i,j) == 1

                        % Dual variables (reciprocity constraints) update:
                        
                        % Aux vector
                        sl.c_tr{i,j}(:,t+1) = s.p_tr{i,j}(:,t+1) + s.p_tr{j,i}(:,t+1);
                        
                        
                        % dual ascent
                        s.mu_tr{i,j}(:,t+1) = s.mu_tr_tilde{i,j}(:,t) + np.beta_tr(i,j)*(sl.c_tr_tilde{i,j}(:,t)-sl.c_tr{i,j}(:,t));
                        
                        % Extra: compute norm of residual
                        res = sl.c_tr{i,j}(:,t+1);
                        res2a(t) = norm([res2(t);res],inf);

                    end
                end
            end
            
            % DSO
            for y=1:np.b
                nablaL_DSO = [zeros(2*np.h,1); -s.mu_tg(:,t)-s.mu_pb{y}(:,t);kron(ones(length(np.B{y}),1),-s.mu_pb{y}(:,t));zeros(np.h*length(np.B{y}),1)];
                nablaL_tilde_DSO = [zeros(2*np.h,1); -s.mu_tg_tilde(:,t)-s.mu_pb_tilde{y}(:,t);kron(ones(length(np.B{y}),1),-s.mu_pb_tilde{y}(:,t));zeros(np.h*length(np.B{y}),1)];
                
                s.u_DSO{y}(:,t+1) = s.u_DSO_tilde{y}(:,t) - np.gamma_mg*(nablaL_tilde_DSO - nablaL_DSO ) ;
                
                s.p_tg{y}(:,t+1) = s.u_DSO{y}(np.h*2+1:np.h*3,t+1);
                for zz = 1:length(np.B{y})
                    z = np.B{y}(zz);
                    s.p_l{y,z}(:,t+1) =  s.u_DSO{y}(np.h*(3+zz-1)+1:np.h*(3+zz),t+1);
                    s.q_l{y,z}(:,t+1) =  s.u_DSO{y}(np.h*(3+length(np.B{y})+zz-1)+1:np.h*(3+length(np.B{y})+zz),t+1);
                end
            end
            
            % Dual variable (grid constraints)
            sl.b_DSO(:,t+1) = kron([1;-1],s.sigma_mg(:,t+1))-[np.pmg_max*ones(np.h,1);-np.pmg_min*ones(np.h,1)];
            s.lambda_mg(:,t+1) = max(0, s.lambda_mg_tilde(:,t) + np.gamma_mg*( sl.b_DSO_tilde(:,t)-sl.b_DSO(:,t) ));        

            % Dual variable (local power balance of busses)
            res_dso2(t)=0;
            dres_dso2(t)=0;
            for y=1:np.b
                sl.c_pb{y}(:,t+1) = - s.p_tg{y}(:,t+1);
                for ii = 1:length(np.Pasag_b{y})
                    i = np.Pasag_b{y}(ii);
                    sl.c_pb{y}(:,t+1) = sl.c_pb{y}(:,t+1) + np.Pd(i+np.n,:)';
                end
                for ii = 1:length(np.N_b{y})
                    i = np.N_b{y}(ii);
                    sl.c_pb{y}(:,t+1) = sl.c_pb{y}(:,t+1) + sl.b{i}(:,t+1);
                end

                for zz = 1:length(np.B{y})
                    z = np.B{y}(zz);
                    sl.c_pb{y}(:,t+1) = sl.c_pb{y}(:,t+1) - s.p_l{y,z}(:,t+1);
                end
                s.mu_pb{y}(:,t+1) = s.mu_pb_tilde{y}(:,t) + np.beta_pb(y)*(sl.c_pb_tilde{y}(:,t)-sl.c_pb{y}(:,t) ) ;
                
                % error
                res_dso2(t) = norm([res_dso2(t);sl.c_pb{y}(:,t)],inf);
                
                
                res1(np.n+y,t+1)= norm(s.u_DSO{y}(:,t+1)-s.u_DSO{y}(:,t),inf);
            end

            % Dual variable (trading with the main grid)
            s.sigma_tg(:,t+1) = zeros(np.h,1) ;
            for yy = 1:length(np.B_mg)
                y = np.B_mg(yy);
                s.sigma_tg(:,t+1) = s.sigma_tg(:,t+1) + s.p_tg{y}(:,t+1);                
            end
            sl.c_tg(:, t+1) = s.sigma_mg(:,t+1) - s.sigma_tg(:,t+1);
            
            s.mu_tg(:,t+1) = s.mu_tg_tilde(:,t) + np.beta_tg*(sl.c_tg_tilde(:,t)-sl.c_tg(:,t));
            
            res_gridTrading = norm(sl.c_tg(:,t),inf);
            
            % Extra: stopping criterion
            %s.comp_time(t) = toc;
            s.res(t) = (res2(t));           
            s.res1(t) = norm(res1(:,t+1),inf);
            
            error_v(:,t) = [s.res(t);res_dso2(t);s.res1(t);res_gridTrading];
            
            error_v(:,t)
            error(t) = norm(error_v(:,t),inf);
            
            temp.error_v = error_v;
            temp.error = error;
            %temp.time = s.comp_time;
            if error(t) < np.er_max
                break
            end
            

            t= t+1;
            %save(['p2p_temp',date],'temp')
            %toc
        end
        
        s.error = error;
        s.error_v = error_v;
end