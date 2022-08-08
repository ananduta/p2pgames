function [s,sl,np] = sd_alg_eP2P_l1_DSO2(np)
%% Distributed algorithm
% extended P2P market
% W. Ananduta
% 12/11/2020

        %% INITIALIZATION
        % Assigning parameters of the algorithm
        %np = alg_param(np);
        %np = alg_param_37b_l2(np); 
        np = alg_param_37b_f(np);
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
            tic
            t
            % 1) Strategy update
            s.sigma_mg(:,t+1) = np.sumPd;
            
            
           for i=1:np.n
                
                % primal update of prosumer (quadratic prog.)
                tic
                s = loc_opt_qprog_l1(np,s,t,i);
                s.comp_p(i) = toc;
                % local load imbalance of prosumer i
                sl.b{i}(:,t+1) = np.Pd(i,1:np.h)' - s.p_di{i}(:,t+1) - s.p_ds{i}(:,t+1) + s.p_ch{i}(:,t+1);

                % forward p_mg to DSO
                s.sigma_mg(:,t+1) = s.sigma_mg(:,t+1) + s.p_mg{i}(:,t+1);
                
                % compute error

                res1(i,t+1)= norm(s.u{i}(:,t+1)-s.u{i}(:,t),inf);
            end

            res2(t) = 0; % Extra: compute norm of residual
            dres2(t) = 0;

            % dual variables (reciprocity constraint) Prosumers
            for i=1:np.n

                for j=1:np.n
                    if np.Adj(i,j) == 1

                        % Dual variables (reciprocity constraints) update:
                        
                        % Aux vector
                        sl.c_tr{i,j}(:,t+1) = s.p_tr{i,j}(:,t+1) + s.p_tr{j,i}(:,t+1);
                        
                        % Reflected dual ascent
                        s.mu_tr{i,j}(:,t+1) = s.mu_tr{i,j}(:,t) + np.beta_tr(i,j)*(2*sl.c_tr{i,j}(:,t+1) - sl.c_tr{i,j}(:,t));

                        % Extra: compute norm of residual
                        res = sl.c_tr{i,j}(:,t+1);
                        res2(t) = norm([res2(t);res],inf);

                    end
                end
            end


            % DSO=========================================================================================================================================
            % Primal update 
            
            for y=1:np.b
                % create concatenations of p^l_(y,z) (k) and q^l_(y,z) (k)
                pl = s.u_DSO{y}(np.h*3+1:np.h*(3+length(np.B{y})),t);
                ql = s.u_DSO{y}(np.h*(3+length(np.B{y}))+1:end,t);
                % add the concatenations into a_DSO{y}
                a_DSO{y} = [zeros(2*np.h,1); -s.mu_tg(:,t)-s.mu_pb{y}(:,t);kron(ones(length(np.B{y}),1),-s.mu_pb{y}(:,t))+pl;zeros(np.h*length(np.B{y}),1)+ql]; 
            end
            tic
            s = projDSO2(s,a_DSO,np,t);                                         
            s.comp_DNO = toc;                       
            % Dual variable (grid constraints)
            sl.b_DSO(:,t+1) = kron([1;-1],2*s.sigma_mg(:,t+1)-s.sigma_mg(:,t))-[np.pmg_max*ones(np.h,1);-np.pmg_min*ones(np.h,1)];
            s.lambda_mg(:,t+1) = max(0, s.lambda_mg(:,t) + np.gamma_mg*sl.b_DSO(:,t+1));        

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
                s.mu_pb{y}(:,t+1) = s.mu_pb{y}(:,t) + np.beta_pb(y)*(2*sl.c_pb{y}(:,t+1) - sl.c_pb{y}(:,t)) ;
                
                % error
                res_dso2(t) = norm([res_dso2(t);sl.c_pb{y}(:,t+1)],inf);
                
                
                res1(np.n+y,t+1)= norm(s.u_DSO{y}(:,t+1)-s.u_DSO{y}(:,t),inf);
            end

            % Dual variable (trading with the main grid)
            s.sigma_tg(:,t+1) = zeros(np.h,1) ;
            for yy = 1:length(np.B_mg)
                y = np.B_mg(yy);
                s.sigma_tg(:,t+1) = s.sigma_tg(:,t+1) + s.p_tg{y}(:,t+1);                
            end
            sl.c_tg(:, t+1) = s.sigma_mg(:,t+1) - s.sigma_tg(:,t+1);
            s.mu_tg(:,t+1) = s.mu_tg(:,t) + np.beta_tg*(2*sl.c_tg(:,t+1)-sl.c_tg(:,t));
            
            
            s.comp_time(t) = toc;
            
            res_gridTrading = norm(sl.c_tg(:,t+1),inf);
            
            
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