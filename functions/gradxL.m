function nablaL = gradxL(np,pVar,dVar,i)
    
    
    nablaFi = np.d_mg*(pVar.sigma_mg-pVar.p_mg);
    
    c1 =[];
    
    for h = 1:np.h
        c_1 =  [0; 0; 0; nablaFi(h,1); zeros(2*sum(np.Adj(i,:)),1)]; % Nash eq.
        %c_1 = [0; 0; 0; np.d_mg*(s.sigma_mg(h,t)); zeros(2*sum(np.Adj(i,:)),1)];            % Wadrop eq.
        c1 = [c1;c_1];
    end
    
    % linear term associated with p_di
    
    c2 = np.Sdi{i}'*-dVar.mu_pb;
    
    % linear term associated with p_ch
    c3a = np.Sch{i}'*dVar.mu_pb;
    
    % linear term associated with p_ds
    c3b = np.Sds{i}'*-dVar.mu_pb;
    
    c3 = c3a+c3b;
    
    % linear term associated with grid constraints
    c4 = dVar.lambda_mg'*[np.Smg{i};-np.Smg{i}]+dVar.mu_tg'*np.Smg{i};
    c4 = c4';
    
    % linear term associated with reciprocity constraints
    c5 = zeros(np.h*(np.no_localDecision+2*sum(np.Adj(i,:))),1);
    for j=1:np.n
        if np.Adj(i,j) == 1
            c5 = c5 + np.Str{i,j}'*dVar.mu_tr{i,j};
        end     
    end
    %c3 = c3;
    
    % linear term associated with proximal term
    c6 = -np.A{i}*pVar.u;
    
    nablaL = c1+ c2 + c3 + c4 + c5 + c6;


end