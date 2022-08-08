function [s,o] = com_cost(np,o)
    np = build_mat_exP2P_tr1(np);
    
    for i = 1:np.n
        Qh = np.Qh{i};
        ch = np.ch{i};
        
        
        % Compute potential function
        negQ = diag([0 0 0 np.d_mg zeros(1,length(np.N{i})) zeros(1,length(np.N{i}))]);
        negQ = kron(eye(np.h),negQ);
        Qloc = Qh - negQ;
        floc = o.u{i}'*Qloc*o.u{i} + ch'*o.u{i};
        
        D = np.Smg{i}'*np.d_mg*np.Smg{i};
        
        o.P(i) = 0.5*(o.J(i) + floc + o.u{i}'*D*o.u{i});
        
    end
    
    
    % Compute cost
        
    
    o.Pt = sum(o.P);
    o.iter = s.t;
end