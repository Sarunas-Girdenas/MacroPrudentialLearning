classdef SteadyStateSolution
   
    % steady state solution for learning paper
    
    properties(GetAccess = 'public', SetAccess = 'public')
        
        epsilon;
        v;
        gamma;
        beta;
        m_star_st;
        j;
        eta;
        theta;
        m;
        r_pi;
        chi; %problematic
       
        
    end
    
    methods
       
        function f = solveSteadyState(obj)
            
            % here we use ratio_hh to compute all the steady state values
            
            ratio_hh = fzero(@find_h,0);
            h_prim_st = 1/(1+ratio_hh);
            h_st      = 1-h_prim_st;
            L_st      = ((1+(1-obj.v)*(obj.epsilon-1))/((1-obj.v)*(obj.epsilon-1)-obj.m_star_st*obj.j*ratio_hh*(1-obj.v)*(obj.epsilon-1)))^(-1/obj.eta);
            Y_st      = h_st^obj.v*L_st^(1-obj.v);
            R_st      = 1/obj.beta;
            c_prim_st = (1-obj.v)*(obj.epsilon-1)*Y_st/obj.epsilon/L_st^(-obj.eta);
            c_st      = Y_st-c_prim_st;
            w_st      = (1-obj.v)*(obj.epsilon-1)/obj.epsilon*Y_st/L_st;
            b_st      = obj.m_star_st*obj.beta*obj.j*c_prim_st*h_st/(1-obj.beta)/h_prim_st;
            b_prim_st = -b_st;
            q_st      = obj.j/(1-obj.beta)*c_prim_st/h_prim_st;
            lambda_st = 1/c_st-obj.gamma*R_st/c_st;
            X_st      = obj.epsilon/(obj.epsilon-1);
            Md_st     = obj.chi*c_prim_st/(1-obj.beta);
            Ms_st     = Md_st;
            
            f = struct('ratio_hh',ratio_hh,'h_prim_st',h_prim_st,'h_st',h_st,'L_st',L_st,'Y_st',Y_st,'R_st',R_st,...
                'c_prim_st',c_prim_st,'c_st',c_st,'w_st',w_st,'b_st',b_st,'b_prim_st',b_prim_st,'q_st',q_st,...
                'lambda_st',lambda_st,'X_st',X_st,'Md_st',Md_st,'Ms_st', Ms_st);
            
        end
            
        function f = find_h(obj)
            
            % function that gives ratio_hh when evaluated at zero
                 
            f = (1-obj.gamma-obj.m_star_st*obj.beta*(1-obj.gamma/obj.beta))*obj.j...
            *(1+(1-obj.v)*(obj.epsilon-1))*x/(1-obj.beta)/((obj.epsilon-1)-(obj.epsilon-1)*obj.m_star_st*obj.j*obj.x)-obj.gamma*obj.v;

        end         
        
    end
    
   
    
    
end