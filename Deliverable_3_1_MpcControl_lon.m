classdef MpcControl_lon < MpcControlBase
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   x0           - initial state (estimate)
            %   V_ref, u_ref - reference state/input
            %   d_est        - disturbance estimate
            %   x0other      - initial state of other car
            % OUTPUTS
            %   u0           - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(mpc.H/mpc.Ts); % Horizon steps
            N = N_segs + 1;              % Last index in 1-based Matlab indexing

            [nx, nu] = size(mpc.B);
            % Targets
            V_ref = sdpvar(1);
            u_ref = sdpvar(1);

            % Disturbance estimate (Ignore this before Todo 4.1)
            d_est = sdpvar(1);

            % Initial states
            x0 = sdpvar(nx, 1);
            x0other = sdpvar(nx, 1); % (Ignore this before Todo 5.1)

            % Input to apply to the system
            u0 = sdpvar(nu, 1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D
            %       are the DISCRETE-TIME MODEL of your system.
            %       You can find the linearization steady-state
            %       in mpc.xs and mpc.us.
            
            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            X = sdpvar(1,N); % state trajectory
            U = sdpvar(1,N); % input trajectory 

            A = mpc.A(2,2);
            B = mpc.B(2,1);
            xs = mpc.xs(2);
            us = mpc.us;
            f_xs_us = mpc.f_xs_us(2);
            Q = 50;
            R = 5;
            % Constraints 
            % input constraints: 
            M = [1;-1];
            m = [1;1];
            % State constraints:
            F = -1;
            f = 0;
            % LQR unconstrained: 
            [K, Qf, ~] = dlqr(A,B,Q,R);
            K = -K; 
            % Maximal invariant set:
            P1 = [F; M*K];
            p2 = [f;m+M*(K*xs-us)];

            

            Xf = polytope(P1,p2);
            A_cl = A+B*K;
            double(Xf)
            % Compute the maximal invariant set
            while 1
                Xf_prev = Xf;
                [T,t] = double(Xf);
                preS_Xf = polytope(T*A_cl,t);
                Xf = intersect(Xf,preS_Xf);
                if isequal(Xf,Xf_prev)
                    break
                end
            end
            [Hf,hf] = double(Xf); % get the matrix information from the terminal set
            obj = 0;
            con = [];
            % figure;
            % plot(polytope(Hf, hf),'g');
            % title('Terminal Invariant Set for Longitudinal Subsystem');
            % xlabel('Velocity (m/s)');
            % ylabel('Throttle Input');
            

            % Initial condition
            con = con + (X(:, 1) == x0(2));


            for i = 1:N-1
                % system dynamics
                con = con + (X(:,i+1) == A*(X(:,i)-xs) + B*(U(:,i)-us) + f_xs_us);
                % input constraints
                con = con + (M*U(:,i) <= m); 
                % state constraints
                con = con + (F*X(:,i) <= f);
                % compute cost
                obj = obj + ((X(:,i)-V_ref)'*Q*(X(:,i)-V_ref) + (U(:,i)-u_ref)'*R*(U(:,i)-u_ref));  
            end
            % Terminal Set Constraints
            con = con + (Hf*X(:,N) <= hf);
            % Terminal cost
            obj = obj +((X(:,N)-V_ref)'*Qf*(X(:,N)-V_ref)); 

            % Replace this line and set u0 to be the input that you
            % want applied to the system. Note that u0 is applied directly
            % to the nonlinear system. You need to take care of any 
            % offsets resulting from the linearization.
            % If you want to use the delta formulation make sure to
            % substract mpc.xs/mpc.us accordingly.
            con = con + ( u0 == U(:,1) );

            % Pass here YALMIP sdpvars which you want to debug. You can
            % then access them when calling your mpc controller like
            % \\[u, X, U] = mpc_lon.get_u(x0, ref);
            % with debugVars = {X_var, U_var};
            debugVars = {X,U};
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {x0, V_ref, u_ref, d_est, x0other}, {u0, debugVars{:}});
        end
        
        % Computes the steady state target which is passed to the
        % controller
        function [Vs_ref, us_ref] = compute_steady_state_target(mpc, ref, d_est)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            %   d_est  - disturbance estimate (Ignore before Todo 4.1)
            % OUTPUTS
            %   Vs_ref, us_ref - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Steady-state subsystem
            A = mpc.A(2, 2);
            B = mpc.B(2, 1);

            % Subsystem linearization steady-state
            xs = mpc.xs(2);
            us = mpc.us;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            Vs_ref = ref;
            us_ref = ((1-A)/B)*(ref-xs)+us;
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
