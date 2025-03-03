classdef NmpcControl < handle

    properties
        % The NMPC problem
        opti

        % Problem parameters
        x0, ref, x0other

        % Most recent problem solution
        sol

        % The input that you want to apply to the system
        u0

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add any variables you would like to read to debug here
        % and then store them in the NmpcControl function below.
        % e.g., you could place X here and then add obj.X = X
        % in the NmpcControl function below.
        % 
        % After solving the problem, you can then read these variables 
        % to debug via
        %   nmpc.sol.value(nmpc.X)
        % 
        X, U
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end

    methods
        function obj = NmpcControl(car, H)

            import casadi.*

            N_segs = ceil(H/car.Ts); % Horizon steps
            N = N_segs + 1;          % Last index in 1-based Matlab indexing

            nx = 4;
            nu = 2;

            % Define the NMPC optimization problem
            opti = casadi.Opti();
            
            % Parameters (symbolic)
            obj.x0 = opti.parameter(nx, 1);       % initial state
            obj.ref = opti.parameter(2, 1);       % target y, velocity
            obj.x0other = opti.parameter(nx, 1);  % initial state of other car

            % SET THIS VALUE TO BE YOUR CONTROL INPUT
            obj.u0 = opti.variable(nu, 1);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE

            % Define your problem using the opti object created above
            f= @(x,u) car.f(x,u);
            Ts=car.Ts;
            obj.X=opti.variable(nx,N);
            obj.U=opti.variable(nu,N_segs);
            y_ref = obj.ref(1);
            v_ref = obj.ref(2);
            Ref=[0;y_ref;0;v_ref];

            %Qt = 100; % Terminal weight
            
            Q=diag([0,30,0,30]);
            R=diag([5;5]);

            % Initialize cost
            cost = 0;
            
            % Loop over the prediction horizon
            for k = 1:N_segs
                cost= cost + (obj.X(:,k) - Ref)'*Q*(obj.X(:,k) - Ref) + (obj.U(:,k))'*R*obj.U(:, k);
               
                %RK4 constraints
                opti.subject_to(obj.X(:,k+1) == RK4(obj.X(:,k), obj.U(:,k),Ts,f));

            end
            
            %cost = cost + Qt * (obj.X(2, N) - y_ref)^2 + Qt * (obj.X(4, N) - v_ref)^2;
            
            % change this line accordingly
            opti.subject_to( obj.u0 == obj.U(:,1));



            %Input constraints
            opti.subject_to( -0.5236 <= obj.U(1,:) <= 0.5236);
            opti.subject_to( -1 <= obj.U(2,:) <= 1);

            %State constraints
            opti.subject_to( -0.5 <= obj.X(2,:) <= 3.5);
            opti.subject_to( -0.0873 <= obj.X(3,:) <= 0.0873);

            opti.subject_to(obj.X(:,1)==obj.x0);


            opti.minimize(cost);

            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Store the defined problem to solve in get_u
            obj.opti = opti;

            % Setup solver
            options = struct;
            options.ipopt.print_level = 0;
            options.print_time = 0;
            options.expand = true;
            obj.opti.solver('ipopt', options);
        end

        function u = get_u(obj, x0, ref, x0other)

            if nargin < 4
                x0other = zeros(4, 1);
            end

            % Compute solution from x0
            obj.solve(x0(1:4), ref, x0other(1:4));
            % disp(obj.sol.value(obj.U(:,1)));
            u = obj.sol.value(obj.u0);
        end

        function solve(obj, x0, ref, x0other)

            % Pass parameter values
            obj.opti.set_value(obj.x0, x0);
            obj.opti.set_value(obj.ref, ref);
            obj.opti.set_value(obj.x0other, x0other);

            obj.sol = obj.opti.solve();   % actual solve
            % Set warm start for next solve
            obj.opti.set_initial(obj.sol.value_variables());
            obj.opti.set_initial(obj.opti.lam_g, obj.sol.value(obj.opti.lam_g));
        end
    end
end


function [x_next] = RK4(X,U,h,f)
%
% Inputs : 
%    X, U current state and input
%    h    sample period
%    f    continuous time dynamics f(x,u)
% Returns
%    State h seconds in the future
%

% Runge-Kutta 4 integration
% write your function here
   k1 = f(X,         U);
   k2 = f(X+h/2*k1, U);
   k3 = f(X+h/2*k2, U);
   k4 = f(X+h*k3,   U);
   x_next = X + h/6*(k1+2*k2+2*k3+k4);
end


