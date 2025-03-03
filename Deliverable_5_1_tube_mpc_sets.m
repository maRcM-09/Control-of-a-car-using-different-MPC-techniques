Ts = 1/10; % Sample time
car = Car(Ts);
[xs, us] = car.steady_state(120 / 3.6);
sys = car.linearize(xs, us);
[sys_lon, sys_lat] = car.decompose(sys);
[~, A, B, ~, ~] = Car.c2d_with_offset(sys_lon, Ts);
F_w = [1;-1];
f_w = [0.5;0.5];
W = B*Polyhedron(F_w,f_w);
Q = diag([50,50]);
R = 5;
%% Design stabilizing controller K such that ||A-BK|| < 1
[K,Qf,P] = dlqr(A,-B,Q,R);
K = -K;
if norm(max(abs(eig(A-B*K)))) < 1
    disp("Good controller")
    disp(max(abs(eig(A-B*K))))
else
    disp("Not stabilizing")
    disp(norm(A-B*K))
end

%% compute minimal robust invariant set 
Eps = W;% Omega = {0}
A_cl = A-B*K;
k = 1;
while true
    % Compute A^i * W
    Eps = Eps + (A_cl^k)*W;
    Eps = minHRep(Eps);
    % Check for convergence: Omega_(i+1) == Omega_i
    if norm(A_cl^k) < 1e-2
        break;
    end
    k = k+1;
end
figure(1)
plot(Eps);

%% Caculate the tightenen constraints
% tightening state constraints
F_x = [-1 0];
x_safe = 10;
f_x = -6+x_safe; 
X = Polyhedron(F_x, f_x);
X_m_E = X - Eps;
% tightening input set
F_u = [1; -1];
f_u = [1; 1];
U = Polyhedron(F_u, f_u);
U_m_KE = U - K * Eps;
% save these variables in the following 
F_x_t = X_m_E.A;
f_x_t = X_m_E.b;
M_u_t = U_m_KE.A;
m_u_t = U_m_KE.b;

%% Terminal Set computation
Xf = Polyhedron([F_x_t; M_u_t*K],[f_x_t; m_u_t]);% we want to regulate relative dynamics here
figure(2)
while true
    % save set from previous iteration
    Xf_prev = Xf;
    Xf_preS = Polyhedron(Xf.A*A_cl,Xf.b);
    Xf = intersect(Xf, Xf_preS);
    if abs(Xf.volume - Xf_prev.volume) < 1e-10
        break
    end
end
plot(Xf);
hold off;
%% Save the data in current directory 
Xf_A = Xf.A;
Xf_b = Xf.b;
Eps_A = Eps.A;
Eps_b = Eps.b;
save('tube_mpc_data.mat', 'F_x_t', 'M_u_t', 'f_x_t', 'm_u_t', 'Xf_A', 'Xf_b', 'Qf', 'Eps_A', 'Eps_b', 'K', 'Q', 'R');
