addpath D:\Charbel\EPFL\'ME 425 - Model Predictive Control'\Matlab\casadi-3.6.6-windows64-matlab2018b
%%
Ts=0.1;
car = Car(Ts);

H = 5;
mpc = NmpcControl(car, H);
x0 = [0 0 0 80/3.6]';
ref = [3 100/3.6]';
u = mpc.get_u(x0, ref);

params = {};
params.Tf = 15;
params.myCar.model = car;
params.myCar.x0 = [0 0 0 80/3.6]';
params.myCar.u = @mpc.get_u;
ref1 = [0 80/3.6]';
ref2 = [3 100/3.6]';
params.myCar.ref = car.ref_step(ref1, ref2, 2);
result = simulate(params);
visualization(car, result);