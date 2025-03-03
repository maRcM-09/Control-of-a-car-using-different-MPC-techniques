addpath D:\Charbel\EPFL\'ME 425 - Model Predictive Control'\Matlab\casadi-3.6.6-windows64-matlab2018b
%%
Ts=0.1;
car = Car(Ts);
H = 8;
mpc = NmpcControlOvertake(car,H);
x0_ego = [0 0 0 80/3.6]';
x0_other = [20 0 0 80/3.6]';
ref1 = [0 80/3.6]';
ref2 = [0 100/3.6]';
params = {};
params.Tf = 15;
params.myCar.model = car;
params.myCar.x0 = x0_ego;
params.myCar.u = @mpc.get_u;
params.myCar.ref = car.ref_step(ref1, ref2, 1);
params.otherCar.model = car;
params.otherCar.x0 = x0_other;
params.otherCar.u = car.u_const(80/3.6);
result = simulate(params);
visualization(car, result);