Desired_Wp = [5;5;5];

model = multirotor;
s = state(model);

s(end) = 0;

Time_Of_Sim = [0 20];

[x, y, z, vx, vy, vz, yaw, pitch, roll, p, q, r, t] = ode45(@(state) dState(state), Time_Of_Sim, s)
dState(s)

