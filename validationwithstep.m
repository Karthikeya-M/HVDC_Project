close all;clc  % Close open windows, clear the workspace and the command window
set(0,'defaulttextinterpreter','latex') % Font settings for good-looking plots
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',14)
format compact; format long % Number format to be shown in the command window
format short G

%% System parameters
wn = 2*pi*50; % Base angular frequency [rad/s]
% Coupling matrix resulting from the abc to dq frame transformation (T) of diff. eqs.
Wpu = [0 1;-1 0]; % = T*d(T^-1)/dt, T is the abc-to-dq matrix rotation

% Small decoupling capacitance
c_t = 0.05; % [pu]

% Converter AC-side parameters
l_f = 0.08; % Filter inductor [pu]
r_f = 0.01*l_f; % 1 percent losses

l_t = 0.15; % Transformer leackage impedance [pu]
r_t = l_t/20; % Transformer X/R ratio is usually larger than 10

% Thevenin equivalent
% Inductance CANNOT be set to zero to obtain the VSC admittance
% Instead, a very high SCR value can be used
XR = 10; % A typical X/R ratio in transmission systems of 10 is used
SCR = 5; % Short-circuit ratio at the HV point of connection
l_g = 1/(SCR*sqrt(1+1/XR^2)); % Grid reactance [pu]
r_g = l_g/XR; % Grid resistance [pu]

% Current controller (tuning via IMC)
w_CL = 1e3; % [rad/s]
kp_CL = l_f*w_CL/wn; % Internal model control (IMC) for a 1st order response
ki_CL = r_f*w_CL;
wf_v = 5*w_CL;  % Feedforward voltages filter cut-off [rad/s]
wf_i = 5*w_CL;  % Measurement currents filter cut-off [rad/s]

% PLL parameters (tuning via IMC)
z_pll = 1.25; % PLL damping ratio
w_pll = 0.05*w_CL; % PLL designed BW [rad/s]
ki_pll = w_pll^2/wn; 
kp_pll = 2*z_pll*w_pll/wn; % Proportional and integral gains
wf_PLL = 2*pi*50; % Voltage PF [rad/s]

% Outer controllers (tuning via IMC)
w_P = w_CL/25; % Power control bandwidth [rad/s]
kp_P = w_P/w_CL; % Pac control parameters
ki_P = w_P;
kp_Q = w_P/w_CL; % Q control parameters
ki_Q = w_P;

% % Grid-supporting functions
% w_f = w_P/10; % Frequency support BW [rad/s]
% kp_f = 5; % Proportional gain: dp/dw i.e. inverse of droop [pu/pu]
% w_v = w_P/10; % Vac support BW [rad/s]
% kp_v = 5; % Proportional gain: dq/dv i.e. inverse of droop [pu/pu]

%% Declaration of symbolic variables
% Electrical variables
w = sym('w','real'); % Infinite bus (IB) frequency [pu]

vd_th=sym('vd_th','real'); % Voltages at the IB
vq_th=sym('vq_th','real');
v_th = [vd_th;vq_th]; % Merge into a vector

id_pcc=sym('id_pcc','real'); % Currents at the PCC
iq_pcc=sym('iq_pcc','real');
i_pcc = [id_pcc;iq_pcc]; % Merge into a vector

id_t=sym('id_t','real'); % L-filter current
iq_t=sym('iq_t','real');
i_t = [id_t;iq_t]; % Vector form

vd_t=sym('vd_t','real'); % Capacitor voltage
vq_t=sym('vq_t','real');
v_t = [vd_t;vq_t]; % Vector form

vd_m=sym('vd_m','real'); % VSC modullated voltage
vq_m=sym('vq_m','real');
v_m = [vd_m;vq_m]; % Merge into a vector

% Control variables
id_f=sym('id_f','real'); % Filtered currents
iq_f=sym('iq_f','real');
i_f = [id_f;iq_f];

vd_f=sym('vd_f','real'); % L-filter output voltages filtered
vq_f=sym('vq_f','real');
v_f = [vd_f;vq_f];

% Voltage and frequency support
dqref_f=sym('dqref_f','real'); % Filtered delta reactive power reference [pu]
dpref_f=sym('dpref_f','real'); % Filtered delta active power reference [pu]

e_P=sym('e_P','real'); % Epsilon power control
e_Q=sym('e_Q','real'); % Epsilon = integrator states

e_pll = sym('e_pll','real'); % Epsilon PLL
theta_vsc = sym('theta_vsc','real'); % Angle between the grid and the VSC frames
vq_pll=sym('vq_pll','real'); % Filtered q-axis voltage into the PLL

ed_CL=sym('ed_CL','real'); % Epsilon current control (d-axis)
eq_CL=sym('eq_CL','real'); % Epsilon current control (q-axis)
e_CL = [ed_CL;eq_CL]; % Merge into a vector

% Reference values: grid-support and PQ-controls
v_ref=sym('v_ref','real'); % Reference voltage magnitude [pu]
w_ref=sym('w_ref','real'); % Frequency setpoint [pu]
p_ref=sym('p_ref','real'); % Active power reference [pu]
q_ref=sym('q_ref','real'); % Reactive power reference [pu]

%% System equations
% Series RL dynamics: filter + transformer
di_t_dt = wn/(l_f+l_t)*(v_m - v_t - (r_f+r_t)*i_t - (l_f+l_t)*Wpu*w*i_t);
% Shunt capacitor
dv_t_dt = wn/c_t*(i_t - i_pcc - c_t*w*Wpu*v_t); % Shunt capacitance
% Series RL: grid impedance
di_pcc_dt = wn/l_g*(v_t - v_th - r_g*i_pcc - l_g*Wpu*w*i_pcc); 

% 2L-VSC control
% Frame rotation from grid-frame to the VSC-frame
T = [cos(theta_vsc), -sin(theta_vsc);
    sin(theta_vsc),  cos(theta_vsc)];

% Control point, Point of Common Coupling (PCC), is at the transformer LV-side
ic_t = T*i_t; % Current in the converter frame
vc_t = simplify(T * (v_t + r_t*i_t + l_t/wn*di_t_dt + l_t*w*Wpu*i_t) ); % Filter output voltage (transformer LV)

% LPF v&i
% wf_v = 2*pi/2e-4;
% wf_i = 2*pi/2e-4;
dv_f_dt = wf_v*(vc_t - v_f);
di_f_dt = wf_i*(ic_t - i_f);

% PLL: drives the q-axis voltage to zero
dvq_pll_dt = wf_PLL*(v_f(2) - vq_pll); % LPF the q-axis voltage
de_pll_dt = ki_pll*(0 - vq_pll); % PLL integrator
w_pll = kp_pll*(0 - vq_pll) + e_pll; % PLL frequency
dtheta_vsc_dt = wn*(w_pll - w); % Angular difference between frames

% Voltage support: creates a delta in the reference reactive power
% proportional to the voltage magnitude difference between the reference
% and the measurement. Idem with frequency support and active power.
% ddqref_f_dt = w_v*(kp_v*(v_ref - sqrt(v_f.'*v_f)) - dqref_f); % Voltage support
% ddpref_f_dt = w_f*(kp_f*(w_ref - w_pll) - dpref_f); % Frequency support

% Power control
p_m = v_f.'*i_f; % Non-linear power computation in the dq-frame
q_m = v_f.'*Wpu*i_f; % In per unit there is no scaling factor, i.e. 3/2
de_P_dt = ki_P*(p_ref - p_m); % Integral control
de_Q_dt = ki_Q*(q_ref - q_m);

i_ref = [kp_P*(p_ref  - p_m) + e_P;
    kp_Q*(q_ref - q_m) + e_Q;]; % Current references

% Current control
de_CL_dt = ki_CL*(i_ref - i_f); % Integral action
vc_m = v_f + kp_CL*(i_ref - i_f) + e_CL + w_pll*l_f*Wpu*i_f; % PI + decoupling and FF

% Modullation and inverse rotation to the grid reference frame
v_m = simplify(T\vc_m); % DC voltage is considered constant

%% Steady-state operating point
% Setpoints and inputs [pu]
 v_th0 = [0.982363188981222;  0]; % Thevenin voltage if using a non-zero grid impedance, otherwise HV-trafo voltage
%v_th0 = [1.0;  0];
p_ref0 = 1.0; % Active power setpoint
q_ref0 = 0; % Reactive power setpoint
%v_ref0 = 1.00; % Voltage magnitude setpoint
%w_ref0 = 1.00; % Frequency setpoint
w0 = 1.0; % Steady-state grid frequency

op = [v_th;p_ref;q_ref;w];      % Initial op point variables
op0 = [v_th0;p_ref0;q_ref0;w0]; % Initial operating point values

% Differential equations: dx/dt
f = [di_pcc_dt;di_t_dt;dv_t_dt;dvq_pll_dt;de_pll_dt;...
    dtheta_vsc_dt;de_CL_dt;di_f_dt;dv_f_dt;de_P_dt;de_Q_dt];
f = subs(f,[vd_m;vq_m],v_m); % Substitue VSC output variables to close the loop
f = subs(f,op,op0); % Substitue the setpoints and inputs

x_sol = symvar(f); % Variables to be found
x_sol_str = string(x_sol);

% Solution range for the unknowns as [min_value, max_value]
x00 = 1.3*[-ones(length(x_sol),1),ones(length(x_sol),1)]; % Search range for solution is around +-1 pu
x00(x_sol_str == "theta_vsc",1:2) = [-pi,pi]/2; % Wrap angle between -90 and 90 (static stability)
x00(x_sol_str == "vd_f",1:2) = [0.7,1.3]; % VSC d-axis voltage aligned

% Solve the non-linear system of equations
x00 = vpasolve(f,x_sol,x00);
x00 = struct2cell(x00); % Write the solution into a vector format
x00 = double([x00{:}])';

% Retrieve the solution
e_pll0 = x00(x_sol_str == "e_pll");
theta_vsc0 = x00(x_sol_str == "theta_vsc");
vq_pll0 = x00(x_sol_str == "vq_pll");
e_CL0 = [x00(x_sol_str == "ed_CL");x00(x_sol_str == "eq_CL")];
i_pcc0 = [x00(x_sol_str == "id_pcc");x00(x_sol_str == "iq_pcc")];
e_P0 = x00(x_sol_str == "e_P");
e_Q0 = x00(x_sol_str == "e_Q");
v_f0 = [x00(x_sol_str == "vd_f");x00(x_sol_str == "vq_f")];
i_f0 = [x00(x_sol_str == "id_f");x00(x_sol_str == "iq_f")];

% dqref_f0 = x00(x_sol_str == "dqref_f");
% dpref_f0 = x00(x_sol_str == "dpref_f");
v_mag_LV0 = sqrt(v_f0.'*v_f0); % Voltage magnitude at the LV-side
v_HV = v_th + r_g*i_pcc + l_g*w*Wpu*i_pcc; % Steady-state voltage magnitude at the HV-side
% Get the numeric value (use as input v_th0 in case of high grid impedance)
v_mag_HV0 = double(subs(sqrt(v_HV.'*v_HV),[x_sol.';op],[x00;op0])); 

% State-space form
%Define the derivatives
dx_dt = subs([di_t_dt;dvq_pll_dt;de_pll_dt;...
    dtheta_vsc_dt;de_CL_dt;di_f_dt;dv_f_dt;de_P_dt;de_Q_dt;dv_t_dt;di_pcc_dt],[vd_m;vq_m],v_m); % Substitue VSC output variables to close the loop
x = [i_t;vq_pll;e_pll;theta_vsc;e_CL;i_f;v_f;e_P;e_Q;v_t;i_pcc]; % State variables
x0 = double(subs(x,[x_sol.';op],[x00;op0]));

% Inputs
u = [p_ref;q_ref]; % Both inputs and outputs are in the global frame so no need for additional frame rotations
u0 = double(subs(u,[x_sol.';op],[x00;op0]));

% Outputs
y = [i_f;e_P;e_Q;vd_f;vq_f]; % Both inputs and outputs are in the global frame
y0 = double(subs(y,[x_sol.';op],[x00;op0]));

% syms_left = symvar(dx_dt);
% disp(setdiff(syms_left, [x;u;op]));

% Linearize around the operating point
A = double(subs(jacobian(dx_dt, x),[x;u;op],[x0;u0;op0]));
B = double(subs(jacobian(dx_dt, u),[x;u;op],[x0;u0;op0]));
C = double(subs(jacobian(y, x),[x;u;op],[x0;u0;op0]));
D = double(subs(jacobian(y, u),[x;u;op],[x0;u0;op0]));

%Build the state-space model
model = ss(A,B,C,D,'StateName',string(x),'InputName',string(u),'OutputName',string(y));

%% --- Time-domain simulation setup (with step change) ---
Tfinal = 5;               % total time [s]
dt = 1e-5;                % time step [s]
t = 0:dt:Tfinal;          % time vector

% Initialize input matrix [p_ref; q_ref]
u = zeros(2, length(t));  
u(1,:) = 1.0;   % initial p_ref = 1.0
u(2,:) = 0.0;   % initial q_ref = 0.0

% Apply step changes at 2 seconds
u(1, t >= 2) = 1.0 - 0.03;  % Step of +0.01 in p_ref
u(2, t >= 2) = 0.0 + 0.05;  % Step of +0.01 in q_ref

%% --- Simulate linearized system response ---
[y, t_out, x_out] = lsim(model, u', t,x0);

%% --- Quick visualization of step inputs ---
figure;
plot(t_out, u(1,:), 'b', 'LineWidth', 1.4); hold on;
plot(t_out, u(2,:), 'r', 'LineWidth', 1.4);
xlabel('Time [s]');
ylabel('Input [pu]');
title('Step Inputs: $p_{ref}$ and $q_{ref}$','Interpreter','latex');
legend('$p_{ref}$','$q_{ref}$','Interpreter','latex','Location','best');
grid on;

%% --- Load PSCAD data ---
data = readmatrix('var.txt');  % Columns: Domain, e_P, e_Q, id, iq, vd, vq, ed_CL, eq_CL

t_txt  = data(:,1);  % Time [s]
eP_txt = data(:,2);  % e_P
eQ_txt = data(:,3);  % e_Q
id_txt = data(:,4);  % i_d,pcc
iq_txt = data(:,5);  % i_q,pcc
vd_txt = data(:,6);  % v_d,f
vq_txt = data(:,7);  % v_q,f
% edCL_txt = data(:,8);  % e_d,CL
% eqCL_txt = data(:,9);  % e_q,CL


%% --- Extract MATLAB simulation outputs ---
% Note: adjust indices if 'OutputName' order differs
% (use get(model,'OutputName') to verify)
id_m = y(:,1);
iq_m = y(:,2);
eP_m = y(:,3);
eQ_m = y(:,4);
vd_m = y(:,5);
vq_m = y(:,6);
% edCL_m = y(:,7);
% eqCL_m = y(:,8);


%% --- Figure 1: Power control integrators (e_P, e_Q) ---
figure;
subplot(2,1,1);
plot(t_out, eP_m, 'b', 'LineWidth', 1.3); hold on;
plot(t_txt, eP_txt, 'r--', 'LineWidth', 1.3);
xlabel('Time [s]'); ylabel('$e_P$ [pu]','Interpreter','latex');
title('Active Power Controller Integral ($e_P$)','Interpreter','latex');
legend('MATLAB','PSCAD','Location','best'); grid on;

subplot(2,1,2);
plot(t_out, eQ_m, 'b', 'LineWidth', 1.3); hold on;
plot(t_txt, eQ_txt, 'r--', 'LineWidth', 1.3);
xlabel('Time [s]'); ylabel('$e_Q$ [pu]','Interpreter','latex');
title('Reactive Power Controller Integral ($e_Q$)','Interpreter','latex');
legend('MATLAB','PSCAD','Location','best'); grid on;
sgtitle('Comparison of Power Controller Integrators','Interpreter','latex');

%% --- Figure 2: PCC Currents (id, iq) ---
figure;
subplot(2,1,1);
plot(t_out, id_m, 'b', 'LineWidth', 1.3); hold on;
plot(t_txt, id_txt, 'r--', 'LineWidth', 1.3);
xlabel('Time [s]'); ylabel('$i_d$ [pu]','Interpreter','latex');
title('Filtered $d$-axis Current','Interpreter','latex');
legend('MATLAB','PSCAD','Location','best'); grid on;

subplot(2,1,2);
plot(t_out, iq_m, 'b', 'LineWidth', 1.3); hold on;
plot(t_txt, iq_txt, 'r--', 'LineWidth', 1.3);
xlabel('Time [s]'); ylabel('$i_q$ [pu]','Interpreter','latex');
title('Filtered $q$-axis Current','Interpreter','latex');
legend('MATLAB','PSCAD','Location','best'); grid on;
sgtitle('Filtered Currents Comparison','Interpreter','latex');

% %% --- Figure 3: Filtered Voltages (vd_f, vq_f) ---
% figure;
% subplot(2,1,1);
% plot(t_out, vd_m, 'b', 'LineWidth', 1.3); hold on;
% plot(t_txt, vd_txt, 'r--', 'LineWidth', 1.3);
% xlabel('Time [s]'); ylabel('$v_d$ [pu]','Interpreter','latex');
% title('Filtered $d$-axis Voltage','Interpreter','latex');
% legend('MATLAB','PSCAD','Location','best'); grid on;
% 
% subplot(2,1,2);
% plot(t_out, vq_m, 'b', 'LineWidth', 1.3); hold on;
% plot(t_txt, vq_txt, 'r--', 'LineWidth', 1.3);
% xlabel('Time [s]'); ylabel('$v_q$ [pu]','Interpreter','latex');
% title('Filtered $q$-axis Voltage','Interpreter','latex');
% legend('MATLAB','PSCAD','Location','best'); grid on;
% sgtitle('Filtered Voltages Comparison','Interpreter','latex');

% Extract steady-state value for vd_f from x00
vd_f0 = x00(x_sol_str == "vd_f"); % This is your equilibrium value, e.g., 0.99

% Add the steady-state value to the lsim output
vd_m_absolute = vd_m + vd_f0;

% Now plot the absolute value against PSCAD data
figure;
subplot(2,1,1);
plot(t_out, vd_m_absolute, 'b', 'LineWidth', 1.3); hold on;
plot(t_txt, vd_txt, 'r--', 'LineWidth', 1.3);
xlabel('Time [s]'); ylabel('$v_d$ [pu]','Interpreter','latex');
title('Filtered $d$-axis Voltage','Interpreter','latex');
legend('MATLAB','PSCAD','Location','best'); grid on;

% Repeat similar process for vq_m if needed:
vq_f0 = x00(x_sol_str == "vq_f");
vq_m_absolute = vq_m + vq_f0;

subplot(2,1,2);
plot(t_out, vq_m_absolute, 'b', 'LineWidth', 1.3); hold on;
plot(t_txt, vq_txt, 'r--', 'LineWidth', 1.3);
xlabel('Time [s]'); ylabel('$v_q$ [pu]','Interpreter','latex');
title('Filtered $q$-axis Voltage','Interpreter','latex');
legend('MATLAB','PSCAD','Location','best'); grid on;
sgtitle('Filtered Voltages Comparison','Interpreter','latex');

% id_f0   = x00(x_sol_str == "id_f");   % Filtered d-axis current
% id_m = y(:,1);
% id_m_absolute  = id_m + id_f0;
% figure;
% subplot(2,1,1);
% plot(t_out, id_m_absolute, 'b', 'LineWidth', 1.3); hold on;
% plot(t_txt, id_txt, 'r--', 'LineWidth', 1.3);
% xlabel('Time [s]');
% ylabel('$i_d$ [pu]','Interpreter','latex');
% title('Filtered $d$-axis Current (Absolute)','Interpreter','latex');
% legend('MATLAB (Absolute)','PSCAD','Location','best');
% grid on;
% %% --- Figure 4: Inner current controller voltages (e_dCL, e_qCL) ---
% figure;
% subplot(2,1,1);
% plot(t_out, edCL_m, 'b', 'LineWidth', 1.3); hold on;
% plot(t_txt, edCL_txt, 'r--', 'LineWidth', 1.3);
% xlabel('Time [s]'); ylabel('$e_{d,CL}$ [pu]','Interpreter','latex');
% title('Current Controller $d$-axis Voltage Command','Interpreter','latex');
% legend('MATLAB','PSCAD','Location','best'); grid on;
% 
% subplot(2,1,2);
% plot(t_out, eqCL_m, 'b', 'LineWidth', 1.3); hold on;
% plot(t_txt, eqCL_txt, 'r--', 'LineWidth', 1.3);
% xlabel('Time [s]'); ylabel('$e_{q,CL}$ [pu]','Interpreter','latex');
% title('Current Controller $q$-axis Voltage Command','Interpreter','latex');
% legend('MATLAB','PSCAD','Location','best'); grid on;
% sgtitle('Inner Current Controller Voltage Commands','Interpreter','latex');

% %% Small-signal analysis
% % Eigenvalues
% names = string(model.StateName);
% [V,eigens] = eig(model.A); % Right eigenvectors (V) and eigenvalues, A = V*eigens*inv(V);
% [~,ind] = sort(abs(imag(diag(eigens))),'descend'); % Sort by frequency
% eigens = eigens(ind,ind);
% V = V(:,ind); % End sorting
% W = inv(V); % Left (W) eigenvectors | lambda = W*model.A*V
% pf = V.*W.'; % Participation factors (classic computation)
% 
% modes = diag(eigens);
% freq = abs(imag(modes))/(2*pi);
% damp = -real(modes)./abs(modes);
% 
% if sum(real(modes)>=0)>0
%     fprintf('\n Unstable system !! \n')
% else
%     fprintf('\n Stable system \n')
% end
% 
% figure('Name','Damping and frequency','NumberTitle','off')
% subplot(2,1,1)
% box on;hold on;grid on
% bar(1:1:length(modes),freq,'b');
% set(gca, 'yScale', 'log')
% ylim(10.^[floor(log10(min(freq(freq>1e-5)))), ceil(log10(max(freq)))])
% yticks([1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5]);
% xlabel('Modes'); ylabel("Frequency [Hz]");
% xticks(1:1:length(modes))
% subplot(2,1,2)
% box on;hold on;grid on
% bar(1:1:length(modes),damp,'b')
% xlabel('Modes'); ylabel("Damping ratio [pu]");
% xticks(1:1:length(modes))
% sgtitle(strcat("System modes: VSC + grid equivalent with $|Z_g| =$ ",num2str(sqrt(r_g^2+l_g^2),'%.2f')," pu"))
% set(gcf,'position',[50,50,1000,720])
