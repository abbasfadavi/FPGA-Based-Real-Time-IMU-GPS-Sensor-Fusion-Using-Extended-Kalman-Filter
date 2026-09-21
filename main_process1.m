function [state, P, orient_block, pos_block] = main_process1(state, P, accel_block, gyro_block, mag_block, gps_pos, gps_vel)
numSamples = 160;
orient_block = zeros(numSamples, 4, 'single');
pos_block    = zeros(numSamples, 3, 'single');

for imuSample = 1:numSamples
    accel = accel_block(imuSample, :);
    gyro  = gyro_block(imuSample, :);
    [state, P, fusedOrient_vec, fusedPos] = imu_predict(state, P, accel, gyro);
    orient_block(imuSample, :) = fusedOrient_vec;
    pos_block(imuSample, :) = fusedPos;
end

[state, P] = gps_mag_update(state, P, gps_pos, gps_vel, mag_block(end,:));
end
%% ------------------------------------------------------------------------
function [state, P, fusedOrient_vec, fusedPos] = imu_predict(state, P, accel, gyro)
Ts = 0.00625;
state = single(state);
P = single(P);
g_ned = single([0; 0; 9.80665]);
Ts = single(Ts);
% IMU prediction step (single precision)
q     = state(1:4);
pos   = state(5:7);
vel   = state(8:10);
b_g   = state(11:13);
b_a   = state(14:16);
omega = gyro.' - b_g;
dtheta = omega * Ts;

half_theta = single(0.5) * dtheta;
half_norm2 = sum(half_theta.^2, 1);
dq_w   = single(1) - single(0.5) * half_norm2;

% Quaternion multiplication (update)
q_new = [ dq_w   * q(1) - half_theta(1)*q(2) - half_theta(2)*q(3) - half_theta(3)*q(4);
          dq_w   * q(2) + half_theta(1)*q(1) + half_theta(3)*q(3) - half_theta(2)*q(4);
          dq_w   * q(3) - half_theta(3)*q(2) + half_theta(2)*q(1) + half_theta(1)*q(4);
          dq_w   * q(4) + half_theta(3)*q(1) - half_theta(1)*q(3) + half_theta(2)*q(2)];

% Normalize quaternion
q_norm = sqrt(sum(q_new.^2, 1));   
q_new = q_new / q_norm;
% Rotation matrix from quaternion
q0 = q_new(1);
q1 = q_new(2);
q2 = q_new(3);
q3 = q_new(4);
Rb2n = [ q0^2+q1^2-q2^2-q3^2, 2*(q1*q2 - q0*q3),   2*(q1*q3 + q0*q2);
               2*(q1*q2 + q0*q3),   q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
               2*(q1*q3 - q0*q2),   2*(q2*q3 + q0*q1),   q0^2 - q1^2 - q2^2 + q3^2 ];

% Acceleration in body frame minus bias
f_body = accel.' - b_a;
% Transform to navigation frame and add gravity
f_ned  = Rb2n * f_body;
a_ned = f_ned + g_ned;

% Velocity and position update (trapezoidal)
vel_new = vel + a_ned * Ts;
pos_new = pos + (vel + single(0.5)*(vel_new - vel)) * Ts;

% Update state
state(1:4) = q_new;
state(5:7)  = pos_new;
state(8:10) = vel_new;

% Diagonal covariance increase (process noise)
diag_increase = single([1e-6 * ones(4,1);   % quat
                        2e-2 * ones(3,1);   % pos
                        5e-2 * ones(3,1);   % vel
                        1e-4 * ones(3,1);   % gyro bias
                        5e-3 * ones(3,1);   % acc bias
                        1e-6 * ones(6,1)]); % mag parts
P = P + diag(diag_increase);

% Outputs
fusedOrient_vec = q_new.';   % row vector 1x4
fusedPos = pos_new.';

end
%% ------------------------------------------------------------------------
function [state, P] = gps_mag_update(state, P, gps_pos, gps_vel, mag)
Rpos = single(5.169);      
Rvel = single(0.0051);     
Rmag = single(0.0862); 
% GPS update
H_gps = zeros(6,22,'single');
H_gps(1:3, 5:7) = eye(3,'single');
H_gps(4:6, 8:10) = eye(3,'single');
R_gps = blkdiag(Rpos * eye(3,'single'), Rvel * eye(3,'single'));
w1 = H_gps * P;
w2 = w1 * H_gps';
S = w2' + R_gps;
X = P * H_gps';          
Y = X';                 
Kt = chol_solve(S, Y);  
K = Kt';
w3 = K * H_gps;
w3 = eye(22,'single') - w3;
w4 = w3 * P * w3';
w5 = K * R_gps * K';
P = w4 + w5;
%%
nedPos_meas = gps_pos2ned_simple(gps_pos); 
residual_pos = nedPos_meas(:) - state(5:7);

residual_vel = gps_vel(:) - state(8:10);
z_gps = [residual_pos; residual_vel];
w6 = K * z_gps; 
state = state + w6;
% Magnetometer update with analytical Jacobian
H_mag = analytic_mag_jacobian(state);
w7 = H_mag * P;
w8 = w7 * H_mag';
SS = w8 + Rmag * eye(3,'single');
XX = P * H_mag';          
YY = XX';   
Ktt = chol_solve(SS,YY);
KK = Ktt';
q = state(1:4);
m_ned = state(17:19);
b_mag = state(20:22);
Rb2n = quat2rotmat(q);
R_ned2body = Rb2n';

pred_mag = R_ned2body * m_ned + b_mag;
z_mag = mag(:) - pred_mag;
w9 = KK * z_mag;
state = state + w9;
%P = (eye(22,'single') - KK * H_mag) * P * (eye(22,'single') - KK * H_mag)' + KK * (Rmag * eye(3,'single')) * KK';
% Normalize quaternion
state(1:4) = state(1:4) / sqrt(sum(state(1:4).^2,1));

end
%% ------------------------------------------------------------------------
function X = chol_solve(S, Y)
L = chol(S, 'lower');     
Z = L \ Y;
X = L' \ Z;
end
%% ------------------------------------------------------------------------
function ned = gps_pos2ned_simple(gps_pos)
refloc = single([42.2825 -72.3430 53.0352]);
% Convert gps_pos to NED using spherical Earth approximation (single)
R = single(6371000);       % Earth radius [m]

lat0 = deg2rad(single(refloc(1)));
lon0 = deg2rad(single(refloc(2)));
alt0 = single(refloc(3));

lat = deg2rad(single(gps_pos(1)));
lon = deg2rad(single(gps_pos(2)));
alt = single(gps_pos(3));

dy = (lat - lat0) * R;          % north
dx = (lon - lon0) * cos(lat0) * R;  % east
dz = alt0 - alt;                 % down
 
ned = [dy, dx, dz];  
end
%% ------------------------------------------------------------------------
function H = analytic_mag_jacobian(x)
% Analytical Jacobian of magnetometer measurement (single)
q = x(1:4);
m_ned = x(17:19);

Rb2n = quat2rotmat(q);
R_ned2body = Rb2n';
q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

% Derivatives of Rb2n w.r.t. quaternion components
% ∂R/∂q0
dR_q0 = [ 2*q0, -2*q3,  2*q2;
          2*q3,  2*q0, -2*q1;
         -2*q2,  2*q1,  2*q0 ];
% ∂R/∂q1
dR_q1 = [ 2*q1,  2*q2,  2*q3;
          2*q2, -2*q1, -2*q0;
          2*q3,  2*q0, -2*q1 ];
% ∂R/∂q2
dR_q2 = [ -2*q2,  2*q1,  2*q0;
           2*q1,  2*q2,  2*q3;
          -2*q0,  2*q3, -2*q2 ];
% ∂R/∂q3
dR_q3 = [ -2*q3, -2*q0,  2*q1;
           2*q0, -2*q3,  2*q2;
           2*q1,  2*q2,  2*q3 ];

% Derivatives of R_ned2body are transposes
dR_ned_q0 = dR_q0';
dR_ned_q1 = dR_q1';
dR_ned_q2 = dR_q2';
dR_ned_q3 = dR_q3';

% Columns for quaternion part: (dR_ned/dq_i) * m_ned
col1 = dR_ned_q0 * m_ned(:);
col2 = dR_ned_q1 * m_ned(:);
col3 = dR_ned_q2 * m_ned(:);
col4 = dR_ned_q3 * m_ned(:);

H = zeros(3,22,'single');
H(:,1:4) = [col1, col2, col3, col4];
H(:,17:19) = R_ned2body;      % ∂/∂m_ned
H(:,20:22) = eye(3,'single'); % ∂/∂b_mag
end
%% ------------------------------------------------------------------------
function R = quat2rotmat(q)
% Quaternion to rotation matrix (single)
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
R = [ q0^2+q1^2-q2^2-q3^2, 2*(q1*q2 - q0*q3),   2*(q1*q3 + q0*q2);
      2*(q1*q2 + q0*q3),   q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
      2*(q1*q3 - q0*q2),   2*(q2*q3 + q0*q1),   q0^2 - q1^2 - q2^2 + q3^2 ];
end