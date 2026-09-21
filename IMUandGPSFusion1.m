%% IMU and GPS Fusion for Inertial Navigation 
clc, clear, close all
format longG
rng(0)
%% 
imuFs = 160;
gpsFs = 1;
refloc = [42.2825 -72.3430 53.0352];
imuSamplesPerGPS = imuFs / gpsFs;
Ts = 1 / imuFs;
load LoggedQuadcopter.mat trajData;
trajOrient = trajData.Orientation;
trajVel    = trajData.Velocity;
trajPos    = trajData.Position;
trajAcc    = trajData.Acceleration;
trajAngVel = trajData.AngularVelocity;
gps = gpsSensor('UpdateRate', gpsFs, 'ReferenceLocation', refloc);
gps.DecayFactor = 0.5;
gps.HorizontalPositionAccuracy = 1.6;
gps.VerticalPositionAccuracy   = 1.6;
gps.VelocityAccuracy = 0.1;
imu = imuSensor('accel-gyro-mag', 'SampleRate', imuFs);
imu.MagneticField = [19.5281 -5.0741 48.0067];
fusionfilt = insfilterMARG;
fusionfilt.IMUSampleRate = imuFs;
fusionfilt.ReferenceLocation = refloc;
initstate = zeros(22,1);
initstate(1:4) = compact(meanrot(trajOrient(1:100)));
initstate(5:7)  = mean(trajPos(1:100,:),1);
initstate(8:10) = mean(trajVel(1:100,:),1);
initstate(11:13) = imu.Gyroscope.ConstantBias / imuFs;
initstate(14:16) = imu.Accelerometer.ConstantBias / imuFs;
initstate(17:19) = imu.MagneticField;
initstate(20:22) = imu.Magnetometer.ConstantBias;
fusionfilt.State = initstate;    
initP = 1e-9 * eye(22);
fusionfilt.StateCovariance = initP;
secondsToSimulate = 50;
loopBound = secondsToSimulate * imuFs;
loopBound = floor(loopBound / imuFs) * imuFs;
load mat_accel
load mat_gyro
load mat_mag
load mat_lla
load mat_gpsvel
pqorient_mat = zeros(160*50, 4); 
pqpos = zeros(160*50, 3);
state = single(initstate);
P = single(initP);
inp = [];
outp = [];
for block = 1:50
    start_idx   = (block-1) * imuSamplesPerGPS + 1;
    end_idx     = block * imuSamplesPerGPS;
    accel_block = single(mat_accel(start_idx:end_idx, :));
    gyro_block  = single(mat_gyro(start_idx:end_idx, :));
    mag_block   = single(mat_mag(start_idx:end_idx, :));
    gps_pos     = single(mat_lla(block, :));
    gps_vel     = single(mat_gpsvel(block, :));
    %%
    inp = [inp state' reshape(P',1,[]) reshape(accel_block',1,[]) reshape(gyro_block',1,[]) reshape(mag_block',1,[]) gps_pos gps_vel]; 
    [state, P, orient_block, pos_block] = main_process1(state,P,accel_block,gyro_block,mag_block,gps_pos,gps_vel);
    outp = [outp state' reshape(P',1,[]) reshape(single(orient_block)',1,[]) reshape(single(pos_block)',1,[])];
    pqorient_mat(start_idx:end_idx, :) = orient_block;
    pqpos(start_idx:end_idx, :) = pos_block;
end
pqorient = quaternion(pqorient_mat);
%% result
posd  = pqpos(1:loopBound,:) - trajPos(1:loopBound,:);
quatd = rad2deg(dist(pqorient(1:loopBound), trajOrient(1:loopBound)));
msep  = sqrt(mean(posd.^2,1));
fprintf('X =  %.2f m\n', msep(1));
fprintf('Y =  %.2f m\n', msep(2));
fprintf('Z =  %.2f m\n', msep(3));
fprintf('angle = %.2f degrees \n', sqrt(mean(quatd.^2)));
%% save
save_file(inp,'inp.bin');
save_file(outp,'outp.bin');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_file(x,name)
fid = fopen(name, 'wb');
fwrite(fid,x, 'float32');
fclose(fid);
end


