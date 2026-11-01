function writeImuCsv(fileName, imu)
%WRITEIMUCSV Write IMU measurements (time, f [m/s^2], w [rad/s], body frame).
fid = fopen(fileName, 'w');
fprintf(fid, 'Time [s],fx [m/s^2],fy [m/s^2],fz [m/s^2],wx [rad/s],wy [rad/s],wz [rad/s]\n');
fprintf(fid, '%.6f,%.9f,%.9f,%.9f,%.12e,%.12e,%.12e\n', [imu.t imu.f imu.w]');
fclose(fid);
end
