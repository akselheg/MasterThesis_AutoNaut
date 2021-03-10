function data = load_data(path,id,names)
    save_directory = pwd;
    % Check if data already exists
    if exist([save_directory '/' id '_data.mat']) == 2
        load([id '_data.mat']);
    else
        data = struct;
        for i=1:length(names)
            name = names{i};
            switch name
                case 'AngularVelocity'
                    AngularVelocity = read_csv_file('AngularVelocity.csv', path);

%                     data.angvelo_sampl_freq = round(10.^4*(AngularVelocity.timestamp_secondsSince01_01_1970_(2,:)-AngularVelocity.timestamp_secondsSince01_01_1970_(1,:)))/10.^4;
                    data.angvelo_sampl_freq = mean(diff(AngularVelocity.timestamp_secondsSince01_01_1970_));
                    data.angvelo_sampl_freq_hz = 1/data.angvelo_sampl_freq;

                    angvelo_timest_original = AngularVelocity.timestamp_secondsSince01_01_1970_-data.t0;

                    av_mask = all(char(AngularVelocity.entity{:}) == 'STB - ADIS',2);
                    gps_r_mask = all(char(AngularVelocity.entity{:}) == 'GPS       ',2);
                    data.p = AngularVelocity.x_rad_s_(av_mask);
                    data.q = AngularVelocity.y_rad_s_(av_mask);
                    data.r = AngularVelocity.z_rad_s_(av_mask);
                    data.r_gps = AngularVelocity.z_rad_s_(gps_r_mask);
                    data.angvelo_timest_original = angvelo_timest_original(av_mask);

                case 'EulerAngles'
                    EulerAngles = read_csv_file('EulerAngles.csv', path);
                    data.eu_sampl_freq = round(10.^2*(EulerAngles.timestamp_secondsSince01_01_1970_(2,:)-EulerAngles.timestamp_secondsSince01_01_1970_(1,:)))/10.^2;
%                     data.eu_sampl_freq = mean(diff(EulerAngles.timestamp_secondsSince01_01_1970_));
                    data.eu_sampl_freq_hz = 1/data.eu_sampl_freq;
                    data.eu_timest_original = round(EulerAngles.timestamp_secondsSince01_01_1970_-data.t0,2);

                    data.phi = EulerAngles.phi_rad_;
                    data.theta = EulerAngles.theta_rad_;
                    data.psi = EulerAngles.psi_rad_;
                    data.psi_magnetic = EulerAngles.psi_magnetic_rad_;

                case 'GpsFix' 
                    GpsFix = read_csv_file('GpsFix.csv', path);
                    data.t0 = GpsFix.timestamp_secondsSince01_01_1970_(1);
                    data.gps_sampl_freq = round(10.^2*(GpsFix.timestamp_secondsSince01_01_1970_(2,:)-GpsFix.timestamp_secondsSince01_01_1970_(1,:)))/10.^2;
                    data.gps_sampl_freq_hz = 1/data.gps_sampl_freq;
                    data.gps_timest_original = round(GpsFix.timestamp_secondsSince01_01_1970_-GpsFix.timestamp_secondsSince01_01_1970_(1),2);

                    data.cog = GpsFix.cog_rad_;
                    data.sog = GpsFix.sog_m_s_;
                    data.lat = GpsFix.lat_rad_;
                    data.lon = GpsFix.lon_rad_;
                    data.h = GpsFix.height_m_;
                    
                case 'Acceleration'
                    Acceleration = read_csv_file('Acceleration.csv', path);
%                     data.acc_sampl_freq = round(10.^4*(Acceleration.timestamp_secondsSince01_01_1970_(2,:)-Acceleration.timestamp_secondsSince01_01_1970_(1,:)))/10.^4;
                    data.acc_sampl_freq = mean(diff(Acceleration.timestamp_secondsSince01_01_1970_));
                    data.acc_sampl_freq_hz = 1/data.acc_sampl_freq;
                    data.acc_x = Acceleration.x_m_s_s_;
                    data.acc_y = Acceleration.y_m_s_s_;
                    data.acc_z = Acceleration.z_m_s_s_;
                    data.acc_timest_original = Acceleration.timestamp_secondsSince01_01_1970_-data.t0;
                case 'SetServoPosition'
                    SetServoPosition = read_csv_file('SetServoPosition.csv', path);
                    data.servo_sample_freq = mean(diff(SetServoPosition.timestamp_secondsSince01_01_1970_-SetServoPosition.timestamp_secondsSince01_01_1970_(1)));
                    data.servo_sample_freq_hz = 1/data.servo_sample_freq;
                    data.servo_timest_original = SetServoPosition.timestamp_secondsSince01_01_1970_-data.t0;
                    data.delta = SetServoPosition.value_rad_;
                case 'Heave'
                    Heave = read_csv_file('Heave',path);

                    data.heave_sampl_freq = round(10.^4*(Heave.timestamp_secondsSince01_01_1970_(2,:)-Heave.timestamp_secondsSince01_01_1970_(1,:)))/10.^4;
                    data.heave_sampl_freq_hz = 1/data.heave_sampl_freq;

                    heave_timest_original = Heave.timestamp_secondsSince01_01_1970_-Heave.timestamp_secondsSince01_01_1970_(1);

                    heave_mask = all(char(Heave.entity{:}) == 'GPS      ',2);
                    data.heave = Heave.value_m_(heave_mask);
                    data.heave_timest_original = heave_timest_original(heave_mask);
                case 'AbsoluteWind'
                    AbsoluteWind = read_csv_file('AbsoluteWind',path);

                    data.wind_sampl_freq = round(10.^4*(AbsoluteWind.timestamp_secondsSince01_01_1970_(2,:)-AbsoluteWind.timestamp_secondsSince01_01_1970_(1,:)))/10.^4;
                    data.widn_sampl_freq_hz = 1/data.wind_sampl_freq;
                    data.wind_timest_original = AbsoluteWind.timestamp_secondsSince01_01_1970_-AbsoluteWind.timestamp_secondsSince01_01_1970_(1);
                    % Upsample to 2Hz
                    t = [data.wind_timest_original(1):0.5:data.wind_timest_original(end)];
                    data.wind_dir = interp1(data.wind_timest_original,AbsoluteWind.dir___,t);
                    data.wind_speed = interp1(data.wind_timest_original,AbsoluteWind.speed_m_s_,t);
                    data.wind_timest_original = t;
    %                 data.wind_dir = AbsoluteWind.dir___;
    %                 data.wind_speed = AbsoluteWind.speed_m_s_;
                case 'DesiredHeading'
                    DesiredHeading = read_csv_file('DesiredHeading',path);
                    data.chi_d_sample_freq = mean(diff(DesiredHeading.timestamp_secondsSince01_01_1970_-DesiredHeading.timestamp_secondsSince01_01_1970_(1)));
                    data.chi_d_sample_freq_hz = 1/data.chi_d_sample_freq;
                    data.chi_d_timest_original = DesiredHeading.timestamp_secondsSince01_01_1970_-data.t0;
                    data.chi_d = DesiredHeading.value_rad_;
                    % Upsample to 2Hz
%                     t = [data.desired_heading_timest_original(1):0.5:data.desired_heading_timest_original(end)];
%                     data.desired_heading = interp1(data.desired_heading_timest_original,DesiredHeading.value_rad_,t);
%                     data.desired_heading_timest_original = t;
            end     
        end
        save([id '_data.mat'],'data')
    end
end

