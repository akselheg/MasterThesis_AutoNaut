function loadCSVtoMat(file)
    switch file
        case 'AngularVelocity.csv'
            data  = read_csv_file('AngularVelocity.csv');
            AngularVelocity.timestamp = data.timestamp_secondsSince01_01_1970_;
            av_mask = all(char(data.entity{:}) == 'STB - ADIS',2);
            gps_r_mask = all(char(data.entity{:}) == 'GPS       ',2);
            AngularVelocity.x = data.x_rad_s_(av_mask);
            AngularVelocity.y = data.y_rad_s_(av_mask);
            AngularVelocity.z = data.z_rad_s_(av_mask);
            period = mean(diff(data.timestamp_secondsSince01_01_1970_(av_mask)));
            AngularVelocity.timestamp = data.timestamp_secondsSince01_01_1970_;
            AngularVelocity.freq = 1/period;
            AngularVelocity.z_gps = data.z_rad_s_(gps_r_mask);
            save('AngularVelocity.mat', 'AngularVelocity','-v7.3')
        case 'Acceleration.csv'
            data  = read_csv_file('Acceleration.csv');
            Acceleration.timestamp = data.timestamp_secondsSince01_01_1970_;
            Acceleration.x = data.x_m_s_s_;
            Acceleration.y = data.y_m_s_s_;
            Acceleration.z = data.z_m_s_s_;
            period = mean(diff(data.timestamp_secondsSince01_01_1970_));
            Acceleration.freq = 1/period;
            save('Acceleration.mat', 'Acceleration','-v7.3')
    end
end