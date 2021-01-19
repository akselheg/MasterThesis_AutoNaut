%%
clearvars;close all;clc;


startDate = [2020, 07, 01];
endDate = [2020, 07, 02];
%place = "vestlandet";
place = "midtnorge";
ch = char(place);
ch = string(ch(1));
% https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800vhf/mywavewam800_vestlandet.an.2019052818.nc
% https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800vhf/mywavewam800_midtnorge.an.2020021918.nc
                  %     /thredds/dodsC/fou-hi/mywavewam800mhf/mywavewam800_midtnorge.an.2020022018.nc
startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    previousDate = daysadd(startDateString,i-1); 
    numberPreviousDate = yyyymmdd(datetime(datestr(previousDate)));
    if i == 0
        waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place +".an." + string(numberPreviousDate) + "18.nc";
        waveInfo = ncinfo(waveData);
        
        latitudeMapWave = ncread(waveData,'latitude');
        longitudeMapWave = ncread(waveData,'longitude');
        waveDimentions = [size(latitudeMapWave) , 24*(numDays+1) + 18];
        
        waveSize = zeros(waveDimentions);
        waveDir = zeros(waveDimentions);
        windDir = zeros(waveDimentions); 
        waveHZ = zeros(waveDimentions); 
        windSpeed = zeros(waveDimentions); 
        
        waveSize(:,:,1:18) = ncread(waveData,'hs', [1 1 7], [inf inf inf]);
        waveDir(:,:,1:18) = ncread(waveData,'thq', [1 1 7], [inf inf inf]); % I believe thq is total wave dir
        windDir(:,:,1:18) = ncread(waveData,'dd', [1 1 7], [inf inf inf]);
        waveHZ(:,:,1:18) = ncread(waveData,'tp', [1 1 7], [inf inf inf]);
        windSpeed(:,:,1:18) = ncread(waveData,'ff', [1 1 7], [inf inf inf]);
        
    end
    waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place + ".an." + string(numberDate) + "18.nc";
    
    waveSize(:,:,24*i+7:24*i+30) = ncread(waveData,'hs');
    waveDir(:,:,24*i+19:24*i+42) = ncread(waveData,'thq');  % I believe thq is total wave dir 
    windDir(:,:,24*i+19:24*i+42) = ncread(waveData,'dd');
    waveHZ(:,:,24*i+19:24*i+42) = ncread(waveData,'tp');
    windSpeed(:,:,24*i+19:24*i+42) = ncread(waveData,'ff');
end
startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "weatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'windSpeed','waveHZ','waveSize','waveDir', 'latitudeMapWave', 'longitudeMapWave','windDir', 'waveInfo','-v7.3');
pause(10);
disp('Done');
%%
clearvars;close all;
pause(3);

startDate = [2020, 07, 03];
endDate = [2020, 07, 04];
place = "midtnorge";
ch = char(place);
ch = string(ch(1));

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    previousDate = daysadd(startDateString,i-1); 
    numberPreviousDate = yyyymmdd(datetime(datestr(previousDate)));
    if i == 0
        waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place +".an." + string(numberPreviousDate) + "18.nc";
        waveInfo = ncinfo(waveData);
        
        latitudeMapWave = ncread(waveData,'latitude');
        longitudeMapWave = ncread(waveData,'longitude');
        waveDimentions = [size(latitudeMapWave) , 24*(numDays+1) + 18];
        
        waveSize = zeros(waveDimentions);
        waveDir = zeros(waveDimentions);
        windDir = zeros(waveDimentions); 
        waveHZ = zeros(waveDimentions); 
        windSpeed = zeros(waveDimentions); 
        
        waveSize(:,:,1:18) = ncread(waveData,'hs', [1 1 7], [inf inf inf]);
        waveDir(:,:,1:18) = ncread(waveData,'thq', [1 1 7], [inf inf inf]); % I believe thq is total wave dir
        windDir(:,:,1:18) = ncread(waveData,'dd', [1 1 7], [inf inf inf]);
        waveHZ(:,:,1:18) = ncread(waveData,'tp', [1 1 7], [inf inf inf]);
        windSpeed(:,:,1:18) = ncread(waveData,'ff', [1 1 7], [inf inf inf]);
        
    end
    waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place + ".an." + string(numberDate) + "18.nc";
    
    waveSize(:,:,24*i+7:24*i+30) = ncread(waveData,'hs');
    waveDir(:,:,24*i+19:24*i+42) = ncread(waveData,'thq');  % I believe thq is total wave dir 
    windDir(:,:,24*i+19:24*i+42) = ncread(waveData,'dd');
    waveHZ(:,:,24*i+19:24*i+42) = ncread(waveData,'tp');
    windSpeed(:,:,24*i+19:24*i+42) = ncread(waveData,'ff');
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "weatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'windSpeed','waveHZ','waveSize','waveDir', 'latitudeMapWave', 'longitudeMapWave','windDir', 'waveInfo','-v7.3');
pause(10);
disp('Done2');
%% 
clearvars;close all;
pause(3);

startDate = [2020, 07, 05];
endDate = [2020, 07, 05];
place = "midtnorge";
ch = char(place);
ch = string(ch(1));

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    previousDate = daysadd(startDateString,i-1); 
    numberPreviousDate = yyyymmdd(datetime(datestr(previousDate)));
    if i == 0
        waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place +".an." + string(numberPreviousDate) + "18.nc";
        waveInfo = ncinfo(waveData);
        
        latitudeMapWave = ncread(waveData,'latitude');
        longitudeMapWave = ncread(waveData,'longitude');
        waveDimentions = [size(latitudeMapWave) , 24*(numDays+1) + 18];
        
        waveSize = zeros(waveDimentions);
        waveDir = zeros(waveDimentions);
        windDir = zeros(waveDimentions); 
        waveHZ = zeros(waveDimentions); 
        windSpeed = zeros(waveDimentions); 
        
        waveSize(:,:,1:18) = ncread(waveData,'hs', [1 1 7], [inf inf inf]);
        waveDir(:,:,1:18) = ncread(waveData,'thq', [1 1 7], [inf inf inf]); % I believe thq is total wave dir
        windDir(:,:,1:18) = ncread(waveData,'dd', [1 1 7], [inf inf inf]);
        waveHZ(:,:,1:18) = ncread(waveData,'tp', [1 1 7], [inf inf inf]);
        windSpeed(:,:,1:18) = ncread(waveData,'ff', [1 1 7], [inf inf inf]);
        
    end
    waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place + ".an." + string(numberDate) + "18.nc";
    
    waveSize(:,:,24*i+7:24*i+30) = ncread(waveData,'hs');
    waveDir(:,:,24*i+19:24*i+42) = ncread(waveData,'thq');  % I believe thq is total wave dir 
    windDir(:,:,24*i+19:24*i+42) = ncread(waveData,'dd');
    waveHZ(:,:,24*i+19:24*i+42) = ncread(waveData,'tp');
    windSpeed(:,:,24*i+19:24*i+42) = ncread(waveData,'ff');
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "weatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'windSpeed','waveHZ','waveSize','waveDir', 'latitudeMapWave', 'longitudeMapWave','windDir', 'waveInfo','-v7.3');
pause(10);
disp('Done3');
%%
clearvars;close all;
pause(3)

startDate = [2020, 07, 06];
endDate = [2020, 07, 06];
place = "midtnorge";
ch = char(place);
ch = string(ch(1));

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    previousDate = daysadd(startDateString,i-1); 
    numberPreviousDate = yyyymmdd(datetime(datestr(previousDate)));
    if i == 0
        waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place +".an." + string(numberPreviousDate) + "18.nc";
        waveInfo = ncinfo(waveData);
        
        latitudeMapWave = ncread(waveData,'latitude');
        longitudeMapWave = ncread(waveData,'longitude');
        waveDimentions = [size(latitudeMapWave) , 24*(numDays+1) + 18];
        
        waveSize = zeros(waveDimentions);
        waveDir = zeros(waveDimentions);
        windDir = zeros(waveDimentions); 
        waveHZ = zeros(waveDimentions); 
        windSpeed = zeros(waveDimentions); 
        
        waveSize(:,:,1:18) = ncread(waveData,'hs', [1 1 7], [inf inf inf]);
        waveDir(:,:,1:18) = ncread(waveData,'thq', [1 1 7], [inf inf inf]); % I believe thq is total wave dir
        windDir(:,:,1:18) = ncread(waveData,'dd', [1 1 7], [inf inf inf]);
        waveHZ(:,:,1:18) = ncread(waveData,'tp', [1 1 7], [inf inf inf]);
        windSpeed(:,:,1:18) = ncread(waveData,'ff', [1 1 7], [inf inf inf]);
        
    end
    waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place + ".an." + string(numberDate) + "18.nc";
    
    waveSize(:,:,24*i+7:24*i+30) = ncread(waveData,'hs');
    waveDir(:,:,24*i+19:24*i+42) = ncread(waveData,'thq');  % I believe thq is total wave dir 
    windDir(:,:,24*i+19:24*i+42) = ncread(waveData,'dd');
    waveHZ(:,:,24*i+19:24*i+42) = ncread(waveData,'tp');
    windSpeed(:,:,24*i+19:24*i+42) = ncread(waveData,'ff');
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "weatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'windSpeed','waveHZ','waveSize','waveDir', 'latitudeMapWave', 'longitudeMapWave','windDir', 'waveInfo','-v7.3');
pause(10);
disp('Done4');
%%
clearvars;close all;
pause(3);

startDate = [2020, 07, 09];
endDate = [2020, 07, 09];
place = "midtnorge";
ch = char(place);
ch = string(ch(1));

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    previousDate = daysadd(startDateString,i-1); 
    numberPreviousDate = yyyymmdd(datetime(datestr(previousDate)));
    if i == 0
        waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place +".an." + string(numberPreviousDate) + "18.nc";
        waveInfo = ncinfo(waveData);
        
        latitudeMapWave = ncread(waveData,'latitude');
        longitudeMapWave = ncread(waveData,'longitude');
        waveDimentions = [size(latitudeMapWave) , 24*(numDays+1) + 18];
        
        waveSize = zeros(waveDimentions);
        waveDir = zeros(waveDimentions);
        windDir = zeros(waveDimentions); 
        waveHZ = zeros(waveDimentions); 
        windSpeed = zeros(waveDimentions); 
        
        waveSize(:,:,1:18) = ncread(waveData,'hs', [1 1 7], [inf inf inf]);
        waveDir(:,:,1:18) = ncread(waveData,'thq', [1 1 7], [inf inf inf]); % I believe thq is total wave dir
        windDir(:,:,1:18) = ncread(waveData,'dd', [1 1 7], [inf inf inf]);
        waveHZ(:,:,1:18) = ncread(waveData,'tp', [1 1 7], [inf inf inf]);
        windSpeed(:,:,1:18) = ncread(waveData,'ff', [1 1 7], [inf inf inf]);
        
    end
    waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place + ".an." + string(numberDate) + "18.nc";
    
    waveSize(:,:,24*i+7:24*i+30) = ncread(waveData,'hs');
    waveDir(:,:,24*i+19:24*i+42) = ncread(waveData,'thq');  % I believe thq is total wave dir 
    windDir(:,:,24*i+19:24*i+42) = ncread(waveData,'dd');
    waveHZ(:,:,24*i+19:24*i+42) = ncread(waveData,'tp');
    windSpeed(:,:,24*i+19:24*i+42) = ncread(waveData,'ff');
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "weatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'windSpeed','waveHZ','waveSize','waveDir', 'latitudeMapWave', 'longitudeMapWave','windDir', 'waveInfo','-v7.3');
pause(10);
disp('Done5');
%%
clearvars;close all;
pause(3);

startDate = [2020, 05, 28];
endDate = [2020, 05, 28];
place = "vestlandet";
ch = char(place);
ch = string(ch(1));

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    previousDate = daysadd(startDateString,i-1); 
    numberPreviousDate = yyyymmdd(datetime(datestr(previousDate)));
    if i == 0
        waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place +".an." + string(numberPreviousDate) + "18.nc";
        waveInfo = ncinfo(waveData);
        
        latitudeMapWave = ncread(waveData,'latitude');
        longitudeMapWave = ncread(waveData,'longitude');
        waveDimentions = [size(latitudeMapWave) , 24*(numDays+1) + 18];
        
        waveSize = zeros(waveDimentions);
        waveDir = zeros(waveDimentions);
        windDir = zeros(waveDimentions); 
        waveHZ = zeros(waveDimentions); 
        windSpeed = zeros(waveDimentions); 
        
        waveSize(:,:,1:18) = ncread(waveData,'hs', [1 1 7], [inf inf inf]);
        waveDir(:,:,1:18) = ncread(waveData,'thq', [1 1 7], [inf inf inf]); % I believe thq is total wave dir
        windDir(:,:,1:18) = ncread(waveData,'dd', [1 1 7], [inf inf inf]);
        waveHZ(:,:,1:18) = ncread(waveData,'tp', [1 1 7], [inf inf inf]);
        windSpeed(:,:,1:18) = ncread(waveData,'ff', [1 1 7], [inf inf inf]);
        
    end
    waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place + ".an." + string(numberDate) + "18.nc";
    
    waveSize(:,:,24*i+7:24*i+30) = ncread(waveData,'hs');
    waveDir(:,:,24*i+19:24*i+42) = ncread(waveData,'thq');  % I believe thq is total wave dir 
    windDir(:,:,24*i+19:24*i+42) = ncread(waveData,'dd');
    waveHZ(:,:,24*i+19:24*i+42) = ncread(waveData,'tp');
    windSpeed(:,:,24*i+19:24*i+42) = ncread(waveData,'ff');
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "weatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'windSpeed','waveHZ','waveSize','waveDir', 'latitudeMapWave', 'longitudeMapWave','windDir', 'waveInfo','-v7.3');
pause(10);
disp('Done6');
%%
clearvars;close all;
pause(3);

startDate = [2020, 02, 20];
endDate = [2020, 02, 20];
place = "midtnorge";
ch = char(place);
ch = string(ch(1));

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    previousDate = daysadd(startDateString,i-1); 
    numberPreviousDate = yyyymmdd(datetime(datestr(previousDate)));
    if i == 0
        waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place +".an." + string(numberPreviousDate) + "18.nc";
        waveInfo = ncinfo(waveData);
        
        latitudeMapWave = ncread(waveData,'latitude');
        longitudeMapWave = ncread(waveData,'longitude');
        waveDimentions = [size(latitudeMapWave) , 24*(numDays+1) + 18];
        
        waveSize = zeros(waveDimentions);
        waveDir = zeros(waveDimentions);
        windDir = zeros(waveDimentions); 
        waveHZ = zeros(waveDimentions); 
        windSpeed = zeros(waveDimentions); 
        
        waveSize(:,:,1:18) = ncread(waveData,'hs', [1 1 7], [inf inf inf]);
        waveDir(:,:,1:18) = ncread(waveData,'thq', [1 1 7], [inf inf inf]); % I believe thq is total wave dir
        windDir(:,:,1:18) = ncread(waveData,'dd', [1 1 7], [inf inf inf]);
        waveHZ(:,:,1:18) = ncread(waveData,'tp', [1 1 7], [inf inf inf]);
        windSpeed(:,:,1:18) = ncread(waveData,'ff', [1 1 7], [inf inf inf]);
        
    end
    waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place + ".an." + string(numberDate) + "18.nc";
    
    waveSize(:,:,24*i+7:24*i+30) = ncread(waveData,'hs');
    waveDir(:,:,24*i+19:24*i+42) = ncread(waveData,'thq');  % I believe thq is total wave dir 
    windDir(:,:,24*i+19:24*i+42) = ncread(waveData,'dd');
    waveHZ(:,:,24*i+19:24*i+42) = ncread(waveData,'tp');
    windSpeed(:,:,24*i+19:24*i+42) = ncread(waveData,'ff');
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "weatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'windSpeed','waveHZ','waveSize','waveDir', 'latitudeMapWave', 'longitudeMapWave','windDir', 'waveInfo','-v7.3');
pause(10);
disp('DoneDone');