clearvars;close all;clc;


startDate = [2020, 07, 01];
endDate = [2020, 07, 03];

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    windCurrentData =  "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/NorKyst-800m_ZDEPTHS_his.an." + string(numberDate) + "00.nc";
    if i == 0
        CurrentInfo = ncinfo(windCurrentData);
        latitudeCurrentMap = ncread(windCurrentData,'lat');
        longitudeCurrentMap = ncread(windCurrentData,'lon');
        windCurrentDimentions = [size(latitudeCurrentMap), 24*(numDays+1)];
        currentEast = zeros(windCurrentDimentions);
        currentNorth = zeros(windCurrentDimentions);
    end
    for curhour = 0 : 23
        currentNorth(:,:,24*i+curhour+1) = ncread(windCurrentData,'v_northward',[endDate1 1 1 1+curhour],[inf inf 1 1]);
        currentEast(:,:,24*i+curhour+1) = ncread(windCurrentData,'u_eastward',[1 1 1 1+curhour],[inf inf 1 1]);
    end
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "currentweatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'currentNorth', 'currentEast','longitudeCurrentMap', 'latitudeCurrentMap', 'CurrentInfo','-v7.3');
disp('Done')
%%
clearvars;close all;clc;


startDate = [2020, 07, 03];
endDate = [2020, 07, 04];

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    windCurrentData =  "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/NorKyst-800m_ZDEPTHS_his.an." + string(numberDate) + "00.nc";
    if i == 0
        CurrentInfo = ncinfo(windCurrentData);
        latitudeCurrentMap = ncread(windCurrentData,'lat');
        longitudeCurrentMap = ncread(windCurrentData,'lon');
        windCurrentDimentions = [size(latitudeCurrentMap), 24*(numDays+1)];
        currentEast = zeros(windCurrentDimentions);
        currentNorth = zeros(windCurrentDimentions);
    end
    for curhour = 0 : 23
        currentNorth(:,:,24*i+curhour+1) = ncread(windCurrentData,'v_northward',[1 1 1 1+curhour],[inf inf 1 1]);
        currentEast(:,:,24*i+curhour+1) = ncread(windCurrentData,'u_eastward',[1 1 1 1+curhour],[inf inf 1 1]);
    end
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "currentweatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'currentNorth', 'currentEast','longitudeCurrentMap', 'latitudeCurrentMap', 'CurrentInfo','-v7.3');
disp('Done')
%%
clearvars;close all;clc;


startDate = [2020, 07, 04];
endDate = [2020, 07, 05];

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    windCurrentData =  "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/NorKyst-800m_ZDEPTHS_his.an." + string(numberDate) + "00.nc";
    if i == 0
        CurrentInfo = ncinfo(windCurrentData);
        latitudeCurrentMap = ncread(windCurrentData,'lat');
        longitudeCurrentMap = ncread(windCurrentData,'lon');
        windCurrentDimentions = [size(latitudeCurrentMap), 24*(numDays+1)];
        currentEast = zeros(windCurrentDimentions);
        currentNorth = zeros(windCurrentDimentions);
    end
    for curhour = 0 : 23
        currentNorth(:,:,24*i+curhour+1) = ncread(windCurrentData,'v_northward',[1 1 1 1+curhour],[inf inf 1 1]);
        currentEast(:,:,24*i+curhour+1) = ncread(windCurrentData,'u_eastward',[1 1 1 1+curhour],[inf inf 1 1]);
    end
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "currentweatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'currentNorth', 'currentEast','longitudeCurrentMap', 'latitudeCurrentMap', 'CurrentInfo','-v7.3');
disp('Done')
%%
clearvars;close all;clc;


startDate = [2020, 07, 05];
endDate = [2020, 07, 05];

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    windCurrentData =  "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/NorKyst-800m_ZDEPTHS_his.an." + string(numberDate) + "00.nc";
    if i == 0
        CurrentInfo = ncinfo(windCurrentData);
        latitudeCurrentMap = ncread(windCurrentData,'lat');
        longitudeCurrentMap = ncread(windCurrentData,'lon');
        windCurrentDimentions = [size(latitudeCurrentMap), 24*(numDays+1)];
        currentEast = zeros(windCurrentDimentions);
        currentNorth = zeros(windCurrentDimentions);
    end
    for curhour = 0 : 23
        currentNorth(:,:,24*i+curhour+1) = ncread(windCurrentData,'v_northward',[1 1 1 1+curhour],[inf inf 1 1]);
        currentEast(:,:,24*i+curhour+1) = ncread(windCurrentData,'u_eastward',[1 1 1 1+curhour],[inf inf 1 1]);
    end
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "currentweatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'currentNorth', 'currentEast','longitudeCurrentMap', 'latitudeCurrentMap', 'CurrentInfo','-v7.3');
disp('Done')
%%
clearvars;close all;clc;


startDate = [2020, 07, 06];
endDate = [2020, 07, 06];

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    windCurrentData =  "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/NorKyst-800m_ZDEPTHS_his.an." + string(numberDate) + "00.nc";
    if i == 0
        CurrentInfo = ncinfo(windCurrentData);
        latitudeCurrentMap = ncread(windCurrentData,'lat');
        longitudeCurrentMap = ncread(windCurrentData,'lon');
        windCurrentDimentions = [size(latitudeCurrentMap), 24*(numDays+1)];
        currentEast = zeros(windCurrentDimentions);
        currentNorth = zeros(windCurrentDimentions);
    end
    for curhour = 0 : 23
        currentNorth(:,:,24*i+curhour+1) = ncread(windCurrentData,'v_northward',[1 1 1 1+curhour],[inf inf 1 1]);
        currentEast(:,:,24*i+curhour+1) = ncread(windCurrentData,'u_eastward',[1 1 1 1+curhour],[inf inf 1 1]);
    end
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "currentweatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'currentNorth', 'currentEast','longitudeCurrentMap', 'latitudeCurrentMap', 'CurrentInfo','-v7.3');
disp('done')
%%
clearvars;close all;clc;


startDate = [2020, 07, 09];
endDate = [2020, 07, 09];

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    windCurrentData =  "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/NorKyst-800m_ZDEPTHS_his.an." + string(numberDate) + "00.nc";
    if i == 0
        CurrentInfo = ncinfo(windCurrentData);
        latitudeCurrentMap = ncread(windCurrentData,'lat');
        longitudeCurrentMap = ncread(windCurrentData,'lon');
        windCurrentDimentions = [size(latitudeCurrentMap), 24*(numDays+1)];
        currentEast = zeros(windCurrentDimentions);
        currentNorth = zeros(windCurrentDimentions);
    end
    for curhour = 0 : 23
        currentNorth(:,:,24*i+curhour+1) = ncread(windCurrentData,'v_northward',[1 1 1 1+curhour],[inf inf 1 1]);
        currentEast(:,:,24*i+curhour+1) = ncread(windCurrentData,'u_eastward',[1 1 1 1+curhour],[inf inf 1 1]);
    end
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "currentweatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'currentNorth', 'currentEast','longitudeCurrentMap', 'latitudeCurrentMap', 'CurrentInfo','-v7.3');
%%
clearvars;close all;clc;


startDate = [2020, 02, 20];
endDate = [2020, 02, 20];

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    windCurrentData =  "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/NorKyst-800m_ZDEPTHS_his.an." + string(numberDate) + "00.nc";
    if i == 0
        CurrentInfo = ncinfo(windCurrentData);
        latitudeCurrentMap = ncread(windCurrentData,'lat');
        longitudeCurrentMap = ncread(windCurrentData,'lon');
        windCurrentDimentions = [size(latitudeCurrentMap), 24*(numDays+1)];
        currentEast = zeros(windCurrentDimentions);
        currentNorth = zeros(windCurrentDimentions);
    end
    for curhour = 0 : 23
        currentNorth(:,:,24*i+curhour+1) = ncread(windCurrentData,'v_northward',[1 1 1 1+curhour],[inf inf 1 1]);
        currentEast(:,:,24*i+curhour+1) = ncread(windCurrentData,'u_eastward',[1 1 1 1+curhour],[inf inf 1 1]);
    end
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "currentweatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'currentNorth', 'currentEast','longitudeCurrentMap', 'latitudeCurrentMap', 'CurrentInfo','-v7.3');
disp('Done')
%%
clearvars;close all;clc;


startDate = [2020, 05, 28];
endDate = [2020, 05, 29];

startDateString = string(startDate(1)) + "/" + string(startDate(2)) + "/" + string(startDate(3));
endDateString = string(endDate(1)) + "/" + string(endDate(2)) + "/" + string(endDate(3));
numDays = daysact(startDateString, endDateString);
for i = 0:numDays
    date = daysadd(startDateString,i);
    str = datestr(date);
    test = datetime(str);
    numberDate = yyyymmdd(datetime(datestr(date)));
    windCurrentData =  "https://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-1h/NorKyst-800m_ZDEPTHS_his.an." + string(numberDate) + "00.nc";
    if i == 0
        CurrentInfo = ncinfo(windCurrentData);
        latitudeCurrentMap = ncread(windCurrentData,'lat');
        longitudeCurrentMap = ncread(windCurrentData,'lon');
        windCurrentDimentions = [size(latitudeCurrentMap), 24*(numDays+1)];
        currentEast = zeros(windCurrentDimentions);
        currentNorth = zeros(windCurrentDimentions);
    end
    for curhour = 0 : 23
        currentNorth(:,:,24*i+curhour+1) = ncread(windCurrentData,'v_northward',[1 1 1 1+curhour],[inf inf 1 1]);
        currentEast(:,:,24*i+curhour+1) = ncread(windCurrentData,'u_eastward',[1 1 1 1+curhour],[inf inf 1 1]);
    end
end

startDateSaveFormat = string(startDate(1)) + "-" + string(startDate(2)) + "-" + string(startDate(3));
endDateSaveFormat =  string(endDate(1)) + "-" + string(endDate(2)) + "-" + string(endDate(3));
savename = "currentweatherData_"+ startDateSaveFormat + "_" + endDateSaveFormat + ".mat"; 
save(savename, 'currentNorth', 'currentEast','longitudeCurrentMap', 'latitudeCurrentMap', 'CurrentInfo','-v7.3');
%%