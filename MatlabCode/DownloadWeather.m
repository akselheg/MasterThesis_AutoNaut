

a = TryDownloadWeather(1,1,1,1,1,1)


function wind = TryDownloadWeather(x0,y0,x1,y1, day,time)
    wind = 1;
    waveData = "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800"+ ch +"hf/mywavewam800_" + place +".an." + string(numberPreviousDate) + "18.nc";
    waveInfo = ncinfo(waveData);    
    latitudeMapWave = ncread(waveData,'latitude');
end