function [output1,output2] = runUWG(CL_EPW_PATH,CL_EPW,CL_XML_PATH,CL_XML)
    % =========================================================================
    %  THE URBAN WEATHER GENERATOR
    % =========================================================================
    % Orig. Author: B. Bueno & edited by A. Nakano & Lingfu Zhang
    % Last modified by Joseph Yang (joeyang@mit.edu) - May, 2016
    % 
    % Notes
    % a. When compiling, add 'z_meso.txt', 'RefDOE.mat','SchDef.m' to the list of files
    % b. Program description can be found in the following papers & websites
    %   - Joseph Yang's Master Thesis (2016)
    %   - https://github.com/hansukyang/UWG_Matlab
    % =========================================================================

    close all;
    ver = 4.1;
    % 4.1 (beta) updates
    %   - changed EPW output to rural wind speed again
    %   - compatibility with xml re-established
    %   - read in 'initialize.m' file for Matlab 
    
    % =========================================================================
    % Section 1 - Definitions for constants / other parameters
    % =========================================================================

    min_thickness = 0.01;   % Minimum layer thickness (to prevent crashing) (m)
    max_thickness = 0.05;   % Maximum layer thickness (m)
    soilTcond = 1;          % http://web.mit.edu/parmstr/Public/NRCan/nrcc29118.pdf (Figly & Snodgrass)
    soilvolHeat = 2e6;      % http://www.europment.org/library/2013/venice/bypaper/MFHEEF/MFHEEF-21.pdf (average taken from Table 1)
    soil = Material(soilTcond,soilvolHeat); % Soil material used for soil-depth padding

    
    % Physical parameters (moved here from Param.m)
    g = 9.81;               % gravity
    cp = 1004.;             % heat capacity for air (J/kg.K)
    vk = 0.40;              % von karman constant
    r = 287.;               % gas constant
    rv = 461.5;             %
    lv = 2.26e6;            % latent heat of evaporation
    sigma = 5.67e-08 ;      % Stefan Boltzmann constant
    waterDens = 1000;       % water density (kg/m^3)
    lvtt = 2.5008e6;        %
    tt = 273.16;            %
    estt = 611.14;          %
    cl = 4.218e3;           %
    cpv = 1846.1;           %
    b = 9.4;                % Coefficients derived by Louis (1979)
    cm = 7.4;               %
    colburn = (0.713/0.621)^(2./3.); % (Pr/Sc)^(2/3) for Colburn analogy in water evaporation
    
    % Site-specific parameters    
    wgmax = 0.005;          % maximum film water depth on horizontal surfaces (m)


    % =========================================================================
    % Section 2 - Read EPW file
    % =========================================================================
    %try
    climate_data = strcat(CL_EPW_PATH,'\',CL_EPW);
    fullyScripted = 1;
    
    disp(['Rural weather file selected: ',climate_data]);
    epwid = fopen(climate_data);
    C = importdata(climate_data, ',', 8);
    
    % Read header lines (1 to 8) from EPW and ensure TMY2 format
    % Note that TMY3 format is not compatible with the current version of UWG
    header = C.textdata(1:8,1);
    TMY3 = strfind(header(1), 'TMY3');
    
    if ~isempty(TMY3{1})
        disp('UWG unable to run: UWG requires TMY2 format for weather data');
        new_climate_file = 0;
        return;
    end
    
    % Read Lat, Long (line 1 of EPW)
    line1 = strsplit(header{1},',');
    lat = str2double(line1{7});
    lon = str2double(line1{8});
    GMT = str2double(line1{9}); 

    % Read in soil temperature data (assumes this is always there)
    soildata = strsplit(header{4},',');
    n_soil = str2double(soildata{2});
    Tsoil = zeros(n_soil,12);
    depth = zeros(n_soil,1);
    
    % Read monthly data for each layer of soil from EPW file
    for i = 1:n_soil
        depth(i) = str2double(soildata{3 + (i-1)*13});
        % Monthly data
        for j = 1:12
            Tsoil(i,j) = str2double(soildata{3+(i-1)*13+j})+273.15;
        end
    end
    
    % Read weather data from EPW for each time step in weather file
    i = 1;
    
    readin = fgetl(epwid);
    
    while (readin ~= -1)
        epwinput.values(i,:) = textscan(readin, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %*[^\n]', 'delimiter',',','MultipleDelimsAsOne',1);
        i=i+1;
        readin = fgetl(epwid);
    end
    epwinput.values(1:8,:) = []; % Ignore the header lines (first 8)
    
    
    fclose all;
    disp("finish");
    
    xml_location = strcat(CL_XML_PATH,'\',CL_XML);
    % Input files for UWG - note that soil layer buffering and 
    % layer thickness control are only performed for XML. (should update) 
    [~,~,ext] = fileparts(xml_location);
    
    %-----------------
    if strcmp(ext,'.m')
        % Run matlab script to generate UCM, UBL, etc.
        run(xml_location);
        nightStart = 18;        % arbitrary values (not used for XLSM)
        nightEnd = 8;

        simTime = SimParam(dtSim,dtWeather,Month,Day,nDay);
        weather = Weather(climate_data,simTime.timeInitial,simTime.timeFinal);
        forcIP = Forcing(weather.staTemp,weather); 
        forc = Forcing;

        % Road (Assume 0.5m of asphalt)
        emis = 0.93;
        asphalt = Material (1.0,1.6e6);
        thickness = 0.05 * ones (ceil(d_road/0.05),1);
        road = Element(alb_road,emis,thickness,[asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;asphalt;],0,293,1);
        road.vegCoverage = min(vegCover/(1-bldDensity),1);

        rural = road;
        rural.vegCoverage = rurVegCover;
        T_init = weather.staTemp(1);
        H_init = weather.staHum(1);

        geoParam = Param(h_ubl1,h_ubl2,h_ref,h_temp,h_wind,c_circ,maxDay,maxNight,...
            latTree,latGrss,albVeg,vegStart,vegEnd,nightStart,nightEnd,windMin,wgmax,c_exch,maxdx,...
            g, cp, vk, r, rv, lv, pi(), sigma, waterDens, lvtt, tt, estt, cl, cpv, b, cm, colburn);
        UBL = UBLDef('C',charLength,weather.staTemp(1),maxdx,geoParam.dayBLHeight,geoParam.nightBLHeight); 

        % Define BEM for each DOE type (read the fraction)
        load ('RefDOE.mat');

        % Define building energy models
        k = 0;
        r_glaze = 0;
        SHGC = 0;
        alb_wall = 0;
        area = bld*charLength^2*bldDensity*bldHeight/h_floor;  % building floor area
        for i = 1:16
            for j = 1:3
                if bld(i,j) > 0
                    k = k + 1;
                    BEM(k) = refBEM(i,j,zone);
                    BEM(k).frac = bld(i,j);
                    BEM(k).fl_area = area(i,j);
                    r_glaze = r_glaze + BEM(k).frac * BEM(k).building.glazingRatio;
                    SHGC = SHGC + BEM(k).frac * BEM(k).building.shgc;
                    alb_wall = alb_wall + BEM(k).frac * BEM(k).wall.albedo;
                    BEM(k).Qocc = BEM(k).Qocc;
                    Sch(k) = Schedule(i,j,zone);                    
                end
            end 
        end        
        
        UCM = UCMDef(bldHeight,bldDensity,verToHor,treeCoverage,...
            sensAnth,latAnth,T_init,H_init,weather.staUmod(1),geoParam,r_glaze,SHGC,alb_wall,road); 
        UCM.h_mix = h_mix;
        
        % Reference site class (also include VDM)
        RSM = RSMDef(lat,lon,GMT,h_obs,weather.staTemp(1),weather.staPres(1),geoParam);
        USM = RSMDef(lat,lon,GMT,bldHeight/10,weather.staTemp(1),weather.staPres(1),geoParam);
        
        % For .m file, assume the soil depth is close to one of the ground
        % soil depth specified in EPW (0.5, 1.0, 2.0)
        for i = 1:n_soil
            if sum(road.layerThickness) <= depth(i)
                soilindex1 = i;
                break;
            end
        end

        % Same for rural road
        for i = 1:n_soil
            if sum(rural.layerThickness) <= depth(i)
                soilindex2 = i;
                break;
            end
        end
        
    % =========================================================================
    % Section 6 - HVAC Autosizing (unlimited cooling & heating)
    % =========================================================================
    for j = 1:numel(BEM)
        if autosize
            BEM(j).building.coolCap = 9999;
            BEM(j).building.heatCap = 9999;  
        end
    end

    % =========================================================================
    % Section 7 - UWG main section
    % =========================================================================
    N = simTime.days * 24;
    n = 0;
    ph = simTime.dt/3600;       % per hour

    % Data dump variables
    time = transpose(1:1:simTime.days*24);
    WeatherData (N,1) = Forcing;
    UCMData (N,1) = UCMDef;
    UBLData (N,1) = UBLDef;
    RSMData (N,1) = RSMDef;
    USMData (N,1) = RSMDef;

    bTemp = zeros (N,numel(BEM));
    bRHum = zeros (N,numel(BEM));
    bPelec = zeros (N,numel(BEM));
    bQgas = zeros (N,numel(BEM));
    bPequip = zeros (N,numel(BEM));
    bPlight = zeros (N,numel(BEM));
    bQocc = zeros (N,numel(BEM));
    bFluxMass = zeros (N,numel(BEM));
    bFluxRoof = zeros(N,numel(BEM));
    bFluxWall = zeros (N,numel(BEM));
    bFluxSolar = zeros (N,numel(BEM));
    bFluxWindow = zeros (N,numel(BEM));
    bFluxInfil = zeros (N,numel(BEM));
    bFluxVent = zeros (N,numel(BEM));
    bCoolConsump = zeros (N,numel(BEM));
    bHeatConsump = zeros (N,numel(BEM));
    bCoolDemand = zeros (N,numel(BEM));
    bHeatDemand = zeros (N,numel(BEM));
    bTwallext = zeros (N,numel(BEM));
    bTroofext = zeros (N,numel(BEM));
    bTwallin = zeros (N,numel(BEM));
    bTroofin = zeros (N,numel(BEM));  
    bTmassin = zeros (N,numel(BEM));  
    bCOP = zeros(N,numel(BEM));
    bVent = zeros (N,numel(BEM));

    for it=1:(simTime.nt-1)

        % Update water temperature (estimated)
        if n_soil == 0
            forc.deepTemp = mean([forcIP.temp]);            % for BUBBLE/CAPITOUL/Singapore only
            forc.waterTemp = mean([forcIP.temp]) - 10;      % for BUBBLE/CAPITOUL/Singapore only
        else
            forc.deepTemp = Tsoil(soilindex1,simTime.month);
            forc.waterTemp = Tsoil(3,simTime.month);
        end
        
        % There's probably a better way to update the weather...
        simTime = UpdateDate(simTime);
        forc.infra = forcIP.infra(ceil(it*ph));       
        forc.wind = max(forcIP.wind(ceil(it*ph)),geoParam.windMin);     
        forc.uDir = forcIP.uDir(ceil(it*ph));
        forc.hum = forcIP.hum(ceil(it*ph));
        forc.pres = forcIP.pres(ceil(it*ph));
        forc.temp = forcIP.temp(ceil(it*ph));
        forc.rHum = forcIP.rHum(ceil(it*ph));
        forc.prec = forcIP.prec(ceil(it*ph));
        forc.dir = forcIP.dir(ceil(it*ph));
        forc.dif = forcIP.dif(ceil(it*ph));
        UCM.canHum = forc.hum;      % Canyon humidity (absolute) same as rural
        
        % Update solar flux
        [rural,UCM,BEM] = SolarCalcs(UCM,BEM,simTime,RSM,forc,geoParam,rural);
            
        % Update buildling & traffic schedule
        if strcmp(ext,'.xlsm') || strcmp(ext,'.m')

            % Assign day type (1 = weekday, 2 = sat, 3 = sun/other)
            if mod (simTime.julian,7) == 0      % Sunday
                dayType = 3;
            elseif mod (simTime.julian,7) == 6  % Saturday
                dayType = 2;
            else                                % Weekday
                dayType = 1;
            end

            % Update anthropogenic heat load for each hour (building & UCM)
            UCM.sensAnthrop = sensAnth*(SchTraffic(dayType,simTime.hourDay+1));
            
            for i = 1:numel(BEM)
                
                % Set temperature
                BEM(i).building.coolSetpointDay = Sch(i).Cool(dayType,simTime.hourDay+1) + 273.15;
                BEM(i).building.coolSetpointNight = BEM(i).building.coolSetpointDay;
                BEM(i).building.heatSetpointDay = Sch(i).Heat(dayType,simTime.hourDay+1) + 273.15;
                BEM(i).building.heatSetpointNight = BEM(i).building.heatSetpointDay;
                disp("firstbem elec");
                disp(BEM(i).Elec);
                disp(Sch(i).Qelec(i));
                disp("--fin--");
                % Internal Heat Load Schedule (W/m^2 of floor area for Q)
                BEM(i).Elec = Sch(i).Qelec*Sch(i).Elec(dayType,simTime.hourDay+1);
                BEM(i).Light = Sch(i).Qlight*Sch(i).Light(dayType,simTime.hourDay+1);
                BEM(i).Nocc = Sch(i).Nocc*Sch(i).Occ(dayType,simTime.hourDay+1);
                BEM(i).Qocc = sensOcc*(1-LatFOcc)*BEM(i).Nocc;

                % SWH and ventilation schedule
                BEM(i).SWH = Sch(i).Vswh*Sch(i).SWH(dayType,simTime.hourDay+1);     % litres per hour / m^2 of floor space
                BEM(i).building.vent = Sch(i).Vent;                                 % m^3/s/m^2 of floor
                BEM(i).Gas = Sch(i).Qgas * Sch(i).Gas(dayType,simTime.hourDay+1);   % Gas Equip Schedule, per m^2 of floor

                % This is quite messy, should update
                intHeat = BEM(i).Light+BEM(i).Elec+BEM(i).Qocc;
                BEM(i).building.intHeatDay = intHeat;
                BEM(i).building.intHeatNight = intHeat;
                BEM(i).building.intHeatFRad = (RadFLight *BEM(i).Light + RadFEquip*BEM(i).Elec)/intHeat;
                BEM(i).building.intHeatFLat = LatFOcc*sensOcc*BEM(i).Nocc/intHeat;
                
                BEM(i).T_wallex = BEM(i).wall.layerTemp(1);
                BEM(i).T_wallin = BEM(i).wall.layerTemp(end);
                BEM(i).T_roofex = BEM(i).roof.layerTemp(1);
                BEM(i).T_roofin = BEM(i).roof.layerTemp(end);
            end
        end
    end
    output1 = 1;
    output2 = 2;
end