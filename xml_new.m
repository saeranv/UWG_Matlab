function [new_climate_file] = xml_new(CL_EPW_PATH,CL_EPW,CL_XML_PATH,CL_XML,CL_RE_PATH,CL_RE)
    % =========================================================================
    %  THE URBAN WEATHER GENERATOR
    % =========================================================================
    % Orig. Author: B. Bueno & edited by A. Nakano & Lingfu Zhang
    % Last modified by Joseph Yang (joeyang@mit.edu)
    % latest modification 2015-11-18
    % Changes
    % 1. Minor typo fixed ('Wa11' to 'Wall')
    % 2. Shortened code and added a section to break up single material
    % 3. Re-added section to break up thick layers (>5cm)
    % 4. (TBD) Modify functionality to either read from XML or matlab file
    % =========================================================================
    
    % =========================================================================
    % Section 1 - Definitions for constants / other parameters
    % =========================================================================
    % echo off;
    sim_dt = 300;           % Simulation time step (s)
    weather_dt = 3600;      % Weather data time step (EPW) (s)
    fullyScripted = 1;
    autosize = 1;           % HAVC autosizing

    % Physical parameters
    g = 9.81;               % gravity
    cp = 1004.;             % heat capacity for air (constant pressure)
    vk = 0.40;              % von karman constant
    r = 287.;               % gas constant
    rv = 461.5;
    lv = 2.26e6;            % latent heat of evaporation
    sigma = 5.67e-08 ;      % Stefan Boltzmann constant
    waterDens = 1000;       % water density
    lvtt = 2.5008e6;
    tt = 273.16;
    estt = 611.14;
    cl = 4.218e3;
    cpv = 1846.1;
    b   = 9.4;              % Coefficients derived by Louis (1979)
    cm  = 7.4; 
    colburn = (0.713/0.621)^(2./3.); % (Pr/Sc)^(2/3) for Colburn analogy in water evaporation

    % =========================================================================
    % Section 2 - Read EPW file
    % =========================================================================
    try
        climate_data = strcat(CL_EPW_PATH,'\',CL_EPW);
        epwPathName = CL_EPW_PATH;
        epwFileName = CL_EPW;
        [pathstr,name,ext] = fileparts(climate_data);
        epwFileExt = ext;
    catch
        [epwFileName,epwPathName] = uigetfile('.epw','Select Rural EnergyPlus Weather File');
        climate_data = strcat(epwPathName,epwFileName);
        epwPathName = epwPathName(1:end-1);
    end

    disp(['Rural weather file selected: ',climate_data])
    epwid = fopen(climate_data);

    delimiterIn = ',';
    headerlinesIn = 8;
    C = importdata(climate_data, delimiterIn, headerlinesIn);

    for i = 1:8
        header(i) = C.textdata(i,1);
    end

    i = 1;
    while 1
        readin = fgetl(epwid);   
        if readin == -1 
            break;
        end
        epwinput.values(i,:) = textscan(readin, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %*[^\n]', 'delimiter',',','MultipleDelimsAsOne',1);
        i=i+1;
    end
    epwinput.values(1:8,:) = [];
    fclose all;

    % =========================================================================
    % Section 3 - Read XML Input File
    % =========================================================================
    try
        xml_location = strcat(CL_XML_PATH,'\',CL_XML);
    catch
        [FileName,PathName] = uigetfile('.xml','Select Urban Parameter Input file');
        xml_location = strcat(PathName,FileName);
    end

    xml_input = xml_read(xml_location);
    [pathstr,name,ext] = fileparts(xml_location);
    disp(['Urban Parameter file selected: ',xml_location]);

    % Re-naming for file readability
    xmlTyp(1) = xml_input.typology1;
    xmlTyp(2) = xml_input.typology2;
    xmlTyp(3) = xml_input.typology3;
    xmlTyp(4) = xml_input.typology4;
    building(1) = xmlTyp(1).building;
    building(2) = xmlTyp(2).building;
    building(3) = xmlTyp(3).building;
    building(4) = xmlTyp(4).building;
    xmlParam = xml_input.parameter;
    xmlUArea = xml_input.urbanArea;
    xmlRSite = xml_input.referenceSite;

    % =========================================================================
    % Section 4 - Get save location for the new EPW file
    % =========================================================================
    try
        new_climate_file = strcat(CL_RE_PATH,'\',CL_RE);
        newPathName = CL_RE_PATH;
        newFileName_withExt = CL_RE;    
    catch
        [newFileName_withExt,newPathName] = uiputfile('.epw','Select save location');
        new_climate_file = strcat(newPathName,newFileName_withExt);    
        newPathName = newPathName(1:end-1);
    end

    [new_pathstr,new_name,new_ext] = fileparts(new_climate_file);
    newFileName = new_name;
    disp(['Save location selected: ',new_climate_file]);

    % =========================================================================
    % Section 5 - Read XML file for urban area definition
    % =========================================================================
    % Reconstruct xml Materials (removed for now - unsure regarding methodology)

    % Simulation paramters
    simParam = SimParam(sim_dt,weather_dt,xmlParam.simuStartMonth,...
        xmlParam.simuStartDay, xmlParam.simuDuration);
    weather = Weather(climate_data,simParam.timeInitial,simParam.timeFinal);
    nightStart = mean ([building.nightSetStart]);
    nightEnd = mean ([building.nightSetEnd]);
    parameter = Param(xmlUArea.daytimeBLHeight,xmlUArea.nighttimeBLHeight,...
        xmlUArea.refHeight,xmlParam.tempHeight,xmlParam.windHeight,...
        1.2,200,50,xmlUArea.treeLatent,xmlUArea.grassLatent,xmlUArea.vegAlbedo,...
        xml_input.urbanArea.vegStart,xml_input.urbanArea.vegEnd,...
        nightStart,nightEnd,0.1,10,0.005,0.3,500,...
        g, cp, vk, r, rv, lv, pi(), sigma, waterDens, lvtt, tt, estt, cl,...
        cpv, b, cm, colburn);

    % Define Road Element
    roadMat = [];
    newthickness = [];
    materials = xmlUArea.urbanRoad.materials;
    urbanRoad = xmlUArea.urbanRoad;
    k = xmlUArea.urbanRoad.materials.thermalConductivity;
    Vhc = xmlUArea.urbanRoad.materials.volumetricHeatCapacity;
    if numel(materials.thickness)>1
        for j = 1:numel(materials.thickness)
            % Break up each layer that's more than 5cm thick
            if materials.thickness(j) > 0.05
                nlayers = ceil(materials.thickness(i)/0.05);
                for l = 1:nlayers
                    roadMat = [roadMat Material(k{j},Vhc{j})];
                    newthickness = [newthickness materials.thickness(j)/nlayers];
                end
            else
                roadMat = [roadMat Material(k{j},Vhc{j})];
                newthickness = [newthickness materials.thickness(j)];
            end
        end
    else
        % Divide single layer into two (UWG assumes at least 2 layers)
        newthickness = [materials.thickness/2 materials.thickness/2];
        roadMat = [Material(k,Vhc) Material(k,Vhc)];
    end
    road = Element(urbanRoad.albedo,urbanRoad.emissivity,newthickness,roadMat,...
        urbanRoad.vegetationCoverage,urbanRoad.initialTemperature + 273.15,xmlUArea.urbanRoad.inclination);

    % Define Rural Element
    ruralMat = [];
    newthickness = [];
    materials = xmlRSite.ruralRoad.materials;
    ruralRoad = xmlRSite.ruralRoad;
    k = xmlUArea.urbanRoad.materials.thermalConductivity;
    Vhc = xmlUArea.urbanRoad.materials.volumetricHeatCapacity;
    if numel(materials.thickness)>1
        for j = 1:numel(materials.thickness)
            % Break up each layer that's more than 5cm thick
            if materials.thickness(j) > 0.05
                nlayers = ceil(materials.thickness(j)/0.05);
                for l = 1:nlayers
                    ruralMat = [ruralMat Material(k{j},Vhc{j})];
                    newthickness = [newthickness materials.thickness(j)/nlayers];
                end
            else
                ruralMat = [wallMat Material(k{j},Vhc{j})];
                newthickness = [newthickness materials.thickness(j)];
            end
        end
    else
        newthickness = [materials.thickness/2 materials.thickness/2];    
        ruralMat = [Material(k,Vhc) Material(k,Vhc)];
    end
    rural = Element(ruralRoad.albedo,ruralRoad.emissivity,newthickness,...
        ruralMat,ruralRoad.vegetationCoverage,ruralRoad.initialTemperature + 273.15,ruralRoad.inclination);

    for i = 1:4
        % Define Wall
        wallMat = [];
        newthickness = [];
        materials = xmlTyp(i).construction.wall.materials;
        xwall = xmlTyp(i).construction.wall;
        k = materials.thermalConductivity;
        Vhc = materials.volumetricHeatCapacity;

        if numel(materials.thickness) > 1
            for j = 1:numel(materials.thickness)
                % Break up each layer that's more than 5cm thick
                if materials.thickness(j) > 0.05
                    nlayers = ceil(materials.thickness(j)/0.05);
                    for l = 1:nlayers
                        wallMat = [wallMat Material(k{j},Vhc{j})];
                        newthickness = [newthickness materials.thickness(j)/nlayers];
                    end
                else
                    wallMat = [wallMat Material(k{j},Vhc{j})];
                    newthickness = [newthickness materials.thickness(j)];
                end
            end
        else
            newthickness = [materials.thickness/2 materials.thickness/2];
            wallMat = [Material(k, Vhc) Material(k, Vhc)];
        end
        wall(i) = Element(xwall.albedo,xwall.emissivity,newthickness,wallMat,...
            xwall.vegetationCoverage,xwall.initialTemperature + 273.15,xwall.inclination);

        % Define Roof element
        roofMat = [];
        newthickness = [];
        materials = xmlTyp(i).construction.roof.materials;
        xroof = xmlTyp(i).construction.roof;
        k = materials.thermalConductivity;
        Vhc = materials.volumetricHeatCapacity;

        if numel(materials.thickness) > 1
            for j = 1:numel(materials.thickness)
                % Break up each layer that's more than 5cm thick
                if materials.thickness(j) > 0.05
                    nlayers = ceil(materials.thickness(j)/0.05);
                    for l = 1:nlayers
                        roofMat = [roofMat Material(k{j},Vhc{j})];
                        newthickness = [newthickness materials.thickness(j)/nlayers];
                    end
                else
                    roofMat = [roofMat Material(k{j},Vhc{j})];
                    newthickness = [newthickness materials.thickness(j)];
                end
            end
        else
            newthickness = [materials.thickness/2 materials.thickness/2];
            roofMat = [Material(k,Vhc) Material(k,Vhc)];
        end
        roof(i) = Element(xroof.albedo,xroof.emissivity,newthickness,roofMat,...
            xroof.vegetationCoverage,xroof.initialTemperature + 273.15,xroof.inclination);

        % Define Mass
        massMat = [];
        newthickness = [];
        materials = xmlTyp(i).construction.mass.materials;
        xmass = xmlTyp(i).construction.mass;
        k = materials.thermalConductivity;
        Vhc = materials.volumetricHeatCapacity;

        if numel(materials.thickness) > 1
            for j = 1:numel(materials.thickness)
                % Break up each layer that's more than 5cm thick
                if materials.thickness(j) > 0.05
                    nlayers = ceil(materials.thickness(j)/0.05);
                    for l = 1:nlayers
                        massMat = [massMat Material(k{j},Vhc{j})];
                        newthickness = [newthickness materials.thickness(j)/nlayers];
                    end
                else
                    massMat = [massMat Material(k{j},Vhc{j})];
                    newthickness = [newthickness materials.thickness(j)];
                end
            end
        else
            newthickness = [materials.thickness/2 materials.thickness/2];
            massMat = [Material(k,Vhc) Material(k,Vhc)];
        end
        mass(i) = Element(xmass.albedo,xmass.emissivity,newthickness,massMat,...
            xmass.vegetationCoverage,xmass.initialTemperature + 273.15,xmass.inclination);

        % Define building typology
        typology(i) = Building(xmlTyp(i).building.floorHeight,...
            xmlTyp(i).building.nightInternalGains,...
            xmlTyp(i).building.dayInternalGains,...
            xmlTyp(i).building.radiantFraction,...
            xmlTyp(i).building.latentFraction,...
            xmlTyp(i).building.infiltration,...
            xmlTyp(i).building.ventilation,...
            xmlTyp(i).construction.glazing.glazingRatio,...
            xmlTyp(i).construction.glazing.windowUvalue,...
            xmlTyp(i).construction.glazing.windowSHGC,...
            xmlTyp(i).building.coolingSystemType,...
            xmlTyp(i).building.coolingCOP,...
            xmlTyp(i).building.heatReleasedToCanyon,...
            xmlTyp(i).building.daytimeCoolingSetPoint + 273.15,...
            xmlTyp(i).building.nighttimeCoolingSetPoint + 273.15,...
            xmlTyp(i).building.daytimeHeatingSetPoint + 273.15,...
            xmlTyp(i).building.nighttimeHeatingSetPoint + 273.15,...
            xmlTyp(i).building.coolingCapacity,...
            xmlTyp(i).building.heatingEfficiency,...
            xmlTyp(i).building.initialT + 273.15);

        % Urban Configuration [building,mass,wall,roof,road]
        urbanConf(i) = UrbanConf(typology(i),mass(i),wall(i),roof(i),road);
    end 

    % Reference site class
    refSite = ReferenceSite(xmlRSite.latitude,xmlRSite.longitude,xmlRSite.averageObstacleHeight,...
        weather.staTemp(1),weather.staPres(1),parameter);

    % Urban Usage [Fraction of urban configurations,urban configurations]
    typDist = [xmlTyp(1).dist xmlTyp(2).dist xmlTyp(3).dist xmlTyp(4).dist]/100;
    urbanUsage = UrbanUsage(typDist,[urbanConf(1), urbanConf(2), urbanConf(3), urbanConf(4)]);   

    % Define dominant typology
    [M, I] = max(typDist);
    dominantTypology = typology(I);
    dominantWall = wall(I);
    disp(dominantTypology);

    urbanArea = UrbanArea(xmlUArea.averageBuildingHeight,xmlUArea.siteCoverageRatio,... 
        xmlUArea.facadeToSiteRatio,xmlUArea.treeCoverage, xmlUArea.nonBldgSensibleHeat,...               
        xmlUArea.nonBldgLatentAnthropogenicHeat,weather.staTemp(1),weather.staHum(1),...
        weather.staUmod(1),parameter,dominantTypology,dominantWall,road,rural); 
    ublVars = UblVars('C',xml_input.urbanArea.charLength,weather.staTemp(1),parameter.maxdx); 
    
    save ('UWGdump.mat');

    % =========================================================================
    % Section 6 - HVAC Autosizing (if needed)
    % =========================================================================
    Fc = fopen('coolCap.txt','r+');
    for i = 1:numel(urbanArea)
        for j = 1:numel(urbanUsage(i).urbanConf)
            if autosize
                coolCap = Autosize(urbanArea(i),ublVars(i),...
                    urbanUsage(i).urbanConf(j),climate_data,rural,refSite,parameter);
                fprintf(Fc,'%1.3f\n',coolCap);
                urbanUsage(i).urbanConf(j).building.coolCap = coolCap;
            else
                urbanUsage(i).urbanConf(j).building.coolCap = fscanf(Fc,'%f\n',1);
            end
        end
    end
    fclose(Fc);

    % =========================================================================
    % Section 7 - UWG main section
    % =========================================================================
    forc = Forcing(weather.staTemp);
    sensHeat = zeros(numel(urbanArea),1);

    timeCount = 0.;
    Can_Tdb = [];
    Can_hum = [];
    for it=1:simParam.nt
        timeCount=timeCount+simParam.dt;
        simParam = UpdateDate(simParam);
        % print file
        if eq(mod(timeCount,simParam.timePrint),0)
            for i = 1:numel(urbanArea)
                Can_Tdb = [Can_Tdb (urbanArea(i).canTemp-273.15)];
                Can_hum = [Can_hum urbanArea(i).canHum];
                radTemp = [urbanUsage(1).urbanConf(1).wall.layerTemp(1) forc.skyTemp  urbanUsage(1).urbanConf(1).road.layerTemp(1)];
                radTemp_calc = [radTemp, urbanArea(1).bldHeight urbanArea(1).canWidth];
                filename = strcat(newPathName,'\','Trad_',newFileName,'.csv');
                dlmwrite(filename,radTemp_calc,'delimiter',',','-append');
            end 
        end

        % start of Bruno's code
        % read forcing
        if eq(mod(timeCount,simParam.timeForcing),0) || eq(timeCount,simParam.dt)
            if le(forc.itfor,simParam.timeMax/simParam.timeForcing)
                forc = ReadForcing(forc,weather,parameter);
                % solar calculations
                [rural,urbanArea,urbanUsage] = SolarCalcs(urbanArea,urbanUsage,...
                  simParam,refSite,forc,parameter,rural);
                [first,second,third] = PriorityIndex(forc.uDir,ublVars);
            end
        end
        % rural heat fluxes
        rural.infra = forc.infra-parameter.sigma*rural.layerTemp(1)^4.;
        rural = SurfFlux( rural,forc,parameter,simParam,forc.hum,forc.temp,...
            forc.wind,2,0.);
        % vertical profiles of meteorological variables at the rural site
        refSite = VerticalDifussionModel( refSite,forc,...
            rural,parameter,simParam );
        for i = 1:numel(urbanArea)
            % urban heat fluxes
            [urbanArea(i),ublVars(i),urbanUsage(i),forc] = UrbFlux(urbanArea(i),ublVars(i),...
            urbanUsage(i),forc,parameter,simParam,refSite);
            % urban canyon temperature and humidity
            urbanArea(i) = UrbThermal(urbanArea(i),ublVars(i).ublTemp,urbanUsage(i),forc,parameter,simParam);
        end
        % urban boundary layer temperature
        for i=1:numel(urbanArea)
            sensHeat(i) = urbanArea(i).sensHeat;
        end
        ublVars = UrbanBoundaryLayerModel(ublVars,sensHeat,refSite,rural,forc,...
            parameter,simParam,first,second,third);
        % print progress
        if eq(mod(timeCount,3600.*24),0)
          progressbar(timeCount/(3600.*24)/365);
        end
    end

    % =========================================================================
    % Section 8 - Writing new EPW file
    % =========================================================================
    disp('Calculating new Temperature and humidity values')
    for iJ = 1:size(Can_Tdb,2)
        epwinput.values{iJ,7}{1,1} = num2str(Can_Tdb(iJ),'%0.1f'); % dry bulb temperature  [�C]
        [Tdb, w, Can_phi(iJ), h, Can_Tdp(iJ), v] = Psychrometrics(Can_Tdb(iJ), Can_hum(iJ), str2num(epwinput.values{iJ,10}{1,1}));
        epwinput.values{iJ,8}{1,1} = num2str(Can_Tdp(iJ),'%0.1f'); % dew point temperature [�C]
        epwinput.values{iJ,9}{1,1} = num2str(Can_phi(iJ),'%0.0f'); % relative humidity     [%]
    end
    progressbar(360/365);
    disp('writing new EPW file');

    new_climate_file = strcat(newPathName,'\',newFileName,'.epw');
    epwnewid = fopen(new_climate_file,'w');
    for i = 1:8
        fprintf(epwnewid,'%s\r\n',header{i});
    end
    for i = 1:size(epwinput.values,1)
        printme = [];
        for e = 1:34
            printme = [printme epwinput.values{i,e}{1,1} ','];
        end
        printme = [printme epwinput.values{i,e}{1,1}];
        fprintf(epwnewid,'%s\r\n',printme);
    end

    progressbar(1);
    if fullyScripted
        disp('Inputs scripted, supressing pop-up notification...');
    else
        h = msgbox('Urban Weather Generation Complete','UWG 2.0','help');
    end
    disp(['New climate file generated: ',new_climate_file]);

    save ('UWGdump.mat');
    fclose all;
end