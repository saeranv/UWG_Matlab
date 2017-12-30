function [new_climate_file] = xml_inputs_outputTrad4(CL_EPW_PATH,CL_EPW,CL_XML_PATH,CL_XML,CL_RE_PATH,CL_RE,IN_MON,IN_DAY,IN_DUR)
% =========================================================================
%  THE URBAN WEATHER GENERATOR
% % Generate new EPW file
% =========================================================================
% Author: B. Bueno
% Packaged by J. Sullivan-Fedock
% Modified by A. Nakano & Lingfu Zhang
% latest modification 2014-10-29
% -------------------------------------------------------------------------

fullyScripted = 1;
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
    fullyScripted = 0;
end
disp(['Rural weather file selected: ',climate_data])

%% Read EPW file
epwid = fopen(climate_data);

delimiterIn = ',';
headerlinesIn = 8;
C = importdata(climate_data, delimiterIn, headerlinesIn);

for i = 1:8
    header(i) = C.textdata(i,1);
end
fclose all;
epwid = fopen(climate_data);

i = 1;
while 1
    
    readin = fgetl(epwid);   
    if readin == -1 break;end
    % epwinput.values(i,:) = textscan(readin, '%f %f %f %f %f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %*[^\n]', 'delimiter',',','MultipleDelimsAsOne',1);
    epwinput.values(i,:) = textscan(readin, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %*[^\n]', 'delimiter',',','MultipleDelimsAsOne',1);
    i=i+1;
end

epwinput.values(1:8,:) = [];
fclose all;

%% Import Simulation Parameters
% XML import
try
    xml_location = strcat(CL_XML_PATH,'\',CL_XML);
    
catch
    [FileName,PathName] = uigetfile('.xml','Select Urban Parameter Input file');
    xml_location = strcat(PathName,FileName);
   
    fullyScripted = 0;
end
xml_input = xml_read(xml_location);
disp(['Urban Parameter file selected: ',xml_location]);

[pathstr,name,ext] = fileparts(xml_location);
xmlFileName = name;

% Loop around typologies 1 - 4
%typology = {'typology2','typology2', 'typology3', 'typology4'};

% Class housecleaning
%wall
%%
%typ1
if isa(xml_input.typology1.construction.wall.materials.names,'char')
    xml_input.typology1.construction.wall.materials.names = {xml_input.typology1.construction.wall.materials.names};
end
if isa(xml_input.typology1.construction.wall.materials.thermalConductivity,'double')
    xml_input.typology1.construction.wall.materials.thermalConductivity = {xml_input.typology1.construction.wall.materials.thermalConductivity};
end
if isa(xml_input.typology1.construction.wall.materials.volumetricHeatCapacity,'double')
    xml_input.typology1.construction.wall.materials.volumetricHeatCapacity = {xml_input.typology1.construction.wall.materials.volumetricHeatCapacity};
end

%roof
if isa(xml_input.typology1.construction.roof.materials.names,'char')
    xml_input.typology1.construction.roof.materials.names = {xml_input.typology1.construction.roof.materials.names};
end
if isa(xml_input.typology1.construction.roof.materials.thermalConductivity,'double')
    xml_input.typology1.construction.roof.materials.thermalConductivity = {xml_input.typology1.construction.roof.materials.thermalConductivity};
end
if isa(xml_input.typology1.construction.roof.materials.volumetricHeatCapacity,'double')
    xml_input.typology1.construction.roof.materials.volumetricHeatCapacity = {xml_input.typology1.construction.roof.materials.volumetricHeatCapacity};
end

%Mass
if isa(xml_input.typology1.construction.mass.materials.names,'char')
    xml_input.typology1.construction.mass.materials.names = {xml_input.typology1.construction.mass.materials.names};
end
if isa(xml_input.typology1.construction.mass.materials.thermalConductivity,'double')
    xml_input.typology1.construction.mass.materials.thermalConductivity = {xml_input.typology1.construction.mass.materials.thermalConductivity};
end
if isa(xml_input.typology1.construction.mass.materials.volumetricHeatCapacity,'double')
    xml_input.typology1.construction.mass.materials.volumetricHeatCapacity = {xml_input.typology1.construction.mass.materials.volumetricHeatCapacity};
end

%%
%typ2
%wall
if isa(xml_input.typology2.construction.wall.materials.names,'char')
    xml_input.typology2.construction.wall.materials.names = {xml_input.typology2.construction.wall.materials.names};
end
if isa(xml_input.typology2.construction.wall.materials.thermalConductivity,'double')
    xml_input.typology2.construction.wall.materials.thermalConductivity = {xml_input.typology2.construction.wall.materials.thermalConductivity};
end
if isa(xml_input.typology2.construction.wall.materials.volumetricHeatCapacity,'double')
    xml_input.typology2.construction.wall.materials.volumetricHeatCapacity = {xml_input.typology2.construction.wall.materials.volumetricHeatCapacity};
end

%roof
if isa(xml_input.typology1.construction.roof.materials.names,'char')
    xml_input.typology2.construction.roof.materials.names = {xml_input.typology2.construction.roof.materials.names};
end
if isa(xml_input.typology2.construction.roof.materials.thermalConductivity,'double')
    xml_input.typology2.construction.roof.materials.thermalConductivity = {xml_input.typology2.construction.roof.materials.thermalConductivity};
end
if isa(xml_input.typology2.construction.roof.materials.volumetricHeatCapacity,'double')
    xml_input.typology2.construction.roof.materials.volumetricHeatCapacity = {xml_input.typology2.construction.roof.materials.volumetricHeatCapacity};
end

%Mass
if isa(xml_input.typology2.construction.mass.materials.names,'char')
    xml_input.typology2.construction.mass.materials.names = {xml_input.typology2.construction.mass.materials.names};
end
if isa(xml_input.typology2.construction.mass.materials.thermalConductivity,'double')
    xml_input.typology2.construction.mass.materials.thermalConductivity = {xml_input.typology2.construction.mass.materials.thermalConductivity};
end
if isa(xml_input.typology2.construction.mass.materials.volumetricHeatCapacity,'double')
    xml_input.typology2.construction.mass.materials.volumetricHeatCapacity = {xml_input.typology2.construction.mass.materials.volumetricHeatCapacity};
end

%%
%typ3
if isa(xml_input.typology3.construction.wall.materials.names,'char')
    xml_input.typology{i}.construction.wall.materials.names = {xml_input.typology{i}.construction.wall.materials.names};
end
if isa(xml_input.typology3.construction.wall.materials.thermalConductivity,'double')
    xml_input.typology3.construction.wall.materials.thermalConductivity = {xml_input.typology3.construction.wall.materials.thermalConductivity};
end
if isa(xml_input.typology3.construction.wall.materials.volumetricHeatCapacity,'double')
    xml_input.typology3.construction.wall.materials.volumetricHeatCapacity = {xml_input.typology3.construction.wall.materials.volumetricHeatCapacity};
end

%roof
if isa(xml_input.typology3.construction.roof.materials.names,'char')
    xml_input.typology3.construction.roof.materials.names = {xml_input.typology3.construction.roof.materials.names};
end
if isa(xml_input.typology3.construction.roof.materials.thermalConductivity,'double')
    xml_input.typology3.construction.roof.materials.thermalConductivity = {xml_input.typology3.construction.roof.materials.thermalConductivity};
end
if isa(xml_input.typology3.construction.roof.materials.volumetricHeatCapacity,'double')
    xml_input.typology3.construction.roof.materials.volumetricHeatCapacity = {xml_input.typology3.construction.roof.materials.volumetricHeatCapacity};
end

%Mass
if isa(xml_input.typology3.construction.mass.materials.names,'char')
    xml_input.typology3.construction.mass.materials.names = {xml_input.typology3.construction.mass.materials.names};
end
if isa(xml_input.typology3.construction.mass.materials.thermalConductivity,'double')
    xml_input.typology3.construction.mass.materials.thermalConductivity = {xml_input.typology3.construction.mass.materials.thermalConductivity};
end
if isa(xml_input.typology3.construction.mass.materials.volumetricHeatCapacity,'double')
    xml_input.typology3.construction.mass.materials.volumetricHeatCapacity = {xml_input.typology3.construction.mass.materials.volumetricHeatCapacity};
end

%%
%typ4
if isa(xml_input.typology4.construction.wall.materials.names,'char')
    xml_input.typology{i}.construction.wall.materials.names = {xml_input.typology{i}.construction.wall.materials.names};
end
if isa(xml_input.typology4.construction.wall.materials.thermalConductivity,'double')
    xml_input.typology4.construction.wall.materials.thermalConductivity = {xml_input.typology4.construction.wall.materials.thermalConductivity};
end
if isa(xml_input.typology4.construction.wall.materials.volumetricHeatCapacity,'double')
    xml_input.typology4.construction.wall.materials.volumetricHeatCapacity = {xml_input.typology4.construction.wall.materials.volumetricHeatCapacity};
end

%roof
if isa(xml_input.typology4.construction.roof.materials.names,'char')
    xml_input.typology4.construction.roof.materials.names = {xml_input.typology4.construction.roof.materials.names};
end
if isa(xml_input.typology4.construction.roof.materials.thermalConductivity,'double')
    xml_input.typology4.construction.roof.materials.thermalConductivity = {xml_input.typology4.construction.roof.materials.thermalConductivity};
end
if isa(xml_input.typology4.construction.roof.materials.volumetricHeatCapacity,'double')
    xml_input.typology4.construction.roof.materials.volumetricHeatCapacity = {xml_input.typology4.construction.roof.materials.volumetricHeatCapacity};
end

%Mass
if isa(xml_input.typology4.construction.mass.materials.names,'char')
    xml_input.typology4.construction.mass.materials.names = {xml_input.typology4.construction.mass.materials.names};
end
if isa(xml_input.typology4.construction.mass.materials.thermalConductivity,'double')
    xml_input.typology4.construction.mass.materials.thermalConductivity = {xml_input.typology4.construction.mass.materials.thermalConductivity};
end
if isa(xml_input.typology4.construction.mass.materials.volumetricHeatCapacity,'double')
    xml_input.typology4.construction.mass.materials.volumetricHeatCapacity = {xml_input.typology4.construction.mass.materials.volumetricHeatCapacity};
end

% road
if isa(xml_input.urbanArea.urbanRoad.materials.names,'char')
    xml_input.urbanArea.urbanRoad.materials.names = {xml_input.urbanArea.urbanRoad.materials.names};
end
if isa(xml_input.urbanArea.urbanRoad.materials.thermalConductivity,'double')
    xml_input.urbanArea.urbanRoad.materials.thermalConductivity = {xml_input.urbanArea.urbanRoad.materials.thermalConductivity};
end
if isa(xml_input.urbanArea.urbanRoad.materials.volumetricHeatCapacity,'double')
    xml_input.urbanArea.urbanRoad.materials.volumetricHeatCapacity = {xml_input.urbanArea.urbanRoad.materials.volumetricHeatCapacity};
end

% Rural
if isa(xml_input.referenceSite.ruralRoad.materials.names,'char')
    xml_input.referenceSite.ruralRoad.materials.names = {xml_input.referenceSite.ruralRoad.materials.names};
end
if isa(xml_input.referenceSite.ruralRoad.materials.thermalConductivity,'double')
    xml_input.referenceSite.ruralRoad.materials.thermalConductivity = {xml_input.referenceSite.ruralRoad.materials.thermalConductivity};
end
if isa(xml_input.referenceSite.ruralRoad.materials.volumetricHeatCapacity,'double')
    xml_input.referenceSite.ruralRoad.materials.volumetricHeatCapacity = {xml_input.referenceSite.ruralRoad.materials.volumetricHeatCapacity};
end

%% Reconstruct xml Materials

% Wall
% Break wall into pieces
for i = 1:size(xml_input.typology1.construction.wall.materials.thickness,2)
    wallGrid_t{1,i} = xml_input.typology1.construction.wall.materials.thickness(i);
    wallGrid_n{1,i} = xml_input.typology1.construction.wall.materials.names(i);
    wallGrid_c{1,i} = xml_input.typology1.construction.wall.materials.thermalConductivity(i);
    wallGrid_h{1,i} = xml_input.typology1.construction.wall.materials.volumetricHeatCapacity(i);
end

for i = 1:size(wallGrid_t,2)
    clear Lr;
    Lr = wallGrid_t{i};
    clear Lm;
    if Lr <= 0.05;
                    Lm(1)=Lr;
    elseif Lr <= 0.1;
                    Lm(1)=Lr/2;
                    Lm(2)=Lr/2;
    elseif Lr <= 0.2;
                    Lm(1)=Lr/4;
                    Lm(2)=Lr/2;
                    Lm(3)=Lr/4;
    else Lr <= 0.3;
                    Lm(1)=Lr/6;
                    Lm(2)=Lr/3;
                    Lm(3)=Lr/3;
                    Lm(4)=Lr/6;
    end
    % Reassamble wall
    wallGrid_t{i} = Lm;
    clear names;
    clear thermalConductivity;
    clear volumetricHeatCapacity;
    names = wallGrid_n{1,i};
    thermalConductivity = wallGrid_c{1,i};
    volumetricHeatCapacity = wallGrid_h{1,i};

    for j = 1:size(Lm,2)
        wallGrid_n{j,i} = names;
        wallGrid_c{j,i} = thermalConductivity;
        wallGrid_h{j,i} = volumetricHeatCapacity;
    end
end
xml_input.typology1.construction.wall.materials.thickness = [];
for i = 1:size(wallGrid_t,2)
    xml_input.typology1.construction.wall.materials.thickness = [xml_input.typology1.construction.wall.materials.thickness wallGrid_t{i}];
end
xml_input.typology1.construction.wall.materials.names = [];
for j = 1:size(wallGrid_n,2)
    for i = 1:size(wallGrid_n,1)
        xml_input.typology1.construction.wall.materials.names = [xml_input.typology1.construction.wall.materials.names wallGrid_n{i,j}];
    end
end
xml_input.typology1.construction.wall.materials.thermalConductivity = [];
for j = 1:size(wallGrid_c,2)
    for i = 1:size(wallGrid_c,1)
        xml_input.typology1.construction.wall.materials.thermalConductivity = [xml_input.typology1.construction.wall.materials.thermalConductivity wallGrid_c{i,j}];
    end
end
xml_input.typology1.construction.wall.materials.volumetricHeatCapacity = [];
for j = 1:size(wallGrid_h,2)
    for i = 1:size(wallGrid_h,1)
        xml_input.typology1.construction.wall.materials.volumetricHeatCapacity = [xml_input.typology1.construction.wall.materials.volumetricHeatCapacity wallGrid_h{i,j}];
    end
end
xml_input.typology1.construction.wall.materials.thickness = xml_input.typology1.construction.wall.materials.thickness';

%%
%typ2
for i = 1:size(xml_input.typology2.construction.wall.materials.thickness,2)
    wallGrid_t2{1,i} = xml_input.typology2.construction.wall.materials.thickness(i);
    wallGrid_n2{1,i} = xml_input.typology2.construction.wall.materials.names(i);
    wallGrid_c2{1,i} = xml_input.typology2.construction.wall.materials.thermalConductivity(i);
    wallGrid_h2{1,i} = xml_input.typology2.construction.wall.materials.volumetricHeatCapacity(i);
end
for i = 1:size(wallGrid_t2,2)
    clear Lr;
    Lr = wallGrid_t2{i};
    clear Lm;
    if Lr <= 0.05;
                        Lm(1)=Lr;
        elseif Lr <= 0.1;
                        Lm(1)=Lr/2;
                        Lm(2)=Lr/2;
        elseif Lr <= 0.2;
                        Lm(1)=Lr/4;
                        Lm(2)=Lr/2;
                        Lm(3)=Lr/4;
        else Lr <= 0.3;
                        Lm(1)=Lr/6;
                        Lm(2)=Lr/3;
                        Lm(3)=Lr/3;
                        Lm(4)=Lr/6;
    end
    % Reassamble wall
    wallGrid_t2{i} = Lm;
    
    clear names;
    clear thermalConductivity;
    clear volumetricHeatCapacity;
    names = wallGrid_n2{1,i};
    thermalConductivity = wallGrid_c2{1,i};
    volumetricHeatCapacity = wallGrid_h2{1,i};

    for j = 1:size(Lm,2)
        wallGrid_n2{j,i} = names;
        wallGrid_c2{j,i} = thermalConductivity;
        wallGrid_h2{j,i} = volumetricHeatCapacity;
    end
end

xml_input.typology2.construction.wall.materials.thickness = [];
for i = 1:size(wallGrid_t2,2)
    xml_input.typology2.construction.wall.materials.thickness = [xml_input.typology2.construction.wall.materials.thickness wallGrid_t2{i}];
end
xml_input.typology2.construction.wall.materials.names = [];
for j = 1:size(wallGrid_n2,2)
    for i = 1:size(wallGrid_n2,1)
        xml_input.typology2.construction.wall.materials.names = [xml_input.typology2.construction.wall.materials.names wallGrid_n2{i,j}];
    end
end
xml_input.typology2.construction.wall.materials.thermalConductivity = [];
for j = 1:size(wallGrid_c2,2)
    for i = 1:size(wallGrid_c2,1)
        xml_input.typology2.construction.wall.materials.thermalConductivity = [xml_input.typology2.construction.wall.materials.thermalConductivity wallGrid_c2{i,j}];
    end
end
xml_input.typology2.construction.wall.materials.volumetricHeatCapacity = [];
for j = 1:size(wallGrid_h2,2)
    for i = 1:size(wallGrid_h2,1)
        xml_input.typology2.construction.wall.materials.volumetricHeatCapacity = [xml_input.typology2.construction.wall.materials.volumetricHeatCapacity wallGrid_h2{i,j}];
    end
end
xml_input.typology2.construction.wall.materials.thickness = xml_input.typology2.construction.wall.materials.thickness';
%%
%typ3
for i = 1:size(xml_input.typology3.construction.wall.materials.thickness,3)
    wallGrid_t3{1,i} = xml_input.typology3.construction.wall.materials.thickness(i);
    wallGrid_n3{1,i} = xml_input.typology3.construction.wall.materials.names(i);
    wallGrid_c3{1,i} = xml_input.typology3.construction.wall.materials.thermalConductivity(i);
    wallGrid_h3{1,i} = xml_input.typology3.construction.wall.materials.volumetricHeatCapacity(i);
end
for i = 1:size(wallGrid_t3,2)
    clear Lr;
    Lr = wallGrid_t3{i};
    clear Lm;
    if Lr <= 0.05;
                        Lm(1)=Lr;
        elseif Lr <= 0.1;
                        Lm(1)=Lr/2;
                        Lm(2)=Lr/2;
        elseif Lr <= 0.2;
                        Lm(1)=Lr/4;
                        Lm(2)=Lr/2;
                        Lm(3)=Lr/4;
        else Lr <= 0.3;
                        Lm(1)=Lr/6;
                        Lm(2)=Lr/3;
                        Lm(3)=Lr/3;
                        Lm(4)=Lr/6;
    end
    % Reassamble wall
    wallGrid_t3{i} = Lm;
    
    clear names;
    clear thermalConductivity;
    clear volumetricHeatCapacity;
    names = wallGrid_n3{1,i};
    thermalConductivity = wallGrid_c3{1,i};
    volumetricHeatCapacity = wallGrid_h3{1,i};

    for j = 1:size(Lm,2)
        wallGrid_n3{j,i} = names;
        wallGrid_c3{j,i} = thermalConductivity;
        wallGrid_h3{j,i} = volumetricHeatCapacity;
    end
end

xml_input.typology3.construction.wall.materials.thickness = [];
for i = 1:size(wallGrid_t3,2)
    xml_input.typology3.construction.wall.materials.thickness = [xml_input.typology3.construction.wall.materials.thickness wallGrid_t3{i}];
end
xml_input.typology3.construction.wall.materials.names = [];
for j = 1:size(wallGrid_n3,2)
    for i = 1:size(wallGrid_n3,1)
        xml_input.typology3.construction.wall.materials.names = [xml_input.typology3.construction.wall.materials.names wallGrid_n3{i,j}];
    end
end
xml_input.typology3.construction.wall.materials.thermalConductivity = [];
for j = 1:size(wallGrid_c3,2)
    for i = 1:size(wallGrid_c3,1)
        xml_input.typology3.construction.wall.materials.thermalConductivity = [xml_input.typology3.construction.wall.materials.thermalConductivity wallGrid_c3{i,j}];
    end
end
xml_input.typology3.construction.wall.materials.volumetricHeatCapacity = [];
for j = 1:size(wallGrid_h3,2)
    for i = 1:size(wallGrid_h3,1)
        xml_input.typology3.construction.wall.materials.volumetricHeatCapacity = [xml_input.typology3.construction.wall.materials.volumetricHeatCapacity wallGrid_h3{i,j}];
    end
end
xml_input.typology3.construction.wall.materials.thickness = xml_input.typology3.construction.wall.materials.thickness';
%%
%typ4
for i = 1:size(xml_input.typology4.construction.wall.materials.thickness,2)
    wallGrid_t4{1,i} = xml_input.typology4.construction.wall.materials.thickness(i);
    wallGrid_n4{1,i} = xml_input.typology4.construction.wall.materials.names(i);
    wallGrid_c4{1,i} = xml_input.typology4.construction.wall.materials.thermalConductivity(i);
    wallGrid_h4{1,i} = xml_input.typology4.construction.wall.materials.volumetricHeatCapacity(i);
end
for i = 1:size(wallGrid_t4,2)
    clear Lr;
    Lr = wallGrid_t4{i};
    clear Lm;
    if Lr <= 0.05;
                        Lm(1)=Lr;
        elseif Lr <= 0.1;
                        Lm(1)=Lr/2;
                        Lm(2)=Lr/2;
        elseif Lr <= 0.2;
                        Lm(1)=Lr/4;
                        Lm(2)=Lr/2;
                        Lm(3)=Lr/4;
        else Lr <= 0.3;
                        Lm(1)=Lr/6;
                        Lm(2)=Lr/3;
                        Lm(3)=Lr/3;
                        Lm(4)=Lr/6;
    end
    % Reassamble wall
    wallGrid_t4{i} = Lm;
    
    clear names;
    clear thermalConductivity;
    clear volumetricHeatCapacity;
    names = wallGrid_n4{1,i};
    thermalConductivity = wallGrid_c4{1,i};
    volumetricHeatCapacity = wallGrid_h4{1,i};

    for j = 1:size(Lm,2)
        wallGrid_n4{j,i} = names;
        wallGrid_c4{j,i} = thermalConductivity;
        wallGrid_h4{j,i} = volumetricHeatCapacity;
    end
end

xml_input.typology4.construction.wall.materials.thickness = [];
for i = 1:size(wallGrid_t4,2)
    xml_input.typology4.construction.wall.materials.thickness = [xml_input.typology4.construction.wall.materials.thickness wallGrid_t4{i}];
end
xml_input.typology4.construction.wall.materials.names = [];
for j = 1:size(wallGrid_n4,2)
    for i = 1:size(wallGrid_n4,1)
        xml_input.typology4.construction.wall.materials.names = [xml_input.typology4.construction.wall.materials.names wallGrid_n4{i,j}];
    end
end
xml_input.typology4.construction.wall.materials.thermalConductivity = [];
for j = 1:size(wallGrid_c4,2)
    for i = 1:size(wallGrid_c4,1)
        xml_input.typology4.construction.wall.materials.thermalConductivity = [xml_input.typology4.construction.wall.materials.thermalConductivity wallGrid_c4{i,j}];
    end
end
xml_input.typology4.construction.wall.materials.volumetricHeatCapacity = [];
for j = 1:size(wallGrid_h4,2)
    for i = 1:size(wallGrid_h4,1)
        xml_input.typology4.construction.wall.materials.volumetricHeatCapacity = [xml_input.typology4.construction.wall.materials.volumetricHeatCapacity wallGrid_h4{i,j}];
    end
end
xml_input.typology4.construction.wall.materials.thickness = xml_input.typology4.construction.wall.materials.thickness';
%% Roof
for i = 1:size(xml_input.typology1.construction.roof.materials.thickness,2)
    roofGrid_t{1,i} = xml_input.typology1.construction.roof.materials.thickness(i);
    roofGrid_n{1,i} = xml_input.typology1.construction.roof.materials.names(i);
    roofGrid_c{1,i} = xml_input.typology1.construction.roof.materials.thermalConductivity(i);
    roofGrid_h{1,i} = xml_input.typology1.construction.roof.materials.volumetricHeatCapacity(i);
end
for i = 1:size(roofGrid_t,2)
    clear Lr;
    Lr = roofGrid_t{i};
    clear Lm;
    if Lr <= 0.05;
                    Lm(1)=Lr;
    elseif Lr <= 0.1;
                    Lm(1)=Lr/2;
                    Lm(2)=Lr/2;
    elseif Lr <= 0.2;
                    Lm(1)=Lr/4;
                    Lm(2)=Lr/2;
                    Lm(3)=Lr/4;
    else Lr <= 0.3;
                    Lm(1)=Lr/6;
                    Lm(2)=Lr/3;
                    Lm(3)=Lr/3;
                    Lm(4)=Lr/6;
    end

    %Reassemble roof
    roofGrid_t{i} = Lm;
    clear names;
    clear thermalConductivity;
    clear volumetricHeatCapacity;
    names = roofGrid_n{1,i};
    thermalConductivity = roofGrid_c{1,i};
    volumetricHeatCapacity = roofGrid_h{1,i};

    for j = 1:size(Lm,2)
        roofGrid_n{j,i} = names;
        roofGrid_c{j,i} = thermalConductivity;
        roofGrid_h{j,i} = volumetricHeatCapacity;
    end
end

xml_input.typology1.construction.roof.materials.thickness = [];
for i = 1:size(roofGrid_t,2)
    xml_input.typology1.construction.roof.materials.thickness = [xml_input.typology1.construction.roof.materials.thickness roofGrid_t{i}];
end
xml_input.typology1.construction.roof.materials.names = [];
for j = 1:size(roofGrid_n,2)
    for i = 1:size(roofGrid_n,1)
        xml_input.typology1.construction.roof.materials.names = [xml_input.typology1.construction.roof.materials.names roofGrid_n{i,j}];
    end
end
xml_input.typology1.construction.roof.materials.thermalConductivity = [];
for j = 1:size(roofGrid_c,2)
    for i = 1:size(roofGrid_c,1)
        xml_input.typology1.construction.roof.materials.thermalConductivity = [xml_input.typology1.construction.roof.materials.thermalConductivity roofGrid_c{i,j}];
    end
end
xml_input.typology1.construction.roof.materials.volumetricHeatCapacity = [];
for j = 1:size(roofGrid_h,2)
    for i = 1:size(roofGrid_h,1)
        xml_input.typology1.construction.roof.materials.volumetricHeatCapacity = [xml_input.typology1.construction.roof.materials.volumetricHeatCapacity roofGrid_h{i,j}];
    end
end
xml_input.typology1.construction.roof.materials.thickness = xml_input.typology1.construction.roof.materials.thickness';

%%


for i = 1:size(xml_input.typology2.construction.roof.materials.thickness,2)
    roofGrid_t2{1,i} = xml_input.typology2.construction.roof.materials.thickness(i);
    roofGrid_n2{1,i} = xml_input.typology2.construction.roof.materials.names(i);
    roofGrid_c2{1,i} = xml_input.typology2.construction.roof.materials.thermalConductivity(i);
    roofGrid_h2{1,i} = xml_input.typology2.construction.roof.materials.volumetricHeatCapacity(i);
end
for i = 1:size(roofGrid_t2,2)
    clear Lr;
    Lr = roofGrid_t2{i};
    clear Lm;
    if Lr <= 0.05;
                    Lm(1)=Lr;
    elseif Lr <= 0.1;
                    Lm(1)=Lr/2;
                    Lm(2)=Lr/2;
    elseif Lr <= 0.2;
                    Lm(1)=Lr/4;
                    Lm(2)=Lr/2;
                    Lm(3)=Lr/4;
    else Lr <= 0.3;
                    Lm(1)=Lr/6;
                    Lm(2)=Lr/3;
                    Lm(3)=Lr/3;
                    Lm(4)=Lr/6;
    end

    % Reassemble roof
    roofGrid_t2{i} = Lm;
    clear names;
    clear thermalConductivity;
    clear volumetricHeatCapacity;
    names = roofGrid_n2{1,i};
    thermalConductivity = roofGrid_c2{1,i};
    volumetricHeatCapacity = roofGrid_h2{1,i};

    for j = 1:size(Lm,2)
        roofGrid_n2{j,i} = names;
        roofGrid_c2{j,i} = thermalConductivity;
        roofGrid_h2{j,i} = volumetricHeatCapacity;
    end
end

xml_input.typology2.construction.roof.materials.thickness = [];
for i = 1:size(roofGrid_t2,2)
    xml_input.typology2.construction.roof.materials.thickness = [xml_input.typology2.construction.roof.materials.thickness roofGrid_t2{i}];
end
xml_input.typology2.construction.roof.materials.names = [];
for j = 1:size(roofGrid_n2,2)
    for i = 1:size(roofGrid_n2,1)
        xml_input.typology2.construction.roof.materials.names = [xml_input.typology2.construction.roof.materials.names roofGrid_n2{i,j}];
    end
end
xml_input.typology2.construction.roof.materials.thermalConductivity = [];
for j = 1:size(roofGrid_c2,2)
    for i = 1:size(roofGrid_c2,1)
        xml_input.typology2.construction.roof.materials.thermalConductivity = [xml_input.typology2.construction.roof.materials.thermalConductivity roofGrid_c2{i,j}];
    end
end
xml_input.typology2.construction.roof.materials.volumetricHeatCapacity = [];
for j = 1:size(roofGrid_h2,2)
    for i = 1:size(roofGrid_h2,1)
        xml_input.typology2.construction.roof.materials.volumetricHeatCapacity = [xml_input.typology2.construction.roof.materials.volumetricHeatCapacity roofGrid_h2{i,j}];
    end
end
xml_input.typology2.construction.roof.materials.thickness = xml_input.typology2.construction.roof.materials.thickness';
%%
%typ3
for i = 1:size(xml_input.typology3.construction.roof.materials.thickness,2)
    roofGrid_t3{1,i} = xml_input.typology3.construction.roof.materials.thickness(i);
    roofGrid_n3{1,i} = xml_input.typology3.construction.roof.materials.names(i);
    roofGrid_c3{1,i} = xml_input.typology3.construction.roof.materials.thermalConductivity(i);
    roofGrid_h3{1,i} = xml_input.typology3.construction.roof.materials.volumetricHeatCapacity(i);
end
for i = 1:size(roofGrid_t3,2)
    clear Lr;
    Lr = roofGrid_t3{i};
    clear Lm;
    if Lr <= 0.05;
                    Lm(1)=Lr;
    elseif Lr <= 0.1;
                    Lm(1)=Lr/2;
                    Lm(2)=Lr/2;
    elseif Lr <= 0.2;
                    Lm(1)=Lr/4;
                    Lm(2)=Lr/2;
                    Lm(3)=Lr/4;
    else Lr <= 0.3;
                    Lm(1)=Lr/6;
                    Lm(2)=Lr/3;
                    Lm(3)=Lr/3;
                    Lm(4)=Lr/6;
    end

    %Reassemble roof
    roofGrid_t3{i} = Lm;
    clear names;
    clear thermalConductivity;
    clear volumetricHeatCapacity;
    names = roofGrid_n3{1,i};
    thermalConductivity = roofGrid_c3{1,i};
    volumetricHeatCapacity = roofGrid_h3{1,i};

    for j = 1:size(Lm,2)
        roofGrid_n3{j,i} = names;
        roofGrid_c3{j,i} = thermalConductivity;
        roofGrid_h3{j,i} = volumetricHeatCapacity;
    end
end

xml_input.typology3.construction.roof.materials.thickness = [];
for i = 1:size(roofGrid_t3,2)
    xml_input.typology3.construction.roof.materials.thickness = [xml_input.typology3.construction.roof.materials.thickness roofGrid_t3{i}];
end
xml_input.typology3.construction.roof.materials.names = [];
for j = 1:size(roofGrid_n3,2)
    for i = 1:size(roofGrid_n3,1)
        xml_input.typology3.construction.roof.materials.names = [xml_input.typology3.construction.roof.materials.names roofGrid_n3{i,j}];
    end
end
xml_input.typology3.construction.roof.materials.thermalConductivity = [];
for j = 1:size(roofGrid_c3,2)
    for i = 1:size(roofGrid_c3,1)
        xml_input.typology3.construction.roof.materials.thermalConductivity = [xml_input.typology3.construction.roof.materials.thermalConductivity roofGrid_c3{i,j}];
    end
end
xml_input.typology3.construction.roof.materials.volumetricHeatCapacity = [];
for j = 1:size(roofGrid_h3,2)
    for i = 1:size(roofGrid_h3,1)
        xml_input.typology3.construction.roof.materials.volumetricHeatCapacity = [xml_input.typology3.construction.roof.materials.volumetricHeatCapacity roofGrid_h3{i,j}];
    end
end
xml_input.typology3.construction.roof.materials.thickness = xml_input.typology3.construction.roof.materials.thickness';

%%
%typ4
for i = 1:size(xml_input.typology4.construction.roof.materials.thickness,2)
    roofGrid_t4{1,i} = xml_input.typology4.construction.roof.materials.thickness(i);
    roofGrid_n4{1,i} = xml_input.typology4.construction.roof.materials.names(i);
    roofGrid_c4{1,i} = xml_input.typology4.construction.roof.materials.thermalConductivity(i);
    roofGrid_h4{1,i} = xml_input.typology4.construction.roof.materials.volumetricHeatCapacity(i);
end
for i = 1:size(roofGrid_t4,2)
    clear Lr;
    Lr = roofGrid_t4{i};
    clear Lm;
    if Lr <= 0.05;
                    Lm(1)=Lr;
    elseif Lr <= 0.1;
                    Lm(1)=Lr/2;
                    Lm(2)=Lr/2;
    elseif Lr <= 0.2;
                    Lm(1)=Lr/4;
                    Lm(2)=Lr/2;
                    Lm(3)=Lr/4;
    else Lr <= 0.3;
                    Lm(1)=Lr/6;
                    Lm(2)=Lr/3;
                    Lm(3)=Lr/3;
                    Lm(4)=Lr/6;
    end

    %Reassemble roof
    roofGrid_t4{i} = Lm;
    clear names;
    clear thermalConductivity;
    clear volumetricHeatCapacity;
    names = roofGrid_n4{1,i};
    thermalConductivity = roofGrid_c4{1,i};
    volumetricHeatCapacity = roofGrid_h4{1,i};

    for j = 1:size(Lm,2)
        roofGrid_n4{j,i} = names;
        roofGrid_c4{j,i} = thermalConductivity;
        roofGrid_h4{j,i} = volumetricHeatCapacity;
    end
end

xml_input.typology4.construction.roof.materials.thickness = [];
for i = 1:size(roofGrid_t4,2)
    xml_input.typology4.construction.roof.materials.thickness = [xml_input.typology4.construction.roof.materials.thickness roofGrid_t4{i}];
end
xml_input.typology4.construction.roof.materials.names = [];
for j = 1:size(roofGrid_n4,2)
    for i = 1:size(roofGrid_n4,1)
        xml_input.typology4.construction.roof.materials.names = [xml_input.typology4.construction.roof.materials.names roofGrid_n4{i,j}];
    end
end
xml_input.typology4.construction.roof.materials.thermalConductivity = [];
for j = 1:size(roofGrid_c4,2)
    for i = 1:size(roofGrid_c4,1)
        xml_input.typology4.construction.roof.materials.thermalConductivity = [xml_input.typology4.construction.roof.materials.thermalConductivity roofGrid_c4{i,j}];
    end
end
xml_input.typology4.construction.roof.materials.volumetricHeatCapacity = [];
for j = 1:size(roofGrid_h4,2)
    for i = 1:size(roofGrid_h4,1)
        xml_input.typology4.construction.roof.materials.volumetricHeatCapacity = [xml_input.typology4.construction.roof.materials.volumetricHeatCapacity roofGrid_h4{i,j}];
    end
end
xml_input.typology4.construction.roof.materials.thickness = xml_input.typology4.construction.roof.materials.thickness';

%% Mass
%typ1
% Break mass into pieces
for i902 = 1:size(xml_input.typology1.construction.mass.materials.thickness,2)
	massGrid_t{1,i902} = xml_input.typology1.construction.mass.materials.thickness(i902);
	massGrid_n{1,i902} = xml_input.typology1.construction.mass.materials.names(i902);
	massGrid_c{1,i902} = xml_input.typology1.construction.mass.materials.thermalConductivity(i902);
	massGrid_h{1,i902} = xml_input.typology1.construction.mass.materials.volumetricHeatCapacity(i902);
end
for i903 = 1:size(massGrid_t,2)
clear Lr;
Lr = massGrid_t{i903};
clear Lm;
if Lr <= 0.05;
				Lm(1)=Lr;
elseif Lr <= 0.1;
				Lm(1)=Lr/2;
				Lm(2)=Lr/2;
elseif Lr <= 0.2;
				Lm(1)=Lr/4;
				Lm(2)=Lr/2;
				Lm(3)=Lr/4;
else Lr <= 0.3;
				Lm(1)=Lr/6;
				Lm(2)=Lr/3;
				Lm(3)=Lr/3;
				Lm(4)=Lr/6;
end

% Reassamble mass
massGrid_t{i903} = Lm;
clear names;
clear thermalConductivity;
clear volumetricHeatCapacity;
names = massGrid_n{1,i903};
thermalConductivity = massGrid_c{1,i903};
volumetricHeatCapacity = massGrid_h{1,i903};

for j902 = 1:size(Lm,2)
	massGrid_n{j902,i903} = names;
	massGrid_c{j902,i903} = thermalConductivity;
	massGrid_h{j902,i903} = volumetricHeatCapacity;
end
end
xml_input.typology1.construction.mass.materials.thickness = [];
for i904 = 1:size(massGrid_t,2)
	xml_input.typology1.construction.mass.materials.thickness = [xml_input.typology1.construction.mass.materials.thickness massGrid_t{i904}];
end
xml_input.typology1.construction.mass.materials.names = [];
for j905 = 1:size(massGrid_n,2)
	for i905 = 1:size(massGrid_n,1)
		xml_input.typology1.construction.mass.materials.names = [xml_input.typology1.construction.mass.materials.names massGrid_n{i905,j905}];
	end
end
xml_input.typology1.construction.mass.materials.thermalConductivity = [];
for j906 = 1:size(massGrid_c,2)
	for i906 = 1:size(massGrid_c,1)
		xml_input.typology1.construction.mass.materials.thermalConductivity = [xml_input.typology1.construction.mass.materials.thermalConductivity massGrid_c{i906,j906}];
	end
end
xml_input.typology1.construction.mass.materials.volumetricHeatCapacity = [];
for j907 = 1:size(massGrid_h,2)
	for i907 = 1:size(massGrid_h,1)
		xml_input.typology1.construction.mass.materials.volumetricHeatCapacity = [xml_input.typology1.construction.mass.materials.volumetricHeatCapacity massGrid_h{i907,j907}];
	end
end
xml_input.typology1.construction.mass.materials.thickness = xml_input.typology1.construction.mass.materials.thickness';

%%typ2
% Break mass into pieces
for i902 = 1:size(xml_input.typology2.construction.mass.materials.thickness,2)
	massGrid_t2{1,i902} = xml_input.typology2.construction.mass.materials.thickness(i902);
	massGrid_n2{1,i902} = xml_input.typology2.construction.mass.materials.names(i902);
	massGrid_c2{1,i902} = xml_input.typology2.construction.mass.materials.thermalConductivity(i902);
	massGrid_h2{1,i902} = xml_input.typology2.construction.mass.materials.volumetricHeatCapacity(i902);
end
for i903 = 1:size(massGrid_t2,2)
clear Lr;
Lr = massGrid_t2{i903};
clear Lm;
if Lr <= 0.05;
				Lm(1)=Lr;
elseif Lr <= 0.1;
				Lm(1)=Lr/2;
				Lm(2)=Lr/2;
elseif Lr <= 0.2;
				Lm(1)=Lr/4;
				Lm(2)=Lr/2;
				Lm(3)=Lr/4;
else Lr <= 0.3;
				Lm(1)=Lr/6;
				Lm(2)=Lr/3;
				Lm(3)=Lr/3;
				Lm(4)=Lr/6;
end

%Reassemble mass
massGrid_t2{i903} = Lm;
clear names;
clear thermalConductivity;
clear volumetricHeatCapacity;
names = massGrid_n2{1,i903};
thermalConductivity = massGrid_c2{1,i903};
volumetricHeatCapacity = massGrid_h2{1,i903};

for j902 = 1:size(Lm,2)
	massGrid_n2{j902,i903} = names;
	massGrid_c2{j902,i903} = thermalConductivity;
	massGrid_h2{j902,i903} = volumetricHeatCapacity;
end
end
xml_input.typology2.construction.mass.materials.thickness = [];
for i904 = 1:size(massGrid_t2,2)
	xml_input.typology2.construction.mass.materials.thickness = [xml_input.typology2.construction.mass.materials.thickness massGrid_t2{i904}];
end
xml_input.typology2.construction.mass.materials.names = [];
for j905 = 1:size(massGrid_n2,2)
	for i905 = 1:size(massGrid_n2,1)
		xml_input.typology2.construction.mass.materials.names = [xml_input.typology2.construction.mass.materials.names massGrid_n2{i905,j905}];
	end
end
xml_input.typology2.construction.mass.materials.thermalConductivity = [];
for j906 = 1:size(massGrid_c2,2)
	for i906 = 1:size(massGrid_c2,1)
		xml_input.typology2.construction.mass.materials.thermalConductivity = [xml_input.typology2.construction.mass.materials.thermalConductivity massGrid_c2{i906,j906}];
	end
end
xml_input.typology2.construction.mass.materials.volumetricHeatCapacity = [];
for j907 = 1:size(massGrid_h2,2)
	for i907 = 1:size(massGrid_h2,1)
		xml_input.typology2.construction.mass.materials.volumetricHeatCapacity = [xml_input.typology2.construction.mass.materials.volumetricHeatCapacity massGrid_h2{i907,j907}];
	end
end
xml_input.typology2.construction.mass.materials.thickness = xml_input.typology2.construction.mass.materials.thickness';

%%typ3
% Break mass into pieces
for i902 = 1:size(xml_input.typology3.construction.mass.materials.thickness,2)
	massGrid_t3{1,i902} = xml_input.typology3.construction.mass.materials.thickness(i902);
	massGrid_n3{1,i902} = xml_input.typology3.construction.mass.materials.names(i902);
	massGrid_c3{1,i902} = xml_input.typology3.construction.mass.materials.thermalConductivity(i902);
	massGrid_h3{1,i902} = xml_input.typology3.construction.mass.materials.volumetricHeatCapacity(i902);
end
for i903 = 1:size(massGrid_t3,2)
clear Lr;
Lr = massGrid_t3{i903};
clear Lm;
if Lr <= 0.05;
				Lm(1)=Lr;
elseif Lr <= 0.1;
				Lm(1)=Lr/2;
				Lm(2)=Lr/2;
elseif Lr <= 0.2;
				Lm(1)=Lr/4;
				Lm(2)=Lr/2;
				Lm(3)=Lr/4;
else Lr <= 0.3;
				Lm(1)=Lr/6;
				Lm(2)=Lr/3;
				Lm(3)=Lr/3;
				Lm(4)=Lr/6;
end

%Reassemble mass
massGrid_t3{i903} = Lm;
clear names;
clear thermalConductivity;
clear volumetricHeatCapacity;
names = massGrid_n3{1,i903};
thermalConductivity = massGrid_c3{1,i903};
volumetricHeatCapacity = massGrid_h3{1,i903};

for j902 = 1:size(Lm,2)
	massGrid_n3{j902,i903} = names;
	massGrid_c3{j902,i903} = thermalConductivity;
	massGrid_h3{j902,i903} = volumetricHeatCapacity;
end
end
xml_input.typology3.construction.mass.materials.thickness = [];
for i904 = 1:size(massGrid_t3,2)
	xml_input.typology3.construction.mass.materials.thickness = [xml_input.typology3.construction.mass.materials.thickness massGrid_t3{i904}];
end
xml_input.typology3.construction.mass.materials.names = [];
for j905 = 1:size(massGrid_n3,2)
	for i905 = 1:size(massGrid_n3,1)
		xml_input.typology3.construction.mass.materials.names = [xml_input.typology3.construction.mass.materials.names massGrid_n3{i905,j905}];
	end
end
xml_input.typology3.construction.mass.materials.thermalConductivity = [];
for j906 = 1:size(massGrid_c3,2)
	for i906 = 1:size(massGrid_c3,1)
		xml_input.typology3.construction.mass.materials.thermalConductivity = [xml_input.typology3.construction.mass.materials.thermalConductivity massGrid_c3{i906,j906}];
	end
end
xml_input.typology3.construction.mass.materials.volumetricHeatCapacity = [];
for j907 = 1:size(massGrid_h3,2)
	for i907 = 1:size(massGrid_h3,1)
		xml_input.typology3.construction.mass.materials.volumetricHeatCapacity = [xml_input.typology3.construction.mass.materials.volumetricHeatCapacity massGrid_h3{i907,j907}];
	end
end
xml_input.typology3.construction.mass.materials.thickness = xml_input.typology3.construction.mass.materials.thickness';

%%typ4
% Break mass into pieces
for i902 = 1:size(xml_input.typology4.construction.mass.materials.thickness,2)
	massGrid_t4{1,i902} = xml_input.typology4.construction.mass.materials.thickness(i902);
	massGrid_n4{1,i902} = xml_input.typology4.construction.mass.materials.names(i902);
	massGrid_c4{1,i902} = xml_input.typology4.construction.mass.materials.thermalConductivity(i902);
	massGrid_h4{1,i902} = xml_input.typology4.construction.mass.materials.volumetricHeatCapacity(i902);
end
for i903 = 1:size(massGrid_t4,2)
clear Lr;
Lr = massGrid_t4{i903};
clear Lm;
if Lr <= 0.05;
				Lm(1)=Lr;
elseif Lr <= 0.1;
				Lm(1)=Lr/2;
				Lm(2)=Lr/2;
elseif Lr <= 0.2;
				Lm(1)=Lr/4;
				Lm(2)=Lr/2;
				Lm(3)=Lr/4;
else Lr <= 0.3;
				Lm(1)=Lr/6;
				Lm(2)=Lr/3;
				Lm(3)=Lr/3;
				Lm(4)=Lr/6;
end

%Reassemble mass
massGrid_t4{i903} = Lm;
clear names;
clear thermalConductivity;
clear volumetricHeatCapacity;
names = massGrid_n4{1,i903};
thermalConductivity = massGrid_c4{1,i903};
volumetricHeatCapacity = massGrid_h4{1,i903};

for j902 = 1:size(Lm,2)
	massGrid_n4{j902,i903} = names;
	massGrid_c4{j902,i903} = thermalConductivity;
	massGrid_h4{j902,i903} = volumetricHeatCapacity;
end
end
xml_input.typology4.construction.mass.materials.thickness = [];
for i904 = 1:size(massGrid_t4,2)
	xml_input.typology4.construction.mass.materials.thickness = [xml_input.typology4.construction.mass.materials.thickness massGrid_t4{i904}];
end
xml_input.typology4.construction.mass.materials.names = [];
for j905 = 1:size(massGrid_n4,2)
	for i905 = 1:size(massGrid_n4,1)
		xml_input.typology4.construction.mass.materials.names = [xml_input.typology4.construction.mass.materials.names massGrid_n4{i905,j905}];
	end
end
xml_input.typology4.construction.mass.materials.thermalConductivity = [];
for j906 = 1:size(massGrid_c4,2)
	for i906 = 1:size(massGrid_c4,1)
		xml_input.typology4.construction.mass.materials.thermalConductivity = [xml_input.typology4.construction.mass.materials.thermalConductivity massGrid_c4{i906,j906}];
	end
end
xml_input.typology4.construction.mass.materials.volumetricHeatCapacity = [];
for j907 = 1:size(massGrid_h4,2)
	for i907 = 1:size(massGrid_h4,1)
		xml_input.typology4.construction.mass.materials.volumetricHeatCapacity = [xml_input.typology4.construction.mass.materials.volumetricHeatCapacity massGrid_h4{i907,j907}];
	end
end
xml_input.typology4.construction.mass.materials.thickness = xml_input.typology4.construction.mass.materials.thickness';

%% Road
% Break road into pieces
for i902 = 1:size(xml_input.urbanArea.urbanRoad.materials.thickness,2)
	roadGrid_t{1,i902} = xml_input.urbanArea.urbanRoad.materials.thickness(i902);
	roadGrid_n{1,i902} = xml_input.urbanArea.urbanRoad.materials.names(i902);
	roadGrid_c{1,i902} = xml_input.urbanArea.urbanRoad.materials.thermalConductivity(i902);
	roadGrid_h{1,i902} = xml_input.urbanArea.urbanRoad.materials.volumetricHeatCapacity(i902);
end
for i903 = 1:size(roadGrid_t,2)
	clear Lr;
	Lr = roadGrid_t{i903};
	clear Lm;
	if Lr <= 0.05;
					Lm(1)=Lr;
	elseif Lr <= 0.1;
					Lm(1)=Lr/2;
					Lm(2)=Lr/2;
	elseif Lr <= 0.2;
					Lm(1)=Lr/4;
					Lm(2)=Lr/2;
					Lm(3)=Lr/4;
	else Lr <= 0.3;
					Lm(1)=Lr/6;
					Lm(2)=Lr/3;
					Lm(3)=Lr/3;
					Lm(4)=Lr/6;
	end

	% Reassemble road
	roadGrid_t{i903} = Lm;
	clear names;
	clear thermalConductivity;
	clear volumetricHeatCapacity;
	names = roadGrid_n{1,i903};
	thermalConductivity = roadGrid_c{1,i903};
	volumetricHeatCapacity = roadGrid_h{1,i903};

	for j902 = 1:size(Lm,2)
		roadGrid_n{j902,i903} = names;
		roadGrid_c{j902,i903} = thermalConductivity;
		roadGrid_h{j902,i903} = volumetricHeatCapacity;
	end
end
xml_input.urbanArea.urbanRoad.materials.thickness = [];
for i904 = 1:size(roadGrid_t,2)
	xml_input.urbanArea.urbanRoad.materials.thickness = [xml_input.urbanArea.urbanRoad.materials.thickness roadGrid_t{i904}];
end
xml_input.urbanArea.urbanRoad.materials.names = [];
for j905 = 1:size(roadGrid_n,2)
	for i905 = 1:size(roadGrid_n,1)
		xml_input.urbanArea.urbanRoad.materials.names = [xml_input.urbanArea.urbanRoad.materials.names roadGrid_n{i905,j905}];
	end
end
xml_input.urbanArea.urbanRoad.materials.thermalConductivity = [];
for j906 = 1:size(roadGrid_c,2)
	for i906 = 1:size(roadGrid_c,1)
		xml_input.urbanArea.urbanRoad.materials.thermalConductivity = [xml_input.urbanArea.urbanRoad.materials.thermalConductivity roadGrid_c{i906,j906}];
	end
end
xml_input.urbanArea.urbanRoad.materials.volumetricHeatCapacity = [];
for j907 = 1:size(roadGrid_h,2)
	for i907 = 1:size(roadGrid_h,1)
		xml_input.urbanArea.urbanRoad.materials.volumetricHeatCapacity = [xml_input.urbanArea.urbanRoad.materials.volumetricHeatCapacity roadGrid_h{i907,j907}];
	end
end
xml_input.urbanArea.urbanRoad.materials.thickness = xml_input.urbanArea.urbanRoad.materials.thickness';

%% Rural Road
% Break rural into pieces
for i902 = 1:size(xml_input.referenceSite.ruralRoad.materials.thickness,2)
	ruralGrid_t{1,i902} = xml_input.referenceSite.ruralRoad.materials.thickness(i902);
	ruralGrid_n{1,i902} = xml_input.referenceSite.ruralRoad.materials.names(i902);
	ruralGrid_c{1,i902} = xml_input.referenceSite.ruralRoad.materials.thermalConductivity(i902);
	ruralGrid_h{1,i902} = xml_input.referenceSite.ruralRoad.materials.volumetricHeatCapacity(i902);
end
for i903 = 1:size(ruralGrid_t,2)
clear Lr;
Lr = ruralGrid_t{i903};
clear Lm;
if Lr <= 0.05;
				Lm(1)=Lr;
elseif Lr <= 0.1;
				Lm(1)=Lr/2;
				Lm(2)=Lr/2;
elseif Lr <= 0.2;
				Lm(1)=Lr/4;
				Lm(2)=Lr/2;
				Lm(3)=Lr/4;
else Lr <= 0.3;
				Lm(1)=Lr/6;
				Lm(2)=Lr/3;
				Lm(3)=Lr/3;
				Lm(4)=Lr/6;
end

% Reassemble rural
ruralGrid_t{i903} = Lm;
clear names;
clear thermalConductivity;
clear volumetricHeatCapacity;
names = ruralGrid_n{1,i903};
thermalConductivity = ruralGrid_c{1,i903};
volumetricHeatCapacity = ruralGrid_h{1,i903};

for j902 = 1:size(Lm,2)
	ruralGrid_n{j902,i903} = names;
	ruralGrid_c{j902,i903} = thermalConductivity;
	ruralGrid_h{j902,i903} = volumetricHeatCapacity;
end
end
xml_input.referenceSite.ruralRoad.materials.thickness = [];
for i904 = 1:size(ruralGrid_t,2)
	xml_input.referenceSite.ruralRoad.materials.thickness = [xml_input.referenceSite.ruralRoad.materials.thickness ruralGrid_t{i904}];
end
xml_input.referenceSite.ruralRoad.materials.names = [];
for j905 = 1:size(ruralGrid_n,2)
	for i905 = 1:size(ruralGrid_n,1)
		xml_input.referenceSite.ruralRoad.materials.names = [xml_input.referenceSite.ruralRoad.materials.names ruralGrid_n{i905,j905}];
	end
end
xml_input.referenceSite.ruralRoad.materials.thermalConductivity = [];
for j906 = 1:size(ruralGrid_c,2)
	for i906 = 1:size(ruralGrid_c,1)
		xml_input.referenceSite.ruralRoad.materials.thermalConductivity = [xml_input.referenceSite.ruralRoad.materials.thermalConductivity ruralGrid_c{i906,j906}];
	end
end
xml_input.referenceSite.ruralRoad.materials.volumetricHeatCapacity = [];
for j907 = 1:size(ruralGrid_h,2)
	for i907 = 1:size(ruralGrid_h,1)
		xml_input.referenceSite.ruralRoad.materials.volumetricHeatCapacity = [xml_input.referenceSite.ruralRoad.materials.volumetricHeatCapacity ruralGrid_h{i907,j907}];
	end
end
xml_input.referenceSite.ruralRoad.materials.thickness = xml_input.referenceSite.ruralRoad.materials.thickness';

%% 
%IN_MON2=str2num(IN_MON);
%IN_DAY2=str2num(IN_DAY);
%IN_DUR2=str2num(IN_DUR);
% Simulation paramters
simParam = SimParam(...
    300.,...            % Simulation time-step
    3600,...            % Weather data time-step
    1,...               % Begin month
    1,...               % Begin day of the month
    365);

%% for production!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%    IN_MON2,...        % Begin month
%    IN_DAY2,...        % Begin day of the month
%    IN_DUR2);
%

   % xml_input.parameter.simuDuration);               % Number of days of simulation
weather = Weather(climate_data,...
    simParam.timeInitial,simParam.timeFinal);
%    xml_input.simParam.weatherDataTimeStep * 60,...  % Weather data time-step

%% Param
nightSetStart =  0.25 * (xml_input.typology2.building.nightSetStart +  xml_input.typology2.building.nightSetStart +  xml_input.typology3.building.nightSetStart +  xml_input.typology4.building.nightSetStart);
nightSetEnd =  0.25 * (xml_input.typology2.building.nightSetEnd +  xml_input.typology2.building.nightSetEnd +  xml_input.typology3.building.nightSetEnd +  xml_input.typology4.building.nightSetEnd);

parameter = Param(xml_input.urbanArea.daytimeBLHeight,...
    xml_input.urbanArea.nighttimeBLHeight,...
    xml_input.urbanArea.refHeight,...
    xml_input.parameter.tempHeight,...
    xml_input.parameter.windHeight,...
    1.2,...
    200,...
    50,...
    xml_input.urbanArea.treeLatent,...
    xml_input.urbanArea.grassLatent,...
    xml_input.urbanArea.vegAlbedo,...
    xml_input.urbanArea.vegStart,...
    xml_input.urbanArea.vegEnd,...
    nightSetStart,...
    nightSetEnd,...
    0.1,...
    10,...
    0.005,...
    0.3,...
    500,...
    9.81, 1004, 0.40, 287, 461.5, 2260000, 3.141592653, 0.0000000567, 1000, 2500800, 273.16, 611.14, 4218,...
    1846.1, 9.4, 7.4, 1.09647471147);

%%
% Define Wall
wallMat1 = [];
wallMat2 = [];
wallMat3 = [];
wallMat4 = [];

%typ1
for j = 1:size(xml_input.typology1.construction.wall.materials.names,2)
	wallMat1 = [wallMat1 Material(xml_input.typology1.construction.wall.materials.thermalConductivity{j},xml_input.typology1.construction.wall.materials.volumetricHeatCapacity{j})];
end
wall1 = Element(xml_input.typology1.construction.wall.albedo,...
xml_input.typology1.construction.wall.emissivity,...
xml_input.typology1.construction.wall.materials.thickness,...
wallMat1,...
xml_input.typology1.construction.wall.vegetationCoverage,...
xml_input.typology1.construction.wall.initialTemperature + 273.15,...
xml_input.typology1.construction.wall.inclination);

%typ2
for j = 1:size(xml_input.typology2.construction.wall.materials.names,2)
	wallMat2 = [wallMat2 Material(xml_input.typology2.construction.wall.materials.thermalConductivity{j},xml_input.typology2.construction.wall.materials.volumetricHeatCapacity{j})];
end
wall2 = Element(xml_input.typology2.construction.wall.albedo,...
xml_input.typology2.construction.wall.emissivity,...
xml_input.typology2.construction.wall.materials.thickness,...
wallMat2,...
xml_input.typology2.construction.wall.vegetationCoverage,...
xml_input.typology2.construction.wall.initialTemperature + 273.15,...
xml_input.typology2.construction.wall.inclination);

%typ3
for j = 1:size(xml_input.typology3.construction.wall.materials.names,2)
	wallMat3 = [wallMat3 Material(xml_input.typology3.construction.wall.materials.thermalConductivity{j},xml_input.typology3.construction.wall.materials.volumetricHeatCapacity{j})];
end
wall3 = Element(xml_input.typology3.construction.wall.albedo,...
xml_input.typology3.construction.wall.emissivity,...
xml_input.typology3.construction.wall.materials.thickness,...
wallMat3,...
xml_input.typology3.construction.wall.vegetationCoverage,...
xml_input.typology3.construction.wall.initialTemperature + 273.15,...
xml_input.typology3.construction.wall.inclination);

%typ4
for j = 1:size(xml_input.typology4.construction.wall.materials.names,2)
	wallMat4 = [wallMat4 Material(xml_input.typology4.construction.wall.materials.thermalConductivity{j},xml_input.typology4.construction.wall.materials.volumetricHeatCapacity{j})];
end
wall4 = Element(xml_input.typology4.construction.wall.albedo,...
xml_input.typology4.construction.wall.emissivity,...
xml_input.typology4.construction.wall.materials.thickness,...
wallMat4,...
xml_input.typology4.construction.wall.vegetationCoverage,...
xml_input.typology4.construction.wall.initialTemperature + 273.15,...
xml_input.typology4.construction.wall.inclination);

%% Define Roof
roofMat1 = [];
roofMat2 = [];
roofMat3 = [];
roofMat4 = [];

%typ1
for j = 1:size(xml_input.typology1.construction.roof.materials.names,2)
	roofMat1 = [roofMat1 Material(xml_input.typology1.construction.roof.materials.thermalConductivity{j},xml_input.typology1.construction.roof.materials.volumetricHeatCapacity{j})];
end
roof1 = Element(xml_input.typology1.construction.roof.albedo,...
xml_input.typology1.construction.roof.emissivity,...
xml_input.typology1.construction.roof.materials.thickness,...
roofMat1,...
xml_input.typology1.construction.roof.vegetationCoverage,...
xml_input.typology1.construction.roof.initialTemperature + 273.15,...
xml_input.typology1.construction.roof.inclination);

%typ2
for j = 1:size(xml_input.typology2.construction.roof.materials.names,2)
	roofMat2 = [roofMat2 Material(xml_input.typology2.construction.roof.materials.thermalConductivity{j},xml_input.typology2.construction.roof.materials.volumetricHeatCapacity{j})];
end
roof2 = Element(xml_input.typology2.construction.roof.albedo,...
xml_input.typology2.construction.roof.emissivity,...
xml_input.typology2.construction.roof.materials.thickness,...
roofMat2,...
xml_input.typology2.construction.roof.vegetationCoverage,...
xml_input.typology2.construction.roof.initialTemperature + 273.15,...
xml_input.typology2.construction.roof.inclination);

%typ3
for j = 1:size(xml_input.typology3.construction.roof.materials.names,2)
	roofMat3 = [roofMat3 Material(xml_input.typology3.construction.roof.materials.thermalConductivity{j},xml_input.typology3.construction.roof.materials.volumetricHeatCapacity{j})];
end
roof3 = Element(xml_input.typology3.construction.roof.albedo,...
xml_input.typology3.construction.roof.emissivity,...
xml_input.typology3.construction.roof.materials.thickness,...
roofMat3,...
xml_input.typology3.construction.roof.vegetationCoverage,...
xml_input.typology3.construction.roof.initialTemperature + 273.15,...
xml_input.typology3.construction.roof.inclination);

%typ4
for j = 1:size(xml_input.typology4.construction.roof.materials.names,2)
	roofMat4 = [roofMat4 Material(xml_input.typology4.construction.roof.materials.thermalConductivity{j},xml_input.typology4.construction.roof.materials.volumetricHeatCapacity{j})];
end
roof4 = Element(xml_input.typology4.construction.roof.albedo,...
xml_input.typology4.construction.roof.emissivity,...
xml_input.typology4.construction.roof.materials.thickness,...
roofMat4,...
xml_input.typology4.construction.roof.vegetationCoverage,...
xml_input.typology4.construction.roof.initialTemperature + 273.15,...
xml_input.typology4.construction.roof.inclination);

% Define Mass
%tup1
massMat1 = [];
massMat2 = [];
massMat3 = [];
massMat4 = [];

for j = 1:size(xml_input.typology1.construction.mass.materials.names,2)
	massMat1 = [massMat1 Material(xml_input.typology1.construction.mass.materials.thermalConductivity{j},xml_input.typology1.construction.mass.materials.volumetricHeatCapacity{j})];
end
mass1 = Element(xml_input.typology1.construction.mass.albedo,...
xml_input.typology1.construction.mass.emissivity,...
xml_input.typology1.construction.mass.materials.thickness,...
massMat1,...
xml_input.typology1.construction.mass.vegetationCoverage,...
xml_input.typology1.construction.mass.initialTemperature + 273.15,...
xml_input.typology1.construction.mass.inclination);

%typ2
for j = 1:size(xml_input.typology2.construction.mass.materials.names,2)
	massMat2 = [massMat2 Material(xml_input.typology2.construction.mass.materials.thermalConductivity{j},xml_input.typology2.construction.mass.materials.volumetricHeatCapacity{j})];
end
mass2 = Element(xml_input.typology2.construction.mass.albedo,...
xml_input.typology2.construction.mass.emissivity,...
xml_input.typology2.construction.mass.materials.thickness,...
massMat2,...
xml_input.typology2.construction.mass.vegetationCoverage,...
xml_input.typology2.construction.mass.initialTemperature + 273.15,...
xml_input.typology2.construction.mass.inclination);

%typ3
for j = 1:size(xml_input.typology3.construction.mass.materials.names,2)
	massMat3 = [massMat3 Material(xml_input.typology3.construction.mass.materials.thermalConductivity{j},xml_input.typology3.construction.mass.materials.volumetricHeatCapacity{j})];
end
mass3 = Element(xml_input.typology3.construction.mass.albedo,...
xml_input.typology3.construction.mass.emissivity,...
xml_input.typology3.construction.mass.materials.thickness,...
massMat3,...
xml_input.typology3.construction.mass.vegetationCoverage,...
xml_input.typology3.construction.mass.initialTemperature + 273.15,...
xml_input.typology3.construction.mass.inclination);
	
%typ4
for j = 1:size(xml_input.typology4.construction.mass.materials.names,2)
	massMat4 = [massMat4 Material(xml_input.typology4.construction.mass.materials.thermalConductivity{j},xml_input.typology4.construction.mass.materials.volumetricHeatCapacity{j})];
end
mass4 = Element(xml_input.typology4.construction.mass.albedo,...
xml_input.typology4.construction.mass.emissivity,...
xml_input.typology4.construction.mass.materials.thickness,...
massMat4,...
xml_input.typology4.construction.mass.vegetationCoverage,...
xml_input.typology4.construction.mass.initialTemperature + 273.15,...
xml_input.typology4.construction.mass.inclination);

% Define Road
roadMat = [];
for j = 1:size(xml_input.urbanArea.urbanRoad.materials.names,2)
	roadMat = [roadMat Material(xml_input.urbanArea.urbanRoad.materials.thermalConductivity{j},xml_input.urbanArea.urbanRoad.materials.volumetricHeatCapacity{j})];
end
road = Element(xml_input.urbanArea.urbanRoad.albedo,...
xml_input.urbanArea.urbanRoad.emissivity,...
xml_input.urbanArea.urbanRoad.materials.thickness,...
roadMat,...
xml_input.urbanArea.urbanRoad.vegetationCoverage,...
xml_input.urbanArea.urbanRoad.initialTemperature + 273.15,...
xml_input.urbanArea.urbanRoad.inclination);

% Define Rural
ruralMat = [];
for j = 1:size(xml_input.referenceSite.ruralRoad.materials.names,2)
	ruralMat = [ruralMat Material(xml_input.referenceSite.ruralRoad.materials.thermalConductivity{j},xml_input.referenceSite.ruralRoad.materials.volumetricHeatCapacity{j})];
end
rural = Element(xml_input.referenceSite.ruralRoad.albedo,...
xml_input.referenceSite.ruralRoad.emissivity,...
xml_input.referenceSite.ruralRoad.materials.thickness,...
ruralMat,...
xml_input.referenceSite.ruralRoad.vegetationCoverage,...
xml_input.referenceSite.ruralRoad.initialTemperature + 273.15,...
xml_input.referenceSite.ruralRoad.inclination);

%% Define Buildings
%typ1
typology1 = Building(xml_input.typology1.building.floorHeight,...
xml_input.typology1.building.nightInternalGains,...
xml_input.typology1.building.dayInternalGains,...
xml_input.typology1.building.radiantFraction,...
xml_input.typology1.building.latentFraction,...
xml_input.typology1.building.infiltration,...
xml_input.typology1.building.ventilation,...
xml_input.typology1.construction.glazing.glazingRatio,...
xml_input.typology1.construction.glazing.windowUvalue,...
xml_input.typology1.construction.glazing.windowSHGC,...
xml_input.typology1.building.coolingSystemType,...
xml_input.typology1.building.coolingCOP,...
xml_input.typology1.building.heatReleasedToCanyon,...
xml_input.typology1.building.daytimeCoolingSetPoint + 273.15,...
xml_input.typology1.building.nighttimeCoolingSetPoint + 273.15,...
xml_input.typology1.building.daytimeHeatingSetPoint + 273.15,...
xml_input.typology1.building.nighttimeHeatingSetPoint + 273.15,...
xml_input.typology1.building.coolingCapacity,...
xml_input.typology1.building.heatingEfficiency,...
xml_input.typology1.building.initialT + 273.15);

%typ2
typology2 = Building(xml_input.typology2.building.floorHeight,...
xml_input.typology2.building.nightInternalGains,...
xml_input.typology2.building.dayInternalGains,...
xml_input.typology2.building.radiantFraction,...
xml_input.typology2.building.latentFraction,...
xml_input.typology2.building.infiltration,...
xml_input.typology2.building.ventilation,...
xml_input.typology2.construction.glazing.glazingRatio,...
xml_input.typology2.construction.glazing.windowUvalue,...
xml_input.typology2.construction.glazing.windowSHGC,...
xml_input.typology2.building.coolingSystemType,...
xml_input.typology2.building.coolingCOP,...
xml_input.typology2.building.heatReleasedToCanyon,...
xml_input.typology2.building.daytimeCoolingSetPoint + 273.15,...
xml_input.typology2.building.nighttimeCoolingSetPoint + 273.15,...
xml_input.typology2.building.daytimeHeatingSetPoint + 273.15,...
xml_input.typology2.building.nighttimeHeatingSetPoint + 273.15,...
xml_input.typology2.building.coolingCapacity,...
xml_input.typology2.building.heatingEfficiency,...
xml_input.typology2.building.initialT + 273.15);

%typ3
typology3 = Building(xml_input.typology3.building.floorHeight,...
xml_input.typology3.building.nightInternalGains,...
xml_input.typology3.building.dayInternalGains,...
xml_input.typology3.building.radiantFraction,...
xml_input.typology3.building.latentFraction,...
xml_input.typology3.building.infiltration,...
xml_input.typology3.building.ventilation,...
xml_input.typology3.construction.glazing.glazingRatio,...
xml_input.typology3.construction.glazing.windowUvalue,...
xml_input.typology3.construction.glazing.windowSHGC,...
xml_input.typology3.building.coolingSystemType,...
xml_input.typology3.building.coolingCOP,...
xml_input.typology3.building.heatReleasedToCanyon,...
xml_input.typology3.building.daytimeCoolingSetPoint + 273.15,...
xml_input.typology3.building.nighttimeCoolingSetPoint + 273.15,...
xml_input.typology3.building.daytimeHeatingSetPoint + 273.15,...
xml_input.typology3.building.nighttimeHeatingSetPoint + 273.15,...
xml_input.typology3.building.coolingCapacity,...
xml_input.typology3.building.heatingEfficiency,...
xml_input.typology3.building.initialT + 273.15);

%typ4
typology4 = Building(xml_input.typology4.building.floorHeight,...
xml_input.typology4.building.nightInternalGains,...
xml_input.typology4.building.dayInternalGains,...
xml_input.typology4.building.radiantFraction,...
xml_input.typology4.building.latentFraction,...
xml_input.typology4.building.infiltration,...
xml_input.typology4.building.ventilation,...
xml_input.typology4.construction.glazing.glazingRatio,...
xml_input.typology4.construction.glazing.windowUvalue,...
xml_input.typology4.construction.glazing.windowSHGC,...
xml_input.typology4.building.coolingSystemType,...
xml_input.typology4.building.coolingCOP,...
xml_input.typology4.building.heatReleasedToCanyon,...
xml_input.typology4.building.daytimeCoolingSetPoint + 273.15,...
xml_input.typology4.building.nighttimeCoolingSetPoint + 273.15,...
xml_input.typology4.building.daytimeHeatingSetPoint + 273.15,...
xml_input.typology4.building.nighttimeHeatingSetPoint + 273.15,...
xml_input.typology4.building.coolingCapacity,...
xml_input.typology4.building.heatingEfficiency,...
xml_input.typology4.building.initialT + 273.15);

% Urban Configuration [building,mass,wall,roof,road]
urbanConf1 = UrbanConf(typology1,mass1,wall1,roof1,road);
urbanConf2 = UrbanConf(typology2,mass2,wall2,roof2,road);
urbanConf3 = UrbanConf(typology3,mass3,wall3,roof3,road);    
urbanConf4 = UrbanConf(typology4,mass4,wall4,roof4,road);    

%find distro
typ1_distro = xml_input.typology1.dist;
typ2_distro = xml_input.typology2.dist;
typ3_distro = xml_input.typology3.dist;
typ4_distro = xml_input.typology4.dist;

% Urban Usage [Fraction of urban configurations,urban configurations]
urbanUsage = UrbanUsage([typ1_distro/100, typ2_distro/100, typ3_distro/100, typ4_distro/100],[urbanConf1, urbanConf2, urbanConf3, urbanConf4]);
   
%% Reference site
refSite = ReferenceSite(xml_input.referenceSite.latitude,...
    xml_input.referenceSite.longitude,...
    xml_input.referenceSite.averageObstacleHeight,...
    weather.staTemp(1),weather.staPres(1),parameter);

%% Urban areas

% %weightedAverage
% floorHeight_avg = (typ1_distro/100)* xml_input.typology1.building.floorHeight...
%     + (typ2_distro/100)* xml_input.typology2.building.floorHeight...
%     + (typ3_distro/100)* xml_input.typology3.building.floorHeight...
%     + (typ4_distro/100)* xml_input.typology4.building.floorHeight;
% nightInternalGains_avg = (typ1_distro/100)* xml_input.typology1.building.nightInternalGains...
%     + (typ2_distro/100)* xml_input.typology2.building.nightInternalGains...
%     + (typ3_distro/100)* xml_input.typology3.building.nightInternalGains...
%     + (typ4_distro/100)* xml_input.typology4.building.nightInternalGains;
% dayInternalGains_avg = (typ1_distro/100)* xml_input.typology1.building.dayInternalGains...
%     + (typ2_distro/100)* xml_input.typology2.building.dayInternalGains...
%     + (typ3_distro/100)* xml_input.typology3.building.dayInternalGains...
%     + (typ4_distro/100)* xml_input.typology4.building.dayInternalGains;
% radiantFraction_avg = (typ1_distro/100)* xml_input.typology1.building.radiantFraction...
%     + (typ2_distro/100)* xml_input.typology2.building.radiantFraction...
%     + (typ3_distro/100)* xml_input.typology3.building.radiantFraction...
%     + (typ4_distro/100)* xml_input.typology4.building.radiantFraction;
% latentFraction_avg = (typ1_distro/100)* xml_input.typology1.building.latentFraction...
%     + (typ2_distro/100)* xml_input.typology2.building.latentFraction...
%     + (typ3_distro/100)* xml_input.typology3.building.latentFraction...
%     + (typ4_distro/100)* xml_input.typology4.building.latentFraction;
% infiltration_avg = (typ1_distro/100)* xml_input.typology1.building.infiltration...
%     + (typ2_distro/100)* xml_input.typology2.building.infiltration...
%     + (typ3_distro/100)* xml_input.typology3.building.infiltration...
%     + (typ4_distro/100)* xml_input.typology4.building.infiltration;
% glazingRatio_avg = (typ1_distro/100)* xml_input.typology1.construction.glazing.glazingRatio...
%     + (typ2_distro/100)* xml_input.typology2.construction.glazing.glazingRatio...
%     + (typ3_distro/100)* xml_input.typology3.construction.glazing.glazingRatio...
%     + (typ4_distro/100)* xml_input.typology4.construction.glazing.glazingRatio;
% windowUvalue_avg = (typ1_distro/100)* xml_input.typology1.construction.glazing.windowUvalue...
%     + (typ2_distro/100)* xml_input.typology2.construction.glazing.windowUvalue...
%     + (typ3_distro/100)* xml_input.typology3.construction.glazing.windowUvalue...
%     + (typ4_distro/100)* xml_input.typology4.construction.glazing.windowUvalue;
% windowSHGC_avg = (typ1_distro/100)* xml_input.typology1.construction.glazing.windowSHGC...
%     + (typ2_distro/100)* xml_input.typology2.construction.glazing.windowSHGC...
%     + (typ3_distro/100)* xml_input.typology3.construction.glazing.windowSHGC...
%     + (typ4_distro/100)* xml_input.typology4.construction.glazing.windowSHGC;
% coolingSystemType_avg = (typ1_distro/100)* xml_input.typology1.building.coolingSystemType...
%     + (typ2_distro/100)* xml_input.typology2.building.coolingSystemType...
%     + (typ3_distro/100)* xml_input.typology3.building.coolingSystemType...
%     + (typ4_distro/100)* xml_input.typology4.building.coolingSystemType;
% 
% if (typ1_distro > typ2_distro && typ1_distro > typ3_distro && typ1_distro > typ4_distro)
%     coolingSystemType_avg = xml_input.typology1.building.coolingSystemType;
% elseif (typ2_distro > typ1_distro && typ2_distro > typ3_distro && typ2_distro > typ4_distro)
%     coolingSystemType_avg = xml_input.typology2.building.coolingSystemType;
% elseif (typ3_distro > typ1_distro && typ3_distro > typ2_distro && typ3_distro > typ4_distro)
%     coolingSystemType_avg = xml_input.typology3.building.coolingSystemType;
% elseif (typ4_distro > typ1_distro && typ4_distro > typ2_distro && typ4_distro > typ3_distro)
%     coolingSystemType_avg = xml_input.typology4.building.coolingSystemType;
% else coolingSystemType_avg = xml_input.typology1.building.coolingSystemType;
% end 
% 
% coolingCOP_avg = (typ1_distro/100)* xml_input.typology1.building.coolingCOP...
%     + (typ2_distro/100)* xml_input.typology2.building.coolingCOP...
%     + (typ3_distro/100)* xml_input.typology3.building.coolingCOP...
%     + (typ4_distro/100)* xml_input.typology4.building.coolingCOP;
% heatReleasedToCanyon_avg = (typ1_distro/100)* xml_input.typology1.building.heatReleasedToCanyon...
%     + (typ2_distro/100)* xml_input.typology2.building.heatReleasedToCanyon...
%     + (typ3_distro/100)* xml_input.typology3.building.heatReleasedToCanyon...
%     + (typ4_distro/100)* xml_input.typology4.building.heatReleasedToCanyon;
% daytimeCoolingSetPoint_avg = 273.15 + (typ1_distro/100)* xml_input.typology1.building.daytimeCoolingSetPoint...
%     + (typ2_distro/100)* xml_input.typology2.building.daytimeCoolingSetPoint...
%     + (typ3_distro/100)* xml_input.typology3.building.daytimeCoolingSetPoint...
%     + (typ4_distro/100)* xml_input.typology4.building.daytimeCoolingSetPoint;
% nighttimeCoolingSetPoint_avg = 273.15 + (typ1_distro/100)* xml_input.typology1.building.nighttimeCoolingSetPoint...
%     + (typ2_distro/100)* xml_input.typology2.building.nighttimeCoolingSetPoint...
%     + (typ3_distro/100)* xml_input.typology3.building.nighttimeCoolingSetPoint...
%     + (typ4_distro/100)* xml_input.typology4.building.nighttimeCoolingSetPoint;
% daytimeHeatingSetPoint_avg = 273.15 + (typ1_distro/100)* xml_input.typology1.building.daytimeHeatingSetPoint...
%     + (typ2_distro/100)* xml_input.typology2.building.daytimeHeatingSetPoint...
%     + (typ3_distro/100)* xml_input.typology3.building.daytimeHeatingSetPoint...
%     + (typ4_distro/100)* xml_input.typology4.building.daytimeHeatingSetPoint;
% nighttimeHeatingSetPoint_avg = 273.15 + (typ1_distro/100)* xml_input.typology1.building.nighttimeHeatingSetPoint...
%     + (typ2_distro/100)* xml_input.typology2.building.nighttimeHeatingSetPoint...
%     + (typ3_distro/100)* xml_input.typology3.building.nighttimeHeatingSetPoint...
%     + (typ4_distro/100)* xml_input.typology4.building.nighttimeHeatingSetPoint;
% coolingCapacity_avg = (typ1_distro/100)* xml_input.typology1.building.coolingCapacity...
%     + (typ2_distro/100)* xml_input.typology2.building.coolingCapacity...
%     + (typ3_distro/100)* xml_input.typology3.building.coolingCapacity...
%     + (typ4_distro/100)* xml_input.typology4.building.coolingCapacity;
% heatingEfficiency_avg = (typ1_distro/100)* xml_input.typology1.building.heatingEfficiency...
%     + (typ2_distro/100)* xml_input.typology2.building.heatingEfficiency...
%     + (typ3_distro/100)* xml_input.typology3.building.heatingEfficiency...
%     + (typ4_distro/100)* xml_input.typology4.building.heatingEfficiency;
% initialT_avg = 273.15 + (typ1_distro/100)* xml_input.typology1.building.initialT...
%     + (typ2_distro/100)* xml_input.typology2.building.initialT...
%     + (typ3_distro/100)* xml_input.typology3.building.initialT...
%     + (typ4_distro/100)* xml_input.typology4.building.initialT;
% 
% typology_avg = Building(floorHeight_avg,...
% nightInternalGains_avg,...
% dayInternalGains_avg,...
% radiantFraction_avg,...
% latentFraction_avg,...
% infiltration_avg,...
% ventilation_avg,...
% glazingRatio_avg,...
% windowUvalue_avg,...
% windowSHGC_avg,...
% coolingSystemType_avg,...
% coolingCOP_avg,...
% heatReleasedToCanyon_avg,...
% xdaytimeCoolingSetPoint_avg,...
% nighttimeCoolingSetPoint_avg,...
% daytimeHeatingSetPoint_avg,...
% nighttimeHeatingSetPoint_avg,...
% coolingCapacity_avg,...
% heatingEfficiency_avg,...
% initialT_avg);

if (typ1_distro > typ2_distro && typ1_distro > typ3_distro && typ1_distro > typ4_distro)
    dominantTypology = typology1;
    dominantWall = wa111;
elseif (typ2_distro > typ1_distro && typ2_distro > typ3_distro && typ2_distro > typ4_distro)
    dominantTypology = typology2;
    dominantWall = wa112;
elseif (typ3_distro > typ1_distro && typ3_distro > typ2_distro && typ3_distro > typ4_distro)
    dominantTypology = typology3;
    dominantWall = wa113;
elseif (typ4_distro > typ1_distro && typ4_distro > typ2_distro && typ4_distro > typ3_distro)
    dominantTypology = typology4;
    dominantWall = wa114;
else
    dominantTypology = typology1;
    dominantWall = wall1;
end 

disp(dominantTypology);

urbanArea = UrbanArea(xml_input.urbanArea.averageBuildingHeight,...               
    xml_input.urbanArea.siteCoverageRatio,...              
    xml_input.urbanArea.facadeToSiteRatio,...              
    xml_input.urbanArea.treeCoverage,...              
    xml_input.urbanArea.nonBldgSensibleHeat,...               
    xml_input.urbanArea.nonBldgLatentAnthropogenicHeat,...               
    weather.staTemp(1),weather.staHum(1),weather.staUmod(1),parameter,...
    dominantTypology,dominantWall,road,rural); 

%% UblVars
ublVars = UblVars('C',...
    xml_input.urbanArea.charLength,...
    weather.staTemp(1),parameter.maxdx); 

%% HVAC autosize
autosize = 1;
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
%% Run UWG

forc = Forcing(weather.staTemp);
sensHeat = zeros(numel(urbanArea),1);

disp '==========================='
disp 'Start UWG simulation '
disp '==========================='
% John's code, no printfile section
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
            
            %output values for radTemp
            %radTemp_calc = [urbanArea(1).canTemp forc.skyTemp urbanUsage(1).urbanConf(1).wall.layerTemp(1) urbanUsage(1).urbanConf(1).road.layerTemp(1) urbanArea(1).canWidth];
            radTemp = [urbanUsage(1).urbanConf(1).wall.layerTemp(1) forc.skyTemp  urbanUsage(1).urbanConf(1).road.layerTemp(1)];
            radTemp_calc = [radTemp, urbanArea(1).bldHeight urbanArea(1).canWidth];
            filename = strcat(epwPathName,'\','Trad_',xmlFileName,'.csv');
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
    % print day
    if eq(mod(timeCount,3600.*24),0)
      varname = timeCount/(3600.*24);
      fprintf('DAY %g\n', varname);
      
% end of bruno's code
    
      progressbar(varname/365)
    end
end
%% Write modified values to epwinput.value

disp('Calculating new Temperature and humidity values')
for iJ = 1:size(Can_Tdb,2)
    epwinput.values{iJ,7}{1,1} = num2str(Can_Tdb(iJ),'%0.1f'); % dry bulb temperature  [�C]
    
    [Tdb, w, Can_phi(iJ), h, Can_Tdp(iJ), v] = Psychrometrics(Can_Tdb(iJ), Can_hum(iJ), str2num(epwinput.values{iJ,10}{1,1}));
    
    epwinput.values{iJ,8}{1,1} = num2str(Can_Tdp(iJ),'%0.1f'); % dew point temperature [�C]
    epwinput.values{iJ,9}{1,1} = num2str(Can_phi(iJ),'%0.0f'); % relative humidity     [%]
end

    progressbar(360/365)
disp '==========================='
disp '    UWG ended correctly    '
disp '==========================='


%% Write new EPW file

%from old xml code /needed
disp('writing new EPW file')

%%FOR PRODUCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%new_climate_file = strcat(CL_RE_PATH,'\',xmlFileName,'.epw');
%%END FOR PRODUCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_climate_file = strcat(epwPathName,'\',xmlFileName,'.epw');
if exist(new_climate_file)
delete(new_climate_file)
end

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
%end of old xml code

progressbar(1)
if fullyScripted
    disp('Inputs scripted, supressing pop-up notification...')
else
    h = msgbox('Urban Weather Generation Complete','UWG 2.0','help');
end
disp(['New climate file generated: ',new_climate_file])

fclose all;
end