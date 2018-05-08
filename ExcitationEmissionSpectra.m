delete(instrfindall) 
clear all
addpath('C:\Users\FontenotRS\Desktop\Excitation Emission\')
%% Variables for the Measurement

% PMT Settings

voltage = 'D'; %Default PMT Voltage Setting

% Number of Scans
NumRuns = 1; %This is the number of runs to perform on one excitation run

%% Excitation Monochromator Settings

% SPEX 1680 Double Monochromator Settings 
MinWavelength = 300; %Starting Wavelength
SPEX_Wavelength = 400.2; %Wavelength on the SPEX 1680 with adjustment
End_Wavelength = 700; %This sets the monochromator to the CdSe peak (nm)
Inc_Wavelength = 0.25; %Step size in nm and must be in increments of 1/50 or 0.05

%% Emission Monochromator Settings

Emission_Wavelength = 636; %This sets the single monochromator to the emission you want to observe
Ent_Slit_Width = 9; %Unit is in mm and must be greater than 0 and less than 11.



%% Checking inputs

if MinWavelength > 1100 || MinWavelength < 100
    error('Starting wavelength is outside of the SPEX range. Choose between 100 and 1100 nm!')
end

if End_Wavelength > 1100 || End_Wavelength < 100
    error('Ending wavelength is outside of the SPEX range. Choose between 100 and 1100 nm!')
end

if SPEX_Wavelength > 1100 || SPEX_Wavelength < 100
    error('SPEX wavelength is outside of the SPEX range. Choose between 100 and 1100 nm!')
end

if   ~mod(Inc_Wavelength,0.05)~=1
     error('The step size must be divisible by 0.05 or 1/50. Choose a different step size!')
end

if Ent_Slit_Width > 11.0 || Ent_Slit_Width <=0
    error('The slit width must be greater than 0 and less than 11')
end

if Ent_Slit_Width*1.23 > 11
    Exit_Slit_Width = 11.0
else
    Exit_Slit_Width = Ent_Slit_Width*1.23  %This lets more light in without losing the resolution
end

if Emission_Wavelength > 1100 || Emission_Wavelength < 0
    error('Emission wavelength is outside of the SPEX 270M range. Choose between 100 and 1100 nm!')
end

%% Setting up the locations to Save the Files

dname = 'C:\Users\FontenotRS\Desktop\Excitation Results\'; %This sets the directory where the files to be analyzed are created.
formatOut = 'mm-dd-yy';
DateVector = datetime('now');
Date = datestr(DateVector, formatOut);
exist Date;
test = ans;
dash = '\';
if test~=7
    mkdir(dname,Date); % This creates the days folder
end

WorkingDir = fullfile(dname,Date); %This checks to see if there are any files
cd(WorkingDir);
listOfRuns = ls('Run*');

if isempty(listOfRuns)
    len = 1;
    CurrentNum = 1;
else
    len = length(listOfRuns(:,1)); %This counts how many folders are currently in the folder
    for a=1:1:len
    LastRun = listOfRuns(a,:); %This calls up the last Run folder
    len2 = length(LastRun)-4;
    if len2 == 1
        num(a) = str2num(LastRun(5));
    else
        num(a) = str2num(LastRun(5:4+len2));
    end
    CurrentNum = max(num) + 1;
    end
end

%This finds the highest number for the folder and takes into account
%Matlabs stupid sorting


NewDirFolder = sprintf('Run %d', CurrentNum);

mkdir(WorkingDir, NewDirFolder)

SavingDir = fullfile(dname,Date,NewDirFolder);
cd(SavingDir) %This puts the program in the desired directory to save the files

%Creating the run file
fid1=fopen(sprintf('Run %d Excitation Spectra.txt', CurrentNum), 'w');
fid2=fopen(sprintf('Run %d Experimental Parameters.txt', CurrentNum), 'w');
fprintf(fid2, 'This is an excitation spectra run. \r\n');
fprintf(fid2, 'The voltage for this test is %s V. \r\n',voltage);
fprintf(fid2, 'The starting wavelength for this test is %d nm. \r\n',MinWavelength); 
fprintf(fid2, 'The ending wavelength is %d nm. \r\n',End_Wavelength); 
fprintf(fid2, 'The step size is %d nm. \r\n',Inc_Wavelength); 
fprintf(fid2, 'The emission wavelength is %0.4f nm. \r\n',Emission_Wavelength);
fprintf(fid2, 'The enterance slit width is %0.4f nm. \r\n',Ent_Slit_Width);
fprintf(fid2, 'The exit slit width is %0.4f nm. \r\n',Exit_Slit_Width);
fprintf(fid1, '%s \t %s\r\n','Wavelength (nm)', 'Counts'); 


%% Setting up communications with the SPEX Single and Double Monochromators

%Double Monochromator
SPEX=gpib('ni',0,1); %creates the gpib object for the SPEX
fopen(SPEX)
SPEX.EOIMode = 'off';
SPEX.EOSMode = 'none';
fwrite(SPEX, uint8([32]));
check = char(fread(SPEX,1));

%Single Monochromator
SPEX_Em=gpib('ni',1,3); %creates the gpib object for the SPEX
fopen(SPEX_Em)
SPEX_Em.EOIMode = 'off';
SPEX_Em.EOSMode = 'none';
fwrite(SPEX_Em, uint8([32]));
check = char(fread(SPEX_Em,1));

if check == 'B'
    fwrite(SPEX,uint8([79; 50; 48; 48; 48; 0])) %Mean O2000\0
    char(fread(SPEX,1))
    fwrite(SPEX_Em,uint8([79; 50; 48; 48; 48; 0])) %Mean O2000\0
    char(fread(SPEX_Em,1))
    pause(2)
    fwrite(SPEX,uint8([65])) %A (initialze motor)
    fwrite(SPEX_Em,uint8([65])) %A (initialze motor)
        h = waitbar(10,'Initializing Motor...');
        steps = 12;
        for step = 1:steps
            pause(10)
            waitbar(step / steps)
        end
    close(h)
    char(fread(SPEX,1))
    char(fread(SPEX_Em,1))
    disp('SPEX Monochromator successfully initialized')
elseif check == 'o'
    set(handles.Messages, 'string', 'Monochromator is not in the correct state. Restart!')
    error('Monochromator is not in the correct state. Restart!')
elseif check == 'F'  
    disp( 'Monochromator is not in the correct state. Restart!')
    error('Monochromator is not in the correct state. Restart!')
else 
    disp('Monochromator in weird state. Restart!')
    error('Monochromator in weird state. Restart!')
end

SPEX.EOSMode = 'read';
SPEX.EOSCharCode = 13;

%% Setting up communications with the PMT

PMT=serial('COM7'); %creates the serial object in this case PMT tube

% This writes the correct PMT settings specified by the manufacturer.

PMT.BaudRate=9600; %The following are the manufacturer specifications
PMT.DataBits=8; 
PMT.StopBits=1; 
PMT.Terminator='CR';
PMT.Parity='none';
PMT.FlowControl='none';

% Communicating with the PMT

fopen(PMT);
pause(0.5)
PMT.ReadAsyncMode='continuous';
disp('PMT setup successful');

%% Setting the voltage

if voltage == 'D'
    fwrite(PMT,uint8(hex2dec(['44'; '0D'])))
    pause(5)
    PMTErrorCheck(1,PMT)
    disp('Voltage Set to Default')
elseif voltage >= 0 & voltage <= 1200
    fwrite(PMT,uint8([86; bitshift(voltage,-8); bitand(voltage,255); 13])) %set voltage to #
    PMTErrorCheck(1,PMT)
    sprintf('Voltage Set to %d V.', voltage)
else
    fwrite(PMT,uint8([86; bitshift(0,-8); bitand(0,255); 13]))
    error('Invalid Voltage! Program shutdown and voltage set to 0')
end

%% Moving the Single Monochromator to the Desired Wavelength

PauseTime = Inc_Wavelength*0.1;
StepsToMove = double(transpose(int2str(((MinWavelength-SPEX_Wavelength)*50))));
fwrite(SPEX,uint8([70; 48; 44; StepsToMove; 13])); %Moves the grating!
char(fread(SPEX,1));
SPEX_Wavelength = MinWavelength;

ptime(120);
MaxWavelength = End_Wavelength;
Wavelength_Plot=zeros(length(1:1:NumRuns),length(MinWavelength:Inc_Wavelength:MaxWavelength));

%% Setting up the Slits

SPEX_Em.EOSMode = 'read';
SPEX_Em.EOSCharCode = 13;

%Checking the status of the motor
fwrite(SPEX_Em, uint8([69]))
busy = fread(SPEX_Em,2)
 while busy(2) == 113
     disp('Motor is busy')
     ptime(30)
     fwrite(SPEX_Em, uint8([69]))
     busy = fread(SPEX_Em,1);
 end

%Reading the Position of the Enterance slit 
fwrite(SPEX_Em,uint8([106;48;44;48;13]))
sline = 1; i=1;
while sline ~= 13
    sline (i) = fread(SPEX_Em,1);
    i=i+1;
end
startSlitEntPosition = str2num(char(sline(2:length(sline)-1)'));

%Moving the Enterance slit to the correct opening
SlitStepsToMove = double(transpose(int2str((Ent_Slit_Width*80 - startSlitEntPosition))))
fwrite(SPEX_Em, uint8([107;48;44;48;44; SlitStepsToMove; 13]))
ptime(60)
fread(SPEX_Em,1)

%Reading the Position of the Exit Slit
fwrite(SPEX_Em,uint8([106;50;44;48;13]))
sline = 1; i=1;
while sline ~= 13
    sline(i) = fread(SPEX_Em,1);
    i=i+1;
end

%Moving the Exit slit the correct opening
startSlitExitPosition = str2num(char(sline(2:length(sline)-1)'))
ExitSlitStepsToMove = double(transpose(int2str((Exit_Slit_Width*80 - startSlitExitPosition))))
fwrite(SPEX_Em, uint8([107;48;44;51;44; ExitSlitStepsToMove; 13]))
ptime(60)
fread(SPEX_Em,1)

%% Setting the SPEX 270M emission Wavlength 

fwrite(SPEX, uint8([72;48;13]))
pause(1)
tline = fread(SPEX)
Current_270M_Wavelength = 0.03125*str2num(char(tline(2:length(tline)-1)'))-7.8415;
%The 7.5729 takes into the account the callibration error of the 270M. 

StepsToMove270M=double(transpose(int2str(((Emission_Wavelength-Current_270M_Wavelength)*32))));
fwrite(SPEX,uint8([70; 48; 44; StepsToMove270M; 13])); %Moves the grating!
pause(120)
sprintf('The Monochromator has moved to %d nm.', Current_270M_Wavelength)
Current_270M_Wavelength = Emission_Wavelength;

%% Excitation Scans start here
N=1;
for N=1:1:NumRuns
    start = 1;
    SatCheck = 0;
 for Wavelength=MinWavelength:Inc_Wavelength:MaxWavelength
    StepsToMoveInc = double(transpose(int2str(((Inc_Wavelength)*50))));
    fwrite(PMT,uint8(hex2dec(['53'; '0D'])));
    pause(1.1)
    sol=fread(PMT,PMT.BytesAvailable);
    if sol(1)==255 
        SatCheck=SatCheck+1;
    end
    if SatCheck==3
        fwrite(PMT,uint8([86; bitshift(0,-8); bitand(0,255); 13]))
        msgbox('Program shutdown due to too much light! Reduce Integration time or light!')
        error('Too much light! 3 saturation readings!')
    end
    counts(N,start)= sol(1)*16777216 +sol(2)*65536 + sol(3)*256 + sol(4);
    fprintf(fid1, '%0.3f \t %0.3f\r\n',[Wavelength; counts(start)]);
    Wavelength_Plot(N,start)=Wavelength;
    
    plot(Wavelength_Plot(1,start), sum(counts(1:1:N,start),1), '*r', 'LineWidth', 1.5);
    xlabel('\bf Wavelength \rm (nm)')
    ylabel('\bf Counts')
    start=start+1;
    drawnow;

    fwrite(SPEX,uint8([70; 48; 44; StepsToMoveInc; 13])); %Moves the grating!
    char(fread(SPEX,1));
    pause(PauseTime);
    fwrite(SPEX, uint8([32]))
    char(fread(SPEX,1))
    SPEX_Wavelength = Wavelength;
    
 end

fname = sprintf('Experiment %d Run %d of %d.mat', CurrentNum, N, NumRuns)
ExcSavedWavelength = Wavelength_Plot(N,:);
SavedCounts = counts(N,:);
save(fname, 'ExcSavedWavelength', 'SavedCounts') 

N=N+1;
StepsToMoveNew = double(transpose(int2str(((MinWavelength-MaxWavelength)*50-50*Inc_Wavelength))));
fwrite(SPEX,uint8([70; 48; 44; StepsToMoveNew; 13]));
char(fread(SPEX,1));
ptime(120)
clear StepsToMoveInc;
cla
end

hold off
% figure
% plot(Wavelength_Plot(1,:), sum(counts,1), '-r', 'LineWidth', 1.5)
% xlabel('\bf Wavelength \rm (nm)')
% ylabel('\bf Counts')
% ax=gca;
% set(ax, 'XMinorTick','on', 'YMinorTick','on')
% set(gca, 'YTickLabel', num2str(transpose(get(gca,'YTick'))))
% saveas(gcf,'Varying_Wavelength', 'tiffn')
% saveas(gcf,'Varying_Wavelength', 'jpg')
% saveas(gcf,'Varying_Wavelength', 'fig')
% saveas(gcf,'Varying_Wavelength', 'm')

