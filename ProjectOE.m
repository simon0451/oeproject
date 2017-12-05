%%%% Costal Engineering Project %%%%%

%Deliverable List
%Dune Profile for the System pre and post storm

dat=csvread('Bathymetry.csv',1,0);
xcord=dat(:,1);
zmarch=dat(:,3);
zjune=dat(:,2);

figure(1)
plot(xcord,zmarch,xcord,zjune)
ylabel('Depth [m]')
xlabel('Offshore Distance [m]')
legend('Bathymetry (March 16, 2012)','----------------(June 12, 2012)')


% %2a Wind Conditions

filename = 'E:\OE\Project\WindSpeed2012.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%s%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
date_time = dataArray{:, 1};
WindSpeedms = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

figure(2)
[U,X]=hist(WindSpeedms,0:1:50);
bar(X,U/length(WindSpeedms)*100)
xlim([0 25])
xlabel('Wind Speed [m/s]')
ylabel('Occurence Percentage')

filename = 'E:\OE\Project\WindDirection2012.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%s%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
date_time = dataArray{:, 1};
WindDirectiondegree = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

figure(3)
[Deg,X1]=hist(WindDirectiondegree,0:10:360);
bar(X1,Deg/length(WindDirectiondegree)*100)

%storm wind data

filename = 'E:\OE\Project\WindConditions.csv';
delimiter = ',';
formatSpec = '%*s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
rawData = dataArray{1};
for row=1:size(rawData, 1);
    % Create a regular expression to detect and remove non-numeric prefixes and
    % suffixes.
    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
    try
        result = regexp(rawData{row}, regexstr, 'names');
        numbers = result.numbers;
        
        % Detected commas in non-thousand locations.
        invalidThousandsSeparator = false;
        if any(numbers==',');
            thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
            if isempty(regexp(thousandsRegExp, ',', 'once'));
                numbers = NaN;
                invalidThousandsSeparator = true;
            end
        end
        % Convert numeric strings to numbers.
        if ~invalidThousandsSeparator;
            numbers = textscan(strrep(numbers, ',', ''), '%f');
            numericData(row, 1) = numbers{1};
            raw{row, 1} = numbers{1};
        end
    catch me
    end
end
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells
WindSpeedms1 = cell2mat(raw(:, 1));
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;

figure(4)
[U1,X2]=hist(WindSpeedms1,0:1:50);
bar(X2,U1/length(WindSpeedms1)*100)
xlim([0 15])
xlabel('Wind Speed [m/s]')
ylabel('Occurence Percentage')

filename = 'E:\OE\Project\WindDirectionMarch15June15.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%s%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
date_time1 = dataArray{:, 1};
WindDirectiondegree1 = dataArray{:, 2};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

figure(5)
[U2,X3]=hist(WindDirectiondegree1,0:10:360);
bar(X3,U2/length(WindDirectiondegree1))
xlabel('Wind Direction [^o]')
ylabel('Frequency of Occurence')

%2b Annual Wave Conditions
Jan=importdata('January2012.txt');
HJan=Jan.data(:,6);
TJan=Jan.data(:,9);
Feb=importdata('February2012.txt');
HFeb=Feb.data(:,6);
TFeb=Feb.data(:,9);
Mar=importdata('March2012.txt');
HMar=Mar.data(:,6);
TMar=Mar.data(:,9);
May=importdata('May2012.txt');
HMay=May.data(:,6);
TMay=May.data(:,9);
Jun=importdata('June2012.txt');
HJun=Jun.data(:,6);
TJun=Jun.data(:,9);
Dec=importdata('December2012.txt');
HDec=Dec.data(:,6);
TDec=Dec.data(:,9);

H=[HJan; HFeb; HMar; HMay; HJun; HDec];
T=[TJan; TFeb; TMar; TMay; TJun; TDec];
H(isnan(H(:,1)),:) = [];
T(isnan(T(:,1)),:) = [];

figure(6)
histnorm(H)
xlabel('Normalized Wave Height [m]')
ylabel('Occurence Percentage [%]')

figure(7)
histnorm(T)
xlabel('Normalized Wave Period [s]')
ylabel('Occurence Probability')




% %% Sediment
% clear all
% close all
% %% data from the frf army study at duck beach in 1994
% %% this is at a cross-shore distance of 69.4m, which is roughly where the frf army study definined the onshore beach right before the waterline
% seive_size=[1.414 1.189 1 .841 .707 .595 .5 .42 .354 .297 .25 .21 .177 .149 .125 .105]; % in mm 
% location_percent_in_each_seive=[0 .2 .29 .84 3.73 11.39 30.5 15.42 10.66 7.52 9.62 4.67 3.24 1.43 .34 .15]; % percentage of the total weight at each seive size
% % percentfiner=[0 .2 .49 1.33 5.06 16.45 46.95 62.37 73.03 80.55 90.17 94.84 98.08 99.51 99.85 100]; % percent finer
% percentfiner=[100 99.85 99.51 98.08 94.84 90.17 80.55 73.03 62.37 46.95 16.45 5.06 1.33 .49 .2 0];
% figure
% plot(seive_size,percentfiner)
% ylabel('% Finer (cumulative)')
% xlabel('Seive Size [mm]')
% title('CDF')
% grid on
% figure
% plot(seive_size,location_percent_in_each_seive)
% ylabel('% of total weight')
% xlabel('Seive Size [mm]')
% grid on
% %%
% d50=.308275;
% g=9.81; 
% u_mean=mean(us);
% u_given=(us);
% u_prime=u_given-u_mean;
% v_mean=mean(vs);
% v_given=(vs);
% v_prime=v_given-v_mean;
% urms=sqrt(mean(u_prime.^2));
% vrms=sqrt(mean(v_prime.^2));
% u_o=sqrt(2)*sqrt(urms^2+vrms^2);
% sigma_press=std(ps);
% p_mean=mean(ps);
% ds=d_50/1000;
% Cd=.5;
% k=.4;
% S=2.54;
% theta_crit=.04; %% might need to be changed
% kw=.0046/((S-1)*g);
% kc=.0053/((S-1)*g);
% eb=.135;
% phaseR=33; 
% es=.015;
% Hs=4*sigma_press; %%% this might need to be set to the mean Hs from the data we have, unsure though
% T=5; %%% this needs to be set to the mean(period) from the data we have
% d_o=Hs*(T/(2*pi))*sqrt(g/p_mean);
% 
% fs=exp(((5.213*((2.5*ds)/(d_o/2))^.194)-5.977));
% tau=(.5)*fs*rho_water*(((d_o/2)*((2*pi)/T))^2);
% Ustar=sqrt(tau/rho_water);
% ws=sqrt((4/3)*(g*ds/Cd)*(S-1));
% Ro=ws/(k*Ustar);
% theta=(tau/(rho_water*(S-1)*g*ds));
% q_b_mpm=(8*((theta-theta_crit)^1.5)*ds*sqrt(g*ds*(S-1)))*(3600*24); % bedload  [in days]
% q_b_eng=abs((kw*((eb*(u_o^3))/tan(phaseR))+kc*((eb*((sqrt(u_mean^2 + u_o^2))^2)*u_mean)/tan(phaseR)))*(3600*24));
% q_s_cem=((7*(theta-theta_crit)*sqrt(theta)))*ws*ds*ds*sqrt(g*ds*(S-1))*(3600*24);
% q_s_eng=(kw*((es*(u_o^4))/ws)+kc*((es*((sqrt(u_mean^2 + u_o^2))^3)*u_mean)
