clc ; clear ; close all
% This code takes the raw data from the GPPH experiment and plots just the 
% raw data with no calibration, adjusting, or filtering

% DataPlotter_v71

%% Collecting Raw Data
load('runID_parametricstudy.mat')
PressLabels = {'P11','P12','P21','P22'} ;
Ln = 3 ;
font = 20 ;
yellow = [0.831, 0.686, 0.216] ; 
orange = [0.8, 0.333, 0] ; 
maroon = [0.5 0 0] ; 
purple = [0.7, 0.3, 0.8] ; 
teal = [0.784,0.635,0.784] ;
Colors = [yellow; orange; maroon; purple; teal] ;
opacity = 1 ;
Lines = ["-", "-.", ":"] ;

%% Load directory
% cd runs_parametricstudy
cd 'G:\Shared drives\Hydroelasticity Laboratory - Gilbert Research Group\Personnel\Z_Past Students\Past Students with Papers to Finish\Van Erem\Testing\Parametric Study\runs_parametricstudy'
lastrun = 162 ;

for ii = 2
    run = Runs{ii} ;
    preimpacttime = 0 ;
    endimpacttime = 0.4 ;
    endpressurewindow = 0.15 ;
     if Trim(ii) == 0 % deg
         keelentry = 0.241 ; % m
         LA0startpos = 0 ; % m
         LA1startpos = 0 ; % m
         Sw(ii) = 0.419 ; % m2
     elseif Trim(ii) == -5 % deg
         keelentry = 0.241 - 0.0555 ; % m
         LA0startpos = 0.0555 ; % m
         LA1startpos = 0 ; % m
         Sw(ii) = 0.323 ; % m2
     elseif Trim(ii) == -10 % deg
         keelentry = 0.241 - 0.1128 ; % m
         LA0startpos = 0.1128 ; % m
         LA1startpos = 0 ; % m
         Sw(ii) = 0.196 ; % m2
     elseif Trim(ii) == -15 % deg
         keelentry = 0.241 - 0.170148; % m
         LA0startpos = 0.170148; % m
         LA1startpos = 0; % m
         Sw(ii) = 0.067; % m2
     end 

for i = 1 %:length(run)

    [LATrial, averagetrim] = GetLAData(run(i), keelentry, LA0startpos, LA1startpos, Heave(ii));
    [ForceTrial] = GetForceData(run(i), keelentry, LA0startpos, LA1startpos, Trim(ii), Heave(ii), Surge(ii)) ;
    [PressureTrial] = GetPressureData(run(i), keelentry, LA0startpos, LA1startpos, Heave(ii));
    [PotandIncTrial] = GetPotandIncData(run(i), keelentry, LA0startpos, LA1startpos, Heave(ii));
    [CarTrial] = GetCarData(run(i), keelentry, LA0startpos, LA1startpos, Heave(ii));

    fprintf('Running test condition: %d\n', ii)
    fprintf('Showing data for run # %d\n', run(i))

%% Plots
% Linear Actuator Data
    % Position
        figure(1)
        subplot(2,2,1)
        hold on
        plot(LATrial(:, 1), LATrial(:, 2:5),'LineWidth',Ln)
        
        xlabel("Time [s]", 'FontSize', font)
        ylabel("Position [m]", 'FontSize', font)
        
        set(gca,'FontSize', font);
        grid on

    % Velocity
        figure (1)
        subplot(2,2,2)
        hold on
        plot(LATrial(:, 1), LATrial(:, 6:9),'LineWidth', Ln)
        
        xlabel("Time [s]", 'FontSize',font)
        ylabel("Velocity [m/s]", 'FontSize',font)
        
        set(gca,'FontSize',font);
        grid on

    % Carriage Position
        figure(1)
        subplot(2,2,3)
        hold on 
        plot(CarTrial(:, 1), CarTrial(:, 2),'LineWidth',Ln)
       
        xlabel("Time [s]", 'FontSize',font)
        ylabel("Tank Position [m]", 'FontSize',font)
        
        set(gca,'FontSize',font);
        grid on

    % Carriage Velocity
        figure(1), subplot(2,2,4)
        hold on 
        plot(CarTrial(:, 1), CarTrial(:, 3),'LineWidth',Ln)
        %title("Carriage Velocity Data", 'FontSize',font)
        xlabel("Time [s]", 'FontSize',font)
        ylabel("Velocity [m/s]", 'FontSize',font)
        ylim([0 8])
        
        set(gca,'FontSize',font);
        grid on
        hold off

% Steady State Windows    
      figure(2)
      hold on
      plot(CarTrial(:, 1), CarTrial(:, 3)/Surge(ii), 'LineWidth', Ln)
      plot(LATrial(:, 1), LATrial(:, 8:9)*(Surge(ii)/Heave(ii))/Surge(ii), 'LineWidth', Ln)
      title('Steady State Window Comparison', 'FontSize', font)
      xlabel('Time [s]', 'FontSize', font)
      ylabel('Velocity', 'FontSize', font)
      legend('Carriage', 'Actuator 0', 'Actuator 1', 'FontSize', font)
      set(gca, 'FontSize', font) ;
      grid on
      hold off

% Force Data
        figure(9), subplot(3, 1, 1)
        hold on
        plot(ForceTrial(:, 1), ForceTrial(:, 2), 'LineWidth', Ln)
        xlabel('Time [s]','FontSize', font)
        ylabel('Force [N]', 'FontSize', font)
        legend('x-axis', 'location', 'best')
        
        set(gca, 'FontSize', font) ;
        grid on

        figure(9), subplot(3, 1, 2)
        plot(ForceTrial(:, 1), ForceTrial(:, 3), 'LineWidth', Ln)
        xlabel('Time [s]','FontSize', font)
        ylabel('Force [N]', 'FontSize', font)
        legend('y-axis', 'location', 'best')
        
        set(gca, 'FontSize', font) ;
        grid on

        figure(9), subplot(3, 1, 3)
        plot(ForceTrial(:, 1), ForceTrial(:, 4), 'LineWidth', Ln)
        xlabel('Time [s]','FontSize', font)
        ylabel('Force [N]', 'FontSize', font)
        legend('z-axis', 'location', 'best')
        
        set(gca, 'FontSize', font) ;
        grid on
        hold off
        figure(8)
        plot(ForceTrial(:, 1), ForceTrial(:, 2:4), 'LineWidth', Ln) % forcetime, fxforce, fyforce, fzforce
        xlabel('Time [s]', 'FontSize', font)
        ylabel('Force [N]', 'FontSize', font)
        
        legend('x-axis', 'y-axis', 'z-axis', 'location', 'best') 
        set(gca, 'FontSize', font) ;
        grid on

    % Total Force Data
        figure(3), subplot(3,1,3)
        hold on 
        plot(ForceTrial(:, 1), ForceTrial(:, 6),'LineWidth',Ln) % forcetime, totalzforce
        
        xlabel("Time [s]", 'FontSize',font)
        ylabel("Total Force [N]", 'FontSize',font)
        
        set(gca,'FontSize',font);
        grid on

    % Midship Force Data
        figure(3), subplot(3,1,1)
        hold on 
        plot(ForceTrial(:, 1), ForceTrial(:, 4),'LineWidth',Ln) % forcetime, fzforce
        xlabel("Time [s]", 'FontSize',font)
        ylabel("Midship Force [N]", 'FontSize',font)
        
        set(gca,'FontSize',font);
        grid on

    % Stern Force Data
        figure(3),subplot(3,1,2)
        hold on
        plot(ForceTrial(:, 1), ForceTrial(:, 5),'LineWidth',Ln) % forcetime, ozforce
        xlabel("Time [s]", 'FontSize',font)
        ylabel("Stern Force [N]", 'FontSize',font)
        
        set(gca,'FontSize',font);
        grid on
        hold off


% Pressure Data
    figure(4)
    title(['Pressure Data Cond ', num2str(Cond(ii))], 'FontSize',font)
    opacity = 1;
        % Individual Pressure Sensor Data
            for j = 2:5 % p11pressure, p12pressure, p21pressure, p22pressure
                figure(4), 
                subplot(3, 2, j-1)
                hold on 
                
                individualpress_plot = plot(PressureTrial(:, 1), PressureTrial(:, j)) ;
                set(individualpress_plot, 'LineWidth' ,Ln, 'Color', [Colors(j-1, :), opacity])
                xlabel("Time [s]", 'FontSize',font)
                
                ylabel(PressLabels{j - 1}, 'FontSize', font)
                xlim([-0.25 0.25])
                set(gca,'FontSize', font) ;
                grid on
            end
    % Total Pressure Data
        figure(4), subplot(3, 2, [5,6])
        hold on
        markers = {'^', 'o', 's', 'd'};
        for k = 1:length(PressLabels)
            
            pressplots = plot(PressureTrial(:, 1), PressureTrial(:, k + 1) ) ;
            set(pressplots, 'LineWidth',Ln, 'Color', [Colors(k, :) opacity]);
            
        end
        hold off
        xlabel("Time [s]", 'FontSize',font)
        ylabel('nondim P', 'FontSize', font)
        xlim([-0.25 0.25])
        legend('P11', 'P12', 'P21', 'P22', 'Location', 'eastoutside')
        set(gca,'FontSize',font);
        grid on

% Potentiometer and Inclinometer Data
    % String Pot
        figure(6)
        subplot(2,1,1)
        hold on
        plot(PotandIncTrial(:, 1), PotandIncTrial(:, 2),'LineWidth',Ln) % PotandInctime, Pot
        plot(LATrial(:,1), LATrial(:, 4),'LineWidth',Ln) % LAtime, LAmeaspos0
        % title('String Potentiometer Data')
        xlabel('Time [s]', 'FontSize',font)
        ylabel('Position (String Pot) [m]', 'FontSize',font)
        set(gca,'FontSize',font);
        grid on
        hold off
    % Inclinometer
        figure(6)
        subplot(2,1,2)
        hold on
        plot(PotandIncTrial(:, 1), PotandIncTrial(:, 3),'LineWidth',Ln) % PotandInctime
        plot(LATrial(:, 1), LATrial(:, 12),'LineWidth',Ln) % trim 
        xlabel('Time [s]', 'FontSize',font)
        ylabel('Trim Angle [deg]', 'FontSize',font)
        legend('Inclinometer','Measured LA', 'Position',[0.607142857142857 0.414695943674116 0.285714285714285 0.0939177101967799],'FontSize',font)
        set(gca,'FontSize',font);
        grid on
        hold off        
end
end

%% FUNCTION FILES %%
function [func, averagetrim] = GetLAData(Run, keelentry, LA0startpos, LA1startpos, Heave)
    if exist(Run + "LAData.csv ", 'file') ~= 0

    [LA] = readtable(Run + "LAData.csv") ;
        if size(LA, 1) ~= 0
    endrun = 300 ;
    %endrun = find(LA.Var9(11:end, :), 1, 'last'); % currently only accounts for actuator 0, will not work for async runs
    LAcycletime = LA.Var2(11:endrun, :);            % Cycle Time
    LAexppos0 = keelentry - LA.Var5(11:endrun,:);   % Expected Position 0
    LAmeaspos0 = keelentry - LA.Var6(11:endrun,:);  % Measured Position 0
    LAexppos1 = keelentry - LA.Var7(11:endrun,:);   % Expected Position 1
    LAmeaspos1 = keelentry - LA.Var8(11:endrun,:);  % Measured Position 1
    LAexpvel0 = LA.Var9(11:endrun, :);              % Expected Velocity 0
    LAmeasvel0 = LA.Var10(11:endrun, :);            % Measured Velocity 0
    LAexpvel1 = LA.Var11(11:endrun, :);             % Expected Velocity 1
    LAmeasvel1 = LA.Var12(11:endrun, :);            % Measured Velocity 1
  
    % Calculating trim 
    distbetweenLA = 0.64; % m
    distfromhome0 = LA0startpos + LAmeaspos0; 
    distfromhome1 = LA1startpos + LAmeaspos1;
    trim = atand((distfromhome1 - distfromhome0) ./ distbetweenLA);
    averagetrim = mean(trim);

    func = zeros(length(LAcycletime), 18);
    func(1:length(LAcycletime), 1) = LAcycletime;
    func(1:length(LAexppos0), 2) = LAexppos0;
    func(1:length(LAexppos1), 3) = LAexppos1;
    func(1:length(LAmeaspos0), 4) = LAmeaspos0;
    func(1:length(LAmeaspos1), 5) = LAmeaspos1;
    func(1:length(LAexpvel0), 6) = LAexpvel0;
    func(1:length(LAexpvel1), 7) = LAexpvel1;
    func(1:length(LAmeasvel0), 8) = LAmeasvel0;
    func(1:length(LAmeasvel1), 9) = LAmeasvel1;
    %func(1, 10) = impactcycletime;
    %func(1, 11) = endsteadystatetime;
    func(1:length(trim), 12) = trim;
    %func(1:length(LAmeasacc0), 13) = LAmeasacc0;
    %func(1:length(LAmeasacc1), 14) = LAmeasacc1;
    func(1, 15) = averagetrim;
    %func(1:length(LAtime_isolated), 16) = LAtime_isolated;
    %func(1:length(LAtime_isolated), 17) = LAmeaspos0_isolated;
    %func(1,18) = submergence_peak_time;
        else
            func = zeros(290, 18) ;
            averagetrim = zeros(1, 1) ;
        end

    else
        func = zeros(290, 18) ;
        averagetrim = zeros(1, 1) ;
    end
end

% Carriage Data
function func = GetCarData(Run, keelentry, LA0startpos, LA1startpos, Heave)
    if exist(Run + "Car.csv", 'file') ~= 0
    [car] = readtable(Run + "Car.csv");
        if size(car, 1) ~= 0
        if Run ~= 71
    carstart = find(car.ActualVel(:, :) > 0.01, 1, 'first');
    carend = carstart + 1000;
    
    % find LA cycle time at impact event
    LAData = GetLAData(Run, keelentry, LA0startpos, LA1startpos, Heave); 
        if ~ismember(Run, [2 99 30 32 34 66 68 69])
            impactcycletime = LAData(1, 10);
        elseif Run == 2
            impactcycletime = 5.6050e03 ;
        elseif Run == 99
            impactcycletime = 1.0014e03 ;
        elseif Run == 30
            impactcycletime = 2.9725e03 ;
        elseif Run == 32
            impactcycletime = 3.9725e03 ;
        elseif Run == 34
            impactcycletime = 3.6785e03 ;
        elseif Run == 66
            impactcycletime = 7.6230e03 ;
        elseif Run == 68
            impactcycletime = 9.2357e03 ;
        elseif Run == 69
            impactcycletime = 9.5495e03 ;
        end
    
    cartime = car.CycleTime(carstart:carend, :) - impactcycletime;
    carpos = car.ActualPos(carstart:carend, :);
    carvel = car.ActualVel(carstart:carend, :);
    
    func = zeros(length(cartime), 3);
    func(1:length(cartime), 1) = cartime;
    func(1:length(carpos), 2) = carpos;
    func(1:length(carvel), 3) = carvel;
        else
            func = zeros(1001, 3) ;
        end
        else
        func = zeros(1001, 3);    
        end
    else
        func = zeros(1001, 3);    
    end
end

% Force Data
function [func, Vn, nondimfactor_force] = GetForceData(Run, keelentry, LA0startpos, LA1startpos, trim, Heave, Surge)
    if exist(Run + "Force.csv", 'file') ~= 0
    [force] = readtable(Run + "Force.csv");
        if size(force, 1) ~= 0
    % find LA cycle time at impact event
    LAData = GetLAData(Run, keelentry, LA0startpos, LA1startpos, Heave);
        if ~ismember(Run, [99 30 32 66 68 69 71])
            impactcycletime = LAData(1, 10);
            endsteadystatetime = LAData(1, 11); 
        elseif Run == 99
            impactcycletime = 1.0014e03 ;
            endsteadystatetime = 0.0832 ;
        elseif Run == 66
            impactcycletime = 7.6230e03 ;
            endsteadystatetime = 0.2880 ;
        elseif Run == 68
            impactcycletime = 9.2357e03 ;
            endsteadystatetime = 0.3000 ;
        end
    
    forcecycletimestart = force.XValuesGuaranteedValidOnlyForCDAQ1Mod2_ai0(11, :);
    forcetimesync = forcecycletimestart - impactcycletime;
    forcetime = force.Notes(11:end, :) + forcetimesync;
    % futek forces
        fxvoltage = force.Var3(11:end,:) ; 
        fyvoltage = force.Var4(11:end, :) ;
        fzvoltage = force.Var5(11:end, :) ;
    % omega force
        ozvoltage = force.Var6(11:end, :) ;

    totalzforce = fzvoltage + ozvoltage ;
    fxmax = max(fxvoltage) ;
    fymax = max(fyvoltage) ;
    fzmax = max(fzvoltage) ;
    ozmax = max(ozvoltage) ;

    func(1:length(forcetime), 1) = forcetime;
    func(1:length(forcetime), 2) = fxvoltage;
    func(1:length(forcetime), 3) = fyvoltage;
    func(1:length(forcetime), 4) = fzvoltage;
    func(1:length(forcetime), 5) = ozvoltage;
    func(1:length(forcetime), 6) = totalzforce;
    func(1, 7) = fxmax; 
    func(1, 8) = fymax; 
    func(1, 9) = fzmax; 
    func(1, 10) = ozmax;
    else
            func = zeros(200000, 31);
            Vn = zeros(1, 1) ;
            nondimfactor_force = zeros(1, 1) ;
        end
    else
        func = zeros(200000, 31);    
    end
end

% Pressure Data
function [func, p11pressuretime_rise_revised] = GetPressureData(Run, keelentry, LA0startpos, LA1startpos, Heave)
    if exist(Run + "ACandPressure.csv", 'file') ~= 0
    [ACandPressure] = readtable(Run + "ACandPressure.csv");
        if size(ACandPressure, 1) ~= 0
    % find LA cycle time at impact event
    LAData = GetLAData(Run, keelentry, LA0startpos, LA1startpos, Heave);
    if ~ismember(Run, [99 30 32 66 68 69 71])
        impactcycletime = LAData(1, 10);
        endsteadystatetime = LAData(1, 11); 
    elseif Run == 99
        impactcycletime = 1.0014e03 ;
        endsteadystatetime = 0.0832 ;
    %elseif Run == 68
        %impactcycletime = 7.6230e03 ;
        %endsteadystatetime = 0.2880 ;
    end
    pressurecycletimestart = ACandPressure.XValuesGuaranteedValidOnlyForCDAQ1Mod3_ai0(11, :);
    pressuretimesync = pressurecycletimestart - impactcycletime;
    pressuretime = ACandPressure.Notes(11:end, :) + pressuretimesync;

    p11voltage = ACandPressure.Var6(11:end, :) ;
    p12voltage = ACandPressure.Var7(11:end, :) ;
    p21voltage = ACandPressure.Var8(11:end, :) ;
    p22voltage = ACandPressure.Var9(11:end, :) ;

    func(1:length(pressuretime), 1) = pressuretime;
    func(1:length(pressuretime), 2) = p11voltage ;
    func(1:length(pressuretime), 3) = p12voltage ;
    func(1:length(pressuretime), 4) = p21voltage ;
    func(1:length(pressuretime), 5) = p22voltage ;
    else
            func = zeros(200000, 15);
            p11pressuretime_rise_revised = zeros(1, 1) ;
        end
    else
        func = zeros(200000, 15);
    end
end

% Potentiometer and Inclinometer Data
function func = GetPotandIncData(Run, keelentry, LA0startpos, LA1startpos, Heave)
    if exist(Run + "PotandInc.csv", 'file') ~= 0 
    [PotandInc] = readtable(Run + "PotandInc.csv");
        if size(PotandInc, 1) ~= 0
    % find LA cycle time at impact event
    LAData = GetLAData(Run, keelentry, LA0startpos, LA1startpos, Heave); 
        if ~ismember(Run, [99 30 32 66 68 69 71])
            impactcycletime = LAData(1, 10);
            endsteadystatetime = LAData(1,11);
        elseif Run == 99
            impactcycletime = 1.0014e03 ;
            endsteadystatetime = 0.0832 ;
        elseif Run == 66
            impactcycletime = 7.6230e03 ;
            endsteadystatetime = 0.2880 ;
        end

    PotandInccycletimestart = PotandInc.XValuesGuaranteedValidOnlyForCDAQ1Mod1_ai0(11, :);
    PotandInctimesync = PotandInccycletimestart - impactcycletime;
    PotandInctime = PotandInc.Notes(11:end, :) + PotandInctimesync;
    
    Potvoltage = PotandInc.Var3(11:end, :) ; 
    Incvoltage = PotandInc.Var4(11:end, :) ;

    averagetrim = mean(Incvoltage);
    
    func = zeros(length(PotandInctime), 4);
    func(1:length(PotandInctime), 1) = PotandInctime;
    func(1:length(Potvoltage), 2) = Potvoltage;
    func(1:length(Incvoltage), 3) = Incvoltage;
    func(1, 4) = averagetrim;
        else
            func = zeros(10000, 4) ;
        end
    else
        func = zeros(10000, 4) ;
    end
end