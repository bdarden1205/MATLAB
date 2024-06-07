clc ; clear ; close all

% This code collects the raw data of the carriage, linear actuators, 
% load cells, accelerometers, pressure sensors, string potentiometer, and inclinometer
% and outputs the data in a graphical format. This code applies the conversion  
% factors of the sensors, offsets the data from the ZEROS, and uses the cycle time 
% to synchonize the time of the data

% DataPlotter_v7 (Darden version)

%% Collecting Raw Data
load('runID_parametricstudy.mat')
PressLabels = {'P11','P12','P21','P22'} ;
Ln = 3 ;
font = 20 ;
paperpos = [0 0 3 2.62] ;
papersize = [3 2.62] ;
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
    [ForceTrial, Vn(ii), nondimfactor_force(ii)] = GetForceData(run(i), keelentry, LA0startpos, LA1startpos, Trim(ii), Heave(ii), Surge(ii));
    [HydrostaticForceTrial] = GetHydroStaticForce(run(i), keelentry, LA0startpos, LA1startpos, Heave(ii));
    [PressureTrial, p11pressuretime_rise_revised(ii,i)] = GetPressureData(run(i), keelentry, LA0startpos, LA1startpos, Heave(ii));
    % [AccelerometerTrial] = GetAccelerometerData(run(i), keelentry, LA0startpos, LA1startpos);
    [PotandIncTrial] = GetPotandIncData(run(i), keelentry, LA0startpos, LA1startpos, Heave(ii));
    [CarTrial] = GetCarData(run(i), keelentry, LA0startpos, LA1startpos, Heave(ii));

fprintf('Running test condition: %d\n', ii)
fprintf('Showing data for run # %d\n', run(i))
% Z force 
    peaktime = ForceTrial(1, 12);
    peakvalue = ForceTrial(1, 11);
    peakvalue_nondim = ForceTrial(1, 26); %22 is umd, 23 is nominal disp, 24 is wetted surface, 25 is specific disp, 26 is test
    %fprintf('The peak force time is: %d\n', peaktime);
    %fprintf('The peak force value is: %d\n', peakvalue);
    %fprintf('For case: %d\n', ii);
    %fprintf('On run: %d\n', run(i));
    peaktotalzforce(ii, i) = peakvalue;
    peaktotalzforce_nondim(ii, i) = peakvalue_nondim;
    
    totalzforce(:, i) = ForceTrial(1:2000, 31);
    totalzforcetime(i, :) = ForceTrial(1:2000, 27);
    % totalzforce(ii, :) = ForceTrial(:, 6);
        
    totalzforce_nondim_met2(i, :) = ForceTrial(:, 17) ; % UMD Nondim
    totalzforce_nondim_met3(i, :) = ForceTrial(:, 18) ; % Nominal Disp
    totalzforce_nondim_met4(i, :) = ForceTrial(:, 19) ; % Wetted Surface 
    totalzforce_nondim_met5(i, :) = ForceTrial(:, 20) ; % Specific Disp
    totalzforce_nondim_met6(i, :) = ForceTrial(:, 21) ; % Test

% Pressures
    rho = 1000 ; % density of water
    VT(ii) = sqrt(Surge(ii)^2 + Heave(ii)^2) ; % Resultant velocity
    nondimfactor_pressure(ii) = rho * VT(ii)^2 ; % dynamic pressure
    
    nondimfactor_force(ii) = rho * Sw(ii) * VT(ii)^2 ; % factor used to nondimensionalize force
    
    peakpressure_p11(ii, i) = PressureTrial(1,6) ; % p11max_mag 
    peakpressure_p12(ii, i) = PressureTrial(2,6) ; % p12max_mag
    peakpressure_p21(ii, i) = PressureTrial(3,6) ; % p21max_mag
    peakpressure_p22(ii, i) = PressureTrial(4,6) ; % p22max_mag
    
    peakpressuretime_p11(ii, i) = PressureTrial(6,6) ; % p11max_time
    peakpressuretime_p12(ii, i) = PressureTrial(7,6) ; % p12max_time
    
    maxpressure_mags(ii, :) = [PressureTrial(1,6), PressureTrial(2,6), PressureTrial(3,6), PressureTrial(4,6)]; % changes each iteration, does not store values between trials
    maxpressure_times(ii, :) = [PressureTrial(6,6), PressureTrial(7,6), PressureTrial(8,6), PressureTrial(9,6)]; % changes each iteration, does not store values between trials

    peakpressure_p11_nondim(ii, i) = 1000*PressureTrial(1,6)./nondimfactor_pressure(ii); 
    peakpressure_p12_nondim(ii, i) = 1000*PressureTrial(2,6)./nondimfactor_pressure(ii);
    peakpressure_p21_nondim(ii, i) = 1000*PressureTrial(3,6)./nondimfactor_pressure(ii);
    peakpressure_p22_nondim(ii, i) = 1000*PressureTrial(4,6)./nondimfactor_pressure(ii);
    maxpressure_mags_nondim(ii, :) = [peakpressure_p11_nondim(ii, i), peakpressure_p12_nondim(ii, i), peakpressure_p21_nondim(ii, i), peakpressure_p22_nondim(ii, i)];
    
    p11pressure(i, :) = PressureTrial(1:20000, 8)  ; % p11pressure_isolated
    p12pressure(i, :) = PressureTrial(1:20000, 9)  ; % p12pressure_isolated
    p21pressure(i, :) = PressureTrial(1:20000, 10) ; % p21pressure_isolated
    p22pressure(i, :) = PressureTrial(1:20000, 11) ; % p22pressure_isolated
    
    p11pressuretime_revised(i,:) = PressureTrial(1:20000, 12) ; % p11pressuretime_revised
    p11pressure_revised(i,:) = PressureTrial(1:20000, 13)     ; % p11pressure_revised
    p12pressuretime_revised(i,:) = PressureTrial(1:20000, 14) ; % p12pressuretime_revised
    p12pressure_revised(i,:) = PressureTrial(1:20000, 15)     ; % p12pressure_revised
    
    [p11pressure_peak_revised(ii,i), p11_peakindex_revised(ii,i)] = max(p11pressure_revised(i,4000:7000));
    [p12pressure_peak_revised(ii,i), p12_peakindex_revised(ii,i)] = max(p12pressure_revised(i,4000:6500));
    p11pressuretime_peak_revised(ii,i) = p11pressuretime_revised(i, p11_peakindex_revised(ii,i)+4000);
    p12pressuretime_peak_revised(ii,i) = p12pressuretime_revised(i, p12_peakindex_revised(ii,i)+4000);
    
    p11pressure_peak_nondim_revised(ii,i) = 1000*p11pressure_peak_revised(ii,i)./nondimfactor_pressure(ii);
    p12pressure_peak_nondim_revised(ii,i) = 1000*p12pressure_peak_revised(ii,i)./nondimfactor_pressure(ii);
    
% Position 
    Pos_measured_LA0(i,:) = LATrial(1:100, 17) ; % LAmeaspos0_isolated 
    endsteadystatetime(i) = LATrial(1, 11) ; % endsteadystatetime

% Peak Submergence 
    submergence_peak_time(ii, i) = LATrial(1,18) ; % submergence_peak_time


%% Plots
% Linear Actuator Data
    % Position
        figure(1)
        subplot(2,2,1)
        hold on
        plot(LATrial(:, 1), LATrial(:, 2:5),'LineWidth',Ln)
        %title("Actuator Position Data", 'FontSize',font)
        xlabel("Time [s]", 'FontSize', font)
        ylabel("Position [m]", 'FontSize', font)
        %xlim([preimpacttime endimpacttime])
        %legend('Expected Pos 0', 'Expected Pos 1','Measured Pos 0','Measured Pos 1', 'FontSize',font)
        set(gca,'FontSize', font);
        grid on

    % Velocity
        figure (1)
        subplot(2,2,2)
        hold on
        plot(LATrial(:, 1), LATrial(:, 6:9),'LineWidth', Ln)
        %title("Actuator Velocity Data", 'FontSize',font)
        xlabel("Time [s]", 'FontSize',font)
        ylabel("Velocity [m/s]", 'FontSize',font)
        %xlim([preimpacttime endimpacttime])
        %legend('Expected Vel 0', 'Expected Vel 1','Measured Vel 0','Measured Vel 1', 'FontSize',font)
        set(gca,'FontSize',font);
        grid on

    % Carriage Position
        figure(1)
        subplot(2,2,3)
        hold on 
        plot(CarTrial(:, 1), CarTrial(:, 2),'LineWidth',Ln)
        %title("Carriage Position Data", 'FontSize',font)
        xlabel("Time [s]", 'FontSize',font)
        ylabel("Tank Position [m]", 'FontSize',font)
        %xlim([preimpacttime endimpacttime])
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
        %xlim([preimpacttime endimpacttime])
        set(gca,'FontSize',font);
        grid on
        hold off

% Steady State Windows
%     figure(2)
%     hold on
%     plot(CarTrial(:, 1), CarTrial(:, 3),'LineWidth',Ln)
%     plot(LATrial(:, 1), LATrial(:, 8),'LineWidth',Ln)
%     plot(LATrial(:, 1), LATrial(:, 9),'LineWidth',Ln)
%     xlabel("Time [s]", 'FontSize',font)
%     ylabel("Velocity [m/s]", 'FontSize',font)
%     %xlim([preimpacttime endimpacttime])
%     legend("Carriage", "Actuator 0", "Actuator 1", 'FontSize',font)
%     set(gca,'FontSize',font);
%     grid on
%     hold off
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
    % Futek Force Data
        figure(8), subplot(2, 1, 2)
        hold on
        plot(ForceTrial(:, 1), ForceTrial(:, 2:4), 'LineWidth', Ln) % forcetime, fxforce, fyforce, fzforce
        xlabel('Time [s]', 'FontSize', font)
        ylabel('Force [N]', 'FontSize', font)
        % xlim([preimpacttime endimpacttime])
        xlim([-2 3])
        legend('x-axis', 'y-axis', 'z-axis', 'location', 'best') 
        set(gca, 'FontSize', font) ;
        grid on
        
        ResultantForce = sqrt(ForceTrial(:, 2).^2 + ForceTrial(:, 3).^2 + ForceTrial(:, 4).^2) ;
        figure(8), subplot(2, 1, 1)
        hold on
        plot(ForceTrial(:, 1), ResultantForce, 'LineWidth', Ln, 'LineStyle', '--') % forcetime, fzforce
        plot(ForceTrial(:, 1), abs(ForceTrial(:, 4)), 'LineWidth', Ln, 'LineStyle', ':')
        xlabel('Time [s]', 'FontSize', font)
        ylabel('Force [N]', 'FontSize', font)
        legend('Resultant', 'z-axis', 'location', 'best')
        xlim([-2 3])
        set(gca, 'FontSize', font) ;
        grid on
        hold off
        hold off

        figure(9), subplot(3, 1, 1)
        hold on
        plot(ForceTrial(:, 1), ForceTrial(:, 2), 'LineWidth', Ln)
        xlabel('Time [s]','FontSize', font)
        ylabel('Force [N]', 'FontSize', font)
        legend('x-axis', 'location', 'best')
        xlim([-2 3])
        set(gca, 'FontSize', font) ;
        grid on

        figure(9), subplot(3, 1, 2)
        plot(ForceTrial(:, 1), ForceTrial(:, 3), 'LineWidth', Ln)
        xlabel('Time [s]','FontSize', font)
        ylabel('Force [N]', 'FontSize', font)
        legend('y-axis', 'location', 'best')
        xlim([-2 3])
        set(gca, 'FontSize', font) ;
        grid on

        figure(9), subplot(3, 1, 3)
        plot(ForceTrial(:, 1), ForceTrial(:, 4), 'LineWidth', Ln)
        xlabel('Time [s]','FontSize', font)
        ylabel('Force [N]', 'FontSize', font)
        legend('z-axis', 'location', 'best')
        xlim([-2 3])
        set(gca, 'FontSize', font) ;
        grid on
        hold off
%         figure(8)
%         plot(ForceTrial(:, 1), ForceTrial(:, 2:4), 'LineWidth', Ln) % forcetime, fxforce, fyforce, fzforce
%         xlabel('Time [s]', 'FontSize', font)
%         ylabel('Force [N]', 'FontSize', font)
%         % xlim([preimpacttime endimpacttime])
%         xlim([-2 3])
%         legend('x-axis', 'y-axis', 'z-axis', 'location', 'best') 
%         set(gca, 'FontSize', font) ;
%         grid on

    % Total Force Data
        figure(3), subplot(3,1,3)
        hold on 
        plot(ForceTrial(:, 1), ForceTrial(:, 6),'LineWidth',Ln) % forcetime, totalzforce
        %plot(ForceTrial(:, 27), ForceTrial(:, 31),'LineWidth',Ln) % forcetime_isolated, totalzforce_isolated
        %scatter(ForceTrial(1, 12), ForceTrial(1, 11), 'Marker', 'd','LineWidth', 5) % totalxpeaktime, totalzpeak
        xlabel("Time [s]", 'FontSize',font)
        ylabel("Total Force [N]", 'FontSize',font)
        % legend('Total', 'Isolated', 'location', 'best')
        % xlim([preimpacttime endimpacttime])
        xlim([-2 3])
        set(gca,'FontSize',font);
        grid on

    % Midship Force Data
        figure(3), subplot(3,1,1)
        hold on 
        plot(ForceTrial(:, 1), ForceTrial(:, 4),'LineWidth',Ln) % forcetime, fzforce
        xlabel("Time [s]", 'FontSize',font)
        ylabel("Midship Force [N]", 'FontSize',font)
        % xlim([preimpacttime endimpacttime])
        xlim([-2 3])
        set(gca,'FontSize',font);
        grid on

    % Stern Force Data
        figure(3),subplot(3,1,2)
        hold on
        plot(ForceTrial(:, 1), ForceTrial(:, 5),'LineWidth',Ln) % forcetime, ozforce
        xlabel("Time [s]", 'FontSize',font)
        ylabel("Stern Force [N]", 'FontSize',font)
        % xlim([preimpacttime endimpacttime])
        xlim([-2 3])
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
                %individualpress_plot = plot(PressureTrial(:, 1), 1000*PressureTrial(:, j)./nondimfactor_pressure);
                individualpress_plot = plot(PressureTrial(:, 1), PressureTrial(:, j)) ;
                set(individualpress_plot, 'LineWidth' ,Ln, 'Color', [Colors(j-1, :), opacity])
                xlabel("Time [s]", 'FontSize',font)
                %ylabel(strcat(PressLabels{j-1},'/rho*V_T^2'), 'FontSize',font)
                %ylabel(strcat('nondim', PressLabels{j - 1}), 'FontSize', font)
                ylabel(PressLabels{j - 1}, 'FontSize', font)
                % xlim([preimpacttime endimpacttime])
                xlim([-0.20 0.20])
                set(gca,'FontSize', font) ;
                grid on
            end
    % Total Pressure Data
        figure(4), subplot(3, 2, [5,6])
        hold on
        markers = {'^', 'o', 's', 'd'};
        for k = 1:length(PressLabels)
            %pressplots = plot(PressureTrial(:,1), 1000*PressureTrial(:, k+1)./nondimfactor_pressure);
            pressplots = plot(PressureTrial(:, 1), PressureTrial(:, k + 1) ) ;
            set(pressplots, 'LineWidth',Ln, 'Color', [Colors(k, :) opacity]);
            % scatter(maxpressure_times(ii, k), maxpressure_mags(ii, k), 'Marker', markers{k}, 'LineWidth', 5);
        end
        hold off
        xlabel("Time [s]", 'FontSize',font)
        %ylabel('P / rho*V_T^2', 'FontSize',font)
        ylabel('nondim P', 'FontSize', font)
        legend('P11', 'P12', 'P21', 'P22', 'Location', 'eastoutside')
        % xlim([preimpacttime endimpacttime])
        xlim([-0.20 0.20])
        set(gca,'FontSize',font);
        grid on

% Acceleration Data
%     figure(5), subplot(4,1,4)
%     hold on
%     plot(LATrial(:,1), LATrial(:, 13),'LineWidth',Ln)
%     plot(LATrial(:,1), LATrial(:, 14),'LineWidth',Ln)
%     xlabel("Time [s]", 'FontSize',font)
%     ylabel(strcat('Acceleration [m/s^2]'), 'FontSize',font)
%     % xlim([preimpacttime endimpacttime])
%     set(gca,'FontSize',font);
%     grid on
%     hold off

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
        % xlim([preimpacttime endimpacttime])
        xlim([-2 3])
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
        % xlim([preimpacttime endimpacttime])
        xlim([-2 3]) 
        legend('Inclinometer','Measured LA', 'Position',[0.607142857142857 0.414695943674116 0.285714285714285 0.0939177101967799],'FontSize',font)
        set(gca,'FontSize',font);
        grid on
        hold off
% Force and Position Data + Nondim Force
    % LA Position
%         figure(7)
%         subplot(2,1,1)
%         hold on 
%         yyaxis right
%         plot(LATrial(:, 1), LATrial(:, 2:5),'LineWidth',Ln,'Color',[0.5
%         0.5 0.5]) % LAtime, LAexppos0, LAexppos1, LAmeaspos0, LAmeaspos1
%         title(strcat('Cond ',num2str(Cond(ii))), 'FontSize',font)
%         xlabel("Time [s]", 'FontSize',font)
%         ylabel("Position [m]", 'FontSize',font)
%         set(gca,'FontSize',font);
%         grid on
%     % Force
%         figure(7)
%         subplot(2,1,1)
%         hold on 
%         yyaxis left
%         plot(ForceTrial(:, 1), ForceTrial(:, 4)+ForceTrial(:,
%         5),'LineWidth',Ln) % forcetime, fzforce + ozforce
%         % scatter(ForceTrial(1, 12), ForceTrial(1, 11), 'Marker', 'd', 'LineWidth', 5)
%         endsteadystate_line = plot([endsteadystatetime(end), endsteadystatetime(end)], [-2, 2], 'LineStyle', ':', 'LineWidth', 2.5, 'Color', yellow,  'Marker', 'none'); % end steady state, LATrial(1,11)
%         % scatter(HydrostaticForceTrial(:, 1), HydrostaticForceTrial(:, 2), 'Marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g')
%         xlabel("Time [s]", 'FontSize',font)
%         ylabel("Total Z-Force [N]", 'FontSize',font)
%         xlim([preimpacttime endimpacttime])
%         set(gca,'FontSize',font);
%         grid on
    % Nondimensionalized Force
%         figure(7), subplot(2,1,2)
%         hold on 
%         % scatter(HydrostaticForceTrial(:, 1), ForceTrial(1:length(HydrostaticForceTrial(:, 1)), 13), 'Marker', 'o', 'LineWidth', 2, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g')
%         % scatter(HydrostaticForceTrial(:, 1), ForceTrial(1:length(HydrostaticForceTrial(:, 1)), 14), 'Marker', 's', 'LineWidth', 2, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')
%         % scatter(HydrostaticForceTrial(:, 1), ForceTrial(1:length(HydrostaticForceTrial(:, 1)), 15), 'Marker', '^', 'LineWidth', 2, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
%         plot(ForceTrial(:, 1), ForceTrial(:, 17), 'b')
%         plot(ForceTrial(:, 1), ForceTrial(:, 18), 'r')
%         plot(ForceTrial(:, 1), ForceTrial(:, 19), '-', 'Color', [1, 0.5, 0]) % orange color
%         plot(ForceTrial(:, 1), ForceTrial(:, 20))
%         plot(ForceTrial(:, 1), ForceTrial(:, 21))
%         % ylim([-1 8])
%         xlabel('Time [s]', 'FontSize',font)
%         ylabel('Nondimensional Force', 'FontSize',font)
%         legend('Max L,B Method', 'Disp Method', 'Wet Surf Method', 'Specific Disp', 'test', 'Location', 'northeast')
%         % legend('Disp(dynamic)', 'UMD', 'Disp', 'Location', 'northwest')
%         set(gca,'FontSize',font);
%         grid on

end 
% put break point here ^^^ when looking at individual trials, comment out when grouping conditions 
% uncomment this when in individual run mode
end 
% uncomment this to display all of the runs together



%% Functions %%

% Linear Actuator Data
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
    
    % Finding impact window, only applicable when motion parameters of LA 0 =
    % LA 1 since only accounting for LA 0
    preimpactindex = find(LAexppos0 < 0.001, 1, 'first') + 1;
    impactcycletime = LAcycletime(preimpactindex);
    
    % Making the slam event the time origin
    LAtime = LAcycletime - impactcycletime;
    
    % Calculating Acceleration
    LAmeasacc0 = zeros(length(LAmeasvel0), 1);
    LAmeasacc1 = zeros(length(LAmeasvel1), 1);
    
    for i = 2 : length(LAmeasvel0)
    LAmeasacc0(i, 1) = (LAmeasvel0(i) - LAmeasvel0(i-1)) / (LAtime(i) - LAtime(i-1));
    end
    
    for i = 2 : length(LAmeasvel1)
    LAmeasacc1(i, 1) = (LAmeasvel1(i) - LAmeasvel1(i-1)) / (LAtime(i) - LAtime(i-1));
    end
    
    % Prescribed Acceleration
    if Heave == 0.25
        prescibedacc = -5;
    else
        prescibedacc = -15;
    end
    
    
    % Cutting off data at the end of steady state
    % endsteadystateindex = find(LAmeasacc0(impactindex:end) <= 0.8*prescribedacceleration, 1); % Multiply by coeff to avoid getting second peak
    % if Run == 3 || Run == 23 || Run == 93 
    %     endsteadystateindex = find(LAmeasacc1(impactindex:end) <= 0.85*prescribedacceleration, 1);
    % end
    %endsteadystateindex = find(LAmeasacc0 <= -4, 1);
    endsteadystateindex = find(LAexpvel0 == max(LAexpvel0), 1, 'last') + 3; % adding 3 to account for PID delay
    % endsteadystateindex = find(LAmeasacc0 == min(LAmeasacc0), 1);
    % if Run == 3 || Run == 23 || Run == 93 || Run == 128 
    %     endsteadystateindex = find(LAmeasacc1 == min(LAmeasacc1), 1);
    % end
    
    % endsteadystatetime = LAtime(endsteadystateindex + impactindex);
    endsteadystatetime = LAtime(endsteadystateindex);
    
    % Finding the time of max bouyancy
    [submergence_peak, submergence_peak_time_index] = min(LAexppos0);
    submergence_peak_time = LAtime(submergence_peak_time_index+3);
    
    % Calculating trim 
    distbetweenLA = 0.64; % m
    distfromhome0 = LA0startpos + LAmeaspos0; 
    distfromhome1 = LA1startpos + LAmeaspos1;
    trim = atand((distfromhome1 - distfromhome0) ./ distbetweenLA);
    averagetrim = mean(trim);
    
    % % Isolating the slam event window so that profiles have a consistent index
    % when averaging
    endisolationwindowindex = preimpactindex + 100; 
    LAtime_isolated = LAtime(preimpactindex:endisolationwindowindex);
    LAmeaspos0_isolated = LAmeaspos0(preimpactindex:endisolationwindowindex);
    
    
    func = zeros(length(LAtime), 18);
    func(1:length(LAtime), 1) = LAtime;
    func(1:length(LAexppos0), 2) = LAexppos0;
    func(1:length(LAexppos1), 3) = LAexppos1;
    func(1:length(LAmeaspos0), 4) = LAmeaspos0;
    func(1:length(LAmeaspos1), 5) = LAmeaspos1;
    func(1:length(LAexpvel0), 6) = LAexpvel0;
    func(1:length(LAexpvel1), 7) = LAexpvel1;
    func(1:length(LAmeasvel0), 8) = LAmeasvel0;
    func(1:length(LAmeasvel1), 9) = LAmeasvel1;
    func(1, 10) = impactcycletime;
    func(1, 11) = endsteadystatetime;
    func(1:length(trim), 12) = trim;
    func(1:length(LAmeasacc0), 13) = LAmeasacc0;
    func(1:length(LAmeasacc1), 14) = LAmeasacc1;
    func(1, 15) = averagetrim;
    func(1:length(LAtime_isolated), 16) = LAtime_isolated;
    func(1:length(LAtime_isolated), 17) = LAmeaspos0_isolated;
    func(1,18) = submergence_peak_time;
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

    % get zeros 
    offset = ForceZeros(Run);
    
    forcecycletimestart = force.XValuesGuaranteedValidOnlyForCDAQ1Mod2_ai0(11, :);
    forcetimesync = forcecycletimestart - impactcycletime;
    forcetime = force.Notes(11:end, :) + forcetimesync;
    % futek forces
        fxvoltage = force.Var3(11:end, :) - offset(1); 
        fyvoltage = force.Var4(11:end, :) - offset(2);
        fzvoltage = force.Var5(11:end, :) - offset(3);
    % omega force
        ozvoltage = force.Var6(11:end, :) - offset(4);
    
    % find preimpacttime index
    % preimpactindex =  find(forcetime > 0, 1, 'first');
    
    % Force time correction
    % expectedtimestep = 0.0005; 
    % actualtimestep = 0.000488;
    %forcetime = forcetime + ((expectedtimestep - actualtimestep) * preimpactindex);
    % forcetime = forcetime - 0.1;
    % 
    % Recalculate impact window index
    preimpactindex =  find(forcetime > 0, 1, 'first');
    endsteadystateindex = find(forcetime > endsteadystatetime, 1, 'first');
    
    % Calculating the calibration factor
    % signal conditoner gains
    fxgain = 1665;
    fygain = 1665;
    fzgain = 1665;
    
    % For Parametric Study
    if Run == 49 || Run == 81 || Run == 82
        fzgain = 509;
    end
    
    % For APSDFD
    % fzgain = 5000;
    % if Run > 64
    %     fzgain = 1665;
    % end
    
    %signal conditioner excitation
    fxexc = 10; %V
    fyexc = 10;
    fzexc = 10;
    
    % calibration sheet data from August 2023
    fxload = [0, 50, 100, 150, 200, 250]; %lbs
    fyload = [0, 50, 100, 150, 200, 250];
    fzload = [0, 100, 200, 300, 400, 500];
    
    fxoutput = [0, -0.3999, -0.7997, -1.1998, -1.6001, -2.006]; %mV/V
    fyoutput = [0, -0.3951, -0.7900, -1.1848, -1.5797, -1.9749];
    fzoutput = [0, -0.1993, -0.3987, -0.5981, -0.7977, -0.9972];
    
    fxcalfactor = calibrationfactor(fxgain, fxexc, fxload, fxoutput);
    fycalfactor = calibrationfactor(fygain, fyexc, fyload, fyoutput);
    fzcalfactor = calibrationfactor(fzgain, fzexc, fzload, fzoutput);
    ozcalfactor = 17.2939; % lbs/V, drawn from calibration data
    if Run == 49 || Run == 81 || Run == 82
        ozcalfactor = 17.2939 / 2;
    end
    
    fxforce = fxvoltage .* (fxcalfactor(1) * 4.44822); % lbs to N
    fyforce = fyvoltage .* (fycalfactor(1) * 4.44822);
    fzforce = fzvoltage .* (fzcalfactor(1) * 4.44822);
    ozforce = ozvoltage .* (ozcalfactor * 4.44822);
    totalzforce = fzforce + ozforce;
    
    if preimpactindex < endsteadystateindex
    [fxmax, maxfxindex] = max(fxforce(preimpactindex:endsteadystateindex));
    [fymax, maxfyindex] = max(fyforce(preimpactindex:endsteadystateindex));
    [fzmax, maxfzindex] = max(fzforce(preimpactindex:endsteadystateindex));
    [ozmax, maxozindex] = max(ozforce(preimpactindex:endsteadystateindex));
    [totalzmax, maxtzindex] = max(totalzforce(preimpactindex:endsteadystateindex));
    
    [fxmin, minfxindex] = min(fxforce(preimpactindex:endsteadystateindex));
    [fymin, minfyindex] = min(fyforce(preimpactindex:endsteadystateindex));
    [fzmin, minfzindex] = min(fzforce(preimpactindex:endsteadystateindex));
    [ozmin, minozindex] = min(ozforce(preimpactindex:endsteadystateindex));
    [totalzmin, mintzindex] = min(totalzforce(preimpactindex:endsteadystateindex));
    
        if abs(totalzmax) > abs(totalzmin)
            totalzpeak = totalzmax;
            impactduration = forcetime(preimpactindex:endsteadystateindex);
            totalzpeaktime = impactduration(maxtzindex);
        else 
            totalzpeak = totalzmin;
            impactduration = forcetime(preimpactindex:endsteadystateindex);
            totalzpeaktime = impactduration(mintzindex);
        end
    
    else 
        fxmax = 0; 
        fymax = 0;
        fzmax = 0; 
        ozmax = 0; 
        totalzpeak = 0;
        totalzpeaktime = 0;
        disp('This run does not have sufficient steady state window')
        disp(Run)
    end
    % OLD MET 1
    % Nondim method 1 - dynamic displacement
    % Get hydrostatic data to nondimensionalize
    Hydrostatics = GetHydroStaticForce(Run, keelentry, LA0startpos, LA1startpos, Heave);
    statictime = Hydrostatics(:,1);
    staticforce = Hydrostatics(:,2);
    ii = 1;
    representativeforcedata = zeros(length(staticforce), 1);
    for i = 1:length(forcetime)
        if forcetime(i) > statictime(ii)
            representativeforcedata(ii) = totalzforce(i);
            ii = ii + 1;
        end
        if ii > length(statictime)
            break
        end
    end 
    met1nondim = representativeforcedata ./ staticforce;
    
    % Nondim method 2 - Max length and beam
    % normal velocity
    R = [cosd(trim), -sind(trim); sind(trim), cosd(trim)]; %rotation matrix
    V = R * [-Heave; Surge];
    Vn = V(1); 
    V_T = sqrt(Heave^2 + Surge^2);
    
    % max length and beam of GPPH, m
    B = 0.36; 
    L = 1.22;
    
    % density of water
    rho = 1000; %kg / m3
    
    % nondimensionalization (is a ridiculously long word)
    % met2nondimscatter = representativeforcedata ./ (rho * Vn * L * B);
    met2nondimpeak = totalzpeak ./ (rho * Vn^2 * L * B);
    met2nondimplot = totalzforce ./ (rho * Vn^2 * L * B);
    
    % method 3 - neutral buoyancy
    displacement = 28.858 * 4.448;
    % met3nondimscatter = representativeforcedata ./ disp;
    met3nondimpeak = totalzpeak ./ displacement;
    met3nondimplot = totalzforce ./ displacement;
    
    % method 4 - dynamic pressure and wetted surface
    if trim == 0
        Sw = 0.419; % m2
    elseif trim == -5
        Sw = 0.323; % m2
    elseif trim == -10
        Sw = 0.196; % m2
    elseif trim == -15
        Sw = 0.067; % m2
    end
    
    met4nondimpeak = totalzpeak ./ (rho * Vn^2 * Sw);
    met4nondimplot = totalzforce ./ (rho * Vn^2 * Sw);
    
    % method 5 - specific case displacement 
    if trim == 0
        displacement_alt = 41.515 * 4.448; % N
    elseif trim == -5
        displacement_alt = 22.358 * 4.448; % N
    elseif trim == -10
        displacement_alt = 12.603 * 4.448; % N
    elseif trim == -15
        displacement_alt = 2.742 * 4.448; % N
    end
    
    met5nondimpeak = totalzpeak ./ displacement_alt;
    met5nondimplot = totalzforce ./ displacement_alt;
    
    % method 6 - test
    if trim == 0
        Sw = 0.419; % m2
    elseif trim == -5
        Sw = 0.323; % m2
    elseif trim == -10
        Sw = 0.196; % m2
    elseif trim == -15
        Sw = 0.067; % m2
    end
    
    nondimfactor_force = (rho * V_T^2 * Sw);
    met6nondimpeak = totalzpeak ./ (rho * V_T^2 * Sw);
    met6nondimplot = totalzforce ./ (rho * V_T^2 * Sw);
    
    % Isolating the slam event window so that forces have a consistent index
    % when averaging
    endisolationwindowindex = preimpactindex + 2000; 
    forcetime_isolated = forcetime(preimpactindex:endisolationwindowindex);
    fxforce_isolated = fxforce(preimpactindex:endisolationwindowindex);
    fyforce_isolated = fyforce(preimpactindex:endisolationwindowindex);
    fzforce_isolated = fzforce(preimpactindex:endisolationwindowindex);
    totalzforce_isolated = totalzforce(preimpactindex:endisolationwindowindex);
    
    func = zeros(length(forcetime), 31);
    func(1:length(forcetime), 1) = forcetime;
    func(1:length(forcetime), 2) = fxforce;
    func(1:length(forcetime), 3) = fyforce;
    func(1:length(forcetime), 4) = fzforce;
    func(1:length(forcetime), 5) = ozforce;
    func(1:length(forcetime), 6) = totalzforce;
    func(1, 7) = fxmax; 
    func(1, 8) = fymax; 
    func(1, 9) = fzmax; 
    func(1, 10) = ozmax;
    func(1, 11) = totalzpeak;
    func(1, 12) = totalzpeaktime;
    func(1:length(representativeforcedata), 13) = met1nondim;
    % func(1:length(representativeforcedata), 14) = met2nondimscatter;
    % func(1:length(representativeforcedata), 15) = met3nondimscatter;
    func(1:length(representativeforcedata), 16) = representativeforcedata;
    func(1:length(forcetime), 17) = met2nondimplot;
    func(1:length(forcetime), 18) = met3nondimplot;
    func(1:length(forcetime), 19) = met4nondimplot;
    func(1:length(forcetime), 20) = met5nondimplot;
    func(1:length(forcetime), 21) = met6nondimplot;
    func(1, 22) = met2nondimpeak;
    func(1, 23) = met3nondimpeak;
    func(1, 24) = met4nondimpeak;
    func(1, 25) = met5nondimpeak;
    func(1, 26) = met6nondimpeak;
    func(1:length(forcetime_isolated), 27) = forcetime_isolated;
    func(1:length(forcetime_isolated), 28) = fxforce_isolated;
    func(1:length(forcetime_isolated), 29) = fyforce_isolated;
    func(1:length(forcetime_isolated), 30) = fzforce_isolated;
    func(1:length(forcetime_isolated), 31) = totalzforce_isolated;
        else
            func = zeros(200000, 31);
            Vn = zeros(1, 1) ;
            nondimfactor_force = zeros(1, 1) ;
        end
    else
        func = zeros(200000, 31);    
    end
end

% Hydrostatic Force
function func = GetHydroStaticForce(Run, keelentry, LA0startpos, LA1startpos, Heave) 
    % Get actuator position and time data
    LAData = GetLAData(Run, keelentry, LA0startpos, LA1startpos, Heave); 
    LAtime = LAData(:, 1);
    LAmeaspos0 = LAData(:, 4);
    
    % hydrostatic forces taken from Rhino model
    staticpos = LAmeaspos0(103:121);
    % staticforce = [0.000, 0.212, 0.870, 1.985, 3.558, 5.587, 8.073, 11.009, 14.391, 18.212] * 4.448; %lbf to N;
    % staticforce = [0.000, 0.642, 2.629, 5.973, 10.666, 16.702, 23.576, 30.820, 38.526, 46.259] * 4.448; %lbf to N
    staticforce = [
    0.003
    0.544
    1.802
    2.927
    5.451
    8.275
    12.960
    17.496
    23.072
    28.024
    30.996
    36.891
    42.939
    45.547
    42.568
    51.982
    51.600
    51.408
    55.832] * 4.448; 
    
    statictime = zeros(1, length(staticforce));
    
    % find when hydrostatic forces occur in LAtime
    ii = 1; 
    for i = 1:length(LAmeaspos0)
        if LAmeaspos0(i) == staticpos(ii)
            statictime(ii) = LAtime(i);
            ii = ii + 1;
            
        end
        if ii > length(staticpos)
            break
        end
    end 
    
    % if trim == 0 
    %     staticpos = 1:1:4; 
    %     staticforce = 1:1:4;
    % elseif trim == 5
    %     staticpos = 1:1:4; 
    %     staticforce = 1:1:4;
    % end
    
    func = zeros(length(statictime), 2);
    func(1:length(statictime), 1) = statictime;
    func(1:length(statictime), 2) = staticforce;
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
    
    % get zeros
    offset = PressureZeros(Run);
    
    Gain = 10;
    
    p11voltage = ACandPressure.Var6(11:end, :) ./ Gain - offset(1);
    p12voltage = ACandPressure.Var7(11:end, :) ./ Gain - offset(2);
    p21voltage = ACandPressure.Var8(11:end, :) ./ Gain - offset(3);
    p22voltage = ACandPressure.Var9(11:end, :) ./ Gain - offset(4);
    
    
    %calfactors = [104.15, 100.48, 102.45, 104.60] / 1000; % V / psi [P11 P12 P21 P22]
    calfactors = [15.106, 15.272, 14.859, 15.171] / (1000); % V / kPa [P11 P12 P21 P22]
    
    p11pressure = p11voltage ./ calfactors(1); 
    p12pressure = p12voltage ./ calfactors(2); 
    p21pressure = p21voltage ./ calfactors(3); 
    p22pressure = p22voltage ./ calfactors(4);
    
    % Finding the time of the peaks
    preimpactindex = find(pressuretime>0.001, 1, 'first');
    endimpactindex = find(pressuretime>0.1, 1, 'first');
    
    % Second offset because of pressure sensor drift during run
    p11drift = p11pressure(preimpactindex);
    p12drift = p12pressure(preimpactindex);
    p21drift = p21pressure(preimpactindex);
    p22drift = p22pressure(preimpactindex);
    
    p11pressure = p11pressure - p11drift;
    p12pressure = p12pressure - p12drift;
    p21pressure = p21pressure - p21drift;
    p22pressure = p22pressure - p22drift;
    
    % Peak Pressure Magnitude
    [p11max_mag, p11max_index] = max(p11pressure(preimpactindex:endimpactindex));
    [p12max_mag, p12max_index] = max(p12pressure(preimpactindex:endimpactindex));
    [p21max_mag, p21max_index] = max(p21pressure(preimpactindex:endimpactindex));
    [p22max_mag, p22max_index] = max(p22pressure(preimpactindex:endimpactindex));
    
    pressuremax_mag = max([p11max_mag, p12max_mag, p21max_mag, p22max_mag]); 
    
    % Peak Pressure Time
    timewindow = pressuretime(preimpactindex:endimpactindex);
    p11max_time = timewindow(p11max_index);
    p12max_time = timewindow(p12max_index);
    p21max_time = timewindow(p21max_index);
    p22max_time = timewindow(p22max_index);
    
    % Isolating the slam event window so that pressures have a consistent index
    % when averaging
    endisolationwindowindex = preimpactindex + 20000 - 1000;
    pressuretime_isolated = pressuretime(preimpactindex-1000:endisolationwindowindex);
    p11pressure_isolated = p11pressure(preimpactindex-1000:endisolationwindowindex);
    p12pressure_isolated = p12pressure(preimpactindex-1000:endisolationwindowindex);
    p21pressure_isolated = p21pressure(preimpactindex-1000:endisolationwindowindex);
    p22pressure_isolated = p22pressure(preimpactindex-1000:endisolationwindowindex);
    
    pressuretolerance_p11 = 0.5;
    pressuretolerance_p12 = 0.2;
    risingslopeindex_p11 = find(p11pressure_isolated > pressuretolerance_p11, 1, 'first');
    risingslopeindex_p11_global = risingslopeindex_p11 + preimpactindex-1000;
    risingslopeindex_p12 = find(p12pressure_isolated > pressuretolerance_p12, 1, 'first');
    risingslopeindex_p12_global = risingslopeindex_p12 + preimpactindex-1000;
    % risingslopeindex_p21 = find(p21pressure > pressuretolerance, 1, 'first');
    % risingslopeindex_p22 = find(p22pressure > pressuretolerance, 1, 'first');
    endisolationwindowindex_p11 = risingslopeindex_p11 + 10000;
    endisolationwindowindex_p11_global = risingslopeindex_p11_global + 20000;
    endisolationwindowindex_p12 = risingslopeindex_p12 + 10000;
    endisolationwindowindex_p12_global = risingslopeindex_p12_global + 20000;
    shift = 500;
    shift_global = 5000;
    %P11
    % p11pressuretime_revised = pressuretime_isolated(risingslopeindex_p11-shift:endisolationwindowindex_p11-shift);
    p11pressuretime_revised = pressuretime(risingslopeindex_p11_global-shift_global : endisolationwindowindex_p11_global-shift_global);
    p11pressuretime_rise_revised =  pressuretime_isolated(risingslopeindex_p11);
    % p11pressure_revised = p11pressure_isolated(risingslopeindex_p11-shift:endisolationwindowindex_p11-shift);
    p11pressure_revised = p11pressure(risingslopeindex_p11_global - shift_global : endisolationwindowindex_p11_global-shift_global);
    
    %P12
    % p12pressuretime_revised = pressuretime_isolated(risingslopeindex_p12-shift:endisolationwindowindex_p12-shift);
    p12pressuretime_revised = pressuretime(risingslopeindex_p12_global-shift_global:endisolationwindowindex_p12_global-shift_global);
    % p12pressure_revised = p12pressure_isolated(risingslopeindex_p12-shift:endisolationwindowindex_p12-shift);
    p12pressure_revised = p12pressure(risingslopeindex_p12_global - shift_global : endisolationwindowindex_p12_global-shift_global);
    
    func = zeros(length(pressuretime), 15);
    func(1:length(pressuretime), 1) = pressuretime;
    func(1:length(pressuretime), 2) = p11pressure;
    func(1:length(pressuretime), 3) = p12pressure;
    func(1:length(pressuretime), 4) = p21pressure;
    func(1:length(pressuretime), 5) = p22pressure;
    func(1,6) = p11max_mag; 
    func(2,6) = p12max_mag;
    func(3,6) = p21max_mag;
    func(4,6) = p22max_mag;
    func(5,6) = pressuremax_mag;
    func(6,6) = p11max_time;
    func(7,6) = p12max_time;
    func(8,6) = p21max_time;
    func(9,6) = p22max_time;
    func(10,6) = preimpactindex;
    func(1:length(pressuretime_isolated), 7) = pressuretime_isolated;
    func(1:length(pressuretime_isolated), 8) = p11pressure_isolated;
    func(1:length(pressuretime_isolated), 9) = p12pressure_isolated;
    func(1:length(pressuretime_isolated), 10) = p21pressure_isolated;
    func(1:length(pressuretime_isolated), 11) = p22pressure_isolated;
    func(1:length(p11pressuretime_revised), 12) = p11pressuretime_revised;
    func(1:length(p11pressuretime_revised), 13) = p11pressure_revised;
    func(1:length(p11pressuretime_revised), 14) = p12pressuretime_revised;
    func(1:length(p11pressuretime_revised), 15) = p12pressure_revised;
        else
            func = zeros(200000, 15);
            p11pressuretime_rise_revised = zeros(1, 1) ;
        end
    else
        func = zeros(200000, 15);
    end
end

% Pressure Shift
function func = ShiftP21andP22(Run, keelentry, LA0startpos, LA1startpos, Heave, p11pressuretime_rise_avg_revised)
    [PressureData, p11pressuretime_rise_revised] = GetPressureData(Run, keelentry, LA0startpos, LA1startpos, Heave);
    % pressuretime_isolated = PressureData(1:20000, 7);
    % p21pressure_isolated = PressureData(1:20000, 10);
    % p22pressure_isolated = PressureData(1:20000, 11);
    % p11pressuretime_revised = PressureData(1:20000, 12);
    pressuretime = PressureData(1:100000, 1);
    p21pressure = PressureData(1:100000, 4);
    p22pressure = PressureData(1:100000, 5);
    
    deltat = p11pressuretime_rise_avg_revised - p11pressuretime_rise_revised;
    % pressuretime_isolated_aligned = pressuretime_isolated + deltat;
    % pressuretime_revised_aligned = p11pressuretime_revised + deltat;
    pressuretime_aligned = pressuretime + deltat;
    % shiftindex = find(pressuretime_isolated_aligned > 0.001, 1, 'first');
    shiftindex_global = find(pressuretime_aligned > 0.001, 1, 'first');
    % startwindow = shiftindex; 
    startwindow_global = shiftindex_global;
    % endwindow = startwindow + 10000;
    endwindow_global = startwindow_global + 25000;
    
    % pressuretime_isolated_shifted = pressuretime_isolated_aligned(startwindow:endwindow);
    % p21pressure_isolated_shifted = p21pressure_isolated(startwindow:endwindow);
    % p22pressure_isolated_shifted = p22pressure_isolated(startwindow:endwindow);
    
    pressuretime_shifted = pressuretime_aligned(startwindow_global:endwindow_global);
    p21pressure_shifted = p21pressure(startwindow_global:endwindow_global);
    p22pressure_shifted = p22pressure(startwindow_global:endwindow_global);
    
    % func = zeros(length(p21pressure_isolated_shifted), 3);
    % func(1:length(p21pressure_isolated_shifted), 1) = pressuretime_isolated_shifted(1:length(p21pressure_isolated_shifted));
    % func(1:length(p21pressure_isolated_shifted), 2) = p21pressure_isolated_shifted;
    % func(1:length(p21pressure_isolated_shifted), 3) = p22pressure_isolated_shifted;
    func = zeros(length(p21pressure_shifted), 3);
    func(1:length(p21pressure_shifted), 1) = pressuretime_shifted;
    func(1:length(p21pressure_shifted), 2) = p21pressure_shifted;
    func(1:length(p21pressure_shifted), 3) = p22pressure_shifted;
end 

% Accelerometer Data
function func = GetAccelerometerData(Run, keelentry, LA0startpos, LA1startpos, Heave)
    if exist(Run + "ACandPressure.csv", 'file') ~= 0
    [ACandPressure] = readtable(Run + "ACandPressure.csv");
        if size(ACandPressure, 1) ~= 0
    % find LA cycle time at impact event
    LAData = GetLAData(Run, keelentry, LA0startpos, LA1startpos, Heave); 
        if Run ~= 99
            impactcycletime = LAData(1, 10);
            endsteadystatetime = LAData(1,11);
        else
            impactcycletime = 1.0014e03 ;
            endsteadystatetime = 0.0832 ;
        end
    
    accelerometercycletimestart = ACandPressure.XValuesGuaranteedValidOnlyForCDAQ1Mod3_ai0(11, :);
    accelerometertimesync = accelerometercycletimestart - impactcycletime;
    accelerometertime = ACandPressure.Notes(11:end, :) + accelerometertimesync;
    
    % get zeros
    offset = AccelerometerZeros(Run); 
    
    a1voltage = ACandPressure.Var3(11:end, :) - offset(1);
    a2voltage = ACandPressure.Var4(11:end, :) - offset(2);
    a3voltage = ACandPressure.Var5(11:end, :) - offset(3);
    
    calfactors = [9.37, 9.82, 10.01] / 1000; % V / m/s2
    
    a1acceleration = a1voltage ./ calfactors(1); 
    a2acceleration = a2voltage ./ calfactors(2); 
    a3acceleration = a3voltage ./ calfactors(3); 
    
    func = zeros(length(accelerometertime), 4);
    func(1:length(accelerometertime), 1) = accelerometertime;
    func(1:length(accelerometertime), 2) = a1acceleration;
    func(1:length(accelerometertime), 3) = a2acceleration;
    func(1:length(accelerometertime), 4) = a3acceleration;
        else
            func = zeros(290, 4) ;
        end
    else
        func = zeros(290, 4) ;
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
    % get zeros
    offset = PotandIncZeros(Run); 
    
    % sync time
    PotandInccycletimestart = PotandInc.XValuesGuaranteedValidOnlyForCDAQ1Mod1_ai0(11, :);
    PotandInctimesync = PotandInccycletimestart - impactcycletime;
    PotandInctime = PotandInc.Notes(11:end, :) + PotandInctimesync;
    
    Potvoltage = PotandInc.Var3(11:end, :) - offset(1);
    Incvoltage = PotandInc.Var4(11:end, :) - offset(2);
    
    calfactors = [8.7199/39.37, 22.2341]; % m / V, deg / V
    
    Pot = keelentry - Potvoltage .* calfactors(1); 
    Inc = Incvoltage .* calfactors(2); 
    
    % Design the Butterworth filter 
    n_order = 4 ;              % Filter order // only ranges 1-10
    cutoff_freq = 100 ;        % Cutoff frequency in Hz // use fft to figure out

    sampling_rate = 1 / (PotandInctime(2) - PotandInctime(1)) ;  %% Not sure if this is correct // 
    % the sampling rate of the Butterworth filter is 1Hz at each the time T = 1/1000Hz 
    % (1000Hz = sampling rate of the input signal)

    norm_cutoff = cutoff_freq / (sampling_rate / 2); % Normalized cutoff
    [b, a] = butter(n_order, norm_cutoff, 'low');

    Pot = filtfilt(b, a, Pot) ;
    Inc = filtfilt(b, a, Inc) ;

    averagetrim = mean(Inc);
    
    func = zeros(length(PotandInctime), 4);
    func(1:length(PotandInctime), 1) = PotandInctime;
    func(1:length(Pot), 2) = Pot;
    func(1:length(Inc), 3) = Inc;
    func(1, 4) = averagetrim;
        else
            func = zeros(10000, 4) ;
        end
    else
        func = zeros(10000, 4) ;
    end
end

% Force Zeros
function func = ForceZeros(Run)
    [force0s] = readtable(Run + "Force_ZEROS.csv");
    fx0s = mean(force0s.Var2(11:end, :));
    fy0s = mean(force0s.Var3(11:end, :));
    fz0s = mean(force0s.Var4(11:end, :));
    oz0s = mean(force0s.Var5(11:end, :));
    
    func = zeros(4, 1);
    func(1) = fx0s; 
    func(2) = fy0s; 
    func(3) = fz0s; 
    func(4) = oz0s;
end

% Pressure Zeros
function func = PressureZeros(Run)
    [p0s] = readtable(Run + "ACandPressure_ZEROS.csv");
    p110s = mean(p0s.Var5(11:end, :));
    p120s = mean(p0s.Var6(11:end, :));
    p210s = mean(p0s.Var9(11:end, :));
    p220s = mean(p0s.Var8(11:end, :));
    
    func = zeros(4, 1);
    func(1) = p110s; 
    func(2) = p120s; 
    func(3) = p210s; 
    func(4) = p220s;
end

% Accelerometer Zeros
function func = AccelerometerZeros(Run)
    [ac0s] = readtable(Run + "ACandPressure_ZEROS.csv");
    ac10s = mean(ac0s.Var2(11:end, :));
    ac20s = mean(ac0s.Var3(11:end, :));
    ac30s = mean(ac0s.Var4(11:end, :));
    
    func = zeros(3, 1);
    func(1) = ac10s; 
    func(2) = ac20s; 
    func(3) = ac30s; 
end

% Potentiometer and Inclinometer Zeros
function func = PotandIncZeros(Run)
    [PI0s] = readtable(Run + "PotandInc_ZEROS.csv");
    P0s = mean(PI0s.Var2(11:end, :));
    I0s = 2.5107; % Taken from 0 deg case calibration of calibration data  %mean(PI0s.Var3(11:end, :)); 
    func = zeros(2, 1);
    func(1) = P0s; 
    func(2) = I0s; 
end

% Calibration Factors
function func = calibrationfactor(gain, exc, load, output)
    Vout = output .* exc .* gain ./ 1000; % V
    P = polyfit(Vout, load, 1) ;
    CalFactor = P(1) ; % Load Unit / V
    ZeroCheck = P(2) ;
    func(1) = CalFactor ; 
    func(2) = ZeroCheck ;
end
