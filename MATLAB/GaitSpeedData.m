clear ; clc ; close all

%% Gait Speed Data Slow Trial 

% Data Extraction
Slow1 = readmatrix("Slow1.csv") ; Slow2 = readmatrix("Slow2.xlsx") ; Slow3 = readmatrix("Slow3.xlsx") ;
Fast1 = readmatrix("Fast1.xlsx") ; Fast2 = readmatrix("Fast2.xlsx") ; Fast3 = readmatrix("Fast3.xlsx") ;
Slow1 = Slow1(7:end,:) ; Slow2 = Slow2(7:end, :) ; Slow3 = Slow3(7:end, :) ;
Fast1 = Fast1(7:end,:) ; Fast2 = Fast2(7:end, :) ; Fast3 = Fast3(7:end, :) ;
Slow1(:,2) = Slow1(:,2) - 3.74 ; Slow2(:,2) = Slow2(:,2) - 4.12 ; Slow3(:,2) = Slow3(:,2) - 4.92 ;
Fast1(:,2) = Fast1(:,2) - 3.67 ; Fast2(:,2) = Fast2(:,2) - 3.03 ; Fast3(:,2) = Fast3(:,2) - 3.01 ;
stzero1 = Slow1(:,2) ; stzero2 = Slow2(:,2) ; stzero3 = Slow3(:,2) ;           % Time where gait is stable is t = 0
ftzero1 = Fast1(:,2) ; ftzero2 = Fast2(:,2) ; ftzero3 = Fast3(:,2) ;
sSter1 = Slow1(:,6:8) ; sSter2 = Slow2(:,6:8) ; sSter3 = Slow3(:,6:8) ;        % Sternum 
fSter1 = Fast1(:,6:8) ; fSter2 = Fast2(:,6:8) ; fSter3 = Fast3(:,6:8) ;
sLGRT1= Slow1(:,21:23) ; sLGRT2 = Slow2(:,21:23) ; sLGRT3= Slow3(:,21:23) ;    % Left Greater Trochanter
fLGRT1= Fast1(:,21:23) ; fLGRT2 = Fast2(:,21:23) ; fLGRT3= Fast3(:,21:23) ;
sRGRT1 = Slow1(:,42:44) ; sRGRT2 = Slow2(:,42:44) ; sRGRT3 = Slow3(:,42:44) ;  % Right Greater Trochanter
fRGRT1 = Fast1(:,42:44) ; fRGRT2 = Fast2(:,42:44) ; fRGRT3 = Fast3(:,42:44) ;
sLLK1 = Slow1(:,24:26) ; sLLK2 = Slow2(:,24:26) ; sLLK3 = Slow3(:,24:26) ;     % Left Lateral Knee
fLLK1 = Fast1(:,24:26) ; fLLK2 = Fast2(:,24:26) ; fLLK3 = Fast3(:,24:26) ;
sRLK1 = Slow1(:,45:47) ; sRLK2 = Slow2(:,45:47) ; sRLK3 = Slow3(:,45:47) ;     % Right Lateral Knee
fRLK1 = Fast1(:,45:47) ; fRLK2 = Fast2(:,45:47) ; fRLK3 = Fast3(:,45:47) ;
sLLA1 = Slow1(:,30:32) ; sLLA2 = Slow2(:,30:32) ; sLLA3 = Slow3(:,30:32) ;     % Left Lateral Ankle
fLLA1 = Fast1(:,30:32) ; fLLA2 = Fast2(:,30:32) ; fLLA3 = Fast3(:,30:32) ;
sRLA1 = Slow1(:,51:53) ; sRLA2 = Slow2(:,51:53) ; sRLA3 = Slow3(:,51:53) ;     % Right Lateral Ankle
fRLA1 = Fast1(:,51:53) ; fRLA2 = Fast2(:,51:53) ; fRLA3 = Fast3(:,51:53) ;
sLH1 = Slow1(:,36:38) ; sLH2 = Slow2(:,36:38) ; sLH3 = Slow3(:,36:38) ;        % Left Heel
fLH1 = Fast1(:,36:38) ; fLH2 = Fast2(:,36:38) ; fLH3 = Fast3(:,36:38) ;
sRH1 = Slow1(:,57:59) ; sRH2 = Slow2(:,57:59) ; sRH3 = Slow3(:,57:59) ;        % Right Heel
fRH1 = Fast1(:,57:59) ; fRH2 = Fast2(:,57:59) ; fRH3 = Fast3(:,57:59) ;
 

%% i) Stride Length from Heel Markers
%{
% Left Side Slow Trial Stride Calculations
sLSL1 = stridelength(sLH1(stzero1 >= 0,:), stzero1(stzero1 >= 0)) ; 
sLSL2 = stridelength(sLH2(stzero2 >= 0,:), stzero2(stzero2 >= 0)) ; 
sLSL3 = stridelength(sLH3(stzero3 >= 0,:), stzero3(stzero3 >= 0)) ; 

% Left Side Fast Trial Stride Calculations
fLSL1 = stridelength(fLH1(ftzero1 >= 0,:), ftzero1(ftzero1 >= 0)) ; 
fLSL2 = stridelength(fLH2(ftzero2 >= 0,:), ftzero2(ftzero2 >= 0)) ; 
fLSL3 = stridelength(fLH3(ftzero3 >= 0,:), ftzero3(ftzero3 >= 0)) ; 

% Right Side Slow Trial Stride Calculations
sRSL1 = stridelength(sRH1(stzero1 >= 0,:), stzero1(stzero1 >= 0)) ; 
sRSL2 = stridelength(sRH2(stzero2 >= 0,:), stzero2(stzero2 >= 0)) ; 
sRSL3 = stridelength(sRH3(stzero3 >= 0,:), stzero3(stzero3 >= 0)) ; 

% Right Side Fast Trial Stride Calculations
fRSL1 = stridelength(fRH1(ftzero1 >= 0,:), ftzero1(ftzero1 >= 0)) ; 
fRSL2 = stridelength(fRH2(ftzero2 >= 0,:), ftzero2(ftzero2 >= 0)) ; 
fRSL3 = stridelength(fRH3(ftzero3 >= 0,:), ftzero3(ftzero3 >= 0)) ; 


% Left Side Stride
figure()
plot(stzero1(stzero1 >= 0), sLH1(stzero1 >= 0,1), ":.", 'linewidth', 1.25) ;
hold on
plot(stzero2(stzero2 >= 0), sLH2(stzero2 >= 0,1), ":.", 'linewidth', 1.25) ;
plot(stzero3(stzero3 >= 0), sLH3(stzero3 >= 0,1), ":.", 'linewidth', 1.25) ;
plot(ftzero1(ftzero1 >= 0), fLH1(ftzero1 >= 0,1), 'linewidth', 1.25) ;
plot(ftzero2(ftzero2 >= 0), fLH2(ftzero2 >= 0,1), 'linewidth', 1.25) ;
plot(ftzero3(ftzero3 >= 0), fLH3(ftzero3 >= 0,1), 'linewidth', 1.25) ;
grid on
xlabel("Time [s]") ;
ylabel('Heel Displacemet [mm]') ;
title('Left Heel Displacement over Time') ;
s1 = sprintf('Slow Trial 1 StrL Xavg = %5f ft', sLSL1) ; s2 = sprintf('Slow Trial 2 StrL Xavg = %5f ft', sLSL2) ;
s3 = sprintf('Slow Trial 3 StrL Xavg = %5f ft', sLSL3) ;
f1 = sprintf('Fast Trial 1 StrL Xavg = %5f ft', fLSL1) ; f2 = sprintf('Fast Trial 2 StrL Xavg = %5f ft', fLSL2) ;
f3 = sprintf('Fast Trial 3 StrL Xavg = %5f ft', fLSL3) ;
legend(s1, s2, s3, f1, f2, f3, 'location', 'best') ;

% Right Side Stride
figure()
plot(stzero1(stzero1 >= 0), sRH1(stzero1 >= 0,1), ":.", 'linewidth', 1.25) ;
hold on
plot(stzero2(stzero2 >= 0), sRH2(stzero2 >= 0,1), ":.", 'linewidth', 1.25) ;
plot(stzero3(stzero3 >= 0), sRH3(stzero3 >= 0,1), ":.", 'linewidth', 1.25) ;
plot(ftzero1(ftzero1 >= 0), fRH1(ftzero1 >= 0,1), 'linewidth', 1.25) ;
plot(ftzero2(ftzero2 >= 0), fRH2(ftzero2 >= 0,1), 'linewidth', 1.25) ;
plot(ftzero3(ftzero3 >= 0), fRH3(ftzero3 >= 0,1), 'linewidth', 1.25) ;
grid on
xlabel("Time [s]") ;
ylabel('Heel Displacemet [mm]') ;
title('Right Heel Displacement over Time') ;
s1 = sprintf('Slow Trial 1 StrL Xavg = %5f ft', sRSL1) ; s2 = sprintf('Slow Trial 2 StrL Xavg = %5f ft', sRSL2) ;
s3 = sprintf('Slow Trial 3 StrL Xavg = %5f ft', sRSL3) ;
f1 = sprintf('Fast Trial 1 StrL Xavg = %5f ft', fRSL1) ; f2 = sprintf('Fast Trial 2 StrL Xavg = %5f ft', fRSL2) ;
f3 = sprintf('Fast Trial 3 StrL Xavg = %5f ft', fRSL3) ;
legend(s1, s2, s3, f1, f2, f3, 'location', 'best') ;
hold off




% Comparing the Two Tests
% Left Side Comparison
sMLSL = mean([sLSL1 sLSL2 sLSL3]) ; sSDLSL = std([sLSL1 sLSL2 sLSL3]) ;
fMLSL = mean([fLSL1 fLSL2 fLSL3]) ; fSDLSL = std([fLSL1 fLSL2 fLSL3]) ;
errhigh = [sSDLSL, fSDLSL] ; errlow = [sSDLSL, fSDLSL] ;
Y = [sMLSL fMLSL]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Average Left Leg Stride Length') ;
ylabel('Stride Length [ft]') ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off

% Right Side Comparison
sMRSL = mean([sRSL1 sRSL2 sRSL3]) ; sSDRSL = std([sRSL1 sRSL2 sRSL3]) ;
fMRSL = mean([fRSL1 fRSL2 fRSL3]) ; fSDRSL = std([fRSL1 fRSL2 fRSL3]) ;
errhigh = [sSDRSL, fSDRSL] ; errlow = [sSDRSL, fSDRSL] ;
Y = [sMRSL fMRSL]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Average Right Leg Stride Length') ;
ylabel('Stride Length [ft]') ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off

%}

%% ii) Forward Walking Speed 
%{
% Slow Trials
svster1 = finitediffv(sSter1) ; svster1 = svster1/1000 ; svster1 = svster1(stzero1 >= 0) ;
svmax1 = findpeaks(svster1, "MinPeakProminence", 0.07) ; svmax1 = max(svmax1) ;
svavg1 = mean(svster1, 'omitnan') ;
svster2 = finitediffv(sSter2) ; svster2 = svster2/1000 ; svster2 = svster2(stzero2 >= 0) ;
svmax2 = findpeaks(svster2, "MinPeakProminence", 0.07) ; svmax2 = max(svmax2) ;
svavg2 = mean(svster2, 'omitnan') ;
svster3 = finitediffv(sSter3) ; svster3 = svster3/1000 ; svster3 = svster3(stzero3 >= 0) ;
svmax3 = findpeaks(svster3, "MinPeakProminence", 0.07) ; svmax3 = max(svmax3) ;
svavg3 = mean(svster3, 'omitnan') ;

% Fast Trials
fvster1 = finitediffv(fSter1) ; fvster1 = fvster1/1000 ; fvster1 = fvster1(ftzero1 >= 0) ;
fvmax1 = findpeaks(fvster1, "MinPeakProminence", 0.07) ; fvmax1 = max(fvmax1) ;
fvavg1 = mean(fvster1, 'omitnan') ;
fvster2 = finitediffv(fSter2) ; fvster2 = fvster2/1000 ; fvster2 = fvster2(ftzero2 >= 0) ;
fvmax2 = findpeaks(fvster2, "MinPeakProminence", 0.07) ; fvmax2 = max(fvmax2) ;
fvavg2 = mean(fvster2, 'omitnan') ;
fvster3 = finitediffv(fSter3) ; fvster3 = fvster3/1000 ; fvster3 = fvster3(ftzero3 >= 0) ;
fvmax3 = findpeaks(fvster3, "MinPeakProminence", 0.07) ; fvmax3 = max(fvmax3) ;
fvavg3 = mean(fvster3, 'omitnan') ;

% Walking Speed
figure()
plot(stzero1(stzero1 >= 0), svster1, ":.", 'linewidth', 1.5) ;
hold on 
plot(stzero2(stzero2 >= 0), svster2, ":.", "linewidth", 1.5) ;
plot(stzero3(stzero3 >= 0), svster3, ":.", "linewidth", 1.5) ;
plot(ftzero1(ftzero1 >= 0), fvster1, 'linewidth', 1.5) ;
plot(ftzero2(ftzero2 >= 0), fvster2, 'linewidth', 1.5) ;
plot(ftzero3(ftzero3 >= 0), fvster3, 'linewidth', 1.5) ;
title("Forward Walking Speed") ;
xlabel("Time [s]") ;
ylabel("Walking Speed [m/s]") ;
grid on
ylim([0 2]) ;
xlim([0 4]) ;
sv1 = sprintf("Slow Trial 1 Vavg = %5f m/s", svavg1) ; sv2 = sprintf("Slow Trial 2 Vavg = %5f m/s", svavg2) ; 
sv3 = sprintf("Slow Trial 3 Vavg = %5f m/s", svavg3) ;
fv1 = sprintf("Fast Trial 1 Vavg = %5f m/s", fvavg1) ; fv2 = sprintf("Fast Trial 2 Vavg = %5f m/s", fvavg2) ;
fv3 = sprintf("Fast Trial 3 Vavg = %5f m/s", fvavg3) ;
legend(sv1, sv2, sv3, fv1, fv2, fv3, "location", 'best') ; 


% Comparing the Two Tests
sMLSV = mean([svavg1 svavg2 svavg3]) ; sSDLSV = std([svavg1 svavg2 svavg3]) ;
fMLSV = mean([fvavg1 fvavg2 fvavg3]) ; fSDLSV = std([fvavg1 fvavg2 fvavg3]) ;
errhigh = [sSDLSV, fSDLSV] ; errlow = [sSDLSV, fSDLSV] ;
Y = [sMLSV fMLSV]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Average Walking Speed') ;
ylabel('Walking Speed [m/s]') ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off


%}

%% iii) Heel Forward Linear Velocity
%{
% Left Heel Slow Velocity Calculations
[svLH1, svL1] = finitediffv(sLH1) ; svLH1 = svLH1/1000 ; svLH1 = svLH1(stzero1 >= 0) ;         % Trial 1
svL1 = svL1/1000 ; svL1 = svL1(stzero1 >= 0,:) ; svxyzavg = mean(svL1(:,1:3),'omitnan') ; 
svxavg = svxyzavg(:,1) ; svyavg = svxyzavg(:,2) ; svzavg = svxyzavg(:,3) ;
svLHavg1 = mean(svLH1, 'omitnan') ; 
[svLH2, svL2] = finitediffv(sLH2) ; svLH2 = svLH2/1000 ; svLH2 = svLH2(stzero2 >= 0) ;             % Trial 2
svLHavg2 = mean(svLH2, 'omitnan') ; svL2 = svL2/1000 ; svL2 = svL2(stzero2 >= 0,:) ;
[svLH3, svL3] = finitediffv(sLH3) ; svLH3 = svLH3/1000 ; svLH3 = svLH3(stzero3 >= 0) ;             % Trial 3
svLHavg3 = mean(svLH3, 'omitnan') ; svL3 = svL3/1000 ; svL3 = svL3(stzero3 >= 0,:) ;

% Left Heel Fast Velocity Calculations
[fvLH1, fvL1] = finitediffv(fLH1) ; fvLH1 = fvLH1/1000 ; fvLH1 = fvLH1(ftzero1 >= 0) ;         % Trial 1
fvL1 = fvL1/1000 ; fvL1 = fvL1(ftzero1 >= 0,:) ; fvxyzavg = mean(fvL1(:,1:3),'omitnan') ; 
fvxavg = fvxyzavg(:,1) ; fvyavg = fvxyzavg(:,2) ; fvzavg = fvxyzavg(:,3) ;
fvLHavg1 = mean(fvLH1, 'omitnan') ; 
[fvLH2, fvL2] = finitediffv(fLH2) ; fvLH2 = fvLH2/1000 ; fvLH2 = fvLH2(ftzero2 >= 0) ;             % Trial 2
fvLHavg2 = mean(fvLH2, 'omitnan') ; fvL2 = fvL2/1000 ; fvL2 = fvL2(ftzero2 >= 0,:) ;
[fvLH3, fvL3] = finitediffv(fLH3) ; fvLH3 = fvLH3/1000 ; fvLH3 = fvLH3(ftzero3 >= 0) ;             % Trial 3
fvLHavg3 = mean(fvLH3, 'omitnan') ; fvL3 = fvL3/1000 ; fvL3 = fvL3(ftzero3 >= 0,:) ;

% Right Heel Slow Velocity Calculations
[svRH1, svR1] = finitediffv(sRH1) ; svRH1 = svRH1/1000 ; svRH1 = svRH1(stzero1 >= 0) ;             % Trial 1
svRHavg1 = mean(svRH1, 'omitnan') ; svR1 = svR1/1000 ; svR1 = svR1(stzero1 >= 0,:) ;
[svRH2, svR2] = finitediffv(sRH2) ; svRH2 = svRH2/1000 ; svRH2 = svRH2(stzero2 >= 0) ;             % Trial 2
svRHavg2 = mean(svRH2, 'omitnan') ; svR2 = svR2/1000 ; svR2 = svR2(stzero2 >= 0,:) ;
[svRH3, svR3] = finitediffv(sRH3) ; svRH3 = svRH3/1000 ; svRH3 = svRH3(stzero3 >= 0) ;             % Trial 3
svRHavg3 = mean(svRH3, 'omitnan') ; svR3 = svR3/1000 ; svR3 = svR3(stzero3 >= 0,:) ;

% Right Heel Fast Velocity Calculations
[fvRH1, fvR1] = finitediffv(fRH1) ; fvRH1 = fvRH1/1000 ; fvRH1 = fvRH1(ftzero1 >= 0) ;             % Trial 1
fvRHavg1 = mean(fvRH1, 'omitnan') ; fvR1 = fvR1/1000 ; fvR1 = fvR1(ftzero1 >= 0,:) ;
[fvRH2, fvR2] = finitediffv(fRH2) ; fvRH2 = fvRH2/1000 ; fvRH2 = fvRH2(ftzero2 >= 0) ;             % Trial 2
fvRHavg2 = mean(fvRH2, 'omitnan') ; fvR2 = fvR2/1000 ; fvR2 = fvR2(ftzero2 >= 0,:) ;
[fvRH3, fvR3] = finitediffv(fRH3) ; fvRH3 = fvRH3/1000 ; fvRH3 = fvRH3(ftzero3 >= 0) ;             % Trial 3
fvRHavg3 = mean(fvRH3, 'omitnan') ; fvR3 = fvR3/1000 ; fvR3 = fvR3(ftzero3 >= 0,:) ; 


% Slow Trial 1 Showing Resultant and Components
figure()
plot(stzero1(stzero1 >= 0), svLH1, 'linewidth', 1.25) ;
hold on
plot(stzero1(stzero1 >= 0), svL1, '--', 'linewidth', 1.25) ;
title('Left Heel Linear Velocity') ;
xlabel('Time [s]') ;
ylabel('Linear Velocity [m/s]') ;
sv1 = sprintf('Resultant Velocity Vavg = %5f m/s', svLHavg1) ; sv2 = sprintf('Heel x Velocity Vxavg = %5f m/s', svxavg) ;
sv3 = sprintf('Heel y Velocity Vyavg = %5f m/s', svyavg) ; sv4 = sprintf('Heel z Velocity Vzavg = %5f m/s', svzavg) ;
legend(sv1, sv2, sv3, sv4, 'location', 'best') ;
hold off

% Fast Trial 1 Showing Resultant and Components
figure()
plot(ftzero1(ftzero1 >= 0), fvLH1, 'linewidth', 1.25) ;
hold on
plot(ftzero1(ftzero1 >= 0), fvL1, '--', 'linewidth', 1.25) ;
title('Left Heel Linear Velocity') ;
xlabel('Time [s]') ;
ylabel('Linear Velocity [m/s]') ;
fv1 = sprintf('Resultant Velocity Vavg = %5f m/s', fvLHavg1) ; fv2 = sprintf('Heel x Velocity Vxavg = %5f m/s', fvxavg) ;
fv3 = sprintf('Heel y Velocity Vyavg = %5f m/s', fvyavg) ; fv4 = sprintf('Heel z Velocity Vzavg = %5f m/s', fvzavg) ;
legend(fv1, fv2, fv3, fv4, 'location', 'best') ;
hold off

% Left Heel Resultant Velocity
figure()
plot(stzero1(stzero1 >= 0), svLH1, ":.", 'linewidth', 1.25) ;
hold on
plot(stzero2(stzero2 >= 0), svLH2, ":.", 'linewidth', 1.25) ;
plot(stzero3(stzero3 >= 0), svLH3, ":.", 'linewidth', 1.25) ;
plot(ftzero1(ftzero1 >= 0), fvLH1, 'linewidth', 1.25) ;
plot(ftzero2(ftzero2 >= 0), fvLH2, 'linewidth', 1.25) ;
plot(ftzero3(ftzero3 >= 0), fvLH3, 'linewidth', 1.25) ;
title("Left Heel Resultant Velocity") ;
xlabel('Time [s]') ;
ylabel('Heel Velocity [m/s]') ;
sv1 = sprintf("Slow Trial 1 Vavg = %5f m/s", svLHavg1) ; sv2 = sprintf("Slow Trial 2 Vavg = %5f m/s", svLHavg2) ;
sv3 = sprintf("Slow Trial 3 Vavg = %5f m/s", svLHavg3) ;
fv1 = sprintf("Fast Trial 1 Vavg = %5f m/s", fvLHavg1) ; fv2 = sprintf("Fast Trial 2 Vavg = %5f m/s", fvLHavg2) ;
fv3 = sprintf("Fast Trial 3 Vavg = %5f m/s", fvLHavg3) ;
legend(sv1, sv2, sv3, fv1, fv2, fv3, 'location', 'best') ;
hold off

% Right Heel Resultant Velocity
figure()
plot(stzero1(stzero1 >= 0), svRH1, ":.", 'linewidth', 1.25) ;
hold on
plot(stzero2(stzero2 >= 0), svRH2, ":.", 'linewidth', 1.25) ;
plot(stzero3(stzero3 >= 0), svRH3, ":.", 'linewidth', 1.25) ;
plot(ftzero1(ftzero1 >= 0), fvRH1, 'linewidth', 1.25) ;
plot(ftzero2(ftzero2 >= 0), fvRH2, 'linewidth', 1.25) ;
plot(ftzero3(ftzero3 >= 0), fvRH3, 'linewidth', 1.25) ;
title("Right Heel Resultant Velocity") ;
xlabel('Time [s]') ;
ylabel('Heel Velocity [m/s]') ;
sv1 = sprintf("Slow Trial 1 Vavg = %5f m/s", svRHavg1) ; sv2 = sprintf("Slow Trial 2 Vavg = %5f m/s", svRHavg2) ;
sv3 = sprintf("Slow Trial 3 Vavg = %5f m/s", svRHavg3) ;
fv1 = sprintf("Fast Trial 1 Vavg = %5f m/s", fvRHavg1) ; fv2 = sprintf("Fast Trial 2 Vavg = %5f m/s", fvRHavg2) ;
fv3 = sprintf("Fast Trial 3 Vavg = %5f mls", fvRHavg3) ;
legend(sv1, sv2, sv3, fv1, fv2, fv3, 'location', 'best') ;
hold off

% Comparing the Two Tests
% Left Side Comparison
sMvLH = mean([svLHavg1 svLHavg2 svLHavg3]) ; sSDvLH = std([svLHavg1 svLHavg2 svLHavg3]) ;
fMvLH = mean([fvLHavg1 fvLHavg2 fvLHavg3]) ; fSDvLH = std([fvLHavg1 fvLHavg2 fvLHavg3]) ;
errhigh = [sSDvLH, fSDvLH] ; errlow = [sSDvLH, fSDvLH] ;
Y = [sMvLH fMvLH]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Average Left Heel Resultant Velocity') ;
ylabel('Resultant Velocity [m/s]') ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off

% Right Side Comparison
sMvRH = mean([svRHavg1 svRHavg2 svRHavg3]) ; sSDvRH = std([svRHavg1 svRHavg2 svRHavg3]) ;
fMvRH = mean([fvRHavg1 fvRHavg2 fvRHavg3]) ; fSDvRH = std([fvRHavg1 fvRHavg2 fvRHavg3]) ;
errhigh = [sSDvRH, fSDvRH] ; errlow = [sSDvRH, fSDvRH] ;
Y = [sMvRH fMvRH]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Average Right Heel Resultant Velocity') ;
ylabel('Resultant Velocity [m/s]') ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off

%}

%% iv) Heel Forward Linear Acceleration
%{
% Left Heel Slow Acceleration Calculations
[saLH1, saxyz] = finitediffa(svL1) ;                                                            % Trial 1
safLH1 = filtera(saLH1) ; safx = filtera(saxyz(:,1)) ; safy = filtera(saxyz(:,2)) ;
safz = filtera(saxyz(:,3)) ; safxyz = [safx, safy, safz] ;
saxavg = mean(safx,'omitnan') ; sayavg = mean(safy,'omitnan') ; sazavg = mean(safz,'omitnan') ;
saLHavg1 = mean(safLH1(safLH1 > 0)) ; 
[saLH2, ~] = finitediffa(svL2) ;                                                                % Trial 2
safLH2 = filtera(saLH2) ;
saLHavg2 = mean(safLH2(safLH2 > 0)) ; 
[saLH3, ~] = finitediffa(svL3) ;                                                                % Trial 3
safLH3 = filtera(saLH3) ;
saLHavg3 = mean(safLH3(safLH3 > 0)) ; 

% Left Heel Fast Acceleration Calculations
[faLH1, faxyz] = finitediffa(fvL1) ;                                                            % Trial 1 
fafLH1 = filtera(faLH1) ; fafx = filtera(faxyz(:,1)) ; fafy = filtera(faxyz(:,2)) ;
fafz = filtera(faxyz(:,3)) ; fafxyz = [fafx, fafy, fafz] ;
faxavg = mean(fafx,'omitnan') ; fayavg = mean(fafy,'omitnan') ; fazavg = mean(fafz,'omitnan') ;
faLHavg1 = mean(fafLH1(fafLH1 > 0)) ; 
[faLH2, ~] = finitediffa(fvL2) ;                                                                % Trial 2
fafLH2 = filtera(faLH2) ;
faLHavg2 = mean(faLH2(faLH2 > 0)) ; 
[faLH3, ~] = finitediffa(fvL3) ;                                                                % Trial 3
fafLH3 = filtera(faLH3) ;
faLHavg3 = mean(fafLH3(fafLH3 > 0)) ;

% Right Heel Slow Acceleration Calculations
[saRH1, ~] = finitediffa(svR1) ;                                                                % Trial 1
safRH1 = filtera(saRH1) ;
saRHavg1 = mean(safRH1(safRH1 > 0)) ; 
[saRH2, ~] = finitediffa(svR2) ;                                                                % Trial 2
safRH2 = filtera(saRH2) ;
saRHavg2 = mean(safRH2(safRH2 > 0)) ; 
[saRH3, ~] = finitediffa(svR3) ;                                                                % Trial 3
safRH3 = filtera(saRH3) ;
saRHavg3 = mean(safRH3(safRH3 > 0)) ; 

% Right Heel Fast Acceleration Calculations
[faRH1, ~] = finitediffa(fvR1) ;                                                                % Trial 1
fafRH1 = filtera(faRH1) ;
faRHavg1 = mean(fafRH1(fafRH1 > 0)) ; 
[faRH2, ~] = finitediffa(fvR2) ;                                                                % Trial 2
fafRH2 = filtera(faRH2) ;
faRHavg2 = mean(fafRH2(fafRH2 > 0)) ; 
[faRH3, ~] = finitediffa(fvR3) ;                                                                % Trial 3
fafRH3 = filtera(faRH3) ;
faRHavg3 = mean(fafRH3(fafRH3 > 0)) ; 

% Slow Trial 1 Showing Resultant and Components
figure()
plot(stzero1(stzero1 >= 0), safLH1, 'linewidth', 1.25) ;
hold on
plot(stzero1(stzero1 >= 0), safxyz, '--', 'linewidth', 1.25) ;
title('Left Heel Linear Acceleration') ;
xlabel('Time [s]') ;
ylabel('Linear Acceleration [m/s^2]') ;
svL1 = sprintf('Resultant Acceleration Aavg = %5f m/s^2', saLHavg1) ; svL2 = sprintf('Heel x Acceleration Axavg = %5f m/s^2', saxavg) ;
svL3 = sprintf('Heel y Acceleration Ayavg = %5f m/s^2', sayavg) ; sv4 = sprintf('Heel z Acceleration Azavg = %5f m/s^2', sazavg) ;
legend(svL1, svL2, svL3, sv4, 'location', 'best') ;
hold off

% Fast Trial 1 Showing Resultant and Components
figure()
plot(ftzero1(ftzero1 >= 0), fafLH1, 'linewidth', 1.25) ;
hold on
plot(ftzero1(ftzero1 >= 0), fafxyz, '--', 'linewidth', 1.25) ;
title('Left Heel Linear Acceleration') ;
xlabel('Time [s]') ;
ylabel('Linear Acceleration [m/s^2]') ;
svL1 = sprintf('Resultant Acceleration Aavg = %5f m/s^2', faLHavg1) ; svL2 = sprintf('Heel x Acceleration Axavg = %5f m/s^2', faxavg) ;
svL3 = sprintf('Heel y Acceleration Ayavg = %5f m/s^2', fayavg) ; sv4 = sprintf('Heel z Acceleration Azavg = %5f m/s^2', fazavg) ;
legend(svL1, svL2, svL3, sv4, 'location', 'best') ;
hold off

% Left Heel Resultant Acceleration
 figure()
 plot(stzero1(stzero1 >= 0), safLH1, ":.", "linewidth", 1.25) ;
 hold on
 plot(stzero2(stzero2 >= 0), safLH2, ":.", 'linewidth', 1.25) ;
 plot(stzero3(stzero3 >= 0), safLH3, ":.", 'linewidth', 1.25) ;
 plot(ftzero1(ftzero1 >= 0), fafLH1, 'linewidth', 1.25) ;
 plot(ftzero2(ftzero2 >= 0), fafLH2, 'linewidth', 1.25) ;
 plot(ftzero3(ftzero3 >= 0), fafLH3, 'linewidth', 1.25) ;
 title("Filtered Left Heel Resultant Acceleration") ;
 xlabel("Time [s]") ;
 ylabel("Linear Acceleration [m/s^2]") ;
 sa1 = sprintf("Slow Trial 1 Aavg = %5f m/s^2", saLHavg1) ; sa2 = sprintf("Slow Trial 2 Aavg = %5f m/s^2", saLHavg2) ; 
 sa3 = sprintf("Slow Trial 3 Aavg = %5f m/s^2", saLHavg3) ;
 fa1 = sprintf("Fast Trial 1 Aavg = %5f m/s^2", faLHavg1) ; fa2 = sprintf("Fast Trial 2 Aavg = %5f m/s^2", faLHavg2) ;
 fa3 = sprintf("Fast Trial 3 Aavg = %5f m/s^2", faLHavg3) ;
 legend(sa1, sa2, sa3, fa1, fa2, fa3, 'location', 'best') ;
 

% Right Heel Resultant Acceleration
 figure()
 plot(stzero1(stzero1 >= 0), safRH1, ":.", "linewidth", 1.25) ;
 hold on
 plot(stzero2(stzero2 >= 0), safRH2, ":.", 'linewidth', 1.25) ;
 plot(stzero3(stzero3 >= 0), safRH3, ":.", 'linewidth', 1.25) ;
 plot(ftzero1(ftzero1 >= 0), fafRH1, 'linewidth', 1.25) ;
 plot(ftzero2(ftzero2 >= 0), fafRH2, 'linewidth', 1.25) ;
 plot(ftzero3(ftzero3 >= 0), fafRH3, 'linewidth', 1.25) ;
 title("Filtered Right Heel Resultant Acceleration") ;
 xlabel("Time [s]") ;
 ylabel("Linear Acceleration [m/s^2]") ;
 sa1 = sprintf("Slow Trial 1 Aavg = %5f m/s^2", saRHavg1) ; sa2 = sprintf("Slow Trial 2 Aavg = %5f m/s^2", saRHavg2) ; 
 sa3 = sprintf("Slow Trial 3 Aavg = %5f m/s^2", saRHavg3) ;
 fa1 = sprintf("Fast Trial 1 Aavg = %5f m/s^2", faRHavg1) ; fa2 = sprintf("Fast Trial 2 Aavg = %5f m/s^2", faRHavg2) ;
 fa3 = sprintf("Fast Trial 3 Aavg = %5f m/s^2", faRHavg3) ;
 legend(sa1, sa2, sa3, fa1, fa2, fa3, 'location', 'best') ;


% Comparing the Two Tests 
% Left Side Comparison
sMaLH = mean([saLHavg1 saLHavg2 saLHavg3]) ; sSDaLH = std([saLHavg1 saLHavg2 saLHavg3]) ;
fMaLH = mean([faLHavg1 faLHavg2 faLHavg3]) ; fSDaLH = std([faLHavg1 faLHavg2 faLHavg3]) ;
errhigh = [sSDaLH, fSDaLH] ; errlow = [sSDaLH, fSDaLH] ;
Y = [sMaLH fMaLH]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Average Left Heel Resultant Acceleration') ;
ylabel('Resultant Acceleration [m/s^2]') ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off

% Right Side Comparison
sMaRH = mean([saRHavg1 saRHavg2 saRHavg3]) ; sSDaRH = std([saRHavg1 saRHavg2 saRHavg3]) ;
fMaRH = mean([faRHavg1 faRHavg2 faRHavg3]) ; fSDaRH = std([faRHavg1 faRHavg2 faRHavg3]) ;
errhigh = [sSDaRH, fSDaRH] ; errlow = [sSDaRH, fSDaRH] ;
Y = [sMaRH fMaRH]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Average Right Heel Resultant Acceleration') ;
ylabel('Resultant Acceleration [m/s^2]') ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off

%}

%% v) Knee Relative Joint Angle
%{
% Slow Left Side Knee Angle Calculations
[sltheta1,~,~] = angles(sLGRT1(stzero1 >= 0, :), sLLK1(stzero1 >= 0, :), sLLA1(stzero1 >= 0, :)) ;
sl1max = max(sltheta1) ;
[sltheta2,~,~] = angles(sLGRT2(stzero2 >= 0, :), sLLK2(stzero2 >= 0, :), sLLA2(stzero2 >= 0, :)) ;
sl2max = max(sltheta2) ;
[sltheta3,~,~] = angles(sLGRT3(stzero3 >= 0, :), sLLK3(stzero3 >= 0, :), sLLA3(stzero3 >= 0, :)) ;
sl3max = max(sltheta3) ;

% Fast Left Side Knee Angle Calculations
[fltheta1,~,~] = angles(fLGRT1(ftzero1 >= 0, :), fLLK1(ftzero1 >= 0, :), fLLA1(ftzero1 >= 0, :)) ;
fl1max = max(fltheta1) ;
[fltheta2,~,~] = angles(fLGRT2(ftzero2 >= 0, :), fLLK2(ftzero2 >= 0, :), fLLA2(ftzero2 >= 0, :)) ;
fl2max = max(fltheta2) ;
[fltheta3,~,~] = angles(fLGRT3(ftzero3 >= 0, :), fLLK3(ftzero3 >= 0, :), fLLA3(ftzero3 >= 0, :)) ;
fl3max = max(fltheta3) ;

% Slow Right Side Knee Angle Calculations
[srtheta1,~,~] = angles(sRGRT1(stzero1 >= 0, :), sRLK1(stzero1 >= 0, :), sRLA1(stzero1 >= 0, :)) ;
sr1max = max(srtheta1) ;
[srtheta2,~,~] = angles(sRGRT2(stzero2 >= 0, :), sRLK2(stzero2 >= 0, :), sRLA2(stzero2 >= 0, :)) ;
sr2max = max(srtheta2) ;
[srtheta3,~,~] = angles(sRGRT3(stzero3 >= 0, :), sRLK3(stzero3 >= 0, :), sRLA3(stzero3 >= 0, :)) ;
sr3max = max(srtheta3) ;

% Fast Right Side Knee Angle Calculations
[frtheta1,~,~] = angles(fRGRT1(ftzero1 >= 0, :), fRLK1(ftzero1 >= 0, :), fRLA1(ftzero1 >= 0, :)) ;
fr1max = max(frtheta1) ;
[frtheta2,~,~] = angles(fRGRT2(ftzero2 >= 0, :), fRLK2(ftzero2 >= 0, :), fRLA2(ftzero2 >= 0, :)) ;
fr2max = max(frtheta2) ;
[frtheta3,~,~] = angles(fRGRT3(ftzero3 >= 0, :), fRLK3(ftzero3 >= 0, :), fRLA3(ftzero3 >= 0, :)) ;
fr3max = max(frtheta3) ;

% Left Side Knee Angle
figure()
plot(stzero1(stzero1 >= 0), sltheta1, ":.", "linewidth", 1.25) ;
hold on
plot(stzero2(stzero2 >= 0), sltheta2, ":.", "linewidth", 1.25) ;
plot(stzero3(stzero3 >= 0), sltheta3, ":.", "linewidth", 1.25) ;
plot(ftzero1(ftzero1 >= 0), fltheta1, 'linewidth', 1.25) ;
plot(ftzero2(ftzero2 >= 0), fltheta2, 'linewidth', 1.25) ;
plot(ftzero3(ftzero3 >= 0), fltheta3, 'linewidth', 1.25) ;
title("Left Knee Relative Joint Angle") ;
xlabel("Time [s]") ;
ylabel("Joint Angle [deg]") ;
ylim([100 200]) ;
st1 = sprintf("Slow Trial 1 theta_{max} = %5f", sl1max) ; st2 = sprintf("Slow Trial 2 theta_{max} = %5f", sl2max) ;
st3 = sprintf("Slow Trial 3 theta_{max} = %5f", sl3max) ;
ft1 = sprintf("Fast Trial 1 theta_{max} = %5f", fl1max) ; ft2 = sprintf("Fast Trial 2 theta_{max} = %5f", fl3max) ;
ft3 = sprintf("Fast Trial 3 theta_{max} = %5f", fl3max) ;
legend(st1, st2, st3, ft1, ft2, ft3, "location", "best") ;
grid on


% Right Side Knee Angle
figure()
plot(stzero1(stzero1 >= 0), srtheta1, ":.", "linewidth", 1.25) ;
hold on
plot(stzero2(stzero2 >= 0), srtheta2, ":.", "linewidth", 1.25) ;
plot(stzero3(stzero3 >= 0), srtheta3, ":.", "linewidth", 1.25) ;
plot(ftzero1(ftzero1 >= 0), frtheta1, 'linewidth', 1.25) ;
plot(ftzero2(ftzero2 >= 0), frtheta2, 'linewidth', 1.25) ;
plot(ftzero3(ftzero3 >= 0), frtheta3, 'linewidth', 1.25) ;
title("Right Knee Relative Joint Angle") ;
xlabel("Time [s]") ;
ylabel("Joint Angle [deg]") ;
ylim([100 200]) ;
st1 = sprintf("Slow Trial 1 theta_{max} = %5f", sr1max) ; st2 = sprintf("Slow Trial 2 theta_{max} = %5f", sr2max) ;
st3 = sprintf("Slow Trial 3 theta_{max} = %5f", sr3max) ;
ft1 = sprintf("Fast Trial 1 theta_{max} = %5f", fr1max) ; ft2 = sprintf("Fast Trial 2 theta_{max} = %5f", fr3max) ;
ft3 = sprintf("Fast Trial 3 theta_{max} = %5f", fr3max) ;
legend(st1, st2, st3, ft1, ft2, ft3, "location", "best") ;
grid on

% Comparing the Two Tests 
% Left Side Comparison
sMla = mean([sl1max sl2max sl3max]) ; sSDla = std([sl1max sl2max sl3max]) ;
fMla = mean([fl1max fl2max fl3max]) ; fSDla = std([fl1max fl2max fl3max]) ;
errhigh = [sSDla, fSDla] ; errlow = [sSDla, fSDla] ;
Y = [sMla fMla]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Max Left Relative Knee Angle') ;
ylabel('Relative Knee Angle [deg]') ;
ylim([180 190]) ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off

% Right Side Comparison
sMra = mean([sr1max sr2max sr3max]) ; sSDra = std([sr1max sr2max sr3max]) ;
fMra = mean([fr1max fr2max fr3max]) ; fSDra = std([fr1max fr2max fr3max]) ;
errhigh = [sSDra, fSDra] ; errlow = [sSDra, fSDra] ;
Y = [sMra fMra]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Max Right Relative Knee Angle') ;
ylabel('Relative Knee Angle [deg]') ;
ylim([180 190]) ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off

%}

%% vi) Shank Angular Velocity
%{
% Left Slow Shank Angular Velocity Calculations
[~,~,slthetas1] = angles(sLGRT1(:,:), sLLK1(:,:), sLLA1(:,:)) ;
slav1 = finitediffav(slthetas1) ; slav1 = slav1(stzero1 >= 0) ; 
[~,~,slthetas2] = angles(sLGRT2(:,:), sLLK2(:,:), sLLA2(:,:)) ;
slav2 = finitediffav(slthetas2) ; slav2 = slav2(stzero2 >= 0) ; 
[~,~,slthetas3] = angles(sLGRT3(:,:), sLLK3(:,:), sLLA3(:,:)) ;
slav3 = finitediffav(slthetas3) ; slav3 = slav3(stzero3 >= 0) ;

% Left Fast Shank Angular Velocity Calculations
[~,~,flthetas1] = angles(fLGRT1(:,:), fLLK1(:,:), fLLA1(:,:)) ;
flav1 = finitediffav(flthetas1) ; flav1 = flav1(ftzero1 >= 0) ; 
[~,~,flthetas2] = angles(fLGRT2(:,:), fLLK2(:,:), fLLA2(:,:)) ;
flav2 = finitediffav(flthetas2) ; flav2 = flav2(ftzero2 >= 0) ; 
[~,~,flthetas3] = angles(fLGRT3(:,:), fLLK3(:,:), fLLA3(:,:)) ;
flav3 = finitediffav(flthetas3) ; flav3 = flav3(ftzero3 >= 0) ;

% Right Slow Shank Angular Velocity Calculations
[~,~,srthetas1] = angles(sRGRT1(:,:), sRLK1(:,:), sRLA1(:,:)) ;
srav1 = finitediffav(srthetas1) ; srav1 = srav1(stzero1 >= 0) ; 
[~,~,srthetas2] = angles(sRGRT2(:,:), sRLK2(:,:), sRLA2(:,:)) ;
srav2 = finitediffav(srthetas2) ; srav2 = srav2(stzero2 >= 0) ; 
[~,~,srthetas3] = angles(sRGRT3(:,:), sRLK3(:,:), sRLA3(:,:)) ;
srav3 = finitediffav(srthetas3) ; srav3 = srav3(stzero3 >= 0) ;

% Right Fast Shank Angular Velocity Calculations
[~,~,frthetas1] = angles(fRGRT1(:,:), fRLK1(:,:), fRLA1(:,:)) ;
frav1 = finitediffav(frthetas1) ; frav1 = frav1(ftzero1 >= 0) ; 
[~,~,frthetas2] = angles(fRGRT2(:,:), fRLK2(:,:), fRLA2(:,:)) ;
frav2 = finitediffav(frthetas2) ; frav2 = frav2(ftzero2 >= 0) ; 
[~,~,frthetas3] = angles(fRGRT3(:,:), fRLK3(:,:), fRLA3(:,:)) ;
frav3 = finitediffav(frthetas3) ; frav3 = frav3(ftzero3 >= 0) ;


% Left Side Shank Angular Velocity
figure()
plot(stzero1(stzero1 >= 0), slav1*pi/180, ":.", 'linewidth', 1.25) ;
hold on
plot(stzero2(stzero2 >= 0), slav2*pi/180, ":.", 'linewidth', 1.25) ;
plot(stzero3(stzero3 >= 0), slav3*pi/180, ":.", "linewidth", 1.25) ;
plot(ftzero1(ftzero1 >= 0), flav1*pi/180, 'linewidth', 1.25) ;
plot(ftzero2(ftzero2 >= 0), flav2*pi/180, 'linewidth', 1.25) ;
plot(ftzero3(ftzero3 >= 0), flav3*pi/180, 'linewidth', 1.25) ;
title("Left Shank Angular Velocity")
xlabel("Time [s]") ;
ylabel("Angular Velocity [rad/s]") ;
sav1 = sprintf("Slow Trial 1 omega_{max} = %5f rad/s", max(slav1*pi/180)) ; sav2 = sprintf("Slow Trial 2 omega_{max} = %5f rad/s", max(slav2*pi/180)) ;
sav3 = sprintf("Slow Trial 3 omega_{max} = %5f rad/s", max(slav3*pi/180)) ;
fav1 = sprintf("Fast Trial 1 omega_{max} = %5f rad/s", max(flav1*pi/180)) ; fav2 = sprintf("Fast Trial 2 omega_{max} = %5f rad/s", max(flav2*pi/180)) ; 
fav3 = sprintf("Fast Trial 3 omega_{max} = %5f rad/s", max(flav3*pi/180)) ; 
legend(sav1, sav2, sav3, fav1, fav2, fav3, 'location', 'best') ;
grid on
hold off

% Right Side Shank Angular Velocity
figure()
plot(stzero1(stzero1 >= 0), srav1*pi/180, ":.", 'linewidth', 1.25) ;
hold on
plot(stzero2(stzero2 >= 0), srav2*pi/180, ":.", 'linewidth', 1.25) ;
plot(stzero3(stzero3 >= 0), srav3*pi/180, ":.", "linewidth", 1.25) ;
plot(ftzero1(ftzero1 >= 0), frav1*pi/180, 'linewidth', 1.25) ;
plot(ftzero2(ftzero2 >= 0), frav2*pi/180, 'linewidth', 1.25) ;
plot(ftzero3(ftzero3 >= 0), frav3*pi/180, 'linewidth', 1.25) ;
title("Right Shank Angular Velocity")
xlabel("Time [s]") ;
ylabel("Angular Velocity [rad/s]") ;
sav1 = sprintf("Slow Trial 1 omega_{max} = %5f rad/s", max(srav1*pi/180)) ; sav2 = sprintf("Slow Trial 2 omega_{max} = %5f rad/s", max(srav2*pi/180)) ;
sav3 = sprintf("Slow Trial 3 omega_{max} = %5f rad/s", max(srav3*pi/180)) ;
fav1 = sprintf("Fast Trial 1 omega_{max} = %5f rad/s", max(frav1*pi/180)) ; fav2 = sprintf("Fast Trial 2 omega_{max} = %5f rad/s", max(frav2*pi/180)) ; 
fav3 = sprintf("Fast Trial 3 omega_{max} = %5f rad/s", max(frav3*pi/180)) ; 
legend(sav1, sav2, sav3, fav1, fav2, fav3, 'location', 'best') ;
grid on
hold off


% Comparing the Two Tests 
% Left Side Comparison
sMlav = mean([max(slav1*pi/180) max(slav2*pi/180) max(slav3*pi/180)]) ; sSDlav = std([max(slav1*pi/180) max(slav2*pi/180) max(slav3*pi/180)]) ;
fMlav = mean([max(flav1*pi/180) max(flav2*pi/180) max(flav3*pi/180)]) ; fSDlav = std([max(flav1*pi/180) max(flav2*pi/180) max(flav3*pi/180)]) ;
errhigh = [sSDlav, fSDlav] ; errlow = [sSDlav, fSDlav] ;
Y = [sMlav fMlav]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Max Left Shank Angular Velocity') ;
ylabel('Angular Velocity [rad/s]') ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off

% Right Side Comparison
sMrav = mean([max(srav1*pi/180) max(srav2*pi/180) max(srav3*pi/180)]) ; sSDrav = std([max(srav1*pi/180) max(srav2*pi/180) max(srav3*pi/180)]) ;
fMrav = mean([max(frav1*pi/180) max(frav2*pi/180) max(frav3*pi/180)]) ; fSDrav = std([max(frav1*pi/180) max(frav2*pi/180) max(frav3*pi/180)]) ;
errhigh = [sSDrav, fSDrav] ; errlow = [sSDrav, fSDrav] ;
Y = [sMrav fMrav]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Max Right Shank Angular Velocity') ;
ylabel('Angular Velocity [rad/s]') ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off

%}

%% vii) Shank Angular Acceleration
%{
% Left Slow Shank Angular Acceleration Calculations
[~,~,slthetas1] = angles(sLGRT1(:,:), sLLK1(:,:), sLLA1(:,:)) ;
slaa1 = finitediffaa(slthetas1) ; slaa1 = slaa1(stzero1 >= 0) ; 
slfaa1 = filtera(slaa1) ;
[~,~,slthetas2] = angles(sLGRT2(:,:), sLLK2(:,:), sLLA2(:,:)) ;
slaa2 = finitediffaa(slthetas2) ; slaa2 = slaa2(stzero2 >= 0) ; 
slfaa2 = filtera(slaa2) ;
[~,~,slthetas3] = angles(sLGRT3(:,:), sLLK3(:,:), sLLA3(:,:)) ;
slaa3 = finitediffaa(slthetas3) ; slaa3 = slaa3(stzero3 >= 0) ;
slfaa3 = filtera(slaa3) ;

% Left Fast Shank Angular Acceleration Calculations
[~,~,flthetas1] = angles(fLGRT1(:,:), fLLK1(:,:), fLLA1(:,:)) ;
flaa1 = finitediffaa(flthetas1) ; flaa1 = flaa1(ftzero1 >= 0) ; 
flfaa1 = filtera(flaa1) ;
[~,~,flthetas2] = angles(fLGRT2(:,:), fLLK2(:,:), fLLA2(:,:)) ;
flaa2 = finitediffaa(flthetas2) ; flaa2 = flaa2(ftzero2 >= 0) ; 
flfaa2 = filtera(flaa2) ;
[~,~,flthetas3] = angles(fLGRT3(:,:), fLLK3(:,:), fLLA3(:,:)) ;
flaa3 = finitediffaa(flthetas3) ; flaa3 = flaa3(ftzero3 >= 0) ;
flfaa3 = filtera(flaa3) ;

% Right Slow Shank Angular Acceleration Calculations
[~,~,srthetas1] = angles(sRGRT1(:,:), sRLK1(:,:), sRLA1(:,:)) ;
sraa1 = finitediffaa(srthetas1) ; sraa1 = sraa1(stzero1 >= 0) ; 
srfaa1 = filtera(sraa1) ;
[~,~,srthetas2] = angles(sRGRT2(:,:), sRLK2(:,:), sRLA2(:,:)) ;
sraa2 = finitediffaa(srthetas2) ; sraa2 = sraa2(stzero2 >= 0) ; 
srfaa2 = filtera(sraa2) ;
[~,~,srthetas3] = angles(sRGRT3(:,:), sRLK3(:,:), sRLA3(:,:)) ;
sraa3 = finitediffaa(srthetas3) ; sraa3 = sraa3(stzero3 >= 0) ;
srfaa3 = filtera(sraa3) ;

% Right Fast Shank Angular Accleration Calculations
[~,~,frthetas1] = angles(fRGRT1(:,:), fRLK1(:,:), fRLA1(:,:)) ;
fraa1 = finitediffaa(frthetas1) ; fraa1 = fraa1(ftzero1 >= 0) ; 
frfaa1 = filtera(fraa1) ;
[~,~,frthetas2] = angles(fRGRT2(:,:), fRLK2(:,:), fRLA2(:,:)) ;
fraa2 = finitediffaa(frthetas2) ; fraa2 = fraa2(ftzero2 >= 0) ; 
frfaa2 = filtera(fraa2) ;
[~,~,frthetas3] = angles(fRGRT3(:,:), fRLK3(:,:), fRLA3(:,:)) ;
fraa3 = finitediffaa(frthetas3) ; fraa3 = fraa3(ftzero3 >= 0) ;
frfaa3 = filtera(fraa3) ;


% Left Side Shank Angular Acceleration
figure()
plot(stzero1(stzero1 >= 0), slfaa1*pi/180, ":.", 'linewidth', 1.25) ;
hold on
plot(stzero2(stzero2 >= 0), slfaa2*pi/180, ":.", 'linewidth', 1.25) ;
plot(stzero3(stzero3 >= 0), slfaa3*pi/180, ":.", "linewidth", 1.25) ;
plot(ftzero1(ftzero1 >= 0), flfaa1*pi/180, 'linewidth', 1.25) ;
plot(ftzero2(ftzero2 >= 0), flfaa2*pi/180, 'linewidth', 1.25) ;
plot(ftzero3(ftzero3 >= 0), flfaa3*pi/180, 'linewidth', 1.25) ;
title("Left Shank Angular Acceleration")
xlabel("Time [s]") ;
ylabel("Angular Acceleration [rad/s^2]") ;
saa1 = sprintf("Slow Trial 1 alpha_{max} = %5f rad/s^2", max(slaa1*pi/180)) ; saa2 = sprintf("Slow Trial 2 alpha_{max} = %5f rad/s^2", max(slaa2*pi/180)) ;
saa3 = sprintf("Slow Trial 3 alpha_{max} = %5f rad/s^2", max(slaa3*pi/180)) ;
faa1 = sprintf("Fast Trial 1 alpha_{max} = %5f rad/s^2", max(flaa1*pi/180)) ; faa2 = sprintf("Fast Trial 2 alpha_{max} = %5f rad/s^2", max(flaa2*pi/180)) ; 
faa3 = sprintf("Fast Trial 3 alpha_{max} = %5f rad/s^2", max(flaa3*pi/180)) ; 
legend(saa1, saa2, saa3, faa1, faa2, faa3, 'location', 'best') ;
grid on
hold off

% Right Side Shank Angular Acceleration
figure()
plot(stzero1(stzero1 >= 0), srfaa1*pi/180, ":.", 'linewidth', 1.25) ;
hold on
plot(stzero2(stzero2 >= 0), srfaa2*pi/180, ":.", 'linewidth', 1.25) ;
plot(stzero3(stzero3 >= 0), srfaa3*pi/180, ":.", "linewidth", 1.25) ;
plot(ftzero1(ftzero1 >= 0), frfaa1*pi/180, 'linewidth', 1.25) ;
plot(ftzero2(ftzero2 >= 0), frfaa2*pi/180, 'linewidth', 1.25) ;
plot(ftzero3(ftzero3 >= 0), frfaa3*pi/180, 'linewidth', 1.25) ;
title("Right Shank Angular Acceleration")
xlabel("Time [s]") ;
ylabel("Angular Acceleration [rad/s^2]") ;
saa1 = sprintf("Slow Trial 1 alpha_{max} = %5f rad/s^2", max(sraa1*pi/180)) ; saa2 = sprintf("Slow Trial 2 omega_{max} = %5f rad/s^2", max(sraa2*pi/180)) ;
saa3 = sprintf("Slow Trial 3 alpha_{max} = %5f rad/s^2", max(sraa3*pi/180)) ;
faa1 = sprintf("Fast Trial 1 alpha_{max} = %5f rad/s^2", max(fraa1*pi/180)) ; faa2 = sprintf("Fast Trial 2 omega_{max} = %5f rad/s^2", max(fraa2*pi/180)) ; 
faa3 = sprintf("Fast Trial 3 alpha_{max} = %5f rad/s^2", max(fraa3*pi/180)) ; 
legend(saa1, saa2, saa3, faa1, faa2, faa3, 'location', 'best') ;
grid on
hold off


% Comparing the Two Tests 
% Left Side Comparison
sMlaa = mean([max(slaa1*pi/180) max(slaa2*pi/180) max(slaa3*pi/180)]) ; sSDlaa = std([max(slaa1*pi/180) max(slaa2*pi/180) max(slaa3*pi/180)]) ;
fMlaa = mean([max(flaa1*pi/180) max(flaa2*pi/180) max(flaa3*pi/180)]) ; fSDlaa = std([max(flaa1*pi/180) max(flaa2*pi/180) max(flaa3*pi/180)]) ;
errhigh = [sSDlaa, fSDlaa] ; errlow = [sSDlaa, fSDlaa] ;
Y = [sMlaa fMlaa]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Max Left Shank Angular Acceleration') ;
ylabel('Angular Velocity [rad/s^2]') ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off

% Right Side Comparison
sMraa = mean([max(sraa1*pi/180) max(sraa2*pi/180) max(sraa3*pi/180)]) ; sSDraa = std([max(sraa1*pi/180) max(sraa2*pi/180) max(sraa3*pi/180)]) ;
fMraa = mean([max(fraa1*pi/180) max(fraa2*pi/180) max(fraa3*pi/180)]) ; fSDraa = std([max(fraa1*pi/180) max(fraa2*pi/180) max(fraa3*pi/180)]) ;
errhigh = [sSDraa, fSDraa] ; errlow = [sSDraa, fSDraa] ;
Y = [sMraa fMraa]' ; X = categorical(["Slow Trials", "Fast Trials"]) ;
figure() 
bar(X, Y) ;
hold on 
title('Max Right Shank Angular Acceleration') ;
ylabel('Angular Velocity [rad/s^2]') ;
er = errorbar(X,Y,errlow, errhigh) ;
er.Color = [0 0 0] ; er.LineStyle = 'none' ;
hold off



%}



%% Function Files
% Stride Length
function sl = stridelength(z,t)
    z1 = z(:,1) ; z3 = z(:,3) ; ts = 0.01 ;
    [~, idx] = findpeaks(-z3, t, 'MinPeakProminence', 9.0) ;
   
    idx = idx/ts ; idx = round(idx,0) ;
 
    sl = z1(idx) ;
    test = sl ; 
    ii = 1 ;
    while (1) 
                sl(ii) = test(ii+1) - test(ii) ;
                ii = ii + 1 ;
                if ii == length(test)
                    sl = sl(1:ii-1) ;
                    break
                end
    end
    sl = mean(sl) ;
    sl = sl/(10*2.54*12) ;
    end
    %}   

% Velocity Finite Difference
function [v, v1] = finitediffv(x)
    n = length(x) ;
    v1 = zeros(n,3) ;
    v = zeros(n,1) ;
    ts = 0.01 ;
    for i = 2:n-1
        v1(i,1) = (x(i+1,1)-x(i-1,1))/(2*ts) ;
        v1(i,2) = (x(i+1,2)-x(i-1,2))/(2*ts) ;
        v1(i,3) = (x(i+1,3)-x(i-1,3))/(2*ts) ;
    end
    for i = 1:length(v)
    v(i,1) = sqrt(v1(i,1)^2 + v1(i,2)^2 + v1(i,3)^2) ;
    end
end

% Acceleration Finite Difference
function [a, a1] = finitediffa(V)
    s = length(V) ;
    a = zeros(s,1) ;
    a1 = zeros(s,3) ;
    ts = 0.01 ;
    for i = 3:s-2
        a1(i,1) = (V(i+1,1)-V(i-1,1))/(2*ts) ;
        a1(i,2) = (V(i+1,2)-V(i-1,2))/(2*ts) ;
        a1(i,3) = (V(i+1,3)-V(i-1,3))/(2*ts) ;
    end
    for i = 1:length(a)
    a(i,1) = sqrt(a1(i,1)^2 + a1(i,2)^2 + a1(i,3)^2) ;
    end
end

% Acceleration Filtering
function af = filtera(A)
    nn = length(A) ;
    af = zeros(nn,1) ;
    for i = 11:nn-11
        af(i) = mean(A(i-10:i+10)) ;
    end
end

% Joint Angles
function [thetak, thetat, thetas] = angles(hip,knee,ankle)
    nnn = length(hip) ;
    thetak = zeros(nnn,1) ; thetat = zeros(nnn,1) ; thetas = zeros(nnn,1) ;
    for i = 1:nnn
        ax(i,1) = (hip(i,1) - knee(i,1)) ; az(i,1) = (hip(i,3) - knee(i,3)) ;
        a(i,1) = sqrt(ax(i,1)^2 + az(i,1)^2) ;
        bx(i,1) = (knee(i,1) - ankle(i,1)) ; bz(i,1) = (knee(i,3) - ankle(i,3)) ;
        b(i,1) = sqrt(bx(i,1)^2 + bz(i,1)^2) ;
        cx(i,1) = (hip(i,1) - ankle(i,1)) ; cz(i,1) = (hip(i,3) - ankle(i,3))  ;
        c(i,1) = sqrt(cx(i,1)^2 + cz(i,1)^2) ;

        thetak(i,1) = acosd((a(i,1)^2 + b(i,1)^2 - c(i,1)^2)/(2*a(i,1)*b(i,1))) ;
        if hip(i,1) <= knee(i,1)
            thetat(i,1) = 180 - atand(az(i,1)/-ax(i,1)) ;
        elseif hip(i,1) > knee(i,1)
            thetat(i,1) = atand(az(i,1)/ax(i,1)) ;
        end
        if knee(i,1) <= ankle(i,1)
            thetas(i,1) = 180 - atand(bz(i,1)/-bx(i,1)) ;
        elseif knee(i,1) > ankle(i,1)
            thetas(i,1) = atand(bz(i,1)/bx(i,1)) ;
        end
        if thetat(i,1) <= thetas(i,1) && ankle(i,1) >= knee(i,1)
            thetak(i,1) = 360 - thetak(i,1) ;
        end

    end
    %plot(thetak) ;
    
    % thigh angle < shank angle & anklex > kneex
end

% Angular Velocity Finite Difference
function av = finitediffav(theta)
    n = length(theta) ;
    v = zeros(n,1) ;
    ts = 0.01 ;
    for i = 2:n-1
        if theta(i+1) == 0
            v(i-1,1) = 0 ;
        end
        v(i,1) = (theta(i+1,1)-theta(i-1,1))/(2*ts) ;
    end
    av = v ;
end

% Angular Acceleration Finite Difference
function aa = finitediffaa(theta)
    n = length(theta) ;
    a = zeros(n,1) ;
    ts = 0.01 ;
    for i = 2:n-1
        if theta(i+1) == 0
            a(i-1,1) = 0 ;
        end
        a(i,1) = (theta(i+1,1)-2*theta(i,1)+theta(i-1,1))/(ts^2) ;
    end
    aa = a ;
end