clear ; clc ; close all

%% Gait Speed Data Slow Trial 1

% Data Extraction
Slow1 = readmatrix("Slow1.csv") ;
Slow1 = Slow1(7:end,:) ;
Slow1(:,2) = Slow1(:,2) - 3.74 ;
tzero = Slow1(:,2) ;            % Time where gait is stable is t = 0
C7 = Slow1(:,3:5) ;             % 7th Cervical Vertebra
Ster = Slow1(:,6:8) ;           % Sternum 
LASIS = Slow1(:,9:11) ;         % Left Anterior Superior Illiac Spine
RASIS = Slow1(:,12:14) ;        % Right Anterior Superior Illiac Spine
LPSIS = Slow1(:,15:17) ;        % Left Posterior Superior Illiac Spine
RPSIS = Slow1(:,18:20) ;        % Right Posterior Superior Illiac Spine
LGRT = Slow1(:,21:23) ;         % Left Greater Trochanter
RGRT = Slow1(:,42:44) ;         % Right Greater Trochanter
LLK = Slow1(:,24:26) ;          % Left Lateral Knee
RLK = Slow1(:,45:47) ;          % Right Lateral Knee
LMK = Slow1(:,27:29) ;          % Left Medial Knee
RMK = Slow1(:,48:50) ;          % Right Medial Knee
LLA = Slow1(:,30:32) ;          % Left Lateral Ankle
RLA = Slow1(:,51:53) ;           % Right Lateral Ankle
LMA = Slow1(:,33:35) ;          % Left Medial Ankle
RMA = Slow1(:,54:56) ;          % Right Medial Ankle
LH = Slow1(:,36:38) ;           % Left Heel
RH = Slow1(:,57:59) ;           % Right Heel
LT = Slow1(:,39:41) ;           % Left Toe
RT = Slow1(:,60:62) ;           % Right Toe

%% i) Stride Length from Heel Markers
%{
figure(1)
plot(tzero(tzero >= 0), LH(tzero >= 0)) ;
xlabel("Time [s]") ;
ylabel('Heel Displacemet [mm]') ;




%}

%% Forward Walking Speed 
%{
vster = finitediffv(Ster(tzero >= 0)) ; vster = vster/1000 ;
vlasis = finitediffv(LASIS(tzero >= 0)) ; vlasis = vlasis/1000 ;
vrasis = finitediffv(RASIS(tzero >= 0)) ; vrasis = vrasis/1000 ;

figure(2)
%plot(tzero(tzero >= 0), Ster(tzero >= 0)) ;
%hold on 
plot(tzero(tzero >= 0), vster) ;
hold on 
plot(tzero(tzero >= 0), vlasis) ;
plot(tzero(tzero >= 0), vrasis) ;
title("Forward Walking Speed") ;
xlabel("Time [s]") ;
ylabel("Walking Speed [m/s]") ;
legend("Sternum", "LASIS", "RASIS") ;
%}

%% Heel Forward Linear Velocity

vLH = finitediffv(LH(tzero >= 0)) ; vLH = vLH/1000 ;
vRH = finitediffv(RH(tzero >= 0)) ; vRH = vRH/1000 ;
figure(3)
plot(tzero(tzero >= 0), vLH) ;
hold on
plot(tzero(tzero >= 0), vRH) ;
ylim([-1 4]) ;
title("Heel Forward Velocity") ;
xlabel('Time [s]') ;
ylabel('Heel Velocity [m/s]') ;
legend('Left Heel', 'Right Heel', 'location', 'best') ;
%}

%% Heel Forward Linear Acceleration

aLH = finitediffa(vLH) ; 
afLH = filtera(aLH) ;
aRH = finitediffa(vRH) ;
afRH = filtera(aRH) ;
figure(4) 
plot(tzero(tzero >= 0), aLH, "linewidth", 1.25) ;
hold on
plot(tzero(tzero >= 0), afLH, 'linewidth', 1.0) ;
title('Left Heel Forward Acceleration') ;
xlabel('Time [s]') ;
ylabel('Left Heel Acceleration [m/s^2]') ;
legend('Raw Data', 'Moving Avg') ;
hold off
figure(5)
plot(tzero(tzero >= 0), aRH, 'linewidth', 1.25) ;
hold on
title('Right Heel Forward Acceleration') ;
xlabel('Time [s]') ;
ylabel('Right Heel Acceleration [m/s^2]') ;
plot(tzero(tzero >= 0), afRH, 'linewidth', 1.0) ;
legend("Raw Data", 'Moving Avg') ;
%}
%% Knee Relative Joint Angle

%% Shank Angular Velocity

%% Shank Angular Acceleration







%% Function Files
function v = finitediffv(x)
    n = length(x) ;
    v = zeros(n,1) ;
    ts = 0.01 ;
    for i = 2:n-1
        v(i) = (x(i+1)-x(i-1))/(2*ts) ;
    end
end
function a = finitediffa(V)
    s = length(V) ;
    a = zeros(s,1) ;
    ts = 0.01 ;
    for i = 3:s-2
        a(i) = (V(i+1)-V(i-1))/(2*ts) ;
    end
end
function af = filtera(A)
    nn = length(A) ;
    af = zeros(nn,1) ;
    for i = 11:nn-11
        af(i) = mean(A(i-10:i+10)) ;
    end
end