clear all;
close all;
clc;

%% Plot BEM and Lub Thry Results
close all;
bemdir = '/Users/JMB/Documents/code/bem3d_ves_tube/blood95/postprocess/profiles_bem/';
lubdir = '/Users/JMB/Documents/code/bem3d_ves_tube/blood95/postprocess/profiles_lubthry/';

% SHAPE: vred = 0.95, conf = 0.90, Ca = 0.16
figure;
bemdata = load(strcat(bemdir,'shape_conf90_v1eb4.dat'));
lubdata = load(strcat(lubdir,'shape_v95_conf90_Ca016.dat'));
% mirror profile for lubrication theory
lubdata = [lubdata(:,:);
    flipud(lubdata(:,1)),flipud(-lubdata(:,2))]
% recenter
bemdata(1,:) = bemdata(1,:) - mean(bemdata(1,:));
bemdata(2,:) = bemdata(2,:) - mean(bemdata(2,:));
bemdata(3,:) = bemdata(3,:) - mean(bemdata(3,:));
lubdata(:,1) = lubdata(:,1) - mean(lubdata(:,1));
lubdata(:,2) = lubdata(:,2) - mean(lubdata(:,2));

plot(bemdata(1,:),bemdata(2,:), 'b.'); hold on;
plot(lubdata(:,1),lubdata(:,2), 'k-'); hold on;
xlabel('x');
ylabel('r');
title('v = 0.95, conf = 0.90, Ca = 0.16')

% TENSION: vred = 0.95, conf = 0.90, Ca = 0.16
figure;
bemdata = load(strcat(bemdir,'tension_conf90_v1eb4.dat'));
lubdata = load(strcat(lubdir,'tension_v95_conf90_Ca016.dat'));
% recenter
bemdata(1,:) = bemdata(1,:) - mean(bemdata(1,:));
bemdata(2,:) = bemdata(2,:) - mean(bemdata(2,:));
lubdata(:,1) = lubdata(:,1) - mean(lubdata(:,1));
lubdata(:,2) = lubdata(:,2) - mean(lubdata(:,2));

plot(bemdata(1,:),bemdata(2,:), 'b.'); hold on;
plot(lubdata(:,1),lubdata(:,2), 'k-'); hold on;
xlabel('x');
ylabel('$\sigma$');
title('v = 0.95, conf = 0.90, Ca = 0.16')


%% Plot
close all;
figure; pL = plot(log(Ca1), L1, 'b-o', log(Ca2), L2, 'r-o',...
	log(CaLub1), LLub1, 'b:o', log(CaLub2), LLub2, 'r:o');
xlabel('log Ca');
ylabel('$L$');
legend('3D BEM, v = 0.95, conf = 0.90',...
    '3D BEM, v = 0.95, conf = 0.95',...
    'lub thry, v = 0.95, conf = 0.90',...
    'lub thry, v = 0.95, conf = 0.95', 'Location', 'NorthOutside');
legend boxoff;
figure; pR = plot(log(Ca1), Rproj1, 'b-o', log(Ca2), Rproj2, 'r-o',...
    log(CaLub1), RprojLub1, 'b:o', log(CaLub2), RprojLub2, 'r:o');
xlabel('log Ca');
ylabel('$R_{proj}$');
legend('3D BEM, v = 0.95, conf = 0.90',...
    '3D BEM, v = 0.95, conf = 0.95',...
    'lub thry, v = 0.95, conf = 0.90',...
    'lub thry, v = 0.95, conf = 0.95', 'Location', 'NorthOutside');
legend boxoff;
figure; pS = plot(log(Ca1), Sigavg1, 'b-o', log(Ca2), Sigavg2, 'r-o',...
    log(CaLub1), SigavgLub1, 'b:o', log(CaLub2), SigavgLub2, 'r:o',...
    [-7,1], [0,0], 'k-');
xlabel('log Ca');
ylabel('$\sigma$');
legend('3D BEM, v = 0.95, conf = 0.90',...
    '3D BEM, v = 0.95, conf = 0.95',...
    'lub thry, v = 0.95, conf = 0.90',...
    'lub thry, v = 0.95, conf = 0.95', 'Location', 'NorthOutside');
legend boxoff;