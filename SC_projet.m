%% SC1 
clear all; clc;close all;
load("data_1v_4-09_100hz.mat", "Vm", "servo", "omega_c", "tsimu");

%Variables
ms  = 0.064;
Js  = 4.129e-6;
rs  = 0.0127;
g   = 9.81;
km  = 0.0076776;
kt  = 0.0076830;
nm  = 0.69;
Jm  = 3.9001e-7;
Jeq = 0.0017728; % Jeq = Jm + N^2*Jc
ng  = 0.9;
Kg  = 70;
rarm = 0.0254;
L   = 0.4254;
N = ng*Kg;

omega_c_rad = deg2rad(omega_c);
A_moy = mean(omega_c_rad(114:end));

gn = log((A_moy - omega_c_rad(101:108))./omega_c_rad(101:108));
c = [tsimu(101:108) 1.^tsimu(101:108)];
a = pinv(c)*real(gn);
alpha = a(1);
beta = a(2);
lambda = -alpha;
t0 = beta./lambda;
omega_c_lisse = A_moy./(1 + exp(a(1).*tsimu(1:end-1) + a(2)));

%----Méthode des moindres carrés----%
Y = omega_c_lisse(1:end-1);
X = [Vm(1:end-2) diff(omega_c_lisse)./diff(tsimu(1:end-1))];

R = X' * X;
P = X' * Y;

Rinv = inv(R);
A_moindres_carres = Rinv*P;

K = A_moindres_carres(1);
tau = -A_moindres_carres(2); % Verifier le signe de Tau --> Équivalent a  :1 - alpha_sortie1 et comme A(2) négatif alors tau positif 
% Calculer Rm et Beq à partir de la fonction de transfert standard
Beq = (Jeq./tau) * ((nm*kt*Kg*ng./K)-nm*kt*km*(Kg.^2)*ng)./(nm*kt*Kg*ng./K);
Rm = (1./Beq)*(Kg*ng*kt*nm./K-nm*kt*km*(Kg.^2)*ng);

%% SC-2 R^2 et erreur RMS
num = [K];
den = [tau 1];
tfMC = tf(num, den);
ySim = lsim(tfMC, Vm, tsimu); % systeme obtenu analytiquement

E_squared = sum((ySim(1:end-1) - omega_c_rad).^2); 
N = length(tsimu);
RMSE = sqrt(1/N*E_squared);

y_moy = 1/N*sum(deg2rad(omega_c));
R_squared = sum((ySim(1:end-1) - y_moy).^2)/(sum((deg2rad(omega_c) - y_moy).^2));

%% SC-3 Reponse de la fct transfert Gcm(s) à l'échelon Vm et donnees experimentales 

num = [(Kg*ng*nm*kt)./(Rm*Beq + nm*kt*km*Kg.^2*ng)];
den = [(Rm*Jeq)./(Rm*Beq+nm*kt*km*Kg.^2*ng) 1];
FT = tf(num, den); % des moindres carrés 
[y_Gcm, x_Gcm] = lsim(FT, Vm,tsimu);

figure
plot(x_Gcm, y_Gcm);
hold on; 
plot(tsimu(1:end-1), deg2rad(omega_c));
title("Reponses à Vm de la FCT G_c_m(s) comparés aux données expérimentales")
xlabel('Temps (secondes)')
ylabel('Vitesse (rad/s)')
legend('FT G_C_M(s)', 'Données expérimentales')






