clc; close all; clear all;

% Chargement des donnes du banc d'essai
load("data_1v_4-09_100hz.mat", "Vm", "servo", "omega_c", "tsimu");

% Variables
ms  = 0.064;
Js  = 4.1290e-6; 
rs  = 0.0127;
g   = 9.81;
km  = 0.0076776;
kt  = 0.0076830;
nm  = 0.69;
Jm  = 3.9001e-7;
Jeq = 0.0017728; % Jeq = Jm*kg.^2*ng + Jc
ng  = 0.9;
Kg  = 70;
rarm= 0.0254;
L   = 0.4254;
Rm  = 3.6762;
Beq = 0.0073; % Beq = Bm*kg.^2*ng + Bc
N = ng*Kg;

%% SM
kbb = 5*g*rarm/(L*7);

% SM-3, valeur pour Simulink - Gcm.slx - Gsc.slx
K_SM_3 = nm*kt*Kg*ng;
Gcm_den_s1_SM_3 = 1/(Rm*Beq + nm*kt*km*Kg.^2*ng);
Gcm_den_s2_SM_3 = 1/(Rm*Jeq);
Time_Vm = [tsimu Vm];
Time_servo = [tsimu servo];

% SM-5, matrice 4x4 sans la variable d'etat im
coefA = -(Rm*Beq+nm*kt*km*Kg.^2*ng)/(Rm*Jeq);
A1 = [0 1 0 0 ;
      0 0 kbb 0;
      0 0 0 1;
      0 0 0 coefA];
	  
B1 = [0; 0; 0; (nm*kt*Kg*ng)/(Rm*Jeq)];

C1 = [1 0 0 0;
      0 0 1 0];

D1 = [0; 0];

% valeurs propres
Vp = eig(A1);

% SM-6
% FTBO entre la tension en entrée (Vm) et l'angle du moteur (oc)
Gcm_num = (nm*kt*Kg*ng)/(Rm*Jeq);
Gcm_den = [1 (Rm*Beq + nm*kt*km*(Kg.^2)*ng)/(Rm*Jeq) 0]; % division par (Rm*Jeq) pour avoir forme standard
Gcm = numden2system(Gcm_num, Gcm_den);
% poles de la FTBO
%disp(Gcm.p) 

% FTBO entre l'angle du moteur (oc) et la position de la charge (x) - GOOD
Gsc_num = (5*g*rarm)/(L*7);
Gsc_den = [1 0 0]; % division par (L*7) pour avoir forme standard
Gsc = numden2system(Gsc_num, Gsc_den);
% poles de la FTBO
%disp(Gsc.p) 

% SM-7, FTBO de (x)/(Vm) - GOOD 
Gsm_num = (5*g*rarm*nm*kt*Kg*ng)/(Rm*Jeq*7*L);
Gsm_den = [1 (((Rm*Beq+nm*kt*km*Kg.^2*ng)*7*L)/(Rm*Jeq*7*L)) 0 0 0]; % division par (Rm*Jeq*5*L) pour avoir forme standard
Gsm = numden2system(Gsm_num, Gsm_den);
% poles de la FTBO
%disp(Gsm.p)


%% SC 
% Fait de son cote pour eviter de saturer le workspace
% Valeur tirées de SC
% LES VALEURS NE CHANGERONT PAS, Rm et Beq est du côté moteur/électrique
Rm  = 3.6267;
Beq = 0.0073;


%% SI
% Lieu des racines SI-1 (à garder en commentaire si pas utiliser)
%figure
%rlocus(Gcm) 
%title('Lieu des racines de la fonction G_c_m')	

% SI-2
Kcrit = 4.94; % a)
ts_b = 4/16; % b) Valeur 16 = -zeta*wn pris appartir du rlocus  
ts_c = 4./(16.019);  % c) À partir de la règle 7 Rlocus a la main

% SI-2 d) Valeur TS avec comparaison avec forme standard
%Gcm_BF = feedback(Kcrit*Gcm,1);
%denGcm_BF = [Gcm_BF.Denominator{1, 1}(1) Gcm_BF.Denominator{1, 1}(2) Gcm_BF.Denominator{1, 1}(3)]/Gcm_BF.Denominator{1, 1}(1);
%numGcm_BF = Gcm_BF.Numerator{1, 1}(3)/Gcm_BF.Denominator{1, 1}(1);
%tf_ = tf(numGcm_BF,denGcm_BF);
%wn = sqrt(denGcm_BF(3));
%zeta = denGcm_BF(2)/(2*wn);
%ts_d = 4/wn*zeta;

% SI-2 F)
i = (Kcrit:0.1:100);
p_ts_cst = rlocus(Gcm.TF,i);
phi_SI2 = pi - angle(p_ts_cst(1,:));
Mp_SI2 = 100*exp(-pi./tan(phi_SI2));

%figure(2)
%plot(phi,Mp)
%title('MP-phi pour gain avec ts identique')
%xlabel('phi (rad)')
%ylabel('MP(%)')

%SI-3
Kint_Rac = 7.72; % trouver avec rlocus à damping = 0.8
zeta_SI3 = 0.8;

Kint_Calc = ((Beq*Rm+nm*kt*Kg*Kg*ng*km)/(Rm*Jeq*2*zeta_SI3)).^2 * (Rm*Jeq)/(nm*kt*Kg*ng);
Gcm_Kint = numden2system(Kint_Calc*Gcm.num,Gcm.den);
%disp(Gcm_Kint.BF);

% SI-4
%disp(Gcm_Kint.p);

figure
rlocus(Gcm.TF)
hold on
plot(real(roots(Gcm_Kint.BF.Denominator{1,1})),  imag(roots(Gcm_Kint.BF.Denominator{1,1})),'Diamond','Color','r')
hold on
plot(real(roots(Gcm_Kint.BF.Denominator{1,1})),  -imag(roots(Gcm_Kint.BF.Denominator{1,1})),'Diamond','Color','r')
title('Lieu des racines de la fonction G_c_m')
legend('G_c_m','Pôles désirés')

% SI-5
figure
margin(Gcm.TF)
title('Diagramme de Bode de G_c_m')

% SI-6
zeta_SI6 = 0.5*sqrt(tand(Gcm_Kint.Pm)*sind(Gcm_Kint.Pm));
PM_Calc = atand((2*zeta_SI6)/(sqrt(sqrt(1+4*zeta_SI6.^4)-2*zeta_SI6.^2)));

% SI-7
Aint = (A1 - B1*Kint_Calc*C1(2,:));
Bint = B1*Kint_Calc;
Cint = C1;
Dint = D1;
vp_Aint = eig(Aint);

% SI-8
[numInt, denInt] = ss2tf(Aint,Bint,Cint,Dint,1);
Gsm_int = numden2system(numInt(1,:),denInt);
Gcm_int = numden2system(numInt(2,:),denInt);
%poles_FT_int = roots(denInt); % Valeur propre et poles des FT sont identiques

% SI-9 
%La FT Gsm_int est classe 2 et la FT Gsm etait classe 3
Gsm_intTF = tf(Gsm_int.num, Gsm_int.den); % AJOUT
% SI-10
figure()
rlocus(Gsm_int.TF)
title('Lieu des racines de la G_s_m de la boucle interne')

%% SE
% SE-1
%flag_in = 1; % VERIFIER SI LES 2 SIGNAUX D'ENTREE DE SIMULINK SONT PERTINENTS
out = sim('SE_1');
t = [out.x_sphere.Time];
u = ones(size(t));

% SE-2
y = lsim(Gsm_int.TF ,Time_Vm(:,2),Time_Vm(:,1));

figure()
plot(Time_Vm(:,1),y)
title('FT Gsm-int entre Vm - version linéaire')
ylabel('Position (m)')
xlabel('Time(s)')

% Entre en echelon
y_lin_ech = lsim(Gsm_int.BF, u, t);

figure()
plot(t, y_lin_ech)
hold on
plot(t, out.x_sphere.Data)
title('Réponse à l échelon du système')
xlabel('temps (s)')
ylabel('Position (m)')
legend('Linéaire', 'Non-Linéaire')
%plot(t, 0.98*y_lin_ech(end)*[1:1], 'r--', 'linewidth', 2);
%plot(t, 1.02*y_lin_ech(end)*[1:1], 'r--', 'linewidth', 2);

%Erreur_lin_nlin = y_lin_ech-out.x_sphere;
%figure()
%plot(t, (Erreur_lin_nlin)
%title('Erreur entre Simulink (non-lin) et lsim (lin)')

%% Conception AvPh par critères temporels (lieu des racines)
% Conception initiale
Mp_ini_t = 5;
ts_ini_t = 4;
tr_ini_t = 2;
tp_ini_t = 3; 

phi_t = atan(-pi./(log(Mp_ini_t./100)));
zeta_t = cos(phi_t);

% calcul des omega n --> prendre la plus grande valeur
wn_ts_t = 4./(zeta_t*ts_ini_t);
wn_tr_t = (1+1.1*zeta_t+1.4*zeta_t.^2)./tr_ini_t;
wn_tp_t = pi./(tp_ini_t*sqrt(1-zeta_t.^2));
wn_t = max([wn_ts_t wn_tr_t wn_tp_t]); % omega maximum calculé ici

% pôles désirés 
s_des_t = -zeta_t*wn_t + 1i*(wn_t*sqrt(1-zeta_t.^2));

% Lieu des racines pour voir si un gain K peut ajuster (garder commenter)
figure
rlocus(Gsm_int.TF)
hold on;
plot(real(s_des_t), imag(s_des_t),'p');
hold on;
plot(real(s_des_t), -imag(s_des_t),'p')
title('Lieu des racines de la G_s_m')
% CTRL + MAJ + R enleve les commentaires
% UN GAIN K NE COMPENSE PAS 

pol_t  = polyval(Gsm_int.den, s_des_t);
ph_t   = -rad2deg((2*pi+angle(pol_t))); 
dphi_t = -180 - ph_t + 7.6; 

alpha_t = 180 - rad2deg(phi_t);
zphi_t  = (alpha_t + dphi_t)./2;
pphi_t  = (alpha_t - dphi_t)./2;

z_t = real(s_des_t) - imag(s_des_t)./tand(zphi_t);
p_t = real(s_des_t) - imag(s_des_t)./tand(pphi_t);

Ka_t = 1./abs(((s_des_t - z_t)./(s_des_t - p_t))*polyval(Gsm_int.num, s_des_t)./polyval(Gsm_int.den, s_des_t));

% Fonction de transfert Compensateur AvPh
numGa_t = 1.19*Ka_t*[1 -z_t];
denGa_t = [1 -p_t];
ftGa_t  = numden2system(numGa_t, denGa_t);

% Fonction de transfert Gsm(s) * Ga(s)
Gext_t = numden2system(conv(Gsm.num, ftGa_t.num),conv(Gsm.den, ftGa_t.den));

% Lieu des racines avant tune
figure
rlocus(Gext_t.TF)
hold on;
plot(real(s_des_t), imag(s_des_t),'p');
hold on;
plot(real(s_des_t), -imag(s_des_t),'p')
title('Lieu des racines de la G_s_m(s)*G_a(s)')


% critères de performances finales 
figure
margin(Gext_t.TF) % Première observation : pas assez de phase manque 3 deg

info = stepinfo(Gext_t.BF);

%%%%%%%%%%%%%%% Temporelle section en haut 

%%%% Frequentielle %%%%
% Définition de PM et BW
PM_f  = deg2rad(45);  
BW_f  = 2.3;          
err_f = 0.005;       

% Calcul du coefficient d'amortissement zeta basé sur la marge de phase
zeta_f = 0.5 * sqrt(tan(PM_f) * sin(PM_f));

% Calcul du numérateur et du dénominateur pour la fréquence à laquelle la marge est mesurée
wgNum_f = sqrt(sqrt(1 + 4 * zeta_f^4) - (2 * zeta_f^2));
wgDen_f = sqrt((1 - 2 * zeta_f^2) + sqrt(4 * zeta_f^4 - 4 * zeta_f^2 + 2));
wg_f = BW_f * (wgNum_f / wgDen_f);    
%+ 0.96;  % Fréquence de croisement de la bande passante

% Calcul du gain requis pour obtenir la bande passante désirée
[mag_f, phase_f] = bode(Gsm.TF, wg_f);  
Kdes_f = 1 / mag_f  ;             

% Création d'une nouvelle fonction de transfert avec le gain correctif
Gext_temp_f = numden2system(Kdes_f*Gsm.num, Gsm.den);

% Calcul de la marge de phase souhaitée en degrés
PMdes_f = rad2deg(PM_f);  % normal de reprendre PM initial aulieu de celui de G_ext_f?

% Calcul du déphasage nécessaire pour atteindre la marge de phase désirée
deltaPhi_f = PMdes_f - Gext_temp_f.Pm + 5; %- 2.8;

% Calcul de alpha pour le correcteur d'avance de phase
alpha_f = (1 - sind(deltaPhi_f)) / (1 + sind(deltaPhi_f));

% Calcul du paramètre T du correcteur d'avance de phase
T_f = 1 / (wg_f * sqrt(alpha_f));
z_f = -1 / T_f;
p_f = -1 / (alpha_f * T_f);  

% Calcul du gain du correcteur d'avance de phase
Ka_f = (Kdes_f / sqrt(alpha_f)); %* 1.21

% Définition des coefficients du numérateur et du dénominateur du correcteur
numAvPh_f = Ka_f * [1 1/T_f];
denAvPh_f = [1 1 /(alpha_f*T_f)];

% Création de la fonction de transfert du correcteur d'avance de phase
ftGa_f = numden2system(numAvPh_f, denAvPh_f);

% Création de la fonction de transfert totale (système corrigé)
Gext_f = numden2system(conv(ftGa_f.num, Gsm_int.num),conv(ftGa_f.den, Gsm_int.den));

figure;
margin(Gext_f.TF)

figure % Cas de la réponse à l’échelon 
t = [1:0.01:20]';
u = ones(size(t));  % Échelon unitaire 
% ou FTBF = feedback(FTBO,1) 
ybf = lsim(Gext_f.BF, u, t); 
plot(t,ybf,'b', 'linewidth', 2) 
grid on 
hold on 
plot([t(1); t(end)], 0.98*ybf(end)*[1;1], 'r', 'linewidth', 2) 
plot([t(1); t(end)], 1.02*ybf(end)*[1;1], 'r', 'linewidth', 2) 

figure;
u2 = 6*ones(size(t));
ybf2 = lsim(Gext_f.BF, u2, t); 
plot(t,ybf2,'b', 'linewidth', 2) 


% Définir la valeur de stabilisation finale
y_final = ybf(end);

% Calcul des limites de 2%
lim_sup = 1.02 * y_final;
lim_inf = 0.98 * y_final;

% Déterminer le temps où la réponse reste dans les limites
idx_stable = find(ybf >= lim_inf & ybf <= lim_sup);

% Vérifier à partir de quel indice la réponse reste stable

for i = 1:length(idx_stable)
    if all(ybf(idx_stable(i):end) >= lim_inf & ybf(idx_stable(i):end) <= lim_sup)
        t_stab = t(idx_stable(i)); 
        break;
    end
end

% Affichage du temps de stabilisation
fprintf('Ts (2%%) : %.2f secondes\n', t_stab);

% Ajout d'une ligne verticale sur le graphique
if ~isnan(t_stab)
    hold on;
    plot([t_stab, t_stab], [0, y_final], '--g', 'LineWidth', 2);
end