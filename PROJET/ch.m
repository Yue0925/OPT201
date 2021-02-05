% Prog principal
global A B L
global NB_SIMUL = 0;  % nombre d'appel de simulateur au début de k-ième itération

%% TP2 cas tests
A = 1;
B = -1;
L = [0.7, 0.5, 0.3, 0.2, 0.5]'; % longeur des barres fixées

%xy = [0.2   0.4   0.6   0.8   ...
%      -1.0   -1.5   -1.5   -1.3]'; % cas test 2a ok

%xy = [0.2   0.4   0.6   0.8   ...
%      1.0   1.5   1.5   1.3]'; % cas test 2b ok
      
%xy = [0.2   0.4   0.6   0.8   ...
%      -1.0   -1.5   1.5   -1.3]'; % cas test 2c ok

xy = [0.2   0.4   0.6   0.8   ...
      1.0   -1.2   1.5   -1.3]'; % cas test 2d atteint max itération avec pas unitaire; 17 itérations avec recherche lineair

%% TP3 cas tests 3a
%L = [0.6, 0.6]';
%A = 1;
%B = 0;
%xy = [0.5 0.4]';


%% TP3 cas tests 3b
%L = [2, 1]';
%A = 1;
%B = 0;
%xy = [0.5 0.3]';


%% TP3 cas tests 3c
%L = [2, 1]';
%A = 0;
%B = -1;
%xy = [0.3 0.3]';


% calculer lambda initial par moindre carrées
lm = [];
[e, c, g, a, hl, indic] = chs(4, xy, lm);
lm = -a'\g; 
%lm % vérifié ok


% vérifier le simulateur (déjà fait pour TP1)
if 0 % 1 ou 0
  printf("VERIFICATION LE SIMULATEUR \n");
  for indic = 1:5
    [e, c, g, a, hl, indic] = chs(indic, xy, lm);
    printf("\n\n\n");
  endfor
endif


options.tol = [1e-6 1e-6]; %précision
options.maxit = 1000; %nb max itérations
options.rl = 0; % 0: avec recherche lineair; 1: le pas unitaire
options.verb = 1; %1: les sorties de chaque itération; 2: détailées de la boucle de RL sur la recherche de alpha

[x, lm, info] = sqp(@chs, xy, lm, options);
info

