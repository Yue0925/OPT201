% Prog principal
global A B L

A = 1;
B = -1;

L = [0.7, 0.5, 0.3, 0.2, 0.5]'; % longeur des barres fixées

xy = [0.2   0.4   0.6   0.8   ...
      -1.0   -1.5   -1.5   -1.3]'; % cas test 2a ok

%xy = [0.2   0.4   0.6   0.8   ...
%      1.0   1.5   1.5   1.3]'; % cas test 2b ok
      
%xy = [0.2   0.4   0.6   0.8   ...
%      -1.0   -1.5   1.5   -1.3]'; % cas test 2c ok

%xy = [0.2   0.4   0.6   0.8   ...
%      1.0   -1.2   1.5   -1.3]'; % cas test 2d non converge atteint max itération
      

lm = [0.01, 0.01, 0.01, 0.01, 0.01]'; 
[e, c, g, a, hl, indic] = chs(4, xy, lm);
% calculer lambda par moindre carrées
lm = -a'\g; 
%lm % vérifié ok


% vérifier le simulateur (fait pour tp1)
if 0 % 1 ou 0
  printf("VERIFICATION LE SIMULATEUR \n");
  for indic = 1:5
    [e, c, g, a, hl, indic] = chs(indic, xy, lm);
    printf("\n\n\n");
  endfor
endif


options.tol = [1e-8 1e-8]; %précision
options.maxit = 500; %nb max itérations


[x, lm, info] = sqp(@chs, xy, lm, options);
info

