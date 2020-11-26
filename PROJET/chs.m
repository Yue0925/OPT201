function [e, c, g, a, hl, indic] = chs(indic,xy,lm)
  
  global EXIT_SUCESS = 0; % sortie normale 
  global EXIT_FAILURE = 1; % paramètre(s) d’entrée non correct(s)

  % initialisation des output de fonction
  e = [];
  c = [];
  g = [];
  a = [];
  hl = [];
  
  % pré-calculs et vérifier les données en entrées
  [nn, nb, x, y, x_complet, y_complet] = pre_calcul(xy);
  
  % vérification la bonne entrées
  if mod(length(xy),2)==1
    indic = EXIT_FAILURE; % length(xy) impaire
    return
  endif
  
  if length(lm) != nb
    indic = EXIT_FAILURE; % longueur de lm différent de nb
    return
  endif 
  
  % pilote indic
  switch (indic)
    
    case 1 % tracer la chaîne
      printf("in case %d\n", indic);
      plot(x_complet, y_complet, '-b') % tracé de la chaîne
      indic = EXIT_SUCESS;
      return
      
    case 2 % calcul de e, c
      printf("in case %d\n", indic);
      [e] = cal_e(xy);
      [c] = cal_c(xy);
      indic = EXIT_SUCESS;
      return
      
    case 4 % calcul de e, c, g, a
      printf("in case %d\n", indic);
      [e] = cal_e(xy);
      [c] = cal_c(xy);
      [g] = cal_g(xy);
      verify_gradient(xy);
      [a] = cal_a(xy);
      verify_Jacobian(xy);
      indic = EXIT_SUCESS;
      return

    case 5 % calcul de hl hessien de lagrangien
      printf("in case %d\n", indic);
      [hl] = calcul_hl(xy,lm);  % calcul de hl
      %verify_Hess_Lag(xy, lm);
      indic = EXIT_SUCESS;
      return
      
    otherwise
      printf("ERROR invalid indic number %d\n ", indic);
      indic = EXIT_FAILURE;
      return
   endswitch
endfunction



% pre-calcul / traitements les donées en entrées
function [nn, nb, x, y, x_complet, y_complet] = pre_calcul(xy)
  global A B EXIT_FAILURE
  
  nn = length(xy)/2; % nombre noeuds(sauf extrémités)
  nb = nn+1; % nombre barres
  
  x = xy(1:nn); 
  y = xy(nn+1 : end);
  
  x_complet = [0;x;A]; % les abscisses en ajoutant les deux extrémités
  y_complet = [0;y;B]; % les ordonnées en ajoutant les deux extrémités
endfunction



% fonction qui calcule e énergie potentielle
function [e] = cal_e(xy)
  global L
  [nn, nb, x, y, x_complet, y_complet] = pre_calcul(xy);
  e = 0;
  
  for i = 2:nb+1
    e = e + L(i-1) * (y_complet(i) + y_complet(i-1))/2; % énergie
  endfor
  return
endfunction



% fonction qui calcule c les contraintes
function [c] = cal_c(xy)
  global L
  [nn, nb, x, y, x_complet, y_complet] = pre_calcul(xy);
  c = zeros(nb, 1);

  for i = 2:nb+1
    l_i = sqrt( (x_complet(i) - x_complet(i-1))^2 + (y_complet(i) - y_complet(i-1))^2 );
    c(i-1) = l_i^2 - L(i-1)^2; % contrainte
  endfor
  return
endfunction



% fonction qui calcule g le gradient
function [g] = cal_g(xy)
  global L
  [nn, nb, x, y, x_complet, y_complet] = pre_calcul(xy);
  
  g = zeros(2*nn, 1); % gradient de xy
  for i =1:nb-1
    g(i+nn) = (L(i) + L(i+1))/2;
  endfor
  return
endfunction



% fonction qui calcule a Jacobians des contraintes
function [a] = cal_a(xy)
  global A B
  [nn, nb, x, y, x_complet, y_complet] = pre_calcul(xy);
  a = sparse(nb,2*nn);
  
  diff2 = [2*xy(1)]; %2(x_1-x_0)
  for i = 2:nn
    diff2(end+1) = 2*(xy(i) - xy(i-1));
  endfor
  diff2(end+1) = 2*(A-xy(nn)); % 2(A-x_n)
  
  diff2(end+1) = 2*xy(nn+1); % 2(y_1 - y_0)
  for i = nn+2:2*nn
    diff2(end+1) = 2*(xy(i) - xy(i-1));
  endfor
  diff2(end+1) = 2*(B-xy(2*nn));

  a(1,1) = diff2(1);
  a(1,1+nn) = diff2(1+nb);
  for i=2:nn;
    a(i,i) = diff2(i);
    a(i,i-1) = -diff2(i);
    a(i,nn+i) = diff2(nb+i);
    a(i,nn+i-1) = -diff2(nb+i);
  endfor
  a(nb,nn) = -diff2(nb);
  a(nb,2*nn) = -diff2(2*nb);
  return
endfunction



% calcul de la matrice hessien de lagrangien
function [hl] = calcul_hl(xy,lm)
  [nn, nb, x, y, x_complet, y_complet] = pre_calcul(xy);
  hl = sparse(2*nn,2*nn); 
  
  hl(1,1) = 2*(lm(1)+lm(2));
  hl(nn+1,nn+1) = 2*(lm(1)+lm(2));
  for i = 2:nn
    hl(i,i) = 2*(lm(i)+lm(i+1));
    hl(i,i-1) = -2*lm(i);
    hl(i-1,i) = -2*lm(i);
    hl(i+nn,i+nn) = 2*(lm(i)+lm(i+1));
    hl(i+nn,i-1+nn) = -2*lm(i);
    hl(i-1+nn,i+nn) = -2*lm(i);
  endfor
  return
endfunction



% calcul dérivées de i-ième composant
function [df, t_i] = cal_df(phi, i, xy)
  t_i = sqrt(eps) * max(1, abs(xy(i)));
  e_i = double(1:length(xy) == i); % vector ligne
  df = (phi(xy + (t_i*e_i)') - phi(xy - (t_i*e_i)') ) / (2 * t_i); % d'abord transformer en vec colonne
  return
endfunction



% Vérification gradient affichage un composant par un
function verify_gradient(xy)
  printf("Vérification vecteur gradient \n");
  printf("i \t pas \t \t f'(i) \t \t DF \t \t erreur \t \t < 10e-8? \n");
  g = cal_g(xy); % gradient
  for i=1:length(xy)
    [df, t_i] = cal_df(@cal_e, i, xy);
    printf("%d \t %e \t", i, t_i);
    printf("%e \t %e \t ", g(i), df);
    erreur = abs(df - g(i)); % erreur absolue
    if abs(g(i))!=0
      erreur = erreur/abs(g(i)); % erreur relaive
      printf("%e rel \t", erreur);
    else
      printf("%e abs \t", erreur);
    endif
    if erreur < 10e-8
      printf("TRUE \n");
    else
      printf("FALSE \n");
    endif
  endfor
endfunction



% calcul Jacobians
function [jacobian] = cal_Jacob(phi, xy)
  lin = length(phi(xy));
  col = length(xy);
  jacobian = sparse(lin, col);

  for i=1:col
    [jacobian(1:end, i), t_i] = cal_df(phi, i, xy);
  endfor
  return
endfunction



function verify_Jacobian(xy)
  printf("Vérification Jacobians \n");
  [jacobian] = cal_Jacob(@cal_c, xy);
  [a] = cal_a(xy);
  erreurs = abs(jacobian - a) ./ abs(a);
  %erreurs
  printf("erreurs relatives < 10e-8? (1:TRUE; 0:FALSE) \n");
  erreurs < 10e-8
endfunction



global I=0;

% calcul hessien de lagrangien de i-ème composant
% on vérifie pour chaque composant de Jacob son n composants (en supposant dérivées premiers correctes)
function [hess_lag] = cal_Hess_Lag(xy, lm)
  global I
  nn = length(xy);
  hess_lag = sparse(nn,nn);
  
  for j=1:length(lm)
    I = j;
    tmp = sparse(nn,nn);
    for i=1:nn
      [HL_compo, t_i] = cal_df(@tmp, i, xy);
      tmp(1:end,i) = HL_compo;
    endfor
    hess_lag = hess_lag + tmp *lm(j);
  endfor
  %hess_lag
  return
endfunction



%i-ème ligne du Jacobians
function [comp_i] = tmp(xy)
  global I
  [jacobian] = cal_Jacob(@cal_c, xy);
  comp_i = jacobian(I,1:end)';
endfunction



function verify_Hess_Lag(xy, lm)
  [hess_lag] = cal_Hess_Lag(xy, lm);
  [hl] = calcul_hl(xy,lm);
  erreur = abs(hess_lag - hl)/abs(hl);
  %hess_lag
  erreur
  erreur < 10e-8
endfunction


