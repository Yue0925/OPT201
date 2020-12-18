function [x, lm, info] = sqp(simul, x, lm, options)
  % initializer les valeurs sorties
  info.status = 0; %sortie normale
  info.niter = 0; %nb itérations
  
  n = length(x);
  m = length(lm);
  
  % vérification vars entrées
  if n==0
    error("ERROR: valeur initial x est vide. \n");
    info.status = 1; %inconsistance des arguments d’entrée
    return
  endif

  if options.maxit<1
    error("ERROR: le maximum d'itération est inférieur à 1 ici. \n");
    info.status = 1; %inconsistance des arguments d’entrée
    return
  endif
  
  if options.tol(1)<0 || options.tol(2)<0 || options.tol(1) >1 || options.tol(2)>1
    error("ERROR: les seuils de tolérance ne sont pas compris entre 0 et 1. \n");
    info.status = 1;
    return
  endif
  
  if m==0
    [e, c, g, a] = simul(4, x, lm);
    lm = -a'\g; %estimer lm moindres carrées
  endif
  
  % le pas initial
  printf("\n k=%d \n", info.niter);
  x
  lm
  [e] = simul(2, x, lm); % calcul e
  e
  figure(1); clf(1) ;
  simul(1, x, lm);
  
  x_h=[x]; % historique des x
  lm_h=[lm]; % historique des lm
  F_h = [ ];
  
  % algo méthode de Newton
  while 1
    
    % sortir si max itération atteint
    if info.niter >= options.maxit
      info.status = 2;
      break
    endif
    
    % calculer la nouvelle itération
    info.niter = info.niter + 1;
    [e, c, g, a, hl, indic] = simul(5, x, lm); % calcul hessien de lagrangien
    [e, c, g, a] = simul(4, x, lm); % calcul c, g, a
    A = [hl a'; a sparse(m, m)];
    b = [g' c']';
    res = -A\b;
    d = res(1:n);
    x = x + d;
    lm = res(n+1:end);
    grad_lag = g + a'*lm;
    [e] = simul(2, x, lm); % calcul e
    
    x_h=[x_h , x] ;
    lm_h=[lm_h , lm] ;
    
    % afficher les valeurs intermédiaires
    printf("\n k=%d \n", info.niter);
    x
    lm
    e
    printf("norm sup de gradient lagrangien "); 
    norm(grad_lag, Inf)
    printf("norm sup de contraintes égalités "); 
    norm(c, Inf)
    figure(info.niter); clf(info.niter) ;
    simul(1, x, lm);
    
    F = max( [norm(grad_lag, Inf), norm(c, Inf)] );
    F_h=[F_h , F] ;

    % sortir en cas de test d'arrêt vérifié
    if norm(grad_lag, Inf) <= options.tol(1) && norm(c, Inf) <= options.tol(2)
      info.status = 0; % terminaison normale
      printf("On trouve un point critique, voici son hessien lagrangien "); 
      full(hl)
      break
    endif

endwhile

  % calcul convergence façon 1
  sig = [ ]; 
  for i = 2:info.niter-1 %k>=1
    sig = [sig, norm(x_h(1:end, i+1) - x, 1)/(norm(x_h(1:end, i) - x, 1)^2)];
  endfor

  % calcul convergence façon 2
  sig2 = [ ]; 
  for i = 1:info.niter-1 %k>=1
    sig2 = [sig2, abs(F_h(i+1))/(abs(F_h(i))^2)];
  endfor  
  
  printf("Affichage les résultats historiques:\n "); 
  printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n "); 
  x_h
  lm_h
  F_h
  sig
  sig2

  % matrice définie-positive
  if issymmetric(hl)
    d = eig(hl);
    if all(d > 0)
      printf("La matrice hessien lagrangien est definie positive.\n");
      return
    endif
    if all(d >= 0)
      printf("La matrice hessien lagrangien est semi-definie positive.\n");
      return
    endif
    if all(d < 0)
      printf("La matrice hessien lagrangien est definie négative.\n");
      return
    endif
    if all(d <= 0)
      printf("La matrice hessien lagrangien est semi-definie négative.\n");
      return
    endif
    printf("Cette solution est un point selle.\n");
  endif
  
  return
endfunction
