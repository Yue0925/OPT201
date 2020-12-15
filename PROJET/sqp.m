function [x, lm, info] = sqp(simul, x, lm, options)
  % initializer les valeurs sorties
  info.status = 0; %sortie normale
  info.niter = 0; %nb itérations
  
  n = length(x);
  m = length(lm);
  
  % vérification vars entrées
  if n==0
    printf("ERROR: valeur initial x est vide. \n");
    info.status = 1; %inconsistance des arguments d’entrée
    return
  endif
  
  if options.tol(1)<0 || options.tol(2)<0 || options.tol(1) >1 || options.tol(2)>1
    printf("ERROR: les seuils de tolérance ne sont pas compris entre 0 et 1. \n");
    info.status = 1;
    return
  endif
  
  if m==0
    [e, c, g, a] = simul(4, x, lm);
    lm = -a'\g; %estimer lm moindre carrées
  endif
  
  % le pas initial
  printf("\n k=%d \n", info.niter);
  x
  lm
  figure(1); clf(1) ;
  simul(1, x, lm);
  %pause;
  
  % algo méthode de Newton
  while 1
    
    % sortir si max itération atteint
    if info.niter > options.maxit
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
    
    % afficher les valeurs intermédiaires
    printf("\n k=%d \n", info.niter);
    x
    lm
    [e] = simul(2, x, lm) % calcul e
    printf("norm sup de gradient lagrangien "); 
    norm(grad_lag, Inf)
    printf("norm sup de contraintes égalités "); 
    norm(c, Inf)
    figure(info.niter); clf(info.niter) ;
    simul(1, x, lm);
    %pause; %taper un tab pour avancer

    % sortir en cas de test d'arrêt vérifié
    if norm(grad_lag, Inf) <= options.tol(1) && norm(c, Inf) <= options.tol(2)
      info.status = 0; % terminaison normale
      break
    endif

  endwhile
  
  return
endfunction
