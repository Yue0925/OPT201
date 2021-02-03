function [x, lm, info] = sqp(simul, x, lm, options)
  global NB_SIMUL
  %% initializer les valeurs sorties
  info.status = 0; %sortie normale
  info.niter = 0; %nb itérations
  
  n = length(x);
  m = length(lm);
  
  %% vérification vars entrées
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
  
  if m==0 % en cas de lambda est vide
    [e, c, g, a] = simul(4, x, lm);
    lm = -a'\g; %estimer lm moindres carrées
  endif
  
  %% les données initiales
  printf("\n Afficher les données initiales\n");
  x
  lm
  figure(1); clf(1) ;
  simul(1, x, lm); % tracer la chaîne
  if options.verb == 1
    printf('---------------------------------------------------------------------------------------------------------\n');
    printf("iter \t\t |gl| \t\t |ce| \t\t |x| \t\t |lm| \t\t alpha \t\t phi \n");
  endif
  
  x_h=[x]; % historique des x
  lm_h=[lm]; % historique des lm
  F_h = [ ]; % historique des norme infini de F
  
  %% algo méthode de Newton
  while 1
    %% sortir si max itération atteint
    if info.niter >= options.maxit
      info.status = 2;
      break
    endif
    
    %% calculer la nouvelle itération
    info.niter = info.niter + 1;
    NB_SIMUL = 0;
    [grad_lag, F_z, grad_F] = calcul_F(x, lm, simul);
    p = -grad_F\F_z; % (la direction de descent) résoudre système
    %p = -(grad_F' * grad_F)\(grad_F' * F_z); % mais en cas de matrice grad_F n'est pas inversible
    d = p(1:n); %d_k
    mu = p(n+1:end); %mu_k
 
    omega = 10e-4;
    i = 0;
    alpha = (1/2)^i; 
    phi_z = phi(x, lm, simul); % phi(z_k)
    pente = -2 * phi_z;  
    
    if options.rl % 1: pas unitaire
      x = x + d;
      lm = lm + mu;
    else % 0: recherche lineair d'Armijo
      if options.verb == 2
        printf("-------------------------------------------------------------\n");
        fprintf('iter %d, simul %d, phi %e, pente %e \n', info.niter, NB_SIMUL, phi_z, pente);
        fprintf("\t recherche linéair d'Armijo: \n");
        printf("\t\t alpha \t phip-phi \t DF(phi) \n");
        phip_phi = phi(x + alpha*d, lm + alpha*mu, simul)-phi_z; % affichage de i=0
        fprintf('\t\t %e \t %e \t %e \n', alpha, phip_phi, phip_phi/alpha); % affichage de i=0
      endif
      
      while phi(x + alpha*d, lm + alpha*mu, simul) > phi_z + omega*alpha*pente
        i = i+1;
        alpha = (1/2)^i;
        if options.verb == 2
          phip_phi = phi(x + alpha*d, lm + alpha*mu, simul)-phi_z;
          fprintf('\t\t %e \t %e \t %e \n', alpha, phip_phi, phip_phi/alpha);
        endif
      endwhile
      
      x = x + alpha*d;
      lm = lm + alpha*mu;
    endif
    
    %% après avoir obtenu le nouveau z_k:=(x_k, lambda_k)
    [e, c, g, a, hl, indic] = simul(5, x, lm); % calcul hessien de lagrangien
    [e, c, g, a] = simul(4, x, lm); % calcul c, g, a
    grad_lag = g + a'*lm;
    x_h=[x_h , x] ;
    lm_h=[lm_h , lm] ;
    F = max( [norm(grad_lag, Inf), norm(c, Inf)] );
    F_h=[F_h , F] ;
    
    if options.verb == 1
      format_sortie = '%d \t %e \t %e \t %e \t %e \t %e \t %e \n';
      fprintf(format_sortie, info.niter, norm(grad_lag, Inf), norm(c, Inf), norm(x, Inf), norm(lm, Inf), alpha, phi(x, lm, simul));
    endif
    if options.verb == 2 && !options.rl
      fprintf('\t\t |gl| = %e, \t |ce| = %e \n', norm(grad_lag, Inf), norm(c, Inf));
    endif

    %% sortir en cas de test d'arrêt vérifié
    if norm(grad_lag, Inf) <= options.tol(1) && norm(c, Inf) <= options.tol(2)
      info.status = 0; % terminaison normale
      printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n "); 
      printf("On trouve un point stationnaire z_*, voici son hessien lagrangien "); 
      full(hl)
      break
    endif
endwhile

  figure(2); clf(2) ;
  simul(1, x, lm); % tracer la chaîne

  %% calcul vitesse de convergence façon 2
  sig2 = [ ]; 
  for i = 1:info.niter-1 %k>=1
    sig2 = [sig2, abs(F_h(i+1))/(abs(F_h(i))^2)];
  endfor  
  
  printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n "); 
  printf("Affichage les résultats historiques:\n "); 
  x_h
  lm_h 
  F_h
  
  %% examiner le nature de solution
  a
  N = null(a) % noyeau de c'(x)
  hl_reduit = N'*hl*N % hessien réduite de lagrangien
  d = eig(hl_reduit) % vecteur des valeurs propres de hessien réduite
  eig(hl)

  if all(d > 0)
    printf("La matrice hessien réduite lagrangien est definie positive.\n");
    return
  endif
  if all(d >= 0)
    printf("La matrice hessien réduite lagrangien est semi-definie positive.\n");
    return
  endif
  if all(d < 0)
    printf("La matrice hessien réduite lagrangien est definie négative.\n");
    return
  endif
  if all(d <= 0)
    printf("La matrice hessien réduite lagrangien est semi-definie négative.\n");
    return
  endif
  printf("CONCLUSION: cette solution est un point selle.\n");
  
  return
endfunction


%% fonction calcule vecteur de F, matrice de F' et le gradient lagrangean
function [grad_lag, F_z, grad_F] = calcul_F(x, lm, simul)
  [e, c, g, a, hl, indic] = simul(5, x, lm); % calcul hessien de lagrangien
  [e, c, g, a] = simul(4, x, lm); % calcul c, g, a
  m = length(lm);
  grad_lag = g + a'*lm;
  F_z = [grad_lag' c']';
  grad_F = [hl a'; a sparse(m, m)];
  return
endfunction


%% la fonction de moindres-carrées
function phi_z = phi(x, lm, simul)
  [grad_lag, F_z, grad_F] = calcul_F(x, lm, simul);
  phi_z = 1/2 * norm(F_z, 2)^2;
  return
endfunction
