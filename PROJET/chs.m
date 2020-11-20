function [e, c, g, a, hl, indic] = chs(indic,xy,lm)
  global A B L
  
  EXIT_SUCESS = 0; % sortie normale 
  EXIT_FAILURE = 1; % paramètre(s) d’entrée non correct(s)

  % initialisation des output de fonction
  e = [];
  c = [];
  g = [];
  a = [];
  hl = [];

  nn = length(xy)/2; % nombre noeuds(sauf extrémités)
  
  % vérification la bonne entrées
  if nn*2 ~=length(xy)
    indic = EXIT_FAILURE; % length(xy) impaire
    return
  endif
  
  nb = nn+1; % nombre barres
  x = xy(1:nn); 
  y = xy(nn+1 : end);
  
  x_complet = [0;x;A]; % les coordonées en ajoutant deux extrémités
  y_complet = [0;y;B];
  
  % pilote indic
  switch (indic)
    
    case 1 % tracer la chaîne
      printf("in case %d\n", indic);
      plot(x_complet, y_complet, '-b')
      indic = EXIT_SUCESS;
      return
      
    case 2 % calcul de e, c
      printf("in case %d\n", indic);
      [e, c] = cal_e_c(c,nb,x_complet,y_complet);
      indic = EXIT_SUCESS;
      return
      
    case 4 % calcul de e, c, g, a
      printf("in case %d\n", indic);
      [e, c] = cal_e_c(c,nb,x_complet,y_complet);
      [g, a] = cal_g_a(nn,nb,xy);
      indic = EXIT_SUCESS;
      return
      
    otherwise
      printf("ERROR invalid indic number %d\n ", indic);
      indic = EXIT_FAILURE;
      return
   endswitch

endfunction


function [e, c] = cal_e_c(c,nb,x_complet,y_complet)
  global A B L
  e = 0;
  for i = 2:nb+1
    l_i = sqrt( (x_complet(i) - x_complet(i-1))^2 + (y_complet(i) - y_complet(i-1))^2 );
    e = e + l_i * (y_complet(i) + y_complet(i-1))/2; % énergie
    c(end+1) = l_i^2 - L(i-1)^2; % contrainte
  endfor
  return
endfunction


function [g, a] = cal_g_a(nn,nb,xy)
  global A B L
  g = zeros(2*nn, 1); % gradient de xy
  for i =1:nb-1
    g(i+nn) = (L(i) + L(i+1))/2;
  endfor
  
  a = sparse(nb,2*nn);
  
  diff2 = [2*xy(1)]; %2(x_1-x_0)
  for i = 2:nn
    diff2(end+1) = 2*(xy(i) - xy(i-1));
  endfor
  diff2(end+1) = 2*(xy(nn)-A); % 2(x_n-A)
  
  diff2(end+1) = 2*xy(nn+1); % 2(y_1 - y_0)
  for i = nn+2:2*nn
    diff2(end+1) = 2*(xy(i) - xy(i-1));
  endfor
  diff2(end+1) = 2*(xy(2*nn)-B);

  for i=1:nn;
    a(i,i) = diff2(i);
    a(i,nn+i) = diff2(nb+i);
  endfor
  a(nb,nn) = diff2(nb);
  a(nb,2*nn) = diff2(2*nb);
  return
endfunction

