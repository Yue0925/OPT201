function [e, c, g, a, hl, indic] = chs(indic,xy,lm)
  global A B L
  % initialisation des output de fonction
  e = [];
  c = [];
  g = [];
  a = [];
  hl = [];

  nn = length(xy)/2; % nombre noeuds(sauf extrémités)
  
  % vérification la bonne entrées
  if nn*2 ~=length(xy)
    indic = -1; % length(xy) impaire
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
      
    case 2 % calcul de e, c
      printf("in case %d\n", indic);
      [e, c] = cal_e_c(0,c,nb,x_complet,y_complet);
      
    case 4 % calcul g, a
      printf("in case %d\n", indic);
      [e, c] = cal_e_c(0,c,nb,x_complet,y_complet);
      g = zeros(2*nn, 1); % gradient de xy
      for i =1:nb-1
        g(i+nn) = (L(i) + L(i+1))/2;
      endfor
      a = sparse(nb,2*nn);
      for i=1:nb-1;
        a(i,i) = 2*x(i);
        a(i, i+nn) = 2*y(i);
      endfor
      
    otherwise
      printf("ERROR invalid indic number %d\n ", indic);
   endswitch

endfunction


function [e, c] = cal_e_c(e,c,nb,x_complet,y_complet)
  global A B L
  for i = 2:nb+1
    l_i = sqrt( (x_complet(i) - x_complet(i-1))^2 + (y_complet(i) - y_complet(i-1))^2 );
    e = e + l_i * (y_complet(i) + y_complet(i-1))/2; % énergie
    c(end+1) = l_i^2 - L(i-1)^2; % contrainte
  endfor
  return
endfunction
