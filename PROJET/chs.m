1;

function [e, c, g, a, hl, indic] = chs(indic,xy,lm)  
  
  global A B L
  
  % initialisation des output de fonction
  e = [];
  c = [];
  g = [];
  a = [];
  hl = [];

  nn = length(xy)/2; % nombre noeuds (sauf exétrimités)
  
  % vérification la bonne entrées
  if nn*2 ~=length(xy)
    indic = -1; % length(xy) impaire
    return
  endif
  
  nb = nn+1; % nombre barres
  x = xy(1:nn); 
  y = xy(nn+1 : end);
  
  % tracer la chaîne
  if indic == 1
    plot([0;x;A], [0;y;B], '-b');
    return
  endif
  
  % calcul de e, c

endfunction