global A B L

A = 1;
B = -1;
L = [0.7, 0.5, 0.3, 0.2, 0.5]';

xy = [0.2   0.4   0.6   0.8   ...
      -1.0   -1.5   -1.5   -1.3]';
      
lm = [0.01, 0.01, 0.01, 0.01, 0.01]';
[e, c, g, a, hl, indic] = chs(4, xy, lm);
lm = -a'\g;
%lm % vérifié ok

for indic = 1:5
  [e, c, g, a, hl, indic] = chs(indic, xy, lm)
  full(a)
  printf("\n");
endfor
