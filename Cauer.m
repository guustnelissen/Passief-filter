syms s

den = 0.1126*s^4 + 0.1126*s^2;

den = 2*s^6 + 6.804*s^4 + 6.804*s^2 + 2;

z11 = (0.02622*s^5 + 0.06032*s^3 + 0.02622*s) / den;

z22 = (1.789*s^5 + 3.759*s^3 + 1.789*s) / den;

z12 = (0.03782*s^3) / den;


%% Polen en nullen van Z11
[numZ, denZ] = numden(z11);

numZ_poly = sym2poly(numZ);
denZ_poly = sym2poly(denZ);

zeros_z11 = roots(numZ_poly);
poles_z11 = roots(denZ_poly);

disp("Nullen van z11:")
disp(zeros_z11)

disp("Polen van z11:")
disp(poles_z11)
disp('ook nul in oneindig voor z11')

%% transmissienullen

display('tweevoudige nul in oneindig')
display('drievoudige nul in nul')


%% Eerste stap
y11 = simplify(1/z11);

A = limit(y11/s, s, inf);

Y_star = simplify(y11 - A*s);

disp("A = lim s->inf y11/s =")
disp(vpa(A, 10))

disp("Y_star = y11 - A*s =")
disp(Y_star)
disp('= capaciteit in parallel')

[numYs, denYs] = numden(Y_star);
numYs_poly = sym2poly(numYs);
denYs_poly = sym2poly(denYs);
zeros_Ys = roots(numYs_poly);
poles_Ys = roots(denYs_poly);
disp("Nullen van Ys:")
disp(zeros_Ys)
disp("Polen van Ys:")
disp(poles_Ys)

%% Tweede stap
Z_star = simplify(1/Y_star);

B = limit(Z_star/s, s, inf);

Z_hat = simplify(Z_star - B*s);

disp("B = lim s->inf Z_star/s =")
disp(vpa(B, 10))

disp("Z_hat = Z_star - B*s =")
disp(Z_hat)
disp('= spoel in serie')

[numZh, denZh] = numden(Z_hat);
numZh_poly = sym2poly(numZh);
denZh_poly = sym2poly(denZh);
zeros_Zh = roots(numZh_poly);
poles_Zh = roots(denZh_poly);
disp("Nullen van Zh:")
disp(zeros_Zh)
disp("Polen van Zh:")
disp(poles_Zh)

%% Derde stap
Y_hat = simplify(1/Z_hat);

C = limit(Y_hat*s, s, 0);

Y_tri = simplify(Y_hat - C/s);

disp("C = lim s->0 Y_hat*s =")
disp(vpa(C, 10))
disp('1/C = ')
disp(vpa(1/C, 10))

disp("Y_tri = Y_hat - C/s =")
disp(Y_tri)
disp('= spoel in parallel')

[numYt, denYt] = numden(Y_tri);
numYt_poly = sym2poly(numYt);
denYt_poly = sym2poly(denYt);
zeros_Yt = roots(numYt_poly);
poles_Yt = roots(denYt_poly);
disp("Nullen van Yt:")
disp(zeros_Yt)
disp("Polen van Yt:")
disp(poles_Yt)

%% Vierde stap
Z_tri = simplify(1/Y_tri);

D = limit(Z_tri*s, s, 0);

Z_prime = simplify(Z_tri - D/s);

disp("D = lim s->0 Z_tri*s =")
disp(vpa(D, 10))
disp('1/D = ')
disp(vpa(1/D, 10))

disp("Z_prime = Z_tri -D/s =")
disp(Z_prime)
disp('= capaciteit in serie')

[numZp, denZp] = numden(Z_prime);
numZp_poly = sym2poly(numZp);
denZp_poly = sym2poly(denZp);
zeros_Zp = roots(numZp_poly);
poles_Zp = roots(denZp_poly);
disp("Nullen van Zp:")
disp(zeros_Zp)
disp("Polen van Zp:")
disp(poles_Zp)

%% Vijfde stap: Bereken parallel LC uit Z_prime

[numZp, denZp] = numden(Z_prime);
numZp = expand(numZp);
denZp = expand(denZp);

% Voor een parallel LC geldt:
% Z = (L*s)/(L*C*s^2 + 1)
%
% Algemene vorm:
% Z_prime = (a*s)/(b*s^2 + c)
%
% Dan:
% L = a/c
% C = b/a

a = simplify(subs(diff(numZp, s), s, 0));        % coëfficiënt van s in teller
b = simplify(subs(diff(denZp, s, 2)/2, s, 0));   % coëfficiënt van s^2 in noemer
c = simplify(subs(denZp, s, 0));                 % constante term in noemer

L_lc = simplify(a/c);
C_lc = simplify(b/a);

disp("L van parallel LC =")
disp(vpa(L_lc, 10))

disp("C van parallel LC =")
disp(vpa(C_lc, 10))

% Controle
Z_check = simplify(1 / (s*C_lc + 1/(s*L_lc)));

disp("Controle Z_prime - Z_check =")
disp(simplify(Z_prime - Z_check))

%% Overzicht gevonden Cauer-elementen

waarden = vpa([A; B; 1/C; 1/D; L_lc; C_lc], 10);

namen = ["A"; "B"; "1/C"; "1/D"; "L_lc"; "C_lc"];

betekenis = [
    "capaciteit in parallel"
    "spoel in serie"
    "spoel in parallel"
    "capaciteit in serie"
    "spoel van parallel LC"
    "capaciteit van parallel LC"
    ];

T = table(namen, waarden, betekenis, ...
    'VariableNames', {'Element', 'Waarde', 'Voorstelling'});

disp(T)
disp('L_lc en C_lc in parallel met elkaar, en dan in serie met de rest')
