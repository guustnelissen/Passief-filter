clear; clc;

syms s

%% Gegeven admittanties
den = 0.1259*s^5 + 0.2895*s^3 + 0.1259*s;

Y_n1 = (9.6*s^6 + 32.66*s^4 + 32.66*s^2 + 9.6) / den;
Y_2  = (0.1126*s^4 + 0.1126*s^2) / den;

Y_12 = -0.1815 / den

%% Polen en nullen van Y_n1
[numY, denY] = numden(Y_n1);

numY_poly = sym2poly(numY);
denY_poly = sym2poly(denY);

zeros_Yn1 = roots(numY_poly);
poles_Yn1 = roots(denY_poly);

disp("Nullen van Y_n1:")
disp(zeros_Yn1)

disp("Polen van Y_n1:")
disp(poles_Yn1)
disp('ook pool in oneindig voor Y_n1')

%% transmissienullen

display('vijfvoudige nul in oneindig')


%% Eerste Cauer-stap: spoel naar massa
% Y_n1 = A*s + Y'
A = limit(Y_n1/s, s, inf);

Y_prime = simplify(Y_n1 - A*s);

disp("A = lim s->inf Y_n1/s =")
disp(vpa(A, 10))
disp(' = capaciteit in parallel')

disp("Y_prime = Y_n1 - A*s =")
disp(Y_prime)

[numYp, denYp] = numden(Y_prime);
numYp_poly = sym2poly(numYp);
zeros_Yp = roots(numYp_poly);
disp("Nullen van Yp:")
disp(zeros_Yp)
disp('polen zijn hetzelfde')

%% Tweede stap: inverteren naar impedantie
Z_prime = simplify(1/Y_prime);

% B = lim s->inf Z_prime / s
B = limit(Z_prime / s, s, inf);

Z_hat = simplify(Z_prime - B*s);

disp("B = lim s->inf s*(1/Y') =")
disp(vpa(B, 10))
disp('= spoel in serie')

disp("Z_hat = 1/Y' - B*s =")
disp(Z_hat)

[numZh, denZh] = numden(Z_hat);
numZh_poly = sym2poly(numZh);
zeros_Zh = roots(numZh_poly);
disp("Nullen van Zh:")
disp(zeros_Zh)
disp('polen van Zh zijn de nullen van Yp')

%% Derde stap: opnieuw inverteren
Y_hat = simplify(1/Z_hat);

C = limit(Y_hat/s, s, inf);

Y_star = simplify(Y_hat - C*s);

disp("C = lim s->inf Y_hat/s =")
disp(vpa(C, 10))

disp("Y_star = Y_hat - C*s =")
disp(Y_star)
disp('= capaciteit in parallel')

[numYs, denYs] = numden(Y_star);
numYs_poly = sym2poly(numYs);
zeros_Ys = roots(numYs_poly);
disp("Nullen van Ys:")
disp(zeros_Ys)
disp('polen van Ys zijn hetzelfde als de nullen van Yh')

%% vierde stap

Z_star = simplify(1/Y_star);

% D = lim s->inf Z_star / s
D = limit(Z_star / s, s, inf);

Z_tri = simplify(Z_star - D*s);

disp("D = lim s->inf Z_star/s =")
disp(vpa(D, 10))
disp('= spoel in serie')

disp("Z_tri = Z_star - D*s =")
disp(Z_tri)

[numZt, denZt] = numden(Z_tri);
numZt_poly = sym2poly(numZt);
zeros_Zt = roots(numZt_poly);
disp("Nullen van Z_tri:")
disp(zeros_Zt)
disp('polen van Z_tri zijn de nullen van Y_star')

%% vijfde stap

Y_tri = simplify(1/Z_tri);

E = limit(Y_tri/s, s, inf);

Y_sq = simplify(Y_tri - E*s);

disp("E = lim s->inf Y_tri/s =")
disp(vpa(E, 10))
disp('= capaciteit in parallel')


disp("Y_sq = Y_tri - E*s =")
disp(Y_sq)

[numYsq, denYsq] = numden(Y_sq);
numYsq_poly = sym2poly(numYsq);
zeros_Ysq = roots(numYsq_poly);
disp("Nullen van Ysq:")
disp(zeros_Ys)
disp('polen van Ysq zijn hetzelfde als de nullen van Yt')


%% Overzicht gevonden Cauer-elementen

waarden = vpa([A; B; C; D; E], 10);

namen = ["A"; "B"; "C"; "D"; "E"];
betekenis = [
    "capaciteit in parallel"
    "spoel in serie"
    "capaciteit in parallel"
    "spoel in serie"
    "capaciteit in parallel"
    ];

T = table(namen, waarden, betekenis, ...
    'VariableNames', {'Element', 'Waarde', 'Voorstelling'});

disp(T)