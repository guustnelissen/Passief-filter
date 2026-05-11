num = [0.09077 0 0 0];   % 0.09077*s^3
den = [1 0.9574 3.458 2.024 3.458 0.9574 1];
Z12 = 1;
z12 = 100;

T = tf(num, den);
maxIter = 100;
iter = 0;
ratio_val = z12/Z12

tol = 1e-6;

while abs(ratio_val - 1) > tol  && iter < maxIter
    %% Berekenen van T(-s)
    % Haal de teller en noemer van de geschaalde transferfunctie T op
    [numT, denT] = tfdata(T, 'v');
    
    % Bepaal de lengte voor de tekenvector
    len_numT = length(numT);
    len_denT = length(denT);
    
    % Creëer de vectoren die de tekens van de oneven machten omdraaien
    tekens_numT = (-1).^( (len_numT-1) : -1 : 0 );
    tekens_denT = (-1).^( (len_denT-1) : -1 : 0 );
    
    % Pas de transformatie s -> -s toe
    numT_min_s = numT .* tekens_numT;
    denT_min_s = denT .* tekens_denT;
    
    % Maak de transferfunctie T(-s) aan
    T_min_s = tf(numT_min_s, denT_min_s);
    
    
    %% (6.25) (6.26)
    %T(s)*T(-s)
    T_kwad_bd = T * T_min_s;
    
    Ro_kwad = 1 - (4*(R1/R2) * T_kwad_bd); %het lijkt alsof teller en noemer hetzelfde is maar dit is niet, matlab rondt te veel af
    %% n en m berekenen
    %n en m berkenen kan ook al sneller uit T(s) met formule (6.25 (a))
    % 1. Haal de stabiele noemer (m + n)
    [~, den_coeffs] = tfdata(Ro_kwad, 'v');
    all_poles = roots(den_coeffs);
    stable_poles = all_poles(real(all_poles) < -1e-5); %pak de neg polen in het LHV
    D_s = poly(stable_poles); % Dit is m + n
    
    % 2. Splits m en n
    m_coeffs = zeros(size(D_s)); %initiatie
    n_coeffs = zeros(size(D_s));
    
    % Even indices (s^0, s^2...) en Oneven indices (s^1, s^3...)
    % Let op: MATLAB indexeert van hoog naar laag [s^3, s^2, s^1, s^0]
    indices = length(D_s)-1:-1:0;
    m_mask = mod(indices, 2) == 0; %is het getal deelbaar door 2?
    n_mask = mod(indices, 2) ~= 0; %is er een restwaarde?
    
    m_coeffs(m_mask) = D_s(m_mask);
    n_coeffs(n_mask) = D_s(n_mask);
    
    m = tf(m_coeffs, 1);
    n = tf(n_coeffs, 1);
    n2_m2 = m^2 - n^2; %dit is inderdaad de noemer van ro_kwad, :)
    
    %% nr en mr bereken
    
    [num_coeffs, ~] = tfdata(Ro_kwad, 'v');
    all_zeros = roots(num_coeffs);
    stable_zeros = all_zeros(real(all_zeros) < 0); %We kiezen hier even ez de linkse nullen. 
    % Dit kunnen ook andere zijn, dit gaan mss ook zo moeten want je gaat mss niet de goede K factor vinden met de nullen dat je hebt gekozen
    F_s = poly(stable_zeros); % Dit is mr + nr
    
    % 2. Splits m en n
    mr_coeffs = zeros(size(F_s)); %initiatie
    nr_coeffs = zeros(size(F_s));
    
    % Even indices (s^0, s^2...) en Oneven indices (s^1, s^3...)
    % Let op: MATLAB indexeert van hoog naar laag [s^3, s^2, s^1, s^0]
    indices = length(F_s)-1:-1:0;
    mr_mask = mod(indices, 2) == 0; %is het getal deelbaar door 2?
    nr_mask = mod(indices, 2) ~= 0; %is er een restwaarde?
    
    mr_coeffs(mr_mask) = F_s(mr_mask);
    nr_coeffs(nr_mask) = F_s(nr_mask);
    
    mr = tf(mr_coeffs, 1);
    nr = tf(nr_coeffs, 1);
    nr2_mr2 = mr^2 - nr^2; %dit is inderdaad de teller van ro_kwad, :)
    

    %% N12 door (6.25 (a)) om te vormen
    fef =  2* sqrt(R1/R2) * T;
    [num_coeffs, ~] = tfdata(fef, 'v');
    [~, den_coeffs] = tfdata(fef, 'v');
    N12 = num_coeffs;
    N12_pol = tf(N12,1);
    %als check:
    %den_coeffs = n + m ; KLOPT :)
    
    %% Z-parameters berekenen met Table 6.2
    z11 = (R1*(n-nr))/(m+mr);
    z22 = (R2*(n+nr))/(m+mr);
    z12 = sqrt(R1*R2)*(N12_pol/(m+mr));
    zpk(z11);
    zpk(z22);
    zpk(z12);

    z11_tf = minreal((R1*(n-nr))/(m+mr));
    z22_tf = minreal((R2*(n+nr))/(m+mr));
    z12_tf = minreal(sqrt(R1*R2)*(N12_pol/(m+mr)));

    syms s

    % Convert calculated tf objects to symbolic expressions
    z11_sym = tf_to_sym(z11, s);
    z22_sym = tf_to_sym(z22, s);
    z12_sym = tf_to_sym(z12, s);

     %% Polen en nullen van Z11
    [numZ, denZ] = numden(z11_sym);
    
    numZ_poly = sym2poly(numZ);
    denZ_poly = sym2poly(denZ);
    
    zeros_z11 = roots(numZ_poly);
    poles_z11 = roots(denZ_poly);
    
    %% transmissienullen
    
    %% Eerste stap
    y11 = simplify(1/z11_sym);
    
    A = limit(y11/s, s, inf);
    
    Y_star = simplify(y11 - A*s);
    
    [numYs, denYs] = numden(Y_star);
    numYs_poly = sym2poly(numYs);
    denYs_poly = sym2poly(denYs);
    zeros_Ys = roots(numYs_poly);
    poles_Ys = roots(denYs_poly);
    
    %% Tweede stap
    Z_star = simplify(1/Y_star);
    
    B = limit(Z_star/s, s, inf);
    
    Z_hat = simplify(Z_star - B*s);
    
    
    [numZh, denZh] = numden(Z_hat);
    numZh_poly = sym2poly(numZh);
    denZh_poly = sym2poly(denZh);
    zeros_Zh = roots(numZh_poly);
    poles_Zh = roots(denZh_poly);
    
    %% Derde stap
    Y_hat = simplify(1/Z_hat);
    
    C = limit(Y_hat*s, s, 0);
    
    Y_tri = simplify(Y_hat - C/s);
    
    
    [numYt, denYt] = numden(Y_tri);
    numYt_poly = sym2poly(numYt);
    denYt_poly = sym2poly(denYt);
    zeros_Yt = roots(numYt_poly);
    poles_Yt = roots(denYt_poly);
    
    %% Vierde stap
    Z_tri = simplify(1/Y_tri);
    
    D = limit(Z_tri*s, s, 0);
    
    Z_prime = simplify(Z_tri - D/s);
    
    
    [numZp, denZp] = numden(Z_prime);
    numZp_poly = sym2poly(numZp);
    denZp_poly = sym2poly(denZp);
    zeros_Zp = roots(numZp_poly);
    poles_Zp = roots(denZp_poly);
    
    %% Vijfde stap
    Y_prime = simplify(1/Z_prime);
    
    E = limit(Y_prime*s, s, 0);
    
    Y_cir = simplify(Y_prime - E/s);
    
    [numYc, denYc] = numden(Y_cir);
    numYc_poly = sym2poly(numYc);
    denYc_poly = sym2poly(denYc);
    zeros_Yc = roots(numYc_poly);
    poles_Yc = roots(denYc_poly);
    
    %% Zesde stap
    %hiermee wordt de laatste transmitienullen gerealiseerd
    %overgebleven component realisere
    
    F = Y_cir / s;

    %% ABCD matrices berekenen om zo te checken of Z-param juist zijn

    % Terug van symbolisch naar numerieke componentwaarden

    C1 = sym_to_double(A);      % capaciteit in parallel
    L1 = sym_to_double(B);      % spoel in serie
    L2 = 1/sym_to_double(C);    % spoel in parallel
    C2 = 1/sym_to_double(D);    % capaciteit in serie
    L3 = 1/sym_to_double(E);    % spoel van parallel LC
    C3 = sym_to_double(F);      % capaciteit van parallel LC

    s = tf('s');
    M_C1shunt = [1, 0; -C1*s, 1];
    M_L1serie = [1, -L1*s; 0, 1];
    M_L2shunt = [1, 0; -1/(L2*s), 1];
    M_C2serie = [1, -1/(C2*s); 0, 1];
    M_L3shunt = [1, 0; -1/(L3*s), 1];
    M_C3shunt = [1, 0; -C3*s, 1];
    
    ABCD = M_C3shunt * M_L3shunt * M_C2serie * M_L2shunt * M_L1serie * M_C1shunt;
    ABCD_final = ABCD;
    
    % Z-parameters extraheren uit de ABCD matrix
    % z11 = A/C, z12 = 1/C, z22 = D/C
    A = remove_common_factors(ABCD_final(1,1));
    B = remove_common_factors(ABCD_final(1,2));
    C = remove_common_factors((zpk(ABCD_final(2,1))));
    D = remove_common_factors(ABCD_final(2,2));
    
    Z11 = remove_common_factors(zpk(-D/C));
    Z22 = remove_common_factors(zpk(-A/C));
    Z12 = remove_common_factors(zpk(-1/C));

    Z11_tf = tf(Z11);
    Z22_tf = tf(Z22);
    Z12_tf = tf(Z12);
    
    %%
    
    ratio_sys = remove_common_factors(zpk(z12/Z12));
    ratio_val = dcgain(ratio_sys);

    
    disp('Ratio z12/Z12 is:');
    disp(ratio_val);
    
    if ratio_val > 1
        T = 1.1*T;
    elseif ratio_val < 1
        T = 0.9*T;
    end
end


function Zsym = tf_to_sym(Ztf, s)
[num, den] = tfdata(Ztf, 'v');

num_sym = poly2sym(num, s);
den_sym = poly2sym(den, s);

Zsym = simplify(num_sym / den_sym);
end

function x = sym_to_double(xsym)
xsym = simplify(xsym);

if ~isempty(symvar(xsym))
    error('Deze expressie bevat nog symbolische variabelen en kan niet naar double omgezet worden.');
end

x = double(vpa(xsym));
end