num = [0.09077 0 0 0];   % 0.09077*s^3
den = [1 0.9574 3.458 2.024 3.458 0.9574 1];
Z12 = 1;
z12 = 100;
Q = 30;

T0 = tf(num, den);

% Nullen, polen en gain ophalen
[z, p, k] = zpkdata(T0, 'v');

% Polen naar rechts verschuiven met factor 1/Q
% Interpretatie: reële deel wordt 1/Q dichter bij de imaginaire as gebracht
p_new = real(p) * (1 - 1/Q) + 1i*imag(p);

% Nieuwe transferfunctie maken met dezelfde nullen en gain
T0new_zpk = zpk(z, p_new, k);
T0new = tf(T0new_zpk);


T0 = T0new;

factors = linspace(1, 0.1, 5);
maxIter = length(factors);
iter = 0;
ratio_val = z12/Z12;
R2 = 1;
R1 = 50/240;

tol = 1e-6;


fig1 = figure;
hold on;
grid on;
xlabel('factor');
ylabel('Z12/z12');
title('K-waarde van Z12/z12 in functie van factor');

paar_matrix = [
    1 2 3
    1 2 4
    1 3 4
    1 4 5
    2 3 6
    2 5 6
    3 5 6
    4 5 6
    ];

for j = 1:size(paar_matrix, 1)

    iter = 0
    factor_vec = NaN(1, maxIter);
    K_vec      = NaN(1, maxIter);
    T = T0;

    while iter < maxIter
        iter = iter + 1;
    
        factor = factors(iter);
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
        
        %% nr en mr berekenen — alle geldige combinaties
    
        [num_coeffs, ~] = tfdata(Ro_kwad, 'v');
        all_zeros = roots(num_coeffs);
        n_zeros = length(all_zeros);
        
        % Sorteer nulpunten van meest negatief reëel deel naar meest positief
        [~, idx] = sort(real(all_zeros), 'ascend');
        all_zeros_sorted = all_zeros(idx);
        
        % Groepeer telkens 2 opeenvolgende nulpunten tot één nulpaar
        nulparen = reshape(all_zeros_sorted, 2, []).';
        
        % nulparen is nu een matrix:
        % rij 1 = nulpaar 1, kleinste reëel deel
        % rij 2 = nulpaar 2
        % ...
        % rij 6 = nulpaar 6, grootste reëel deel
    
        gekozen_paren = paar_matrix(j,:);
    
        chosen_zeros = nulparen(gekozen_paren, :).';
        chosen_zeros = chosen_zeros(:);
    
        F_s = real(poly(chosen_zeros));   % Dit is mr + nr
        
        
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
    
        z11_tf = minreal((R1*(n-nr))/(m+mr), 1e-8);
        z22_tf = minreal((R2*(n+nr))/(m+mr), 1e-8);
        z12_tf = minreal(sqrt(R1*R2)*(N12_pol/(m+mr)), 1e-8);
    
        syms s
    
        % Convert calculated tf objects to symbolic expressions
        z11_sym = tf_to_sym(z11_tf, s);
        z22_sym = tf_to_sym(z22_tf, s);
        z12_sym = tf_to_sym(z12_tf, s);
    
         %% Polen en nullen van Z11
        [numZ, denZ] = numden(z11_sym);
        
        numZ_poly = sym2poly(numZ);
        denZ_poly = sym2poly(denZ);
        
        zeros_z11 = roots(numZ_poly);
        poles_z11 = roots(denZ_poly);
        
        %% transmissienullen
        
        %% Cauersynthese met opschonen na elke afsplitsing
        
        cauerTol = 1e-10;
        
        % Zorg dat de startfunctie zelf ook proper is
        z11_sym = ratclean(z11_sym, s, cauerTol);
        
        %% Eerste stap: shunt-C aan ingang
        y11 = ratclean(1/z11_sym, s, cauerTol);
        
        A = scalarclean(limit(y11/s, s, inf), cauerTol);
        
        Y_star = ratclean(y11 - A*s, s, cauerTol);
        
        [numYs, denYs] = numden(Y_star);
        numYs_poly = sym2poly(numYs);
        denYs_poly = sym2poly(denYs);
        zeros_Ys = roots(numYs_poly);
        poles_Ys = roots(denYs_poly);
        
        %% Tweede stap: serie-L
        Z_star = ratclean(1/Y_star, s, cauerTol);
        
        B = scalarclean(limit(Z_star/s, s, inf), cauerTol);
        
        Z_hat = ratclean(Z_star - B*s, s, cauerTol);
        
        [numZh, denZh] = numden(Z_hat);
        numZh_poly = sym2poly(numZh);
        denZh_poly = sym2poly(denZh);
        zeros_Zh = roots(numZh_poly);
        poles_Zh = roots(denZh_poly);
        
        %% Derde stap: shunt-L
        Y_hat = ratclean(1/Z_hat, s, cauerTol);
        
        C = scalarclean(limit(Y_hat*s, s, 0), cauerTol);
        
        Y_tri = ratclean(Y_hat - C/s, s, cauerTol);
        
        [numYt, denYt] = numden(Y_tri);
        numYt_poly = sym2poly(numYt);
        denYt_poly = sym2poly(denYt);
        zeros_Yt = roots(numYt_poly);
        poles_Yt = roots(denYt_poly);
        
        %% Vierde stap: serie-C
        Z_tri = ratclean(1/Y_tri, s, cauerTol);
        
        D = scalarclean(limit(Z_tri*s, s, 0), cauerTol);
        
        Z_prime = ratclean(Z_tri - D/s, s, cauerTol);
        
        [numZp, denZp] = numden(Z_prime);
        numZp_poly = sym2poly(numZp);
        denZp_poly = sym2poly(denZp);
        zeros_Zp = roots(numZp_poly);
        poles_Zp = roots(denZp_poly);
        
        %% Vijfde stap: shunt-L van parallel-LC
        Y_prime = ratclean(1/Z_prime, s, cauerTol);
        
        E = scalarclean(limit(Y_prime*s, s, 0), cauerTol);
        
        Y_cir = ratclean(Y_prime - E/s, s, cauerTol);
        
        [numYc, denYc] = numden(Y_cir);
        numYc_poly = sym2poly(numYc);
        denYc_poly = sym2poly(denYc);
        zeros_Yc = roots(numYc_poly);
        poles_Yc = roots(denYc_poly);
        
        %% Zesde stap: shunt-C van parallel-LC
        F = ratclean(Y_cir / s, s, cauerTol);
        
        if ~isempty(symvar(F))
            error('F is nog afhankelijk van s. De afsplitsing is niet volledig gelukt of de topologie past niet.');
        end
    
        % %% Debug print A t/m F
        % fprintf('\n================ DEBUG A t/m F ================\n');
        % 
        % vars_names = {'A','B','C','D','E','F'};
        % vars_values = {A, B, C, D, E, F};
        % 
        % for k = 1:length(vars_names)
        %     fprintf('\n%s = \n', vars_names{k});
        %     disp(vpa(simplify(vars_values{k}), 8));
        % 
        %     fprintf('symvar(%s) = ', vars_names{k});
        %     disp(symvar(vars_values{k}));
        % 
        %     try
        %         fprintf('double(%s) = %.12g\n', vars_names{k}, sym_to_double(vars_values{k}));
        %     catch ME
        %         fprintf('double(%s) lukt niet: %s\n', vars_names{k}, ME.message);
        %     end
        % end
        % 
        % fprintf('\n================================================\n');
    
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
       redTol = 1e-12;
    
        A = remove_common_factors(ABCD_final(1,1), Tol=redTol);
        B = remove_common_factors(ABCD_final(1,2), Tol=redTol);
        C = remove_common_factors(zpk(ABCD_final(2,1)), Tol=redTol);
        D = remove_common_factors(zpk(ABCD_final(2,2)), Tol=redTol);
        
        Z11 = remove_common_factors(zpk(-D/C), Tol=redTol);
        Z22 = remove_common_factors(zpk(-A/C), Tol=redTol);
        Z12 = remove_common_factors(zpk(-1/C), Tol=redTol);
        
    
        Z11_tf = tf(Z11);
        Z22_tf = tf(Z22);
        Z12_tf = tf(Z12);
        
        %% iteratie
    
        zpk_z12 = remove_common_factors(zpk(z12));
        zpk_Z12 = remove_common_factors(zpk(Z12));
        ratio_Z12_z12 = zpk_Z12 / zpk_z12;
    
        ratio_val = dcgain(ratio_Z12_z12)
        
        % kleine imaginaire numerieke rommel verwijderen
        if abs(imag(ratio_val)) < 1e-10
            ratio_val = real(ratio_val);
        end
        
        
        factor_vec(iter) = factor;
        K_vec(iter)      = ratio_val;
        T = factor * T0;
    
    end

    factor_vec = factor_vec(1:iter);
    K_vec      = K_vec(1:iter);

    figure(fig1);
    plot(factor_vec, K_vec, 'o-', ...
         'DisplayName', sprintf('combi %d: [%d %d %d]', ...
         j, paar_matrix(j,1), paar_matrix(j,2), paar_matrix(j,3)));

    legend show;
end









function Zsym = tf_to_sym(Ztf, s)
    tfTol = 1e-10;

    Ztf = minreal(Ztf, 1e-8);

    [num, den] = tfdata(Ztf, 'v');

    scale = max([1, abs(num), abs(den)]);
    absTol = tfTol * scale;

    num = clean_poly(num, absTol);
    den = clean_poly(den, absTol);

    num = strip_leading_zeros(num, absTol);
    den = strip_leading_zeros(den, absTol);

    num_sym = poly2sym(sym(num), s);
    den_sym = poly2sym(sym(den), s);

    Zsym = ratclean(num_sym / den_sym, s, tfTol);
end

function x = sym_to_double(xsym)
xsym = simplify(xsym);

if ~isempty(symvar(xsym))
    error('Deze expressie bevat nog symbolische variabelen en kan niet naar double omgezet worden.');
end

x = double(vpa(xsym));
end



function expr_clean = ratclean(expr, s, relTol)
% RATCLEAN ruimt numerieke rommel op in een symbolische rationale functie.
% Gebruikt clean_poly, poly_gcd en strip_leading_zeros.

expr = simplify(expr);
[num, den] = numden(expr);

num = expand(num);
den = expand(den);

num_c = sym2poly(num);
den_c = sym2poly(den);

num_c = force_real_poly(num_c, relTol);
den_c = force_real_poly(den_c, relTol);

scale = max([1, abs(num_c), abs(den_c)]);
absTol = relTol * scale;

num_c = clean_poly(num_c, absTol);
den_c = clean_poly(den_c, absTol);

num_c = strip_leading_zeros(num_c, absTol);
den_c = strip_leading_zeros(den_c, absTol);

if isempty(num_c) || all(abs(num_c) < absTol)
    expr_clean = sym(0);
    return;
end

gcd_c = poly_gcd(num_c, den_c, absTol);

[num_c, rem_num] = deconv(num_c, gcd_c);
[den_c, rem_den] = deconv(den_c, gcd_c);

num_c = clean_poly(num_c, absTol);
den_c = clean_poly(den_c, absTol);
rem_num = clean_poly(rem_num, absTol);
rem_den = clean_poly(rem_den, absTol);

num_c = strip_leading_zeros(num_c, absTol);
den_c = strip_leading_zeros(den_c, absTol);

% Normaliseer niet naar monisch, want de gain moet behouden blijven.
expr_clean = simplify(poly2sym(sym(num_c), s) / poly2sym(sym(den_c), s));
end


function x_clean = scalarclean(x, tol)
% SCALARCLEAN maakt een symbolische scalar numeriek proper.

x = simplify(x);

if ~isempty(symvar(x))
    error('scalarclean: waarde bevat nog symbolische variabelen.');
end

xd = double(vpa(x, 16));

if abs(xd) < tol
    xd = 0;
end

x_clean = sym(xd);
end


function p = force_real_poly(p, tol)
% FORCE_REAL_POLY verwijdert verwaarloosbare imaginaire delen.

p = double(p);

scale = max([1, abs(real(p))]);

if max(abs(imag(p))) < tol * scale
    p = real(p);
end
end


