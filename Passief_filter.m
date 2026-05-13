fs1 = 460 * 10^3;
fd1 = 1260 * 10^3;
fd2 = 1740 * 10^3;
fs2 = 8500 * 10^3;
RL = 240;
RE = 50;
B = (fd2-fd1)/sqrt(fd2*fd1);
Ad = -0.40; %in dB
As = -45.0; %in dB

prod_d = fd1 * fd2;
prod_s = fs1 * fs2;

fs2_new = zoek_fs2(fs1, prod_d, fs2);
fprintf('B   = %.2f\n', B);
fprintf('fs2_new    = %.2f\n', fs2_new);
fprintf('prod_s_new = %.2f\n', fs1 * fs2_new);
fprintf('prod_d     = %.2f\n', prod_d);

function fs2_new = zoek_fs2(fs1, prod_d, fs2_huidige)
    % Basisgeval: is het al gelijk genoeg?
    if abs(fs1 * fs2_huidige - prod_d) < 1e-6
        fs2_new = fs2_huidige;
        disp('Gelijk');
    else
        % Niet gelijk: bereken betere waarde en probeer opnieuw
        fs2_new = zoek_fs2(fs1, prod_d, prod_d / fs1);
    end
end

f_center = sqrt(prod_d);
display(f_center)
%% Frequentie normalisatie

fs1_n = fs1/f_center
fs2_n = fs2_new/f_center
fd1_n = fd1/f_center
fd2_n = fd2/f_center
f_center_n = f_center/f_center;

%Weerstandsnormalisatie
RL_norm = RL/RL % 1
RE_norm = RE/RL
%BD naar LD
fs = abs((fs1_n - (1/fs1_n))/B)
fss = abs((fs2_n - (1/fs2_n))/B); %ff testen of deze gelijk zijn (is dus wel)
%de abs is nodig omdat fs_n kleiner is dan 1 en dan krijgen we een neg
%resultaat maar de abs-waarde is gelijk

fd = abs((fd1_n - (1/fd1_n))/B)
fdd = abs((fd2_n - (1/fd2_n))/B);
%orde berekenen
k = fd/fs
k = (fd2_n - fd1_n)/(fs2_n - fs1_n);
rimpel = sqrt(10^(-Ad/10)-1)
n = (log(sqrt(10^(-As/10)-1)/rimpel)) / log(1/k)
orde = ceil(n)
%% Polen en transferfunctie berekenen voor cauer synthese
close all
polen = zeros(1,orde);
m = 0:(orde-1);
    polen(m+1) = ((1i*fd)/rimpel^(1/orde)) * exp(pi*1i*(2*m+1)/(2*orde));

    
disp('De berekende polen zijn:');
disp(polen.');   %met .' ga je van een rij-vector naar een kolom vector

figure;
zplane([], polen.'); % Toont de polen op de cirkel
title('Pool-Nulpunten diagram (Butterworth Orde 3)');
grid on

TF_teller = real(prod(-polen)); %dit kan enkel bij transferfuncties van laagdoorlaat filters

TF_noemer = real(poly(polen));

H = tf(TF_teller, TF_noemer); %Dit is voor de laagdoorlaat
display(H)

%% Transformatie van Laagdoorlaat (H) naar Banddoorlaat (H_bp) 

% Gebruik MATLAB's lp2bp functie. Deze past exact de substitutie 
% p = (q^2 + 1)/(B*q) toe op de polynomen van de transferfunctie.
[TF_teller_bp, TF_noemer_bp] = lp2bp(TF_teller, TF_noemer, f_center_n, B);

% Maak het nieuwe Transfer Function object aan
H_bp = tf(TF_teller_bp, TF_noemer_bp);
disp('De genormaliseerde banddoorlaat transferfunctie H_bp(s) is:');
display(H_bp)

% Plot het resultaat om te zien hoe de polen gesplitst zijn en
% de nulpunten in de oorsprong zijn toegevoegd
figure
pzmap(H_bp);
title('Pool-Nulpunten diagram (Banddoorlaat)');
grid on;

% Bode plot om de frequentierespons te controleren
figure
bode(H_bp);
title('Bode plot van het Banddoorlaat Filter (Genormaliseerd)');
grid on;

%% Schaling H met spanningsdeler (6.22)
R2 = RL_norm; %1
R1 = RE_norm; %Rn = 0.23

T = (R2/(R1 + R2)) * H_bp;
display(T)
figure
bode(T);

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

disp('De transferfunctie T(-s) is:');
display(T_min_s)


%%
T_kwad_bd = T * T_min_s;
display(T_kwad_bd)

Ro_kwad = 1 - (4*(R1/R2) * T_kwad_bd) %het lijkt alsof teller en noemer hetzelfde is maar dit is niet, matlab rondt te veel af

%[num_ro2, denum_ro2] = tfdata(Ro_kwad, 'v');
%poles = roots(denum_ro2);
figure
pzplot(Ro_kwad);

%% n en m berekenen
%n en m berkenen kan ook al sneller uit T(s) met formule (6.25 (a))
% 1. Haal de stabiele noemer (m + n)
[~, den_coeffs] = tfdata(Ro_kwad, 'v');
all_poles = roots(den_coeffs);
stable_poles = all_poles(real(all_poles) < -1e-5); %pak de neg polen in het LHV want stabiel
D_s = poly(stable_poles) % Dit is m + n
figure;
plot(real(stable_poles), imag(stable_poles), 'x')
grid on
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

m = tf(m_coeffs, 1)
n = tf(n_coeffs, 1)
n2_m2 = m^2 - n^2 %dit is inderdaad de noemer van ro_kwad, :)
%%
function ok = is_geldig(zeros_vec, tol)
    if nargin < 2, tol = 1e-6; end
    remaining = zeros_vec(:);
    ok = true;

    while ~isempty(remaining)
        z = remaining(1);
        remaining(1) = [];

        if abs(imag(z)) < tol
            continue  % reële nul, geen partner nodig
        end

        % Check 1: is de geconjugeerde aanwezig?
        diffs_conj = abs(remaining - conj(z));
        [minVal, idx] = min(diffs_conj);
        if minVal > tol
            ok = false;
            return  % geen geconjugeerde gevonden
        end
        remaining(idx) = [];  % geconjugeerde gevonden, verwijder

        % Check 2: zit het gespiegeld paar (+a+jb en -a+jb) er NIET in?
        for j = 1:length(remaining)
            zj = remaining(j);
            if abs(real(z) + real(zj)) < tol && abs(imag(z) - imag(zj)) < tol
                ok = false;
                return  % gespiegeld paar gevonden → ongeldig
            end
        end
    end
end
%% nr en mr berekenen — alle geldige combinaties

[num_coeffs, ~] = tfdata(Ro_kwad, 'v');
all_zeros = roots(num_coeffs);
n_zeros = length(all_zeros);
k = 6;

combo_indices = nchoosek(1:n_zeros, k);
n_combos = size(combo_indices, 1);
fprintf('Totaal aantal combinaties: %d\n', n_combos);

geldige_combos = {};

for i = 1:n_combos
    kandidaat = all_zeros(combo_indices(i,:));

    if is_geldig(kandidaat)  % <-- bevat nu BEIDE checks
        geldige_combos{end+1} = kandidaat;
        fprintf('\nGeldige combinatie %d: index [%s]\n', ...
            length(geldige_combos), num2str(combo_indices(i,:)));
        disp(kandidaat);

        % F_s berekenen voor deze combinatie
        F_s = real(poly(kandidaat));  % real() om numerieke ruis weg te halen

        % Splits mr en nr
        mr_coeffs = zeros(size(F_s));
        nr_coeffs = zeros(size(F_s));
        indices = length(F_s)-1:-1:0;
        mr_coeffs(mod(indices,2)==0) = F_s(mod(indices,2)==0);
        nr_coeffs(mod(indices,2)~=0) = F_s(mod(indices,2)~=0);
        % Strip leading zeros
        mr_coeffs_clean = mr_coeffs(find(mr_coeffs, 1) : end);
        nr_coeffs_clean = nr_coeffs(find(nr_coeffs, 1) : end);
        
        mr = tf(mr_coeffs_clean, 1);
        nr = tf(nr_coeffs_clean, 1);

        display(mr)
        display(nr)
    end
end

fprintf('\nAantal geldige combinaties: %d\n', length(geldige_combos));
nr2_mr2 = mr^2 - nr^2; %dit is inderdaad de teller van ro_kwad, :)


%% N12 door (6.25 (a)) om te vormen
fef =  2* sqrt(R1/R2) * T;
[num_coeffs, ~] = tfdata(fef, 'v');
[~, den_coeffs] = tfdata(fef, 'v');
N12 = num_coeffs
N12_pol = tf(N12,1)
%als check:
%den_coeffs = n + m ; KLOPT :)

%% Z-parameters berekenen met Table 6.2
%N12 is odd
z11 = (R1*(n-nr))/(m+mr)
z22 = (R2*(n+nr))/(m+mr)
z12 = sqrt(R1*R2)*(N12_pol/(m+mr))
zpk(z11)
zpk(z22)
zpk(z12)
%% ABCD matrices berekenen om zo te checken of Z-param juist zijn
Rl = 1;
Rs = RE_norm;
C1 = 76.27765065;
L1 = 0.01190231866;
L2 = 0.001207681337;
C2 = 912.0502424;
L3 = 0.0003295146798;
C3 = 3034.766162;

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
A = remove_common_factors(ABCD_final(1,1))
B = remove_common_factors(ABCD_final(1,2))
C = remove_common_factors((zpk(ABCD_final(2,1))))
D = remove_common_factors(ABCD_final(2,2))

z11 = remove_common_factors(zpk(-D/C))
z22 = remove_common_factors(zpk(-A/C))
z12 = remove_common_factors(zpk(-1/C))