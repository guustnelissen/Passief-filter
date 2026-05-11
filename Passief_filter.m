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

%% Transformatie van LPP polen naar BPF polen

% Initialiseer de array voor de BPF polen
% Omdat elke LPP pool 2 BPF polen genereert, hebben we 2*orde polen
polen_bpf = zeros(1, 2*orde);

% Omega_0 is de genormaliseerde centrumfrequentie (dit is altijd 1)
omega_0_norm = 1;

% Loop door elke LPP pool en bereken de twee bijbehorende BPF polen
for i = 1:orde
    p_m = polen(i);
    
    % De kwadratische vergelijking is: s^2 - (B * p_m) * s + omega_0^2 = 0
    % We gebruiken de abc-formule om s te vinden: 
    % s = [ (B*p_m) +/- sqrt( (B*p_m)^2 - 4*omega_0^2 ) ] / 2
    
    term_b = B * p_m;
    discriminant = term_b^2 - 4 * omega_0_norm^2;
    
    % Bereken de twee BPF polen voor deze ene LPP pool
    s1 = (term_b + sqrt(discriminant)) / 2;
    s2 = (term_b - sqrt(discriminant)) / 2;
    
    % Sla ze op in de nieuwe array
    polen_bpf(2*i - 1) = s1;
    polen_bpf(2*i) = s2;
end

disp('De getransformeerde polen voor de Banddoorlaat (BPF) zijn:');
disp(polen_bpf.');

% Transferfunctie van de Banddoorlaat opstellen

% 1. Nulpunten (Zeros) bepalen
% Een LPP (zoals Butterworth) heeft geen nulpunten. Bij transformatie naar
% banddoorlaat ontstaan er 'orde' (n) nulpunten op de oorsprong (s = 0).
nulpunten_bpf = zeros(orde, 1); 

% Zorg dat polen en nulpunten kolomvectoren zijn voor MATLAB's zp2tf functie
polen_bpf_kolom = polen_bpf.';

% 2. Versterkingsfactor (Gain K) bepalen
% De theoretische gain voor een getransformeerd all-pole filter is:
% K = (Bandbreedte^orde) * (Product van de absolute waarden van de LPP polen)
gain_K = (B^orde) * prod(abs(polen));

% 3. Omzetten van Zeros, Poles en Gain naar Polynomen (teller en noemer)
[teller, noemer] = zp2tf(nulpunten_bpf, polen_bpf_kolom, gain_K);

% 4. Maak het transferfunctie object aan
% (Zorg dat je de Control System Toolbox geïnstalleerd hebt voor de 'tf' functie)
H_bpf = tf(teller, noemer);

disp('De transferfunctie van het Banddoorlaatfilter is:');
display(H_bpf)
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

%% (6.25) (6.26)
%T(s)*T(-s)
% Index:         s^3    s^2     s^1     s^0
noemer_s     = [ 1,     2.953,  4.361,  3.219];
noemer_min_s = [-1,     2.953, -4.361,  3.219];

noemer_kwad = conv(noemer_s, noemer_min_s);

disp(noemer_kwad)

teller_kwad = 2.664*2.664;

T_kwad = tf(teller_kwad, noemer_kwad);

Ro_kwad = 1 - 4*(R1/R2) * T_kwad
%%
T_kwad_bd = T * T_min_s;
display(T_kwad_bd)

Ro_kwad = 1 - (4*(R1/R2) * T_kwad_bd) %het lijkt alsof teller en noemer hetzelfde is maar dit is niet, matlab rondt te veel af
%% n en m berekenen
%n en m berkenen kan ook al sneller uit T(s) met formule (6.25 (a))
% 1. Haal de stabiele noemer (m + n)
[~, den_coeffs] = tfdata(Ro_kwad, 'v');
all_poles = roots(den_coeffs);
stable_poles = all_poles(real(all_poles) < -1e-5); %pak de neg polen in het LHV
D_s = poly(stable_poles) % Dit is m + n

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

%% nr en mr bereken

[num_coeffs, ~] = tfdata(Ro_kwad, 'v');
all_zeros = roots(num_coeffs);
stable_zeros = all_zeros(real(all_zeros) < 0) %We kiezen hier even ez de linkse nullen. 
% Dit kunnen ook andere zijn, dit gaan mss ook zo moeten want je gaat mss niet de goede K factor vinden met de nullen dat je hebt gekozen
F_s = poly(stable_zeros) % Dit is mr + nr

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

mr = tf(mr_coeffs, 1)
nr = tf(nr_coeffs, 1)
nr2_mr2 = mr^2 - nr^2; %dit is inderdaad de teller van ro_kwad, :)

%% Z_in en z11 bepalen (Case B voor Banddoorlaat)
rho = tf(F_s,D_s);
syms s

% Zet de gevonden polynomen om naar symbolische vorm (makkelijker voor Cauer)
m_sym = poly2sym(m_coeffs, s);
n_sym = poly2sym(n_coeffs, s);
mr_sym = poly2sym(mr_coeffs, s);
nr_sym = poly2sym(nr_coeffs, s);

R = 50; % Je bronweerstand

% Bereken m1, n1, m2, n2 volgens Case B (zie Uitleg 3de categorie.pdf, Tabel 6-1)
% z11 moet een oneven/even functie zijn
m1 = simplify(m_sym - mr_sym); % Of m + mr, afhankelijk van het ± teken in je boek
n1 = simplify(n_sym + nr_sym); % Of n - nr
m2 = simplify(m_sym + mr_sym); %even check :)
n2 = simplify(n_sym - nr_sym); %oneven check :)

% Controleer of z11 (n1/m2) de juiste nullen in de oorsprong heeft voor een banddoorlaat
z11 = simplify((n1) / (m2));

disp('De z11 voor de Cauer synthese is:');
disp(vpa(z11, 5));

z12 = sqrt(R) * (sqrt(-(m1*m2 - n1*n2)) / m2);
disp('De z12 voor de Cauer synthese is:');
disp(vpa(z12, 5));

% --- Betrouwbare methode voor Z12 layout ---

% 1. Definieer de teller en noemer symbolisch
teller_z12_kwadraat = simplify(R * -(m1*m2 - n1*n2));
noemer_z12_sym = m2;

% 2. Pak de coëfficiënten van de noemer (dit gaat meestal goed)
den_coeffs = double(poly2sym2poly(noemer_z12_sym, s));

% 3. Voor de teller: we berekenen de coëfficiënten van de term ONDER de wortel
teller_poly_onder_wortel = double(poly2sym2poly(teller_z12_kwadraat, s));

% 4. Omdat we weten dat voor een banddoorlaat z12 de vorm K*s^n / noemer heeft,
% zoeken we de enige coëfficiënt in de teller die NIET nul is (of de grootste).
K_val = sqrt(max(abs(teller_poly_onder_wortel))); 

% 5. Bouw de teller vector handmatig op basis van de graad
% Voor een 6de orde banddoorlaat zit de s^3 term in het midden
num_coeffs = zeros(1, length(den_coeffs));
macht_van_s = 3; % Pas dit aan als je teller een andere macht van s heeft (bijv. s^1 of s^2)
num_coeffs(end - macht_van_s) = K_val;

% 6. Maak de TF aan
Z12_pretty = tf(num_coeffs, den_coeffs);

% Zorg dat de vectoren exact als 'double' rijen worden doorgegeven
Z12_final = tf(double(num_coeffs), double(den_coeffs));

% Forceer de weergave
fprintf('\nZ12 weergave:\n');
Z12_final  % <--- GEEN punt-komma hier!

% --- Hulpfunctie (plaats deze onderaan je script of voer dit uit) ---
function p = poly2sym2poly(sym_expr, var)
    % Forceert een symbolische expressie naar een numerieke vector
    c = coeffs(expand(vpa(sym_expr, 8)), var, 'All');
    p = double(c);
end

% --- Z11 omzetten naar mooie TF layout ---

% 1. Haal de teller (n1) en noemer (m2) op uit de symbolische z11
% Gebruik vpa om breuken om te zetten naar decimalen voor de weergave
num_z11_sym = n1;
den_z11_sym = m2;

% 2. Gebruik de hulpfunctie om de coëfficiënten naar numerieke vectoren te halen
% Dit zorgt ervoor dat kleine symbolische restjes verdwijnen
num_z11_vec = double(poly2sym2poly(num_z11_sym, s));
den_z11_vec = double(poly2sym2poly(den_z11_sym, s));

% 3. Maak het Transfer Function object aan
Z11_final = tf(num_z11_vec, den_z11_vec);

% 4. Forceer de visuele weergave (zoals in image_29c8af.png)
fprintf('\nDe z11 voor de Cauer synthese (TF weergave):\n');
Z11_final % <--- GEEN punt-komma voor de layout
%% N12 door (6.25 (a)) om te vormen
fef =  2* sqrt(R1/R2) * T;
[num_coeffs, ~] = tfdata(fef, 'v');
[~, den_coeffs] = tfdata(fef, 'v');
N12 = num_coeffs
N12_pol = tf(N12,1)
%als check:
%den_coeffs = n + m ; KLOPT :)

%% Y-parameters berekenen met Table 6.2
%N12 is odd
y11 = (1/R1) * (n + nr) / (m - mr) %teller en n
y22 = (1/R2) * (n - nr) / (m - mr)
y12 = -(1/sqrt(R1*R2)) * (N12_pol / (m-mr))

[num_coeffs, ~] = tfdata(y11, 'v');
[~, den_coeffs] = tfdata(y11, 'v');
teller_y11 = num_coeffs;
noemer_y11 = den_coeffs;
zeros_y11 = round(roots(teller_y11),3) 
polen_y11 = round(roots(noemer_y11),3)

z11 = (R1*(n-nr))/(m+mr);
z22 = (R2*(n+nr))/(m+mr);
z12 = sqrt(R1*R2)*(N12_pol/(m+mr));
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
ABCD_final = ABCD

% Z-parameters extraheren uit de ABCD matrix
% z11 = A/C, z12 = 1/C, z22 = D/C
A = remove_common_factors(ABCD_final(1,1))
B = remove_common_factors(ABCD_final(1,2))
C = remove_common_factors((zpk(ABCD_final(2,1))))
D = remove_common_factors(ABCD_final(2,2))

z11 = remove_common_factors(zpk(-D/C))
z22 = remove_common_factors(zpk(-A/C))
z12 = remove_common_factors(zpk(-1/C))
