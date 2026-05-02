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
zplane([], polen.'); % Toont de polen op de ellips
title('Pool-Nulpunten diagram (Chebyshev Orde 5)');
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
H_bpf
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

Ro_kwad = 1 - 4*(R1/R2) * T_kwad_bd
%% n en m berekenen
%n en m berkenen kan ook al sneller uit T(s) met formule (6.25 (a))
% 1. Haal de stabiele noemer (m + n)
[~, den_coeffs] = tfdata(Ro_kwad, 'v');
all_poles = roots(den_coeffs);
stable_poles = all_poles(real(all_poles) < -1e-5); %pak de neg polen
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
n2_m2 = m^2 - n^2 %dit is inderdaad de noemer van ro_kwad, :)

%% nr en mr bereken

[num_coeffs, ~] = tfdata(Ro_kwad, 'v');
all_zeros = roots(num_coeffs);
stable_zeros = all_zeros(real(all_zeros) < 0) %We kiezen hier even ez de linkse nullen. 
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
%als check:
%den_coeffs = n + m ; KLOPT :)

%% Y-parameters berekenen met Table 6.2

y11 = (1/R1) * (m + mr) / (n - nr) %teller en n
y22 = (1/R2) * (m - mr) / (n - nr)
y12 = -(1/sqrt(R1*R2)) * (N12 / (n-nr))

[num_coeffs, ~] = tfdata(y11, 'v');
[~, den_coeffs] = tfdata(y11, 'v');
teller_y11 = num_coeffs;
noemer_y11 = den_coeffs;
zeros_y11 = roots(teller_y11) 
polen_y11 = roots(noemer_y11)

