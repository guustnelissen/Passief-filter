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

fs = abs((fs1_n - (1/fs1_n))/B)
fss = abs((fs2_n - (1/fs2_n))/B); %ff testen of deze gelijk zijn (is dus wel)
%de abs is nodig omdat fs_n kleiner is dan 1 en dan krijgen we een neg
%resultaat maar de abs-waarde is gelijk

fd = abs((fd1_n - (1/fd1_n))/B);
fdd = abs((fd2_n - (1/fd2_n))/B);

k = fd/fs
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

TF_teller = real(prod(-polen));

TF_noemer = real(poly(polen));

H = tf(TF_teller, TF_noemer); %Dit is voor de laagdoorlaat
display(H)
%% Schaling H met spanningsdeler (6.22)
R2 = RL_norm; %1
R1 = RE_norm; %Rn = 0.23

T = (R2/(R1 + R2)) * H;
display(T)
figure
bode(T);
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
