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

fs1_n = fs1/f_center;
fs2_n = fs2_new/f_center;
fd1_n = fd1/f_center;
fd2_n = fd2/f_center;
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

H = tf(TF_teller, TF_noemer); %nog te schalen!!!
display(H)
%% Schaling H met spanningsdeler (6.22)
R2 = RL_norm; %1
R1 = RE_norm; %Rn = 0.23

T = (R2/(R1 + R2)) * H;
display(T)

%% (6.25) (6.26)
%T(s)*T(-s)
% Index:         s^3    s^2     s^1     s^0
noemer_s     = [ 1,     2.953,  4.361,  3.219];
noemer_min_s = [-1,     2.953, -4.361,  3.219];

noemer_kwad = conv(noemer_s, noemer_min_s);

disp(noemer_kwad)

teller_kwad = 2.664*2.664;

T_kwad = tf(teller_kwad, noemer_kwad);

Ro_kwad = 1 - 4*(R1/R2) * T_kwad;

