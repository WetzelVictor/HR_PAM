close all;
clearvars;
clc;

%% Initial signal and spectral analysis

[x,fs] = audioread('../Mesures/Table2/Jeu2_dAddarioXL/Repetabilite/Repet_P2_Plectre_Down_Sol_15.wav');
Nx = length(x);
w = hann(Nx);
y = x.*w;

[sigMax,indMax] = max(abs(x));
noise = y(1:(indMax-1));
noiseFloor = rms(noise);
tail = y((indMax+(0.5*fs)):end);
Ntail = length(tail);

snr = db(sigMax/(noiseFloor));

Nfft = 2^nextpow2(Ntail);
tailSpec = fft(tail,Nfft);
dFfft = fs/Nfft;
fSpec = 0:dFfft:(fs-(1/Nfft));
tailSpec = tailSpec(1:(0.5*Nfft));
fSpec = fSpec(1:(0.5*Nfft));

f1 = frequence(tail,fs,190,200);
Bvect = 1e-7:1e-7:1e-4;
Nm = 25;
L_s = 0.655; % Vibrating length of string
dens_s = mean([7.75e3,8.05e3]); % Density of steel in kg/m^3
% dens_s = 1.15e3; % Density of nylon in kg/m^3
d_s = 4.4e-4; % d'Addario XL Steel G (4th) string
% d_s = 8.5e-4; % Alliance Savarez Nylong G (4th) string
rho_s = dens_s*(0.25*pi*d_s^2); % Linear mass density of string
I_s = (pi*d_s^4)/64; % Cylinder moment of inertia for string
E_s = 2e11; % Most used value for Young modulus of steel
% E_s = 2e9; % Most used value (low) for Young modulus of Nylon
T_s = (rho_s*(2*L_s*f1)^2)-(((pi^2)*E_s*I_s)/(L_s^2)); % Required tension
alpha = 1e-2;

cond = 1;
while cond == 1
    E = zeros(length(Bvect),1);
    for k = 1:length(Bvect)
        E(k) = abs(tailSpec(floor(f1/dFfft)+1))^2;
        for n = 1:Nm
            fn = (n/(2*L_s))*sqrt((T_s/rho_s)*(1+(Bvect(k)*(n^2))));
            fnMax = fn*(1+alpha);
            fnMin = fn*(1-alpha);
            [Emax,maxInd] = max(abs(tailSpec((floor(fnMin/dFfft)+1):(floor(fnMax/dFfft)+1))).^2);
            E(k) = E(k)+Emax;
        end
    end
    bInd = find(E==max(E));
    b = mean(Bvect(bInd));
    prevVals = [T_s,E_s];
    E_s = (b*prevVals(1)*(L_s^2))/((pi^2)*I_s);
    T_s = (rho_s*(2*L_s*f1)^2)-(((pi^2)*E_s*I_s)/(L_s^2));
    errE = sqrt((prevVals(2)-E_s)^2);
    errT = sqrt((prevVals(1)-T_s)^2);
    if (errE<=eps)||(errT<=eps)
        cond = 0;
    end
end

figure;
plot(fSpec,db(tailSpec));
xlim([0,5000]);

%% Extraction of HR damping coefficients by mode

load('../Traitements/repetabilite/plectre/P2_Down_Sol_15.mat');

fMes = fs*f;
fQ = zeros((size(fMes,1)*size(fMes,2)),1);
Qmes = zeros(size(delta));
for k = 1:size(fMes,2)
    fk = fMes(:,k);
    Qk = (-pi*fk)./(fs*delta(:,k));
    fQ(((k-1)*size(fMes,1)+1):(k*size(fMes,1))) = fk;
    Qmes(:,k) = Qk;
end

f1Max = f1*(1+alpha);
f1Min = f1*(1-alpha);
c = reshape(a,size(fQ));
fQ1 = fQ(fQ>=f1Min);
c = c(fQ>=f1Min);
f1mes = fQ1(fQ1<=f1Max);
c = c(fQ1<=f1Max);
[~,f1Ind] = max(abs(c));
f1 = f1mes(f1Ind);
% if not(isempty(f1mes))
%     [~,closest] = min(abs(f1mes-f1));
%     f1 = f1mes(closest);
% end

T_s = (rho_s*(2*L_s*f1)^2)-(((pi^2)*E_s*I_s)/(L_s^2)); % Required tension
B = ((pi^2)*E_s*I_s)/(T_s*(L_s^2));

Qmoy = mean(mean(Qmes));
Qstd = std(std(Qmes));
fQn = zeros(Nm,size(fMes,2));
Qn = zeros(Nm,size(fMes,2));
for n=1:Nm
    fn = (n/(2*L_s))*sqrt((T_s/rho_s)*(1+(B*(n^2))));
    fMax = fn*(1+alpha);
    fMin = fn*(1-alpha);
    for k = 1:size(fMes,2)
        fk = fMes(:,k);
        Qk = Qmes(:,k);
        ak = a(:,k);
        fkn = fk(fk>=fMin);
        Qk = Qk(fk>=fMin);
        ak = ak(fk>=fMin);
        fkn = fkn(fkn<=fMax);
        Qk = Qk(fkn<=fMax);
        ak = ak(fkn<=fMax);
        if not(isempty(fkn))
%             [~,closest] = min(abs(fkn-fn));
%             fkn = fkn(closest);
%             fkn_ind = find(fk==fkn);
%             Qkn = Qk(fkn_ind);
            fknInd = find(ak==max(ak));
            Qkn = mean(Qk(fknInd));
            if (Qkn>=Qmoy)&&(Qkn<=(Qmoy+Qstd))
                fQn(n,k) = mean(fkn(fknInd));
                Qn(n,k) = Qkn;
            end
        end
    end
end
fQn(Qn==0) = NaN;
Qn(Qn==0) = NaN;

% All extracted damping coefficients (comment out next section)
% fQn = reshape(fQn,[(Nm*size(fMes,2)),1]);
% Qn = reshape(Qn,[(Nm*size(fMes,2)),1]);
% nnnInd = find(not(isnan(Qn)));
% Qn = Qn(nnnInd);
% fQn = fQn(nnnInd);

% Mean damping per mode (comment out previous section)
QnMoy = nanmean(Qn,2);
fQnMoy = nanmean(fQn,2);

save('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_15','QnMoy','fQnMoy','E_s');

figure;
semilogx(fQnMoy,QnMoy,'*');
xlabel('Fréquence [log_{10}(Hz)]');
ylabel('Q');
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 16;
% fig = gcf;
% pos = fig.PaperPosition;
% fig.PaperSize = [pos(3),pos(4)];
% print('PickNylonG_DownP2_Qn','-dpdf');

%%

f1 = 196;
L_s = 0.655; % Vibrating length of string
d_sSt = 4.4e-4; % d'Addario XL Steel G (4th) string
% d_sNy = 8.5e-4; % Alliance Savarez Nylong G (4th) string
dens_sSt = mean([7.75e3,8.05e3]); % Density of steel in kg/m^3
% dens_sNy = 1.15e3; % Density of nylon in kg/m^3
rho_sSt = dens_sSt*(0.25*pi*d_sSt^2); % Linear mass density of string
% rho_sNy = dens_sNy*(0.25*pi*d_sNy^2); % Linear mass density of string
I_sSt = (pi*d_sSt^4)/64; % Cylinder moment of inertia for string
% I_sNy = (pi*d_sNy^4)/64; % Cylinder moment of inertia for string
eta = 1.85e-5; % Air viscosity coefficient at 300K
rho_a = 1.184; % Air density at ~300K

lSt1 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_1');
lSt2 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_2');
lSt3 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_3');
lSt4 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_4');
lSt5 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_5');
lSt6 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_6');
lSt7 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_7');
lSt8 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_8');
lSt9 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_9');
lSt10 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_10');
lSt11 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_11');
lSt12 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_12');
lSt13 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_13');
lSt14 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_14');
lSt15 = load('Table2/Repetabilite/Steel_PickG_DownP2_QnMoy_15');

QmoySt = [lSt1.QnMoy,lSt2.QnMoy,lSt3.QnMoy,lSt4.QnMoy,lSt5.QnMoy,lSt6.QnMoy,lSt7.QnMoy,...
    lSt8.QnMoy,lSt9.QnMoy,lSt10.QnMoy,lSt11.QnMoy,lSt12.QnMoy,lSt13.QnMoy,lSt14.QnMoy,lSt15.QnMoy];
fQmoySt = [lSt1.fQnMoy,lSt2.fQnMoy,lSt3.fQnMoy,lSt4.fQnMoy,lSt5.fQnMoy,lSt6.fQnMoy,lSt7.fQnMoy,...
    lSt8.fQnMoy,lSt9.fQnMoy,lSt10.fQnMoy,lSt11.fQnMoy,lSt12.fQnMoy,lSt13.fQnMoy,lSt14.fQnMoy,lSt15.fQnMoy];
errQmesSt = std(QmoySt,0,2,'omitnan');
QmesSt = nanmean(QmoySt,2);
fQmesSt = nanmean(fQmoySt,2);

% lNy1 = load('Table1/Nylon/WireG_0degP1_QnMoy');
% lNy2 = load('Table1/Nylon/WireG_0degP2_QnMoy');
% lNy3 = load('Table1/Nylon/WireG_0degP3_QnMoy');
% lNy4 = load('Table1/Nylon/WireG_45degP1_QnMoy');
% lNy5 = load('Table1/Nylon/WireG_45degP2_QnMoy');
% lNy6 = load('Table1/Nylon/WireG_45degP3_QnMoy');
% lNy7 = load('Table1/Nylon/WireG_90degP1_QnMoy');
% lNy8 = load('Table1/Nylon/WireG_90degP2_QnMoy');
% lNy9 = load('Table1/Nylon/WireG_90degP3_QnMoy');
% lNy10 = load('Table1/Nylon/PickG_DownP1_QnMoy');
% lNy11 = load('Table1/Nylon/PickG_DownP2_QnMoy');
% lNy12 = load('Table1/Nylon/PickG_DownP3_QnMoy');
% lNy13 = load('Table1/Nylon/PickG_UpP1_QnMoy');
% lNy14 = load('Table1/Nylon/PickG_UpP2_QnMoy');
% lNy15 = load('Table1/Nylon/PickG_UpP3_QnMoy');
% 
% QmoyNy = [lNy1.QnMoy,lNy2.QnMoy,lNy3.QnMoy,lNy4.QnMoy,lNy5.QnMoy,lNy6.QnMoy,lNy7.QnMoy,...
%     lNy8.QnMoy,lNy9.QnMoy,lNy10.QnMoy,lNy11.QnMoy,lNy12.QnMoy,lNy13.QnMoy,lNy14.QnMoy,lNy15.QnMoy];
% fQmoyNy = [lNy1.fQnMoy,lNy2.fQnMoy,lNy3.fQnMoy,lNy4.fQnMoy,lNy5.fQnMoy,lNy6.fQnMoy,lNy7.fQnMoy,...
%     lNy8.fQnMoy,lNy9.fQnMoy,lNy10.fQnMoy,lNy11.fQnMoy,lNy12.fQnMoy,lNy13.fQnMoy,lNy14.fQnMoy,lNy15.fQnMoy];
% QmesNy = nanmean(QmoyNy,2);
% fQmesNy = nanmean(fQmoyNy,2);

E_sSt = [lSt1.E_s,lSt2.E_s,lSt3.E_s,lSt4.E_s,lSt5.E_s,lSt6.E_s,lSt7.E_s,...
    lSt8.E_s,lSt9.E_s,lSt10.E_s,lSt11.E_s,lSt12.E_s,lSt13.E_s,lSt14.E_s,lSt15.E_s];
errE_sSt = std(E_sSt);
E_sSt = mean(E_sSt);

% E_sNy = [lNy1.E_s,lNy2.E_s,lNy3.E_s,lNy4.E_s,lNy5.E_s,lNy6.E_s,lNy7.E_s,...
%     lNy8.E_s,lNy9.E_s,lNy10.E_s,lNy11.E_s,lNy12.E_s,lNy13.E_s,lNy14.E_s,lNy15.E_s];
% E_sNy = mean(E_sNy);

% T_sSt = (rho_sSt*(2*L_s*f1)^2)-(((pi^2)*E_sSt*I_sSt)/(L_s^2));
% T_sNy = (rho_sNy*(2*L_s*f1)^2)-(((pi^2)*E_sNy*I_sNy)/(L_s^2));

% Qth = 3e4; % Thermo-elastic dissipation factor
% del_veSt = 6e-3; % Visco-elastic dissipation factor
% del_veNy = 5e-2;

% fSyn = round(round(min(fQmesSt)):round(max(fQmesSt)));
% omega = 2*pi*fSyn;
% Mst = (d_sSt/4)*sqrt(omega/eta);
% Rst = (pi^2)*rho_a*fSyn*((d_sSt^2)/2).*((sqrt(2)./Mst)+(1./(2*(Mst.^2))));
% QsynSt = 1./((Rst./(omega*rho_sSt))+(((omega.^2)*rho_sSt*E_sSt*I_sSt*del_veSt)/(T_sSt^2))+(1/Qth));

% Mny = (d_sNy/4)*sqrt(omega/eta);
% Rny = (pi^2)*rho_a*fSyn*((d_sNy^2)/2).*((sqrt(2)./Mny)+(1./(2*(Mny.^2))));
% QsynNy = 1./((Rny./(omega*rho_sNy))+(((omega.^2)*rho_sNy*E_sNy*I_sNy*del_veNy)/(T_sNy^2))+(1/Qth));

figure;
ax = gca;
errplt = errorbar(fQmesSt,QmesSt,errQmesSt,'*','LineWidth',2);
errplt.MarkerEdgeColor = [0.8500,0.3250,0.0980];
errplt.MarkerFaceColor = [0.8500,0.3250,0.0980];
xlabel('Fréquence [log_{10}(Hz)]');
ylabel('Q');
ax.XScale = 'log';
xlim([min(fQmesSt),max(fQmesSt)]);
ax.LineWidth = 2;
ax.FontSize = 16;
fig = gcf;
pos = fig.PaperPosition;
fig.PaperSize = [pos(3),pos(4)];
print('Table2/RepetabiliteQ','-dpdf');