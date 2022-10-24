% USE CASE 1: ECG
% plotto le due distribuzoni di dati forniti
close all
load('data.mat');
N=length(data1); %trovo quanti dati sono
figure % apro unna finestra grafica
tmax=N/fs; %trovo il tempo di durata della ECG

% plotto i dati raw

subplot(2,1,1)
plot(time,data1);
xlabel("time[sec]");
ylabel("amplitude"); 
title("Raw Data1");

subplot(2,1,2)
plot(time ,data2);
xlabel("time[sec]");
ylabel("amplitude"); 
title("Raw Data2");
print('DataRaw','-dpng');

% detrendo

k =0.5; % [s]
nrfin=tmax/k;
newdata1=lineardetrend(data1,N,nrfin);
newdata2=lineardetrend(data2,N,nrfin);

%plotto i dati detrendati 

figure
subplot(2,1,1)
plot(time,newdata1);
xlabel("time[sec]");
ylabel("amplitude"); 
title("Detrended Data1");

subplot(2,1,2)
plot(time,newdata2);
xlabel("time[sec]");
ylabel("amplitude"); 
title("Detrended Data2");
print('DataDetrend','-dpng');

%% Find Peaks
% identifichiamo i picchi significativi del segnale

m = 6;  % valore del moltiplicando
S=find_thr(newdata1,m); %%  calcolo valore di soglia

[values2, Maxindex2] = find_ecg_peak(newdata2,S);  
[values1, Maxindex1] = find_ecg_peak(newdata1,S);  %% vettore dei picchi

% plotto il risultato

figure;
subplot(2,1,1)
plot(time,newdata1);
xlabel("time[sec]");
ylabel("amplitude"); 
title("FindpeaksData");
hold on;
scatter(Maxindex1/double(fs),values1,'filled');


subplot(2,1,2)
plot(time,newdata2);
xlabel("time[sec]");
ylabel("amplitude"); 
title("FindpeaksData");
hold on;
scatter(Maxindex2/double(fs),values2,'filled');
print('DataPeaks','-dpng');

%% Troviamo la frequenza cardiaca
bpm1=compute_heart_rate(Maxindex1,tmax);
bpm2=compute_heart_rate(Maxindex2,tmax);

%% intervallo RR

[meanRR1,RR1] = RRInterval(Maxindex1,time); %stimo la durata dell'intervallo RR dei singoli campioni RR1 e RR2
[meanRR2,RR2] = RRInterval(Maxindex2,time);

disp('intervallo RR medio 1 [sec]')
disp(meanRR1)
disp('intervallo RR medio 2 [sec]')
disp(meanRR2)
diffRR = meanRR2 - meanRR1;
disp('differenza medie intervalli RR [sec]')
disp(diffRR)
meanRR = (meanRR2 + meanRR1)/2;

%% QRS Complex
% funzione per il calcolo della durata del complesso QRS e per la
% visualizzazione del complesso QRS medio.

[meanQRS1,timeQRS1,matrixQRS1,Numfin1,Nm1] = QRSInterval(Maxindex1,meanRR1,values1,tmax,data1);

grayColor = [.7 .7 .7];
figure('Name','QRS Complex 1');
plot(timeQRS1, matrixQRS1, ':', 'Color', grayColor);
hold on
legend('Singles QRS Complexes 1');
plot(timeQRS1,(meanQRS1)', 'k','DisplayName', 'Mean QRS Complex 1', 'LineWidth', 3);
hold off % plotto il complesso medio QRS 
print('DataQRS1','-dpng'); %stampo l'immagine per aggiungerla alla relazione

% Note: la forma è corretta, individuo i vari picchi P, Q, R, S e T;
        % non so bene che intervallo M prendere
        % manca la rappresentazione del complesso QRS medio della seconda
        % ECG
        

[meanQRS2,timeQRS2,matrixQRS2,Numfin2,Nm2] = QRSInterval(Maxindex2,meanRR2,values2,tmax,data2);

figure('Name','QRS Complex 2');
plot(timeQRS2,matrixQRS2, ':', 'Color', grayColor);
hold on
legend('Singles QRS Complexes 2');
plot(timeQRS2,(meanQRS2)', 'k','DisplayName', 'Mean QRS Complex 2', 'LineWidth', 3);
hold off
print('DataQRS2','-dpng');

%% QRS Complex Duration
% trovo i picchi con findpeaks (senza ribaltare) per ogni finestra di campioni

k1 = 5; % campioni per finestra detrend
k2 = 3;

K1 = 1.73;   % rapporto tra QRS medio e RS medio nel primo e secondo caso ottenuto da grafici QRS medi 
K2 = 1.82;   % dalla durata RS posso stimare quella del complesso totale 

[meandurQRS1,durQRS1] = QRSComplexDuration(Numfin1,Nm1,matrixQRS1,timeQRS1,k1,K1);
[meandurQRS2,durQRS2] =  QRSComplexDuration(Numfin2,Nm2,matrixQRS2,timeQRS2,k2,K2);

disp('Mean QRS Complex Duration 1 [sec] 1')
disp(meandurQRS1)
disp('Mean QRS Complex Duration 2 [sec] 1')
disp(meandurQRS2)

 %% Histogram
 
 % istogramma intervallo RR
n_binsRR = 100; 
HistogramRR(n_binsRR,RR1,RR2);

 % istogramma QRS duration
n_binsQRS = 50; 
HistogramQRS(n_binsQRS,durQRS1,durQRS2);

% written by Davide, Edo, Anna, Matte