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

%% Troviamo la frequenza cardiaca
bpm1=compute_heart_rate(Maxindex1,tmax);
bpm2=compute_heart_rate(Maxindex2,tmax);

%% intervallo RR
RR1 = rand([1 length(Maxindex1)]);
for i=1:length(Maxindex1)
    if (i < length(Maxindex1))
       RR1(i) = time(Maxindex1(i+1)) - time(Maxindex1(i));
    else
        RR1(i) = 0;
    end
end
meanRR1 = mean(RR1);
RR2 = rand([1 length(Maxindex2)]);
for i=1:length(Maxindex2)
    if (i < length(Maxindex2))
       RR2(i) = time(Maxindex2(i + 1)) - time(Maxindex2(i));
    else
        RR2(i) = 0;
    end
end
meanRR2 = mean(RR2);
disp('intervallo RR medio 1')
disp(meanRR1)
disp('intervallo RR medio 2')
disp(meanRR2)
diffRR = meanRR2 - meanRR1;
disp('differenza medie intervalli RR')
disp(diffRR)
meanRR = (meanRR2 + meanRR1)/2;
%% QRS Complex
% funzione per il calcolo della durata del complesso QRS e per la
% visualizzazione del complesso QRS medio.
M = meanRR;   % intervallo temporale in sec finestra per QRS (0.08-0.12 sec)
Numfin1 = length(values1);  % numero finestre totale
Numfin2 = length(values2);
Nm = tmax/M; % numero campioni per ogni finestra 
matrixQRS1 = rand([Numfin1 Nm]);
for i=1:Numfin1
    for j=1:Nm
            matrixQRS1(i,j) = data1(Maxindex1(i) - Nm/2 + j);
    end
end
meanQRS1 = mean(matrixQRS1);
timeQRS = linspace((-M/2),(M/2),Nm);
f7 = figure('Name','QRS Complex 1');
plot(timeQRS, matrixQRS1, ':');
hold on
legend('Singles QRS Complexes 1');
plot(timeQRS,(meanQRS1)', 'k','DisplayName', 'Mean QRS Complex 1');
hold off  % plotto il complesso medio QRS 
% Note: la forma è corretta, individuo i vari picchi P, Q, R, S e T;
        % non so bene che intervallo M prendere
        % manca la rappresentazione del complesso QRS medio della seconda
        % ECG
matrixQRS2 = rand([Numfin2 Nm]);
for i=1:Numfin2
    for j=1:Nm
            matrixQRS2(i,j) = data2(Maxindex2(i) - Nm/2 + j);
    end
end
meanQRS2 = mean(matrixQRS2);
f8 = figure('Name','QRS Complex 2');
plot(timeQRS,matrixQRS2, ':');
hold on
legend('Singles QRS Complexes 2');
plot(timeQRS,(meanQRS2)', 'k','DisplayName', 'Mean QRS Complex 2');
hold off
%% QRS Complex Duration
 % trovo intervallo tra R e S e lo stimo pari a metà della durata QRS
 % QRS = RS*2; stima durata QRS
 %% Histogram
 % istogramma