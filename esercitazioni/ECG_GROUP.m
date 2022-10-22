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
M1 = meanRR1;   % intervallo temporale in sec finestra 
Numfin1 = length(values1);  % numero finestre totale
Numfin2 = length(values2);
Nm1 = T/M1; % numero campioni per ogni finestra 
matrixQRS1 = rand([Numfin1 Nm1]);
grayColor = [.7 .7 .7];
for i=1:Numfin1
    for j=1:Nm1
            matrixQRS1(i,j) = data1(Maxindex1(i) - Nm1/2 + j);
    end
end
meanQRS1 = mean(matrixQRS1);
timeQRS1 = linspace((-M1/2),(M1/2),Nm1);
f7 = figure('Name','QRS Complex 1');
plot(timeQRS1, matrixQRS1, ':', 'Color', grayColor);
hold on
legend('Singles QRS Complexes 1');
plot(timeQRS1,(meanQRS1)', 'k','DisplayName', 'Mean QRS Complex 1', 'LineWidth', 3);
hold off  % plotto il complesso medio QRS 
% Note: la forma è corretta, individuo i vari picchi P, Q, R, S e T;
        % non so bene che intervallo M prendere
        % manca la rappresentazione del complesso QRS medio della seconda
        % ECG
M2 = meanRR2;
Nm2 = T/M2;
matrixQRS2 = rand([Numfin2 Nm2]);
for i=1:Numfin2
    for j=1:Nm2
            matrixQRS2(i,j) = data2(Maxindex2(i) - Nm2/2 + j);
    end
end
meanQRS2 = mean(matrixQRS2);
timeQRS2 = linspace((-M2/2),(M2/2),Nm2);
f8 = figure('Name','QRS Complex 2');
plot(timeQRS2,matrixQRS2, ':', 'Color', grayColor);
hold on
legend('Singles QRS Complexes 2');
plot(timeQRS2,(meanQRS2)', 'k','DisplayName', 'Mean QRS Complex 2', 'LineWidth', 3);
hold off
%% QRS Complex Duration
 % trovo intervallo tra R e S e lo stimo pari a metà della durata QRS
 % QRS = RS*2; stima durata QRS
MaxQRS1 = max(meanQRS1);
MinQRS1 = min(meanQRS1);
for i=1:Nm1
    if(meanQRS1(i) == MaxQRS1)
        Imax1 = i;
    end
end
for i=1:Nm1
    if(meanQRS1(i) == MinQRS1)
        Imin1 = i;
    end
end
durQRS1 = (timeQRS1(Imin1) - timeQRS1(Imax1))*2;
disp('Mean Duration QRS Complex 1 [sec]')
disp(durQRS1)
MaxQRS2 = max(meanQRS2);
MinQRS2 = min(meanQRS2);
for i=1:Nm2
    if(meanQRS2(i) == MaxQRS2)
        Imax2 = i;
    end
end
for i=1:Nm2
    if(meanQRS2(i) == MinQRS2)
        Imin2 = i;
    end
end
durQRS2 = (timeQRS2(Imin2) - timeQRS2(Imax2))*2;
disp('Mean Duration QRS Complex 2 [sec]')
disp(durQRS2)

 %% Histogram
 % istogramma
 
