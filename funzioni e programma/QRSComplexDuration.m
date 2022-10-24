function [meandurQRS,durQRS] = QRSComplexDuration(Numfin,Nm,matrixQRS,timeQRS,k,K)

 % calcolo la durata temporale tra il secondo picco prima e dopo rispetto al picco R e
 % lo stimo circa pari alla QRS duration
 
durQRS = zeros([1,Numfin]);

% detrend
newmatrixQRS = zeros([Numfin,Nm]);

for i=1:Numfin
    Nk = round(numel(matrixQRS(i,:))/k);
    D = detrend(reshape(matrixQRS(i,:),Nk,[]));
    newmatrixQRS(i,:) = reshape(D,[1,numel(matrixQRS(i,:))]);
end

% abs
ABS = zeros([Numfin,Nm]);

for i=1:Numfin
    ABS(i,:) = abs(newmatrixQRS(i,:));
end

% findpeaks (selezione solo picchi R e minimi S)

sogliaQRS = 0.013;   % rapporto tra QRS medio e RS medio nel primo e secondo caso ottenuto da grafici QRS medi 
                     % dalla durata RS posso stimare quella del complesso totale 

for i=1:Numfin
    [picchi,indici] = findpeaks(ABS(i,:), 'MinPeakHeight',sogliaQRS,'NPeaks',2,'MinPeakDistance',8);
    durQRS(1,i) = (timeQRS(indici(2)) - timeQRS(indici(1)))*K;
end
% volendo si può aggiungere un controllo in più con un ciclo for che prende
% solo le durate minori di 0,12 sec e ne fa la media, in modo da eliminare
% gli effetti del rumore ed eventuali errori 

% meandurQRS1 = mean(durQRS1);
% meandurQRS2 = mean(durQRS2);

sumQRS = 0;
cont = 0;
for i=1:Numfin
    if(durQRS(1,i) <= 0.12)
        tmt = sumQRS;
        sumQRS = tmt + durQRS(1,i);
        contt = cont;
        cont = contt + 1;
    end
end

meandurQRS = sumQRS/cont;

end

