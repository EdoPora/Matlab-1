function [meanQRS,timeQRS,matrixQRS,Numfin,Nm] = QRSInterval(Maxindex,meanRR,values,tmax,data)

M = meanRR;   % intervallo temporale in sec finestra 
Numfin = length(values);  % numero finestre totale
Nm = tmax/M; % numero campioni per ogni finestra 
matrixQRS = rand([Numfin Nm]);

for i=1:Numfin
    for j=1:Nm
            matrixQRS(i,j) = data(Maxindex(i) - Nm/2 + j);
    end
end

meanQRS = mean(matrixQRS);
timeQRS = linspace((-M/2),(M/2),Nm);

end

