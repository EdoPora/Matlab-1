function [newdata] = lineardetrend(data,N,nrfin)

Nk=N/nrfin;
newmatrix = reshape(data,[Nk,nrfin]);   
D = detrend(newmatrix);  
newdata = reshape(D,[1,N]);

end

