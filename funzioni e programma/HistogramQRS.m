function [] = HistogramQRS(n_bins,durQRS1,durQRS2)

min_minQRS = min([durQRS1,durQRS2]);
max_maxQRS = max([durQRS1,durQRS2]);

bin_edgesQRS = linspace(min_minQRS, max_maxQRS, n_bins);  

% bin_centersRR = bin_edgesRR(1:end-1) + mean(diff((bin_edgesRR))); % voglio un vettore con i centri dei vari bins, in modo che il valore del singolo bin cada centralmente
figure('Name','RR interval Histogram');
histogram(durQRS1,bin_edgesQRS);
hold on
histogram(durQRS2, bin_edgesQRS,'FaceColor','r');
legend('PDF RR interval 1', 'PDF RR interval 2');
print('HistogramQRS','-dpng');

end