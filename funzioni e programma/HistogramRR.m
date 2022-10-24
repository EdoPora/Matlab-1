function [] = HistogramRR(n_bins,RR1,RR2)

min_minRR = min([RR1,RR2]);
max_maxRR = max([RR1,RR2]);

bin_edgesRR = linspace(min_minRR, max_maxRR, n_bins);  

% bin_centersRR = bin_edgesRR(1:end-1) + mean(diff((bin_edgesRR))); % voglio un vettore con i centri dei vari bins, in modo che il valore del singolo bin cada centralmente
figure('Name','RR interval Histogram');
histogram(RR1,bin_edgesRR);
hold on
histogram(RR2, bin_edgesRR,'FaceColor','r');
legend('PDF RR interval 1', 'PDF RR interval 2');
print('HistogramRR','-dpng');

end

