function [bpm] = compute_heart_rate(Maxindex,tmax)

fs=length(Maxindex)/double(tmax);
bpm=fs*60;

end