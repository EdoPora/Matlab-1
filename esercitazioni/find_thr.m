function [S] = find_thr(newdata,m)

S=mean(newdata)+m*std(newdata);
end