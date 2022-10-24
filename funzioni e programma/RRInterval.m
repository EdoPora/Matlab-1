function [meanRR,RR] =RRInterval(Maxindex,time)

RR = rand([1 length(Maxindex)]);

for i=1:length(Maxindex)
    if (i < length(Maxindex))
       RR(i) = time(Maxindex(i+1)) - time(Maxindex(i));
    else
        RR(i) = 0;
    end
end

meanRR = mean(RR);

end

