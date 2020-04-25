% weight
function w=weight(p,q,mc,g)
    w = 1/(1+exp(-g*(abs(p-q)-mc)));
end

