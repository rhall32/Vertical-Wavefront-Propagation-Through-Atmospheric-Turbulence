function m=findClosestMultiples(n,actual)
    mult=[];
    for i = 1:n
        if mod(n,i)==0
            mult=[mult,i];
        end
    end
    [a, m] = min(abs(mult-actual));
    m= mult(m);
end

            
            