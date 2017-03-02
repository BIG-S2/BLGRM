function [aq] = trace1(A)
    n=size(A,3);
    aq=0;
    for i=1:n
        aq=aq+trace(A(:,:,i));
    end
    
end
