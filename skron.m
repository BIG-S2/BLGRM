function [AB] = skron(A,B)

    AI=skronAI(A);
    BI=skronAI(B);
    ABI=skronAI(A*B);
    AB=2*AI*BI-ABI;
end

    
    