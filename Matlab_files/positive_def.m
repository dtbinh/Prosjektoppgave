function y = positive_def(A)
if size(A,1) ~= size(A,2)
    y = 0;
else
    y = 1;
    n = size(A,1);
    for i = 1:n
        if det(A(1:i,1:i)) < 0
            y = 0;
            break;
        end
    end
end


