function mat = A_func(A)
%EX1_FUNC Summary of this function goes here
%   Detailed explanation goes here
% Takes a 2D matrix as input and prints where element at a location is greater or lesser than 0.5
nr = size(A,1);
nc = size(A,2);
mat = zeros(nc*nr,3);
iter = 1;
for r=1:nr
    for c=1:nc
        if (A(r,c)>0.5)
            sprintf('Number at %s row and %s column: %f is greater than 0.5',iptnum2ordinal(r),iptnum2ordinal(c),A(r,c))
            mat(iter,1) = r;
            mat(iter,2) = c;
            mat(iter,3) = 1;
            iter=iter+1;
        else
            sprintf('Number at %s row and %s column: %f is lesser than 0.5',iptnum2ordinal(r),iptnum2ordinal(c),A(r,c))
            mat(iter,1) = r;
            mat(iter,2) = c;
            mat(iter,3) = 0;
            iter=iter+1;
        end
    end
end
end

