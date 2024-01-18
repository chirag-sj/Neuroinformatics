A = rand(4,8);

nr = size(A,1);
nc = size(A,2);
for r=1:nr
    for c=1:nc
        if (A(r,c)>0.5)
            sprintf('Number at %s row and %s column: %f is greater than 0.5',iptnum2ordinal(r),iptnum2ordinal(c),A(r,c))
        else
            sprintf('Number at %s row and %s column: %f is lesser than 0.5',iptnum2ordinal(r),iptnum2ordinal(c),A(r,c))
        end
    end
end