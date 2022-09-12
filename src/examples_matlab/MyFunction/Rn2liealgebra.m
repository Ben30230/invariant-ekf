function output = Rn2liealgebra(xi)
xi=xi(:);
N=length(xi)/3+2;
output=zeros(N);
for i=3:N
    if i==3
        output(1:3,1:3)=axis2skew(xi(1:3));
    else
        output(1:3,i)=xi(((i-3)*3+1):((i-3)*3+3));
    end
end
end