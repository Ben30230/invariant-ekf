function [R,v,p] = generateIMUdata(w,bg,a,ba,R0,v0,p0,deltat,gravity)
%
n=length(w);
R = zeros(3,3,n+1);
v = zeros(3,n+1);
p = zeros(3,n+1);

R(:,:,1) = R0;
v(:,1) = v0;
p(:,1) = p0;

for i=1:n
    R(:,:,i+1) = R(:,:,i) *Exp((w(:,i)-bg)*deltat);
    v(:,i+1) = v(:,i) + gravity*deltat + R(:,:,i)*(a(:,i)-ba)*deltat;
    p(:,i+1) = p(:,i) + v(:,i)*deltat + 0.5 * (gravity+R(:,:,i)*(a(:,i)-ba))*deltat^2;
end

end