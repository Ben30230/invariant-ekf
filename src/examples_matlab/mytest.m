%% For optimization
clc,clear
a=[0 0 1]';
b=[1 0 0]';
R_true = angvec2r(pi/2,[0 1 0])
% Rop = Exp(rand(3,1));
Rop = eye(3);
alpha = 0.1;
for i=1:1000
    u = (Rop*a-b);
    J = -axis2skew(Rop*a);
%     rank(J'*J)
%     phi = J'*J\J'*u;
    phi = -alpha *(-a'*Rop'*axis2skew(b))';
    Rop = Exp(phi)*Rop;
    if norm(phi)<1e-6
        break;
    end
end
Rop

%%
clc,clear
a=[0 0 1]';
b=[1 0 0]';
Rop = eye(3);
alpha = 0.1;
for i=1:1000
    delta = -(Rop*a-b)'*axis2skew(Rop*a);
    if abs(delta-(-a'*Rop'*axis2skew(b))')>1e-6
        error("!")
    end
    phi = -alpha *delta';
    Rop = Exp(phi)*Rop;
    if norm(phi)<1e-6
        break;
    end
end
Rop

%%
clc,clear
a=[0 0 1]';
b=[1 0 0]';
Rop = eye(3);
alpha = 0.1;
for i=1:1000
    u = (Rop*a-b);
    J = -axis2skew(Rop*a);
    phi = lsqr(J'*J,-J'*u);
    Rop = Exp(phi)*Rop;
    if norm(phi)<1e-6
        break;
    end
end
Rop

