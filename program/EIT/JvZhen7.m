function [Re12,Im12] = JvZhen7(Om, Oc, Oa,  dp, dc, gme, gmm, g23, G21, G31, G32,theta_m,theta_p,theta_c)
%UNTITLED Summary of this function goes he re
%   Detailed explanation goes here
Op3 = Om*exp(1i*theta_p); Op4 = Om*exp(-1i*theta_p);
Oc1 = Oc*exp(1i*theta_c); Oc2 = Oc*exp(-1i*theta_c);


A=zeros(9,9);
B=zeros(9,1);
x=zeros(9,1);



 A =[                 0,               0, Oa*sin(theta_m),                  0,             -g23,  Oa*cos(theta_m),                Op4,               dc;...
 2*Oa*sin(theta_m), Oa*sin(theta_m),               0,               -gmm,                0,             -Oc2,            dc - dp,              Op3;...
                 0,               0,            -gme,                  0, -Oa*sin(theta_m),              -dp,               -Oc1, -Oa*cos(theta_m);...
              -G31,       G21 - G31,               0, -2*Oa*sin(theta_m),                0,      - Op3 - Op4, -2*Oa*cos(theta_m),                0;...
              -G32,     - G21 - G32,               0,                  0,                0,        Op3 + Op4,                  0,      - Oc1 - Oc2;...
               Oc2,           2*Oc2, Oa*cos(theta_m),               -Op4,              -dc, -Oa*sin(theta_m),                  0,             -g23;...
 2*Oa*cos(theta_m), Oa*cos(theta_m),             Oc2,            dp - dc,             -Op3,                0,               -gmm,                0;...
               Op3,            -Op3,              dp,                Oc1, -Oa*cos(theta_m),             -gme,                  0,  Oa*sin(theta_m)];
 
    B=[0; Oa*sin(theta_m); 0; -G31; -G32; Oc2; Oa*cos(theta_m); 0];
x=A\B;
alpha3=x(3)+1i*x(6);
Re12=real(alpha3);
Im12=imag(alpha3);

end