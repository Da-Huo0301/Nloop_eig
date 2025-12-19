 % THE NUMERICAL RESULT OF THE LADDER-STYLE ELECTROMAGNETICALLY INDUCED
% TRANSPARENCY.
% clear;clc;
tic
G21 = 3*2*pi*1e6;
G23 = 3*2*pi*1e6;
G43 = 3*2*pi*1e6;
g21=(G23+G21)/2;
g31 = 0.001*2*pi*1e3;
g23=(G23+G21)/2;
g24=(G23+G21+G43)/2;
g43=G43/2;
g41=G43/2;
Dc =0*2*pi*1e6;
Dr = 0*2*pi*1e6;
Dpr = (-5.0:0.005:5.0)*2*pi*1e6;
theta_p = 0;
theta_c = 0;
theta_m = 3*pi/2;
Oc =5*2*pi*1e6;
Op =0.01*2*pi*1e6;
Od = 5*2*pi*1e6;
Om =0.5*2*pi*1e6;

N0=4*1e10*1e6;
DD=3.58*10^(-29);
e0=8.854188*10^(-12);
Hbar=1.054572*10^(-34);
yet=N0*DD^2/(2*e0*Hbar);
   
epslon=yet;


Op1 = Op*exp(1i*theta_p); Op2 = Op*exp(-1i*theta_p);
Oc1 = Oc*exp(1i*theta_c); Oc2 = Oc*exp(-1i*theta_c);
Od1 = Od; Od2 = Od;

x = zeros(15,1);
re = zeros('like',Dpr);
im = zeros('like',Dpr);
index = 1;
for Dp = Dpr
    % A and b are obtaied from "ladder_eit_sym.m". 
   A =[               0,           0,                0,             -g21, -Om*sin(theta_m),       0,                  0,                0,                0,             -Dp,  Om*cos(theta_m),       0,               -Oc1,                0,               0;...
               0,           0,                0,  Om*sin(theta_m),             -g23,       0,                  0,                0,                0, Om*cos(theta_m),              -Dc,     Od2,                Op1,                0,               0;...
 Om*sin(theta_m),           0, -Om*sin(theta_m),                0,                0,       0,               -g31,                0,                0,            -Oc2,             -Op1,       0,            Dc - Dp,                0,            -Od1;...
               0,           0,                0,                0,                0,    -g24,                  0,                0,                0,               0,              Od1, Dr - Dc,                  0,             -Oc1,             Op1;...
               0,           0,                0,                0,                0,       0,                  0, -Om*sin(theta_m),             -g41,               0,                0,    -Op1,               -Od2, -Om*cos(theta_m),    Dc - Dp - Dr;...
               0,           0,                0,                0,                0,       0,                  0,             -g43,  Om*sin(theta_m),               0,                0,    -Oc1,                  0,               Dr, Om*cos(theta_m);...
               0,         G21,                0,                0,                0,       0, -2*Om*sin(theta_m),                0,                0,     - Op1 - Op2,                0,       0, -2*Om*cos(theta_m),                0,               0;...
               0, - G21 - G23,                0,                0,                0,       0,                  0,                0,                0,       Op1 + Op2,        Oc1 + Oc2,       0,                  0,                0,               0;...
            -G43,   G23 - G43,             -G43,                0,                0,       0,  2*Om*sin(theta_m),                0,                0,               0,      - Oc1 - Oc2,       0,  2*Om*cos(theta_m),        Od1 + Od2,               0;...
             Op1,        -Op1,                0,               Dp, -Om*cos(theta_m),       0,                Oc1,                0,                0,            -g21, -Om*sin(theta_m),       0,                  0,                0,               0;...
               0,        -Oc1,              Oc1, -Om*cos(theta_m),               Dc,    -Od2,                Op1,                0,                0, Om*sin(theta_m),             -g23,       0,                  0,                0,               0;...
 Om*cos(theta_m),           0, -Om*cos(theta_m),              Oc2,             -Op1,       0,            Dp - Dc,                0,              Od1,               0,                0,       0,               -g31,                0,               0;...
               0,           0,                0,                0,             -Od1, Dc - Dr,                  0,              Oc1,              Op1,               0,                0,    -g24,                  0,                0,               0;...
               0,           0,                0,                0,                0,    -Op1,                Od2, -Om*cos(theta_m),     Dp - Dc + Dr,               0,                0,       0,                  0,  Om*sin(theta_m),            -g41;...
             Od2,         Od2,            2*Od2,                0,                0,    -Oc1,                  0,               Dr, -Om*cos(theta_m),               0,                0,       0,                  0,              g43, Om*sin(theta_m)];
 
    b =[0; 0; 0; 0; 0; 0; 0; 0;-G43; 0; 0; 0; 0; 0; Od2];
    x = A\b;
    ddp(index)=Dp;
    rho_21 = (x(4) + 1i*x(10))*epslon;
    re3(index) = real(rho_21)/Op;
    im3(index) = imag(rho_21)/Op;
    index = index + 1;
    
end
subplot(2,4,1)
plot(Dpr/(2*pi*1e6),im0,'b-',Dpr/(2*pi*1e6),imag(rou310),'r--');
xlabel('\delta p');
ylabel('\chi^{I}_p');
axis([-5 5,-inf,inf]);
set(gca,'FontSize',20,'Fontname', 'Times New Roman');
subplot(2,4,2)
plot(Dpr/(2*pi*1e6),im1,'b-',Dpr/(2*pi*1e6),imag(rou311),'r--');
xlabel('\delta p');
ylabel('\chi^{I}_p2');
axis([-5 5,-inf,inf]);
set(gca,'FontSize',20,'Fontname', 'Times New Roman');
subplot(2,4,3)
plot(Dpr/(2*pi*1e6),im2,'b-',Dpr/(2*pi*1e6),imag(rou312),'r--');
xlabel('\delta p');
ylabel('\chi^{I}_p3');
axis([-5 5,-inf,inf]);
set(gca,'FontSize',20,'Fontname', 'Times New Roman');
subplot(2,4,4)
plot(Dpr/(2*pi*1e6),im3,'b-',Dpr/(2*pi*1e6),imag(rou313),'r--');
xlabel('\delta p');
ylabel('\chi^{I}_p4');
axis([-5 5,-inf,inf]);
set(gca,'FontSize',20,'Fontname', 'Times New Roman');

subplot(2,4,5)
plot(Dpr/(2*pi*1e6),re0,'b-',Dpr/(2*pi*1e6),real(rou310),'r--');
xlabel('\delta p');
ylabel('\chi^{R}_p');
axis([-5 5,-inf,inf]);
set(gca,'FontSize',20,'Fontname', 'Times New Roman');
subplot(2,4,6)
plot(Dpr/(2*pi*1e6),re1,'b-',Dpr/(2*pi*1e6),real(rou311),'r--');
xlabel('\delta p');
ylabel('\chi^{R}_p2');
axis([-5 5,-inf,inf]);
set(gca,'FontSize',20,'Fontname', 'Times New Roman');
subplot(2,4,7)
plot(Dpr/(2*pi*1e6),re2,'b-',Dpr/(2*pi*1e6),real(rou312),'r--');
xlabel('\delta p');
ylabel('\chi^{R}_p3');
axis([-5 5,-inf,inf]);
set(gca,'FontSize',20,'Fontname', 'Times New Roman');
subplot(2,4,8)
plot(Dpr/(2*pi*1e6),re3,'b-',Dpr/(2*pi*1e6),real(rou313),'r--');
xlabel('\delta p');
ylabel('\chi^{R}_p4');
axis([-5 5,-inf,inf]);
set(gca,'FontSize',20,'Fontname', 'Times New Roman');
hold on
toc