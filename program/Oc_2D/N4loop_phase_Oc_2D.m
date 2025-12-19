% THE NUMERICAL RESULT OF THE LADDER-STYLE ELECTROMAGNETICALLY INDUCED
% TRANSPARENCY.
clear;clc;
tic
% danwei=2*pi*10
G21 = 3*2*pi*1e6;
G23 = 3*2*pi*1e6;
G43 = 3*2*pi*1e6;
g21=(G23+G21)/2;
g31 = 0.0011*2*pi*1e3;
g23=(G23+G21)/2;
g24=(G23+G21+G43)/2;
g43=G43/2;
g41=G43/2;
Dc = 0.0;
Dr0 = 6*2*pi*1e6;
%Dpr = (-10.0:0.01:10.0)*2*pi*1e6;
theta_p = 0;
theta_c = 0;
theta_m =0*pi/2;
%Oc = 3*2*pi*1e6;
Op =0.01*2*pi*1e6;
Od = 5*2*pi*1e6;
Om = 0.4*2*pi*1e6;

Op1 = Op*exp(1i*theta_p); Op2 = Op*exp(-1i*theta_p);

Od1 = Od; Od2 = Od;



y = zeros(15,1);
%re = zeros('like',Dpr);
%im = zeros('like',Dpr);
lambda=780*10^(-9);
% Dc0=0*2*pi*1e6;
a=2*lambda;
N0=4.2*1e10*1e6;
%N0=5.5*1e11*1e6;% 原子密度;%%%%%%%%%%%%%%%%设多少合适？
h=0.1*a;%%%%%%%%%%%%%%半高宽设多少合适？

DD=3.58*10^(-29);
e0=8.854188*10^(-12);
Hbar=1.054572*10^(-34);
yet=N0*DD^2/(2*e0*Hbar);
   
epslon=yet;
Dp=0*2*pi*1e6;

R=6;M=10;lamdp=780*1e-9;
 L=0.6*1e-5;%这个是介质衍射方向长度，自己设定，越长增益越大


cc=0;
for x=-0.5*a:a/1000:0.5*a
    cc=cc+1;
    Dr=Dr0*sin(pi*x/a);
    xx(cc)=x;
    cc1=0;
for Oc=0*2*pi*1e6:10*2*pi*1e6/1000:10*2*pi*1e6 
           cc1=cc1+1;
           OOCC(cc1)=Oc/(2*pi*1e6);
       Oc1= Oc*exp(1i*theta_c); Oc2= Oc*exp(-1i*theta_c); 
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
    y = A\b;
    rho_21 = (y(4) + 1i*y(10))*epslon;
    re(cc,cc1) = real(rho_21)/Op;
    im(cc,cc1) = imag(rho_21)/Op; 
    T1(cc,cc1)=exp(-2*pi*L*im(cc,cc1)/lamdp)*exp(i*2*pi*L*re(cc,cc1)/lamdp);
end
end

 %save matfile1.mat xx re im ddp;
% save re;
% save im;

% figure(1)
% plot(xx,re,'r-',xx,im,'b-');
% hold on
Ip=zeros(3,1001);
for index=1:1001   
for cc=1:1001
T1(cc,index)=exp(-2*pi*L*im(cc,index)/lamdp)*exp(1i*2*pi*L*re(cc,index)/lamdp);
end
for k=1:3
Ip(k,index)=diffraction_1D_2000_4(R,M,T1,index,k);
end
end
figure(1)
plot(OOCC,Ip(1,:),'r-',OOCC,Ip(2,:),'b-',OOCC,Ip(3,:),'k-','LineWidth',2);
%红蓝黑代表123级
set(gca,'FontSize',20,'Fontname', 'Times New Roman');
xlabel('\Omega_c','FontSize',30,'Fontname', 'Times New Roman');
ylabel('I_p(\theta)','FontSize',30,'Fontname', 'Times New Roman');
  

 toc


