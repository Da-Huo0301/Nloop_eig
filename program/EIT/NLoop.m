%% 修改过的MainA文件
%%
tic
clear
%% Fundamental Physical Constants
    %1. Speed of Light(m/s)
    c=2.997*1e8;
    %2. Permittivity of Vacuum(F/m)
    Epsilon=8.85*1e-12;
    %3. Planck's Constant(J.s)
    Hbar=1.054*1e-34;  

%% Parameters for this three-level Ladder Model
% [Rb:|g>=5^2S_1/2|F=2,mF=2>,|e>=5^2P_3/2|F=3,mF=3>,|r>=|60S_1/2>]

    %1. Dipole Monment((C.m)[5^2S_1/2--->5^2P_3/2:F=2-->F'=3])
    Pge=1.73*1e-29;
    %2. Frequency of Probe(Hz)
    wp=2*pi*384*1e12;
    %3. Decay Rate/Natural Line Width(FWHM)(Hz)
    %Gme=38*1e6;Gmr=5*1e3;Gmd=38*1e6; %-------------------- %% 这里添加了Gmd=Gme
    %4. Length of the sample(m)
    L=5*1e-4;
   
    %7. Detuning of control field(Hz)
    dc=0;dr=0; %-------------------------------- %% 这里添加了dd=-3*dc
    %8. Intensity of control field(Hz^2?)
    Oc=1*2*pi*1e6;
    Od=3*2*pi*1e6;
    Ic=(Oc)^2;
    Ir=(Od)^2;%-------------------------- %% 这里添加了Id=Ic
    %9. Peak of atomic density(m^(-3))
    rho_max=5*1e11/1e-6;
    %10. Half-width of atomic density pulse(m)
    sigma_rho=0.7*1e-3;%------------------------------%此处不知道是否需要改？
    %11. one and two-photon laser linewidths
    %dw1=5.7*2*pi*1e4;dw2=11*2*pi*1e4;dw3=11*2*pi*1e4;%% 这里添加了dw3=11*2*pi*1e4
    %12  Maximal intensity of Probe field
    Om=0.01*2*pi*1e6;
    theta_m=0*pi/2;
    theta_p=0;theta_c=0;
    Oa=2*2*pi*1e6;
%% Generating Parameters

   %1. Coherent decay rate
  gme=-3*2*pi*1e6;gmr=-10*2*pi*1e3;gmm=-1*2*pi*1e3;%% 增加了gmd=(Gmd+dw2)/2，并修改gmr=(Gmr+dw3)/2
  G21=-3*2*pi*1e6;G23=-3*2*pi*1e6;G43=-1*2*pi*1e3;G31=-0.1*2*pi*1e3;G32=-3*2*pi*1e6;
  g23=-1*2*pi*1e6;g24=-1*2*pi*1e6;g43=-1*2*pi*1e3;
   
    %6. Blockade radius(m)
    C6=-1.4*1e11*1e-36;
    Rsa=(C6*gme/((Oc)^2+(Od)^2))^(1/6);
    %5. Average atoms within Rydberg dipole blockade sphere
    n_sa= rho_max*(4/3)*pi*Rsa^3;
   %2. Step/interval
   h=2*Rsa;
   %3. Numbers of interval
   NL=0:h:L;steps2=length(NL);
   %4. Atom-field coupling strength
   eta=Pge*sqrt(wp/(2*Hbar*Epsilon*4/3*pi*Rsa^3));%%-------此处(L/2)^2是否应替换为Rsa^3
   %5. Atomic absorption cross-section
   zeta_0=wp*Pge^2/(Hbar*Epsilon*c*gme);
   
   %6. Half-width of the <sigma_rr> distribution
   w=Ic/gme;%--------------------此处不知道是否需要修改？似乎这一行没啥用
   %7. Central of atomic density pulse
   z0=L/2;
   %
   ab=Pge^2*wp*rho_max/(Hbar*Epsilon*c*gme);
%% Initial values 

    steps=40;
    steps1=500;
    dpmin=-15;dpmax=15;
    Dp=linspace(dpmin,dpmax,steps1)*2*pi*1e6;
    Ip=zeros(steps1,steps2,1);Ip(:,1)=Om^2/eta^2;
    g=zeros(steps1,steps2,1);g(:,1)=1;
    phi=zeros(steps1,steps2,1);phi(:,1)=0;%此处增加了相位的初始条件
    Pol=zeros(steps1,steps2,1);
    A=zeros(steps,steps1);B=A;
%%
%--------------------------------------------------------------------------
% Main loop
%--------------------------------------------------------------------------
    for nn=1:steps
    for n=1:steps1
        dp=Dp(n);
        %d1=dp+dc;
        %d2=dp-dd;
        z=0;
            for k=1:steps2-1  
                z=z+h;       
                rho=rho_max*exp(-(z-z0)^2/(2*sigma_rho^2));
                %rho=1.5*1e7/1e-9;
                K=zeta_0*rho;
               
                  Op0=sqrt(n_sa*Ip(n,k)*eta^2*g(n,k));
      [rou44] = JvZhen5(Op0,Oc,Od,Oa,dp,dc,dr,gme,g23,gmm,g24,gmr,g43,G21,G23,G43,theta_m,theta_p,theta_c);
      sigma=rou44;
%        E1=Ic^2*((dp+dr-dc)^2+gmr^2);
%        E2=2*Ic*((dp+dr-dc)^2+gmr^2)*(gmm*gme-dp*(dp-dc));
%        E3=2*Ic*Ir*(gme*gmr+dp*(dp+dr-dc));
%        E4=(dp^2+gme^2)*(((dp+dr-dc)^2+gmr^2)*((dp-dc)^2+gmm^2)+Ir^2);
%         sigma=Ir*Ic*n_sa*Ip(n,k)*eta^2*g(n,k)/(Ir*Ic*n_sa*Ip(n,k)*eta^2*g(n,k)+E1+E2+E3+E4);
       
% Generating random number to realize with Monte Carlo 
                p=rand(1);
                if p<=sigma
                    sigma=1;
                else
                    sigma=0;
                end
    
                sr=w/8*sigma;%-----------------这行代码并未用到
                [re12,im12] = JvZhen6(Om,Oc,Od,dp,dc,dr,gme,g23,gmm,g24,gmr,g43,G21,G23,G43,theta_p,theta_c);
             [Re12,Im12]= JvZhen7(Om, Oc, Oa,  dp, dc, gme, gmm, g23, G21, G31, G32,theta_m,theta_p,theta_c);
              alpha3=Re12+1i*Im12;
                alpha4=re12+1i*im12;
                %s2=1.4*1e11/10*(1e6)^6/(8*Rsa^6)*sigma;
%                 alpha3=-1i*gme/(1i*dp+gme+Ic/(1i*(dp+dc)+gmm));
%                 alpha4=-1i*gme/(1i*dp+gme+Ic/(1i*(dp-dc)+gmm+Ir/(1i*(dp+dr-dc+sr)+gmr)));
                
                alpha=sigma*alpha3+(1-sigma)*alpha4;%------ %%此处增加了alpha3和alpha4，用来表示总的极化率
                
                Pol(n,k+1)=alpha;
                Ip(n,k+1)=Ip(n,k)*(1-imag(alpha)*K*h);  
                 g(n,k+1)=g(n,k)*(1-K*h*imag(alpha-alpha4));
                %g(n,k+1)=g(n,k)*(1-K*h*imag(alpha-(-1i*gme/((1i*dp+gme)+Id/(1i*d2+gmm)+Ic/(1i*(d1-0)+gmr)))));
                phi(n,k+1)=phi(n,k)+K*h*real(alpha)/2;%此处增加了相位的相位满足的传播方程
                
            end
        
        Pn(n)=mean(Pol(n,:));
        An(n)=Ip(n,end)/Ip(n,1);
        Bn(n)=g(n,end)/g(n,1);
        Cn(n)=phi(n,end);
    end
    PP(nn,:)=Pn;
    A(nn,:)=An;
    B(nn,:)=Bn;
    C(nn,:)=Cn;
    end
%%   
    P=mean(PP);
    T=mean(A);
    G=mean(B);
    Phase=mean(C);%此处增加了相位的部分
%% 透射率VS探测失谐
figure(1)

subplot(2,1,1)
hold on
plot(Dp/1e6/(2*pi),T,'b','LineWidth',2)
ylabel('Transmissian I_p(L)/I_p(0)')
axis([dpmin dpmax 0 2])
% title([ ' \Omega_p=',num2str(Om) ]);
hold on
%% 关联函数VS探测失谐
subplot(2,1,2)
hold on
plot(Dp/1e6/(2*pi),G,'b')
ylabel('g(2)_p(L)')
xlabel('{\Delta_p}/2{\pi}[MHz]','LineWidth',2)
axis([dpmin dpmax 0 1.6])

%% 相位VS探测失谐
%此处增加了相位的画图部分
% subplot(3,1,3)
% hold on
% plot(Dp/1e6/(2*pi),Phase,'r')
% ylabel('\phi(L)/\pi','LineWidth',2)
% xlabel('{\Delta_p}/2{\pi}[MHz]','LineWidth',2)

%%
figure(2)
hold on
plot(Dp/1e6/(2*pi),imag(P),'b--',Dp/1e6/(2*pi),real(P),'r')
toc

% Reproduce the paper of PRL: Electromagnetically Induced Transparency with Rydberg Atoms