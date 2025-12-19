function[IP]=diffraction_1D(R,M,T1,index,k)
%BIpX1=zeros(500,3);
%R=4;M=10;
for p=1:501
    sQ=-1.0001+0.004*p;
    Ep1=0;Epa=0;
for o=1:2001
       x=-0.50005+o/2000;
       %TT=T(2,o);%hs画连续介质光栅
       TT=T1(o,index);%hs画分段图时
       Ep1=Ep1+TT*exp(-i*2*pi*x*sQ*R)/2000;
       %Epa=Epa+1*exp(-i*2*pi*x*sQ*R)/2000;
end

    SSN=sin(M*pi*sQ*R)^2/(M*sin(pi*sQ*R))^2;
     Ep2=abs(Ep1)^2;
    % Epb=abs(Epa)^2;
     %IpX=Ep2^2*SSN;
     IpX=Ep2*SSN;
     %IpXa=Epb^2*SSN;    
    % BIpX1(p,1)=sQ;
     BIpX1(p)=IpX;
    % BIpX1(p,3)=IpX/IpXa;     
end
IP=BIpX1(round(k*250/R)+251);%k取1到R的正整数
% CIp=0;
% for r=1:2001
% CIp=CIp+BIpX1(r,2);
% end
% EIpX=BIpX1(1,2)+BIpX1(500,2)+BIpX1(375,2)+BIpX1(250,2)+BIpX1(125,2)+...
%     BIpX1(62,2)+BIpX1(63,2)+BIpX1(187,2)+BIpX1(188,2)+...
%     BIpX1(437,2)+BIpX1(438,2)+BIpX1(312,2)+BIpX1(313,2);
% for u=1:2001
%      DIpX4(u,1)=BIpX1(u,1);
%      DIpX4(u,2)=BIpX1(u,2);%/EIpX; %EIpXCIp
% end
%OpticalDeepth=8*pi*N00*DD^2*L/(e0*Hbar*2*gamma14*lmdap)
% (BIpX1(312,2)-BIpX1(187,2))/(BIpX1(312,2)+BIpX1(187,2))
% (BIpX1(375,2)-BIpX1(125,2))/(BIpX1(375,2)+BIpX1(125,2))
% (BIpX1(438,2)-BIpX1(62,2))/(BIpX1(438,2)+BIpX1(62,2))
% figure
% plot(BIpX1(:,1),BIpX1(:,2),'r-')%衍射图样
% figure
% plot(BIpX1(:,1),BIpX1(:,3),'b--')
%figure
%plot(DIpX4(:,1),DIpX4(:,2),'b-')
%title('diffraction');

