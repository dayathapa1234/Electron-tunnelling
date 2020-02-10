clear all
close all
hBar = 1.054571800e-34;
N = 1000; 
%L = [5,7,15].*1e-9; %Barrier Length nm
L1=0.1e-9;%Barrier
L2=0.1e-9;%Barrier
L3=0.1e-9;%Well
V1 = [10e-19];
V2 = [10]*1e-19;
Maxenergy = 5;

%Particle Properties ------------------------------------------------------
E = linspace(0,Maxenergy, N)*1e-19; %Incident particle energies eV
m = 9.10938356e-31; 

AP1 = 1;
AN3 = 0;
T1 = zeros(1,N);
R1= zeros(1,N);
AP11 = 1;
        AN33 = 0;
        T2 = zeros(1,N);
        R2= zeros(1,N);

    for n = 1:N;
        k1 = sqrt(2*m*E(n))/hBar;
        k2 = sqrt(2*m*(E(n) - V1))/hBar;
        k3 = k1;
        
        D2 = 0.5 * [(1+k2/k1) (1-k2/k1); (1-k2/k1) (1+k2/k1)];
        P2 = [exp(-1i.*k2.*L1) 0; 0 exp(1i.*k2.*L1)];
        D3 = 0.5*[(1+k3/k2) (1-k3/k2); (1-k3/k2) (1+k3/k2)];

        Q2 = D2 * P2 * D3; %Transfer Matrix

        AP3 = 1/Q2(1,1); 
        AP1 = Q2(2,1)/Q2(1,1);

        T1(n) = abs(AP3)^2; %Transmission probability
        R1(n) = abs(AP1)^2;
       
        k4 = sqrt(2*m*T1(n).*E(n))/hBar;
        k5 = sqrt(2*m*T1(n).*(E(n) - V2))/hBar;
        k6 = k4;
        
        D22 = 0.5 * [(1+k5/k4) (1-k5/k4); (1-k5/k4) (1+k5/k4)];
        P22 = [exp(-1i.*k5.*L2) 0; 0 exp(1i.*k5.*L2)];
        D33 = 0.5*[(1+k6/k5) (1-k6/k5); (1-k6/k5) (1+k6/k5)];

        Q22 = D22 * P22 * D33; %Transfer Matrix

        AP33 = 1/Q22(1,1); 
        AP11 = Q22(2,1)/Q22(1,1);

        T2(n) = abs(AP33)^2; %Transmission probability
        R2(n) = abs(AP11)^2;
       end;
  Alp=  T1.*T2;
  Bet= (1-sqrt(R1.*R2)).^2;
  Gam= 4*sqrt(R1.*R2);
  Zta= k4*L3 +(angle(Q2(1,2))+angle(Q22(2,1))-angle(Q2(1,1))-angle(Q22(1,1)))/2;
  TRTD= Alp./(Bet+Gam.*cos(Zta).^2);
  plot(E,TRTD);
 
  

 for i = 1:N;
     Eres(i) = ((hBar^2)*(i^2)*(pi^2))/(2*m*L3^2);
   
     xa(i)=2*pi*i;
     Xaa=sin((1/2)*xa).^2;
 end
   %RowL= ((hBar^2)*k2/m)*((T1.^2)/trapz(T1.^2))
   %dEdZta= gradient(E)/gradient(ZtaRes);
   %Row= ((dEdZta*T1)/2).*((dEdZta*T2)/2);
   %DRow = (((dEdZta*T1)/2)+((dEdZta*T2)/2)).^2;
   %De = 4+(E-Eres);   
   %Tres = (Row)./(DRow./De);
  Al=  T1.*T2;
  Be= (1-sqrt(R1.*R2)).^2;
  Ga= 4*sqrt(R1.*R2);
  
  TRT= Al./(Be+Ga.*Xaa);
   
   

   figure
   plot(E,TRT);
   
