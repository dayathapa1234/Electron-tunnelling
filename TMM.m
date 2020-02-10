clear all
close all
hBar = 1.054571800e-34;
N = 1000; 
L = [0.5,1,2].*1e-9; %Barrier Length nm
%L=0.5e-9;
%Vo=[0;4;6];
Vn = linspace(0,2, N);
V =4.*((((Vn-1).^2)-1)-(((Vn-1).^4)-1)).*1e-19;%2*(exp(-5*Vn)).*1e-19;%((-Vn/2)+2).*1e-19;% ((-Vn/2)+2).*1e-19;%(-(1/8).*((Vn-2).^2).*((Vn-2).^2)+2).*1e-19;%
%V=2e-19;
%Particle Properties ------------------------------------------------------
E = linspace(0,5, N)*1e-19; %Incident particle energies eV
m = 9.10938356e-31; 
figure;
hold on
AP1 = 1;
AN3 = 0;
T = zeros(1,N);
for r = 1:3;
    for n = 1:N;
        k0 = sqrt(2*m*E(n))/hBar;
        k1 = sqrt(2*m*(E(n) - V(n)))/hBar;
        k2 = k0;
        
        D12 = 0.5 * [(1+k1/k0) (1-k1/k0); (1-k1/k0) (1+k1/k0)];
        P2 = [exp(-i.*k1.*L(r)) 0; 0 exp(i.*k1.*L(r))];
        D21 = 0.5*[(1+k2/k1) (1-k2/k1); (1-k2/k1) (1+k2/k1)];

        Q = D12 * P2 * D21; %Transfer Matrix

        AP3 = 1/Q(1,1); 
        AP1 = Q(2,1)*AP3;

        T(n) = abs(AP3)^2; %Transmission probability
    end

    plot((E./1e-19), (T),'LineWidth',3)
   
    ylabel('T')
    xlabel('Electron Energy (eV)')
    xlim([0 2]);
    ylim([0 1]);
   
    [FX] = gradient((E./1e-19),T);
end 
hold on
 title('Exponential Barrier of 2eV')
 plot(n,(V/1e-19));
 hold off
    legend('0.5nm','1nm','2nm')
    set(gca,'fontsize',25)
grid on
figure;
n=[1:1:N];
plot(n,(V/1e-19));


