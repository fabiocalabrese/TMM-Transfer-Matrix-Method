clc
clear all
n0=1; 
ns=1.6; %refractive index of air and glass substrate 
N = 14; % number of layers

lambda = 300:5:1200;
ks = [0.000132, 0.000125, 0.000118, 0.000112, 0.000107, 0.000102];
kto2 = [ 0.880
    0.7740
    0.6621
    0.5875
    0.4755
    0.3672
    0.3006
    0.2139
    0.1452
    0.1092
    0.0681
    0.0481
    0.0269
    0.0176
    0.0114
    0.0060
    0.0040
    0.0026
    0.0017
    0.0009
    0.0006
    0.0004
    0.0003
    0.0002
    0.0001
    0.0001
    0.0001];
% Inizializza tutto il vettore con zeri
ve = zeros(size(lambda));
vt = zeros(size(lambda));
% Inserisci i dati noti all'inizio
ve(1:length(ks)) = ks;
vt(1:length(kto2)) = kto2;
 
M1=[0, 0; 0, 0];
M2=[0, 0; 0, 0];
Md3=[0, 0; 0, 0];
ii=0; 
for j=300:5:1200
   ii=ii+1;
   l(ii)=j;
   ll(ii)=l(ii)/1000; %micrometri
    % calcolare l'indice di rifrazione per ogni lunghezza d'onda
%     n1(ii) = sqrt(((A(1) * l(ii)^2) / (l(ii)^2 - B(1)^2)) ...
%          + ((C(1) * l(ii)^2) / (l(ii)^2 - D(1)^2)) ...
%          + ((E(1) * l(ii)^2) / (l(ii)^2 - F(1)^2)) + 1);
%      n2(ii) = sqrt(((A(2) * l(ii)^2) / (l(ii)^2 - B(2)^2)) + ...
%          ((C(2) * l(ii)^2) / (l(ii)^2 - D(2)^2)) + ...
%          ((E(2) * l(ii)^2) / (l(ii)^2 - F(2)^2)) + 1);
%n1(ii) = 2.4;
%n2(ii) = 1.8;
    % TiO2
    n1(ii)=sqrt(4.99+...
        1/(96.6*(ll(ii)^1.1))+...
        1/(4.6*(ll(ii)^1.95))) +1i*vt(ii);
 %   SiO2
    n2(ii)=sqrt(0.82+...
         (0.7161663*ll(ii)^2/(ll(ii)^2-0.0684043^2))+...
         (0.4079426*ll(ii)^2/(ll(ii)^2-0.1162414^2))+...
         (0.8974794*ll(ii)^2/(ll(ii)^2-9.896161^2))+...
         1/(700*(ll(ii)^6.2))) + 1i*ve(ii);
       % ZrO2         
n3(ii)=sqrt(1+...
    (1.347091*ll(ii)^2/(ll(ii)^2-0.062543^2))+...
     (2.117788*ll(ii)^2/(ll(ii)^2-0.166739^2))+...
     (9.452943*ll(ii)^2/(ll(ii)^2-24.32057^2)));

    
    d1(ii) = 100 ; % Spessore (nm) dello strato 1 
    d2(ii) = 100; % Spessore (nm) dello strato 2 
    dd(ii) = 145; 
    d3 = 100;

    % Calcolo delle matrici per ogni lunghezza d'onda
    Sd(ii) = cos(2 * pi * n1(ii) * dd(ii) / l(ii));   
    Pd(ii) = -(1i * sin(2 * pi * n1(ii) * dd(ii) / l(ii))) / n1(ii); 
    Qd(ii) = -1i * n1(ii) * sin(2 * pi * n1(ii) * dd(ii) / l(ii));   
    Rd(ii) = cos(2 * pi * n1(ii) * dd(ii) / l(ii));
    Md = [Sd(ii), Pd(ii); Qd(ii), Rd(ii)]; % Matrice difetto

    % difetto con n3
    Sd3(ii) = cos(2 * pi * n3(ii) * d3 / l(ii));   
    Pd3(ii) = -(1i * sin(2 * pi * n3(ii) * d3 / l(ii))) / n3(ii); 
    Qd3(ii) = -1i * n3(ii) * sin(2 * pi * n3(ii) * d3 / l(ii));   
    Rd3(ii) = cos(2 * pi * n3(ii) * d3 / l(ii));
    Md3 = [Sd3(ii), Pd3(ii); Qd3(ii), Rd3(ii)]; % Matrice difetto

    S1(ii)=cos(2*pi*n1(ii)*d1(ii)/l(ii));   
    P1(ii)=-(1i*sin(2*pi*n1(ii)*d1(ii)/l(ii)))/n1(ii); 
    Q1(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d1(ii)/l(ii));   
    R1(ii)=cos(2*pi*n1(ii)*d1(ii)/l(ii)); 
    M1=[S1(ii), P1(ii); Q1(ii), R1(ii)]; %matrix 1 for layer 1 

    S2(ii)=cos(2*pi*n2(ii)*d2(ii)/l(ii));  
    P2(ii)=-(1i*sin(2*pi*n2(ii)*d2(ii)/l(ii)))/n2(ii);
    Q2(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d2(ii)/l(ii));
    R2(ii)=cos(2*pi*n2(ii)*d2(ii)/l(ii)); 
    M2 = [S2(ii), P2(ii); Q2(ii), R2(ii)]; %matrix 2 for layer 2 
    
    Minitial = (M1 * M2)^(N/2); % Moltiplicazione elemento per elemento delle matrici
    Mintermediate = Minitial * Md; % Inserimento dello strato difetto
    Mintermediate3 = Minitial * Md3; % inserimento difetto n3
    M = Mintermediate * Minitial; % Completamento della matrice
    M3 = Mintermediate3*Minitial;
    Mclean = (M1 * M2)^N;



    t(ii)=(2*ns)/(((M(1)+(M(3)*n0))*ns)+(M(2)+(M(4)*n0))); %Luigino SPIE
   T(ii)=(n0/ns)*((abs(t(ii)))^2);
  
   t2(ii)=(2*ns)/(((M3(1)+M3(3)*n0)*ns)+(M3(2)+M3(4)*n0)); %Luigino SPIE   
   T3(ii)=(n0/ns)*((abs(t2(ii)))^2);

      t_clean(ii)=(2*ns)/(((Mclean(1)+Mclean(3)*n0)*ns)+(Mclean(2)+Mclean(4)*n0)); %Luigino SPIE   
   Tclean(ii)=(n0/ns)*((abs(t_clean(ii)))^2);



    
end 

figure
hold on 

subplot(3,1,1)
plot(l,Tclean,'k', LineWidth = 1.5)
title(" SiO2/TiO2 Bragg mirrors");
xlabel("Wavelength (nm)");
ylabel("Trasmittance");

subplot(3,1,2)
plot(l,T,'k',LineWidth = 1.5);
title(" SiO2/TiO2 Bragg mirrors with defect d");
xlabel("Wavelength (nm)");
ylabel("Trasmittance");

subplot(3,1,3)
plot(l,T3,'k',LineWidth = 1.5)
title("SiO2/TiO2 Bragg mirrors with defect n3");
xlabel("Wavelength (nm)");
ylabel("Trasmittance");

