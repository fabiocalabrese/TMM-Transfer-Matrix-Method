%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The base code was provided by Francesco Scotognella "from 01HGBUU -
%Materials and Processes for quantum sensing, metrology and qubit devices"
%course (Politecnico di Torino)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fabio Calabrese 
% Modification of the original code are related to the angular dependence,
% starting from the formulas up to the plot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Transfer Matrix Method
%system: air-multilayer-glass
% angular dependence
clear all;clc;
n0=1; %air
Em=1; %E field in glass
ns=1.6; %refractive index glass
Hm=ns; %H field in glass
%n1=2.4;
d1=100; %refractive index and thickness (nm) layer 1
%n2=1.45;
d2=100; %refractive index and thickness (nm) layer 2
C=0.9;%coefficiente perdite
lambda = 300:5:1200;

% Immaginary refractive index of SiO2
ks = [0.000132, 0.000125, 0.000118, 0.000112, 0.000107, 0.000102];
% Immaginary refractive index of GaAs
k = [ ...
    2, 1.9522, 1.9200, 1.9099, 1.9052, 1.9082, 1.9183, 1.9307, 1.9489, 1.9679, ...
    1.9924, 2.0178, 2.0478, 2.0783, 2.1174, 2.1565, 2.2079, 2.2613, 2.2594, 2.2021, ...
    2.1445, 2.0706, 1.9968, 1.9208, 1.8409, 1.7610, 1.6017, 1.3632, 1.1246, 0.9514, ...
    0.8615, 0.7715, 0.6889, 0.6446, 0.6002, 0.5559, 0.5230, 0.4972, 0.4714, 0.4456, ...
    0.4271, 0.4102, 0.3933, 0.3764, 0.3633, 0.3520, 0.3406, 0.3293, 0.3184, 0.3094, ...
    0.3004, 0.2914, 0.2825, 0.2741, 0.2674, 0.2607, 0.2540, 0.2473, 0.2405, 0.2355, ...
    0.2306, 0.2256, 0.2207, 0.2158, 0.2109, 0.2060, 0.2011, 0.1962, 0.1913, 0.1864, ...
    0.1815, 0.1771, 0.1732, 0.1694, 0.1655, 0.1616, 0.1578, 0.1539, 0.1498, 0.1450, ...
    0.1402, 0.1354, 0.1306, 0.1258, 0.1210, 0.1161, 0.1117, 0.1094, 0.1071, 0.1048, ...
    0.1025, 0.1002, 0.0979, 0.0956, 0.0933, 0.0910, 0.0899, 0.0889, 0.0878, 0.0867, ...
    0.0857, 0.0846, 0.0835, 0.0825, 0.0814, 0.0803 ];
% Immaginary refractive index of TiO2
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
vk = zeros(size(lambda));
vk(1:length(k)) = k;
%prova per ordinato
ddd=[83,96,81,89,81,88,86,91,87,86,...
104,117,118,89,100,96,101,91,83,98,87,82,119,98,119,111,81,108,48];
d=ddd*0.99;

ii=0;
Tt=0;


% angle values
thd=[0 3 6 10 13 16 20 23 27 30 33 37 40 ...
 43 47 50 53 56.7 60 63 67 70 73 77 80 83 87 90];

% conversion grad-radiant
th=(thd*pi)./180;

for j = 1:28
    M1=[0, 0; 0, 0];M2=[0, 0; 0, 0];
    M3=[0, 0; 0, 0];M4=[0, 0; 0, 0];
    M5=[0, 0; 0, 0];M6=[0, 0; 0, 0];
    M7=[0, 0; 0, 0];M8=[0, 0; 0, 0];
    M9=[0, 0; 0, 0];M10=[0, 0; 0, 0];
    M11=[0, 0; 0, 0];M12=[0, 0; 0, 0];
    M13=[0, 0; 0, 0];M14=[0, 0; 0, 0];
    M15=[0, 0; 0, 0];M16=[0, 0; 0, 0];
    M17=[0, 0; 0, 0];M18=[0, 0; 0, 0];
    M19=[0, 0; 0, 0];M20=[0, 0; 0, 0];
    M21=[0, 0; 0, 0];M22=[0, 0; 0, 0];
    M23=[0, 0; 0, 0];M24=[0, 0; 0, 0];
    M25=[0, 0; 0, 0];M26=[0, 0; 0, 0];
    M27=[0, 0; 0, 0];M28=[0, 0; 0, 0];
    M41=[0, 0; 0, 0];M42=[0, 0; 0, 0];
    ii = 0;
    for f=300:5:1200
       ii=ii+1;
       l(ii)=f;
       ll(ii)=l(ii)/1000; %micrometri
       %TiO2
       n1(ii)=sqrt(4.99+...
            1/(96.6*(ll(ii)^1.1))+...
            1/(4.6*(ll(ii)^1.95))) +1i*vt(ii);
       % SiO2
       n2(ii)=sqrt(0.82+...
           (0.7161663*ll(ii)^2/(ll(ii)^2-0.0684043^2))+...
           (0.4079426*ll(ii)^2/(ll(ii)^2-0.1162414^2))+...
           (0.8974794*ll(ii)^2/(ll(ii)^2-9.896161^2))+...
           1/(700*(ll(ii)^6.2))) +1i*ve(ii);
       %% for GaAs
       %n2(ii) = sqrt(11.10+0.743*ll(ii)^2/(ll(ii)^2-0.436^2)) ;%+1i*vk(ii);
        %n2(ii) = n2(ii) + 1i*vk(ii);
       %% for Al2O3
%          n1(ii)=sqrt(1+...
%             (1.023798*ll(ii)^2/(ll(ii)^2-0.0614482^2))+...
%             (1.058264*ll(ii)^2/(ll(ii)^2-0.1106997^2))+...
%             (5.280792*ll(ii)^2/(ll(ii)^2-17.92656^2)));
        %% for ZrO2
%          n2(ii)=sqrt(1+...
%        (1.347091*ll(ii)^2/(ll(ii)^2-0.062543^2))+...
%        (2.117788*ll(ii)^2/(ll(ii)^2-0.166739^2))+...
%        (9.452943*ll(ii)^2/(ll(ii)^2-24.32057^2)));

        

    
       %TE - s-polarization
       coth1(j,ii)=sqrt(1-((n0*n0*sin(th(j))*sin(th(j)))/(n1(ii)*n1(ii))));
       p1(j,ii)=n1(ii)*coth1(j,ii);
       coth2(j,ii)=sqrt(1-((n0*n0*sin(th(j))*sin(th(j)))/(n2(ii)*n2(ii))));
       p2(j,ii)=n2(ii)*coth2(j,ii);
    
       %matrix 1 for layer 1
       S1(ii)=cos(2*pi*n1(ii)*d(1)*coth1(j,ii)/l(ii));
       P1(ii)=-(1i*sin(2*pi*n1(ii)*d(1)*coth1(j,ii)/l(ii)))/n1(ii);
       Q1(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d(1)*coth1(j,ii)/l(ii));
       R1(ii)=cos(2*pi*n1(ii)*d(1)*coth1(j,ii)/l(ii));
       M1=[S1(ii), P1(ii); Q1(ii), R1(ii)];
       %matrix 2 for layer 2
       S2(ii)=cos(2*pi*n2(ii)*d(2)*coth2(j,ii)/l(ii));
       P2(ii)=-(1i*sin(2*pi*n2(ii)*d(2)*coth2(j,ii)/l(ii)))/n2(ii);
       Q2(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(2)*coth2(j,ii)/l(ii));
       R2(ii)=cos(2*pi*n2(ii)*d(2)*coth2(j,ii)/l(ii));
       M2=[S2(ii), P2(ii); Q2(ii), R2(ii)];
       %matrix 3 for layer 3
       S3(ii)=cos(2*pi*n1(ii)*d(3)*coth1(j,ii)/l(ii));
       P3(ii)=-(1i*sin(2*pi*n1(ii)*d(3)*coth1(j,ii)/l(ii)))/n1(ii);
       Q3(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d(3)*coth1(j,ii)/l(ii));
       R3(ii)=cos(2*pi*n1(ii)*d(3)*coth1(j,ii)/l(ii));
       M3=[S3(ii), P3(ii); Q3(ii), R3(ii)];
       %matrix 4 for layer 4
       S4(ii)=cos(2*pi*n2(ii)*d(4)*coth2(j,ii)/l(ii));
       P4(ii)=-(1i*sin(2*pi*n2(ii)*d(4)*coth2(j,ii)/l(ii)))/n2(ii);
       Q4(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(4)*coth2(j,ii)/l(ii));
       R4(ii)=cos(2*pi*n2(ii)*d(4)*coth2(j,ii)/l(ii));
       M4=[S4(ii), P4(ii); Q4(ii), R4(ii)];
       %matrix 5 for layer 5
       S5(ii)=cos(2*pi*n1(ii)*d(5)*coth1(j,ii)/l(ii));
       P5(ii)=-(1i*sin(2*pi*n1(ii)*d(5)*coth1(j,ii)/l(ii)))/n1(ii);
       Q5(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d(5)*coth1(j,ii)/l(ii));
       R5(ii)=cos(2*pi*n1(ii)*d(5)*coth1(j,ii)/l(ii));
       M5=[S5(ii), P5(ii); Q5(ii), R5(ii)];
       %matrix 6 for layer 6
       S6(ii)=cos(2*pi*n2(ii)*d(6)*coth2(j,ii)/l(ii));
       P6(ii)=-(1i*sin(2*pi*n2(ii)*d(6)*coth2(j,ii)/l(ii)))/n2(ii);
       Q6(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(6)*coth2(j,ii)/l(ii));
       R6(ii)=cos(2*pi*n2(ii)*d(6)*coth2(j,ii)/l(ii));
       M6=[S6(ii), P6(ii); Q6(ii), R6(ii)];
       % matrix 7 for layer 7
       S7(ii) = cos(2*pi*n1(ii)*d(7)*coth1(j,ii)/l(ii));
       P7(ii) = -(1i*sin(2*pi*n1(ii)*d(7)*coth1(j,ii)/l(ii))) / n1(ii);
       Q7(ii) = -1i*n1(ii)*sin(2*pi*n1(ii)*d(7)*coth1(j,ii)/l(ii));
       R7(ii) = cos(2*pi*n1(ii)*d(7)*coth1(j,ii)/l(ii));
       M7 = [S7(ii), P7(ii); Q7(ii), R7(ii)];
        
       % matrix 8 for layer 8
       S8(ii) = cos(2*pi*n2(ii)*d(8)*coth2(j,ii)/l(ii));
       P8(ii) = -(1i*sin(2*pi*n2(ii)*d(8)*coth2(j,ii)/l(ii))) / n2(ii);
       Q8(ii) = -1i*n2(ii)*sin(2*pi*n2(ii)*d(8)*coth2(j,ii)/l(ii));
       R8(ii) = cos(2*pi*n2(ii)*d(8)*coth2(j,ii)/l(ii));
       M8 = [S8(ii), P8(ii); Q8(ii), R8(ii)];

       % matrix 9 for layer 9
       S9(ii) = cos(2*pi*n1(ii)*d(9)*coth1(j,ii)/l(ii));
       P9(ii) = -(1i*sin(2*pi*n1(ii)*d(9)*coth1(j,ii)/l(ii))) / n1(ii);
       Q9(ii) = -1i*n1(ii)*sin(2*pi*n1(ii)*d(9)*coth1(j,ii)/l(ii));
       R9(ii) = cos(2*pi*n1(ii)*d(9)*coth1(j,ii)/l(ii));
       M9 = [S9(ii), P9(ii); Q9(ii), R9(ii)];

       % matrix 10 for layer 10
       S10(ii)=cos(2*pi*n2(ii)*d(10)*coth2(j,ii)/l(ii));
       P10(ii)=-(1i*sin(2*pi*n2(ii)*d(10)*coth2(j,ii)/l(ii)))/n2(ii);
       Q10(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(10)*coth2(j,ii)/l(ii));
       R10(ii)=cos(2*pi*n2(ii)*d(10)*coth2(j,ii)/l(ii));
       M10=[S10(ii), P10(ii); Q10(ii), R10(ii)];

       % matrix 11 for layer 11
       S11(ii)=cos(2*pi*n1(ii)*d(11)*coth1(j,ii)/l(ii));
       P11(ii)=-(1i*sin(2*pi*n1(ii)*d(11)*coth1(j,ii)/l(ii)))/n1(ii);
       Q11(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d(11)*coth1(j,ii)/l(ii));
       R11(ii)=cos(2*pi*n1(ii)*d(11)*coth1(j,ii)/l(ii));
       M11=[S11(ii), P11(ii); Q11(ii), R11(ii)];

       % matrix 12 for layer 12
       S12(ii)=cos(2*pi*n2(ii)*d(12)*coth2(j,ii)/l(ii));
       P12(ii)=-(1i*sin(2*pi*n2(ii)*d(12)*coth2(j,ii)/l(ii)))/n2(ii);
       Q12(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(12)*coth2(j,ii)/l(ii));
       R12(ii)=cos(2*pi*n2(ii)*d(12)*coth2(j,ii)/l(ii));
       M12=[S12(ii), P12(ii); Q12(ii), R12(ii)];

       % matrix 13 for layer 13
       S13(ii)=cos(2*pi*n1(ii)*d(13)*coth1(j,ii)/l(ii));
       P13(ii)=-(1i*sin(2*pi*n1(ii)*d(13)*coth1(j,ii)/l(ii)))/n1(ii);
       Q13(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d(13)*coth1(j,ii)/l(ii));
       R13(ii)=cos(2*pi*n1(ii)*d(13)*coth1(j,ii)/l(ii));
       M13=[S13(ii), P13(ii); Q13(ii), R13(ii)];

       % matrix 14 for layer 14
       S14(ii)=cos(2*pi*n2(ii)*d(14)*coth2(j,ii)/l(ii));
       P14(ii)=-(1i*sin(2*pi*n2(ii)*d(14)*coth2(j,ii)/l(ii)))/n2(ii);
       Q14(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(14)*coth2(j,ii)/l(ii));
       R14(ii)=cos(2*pi*n2(ii)*d(14)*coth2(j,ii)/l(ii));
       M14=[S14(ii), P14(ii); Q14(ii), R14(ii)];

       % matrix 15 for layer 15
       S15(ii)=cos(2*pi*n1(ii)*d(15)*coth1(j,ii)/l(ii));
       P15(ii)=-(1i*sin(2*pi*n1(ii)*d(15)*coth1(j,ii)/l(ii)))/n1(ii);
       Q15(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d(15)*coth1(j,ii)/l(ii));
       R15(ii)=cos(2*pi*n1(ii)*d(15)*coth1(j,ii)/l(ii));
       M15=[S15(ii), P15(ii); Q15(ii), R15(ii)];

       % matrix 16 for layer 16
       S16(ii)=cos(2*pi*n2(ii)*d(16)*coth2(j,ii)/l(ii));
       P16(ii)=-(1i*sin(2*pi*n2(ii)*d(16)*coth2(j,ii)/l(ii)))/n2(ii);
       Q16(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(16)*coth2(j,ii)/l(ii));
       R16(ii)=cos(2*pi*n2(ii)*d(16)*coth2(j,ii)/l(ii));
       M16=[S16(ii), P16(ii); Q16(ii), R16(ii)];

       % matrix 17 for layer 17
       S17(ii)=cos(2*pi*n1(ii)*d(17)*coth1(j,ii)/l(ii));
       P17(ii)=-(1i*sin(2*pi*n1(ii)*d(17)*coth1(j,ii)/l(ii)))/n1(ii);
       Q17(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d(17)*coth1(j,ii)/l(ii));
       R17(ii)=cos(2*pi*n1(ii)*d(17)*coth1(j,ii)/l(ii));
       M17=[S17(ii), P17(ii); Q17(ii), R17(ii)];

       % matrix 18 for layer 18
       S18(ii)=cos(2*pi*n2(ii)*d(18)*coth2(j,ii)/l(ii));
       P18(ii)=-(1i*sin(2*pi*n2(ii)*d(18)*coth2(j,ii)/l(ii)))/n2(ii);
       Q18(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(18)*coth2(j,ii)/l(ii));
       R18(ii)=cos(2*pi*n2(ii)*d(18)*coth2(j,ii)/l(ii));
       M18=[S18(ii), P18(ii); Q18(ii), R18(ii)];

       % matrix 19 for layer 19
       S19(ii)=cos(2*pi*n1(ii)*d(19)*coth1(j,ii)/l(ii));
       P19(ii)=-(1i*sin(2*pi*n1(ii)*d(19)*coth1(j,ii)/l(ii)))/n1(ii);
       Q19(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d(19)*coth1(j,ii)/l(ii));
       R19(ii)=cos(2*pi*n1(ii)*d(19)*coth1(j,ii)/l(ii));
       M19=[S19(ii), P19(ii); Q19(ii), R19(ii)];

       % matrix 20 for layer 20
       S20(ii)=cos(2*pi*n2(ii)*d(20)*coth2(j,ii)/l(ii));
       P20(ii)=-(1i*sin(2*pi*n2(ii)*d(20)*coth2(j,ii)/l(ii)))/n2(ii);
       Q20(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(20)*coth2(j,ii)/l(ii));
       R20(ii)=cos(2*pi*n2(ii)*d(20)*coth2(j,ii)/l(ii));
       M20=[S20(ii), P20(ii); Q20(ii), R20(ii)];

       % matrix 21 for layer 21
       S21(ii)=cos(2*pi*n1(ii)*d(21)*coth1(j,ii)/l(ii));
       P21(ii)=-(1i*sin(2*pi*n1(ii)*d(21)*coth1(j,ii)/l(ii)))/n1(ii);
       Q21(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d(21)*coth1(j,ii)/l(ii));
       R21(ii)=cos(2*pi*n1(ii)*d(21)*coth1(j,ii)/l(ii));
       M21=[S21(ii), P21(ii); Q21(ii), R21(ii)];

       % matrix 22 for layer 22
       S22(ii)=cos(2*pi*n2(ii)*d(22)*coth2(j,ii)/l(ii));
       P22(ii)=-(1i*sin(2*pi*n2(ii)*d(22)*coth2(j,ii)/l(ii)))/n2(ii);
       Q22(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(22)*coth2(j,ii)/l(ii));
       R22(ii)=cos(2*pi*n2(ii)*d(22)*coth2(j,ii)/l(ii));
       M22=[S22(ii), P22(ii); Q22(ii), R22(ii)];

       % matrix 23 for layer 23
       S23(ii)=cos(2*pi*n1(ii)*d(23)*coth1(j,ii)/l(ii));
       P23(ii)=-(1i*sin(2*pi*n1(ii)*d(23)*coth1(j,ii)/l(ii)))/n1(ii);
       Q23(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d(23)*coth1(j,ii)/l(ii));
       R23(ii)=cos(2*pi*n1(ii)*d(23)*coth1(j,ii)/l(ii));
       M23=[S23(ii), P23(ii); Q23(ii), R23(ii)];

       % matrix 24 for layer 24
       S24(ii)=cos(2*pi*n2(ii)*d(24)*coth2(j,ii)/l(ii));
       P24(ii)=-(1i*sin(2*pi*n2(ii)*d(24)*coth2(j,ii)/l(ii)))/n2(ii);
       Q24(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(24)*coth2(j,ii)/l(ii));
       R24(ii)=cos(2*pi*n2(ii)*d(24)*coth2(j,ii)/l(ii));
       M24=[S24(ii), P24(ii); Q24(ii), R24(ii)];

       % matrix 25 for layer 25
       S25(ii)=cos(2*pi*n1(ii)*d(25)*coth1(j,ii)/l(ii));
       P25(ii)=-(1i*sin(2*pi*n1(ii)*d(25)*coth1(j,ii)/l(ii)))/n1(ii);
       Q25(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d(25)*coth1(j,ii)/l(ii));
       R25(ii)=cos(2*pi*n1(ii)*d(25)*coth1(j,ii)/l(ii));
       M25=[S25(ii), P25(ii); Q25(ii), R25(ii)];

       % matrix 26 for layer 26
       S26(ii)=cos(2*pi*n2(ii)*d(26)*coth2(j,ii)/l(ii));
       P26(ii)=-(1i*sin(2*pi*n2(ii)*d(26)*coth2(j,ii)/l(ii)))/n2(ii);
       Q26(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(26)*coth2(j,ii)/l(ii));
       R26(ii)=cos(2*pi*n2(ii)*d(26)*coth2(j,ii)/l(ii));
       M26=[S26(ii), P26(ii); Q26(ii), R26(ii)];

       % matrix 27 for layer 27
       S27(ii)=cos(2*pi*n1(ii)*d(27)*coth1(j,ii)/l(ii));
       P27(ii)=-(1i*sin(2*pi*n1(ii)*d(27)*coth1(j,ii)/l(ii)))/n1(ii);
       Q27(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d(27)*coth1(j,ii)/l(ii));
       R27(ii)=cos(2*pi*n1(ii)*d(27)*coth1(j,ii)/l(ii));
       M27=[S27(ii), P27(ii); Q27(ii), R27(ii)];

       % matrix 28 for layer 28
       S28(ii)=cos(2*pi*n2(ii)*d(28)*coth2(j,ii)/l(ii));
       P28(ii)=-(1i*sin(2*pi*n2(ii)*d(28)*coth2(j,ii)/l(ii)))/n2(ii);
       Q28(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(28)*coth2(j,ii)/l(ii));
       R28(ii)=cos(2*pi*n2(ii)*d(28)*coth2(j,ii)/l(ii));
       M28=[S28(ii), P28(ii); Q28(ii), R28(ii)];

       % matrix 29 for layer 29
       S29(ii)=cos(2*pi*n2(ii)*d(29)*coth1(j,ii)/l(ii));
       P29(ii)=-(1i*sin(2*pi*n2(ii)*d(29)*coth1(j,ii)/l(ii)))/n2(ii);
       Q29(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d(29)*coth1(j,ii)/l(ii));
       R29(ii)=cos(2*pi*n2(ii)*d(29)*coth1(j,ii)/l(ii));
       M29=[S29(ii), P29(ii); Q29(ii), R29(ii)];

       
       %matrix 41 for layer 41
       S41(ii)=cos(2*pi*n1(ii)*d1*coth1(j,ii)/l(ii));
       P41(ii)=-(1i*sin(2*pi*n1(ii)*d1*coth1(j,ii)/l(ii)))/n1(ii);
       Q41(ii)=-1i*n1(ii)*sin(2*pi*n1(ii)*d1*coth1(j,ii)/l(ii));
       R41(ii)=cos(2*pi*n1(ii)*d1*coth1(j,ii)/l(ii));
       M41=[S41(ii), P41(ii); Q41(ii), R41(ii)];
       %matrix 42 for layer 42
       S42(ii)=cos(2*pi*n2(ii)*d2*coth2(j,ii)/l(ii));
       P42(ii)=-(1i*sin(2*pi*n2(ii)*d2*coth2(j,ii)/l(ii)))/n2(ii);
       Q42(ii)=-1i*n2(ii)*sin(2*pi*n2(ii)*d2*coth2(j,ii)/l(ii));
       R42(ii)=cos(2*pi*n2(ii)*d2*coth2(j,ii)/l(ii));
       M42=[S42(ii), P42(ii); Q42(ii), R42(ii)];
       
       %matrix product and transmission
       M=(M29*M1*M2*M3*M4*M5*M6*M7*M8*M9*M10*M11*M12*M13*M14*M15*M16*M17*M18*M19*M20*M21*M22*M23*M24*M25*M26*M27*M28);

       p0=n0*cos(th(j));
       ps=ns*sqrt(1-((n0*n0*sin(th(j))*sin(th(j)))/(ns*ns)));
       te=2*p0/(((M(1,1)+M(1,2)*ps)*p0)+(M(2,1)+M(2,2)*ps));
       TE(j,ii)=C.*(ps/p0)*(abs(te)^2);%TE transmission

       M2=M29*((M41*M42)^14);

       te2=2*p0/(((M2(1,1)+M2(1,2)*ps)*p0)+(M2(2,1)+M2(2,2)*ps));
       TE2(j,ii)=C.*(ps/p0)*(abs(te2)^2);%TE transmission


       p0M=cos(th(j))/n0;
       psM=(sqrt(1-((n0*n0*sin(th(j))*sin(th(j)))/(ns*ns))))/ns;
       tm=2*p0M/(((M(1,1)+M(1,2)*psM)*p0M)+(M(2,1)+M(2,2)*psM));
       TM(j,ii)=C.*(psM/p0M)*(abs(tm)^2);%TM transmission

       tm2=2*p0M/(((M2(1,1)+M2(1,2)*psM)*p0M)+(M2(2,1)+M2(2,2)*psM));
       TM2(j,ii)=C.*(psM/p0M)*(abs(tm2)^2);%TM transmission



    end
end


hold on
subplot(2,1,1)
pcolor(l,thd,TE);shading interp;
set(gca,'fontsize',18,'FontWeight','bold');
ylabel('Angle (degrees)','FontWeight','bold');%y-axis label
colorbar
subplot(2,1,2)
pcolor(l,thd,TE2);shading interp;
set(gca,'fontsize',18,'FontWeight','bold');
ylabel('Angle (degrees)','FontWeight','bold');%y-axis label
colorbar
sgtitle('TE Modes','FontWeight','bold','FontSize',20)

figure

subplot(2,1,1)
pcolor(l,thd,TM);shading interp;
set(gca,'fontsize',18,'FontWeight','bold');
ylabel('Angle (degrees)','FontWeight','bold');%y-axis label
colorbar

subplot(2,1,2)
pcolor(l,thd,TM2);shading interp;
set(gca,'fontsize',18,'FontWeight','bold');
ylabel('Angle (degrees)','FontWeight','bold');%y-axis label
colorbar
sgtitle('TM Modes','FontWeight','bold','FontSize',20)
hold off