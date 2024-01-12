%PSA Exergy Analysis
clc
clear
close all
format shortEng

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS

%Ranges of Variables Under Consideration
ppu=0.79:0.00001:0.99999; %Product Purity
r=0.01:0.01:0.99;         %Recovery

%Constants
V=1;         %Volume (m^3) of Product
PO=101.325;  %Dead State Pressure (kPa, abs)
P=input('Input Product Pressure (kPa)   ');%Product Pressure (kPa)
T=288.15;    %Temperature (K)
R=8.314;     %Universal Gas Constant (kJ/kmol*K)
yO2_O=0.21;  %Dead State Oxygen Molar Concentration
yN2_O=0.79;  %Dead State Nitrogen Molar Concentration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OTHER CALCULATIONS BASED ON INPUTS

yN2_P=ppu; %The purity represents the molar fraction of N2 in the product
yO2_P=zeros(1,length(yN2_P));
for i=1:length(yN2_P)
    yO2_P(i)=1-yN2_P(i);
end
NP=P*V/(R*T);%kmols in the product mixture

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE VECTORS/MATRICES FOR CALCULATIONS

NN2P=zeros(length(ppu),1);          %Moles of N2 in product (kmol)
NO2P=zeros(length(ppu),1);          %Moles of O2 in product (kmol)
NN2E=zeros(length(ppu),length(r));  %Moles of N2 in exhaust (kmol)
NO2E=zeros(length(ppu),length(r));  %Moles of O2 in exhaust (kmol)
yN2_E=zeros(length(ppu),length(r)); %Molar concentration of N2 in exhaust
yO2_E=zeros(length(ppu),length(r)); %Molar concentration of O2 in exhaust
X=zeros(length(ppu),length(r));  %Exergy (kJ) for product at 800 kPa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXERGY CALCULATIONS

%Variables that only depend on product purity
for i=1:length(ppu)
        NN2P(i)=yN2_P(i)*NP;
        NO2P(i)=NP-NN2P(i);
end

%Variables that depend on product purity and recovery
for i=1:length(ppu)
    for j=1:length(r)
        NN2E(i,j)=NN2P(i)*((1/r(j))-1);
        NO2E(i,j)=((NN2P(i)+NN2E(i,j))/yN2_O)-NP-NN2E(i,j);
        yN2_E(i,j)=NN2E(i,j)/(NN2E(i,j)+NO2E(i,j));
        yO2_E(i,j)=1-yN2_E(i,j);
        X(i,j)=R*T*(NP*((PO/P)-1+log(P/PO))+NN2P(i)*...
            log(yN2_P(i)/yN2_O)+NO2P(i)*log(yO2_P(i)/yO2_O)+NN2E(i,j)*...
            log(yN2_E(i,j)/yN2_O)+NO2E(i,j)*log(yO2_E(i,j)/yO2_O));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot for Reversible Work Input
figure
per=ppu*100;
rec=r*100;
surf(rec,per,X);
shading interp
xlabel('Recovery (%)')
xticks([1 25 50 75 99])
xlim([1 99])
ylabel('Product Purity (%)')
yticks([79 85 90 95 99.999])
ylim([79 99.999])
zlabel('Reversible Work Input (kJ/m^3)')

%Plot for 99.9% Purity
figure
plot(X(20901,:))
xlabel('Recovery (%)')
ylabel('Reversible Work Input (kJ/m^3)')
title('Reversible Work Input vs Recovery at 99.9% Purity')