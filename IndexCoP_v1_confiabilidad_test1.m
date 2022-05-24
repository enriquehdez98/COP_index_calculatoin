%% PROYECTO SEGUNDO PARCIAL
%%% Calculo de indices del CoP %%%
     %Enrique Hernandez Laredo - ehernandezl190@alumno.uaemex.mx
     %Lilyam Lizette Olmos García Rojas - lolmosg195@alumno.uaemex.mx    
warning('off')
close all;
clear all;
clc;
% Para toda la información se entiende que:
%%CoP ML  (medio-lateral)= ML0 = CoPx = Columna 7
%%CoP AP  (anterior-posterior)= AP0 = CoPy = Columna 8
[ML0A,AP0A,tm,Fs] = data_test1();
[M,N]=size(ML0A);
t=0:1/Fs:(length(ML0A(1,:))/(Fs))-(1/Fs);
%% Calculamos la media de AP
APmA=zeros(1,M);
for i=1:M
APmA(i)=mean(AP0A(i,:));
APma=mean(APmA);
end
%% Calculamos la media de ML
MLPmA=zeros(1,M);
for i=1:M
MLmA(i)=mean(ML0A(i,:));
MLma=mean(MLmA);
end
%% Calculamos AP COP medio.
APA=zeros(1,N);
for i=1:M
APA(i,:)=AP0A(i,:)-APmA(i);
end
%% Calculamos ML COP medio.
MLA=zeros(1,6000);
for i=1:M
MLA(i,:)=ML0A(i,:)-MLmA(i);
end
%% Calculamos RD
% La serie de tiempo de distancia resultante (RD)= es la distancia vectorial 
% desde el COP medio a cada par de puntos en las series de tiempo APo y MLo.
RDA=zeros(1,N);
for i=1:M
RDA(i,:)=sqrt(((APA(i,:)).^2)+(MLA(i,:).^2));
end
%% Calculamos MDIST
%La distancia media (MDIST) es la media de la serie de tiempo RD y representa 
%la distancia media desde el COP medio
lim_x=163;

MDISTA=zeros(1,M);
for i=1:M
MDISTA(i)=mean(RDA(i,:));
end
MDIST=[MDISTA(1:lim_x)'];
%% Calculamos MDIST_AP
% % La distancia media_AP es el valor absoluto medio de la serie 
% de tiempo AP y representa la distancia AP media desde el COP medio
MDIST_APA=zeros(1,M);
for i=1:M
MDIST_APA(i)=mean(abs(APA(i,:)));
end
MDIST_AP=[MDIST_APA(1:lim_x)'];

%% Calculamos MDIST_ML
% % La distancia media_ML es el valor absoluto medio de la serie 
% de tiempo AP y representa la distancia AP media desde el COP medio
MDIST_MLA=zeros(1,M);
for i=1:M
MDIST_MLA(i)=mean(abs(MLA(i,:)));
end
MDIST_ML=[MDIST_MLA(1:lim_x)'];

%% Calculamos RDIST
% La distancia rms (RDIST) desde el COP medio es el valor RMS de la serie 
% de tiempo RD
RDISTA=zeros(1,M);
for i=1:M
RDISTA(i)=sqrt(mean((RDA(i,:)).^2));
end
RDIST=[RDISTA(1:lim_x)'];

%%  Calculamos RDIST AP
RDIST_APA=zeros(1,M);
for i=1:M
RDIST_APA(i)=sqrt(mean((APA(i,:)).^2));
end
RDIST_AP=[RDIST_APA(1:lim_x)'];

%%  Calculamos RDIST ML
RDIST_MLA=zeros(1,M);
for i=1:M
RDIST_MLA(i)=sqrt(mean((MLA(i,:)).^2));
end
RDIST_ML=[RDIST_MLA(1:lim_x)'];

%%  Calculamos TOTEX
[n_subject,n_muestras]=size(RDA);
TOTEXA=zeros(1,M);
for i=1:M
    for j=1:n_muestras-1
        TOTEXA(i)=TOTEXA(i)+ sqrt(((((APA(i,j+1))-(APA(i,j))).^2)+(((MLA(i,j+1))-(MLA(i,j))).^2)));
     end
end
TOTEX=[TOTEXA(1:lim_x)'];

%% TOTEXAP
TOTEXAPA=zeros(1,M);
for i=1:M
    for j=1:n_muestras-1 
        TOTEXAPA(i)=TOTEXAPA(i)+ abs((APA(i,j+1))-(APA(i,j)));
    end
end
TOTEXAP=[TOTEXAPA(1:lim_x)'];

%% TOTEXML
TOTEXMLA=zeros(1,M);
for i=1:M
    for j=1:n_muestras-1
        TOTEXMLA(i)=TOTEXMLA(i)+ abs((MLA(i,j+1))-(MLA(i,j)));
    end
end
TOTEXML=[TOTEXMLA(1:lim_x)'];

%% VEL
% VEL_J=zeros(1,M);
VEL_A=zeros(1,M);              
VEL_A=TOTEXA/60;
VEL_=[VEL_A(1:lim_x)'];

%% VELOCIDAD AP
VEL_APA=zeros(1,M);             
VEL_APA=TOTEXAPA/60;
VEL_AP=[VEL_APA(1:lim_x)'];

%% VELOCIDAD ML
VEL_MLA=zeros(1,M);             
VEL_MLA=TOTEXMLA/60;
VEL_ML=[VEL_MLA(1:lim_x)'];

%% RANGE_ML
RANGE_MLA=abs(min(MLA)-max(MLA));
RANGE_ML=[RANGE_MLA(1:lim_x)'];

%% RANGE_AP
RANGE_APA=abs(min(APA)-max(APA));
RANGE_AP=[RANGE_APA(1:lim_x)'];

%% Time-Domain “Area” Measures
%% SRD
SRDA=sqrt(RDISTA.^2 - MDISTA.^2);
SRD=[SRDA(1:lim_x)'];

%%  AREACC
AREACCA=pi.*((MDISTA+(1.645.*SRDA)).^2);
AREACC=[AREACCA(1:lim_x)'];

%% SAPML
SAPMLA=zeros(1,M);
for i=1:M
        SAPMLA(i)=mean(APA(i,:).*MLA(i,:));
end
SAPML=[SAPMLA(1:lim_x)'];

%% AREACE   
AREACEA=6.*pi.*((((RDIST_APA.^2).*(RDIST_MLA.^2))-(SAPMLA.^2)).^0.5);
AREACE=[AREACEA(1:lim_x)'];

%% SW
SWA=zeros(1,M);
for i=1:M
    for j=1:n_muestras-1
        SWA(i)= SWA(i)+abs((APA(i,j+1).*MLA(i,j)) - (MLA(i,j+1).*APA(i,j)));
    end
end
SW=[SWA(1:lim_x)'];

%% AREASW 
AREASWA=SWA/(2*60);
AREASW=[AREASWA(1:lim_x)'];

%% MFREQ
MFREQA=(VEL_A)./(2.*pi.*MDISTA);
MFREQ=[MFREQA(1:lim_x)'];

%% MFREQAP
MFREQAPA=(VEL_APA)./(4.*sqrt(2).*MDIST_APA);
MFREQAP=[MFREQAPA(1:lim_x)'];

%% MFREQML
MFREQMLA=(VEL_MLA)./(4.*sqrt(2).*MDIST_MLA);
MFREQML=[MFREQMLA(1:lim_x)'];

%% dFDCC
dFDCCA=2*(MDISTA+(1.645.*SRDA));
dFDCC=[dFDCCA(1:lim_x)'];

%% FDCC
FDCCA=zeros(1,M);
for i=1:M
       FDCCA(i)=log(n_muestras)/log((n_muestras.*dFDCCA(i))/TOTEXA(i));
end
FDCC=[FDCCA(1:lim_x)'];

%% dFDCE
dFDCEA=sqrt(24.*sqrt(((RDIST_APA.^2).*(RDIST_MLA.^2))-SAPMLA.^2));
dFDCE=[dFDCEA(1:lim_x)'];

%% FDCE
FDCEA=zeros(1,M);
for i=1:M    
    FDCEA(i)=log(n_muestras)/log((n_muestras.*dFDCEA(i))/TOTEXA(i));
end
FDCE=[FDCEA(1:lim_x)'];

%% Densidad espectral de Potencia
% muestra 16=0.15Hz  y 301=5Hz 
% RDA
psRDA=zeros(M,292);
for i=1:M    
    x=RDA(i,:);
    N = length(x);
    xdft = fft(x);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(x):Fs/2;
    psRDA(i,:)=psdx(10:301);
end
%figure,plot(freq(10:301),psRDA(82,:))

% Total Power RD
POWER_RDT=sum(psRDA');
POWERRD_TA=[POWER_RDT'];

% 50% Power RD
POWER_RD50A=zeros(1,M);
POWER_RD50=cumsum(psRDA');
for i=1:M 
    POWER_RD50A(i)=freq( max(find(POWER_RD50(:,i)<=POWERRD_TA(i)*0.5)) );
end
POWER_RD50A=POWER_RD50A';

% 95% Power RD
POWER_RD95A=zeros(1,M);
POWER_RD95=cumsum(psRDA');
for i=1:M 
    POWER_RD95A(i)=freq( max(find(POWER_RD95(:,i)<=POWERRD_TA(i)*0.95)) );
end
POWER_RD95A=POWER_RD95A';

%APA
psAPA=zeros(M,292);
for i=1:M    
    x=APA(i,:);
    N = length(x);
    xdft = fft(x);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(x):Fs/2;
    psAPA(i,:)=psdx(10:301);
end
%figure,plot(freq(10:301),psAPA(82,:))

% Total Power AP
POWER_APT=sum(psAPA');
POWERAP_TA=[POWER_APT'];

% 50% Power AP
POWER_AP50A=zeros(1,M);
POWER_AP50=cumsum(psAPA');
for i=1:M 
    POWER_AP50A(i)=freq( max(find(POWER_AP50(:,i)<=POWERAP_TA(i)*0.5)) );
end
POWER_AP50A=POWER_AP50A';

% 95% Power AP
POWER_AP95A=zeros(1,M);
POWER_AP95=cumsum(psAPA');
for i=1:M 
    POWER_AP95A(i)=freq( max(find(POWER_AP95(:,i)<=POWERAP_TA(i)*0.95)) );
end
POWER_AP95A=POWER_AP95A';

%MLA
psMLA=zeros(M,292);
for i=1:M    
    x=MLA(i,:);
    N = length(x);
    xdft = fft(x);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(x):Fs/2;
    psMLA(i,:)=psdx(10:301);
end
%figure,plot(freq(10:301),psMLA(82,:))

% Total Power ML
POWER_MLT=sum(psMLA');
POWERML_TA=[POWER_MLT'];

% 50% Power ML
POWER_ML50A=zeros(1,M);
POWER_ML50=cumsum(psMLA');
for i=1:M 
    POWER_ML50A(i)=freq( max(find(POWER_ML50(:,i)<=POWERML_TA(i)*0.5)) );
end
POWER_ML50A=POWER_ML50A';

% 95% Power ML
POWER_ML95A=zeros(1,M);
POWER_ML95=cumsum(psMLA');
for i=1:M 
    POWER_ML95A(i)=freq( max(find(POWER_ML95(:,i)<=POWERML_TA(i)*0.95)) );
end
POWER_ML95A=POWER_ML95A';

% Momentos de densidad de potencia RD
u_0RD=POWERRD_TA;
deltaF=freq(2);m=0;k=1;
u_1RD=zeros(M,1);
for i=1:M 
    for m=1:292 
       u_1RD(i)=(((m-1).*deltaF).^(k)).*psRDA(i,m)+u_1RD(i);
    end
end

deltaF=freq(2);m=0;k=2;
u_2RD=zeros(M,1);
for i=1:M 
    for m=1:292   
       u_2RD(i)=(((m-1).*deltaF).^(k)).*psRDA(i,m)+u_2RD(i);
    end
end

CFREQ_RDA=(u_2RD./u_0RD).^(1/2);
FREQD_RDA=abs(((1-u_1RD.^2)./(u_0RD.*u_2RD)).^(1/2));

%%%%%%%%%%%%%%%%% potencia AP
u_0AP=POWERAP_TA;

deltaF=freq(2);m=0;k=1;
u_1AP=zeros(M,1);
for i=1:M 
    for m=1:292  
       u_1AP(i)=(((m-1)*deltaF)^(k))*psAPA(i,m)+u_1AP(i);
    end
end

deltaF=freq(2);m=0;k=2;
u_2AP=zeros(M,1);
for i=1:M 
    for m=1:292 
       u_2AP(i)=(((m-1)*deltaF)^(k))*psAPA(i,m)+u_2AP(i);
    end
end

CFREQ_APA=(u_2AP./u_0AP).^0.5;

FREQD_APA=abs(((1-u_1AP.^2)./(u_0AP.*u_2AP)).^0.5);

%%%%%%%%%%%%%%%%% potencia ML
u_0ML=POWERML_TA;

deltaF=freq(2);m=0;k=1;
u_1ML=zeros(M,1);
for i=1:M 
    for m=1:292  
       u_1ML(i)=(((m-1)*deltaF)^(k))*psMLA(i,m)+u_1ML(i);
    end
end

deltaF=freq(2);m=0;k=2;
u_2ML=zeros(M,1);
for i=1:M 
    for m=1:292 % muestra 16=0.15Hz  y 301=5Hz  
       u_2ML(i)=(((m-1)*deltaF)^(k))*psMLA(i,m)+u_2ML(i);
    end
end

CFREQ_MLA=(u_2ML./u_0ML).^0.5;

FREQD_MLA=abs(((1-u_1ML.^2)./(u_0ML.*u_2ML)).^0.5);
%%
X=[MDIST,MDIST_ML,MDIST_AP,RDIST,RDIST_ML,RDIST_AP,TOTEX,TOTEXML,TOTEXAP,VEL_,VEL_ML,VEL_AP,SRD,AREACC,SAPML,AREACE,AREASW,MFREQ,MFREQML,MFREQAP,FDCC,FDCE,RANGE_ML,RANGE_AP,POWERRD_TA, POWER_RD50A, POWER_RD95A, POWERAP_TA, POWER_AP50A, POWER_AP95A, POWERML_TA, POWER_ML50A, POWER_ML95A, CFREQ_RDA, FREQD_RDA, CFREQ_APA, FREQD_APA, CFREQ_MLA, FREQD_MLA];
X=X/10;