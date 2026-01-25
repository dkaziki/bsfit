load data_single %loads data with single precision
data=double(data_single);

[ndat,nchan]=size(data); 

fs=256; %sampling rate 
segleng=fs; %segments are 1 second long, i.e. frequency resolution is 1 Hz. 
epleng=2*fs; %divides data into 2 seconds epochs 
      % (this is not necessary, one can also set epleng to ndat. 
      % It just illustrates some results below)

segshift=segleng/2; % segments overlap by 50%

para=[]; %use defaults
freqpairs=[14,14]; %calculate cross-bsiepctrum at f1=f2=13Hz. 

%calculate cross-bispectrum
[bsall,nave]=data2bs_event(data,segleng,segshift,epleng,freqpairs,para); 

n=2; %number of assumed sources

% fit the model; 
        [a,d,err,err_all,bsmodel]=bsfit(bsall,n,para);

%%
addpath(genpath('c:/nolte/meth'));
load sa_eeg

A=mkfilt_eloreta(sa.V_medium); %eloreta filter for later use
grid=sa.grid_medium;
grid2=sa.cortex10K.vc;
cortex=sa.cortex10K;


        [ds,dsmag]=data2source(a,A);
        vout1=spatfiltergauss(dsmag(:,1),grid,.5,grid2);
        vout2=spatfiltergauss(dsmag(:,2),grid,.5,grid2);

        [Fout,wall]=moca_ncomp(ds);
        Foutmag=squeeze(sqrt(sum(Fout.^2,2)));

        vout1b=spatfiltergauss(Foutmag(:,1),grid,.5,grid2);
        vout2b=spatfiltergauss(Foutmag(:,2),grid,.5,grid2);

% makes the figure in source space showing mixed and demixed sources
        figure;
        faktor=.75;
        wfaktor=.4;
        para.myviewdir=[0 -1 1];
        axes('Position',[faktor*.3   faktor*.85 wfaktor wfaktor]);
        showsurface(cortex,para,vout1);colorbar off

        axes('Position',[faktor*.6 faktor*.85 wfaktor wfaktor]);
        showsurface(cortex,para,vout2);colorbar off
        text(-45,0,'mixed','fontsize',16,'fontweight','bold')

        axes('Position',[faktor*.3 faktor*.45 wfaktor wfaktor]);
        showsurface(cortex,para,vout1b);colorbar off

        axes('Position',[faktor*.6 faktor*.45 wfaktor wfaktor]);
        showsurface(cortex,para,vout2b);colorbar off
        text(-45,0,'unmixed','fontsize',16,'fontweight','bold')

