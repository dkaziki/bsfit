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
    