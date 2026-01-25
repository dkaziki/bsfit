function [cs,nave]=data2bs_event(data,segleng,segshift,epleng,freqpairs,para)
% calculates cross-bispectral tensors  from data in general for event related
% measurements
%
% usage: [cs,nave]=data2bs_event(data,segleng,segshift,epleng,freqpairs,para);
%
% input: 
% data: ndat times nchan matrix each colum is the time-series in one
%             channel;
% segleng: length of each segment in bins, e.g. segleng=1000;  
% segshift: number of bins by which neighboring segments are shifted;
%           e.g. segshift=segleng/2 makes overlapping segments
% epleng: length of each epoch
% freqpairs: Kx2 matrix for K pairs of frequencies. Pairs of  frequency in bins. 
%            The meaning depends on the 
%            on the physical length of the segments. E.g., if the sampling
%            rate is 1000 Hz and segleng=500, then the segment is 500 ms
%            long, and the frequency resolution is .5 Hz. Then, e.g. 
%            f1=6 corresponds to 10 Hz, because counting starts from zero 
%            and 10 Hz is the 6.th frequency for that frequency resolution.
%            freqpairs=[6,6] would correspond to the coupling between 10
%            and 20 Hz, i.e. in physical units: f1=10Hz, f2=10Hz and f3=f1+f2=20Hz. 
%            
% para: structure which is eventually used later
%
% output: 
% cs: nchan by nchan by nchan by number_of_frequency_pairs 
%      (by number_of_segments) tensor such that 
%  cs(i,j,k,f)=<x(f1)_i*x(f2)_j*conj(x(f1+f2-1)_k)>
%  where f1=freqpairs(f,1) and  f2=freqpairs(f,1),
%  x=fft(data) and the average is over epeochs and segments
% 
% if para.fave=0 then cs contains a fifth argument denoting 
% the   segment. 
% if para.fave=1 or ommited, then cs was averaged over segments. 

% nave: number of averages

[ndat,nchan]=size(data);
[nf,ndum]=size(freqpairs);

fave=1;
mywindow=repmat(hanning(segleng),1,nchan);
if nargin>5
    if isfield(para,'fave')
      fave=para.fave;
     end;
    if isfield(para,'mywindow');
       mywindow=repmat(para.mywindow,1,nchan);
    end
end





nep=floor(ndat/epleng);

nseg=floor((epleng-segleng)/segshift)+1; %total number of segments

if fave==1
 cs=zeros(nchan,nchan,nchan,nf);
% csn=zeros(nchan,nchan,nchan,nf);
else
 cs=zeros(nchan,nchan,nchan,nf,nseg);
% csn=zeros(nchan,nchan,nchan,nf,nseg);
end


csloc=zeros(nchan,nchan,nchan);
%figure;plot(mywindow);
nave=0;
for j=1:nep;
     %disp(j)
    dataep=data((j-1)*epleng+1:j*epleng,:);
    for i=1:nseg; %average over all segments;
      dataloc=dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);
      datalocfft=fft(detrend(dataloc).*mywindow);
      for f=1:nf % for all frequencies
         f1=freqpairs(f,1);f2=freqpairs(f,2);
         xx=transpose(datalocfft(f1,:))*datalocfft(f2,:);
         for k=1:nchan;
            csloc(:,:,k)=xx*conj(datalocfft(f1+f2-1,k));
            %cslocn(:,:,k)=(abs(datalocfft(f1,:)).^2)'*abs(datalocfft(f2,:)).^2*abs((datalocfft(f1+f2-1,k))).^2;
         end    
         if fave==1
            cs(:,:,:,f)=cs(:,:,:,f)+csloc;
            % csn(:,:,:,f)=csn(:,:,:,f)+cslocn;
         else             
            cs(:,:,:,f,i)=cs(:,:,:,f,i)+csloc;
           % csn(:,:,:,f,i)=csn(:,:,:,f,i)+cslocn;
         end
      end
      if fave==1; nave=nave+1;end
    end 
    if ~(fave==1); nave=nave+1;end
end

if fave==1
  for i=1:nchan;for j=1:nchan;
    cs(i,j,:,:)=cs(i,j,:,:)/nave;
  % csn(i,j,:,:)=sqrt(csn(i,j,:,:)/nave);
  end;end
else
  for i=1:nchan;for j=1:nchan;
    cs(i,j,:,:,:)=cs(i,j,:,:,:)/nave;
  %   csn(i,j,:,:,:)=sqrt(csn(i,j,:,:,:)/nave);
  end;end
end


    
    
   
    

return;