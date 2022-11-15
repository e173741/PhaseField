lengthx=0.06 ; lengthy=0.12 ;number=0 ; hpixel=80 ; hwidth=40 ;
thickness=0.06 ; v=0.3 ;displacementrate=-5*10^-6;
total=100 ; totalcracknumber=20 ;
xbound=0.9*lengthx ;
ybound=0.9*lengthy ;

maxcrack=10 ;
mincrack=3 ;
cracknumber=round((maxcrack-mincrack).*rand(totalcracknumber,1)+1) ;

crackmatrix=cell(totalcracknumber+7,1) ;

for i=1:totalcracknumber
 
cord1x=xbound.*rand(cracknumber(i,1),1) ;
cord1y=ybound.*rand(cracknumber(i,1),1) ;
cord2x=xbound.*rand(cracknumber(i,1),1) ;
cord2y=ybound.*rand(cracknumber(i,1),1) ;

for j=1:cracknumber(i,1)
len=((cord2x(j,1)-cord1x(j,1))^2+(cord2y(j,1)-cord1y(j,1))^2)^0.5 ;    
if len>0.01 && len<0.04 && cord1x(j,1)>0.1*lengthx && cord1x(j,1)<0.9*lengthx && cord2x(j,1)>0.1*lengthx && cord2x(j,1)<0.9*lengthx && cord1y(j,1)>0.1*lengthy && cord1y(j,1)<0.9*lengthy && cord2y(j,1)>0.1*lengthy && cord2y(j,1)<0.9*lengthy
crackmatrix{i,1}(j,1)=cord1x(j,1) ;
crackmatrix{i,1}(j,2)=cord1y(j,1) ;
crackmatrix{i,1}(j,3)=cord2x(j,1) ;
crackmatrix{i,1}(j,4)=cord2y(j,1) ;   
    
end   
end
end

crackmatrix=crackmatrix(~cellfun('isempty',crackmatrix)) ;
a=size(crackmatrix,1) ;
crackmatrix{a+1,1}=[0.02 0.06 0.04 0.06];
crackmatrix{a+2,1}=0.001*[30-5*(3^0.5) 65 30+5*(3^0.5) 55] ;
crackmatrix{a+3,1}=0.001*[30-5*(2^0.5) 60+5*(2^0.5) 30+5*(2^0.5) 60-5*(2^0.5)];
crackmatrix{a+4,1}=0.001*[30 70 30 50];
crackmatrix{a+5,1}=0.001*[20 60 40 60;30 50 30+10*(3^0.5) 40] ;
crackmatrix{a+6,1}=0.001*[20 60 40 60;30 50 40 50-10*(3^0.5)] ;
crackmatrix{a+7,1}=0.001*[20 60 40 60;30 50 30 30] ;
parametermatrix=zeros(total,6) ;
uppermodulus=90*10^9 ;
lowermodulus=30*10^9 ;
modulusmat=(uppermodulus-lowermodulus).*rand(total,1)+lowermodulus ;
upperGc1=100 ;
lowerGc1=30 ;
Gc1mat=(upperGc1-lowerGc1).*rand(total,1)+lowerGc1 ;
ratiomat=(7-3).*rand(total,1)+3 ;
Gc2mat=zeros(6,1) ;
for i=1:total
Gc2mat(i,1)=ratiomat(i,1)*Gc1mat(i,1) ;
end
lowerft=5*10^6 ;
upperft=15*10^6 ;
ftmat=(upperft-lowerft).*rand(total,1)+lowerft ;
cohesionupper=75*10^6 ;
cohesionlower=30*10^6 ;
cohesionmat=(cohesionupper-cohesionlower).*rand(total,1)+cohesionlower ;
frictionupper=50 ;
frictionlower=30 ;
frictionmat=(frictionupper-frictionlower).*rand(total,1)+frictionlower ;
parametermatrix(:,1)=modulusmat ;
parametermatrix(:,2)=Gc1mat ;
parametermatrix(:,3)=Gc2mat ;
parametermatrix(:,4)=ftmat;
parametermatrix(:,5)=cohesionmat ;
parametermatrix(:,6)=frictionmat ;
lc1=zeros(total,1) ;
for i=1:total
lc1(i,1)=27*modulusmat(i,1)*Gc1mat(i,1)/(256*ftmat(i,1)^2) ;
end
sumb=0;
for i=1:total
if lc1(i,1)<0.0028
    parametermatrix(i,:)=0 ;
sumb=sumb+1 ;
lc1(i,1)=0 ;
end
if lc1(i,1)>0.0035
      parametermatrix(i,:)=0 ;
sumb=sumb+1 ;
lc1(i,1)=0;
end
end

parametermatrix(all(~parametermatrix,2),:) = [];
lc1(all(~lc1,2),:) = [];
totalcracknumber=size(crackmatrix,1) ;
total=total-sumb ;
allphase=cell(total*3,8) ;
for i=1:total
modulus=parametermatrix(i,1) ;
Gc1=parametermatrix(i,2) ;
Gc2=parametermatrix(i,3) ;
ft=parametermatrix(i,4) ;
cohesion=parametermatrix(i,5) ;
friction=parametermatrix(i,6) ;
randoms=floor((total-1).*rand(3,1)+1) ;
for j=1:totalcracknumber
crack=crackmatrix{randoms(j,1),1} ;
[CNNinput,recordphase,t,cord,rec,phase,center,displace]=phasefieldsimulator1(lengthx,lengthy,thickness,modulus,Gc1,Gc2,ft,displacementrate,cohesion,friction,v,crack);
outputfile1=['phasefieldsimulatorPINN_',num2str(number),'.mat'] ;
allphase{number,1}=recordphase ;
allphase{number,2}=CNNinput ;
allphase{number,3}=t ;
allphase{number,4}=cord ;
allphase{number,5}=rec ;
allphase{number,6}=phase ;
allphase{number,7}=center ;
allphase{number,8}=displace ;
save(outputfile1,'recordphase','CNNinput','t','cord','rec','phase','center','displace') ;
number=number+1 ;
end
end
