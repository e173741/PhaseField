function [CNNdataout1,CNNdatain1,cordsq]=datapre(recordphase,hpixel,hwidth,center,rec,cord)
lengthy=0.12 ;lengthx=0.06 ;
elsizey=lengthy/hpixel ; elsizex=lengthx/hwidth ;

xm=(elsizex/2):elsizex:lengthx-elsizex/2 ;
ym=(elsizey/2):elsizey:lengthy-elsizey/2 ;
[X,Y]=meshgrid(xm,ym) ;
cordsq=zeros((hpixel)*(hwidth),2) ;

amk=-1 ;
for i=0:size(cordsq,1)-1
  cordsq(i+1,1)=X(1,mod(i,size(Y,2))+1) ;

end

for i=0:size(cordsq,1)-1
  if mod(i,size(Y,2))==0
      amk=amk+1;
  end
    azs=mod(amk,size(Y,1)) ;  anm=mod(i,size(Y,2))+1 ;
cordsq(i+1,2)=Y(azs+1,anm) ;

end
tot=size(cordsq,1) ;
indexarray=zeros(tot,10) ;

r=((elsizex^2+elsizey^2)^0.5)/2 ;
avg=zeros(tot,1) ;
count=1 ;
elementnumber=size(center,1) ;
for i=1:tot 
    for j=1:elementnumber
if ((center(j,1)-cordsq(i,1))^2+(center(j,2)-cordsq(i,2))^2)^0.5<r
indexarray(i,count)=j ;
 count=count+1 ;   
end
    end
avg(i,1)=count-1 ;
count=1;
end

tot=size(cordsq,1) ;
indexarray1=zeros(tot,10) ;

r=((elsizex^2+elsizey^2)^0.5)/2 ;
avg1=avg ;
count=1 ;
cordsize=size(cord,1);
for i=1:tot 
    for j=1:cordsize
if ((cord(j,1)-cordsq(i,1))^2+(cord(j,2)-cordsq(i,2))^2)^0.5<r
indexarray1(i,count)=j ;
 count=count+1 ;   
end
    end
avg1(i,1)=count ;    
count=1;
end
indexarray1( all(~indexarray1,2), : ) = [];
indexarray( all(~indexarray,2), : ) = [];

z=hwidth+1 ;
imageindex=zeros(tot,2) ;
con=0 ;
for i=1:tot+1000
 if mod(i,z)>1   
imageindex(i,1)=mod(i,z) ;
 end
if mod(i,z)==0 
    imageindex(i,1)=1 ;
    con=con+1;
end
imageindex(i,2)=hpixel-con  ;
end
imageindex(1,1)=1;
delarr=zeros(tot+1000,1);

for i=1:tot+1000
    for j=1:2
if imageindex(i,j)==0 || imageindex(i,j)<0
    delarr(i,1)=i ;
end
    end
end

delarr( all(~delarr,2), : ) = [];
imageindex(delarr(:,1),:)=[] ;
fihist=zeros(hpixel,hwidth) ;
Hhist=fihist ;
stress=cell(3,1) ;
for im=1:3
for i=1:hpixel
    for j=1:hwidth
stress{im,1}(i,j)=0;
      
    end
end
end
rec( all(~rec,2), : ) = [];
inc=size(rec,1) ;
CNNdatain1=cell(inc,3) ;
strainall=stress ;
tot=size(indexarray,1) ;
ninc=inc ;
CNNdataout1=cell(ninc-2,3) ;


 for am=1:ninc-1
sum1=0 ;
sum2=0 ;
for i=1:tot
    for j=1:3
    for im=1:avg(i,1)    
  if indexarray(i,im)>0
        sum1=sum1+recordphase{am,2}(indexarray(i,im),j) ;
        sum2=sum2+recordphase{am,3}(indexarray(i,im),j) ;
  end 
    end
    stress{j,1}(imageindex(i,2),imageindex(i,1))=sum1/avg(i,1) ;
    strainall{j,1}(imageindex(i,2),imageindex(i,1))=sum2/avg(i,1) ;
    sum1=0;
    sum2=0;
    end
end
sum3=0 ;
for i=1:tot
    for im=1:avg(i,1)    
  if indexarray(i,im)>0
   sum3=sum3+recordphase{am,1}(indexarray(i,im),1) ;
  end 
    end
    fihist(imageindex(i,2),imageindex(i,1))=sum3/avg(i,1) ;
    sum3=0;
end
sum4=0 ;
for i=1:tot
    for im=1:avg(i,1)    
  if indexarray(i,im)>0
   sum4=sum4+recordphase{am,4}(indexarray(i,im),1) ;
  end 
    end
    Hhist(imageindex(i,2),imageindex(i,1))=sum4/avg(i,1) ;
    sum4=0;
end
CNNdatain1{am,1}=fihist ;
for i=2:3
CNNdatain1{am,i}=stress{i-1,1} ;
end
 end

for amk=1:ninc-1
    
CNNdataout1{amk+1,1}=CNNdatain1{amk,1} ;
CNNdataout1{amk+1,2}=CNNdatain1{amk,2} ;
CNNdataout1{amk+1,3}=CNNdatain1{amk,3} ;
end


end


