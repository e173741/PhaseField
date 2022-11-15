function [CNNinput,recordphase,t,cord,rec,phase,center,displace]=phasefieldsimulator1(lengthx,lengthy,thickness,modulus,Gc1,Gc2,ft,displacementrate,cohesion,friction,v,crack);
lc=27*modulus*Gc1/(256*ft^2) ;
elementsize=lc/2 ;
fd=@(p) drectangle(p,0,lengthx,0,lengthy);
[p,t]=distmesh2d(fd,@huniform,elementsize,[0,0;lengthx,lengthy],[0,0;lengthx,0;0,lengthy;lengthx,lengthy]);
%%%%%%support=[0 0 0.45*lengthx 0 0 1;0.45*lengthx 0 0.55*lengthy 0 1 1;0.55*lengthx 0 lengthx 0 0 1] ;
support=[0 0 lengthx 0 1 1] ;
loading=[0.01*lengthx lengthy 0.99*lengthx lengthy 0 -100] ;
mark=[0.5*lengthx lengthy 2] ;
ninc=100 ;
kb=0;  tolerance=0.45; 
cordsize=size(p,1) ;  elementnumber=size(t,1) ;
cord=p;
total=cordsize ;
F=zeros(total,2) ;
sf=size(support,1) ;
Es=modulus/((1-2*v)*(1+v)); G=0.5*modulus/(1+v) ;
E=[Es*(1-v) Es*v 0 ;Es*v Es*(1-v) 0 ; 0 0 G]; 
ind=zeros(cordsize,1);
elementsize=zeros(elementnumber,1) ;
for i=1:elementnumber
elementsize(i,1)=((cord(t(i,2),1)-cord(t(i,1),1))^2+(cord(t(i,2),2)-cord(t(i,1),2))^2)^0.5 ;
end
elsize1=mean(elementsize) ;
angle=zeros(sf,1);
for i=1:sf
leng=((support(i,4)-support(i,2))^2+(support(i,3)-support(i,1))^2)^0.5 ;
angle(i,1)=abs((support(i,4)-support(i,2))/leng) ;
end
sf1=size(loading,1); 
angle1=zeros(sf1,1) ;
for i=1:sf1
leng1=((loading(i,4)-loading(i,2))^2+(loading(i,3)-loading(i,1))^2)^0.5 ;
angle1(i,1)=abs((loading(i,4)-loading(i,2))/leng1) ;
end


for i=1:sf
for j=1:total
    if round(angle(i,1))==0 && support(i,1)<=cord(j,1) && cord(j,1)<=support(i,3) && support(i,2)-tolerance*elementsize(j,1)<cord(j,2) && support(i,4)+tolerance*elementsize(j,1)>cord(j,2)
    F(j,1)=support(i,5) ;
    F(j,2)=support(i,6) ;
    ind(j,1)=j ;
    end  
   if round(angle(i,1))==1 && support(i,2)<=cord(j,2) && cord(j,2)<=support(i,4) && support(i,1)-tolerance*elementsize(j,1)<cord(j,1) && support(i,3)+tolerance*elementsize(j,1)>cord(j,1)
    F(j,1)=support(i,5) ;
    F(j,2)=support(i,6) ;
    ind(j,1)=j ;
   end
   
end
end

B=1 ;
 for i=1:total
         for j=1:2
      if F(i,j)==1
          F(i,j)=0 ;
      else         
          F(i,j)=B ;
           B=B+1 ;
      end
         end
 end
 
 lm=zeros(elementnumber,6) ;
 
 for i=1:elementnumber
     for j=1:2
         lm(i,j)=F(t(i,1),j) ;
         lm(i,j+2)=F(t(i,2),j) ;
         lm(i,j+4)=F(t(i,3),j) ;
     end
 end

 sigma=cell(elementnumber,1);
 fuint=cell(elementnumber,1); 
 fdintpos=fuint ;fdintneg=fuint  ;
 d=cell(elementnumber,1) ;
 fidpos=cell(elementnumber,1);
  fidneg=cell(elementnumber,1);
 strain=sigma; 
 straine=sigma ;
 strainp=sigma ;
 sigmatrial=sigma;
 S=sigma ;
for i=1:elementnumber
     for sez=1:3
        sigma{i,1}(sez,1)=0;
          strain{i,1}(sez,1)=0;
          straine{i,1}(sez,1)=0;
          strainp{i,1}(sez,1)=0;
          sigmatrial{i,1}(sez,1)=0;
          S{i,1}(sez,1)=0;
          fdintpos{i,1}(sez,1)=0 ;
          fdintneg{i,1}(sez,1)=0 ;
         fidpos{i,1}(sez,1)=0 ;
         fidneg{i,1}(sez,1)=0 ;
     end
end
for i=1:elementnumber
     for sez=1:6
     d{i,1}(sez,1)=0 ;   
   fuint{i,1}(sez,1)=0;
   
     end
end

ne=0; load=size(loading,1) ;
numb=zeros(load,1) ;
indi=cell(load,1); 
for i=1:load
     for sez=1:cordsize
     indi{i,1}(sez,1)=0 ;   
    
     end
end



for i=1:load
    for j=1:cordsize
 
       distance=abs((loading(i,4)-loading(i,2))*cord(j,1)-(loading(i,3)-loading(i,1))*cord(j,2)+loading(i,3)*loading(i,2)-loading(i,4)*loading(i,1))/((loading(i,4)-loading(i,2))^2+(loading(i,3)-loading(i,1))^2)^0.5 ;
      if distance<elsize1*2 && round(angle1(i,1))==0 && loading(i,1)<=cord(j,1) && loading(i,3)>=cord(j,1)
          ne=ne+1; 
          indi{i,1}(j,1)=j ;
      end
      
      if distance<elsize1*1.5 && round(angle1(i,1))==1 && loading(i,2)<=cord(j,2) && loading(i,4)>=cord(j,2)
          ne=ne+1; 
          indi{i,1}(j,1)=j ;
      end
         
    end
    numb(i,1)=ne ;
    ne=0;
end


fivpos=zeros(cordsize,1) ;
center=zeros(elementnumber,2) ;

for i=1:elementnumber
center(i,1)=(cord(t(i,1),1)+cord(t(i,2),1)+cord(t(i,3),1))/3 ;
center(i,2)=(cord(t(i,1),2)+cord(t(i,2),2)+cord(t(i,3),2))/3 ;
end


        disp=zeros(size(mark,1),1) ;
    for i=1:size(mark,1)
        for j=1:cordsize
     dia=((mark(i,1)-cord(j,1))^2+(mark(i,2)-cord(j,2))^2)^0.5 ;    
     if dia<elsize1
    disp(i,1)=j ;
         
     end
        end
    end
    
  
yuk=zeros(cordsize,3);
yuk1=yuk ;
for i=1:1
for j=1:cordsize
if indi{i,1}(j,1)>0
    yuk(indi{i,1}(j,1),1)=j ;
    yuk(indi{i,1}(j,1),2)=loading(i,5)/numb(i,1) ;
    yuk(indi{i,1}(j,1),3)=loading(i,6)/numb(i,1) ;
end
end
end

for i=2:load
for j=1:cordsize
if indi{i,1}(j,1)>0
    yuk1(indi{i,1}(j,1),1)=j ;
    yuk1(indi{i,1}(j,1),2)=loading(i,5)/numb(i,1) ;
    yuk1(indi{i,1}(j,1),3)=loading(i,6)/numb(i,1) ;
end
end
end
momenty=0 ;
for j=1:cordsize
if indi{1,1}(j,1)>0
momenty=momenty+(cord(indi{1,1}(j,1),1)-0.5*lengthx)*yuk(indi{1,1}(j,1),3) ;
end
end
for j=1:cordsize
if momenty<0
if indi{1,1}(j,1)>0 && cord(indi{1,1}(j,1),1)-0.5*lengthx<0 
yuk(indi{1,1}(j,1),3)=yuk(indi{1,1}(j,1),3)-momenty/(cord(indi{1,1}(j,1),1)-0.5*lengthx) ;
break
end
end
end
for j=1:cordsize
if momenty>0
if indi{1,1}(j,1)>0 && cord(indi{1,1}(j,1),1)-0.5*lengthx>0 
yuk(indi{1,1}(j,1),3)=yuk(indi{1,1}(j,1),3)-momenty/(cord(indi{1,1}(j,1),1)-0.5*lengthx) ;
break
end
end
end


kdpos=cell(elementnumber,1) ;
Bnu=kdpos ; ku=kdpos ;Bnd=Bnu ;

Nd=[1/3 1/3 1/3] ;
k=cell(elementnumber,1) ;
fipos=zeros(elementnumber,1) ;
fineg=zeros(elementnumber,1) ;
L=zeros(B-1,1); D=L ;  J=zeros(elementnumber,1) ;
L1=zeros(B-1,1); 
 inc=1 ;
nodalresponse=zeros(B-1,1) ;
nodalresponsedpos=zeros(cordsize,1); 
nodalresponsedneg=zeros(cordsize,1); 
inctot=1 ;
increment=zeros(ninc,1); 
avg=zeros(ninc,1) ;
nodalresponse1=nodalresponse;
incdif=1; ratio=1 ; rec=zeros(ninc,2);
countin=0; 
deltaUf=nodalresponse ; 
error=zeros(500,1); err=zeros(ninc,1); mov=struct('cdata',[],'colormap',[]);g=zeros(B-1,1);
g(F(disp,mark(1,3)),1)=1 ;
Y=L;
H=zeros(elementnumber,1) ; 
H1=H ;A=H ; H1pos=H; H1neg=H; Hpos=H ; Hneg=H ; Y1=Y ;

crack( all(~crack,2), : ) = [];

numberofcrack=size(crack,1) ;
cracklength=zeros(numberofcrack,1);

for i=1:numberofcrack
    cracklength(i,1)=((crack(i,4)-crack(i,2))^2+(crack(i,3)-crack(i,1))^2)^0.5 ;
end

for i=1:numberofcrack  
   xmin=min(crack(i,1),crack(i,3)) ; xmax=max(crack(i,1),crack(i,3)) ; ymin=min(crack(i,2),crack(i,4)) ; ymax=max(crack(i,2),crack(i,4)) ;
    for j=1:elementnumber
        if xmin-elsize1<center(j,1) && xmax+elsize1>center(j,1) && ymax+elsize1>center(j,2) && ymin-elsize1<center(j,2) 
                 distance=abs((crack(i,4)-crack(i,2))*center(j,1)-(crack(i,3)-crack(i,1))*center(j,2)+crack(i,3)*crack(i,2)-crack(i,4)*crack(i,1))/cracklength(i,1) ;
         Hpos(j,1)=10000*(Gc1/(2*lc))*exp(-(distance/(lc/10))^2) ; 
        end
       
    end
end


for i=1:cordsize
for j=1:2
 if yuk(i,1)~=0
    W=F(i,j);
if W~=0  
Y(W)=(yuk(i,j+1));
end
 end
end
end

for i=1:cordsize
for j=1:2
 if yuk1(i,1)~=0
    W=F(i,j);
if W~=0  
Y1(W)=(yuk1(i,j+1));
end
 end
end
end


lambda=0 ;
non=1; totalY=zeros(B-1,1); maxstrain=0;displacementrate1=displacementrate;Y2=Y;

%%%alfa1=(tan(friction*pi/180))/(9+12*(tan(friction*pi/180))^2)^0.5 ;
%%%k1=(3*cohesion)/(9+12*(tan(friction*pi/180))^2)^0.5 ;

%%%alfa1=(2*sin(friction*pi/180))/((3-sin(friction*pi/180))*3^0.5) ;
%%%k1=6*cohesion*cos(friction*pi/180)/((3-sin(friction*pi/180))*3^0.5) ;
sm=0 ;
Ep=cell(elementnumber,1) ;
modulus1=zeros(elementnumber,1) ;
Ew=cell(elementnumber,1);

lame=(v*modulus)/((1+v)*(1-2*v)) ;
nu=modulus/(2*(1+v)) ;
K0=lame+nu ;

for i=1:elementnumber
Ew{i,1}=E ;
end
index=[1 1 1 1 1 1;1 1 2 2 1 2;1 1 1 2 1 3;2 2 1 1 2 1;2 2 2 2 2 2;2 2 1 2 2 3;1 2 1 1 3 1;1 2 2 2 3 2;1 2 1 2 3 3] ;
iden=[1 1 0;1 1 0;0 0 0] ;
Ppos=zeros(3,3) ;
Pplus=Ppos ;
Pminus=Ppos ;  recordphase=cell(ninc,5) ;phase=cell(ninc,1); displace=phase;
ftrial=zeros(elementnumber,1) ;energypos=ftrial; energyneg=ftrial ;eqp=ftrial;
while ninc>inc

fuint1=fuint ;
fdint1pos=fdintpos ;
   
for i=1:elementnumber
    if inctot<2
x1=cord(t(i,1),1);
x2=cord(t(i,2),1) ;
x3=cord(t(i,3),1) ;
y1=cord(t(i,1),2) ;
y2=cord(t(i,2),2) ;
y3=cord(t(i,3),2) ;
b1=y2-y3 ; b2=y3-y1 ; b3=y1-y2 ; c1=x3-x2; c2=x1-x3 ;c3=x2-x1 ;
Ax=[1 x1 y1;1 x2 y2;1 x3 y3]; 
A(i,1)=0.5*det(Ax) ;
Bnu{i,1}=(1/(2*A(i,1)))*[b1 0 b2 0 b3 0;0 c1 0 c2 0 c3;c1 b1 c2 b2 c3 b3];
Bnd{i,1}=(1/(2*A(i,1)))*[b1 b2 b3;c1 c2 c3] ;
J(i,1)=(x1-x3)*(y2-y3)-(x2-x3)*(y1-y3) ; 
    end
fipos(i,1)=Nd*fidpos{i,1} ;
fineg(i,1)=Nd*fidneg{i,1} ;
strain{i,1}=strain{i,1}+Bnu{i,1}*d{i,1} ;

%%% burayi cehck etmemiz lazim birde plastic energy functionini tam ne oldugunu bilmiyoruz onada bakmak lazim
strainmatrix(1,1)=strain{i,1}(1,1) ;
strainmatrix(2,2)=strain{i,1}(2,1) ;
strainmatrix(1,2)=0.5*strain{i,1}(3,1) ;
strainmatrix(2,1)=strainmatrix(1,2) ;
devstrainmatrix=strainmatrix-0.5*(strainmatrix(1,1)+strainmatrix(2,2))*[1 0;0 1] ;
sigmapositive=K0*max(0,strainmatrix(1,1)+strainmatrix(2,2))*eye(2)+2*nu*(devstrainmatrix) ;
sigmanegative=K0*min(0,strainmatrix(1,1)+strainmatrix(2,2))*eye(2) ;
sigmatot=sigmapositive+sigmanegative ;
sigmaprinciple=eig(sigmatot) ;
energypos(i,1)=0.5*K0*((max(0,strainmatrix(1,1)+strainmatrix(2,2)))^2)+nu*(strainmatrix(1,1)^2+strainmatrix(2,2)^2+2*strainmatrix(1,2)^2) ;
energyneg(i,1)=0.5*K0*((min(0,strainmatrix(1,1)+strainmatrix(2,2)))^2) ;
sigma11=max(sigmaprinciple(2),sigmaprinciple(1)) ;
sigma33=min(sigmaprinciple(2),sigmaprinciple(1)) ;
Epositive=K0*heaviside(strainmatrix(1,1)+strainmatrix(2,2))*iden+2*nu*[0.5 -0.5 0;-0.5 0.5 0;0 0 0.5] ;
Enegative=K0*heaviside(-(strainmatrix(1,1)+strainmatrix(2,2)))*iden ;
sigmapos(1,1)=sigmapositive(1,1); sigmapos(2,1)=sigmapositive(2,2); sigmapos(3,1)=sigmapositive(1,2) ;
sigmaneg(1,1)=sigmanegative(1,1); sigmaneg(2,1)=sigmanegative(2,2); sigmaneg(3,1)=sigmanegative(1,2) ;
sigma1=((1-fipos(i,1))^2+kb)*(sigmapos+sigmaneg) ;
sigma{i,1}=sigma1 ;
wc=(max(0,((0.5*(sigma11-sigma33)/cos(friction*pi/180))+(0.5*(sigma11+sigma33)*tan(friction*pi/180)))-cohesion))^2 ;
if inctot>3
Ew{i,1}=((1-fipos(i,1))^2+kb)*(Epositive+Enegative) ;
end
H1pos(i,1)=Hpos(i,1) ;
Hpos(i,1)=(wc/(2*G))   ; 
if H1pos(i,1)>=Hpos(i,1)
    Hpos(i,1)=H1pos(i,1) ;
 
end
fuint{i,1}=A(i,1)*thickness*transpose(Bnu{i,1})*sigma1 ; 
ku{i,1}=A(i,1)*thickness*transpose(Bnu{i,1})*Ew{i,1}*Bnu{i,1} ; 
end
 

for i=1:elementnumber
kdpos{i,1}=((Gc2/lc+2*Hpos(i,1))*transpose(Nd)*Nd)*A(i,1)*thickness+(Gc2*lc*transpose(Bnd{i,1})*Bnd{i,1})*A(i,1)*thickness ; %%%% burdada sorun var plasticity koymadigimiz icin  damage 1 den buyuk cikiyor. 
fdintpos{i,1}=(((Gc2/lc)*fipos(i,1)-2*(1-fipos(i,1))*Hpos(i,1))*transpose(Nd))*thickness*A(i,1)+(Gc2*lc*transpose(Bnd{i,1})*Bnd{i,1}*fidpos{i,1})*thickness*A(i,1) ; %%% hata var burda
end


Ku=zeros(B-1,B-1) ;

  for a=1:elementnumber
    for i=1:6
       for j=1:6 
       P=lm(a,i) ;
        Q=lm(a,j) ;
            if P~=0 && Q~=0  
                Ku(P,Q)=Ku(P,Q)+ku{a,1}(i,j);
                
            end
        end
    end
  end
Kdpos=zeros(cordsize,cordsize); 
   for a=1:elementnumber
    for i=1:3
       for j=1:3 
       P=t(a,i) ;
        Q=t(a,j) ;
            if P~=0 && Q~=0  
              Kdpos(P,Q)=Kdpos(P,Q)+kdpos{a,1}(i,j);                  
            end
        end
    end
  end 

for i=1:elementnumber
    for j=1:6
P=lm(i,j) ;
if P~=0
nodalresponse(P,1)=nodalresponse(P,1)+fuint{i,1}(j,1)-fuint1{i,1}(j,1) ;
end
    end
end

if non==1  %%%%% papaeia bi daha bak cok major bir hata yapiyorsun ordan sikinti oluiyor. 
Uft=Ku\Y ;
deltalambda=displacementrate/(transpose(g)*Uft) ;
deltaUf=deltalambda*Uft ;
D=D+deltaUf ;
lambda=lambda+deltalambda ;
end 
if non>1
deltaUfr=Ku\(lambda*Y-nodalresponse);
Uft=Ku\Y ;
deltalambda=-(transpose(g)*deltaUfr)/(transpose(g)*Uft) ;
lambda=lambda+deltalambda ;
deltaUf=deltaUfr+deltalambda*Uft ;
D=D+deltaUf ;
L=lambda*Y ;
end
  


for i=1:elementnumber
    for j=1:3
P=t(i,j) ;
if P~=0
nodalresponsedpos(P,1)=nodalresponsedpos(P,1)+fdintpos{i,1}(j,1)-fdint1pos{i,1}(j,1) ;
end
    end
end  


for  iz=1:elementnumber
     for iw=1:6
        if  lm(iz,iw)==0  
d{iz,1}(iw,1)=0;
    else
d{iz,1}(iw,1)=deltaUf(lm(iz,iw),1);
        end
     end
end

er1=norm(L-nodalresponse)/norm(L) ;
error(inctot,1)=er1 ;

if er1<10^-10 && non>1  %%% phase ile ilgili bir kriter ekleyebiliriz
    increment(inc,1)=incdif ;
    avg(inc,1)=inctot/inc ;
    rec(inc,2)=sum(L)/(lengthx*thickness) ;
    rec(inc,1)=D(F(disp,mark(1,3)),1)/lengthy;
    err(inc,1)=er1 ;
    incdif=0 ;
    
deltafivpos=Kdpos\-nodalresponsedpos;  
fivpos=fivpos+deltafivpos ; 

for  iz=1:elementnumber
     for iw=1:3
         if  t(iz,iw)==0  
fidpos{iz,1}(iw,1)=0;
    else
fidpos{iz,1}(iw,1)=fivpos(t(iz,iw),1);
        end
     end
end
phase{inc,1}=fivpos;

for i=1:cordsize
    for j=1:2
  W=F(i,j) ;
  if W==0
    displace{inc,1}(i,j)=0 ;
  end
  if W~=0  
displace{inc,1}(i,j)=D(W,1) ;
  end
    end
end

recordphase{inc,1}=fipos;
recordphase{inc,4}=Hpos ;
for i=1:elementnumber
    for j=1:3
recordphase{inc,2}(i,j)=strain{i,1}(j,1) ;
recordphase{inc,3}(i,j)=sigma{i,1}(j,1) ;
    end
end
mov(inc) = getframe;
scatter(rec(:,1),rec(:,2));
title({'force disp ',num2str(incdif);'time ',inc}) ;  

if inc>2 && increment(inc,1)>increment(inc-1,1) && increment(inc,1)>3
    displacementrate=displacementrate1*(3/increment(inc,1));
   
end
if inc>2 && abs(rec(inc,2))<abs(rec(inc-1,2))
    displacementrate=displacementrate1/20;
   Y=Y2*(1/20);
end
inc=inc+1; 
non=1;
end
if er1>10^-10
    non=non+1;
end
if rcond(Ku)<10^-15
    inc=ninc;
end
if inc>10 && abs(sum(L))<0.1*max(abs(rec(:,2))) 
    inc=ninc ;

    
    
end    

inctot=inctot+1 ;
incdif=incdif+1;
end 
  
     figure
patch('Faces',t,'Vertices',cord,'CData',fivpos,'Facecolor','interp','Edgecolor','none');
colormap('jet')
daspect([1 1 1]); 
set(gca,'Visible','off')
  
     figure
patch('Faces',t,'Vertices',cord,'CData',recordphase{5,1},'Facecolor','interp','Edgecolor','none');
colormap('jet')
daspect([1 1 1]); 
set(gca,'Visible','off')

CNNinput{1,1}=Gc1 ;
CNNinput{1,2}=modulus ;
CNNinput{1,3}=Gc2 ;
CNNinput{1,4}=lc ;
CNNinput{1,5}=friction ;
CNNinput{1,6}=cohesion ;
end

