
all=1521 ;
bottom=938 ;
locstrainmax=zeros(all,1) ;
locstrainmin=locstrainmax ;
locstressmax=locstrainmax ;
locstressmin=locstrainmax ;


for i=bottom:all
file1=['phasefieldsimulator_',num2str(i),'.mat'] ;
A=load(file1) ;
A1=A.recordphase ;
A2=A.center ;
rec=A.rec ;
cordsize=size(A2,1) ;
rec(~any(rec,2),:)=[] ;
am=size(rec,1);

strainmax=zeros(am,1) ; strainmin=strainmax ; stressmax=strainmax ;stressmin=stressmax ;
for im=1:am-1
strainmax(im,1)=max(max(A1{im,2}(:,1:3))) ;
strainmin(im,1)=min(min(A1{im,2}(:,1:3))) ;
stressmax(im,1)=max(max(A1{im,3}(:,1:3))) ;
stressmin(im,1)=min(min(A1{im,3}(:,1:3))) ;
end
locstrainmax(i,1)=max(strainmax) ; 
locstrainmin(i,1)=min(strainmin) ;
locstressmax(i,1)=max(stressmax) ; 
locstressmin(i,1)=min(stressmin) ;
end

maxstrain=max(locstrainmax) ;
minstrain=min(locstrainmin) ;
maxstress=max(locstressmax) ;
minstress=min(locstressmin) ;
mstrain=(2/(maxstrain-minstrain)) ;
bstrain=(maxstrain+minstrain)/(minstrain-maxstrain) ;
mstress=(2/(maxstress-minstress)) ;
bstress=(maxstress+minstress)/(minstress-maxstress) ;

backpro=[mstrain bstrain;mstress bstress] ;

writematrix(backpro,'backpro.xlsx') ;

for i=bottom:all
file1=['phasefieldsimulator_',num2str(i),'.mat'] ;
A=load(file1);
A1=A.recordphase ;
A2=A.center ;
rec=A.rec ;
cordsize=size(A2,1) ;
rec(~any(rec,2),:)=[] ;
am=size(rec,1);
Am=cell(am,1) ;
top=1 ;

recordphasenorm=cell(am,4) ;

for iz=1:am-1
    for k1=1:3
    for j1=1:cordsize
       adam=A1{iz,2}(j1,k1) ;
recordphasenorm{iz,2}(j1,k1)=adam*mstrain+bstrain ;
    end
    end
end
for iz=1:am-1
    for k1=1:3
    for j1=1:cordsize
adam=A1{iz,3}(j1,k1) ;
        recordphasenorm{iz,3}(j1,k1)=adam*mstress+bstress ;
    end
    end
end

for is=1:am-1
recordphasenorm{is,1}=A1{is,1}(:,1) ;
end
Hmax=zeros(am-1,1);
for im=1:am-1
Hmax(im,1)=max(A1{im,4}(:,1)) ;
end

maxH=max(Hmax) ;
mh=1/(maxH) ;
for iz=1:am-1
        for j1=1:cordsize
 adam1=A1{iz,4}(j1,1);        
recordphasenorm{iz,4}(j1,1)=adam1*mh ;
        end
end

outputfile1=['recordphasenorm_',num2str(i),'.mat'] ;
save(outputfile1,'recordphasenorm')
end










