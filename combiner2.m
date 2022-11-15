bottom=1 ;
total=1521 ;
sizematrix=zeros(total-bottom+1,2) ;
CNNdatainlast=cell(total-bottom+1,1) ;
CNNdataoutlast=CNNdatainlast ;
CNNinputlast=CNNdatainlast ;

for i=bottom:total
file2=['phasefieldsimulator_',num2str(i),'.mat'] ;
c=load(file2) ;
s = dir(file2);         
rec1=c.rec ;
rec1( ~any(rec1,2), : ) = [];
sizematrix(i-bottom+1,1)=size(rec1,1) ;

if s.bytes>15*10^6
sizematrix(i-bottom+1,2)=1 ;
end
end
sz=0; pos=zeros(total-bottom+1,1) ;
for i=bottom:total
if sizematrix(i-bottom+1,2)>0 && sizematrix(i-bottom+1,1)>=47
sz=sz+1;
pos(i-bottom+1,1)=sz ;
end
end



multiin=zeros(sz,45,80,40,3) ;
multiout=multiin ;
i=1; im=1;
while i<sz
if sizematrix(i,2)>0 && sizematrix(i,1)>=47
file1=['CNNdatainnew_',num2str(i+bottom-1),'.mat'] ;
A2x=load(file1) ;
A2=A2x.CNNdatain1 ;
CNNdatainlast{i+bottom-1,1}=A2 ;
for j=1:45
    for k=1:3

multiout(im,j,1:80,1:40,k)=A2{j,k};
    end
end


file3=['CNNdataoutnew_',num2str(i+bottom-1),'.mat'] ;
A6x=load(file3) ;
A6=A6x.CNNdataout1 ;


A6{1,1}=zeros(80,40) ;
A6{1,2}=A6{1,1} ;
A6{1,3}=A6{1,1} ;


for j=1:45
    for k=1:3
multiin(im,j,1:80,1:40,k)=A6{j,k};
    end
end


im=im+1 ;
end
i=i+1;
end

a=size(multiin);
for iz=1:a(1,1)
    for iz1=1:a(1,2)
        for iz2=1:a(1,3)
            for iz3=1:a(1,4)
                for iz4=1:a(1,5)
                se=isnan(multiin(iz,iz1,iz2,iz3,iz4))  ;
                se1=isnan(multiout(iz,iz1,iz2,iz3,iz4))  ;
      if se==1
          multiin(iz,iz1,iz2,iz3,iz4)=0;
              
                    
      end
        if se1==1
          multiout(iz,iz1,iz2,iz3,iz4)=0;
              
                    
        end
                end
            end
        end
    end
end



save('traininput_v2.mat','multiin', '-v7.3');
save('trainoutput_v2.mat','multiout', '-v7.3');


