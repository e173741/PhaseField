total=6;
nodeinputfeature=cell(total,45) ; 
nodeoutputfeature=nodeinputfeature ;
edgeinputfeature=cell(total,1); 
sender=cell(total,1) ;
receiver=cell(total,1) ;
maximum=zeros(total,1) ;
minumum=maximum ;
nodemax=maximum ;
edgemax=maximum ;
edgefeatures=zeros(total,4) ;
maximum1=zeros(45,1) ;
minumum1=maximum1 ;
for iz=1:total
file1=['phasefieldsimulator_',num2str(iz),'.mat'] ;
A2x=load(file1) ;
A2=A2x.recordphase;
A2 = A2(~any(cellfun('isempty', A2), 2), :);
if size(A2,1)>=45
phase=A2x.phase;
t=A2x.t ;
input=A2x.CNNinput ;
cord=A2x.cord ;
maxnode=max(max(t)) ;
am1=size(t,1) ;
node=zeros(maxnode*6,2) ;
nodemax(iz,1)=maxnode ;

[F,indi]=boundcondition(cord,t) ;

for i=1:am1
      
        node(6*i-5,1)=t(i,1)-1 ;
        node(6*i-5,2)=t(i,2)-1 ;
        node(6*i-4,1)=t(i,1)-1 ;
        node(6*i-4,2)=t(i,3)-1 ;
        node(6*i-3,1)=t(i,2)-1 ;
        node(6*i-3,2)=t(i,1)-1 ;
        node(6*i-2,1)=t(i,2)-1 ;
        node(6*i-2,2)=t(i,3)-1 ;
        node(6*i-1,1)=t(i,3)-1 ;
        node(6*i-1,2)=t(i,1)-1 ;
        node(6*i,1)=t(i,3)-1 ;
        node(6*i,2)=t(i,2)-1 ;
        
end    
node1=unique(node, 'rows') ;
maxmax=size(node1,1);
sender{iz,1}=node1(:,1) ;
receiver{iz,1}=node1(:,2) ;
con=zeros(maxnode,6);
edgemax(iz,1)=maxmax;
k=1 ; tot=zeros(maxnode,1); 
for i=1:am1
    for j=1:3
     for im=1:am1
         for jm=1:3
    if  t(i,j)==t(im,jm)
            con(t(i,j),k)=im ;
            k=k+1 ;
    end
       
         end
     end
     tot(t(i,j),1)=k-1 ;
     k=1 ;
    end
end

a1=0 ; a2=0 ;a3=0 ;


edgefeatures(iz,1)=input{1,2}(1,1) ;
edgefeatures(iz,2)=input{1,4}(1,1) ;
edgefeatures(iz,3)=input{1,3}(1,1) ;
edgefeatures(iz,4)=input{1,1}(1,1) ;

nodestress=cell(45,1) ;
for im=1:45
for i=1:maxnode
    for j=1:tot(i,1)
a1=a1+A2{im,3}(con(i,j),1) ;
a2=a2+A2{im,3}(con(i,j),2) ;
a3=a3+A2{im,3}(con(i,j),3) ;
    end
nodestress{im,1}(i,1)=a1/tot(i,1) ;
nodestress{im,1}(i,2)=a2/tot(i,1) ;
nodestress{im,1}(i,3)=a3/tot(i,1) ;
a1=0; a2=0 ; a3=0 ;
end
end
 
 for i=1:45
     
   maximum1(i,1)=max(max(nodestress{i,1}(:,1),nodestress{i,1}(:,2))) ;
   minumum1(i,1)=min(min(nodestress{i,1}(:,1),nodestress{i,1}(:,2))) ;
 end
 
 maximum(iz,1)=max(maximum1) ;
 minumum(iz,1)=min(minumum1); 
 
 
damage=cell(45,1) ;
for im=1:45
for i=1:maxnode
damage{im,1}(i,1)=phase{im,1}(i,1) ;
end
end

for im=1:44
    for i=1:maxnode
nodeoutputfeature{iz,im}(i,1)=nodestress{im+1,1}(i,1) ;
nodeoutputfeature{iz,im}(i,2)=nodestress{im+1,1}(i,2) ; 
nodeoutputfeature{iz,im}(i,3)=nodestress{im+1,1}(i,3) ;
nodeoutputfeature{iz,im}(i,4)=damage{im+1,1}(i,1) ;

    end
end

for im=1:45
for i=1:maxnode
nodeinputfeature{iz,im}(i,1)=nodestress{im,1}(i,1) ;
nodeinputfeature{iz,im}(i,2)=nodestress{im,1}(i,2) ;
nodeinputfeature{iz,im}(i,3)=nodestress{im,1}(i,3) ;
nodeinputfeature{iz,im}(i,4)=damage{im,1}(i,1) ;
nodeinputfeature{iz,im}(i,5)=im*indi(i,1) ;
nodeinputfeature{iz,im}(i,6)=F(i,1) ;
end
end
end

end


for i=1:total
nodeoutputfeature{i,45}=nodeinputfeature{i,45} ;
end
NO=nodeoutputfeature ;
Ni=nodeinputfeature ;
nodeoutputfeature = nodeoutputfeature(~any(cellfun('isempty', nodeoutputfeature), 2), :);
nodeinputfeature = nodeinputfeature(~any(cellfun('isempty', nodeinputfeature), 2), :);
sender = sender(~any(cellfun('isempty', sender), 2), :);
receiver = receiver(~any(cellfun('isempty', receiver), 2), :);
total=size(nodeoutputfeature,1) ;


nodemax( ~any(nodemax,2), : ) = [];
edgemax( ~any(edgemax,2), : ) = [];
max1=max(maximum) ;
min1=min(minumum) ;

for iz=1:total
for im=1:44
    for i=1:nodemax(iz,1)
nodeoutputfeature{iz,im}(i,1)=(NO{iz,im}(i,1)-min1)/(max1-min1) ;
nodeoutputfeature{iz,im}(i,2)=(NO{iz,im}(i,2)-min1)/(max1-min1)  ; 
nodeoutputfeature{iz,im}(i,3)=(NO{iz,im}(i,3)-min1)/(max1-min1)  ; 

nodeinputfeature{iz,im}(i,1)=(Ni{iz,im}(i,1)-min1)/(max1-min1) ;
nodeinputfeature{iz,im}(i,2)=(Ni{iz,im}(i,2)-min1)/(max1-min1) ; 
nodeinputfeature{iz,im}(i,3)=(Ni{iz,im}(i,3)-min1)/(max1-min1) ;   

    end
end
end

edgefeatures( ~any(edgefeatures,2), : ) = [];

max2=max(edgefeatures(:,1)) ;
min2=min(edgefeatures(:,1)) ;
max3=max(edgefeatures(:,2)) ;
min3=min(edgefeatures(:,2)) ;
max4=max(edgefeatures(:,3)) ;
min4=min(edgefeatures(:,3)) ;
max5=max(edgefeatures(:,4)) ;
min5=min(edgefeatures(:,4)) ;

maxs(1,1)=max1 ; maxs(2,1)=max2; maxs(3,1)=max3;maxs(4,1)=max4;maxs(5,1)=max5 ;

maxs(1,2)=min1 ; maxs(2,2)=min2; maxs(3,2)=min3;maxs(4,2)=min4;maxs(5,2)=min5 ;


for iz=1:total
for i=1:nodemax(iz,1)
edgeinputfeature{iz,1}(i,1)=(edgefeatures(iz,1)-min2)/(max2-min2) ;
edgeinputfeature{iz,1}(i,2)=(edgefeatures(iz,2)-min3)/(max3-min3) ;
edgeinputfeature{iz,1}(i,3)=(edgefeatures(iz,3)-min4)/(max4-min4) ;
edgeinputfeature{iz,1}(i,4)=(edgefeatures(iz,4)-min5)/(max5-min5) ;
end
end
outputfile3=['nodeinput_',num2str(iz),'.mat'] ;
save(outputfile3,'nodeinputfeature','-v7.3')
outputfile4=['nodeoutput_',num2str(iz),'.mat'] ;
save(outputfile4,'nodeoutputfeature','-v7.3')
outputfile5=['edgeinput_',num2str(iz),'.mat'] ;
save(outputfile5,'edgeinputfeature','-v7.3')
outputfile6=['maxs_',num2str(iz),'.mat'] ;
save(outputfile6,'maxs','-v7.3')