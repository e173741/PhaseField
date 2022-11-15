
for iz=1:1000
file1=['phasefieldsimulator_',num2str(i),'.mat'] ;
A2x=load(file1) ;
A2=A2x.recordphase;
phase=A2x.phase;
t=A2x.t ;
input=A2x.CNNinput ;
cord=A2x.cord ;
maxnode=max(max(t)) ;
am1=size(t,1) ;
node=zeros(maxnode*6,2) ;


[F,indi]=boundcondition(cord,t) ;

for i=1:am1
      
        node(6*i-5,1)=t(i,1) ;
        node(6*i-5,2)=t(i,2) ;
        node(6*i-4,1)=t(i,1) ;
        node(6*i-4,2)=t(i,3) ;
        node(6*i-3,1)=t(i,2) ;
        node(6*i-3,2)=t(i,1) ;
        node(6*i-2,1)=t(i,2) ;
        node(6*i-2,2)=t(i,3) ;
        node(6*i-1,1)=t(i,3) ;
        node(6*i-1,2)=t(i,1) ;
        node(6*i,1)=t(i,3) ;
        node(6*i,2)=t(i,2) ;
        
end    
node1=unique(node(:,1:2), 'rows') ;
maxmax=size(node1,1);
nod=cell(4,1) ;
for i=1:4
nod{i,1}=node1+(i-1)*(maxnode)*ones(maxmax,2) ;
end

node2=cat(1,nod{:,1}) ;
nodeinputfeature=cell(4,1) ; 
nodeoutputfeature=zeros(am1,3) ;
con=zeros(maxnode,6); 
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

a1=0 ; a2=0 ;

for im=1:4
for i=1:maxnode
nodeinputfeature{im,1}(i,1)=input{1,2}(1,1)*10^-9 ;
nodeinputfeature{im,1}(i,2)=input{1,4}(1,1) ;
nodeinputfeature{im,1}(i,3)=input{1,3}(1,1)/input{1,1}(1,1) ;
nodeinputfeature{im,1}(i,4)=phase{1,1}(i,1) ;
nodeinputfeature{im,1}(i,5)=im*indi(i,1) ;
nodeinputfeature{im,1}(i,6)=0 ;
nodeinputfeature{im,1}(i,7)=F(i,1) ;
end
end


nodeinputfeature1=cat(1,nodeinputfeature{:,1}) ;

nodestress=cell(4,1) ; k=1;
for im=1:10:40
for i=1:maxnode
    for j=1:tot(i,1)
a1=a1+A2{im,3}(con(i,j),1) ;
a2=a2+A2{im,3}(con(i,j),2) ;
    end
nodestress{k,1}(i,1)=a1/tot(i,1) ;
nodestress{k,1}(i,2)=a2/tot(i,1) ;
a1=0; a2=0 ;
end
k=k+1 ;
end
damage=cell(4,1) ;
k=1;
for im=1:10:40
for i=1:maxnode
damage{k,1}(i,1)=phase{im,1}(i,1) ;
end
k=k+1 ;
end
nodestress1 = cat(1, nodestress{:,1}) ;
maxnode1=size(nodestress1,1) ;
damage1=cat(1,damage{:,1});
for i=1:maxnode1
nodeoutputfeature(i,1)=nodestress1(i,1)*10^-6 ;
nodeoutputfeature(i,2)=nodestress1(i,2)*10^-6 ; 
nodeoutputfeature(i,3)=damage1(i,1) ;
end

outputfile1=['node_',num2str(i),'.mat'] ;
save(outputfile1,'node2','-v7.3')
outputfile2=['input_',num2str(i),'.mat'] ;
save(outputfile2,'nodeinputfeature1','-v7.3')
outputfile3=['output_',num2str(i),'.mat'] ;
save(outputfile3,'nodeoutputfeature','-v7.3')

end
