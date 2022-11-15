
maxi=zeros(750,1) ;
nod1=cell(750,1) ;
nod2=nod1 ;
nod3=nod2 ;
for i=1:750
file1=['node_',num2str(i),'.mat'] ;
file2=['input_',num2str(i),'.mat'] ;
file3=['output_',num2str(i),'.mat'] ;
A1x=load(file1) ;
A1=A1x.node2;
A2x=load(file2) ;
A2=A2x.nodeinputfeature ;
A3x=load(file3);
A3=A3x.nodeoutputfeature ;

con=size(A1,1) ;

maxi(i,1)=max(max(A1)) ;
if i>1
nod1{i,1}=A1+maxi(i-1,1)*ones(con,2) ;
end
if i==1
 nod1{i,1}=A1 ;
end
nod2{i,1}=A2 ;
nod3{i,1}=A3 ;
end

nod1= nod1(~any(cellfun('isempty', nod1), 2), :);
nod2= nod2(~any(cellfun('isempty', nod2), 2), :);
nod3= nod3(~any(cellfun('isempty', nod3), 2), :);

nodetotal=cat(1,nod1{:,1}) ;
inputtotal=cat(1,nod2{:,1}) ;
outputtotal=cat(1,nod3{:,1}) ;

outputfile1=['nodetotal.mat'] ;
save(outputfile1,'nodetotal','-v7.3')
outputfile2=['inputtotal.mat'] ;
save(outputfile2,'inputtotal','-v7.3')
outputfile3=['outputtotal.mat'] ;
save(outputfile3,'outputtotal','-v7.3')




