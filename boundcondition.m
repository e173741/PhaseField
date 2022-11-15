function [F,indi]=boundcondition(cord,t)

lengthx=0.06 ;
lengthy=0.12 ;
support=[0 0 lengthx 0 1 1] ;
loading=[0.01*lengthx lengthy lengthx*0.99 lengthy 0 -1000000] ;
tolerance=0.45 ;
total=size(cord,1) ;
F=zeros(total,1) ;
sf=size(support,1) ;
angle=zeros(sf,1);

 elementnumber=size(t,1) ;
elementsize=zeros(elementnumber,1) ;
for i=1:elementnumber
elementsize(i,1)=((cord(t(i,2),1)-cord(t(i,1),1))^2+(cord(t(i,2),2)-cord(t(i,1),2))^2)^0.5 ;
end

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
    end  
   if round(angle(i,1))==1 && support(i,2)<=cord(j,2) && cord(j,2)<=support(i,4) && support(i,1)-tolerance*elementsize(j,1)<cord(j,1) && support(i,3)+tolerance*elementsize(j,1)>cord(j,1)
    F(j,1)=support(i,5) ;
   end
   
end
end

cordsize=total ;
indi=zeros(cordsize,1);
 load=size(loading,1) ;

for i=1:load
    for j=1:cordsize
 
       distance=abs((loading(i,4)-loading(i,2))*cord(j,1)-(loading(i,3)-loading(i,1))*cord(j,2)+loading(i,3)*loading(i,2)-loading(i,4)*loading(i,1))/((loading(i,4)-loading(i,2))^2+(loading(i,3)-loading(i,1))^2)^0.5 ;
      if distance<elementsize(j,1)*2 && loading(i,1)<cord(j,1) && loading(i,3)>cord(j,1)
          indi(j,1)=1 ;
      end
    end
end



end