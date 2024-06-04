function val= tdoaf(x,xsensor,y)

Nsensor=size(xsensor,2);

val=0;
for i=2:Nsensor
    val=val+(y(i-1)-(norm(x-xsensor(i,:))-norm(x-xsensor(1,:))))^2;
end