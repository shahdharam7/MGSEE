function [m]=change_index(x,r,c)
for i=1:length(x)
    if(mod(x(i),c)==0)
        m(i,1)=floor(x(i)/c);
        m(i,2)=c;
    else
        m(i,1)=floor(x(i)/c)+1;
        m(i,2)=mod(x(i),c);
    end
end
end