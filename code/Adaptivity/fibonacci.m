function [k,list]= fibonacci(n,list)
if n==1
    k=1;
   list  = [1];
elseif n==2
    k=1;
    list = [1,1];
    
else
    k= fibonacci(n-1,list) +fibonacci(n-2,list);
end
    
end