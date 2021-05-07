function [x_se]=DE_to_SE(de,n,h,penultimate_id,a,b,c,d,e,A1,A2)
%% Change DE size to n^3
% Indices to embed 0 in DE n^2
%{
a=zeros(1,n^2);b=zeros(1,n^2);
c=zeros(1,n^2);d=zeros(1,n^2);
ii=0;
for j=1:n
    for i=1:n
        ii=ii+1;
        a(ii)= n^2*(j-1)+ n*(j-1)+i;
        b(ii)= n^2*j -i + 1;
        c(ii)= n^2*(j-1)+ j + n*(i-1);
        d(ii)= n^2*(j-1)+ (n+1)*(i-1)+1;
        
    end
end
e = n^2*(n-1)+(1:n^2);
%}
de_n3 = zeros(n^3,1);
k=1;
for i=1:n^3
     if (ismember(i,a)) || (ismember(i,b)) || (ismember(i,c)) || (ismember(i,d)) || (ismember(i,e))
         
         de_n3(i)=0;    
     else
         
         de_n3(i)= de(k);
         k=k+1;
     end
end

%% DE to SE 
% x_se = (A1|A2)*x_de 
%{
A1=zeros(n^2,n^3);A2=zeros(n^2,n^3);

for i=1:n^2
    
    j=(1:n)+n*(i-1);
    A1(i,j)=1;
end

for i=1:n^2
    g=ceil(i/n);
    k= mod(i,n);
    if k==0
        k=n;
    end
    j=g+(0:n-1)*n+(n^2*(k-1));
    A2(i,j)=1;
end
%}

x_se = (A1 | A2)*de_n3;

% For last edge in network in odd hop length 
if mod(h,2)==1
    x_se(n*(n-1)+penultimate_id)=1;
end
end
