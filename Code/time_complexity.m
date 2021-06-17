function [] = time_complexity()

n=15;
M=8;h=[2 3 4 5 6];
N=(n-1)^2;
L=3;%[2 3 4 5 6];

t1= n*h.*L.^(h-1)+(h-1)*n^3.*(L.^(h-1));
t2= h.*L.^h + h.*(n^3).*L.^h +(L.^(h-1)).*((N*M)+(M.*h)+(M.*h.^2)+(h.^3)); %
t=t2-t1;
display(t1)
display(t2)
figure(12);
subplot(2,2,3)
semilogy(h,t1/1000, 'rs-')%, 'LineWidth', 1);
hold on
semilogy(h,t2/1000, 'gd-')%, 'LineWidth', 1);
axis([0 7 0 10^5])
grid on
legend('n=5 PL-OMP','n=5 L-OMP','n=10 PL-OMP','n=10 L-OMP','n=15 PL-OMP','n=15 L-OMP')%,'n=15 1','n=15 2')
%title('M=8, L=3');
xlabel('h');
%xlabel('L')
ylabel('Time complexity (x 10^3)')
subplot(2,2,4)
semilogy(h,t/1000, 'bs-')%, 'LineWidth', 1);
hold on
legend('n=5','n=10','n=15')%,'n=15')
axis([0 7 0 10^5])
grid on
%title('M=8, L=3');
xlabel('h');
%xlabel('L')
ylabel('Time difference(x 10^3)')
end