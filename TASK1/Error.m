%TASK 1.2
clear all
close all
nelms  = [4, 6, 8, 10];
error_array = zeros(size(nelms,2),1);

figure(2)
hold on
for ielems=1:size(nelms,2)
    [coords,d, Fe] = F1D(nelms(ielems));
    plot(coords, d)
    error_array(ielems) = Fe;
end
hold off
figure(3)
hold on 
plot(log(1./nelms), log(error_array));
hold off