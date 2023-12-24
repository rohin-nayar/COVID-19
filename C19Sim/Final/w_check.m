function [V] = w_check(XY,V,L,N,theta)
%This function ouputs a velocity magnitude that is equal in magnitude and has the opposite direction with an updated theta value to 
%randomise the change in direction.
%Inputs: XY - matrix of x and y coordinates, V - velocty magnitude matrix, L - plot dimension, N - number of individuals, theta - angles array 
for i=1:N
    %check if particles collide with the horizontal walls
    if XY(1,i)>=L-2 || XY(1,i)<=2
        theta = theta+0.5*pi*rand; %Theta values updated by adding a smaller angle
        V(1,i)=-V(1,i); %Velocity Direction Changed and component updated with new theta value
    end
    %check if particles collide with the vertical walls
    if XY(2,i)>=L-2 || XY(2,i)<=2
        theta = theta+0.5*pi*rand;
        V(2,i)=-V(2,i);
    end
end
