%Covid-19 Simulator Program
%Since the Covid-19 Outbreak simulations for tracking covid-19 
%in an environment have proven useful. One of the main purposes of 
%this simulation is to demonstrate one way in which a covid-19 outbreak
%or a similar virus could be simulated with simplifications and
%assumptions made.
%------------------------------------------------------------------------

%Housekeeping
clf
clearvars
close all

%Parameters and Variables
%User Input for the number of individuals to be included in the simulation
N = input('Enter number of individuals for simulation: ');%Number of Individuals
L =1000;       %Width and Height dimensions
dt = 10;       %Timestep
T=10*24*60*60; %Total Simulation time
nstep = T/dt;  %Total number of steps

%Allocating healthy status variables values 
healthy = 0; sick = 1; asymptomatic = 2; recovered = 3;

%Defining empty storage arrays
XY = zeros(2,N); %X coordinates row 1, Y coordinates row 2
A = zeros(3,N);  %Individual status row 1, Status time holder rows 2 & 3
A(1) = sick;     %Initial sick individual
info.timesec = zeros(nstep,1);   %Time structure (seconds)
info.time = zeros(nstep,1);      %Time structure (DD:HH:MM:SS)
info.healthy = zeros(nstep,1);   %Number of healthy individuals at any given time
info.asymptomatic = zeros(nstep,1);  %Number of asymptomatic individuals at any given time
info.sick = zeros(nstep,1);      %Number of sick individuals at any given time
info.recovered = zeros(nstep,1); %Number of recovered individuals at any given time
num.status = zeros(1,N);         %Total number of healthy individuals
num.healthy = zeros(1);        %Toal number of asymptomatic individuals
num.sick = zeros(1);           %Total number of infected individuals
num.asymptomatic = zeros(1);   %Total number of sick individuals
num.recovered = zeros(1);      %Total number of recovered individuals

%Allocating random initial individual positions
for i=1:N
    XY(1,i)=randi([0,L],1,1); %Random X coordinates
    XY(2,i)=randi([0,L],1,1); %Random Y coordinates
end

%Defining initial velocity components and theta values 
for i = 1:N
    vel = 0.1:0.001:0.2; %Velocity range
    vel_magnitude = vel(randi([1,numel(vel)],1,i)); %Generates random velocity magnitudes
    theta = 0:0.01:2*pi; %Theta range
    theta = theta(randi([1,numel(theta)],1,i)); %Theta values
    u_comp = vel_magnitude.*cos(theta);         %u velocity component
    v_comp = vel_magnitude.*sin(theta);         %v velocity component
end
V = [u_comp;v_comp]; %Defined Velocity Matrix

%Initial Plot
for i = 2:N
    if A(i) == 0
    plot(XY(1,i),XY(2,i),'o','MarkerEdgeColor', ... %Healthy Individuals Plotted
        [0 0.6 0.3],'MarkerFaceColor',[0 0.6 0.3]);
    end
    hold on
    plot(XY(1,1),XY(2,2),'o','MarkerEdgeColor', ... %Sick Individual Plotted
        [1,0.2,0.2],'MarkerFaceColor','[1,0.2,0.2]');
    title('Initial Plot')
end
box on;axis square %Box outline and axes defined 
xlabel('Width (m)','FontSize',11); ylabel('Height (m)','FontSize',11); %Axes lables
title('Location Plot','FontSize',12);
ax = gca;
exportgraphics(ax,'Initial_Plot.jpeg') %Saving initial plot
hold off

%Initialise History Line Plot 
figure(3)
axh=0;ayh=N-1;a2xa=0;a2ya=0; %Placeholder Variables to intialise history line plot
a3xs=0;a3ys=1;a4xr=0;a4yr=0;
an1 = animatedline(axh,ayh,'Color',[0 0.6 0.3],'LineWidth',2);  %Healthy indivduals line
an2 = animatedline(a2xa,a2ya,'Color',[1 0.5 0],'LineWidth',2);  %Asyptomatic individual ine
an3 = animatedline(a3xs,a3ys,'Color',[1 0.2 0.2],'LineWidth',2);%Sick individual line
an4 = animatedline(a4xr,a4yr,'Color',[0 0.5 1],'LineWidth',2);  %Recovered individual line

%Loop to Update Plot in time up to nstep
for T=1:nstep
    %Update Positions and Velocities with timestep
    XY(1,:) = XY(1,:) + dt*V(1,:); %X Positions Updated
    XY(2,:) = XY(2,:) + dt*V(2,:); %Y Positions Updated 
    V = w_check(XY,V,L,N,theta);   %Collision with walls Function
 
    for i = 1:N %Loop iterates through all individuals
        
        if A(1,i) == sick
            for j = 1:N
                    r = [XY(1,i)-XY(1,j),XY(2,i)-XY(2,j)]; %Finds array of (x;y)distances
                    dist(j) = abs(sqrt(sum(r.^2))); %Absolute distance calculation
            end
            
            for k = 1:length(dist)       %Loops through distance array for N individuals
                if i~=k && dist(k) <= 2  %Distance between any two individuals is less than 2 m condition
                    p = rand([1,1]);     %Generates random probability
                    if p >= 0.5          %If p >= 0.5 the individual in close contact is infected
                        A(1,k) = asymptomatic; %Individual status change to (infected) asymptomatic
                    end
                end
            end 
            
            if A(1,i) == asymptomatic %If an indvidual is infected and asymptomatic
                for m = 1:N
                    r  = [XY(1,i)-XY(1,m),XY(2,i)-XY(2,m)]; %Finds array of (x;y) distances
                    dist(m) = abs(sqrt(sum(r.^2))); %Absolute distance calculation
                    if (dist(m) <= 2)
                        p = rand([1,1]);     
                        if (p <= 0.3)      %If p <= 0.3 the infected individual will infect healthy individuals 
                            A(1,m) = sick; %Individual status change to (infected) sick
                        end
                    end
                end
            end
            
        end
        
        %Status Change after certain specified period of time
        if A(1,i) == asymptomatic %infected asymptomatic
            A(2,i) = A(2,i)+1; %counter begins to track time
        elseif A(1,i) == sick  %infected sick
            A(3,i) = A(3,i)+1; %counter begins to track time
        end
        if A(2,i) >= 17280 %after 2 days infected asymptomatic individual becomes sick
            A(1,i) = sick;
        end
        if A(3,i) >= 25920 %after 3 days sick individual becomes recovered
            A(1,i) = recovered;
        end
        
        %Allocating the total number of any indiviidual at any time to a structure 
        num.status = A(1,:); %Status of each individual structure
        if num.status(i) == 0
            num.healthy = length(find(A(1,:)==0)); %Ongoing total of healthy individuals
        elseif num.status(i) == 1
            num.sick = length(find(A(1,:)==1)); %Ongoing total of asymptomatic individuals
        elseif num.status(i) == 2
            num.asymptomatic = length(find(A(1,:)==2)); %Ongoing total of sick individuals
        elseif num.status(i) == 3
            num.recovered = length(find(A(1,:)==3));%Ongoing total of recovered individuals
        end
       
    end %End of loop for all indivduals
    
    %Storage for the number of any individual at any given timestep
    info(T).healthy = length(find(A(1,:)==0));      %Stores the total number of healthy individuals
    info(T).asymptomatic = length(find(A(1,:)==2)); %Stores the total number of asymptomatic individuals
    info(T).sick = length(find(A(1,:)==1));         %Stores the total number of sick individuals 
    info(T).recovered = length(find(A(1,:)==3));    %Stores the total number of recovered individuals
    time_temp = datevec(seconds(T));                %Stores time over nsteps (DD:HH:MM:SS)
    info(T).time = time_temp(3:6);
    time_handle = datevec(seconds(T*dt));    %Convert seconds to days and hours
    time = [num2str(time_handle(3)),' days ',num2str(time_handle(4)),' hours']; 
    info(T).timesec = T;                     %Stores time over nsteps (seconds)

    if mod(T,100) == 0 %After short individual the plot is updated
        
        figure(1) %Simulation Figure
        clf
        hold on
         for i=1:N
             if A(1,i) == healthy  %Healthy Individuals Plotted
                 plot(XY(1,i),XY(2,i),'o','color',[0 0.6 0.3], ...
                     'MarkerFaceColor',[0 0.6 0.3])
             elseif A(1,i) == asymptomatic  %Asymptomatic Individuals Plotted
                 plot(XY(1,i),XY(2,i),'o','color',[1 0.5 0], ...
                     'MarkerFaceColor',[1 0.5 0])
             elseif A(1,i) == sick %Sick Individuals Plotted
                 plot(XY(1,i),XY(2,i),'o','color',[1,0.2,0.2], ...
                     'MarkerFaceColor',[1,0.2,0.2])
             elseif A(1,i) == recovered %Immune Individuals Plotted
                 plot(XY(1,i),XY(2,i),'o','color',[0 0.5 1], ...
                     'MarkerFaceColor',[0 0.5 1])
             end           
         end
         %Plot Characteristics
         box on; axis square;
         xlabel('Width (m)','FontSize',11); 
         ylabel('Height (m)','FontSize',11);
         title(['Location Plot: ',time],'FontSize',12);
         pause(0.0001);
         
         %Saving Location Plot Code
         if T == 17600
            saveas(gcf,'Simulation Day 2','jpeg');
         elseif T == 34800
            saveas(gcf,'Simulation Day 4','jpeg');
         elseif T == 52000
            saveas(gcf,'Simulation Day 6','jpeg');
         end
    end
         
    if mod(T,200) ==0
        figure(2); %Bar Plot Figure
        Xbar = categorical({'Healthy','Asymptomatic','Sick','Recovered'});
        Xbar = reordercats(Xbar,{'Healthy','Asymptomatic','Sick','Recovered'});
        b_graph = bar(Xbar, [info(T).healthy info(T).asymptomatic info(T).sick info(T).recovered]);
        b_graph.FaceColor = 'flat';
        b_graph.CData(1,:) = [0 0.6 0.3]; %Healthy Colour
        b_graph.CData(2,:) = [1 0.5 0];   %Infected Colour
        b_graph.CData(3,:) = [1,0.2,0.2]; %Sick Colour
        b_graph.CData(4,:) = [0 0.5 1];   %Immune Colour
        
        %Plot Characteristics 
        ylabel('Number of Individuals','FontSize',11);    %Y Axis Label
        xlabel('Individual Health Status','FontSize',11); %X Axis Label
        title(['Health Status Tracker: ',time],'FontSize',12);
        ylim([0 N]);
        pause(0.0001);
        
        %Saving Bar Plot Code
        if T == 17600
            saveas(gcf,'Bar Graph Day 2','jpeg');
        elseif T == 34800
            saveas(gcf,'Bar Graph Day 4','jpeg');
        elseif T == 52000
            saveas(gcf,'Bar Graph Day 6','jpeg');
        end
        
        figure(3); %History Line Plot Figure
        timeline = info(T).timesec; 
        %Addpoint function adds the points in time for each individual line
        addpoints(an1,timeline/360,info(T).healthy);
        addpoints(an2,timeline/360,info(T).asymptomatic);
        addpoints(an3,timeline/360,info(T).sick);
        addpoints(an4,timeline/360,info(T).recovered);
        an1.LineWidth = 2; an2.LineWidth = 2; %Line Width 2
        an3.LineWidth = 2; an4.LineWidth = 2;
        
        %Plot Characteristics
        box on;
        title(['History of Outbreak: ',time],'FontSize',12);
        ylabel('Number of Individuals','FontSize',11);
        xlabel('Time (Hours)','FontSize',11);
        legend('Healthy','Asymptomatic','Sick','Recovered','Location','northeast','FontSize',8);
        pause(0.00001);
        
        %Saving History Plots Code
        if T == 17600
            saveas(gcf,'History Plot Day 2','jpeg');
        elseif T == 34800
            saveas(gcf,'History Plot Day 4','jpeg');
        elseif T == 52000
            saveas(gcf,'History Plot Day 6','jpeg');
        end 
    end
    
    %Daily Summaries Output to Command Window
    if mod(T,8640) == 0
        d = T/8640; %Time for each day
        diary on
        diary DailySummaries.txt
        disp(['Day: ', num2str(d)]);
        disp(['Number of healthy people: ', num2str(num.healthy)]);
        disp(['Number of asymptomatic people: ', num2str(num.asymptomatic)]);
        disp(['Nnumber of sick people: ', num2str(num.sick)]);
        disp(['Number of recovered people: ', num2str(num.recovered)]);
        diary off
    end
end