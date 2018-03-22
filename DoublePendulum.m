%Written by Lance M. Schaefer
%Written in March 2018

%This code simulates the motion of a double pendulum. It does this by
%numerically solving the equations of motion of a double pendulum. 
%You can solve these eqautions yourself either with the transport theorem
%or by formulating a Lagrangian and differentiating.
%The numerical integration technique is 4th order Runge-Kutta (RK4). 
%The user can input a number of different parameters to control the initial
%conditions of the system. These parameters include the initial angles, the
%initial angular velocities, the mass of each ball, the length of each rod,
%and the effective acceleration due to gravity. As this code runs, it will
%play through an animation of the motion as it actually plots each
%successive figure window. Once it is done running, it will save the
%animation to a .avi file titled DoublePendulum.avi. The user can further
%choose to replay a live animation at any number of frames per second 
%immeditaley once the code terminates

%A few things to note:
% 1) All angles are measured with respect to the vertical axis (i.e both
% the angle between the first rod and ball and the angle between the second
% rod and ball are measured with respect to the vertical axis)
% 2) The coordinate system obeys the right-hand-rule.
% 3) The positve x-axis points to the right and the positive y-axis
% points upward
% 4) This code was written in MATLab version 2016a. If you try to run it in
% a newer version, some functionality may be outdated.
% 5) There is no dissipation of energy in this simulation.
% 6) Remember that RK4 is not perfect. For long-time solutions small errors
% accumulate; and, for crazy initial conditions (typically ones with large 
% initial angular veclocities), you may see divergent solutions with
% infinities popping up


%Do not try to enter values directly in the code; instead, run the code and
%input values into the dialogue box.


%Beginning of code
clear; clc; close all

%initialize the dialogue box
prompt={'Enter the time interval in seconds. Enter it as a vector: ',...
    'Enter the smallest time increment in seconds: ',...
    'Enter the mass of ball 1 and 2 in kilograms. Enter them a as a vector: ',...
    'Enter the gravitational constant in m/s^2: ',...
    'Enter the length of each rod in meters. Enter them as a vector: ',...
    'Enter the initial angles between each rod and the vertical axis in radians. Enter them as a vector: ',...
    'Enter the initial angular velocities of each ball in m/s. Enter them as a vector: ',...
    'Enter the number of frames in the video as a whole number: ',...
    'Enter 1 to play an animation once the code terminates and save to a .avi file; or, enter 0 to just save to a .avi file: '};
dlg_title='Input Parameters';
defaultans={'[0 6]','0.01','[1 1]','9.81','[1 1]','[pi/4 pi]','[0 0]','30','0'};
parameters=inputdlg(prompt,dlg_title,1,defaultans);     

%define to what variables the entries in the dialogue box correspond
time_interval=str2num(parameters{1}); %in seconds
a=time_interval(1); %seconds
b=time_interval(2); %seconds
h=str2num(parameters{2}); %seconds

masses=str2num(parameters{3}); %in kg
m1=masses(1); %kg
m2=masses(2); %kg
g=str2num(parameters{4}); %m/s^2

lengths=str2num(parameters{5}); %in m
L1=lengths(1); %m
L2=lengths(2); %m

angles=str2num(parameters{6});
theta1=angles(1); %radians initial condition at t=0;
theta2=angles(2); %radians initial condition at t=0;

angular_freq=str2num(parameters{7});
omega1=angular_freq(1); %rad/s initial condition at t=0;
omega2=angular_freq(2); %rad/s initial condition at t=0;

%express the four equations of motion governing this double pendulum.

%eqn1
theta1_prime=@(t,theta1,omega1) omega1;

%eqn2
theta2_prime=@(t,theta2,omega2) omega2;

%eqn3
omega1_prime=@(t,theta1,theta2,omega1,omega2) (-g*(2*m1+m2)*sin(theta1)-m2*g*sin(theta1-2*theta2)-...
               2*sin(theta1-theta2)*m2*((omega2^2)*L2+(omega1^2)*L1*cos(theta1-theta2)))/...
               (L1*(2*m1+m2- m2*cos(2*theta1-2*theta2)));           

%eqn4           
omega2_prime=@(t,theta1,theta2,omega1,omega2) (2*sin(theta1-theta2)*((omega1^2)*L1*(m1 + m2)+...
              g*(m1+m2)*cos(theta1) + (omega2^2)*L2*m2*cos(theta1-theta2)))/...
                (L2*(2*m1 + m2-m2*cos(2*theta1-2*theta2)));        
      
            

%begin numerical integration via RK4
DivCheck=0;
inc=1;
for ii=a:h:b-h
    
    time_int(inc)=ii;
    
    %eqn1: solve for theta1 by incrementing time and theta1
    k1=double(theta1_prime(ii,theta1(inc),omega1(inc)));
    k2=double(theta1_prime(ii+h/2,theta1(inc)+((1/2)*(k1*h)),omega1(inc)));
    k3=double(theta1_prime(ii+h/2,theta1(inc)+((1/2)*(k2*h)),omega1(inc)));
    k4=double(theta1_prime(ii+h,theta1(inc)+k3*h,omega1(inc)));
    theta1(inc+1)=theta1(inc)+(1/6)*(k1+2*k2+2*k3+k4)*h;
    clear k1 k2 k3 k4
    
    %eqn2: solve for theta2 by incrementing time and theta2
    k1=double(theta2_prime(ii,theta2(inc),omega2(inc)));
    k2=double(theta2_prime(ii+h/2,theta2(inc)+((1/2)*(k1*h)),omega2(inc)));
    k3=double(theta2_prime(ii+h/2,theta2(inc)+((1/2)*(k2*h)),omega2(inc)));
    k4=double(theta2_prime(ii+h,theta2(inc)+k3*h,omega2(inc)));
    theta2(inc+1)=theta2(inc)+(1/6)*(k1+2*k2+2*k3+k4)*h;
    clear k1 k2 k3 k4
    
    %eqn3: solve for omega1 by incrementing time and omega1
    k1=double(omega1_prime(ii,theta1(inc),theta2(inc),omega1(inc),omega2(inc)));
    k2=double(omega1_prime(ii+h/2,theta1(inc),theta2(inc),omega1(inc)+((1/2)*(k1*h)),omega2(inc)));
    k3=double(omega1_prime(ii+h/2,theta1(inc),theta2(inc),omega1(inc)+((1/2)*(k2*h)),omega2(inc)));
    k4=double(omega1_prime(ii+h,theta1(inc),theta2(inc),omega1(inc)+k3*h,omega2(inc)));
    omega1(inc+1)=omega1(inc)+(1/6)*(k1+2*k2+2*k3+k4)*h;
    clear k1 k2 k3 k4
    
    %eqn4: solve for omega2 by incrementing time and omega2
    k1=double(omega2_prime(ii,theta1(inc),theta2(inc),omega1(inc),omega2(inc)));
    k2=double(omega2_prime(ii+h/2,theta1(inc),theta2(inc),omega1(inc),omega2(inc)+((1/2)*(k1*h))));
    k3=double(omega2_prime(ii+h/2,theta1(inc),theta2(inc),omega1(inc),omega2(inc)+((1/2)*(k2*h))));
    k4=double(omega2_prime(ii+h,theta1(inc),theta2(inc),omega1(inc),omega2(inc)+k3*h));
    omega2(inc+1)=omega2(inc)+(1/6)*(k1+2*k2+2*k3+k4)*h;
    clear k1 k2 k3 k4
    
    inc=inc+1;
end
time_int(inc)=b; %add the last time step to the time vector so it matches the other vectors in length

%position of mass 1:
x1=L1.*sin(theta1);
y1=-L1.*cos(theta1);

%position of mass 2:
x2=x1+L2.*sin(theta2);
y2=y1-L2.*cos(theta2);


%lines between origin (0,0) and mass 1, and between mass 1 and mass 2:
m1=(y1-0)./(x1-0); %slope of all the lines between the origin and ball 1
m2=(y2-y1)./(x2-x1); %slope of all the lines between the ball 1 and ball 2

sidelim=(L1+L2)*1.1; %limit the size of the figure window to 1.1 times the mass length of both rods fully extended

figure(1) %initialize the figure window
F(length(x1)) = struct('cdata',[],'colormap',[]); %setup storage for each figure's snapshot
for ii=1:length(x1)
    p0=plot(0,0,'ks'); %plot the anchor point
    p0.MarkerFaceColor='k';
    p0.MarkerSize=9;
    hold on
    p1=plot(x1(ii),y1(ii),'bo'); %plot ball 1
    p1.MarkerFaceColor='b';
    p1.MarkerSize=10;
    p2=plot(x2(ii),y2(ii),'ro'); %plot ball 2
    p2.MarkerFaceColor='r';
    p2.MarkerSize=10;
    xline1=linspace(0,x1(ii),2); %two x-coords., 
    yline1=m1(ii).*xline1; %two y-coords.
    xline2=linspace(x1(ii),x2(ii),2);
    if m2(ii)>10^14 || m2(ii)<-10^14 %ensure that we don't have a massive slope for a near veritcle rod 2 orientation
        yline2(1)=yline1(2);
        yline2(2)=yline1(2)+L2;
    else
        yintercept=y2(ii)-m2(ii)*x2(ii);
        yline2=m2(ii).*(xline2)+yintercept;
    end
    p3=plot(xline1,yline1,'k-'); %plot rod 1
    p3.LineWidth=1.5;
    p4=plot(xline2,yline2,'k-'); %plot rod 2
    p4.LineWidth=1.5;
    xlim([-sidelim, sidelim]);
    ylim([-sidelim, sidelim]);
    grid on
    xlabel('meters')
    ylabel('meters')
    title(sprintf('theta1=%4.4f degrees, theta2=%4.4f degrees\n omega1=%4.4f rad/s, omega2=%4.4f rad/s',...
         theta1(ii)*180/pi,theta2(ii)*180/pi,omega1(ii),omega2(ii)))
    hold off %make sure to turn the hold off otherwise each successive figure will save the previous position's orientation
    F(ii)=getframe(gcf); %save the figure window to store as a frame in the video
    if theta1(ii)>10^20 || theta2(ii)>10^20
        warndlg('Solution is divergent, simulation terminated','Divergence')
        DivCheck=1;
        break;
    end
        
end

if DivCheck==0 %only if the solution doesn't diverge finish the code
    %write to an .avi file
    v=VideoWriter('DoublePendulum.avi'); %write to a video file
    open(v)
    writeVideo(v,F)
    close(v)

    %if the user wanted to play an imediate animation, play it
    check1=str2num(parameters{9});
    if check1==1
        fps=str2num(parameters{8}); %frames per second
        axes('Position',[0 0 1 1]);
        movie(F,1,fps) %play the movie once
    end
end
