%This code creates a movie file titled 'JuliaSetVideo.avi'
%you can input different real and imaginary constants and it will play 
%through different imaginary constants added on to the real constant 
%(e.g. 0.5+x*i where x ranges from 0.2 to 0.95)
%made by Lance Schaefer

clear;clc;close all

prompt={'Enter the image resolution as a whole number: ',...
    'Enter the x-range as a vector: ','Enter the y-range as a vector: ',...
    'Enter the swept range of imaginary values as a vector: ',...
    'Enter the real part of the Julia-set constant: ',...
    'Enter the number of frames in the video as a whole number: ',...
    'Enter 1 to play an animation and save to a .avi file, enter 0 to just save to a file: '};
dlg_title='Input Parameters';
defaultans={'650','[-1 1]','[-1 1]','[0.35 0.85]','-0.4','200','0'};
parameters=inputdlg(prompt,dlg_title,1,defaultans);           

res=str2num(parameters{1}); %resolution of images
xrange=str2num(parameters{2});
yrange=str2num(parameters{3});
X=linspace(xrange(1),xrange(2),res);
Y=linspace(yrange(1),yrange(2),res);
Bounds=str2num(parameters{4});
lowerImBound=Bounds(1); %lower imaginary bound
upperImBound=Bounds(2); %upper imaginary bound
realConstant=str2num(parameters{5}); %real component of julia set constant
loops=str2num(parameters{6}); %number of different frames, note that the number of frames 
                              %in the movie will actually be double since it loops back to the beginning
check1=str2num(parameters{7}); %if 0: do not play an animation, just save to a video file
                               %if 1: do play an animation and save to a video file

count=0;

F(loops) = struct('cdata',[],'colormap',[]); %setup storage for each snapshot
for ii=linspace(lowerImBound,upperImBound,loops)
   count=count+1;
   [Z,out2]=Julia(res,[-1 1],[-1 1],(ii*1j+realConstant),0,0); %use the Julia set function made by Lance Schaefer
   pcolor(X,Y,Z);
   axis equal
   xlim([xrange(1) xrange(2)])
   ylim([yrange(1) xrange(2)])
   shading flat %set flat shading    
   colormap(hot); %set color scheme
   xlabel('Real Numbers')
   ylabel('Imaginary Numbers')
   title(sprintf('Julia Set Constant: %4.4f+%4.4fj',realConstant,ii))
   F(count)=getframe(gcf);
   
end
fps=32; %frames per second

L=length(F);
for ii=1:L-1 %this whole loop is to extend the length of F in order to reverse the movie back to the begining 
  F(L+ii)=F(L-ii);
end
  
v=VideoWriter('JuliaSetVideo.avi'); %write to a video file
open(v)
writeVideo(v,F)
close(v)

if check1==1; %play a movie now
    axes('Position',[0 0 1 1]);
    movie(F,-3,fps) %negative 3 means to play the movie back and forth 3 times
end
