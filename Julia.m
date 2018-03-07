function [ IterNum,MSplot ] = Julia(res,xrange,yrange,JuliaSetConstant,PlotImage,ExportImage)
%This function plots the Julia set. 
%   The required inputs are the resolution given by an integer number,
%   the xrange (e.g. [-2 1]) and yrange (e.g. [-1 1]).
%   The resolution indicates the number of steps taken along the x and y
%   axis.
%   The output of this function is a matrix of values from 1 to 100
%   corresponding to the number of iterations required for divergence past
%   2. After 100 iterations, the sequence is terminated and taken to
%   converge. The other outputs are plots of the actual Julia set using
%   the pcolor function. One plot has axis, the other doesn't. The one
%   without the axis is used to export an image file if you want one. The
%   other required inputs are the Julia Set Constant which determines how
%   the julia set will look and whether or not you want the function to
%   export a image. The values for Julia Set constant should be real or
%   imaginary. The value for PlotImage and ExportImage should be either 0 or 1.
check1=PlotImage;
check2=ExportImage;
xlow=xrange(1); xhigh=xrange(2);
ylow=yrange(1); yhigh=yrange(2);
xint=res; yint=xint;
[x,y]=meshgrid(linspace(xlow,xhigh,xint),linspace(yhigh,ylow,yint));
xy_grid=x+1j*y;
c_grid=zeros(size(xy_grid));
c_grid(:,:)=JuliaSetConstant;
IterNum=zeros(size(xy_grid));
z=xy_grid;
count=100; %100 is set to the maximum number of iterations
for ii=0:count
   z=z.^2+c_grid;
   IterNum(abs(z)>2 & IterNum==0)=ii; %rewrite over values in IterNum for 
   %which the corresponding element in z results in abs(z)>2. Only rewrite
   %over values in IterNum however if they have not already been assigned a
   %convergence/divergence values (i.e. if the IterNum element is still 0)
end   
IterNum(IterNum==0)=count; %set all values that did not diverge past 2 to 
%the maximum iteration munber
if check1==1
    figure(1) %first figure with axis labels for data
    realAxis=linspace(xlow,xhigh,xint);
    ImAxis=linspace(yhigh,ylow,yint);
    MSplot=pcolor(realAxis,ImAxis,IterNum);
    titleString=sprintf('Julia Set with C = %g+%g*i \n',real(c_grid(1,1)),...
                imag(c_grid(1,1)));
    title(titleString)
    xlabel('Real')
    ylabel('Imaginery')
    axis equal
    xlim([xlow xhigh])
    ylim([ylow yhigh])
    shading flat %set flat shading
    colormap(hot); %set color scheme
else
    MSplot=0;
end
if check2==1
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 50 50];
    fig.PaperPositionMode = 'auto';
    print('export','-dpng','-r0')
end

end
