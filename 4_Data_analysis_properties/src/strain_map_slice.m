function [] = strainFinal05_er_yaoy(...
    disp3D,err3D,xv,yv,zv,dxyz,posDisp,plotPlanes)

pixelsize = 1/(xv(2) - xv(1));

fsize = 12;
n = 2;
densRange = [.02 .04]; % min to max for density

dispRange = [-2 2;
    -2 2;
    -2 2]*1;
strainRange = [-1 1]*.08*1;

if nargin < 8
    plotPlanes = [12:16]; % designed for HEA-16 particle
end
Ngrid = [length(plotPlanes) 1+3+6];

h = figure(1);
clf
set(h,'color','w','outerposition',[10  140  1500  900]);

% Determine z plane locations
zInd = zeros(length(plotPlanes),1);
for a0 = 1:length(plotPlanes)
    sub = posDisp(:,14) == plotPlanes(a0);
    zMean = mean(posDisp(sub,3));
    [~,ind] = min(abs(zMean-zv));
    zInd(a0) = ind;
end

dxyGrid = 1./Ngrid;
dxyGrid(1) = dxyGrid(1)*.9;
ha = cell(Ngrid);
for a0 = 1:Ngrid(1)
    for a1 = 1:Ngrid(2)
        ha{a0,a1} = axes('position',...
            [(a1-1)*dxyGrid(2) (a0-1)*dxyGrid(1) dxyGrid(2) dxyGrid(1)]);
%         axis equal off
        set(gca,'xtick',[],'ytick',[])
        box on
    end
end


% Legends
hleg = cell(4,1);
b = 0.04;
hleg{1} = axes('position',...
    [dxyGrid(2)*0 dxyGrid(1)*length(plotPlanes)+.01 dxyGrid(2)*1 .08]);
hleg{2} = axes('position',...
    [dxyGrid(2)*1+b/2 dxyGrid(1)*length(plotPlanes)+.08 dxyGrid(2)*3-b .015]);
hleg{3} = axes('position',...
    [dxyGrid(2)*4+b/2 dxyGrid(1)*length(plotPlanes)+.08 dxyGrid(2)*6-b .015]);
% hleg{4} = axes('position',...
%     [dxyGrid(2)*10+b/2 dxyGrid(1)*length(plotPlanes)+.08 dxyGrid(2)*3-b .015]);

% Axes
set(h,'currentaxes',hleg{1})
hold on
line([0 0 5],[5 0 0],'linewidth',2,'color','k')
line([-1 0 1],[-1 0 -1]*1.5+5,'linewidth',2,'color','k')
line([-1 0 -1]*1.5+5,[-1 0 1],'linewidth',2,'color','k')
text(0,10,'y','fontsize',fsize+2,'interpreter','latex',...
    'horizontalalign','center')
text(0,7,'[100]','fontsize',fsize,'interpreter','latex',...
    'horizontalalign','center')

text(8,2,'x','fontsize',fsize+2,'interpreter','latex',...
    'horizontalalign','center')
% text(8,-2,'[$0\overline{1}1$]','fontsize',fsize,'interpreter','latex',...
%     'horizontalalign','center')
text(8,-2,'[010]','fontsize',fsize,'interpreter','latex',...
    'horizontalalign','center')

text(13,7+2,'z','fontsize',fsize+2,'interpreter','latex',...
    'horizontalalign','center')
% text(13,7-2,'[011]','fontsize',fsize,'interpreter','latex',...
%     'horizontalalign','center')
text(13,7-2,'[001]','fontsize',fsize,'interpreter','latex',...
    'horizontalalign','center')

t = linspace(0,2*pi,128+1);
plot(8+cos(t),8+sin(t),'linewidth',2,'color','k')
patch(8+cos(t)/2,8+sin(t)/2,'k','edgecolor','none')

hold off
axis equal off
xlim([0 20]-4)
ylim([0 12]-1)
% box on

set(h,'currentaxes',hleg{2})
I = ind2rgb(1:256,cyanToRed);
% v = linspace(dispRange(1,1)*100,dispRange(1,2)*100,256);
v = linspace(dispRange(1,1),dispRange(1,2),256);
% ; v([1 end]) = [];
imagesc(I,'xdata',v)
% axis off
% xlabel('Lattice Displacement [pm]',...
%     'fontsize',fsize,'interpreter','latex')
xlabel('Lattice Displacement [pixel]',...
    'fontsize',fsize,'interpreter','latex')
set(gca,'ytick',[],'tickdir','out')
% xlim(dispRange(1,:)*100)
box off


set(h,'currentaxes',hleg{3})
I = ind2rgb(1:256,jet(256));
v = linspace(strainRange(1)*100,strainRange(2)*100,256);
% ; v([1 end]) = [];
imagesc(I,'xdata',v)
% axis off
xlabel('Infinitesimal Strain [\%]',...
    'fontsize',fsize,'interpreter','latex')
set(gca,'ytick',[],'tickdir','out')
xlim(100*strainRange)
box off

% Density plots
N = size(disp3D);
w = zeros(N(1),N(2),length(plotPlanes));
for a0 = 1:length(plotPlanes)
    disp_temp = disp3D(:,:,zInd(a0),1);
    disp_temp(disp_temp<0) = 0;
    w(:,:,a0) = (disp_temp-densRange(1)) ...
        / (densRange(2) - densRange(1));
end
w(w<0) = 0;
w(w>1) = 1;
mSize = 2;
for a0 = 1:length(plotPlanes)
    set(h,'currentaxes',ha{a0,1})
    imagesc(repmat(rot90(1-w(:,:,a0)/4,1),[1 1 3]),...
        'xdata',xv,'ydata',yv);
    hold on
    sub1 = posDisp(:,14) == plotPlanes(a0) & posDisp(:,7) == 1;
    sub2 = posDisp(:,14) == plotPlanes(a0) & posDisp(:,7) == 0;
    scatter(posDisp(sub1,1),-posDisp(sub1,2),...
        'marker','o','sizedata',mSize,...
        'markeredgecolor','none','markerfacecolor',[0 0 0])
    scatter(posDisp(sub2,1),-posDisp(sub2,2),...
        'marker','o','sizedata',mSize,...
        'markeredgecolor','none','markerfacecolor',[1 0 0])
    line([10,30]/pixelsize,[30,30]/pixelsize,'linewidth',5,'color','k');
    hold off
    axis equal off
    axis padded
end

% Displacement plots
for a0 = 1:length(plotPlanes)
    for a1 = 1:3
        I = min(max((rot90(disp3D(:,:,zInd(a0),a1+1)) ...
            - dispRange(a1,1)) ...
             / (dispRange(a1,2) - dispRange(a1,1)),0),1);
        I = ind2rgb(round(255*I)+1,cyanToRed) ...
            .* repmat(rot90(w(:,:,a0)),[1 1 3]) ...
            + repmat(rot90(1-w(:,:,a0)),[1 1 3]);
        
        set(h,'currentaxes',ha{a0,a1+1})
        
        I_intp = zeros(2^n * (size(I,1)-1) + 1,2^n * (size(I,2)-1) + 1,3);
        for i = 1:3
            I_temp = I(:,:,i);
            I_intp(:,:,i) = interp2(I_temp,n,'spline');
        end
        
        imagesc(I_intp)
        hold on
        hold off
        axis equal off
        axis padded
    end
end



% ii strain
for a0 = 1:length(plotPlanes)
    for a1 = 1:3
        switch a1
            case 1
                I = (circshift(disp3D(:,:,zInd(a0),a1+1),[-1 0]) ...
                    - circshift(disp3D(:,:,zInd(a0),a1+1),[1 0]))/(2*dxyz);
                E = sqrt((circshift(err3D(:,:,zInd(a0),a1+1),[-1 0])).^2 ...
                    + (circshift(err3D(:,:,zInd(a0),a1+1),[1 0])).^2)/(2*dxyz);
            case 2
                I = (circshift(disp3D(:,:,zInd(a0),a1+1),[0 -1]) ...
                    - circshift(disp3D(:,:,zInd(a0),a1+1),[0 1]))/(2*dxyz);
                E = sqrt((circshift(err3D(:,:,zInd(a0),a1+1),[0 -1])).^2 ...
                    + (circshift(err3D(:,:,zInd(a0),a1+1),[0 1])).^2)/(2*dxyz);
            case 3
                if a0 == length(plotPlanes)
                    I = (disp3D(:,:,zInd(a0),a1+1) ...
                        -disp3D(:,:,zInd(a0)-1,a1+1))/(dxyz);
                    E = sqrt((err3D(:,:,zInd(a0),a1+1)).^2 ...
                        +(err3D(:,:,zInd(a0)-1,a1+1)).^2)/(dxyz);
                else
                    I = (disp3D(:,:,zInd(a0)+1,a1+1) ...
                        -disp3D(:,:,zInd(a0)-1,a1+1))/(2*dxyz);
                    E = sqrt((err3D(:,:,zInd(a0)+1,a1+1)).^2 ...
                        +(err3D(:,:,zInd(a0)-1,a1+1)).^2)/2*(dxyz);
                end
        end

        E = ((E) - strainRange(1)) ...
            / (strainRange(2) - strainRange(1));
        E = ind2rgb(round(255*E')+1,jet(256)) ...
            .* repmat((w(:,:,a0))',[1 1 3]) ...
            + repmat((1-w(:,:,a0))',[1 1 3]);
        set(h,'currentaxes',ha{a0,a1+4})
        

        I = ((I) - strainRange(1)) ...
            / (strainRange(2) - strainRange(1));
        I = ind2rgb(round(255*I')+1,jet(256)) ...
            .* repmat((w(:,:,a0))',[1 1 3]) ...
            + repmat((1-w(:,:,a0))',[1 1 3]);
        set(h,'currentaxes',ha{a0,a1+4})

        I_intp = zeros(2^n * (size(I,1)-1) + 1,2^n * (size(I,2)-1) + 1,3);
        for i = 1:3
        I_intp(:,:,i) = interp2(I(:,:,i),n,'spline');
        end
        
        imagesc(I_intp);
        if a1 == 1
            line(4*[20,40],4*[20,20],'linewidth',5,'color','k');
        end
        axis xy
        hold on
        hold off
        axis equal off
        axis padded
    end
end


% ij strain
for a0 = 1:length(plotPlanes)
    for a1 = 1:3
        switch a1
            case 1
                dxy = (circshift(disp3D(:,:,zInd(a0),2),[0 -1]) ...
                    - circshift(disp3D(:,:,zInd(a0),2),[0 1]))/(2*dxyz);
                Exy = sqrt((circshift(err3D(:,:,zInd(a0),2),[0 -1])).^2 ...
                    + (circshift(err3D(:,:,zInd(a0),2),[0 1])).^2)/(2*dxyz);
                
                dyx = (circshift(disp3D(:,:,zInd(a0),3),[-1 0]) ...
                    - circshift(disp3D(:,:,zInd(a0),3),[1 0]))/(2*dxyz);
                Eyx = sqrt((circshift(err3D(:,:,zInd(a0),3),[-1 0])).^2 ...
                    + (circshift(err3D(:,:,zInd(a0),3),[1 0])).^2)/(2*dxyz);
                
                I = (dxy + dyx)/2;
                E = sqrt(Exy.^2 + Eyx.^2)/2;
            case 2
                if a0 == length(plotPlanes)
                    dxz = (circshift(disp3D(:,:,zInd(a0),2),[0 0]) ...
                        - circshift(disp3D(:,:,zInd(a0)-1,2),[0 0]))/(dxyz);
                    Exz = sqrt((circshift(err3D(:,:,zInd(a0),2),[0 0])).^2 ...
                        + (circshift(err3D(:,:,zInd(a0),2)-1,[0 0])).^2)/(dxyz);
                else
                    dxz = (circshift(disp3D(:,:,zInd(a0)+1,2),[0 0]) ...
                        - circshift(disp3D(:,:,zInd(a0)-1,2),[0 0]))/(2*dxyz);
                    Exz = sqrt((circshift(err3D(:,:,zInd(a0)+1,2),[0 0])).^2 ...
                        + (circshift(err3D(:,:,zInd(a0)-1,2),[0 0])).^2)/(2*dxyz);
                end
                dzx = (circshift(disp3D(:,:,zInd(a0),4),[-1 0]) ...
                    - circshift(disp3D(:,:,zInd(a0),4),[1 0]))/(2*dxyz);
                Ezx = sqrt((circshift(err3D(:,:,zInd(a0),4),[-1 0])).^2 ...
                    + (circshift(err3D(:,:,zInd(a0),4),[1 0])).^2)/(2*dxyz);
                
                
                I = (dxz + dzx)/2;
                E = sqrt(Exz.^2 + Ezx.^2)/2;
            case 3
                if a0 == length(plotPlanes)
                    dyz = (circshift(disp3D(:,:,zInd(a0),3),[0 0]) ...
                        - circshift(disp3D(:,:,zInd(a0)-1,3),[0 0]))/(1*dxyz);
                    Eyz = sqrt((circshift(err3D(:,:,zInd(a0),3),[0 0])).^2 ...
                        + (circshift(err3D(:,:,zInd(a0),3)-1,[0 0])).^2)/(1*dxyz);
                else
                    dyz = (circshift(disp3D(:,:,zInd(a0)+1,3),[0 0]) ...
                        - circshift(disp3D(:,:,zInd(a0)-1,3),[0 0]))/(2*dxyz);
                    Eyz = sqrt((circshift(err3D(:,:,zInd(a0)+1,3),[0 0])).^2 ...
                        + (circshift(err3D(:,:,zInd(a0)-1,3),[0 0])).^2)/(2*dxyz);
                end
                dzy = (circshift(disp3D(:,:,zInd(a0),4),[0 -1]) ...
                    - circshift(disp3D(:,:,zInd(a0),4),[0 1]))/(2*dxyz);
                Ezy = sqrt((circshift(err3D(:,:,zInd(a0),4),[0 -1])).^2 ...
                    + (circshift(err3D(:,:,zInd(a0),4),[0 1])).^2)/(2*dxyz);
                
                I = (dyz + dzy)/2;
                E = sqrt(Eyz.^2 + Ezy.^2)/2;
        end
        
        I = (rot90(I) - strainRange(1)) ...
            / (strainRange(2) - strainRange(1));
        I = ind2rgb(round(255*I)+1,jet(256)) ...
            .* repmat(rot90(w(:,:,a0)),[1 1 3]) ...
            + repmat(rot90(1-w(:,:,a0)),[1 1 3]);
        set(h,'currentaxes',ha{a0,a1+7})
        
        
        E = (rot90(E) - strainRange(1)) ...
            / (strainRange(2) - strainRange(1));
        E = ind2rgb(round(255*E)+1,jet(256)) ...
            .* repmat(rot90(w(:,:,a0)),[1 1 3]) ...
            + repmat(rot90(1-w(:,:,a0)),[1 1 3]);
        set(h,'currentaxes',ha{a0,a1+7})      

        I_intp = zeros(2^n * (size(I,1)-1) + 1,2^n * (size(I,2)-1) + 1,3);
        for i = 1:3
        I_intp(:,:,i) = interp2(I(:,:,i),n,'spline');
        end
        
        imagesc(I_intp)
        hold on
        hold off
        axis equal off
        axis padded
    end
end


% Main labels
w = .04;
for a0 = 1:9
    switch a0
        case 1
            t = '$\Delta x$';
        case 2
            t = '$\Delta y$';
        case 3
            t = '$\Delta z$';
        case 4
            t = '$\epsilon_{xx}$';
        case 5
            t = '$\epsilon_{yy}$';
        case 6
            t = '$\epsilon_{zz}$';
        case 7
            t = '$\epsilon_{xy}$';
        case 8
            t = '$\epsilon_{xz}$';
        case 9
            t = '$\epsilon_{yz}$';
            
    end
    
    annotation('textbox',....
        [dxyGrid(2)*(a0+.5)-w/2 .89 w .03],'String',t,...
        'FontSize',fsize+6,'interpreter','latex',...
        'horizontalalign','center','edgecolor','none');
end


end

function [cmap] = cyanToRed(Nvalues)
% Returns a color map with these entries:
% If user doesn't specify number of colours, assume 256 (8 bit)
if nargin == 0
    Nvalues = 256;
end
% coordinate system
x = linspace(0,1,Nvalues)';
% Red
r = -4*(x-1/2).^2 + 1;
r(x>=1/2) = 1;
% % Green
g = -4*(x-1/2).^2 + 1;
% Blue
b = -4*(1/2 - x).^2 + 1;
b(x<=1/2) = 1;
% Cyan map - set green to blue, or maybe average?
g = (b + g)/2;
% Assemble colormap
cmap = [r g b];
% Move to minor gray
cmap = 0.9*cmap + 0.0;
end