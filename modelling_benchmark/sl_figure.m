clear
load plot.mat

et1_kgs   = -arrayfun(@(y) y.qin(1),bcof);
area1     = 1e-2; % hard coding
et1_mmday = et1_kgs/area1*86400;
time_day  = [bcof.tout]/3600/24;
time_nod_day= arrayfun(@(y) y.tout,nod) * c.dayPsec;

x_matrix = reshape(nod(1).terms{x_idx},[inp.nn1,inp.nn2]);
y_matrix = reshape(nod(1).terms{y_idx},[inp.nn1,inp.nn2]);

fig_pos.left   = 0.3;
fig_pos.bottom = 0.7;
fig_pos.length = 0.4;
fig_pos.height = 0.25;

%nt=10;
%a.fig=figure;
a.fs = 10;
a.lw = 2;
a.cz = 8;
fs   = 5; % sampling frequency
% Creating the movie : quality = 100%, no compression is used and the
% timing is adjusted to the sampling frequency of the time interval
qt = 100;
%A number from 0 through 100. Higher quality numbers result in higher video quality and larger file sizes
a.fig = figure;
%set(gcf,'Units','normalized', 'WindowStyle','docked','OuterPosition',[0 0 1 1]);  % maximize the plotting figure
set(gcf,'Units','normalized', 'OuterPosition',[0 0 1 1]);  % maximize the plotting figure
mov           =  VideoWriter('linux.avi');% avifile('pvc1.avi','quality',qt,'compression','indeo5','fps',fs);
mov.FrameRate = 5;
mov.Quality=qt;
open(mov);

for nt=1:round((length(nod)-1)/50):length(nod)-1

    fprintf('plotting the %i out of %i results,\n', nt, length(nod));
    number_of_negative_c=sum(nod(nt).terms{c_idx}<0);
    [max_c,max_c_idx]=max(nod(nt).terms{c_idx});
    fprintf('    max concentration is %f at location  %i, x= %f , y= %f \n',max_c,max_c_idx,  x_nod_mtx(max_c_idx),  y_nod_mtx(max_c_idx));
    if number_of_negative_c>0
          fprintf('    %i negative concentration\n', length(number_of_negative_c));
    end

    %% -------------  sub 1 ET over time  --------------
    a.sub1=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom,...
          fig_pos.length,fig_pos.height]);
    a.plot1=plot(time_day(1:nt),et1_mmday(1:nt),...
             'k-','linewidth',a.lw);hold on
    %a.plot1=plot(eslab(1,:),eslab(2,:),'cx','linewidth',a.lw);
    get(gca,'xtick');
    set(gca,'fontsize',10);
    txt=sprintf('Result at day %.2f',nod(nt).tout*c.dayPsec);
    title(txt);
    xlabel('Time (day)','FontSize',a.fs);
    ylabel('Evt (mm/day)','FontSize',a.fs);
    axis([-0.1 time_day(end) -0.5 21])
    %title('ead profile');

%% -------- contour plot on Saturation ---------
    a.sub5=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-(fig_pos.height+0.07),...
          fig_pos.length+0.1,fig_pos.height]);
		  
    % write pressure and conc in matrix form.
    s_matrix  = reshape(nod(nt).terms{s_idx},[inp.nn1,inp.nn2]);
    contourf(x_matrix,y_matrix,s_matrix,'EdgeColor','none');hold on

%    scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
    color = jet;
    colormap(gca,color);
	caxis([0 1])
    cbsal = colorbar;
    cbsal.Label.String = 'Saturation (-)';
    get(gca,'xtick');
    set(gca,'fontsize',10);
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    title('Saturation')


%% -------- contour plot on concentration ---------
    a.sub6=subplot('position'...
         ,[fig_pos.left,fig_pos.bottom-2*(fig_pos.height+0.07),...
          fig_pos.length+0.1,fig_pos.height]);
    
    % write pressure and conc in matrix form.
    c_matrix = reshape(nod(nt).terms{c_idx},[inp.nn1,inp.nn2]);    
    contourf(x_matrix,y_matrix,c_matrix,'EdgeColor','none');hold on
	
%	scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
    color = jet;
    colormap(gca,color);
	caxis([0 0.264])
    cbsal = colorbar;
    cbsal.Label.String = 'Concentration (-)';
    get(gca,'xtick');
    set(gca,'fontsize',10);
    title('Concentration (kg/kg)');
    xlabel('x (m)','FontSize',a.fs);
    ylabel('Elevation (m)','FontSize',a.fs);
    %axis([10, 40,9,10])
	
	F = getframe(gcf); % save the current figure
    writeVideo(mov,F);% add it as the next frame of the movie

end

close(mov);
close(a.fig);
%figure
% [csal,hcsal] = contourf(x_matrix,y_matrix,p_mtx);
    % hold on
    % set(hcsal,'EdgeColor','none');
    % color = jet;
    % colormap(gca,color);
    % cbsal = colorbar;
    % cbsal.Label.String = 'Pressure (-)';
% scatter (a1(1,(1:f3(1)),1),a1(2,(1:f3(1)),1),2,'filled','w')
% savefig('pressure_contour.fig')						 

figure
[csal,hcsal] = contourf(x_matrix,y_matrix,c_matrix);
    hold on
    set(hcsal,'EdgeColor','none');
    color = jet;
    colormap(gca,color);
    cbsal = colorbar;
	caxis([0 0.246])
    cbsal.Label.String = 'Concentration (-)';
scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
savefig('concentration_contour.fig')						 

figure
[csal,hcsal] = contourf(x_matrix,y_matrix,s_matrix);
    hold on
    set(hcsal,'EdgeColor','none');
    color = jet;
    color = flipud(color);
    colormap(gca,color);
    cbsal = colorbar;
	caxis([0 1])
    cbsal.Label.String = 'Saturation (-)';
scatter(nod(1).terms{x_idx},nod(1).terms{y_idx},2,'filled','w');
savefig('saturation_contour.fig')	
