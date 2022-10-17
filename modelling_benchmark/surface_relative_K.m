%calculate and plor Rayleigh number along soil surface
%can be used for evenly mesh Y of surface layer 
c=ConstantObj();

% layer_thickness_m = 0.01;%the thickness of the surface layer

density_variation = 700; %for sodium chloride d(tho_water)/d(concentration)
porosity          = inp.por(1);
dynamic_viscosity = 1e-3;
diffusivity       = inp.sigmaw;
permeability      = inp.pmax(1);

day_output = [7];%set days need to be plotted			 																  					 																  		
timestep_output = round(day_output*c.secPday/inp.nprint/inp.scalt);

for nt = timestep_output
time_step = nt;
time_day  = [bcof.tout]/3600/24;%second to day
x_matrix = reshape(nod(1).terms{x_idx},[inp.nn1,inp.nn2]);%inp.nn2 is number of nodes in x direction 
y_matrix = reshape(nod(1).terms{y_idx},[inp.nn1,inp.nn2]);
x_ele_matrix = reshape(ele(1).terms{xele_idx},[inp.nn1-1,inp.nn2-1]);
y_ele_matrix = reshape(ele(1).terms{yele_idx},[inp.nn1-1,inp.nn2-1]);

kr_matrix = reshape(ele(time_step).terms{kr_idx},[inp.nn1-1,inp.nn2-1]);

surface_kr= kr_matrix(end,:);
surface_permeability = surface_kr*permeability;
surface_conductivity = surface_permeability*9800*1000;

%plot
a.fs = 10;
a.lw = 1.75; %line width
a.cz = 8; %the size of the marker

figure

semilogy (x_ele_matrix(1,:),surface_conductivity,'k-','linewidth',a.lw);hold off
    grid on
	grid minor
% 	ax = gca;
	ax1.GridAlpha = 0.4;
	ax1.MinorGridAlpha = 0.5;	
    ylabel('Hydraulic conductivity (m/s)','FontSize',a.fs);
    xlim ([0 x_matrix(1,end)])

	figure_name=sprintf('Relative_K_day_%.2f.fig',nod(nt).tout*c.dayPsec);
	saveas(gcf,figure_name)

end


