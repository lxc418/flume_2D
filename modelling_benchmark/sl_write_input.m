clear all
close all
c=ConstantObj();  % all the constants, we suggest to use c.g rather than 9.81 in the script to enhance readerbility;

%% SUTRA.fil

fil=filObj('read_from_file','no');

% the default value for 'fid' is -1, meaning, this will not write to SUTRA.fil if export
fil.terms.('vapinp').('fname')  = 'FLUME.vapinp' ; fil.terms.('vapinp').('fid')   = 50;
fil.terms.('ics').('fname')  = 'FLUME.ics' ; fil.terms.('ics').('fid')   = 55;
fil.terms.('lst').('fname')  = 'FLUME.lst' ; fil.terms.('lst').('fid')   = 60;
fil.terms.('rst').('fname')  = 'FLUME.rst' ; fil.terms.('rst').('fid')   = 66;
fil.terms.('nod').('fname')  = 'FLUME.nod' ; fil.terms.('nod').('fid')   = 30;
fil.terms.('ele').('fname')  = 'FLUME.ele' ; fil.terms.('ele').('fid')   = 40;
fil.terms.('obs').('fname')  = 'FLUME.obs' ; fil.terms.('obs').('fid')   = 70;
fil.terms.('smy').('fname')  = 'FLUME.smy' ; fil.terms.('smy').('fid')   = 80;
fil.terms.('bcof').('fname') = 'FLUME.bcof'; fil.terms.('bcof').('fid')  = 91;
fil.terms.('bcos').('fname') = 'FLUME.bcos'; fil.terms.('bcos').('fid')  = 93;
fil.terms.('bcop').('fname') = 'FLUME.bcop'; fil.terms.('bcop').('fid')  = 92;
fil.terms.('bcou').('fname') = 'FLUME.bcou'; fil.terms.('bcou').('fid')  = 94;


fil.export_to_file();

c_saltwater_kgPkg          = 0.035;
c_freshwater_kgPkg         = 0.0001;
initial_head_aquifer_m     = -2.8;
%initial_pond_water_depth_m = 0.5;
%permeability_silt_m2       = 5.38e-15;
permeability_sand_m2       = 5.67e-11;
pond_radius_m              = 50.1;

%% vapinp file
dx      = 1.0;
dy      = 1;  % vertical direction
dz      = 1;

x_array = 0:dx:100;
y_array = -10:dy:1;
nx      = length(x_array);
ny      = length(y_array);



x0_1=0;
x0_2=1.2;
nx_1 = 60;
x_array = x0_1:(x0_2-x0_1)/nx_1:x0_2;

y0_1=0.;
y0_2=0.;

y1_1=0.1;
y1_2=0.1;
ny_1=10;

y2_1=0.2;
y2_2=0.2;
ny_2=20;

y3_1=0.4;
y3_2=0.3;
ny_3=30;

ny=ny_1+ny_2+ny_3;

if y0_1==y0_2
y_array0(1:nx_1+1) = y0_1;
else
y_array0 = y0_1:(y0_2-y0_1)/nx_1:y0_2;
end
if y1_1==y1_2
y_array1(1:nx_1+1) = y1_1;
else
y_array1 = y1_1:(y1_2-y1_1)/nx_1:y1_2;
end
if y2_1==y2_2
y_array2(1:nx_1+1) = y2_1;
else
y_array2 = y2_1:(y2_2-y2_1)/nx_1:y2_2;
end
if y3_1==y3_2
y_array3(1:nx_1+1) = y3_1;
else
y_array3 = y3_1:(y3_2-y3_1)/nx_1:y3_2;
end

y1_interval=(y_array1(id)-y_array0(id))/ny_1
y2_interval=(y_array2(id)-y_array1(id))/ny_2
y3_interval=(y_array3(id)-y_array2(id))/ny_3

for id=1:nx_1+1

x_nod ( (id-1)*ny+1: id*(ny+1) )= x_array(id);
y_nod ( (id-1)*ny+1: id*(ny+1) )= [y_array0(id):y1_interval:y_array1(id) y_array1(id)+y2_interval:y2_interval:y_array2(id) y_array2(id)+y3_interval:y3_interval:y_array3(id)];

end



nx      = length(x_array);



y_0=0;

ny_1 = 10;
y2_1=0.1;
y2_2=0.1;

ny_2 = 20;
y3_1=0.2;
y3_2=0.2;

ny_3 = 30;
y4_1=0.4;
y4_2=0.3;

y_array = [y_0:(y_1-y_0)/ny_1:x_1,];


nex  = length(x_array)-1;
ney  = length(y_array)-1;

[x_nod_mtx,y_nod_mtx]=meshgrid(x_array,y_array);

keynodes            = zeros(size(x_nod_mtx,1),size(x_nod_mtx,2)+1);
keynodes(:,2:end-1) = (x_nod_mtx(:,1:end-1)+x_nod_mtx(:,2:end))/2;
keynodes(:,1)       = x_nod_mtx(:,1);
keynodes(:,end)     = x_nod_mtx(:,end);
dx_cell_mtx         = diff(keynodes,1,2); % check inbuilding function diff

keynodes            = zeros(size(y_nod_mtx,1)+1,size(y_nod_mtx,2));
keynodes(2:end-1,:) = (y_nod_mtx(1:end-1,:)+y_nod_mtx(2:end,:))/2;
keynodes(1,:)       = y_nod_mtx(1,:);
keynodes(end,:)     = y_nod_mtx(end,:);
dy_cell_mtx         = diff(keynodes,1,1); % check inbuilding function diff


nn       = nx*ny;
ne       = (nx-1)*(ny-1);
sequence = 'yxz';

%dataset 14
ii = (1:nn)';
x_nod_array  = x_nod_mtx(:);
y_nod_array  = y_nod_mtx(:);


%node_index_mtx = flip(reshape(ii,ny,nx));  %note node_index_mtx(52) = 52
node_index_mtx = reshape(ii,ny,nx);  %note node_index_mtx(52) = 52

% the reason of having a gravity compensated matrix is that the node reflects a xy domain with gravity facing down. 
% the original matrix has the same shape as the one with gravity compensation but it is bottom-up

node_index_mtx_gravity_compensated = flip(node_index_mtx);   % this makes a matrix that maps the node position in a xy plane.
y_nod_mtx_gravity_compensated_m      = flip(y_nod_mtx);
x_nod_mtx_gravity_compensated_m      = flip(x_nod_mtx);
dx_cell_mtx_gravity_compensated    = flip(dx_cell_mtx);

idx = 1;

iin1 = zeros(ne,1);
iin2 = zeros(ne,1);
iin3 = zeros(ne,1);
iin4 = zeros(ne,1);
x_ele_array=zeros(ne,1);
y_ele_array=zeros(ne,1);

for j  = 1:nex
    for i = 1:ney
        iin1(idx) = node_index_mtx(i,j);
        iin2(idx) = node_index_mtx(i,j+1);
        iin3(idx) = node_index_mtx(i+1,j+1);
        iin4(idx) = node_index_mtx(i+1,j);
        x_ele_array(idx)= mean([x_nod_array(iin1(idx)),x_nod_array(iin2(idx)),x_nod_array(iin3(idx)),x_nod_array(iin4(idx))] )   ;  % use mean method to get the middle of the cell
        y_ele_array(idx)= mean([y_nod_array(iin1(idx)),y_nod_array(iin2(idx)),y_nod_array(iin3(idx)),y_nod_array(iin4(idx))] )   ;  % use mean method to get the middle of the cell
        idx       = idx+1;

    end
end

x_ele_mtx_m = reshape(x_ele_array, ny-1,nx-1);
y_ele_mtx_m = reshape(y_ele_array, ny-1,nx-1);


x_ele_mtx_gravity_compensated_m = flip(x_ele_mtx_m);
y_ele_mtx_gravity_compensated_m = flip(y_ele_mtx_m);


%% dataset 17
%iqcp = -node_index_mtx_gravity_compensated(1,:)' ; % top cell,negative means evaporation is included
%qinc = dx_cell_mtx_gravity_compensated(1,:)' *dz ;
%uinc = zeros(size(iqcp)) +1e-4;   % evaporating water concentration




%ipbc = node_index_mtx_gravity_compensated(end,:)'; % only the last column. note that the result needs to be in a column
%pbc  = zeros(size(ipbc)) - 7.37325e+03;
%ubc  = zeros(size(ipbc)) + 3.0e-03;


% PART1 is the name
vapinp = vapinpObj('PART1','read_from_file','no');   % setup a empty inpObj

% ##  DATASET 1
vapinp.title1 = '2-D flume simulation for medium sand (2021)';
vapinp.title2 = 'Generating input using sutralab';

% ##  DATASET 2A
%'SUTRA VERSION 2.2 SOLUTE TRANSPORT'
%'2D REGULAR MESH' 105 4001

vapinp.vermin = '2.2';  % note this should be a string not a number;
vapinp.simula = 'SOLUTE';


% ##  DATASET 2B
 %        2d mesh          ==>   ktype(1) = 2
 %        3d mesh          ==>   ktype(1) = 3
 %        irregular mesh   ==>   ktype(2) = 0
 %        layered mesh     ==>   ktype(2) = 1
 %        regular mesh     ==>   ktype(2) = 2
 %        blockwise mesh   ==>   ktype(2) = 3

vapinp.ktype(1)  = 2;  % 2D mesh
vapinp.mshtyp{1} = '2D';
vapinp.mshtyp{2} = 'REGULAR';

vapinp.nn1 = ny;
vapinp.nn2 = nx;


% ##  DATASET 3:  Simulation Control Numbers
vapinp.nn   = nn;
vapinp.ne   = ne;
vapinp.npbc = 0;   %length(ipbc);  revised after dataset 19
vapinp.nubc = 0;
vapinp.nsop = 0;       % dataset 17
vapinp.nsou = 0;
vapinp.nobs = 0;


%%##  DATASET 4:  Simulation Mode Options

vapinp.cunsat = 'UNSATURATED';
vapinp.cssflo = 'TRANSIENT FLOW';
vapinp.csstra = 'TRANSIENT TRANSPORT';
vapinp.cread  = 'COLD' ;
vapinp.istore = 9999;


%%##  DATASET 5:  Numerical Control Parameters
vapinp.up   = 0;
vapinp.gnup = 0.01;
vapinp.gnuu = 0.01;


%  DATASET 6:  Temporal Control and Solution Cycling Data
%  
vapinp.nsch  = 1;
vapinp.npcyc = 1;
vapinp.nucyc = 1;

%DATASET 6:  Temporal Control and Solution Cycling Data

vapinp.schnam = 'TIME_STEPS';
vapinp.schtyp = 'TIME CYCLE';
vapinp.creft  = 'ELAPSED';
vapinp.scalt  = 10.; 
vapinp.ntmax  = 170000; %around 20 days
vapinp.timei  = 0;
vapinp.timel  = 1.e99;
vapinp.timec  = 1.;
vapinp.ntcyc  = 9999;
vapinp.tcmult = 1;
vapinp.tcmin  = 1.e-20;
vapinp.tcmax  = 1.e99;


%##  DATASET 7:  ITERATION AND MATRIX SOLVER CONTROLS
%##  [ITRMAX]        [RPMAX]        [RUMAX]
vapinp.itrmax = 100;
vapinp.rpmax  = 5e+3;
%vapinp.rpmax = 5e-2;  TO200317 too strigent for the first step
vapinp.rumax  = 1.0e-2;
% ##  [CSOLVP]  [ITRMXP]         [TOLP]
vapinp.csolvp = 'ORTHOMIN' ;
vapinp.itrmxp = 2000;
vapinp.tolp   = 1.e-12;

%##  [CSOLVU]  [ITRMXU]         [TOLU]
vapinp.csolvu = 'ORTHOMIN';
vapinp.itrmxu = 2000;
vapinp.tolu   = 1.e-12;


%##  DATASET 8:  Output Controls and Options
%## [NPRINT]  [CNODAL]  [CELMNT]  [CINCID]  [CPANDS]  [CVEL]  [CCORT] [CBUDG]   [CSCRN]  [CPAUSE]
%   2920        'N'        'N'        'N'        'Y'     'Y'        'Y'    'Y'      'Y' 'Y' 'Data Set 8A'

vapinp.nprint = 100;
vapinp.cnodal = 'N';
vapinp.celmnt = 'N';
vapinp.cincid = 'N';
vapinp.cpands = 'N';
vapinp.cvel   = 'N';
vapinp.ccort  = 'N';
vapinp.cbudg  = 'Y';
vapinp.cscrn  = 'N';   % screen output
vapinp.cpause = 'Y';


%## [NCOLPR]    [NCOL]
%     -1000  'N'  'X'  'Y'  'P'  'U'  'S'  '-' 
%## [LCOLPR]    [LCOL]
%     1000 'E'  'X'  'Y'  'VX' 'VY' '-' 
%## [NOBCYC]    [INOB]
%
%##  [NBCFPR]  [NBCSPR]  [NBCPPR]  [NBCUPR]  [CINACT]
%     1000         1000       1000      1000       'Y'
%##

vapinp.ncolpr = -100;
%vapinp.ncol  = 'N'  'X'  'Y'  'P'  'U'  'S'  '-';
vapinp.ncol   = {['N'],['X' ],[ 'Y'  ],['P' ],[ 'U' ],[ 'S' ],[ '-']};

vapinp.lcolpr = 100;
vapinp.lcol   = {[ 'E' ],[ 'X' ],[ 'Y'  ],['VX' ],['VY' ],['-']};


vapinp.nbcfpr = 100;
vapinp.nbcspr = 100;
vapinp.nbcppr = 100;
vapinp.nbcupr = 100;
vapinp.cinact = 'Y';



%##    DATASET 9:  FLUID PROPERTIES
%##     [COMPFL]           [CW]       [SIGMAW]        [RHOW0]       [URHOW0]        [DRWDU]        [VISC0]
%         0.0                1.         3.890D-10         1.0E+3             0.     7.0224E+02         1.0E-3
%##
vapinp.compfl = 1.e-11;
vapinp.cw     = 1.;
vapinp.sigmaw = 1.e-9;
vapinp.rhow0  = 1000;
vapinp.urhow0 = 0.;
vapinp.drwdu  = 700.;
vapinp.visc0  = 0.001;
vapinp.dvidu  = 3.25e-3;

%##      DATASET 10:  SOLID MATRIX PROPERTIES
vapinp.compma = 1.e-7;
vapinp.cs     = 0.;
vapinp.sigmas = 0.;
vapinp.rhos   = 2600.;   %solid density of sodium chloride

%##  DATASET 11:  Adsorption Parameters
%##     [ADSMOD]         [CHI1]         [CHI2]
%#'NONE'
%'FREUNDLICH' 1.D-47 0.05  #less rigid

%vapinp.adsmod = 'FREUNDLICH';
vapinp.adsmod = 'NONE';
vapinp.chi1   = 1.e-46;
vapinp.chi2   = 0.05 ;


%##
%##  DATASET 12:  Production of Energy or Solute Mass
%##     [PRODF0]       [PRODS0]       [PRODF1]       [PRODS1]
%0.             0.             0.             0.
vapinp.prodf0 = 0.;
vapinp.prods0 = 0.;
vapinp.prodf1 = 0.;
vapinp.prods1 = 0.;


%##
%##  DATASET 13:  Orientation of Coordinates to Gravity
%##      [GRAVX]        [GRAVY]        [GRAVZ]
%0.           -9.81          0.
vapinp.gravx  = 0;
vapinp.gravy  = -9.81;
vapinp.gravz  = 0;

%##  DATASET 13B: TIDE FLUCTUATION IN USUBS
%#   TA    -- TIDAL AMPLITUDE(M); 
% #   TP    -- TIDAL PERIOD(S); 
% #   TM    -- MEAN TIDAL LEVEL[M]; 
% #   RHOST -- THE DENSITY FOR TIDE WATER; 
% #   SC    -- SALINITY OF THE SEAWATER; 
% #   ITT   -- ITERRATION CRITERIA FOR BCTIME SHOULD BE LARGE OR EQUAL THAN 2
% ##         [TA]           [TP]           [TM]        [RHOST]        [SC]    [ITT]
            % 0.0       4.32D+4          0.D+0        1025.0D+0    3.57D-02      2

vapinp.ta     = 0.;
vapinp.tp     = 4.32e4;
vapinp.tm     = 0.;
vapinp.rhost  = 1025.;			
vapinp.sc     = 3.57e-2;
vapinp.itt    = 2;				
			

% ##  DATASET 13C: CONTROLLING PARAMETERS
% ##    MET...SWICH OF EVAPORATION. =1 PENMAN EQUATION =2 PDV EQUATION =0 OFF
% ##    MAR...AERODYNAMIC RESISTANCE SWICH =1 CONSTANT =0 OFF
% ##    MSR...SURFACE RESISTANCE SWICH 
% #            MSR=0  SYRFSIS=0;
% #            MSR=1 CAMILO(1986) FORMULAR ;
% #            MSR=2 SUN (1982) FORMULAR;
% #            MSR=3 DAAMEN (1996) FORMULAR;
% #            MSR=4 Van der Griend FORMULAR (1994);
% #            MSR=5 USING HYPERGEOMETRIC FORMULAR
% #            MSR=6 USING APPROXIMATED FORMULAR PROPOSED IN PAPER 2
% #            MSR=9 USING FUNNEL-SHAPED SURFACE RESISTANCE MODEL WITH ONLY TSL 
% #            MSR=10 USING FUNNEL-SHAPED SURFACE RESISTANCE MODEL WITH TSL & NSL (EC=0.5)
% #            MSR=11 USING FUNNEL-SHAPED SURFACE RESISTANCE MODEL WITH TSL & NSL UNDER BROOKS&COREY MODEL
% ##    MSC...SALT RESISTANCE SWICH =0 OFF =1 FUJIMAKI(2006)
% ##    MHT...HEAT BALANCE =0 OFF =1 USING MASS BALANCE EQUATION;=2 ENERGY TRANSPORT EQUATION; =3 USING INITIAL EQUATION
% ##    MVT.. VAPOR TRANSPORT SWICH  ON=1 OFF=0
% ##    MFT...(F)ILM (T)RANSPORT SWICH SEE PAGE 165 SEE PAGE 166 FOR REFERENCE
% #             MFT=0, FILM TRANSPORT DISABLED
% #             MFT=1, FILM TRANSPORT IS ENABLED BASED ON TOKUNAKA (2009) AND ZHANG(2010) 
% #             MFT=2, FILM TRANSPORT IS CALCULATED BY MUALEM (1976) BASED ON TOTAL SATURATION. SEE PAGE 170
% #             MFT=3, FILM TRANSPORT IS ENABLED BASED ON LEBEAU AND KONRAD (2010) SEE PAGE 178 AND 179
% ##    MRK..[NOT USING AT THE MOMENT].(R)ELATIVE HYDRAULIC CONDUCTIVITY CALCULATION METHOD. ONLY APPLY WHEN BROOKS AND COREY METHOD IS USED
% #              MRK=0, USE TRADITIONAL RELATIVE K
% #              MRK=1, USE MUALEM (1976) BASED ON EFFECTIVE SATURATION. 
% #              MRK=2, USE FAYER (1995) BASED RELATIVE K
% ##   [MET] [MAR] [MSR] [MSC]  [MHT]  [MVT]  [MFT]     [MRK]
      % 2      1    9    1      2     1        1         1
vapinp.met    = 2;
vapinp.mar    = 1;
vapinp.msr    = 9;
vapinp.msc    = 1;			
vapinp.mht    = 2;
vapinp.mvt    = 1;				
vapinp.mft    = 1;	
vapinp.mrk    = 1;	
			
      
% ##  DATASET 13D: EVAPORATION SINARIO
% #   WHEN QET=1 THEN USING PENMAN (1948) EQUATION, THIS IS FURTHER MODIFIED BY KONUKCU (2007)
% #   WHEN QET=2 THEN USING PDV (1957) EQUATION.
% ##  THE UNIT OF QET IS (MM/DAY)
% ##  PET SHOULD LESS THAN ZERO ALWAYS. IF NOT, THE MATRIC POENTIAL PARAMETER IN PENMAN EQUATION MIGHT BE LARGE THAN 0! 
% #   UVM -- SOLUBILITY
% #     QET... (MM/DAY) NON OF USE SO FAR
% #     UET... (KG/KG) THE SOLUTE DENSITY OF WATER TAKEN BY EVAPORATION
% #     PET... (PA) WHEN PORE WATER PRESSURE REACHED THE POINT LOWER THAN PET, EVAPORATION IS TAKING PLACE
% #     UVM... (KG/KG) SOLUBILITY
% #     NGT... IF =1 EVAPORATION IS SWICHED OFF DURING NIGHT TIME, IF =0 IT IS ALWAYS ON
% #     ITE... TEMPERATORYLY NOT USING, IT WAS DESIGNED FOR THE NUMBER OF TIME CALL FOR BCTIME
% ##         [QET]         [UET]          [PET]          [UVM]        [NIGHT]     [ITE]
            % 2.0          1.D-5       -0.001D+00      0.264D0           0          1
vapinp.qet    = 2.;
vapinp.uet    = 1.e-5;
vapinp.pet    = -0.001;
vapinp.uvm    = 0.264;			
vapinp.night  = 0;
vapinp.ite    = 1;				
					
            
% ##  DATASET 13E: EVAPORATION PARAMETER
% ##    TMA  -- (A)TMOSPHERIC (T)EMPERATURE ON THE TOP OF THE COLUMN [CELSIUS]                   IT IS USEFULL WHEN MET=2
% ##    TMI  -- (I)NITIAL SOIL TEMPERATURE [CELSIUS]                    IT IS USEFULL WHEN MET=2
% ##    ALF  -- ALBEDO [I]
% ##    RS   -- SHORT WAVE INCOMING RADIATION  [MJ/M2/DAY]
% ##    RH   -- RELATIVE HUMIDITY 0<RH<1                             IT IS USEFULL WHEN MET=2
% ##    AP,BP-- PARAMETER IN WIND FUNCTION                     CURRENTLY NOT USED
% ##    U2   -- DAILY MEAN WIND SPEED AT 2M ABOVE GROUND [KM/DAY]
% ##    TSD  -- (T)EMPERATURE ON THE (SI)DE OF THE COLUMN, THIS IS USEFUL WHEN INSULATION APPLIES ON THE SIDE OF 1-D HEAT TRANSPORT EQ  [CELSIUS]
% ##    SCF  -- (SC)ALING (F)ACTOR BETWEEN THE TOP SURFACE AND SIDE SURFACE OF THE 1-D COLUMN. THIS IS IMPORTANT WHEN USING CUBIC CELLS TO 
% ##              REPRESENT CYLINDRICAL ONES. SEE 2012 PAGE 38,59,68 FOR MORE INFORMATION.
% ##              THE EQUATION OF THE SCF=DELTAX/RADIUS, WHERE DELTAX IS THE DISTANCE BETWEEN THE FIRST AND THE SECOND NODE IN X DIRECTION,
% ##              RADIUS IS RADIUS OF THE PROTOTYPE CYLINDRICAL (LAB)
% ##              SEE PAGE
% ##         [TMA]       [TMI]      [ALF]        [RS]        [RH]        [AP]        [BP]       [U2]       [TSD]    [SCF]
            % 25.D0      25.0       0.2       58.67214    0.52     17.8636        0.044     329.5       24.2D0      0.1818
vapinp.tma   = 35.;
vapinp.tmi   = 25.;
vapinp.alf   = 0.2;
vapinp.rs    = 58.67214;			
vapinp.rh    = 0.52;
vapinp.ap    = 17.8636;					
vapinp.bp    = 0.044;					
vapinp.u2    = 329.5;					
vapinp.tsd   = 24.2;					
vapinp.scf   = 0.1818;					

% ##  DATASET 13F: AERODYNAMIC RESISTANCE TERM
% ##     RAVT... (S/M) AERODYNAMIC RESISTANCE AT THE SOIL SURFACE
% ##     RAVS.. (S/M) AERODYNAMIC RESISTANCE AT THE SIDE OF THE COLUMN, WARNING: BE AWARE OF THE SIZE EFFECT
% ##     SWRAT..  (-)   PARAMETER TO SWICH ON (1.D) OR OFF(0.D0) THE TEMPERATURE CHANGE ON THE SURFACE
% ##         [RAVT]    [RAVS]   [SWRAT]
           % 206.D0     50.D0     0.D0
vapinp.ravt   = 35.;
vapinp.ravs   = 50.;
vapinp.swart  = 0.;
          
% ##  DATASET 13G: SOIL THERMO PROPERTY TERM           
% ##     HSC  -- HEAT SOURCE FROM ATMOSPHERE G [WATT/M2]
% ##     NHER -- NUMER OF NODES THAT CONSIDERS HEAT BALANCE    
% ##     ROUS -- THE DENSITY OF SOLID GRAIN     [KG/M3]
% ##     HCS  -- HEAT CAPACITY OF SOLID GRAIN   [J/KG/K]
% ##      
% ##     [HSC]      [HER]      [ROUS]     [HCS]
      % 9.00D2        1.95D-1      2650.D0   8.75D2
vapinp.hsc   = 900.;
vapinp.her   = 1.95e-1;
vapinp.rous  = 2650.;
vapinp.hcs   = 8.75e+2

% ##  DATASET 13H, PARAMETERS FOR THE SALT RESISTANCE
% ##    AR  FITTING PARAMETER
% ##    BR  FITTING PARAMETER   
% ##          [AR]   [BR]
         % 6.9D-1  -1.04D-0
vapinp.hsc   = 6.9e-1;
vapinp.her   = -1.04;
        
% ##  DATASET 13I SOIL CHARACTERISTIC PARAMETERS
% ##  USING EACH PARAMETERS ARE CONTROLLED BY NREG AND LREG AT DATASET 14 AND DATASET 15
% ##  THE UNIT OF AA1, PHY0 AND PHYB IS [M-1]
% ##  ECTO--- ESSENTRICITY AND TORTUOSITY THE DEFAULT VALUE IS 0.5
% ##  [SWRES1] [AA1]  [VN1]  [SWRES2]    [AA2]   [VN2]  [SWRES3]  [LAM3]   [PHYB3]  [SWRES4]   [LAM4]  [PHYB4]  [PHY0]   [ECTO]
      % 0.06   14.5D0  8.5D0    0.06      15.D0   9.2D0    0.09D0   8.0D0   0.2D0   0.08D0     4.2D0   0.045    5.0D4     0.5D0
vapinp.swres1  = 0.09 
vapinp.aa1     = 4.5   
vapinp.vn1     = 12.
vapinp.swres2  = 0.
vapinp.aa2     = 0.
vapinp.vn2     = 0.
vapinp.swres3  = 0.
vapinp.lam3    = 0.
vapinp.phyb3   = 0.
vapinp.swres4  = 0.
vapinp.lam4    = 0.
vapinp.phyb4   = 0.
vapinp.phy0    = 5.e+4
vapinp.ecto	   = 0.5
	  
	  
% ##  
% ##  DATASET 13J: THERMAL CONDUCTIVITIES OF WATER AND LIQUID
% ##   NTC -- TYPE OF (T)HERMAL (C)ONDUCTIVITIES EQUATION: 
% ##              NTC=1, USING EQUATIONS FROM SUTRA MANUAL, THEN [TCS] AND [TCL] HAS TO BE INPUT 
% ##                     TCS -- THERMAL CONDUCTIVITY OF SOIL  USUALLY BETWEEN 2-4 W/M/K
% ##                     TCL -- THERMAL CONDUCTIVITY OF LIQUID WATER, 0.6      W/M/K
% ##              NTC=2, USING EQUATION FROM CHUNG AND HORTON (1987) THREE PARAMETERS [B1] [B2] [B3] IS REQUIRED 
% ##  [NTC]      [TCS]    [TCL]
% ##    1         3.D0      .6D0
% ##
% ##  [NTC]      [B1]      [B2]        [B3]
      % 2        1.528       -2.406      4.909
vapinp.ntc    = 2.
vapinp.b1     = 1.528
vapinp.b2     = -2.406
vapinp.b3	  = 4.909	  
	  
% ##
% ##  DATASET 13K ENHANCEMENT FACTOR
% ##  SEE CASS (1984) AND SAITO (2006) FOR DETAIL
% ##  SUGGESTED PARAMETERS ARE YA-9.5, YB-3.0, YC-2.6, YD-1.0, YE-4.0,FC--0.001
% ##     FC   -- THE PROPORTION OF CLAY IN THE SOIL  0<FC<1 USED IN
% ##   NEF -- THE SWITCH TO THE ENHANCEMENT FACTOR
% ##            =0, ENHANCEMENT FACTOR EQUALS TO 1
% ##            =1, APPLY ENHANCEMENT FACTOR TO THE TEMPERATURE GRADIENT TERM ONLY
% ##            =2, APPLY ENHANCEMENT FACTOR TO ALL OF THE TERMS
% ##  [NEF]     [YA]   [YB]    [YC]   [YD]   [YE]   [FC]
      % 0       9.5D0  3.0D0   2.6D0  1.0D0 4.0D0  1.D-3
vapinp.nef    = 0.
vapinp.ya     = 9.5
vapinp.yb     = 3.0
vapinp.yc	  = 2.6	  
vapinp.yd     = 1.0
vapinp.ye     = 4.0
vapinp.fc     = 1.e-3
	  
% ##
% ## DATASET 13L PARAMETERS FOR SURFACE RESISTANCE
% ## SEE PAGE 142BLUEBOOK, 98(5) FOR REFERENCE
% ## [TAL]  --  THICKNESS OF AIR LAYER (M) [TAL]
% ## [EC]   --  ECCENTRICITY OF THE ACTIVE PORE, WHEN MSR=6 (DENOTING USING APPROXIMATED SOLUTION GIVEN BY PAPER 2), THIS IS THE EXPOENTIAL
% ##             OF G FUNCTION
% ## [ETR]  --  RESIDUAL EVAPORATION RATE [-] [ETR]
% ##  [TAL]   [EC]   [ETR]     [PSIP]      [CORS]
     % 0.3D-3   5.D-1   1.D0    1000.D0      1.D0
vapinp.tal    = 0.3e-3
vapinp.ec     = 0.5
vapinp.etr    = 1.
vapinp.psip	  = 1000.	  
vapinp.cors   = 1.0	 
	 
% ## DATASET 13M FILM TRANSPORT ALGORITHM 
% ##   THIS ARGORITHM IS BASED ON ZHANG (2010) AND TOKUNAGA (2009)
% ##   ENABLED WHEN MRK=1
% ##   [AGR]  -- AVERAGE GRAIN RADIUS [M]
% ##   [CORF] -- CORRECTION FACTOR FOR THE SATURATED KFILM. SEE ZHANG (2010)
% ##   [PHICM]-- CAPILLARY PRESSURE HEAD THAT SURFACE ROUGHNESS OF THE SOILID GRAIN BECOMES NEGNIGIBLE. (M) IT IS SUGGESTD
% #                AS 1000M FROM TULLER AND OR (2005)
% ##   [ASVL] -- HAMAKER CONSTANT. IT IS SUGGESTED BY TULLER AND OR THAT THIS CONSTANT EQUALS TO -6E-20
% ##    [CORF]             [AGR]     [PHICM]      [ASVL]
        % 1.D0             0.00035      1.E3       -6.0E-20
vapinp.corf    = 1.0
vapinp.agr     = 3.5e-4
vapinp.phicm   = 1000.
vapinp.asvl	   = -6.0e-20	


%##  DATASET 14:  NODEWISE DATA
%%##                              [SCALX] [SCALY] [SCALZ] [PORFAC]
vapinp.scalx  = 1.;
vapinp.scaly  = 1.;
vapinp.scalz  = 1.;
vapinp.porfac = 1.;

%##      [II]    [NRE    G(II)]  [X(II)] [Y(II)] [Z(II)] [POR(II)]
%vapinp.ii   = (1:nn)';
vapinp.nreg = zeros(nn,1)+1;
vapinp.x    = x_nod_array;
vapinp.y    = y_nod_array;
%vapinp.z    = zeros(nn,1)+dz;
vapinp.z    = 2*pi*vapinp.x+0.1;
vapinp.por  = zeros(nn,1)+0.43;


%##                              [PMAXFA]        [PMINFA]        [ANG1FA]        [ALMAXF]        [ALMINF]        [ATMAXF]        [ATMINF]
%'ELEMENT'               1.0000000D+00   1.0000000D+00   1.0000000D+00   2 2 2 2
%##     [L]      [LREG(L)]       [PMAX(L)]       [PMIN(L)]       [ANGLEX(L)]     [ALMAX(L)]      [ALMIN(L)]      [ATMAX(L)]      [ATMIN(L)]
vapinp.pmaxfa = 1.;
vapinp.pminfa = 1.;
vapinp.angfac = 1.;
vapinp.almaxf = 1.;
vapinp.alminf = 1.;
vapinp.atmaxf = 1.;
vapinp.atminf = 1.;
%vapinp.l      = (1:ne)';
vapinp.lreg   = zeros(ne,1)+1;


pmax_mtx_gravity_compensated_m2= zeros(size(x_ele_mtx_gravity_compensated_m)) + permeability_sand_m2 ;   %background permeability

mask_ele_mtx_silt_layer_gravity_compensated = y_ele_mtx_gravity_compensated_m  > -6  ;   % mask matrix, for element matrix 

%pmax_mtx_gravity_compensated_m2 (mask_ele_mtx_silt_layer_gravity_compensated) =  permeability_silt_m2;   % silt layer permeability

pmax_mtx_m2=flip(pmax_mtx_gravity_compensated_m2);

pmax_array_m2 = pmax_mtx_m2(:);



vapinp.pmax   = pmax_array_m2;
vapinp.pmin   = pmax_array_m2;
vapinp.anglex = zeros(ne,1);
vapinp.almax  = zeros(ne,1)+0.5e-0;
vapinp.almin  = zeros(ne,1)+0.5e-0;
vapinp.atmax  = zeros(ne,1)+0.5e-0;
vapinp.atmin  = zeros(ne,1)+0.5e-0;

% note: for SUTRASET when node number is negative, the second input is surface area of the node
% and the third input is the thickness of the cell.

%% DATASET 17   Data for Fluid Source and Sinks

%vapinp.iqcp  = iqcp;
%vapinp.qinc  = qinc;
%vapinp.uinc  = uinc;


% ## DATASET 19:  Data for Specified Pressure Nodes
%###  [IPBC]                [PBC]                [UBC]

mask_nod_mtx_aquifer_boundary_gravity_compensated = and(y_nod_mtx_gravity_compensated_m<-4, x_nod_mtx_gravity_compensated_m>99.99);  % below 4 metre, greater than 200 m away from the centre
ipbc_node_idx_array                           = node_index_mtx_gravity_compensated(mask_nod_mtx_aquifer_boundary_gravity_compensated);
pbc                                           = -(y_nod_mtx_gravity_compensated_m(mask_nod_mtx_aquifer_boundary_gravity_compensated) - initial_head_aquifer_m ) *c.g * (vapinp.rhow0 + vapinp.drwdu * c_saltwater_kgPkg);



%mask_mtx_aquifer_boundary_gravity_compensated_left = and(y_nod_mtx_gravity_compensated_m<-4, x_nod_mtx_gravity_compensated_m<0.01);
%ipbc_node_idx_array_left                           = node_index_mtx_gravity_compensated(mask_mtx_aquifer_boundary_gravity_compensated_left);
%pbc_left                                      = -(y_nod_mtx_gravity_compensated_m(mask_mtx_aquifer_boundary_gravity_compensated) - initial_head_aquifer_m +0.5 ) *c.g * vapinp.rhow0 ;
%
%vapinp.ipbc = [ipbc_node_idx_array; ipbc_node_idx_array_left];
%vapinp.pbc  = [pbc;pbc_left];
%vapinp.ubc  = [zeros(size(pbc))+c_saltwater_kgPkg;zeros(size(pbc_left))+c_freshwater_kgPkg];
%vapinp.npbc = length(vapinp.pbc);



mask_nod_mtx_aquifer_boundary_gravity_compensated_top = and(y_nod_mtx_gravity_compensated_m>-0.1, x_nod_mtx_gravity_compensated_m<pond_radius_m);
ipbc_node_idx_array_top                           = node_index_mtx_gravity_compensated(mask_nod_mtx_aquifer_boundary_gravity_compensated_top);
pbc_top                                      = zeros(size(ipbc_node_idx_array_top));

%vapinp.ipbc = [ipbc_node_idx_array; ipbc_node_idx_array_top];
%vapinp.pbc  = [pbc;pbc_top];
%vapinp.ubc  = [zeros(size(pbc))+c_saltwater_kgPkg;zeros(size(pbc_top))+c_freshwater_kgPkg];
%vapinp.npbc = length(vapinp.pbc);



vapinp.ipbc = ipbc_node_idx_array;
vapinp.pbc  = pbc;
vapinp.ubc  = zeros(size(pbc))+c_saltwater_kgPkg;
vapinp.npbc = length(vapinp.pbc);


%##
%##  DATASET     22:  Ele        ment Incid      ence Data
%##    [LL]      [IIN(1)]        [IIN(2)]        [IIN(3)]        [IIN(4)]

%ne_mtx=reshape(l,ney,nex);
vapinp.iin1 = zeros(vapinp.ne,1);
vapinp.iin2 = zeros(vapinp.ne,1);
vapinp.iin3 = zeros(vapinp.ne,1);
vapinp.iin4 = zeros(vapinp.ne,1);

% DATASET 22, user need to check if the sequence of the node is in clockwise manner as suggested by the manual

%for j  = 1:nex
%    for i = 1:ney
%        vapinp.iin1(idx) = node_index_mtx(i,j);
%        vapinp.iin2(idx) = node_index_mtx(i,j+1);
%        vapinp.iin3(idx) = node_index_mtx(i+1,j+1);
%        vapinp.iin4(idx) = node_index_mtx(i+1,j);
%        idx           = idx+1;
%    end
%end


vapinp.iin1= iin1;
vapinp.iin2= iin2;
vapinp.iin3= iin3;
vapinp.iin4= iin4;

%vapinp.export_to_file();


% setting the initial pressure as hydrostatic, in particular at the sandy aquifer, the silt layer will be overwritten
pm1_mtx_gravity_compensated_pa= - (- initial_head_aquifer_m + y_nod_mtx_gravity_compensated_m)*c.g*c.rhow_pure_water;


mask_nod_mtx_silt_layer_gravity_compensated = y_nod_mtx_gravity_compensated_m  > -6  ;   % mask matrix, for nod matrix 

%set a relatively high suction in the silt layer
pm1_mtx_gravity_compensated_pa ( mask_nod_mtx_silt_layer_gravity_compensated) = -20000;




% put the top cell as zero pressure 
%pm1_mtx_gravity_compensated_pa(mask_nod_mtx_aquifer_boundary_gravity_compensated_top)=initial_pond_water_depth_m*c.rhow_pure_water*c.g;

pm1_mtx_pa=flip(pm1_mtx_gravity_compensated_pa);


% initial solute concentrtion, sandy aquifer has a saline concentration while silt has a fresh concentration
um1_mtx_gravity_compensated_kgPkg= zeros(size(pm1_mtx_gravity_compensated_pa))+ c_saltwater_kgPkg;



um1_mtx_gravity_compensated_kgPkg (mask_nod_mtx_silt_layer_gravity_compensated) = c_freshwater_kgPkg;


um1_mtx_kgPkg=flip(um1_mtx_gravity_compensated_kgPkg);





fprintf('use the original ics file')
% ics file
ics       = icsObj('PART1','read_from_file','no');
ics.tics  = 0;
%ics.cpuni = 'UNIFORM';
%ics.pm1   = -30;
ics.cpuni = 'NONUNIFORM';
ics.pm1   = pm1_mtx_pa(:);
ics.cuuni = 'NONUNIFORM';
ics.um1   = um1_mtx_kgPkg ;
ics.export_to_file();
%
%
%%pm1_mtx=reshape(b.pm1,[vapinp.nn1,vapinp.nn2]);
%%
%%b=icsObj('PART1_cs.csv','inpObj',vapinp);
%%pm1_mtx=reshape(b.pm1,[vapinp.nn1,vapinp.nn2]);
%%um1_mtx=reshape(b.um1,[vapinp.nn1,vapinp.nn2]);
%%vapinp.get_x_nod_mtx;
%%vapinp.get_y_nod_mtx;
%%a=figure;
%%surf(vapinp.x_nod_mtx,vapinp.y_nod_mtx,um1_mtx);
%%savefig(a,'conc.fig');
%%
%%a=figure;
%%surf(vapinp.x_nod_mtx,vapinp.y_nod_mtx,pm1_mtx);
%%savefig(a,'p.fig');
%%
%%saveas(a,'Barchart.png')
%
%
%ics.export_to_file();
