clear all
close all
run('/storage/macondo/s4524462/SutraLab/mfiles/slsetpath.m')																																	   
c=ConstantObj();  % all the constants, we suggest to use c.g rather than 9.81 in the script to enhance readerbility;

%% SUTRA.fil

fil=filObj('read_from_file','no');

% the default value for 'fid' is -1, meaning, this will not write to SUTRA.fil if export
fil.terms.('inp').('fname')  = 'FLUME.inp' ; fil.terms.('inp').('fid')   = 50;
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

%%
%soil properties 
porosity                    = 0.39;
c_saltwater_kgPkg           = 0.035;
% c_freshwater_kgPkg          = 0.0001;
constant_water_table_m      = 0.24;
relative_humidity			= 0.57;
aero_resistance				= 86;
% choose boundary condition from 'right_side', 'left_side', 'bottom','both_side'.
pressure_boundary			='right_side';
boundary_height				=0.05;

initial_temperature_C       = 23.9;
initial_concentration_kgPkg = 0.035;
initial_head_aquifer_m      = 0.24;
permeability_sand_m2        = 1.e-9;

time_scale=5; 
time_cycle=11250; %around 7 days
diffusivity=1e-9;

%% input for meshing 
BLOCK = 1; %input block number (every block must be a quadrilateral and the left and right side
%								should be parallel to the y axis)
x1  = 0;
x2  = 1.2;
nex = 400; %Number of segments along x

y(1,1)=0;
y(1,2)=0;

ney_section(1)=50; %Number of segments between y1 and y2, must be consistente in every block

y(2,1)=0.305;
y(2,2)=0.2;

ney_section(2)=50;

y(3,1)=0.405;
y(3,2)=0.3;


dz=0.01; %default

mesh %mesh must be called in every block
%%
% %BLOCK 2
% BLOCK = 2;
% 
% x1 = 250;
% x2 = 260;
% nex = 10; %Number of segments along x
% 
% 
% y(1,1)=-5;
% y(1,2)=-5;
% 
% ney_section(1)=20;
% 
% y(2,1)=-1;
% y(2,2)=-2;
% 
% ney_section(2)=30;
% 
% y(3,1)=0.4;
% y(3,2)=-1.5;
% 
% ney_section(3)=50;
% 
% y(4,1)=0.9;
% y(4,2)=-0.8;
% 
% dz=0.01;
% mesh
%%
nn       =  nx*ny;
ne       =  nex*ney;
sequence = 'yxz';

%dataset 14
ii = (1:nn)';
x_nod_array  = x_nod_mtx(:);
y_nod_array  = y_nod_mtx(:);


%node_index_mtx = flip(reshape(ii,ny,nx));  %note node_index_mtx(52) = 52
node_index_mtx = reshape(ii,ny,nx);  %note node_index_mtx(52) = 52

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

% the reason of having a gravity compensated matrix is that the node reflects a xy domain with gravity facing down. 
% the original matrix has the same shape as the one with gravity compensation but it is bottom-up

node_index_mtx_gravity_compensated = flip(node_index_mtx);   % this makes a matrix that maps the node position in a xy plane.
y_nod_mtx_gravity_compensated_m    = flip(y_nod_mtx);
x_nod_mtx_gravity_compensated_m    = flip(x_nod_mtx);
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

%% dataset 17 ##Caution: in Sutravap, dateset 17 is used to locate the nodes starting from top layer to bottom
%               iqcp are all node numbers with a negetive sign since the
%               vapor flow not only occurs on the surface. qinc and uinc
%               are row and column where the node is located, respectively.
iqcp = -reshape(node_index_mtx_gravity_compensated',nn,1) ; % top cell,negative means evaporation is included
qinc=ny+1-mod(-iqcp,ny);
qinc(qinc>ny)=1;			   
uinc = (-iqcp + qinc-1)./ny;  



%ipbc = node_index_mtx_gravity_compensated(end,:)'; % only the last column. note that the result needs to be in a column
%pbc  = zeros(size(ipbc)) - 7.37325e+03;
%ubc  = zeros(size(ipbc)) + 3.0e-03;

%% INPUT FROM DATESET1 TO 22
% PART1 is the name
inp = vapinpObj('FLUME','read_from_file','no');   % setup a empty inpObj

% ##  DATASET 1
inp.title1 = '2-D flume simulation for medium sand (2021)';
inp.title2 = 'Generating input using sutralab';

% ##  DATASET 2A
%'SUTRA VERSION 2.2 SOLUTE TRANSPORT'
%'2D REGULAR MESH' 105 4001

inp.vermin = '2.2';  % note this should be a string not a number;
inp.simula = 'SOLUTE';


% ##  DATASET 2B
 %        2d mesh          ==>   ktype(1) = 2
 %        3d mesh          ==>   ktype(1) = 3
 %        irregular mesh   ==>   ktype(2) = 0
 %        layered mesh     ==>   ktype(2) = 1
 %        regular mesh     ==>   ktype(2) = 2
 %        blockwise mesh   ==>   ktype(2) = 3

inp.ktype(1)  = 2;  % 2D mesh
inp.mshtyp{1} = '2D';
inp.mshtyp{2} = 'REGULAR';

inp.nn1 = ny;
inp.nn2 = nx;


% ##  DATASET 3:  Simulation Control Numbers
inp.nn   = nn;
inp.ne   = ne;
inp.npbc = 0;   %length(ipbc);  updated after dataset 19
inp.nubc = 0;
inp.nsop = nn;  %dataset 17
inp.nsou = 0;
inp.nobs = 0;

%##  DATASET 4:  Simulation Mode Options
inp.cunsat = 'UNSATURATED';
inp.cssflo = 'TRANSIENT';
inp.csstra = 'TRANSIENT';
inp.cread  = 'COLD' ;
inp.istore = 9999;


%##  DATASET 5:  Numerical Control Parameters
inp.up   = 0;
inp.gnup = 0.01;
inp.gnuu = 0.01;


%DATASET 6:  Temporal Control and Solution Cycling Data  
inp.nsch  = 1;
inp.npcyc = 1;
inp.nucyc = 1;

%DATASET 6:  Temporal Control and Solution Cycling Data
inp.schnam = 'TIME_STEPS';
inp.schtyp = 'TIME CYCLE';
inp.creft  = 'ELAPSED';
inp.scalt  = time_scale; 
inp.ntmax  = time_cycle; %around 20 days
inp.timei  = 0;
inp.timel  = 1.e99;
inp.timec  = 1.;
inp.ntcyc  = 9999;
inp.tcmult = 1;
inp.tcmin  = 1.e-20;
inp.tcmax  = 1.e99;

%##  DATASET 7:  ITERATION AND MATRIX SOLVER CONTROLS
%##  [ITRMAX]        [RPMAX]        [RUMAX]
inp.itrmax = 100;
inp.rpmax  = 5e+3;
%inp.rpmax = 5e-2;  TO200317 too strigent for the first step
inp.rumax  = 1.0e-2;
% ##  [CSOLVP]  [ITRMXP]         [TOLP]
inp.csolvp = 'ORTHOMIN' ;
inp.itrmxp = 20000;
inp.tolp   = 1.e-12;

%##  [CSOLVU]  [ITRMXU]         [TOLU]
inp.csolvu = 'ORTHOMIN';
inp.itrmxu = 20000;
inp.tolu   = 1.e-12;

%##  DATASET 8:  Output Controls and Options
%## [NPRINT]  [CNODAL]  [CELMNT]  [CINCID]  [CPANDS]  [CVEL]  [CCORT] [CBUDG]   [CSCRN]  [CPAUSE]
%   2920        'N'        'N'        'N'        'Y'     'Y'        'Y'    'Y'      'Y' 'Y' 'Data Set 8A'

inp.nprint = 1000;
inp.cnodal = 'N';
inp.celmnt = 'N';
inp.cincid = 'N';
inp.cpands = 'N';
inp.cvel   = 'N';
inp.ccort  = 'N';
inp.cbudg  = 'Y';
inp.cscrn  = 'N';   % screen output
inp.cpause = 'Y';

%## [NCOLPR]    [NCOL]
%     -1000  'N'  'X'  'Y'  'P'  'U'  'S'  '-' 
%## [LCOLPR]    [LCOL]
%     1000 'E'  'X'  'Y'  'VX' 'VY' '-' 
%## [NOBCYC]    [INOB]
%
%##  [NBCFPR]  [NBCSPR]  [NBCPPR]  [NBCUPR]  [CINACT]
%     1000         1000       1000      1000       'Y'
%##

inp.ncolpr = 1000;
%inp.ncol  = 'N'  'X'  'Y'  'P'  'U'  'S'  '-';
inp.ncol   = {['N'],['X' ],[ 'Y'  ],['P' ],[ 'U' ],[ 'S' ],[ '-']};

inp.lcolpr = 1000;
inp.lcol   = {[ 'E' ],[ 'X' ],[ 'Y'  ],['VX' ],['VY' ],['-']};

inp.nbcfpr = 1000;
inp.nbcspr = 1000;
inp.nbcppr = 1000;
inp.nbcupr = 1000;
inp.cinact = 'Y';

%##    DATASET 9:  FLUID PROPERTIES
%##     [COMPFL]           [CW]       [SIGMAW]        [RHOW0]       [URHOW0]        [DRWDU]        [VISC0]
%         0.0                1.         3.890D-10         1.0E+3             0.     7.0224E+02         1.0E-3
inp.compfl = 1.e-11;
inp.cw     = 1.;
inp.sigmaw = diffusivity;
inp.rhow0  = 1000;
inp.urhow0 = 0.;
inp.drwdu  = 700.;
inp.visc0  = 0.001;
inp.dvidu  = 3.25e-3;

%##      DATASET 10:  SOLID MATRIX PROPERTIES
inp.compma = 1.e-7;
inp.cs     = 0.;
inp.sigmas = 0.;
inp.rhos   = 2600.;   %solid density of sodium chloride

%##  DATASET 11:  Adsorption Parameters
%##     [ADSMOD]         [CHI1]         [CHI2]
%#'NONE'
%'FREUNDLICH' 1.D-47 0.05  #less rigid

inp.adsmod = 'FREUNDLICH';
%inp.adsmod = 'NONE';
inp.chi1   = 1.e-49;
inp.chi2   = 0.05 ;

%##  DATASET 12:  Production of Energy or Solute Mass
%##     [PRODF0]       [PRODS0]       [PRODF1]       [PRODS1]
%0.             0.             0.             0.
inp.prodf0 = 0.;
inp.prods0 = 0.;
inp.prodf1 = 0.;
inp.prods1 = 0.;

%##  DATASET 13:  Orientation of Coordinates to Gravity
%##      [GRAVX]        [GRAVY]        [GRAVZ]
%0.           -9.81          0.
inp.gravx  = 0;
inp.gravy  = -9.81;
inp.gravz  = 0;

%% EVAPORATION RELATED INPUT
%##  DATASET 13B: TIDE FLUCTUATION IN USUBS
%#   TA    -- TIDAL AMPLITUDE(M); 
% #   TP    -- TIDAL PERIOD(S); 
% #   TM    -- MEAN TIDAL LEVEL[M]; 
% #   RHOST -- THE DENSITY FOR TIDE WATER; 
% #   SC    -- SALINITY OF THE SEAWATER; 
% #   ITT   -- ITERRATION CRITERIA FOR BCTIME SHOULD BE LARGE OR EQUAL THAN 2
% ##         [TA]           [TP]           [TM]        [RHOST]        [SC]    [ITT]
            % 0.0       4.32D+4          0.D+0        1025.0D+0    3.57D-02      2

inp.ta     = 0.;
inp.tp     = 4.32e4;
inp.tm     = 0.;
inp.rhost  = 1025.;			
inp.sc     = 3.57e-2;
inp.itt    = 2;				
			

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
inp.met    = 2;
inp.mar    = 1;
inp.msr    = 10;
inp.msc    = 1;			
inp.mht    = 0;
inp.mvt    = 1;				
inp.mft    = 1;	
inp.mrk    = 1;	
			
      
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
% #     ITE... ITE=0 SOLID SALT IS OBTAINED FROM SALT CURVE, ITE=1 SM IS INHERITED FROM *.ICS FILE
% ##         [QET]         [UET]          [PET]          [UVM]        [NIGHT]     [ITE]
            % 2.0          1.D-5       -0.001D+00      0.264D0           0          0
inp.qet    = 2.;
inp.uet    = 1.e-5;
inp.pet    = -0.001;
inp.uvm    = 0.264;			
inp.night  = 0;
inp.ite    = 0;				
					
            
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
inp.tma   = initial_temperature_C;
inp.tmi   = initial_temperature_C;
inp.alf   = 0.2;
inp.rs    = 58.67214;			
inp.rh    = relative_humidity;
inp.ap    = 17.8636;					
inp.bp    = 0.044;					
inp.u2    = 329.5;					
inp.tsd   = 24.2;					
inp.scf   = 0.1818;					

% ##  DATASET 13F: AERODYNAMIC RESISTANCE TERM
% ##     RAVT... (S/M) AERODYNAMIC RESISTANCE AT THE SOIL SURFACE
% ##     RAVS.. (S/M) AERODYNAMIC RESISTANCE AT THE SIDE OF THE COLUMN, WARNING: BE AWARE OF THE SIZE EFFECT
% ##     SWRAT..  (-)   PARAMETER TO SWICH ON (1.D) OR OFF(0.D0) THE TEMPERATURE CHANGE ON THE SURFACE
% ##         [RAVT]    [RAVS]   [SWRAT]
           % 206.D0     50.D0     0.D0
inp.ravt   = aero_resistance;
inp.ravs   = 50.;
inp.swart  = 0.;  
          
% ##  DATASET 13G: SOIL THERMO PROPERTY TERM           
% ##     HSC  -- HEAT SOURCE FROM ATMOSPHERE G [WATT/M2]
% ##     NHER -- NUMER OF NODES THAT CONSIDERS HEAT BALANCE    
% ##     ROUS -- THE DENSITY OF SOLID GRAIN     [KG/M3]
% ##     HCS  -- HEAT CAPACITY OF SOLID GRAIN   [J/KG/K]
% ##      
% ##     [HSC]      [HER]      [ROUS]     [HCS]
      % 9.00D2        1.95D-1      2650.D0   8.75D2
inp.hsc   = 900.;
inp.her   = 1.95e-1;
inp.rous  = 2650.;
inp.hcs   = 8.75e+2;

% ##  DATASET 13H, PARAMETERS FOR THE SALT RESISTANCE
% ##    AR  FITTING PARAMETER
% ##    BR  FITTING PARAMETER   
% ##          [AR]   [BR]
         % 6.9D-1  -1.04D-0
inp.ar   = 6.9e-1;
inp.br   = -1.04;
        
% ##  DATASET 13I SOIL CHARACTERISTIC PARAMETERS
% ##  USING EACH PARAMETERS ARE CONTROLLED BY NREG AND LREG AT DATASET 14 AND DATASET 15
% ##  THE UNIT OF AA1, PHY0 AND PHYB IS [M-1]
% ##  ECTO--- ESSENTRICITY AND TORTUOSITY THE DEFAULT VALUE IS 0.5
% ##  [SWRES1] [AA1]  [VN1]  [SWRES2]    [AA2]   [VN2]  [SWRES3]  [LAM3]   [PHYB3]  [SWRES4]   [LAM4]  [PHYB4]  [PHY0]   [ECTO]
      % 0.06   14.5D0  8.5D0    0.06      15.D0   9.2D0    0.09D0   8.0D0   0.2D0   0.08D0     4.2D0   0.045    5.0D4     0.5D0
inp.swres1  = 0.06; 
inp.aa1     = 15.5;   
inp.vn1     = 8.;
inp.swres2  = 0.;
inp.aa2     = 0.;
inp.vn2     = 0.;
inp.swres3  = 0.;
inp.lam3    = 0.;
inp.phyb3   = 0.;
inp.swres4  = 0.;
inp.lam4    = 0.;
inp.phyb4   = 0.;
inp.phy0    = 5.e+4;
inp.ecto	= 0.5;

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
inp.ntc    = 2.;
inp.b1     = 1.528;
inp.b2     = -2.406;
inp.b3	   = 4.909;	  

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
inp.nef    = 0.;
inp.ya     = 9.5;
inp.yb     = 3.0;
inp.yc	   = 2.6	 ; 
inp.yd     = 1.0;
inp.ye     = 4.0;
inp.fc     = 1.e-3;
	  
% ## DATASET 13L PARAMETERS FOR SURFACE RESISTANCE
% ## SEE PAGE 142BLUEBOOK, 98(5) FOR REFERENCE
% ## [TAL]  --  THICKNESS OF AIR LAYER (M) [TAL]
% ## [EC]   --  ECCENTRICITY OF THE ACTIVE PORE, WHEN MSR=6 (DENOTING USING APPROXIMATED SOLUTION GIVEN BY PAPER 2), THIS IS THE EXPOENTIAL
% ##             OF G FUNCTION
% ## [ETR]  --  RESIDUAL EVAPORATION RATE [-] [ETR]
% ##  [TAL]   [EC]   [ETR]     [PSIP]      [CORS]
     % 0.3D-3   5.D-1   1.D0    1000.D0      1.D0
inp.tal    = 0.3e-3;
inp.ec     = 0.5;
inp.etr    = 1.;
inp.psip   = 1000.;	  
inp.cors   = 1.0	 ;
	 
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
inp.corf    = 1.0;
inp.agr     = 3.5e-4;
inp.phicm   = 1000.;
inp.asvl    = -6.0e-20	;

%% ##  DATASET 14:  NODEWISE DATA
%##                              [SCALX] [SCALY] [SCALZ] [PORFAC]
inp.scalx  = 1.;
inp.scaly  = 1.;
inp.scalz  = 1.;
inp.porfac = 1.;

%##      [II]    [NRE    G(II)]  [X(II)] [Y(II)] [Z(II)] [POR(II)]
inp.ii   = (1:nn)';
inp.nreg = zeros(nn,1)+1;
inp.x    = x_nod_array;
inp.y    = y_nod_array;
inp.z    = zeros(nn,1)+dz;
inp.por  = zeros(nn,1)+porosity;

%% ##  DATASET 15:  ELEMENTWISE DATA
%##             [PMAXFA]        [PMINFA]        [ANG1FA]        [ALMAXF]        [ALMINF]        [ATMAXF]        [ATMINF]
%'ELEMENT'     1.0000000D+00   1.0000000D+00   1.0000000D+00       1               1                1               1
inp.pmaxfa = 1.;
inp.pminfa = 1.;
inp.angfac = 1.;
inp.almaxf = 1.;
inp.alminf = 1.;
inp.atmaxf = 1.;
inp.atminf = 1.;
%##     [L]      [LREG(L)]       [PMAX(L)]       [PMIN(L)]       [ANGLEX(L)]     [ALMAX(L)]      [ALMIN(L)]      [ATMAX(L)]      [ATMIN(L)]

inp.l       =(1:ne)';
inp.lreg    =zeros(ne,1)+1;
inp.pmax    =zeros(ne,1)+permeability_sand_m2;
inp.pmin    =zeros(ne,1)+permeability_sand_m2;
inp.anglex  =zeros(ne,1);
inp.almax   =zeros(ne,1)+0.01;
inp.almin   =zeros(ne,1)+0.01;
inp.atmax   =zeros(ne,1)+0.001;
inp.atmin   =zeros(ne,1)+0.001;

% note: for SUTRASET when node number is negative, the second input is surface area of the node
% and the third input is the thickness of the cell.

%% DATASET 17   Data for Fluid Source and Sinks

inp.iqcp  = iqcp;
inp.qinc  = qinc;
inp.uinc  = uinc;


% DATASET 19:  Data for Specified Pressure Nodes
%###  [IPBC]                [PBC]                [UBC]
if constant_water_table_m ==0
	inp.ipbc = '%';
	inp.pbc  = '%';
	inp.ubc  = '%';
	
	elseif strcmp(pressure_boundary,'bottom')==1
	bottom_nodes                = node_index_mtx_gravity_compensated(ny,:)';  % below 4 metre, greater than 200 m away from the centre
	pbc(1:nx)                   = constant_water_table_m*(c.rhow_pure_water+700*c_saltwater_kgPkg)*c.g;
	ubc(1:nx)                   = c_saltwater_kgPkg;
	inp.ipbc = bottom_nodes;
	inp.pbc  = pbc';
	inp.ubc  = ubc';
	inp.npbc = length(inp.pbc);
	
	elseif strcmp(pressure_boundary,'right_side')==1
	% for i = 1:inp.nn2-1
	% [yminValue(i), yclosestIndex(i)]=min(abs(0.25 - y_ele_matrix(:,i))); %get the node closest to water table by comparing y 
	% x_water_table(i) = x_ele_matrix(yclosestIndex(i),i);
	% y_water_table(i) = y_ele_matrix(yclosestIndex(i),i);
	% end
	right_side_nodes            = node_index_mtx_gravity_compensated(:,end)';  
	right_boundary_nodes		= right_side_nodes;
	right_boundary_nodes(y_nod_mtx(right_side_nodes)>boundary_height) = [];	
	ipbc						= right_boundary_nodes;
	right_boundary_nodes_y 	    = y_nod_mtx(right_boundary_nodes);
	boundary_nodes_y			= right_boundary_nodes_y;
		if boundary_nodes_y < constant_water_table_m
			pbc            	    = (constant_water_table_m-boundary_nodes_y)*(c.rhow_pure_water+700*c_saltwater_kgPkg)*c.g;
		else 
			pbc                 = (constant_water_table_m-boundary_nodes_y)*c.rhow_pure_water*c.g;
		end
	ubc(1:length(pbc))          = c_saltwater_kgPkg;
	inp.ipbc = ipbc';
	inp.pbc  = pbc';
	inp.ubc  = ubc';
	inp.npbc = length(inp.pbc);	
	
	elseif strcmp(pressure_boundary,'left_side')==1
	left_side_nodes             = node_index_mtx_gravity_compensated(:,1)';  
	left_boundary_nodes			= left_side_nodes;
	left_boundary_nodes(y_nod_mtx(left_side_nodes)>boundary_height)=[];	
	% left_side_nodes(1)			= []; %choose if exclude the surface node
	ipbc						= left_boundary_nodes;
	left_boundary_nodes_y 		= y_nod_mtx(left_boundary_nodes);
	boundary_nodes_y			= left_boundary_nodes_y;
		if boundary_nodes_y < constant_water_table_m
			pbc            		= (constant_water_table_m-boundary_nodes_y)*(c.rhow_pure_water+700*c_saltwater_kgPkg)*c.g;
		else 
			pbc             	= (constant_water_table_m-boundary_nodes_y)*c.rhow_pure_water*c.g;
		end
	ubc(1:length(pbc))          = c_saltwater_kgPkg;
	inp.ipbc = ipbc';
	inp.pbc  = pbc';
	inp.ubc  = ubc';
	inp.npbc = length(inp.pbc);	
	
	elseif strcmp(pressure_boundary,'both_side')==1
	left_side_nodes             = node_index_mtx_gravity_compensated(:,1)';  
	left_boundary_nodes			= left_side_nodes;
	left_boundary_nodes(y_nod_mtx(left_side_nodes)>boundary_height)=[];	
	right_side_nodes            = node_index_mtx_gravity_compensated(:,end)';  	
	right_boundary_nodes		= right_side_nodes;
	right_boundary_nodes(y_nod_mtx(right_side_nodes)>boundary_height) = [];	
	ipbc						= [left_boundary_nodes,right_boundary_nodes];
	left_boundary_nodes_y 			= y_nod_mtx(left_boundary_nodes);
	right_boundary_nodes_y 			= y_nod_mtx(right_boundary_nodes);
	boundary_nodes_y				= [left_boundary_nodes_y,right_boundary_nodes_y];
		if boundary_nodes_y < constant_water_table_m
			pbc            	    = (constant_water_table_m-boundary_nodes_y)*(c.rhow_pure_water+700*c_saltwater_kgPkg)*c.g;
		else 
			pbc                 = (constant_water_table_m-boundary_nodes_y)*c.rhow_pure_water*c.g;
		end
	ubc(1:length(pbc))          = c_saltwater_kgPkg;
	inp.ipbc = ipbc';
	inp.pbc  = pbc';
	inp.ubc  = ubc';
	inp.npbc = length(inp.pbc);	
end  

%mask_mtx_aquifer_boundary_gravity_compensated_left = and(y_nod_mtx_gravity_compensated_m<-4, x_nod_mtx_gravity_compensated_m<0.01);
%ipbc_node_idx_array_left                           = node_index_mtx_gravity_compensated(mask_mtx_aquifer_boundary_gravity_compensated_left);
%pbc_left                                      = -(y_nod_mtx_gravity_compensated_m(mask_mtx_aquifer_boundary_gravity_compensated) - initial_head_aquifer_m +0.5 ) *c.g * inp.rhow0 ;
%
%inp.ipbc = [ipbc_node_idx_array; ipbc_node_idx_array_left];
%inp.pbc  = [pbc;pbc_left];
%inp.ubc  = [zeros(size(pbc))+c_saltwater_kgPkg;zeros(size(pbc_left))+c_freshwater_kgPkg];
%inp.npbc = length(inp.pbc);

%##
%##  DATASET     22:  Ele        ment Incid      ence Data
%##    [LL]      [IIN(1)]        [IIN(2)]        [IIN(3)]        [IIN(4)]

%ne_mtx=reshape(l,ney,nex);
inp.iin1 = zeros(inp.ne,1);
inp.iin2 = zeros(inp.ne,1);
inp.iin3 = zeros(inp.ne,1);
inp.iin4 = zeros(inp.ne,1);

% DATASET 22, user need to check if the sequence of the node is in clockwise manner as suggested by the manual

inp.iin1= iin1;
inp.iin2= iin2;
inp.iin3= iin3;
inp.iin4= iin4;

inp.export_to_file();

%% ICS FILE



% setting the initial pressure as hydrostatic, in particular at the sandy aquifer, the silt layer will be overwritten
pm1_mtx_gravity_pa= - (- initial_head_aquifer_m + y_nod_mtx)*c.g*(c.rhow_pure_water+700*initial_concentration_kgPkg);%sutra manual p15


% mask_nod_mtx_silt_layer_gravity_compensated = y_nod_mtx_gravity_compensated_m  > -6  ;   % mask matrix, for nod matrix 

%set a relatively high suction in the silt layer
% pm1_mtx_gravity_compensated_pa ( mask_nod_mtx_silt_layer_gravity_compensated) = -20000;

% put the top cell as zero pressure 
% pm1_mtx_gravity_compensated_pa(mask_nod_mtx_aquifer_boundary_gravity_compensated_top)=initial_pond_water_depth_m*c.rhow_pure_water*c.g;
% pm1_mtx_pa=flip(pm1_mtx_gravity_compensated_pa);


% initial solute concentrtion, sandy aquifer has a saline concentration while silt has a fresh concentration
% um1_mtx_gravity_compensated_kgPkg= zeros(size(pm1_mtx_gravity_compensated_pa))+ c_saltwater_kgPkg;
% um1_mtx_gravity_compensated_kgPkg (mask_nod_mtx_silt_layer_gravity_compensated) = c_freshwater_kgPkg;
% um1_mtx_kgPkg=flip(um1_mtx_gravity_compensated_kgPkg);


fprintf('use the original ics file')
% ics file
ics       = icsObj('FLUME','read_from_file','no');
ics.tics  = 0.0;
ics.cpuni = 'NONUNIFORM';
ics.pm1   = pm1_mtx_gravity_pa;
% ics.cpuni = 'UNIFORM';
% ics.pm1   = 10;
ics.cuuni = 'UNIFORM';
ics.um1   = initial_concentration_kgPkg ;
ics.ctuni = 'UNIFORM';
ics.tm1   = initial_temperature_C + 273.15 ;

ics.export_to_file();
%
%
%%pm1_mtx=reshape(b.pm1,[inp.nn1,inp.nn2]);
%%
%%b=icsObj('PART1_cs.csv','inpObj',inp);
%%pm1_mtx=reshape(b.pm1,[inp.nn1,inp.nn2]);
%%um1_mtx=reshape(b.um1,[inp.nn1,inp.nn2]);
%%inp.get_x_nod_mtx;
%%inp.get_y_nod_mtx;
%%a=figure;
%%surf(inp.x_nod_mtx,inp.y_nod_mtx,um1_mtx);
%%savefig(a,'conc.fig');
%%
%%a=figure;
%%surf(inp.x_nod_mtx,inp.y_nod_mtx,pm1_mtx);
%%savefig(a,'p.fig');
%%
%%saveas(a,'Barchart.png')
%
%
%ics.export_to_file();
