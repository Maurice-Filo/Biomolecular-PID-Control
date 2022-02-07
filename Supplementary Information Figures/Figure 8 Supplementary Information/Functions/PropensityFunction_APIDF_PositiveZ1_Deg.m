function Propensity = PropensityFunction_APIDF_PositiveZ1_Deg(x, Parameters)
% Propensity Function for Cholesterol Network controlled by the APIDF
% controller with degradation inhibition, and Z1 as the positive actuation.

%% Parameters
Ingestion.k1 = Parameters.Ingestion.k1;
ICSmax = Parameters.ICSmax;
ICt = Parameters.ICt;
IS = Parameters.IS;
BileSaltRelease.k1 = Parameters.BileSaltRelease.k1;
BileSaltReturn.k1 = Parameters.BileSaltReturn.k1;
BileSaltExcretion.k1 = Parameters.BileSaltExcretion.k1;
k5 = Parameters.k5;
k6 = Parameters.k6;
k7 = Parameters.k7;
k8 = Parameters.k8;
BCRmax = Parameters.BCRmax;
BCRt = Parameters.BCRt;
BS = Parameters.BS;
HCSmax = Parameters.HCSmax;
HCSt = Parameters.HCSt;
HS = Parameters.HS;
k9 = Parameters.k9;
k10 = Parameters.k10;
k11 = Parameters.k11;
VLDLCholesterolFormation.k1 = Parameters.VLDLCholesterolFormation.k1;
khrs = Parameters.khrs;
HepaticLDLReceptorDegradation.k1 = Parameters.HepaticLDLReceptorDegradation.k1;
VLDLCholesterolReUptake.k1 = Parameters.VLDLCholesterolReUptake.k1;
k15 = Parameters.k15;
IDLCholesterolReUptake.k1 = Parameters.IDLCholesterolReUptake.k1;
k17 = Parameters.k17;
k18 = Parameters.k18;
ReceptorIndependentHepaticUptake.k1 = Parameters.ReceptorIndependentHepaticUptake.k1;
k20 = Parameters.k20;
ReceptorIndependentPeripheralUptake.k1 = Parameters.ReceptorIndependentPeripheralUptake.k1;
kprs = Parameters.kprs;
PeripheralLDLReceptorDegradation.k1 = Parameters.PeripheralLDLReceptorDegradation.k1;
k23 = Parameters.k23;
k24 = Parameters.k24;
PeripheralSteroidProduction.k1 = Parameters.PeripheralSteroidProduction.k1;
k26 = Parameters.k26;
PCSmax = Parameters.PCSmax;
PPCt = Parameters.PPCt;
PCSS = Parameters.PCSS;
k27 = Parameters.k27;
k28 = Parameters.k28;
k29 = Parameters.k29;
% Intake = Parameters.Intake;
% Intestine = Parameters.Intestine;
% HepaticTissue = Parameters.HepaticTissue;
% PeripheralTissue = Parameters.PeripheralTissue;
% Plasma = Parameters.Plasma;
% Excreted = Parameters.Excreted;
mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
k = Parameters.k;
delta_c = Parameters.delta_c;
delta = Parameters.delta;
n = Parameters.n;
beta = Parameters.beta;

%% Extract State Variables
DC = x(1);
IC = x(2);
ICS = x(3);
IBS = x(4);
HBS = x(5);
HFC = x(6);
HCS = x(7);
HCE = x(8);
ACAT = x(9);
CEH = x(10);
HNHDLS = x(11);
HLDLRs = x(12);
HLDLRsS = x(13);
HLDLRD = x(14);
SRB1 = x(15);
PFC = x(16);
PLDLRs = x(17);
PLDLRsS = x(18);
PLDLRD = x(19);
PCE = x(20);
PSS = x(21);
PCS = x(22);
INHDLS = x(23);
NHDL = x(24);
VLDLC = x(25);
IDLC = x(26);
LPL = x(27);
LDLC = x(28);
HSL = x(29);
HDLC = x(30);
LCAT = x(31);
CETP = x(32);
EBS = x(33);
EC = x(34);
Z_1 = x(35);
Z_2 = x(36);

%% Propensities
Propensity = [ ...
    Ingestion.k1*DC; ...
   	ICSmax/(1+power(IC/ICt,IS)); ...
  	BileSaltRelease.k1*HBS; ...
  	BileSaltReturn.k1*IBS; ...
   	BileSaltExcretion.k1*IBS; ...
  	k5*HFC/HBS; ...
  	k6*IC*IBS; ...
  	k7*IC*IBS; ...
  	k8*PFC; ...
   	BCRmax/(1+power(BCRt/HFC,BS)); ...
 	HCSmax/(1+power(HFC/HCSt,HS)); ...
	k9*ACAT*HFC; ...
  	k10*CEH*HCE; ...
  	k11*PFC; ...
  	VLDLCholesterolFormation.k1*HFC; ...
   	khrs*HLDLRsS/HFC; ...
  	HepaticLDLReceptorDegradation.k1*HLDLRs; ...
   	VLDLCholesterolReUptake.k1*VLDLC; ...
  	k15*VLDLC*LPL; ...
  	IDLCholesterolReUptake.k1*IDLC; ...
  	k17*IDLC*HSL; ...
  	k18*LDLC*HLDLRs; ...
   	ReceptorIndependentHepaticUptake.k1*LDLC; ...
   	k20*PLDLRs*LDLC; ...
  	ReceptorIndependentPeripheralUptake.k1*LDLC; ...
   	kprs*PLDLRsS/PFC; ...
  	PeripheralLDLReceptorDegradation.k1*PLDLRs; ...
   	k23*ACAT*PFC; ...
  	k24*CEH*PCE; ...
   	PeripheralSteroidProduction.k1*PFC; ...
	k26*PFC*NHDL*LCAT; ...
   	PCSmax/(1+power(PFC/PPCt,PCSS)); ...
	k27*HDLC*CETP; ...
  	k28*HDLC*CETP; ...
  	k29*HDLC*SRB1; ...
    mu + beta*LDLC; ...
    theta*LDLC; ...
    eta*Z_1*Z_2; ...
    k*Z_1; ...
    delta_c*Z_1; ...
    delta_c*Z_2; ... 
    delta*LDLC^n*IC ...
    ];
end