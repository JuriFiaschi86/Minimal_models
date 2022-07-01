Model`Name = "T13Balpha02g";
Model`NameLaTeX = "T1-3-B (\\alpha = 0) with 2 scalar generations";
Model`Authors = "Simon May";
Model`Date = "2018-04-15";

(*-------------------------------------------*)
(*   particle content                        *)
(*-------------------------------------------*)

(* global symmetries *)
(* discrete ℤ₂ symmetry *)
Global[[1]] = {Z[2], Z2};

(* gauge groups *)
Gauge[[1]] = {B,   U[1], hypercharge, g1, False, 1};
Gauge[[2]] = {WB, SU[2], left,        g2, True,  1};
Gauge[[3]] = {G,  SU[3], color,       g3, False, 1};

(* matter fields *)
(*                   {name, gens, components, Y/2,  SU(2), SU(3), global} *)
(* Standard Model *)
FermionFields[[1]] = {q,    3,    {uL, dL},    1/6, 2,     3,     1};
FermionFields[[2]] = {l,    3,    {vL, eL},   -1/2, 2,     1,     1};
FermionFields[[3]] = {d,    3,    conj[dR],    1/3, 1,    -3,     1};
FermionFields[[4]] = {u,    3,    conj[uR],   -2/3, 1,    -3,     1};
FermionFields[[5]] = {e,    3,    conj[eR],      1, 1,     1,     1};

ScalarFields[[1]]  = {H,    1,    {Hp, H0},    1/2, 2,     1,     1};

(* new fields *)
FermionFields[[6]] = {CPsi, 1, CPsi0, 0, 1, 1, -1};
FermionFields[[7]] = {psip, 1, {psipp, psip0}, 1/2, 2, 1, -1};
FermionFields[[8]] = {psi, 1, {psi0, psim}, -1/2, 2, 1, -1};
ScalarFields[[2]]  = {phi, 2, {{phi0/Sqrt[2], phip}, {conj[phip], -phi0/Sqrt[2]}}, 0, 3, 1, -1};
RealScalars = {phi0};


(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)
NameOfStates = {GaugeES, EWSB};

(* ----- before EWSB ----- *)
DEFINITION[GaugeES][LagrangianInput] = {
	{LagHC,      {AddHC -> True }},
	{LagNoHC,    {AddHC -> False}},
	{LagBSMHC,   {AddHC -> True }},
	{LagBSMNoHC, {AddHC -> False}}
};

(* Standard Model Lagrangian *)
LagHC      = -Yd conj[H].d.q - Ye conj[H].e.l - Yu u.q.H;
LagNoHC    = mu2 conj[H].H - 1/2 λ conj[H].H.conj[H].H;

(* Lagrangian involving the new fields *)
(* can’t include λ3 term with four indices – see
http://stauby.de/sarah_userforum/viewtopic.php?f=4&t=402 *)
LagBSMNoHC = - 1/2 mϕ2 phi.phi \
	 - λ1 conj[H].H.phi.phi (*- λ2 conj[H].phi.phi.H*) (*- λ3 phi.phi.phi.phi*);
LagBSMHC   = - 1/2 mΨ CPsi.CPsi - mψψp psi.psip \
	 - λ4 conj[H].psip.CPsi - λ5 psi.H.CPsi - λ6 l.phi.psip;

(* ----- after EWSB ----- *)
(* gauge sector mixing *)
DEFINITION[EWSB][GaugeSector] = {
	{{VB,     VWB[3]}, {VP, VZ}, ZZ},
	{{VWB[1], VWB[2]}, {VWp, conj[VWp]}, ZW}
};

(* VEVs *)
DEFINITION[EWSB][VEVs] = {
	(* Standard Model Higgs VEV *)
	{H0,
		{v, 1/Sqrt[2]},
		{Ah, \[ImaginaryI]/Sqrt[2]},
		{hh, 1/Sqrt[2]}
	}
};

(* mixing *)
DEFINITION[EWSB][MatterSector] = {
	(* Standard Model mixing *)
	{{{dL}, {conj[dR]}}, {{DL, Vd}, {DR, Ud}}},
	{{{uL}, {conj[uR]}}, {{UL, Vu}, {UR, Uu}}},
	{{{eL}, {conj[eR]}}, {{EL, Ve}, {ER, Ue}}},
	(* mixing of new fields can be added here *)
	{{vL}, {VL, Uneu}},
	{{phip}, {etap, Uetp}},
	{{phi0}, {eta0, Uet0}},
	{{CPsi0, psi0, psip0}, {chi, Uferm}}
};

(* Dirac spinors *)
DEFINITION[EWSB][DiracSpinors] = {
	Fd -> {DL, conj[DR]},
	Fe -> {EL, conj[ER]},
	Fu -> {UL, conj[UR]},
	Fnu -> {VL, conj[VL]},
	(* new fields *)
	Fpsi -> {psim, conj[psipp]},
	Fchi -> {chi, conj[chi]}
};

(* phases, see http://stauby.de/sarah_wiki/index.php?title=Phases *)
DEFINITION[EWSB][Phases] = {
    {psim, phasepsi}
};
