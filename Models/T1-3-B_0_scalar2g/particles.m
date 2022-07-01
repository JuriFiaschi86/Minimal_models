ParticleDefinitions[GaugeES] = {
	(* new fields *)
	{CPsi0, {Description -> "BSM field Ψ0",
			OutputName -> "CPsi0",
			LaTeX -> "\\Psi^0"
		}
	},
	{psipp, {Description -> "BSM field ψ'+",
			OutputName -> "psipp",
			LaTeX -> "\\psi'^+"
		}
	},
	{psip0, {Description -> "BSM field ψ'0",
			OutputName -> "psip0",
			LaTeX -> "\\psi'^0"
		}
	},
	{phip, {Description -> "BSM field ϕ+",
			OutputName -> "phip",
			LaTeX -> "\\phi^+",
			ElectricCharge -> 1
		}
	},
	{phi0, {Description -> "BSM field ϕ0",
			OutputName -> "phi0",
			LaTeX -> "\\phi^0",
			ElectricCharge -> 0
		}
	},
	{psi0, {Description -> "BSM field ψ0",
			OutputName -> "psi0",
			LaTeX -> "\\psi^0"
		}
	},
	{psim, {Description -> "BSM field ψ-",
			OutputName -> "psim",
			LaTeX -> "\\psi^-"
		}
	},

	(* Standard Model *)
	{H0, {
			PDG -> {0},
			Width -> 0,
			Mass -> Automatic,
			FeynArtsNr -> 1,
			LaTeX -> "H^0",
			OutputName -> "H0"
		}
	},
	{Hp, {
			PDG -> {0},
			Width -> 0,
			Mass -> Automatic,
			FeynArtsNr -> 2,
			LaTeX -> "H^+",
			OutputName -> "Hp"
		}
	},

	{VB,  {Description -> "B-Boson"}},
	{VG,  {Description -> "Gluon"}},
	{VWB, {Description -> "W-Bosons"}},
	{gB,  {Description -> "B-Boson Ghost"}},
	{gG,  {Description -> "Gluon Ghost"}},
	{gWB, {Description -> "W-Boson Ghost"}}
};

ParticleDefinitions[EWSB] = {
	(* new fields *)
	(* to be completed manually *)
	{etap, {Description -> "BSM field η+",
			OutputName -> "etp",
			LaTeX -> "\\eta^+",
			PDG -> {921, 922},
			FeynArtsNr -> 921,
			ElectricCharge -> 1
		}
	},
	{eta0, {Description -> "BSM neutral field η0",
			OutputName -> "et0",
			LaTeX -> "\\eta^0",
			PDG -> {931, 932},
			FeynArtsNr -> 931,
			ElectricCharge -> 0
		}
	},

	{Fpsi, {Description -> "BSM charged (-1) Dirac spinor Fψ",
			OutputName -> "psi",
			LaTeX -> "\\psi",
			PDG -> {903},
			FeynArtsNr -> 903,
			ElectricCharge -> -1
		}
	},
	{Fchi, {Description -> "BSM neutral Majorana spinor Fχ",
			OutputName -> "ch",
			LaTeX -> "\\chi",
			PDG -> {911, 912, 913},
			FeynArtsNr -> 911,
			ElectricCharge -> 0
		}
	},

	{Fnu,  {Description -> "Neutrinos",
			Mass -> LesHouches
		}
	},

	(* Standard Model *)
	{hh,   {Description -> "Higgs",
			PDG -> {25},
			PDG.IX -> {101000001}
		}
	},
	{Ah,   {Description -> "Pseudo-Scalar Higgs",
			PDG -> {0},
			PDG.IX -> {0},
			Mass -> {0},
			Width -> {0}
		}
	},
	{Hp,   {Description -> "Charged Higgs",
			PDG -> {0},
			PDG.IX -> {0},
			Width -> {0},
			Mass -> {0},
			LaTeX -> {"H^+", "H^-"},
			OutputName -> {"Hp", "Hm"},
			ElectricCharge -> 1
		}
	},

	{VP,   {Description -> "Photon"}},
	{VZ,   {Description -> "Z-Boson",
			Goldstone -> Ah
		}
	},
	{VG,   {Description -> "Gluon"}},
	{VWp,  {Description -> "W+ - Boson",
			Goldstone -> Hp
		}
	},
	{gP,   {Description -> "Photon Ghost"}},
	{gWp,  {Description -> "Positive W+ - Boson Ghost"}},
	{gWpC, {Description -> "Negative W+ - Boson Ghost"}},
	{gZ,   {Description -> "Z-Boson Ghost"}},
	{gG,   {Description -> "Gluon Ghost"}},

	{Fd,   {Description -> "Down-Quarks"}},
	{Fu,   {Description -> "Up-Quarks"}},
	{Fe,   {Description -> "Leptons"}}
};

WeylFermionAndIndermediate = {
	(* new fields *)
	{phip, {LaTeX -> "\\phi^+"}},
	{phi0, {LaTeX -> "\\phi^0"}},

	{CPsi, {LaTeX -> "\\CPsi"}},
	{psip, {LaTeX -> "\\psi'"}},
	{phi, {LaTeX -> "\\phi"}},

	{chi, {LaTeX -> "\\chi"}},
	(* Standard Model *)
	{H, {
			PDG -> {0},
			Width -> 0,
			Mass -> Automatic,
			LaTeX -> "H",
			OutputName -> ""
		}
	},

	{dR,  {LaTeX -> "d_R"}},
	{eR,  {LaTeX -> "e_R"}},
	{l,   {LaTeX -> "l"}},
	{uR,  {LaTeX -> "u_R"}},
	{q,   {LaTeX -> "q"}},
	{eL,  {LaTeX -> "e_L"}},
	{dL,  {LaTeX -> "d_L"}},
	{uL,  {LaTeX -> "u_L"}},
	{vL,  {LaTeX -> "\\nu_L"}},

	{DR,  {LaTeX -> "D_R"}},
	{ER,  {LaTeX -> "E_R"}},
	{UR,  {LaTeX -> "U_R"}},
	{EL,  {LaTeX -> "E_L"}},
	{DL,  {LaTeX -> "D_L"}},
	{UL,  {LaTeX -> "U_L"}}
};
