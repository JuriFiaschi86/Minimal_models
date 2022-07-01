ParameterDefinitions = {
	(* new parameters *)
	{mϕ2, {Description -> "BSM parameter mϕ^2",
			LaTeX -> "m_\\phi^2",
			OutputName -> "mph2",
			LesHouches -> "MPHI2",
			Real -> True
		}
	},
	{mΨ, {Description -> "BSM parameter mΨ",
			LaTeX -> "m_\\Psi",
			OutputName -> "mCPsi",
			LesHouches -> {"T13B", 12},
			Real -> True (* not required (Real -> False: CP violation) *)
		}
	},
	{mψψp, {Description -> "BSM parameter mψψ'",
			LaTeX -> "m_{\\psi\\psi'}",
			OutputName -> "mpp",
			LesHouches -> {"T13B", 13},
			Real -> True (* not required (Real -> False: CP violation) *)
		}
	},
	{λ1, {Description -> "BSM parameter λ1",
			LaTeX -> "\\lambda_1",
			OutputName -> "lam1",
			LesHouches -> "LAM1",
			Real -> True
		}
	},
	{λ2, {Description -> "BSM parameter λ2",
			LaTeX -> "\\lambda_2",
			OutputName -> "lam2",
			LesHouches -> "LAM2",
			Real -> True
		}
	},
	{λ3, {Description -> "BSM parameter λ3",
			LaTeX -> "\\lambda_3",
			OutputName -> "l3",
			LesHouches -> "LAM3",
			Real -> True
		}
	},
	{λ4, {Description -> "BSM parameter λ4",
			LaTeX -> "\\lambda_4",
			OutputName -> "lam4",
			LesHouches -> {"T13B", 24},
			Real -> True (* not required (Real -> False: CP violation) *)
		}
	},
	{λ5, {Description -> "BSM parameter λ5",
			LaTeX -> "\\lambda_5",
			OutputName -> "lam5",
			LesHouches -> {"T13B", 25},
			Real -> True (* not required (Real -> False: CP violation) *)
		}
	},
	{λ6, {Description -> "BSM parameter λ6",
			LaTeX -> "\\lambda_6",
			OutputName -> "lam6",
			LesHouches -> "LAM6",
			Real -> True (* not required (Real -> False: CP violation) *)
		}
	},
	{phasepsi, {Description -> "BSM phase for fermion component ψ-",
			LaTeX -> "p",
			OutputName -> "ppsi",
			LesHouches -> {"T13B", 30}
		}
	},

	(* mixing matrices must have the same Mathematica symbol and OutputName? See
	http://stauby.de/sarah_userforum/viewtopic.php?f=4&t=403 *)
	{Uetp, {Description -> "BSM charged (+1) scalar mixing matrix",
			LaTeX -> "U_{\\eta^+}",
			OutputName -> "Uetp",
			LesHouches -> "ETAPMIX"
		}
	},
	{Uet0, {Description -> "BSM neutral scalar mixing matrix",
			LaTeX -> "U_{\\eta^0}",
			OutputName -> "Uet0",
			LesHouches -> "ETA0MIX"
		}
	},
	{Uferm, {Description -> "BSM neutral fermion mixing matrix",
			LaTeX -> "U_\\chi",
			OutputName -> "Ufrm",
			LesHouches -> "CHIMIX"
		}
	},

	{Uneu, {Description -> "Neutrino mixing matrix",
			LaTeX -> "U_\\nu",
			OutputName -> "Uneu",
			LesHouches -> "NUMIX"
		}
	},

	(* Standard Model parameters *)
	{g1,        {Description -> "Hypercharge-Coupling"}},
	{g2,        {Description -> "Left-Coupling"}},
	{g3,        {Description -> "Strong-Coupling"}},
	{AlphaS,    {Description -> "Alpha Strong"}},
	{e,         {Description -> "electric charge"}},

	{Gf,        {Description -> "Fermi's constant"}},
	{aEWinv,    {Description -> "inverse weak coupling constant at mZ"}},

	{Yu,        {Description -> "Up-Yukawa-Coupling",
					DependenceNum -> Sqrt[2]/v * {
						{Mass[Fu,1], 0, 0},
						{0, Mass[Fu,2], 0},
						{0, 0, Mass[Fu,3]}
					}
				}
	},
	{Yd,        {Description -> "Down-Yukawa-Coupling",
					DependenceNum -> Sqrt[2]/v * {
						{Mass[Fd,1], 0, 0},
						{0, Mass[Fd,2], 0},
						{0, 0, Mass[Fd,3]}
					}
				}
	},
	{Ye,        {Description -> "Lepton-Yukawa-Coupling",
					DependenceNum -> Sqrt[2]/v * {
						{Mass[Fe,1], 0, 0},
						{0, Mass[Fe,2], 0},
						{0, 0, Mass[Fe,3]}
					}
				}
	},


	{mu2,       {Description -> "SM Mu Parameter"}},
	{\[Lambda], {Description -> "SM Higgs Selfcouplings",
					DependenceNum -> Mass[hh]^2/(v^2)
				}
	},
	{v,         {Description -> "EW-VEV",
					DependenceNum -> Sqrt[4 * Mass[VWp]^2/(g2^2)],
					DependenceSPheno -> None
				}
	},
	{mH2,       {Description -> "SM Higgs Mass Parameter"}},

	{ThetaW,    {Description -> "Weinberg-Angle",
					DependenceNum -> ArcSin[Sqrt[1 - Mass[VWp]^2/Mass[VZ]^2]]
				}
	},

	{ZZ,        {Description -> "Photon-Z Mixing Matrix"}},
	{ZW,        {Description -> "W Mixing Matrix",
					Dependence -> 1/Sqrt[2] {
						{1, 1},
						{\[ImaginaryI], -\[ImaginaryI]}
					}
				}
	},

	{Vu,        {Description -> "Left-Up-Mixing-Matrix"}},
	{Vd,        {Description -> "Left-Down-Mixing-Matrix"}},
	{Uu,        {Description -> "Right-Up-Mixing-Matrix"}},
	{Ud,        {Description -> "Right-Down-Mixing-Matrix"}},
	{Ve,        {Description -> "Left-Lepton-Mixing-Matrix"}},
	{Ue,        {Description -> "Right-Lepton-Mixing-Matrix"}}
};
