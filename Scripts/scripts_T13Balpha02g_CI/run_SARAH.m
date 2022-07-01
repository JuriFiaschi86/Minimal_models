SARAHdirectory = "/home/users/fiaschi/Dokumente/SARAH-4.13.0"
model = "T1-3-B_0_scalar2g"

SetDirectory[SARAHdirectory]

<<"SARAH.m";

Print["SARAH loaded"]

Start[model];

Print["Model loaded"]

(* MakeSPheno[StandardCompiler -> "gfortran", IncludeLoopDecays -> True]; *)
MakeSPheno[];

Print["SPheno source code generated"]

(* CPViolation must be set to True if the model has complex parameters *)
MakeCHep[SLHAinput -> True, CPViolation -> False];

Print["CalcHep model files generated"]

Exit[]
