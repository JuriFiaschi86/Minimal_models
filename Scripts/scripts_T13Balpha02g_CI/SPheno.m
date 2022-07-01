OnlyLowEnergySPheno = True;
(* add Block TREELEVELUNITARITY to SPheno output spectrum file *)
(* AddTreeLevelUnitarityLimits = True; *)

MINPAR = {
    {1, lambdaInput},
    (*{2, mHInput},*)
    {12, mCPsiInput},
    {13, mpsipsipInput},
    {24, lambda4Input},
    {25, lambda5Input}
};

BoundaryLowScaleInput = {
    {mϕ2, LHInput[mϕ2]},
    {mΨ, mCPsiInput},
    {mψψp, mpsipsipInput},
    {λ1, LHInput[λ1]},
    (*{λ2, LHInput[λ2]},*)
    (*{λ3, LHInput[λ3]},*)
    {λ4, lambda4Input},
    {λ5, lambda5Input},
    {λ6, LHInput[λ6]},
    (* Standard Model *)
    {λ, lambdaInput}
    (*{λ, mHInput^2 / vSM^2}*)
};

DEFINITION[MatchingConditions] = {
    {Ye, YeSM},
    {Yd, YdSM},
    {Yu, YuSM},
    {g1, g1SM},
    {g2, g2SM},
    {g3, g3SM},
    {v, vSM}
};

DefaultInputValues = {
    mϕ2[1, 1] -> 510^2,
    mϕ2[2, 2] -> 520^2,
    mCPsiInput -> 610,
    mpsipsipInput -> 720,
    lambda4Input -> 0.014,
    lambda5Input -> 0.015,
    λ6[1, 1] -> 0.016,
    λ6[2, 1] -> 0.016,
    λ6[3, 1] -> 0.016,
    λ6[1, 2] -> 0.016,
    λ6[2, 2] -> 0.016,
    λ6[3, 2] -> 0.016,
    (* Standard Model *)
    (*mHInput -> 125.09,*) (* see PDG *)
    lambdaInput -> 0.2612 (* ≈ 125.09² / 244.74² *)
};

ParametersToSolveTadpoles = {mu2};

ListDecayParticles = {Fu, Fe, Fd, hh, Hp, etap, eta0, Fchi, Fpsi};
ListDecayParticles3B = {{Fu, "Fu.f90"}, {Fe, "Fe.f90"}, {Fd, "Fd.f90"}};

(* see http://stauby.de/sarah_userforum/viewtopic.php?f=4&t=378 *)
(*RenConditionsDecays = {
    {dCosTW, 1/2*Cos[ThetaW] * (PiVWp/(MVWp^2) - PiVZ/(mVZ^2))},
    {dSinTW, -dCosTW/Tan[ThetaW]},
    {dg2, 1/2 * g2 * (
        derPiVPheavy0 + PiVPlightMZ / MVZ^2 - (
            -(PiVWp / MVWp^2) + PiVZ / MVZ^2) / Tan[ThetaW]^2 + (2 * PiVZVP * Tan[ThetaW]
        ) / MVZ^2
    )},
    {dg1, dg2 * Tan[ThetaW] + g2 * dSinTW / Cos[ThetaW] - dCosTW * g2 * Tan[ThetaW] / Cos[ThetaW]}
};*)
