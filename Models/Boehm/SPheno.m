OnlyLowEnergySPheno = True;

MINPAR = {
    {1, lambdaInput},
    (*{2, mHInput},*)
    {11, M2Input},
    {13, M3Input},
    {14, M4Input},
    (*{15, mNInput},*)
    {22, lambda2Input},
    {23, lambda3Input},
    {24, lambda4Input},
    {25, lambda5Input},
    {26, lambda6Input},
    {27, lambda7Input},
    {28, YbsmInput}
};

BoundaryLowScaleInput = {
    {M22, M2Input^2},
    {Mtst, M3Input^2},
    {M42, M4Input^2},
    {MasN, LHInput[MasN]},
    {λ2, lambda2Input},
    {λ3, lambda3Input},
    {λ4, lambda4Input},
    {λ5, lambda5Input},
    {λ6, lambda6Input},
    {λ7, lambda7Input},
    {λ8, LHInput[λ8]},
    {Ybsm, YbsmInput},
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
    M2Input -> 100,
    M3Input -> 150,
    M4Input -> 160,
    lambda2Input -> 0.012,
    lambda3Input -> 0.013,
    lambda4Input -> 0.014,
    lambda5Input -> 0.015,
    lambda6Input -> 0.016,
    lambda7Input -> 0.017,
    YbsmInput -> 0.019,
    λ8[1, 1] -> 0.029,
    λ8[2, 2] -> 0.029,
    λ8[3, 3] -> 0.029,
    MasN[1, 1] -> 20,
    MasN[2, 2] -> 20,
    MasN[3, 3] -> 20,
    (* Standard Model *)
    (*mHInput -> 125.09,*) (* see PDG *)
    lambdaInput -> 0.2612 (* ≈ 125.09² / 244.74² *)
};

ParametersToSolveTadpoles = {mu2};

ListDecayParticles = {Fu, Fe, Fd, hh, Hp, etap, zetaR, zetaI, FmaN};
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
