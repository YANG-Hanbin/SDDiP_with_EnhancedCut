param = (
    verbose             = false,
    MIPGap              = 1e-4,
    TimeLimit           = 3600,
    terminate_threshold = 1e-3,
    MaxIter             = 100,
    θ̲                   = 0.0,
    OPT                 = 0.0,
    tightness           = true,
    numScenarios        = 3,
    LiftIterThreshold   = 10,
    branch_threshold    = 1e-3,
    ## "interval_mid", "exact_point"
    med_method          = "interval_mid",   
    ## :PLC, :SMC, :LC
    cutSelection        = :PLC,             
    ## :SDDPL, :SDDP, :SDDiP
    algorithm           = :SDDPL,   
    T                   = 12, 
    num                 = 10, 
    case                = "case30pwl", # "case_RTS_GMLC", "case30", "case30pwl", "case24_ieee_rts"     
)


param_levelsetmethod = (
    μ             = 0.9,
    λ             = 0.5,
    threshold     = 1e-4,
    nxt_bound     = 1e10,
    MaxIter       = 200,
    verbose       = false,
)

param_PLC = (
    core_point_strategy = "Eps", # "Mid", "Eps"
    δ                   = 1e-3,
    ℓ                   = .0,
)