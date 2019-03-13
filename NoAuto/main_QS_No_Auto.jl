################################################################################
# QS without autoregulation
################################################################################
# coop_Cost: constant for cost of cooperation 
# sig_Cost: constant for cost of signaling 
# lambda: parameter used in the zero-truncated Poisson distribution for
#         generating number of mixing genotypes 
# mu_Cheats: parameter used in the zero-truncated Poisson distribution for
#            generating number of cheaters 
# rep_Num: index number of replications

function main_QS_No_Auto(coop_Cost::Float64,sig_Cost::Float64,lambda::Float64,mu_Cheats::Float64,rep_Num::Int64)

    # set path
    dirPath = "."

    include("$dirPath/fortune_wheel.jl")
    include("$dirPath/mut_parameter.jl")
    include("$dirPath/eval_genotype.jl")
    include("$dirPath/sample_ztp.jl")
    include("$dirPath/randNum_Poisson.jl")

    # reseed random generator Base.Random.GLOBAL_RNG.seed
    srand(rand(RandomDevice(), UInt32, 1000))

    ################################################################################
    # set parameters
    ################################################################################
    # maximum generation
    max_G = 5000
    # population size
    size_Pop = 5000
    # environments
    grid_Size = 100
    # baseline fitness
    baseline = 100.
    # payoff for gene turned 'ON' & cooperation
    coop_Benefit = 1.5

    # mutation rate
    mu_Production = 0.01
    mu_Th_Signal = 0.01

    # maximum cellular density (cells per microliter)
    max_CellDen = 10.^5
    # minimum cellular density
    min_CellDen = 10.^1.5
    # median cellular density
    median_CellDen = median([min_CellDen,max_CellDen])
    # baseline volume
    base_Volume = 10.

    # signal decay rate
    decay_Rate = 10.^-4.

    # initial testing environments
    env_CellDen = collect(linspace(min_CellDen,max_CellDen,grid_Size))

    # maximum cellular production rate
    max_ProRate = 2e-08
    # minimum cellular production rate
    min_ProRate = 0.
    # initial cellular production rate
    init_pro_Rate = 0.5e-08 #mean([min_ProRate,max_ProRate])
    # SD for mutation of cellular production rate
    mu_SD_ProRate = 0.01e-08

    # maximum signal concentration
    max_sigTh = 2e1
    # minimum signal concentration
    min_sigTh = 0.0001e1
    # initial signal threshold
    init_sig_Th = 0.5e1
    # SD for mutation of signal concentration
    mu_SD_sigTh = 0.01e1

    # define cheats parameters
    cheats_proRate = 0.
    cheats_sigTh = 2e1

    ################################################################################
    # initialization
    ################################################################################
    # initialize production rate
    pro_Rate = init_pro_Rate*ones(Float64,size_Pop)
    # initialize signal threshold
    sig_Th = init_sig_Th*ones(Float64,size_Pop)
    # initialize genotype fitness
    fit_Pop = zeros(Float64,size_Pop)
    # initialize cooperation payoff
    coopPayoff_Pop = zeros(Float64,size_Pop)
    # initialize cost for signaling
    sigCost_Pop = zeros(Float64,size_Pop)
    # initialize cost for cooperation
    coopCost_Pop = zeros(Float64,size_Pop)
    # efficiency for cooperation
    coopEff_Pop = zeros(Float64,size_Pop)
    # cheats index
    index_Cheats = zeros(Int64,size_Pop)
    # selection index
    index_Select = zeros(Int64,size_Pop)

    # initialize result saving space
    # all lineages
    fit_Evo = zeros(Float64,max_G)
    pro_Rate_Evo = zeros(Float64,max_G)
    sig_Th_Evo = zeros(Float64,max_G)
    coopPayoff_Evo = zeros(Float64,max_G)
    sigCost_Evo = zeros(Float64,max_G)
    coopCost_Evo = zeros(Float64,max_G)
    coopEff_Evo = zeros(Float64,max_G)
    # non-cheats
    fit_nonCheats_Evo = zeros(Float64,max_G)
    pro_Rate_nonCheats_Evo = zeros(Float64,max_G)
    sig_Th_nonCheats_Evo = zeros(Float64,max_G)
    coopPayoff_nonCheats_Evo = zeros(Float64,max_G)
    sigCost_nonCheats_Evo = zeros(Float64,max_G)
    coopCost_nonCheats_Evo = zeros(Float64,max_G)
    coopEff_nonCheats_Evo = zeros(Float64,max_G)
    # record cheats
    numCheats_Evo = zeros(Int64,max_G)

    ################################################################################
    # Evolution
    ################################################################################
    # initial evaluation
    eval_genotype(fit_Pop,coopPayoff_Pop,coopCost_Pop,coopEff_Pop,sigCost_Pop,pro_Rate,sig_Th,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lambda,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen)
    for g = 1:max_G

        # roulette-wheel selection
        for n = 1:size_Pop
          index_Select[n] = fortune_wheel(fit_Pop)
        end

        # update
        pro_Rate = pro_Rate[index_Select]
        sig_Th = sig_Th[index_Select]
        fit_Pop = fit_Pop[index_Select]
        coopPayoff_Pop = coopPayoff_Pop[index_Select]
        coopCost_Pop = coopCost_Pop[index_Select]
        coopEff_Pop = coopEff_Pop[index_Select]
        sigCost_Pop = sigCost_Pop[index_Select]
        index_Cheats = index_Cheats[index_Select]

        # generate cheats
        if numCheats_Evo[g]<size_Pop
            num_Cheats = randNum_Poisson(mu_Cheats)
            if num_Cheats!=0
                temp_Index = randperm(size_Pop)[1:num_Cheats]
                index_Cheats[temp_Index] = 1
                pro_Rate[temp_Index] = cheats_proRate
                sig_Th[temp_Index] = cheats_sigTh
            end
        end

        # mutate production rate
        mut_parameter(pro_Rate,mu_Production,mu_SD_ProRate,min_ProRate,max_ProRate,index_Cheats,size_Pop)

        # mutate signal threshold
        mut_parameter(sig_Th,mu_Th_Signal,mu_SD_sigTh,min_sigTh,max_sigTh,index_Cheats,size_Pop)

        # record cheats
        numCheats_Evo[g] = sum(index_Cheats)

        # save results
        fit_Evo[g] = mean(fit_Pop)
        pro_Rate_Evo[g] = mean(pro_Rate)
        sig_Th_Evo[g] = mean(sig_Th)
        coopPayoff_Evo[g] = mean(coopPayoff_Pop)
        sigCost_Evo[g] = mean(sigCost_Pop)
        coopCost_Evo[g] = mean(coopCost_Pop)
        coopEff_Evo[g] = mean(coopEff_Pop)

        if numCheats_Evo[g]<size_Pop
            index_nonCheats = (index_Cheats.==0)
            fit_nonCheats_Evo[g] = mean(fit_Pop[index_nonCheats])
            pro_Rate_nonCheats_Evo[g] = mean(pro_Rate[index_nonCheats])
            sig_Th_nonCheats_Evo[g] = mean(sig_Th[index_nonCheats])
            coopPayoff_nonCheats_Evo[g] = mean(coopPayoff_Pop[index_nonCheats])
            sigCost_nonCheats_Evo[g] = mean(sigCost_Pop[index_nonCheats])
            coopCost_nonCheats_Evo[g] = mean(coopCost_Pop[index_nonCheats])
            coopEff_nonCheats_Evo[g] = mean(coopEff_Pop[index_nonCheats])
        else
            fit_nonCheats_Evo[g] = fit_Evo[g]
            pro_Rate_nonCheats_Evo[g] = pro_Rate_Evo[g]
            sig_Th_nonCheats_Evo[g] = sig_Th_Evo[g]
            coopPayoff_nonCheats_Evo[g] = coopPayoff_Evo[g]
            sigCost_nonCheats_Evo[g] = sigCost_Evo[g]
            coopCost_nonCheats_Evo[g] = coopCost_Evo[g]
            coopEff_nonCheats_Evo[g] = coopEff_Evo[g]
        end

        # genotype evaluation
        eval_genotype(fit_Pop,coopPayoff_Pop,coopCost_Pop,coopEff_Pop,sigCost_Pop,pro_Rate,sig_Th,baseline,coop_Benefit,coop_Cost,sig_Cost,size_Pop,lambda,env_CellDen,grid_Size,base_Volume,decay_Rate,median_CellDen)

        fit_Pop[find(fit_Pop.<0)] = 0.
        
    end

end
