#### ------------------------------------ ####
####                TODO                  ####
#### ------------------------------------ ####


# using PyPlot
# using Iterators

@everywhere include("./Tmp.jl")
@everywhere include("./GetData.jl")
@everywhere include("./NFVSimulation.jl")

@everywhere function average_utility(nfv::Vector{Tmp.VNF})
    N = nfv[1].N
    dt = nfv[1].dt
    n = length(nfv)

    x = zeros(Float64, N)
    
    for i = 1:N
        x[i] =  1/n*sum(map(x-> x.u[i], nfv))
    end
    return x
end


#### ------------------------------------ ####
#### Retrieve the data from the .rrd-file ####
#### ------------------------------------ ####

# when simulating all 5 methods, it seems like 7 cores is on the limit

K = 3000 # number of simulations
# U_data = Dict{Int64,Vector{Float64}}()
# U_sota1_data = Dict{Int64,Vector{Float64}}()
# U_sota2_data = Dict{Int64,Vector{Float64}}()
U_data = SharedArray(Float64, K)
U_sota1_data = SharedArray(Float64, K)
U_sota2_data = SharedArray(Float64, K)
U_sota1_ad_data = SharedArray(Float64, K)
U_sota2_ad_data = SharedArray(Float64, K)

# Start = zeros(K)
# t_data = Dict{Int64,Vector{Float64}}()

@everywhere function main()
   
    # Start time used in paper
    # start = 1472468100 + day + 11*hour
    # day = 60*60*24
    # hour = 60*60
    # minute = 60
    # First sample data: 1472558400
    # Last sample date: 1472990400 

    # Generate a random starting time from the sample data
    # START = rand(1472558400:1:1472990400)
    # END = START+hour  # let's run the simulation for 2 hours
    # file = "port-555.rrd"
    dt = 1e-3
    plot_data = true
    peak = 1e7 #let's have the peak at 10 million packets per second
    println("getting input data")
    # @time input, t = GetData.getPacketData(START, END, dt, file, peak);

    # choose the duration for the simulation.
    # available: 1, 2, 6
    
    # for running on our cloud-0x servers it seems like 1hour
    # simulations using 10 cores limits the RAM-memory
    
    @time input, t = GetData.getPacketData(dt,peak);
    println("input data retrieved")
    N = Int64(length(input))

    #### ------------------------------------ ####
    #### Set up the virtual network functions ####
    #### ------------------------------------ ####
    # let's create an array of n network functions
    n = 5
    nfv = Vector{Tmp.VNF}(n)

    # intialize the service chain
    @time initialize!(nfv, input, dt)

    # println("copying the service-chain")
    nfv_sota1 = deepcopy(nfv)
    nfv_sota2 = deepcopy(nfv)
    nfv_sota1_ad = deepcopy(nfv)
    nfv_sota2_ad = deepcopy(nfv)
    
    t_n = Int64(round(N/3000))
    t = collect(t)

    nfv_small = Vector{Tmp.VNF}(n)
    nfv_sota1_small = Vector{Tmp.VNF}(n)
    nfv_sota2_small = Vector{Tmp.VNF}(n)
    nfv_sota1_ad_small = Vector{Tmp.VNF}(n)
    nfv_sota2_ad_small = Vector{Tmp.VNF}(n)
    
    println("simulate state of the art with 'dynamic scaling' ")
#    @time simulate_sota!(nfv_sota1, false)
    @time simulate_sota!(nfv_sota1, false, false)
    @time simulate_sota!(nfv_sota1_ad, false, true)
    println("finished state of the art 1")

    
    println("simulate state of the art with 'dynamic overprovisioning' ")
    # @time simulate_sota!(nfv_sota2, true)
    @time simulate_sota!(nfv_sota2, true, false)
    @time simulate_sota!(nfv_sota2_ad, true, true)
    println("finished state of the art 2")

    
    println("simulate our method")
    @time simulate!(nfv)
    println("finished our method")
    

    println("downsampling the information stored")
    @time t_small = downsample!(nfv, nfv_small, t, t_n)
    @time t_small = downsample!(nfv_sota1, nfv_sota1_small, t, t_n)
    @time t_small = downsample!(nfv_sota2, nfv_sota2_small, t, t_n)
    @time t_small = downsample!(nfv_sota1_ad, nfv_sota1_ad_small, t, t_n)
    @time t_small = downsample!(nfv_sota2_ad, nfv_sota2_ad_small, t, t_n)
    println("finished downsampling")

    # plot our method
    # close("all")
    # @time plot_all(nfv_small, t_small, t_n, false)
    # plot state of the art
    # @time plot_all(nfv_sota_small, t_small, t_n, true)

    
    @time U = average_utility(nfv_small)
    @time U_sota1 = average_utility(nfv_sota1_small)
    @time U_sota2 = average_utility(nfv_sota2_small)
    @time U_sota1_ad = average_utility(nfv_sota1_ad_small)
    @time U_sota2_ad = average_utility(nfv_sota2_ad_small)

    # t_small = t_small - t_small[1]
    # figure("Average Utility", figsize=(10,5))
    # PyPlot.plot(t_small, U, label="U(t)")
    # PyPlot.plot(t_small, U_sota1, label="U_sota1(t)")
    # PyPlot.plot(t_small, U_sota2, label="U_sota2(t)")
    # legend(loc="lower right")

    # figure("Input", figsize=(10,5))
    # PyPlot.plot(t_small, nfv_small[1].r, label="r(t)")
    # legend(loc="lower right")
    return mean(U), mean(U_sota1), mean(U_sota2),  mean(U_sota1_ad), mean(U_sota2_ad)

end
    # U_data[k] = mean(U)
    # U_sota1_data[k] = mean(U_sota1)
    # U_sota2_data[k] = mean(U_sota2)
    # # Start[k] = START
    # # t_data[k] = t_small


@parallel for k in 1:K
    println("starting iteration $k")
    U_data[k], U_sota1_data[k], U_sota2_data[k],  U_sota1_ad_data[k], U_sota2_ad_data[k] = main()
    println("finished iteration $k")
end

# println("simulation is done!!!!")

# println("Our method")
# @show U_data
# @show U_sota1_data
# @show U_sota2_data


