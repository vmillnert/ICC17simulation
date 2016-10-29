# using PyPlot 
# using PGFPlots
# using Iterators


include("./Tmp.jl")
include("./GetData.jl")

using Tmp
#### ---------------------------- ####
#### Initialize the service-chain ####
#### ---------------------------- ####
function initialize!(nfv::Vector{Tmp.VNF}, input::Vector{Float64},
                     dt::Float64)

    println("Setting up the problem")

    n = length(nfv)
    # before initializing them we have to define
    # ns - the nominal service rate
    # time_overhead - the time needed to start/stop a machine
    # Delta - the time_overhead expressed in indicies (i.e. number of dt's it correspond to)
    # D - the local deadline that packets has to pass through the function
    # xilb - lower bound on the machine uncertainties
    # xiub - upper bound on the machine uncertainties

    peak = maximum(input)
    N = Int64(length(input))
    nslb = 1e5 # lower bound for nominal service-rate
    nsub = 2e5 # upper bound for the nominal service-rate
    ns = rand(nslb:1e3:nsub, n) # generate the nominal service-rate randomly

    tlb = 30.0 # lower bound for time-overhead
    tub = 2*60.0  # upper bound for time-overhead
    time_overhead = rand(tlb:1e-2:tub,n) # generate the time-overhead randomly

    # Time-overhead
    Delta = map(x -> Int64(round(x/dt)), time_overhead)

    # Deadline
    D = 0.01*ones(Float64, n)

    xilb = zeros(Float64,n)
    xiub = zeros(Float64,n)
    for i in 1:n
        xilb[i] = rand(linspace(-ns[i]*0.3, 0, 1000), 1)[1]
        xiub[i] = rand(linspace(0, ns[i]*0.3, 1000), 1)[1]
    end

    for i in 1:n

        # estimate how many machines we will need at most (just to make
        # the simulation run faster)
        m_max = Int64(round(peak/ns[i]*10))

        nfv[i] = Tmp.VNF(Int64(i), N, dt, D[i], ns[i], Delta[i], xilb[i], xiub[i], m_max)

        # The number of machines running for the first time_overhead seconds, i.e. for Delta[i] indicies
        nfv[i].mref[1:Delta[i]] = Int64(round(input[1]/ns[i]))*ones(Int64, Delta[i])

        # The number of machines we are starting with at the first time-instance
        nfv[i].m[1:Delta[i]] = nfv[i].mref[1]

        ## THE FOLLOWING IS JUST TO MAKE IT POSSIBLE TO USE A BINARY SEARCH
        ## ALGORITHM FOR FINDING THE DELAY
        # Set up the lower bound for cumulative service-rate
        nfv[i].Slb = Inf16*ones(N)
        nfv[i].Slb[1] = 0.0
        for j = 1:Delta[i]-1
            nfv[i].Slb[j+1]  = nfv[i].Slb[j] + dt*nfv[i].m[j]*(nfv[i].ns+nfv[i].xilb)
        end
        nfv[i].P = Inf16*ones(N)
        nfv[i].P[1] = 0.0
        i==1 ? nfv[i].r = input : nfv[i].r = zeros(Float64,N)

        # ------------------------------------ #
        #   set up the machine uncertainties
        # ------------------------------------ #
        # We assume that each machine will have a certain uncertainty,
        # chosen uniformly at random in the interval [xilb, xiub].
        nfv[i].xi_j = rand(linspace(nfv[i].xilb,nfv[i].xiub, 100), m_max)

        # # We then assume that the machine uncertainty will vary a little
        # # bit around this value over time. Also according to a uniform
        # # distribution withing 1%
        # nfv[i].xi_j = ones(m_max,N) # we assume no machine uncertainty for now...
        # for j in 1:m_max
        #     #nfv[i].xi_j[j,:] = rand(linspace(tmp[j]*0.9, tmp[j]*1.1, 100), N)
        #     nfv[i].xi_j[j,:] = tmp[j]*ones(N)
        # end
        nfv[i].rhat[1] = input[1]
        
    end
    println("Set up of the problem is done")

end



### ----------------------------------------- ###
### Function for simulating the service-chain ###
### ----------------------------------------- ###
function simulate!(nfv::Vector{Tmp.VNF})
    M = 1 # filter constant when computing the derivative used in rhat
    x = 0.0
    N = nfv[1].N
    n = length(nfv)
    for i = 1:N

        # if mod(i, N/10)==0
        #     println("Completed $(Int64(round(i/N*100))) percent")
        # end

        for vnf in nfv
            # @show vnf.id, i
            # check for the expected deadline to see if we should admit
            # packets or not
            vnf.dub[i] <= vnf.D ? vnf.alpha[i]=1 : vnf.alpha[i]=0

            # admit (or not) the new packets
            vnf.p[i] = vnf.alpha[i]*vnf.r[i]

            # compute the numer of machines to use
            if i > vnf.Delta
                vnf.m[i] = vnf.mref[i - vnf.Delta]
            end
            # compute the average machine uncertainty on the running machines
            m = vnf.m[i]
            vnf.xihat[i] = 1/m*sum(vnf.xi_j[1:m])
            # compute the maximum service-rate at this time
            vnf.smax[i] = m*(vnf.ns + vnf.xihat[i])

            # compute the service rate at this time
            vnf.q[i]<=1e-5 && vnf.p[i] <= vnf.smax[i] ? vnf.s[i] = vnf.p[i] : vnf.s[i] = vnf.smax[i]

            # compute rhat, i.e. the expected input at time t+Delta_i
            # This one will depend on whether it is for the first function or not.

            ### TODO ###
            # Have to improve the estimation of incoming input !!!
            
            if i > 1

                if vnf.id == 1
                    vnf.rhat[i] = vnf.r[i] + 1/(1+M)*(vnf.rhat[i-1] -
                                                      vnf.r[i-1]) +
                    vnf.Delta*M/(1+M)*(vnf.r[i]
                                       - vnf.r[i-1])
                else
                    # !!! TODO !!!
                    # COMPUTE RHAT THE OTHER WAY
                    # !!! TODO !!!

                    vnf.rhat[i] =
                    min(nfv[vnf.id-1].mref[i]*(nfv[vnf.id-1].ns+nfv[vnf.id-1].xihat[i]),
                    nfv[vnf.id-1].rhat[i])

                end
                
            end
            
            @assert vnf.rhat[i] < Inf16 "Shit hit the fan! id: $(vnf.id), i:$i, vnf.rhat[i]:$(vnf.rhat[i])"
            
            # compute mref
            x = vnf.rhat[i]/(vnf.ns + vnf.xihat[i])
            if floor(x)*1/x >= x*1/ceil(x)
                vnf.mref[i] = Int64(floor(x))
            else
                vnf.mref[i] = Int64(ceil(x))
            end

            # compute the availability function
            
            vnf.a[i] = vnf.s[i]/vnf.r[i- Int64(round(vnf.d[i]/vnf.dt))]
            # compute the efficiency function
            vnf.e[i] = vnf.s[i]/vnf.smax[i]
            # compute the utility function
            vnf.u[i] = vnf.a[i]*vnf.e[i]

            @assert vnf.a[i] < Inf16 "Shit hit the fan!! i: $i, id: $(vnf.id) "
            
            ## IF WE HAVE MORE THAN ONE FUNCTION IN THE CHAIN:
            if vnf.id < length(nfv)
                nfv[vnf.id+1].r[i] = vnf.s[i]
            end
            
            if i < N
                # update the total amount of admitted packets
                vnf.P[i+1] = vnf.P[i] + vnf.dt*vnf.p[i]
                # update the total amount of served packets
                vnf.S[i+1] = min(vnf.S[i] + vnf.dt*vnf.s[i], vnf.P[i+1])
                # update the queue size
                vnf.q[i+1] = vnf.P[i+1] - vnf.S[i+1]

                # compute the lower bound of the total amount of served
                # packets
                if i < N - vnf.Delta
                    vnf.Slb[i+vnf.Delta] = vnf.Slb[i+vnf.Delta-1] + vnf.dt*vnf.mref[i]*(vnf.ns+vnf.xilb)
                end
                
                # compute the upper bound of the epected delay for packets
                # entering at the next time-instance
                Tmp.delay_ub!(vnf, i+1)
                Tmp.delay!(vnf,i+1)

            end    

        end

    end

end


### ----------------------------------------- ###
### Function for simulating State-of-the-Art  ###
### ----------------------------------------- ###
function simulate_sota!(nfv::Vector{Tmp.VNF}, over_provision::Bool, admission_control::Bool)
    M = 1 # filter constant when computing the derivative used in rhat
    x = 0.0
    N = nfv[1].N
    n = length(nfv)
    for i = 1:N

        # if mod(i, N/10)==0
        #     println("Completed $(Int64(round(i/N*100))) percent")
        # end

        for vnf in nfv

            if admission_control
                # check for the expected deadline to see if we should admit
                # packets or not
                vnf.dub[i] <= vnf.D ? vnf.alpha[i]=1 : vnf.alpha[i]=0
            else
                # No admission control => always admit!
                vnf.alpha[i] = 1
            end
            
            # admit the new packets
            vnf.p[i] = vnf.alpha[i]*vnf.r[i]

            # compute the numer of machines to use
            if i > vnf.Delta
                vnf.m[i] = vnf.mref[i - vnf.Delta]
                @assert vnf.m[i] > 0 "no machines running i:$i, id:$(vnf.id), m:$(vnf.m[i]), vnf.mref[i-Delta]$(vnf.mref[i-vnf.Delta])"

            end
            @assert vnf.m[i] > 0 "no machines running i:$i, id:$(vnf.id), m:$(vnf.m[i])"


            # compute the average machine uncertainty on the running
            # machines
            m = vnf.m[i]
            vnf.xihat[i] = 1/m*sum(vnf.xi_j[1:m])
            # compute the maximum service-rate at this time
            vnf.smax[i] = m*(vnf.ns + vnf.xihat[i])
            
            # compute the service rate at this time
            vnf.q[i]<=1e-5 && vnf.p[i] <= vnf.smax[i] ? vnf.s[i] = vnf.p[i] : vnf.s[i] = vnf.smax[i]


            # ------------------ #
            # Control of service #
            # ------------------ #

            if over_provision
                # over provision by 30%
                vnf.mref[i] = Int64(round(vnf.r[i]/vnf.ns*1.1))
            else
                # If the utilization ("efficiency") is above 70% of max capacity we
                # increase the number of machines
                # If the utilization ("efficiency") is below 30% of max
                # capacity we decrease the number of machines
                if i>1
                    if vnf.e[i-1] > 0.99
                        vnf.mref[i] = vnf.m[i]+1
                    elseif vnf.e[i-1] < 0.95
                        vnf.mref[i] = max(1,vnf.m[i]-1)
                    else
                        vnf.mref[i] = vnf.m[i]
                    end
                end
            end
            
            # Compute the availability function. If the packet miss the
            # deadline the availability is considered to be zero

            if vnf.d[i] > vnf.D
                vnf.a[i] = 0.0
            else
                k = Int64(max(1, i-Int64(round(vnf.d[i]/vnf.dt))))
                if vnf.r[k] == 0
                    vnf.a[i] = 0.0
                else
                    vnf.a[i] = vnf.s[i]/vnf.r[k]
                end
                
                @assert vnf.a[i] < Inf16 "Shit hit the fan!!
                vnf.r[k]:$(vnf.r[k]), vnf.a[i]:$(vnf.a[i]), k:$k, i:
                $i, id: $(vnf.id) vnf.s[i]:$(vnf.s[i])"

            end

            
            # compute the efficiency function
            vnf.e[i] = vnf.s[i]/vnf.smax[i]
            # compute the utility function
            vnf.u[i] = vnf.a[i]*vnf.e[i]
            
            
            
            ## IF WE HAVE MORE THAN ONE FUNCTION IN THE CHAIN:
            if vnf.id < length(nfv)
                nfv[vnf.id+1].r[i] = vnf.s[i]
            end
            
            if i < N
                # update the total amount of admitted packets
                vnf.P[i+1] = vnf.P[i] + vnf.dt*vnf.p[i]
                # update the total amount of served packets
                vnf.S[i+1] = min(vnf.S[i] + vnf.dt*vnf.s[i], vnf.P[i+1])
                # update the queue size
                vnf.q[i+1] = vnf.P[i+1] - vnf.S[i+1]

                # compute the lower bound of the total amount of served
                # packets
                if i < N - vnf.Delta
                    vnf.Slb[i+vnf.Delta] = vnf.Slb[i+vnf.Delta-1] + vnf.dt*vnf.mref[i]*(vnf.ns+vnf.xilb)
                end
                
                # compute the upper bound of the epected delay for packets
                # entering at the next time-instance
                Tmp.delay_ub!(vnf, i+1)
                Tmp.delay!(vnf,i+1)

            end    

        end

    end

end



### -------------------------------------------------------------- ###
### Function for downsampling the data, to make it easier to plot  ###
### -------------------------------------------------------------- ###

function downsample!(nfv::Vector{Tmp.VNF}, nfv_small::Vector{Tmp.VNF},
                     t::Vector{Float64}, n::Int64)


    
    # Remove the beginning and the end of the simulation, i.e. where
    # we don't have information
    j = 1
    for vnf in nfv
        j = max(vnf.Delta,j)
        # for k in 1:N
        #     vnf.d[k] = max(vnf.d[k], 0)
        # end        
    end
    j = j+1
        
    # Downsample the vector to make it easier to plot it
    N = nfv[1].N
    newsize = Int64(round((N-2*j)/n))
    d = zeros(Float64, newsize)
    m = zeros(Int64, newsize)
    mref = zeros(Int64, newsize)
    rhat = zeros(Float64, newsize)
    q = zeros(Float64, newsize)
    s = zeros(Float64, newsize)
    u = zeros(Float64, newsize)
    a = zeros(Float64, newsize)
    e = zeros(Float64, newsize)
    r = zeros(Float64, newsize)
    for vnf in nfv
        
        nfv_small[vnf.id] = Tmp.VNF(Int64(vnf.id), newsize, vnf.dt,
                                    vnf.D, vnf.ns, vnf.Delta, vnf.xilb, vnf.xiub, Int64(length(vnf.xi_j)))

        for k in 1:newsize
            i = Int64(k*n)
            nfv_small[vnf.id].d[k] = vnf.d[i]
            nfv_small[vnf.id].m[k] = vnf.m[i]
            nfv_small[vnf.id].mref[k] = vnf.mref[i]
            nfv_small[vnf.id].rhat[k] = vnf.rhat[i]
            nfv_small[vnf.id].q[k] = vnf.q[i]
            nfv_small[vnf.id].s[k] = vnf.s[i]
            nfv_small[vnf.id].u[k] = vnf.u[i]
            nfv_small[vnf.id].r[k] = vnf.r[i]
            nfv_small[vnf.id].a[k] = vnf.a[i]
            nfv_small[vnf.id].e[k] = vnf.e[i]
        end

    end
    t_small = zeros(Float64, newsize)
    for k = 1:newsize
        i = Int64(k*n)
        t_small[k] = t[i]
    end

    return t_small
    # t = t[j:end-j]
    # t = collect(takenth(t, n))
end



### ------------------------------------------------------------ ###
### Function for plotting the different things in the simulation ###
### ------------------------------------------------------------ ###
# function plot_all(nfv, t::Vector{Float64}, n, sota::Bool)

#     if sota
#         s = "State of the Art"
#     else
#         s = ""
#     end
#     for vnf in nfv
#         id = vnf.id
#         figure("delay $s", figsize=(10,5));
#         PyPlot.plot(t, vnf.d, label=string("d_$id"))
#              legend(loc="upper right")

#         figure("machines $s", figsize=(10,5));
#         PyPlot.plot(t, vnf.m, label=string("m_$id"))
#         legend(loc="upper right")

#         figure("utility function $s", figsize=(10,5));
#         PyPlot.plot(t, vnf.u, label=string("u_$id"))
#         legend(loc="upper right")

#         D = Int64(round(vnf.Delta/n))
#         figure(string("expected and true input $id $s"), figsize=(10,5));
#         PyPlot.plot(t, vnf.r, label=string("r_$id"))
#         PyPlot.plot(t[D:end], vnf.rhat[1:end-D+1], label=string("rhat_$id"))
#         legend(loc="upper right")

#         figure(string("service and inupt rate $id $s"), figsize=(10,5));
#         PyPlot.plot(t,vnf.r, label=string("r_$id"));
#         PyPlot.plot(t,vnf.s, label=string("s_$id"));
#         legend(loc="upper right");
 
#     end
# end

# function save_pgfplot(nfv, t::Vector{Float64}, U::Vector{Float64}, n::Int64, sota::Bool)


#     if sota
#         s="_sota"
#     else
#         s=""
#     end
#     t = t-t[1] # normalize the time
#     t = t/(60*60) # convert it to hours

    
#     # ----------------------- #
#     #          Input          #
#     # ----------------------- #
#     p = Axis(style="width=10cm, height=5cm, scaled ticks=false,
#              yticklabel style={/pgf/number format/fixed,
#                      /pgf/number format/precision=1}", 
#              PGFPlots.Plots.Linear(t, nfv[1].r, mark="none"),
#              xlabel="time (hours)", ylabel="input (pps)", legendPos="north east",
#              xmin=0, xmax=t[end])
#     save("input$s.pdf", p)
    
#     # # ------------------ #
#     # # number of machines #
#     # # ------------------ #
#     p = Axis(style="width=10cm, height=5cm, legend columns=-1, scaled ticks=false",
#              map(x->
#              PGFPlots.Plots.Linear(t, x.m, mark="none",
#              legendentry=string("\$m_$(x.id)(t)\$")), nfv),
#              xlabel="time (hours)", ylabel="machines", legendPos="south east",
#              xmin=0, xmax=t[end], ymin=0)
#     save("machines$s.pdf", p)


#     # # ------------------ #
#     # #       delay        #
#     # # ------------------ #
#     p = Axis(style="width=10cm, height=5cm, scaled ticks=false,
#              yticklabel style={/pgf/number format/fixed,
#                      /pgf/number format/precision=3},
#              try min ticks=5,
#              legend columns=-1",
#              map(x->
#              PGFPlots.Plots.Linear(t, x.d, mark="none",
#              legendentry=string("\$d_$(x.id)(t)\$")), nfv),
#              xlabel="time (hours)", ylabel="delay (s)", legendPos="north east",
#              xmin=0, xmax=t[end],ymin=0.0 )
#     save("delay$s.pdf", p)


#     # ----------------------- #
#     #     average utility     #
#     # ----------------------- #
#     p = Axis(style="width=10cm, height=5cm, scaled ticks=false", 
#              PGFPlots.Plots.Linear(t, U, mark="none",
#                                    legendentry=L"$U(t)$"),
#              xlabel="time (hours)", ylabel="utility", legendPos="south east",
#              ymin=0, ymax=1.01, xmin=0, xmax=t[end])
#     save("average_utility$s.pdf", p)
    

#     for vnf in nfv

#         # ----------------- #
#         # service and input #
#         # ----------------- #
#         # p = Axis(style="width=10cm, height=5cm, scaled ticks=false,
#         #      yticklabel style={/pgf/number format/fixed,
#         #              /pgf/number format/precision=1}", 
#         #          [PGFPlots.Plots.Linear(t, vnf.s, mark="none",
#         #                                 legendentry=string("\$s_$(vnf.id)(t)\$")),
#         #           PGFPlots.Plots.Linear(t, vnf.r, mark="none",
#         #                                 legendentry=string("\$r_$(vnf.id)(t)\$"))],
#         # xlabel="time (hours)", ylabel="packets per second",
#         # xmin=0, xmax=t[end])
#         # save(string("input_service_$(vnf.id)$s.pdf"), p)

#         # ----------------------- #
#         # true and expected input #
#         # ----------------------- #
#         # p = Axis(style="width=10cm, height=5cm, scaled ticks=false,
#         #      yticklabel style={/pgf/number format/fixed,
#         #              /pgf/number format/precision=1}", 
#         #          [PGFPlots.Plots.Linear(t, vnf.r, mark="none",
#         #                                 legendentry=string("\$s_$(vnf.id)(t)\$")),
#         #           PGFPlots.Plots.Linear(t, vnf.rhat, mark="none",
#         #                                 legendentry=string("\$ \\hat{r}_$(vnf.id)(t)\$"))],
#         #          xlabel="time (hours)", ylabel="packets per second")
#         # save(string("true_expected_$(vnf.id)$s.pdf"), p)
        

#         # # ------------------ #
#         # # queue_size  #
#         # # ------------------ #
#         # p = Axis(style="width=10cm, height=5cm, scaled ticks=false", 
#         #          PGFPlots.Plots.Linear(t, vnf.q, mark="none",
#         #                                 legendentry=string("\$q_$(vnf.id)(t)\$")),
#         #          xlabel="seconds", ylabel="packets")
#         # save(string("queue_$(vnf.id).pdf"), p)
        
#     end
# end
    
    
   
