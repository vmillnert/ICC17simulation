# testing to see how to generate packets with headers

module Tmp

# we need a ringbuffer for storing the expected service and expected
# delay
type Buffer
    data::Vector{Float64}
    size::Int64 
    addIndex::Int64
    getIndex::Int64
    Buffer(data::Array{Float64}) = new(data, size(data)[1], 1, 1)
end

# Some functionality to add and get data from the ringbuffer
function add(buf::Buffer,p::Float64)
    buf.data[buf.addIndex] = p
    if buf.addIndex + 1 > buf.size
        buf.addIndex = 1
    else
        buf.addIndex = buf.addIndex + 1
    end
    if buf.addIndex == buf.getIndex
        if buf.getIndex + 1 > buf.size
            buf.getIndex = 1
        else
            buf.getIndex = buf.getIndex + 1
        end
    end
end

function get(buf::Buffer)
    buf.addIndex == buf.getIndex && return null
    p::Float64 = buf.data[buf.getIndex]
    if buf.getIndex + 1 > buf.size
        buf.getIndex = 1
    else
        buf.getIndex = buf.getIndex + 1
    end
    return p
end

# TODO: add initialization for this type
# Special type for a nfv function
type VNF

    id::Int32 # id of the function, i.e. what place it is in the
              # service-chain
    
    N::Int32 # size of the vectors
    dt::Float64 # time-step used for the simulation
    
    # deadline
    D::Float64 # local deadline for this function

    # Input related
    r::Vector{Float64} # input rate
    R::Vector{Float64} # total amount of arrived packets
    rhat::Vector{Float64} # estimated input in \Delta_i time-units

    # admission related
    alpha::Vector{Int64} # addimtance flag (whether we shuold admit
                           # packets into the function or not
    p::Vector{Float64} # rate of admitted packets
    P::Vector{Float64} # total number of admitted packets into the function

    # Buffer and queue related
    q::Vector{Float64} # queue size
    d::Vector{Float64} # delay for the packet exiting the function 
    dub::Vector{Float64} # worst-case delay for the packet entering
                         # the function

    # Service related
    ns::Float64 # nominal service rate 
    s::Vector{Float64} # service rate
    S::Vector{Float64} # number of processed packets
    smax::Vector{Float64} # maximum service-rate
    Smax::Vector{Float64} # cumulative maximum service-rate
    Slb::Vector{Float64} # lower bound on the cumulative service-rate

    # machine related
    Delta::Int64 # machine over-head (measured in indicies of dt)
    m::Vector{Int64} # number of running machines
    mref::Vector{Int64} # reference signal to the service-controller
    xihat::Vector{Float64} # average  machine uncertainty
    xilb::Float64 # lower bound on the machine uncertainty
    xiub::Float64 # upper bound on the machine uncertainty
    xi_j::Vector{Float64} # machine uncertainty for each machine in
                          # the function, we assume that it's constant
                          # but unknown
    
    # Evaluation related u::Vector{Float64} # utility function
    u::Vector{Float64} # utility function
    a::Vector{Float64} # availability function
    e::Vector{Float64} # efficiency function
    Dropped::Vector{Float64} # number of dropped packets
    Missed::Vector{Float64} # number of packets that missed a deadline
    
    VNF(id::Int64, N::Int64, dt::Float64, D::Float64, ns::Float64, Delta::Int64,
    xilb::Float64, xiub::Float64, m_max::Int64) = new(id, N, dt, D, zeros(N), zeros(N),
    zeros(N), zeros(Int64, N), zeros(N), zeros(N), zeros(N), zeros(N),
    zeros(N), ns, zeros(N), zeros(N), zeros(N), zeros(N), zeros(N),
    Delta, zeros(Int64, N), zeros(Int64, N), zeros(N), xilb, xiub, zeros(Float64, m_max),
    zeros(N), zeros(N), zeros(N), zeros(N), zeros(N))

end


# computing the delay for packets at the head of the queue at
# time-instance i
# using binary search
function delay!(vnf::VNF, i)
    dt = vnf.dt
    N = vnf.N
    if i>1
        j = searchsortedfirst(vnf.P, vnf.S[i])
        if j <= N
            vnf.d[i] = max(dt*(i-j), 0)
        else
            vnf.d[i] = 0.0
        end
    else
        vnf.d[i] = 0.0
    end
end

# Compute the upper bound of the expected delay for packets entering
# the function at time-index i
# Uses binary search
function delay_ub!(vnf::VNF, i)
    dt = vnf.dt
    N = length(vnf.s)
    Pi = vnf.P[i]
    Si = vnf.S[i]
    Slbi = vnf.Slb[i]
    # rethink how to compute expected delay for the first time!
    if i < N-vnf.Delta
        j = searchsortedfirst(vnf.Slb, Pi-Si+Slbi)
        if j < N
            vnf.dub[i] = max(dt*(j-i), 0)
        else
            #vnf.dub[i] = dt*vnf.Delta
            vnf.dub[i] = 0.0
        end
    else
        #vnf.dub[i] = dt*vnf.Delta
        vnf.dub[i] = 0.0
    end
end


end
