#=
Simple script to extract the data from a .rrd-file and plot it
Victor Millnert at the Department of Automatic Control
Theory, Lund University
=#
module GetData

using JSON
# using PyPlot
using Interpolations

#=
start and end time in seconds since UNIX epoc
2016-08-29 1472468100
2016-09-05 1473072600
24 hours 86400
=#
# day = 86400
# START=1472468100
# #END  =1473072600
# END = START+day                 

# possible readings from the in-traffic in the .rrd-file
# INOCTETS
# INUCASTPKTS
# INNUCASTPKTS
# INBROADCASTPKTS
# INMULTICASTPKTS

# function getPacketData(t_start::Int64, t_end::Int64, dt_fine::Float64, file_name::String, peak::Float64, plot_data::Bool)
function getPacketData(dt_fine::Float64, peak::Float64)
 
    # READ and PARSE the .rrd-file named "file_name"
    # file = JSON.parse(readstring(`rrdtool xport --json
    #                              --start $t_start --end $t_end
    #                              DEF:inoc=$file_name:INUCASTPKTS:AVERAGE
    #                              XPORT:inoc`));

    # file = JSON.parse(readall(`rrdtool xport --json
    #                              --start $t_start --end $t_end
    #                              DEF:inoc=$file_name:INUCASTPKTS:AVERAGE
    #                              XPORT:inoc`));


    # randomly select one of the 100 input-files
    i = rand(1:1000)
    file = JSON.parsefile("./data_1h/input_$i.txt")
    @show i
    # extract the "data"-fields, containing the number of packets
    data = file["data"];
    dt = file["meta"]["step"];
    t_start = file["meta"]["start"];
    t_end = file["meta"]["end"];
    # println("step size is $(dt) ")

    N = length(data)

    # convert the data-array from vector of {Any,1} to a vector of Float64
    input = Vector{Float64}(N)
    for i in 1:N
        if typeof(data[i][1]) == Float64
            input[i] = data[i][1];
        else
            input[i] = 0.0
        end
    end

    # convert the starting-time and stop-time to a date vector
    t = t_start:dt:t_end;
    if length(t) != N
        t = t_start:dt:t_end-dt;
    end

    itp = interpolate(input, BSpline(Quadratic(Reflect())), OnGrid());
    i_fine = 1:dt_fine/dt:N;
    N_fine = length(i_fine)
    t_fine = t_start:dt_fine:t_end;
    if length(i_fine) != length(t_fine)
        i_fine = 1:dt_fine/dt:N;
        t_fine = t_start:dt_fine:t_end-dt;
    end
    if mod(N_fine/10,10) != 0
        i_fine = i_fine[1:end-1]
        t_fine = t_fine[1:end-1]
    end

    fine = zeros(length(i_fine))
    for i in 1:length(i_fine)
        fine[i] = itp[i_fine[i]]
    end
#    fine = itp[i_fine];

    # scale the input such that the peak is at the specified peak
    scale = peak/maximum(fine)
    fine = fine*scale
    input = input*scale
    
    # if plot_data
    #     # ---------------
    #     # plot the output
    #     # ---------------
    #     t_date = map(x-> Dates.unix2datetime(x), t);
    #     t_fine_date = map(x-> Dates.unix2datetime(x), t_fine);
        
    #     close("all")
    #     # figure("Incoming octets", figsize=(15,10))
    #     # plot(t_date, input)

    #     # figure("Interpolated data", figsize=(15,10))
    #     # plot(t_fine, fine, label="interpolated")
    #     # legend(loc="upper right")

    #     figure("Comparison data", figsize=(15,10))
    #     plot(t_date, input, label="sampled")
    #     plot(t_fine_date, fine, label="interpolated")
    #     legend(loc="upper right")
    # end
    return fine, t_fine
end #end function

# # just to make it possibel to skip specifying if we should plot or
# # not...
# function getPacketData(t_start::Int64, t_end::Int64, dt_fine::Float64, file_name::String, peak::Float64)
#      return getPacketData(t_start, t_end, dt_fine, file_name, peak, false)
# end

end #end module
