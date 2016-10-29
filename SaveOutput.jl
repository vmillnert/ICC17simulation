#=
Simple script to save down the output of the simulation to a text-file
Victor Millnert at the Department of Automatic Control
Theory, Lund University
=#

using JSON


function save_output(U_data, U_sota1_data, U_sota2_data, U_sota1_ad_data, U_sota2_ad_data,
                     filename::ASCIIString)

    # convert the string to JSON format

    s = "{\"U_data\" : $(U_data) , \"U_sota1_data\" : $(U_sota1_data), \"U_sota2_data\" : $(U_sota2_data), \"U_sota1_ad_data\" : $(U_sota1_ad_data), \"U_sota2_ad_data\" : $(U_sota2_ad_data) }";

    f = open("$(filename)", "w")
    write(f, s)
    close(f)

end #end function
