# TODO make more specific by removing just ppm or Ma from a string
# strip_units(s) = ...
"""
Strip units like "_ppm" or "_Ma" from the long name titles
"""
function short_name(s::String)
    if s == "Dy_Yb"
        return "Dy-Yb" # leave the Dy Yb element combination together
    else
        return split(s,"_")[1]
    end
end

function heatmap(M::NamedMatrix; kw...) #may clash when looking at a subarray of a DensityTensor
    heatmap(names(M,1), names(M,2), M; yflip=true, kw...)
end
