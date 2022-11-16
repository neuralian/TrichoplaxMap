using GLMakie
using Colors

# global parameters
worldWide = 1600.0 # microns
pxl_per_um = 0.5      # pixels per micron display
tr_diameter = 1000.0  # um

# derived parameters
tr_location = Point2f(worldWide/2.0, worldWide/2.0)
pixelWide = Int64(round(worldWide*pxl_per_um))

# Trichoplax Utilities
struct Trichoplax
    map::Array{Float64,2}
    radius::Vector{Float64} # nb declaring as array makes value mutable as diameter[]
    location::Vector{Point2f}  
    edge_handle::Poly
    nLabels::Vector{Int64}
    label_handle::Vector{Scatter}  # handles of scatterplots
    label_name::Vector{String}
end


function Trichoplax(radius::Float64, location::Point2f)
    plt_handle =  poly!(Circle(location, radius), 
                  color = RGBA(.8, .8, .8, .1), strokecolor = RGB(.6, .6, .6), strokewidth = .5)
    Trichoplax( zeros(1,1), 
                [radius], [location], plt_handle, 
                [0], [scatter!(Point2f(NaN, NaN))], [""])
end

distance(p::Point2f) = sqrt(p[1]^2 + p[2]^2)
distance(p::Point2f, q::Point2f) = distance(p-q)
distance(p::Point2f, Q::Vector{Point2f}) = [distance(p,q) for q in Q]
# TRUE if any points in Q are within Δ of p
anycloserthan(Δ::Float64, p::Point2f, Q::Vector{Point2f}) = length(findall(distance(p,Q).<Δ))>0


# function draw(trichoplax::Trichoplax)
#     poly!(Circle(trichoplax.location[], trichoplax.radius[]), 
#     color = :white, strokecolor = :black, strokewidth = 1)
# end

function move(trichoplax::Trichoplax, dp::Point2f)

    trichoplax.location[] +=  dp
    trichoplax.edge_handle[1] = Circle(trichoplax.location[], trichoplax.radius[])  # move outline
    for i in 1:trichoplax.nLabels[]   # move labels
        trichoplax.label_handle[i][1] = trichoplax.label_handle[i][1][] .+= dp  
   end
end

function moveto(trichoplax::Trichoplax, dx::Point2f)
    move(trichoplax, p-trichoplax.location[])
end


# sample of n points (type Point2f) from radial intensity function
# intensity(d,r) is expected number of points in a cell of diameter d um at radius r um
# with minimum distance Δ between samples (default=1 cell diameter, ie max 1 sample per cell) 
function label(tr::Trichoplax, labelName::String, n::Int64, d::Float64, intensity::Function, Δ::Float64=-1.0)

    dA = π*(d/2.0)^2  # area element
    if Δ<0.0 Δ = d end   # default min spacing

    count = 0
    P = Array{Point2f,1}(undef, n)  # array of n points, initially undefined
    while count < n
        p = Point2f((2.0*rand(2).-1.0).*tr.radius)  # uniform random point in area containing the trichoplax
        if rand()*dA<intensity(distance(p)) # provisional accept based on density
            if count<1 || !anycloserthan(Δ, p, P[1:count]) 
                count = count + 1    
                P[count] = p           # accept
            end
        end
    end
    trichoplax.nLabels[] = trichoplax.nLabels[] + 1
    trichoplax.label_handle[trichoplax.nLabels[]] = scatter!(ax,tr.location.+P)
    trichoplax.label_name[trichoplax.nLabels[]] = labelName
    return P
end


function uniformBand(r::Float32)
    if (r>100.) & (r<200.) return 1.0
    else return 0.0
    end
end

F = Figure(resolution = (worldWide*pxl_per_um,worldWide*pxl_per_um))
ax = Axis(F[1,1], aspect = 1)
xlims!(0, worldWide)
ylims!(0, worldWide)

# construct world

worldMap = zeros(pixelWide, pixelWide)


# construct and draw Trichoplax
trichoplax  = Trichoplax(tr_diameter/2.0, tr_location)

P = label(trichoplax, "LABEL", 100, 5.0, uniformBand)

#tr_plthandle = draw(trichoplax)

display(F)