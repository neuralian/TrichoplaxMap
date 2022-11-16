using GLMakie
using Colors
using Debugger

# global parameters
worldWide = 1600.0 # microns
pxl_per_um = 0.5      # pixels per micron display
tr_diameter = 1000.0  # um
MAXLABELS = 100

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
                [0], Vector{Scatter}(undef,MAXLABELS), Vector{String}(undef,MAXLABELS) )
             #   [0], [scatter!(Point2f(NaN, NaN))], [""])
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

# sample of n points p::point2f from density(r) defined on (0,rmax), where r is distance of p from origin
# minimum distance Δ between sample points
function radialsample(n::Int64, density::Function,  rmax::Float64, Δ::Float64)
    crowdCount = 0
    overCrowd = 25   # will exit if fail to find room to place a label consecutively this many times

    popCount = 0
    P = Array{Point2f,1}(undef, n)  # array of n points, initially undefined
    while (popCount < n) && (crowdCount < overCrowd)
        p = Point2f((2.0*rand(2).-1.0).*rmax)  # uniform random point in area containing the trichoplax
        if rand()<density(distance(p)) # provisional accept based on density
            if (popCount < 1) || !anycloserthan(Δ, p, P[1:popCount]) 
                popCount += 1    
                P[popCount] = p           # accept
                crowdCount = 0
            else
                crowdCount += 1
            end
        end
    end
    return (P, popCount)
end

# univariate sample from density on [0, rmax] (by rejection)
function bumpsample(n::Int64, density::Function, rmax::Base.Float64)

    x = Vector{Float64}(undef,n)
    count = 0
    while count < n
        candidate = Float32(rand()*rmax)
        if rand() < density(candidate)
            count += 1
            x[count] = candidate
        end
    end

    return x
end


        


# add label (coloured points)
function label(tr::Trichoplax, labelName::String, n::Int64, density::Function, Δ::Float64,
              size::Float64 = 25.0, color::RGBA=RGBA(.8, .4, .2, 1.))

    (P, popCount) = radialsample(n, density, tr.radius[], Δ)
    trichoplax.nLabels[] += 1
    trichoplax.label_handle[trichoplax.nLabels[]] = 
        scatter!(ax,tr.location.+P[1:popCount], markersize = size, color = color)
    trichoplax.label_name[trichoplax.nLabels[]] = labelName
    print("Label: ", labelName, popCount==n ? ". OK: " : ". WARNING: ")
    println(popCount, " of ", n, " points created" )
end

# radial density function r->bumpdensity(r::Float32, ...)
# bump is triangle between r and R with max 1 at midpoint, raised to power q
# q = 1 gives triangle, q>1 concentrates mass near centre, 
# 0<q<1 spreads mass outward,  q = 0 is uniform on (r,R)
function bumpdensity(x::Float32, q::Float64, r::Float64, R::Float64)
    x0 = (r+R)/2.0
    h = (R-r)/2.0
    if (x>r) && (x<=x0) 
        return ((x-r)/h)^q
    elseif (x>x0) && (x<R)
        return ((R-x)/h)^q
    else
        return 0.0f0
    end
end

# compute parameters (clumpmean, clumpsize, nclumps) of random-sized randomly scattered clumps
# a clump is a local bump (radial distribution around some p not at origin)
# clump means have a radial bump distribution, ΔC is min distance between clumps
# clump sizes (s) have a univariate bump distribution
function clumps(x::Float32, tr::Trichoplax, 
                        nClump::Int64, q::Float64,  r::Float64, R::Float64, ΔC::Float64,
                        mean_clumpsize::Float64, range_clumpsize::Float64)

    # clump mean locations from bump distribution
    # nb actual number of clumps (nC) may vary from requested number (nClump) 
    #    because of crowding
    (C, nC) = radialsample(nClump, x->bumpdensity(x, q, r, R), tr.radius[], ΔC)

    # clump sizes from bump distribution
    s = mean_clumpsize-range_clumpsize
    S = mean_clumpsize+range_clumpsize
    (sC, dummy) = bumpsample(nC, x->bumpdensity(x, q, s, S), S)

    return (C, sC, nC)


end



F = Figure(resolution = (worldWide*pxl_per_um,worldWide*pxl_per_um))
ax = Axis(F[1,1], aspect = 1)
xlims!(0, worldWide)
ylims!(0, worldWide)
hidedecorations!(ax)

# construct world

worldMap = zeros(pixelWide, pixelWide)


# construct and draw Trichoplax
trichoplax  = Trichoplax(tr_diameter/2.0, tr_location)

# Trox-2 Jacobs et al. (2004)
Trox2density = x-> bumpdensity(x, 3.0, 440.0, 500.0)
minSeparation = 0.1  # minimum separation between points um
nTrox2_Particles = 10000 # number of label points
Trox2_particlesize = 1.0
Trox2_particlecolor = RGBA(1.0, 0.4, 0.4, 1.0)
label(trichoplax, "Trox_2", nTrox2_Particles, Trox2density, minSeparation, 
     Trox2_particlesize, Trox2_particlecolor)

# trPaxB Hadrys et al. (2005)
#trPaxBdensity = x-> clumpdensity(x, 3.0, 400.0, 480.0)
# minSeparation = 0.1  # minimum separation between points um
# ntrPaxB_Particles = 100 # number of label points
# trPaxB_particlesize = 5.0
# trPaxB_particlecolor = RGBA(8.0, 0.4, 0.8, 1.0)
# label(trichoplax, "trPaxB", ntrPaxB_Particles, trPaxBdensity, minSeparation, 
#        trPaxB_particlesize, trPaxB_particlecolor)

display(F)