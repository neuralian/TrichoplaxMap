using GLMakie
using Colors
using Random
using Infiltrator

Random.seed!(121213)

# global parameters
worldWide = 1600.0 # microns
pxl_per_um = 0.5      # pixels per micron display
tr_diameter = 1000.0  # um
MAXLABELS = 20
MAXCELLTYPES = 20
MAXCELLS = 1000   # maximum number of cells of any type

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
    nCelltypes::Vector{Int64}   # number of cell types
    nCells::Vector{Int64}       # number of cells of each type
    cell_handle::Vector{Vector{Poly}}  # cell_handle[i][j] is handle of jth cell of type i 
    celltype_name::Vector{String}
end


function Trichoplax(radius::Float64, location::Point2f)
    plt_handle =  poly!(Circle(location, radius), 
                  color = RGBA(.6, .6, .6, .5), strokecolor = RGB(.6, .6, .6), strokewidth = .5)
    Trichoplax( zeros(1,1), 
                [radius], [location], plt_handle, 
                [0], Vector{Scatter}(undef,MAXLABELS), Vector{String}(undef,MAXLABELS),
                [0], Vector{Int64}(undef, MAXCELLTYPES),
                [ Vector{Poly}(undef,MAXCELLS) for _ = 1:MAXCELLTYPES], 
                Vector{String}(undef, MAXCELLTYPES) 
            )
             #   [0], [scatter!(Point2f(NaN, NaN))], [""])
end

distance(p::Point2f) = sqrt(p[1]^2 + p[2]^2)
distance(p::Point2f, q::Point2f) = distance(p-q)
distance(p::Point2f, Q::Vector{Point2f}) = [distance(p,q) for q in Q]
# TRUE if any points in Q are within Δ of p
anycloserthan(Δ::Float64, p::Point2f, Q::Vector{Point2f}) = length(findall(distance(p,Q).<Δ))>0

function withinanyclump(p::Point2f, C::Vector{Point2f}, S::Vector{Float64})
    # true if p is within a clump
    for i in 1:length(C)
        if distance(p, C[i])<S[i] return true
        end
    end
    return false
end


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
    print(popCount==n ? "OK: " : "WARNING: ")
    println(popCount, " of ", n, " points created" )
    return P[1:popCount]
end

# sample of n points p::point2f from density(r) defined on (0,rmax), where r is distance of p from origin
# minimum distance Δ between sample points
function radialclumpsample(n::Int64, density::Function,  rmax::Float64, S::Vector{Float64})
    crowdCount = 0
    overCrowd = 25   # will exit if fail to find room to place a label consecutively this many times

    popCount = 0
    P = Array{Point2f,1}(undef, n)  # array of n points, initially undefined
    while (popCount < n) && (crowdCount < overCrowd)
        p = Point2f((2.0*rand(2).-1.0).*rmax)  # uniform random point in area containing the trichoplax
        if rand()<density(distance(p)) # provisional accept based on density
            if (popCount < 1) || !withinanyclump(p, P[1:popCount], 0.75*S[1:popCount]) 
                popCount += 1    
                P[popCount] = p           # accept
                crowdCount = 0
            else
                crowdCount += 1
            end
        end
    end
    print(popCount==n ? "OK: " : "WARNING: ")
    println(popCount, " of ", n, " points created" )
    return P[1:popCount]
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

function rotate(P::Vector{Point2f}, Θ::Float64)
    # rotate points around origin

    Q = Vector{Point2f}(undef,length(P))
    for i in 1:length(P)
        Q[i] = Point2f(P[i][1]*cos(Θ) - P[i][2]*sin(Θ), P[i][1]*sin(Θ)+P[i][2]*cos(Θ) )
    end

    return Q
end




        
# draw label points
function label(trichoplax::Trichoplax, P::Vector{Point2f}, labelName::String, 
              size::Float64 = 25.0, color::RGBA=RGBA(.8, .4, .2, 1.), marker=:circle)

    trichoplax.nLabels[] += 1
    trichoplax.label_handle[trichoplax.nLabels[]] = 
        scatter!(ax,trichoplax.location.+P, markersize = size, color = color)
    trichoplax.label_name[trichoplax.nLabels[]] = labelName
end


# draw cells
# orientation counter-clockwise
function cells(trichoplax::Trichoplax, 
            position::Vector{Point2f}, orientation::Vector{Float64}, shape::Vector{Point2f},
            celltypename::String, size::Float64, color::RGBA, edgecolor::RGBA)

    trichoplax.nCelltypes[] += 1
    trichoplax.celltype_name[ trichoplax.nCelltypes[]] = celltypename

    for j in 1:length(P)
        trichoplax.cell_handle[trichoplax.nCelltypes[]] = 
            poly!(trichoplax.location.+position[j].+rotate(Float32(size)*shape, orientation[j]))
    end

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

# compute parameters (clumpmean, clumpsize) of random-sized randomly scattered clumps
# a clump is a local bump (radial distribution around point C[i])
# clump means have a radial bump distribution with shape exponent q 
# clump sizes have univariate bump distribution with shape exponent qi
# ΔC is min distance between clumps
function make_clumps(nClump::Int64, q::Float64,  qi::Float64, r::Float64, R::Float64, 
                        mean_clumpsize::Float64, range_clumpsize::Float64)

    # clump sizes from bump distribution
    s = mean_clumpsize/4.0 #-range_clumpsize
    S = mean_clumpsize+range_clumpsize
    sC = bumpsample(nClump, x->bumpdensity(x, qi, s, S), S)

    # clump mean locations from bump distribution
    # nb actual number of clumps (nC) may vary from requested number (nClump) 
    #    because of crowding
    C = radialclumpsample(nClump, x->bumpdensity(x, q, r, R), R, sC)

    return (C, sC, qi)

end

# clump density at x defined by parameters clumpmean, clumpsize & q returned by make_clumps
# nb this is not normalized. Individual clumps have max==1 but the max height of superimposed
#    clumps is unknown.  
function clumpdensity(x::Point2f, cmean::Vector{Point{2, Float32}}, s::Vector{Float64}, q::Float64)

    f = 0.0
    for i in 1:length(cmean)
        d = distance(x, cmean[i])
        if d<s[i]    # point is in range of ith clump
            f += (1.0 - d/s[i])^q  # add density due to ith clump 
        end
    end
    return f
end



# sample of points p::Point2f from clump density 
# clump is tuple of clump locations and sizes returned by clumps()
# rmax is trichoplax radius (defining region for candidate selection)
# min distance Δ between sample points
function clumpsample(n::Int64, cmean::Vector{Point{2, Float32}}, s::Vector{Float64}, q::Float64, 
                     rmax::Float64, Δ::Float64)

    crowdCount = 0
    overCrowd = 25   # will exit if fail to find room to place a label consecutively this many times

    popCount = 0
    P = Array{Point2f,1}(undef, n)  # array of n points, initially undefined
    while (popCount < n) && (crowdCount < overCrowd)
        p = Point2f((2.0*rand(2).-1.0).*rmax)  # uniform random point in area containing the trichoplax
        if distance(p)<rmax && rand()<clumpdensity(p, cmean, s, q) # provisional accept based on density
            if (popCount < 1) || !anycloserthan(Δ, p, P[1:popCount]) 
                popCount += 1    
                P[popCount] = p           # accept
                crowdCount = 0
            else
                crowdCount += 1
            end
        end
    end
    return P[1:popCount]    
    

end


F = Figure(resolution = (worldWide*pxl_per_um,worldWide*pxl_per_um))
ax = Axis(F[1,1], aspect = 1)
xlims!(0, worldWide)
ylims!(0, worldWide)
hidedecorations!(ax)

# construct world

worldMap = zeros(pixelWide, pixelWide)


# construct and draw Trichoplax
# scatter!(tr_location, markersize = 20, color = RGB(.1, .4, .6))
trichoplax  = Trichoplax(tr_diameter/2.0, tr_location)


# Trox-2 Jacobs et al. (2004)
Trox2density = x-> bumpdensity(x, 1.5, 425.0, 500.0)
minSeparation = 0.1  # minimum separation between points um
nTrox2_Particles = 15000 # number of label points
Trox2_particle_location = radialsample(nTrox2_Particles, 
                            Trox2density, trichoplax.radius[], minSeparation)
Trox2_particlesize = 2.0
Trox2_particlecolor = RGBA(7.0, 0.4, 0.1, 1.0)
label(trichoplax, Trox2_particle_location, "Trox_2", Trox2_particlesize, Trox2_particlecolor)

# trPaxB Hadrys et al. (2005)
trPaxB_nclumps = 36
trPaxB_clump_distribution_shape =2.0 # q, shape of clump distribution
trPaxB_clump_shape = 1.5  # qi, shape of each clump
trPaxB_r = 410.  # lower bound on clump location
trPaxB_R = 450.  # upper bound on clump location
#trPaxB_Δ = 28.  # min distance between clump means
trPaxB_mean_clumpsize = 50.
trPaxB_range_clumpsize = 40.  # clumpsize is a bump reaching this far above and below mean 
trPaxB_clump = make_clumps(trPaxB_nclumps, trPaxB_clump_distribution_shape, trPaxB_clump_shape, 
                trPaxB_r, trPaxB_R, trPaxB_mean_clumpsize, trPaxB_range_clumpsize)
trPaxBdensity = x-> clumpdensity(x, 3.0, 400.0, 480.0)
minSeparation = 0.1  # minimum separation between points um
ntrPaxB_Particles = 20000 # number of label points
TrPaxB_particle_location = clumpsample(ntrPaxB_Particles, trPaxB_clump[1], trPaxB_clump[2], 
                    trPaxB_clump[3], trichoplax.radius[], minSeparation)
trPaxB_particlesize = 0.25
trPaxB_particlecolor = RGBA(0.7, 0.1, 0.65, 1.0)
label(trichoplax, TrPaxB_particle_location, "trPaxB", trPaxB_particlesize, trPaxB_particlecolor)


crystal_size = 12.0
crystal_colour = RGBA(.85, .9, 1.0, 1.0)
crystal_shape = :circle
label(trichoplax, trPaxB_clump[1], "Crystals", crystal_size, crystal_colour, crystal_shape)


display(F)