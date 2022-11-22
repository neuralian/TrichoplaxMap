using GLMakie
using Colors
using Random, Distributions
using Infiltrator

Random.seed!(121213)

# global parameters

pxl_per_um = 0.5      # pixels per micron display
tr_diameter = 1000.0  # um
worldHigh = 1200.0 # microns
MAXLABELS = 20
MAXCELLTYPES = 20
MAXCELLS = 5000   # maximum number of cells of any type
MAXBAX = 100  # max bacteria per cluster
GLANDCELL_SIZE = 5.f0
TASTERANGE = 12.
MARGIN_WIDTH = 100.  # marginal zone
EAT_THRESHOLD = 5   # for eat decision
# type 2 gland cells happen to be the second cell type to be created
# WARNING: Definition of GLANDCELL_TYPE2_INDEX will have to be changed
#          if cell type creation order changes (cos I'm too lazy to write code to make this happen)
GLANDCELL_TYPE2_INDEX = 2   


# Trichoplax Utilities
struct Trichoplax
    map::Array{Float64,2}
    radius::Vector{Float64} # nb declaring as array makes value mutable as diameter[]
    location::Observable{Point2f}
    edge_handle::Poly
    lipozone_handle::Poly
    nLabels::Vector{Int64}
    label_handle::Vector{Scatter}  # handles of scatterplots
    label_name::Vector{String}
    nCelltypes::Vector{Int64}   # number of cell types
    nCells::Vector{Int64}       # number of cells of each type
    cell_handle::Vector{Vector{Poly}}  # cell_handle[i][j] is handle of jth cell of type i 
    celltype_name::Vector{String}
    cell_location::Vector{Vector{Point2f}}
    taste_map::Observable{Vector{Point2f}}
    eat_trigger::Vector{Int64}     # eat when eat_trigger > EAT_THRESHOLD
    v::Vector{Point2f}      # velocity
end


function Trichoplax(radius::Float64, location::Observable{Point2f}, v0::Point2f=Point2f(0,0))
    plt_handle =  poly!(Circle(location[], radius), 
                  color = RGBA(.8, .7, .7, .25), strokecolor = RGB(.6, .6, .6), strokewidth = .25)
    lipozone_handle = poly!(Circle(location[], radius-MARGIN_WIDTH), 
    color = RGBA(1.0, .75, 1.0, .2) , strokecolor = RGB(.6, .6, .6), strokewidth = .25)
    empty_tastemap = Observable(Vector{Point2f}(undef,0))
    scatter!(empty_tastemap, color = RGBA(0., 1., 0., .35), markersize = 2*TASTERANGE)
    Trichoplax( zeros(1,1), 
                [radius], location, plt_handle, lipozone_handle,
                [0], Vector{Scatter}(undef,MAXLABELS), Vector{String}(undef,MAXLABELS),
                [0], Vector{Int64}(undef, MAXCELLTYPES),
                [ Vector{Poly{Tuple{Vector{Point{2, Float32}}}}}(undef,MAXCELLS) for _ = 1:MAXCELLTYPES], 
                Vector{String}(undef, MAXCELLTYPES),
                [Vector{Point2f}(undef, MAXCELLS) for _ in 1:MAXCELLTYPES],
                empty_tastemap,
                [0],
                [v0]
            )
             #   [0], [scatter!(Point2f(NaN, NaN))], [""])
end

# distance utilities, for placing cells and labels
distance(p::Point2f) = sqrt(p[1]^2 + p[2]^2)
distance(p::Point2f, q::Point2f) = distance(p-q)
distance(p::Point2f, Q::Vector{Point2f}) = [distance(p,q) for q in Q]
# TRUE if any points in Q are within Δ of p
anycloserthan(Δ::Float64, p::Point2f, Q::Vector{Point2f}) = 
               length(Q)>0 ? length(findall(distance(p,Q).<Δ))>0 : false

function anycell_closerthan(Δ::Float64, me::Point2f, whanau::Vector{Point2f}, trichoplax::Trichoplax)

    # whanau are other cells of the same kind, curently under construction
    cellnearby = anycloserthan(Δ, me, whanau)
    if !cellnearby
        for celltype in 1:trichoplax.nCelltypes[]
            if anycloserthan(Δ, me, trichoplax.cell_location[celltype])
                return true
            end
        end
    end
    return cellnearby
end

function withinanyclump(p::Point2f, C::Vector{Point2f}, S::Vector{Float64})
    # true if p is within a clump
    for i in 1:length(C)
        if distance(p, C[i])<S[i] return true
        end
    end
    return false
end


# step trichoplax forward at current velocity
function move(trichoplax::Trichoplax)

    trichoplax.location[] +=  trichoplax.v[]
    trichoplax.edge_handle[1] = Circle(trichoplax.location[], trichoplax.radius[]) 
    trichoplax.lipozone_handle[1] = Circle(trichoplax.location[], trichoplax.radius[]-MARGIN_WIDTH) 
     # move outline
     for i in 1:trichoplax.nLabels[]   # move labels
        trichoplax.label_handle[i][1] = trichoplax.label_handle[i][1][] .+= trichoplax.v[]  
   end
   for i in 1:trichoplax.nCelltypes[]   # move cells
        for j in 1:trichoplax.nCells[i]
            trichoplax.cell_handle[i][j][1] = trichoplax.cell_handle[i][j][1][] .+= trichoplax.v[]  
        end
    end
end

# move to specified point
function moveto(trichoplax::Trichoplax, p::Point2f)
    v = trichoplax.v[]
    trichoplax.v[] =  p-trichoplax.location[]
    move(trichoplax)
    trichoplax.v[] = v
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
# minimum distance Δ between all cells, including other cell types
function radialcellsample(n::Int64, density::Function,  rmax::Float64, Δ::Float64, trichoplax::Trichoplax)
    crowdCount = 0
    overCrowd = 25   # will exit if fail to find room to place a label consecutively this many times

    popCount = 0
    P = Array{Point2f,1}(undef, n)  # array of n points, initially undefined
    while (popCount < n) && (crowdCount < overCrowd)
        p = Point2f((2.0*rand(2).-1.0).*rmax)  # uniform random point in area containing the trichoplax
        if rand()<density(distance(p)) # provisional accept based on density
            if (popCount < 1) || !anycell_closerthan(Δ, p, P[1:popCount], trichoplax) 
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

function rotate(P::Vector{Point2f}, Θ)
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
            celltypename::String, size::Float32, color::RGBA, edgecolor::RGBA, edgewidth::Float64=1.0)
#@infiltrate
    trichoplax.nCelltypes[] += 1
    trichoplax.nCells[trichoplax.nCelltypes[]] = length(position)
    trichoplax.celltype_name[trichoplax.nCelltypes[]] = celltypename
    trichoplax.cell_location[trichoplax.nCelltypes[]] = position

    for j in 1:length(position)
        trichoplax.cell_handle[trichoplax.nCelltypes[]][j] = 
            poly!([trichoplax.location[]].+position[j].+rotate(size*shape, orientation[j]),
            color = color, strokecolor = edgecolor, strokewidth = edgewidth, overdraw = true)
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

# radial density with peak at origin (centre)
function radialdensity(x::Float32, q::Float64, R::Float64)

    if (x<R)
        return ((R-x)/R)^q
    end
    return 0.0
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

# cell shapes
 # crystal cell = trapezoid with origin at centre, "line of sight" axis +x 
#crystal_cell_shape =    [Point2f(-1., 0.8), Point2f(1., 1.25), Point2f(1., -1.25), Point2f(-1., -0.8)]   # trapezoid with origin at centre, "line of sight" axis +x 
function ampulla_shape(height::Float64=1.0, basalwidth::Float64=0.4, apicalwidth::Float64=0.30,
    neckheight::Float64=0.20, necklength::Float64=0.30, neckwidth::Float64=0.20, 
    chamfer::Float64 = 0.05)
# base shape symmetric around vertical axis, apical face up, origin at center of apical face
# vertices listed clockwise from top right
#        
#            16--1
#            /    \
#           15     2
#            |     |
#           14     3
#            \    /
#            13   4
#            |   |
#           12    5  
#           /      \
#          11       6
#           |       |
#          10       7
#           \      /
#           9-----8


[
Point2f(apicalwidth/2.0-chamfer, 0.0), Point2f(apicalwidth/2.0, -chamfer),
Point2f(apicalwidth/2.0, -neckheight+chamfer/2.0), 
Point2f(neckwidth/2.0, -neckheight-chamfer/2.0),
Point2f(neckwidth/2.0, -neckheight-necklength+chamfer/2.0), 
Point2f(basalwidth/2.0-chamfer, -neckheight-necklength-chamfer/2.0),
Point2f(basalwidth/2.0, -height+chamfer), 
Point2f(basalwidth/2.0-chamfer, -height),
Point2f(-basalwidth/2.0+chamfer, -height),
Point2f(-basalwidth/2.0, -height+chamfer),   
Point2f(-basalwidth/2.0+chamfer, -neckheight-necklength-chamfer/2.0),  
Point2f(-neckwidth/2.0, -neckheight-necklength+chamfer/2.0),                     
Point2f(-neckwidth/2.0, -neckheight-chamfer/2.0),
Point2f(-apicalwidth/2.0, -neckheight+chamfer/2.0),       
Point2f(-apicalwidth/2.0, -chamfer), Point2f(-apicalwidth/2.0+chamfer, 0.0)
]
end  






function Trox2()
    # Trox-2 Jacobs et al. (2004); 
    # ALSO TrDLx and TrMnx (patchy in same zone) Monteiro et al. (2006)
     # Dlx is ANTP class. Limb and craniofacial development; axonal sprouting and neuronal migration 
     # Mnx pancreatic insulin cells and motorneurons
     # Monteiro &c also found TrHmx but did not identify spatial distribution
     # Hmx labels sensory placodes in vertebrates
    Trox2density = x-> bumpdensity(x, 1.0, 425.0, 500.0)
    minSeparation = 0.1  # minimum separation between points um
    nTrox2_Particles = 18000 # number of label points
    Trox2_particle_location = radialsample(nTrox2_Particles, 
                            Trox2density, trichoplax.radius[], minSeparation)
    Trox2_particlesize = 2.0
    Trox2_particlecolor = RGBA(0.95, 0.8, 0.25, 1.0)
    label(trichoplax, Trox2_particle_location, "Trox_2", Trox2_particlesize, Trox2_particlecolor)
end


function get_PaxB_clumps()
        # returns clump means

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
    return trPaxB_clump
end


function PaxB(clump::Tuple{Vector{Point{2, Float32}}, Vector{Float64}, Float64})
    # trPaxB Hadrys et al. (2005)
    # returns clump means
    # trPaxB_nclumps = 36
    # trPaxB_clump_distribution_shape =2.0 # q, shape of clump distribution
    # trPaxB_clump_shape = 1.5  # qi, shape of each clump
    # trPaxB_r = 410.  # lower bound on clump location
    # trPaxB_R = 450.  # upper bound on clump location
    # #trPaxB_Δ = 28.  # min distance between clump means
    # trPaxB_mean_clumpsize = 50.
    # trPaxB_range_clumpsize = 40.  # clumpsize is a bump reaching this far above and below mean 
    # trPaxB_clump = make_clumps(trPaxB_nclumps, trPaxB_clump_distribution_shape, trPaxB_clump_shape, 
    #                 trPaxB_r, trPaxB_R, trPaxB_mean_clumpsize, trPaxB_range_clumpsize)
    #trPaxBdensity = x-> clumpdensity(x, 3.0, 400.0, 480.0)
    minSeparation = 0.1  # minimum separation between points um
    ntrPaxB_Particles = 20000 # number of label points
    TrPaxB_particle_location = clumpsample(ntrPaxB_Particles, clump[1], clump[2], 
                        clump[3], trichoplax.radius[], minSeparation)
    trPaxB_particlesize = 0.25
    trPaxB_particlecolor = RGBA(0.7, 0.1, 0.65, 1.0)
    label(trichoplax, TrPaxB_particle_location, "trPaxB", trPaxB_particlesize, trPaxB_particlecolor)

end


function chordinLike()

    TrChrd_density = x-> radialdensity(x, 2.0, 100.0)
    minSeparation = 0.1  # minimum separation between points um
    nTrChrd_Particles = 1000 # number of label points
    TrChrd_particle_location = radialsample(nTrChrd_Particles, 
                            TrChrd_density, trichoplax.radius[], minSeparation)
    TrChrd_particlesize = 2.0
    TrChrd_particlecolor = RGBA(.75, 0.25, 0.0, 1.0)
    label(trichoplax, TrChrd_particle_location, "TrChrd", TrChrd_particlesize, TrChrd_particlecolor)

end


function bmp()
    Trbmp_density = x-> bumpdensity(x, 1.5, 50., 475.)
    minSeparation = 0.1  # minimum separation between points um
    nTrbmp_Particles = 20000 # number of label points
    Trbmp_particle_location = radialsample(nTrbmp_Particles, 
                            Trbmp_density, trichoplax.radius[], minSeparation)
    Trbmp_particlesize = 2.0
    Trbmp_particlecolor = RGBA(0.85, 0.55, 0.0, 1.0)
    label(trichoplax, Trbmp_particle_location, "TrBmp", Trbmp_particlesize, Trbmp_particlecolor)
end



 #[Point2f(-1., 0.8), Point2f(1., 1.25), Point2f(1., -1.25), Point2f(-1., -0.8)] 
function crystals(trichoplax::Trichoplax, location::Vector{Point2f})
    # Crystal cells
    # trapezoid with origin at centre, "line of sight" axis +x 
    xtal_shape = [ Point2f(0.8, -1.0), Point2f(1.25, 1.), Point2f(-1.25, 1.), Point2f(-.8, -1.) ]   
    xtal_wobbleamount = π/12.
    xtal_wobble = rand(Normal(0.0, xtal_wobbleamount), length(location))
    crystal_size = 6.f0
    crystal_colour = RGBA(.85, 0.9, 1.0, 1.)
    outline_colour = RGBA(0.25, 0.25, 1.0, 1.0)
    outline_width = 1.
    #label(trichoplax, trPaxB_clump[1], "Crystals", crystal_size, crystal_colour, crystal_shape)

    cells(trichoplax, location,
        Float64.([ -atan(location[j]...) for j in 1:length(location)].+xtal_wobble), 
        xtal_shape, "crystals", crystal_size, crystal_colour, outline_colour, outline_width)
end



function ampullae(trichoplax::Trichoplax, n::Int64)


    # ampulla cells are default gland cell shape
    ampl_size = 24.0f0
    ampl_shape = ampulla_shape()  
    spacing =  10.0
    lateral_jitter = .1  # lateral shape jitter, sd as proportion of ampl_size
    long_jitter = .35     # length jitter, ... 

    ampl_colour =  RGBA(.1, 1.0, .5, 1.0)  
    outline_colour = RGBA(0., 0., 0., 1.0)
    outline_width = 0.25

    # n random locations around periphery, with at least spacing between, 
    # give up if crowded consecutive attempts fail to find a spot
    count = 0
    crowded = 25
    crowding = 0
    location = Vector{Point2f}(undef,n)
    while count<n && crowding<crowded
        Θ = 2π*rand()
        location[count+1] = trichoplax.radius[]*Point2f(cos(Θ), sin(Θ))
        if !anycloserthan(spacing, location[count+1], location[1:count])
            count += 1
            crowding = 0
        else
            crowding += 1
        end
    end

    trichoplax.nCelltypes[] += 1
    trichoplax.nCells[trichoplax.nCelltypes[]] = length(location)
    trichoplax.celltype_name[trichoplax.nCelltypes[]] = "ampulla"
    trichoplax.cell_location[trichoplax.nCelltypes[]] = location

    for j in 1:length(location)

        # jitter: apical face (1st 2 layers) fixed, 6 deeper layers have random jitter
        # that scales with distance from apical face
        Ljit = rand(truncated(Normal(0.0, lateral_jitter*ampl_size), -lateral_jitter*ampl_size, lateral_jitter*ampl_size), 6).*(1:6)./6.0
        Hjit = rand(truncated(Normal(0.0, long_jitter*GLANDCELL_SIZE), -long_jitter*ampl_size, long_jitter*ampl_size), 6).*(1:6)./6.0
        shape = ampl_size*ampl_shape
       # @infiltrate
        for k in 1:6
            shape[k+2] += Point2f(Ljit[k], Hjit[k])
            shape[15-k] += Point2f(Ljit[k], Hjit[k])
        end

       # @infiltrate
    trichoplax.cell_handle[trichoplax.nCelltypes[]][j] = 
        poly!([trichoplax.location[]].+location[j].+rotate(shape, -atan(location[j]...)),
        color = ampl_colour, strokecolor = outline_colour, strokewidth = outline_width)
    end


end


function glandcell_type2(trichoplax::Trichoplax, n::Int64)
   # Mayorova et al. (2019)

    gland2_colour =  RGBA(.1, 1.0, .5, 1.0)  
    gland2_size = 4.0f0
    outline_colour = RGBA(0., 0., 0., 1.0)
    outline_width = 0.25
    spacing = 10.0

    shape = [Point2f(cos(Θ), sin(Θ)) for Θ in (1:6)*2π/6]
    # nb because area element increases linearly with distance from center, 
    # quadratic radial density gives linear increase in spatial density 
    density = x-> x < (trichoplax.radius[]-spacing) ? (x/trichoplax.radius[])^2 : 0.0
    location = radialcellsample(n, density , trichoplax.radius[], spacing, trichoplax)

    trichoplax.nCelltypes[] += 1
    trichoplax.nCells[trichoplax.nCelltypes[]] = length(location)
    trichoplax.celltype_name[trichoplax.nCelltypes[]] = "Gland_T2"   
    trichoplax.cell_location[trichoplax.nCelltypes[]] = location

    for j in 1:length(location)
        trichoplax.cell_handle[trichoplax.nCelltypes[]][j] = 
        poly!([trichoplax.location[]].+location[j].+rotate(GLANDCELL_SIZE*shape, rand()[]*2π),
        color = gland2_colour, strokecolor = outline_colour, strokewidth = outline_width)
    end

end

function glandcell_type3(trichoplax::Trichoplax, n::Int64)
    # Mayorova et al. (2019)
 
     gland3_colour =  RGBA(.1, .7, .4, 1.0) 
     outline_colour = RGBA(0., 0., 1., 1.0)
     outline_width = 0.25
     spacing = 10.0
     margin_width = 100.
 
     shape = [Point2f(cos(Θ), sin(Θ)) for Θ in (1:6)*2π/6]
     # nb because area element increases linearly with distance from center, 
     # quadratic radial density gives linear increase in spatial density 
     density = x-> x < (trichoplax.radius[]-margin_width) ? 1.0 : 0.0
     location = radialcellsample(n, density , trichoplax.radius[], spacing, trichoplax)

     trichoplax.nCelltypes[] += 1
     trichoplax.nCells[trichoplax.nCelltypes[]] = length(location)
     trichoplax.celltype_name[trichoplax.nCelltypes[]] = "Gland_T3"   
     trichoplax.cell_location[trichoplax.nCelltypes[]] = location

     for j in 1:length(location)
         trichoplax.cell_handle[trichoplax.nCelltypes[]][j] = 
         poly!([trichoplax.location[]].+location[j].+rotate(GLANDCELL_SIZE*shape, rand()[]*2π),
         color = gland3_colour, strokecolor = outline_colour, strokewidth = outline_width)
     end
 
 end

 function glandcell_type1(trichoplax::Trichoplax, n::Int64)
    # Mayorova et al. (2019)
 
     gland1_colour =  RGBA(.2, .8, .75, 1.0) 
     outline_colour = RGBA(0., 0., 1., 1.0)
     outline_width = 0.25
     spacing = 10.0
     margin_width = 100.
 
     shape = [Point2f(cos(Θ), sin(Θ)) for Θ in (1:6)*2π/6]
     # nb because area element increases linearly with distance from center, 
     # quadratic radial density gives linear increase in spatial density 
     location = radialcellsample(n, 
            x->bumpdensity(x, 2.0, trichoplax.radius[]-margin_width, trichoplax.radius[]) , 
            trichoplax.radius[], spacing, trichoplax)

     trichoplax.nCelltypes[] += 1
     trichoplax.nCells[trichoplax.nCelltypes[]] = length(location)
     trichoplax.celltype_name[trichoplax.nCelltypes[]] = "Gland_T3"   
     trichoplax.cell_location[trichoplax.nCelltypes[]] = location

     for j in 1:length(location)
         trichoplax.cell_handle[trichoplax.nCelltypes[]][j] = 
         poly!([trichoplax.location[]].+location[j].+rotate(GLANDCELL_SIZE*shape, rand()[]*2π),
         color = gland1_colour, strokecolor = outline_colour, strokewidth = outline_width)
     end
 
 end

 function lipophil(trichoplax::Trichoplax, n::Int64)
    # Mayorova et al. (2019)
 
     lipo_colour =  RGBA(1.0, .75, 1.0, 1.0) 
     lipo_size = 5.0f0
     outline_colour = RGBA(0., 0., 1., 1.0)
     outline_width = 0.25
     spacing = 10.0
     cell_radius = 1.0
     chamfer = 0.4
 
   
     shape = [Point2f(cell_radius-chamfer, -cell_radius), Point2f(cell_radius, -cell_radius+chamfer),
              Point2f(cell_radius, cell_radius-chamfer),  Point2f(cell_radius-chamfer, cell_radius),
              Point2f(-cell_radius+chamfer, cell_radius), Point2f(-cell_radius, cell_radius-chamfer),
              Point2f(-cell_radius, -cell_radius+chamfer), Point2f(-cell_radius+chamfer, -cell_radius)]
     # nb because area element increases linearly with distance from center, 
     # quadratic radial density gives linear increase in spatial density 
     location = radialcellsample(n,  x-> x < (trichoplax.radius[]-MARGIN_WIDTH) ? 1.0 : 0.0 , 
                trichoplax.radius[], spacing, trichoplax)

     trichoplax.nCelltypes[] += 1
     trichoplax.nCells[trichoplax.nCelltypes[]] = length(location)
     trichoplax.celltype_name[trichoplax.nCelltypes[]] = "Gland_T3"   
     trichoplax.cell_location[trichoplax.nCelltypes[]] = location
 
     for j in 1:length(location)
         trichoplax.cell_handle[trichoplax.nCelltypes[]][j] = 
         poly!([trichoplax.location[]].+location[j].+rotate(lipo_size*shape, rand()[]*2π),
         color = lipo_colour, strokecolor = outline_colour, strokewidth = outline_width)
     end
 
 end

function scatter_bacteria(location::Point2f, nClumps::Int64, nPerClump::Int64)
    # clump means Gaussian scattered around location (x,y)
    # Poisson count per clump, each Gaussian scattered around its mean location

    betweenSD = 100.0
    withinSD = 12.0
    clump = [location .+ Point2f(rand(Normal(0.0, betweenSD), 2)) for _ in 1:nClumps]
    clumpSize = rand(Poisson(nPerClump), nClumps)

    p = Vector{Point2f}(undef, 0)   # empty vector
    for i in 1:nClumps
        for j in 1:clumpSize[i]
            push!(p, clump[i] + Point2f(rand(Normal(0.0, withinSD), 2) ) )
        end
    end
    scatter!(p, markersize = 5.0, color = :red)
 end


# indices of type 2 gland cells within Δ of a bacterium
function taste(trichoplax::Trichoplax, bacteria::Scatter, Δ::Float64)

    bp = bacteria[1][]   # vector of Point2f bacteria locations
    g2p = [trichoplax.location[]] .+ trichoplax.cell_location[GLANDCELL_TYPE2_INDEX]  # vector of gland cell locations
    iTaste = Vector{Int64}(undef, 0)
    for i in 1:length(g2p)
        if any(distance(g2p[i], bp).<Δ)  # any bacterium within Δ of gp2[i]
            push!(iTaste, i)             # ith gland cell tastes a bacterium
        end
    end
   return iTaste
end


function init_maps(trichoplax::Trichoplax)

    scatter!(trichoplax.taste_map, color = RGBA(0., 1., 0., 0.25), markersize = 25.)

end


# proportion of bacteria tasted in lipohil zone to bacteria there and 
# in the marginal zone in movement direction 
function ingested(trichoplax::Trichoplax, iTaste::Vector{Int64})

   if isempty(iTaste)
    println("ERROR in ingested(): No bacteria tasted")
    return -1
   end

    PT = [trichoplax.cell_location[2][i] for i in iTaste]  # tasting cells

    nIngested = sum(distance.(PT).<(trichoplax.radius[]-MARGIN_WIDTH)) # count bacteria in lipophil zone

    nIncoming = 0
    d = distance(trichoplax.v[])
    a = cos(π/4.0)
    for i in 1:length(iTaste)
         # v.p/|v||p| = cosine of angle between tasting cell and movement direction
        if  sum(trichoplax.v[].*PT[i])/(d*distance(PT[i]))>a  && 
                 distance(PT[i])>(trichoplax.radius[]-MARGIN_WIDTH)
            nIncoming += 1
        end
    end

    return nIncoming/(nIngested+nIncoming)
end







video_aspectratio = 16.0/12.0
video_worldHigh = 1600
video_screenpixels = (video_aspectratio*video_worldHigh*pxl_per_um, video_worldHigh*pxl_per_um)

# for video_worldHigh
worldHigh = video_worldHigh
screenpixels = video_screenpixels
worldWide = video_aspectratio*video_worldHigh

# 16:9 aspect
F = Figure(resolution = screenpixels)
ax = Axis(F[1,1], aspect = DataAspect())
xlims!(0, worldWide)
ylims!(0, worldHigh)
hidedecorations!(ax)

# construct world

#worldMap = zeros(pixelWide, pixelWide)

bacteria_location = Point2f(200.0+1.25*tr_diameter, 100.0+0.9*tr_diameter )
bacteria = scatter_bacteria(bacteria_location, 24, 6)

# construct and draw Trichoplax
trichoplax_startpoint = Point2f(200.0+tr_diameter/2.0, 100.0+tr_diameter/2.0)
trichoplax_startpoint = Observable(Point2f(800., 800.))
trichoplax  = Trichoplax(tr_diameter/2.0, trichoplax_startpoint)

PaxB_clump = get_PaxB_clumps()  

showTFs = false
if showTFs
    PaxB(PaxB_clump)
    Trox2()
    chordinLike()
    bmp()
end

showCells = true
if showCells
    # crystal cells at PaxB clump means
    crystals(trichoplax, PaxB_clump[1])

    glandcell_type2(trichoplax, 800)

    glandcell_type3(trichoplax, 100)

    glandcell_type1(trichoplax,100)

    lipophil(trichoplax, 100)

    ampullae(trichoplax, 128)
end



display(F)
Random.seed!(4141)
VIDEO = true
if VIDEO
    trichoplax.v[] = Point2f(2.5,0.0)     # initial velocity
    record(F, "trich1.mkv", 1:200) do i
    #while trichoplax.location[][1]<clumpx
        move(trichoplax)
        iTaste = taste(trichoplax, bacteria, TASTERANGE)

        if isempty(iTaste)
            trichoplax.taste_map[] = Vector{Point2f}(undef,0)
        else
            # set velocity towards mean of tasted bacteria
            r = mean([trichoplax.cell_location[GLANDCELL_TYPE2_INDEX][j] for j in iTaste]) 
            trichoplax.v[] = 0.95*trichoplax.v[] +  0.05*2.5*r/trichoplax.radius[] 

            # display signal release from type 2 cells that taste bacteria
            trichoplax.taste_map[] = [trichoplax.location[]].+
                [trichoplax.cell_location[GLANDCELL_TYPE2_INDEX][j] for j in iTaste]

            # activate lipophil cells if the number of tasting cells in lipophil zone
            # is sufficiently greater than number in leading edge of marginal zone
            # ie there is a diminishing return for continuing to move vs stopping to eat
            if ingested(trichoplax, iTaste) < 0.25
                # for each tasing cell in lipophil zone
                # activate nearby lipophil cells
                trichoplax.eat_trigger[] += 1
                if trichoplax.eat_trigger[]>EAT_THRESHOLD
                   println("EAT!")
                end
            end
        end
        trichoplax.v[] = trichoplax.v[] + Point2f(rand(Normal(0.0, 0.1),2))
    end
end


