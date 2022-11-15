using GLMakie

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
    radius::Array{Float64,1} # nb declaring as array makes value mutable as diameter[]
    location::Array{Point2f,1}  
    edge_handle::Poly
end


function Trichoplax(radius::Float64, location::Point2f)
    plt_handle =  poly!(Circle(location, radius), 
                  color = :white, strokecolor = :black, strokewidth = 1)
    Trichoplax(zeros(1,1), [radius], [location], plt_handle)
end

distance(p::Point2f) = sqrt(p[1]^2 + p[2]^2)
distance(p::Point2f, q::Point2f) = distance(p-q)

function draw(trichoplax::Trichoplax)
    #trichoplax.plot_data_handle = 
    poly!(Circle(trichoplax.location[], trichoplax.radius[]), 
    color = :white, strokecolor = :black, strokewidth = 1)
end

function moveto(trichoplax::Trichoplax, p::Point2f)
   trichoplax.location[] = p
   trichoplax.edge_handle[1] = Circle(p, trichoplax.radius[])  # update plot data
end

function move(trichoplax::Trichoplax, dx::Point2f)
    moveto(trichoplax, trichoplax.location[]+dx)
end




############################




F = Figure(resolution = (worldWide*pxl_per_um,worldWide*pxl_per_um))
ax = Axis(F[1,1], aspect = 1)
xlims!(0, worldWide)
ylims!(0, worldWide)

# construct world

worldMap = zeros(pixelWide, pixelWide)


# construct and draw Trichoplax
trichoplax  = Trichoplax(tr_diameter/2.0, tr_location)




# sample of n points (type Point2f) from density(r) in disc of radius R
# maxDensity is an upper bound on density in the region (for rejection sampling)
# density is 
function radialDensitySample(tr::Trichoplax, n::Int64, radialDensity::Function)

    count = 0
    P = Array{Point2f,1}(undef, n)  # array of n points, initially undefined
    while count < n
        x = (2.0*rand()-1.0)*tr_diameter/2.0  # random x in disc
        y = (2.0*rand()-1.0)*tr_diameter/2.0  # random y in disc
        if (rand()<radialDensity(sqrt(x^2+y^2))) 
            count = count + 1
            P[count] = Point2f(x,y)
        end
    end
    return P
end


function uniformBand(r::Float64)
    if (r>100.) & (r<200.) return 1.0
    else return 0.0
    end
end

#P = radialDensitySample(1000, 500., uniformBand)

#tr_plthandle = draw(trichoplax)

display(F)