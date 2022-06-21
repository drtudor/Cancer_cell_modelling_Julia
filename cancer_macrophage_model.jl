using Agents, Random
using InteractiveDynamics
using GLMakie
using Statistics: mean


mutable struct Macrophage <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    eaten:: Int
end


# Function to create model domain

function create_model(;n_macrophages = 100, density = 0.1, griddims = (100,50), random_seed = 123)
    space = GridSpace(griddims; periodic = false, metric = :euclidean)
    rng = MersenneTwister(random_seed)

    #Domain = ABM(Macrophage, space; scheduler = Schedulers.randomly,rng, 
    #    properties = (cancer_cells = zeros(Int,griddims),cancer_spaces = false(griddims),))

        properties = (
        cancer_space = falses(griddims),
        cancer_cells = zeros(Int, griddims),
    )
    model = ABM(Macrophage, space; properties, rng, scheduler = Schedulers.randomly)
    id = 0 
    for _ in 1:n_macrophages
        id += 1
        macrophage = Macrophage(id, (0,0),0)
        add_agent!(macrophage,model)
    end 

    for p in positions(model)
        if rand(model.rng) < density
            cancer_space = true
            cancer_cells =  cancer_space ? 1 : 0
            model.cancer_cells[p...] = cancer_cells
            model. cancer_space[p...] =  cancer_space
        end 
    end 
    return model
end 

# Initiate a version of the model
Model = create_model()
function agent_step!(macrophage::Macrophage, model)
    macrophage_step!(macrophage,model)
end 


function macrophage_step!(agent, model)
    walk!(agent, rand, model)
    macrophage_eat!(agent,model)
end 

function macrophage_eat!(agent,model)
    if model.cancer_space[agent.pos...] 
        agent.eaten += 1
        model.cancer_space[agent.pos...] = false
        model.cancer_cells[agent.pos...] = 0
    end
end 

function available_food!(agent,model)
    newsite = agent.pos
    neighbors = nearby_positions(agent.post,model,2)
    empty = collect(empty_positions(model))
    if length(empty > 0)
        available_cells = (model.cancer_cells[x,y] for (x,y) in empty)
        max_value = maximum(available_cells)
        if max_value > 0
            cell_sites_inds = findall(x -> x == max_value, collect(available_cells))
            cell_sites = empty[cell_sites_inds]
            for dia in 1:2
                np = nearby_positions(agent.pos,model,dia)
                suitable = intersect(np,cell_sites)
                if length(suitable) > 0
                    newsite = rand(model.rng, suitable)
                    break
                end 
            end 
            newsite != agent.pos && move_agent!(agent, newsite, model)
        end 
    end 
    agent.eaten += 1 
    model.cancer_space[agent.pos...] = false
    model.cancer_cells[agent.pos...] = 0
end 

function cancer_growth_step!(Model)
    for I in findall(isequal(1), Model.cancer_cells)
        for idx in nearby_positions(I.I,Model)
            if Model.cancer_cells[idx...] == 0
                if rand(Model.rng) < 0.01
                    Model.cancer_cells[idx...] = 1
                    Model.cancer_space[idx...] = true 
                else
                    Model.cancer_cells[idx...] = 0
                    Model.cancer_space[idx...] = false
                end 
            end

            
        end 
        
        if rand(Model.rng) < 0.01
            Model.cancer_cells[I] = 2 
            Model.cancer_space[I] = true
        end
        
    end 
    """
    for I in findall(isequal(0), Model.cancer_cells)
        if rand(Model.rng) < 0.1
            Model.cancer_cells[I] = 1
            Model.cancer_spacer[I] = true
        end 
    end 
    """
end 

step!(Model, agent_step!, cancer_growth_step!, 1)

count(t == 1 for t in Model.cancer_cells)
    
Model = create_model(density = 0.01)
#steps = 5

#step!(Model, agent_step!, cancer_growth_step!, steps)

adata, _ = run!(Model, agent_step!, cancer_growth_step!, 1500, adata = [:eaten])
adata[1:10,:]

figure = Figure(resolution = (600, 600))
title_text = Observable("Distribution of cancer cells per macrophage\nStep 1")
figure[1, 1] = Label(figure, title_text; textsize=30, tellwidth=false)
ax = figure[2, 1] = Axis(figure; xlabel="Cancer cells", ylabel="Number of macrophages")
histdata = Observable(adata[adata.step .== 1500, :eaten])
hist!(ax, histdata; bar_position=:step)
#record(figure, "'CC_dist.mp4", 0:20; framerate=3) do i
    #histdata[] = adata[adata.step .== 0, :eaten]
title_text[] = "Distribution of cancer cells per macrophage, step =1500"
xlims!(ax, (0, max(histdata[]...)))
ylims!(ax, (0, 50))
#end



"""
plotkwargs = (
    add_colorbar = false,
    heatarray = :cancer_cells,
    heatkwargs = (
        colorrange = (0,2),
        colormap = cgrad([:white, :black, :grey]; categorical = true),
    ),
)
fig, _ = abmplot(Model; plotkwargs...)
fig


"""
"""
All the code below creates interactive plotting to see how the macrophages eat the cancer cells:
 - They only eat dead (grey) and alive (black) cancer cells, anything which is white is a free space
 - When they have eaten the cancer cell (whether alive or dead), they space they occupied becomes free again, and if there is enough time the cancer can regrown there 

"""
plotkwargs = (
    ac=:green, as = 5, am = :circle,
    heatarray = :cancer_cells,
    heatkwargs = (
        colorrange = (0,2),
        colormap = cgrad([:white, :black, :grey]; categorical = true),
        
    ),
)

params = Dict(
        :n_macrophages => 0:100,
        :density => 0:0.01:1,

)

model = create_model(n_macrophages = 100,density = 0.5)
cancer_count(model) = count(model.cancer_space)
cancer_dead(model) = count(model.cancer_cells == 2)
adata = [(:eaten, mean)]
alabels = ["Eaten (Average)"]
mdata = [cancer_count]
mlabels = ["Number of cancer cells"] 
fig, p = abmexploration(model; 
    agent_step! = agent_step!, model_step! = cancer_growth_step!,params, plotkwargs...,
    adata, alabels,mdata,mlabels
)

fig
save("figure.png", fig)






#macrophage(a) = a.type == :Macrophage
model = create_model(n_macrophages = 50,density = 0.01)

n = 1500
#adata = [(macrophage, macrophage.eaten)]
adata = [(:eaten, mean)]
mdf = run!(model, agent_step!, cancer_growth_step!, n; adata,mdata)
function plot_population_timeseries_cancer(mdf)
    figure = Figure(resolution = (600,600))
    ax = figure[1, 1] = Axis(figure; xlabel = "Step", ylabel = "Number of cancer cells")
    cancercellsl = lines!(ax,mdf.step, mdf.cancer_count, color = :green)
    #eatens = scatter!(ax,mdf[1].step, mdf[1].mean_eaten, color = :green)
    #figure[1, 2] = Legend(figure, [can\cercellsl, eatens], ["Cancer Cells", "Mean number of eaten cancer cells"])
    figure
    save("Cancer_cellls_50_macro.png", figure)

end
function plot_population_timeseries_macro(mdf)
    figure = Figure(resolution = (600, 600))
    ax = figure[1, 1] = Axis(figure; xlabel = "Step", ylabel = "Macrophage cancer cell consumption (mean)")
    cancercellsl = lines!(ax,mdf.step, mdf.mean_eaten, color = :green)
    #eatens = scatter!(ax,mdf[1].step, mdf[1].mean_eaten, color = :green)
    #figure[1, 2] = Legend(figure, [cancercellsl, eatens], ["Cancer Cells", "Mean number of eaten cancer cells"])
    figure
    save("Macrophage_eaten_50_macro.png", figure)

end


plot_population_timeseries_macro(mdf[1])
plot_population_timeseries_cancer(mdf[2])

model = create_model(density = 0.01)

abmvideo(
    "cancer_model.mp4",
    model, agent_step!, cancer_growth_step!;
    
    title = "Cancer model", framerate = 150,frames = 1500,
    adata, alabels,mdata,mlabels,
    plotkwargs...
)


