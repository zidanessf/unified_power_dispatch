include("library/RDDP.jl")
using .RDDP,PowerModels,CSV,Interpolations,DataFrames,Plots,CPLEX
# read input data
N_stage,K_seg = 25,5
RU,RD = 24/(N_stage-1),24/(N_stage-1)
UT,DT = 4 * (N_stage-1)/24, 4 * (N_stage-1)/24
silence()
case = PowerModels.parse_file("datasets/pglib_opf_case118_ieee.m")
# case = PowerModels.parse_file("datasets/case6ww.m")
for gen in keys(case["gen"])
    case["gen"][gen]["pmin"] = max(case["gen"][gen]["pmin"],case["gen"][gen]["pmax"]*0.35)
end
# wind_power,load = DataFrame(),DataFrame()
wind_power = CSV.read("datasets/wind_2020.csv",DataFrame)[1:24,:]
load = CSV.read("datasets/load_data.csv",DataFrame)[1:24,:load]
gens = keys(case["gen"])
loadmax = sum(max(case["gen"][gen]["pmax"],0) for gen in gens)
load = load * loadmax/maximum(load)
load_ipt = extrapolate(interpolate(load,BSpline(Linear())),Flat())
load = [load_ipt(t*24/(N_stage-1)) for t in 2:N_stage]
r = 0.6 * loadmax/maximum(wind_power[!,:wf1])
alpha = 0.3
wind_ipt1 = extrapolate(interpolate(wind_power[!,:wf1],BSpline(Linear())),Flat())
Pw_max = [r*(1+alpha)*wind_ipt1(t*24/(N_stage-1)) for t in 2:N_stage]
Pw_min = [r*(1-alpha)*wind_ipt1(t*24/(N_stage-1)) for t in 2:N_stage]
Pw = [r*wind_ipt1(t*24/(N_stage-1)) for t in 2:N_stage]
ess = [1,2]
re = 0.2 * maximum(Pw)
Emax,Emin,Pmax,E0 = re*[1,1],[0,0],re*[0.5,0.5],[0,0]
msro = @timev RDDP.buildMultiStageRobustModel(
    N_stage = N_stage,
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer
    ),
    MaxIteration = 1000,
    MaxTime = 600,
    Gap = 0.005,
    AutoStageSelection = true
) do ro::JuMP.Model,t
    if t >= 2
        @variable(ro,pw_max,RDDP.Uncertain,lowesr_bound=Pw_min[t-1],upper_bound=Pw_max[t-1])
    end
    # defining operator's strategy
    @variable(ro,Pg[gen in gens],RDDP.State,initial_value = case["gen"][gen]["pmin"],spread = 2:N_stage)
    @variable(ro,pw,lower_bound=0)
    @variable(ro,Pgcost[gens],lower_bound=0)
    @variable(ro,Ses[e in [1,2]],RDDP.State,initial_value = E0[e],spread = 2:N_stage)
    @variable(ro,s1,lower_bound=0)
    @variable(ro,s2,lower_bound=0)
    @variable(ro,Pes[e in [1,2]],lower_bound=-Pmax[e],upper_bound=Pmax[e])
    if t == 1
        @variable(ro,st[gen in gens, t in 2:N_stage],Bin,RDDP.State,initial_value = 1,spread = 1:N_stage)
        @variable(ro,st_aux[gen in gens, t in 2:N_stage],lower_bound=0)
    else
        @variable(ro,st[gen in gens, t in 2:N_stage],RDDP.State,initial_value = 1,spread = 1:N_stage)
    end
    if t >= 2
        @constraint(ro,pw <= pw_max)
        for gen in gens # constraints of thermal units
            @constraint(ro,Pg[gen].out >= st[gen,t].out * case["gen"][gen]["pmin"]) # power output upper limit 
            @constraint(ro,Pg[gen].out <= st[gen,t].out * case["gen"][gen]["pmax"]) # power output lower limit
            if t >= 3
                @constraint(ro,Pg[gen].out - Pg[gen].in <= st[gen,t-1].out * RU * case["gen"][gen]["pmax"] + (1-st[gen,t-1].out) * case["gen"][gen]["pmax"]) # ramp up limit
                @constraint(ro,Pg[gen].in - Pg[gen].out <= st[gen,t].out * RD * case["gen"][gen]["pmax"] + (1-st[gen,t].out) * case["gen"][gen]["pmax"]) # ramp down limit
            end
            if length(case["gen"][gen]["cost"]) == 3 # piecewise linear cost 
                for k = 1:K_seg
                    pk = case["gen"][gen]["pmax"]*k/K_seg
                    costk = (case["gen"][gen]["cost"][2])*pk + (case["gen"][gen]["cost"][1])*pk^2
                    ratek = 2*case["gen"][gen]["cost"][1]*pk
                    @constraint(ro,Pgcost[gen] >= 24/(N_stage-1) * (costk + ratek*(Pg[gen].out - pk)))
                end
            elseif length(case["gen"][gen]["cost"]) == 2
                @constraint(ro,Pgcost[gen] == 24/(N_stage-1) * case["gen"][gen]["cost"][1]*Pg[gen].out)
            elseif length(case["gen"][gen]["cost"]) == 0
                nothing
            else
                error("cost term must be 2 or 3")
            end
        end
        for e in ess # constraints of storages
            @constraint(ro,Ses[e].out == Ses[e].in - 24/(N_stage-1)*Pes[e]) #state transition
            @constraint(ro,Ses[e].out <= Emax[e])
            @constraint(ro,Ses[e].out >= Emin[e])
        end
        for gen in gens
            for s in 2:N_stage
                @constraint(ro,st[gen,s].out == st[gen,s].in)
            end
        end
        @constraint(ro,sum(Pg[gen].out for gen in gens) + sum(Pes) + pw + s1 - s2 == load[t-1]) # power balance
        @objective(ro,Min,sum(Pgcost) + 10000*s1 + 10000*s2 + 2000 * (pw_max - pw))
    else # t = 1
        for s in 2:N_stage
            for gen in gens # minimum stop start time
                if s >= 3
                    if s <= N_stage - UT + 1
                        @constraint(ro,sum(st[gen,k].out for k in s:s+UT-1) >= UT*(st[gen,s].out - st[gen,s-1].out))
                    else
                        @constraint(ro,sum(st[gen,k].out - (st[gen,s].out - st[gen,s-1].out) for k in s:N_stage) >= 0)
                    end
                    if s <= N_stage - DT + 1
                        @constraint(ro,sum(1 - st[gen,k].out for k in s:s+DT-1) >= DT*(st[gen,s-1].out - st[gen,s].out))
                    else
                        @constraint(ro,sum(1 - st[gen,k].out - (st[gen,s-1].out - st[gen,s].out) for k in s:N_stage) >= 0)
                    end
                    @constraint(ro,st_aux[gen,s] >= st[gen,s].out - st[gen,s-1].out)
                end
            end
        end
        @objective(ro,Min,2000 * 24/(N_stage-1) * sum(sum(st_aux)))
    end
end
RDDP.train(msro)
objectives,dist,min_obj,obj_std = RDDP.evaluate(msro)
plot([[value(msro.lower[2][:st][gen,t].out) for t in 2:N_stage] for gen in gens],label=hcat(["G$gen" for gen in gens]...),marker=true)
title!("Unit Commitment")
xlabel!("Time(hour)")
ylabel!("ON-OFF")
# plot([[value(msro.lower[t][:cost_to_go]) for t in 2:N_stage-1],[value(msro.upper[t][:cost_to_go]) for t in 2:N_stage-1]],label=["upper bound" "lower bound"])
# plot([msro.lower[t][:penalty] for t in 1:N_stage])
# TODO: do some comparison with AARUC and independent AARUC + RDDP via matlab xprog