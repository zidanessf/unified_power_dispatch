total_solves = 0
number_of_reject = 0
sceanio_not_change_num = 0
FollowerValuePrev = 0
t_stage1,t_stage2,add_lower_bound_time,add_upper_bound_time,pass_time = 0,0,0,0,0
using CPLEX
Random.seed!(123)
mutable struct MultiStageRobustModel
    upper::Array{JuMP.Model}
    lower::Array{JuMP.Model}
    uncertain::Array{JuMP.Model}
    config::Dict
end
function add_upper_bound(msro,t,n_iter)# TODO: add state history 
    m = msro.upper[t]
    # log history
    for i in 1:length(m[:long_term_states])
        push!(m[:history][i],value(msro.upper[t+1][:long_term_states][i].in))
    end
    # add weight 
    dual_of_obj = @variable(m,lower_bound=0)
    push!(m[:λ_set],dual_of_obj)
    # add constraint
    # x_i - ∑y_k*x_{ki} + τu_i - τl_i == 0 for x_i in states variables
    for i in 1:length(m[:short_term_states])
        set_normalized_coefficient(m[:sum_states][i],dual_of_obj, - value(msro.upper[t+1][:short_term_states][i].in))
    end
    # modify constraint  ∑y_k == 1
    if n_iter == 1
        @constraint(m,sum_y,dual_of_obj == 1)
    end
    set_normalized_coefficient(msro.upper[t][:sum_y],dual_of_obj,1)
    v_upper = objective_value(msro.upper[t+1])
    # push!(m[:history][:cost_to_go],v_upper)
    # v_upper = sum(value(msro.upper[s][:cost_now]) for s in t+1:length(msro.upper))
    # modify objective
    set_normalized_coefficient(msro.upper[t][:upper_bound],dual_of_obj,-v_upper)
end
function add_lower_bound(msro,t)
    lower_cut = AffExpr()
    states_value,cut_slope = [],[]
    m = msro.lower[t]
    for i in 1:length(m[:state])
        if t >= m[:spread][i][1]
            p = dual(FixRef(msro.lower[t+1][:state][i].in))
            push!(states_value,value(msro.lower[t+1][:state][i].in))
            push!(cut_slope,p)
            m[:penalty] = max(m[:penalty],1.01*abs(p))#dynamic method
            add_to_expression!(lower_cut,p * (m[:state][i].out - value(msro.lower[t+1][:state][i].in)))
        end
    end
    v_lower = objective_value(msro.lower[t+1])
    cut_offset = v_lower - cut_slope' * states_value
    # v_lower = sum(value(msro.lower[s][:cost_now]) for s in t+1:length(msro.lower))
    add_to_expression!(lower_cut,v_lower)
    @constraint(m,m[:cost_to_go] >= lower_cut)
    # modify the cuts table
    # push!(m[:cuts_generated],lower_cut)
    # push!(m[:cuts_slopes],cut_slope)
    # push!(m[:cuts_offsets],cut_offset)
    # push!(m[:states_value],states_value)
    # for i in 1:length(m[:cuts_table])
    #     push!(m[:cuts_table][i],m[:cuts_slopes][i]'*states_value + m[:cuts_offsets][i])
    # end
    # I = length(m[:cuts_table]) + 1
    # m[:cuts_table][I] = [cut_slope'*m[:states_value][i] + cut_offset for i in 1:I]
end
function solve_upper(msro,t)
    m = msro.upper[t]
    penalty = msro.lower[t][:penalty]
    gap = 999
    while !isempty(m[:benders_cut])
        con = pop!(m[:benders_cut])
        delete(m,con)
    end
    n = 1
    while true
        optimize!(m)
        if isempty(m[:λ_set])
            break
        end
        lb = value(m[:ρ])
        λ_value = value.(m[:λ_set])
        obj = 0
        dual_λ = zeros(length(m[:λ_set]))
        for i in 1:length(m[:long_term_states])
            tmp = penalty * (value(m[:long_term_states][i].out) - m[:history][i]'*λ_value)
            obj += abs(tmp)
            if tmp >= 0
                dual_λ += - penalty * m[:history][i]
            else
                dual_λ += penalty * m[:history][i]
            end
        end
        cut = AffExpr(obj)
        ub = obj
        gap = abs(ub-lb)/(lb+1e-6)
        # @info("ub = $ub lb = $lb gap = $gap t = $t)")
        # @info(n)
        if gap <= 0.01 || abs(ub-lb) <= 0.01
            break
        end
        n += 1
        for i in 1:length(dual_λ)
            add_to_expression!(cut,dual_λ[i]*(m[:λ_set][i] - λ_value[i]))
        end
        con = @constraint(m,m[:ρ] >= cut)
        push!(m[:benders_cut],con)
    end
end
function solve_max_min(_max::JuMP.Model,_min::JuMP.Model,binder_max,binder_min;lipshitz_constant=100000,mu_upper_bound=1e7,max_iteration=20)
    try
        @assert _max[:mu_tag]
    catch
        @variable(_max,μ,upper_bound=mu_upper_bound)
        @variable(_max,λ[i=1:length(binder_max)],lower_bound=-lipshitz_constant,upper_bound=lipshitz_constant)
        @objective(_max,Max,μ + sum(λ.*binder_max))
        _max[:mu_tag] = true
    end
    # delete exisiting cuts
    try
        while !isempty(_max[:valid_cut])
            con = pop!(_max[:valid_cut])
            delete(_max,con)
        end
    catch #_max[:valid_cut] not exisit
        _max[:valid_cut] = []
    end
    num = 1
    cut_bound = 0
    while true
        Suppressor.@suppress_out optimize!(_max)
        UB = objective_value(_max)
        binder_value = value.(binder_max)
        for i in 1:length(binder_value)
            fix(binder_min[i],binder_value[i];force=true)
        end
        Suppressor.@suppress_out optimize!(_min)
        LB = objective_value(_min)
        @constraint(_max,_max[:μ] + sum(_max[:λ] .* binder_value) <= LB)
        for i in 1:length(binder_value)
            cut_bound = max(cut_bound,abs(dual(FixRef(binder_min[i]))))
        end
        num += 1
        # @info(UB,LB)
        if abs(UB-LB)/LB <= 0.01
            break
        end
        if num >= max_iteration
            @warn("maxmin algorithm did not converge, the current gap is $(abs(UB-LB)/LB)")
            break
        end
        for i in 1:length(binder_value)
            set_upper_bound(_max[:λ][i],min(lipshitz_constant,1.5 * cut_bound))
            set_lower_bound(_max[:λ][i],- min(lipshitz_constant,1.5 * cut_bound))
        end
    end
end
function select_valid_cuts(m::JuMP.Model)
    while !isempty(m[:valid_cuts]) # remove all cuts
        con = pop!(m[:valid_cuts])
        delete(m,con)
    end
    maximum_cut_index = []
    for j in 1:length(m[:cuts_table]) # select non-dominated cuts 
        push!(maximum_cut_index,argmax([m[:cuts_table][i][j] for i in 1:length(m[:cuts_table])]))
    end
    for k in unique(maximum_cut_index) # add non-dominated cuts to model
        con = @constraint(m,m[:cost_to_go] >= m[:cuts_generated][k])
        push!(m[:valid_cuts],con)
    end
end
function ForwardPassPrimal(msro,start,stop)
    global t_stage1,t_stage2,pass_time
    additional = Dict()
    additional[:gaphourly] = []
    # additional[:upper] = []
    gap = 0
    differ = 0
    sceanio_not_change = true
    tt0 = time()
    for t in start:stop #前推步骤
        # fix the previous decision
        # select_valid_cuts(msro.lower[t])
        for i in 1:length(msro.lower[t][:state])
            if t > max(msro.lower[t][:spread][i][1],start)
                fix(msro.upper[t][:state][i].in,value(msro.lower[t-1][:state][i].out))
                fix(msro.lower[t][:state][i].in,value(msro.lower[t-1][:state][i].out))
            end
        end
        # if sceanio_not_change_num >= 20
        #     wst_vertex = msro.upper[t][:wst_vertex]
        # else
        #     wst_vertex = vertex_enumeration_primal(msro,t)
        # end
        t0 = time()
        if msro.uncertain[t][:is_uncertain]
            i = Random.randperm(length(msro.upper[t][:vertex]))[1]
            msro.lower[t][:wst_case] = msro.upper[t][:vertex][i]   
            for i in 1:length(msro.lower[t][:wst_case])
                fix(msro.lower[t][:uncertain][i],msro.lower[t][:wst_case][i]; force=true)
            end
        end
        if t > 1
            t_stage2 += time() - t0
        end
        # printstyled("____________lower problem_____________";color=:green)
        Suppressor.@suppress_out optimize!(msro.lower[t])
        if t < length(msro.lower)
            msro.lower[t][:cost_to_go_value] = value(msro.lower[t][:cost_to_go])
        end
        @assert termination_status(msro.lower[t]) == MOI.OPTIMAL
        global total_solves
        total_solves += 2
        # @info(t1)    
    end
    pass_time += time() - tt0
    additional[:LowerBound] = value(msro.lower[start][:cost_to_go])
    return additional
end
function BackwardPassPrimal(msro,N_ITER,start,stop)
    global t_stage2,add_lower_bound_time,add_upper_bound_time
    for t in [stop-x+start for x in start+1:stop]#回代步骤
        t0 = time()
        # **update the overestimator**
            # **solve the updated upper problem**
        if msro.uncertain[t+1][:is_uncertain]
            wst_obj,wst_case = 0,[],[]
            for wst in msro.upper[t+1][:vertex]
                for i in 1:length(wst)
                    fix(msro.lower[t+1][:uncertain][i],wst[i]; force=true)
                end
                Suppressor.@suppress_out optimize!(msro.lower[t+1])
                @assert termination_status(msro.lower[t+1]) == MOI.OPTIMAL
                if objective_value(msro.lower[t+1]) >= wst_obj
                    wst_obj = objective_value(msro.lower[t+1])
                    wst_case = wst
                end
                msro.lower[t+1][:wst_case] = wst_case
            end
            # wst_vertex = msro.upper[t+1][:wst_vertex]
            for i in 1:length(msro.lower[t+1][:wst_case])
                fix(msro.lower[t+1][:uncertain][i],msro.lower[t+1][:wst_case][i]; force=true)
            end
        end
        # **solve the updated lower problem**
        Suppressor.@suppress_out optimize!(msro.lower[t+1])
        @assert termination_status(msro.lower[t+1]) == MOI.OPTIMAL
        global total_solves
        total_solves += 3
        v_lower = objective_value(msro.lower[t+1])
        if (v_lower > 1.001 * msro.lower[t][:cost_to_go_value]) ||  N_ITER == 1
            add_lower_bound(msro,t)
        end
        # @info(objective_value(msro.lower[t+1]),value(msro.lower[t][:cost_to_go]))
        # add_lower_bound(msro,t)
    end
end
function mature_stage(msro,additional)
    if minimum(additional[:gaphourly][1:end-1]) < 0.01 && maximum(additional[:gaphourly][1:end-1]) < 0.2 
        stop = findfirst(x->x<=0.01,additional[:gaphourly][1:end-1])
    else 
        stop = length(msro.lower)
    end
    # @info("stop = $stop")
    return stop
end
function train(msro)
    n_iter = 0
    global total_solves,sceanio_not_change_num,t_stage1,t_stage2,add_lower_bound_time,add_upper_bound_time,pass_time
    msro.config["train_intraday"] = false
    t_stage1,t_stage2,add_lower_bound_time,add_upper_bound_time,pass_time = 0,0,0,0,0
    total_solves,sceanio_not_change_num,gap_temp,no_improvement = 0,0,9999,0
    tmp_lower_bound,pre_quit = [],false
    stage2value = -99999
    start,stop = 1,length(msro.lower)
    additional = Dict()
    solution_status = DataFrames.DataFrame(LowerBound=[],mean=[],std=[],Time=[],TotalSolves=[]) 
    print_banner(stdout)
    print_iteration_header(stdout)
    t1 = time()
    break_point = 1
    while true
        # for t in 1:length(msro.upper)
        #     set_optimizer(msro.upper[t],CPLEX.Optimizer)
        # end
        n_iter += 1
        additional = ForwardPassPrimal(msro,start,stop)
        # ______________________________________________________________________________________________________________
        additional[:Iteration] = n_iter
        additional[:TotalSolves] = total_solves
        t2 = time()
        additional[:Time] = t2 - t1
        if n_iter == 1
            additional[:mean] = additional[:LowerBound]
            additional[:std] = Inf
        else
            lb_history_10 = [x for x in solution_status[max(break_point,end-10):end,:LowerBound]]
            additional[:mean] = Statistics.mean(vcat(additional[:LowerBound],lb_history_10))
            additional[:std] = Statistics.std(vcat(additional[:LowerBound],lb_history_10))/additional[:mean]
        end
        print_iteration(stdout,additional)
        if n_iter >= 3
            prev_mean = solution_status[end,:mean]
        end
        if n_iter >= 1
            push!(solution_status,additional)
        end
        # quit condition
        if n_iter >= msro.config["MaxIteration"] || additional[:Time] >= msro.config["MaxTime"] || no_improvement >= msro.config["no_improvement"]
            if n_iter >= msro.config["MaxIteration"]
                printstyled("Fail to converge. Maximum iterations reached ";color=:red)
            elseif additional[:Time] >= msro.config["MaxTime"]
                printstyled("Fail to converge. Maximum time reached";color=:red)
            else
                printstyled("Fail to converge. No further improvement";color=:red)
            end
            print("$n_iter iterations in $(additional[:Time]) seconds. \n")
            break
        end
        if n_iter >= 10 && !msro.config["train_intraday"] && abs(prev_mean - additional[:mean])/prev_mean <= msro.config["Gap"]/sqrt(10) && n_iter - break_point >= 5
            printstyled("Converged with stage 1 bounds met. ";color=:green)
            print("$n_iter iterations in $(round(additional[:Time];digits=3)) seconds. \n")
            break
        end
        # smart stage selection
        # if n_iter >= 3 && msro.config["AutoStageSelection"]
        #     if !msro.config["train_intraday"] && abs((prev_lower_bound - additional[:LowerBound])/additional[:LowerBound]) <= 1e-3 && additional[:LowerBound] > 0.1
        #         no_improvement += 1
        #     end
        #     if msro.config["train_intraday"]
        #         if additional[:std] <= msro.config["Gap"]
        #             if abs(additional[:LowerBound] - stage2value)/additional[:LowerBound] <= msro.config["Gap"] && n_iter - break_point >= 5
        #                 printstyled("Converged with stage 2 value not change. ";color=:green)
        #                 print("$n_iter iterations in $(round(additional[:Time];digits=3)) seconds. \n")
        #                 break
        #             else
        #                 stage2value = additional[:LowerBound]
        #             end
        #             msro.config["train_intraday"] = false
        #             break_point = n_iter + 1
        #             printstyled("---------------------training between stage 1 and 2...---------------\n")
        #         end
        #     end
        #     if no_improvement >= 3
        #         msro.config["train_intraday"] = true
        #         no_improvement = 0
        #         state_values = [value(msro.lower[1][:state][i].out) for i in 1:length(msro.lower[1][:state])]
        #         while !isempty(msro.lower[1][:reg_cons])
        #             con = pop!(msro.lower[1][:reg_cons])
        #             delete(msro.lower[1],con)
        #         end
        #         for i in 1:length(msro.lower[1][:state])
        #             con1 = @constraint(msro.lower[1],msro.lower[1][:abs_state_variable][i] >= msro.lower[1][:state][i].out - state_values[i])
        #             con2 = @constraint(msro.lower[1],msro.lower[1][:abs_state_variable][i] >= - msro.lower[1][:state][i].out + state_values[i])
        #             push!(msro.lower[1][:reg_cons],con1)
        #             push!(msro.lower[1][:reg_cons],con2)
        #         end
        #         @objective(msro.lower[1],Min,msro.lower[1][:original_objective] + 20 * n_iter * sum(msro.lower[1][:abs_state_variable]))
        #         break_point = n_iter + 1
        #         printstyled("-------------------------fixing stage 1...---------------------------\n")
        #     end
        #     if msro.config["train_intraday"]
        #         start = 2
        #     else
        #         start = 1
        #     end
        # end
        Suppressor.@suppress_out BackwardPassPrimal(msro,n_iter,start,stop)
    end
    # println("stage 1 solve time: $t_stage1")
    # println("stage 2 solve time: $t_stage2")
    # println("add lower bound time: $add_lower_bound_time")
    # println("add upper bound time: $add_upper_bound_time")
    return solution_status
end
function buildMultiStageRobustModel(creator::Function;
    N_stage::Int,optimizer,
    MaxIteration=100,
    MaxTime=3600,
    Gap=0.01,
    initial_penalty = 1e8,
    use_maxmin_solver=false,
    AutoStageSelection=false)
    msro = Suppressor.@suppress_out MultiStageRobustModel(
        [JuMP.Model(optimizer) for t in 1:N_stage],
        [JuMP.Model(optimizer) for t in 1:N_stage],
        [JuMP.Model(with_optimizer(Ipopt.Optimizer)) for t in 1:N_stage],
        Dict([("MaxIteration",MaxIteration),
        ("MaxTime",MaxTime),
        ("Gap",Gap),
        ("no_improvement",100),
        ("initial_penalty",initial_penalty),
        ("use_maxmin_solver",use_maxmin_solver),
        ("AutoStageSelection",AutoStageSelection)]))
    for t in 1:N_stage # construct upper bound problem
        m = msro.upper[t]
        m[:state] = []
        m[:initial_value] = []
        m[:uncertain] = []
        m[:spread] = []
        m[:long_term_states] = []
        m[:short_term_states] = []
        m[:history] = Dict()
        m[:benders_cut] = []
        m[:λ_set] = []
        # set_parameter(m[:uncertain_problem],"NonConvex",2)
        creator(m,t)
        for i in 1:length(m[:state])
            if m[:spread][i][1] == 1
                push!(m[:long_term_states],m[:state][i])
            else
                push!(m[:short_term_states],m[:state][i])
            end
        end
        for i in 1:length(m[:long_term_states])
            m[:history][i] = []
        end
        for i in 1:length(m[:state])
            if t <= m[:spread][i][1]
                fix(m[:state][i].in,m[:initial_value][i])
            end
        end
        m[:cost_now] = objective_function(m)
        if t < N_stage
            @variable(m,cost_to_go,lower_bound=0)
            @variable(m,τu[1:length(m[:short_term_states])],lower_bound=0)
            @variable(m,τl[1:length(m[:short_term_states])],lower_bound=0)
            @variable(m,ep,lower_bound=0)
            @variable(m,ρ,lower_bound=0)
            @constraint(m,sum_states[i=1:length(m[:short_term_states])],m[:short_term_states][i].out + τu[i] - τl[i] == 0)
            m[:extra_penalty] = @constraint(m,ep == sum(τu[i]*initial_penalty for i in 1:length(m[:short_term_states])) + sum(τl[i]*initial_penalty for i in 1:length(m[:short_term_states])))
            m[:upper_bound] = @constraint(m,cost_to_go >= ep + ρ)
            @objective(m,Min,objective_function(m) + cost_to_go)
            m[:cost_to_go_now] = cost_to_go
        end
        m[:converging] = false
        m[:penalty] = 0.01
    end
    for t in 1:N_stage # construct lower bound problem
        m = msro.lower[t]
        m[:state] = []
        m[:initial_value] = []
        m[:uncertain] = []
        m[:spread] = []
        m[:valid_cuts],m[:cuts_generated],m[:cuts_offsets],m[:cuts_slopes],m[:states_value],m[:cuts_table] = [],[],[],[],[],Dict()
        creator(m,t)
        for i in 1:length(m[:state])
            if t == m[:spread][i][1]
                fix(m[:state][i].in,m[:initial_value][i])
            end
        end
        m[:cost_now] = objective_function(m)
        if t < N_stage
            @variable(m,cost_to_go,lower_bound=0)
            @objective(m,Min,objective_function(m) + cost_to_go)
            m[:cost_to_go_now] = cost_to_go
        end
        if t == 1
            @variable(m,abs_state_variable[1:length(m[:state])],lower_bound=0)
            m[:reg_cons] = []
        end
        m[:original_objective] = objective_function(m)
        m[:converging] = false
        m[:penalty] = 0.01
    end
    for t in 1:N_stage # construct uncertainty parameter problem
        m = msro.uncertain[t]
        m[:state] = []
        m[:initial_value] = []
        m[:uncertain] = []
        m[:spread] = []
        creator(m,t)
        if isempty(m[:uncertain])
            m[:is_uncertain] = false
            continue
        else
            m[:is_uncertain] = true
        end
        for ctype ∈ JuMP.list_of_constraint_types(m)# delete non-uncertain constraints
            for c ∈ JuMP.all_constraints(m,ctype[1],ctype[2])
                try
                    if all(JuMP.normalized_coefficient(c,v) == 0 for v ∈ m[:uncertain])
                        delete(m,c)
                    end
                catch
                end
            end
        end
        for v ∈ JuMP.all_variables(m) # delete non-uncertain variables
            if v ∉ m[:uncertain]
                # fix(v,0;force=true)
                JuMP.delete(m,v)
            end
        end
        m[:valid_cut] = [] # store cuts generated by max-min algorithm
        poly = Polyhedra.polyhedron(m, CDDLib.Library())# compute vertice
        msro.upper[t][:vertex] = [v for v in Polyhedra.points(Polyhedra.vrep(poly))]
    end
    # L1_regularization(msro)
    return msro
end
function L1_regularization(msro::MultiStageRobustModel) # not working
    for t in 1:length(msro.lower)
        m = msro.lower[t]
        @variable(m,z_aux[1:length(m[:state])],lower_bound=0)
        @variable(m,z[1:length(m[:state])],lower_bound=0)
        @objective(m,Min,objective_function(m) + 100000*sum(z_aux))
        for i in 1:length(m[:state])
            @constraint(m,z_aux[i] >= z[i] - m[:state][i].in)
            @constraint(m,z_aux[i] >= - z[i] + m[:state][i].in)
        end
    end
end
function evaluate_under_worst_case(msro,worst_case)
    total_cost = 0
    for t in 1:length(worst_case)
        fixprev(msro,t)
        for (idx,gen) in enumerate(keys(msro.case_dict[:windfarm]))
            vertex = [v for v in worst_case]
            windpower = msro.data[t][:wind_power]
            # @assert abs(vertex[idx]) <= 0.5
            fix(msro.lower[t][:Pw_err][gen],vertex[t][idx]*windpower[gen])
        end
        Suppressor.@suppress_out optimize!(msro.lower[t])
        # println(value(m[t][:cost_now]))
        @assert termination_status(msro.lower[t]) == MOI.OPTIMAL #OPTIMAL
        total_cost += value(msro.lower[t][:cost_now])
    end
    return total_cost
    # ForwardPassPrimal(msro,1,length(msro.lower))
    # return sum(value(msro.lower[t][:cost_now]) for t in 1:length(msro.lower))
end
function evaluate_under_worst_case(leader,follower,worst_case)
    total_cost = 0
    for t in 1:length(worst_case)
        fixprev(leader,t)
        for (idx,gen) in enumerate(keys(leader.case_dict[:windfarm]))
            vertex = [v for v in worst_case]
            windpower = leader.data[t][:wind_power]
            # @assert abs(vertex[idx]) <= 0.5
            fix(leader.intraday[t][:Pw_err][gen],vertex[t][idx]*windpower[gen])
        end
        Suppressor.@suppress_out optimize!(leader.intraday[t])
        # println(value(m[t][:cost_now]))
        @assert termination_status(leader.intraday[t]) == MOI.OPTIMAL #OPTIMAL
        total_cost += value(leader.intraday[t][:cost_now])
        fixprev(follower,value.(leader.intraday[t][:binder]),t)
        Suppressor.@suppress_out optimize!(follower.intraday[t])
        total_cost += value(follower.intraday[t][:cost_now])
    end
    return total_cost
end
function evaluate(msro)
    Random.seed!(1235)
    n = 1000
    T = length(msro.lower)
    objectives = []
    worst_index,tmp_worst = 1,0
    for k in 1:n
        # total_cost[k] = value(dayahead[:cost_now])
        m = msro.lower
        total_cost = objective_value(msro.lower[1]) - value(msro.lower[1][:cost_to_go])
        for t in 2:T
            for i in 1:length(msro.lower[t][:state])
                if t > msro.lower[t][:spread][i][1]
                    fix(msro.lower[t][:state][i].in,value(msro.lower[t-1][:state][i].out))
                end
            end
            if msro.uncertain[t][:is_uncertain]
                i = Random.randperm(length(msro.upper[t][:vertex]))[1]
                msro.lower[t][:wst_case] = msro.upper[t][:vertex][i]
                for i in 1:length(msro.lower[t][:wst_case])
                    fix(msro.lower[t][:uncertain][i],msro.lower[t][:wst_case][i]; force=true)
                end
            end
            Suppressor.@suppress_out optimize!(m[t])
            # println(value(m[t][:cost_now]))
            @assert termination_status(m[t]) == MOI.OPTIMAL #OPTIMAL
            total_cost += value(m[t][:cost_now])
        end
        if total_cost > tmp_worst
            tmp_worst = total_cost
            worst_index = k
        end
        push!(objectives,total_cost)
    end
    dist,min_obj,obj_std = peaksOverThresholdEstimator(objectives,0.95)
    return objectives,dist,min_obj,obj_std
end
function peaksOverThresholdEstimator(data::Array,quantile_α)
    u = Statistics.quantile(data,quantile_α)
    extreme_values = filter(x->x>=u,data)
    nev = (extreme_values .- minimum(extreme_values))/Statistics.std(extreme_values)
    N = length(nev)
    function f!(F,x)
        F[1] = 1/x[1] - (1/x[2] + 1) * 1/N * sum(nev[i]/(1 + x[1]*nev[i]) for i in 1:N)
        F[2] = 1/N * sum(log(1 + x[1]*nev[i]) for i in 1:N) - x[2]
    end
    x = NLsolve.nlsolve(f!,[0.1,0.1],autodiff=:forward)
    x = x.zero
    dist = Distributions.GeneralizedPareto(x[2]/x[1],x[2])
    return dist,minimum(extreme_values),Statistics.std(extreme_values)
end
