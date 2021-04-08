using CPLEX,JuMP
m = Model(with_optimizer(CPLEX.Optimizer))
@variable(m,x[1:10000] >= 0)
@constraint(m,x.<=10)
@constraint(m,sum(x[i]*cos(i) for i in 1:10000) <= 10)
@constraint(m,sum(x[i]*i for i in 1:10000) == 10)
@objective(m,Min,sum(sin(i)*x[i] for i in  1:10000))
printstyled("----------------solve model with variable unfixed----------------";color=:green)
optimize!(m)
printstyled("----------------solve model with variable fixed----------------";color=:green)
set_optimizer(m,CPLEX.Optimizer)
for i in 5000:10000
    fix(x[i],0;force=true)
end
for i in 2500:5000
    @constraint(m,x[i] == x[i+2500])
end
optimize!(m)