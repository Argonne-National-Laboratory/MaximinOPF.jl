
  using JuMP
  using CPLEX, Ipopt 

  mTest = Model(with_optimizer(CPLEX.Optimizer))
  #mTest = Model(with_optimizer(Ipopt.Optimizer))
  @variable(mTest, x)
  @variable(mTest, y)
  @objective(mTest, Min, (x-3)^2 + (y+2)^2)
  println(mTest)
  optimize!(mTest)
  println("Status: ",JuMP.termination_status(mTest))
  println("(x,y)=",JuMP.value(x)," ",JuMP.value(y))
  println("obj val: ",JuMP.objective_value(mTest))
  println("It works!")
