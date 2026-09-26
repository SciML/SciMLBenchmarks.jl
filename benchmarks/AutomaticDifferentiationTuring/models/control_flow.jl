#= 
This model illustrates dynamic control flow inside a model that depends on the
value of a random variable. This will cause problems with ReverseDiff's
compiled tapes, as a tape compiled at a given value of `a` may not be
appropriate for a different value of `a`.

To make sure that the table correctly reflects this issue, the preparation for
the gradient is carried out at a value of `a > 0`, and the gradient is
evaluated at a value of `a < 0`. See `main.jl` for more information.
=#

@model function control_flow()
    a ~ Normal()
    if a > 0
        b ~ Normal()
    else
        b ~ Beta(2, 2)
    end
end

model = control_flow()
