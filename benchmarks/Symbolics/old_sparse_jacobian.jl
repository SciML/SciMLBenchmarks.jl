function old_sparsejacobian(ops::AbstractVector, vars::AbstractVector)
    sp = Symbolics.jacobian_sparsity(ops, vars)
    I,J,_ = findnz(sp)

    exprs = old_sparsejacobian_vals(ops, vars, I, J)

    sparse(I, J, exprs, length(ops), length(vars))
end

function old_sparsejacobian_vals(ops::AbstractVector, vars::AbstractVector, I::AbstractVector, J::AbstractVector; simplify::Bool=false, kwargs...)
    exprs = Num[]
    sizehint!(exprs, length(I))

    for (i,j) in zip(I, J)
        push!(exprs, Num(old_expand_derivatives(Differential(vars[j])(ops[i]), simplify; kwargs...)))
    end
    exprs
end


function old_expand_derivatives(O::SymbolicUtils.Symbolic, simplify=false; throw_no_derivative=false)
    if iscall(O) && isa(operation(O), Differential)
        arg = only(arguments(O))
        arg = old_expand_derivatives(arg, false; throw_no_derivative)
        return old_executediff(operation(O), arg, simplify; throw_no_derivative)
    elseif iscall(O) && isa(operation(O), Integral)
        return operation(O)(old_expand_derivatives(arguments(O)[1]; throw_no_derivative))
    elseif !Symbolics.hasderiv(O)
        return O
    else
        args = map(a->old_expand_derivatives(a, false; throw_no_derivative), arguments(O))
        O1 = operation(O)(args...)
        return simplify ? SymbolicUtils.simplify(O1) : O1
    end
end
function old_expand_derivatives(n::Num, simplify=false; kwargs...)
    Symbolics.wrap(old_expand_derivatives(Symbolics.value(n), simplify; kwargs...))
end

function old_occursin_info(x, expr, fail = true)
    if SymbolicUtils.symtype(expr) <: AbstractArray
        if fail
            error("Differentiation with array expressions is not yet supported")
        else
            return occursin(x, expr)
        end
    end

    # Allow scalarized expressions
    function is_scalar_indexed(ex)
        (iscall(ex) && operation(ex) == getindex && !(SymbolicUtils.symtype(ex) <: AbstractArray)) ||
        (iscall(ex) && (SymbolicUtils.issym(operation(ex)) || iscall(operation(ex))) &&
         is_scalar_indexed(operation(ex)))
    end

    # x[1] == x[1] but not x[2]
    if is_scalar_indexed(x) && is_scalar_indexed(expr) &&
        isequal(first(arguments(x)), first(arguments(expr)))
        return isequal(operation(x), operation(expr)) &&
               isequal(arguments(x), arguments(expr))
    end

    if is_scalar_indexed(x) && is_scalar_indexed(expr) &&
        !occursin(first(arguments(x)), first(arguments(expr)))
        return false
    end

    if is_scalar_indexed(expr) && !is_scalar_indexed(x) && !occursin(x, expr)
        return false
    end

    !iscall(expr) && return isequal(x, expr)
    if isequal(x, expr)
        true
    else
        args = map(a->old_occursin_info(x, a, operation(expr) !== getindex), arguments(expr))
        if all(_isfalse, args)
            return false
        end
        Term{Real}(true, args)
    end
end

function old_occursin_info(x, expr::Sym, fail)
    if SymbolicUtils.symtype(expr) <: AbstractArray && fail
            error("Differentiation of expressions involving arrays and array variables is not yet supported.")
    end
    isequal(x, expr)
end

_isfalse(occ::Bool) = occ === false
_isfalse(occ::SymbolicUtils.Symbolic) = iscall(occ) && _isfalse(operation(occ))

_iszero(x) = false
_isone(x) = false
_iszero(x::Number) = iszero(x)
_isone(x::Number) = isone(x)
_iszero(::SymbolicUtils.Symbolic) = false
_isone(::SymbolicUtils.Symbolic) = false
_iszero(x::Num) = _iszero(value(x))::Bool
_isone(x::Num) = _isone(value(x))::Bool


function old_executediff(D, arg, simplify=false; occurrences=nothing, throw_no_derivative=false)
    if occurrences == nothing
        occurrences = old_occursin_info(D.x, arg)
    end

    _isfalse(occurrences) && return 0
    occurrences isa Bool && return 1 # means it's a `true`

    if !iscall(arg)
        return D(arg) # Cannot expand
    elseif (op = operation(arg); SymbolicUtils.issym(op))
        inner_args = arguments(arg)
        if any(isequal(D.x), inner_args)
            return D(arg) # base case if any argument is directly equal to the i.v.
        else
            return sum(inner_args, init=0) do a
                return old_executediff(Differential(a), arg; throw_no_derivative) *
                old_executediff(D, a; throw_no_derivative)
            end
        end
    elseif op === getindex
        inner_args = arguments(arguments(arg)[1])
        c = 0
        for a in inner_args
            if isequal(a, D.x)
                return D(arg)
            else
                c += Differential(a)(arg) * D(a)
            end
        end
        return old_expand_derivatives(c)
    elseif op === ifelse
        args = arguments(arg)
        O = op(args[1], 
            old_executediff(D, args[2], simplify; occurrences=arguments(occurrences)[2], throw_no_derivative), 
            old_executediff(D, args[3], simplify; occurrences=arguments(occurrences)[3], throw_no_derivative))
        return O
    elseif isa(op, Differential)
        # The recursive expand_derivatives was not able to remove
        # a nested Differential. We can attempt to differentiate the
        # inner expression wrt to the outer iv. And leave the
        # unexpandable Differential outside.
        if isequal(op.x, D.x)
            return D(arg)
        else
            inner = old_executediff(D, arguments(arg)[1], false; throw_no_derivative)
            # if the inner expression is not expandable either, return
            if iscall(inner) && operation(inner) isa Differential
                return D(arg)
            else
                # otherwise give the nested Differential another try
                return old_executediff(op, inner, simplify; throw_no_derivative)
            end
        end
    elseif isa(op, Integral)
        if isa(op.domain.domain, Symbolics.AbstractInterval)
            domain = op.domain.domain
            a, b = Symbolics.DomainSets.endpoints(domain)
            c = 0
            inner_function = arguments(arg)[1]
            if iscall(value(a))
                t1 = SymbolicUtils.substitute(inner_function, Dict(op.domain.variables => value(a)))
                t2 = D(a)
                c -= t1*t2
            end
            if iscall(value(b))
                t1 = SymbolicUtils.substitute(inner_function, Dict(op.domain.variables => value(b)))
                t2 = D(b)
                c += t1*t2
            end
            inner = old_executediff(D, arguments(arg)[1]; throw_no_derivative)
            c += op(inner)
            return Symbolics.value(c)
        end
    end

    inner_args = arguments(arg)
    l = length(inner_args)
    exprs = []
    c = 0

    for i in 1:l
        t2 = old_executediff(D, inner_args[i],false; occurrences=arguments(occurrences)[i], throw_no_derivative)

        x = if _iszero(t2)
            t2
        elseif _isone(t2)
            d = Symbolics.derivative_idx(arg, i)
            if d isa Symbolics.NoDeriv
                throw_no_derivative && error((arg, i))
                D(arg)
            else
                d
            end
        else
            t1 = Symbolics.derivative_idx(arg, i)
            t1 = if t1 isa Symbolics.NoDeriv
                throw_no_derivative && error((arg, i))
                D(arg)
            else
                t1
            end
            t1 * t2
        end

        if _iszero(x)
            continue
        elseif x isa SymbolicUtils.Symbolic
            push!(exprs, x)
        else
            c += x
        end
    end

    if isempty(exprs)
        return c
    elseif length(exprs) == 1
        term = (simplify ? SymbolicUtils.simplify(exprs[1]) : exprs[1])
        return _iszero(c) ? term : c + term
    else
        x = +((!_iszero(c) ? vcat(c, exprs) : exprs)...)
        return simplify ? SymbolicUtils.simplify(x) : x
    end
end
