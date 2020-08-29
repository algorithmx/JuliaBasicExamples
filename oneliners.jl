## 1d numerical integration ===================================

trapezoidal(f_x, h) = (0.5*h) * ( f_x[1] 
                                + 2*sum(f_x[2:(length(f_x)-1)]) 
                                + f_x[end] )


simpson(x,h) = (1//3) * h * ( x[1]
                            + 4sum(x[2:2:length(x)-1])
                            + 2sum(x[3:2:length(x)-2])
                            + x[end] )

## examples

η = 0.01π

x = η.*collect(0:100)

# 22.140692632779267`
trapezoidal0(exp.(x),η)

trapezoidal(exp.(x),η)

simpson(exp.(x),η)

# 0
trapezoidal(cos.(x),η)
simpson(cos.(x),η)

# 2
trapezoidal(sin.(x),η)
simpson(sin.(x),η)

                            
## compute π =============================================
using LinearAlgebra

compute_pi(N) = (4/N)*sum((rand(N).^2 .+ rand(N).^2).<1.0)

##
compute_pi(100_000_000)


## brackets matching

is_bracket_match(s) = sum((s[i]=='(' ? 1 : (s[i]==')' ? -1 : 0)) for i=1:length(s))==0

## test 

s1 = "(((a+b)/c)-d)"
is_bracket_match(s1)
##
s0 = "(((a+b)/c)-d))"
is_bracket_match(s0)


