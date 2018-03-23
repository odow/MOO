julia MOO.jl generate moo.params
julia MOO.jl build moo.params moo.model
julia -p 15 MOO.jl solve moo.model moo.solved.model
julia MOO.jl simulate moo.solved.model moo.results 3.2 70.0 1.0 2500.0
pause
