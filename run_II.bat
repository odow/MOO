julia MOO_II.jl generate moo.params
julia MOO_II.jl build moo.params moo.model
julia -p 15 MOO_II.jl solve moo.model moo.solved.model
julia MOO_II.jl simulate moo.solved.model moo.results 3.2 70.0 1.0 2500.0
pause
