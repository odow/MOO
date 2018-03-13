#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
This file can be used as follows:

generate:
    julia MOO.jl generate [parameterfile]

    Example
        julia MOO.jl generate moo.params

build:
    julia MOO.jl build [parameterfile] [model output]

    Example
        julia MOO.jl build moo.params moo.model

solve:
    julia -p N MOO.jl solve [model output] [solved output]

    Example
        julia -p 6 MOO.jl solve moo.model moo.solved.model

simulate:
    julia MOO.jl simulate solved.model results bcs maintenance lactation pasturecover

    Example
        julia MOO.jl simulate moo.solved.model moo.results 3.2 70 1 2500
=#

using DynamicProgramming, JSON

@everywhere include("ecow.jl")

function buildMOO(params)
    SDPModel(
        stages   = params["stages"],
        sense    = :Max
                ) do sp, t
        cow = ECow.Cow()
        function update_state(
                x_bcs::Float64,
                x_maintenance::Float64,
                x_lactation::Float64,
                u_supplement::Float64,
                u_herbage::Float64,
                stage::Int
            )
            feed_supplement      = ECow.FoodOffer(ECow.PKE(), u_supplement)
            feed_herbage         = ECow.FoodOffer(ECow.Pasture(), u_herbage)
            total_herbage_intake = 0.0
            total_milk_produced  = 0.0
            for day in 1:7
                (b, m, h, s, milk) = ECow.updateday!(
                    cow,
                    7 * (stage - 1) + day,
                    x_bcs,
                    x_maintenance,
                    x_lactation > 0.5,
                    feed_supplement,
                    feed_herbage
                )
                total_herbage_intake += h
                total_milk_produced  += milk
                x_bcs                 = b
                x_maintenance         = m
            end
            x_bcs, x_maintenance, total_herbage_intake, total_milk_produced
        end

        function FEIpenalty(x::Float64)
            if -Inf < x <= 3
                return 0.00
            elseif 3 < x <= 4
                return 0.00 + 0.25 * (x - 3)
            elseif 4 < x <= 5
                return 0.25 + 0.50 * (x - 4)
            elseif 5 < x
                return 0.75 + 1.00 * (x - 5)
            end
        end

        @states!(sp, begin
            # the body condition state (US BCS Units)
            bcs          = params["bcs_discretization"]
            # the maintenance energy state (MJ/day)
            maintenance  = params["maintenance_discretization"]
            # lactation status, 0=dry, 1=milking
            lactation    = [0, 1]
            # average pasture cover kg/Ha
            pasturecover = params["pasturecover_discretization"]
            # stocking rate
            stockingrate = 2.5:0.5:4.5
            # feed in storage
            storage      = 0:200:10_000
        end)

        @controls!(sp, begin
            # Quantity of supplement to feed (kg/cow/day)
            supplement   = params["supplement_discretization"]
            # Quantity of pasture to offer (kg/cow/day)
            herbage      = params["herbage_discretization"]
            # lactation status, 0=dry, 1=milking
            newlactation = [0, 1]
        end)

        dynamics!(sp, (y, x, u, w) -> begin
            (y_bcs, y_maintenance, pasture_eaten, milk_produced) = update_state(
                x[bcs],
                x[maintenance],
                x[lactation],
                u[supplement],
                u[herbage],
                t
            )

            # update states
            y[bcs]          = y_bcs
            y[maintenance]  = y_maintenance
            y[lactation]    = u[newlactation]
            y[pasturecover] = min(
                params["max_cover"],
                x[pasturecover] +
                    params["daysperstage"] * params["grass_growth"][t] -
                    x[stockingrate] * pasture_eaten
            )
            y[stockingrate] = x[stockingrate]
            y[storage]      = x[storage] - u[supplement] * params["daysperstage"] * x[stockingrate]

            # return the stage-objective
            return x[stockingrate] * (
                    params["milk_price"] * milk_produced -
                    # params["supplement_price"] * params["daysperstage"] * u[supplement] -
                    params["daysperstage"] * FEIpenalty(u[supplement])
                ) -
                params["min_penalty"] * max(
                    0.0,
                    params["min_cover"] - y[pasturecover]
                )
        end)

        terminalobjective!(sp, (x) -> begin
            # Set the objective of the final stage
            bcs_penalty = x[stockingrate] * params["bcs_penalty"] * max(
                0.0,
                params["initial_bcs"] - x[bcs]
            )
            cover_penalty = params["min_penalty"] * max(
                0.0,
                params["initial_cover"] - x[pasturecover]
            )
            return -bcs_penalty - cover_penalty
        end)

        constraints!(sp, (x, u, w) -> all([
                # dry-off decision is permanent
                x[lactation] + 0.5 >= u[newlactation],
                # maximum lactation length
                t < 40 || u[lactation] < 0.5,
                # cannot feed supplement we don't have
                u[supplement] * params["daysperstage"] * x[stockingrate] <= x[storage]
            ])
        )

    end
end


if length(ARGS) > 0
    if ARGS[1] == "generate"
        papamoa = 0.8 * [
            50.0, 55.0, 45.0, 41.0, 31.0, 19.0, 19.0, 30.0, 47.0, 74.0, 63.0, 50.0
        ]
        params = Dict(
            # =========== stages ===========
            "stages"                      => 52,
            "daysperstage"                => 7,
            # =========== states ===========
            "bcs_discretization"          => 2.5:0.1:4.0,
            "maintenance_discretization"  => 40.0:5:90.0,
            "pasturecover_discretization" => 1500.0:200.0:3500.0,
            # ========== controls ==========
            "supplement_discretization"   => 0:1.0:6,
            "herbage_discretization"      => 0:5.0:50 ,
            # =========== pasture ==========
            # growth per stage (kg/Ha/Week)
            "grass_growth"                => ECow.monthlytoweekly(
                                                papamoa, Dates.Date(2017,8,1)),
            "max_cover"                   => 3500.0,
            "initial_cover"               => 2500.0,
            "min_cover"                   => 1500.0,
            "min_penalty"                 => 2.0,
            # ============= PKE ============
            "supplement_price"            => 0.5,
            # ============ cows ============
            "stocking_rate"               => 3.0,
            "initial_bcs"                 => 3.2,
            "bcs_penalty"                 => 1000.0,
            # ============ milk ============
            "milk_price"                  => 6.0
        )
        open(ARGS[2], "w") do io
            write(io, JSON.json(params))
        end
    elseif ARGS[1] == "build"
        @assert length(ARGS) == 3
        params = JSON.parsefile(ARGS[2])
        m = buildMOO(params)
        open(ARGS[3], "w") do io
            serialize(io, m)
        end
    elseif ARGS[1] == "solve"
        @assert length(ARGS) == 3
        # load model
        m = open(ARGS[2], "r") do io
            deserialize(io)
        end
        # solve model
        solve(m)
        # save model
        open(ARGS[3], "w") do io
            serialize(io, m)
        end
    elseif ARGS[1] == "simulate"
        @assert length(ARGS) == 7
        # load model
        m = open(ARGS[2], "r") do io
            deserialize(io)
        end
        # Simulate the policy
        results = simulate(m,
            1,
            bcs          = parse(Float64, ARGS[4]),
            maintenance  = parse(Float64, ARGS[5]),
            lactation    = parse(Float64, ARGS[6]),
            pasturecover = parse(Float64, ARGS[7])
        )
        # open(ARGS[3], "w") do io
        #     write(io, JSON.json(results))
        # end

        summary = ECow.summarize(
            parse(Float64, ARGS[4]),
            parse(Float64, ARGS[5]),
            results[:supplement],
            results[:herbage],
            results[:lactation]
        )
        open(ARGS[3], "w") do io
            write(io, JSON.json(summary))
        end
    else
        error("Unknown args $(ARGS)")
    end
end
