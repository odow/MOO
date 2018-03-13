#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using DynamicProgramming

m = SDPModel(
    stages   = params["stages"],
    sense    = :Max
            ) do sp, t

    cow        = ECow.Cow()

    function update_state(
            x_bcs::Float64,
            x_maintenance::Float64,
            x_lactation::Float64,
            u_supplement::Float64,
            u_herbage::Float64,
            stage::Int
        )
        supplement           = ECow.FoodOffer(ECow.PKE(), u_supplement)
        herbage              = ECow.FoodOffer(ECow.Pasture(), u_herbage)
        total_herbage_intake = 0.0
        total_milk_produced  = 0.0
        for day in 1:7
            (b, m, h, s, mlk) = ECow.updateday!(
                cow,
                7 * (week - 1) + day,
                x_bcs,
                x_maintenance,
                x_lactation > 0.5,
                supplement,
                herbage
            )
            total_herbage_intake += h
            total_milk_produced  += milk
            x_bcs = b
            x_maintenance = m
        end
        x_bcs, x_maintenance, total_herbage_intake, total_milk_produced
    end

    @addstates!(sp, begin
        # the body condition state (US BCS Units)
        bcs          = params["bcs_discretization"]
        # the maintenance energy state (MJ/day)
        maintenance  = params["maintenance_discretization"]
        # lactation status, 0=dry, 1=milking
        lactation    = [0, 1]
        # average pasture cover kg/Ha
        pasturecover = params["pasturecover_discretization"]
    end)

    @addcontrols!(sp, begin
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
            x[pasturecover] + params["grass_growth"][t] -
                params["stocking_rate"] * pasture_eaten
        )

        # return the stage-objective
        return params["stocking_rate"] * (
                params["milk_price"] * milk_produced -
                params["supplement_price"] * params["daysperstage"] * u[supplement]
            ) -
            params["min_penalty"] * max(0.0, params["min_cover"] - y[pasturecover])
    end)

    terminalobjective!(sp, (x) -> begin
        # Set the objective of the final stage
        bcs_penalty = params["stocking_rate"] * 1000.0 * max(
            0.0,
            params["initial_bcs"] - x[bcs]
        )
        cover_penalty = params["min_penalty"] * max(
            0.0,
            params["initial_cover"] - x[pasturecover]
        )
        return -bcs_penalty - cover_penalty
    end)

    constraints!(sp, (x, u, w) -> begin
        # dry-off decision is permanent
        x[lactation] + 0.5 >= u[newlactation]
    end)
end

solve(m)

srand(1111)

# Simulate the policy
@time results = simulate(m,
    1,
    bcs          = 3.2,
    maintenance  = 70.0,
    lactation    = 1,
    pasturecover = 2500.0,
    milkprice    = 6.0
)

ECow.saveobject(resultfile, results)
