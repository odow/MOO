#  Copyright 2018, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

module ECow

"""
A type for holding Cow related parameters
"""
struct Cow
    age::Int                                  # Age of cow at start of season (years)
    calfbirthweight::Float64                  # Birth weight of calf (kg)
    initialbcs::Float64                       # BCS at calving (US scale)
    initialliveweight::Float64                # Liveweight at calving (kg)
    conceptiondate::Int                       # Date of conception (days from calving)
    gestationlength::Int                      # Length of pregnancy (days)
    seasonlength::Int                         # Number of days before dried off (days)
    matureliveweight::Float64                 # Liveweight of mature animal (kg)
    wilminkprotein::Tuple{Vararg{Float64, 3}} # Wilmink protein coefficients
    wilminkfat::Tuple{Vararg{Float64, 3}}     # Wilmink fat coefficients
    vetheraniamA::Float64                     # Vetheraniam milk coefficients
    vetheraniamb::Tuple{Vararg{Float64, 3}}   # Vetheraniam milk coefficients
    vetheraniamc::Tuple{Vararg{Float64, 3}}   # Vetheraniam milk coefficients
    bcsperliveweight::Float64                 # The BCS per LW convervsion
    substitutionrate::Float64                 # Rate at which pasture is substituted for supplement

    max_milk_energy::Vector{Float64}            #
    base_lipid_change::Vector{Float64}          #   some storage containers
    base_liveweight_change::Vector{Float64}     #   for precomputed values
    milk_conversion_efficiency::Vector{Float64} #
    energy_pregnancy::Vector{Float64}
    liveweight_growth::Vector{Float64}
end
function Cow(initialbcs, intialliveweight, seasonlength)
    cow = Cow(
        4,                                   # Age of cow at start of season (years)
        40.0,                                # Birth weight of calf (kg)
        initialbcs,                          # BCS at calving (US scale)
        intialliveweight,                    # Liveweight at calving (kg)
        90,                                  # Date of conception (days from calving)
        284,                                 # Length of pregnancy (days)
        seasonlength,                        # Number of days before dried off (days)
        550.0,                               # Liveweight of mature animal (kg)
        (3.072, 1.0, 0.00337),               # Wilmink protein coefficients
        (3.685, 2.62, 0.0049),               # Wilmink fat coefficients
        4.029e-9,                            # Vetheraniam A
        (-9.614e9, 2.642e10, 1.093e9),       #
        (-0.116, -9.907e-4, -5.943),         #
        (5 + initialbcs) / intialliveweight, #
        0.21,                                # Substitution rate Stockdale (2000)
        zeros(365),
        zeros(365),
        zeros(365),
        zeros(365),
        zeros(365),
        zeros(365)
        )
    for day=1:365
        cow.max_milk_energy[day]            = energy(Milk, cow, day, true)
        cow.base_lipid_change[day]          = changeliveweight(Lipid, cow, day)
        cow.milk_conversion_efficiency[day] = kgmsperME(1.0, cow, day)
        cow.energy_pregnancy[day]           = energy(Pregnancy, cow, day)
        cow.liveweight_growth[day]          = changeliveweight(Growth, cow, day)
    end
    cow.base_liveweight_change[1] = changeliveweight(Pregnancy, cow, 1) + cow.base_lipid_change[1] + cow.liveweight_growth[1]
    for day=2:365
        cow.base_liveweight_change[day] = cow.base_liveweight_change[day-1] + changeliveweight(Pregnancy, cow, day) + cow.base_lipid_change[day] + cow.liveweight_growth[day]
    end
    return cow
end
Cow() = Cow(3.2, 460.0, 280)

"""
A type for holding Food related parameters
"""
struct Feed
    digestibility::Float64
    metabolisableenergy::Float64
    neutraldetergentfibre::Float64
    terrainfactor::Float64
    greenforage::Float64
end
# Some common values
Pasture() = Feed(0.85, 10.3, 0.44, 1.0, 1.5)
PKE() = Feed(0.8, 10.3, 0.44, 1.0, 1.5)

struct FoodOffer
    food::Feed
    quantity::Float64
end

"""
From here on down are internal helper methods
"""

const MEDIET = 10.3 # Kludgy hack for now
# Efficiency of ME utilisation for milk synthesis
kl(cow::Cow) = 0.02 * MEDIET + 0.4
# Efficiency of ME utilisation for maintenance
km(cow::Cow) = 0.02 * MEDIET + 0.5
# Efficiency of ME utilisation for fat deposition
function kgain(cow::Cow, islactating::Bool)
    if islactating
        return 0.02 * MEDIET + 0.4
    else
        return 0.042 * MEDIET + 0.006
    end
end

# Adjust the milking potential of cow by age
const age_adj_vals = [0.0, 0.0, 0.75, 0.87, 0.95, 1.0, 1.0, 1.0, 1.0, 0.97, 0.92]
"""
Adjust the milk curve for the cow based on age at end of season
"""
age_adjustment(cow::Cow) = age_adj_vals[cow.age]
"""
Energy content per 1kg of liquid milk
"""
MEperlitre(cow::Cow, t::Int) = (0.376 * fat(cow, t) + 0.209 * protein(cow, t) + 0.948) / kl(cow)
"""
Solids percentage of 1kg of liquid milk
"""
solidspercentage(cow::Cow, t::Int) = (fat(cow, t) + protein(cow, t)) / 100
"""
kgMS produced per MJ of ME set aside for milk production
"""
kgmsperME(ME::Float64, cow::Cow, t::Int) = 1.23 * ME / MEperlitre(cow, t) * solidspercentage(cow, t)
"""
Wilmink fat curve for cow
    % Fat composition (0-100%)
"""
fat(cow::Cow, t::Int) = cow.wilminkfat[1] + cow.wilminkfat[2] * exp(-0.05 * t) + cow.wilminkfat[3] * t
"""
Wilmink protein curve for cow
    % Protein composition (0-100%)
"""
protein(cow::Cow, t::Int) = cow.wilminkprotein[1] + cow.wilminkprotein[2] * exp(-0.05 * t) + cow.wilminkprotein[3] * t

"""
Convert US BCS' (0-5) to NZ (0-10)
"""
NZBCS(USBCS::Float64) = (USBCS - 1.5) / 0.32
"""
Convert NZ BCS' (0-10) to US (0-5)
"""
USBCS(NZBCS::Float64) = 0.32 * NZBCS + 1.5

const Milk = Val{:milk}
"""
Maximum energy available for milk production
output: MJ/day
"""
function energy(::Type{Milk}, cow::Cow, day::Int, islactating::Bool)
    if islactating
        y = 0.
        for i = 1:3
            y += cow.vetheraniamb[i] * exp(day * cow.vetheraniamc[i])
        end
        return  y * age_adjustment(cow) * cow.vetheraniamA / kl(cow)
    else
        return 0.0
    end
end

const Growth = Val{:growth}
"""
The energy required for growth
output: MJ/day
"""
energy(::Type{Growth}, cow::Cow, t::Int, liveweight) = 1.3 * (4.1 + 0.0332 * liveweight - 9e-6 * liveweight^2) / (1 - 0.1475 * cow.liveweight_growth[t]) * cow.liveweight_growth[t]

const Maintenance = Val{:maintenance}
"""
The energy required for maintenance
output: MJ/day
"""
function energy(::Type{Maintenance}, cow::Cow, pasture::Feed, liveweight, actualherbageintake, actualmilkproduction)
    b0 = 1.4 * 0.28 * exp(-0.03 * cow.age) / km(cow)
    b1 = 0.006 * (0.9 - pasture.digestibility) / km(cow)
    # b1 = 0.0025 * (0.9 - feed.D) / km(cow)
    b2 = 0.05 * pasture.terrainfactor / ((3 + pasture.greenforage) * km(cow))
    return b0 * liveweight^0.75 + b1 * actualherbageintake * liveweight + b2 * liveweight + 0.1 * actualmilkproduction
end
const ageliveweightadjustment = [0.0, 0.0, 0.85, 0.92, 0.96, 1.0, 1.0, 1.0]
function changeliveweight(::Type{Growth}, cow::Cow, t::Int)
    return (ageliveweightadjustment[cow.age]* cow.matureliveweight - cow.initialliveweight) / 365
end

const Pregnancy = Val{:pregnancy}
"""
The energy required for pregnancy
output: MJ/day
"""
function energy(::Type{Pregnancy}, cow::Cow, t::Int)
    dp = t - cow.conceptiondate
    if dp < 1
        return 0.
    else
        a = 0.025 * cow.calfbirthweight^2
        b = 10^(151.665 - 151.64 * exp(-5.76e-5 * dp))
        c = 0.0201 * exp(-5.76e-5 * dp)
        return a * b * c / (0.13 * 40)
    end
end
function changeliveweight(::Type{Pregnancy}, cow::Cow, t::Int)
    dp = t - cow.conceptiondate
    if dp < 1
        return 0.
    else
        return cow.calfbirthweight * (18.28 * 0.02 - 0.0000286 * dp) * exp(0.02 * dp - 0.0000143 * dp^2) / 1000
    end
end

function MEperKGLWChange(cow::Cow, bcs::Float64, is_gain::Bool, islactating::Bool)
    # ME per 1kg lw change (MEMJ / KgLW)
    if is_gain
        # More efficient at using energy to gain weight
        return (10.1 + 1.976 * bcs) / kgain(cow, islactating)
    else
        # Less efficient at burning lipid to energy
        return (10.1 + 1.976 * bcs) / kgain(cow, islactating) * 0.84
    end
end

function totalrequiredenergy(cow::Cow, day::Int, maintenance, liveweight, islactating, bcs, herbage::Feed, energy_lipid, energy_milk)
    energy(Growth, cow, day, liveweight) +
    cow.energy_pregnancy[day] +
    energy_milk +
    maintenance +
    energy_lipid
end
function getenergy(::Type{Milk}, cow, day, islactating)
    if islactating
        return cow.max_milk_energy[day]
    else
        return 0.
    end
end

"""
Stage of lactation
add +3 to shift the daily values into the center of the week.
"""
SOL(t::Int) = 0.67 + 0.0972 * (4.0401 * log10((t+3) / 7) - 0.095 * ((t+3) / 7) + 0.095)
"""
Potential Dry Matter Intake is limited by grazing ability, energy requirements and rumen fill
"""
function PotentialDMI(liveweight::Float64, herbage::Feed, requiredenergy::Float64, day::Int)
    return min(DMIgrazing(liveweight, day), DMIenergy(requiredenergy, herbage), DMIrumen(liveweight, herbage, day))
end
"""
Potential Dry Matter Intake due to grazing
"""
DMIgrazing(liveweight::Float64, day::Int) = 0.0375 * liveweight * SOL(day)
"""
Potential Dry Matter Intake due to energy requirements
"""
DMIenergy(requiredenergy::Float64, herbage::Feed) = requiredenergy / herbage.metabolisableenergy
"""
Potential Dry Matter Intake due to rumen size
"""
DMIrumen(liveweight::Float64, herbage::Feed, day::Int) = 0.0165 * liveweight * SOL(day) / herbage.neutraldetergentfibre


"""
Harvesting Efficiency
This is for ryegrass to ground level
fraction [0, 1]
"""
harvestefficiency(pdmi::Float64, q::Float64) = 0.57676 * (pdmi / q)^0.536
# Ryegrass to a 4cm cutting height
# return 83.33 * (potential_DMI / available_pasture)^-0.7
# Lucerne to a 4cm cutting height
# return -0.322 * log(potential_DMI / available_pasture) + 0.7128

"""
Actual herbage intake
kgDM
"""
function herbageintake(cow::Cow, day::Int, liveweight::Float64, herbage::Feed, herbageoffer::Float64, requiredenergy::Float64, supplementation::Float64)
    # Potential Dry Matter Intake
    pdmi = PotentialDMI(liveweight, herbage, requiredenergy, day)
    # Herbage intake without supplement
    no_supplement = harvestefficiency(pdmi, herbageoffer) * herbageoffer
    # Substitution Rate -- N.B. economic answers are v. sensitive to this!
    SR = cow.substitutionrate * no_supplement / liveweight * 100 - 0.18 # Stockdale 2000
    # Actual pasture intake
    return   no_supplement - SR * supplementation
end

# Kg change in body lipid by days in milk
const Lipid = Val{:lipid}
function changeliveweight(::Type{Lipid}, cow::Cow, t::Int)
    # Daily change in lipid mass (kg)

    max_lipid_loss = -1.75
    t_prime = 112.

    # Empty weight of cow
    # empty_weight = cow.constants.LW0 * (1 - 0.1113 / (1 - 0.129 * cow.constants.BCS0))
    BCSstandard = 5
    LWatBCSstandard = cow.initialliveweight / (1 - 0.129 * (BCSstandard - cow.initialbcs))
    kgLWperBCS = 0.129 * LWatBCSstandard
    LWatBCS3 = LWatBCSstandard - 2 * kgLWperBCS
    empty_weight = LWatBCS3 * (1 - 0.15) # the 15% gut fill

    # Mass of body lipid at calving (kg)
    calving = 0.12 * (cow.initialbcs - 0.36) * empty_weight

    # Mass of body lipid at T' (kg)
    prime = 0.26 * empty_weight
    # prime = 0.2988 * empty_weight

    # Desired mass of body lipid at calving (empty kg)
    next_calving = 0.12 * (3.2 - 0.36) * empty_weight
    # next_calving = 0.12 * (3.3 - 0.36) * empty_weight

    # Days from conception
    DFcon = t - cow.conceptiondate

    # Case where next pregnancy occurs after T'
    if t_prime < cow.conceptiondate
        # Rate of change of lipid mass at calving
        change_at_calving = 2 * (prime - calving) / t_prime

        # There is some maximum rate of lipid loss
        if change_at_calving < max_lipid_loss
            # In which case make it the maximum
            change_at_calving = max_lipid_loss

            # And adjust T'
            t_prime = 2 * (prime - calving) / max_lipid_loss
        end
    else
        # Days from conception
        DFcon = t_prime - cow.conceptiondate

        # Cow is aiming for some target at conception
        target = prime + (DFcon)^2 * (next_calving - prime) / cow.gestationlength^2

        dlt = 2 * DFcon * (next_calving - prime) / cow.gestationlength^2

        # Rate of change of lipid mass at calving
        change_at_calving = 2 * (target - calving) / t_prime - dlt

        # There is some maximum rate of lipid loss
        if change_at_calving < max_lipid_loss
            # In which case make it the maximum
            change_at_calving = max_lipid_loss
            warn("Max lipid loss occured")
            # And adjust T'
            t_prime = 2 * (prime - calving) / (max_lipid_loss + dlt)
        end
    end
    if t <= t_prime
        # linear change to zero lipid change at T'
        return change_at_calving * (1 - t / t_prime)
    elseif t_prime < t && t < cow.conceptiondate
        # Day is after T' but still not pregnant so no drive for lipid change
        return 0.
    else
        # Cow driven to reach next_calving body lipid by next calving
        return 2 * (next_calving - prime) * (t - cow.conceptiondate) / cow.gestationlength^2
    end

end

function updateday!(cow::Cow, day::Int, bcs::Float64, maintenance::Float64, islactating::Bool, supplementation::FoodOffer, herbage::FoodOffer)
    # Deterministic change in liveweight
    liveweight = cow.initialliveweight + cow.base_liveweight_change[day] + (bcs - cow.initialbcs) / cow.bcsperliveweight

    energy_milk = getenergy(Milk, cow, day, islactating)
    energy_lipid = cow.base_lipid_change[day] * MEperKGLWChange(cow, bcs, cow.base_lipid_change[day]>0, islactating)

    energy_required = totalrequiredenergy(cow, day, maintenance, liveweight, islactating, bcs, herbage.food, energy_lipid, energy_milk)
    herbage_intake  = herbageintake(cow, day, liveweight, herbage.food, herbage.quantity, energy_required, supplementation.quantity)
    net_energy   = herbage_intake * herbage.food.metabolisableenergy + supplementation.quantity * supplementation.food.metabolisableenergy - energy_required

    allocated_to_milk = 0.0
    if islactating
        max_milk_me = energy_milk                          # Can work out what energy we 'want' to spend on milk (maximum)
        if net_energy  > 0.0
            # Positive energy balance even with max milk
            allocated_to_milk = 0.0
        elseif net_energy  < -(max_milk_me + abs(energy_lipid))
            # We have such a negative energy balance that we can't afford to make any milk
            allocated_to_milk = -max_milk_me
        else
            # Some reduction in milk production to support lipid deposition
            allocated_to_milk = net_energy  * max_milk_me / (max_milk_me + abs(energy_lipid))
        end
    end
    change_bcs = (energy_lipid + net_energy - allocated_to_milk) * cow.bcsperliveweight / MEperKGLWChange(cow, bcs, energy_lipid>0, islactating)

    maintenance = energy(Maintenance, cow, herbage.food, liveweight, herbage_intake, energy_milk + allocated_to_milk)
    milksolids  = (energy_milk + allocated_to_milk) * cow.milk_conversion_efficiency[day]

    return bcs + change_bcs, maintenance, herbage_intake, supplementation.quantity, milksolids
end

function monthlytoweekly(x::Vector{Float64}, date::Dates.Date)
    @assert length(x) == 12 # check we have monthly data
    weekly_growth   = zeros(52)
    number_readings = zeros(Int, 52)
    for i=1:52
        for j=1:7
            weekly_growth[i] += interpolate(x, date)
            number_readings[i] += 1
            date += Dates.Day(1)
        end
    end
    weekly_growth ./ number_readings
end

function interpolate(x::Vector{Float64}, date::Dates.Date, reference_day=15)
    year, month, day = Dates.year(date), Dates.month(date), Dates.day(date)
    left_year, right_year = year, year
    left_month, right_month = if day < reference_day
        max(1, month-1), month
    else
        month, min(12, month+1)
    end
    if day <= reference_day && month == 1
        left_year  = year - 1
        left_month = 12
    elseif day >= reference_day && month==12
        right_year  = year + 1
        right_month = 1
    end
    left_date  = Dates.Date(left_year, left_month, reference_day)
    right_date = Dates.Date(right_year, right_month, reference_day)

    lambda = (date - left_date).value / (right_date - left_date).value

    return (1 - lambda) * x[Dates.month(left_date)] + lambda * x[Dates.month(right_date)]
end

function summarize(bcs, maintenance, u_supplement, u_herbage, u_lactation)
    cow             = ECow.Cow()
    feed_supplement = [ECow.FoodOffer(ECow.PKE(), y) for y in u_supplement]
    feed_herbage    = [ECow.FoodOffer(ECow.Pasture(), y) for y in u_herbage]
    lactation       = [y > 0.5 for y in u_lactation]

    results = Dict{Symbol, Any}(
        # ===== states =====
        :bcs           => Float64[],
        :maintenance   => Float64[],
        :lactation     => Float64[],
        :cover         => Float64[],
        # ==== controls ====
        :supplement    => Float64[],
        :herbage       => Float64[],
        # ===== profit =====
        :milksolids    => Float64[]
    )

    t = 1
    for week in 1:52
        supplement   = feed_supplement[week]
        herbage      = feed_herbage[week]
        is_lactating = lactation[week]
        for day in 1:7
            bcs, maintenance, h, s, ms = updateday!(
                cow,
                t,
                bcs,
                maintenance,
                is_lactating,
                supplement,
                herbage
            )
            t += 1
            push!(results[:bcs], bcs)
            push!(results[:maintenance], maintenance)
            push!(results[:lactation], is_lactating)
            push!(results[:supplement], s)
            push!(results[:herbage], h)
            push!(results[:milksolids], ms)
        end
    end

    results[:final_bcs]        = results[:bcs][end]
    results[:total_herbage]    = sum(results[:herbage])
    results[:total_supplement] = sum(results[:supplement])
    results[:total_milksolids] = sum(results[:milksolids])
    results[:feed_imported]    = results[:total_supplement] / (results[:total_supplement] + results[:total_herbage])

    results
end

end
