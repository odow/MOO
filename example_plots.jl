using ECow

# Initialise our cow
cow = ECow.Cow()
# set our herbage
herbage = ECow.FoodOffer(ECow.Pasture(), 50.0)
# Set the supplement
supplement = ECow.FoodOffer(ECow.PKE(), 1.0)

milkingplan(t,n=365) = vcat(trues(t), falses(n-t))

results = Dict(
    :bcs         => vcat(3.2, zeros(365)),
    :maintenance => vcat(63, zeros(365)),
    :lactation   => milkingplan(40 * 7),
    :supplement  => zeros(365),
    :herbage     => zeros(365),
    :milksolids  => zeros(365)
)

for day in 1:365
    (new_bcs, new_maintenance, actual_herbage_intake, actual_supplement_intake, milksolids) = ECow.updateday!(cow, day, results[:bcs][day], results[:maintenance][day], results[:lactation][day], supplement, herbage)
    results[:bcs][day+1]         = new_bcs
    results[:maintenance][day+1] = new_maintenance
    results[:supplement][day]  = actual_supplement_intake
    results[:herbage][day]     = actual_herbage_intake
    results[:milksolids][day]  = milksolids
end

using Plots
const mm = Plots.mm
const pt = Plots.pt
fntsm = Plots.font("times", 12)
fntlg = Plots.font("times", 14)
default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm,left_margin=10mm,bottom_margin=7.5mm)
default(size=(800,600),top_margin=0mm, right_margin=0mm)
gr()

plot(results[:bcs],
xlabel="Days of the Season",
ylabel="BCS (US Scale)",
label="Actual BCS",
ylims=(2, 3.5),
grid=false)
hline!([3.2],
label="End-of-season target BCS",
linestyle=:dash,
width=2,
legend=:bottomright)
savefig("example_body_condition_score.pdf")
quit()
plot(results[:maintenance],
xlabel="Days of the Season",
ylabel="Maintenance (MJ/Day)",
legend=false,
grid=false)
savefig("example_maintenance.pdf")

plot(results[:herbage],
xlabel="Days of the season",
ylabel="Quantity (kgDM/Day)",
grid=false,
ylims=(0, 15),
label="Herbage",
# legend=:bottomright
)
plot!(results[:supplement], label="Supplement")
savefig("example_feed_intake.pdf")

plot(results[:milksolids],
xlabel="Days of the season",
ylabel="Quantity (kg/Day)",
grid=false,
legend=false)
savefig("example_milk_solids.pdf")




println("The End")
