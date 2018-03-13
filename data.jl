using JSON

params = Dict(
    # =========== stages ===========
    "stages" => 52,
    "daysperstage" => 7,

    # =========== states ===========
    "bcs_discretization" => 2.5:0.025:4.0,
    "maintenance_discretization" => 40.0:5:90.0,
    "pasturecover_discretization" => 1500.0:100.0:3500.0,

    # ========== controls ==========
    "supplement_discretization" => 0:0.5:10,
    "herbage_discretization" => 0:5:50 ,

    # =========== pasture ==========
    # growth per stage (kg/Ha/Week)
    "grass_growth"     => [],
    "initial_cover"    => 2500.0,
    "max_cover"        => 3500.0,
    "min_cover"        => 1500.0,
    "min_penalty"      => 100.0,

    # ============= PKE ============
    "supplement_price" => 0.5,

    # ============ cows ============
    "stocking_rate" => 3.0,
    "initial_bcs" => 3.2,
    "bcs_penalty" => 1000.0,

    # ============ milk ============
    "milk_price" => 6.0,
)
