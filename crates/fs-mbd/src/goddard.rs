#[derive(Debug, Clone)]
pub struct GoddardParams {
    pub chamber_pressure_psi: f64,
    pub fuel_flow_kg_per_sec: f64,
    pub throat_area_cm2: f64,
    pub expansion_ratio: f64,
}

#[derive(Debug, Clone)]
pub struct GoddardResult {
    pub chamber_pressure_psi: f64,
    pub chamber_pressure_pa: f64,
    pub exhaust_velocity_mps: f64,
    pub thrust_newtons: f64,
    pub specific_impulse_sec: f64,
    pub mach_exit: f64,
}

pub fn step_goddard_rocket(params: &GoddardParams) -> GoddardResult {
    let chamber_pressure_pa = params.chamber_pressure_psi * 6894.76;
    let gamma = 1.24; // Combustion products heat capacity ratio
    let chamber_temp_k = (2400.0 + (chamber_pressure_pa / 2.4e6) * 400.0).round();
    let gas_constant_r = 365.0; // J/(kg*K) for gasoline + liquid O2
    let expansion = params.expansion_ratio.max(1.4);

    // Supersonic Mach number at exit via area-Mach relation
    let mach_exit =
        ((2.0 / (gamma - 1.0)) * (params.expansion_ratio.powf(2.0 / (gamma + 1.0)) - 1.0)).sqrt();
    let exhaust_velocity_mps = (((2.0 * gamma) / (gamma - 1.0))
        * gas_constant_r
        * chamber_temp_k
        * (1.0 - 1.0 / expansion.powf(gamma - 1.0)))
    .sqrt()
    .round();
    let thrust_newtons = (params.fuel_flow_kg_per_sec * exhaust_velocity_mps).round();
    let specific_impulse_sec = exhaust_velocity_mps / 9.80665;

    GoddardResult {
        chamber_pressure_psi: params.chamber_pressure_psi,
        chamber_pressure_pa,
        exhaust_velocity_mps,
        thrust_newtons,
        specific_impulse_sec: (specific_impulse_sec * 10.0).round() / 10.0,
        mach_exit: (mach_exit * 100.0).round() / 100.0,
    }
}
