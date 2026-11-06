//! Continuous wall/film/gas balance, independent of the discrete endpoint solve.
use super::fixture::*;
use fs_couple::finite_gas_thermal::wall_film::implicit::WallFilmImplicitControls;
use fs_material::{
    evaporation::gas::wall_film::GasWallFilmInventory, phase::UniformEnthalpyStepInput,
};

pub(super) struct PhysicalBalance {
    pub(super) wall_capacity: f64,
    pub(super) carrier_mass: f64,
    pub(super) total_volatile_mass: f64,
    pub(super) sensible_conductance: f64,
    pub(super) phase_exposure: f64,
}
impl PhysicalBalance {
    // y = [wall temperature, liquid mass, gas temperature]. All inventories
    // remain strictly wet here; dry-out is covered by the separate event test.
    // The fixture's actual original NASA9 species have Cp/R=3.5, while liquid
    // Cp=4200. The original phase law's saturation curve is integrated below.
    fn derivative(&self, y: [f64; 3]) -> [f64; 3] {
        let [tw, film, tg] = y;
        let ra = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028;
        let vapor = self.total_volatile_mass - film;
        assert!(film > 0.0 && vapor > 0.0);
        let pressure = (self.carrier_mass * ra + vapor * R) * tg / GAS_VOLUME;
        let partial = vapor * R * tg / GAS_VOLUME;
        let cp_vapor = 3.5 * R;
        let delta_cp = cp_vapor - CP_LIQUID;
        let latent = LATENT + delta_cp * (tw - 300.0);
        let saturation = 3536.8
            * ((LATENT - delta_cp * 300.0) / R * (1.0 / 300.0 - 1.0 / tw)
                + delta_cp / R * (tw / 300.0).ln())
            .exp();
        let rate = self.phase_exposure * pressure / (R * tg)
            * ((pressure - partial) / (pressure - saturation)).ln();
        let sensible = self.sensible_conductance * (tw - tg);
        let upstream = if rate < 0.0 { tg } else { tw };
        let wall_rate = (-sensible - rate * (latent + cp_vapor * (upstream - tw)))
            / (self.wall_capacity + film * CP_LIQUID);
        // Advected h minus the gas's own u includes pressure-flow energy.
        let gas_rate = (sensible + rate * (cp_vapor * (upstream - tg) + R * tg))
            / (2.5 * (self.carrier_mass * ra + vapor * R));
        [wall_rate, -rate, gas_rate]
    }
    pub(super) fn reference(&self, mut y: [f64; 3], duration: f64) -> [f64; 3] {
        let dt = duration / 4096.0;
        for _ in 0..4096 {
            let a = self.derivative(y);
            let b = self.derivative(std::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
            let c = self.derivative(std::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
            let d = self.derivative(std::array::from_fn(|i| y[i] + dt * c[i]));
            y = std::array::from_fn(|i| y[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
        }
        y
    }
}

#[test]
fn whole_wall_film_and_gas_trajectory_refines_toward_independent_continuous_energy_balance() {
    let source = source();
    let initial_gas = gas(&source, 320.0, 6000.0);
    let binding = binding(&source);
    let curve = curve();
    let body = body(20.0);
    let duration = 0.1;
    let original = input(&curve, 300.0, duration);
    let film_mass = 1e-6;
    let balance = PhysicalBalance {
        wall_capacity: original.mass_kg * CP_WALL,
        carrier_mass: initial_gas.species_mass_kg(&carrier()).unwrap(),
        total_volatile_mass: film_mass + initial_gas.species_mass_kg(&water()).unwrap(),
        sensible_conductance: body.convection_w_per_m2_k * AREA,
        phase_exposure: boundary(&binding).mass_transfer_coefficient_m_s * AREA,
    };
    let reference = balance.reference([300.0, film_mass, 320.0], duration);
    assert!(reference[0] > 300.0 && reference[1] > film_mass && reference[2] < 320.0);
    let mut errors = Vec::new();
    for steps in [1, 2, 4, 8] {
        let mut gas = initial_gas.clone();
        let mut state = original.initial;
        let mut film = GasWallFilmInventory::Wet {
            mass_kg: film_mass,
            specific_enthalpy_j_kg: liquid_h(300.0).unwrap(),
        };
        scope(|cx, _| {
            for _ in 0..steps {
                let step = UniformEnthalpyStepInput {
                    initial: state,
                    duration_s: duration / steps as f64,
                    ..original
                };
                let next = environment(&gas, &source, film, &binding, &body)
                    .advance_implicit(cx, step, WallFilmImplicitControls::default())
                    .unwrap();
                assert!(next.report.phase.transferred_mass_kg < 0.0);
                assert!(next.report.phase.dry_out_time_s.is_none());
                state = next.state;
                film = next.report.next_film;
                gas = next.report.next_gas;
            }
        });
        let actual = [
            state.temperature_k(),
            film.mass_kg(),
            gas.state().state().temperature,
        ];
        let error = (actual[0] - reference[0])
            .abs()
            .max((actual[1] - reference[1]).abs() * LATENT / balance.wall_capacity)
            .max((actual[2] - reference[2]).abs());
        assert!(error.is_finite() && error > 0.0);
        errors.push(error);
        near(
            gas.mass_kg() + film.mass_kg(),
            initial_gas.mass_kg() + film_mass,
            1e-13,
        );
    }
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.8 * pair[0],
            "coupled time refinement stalled: {errors:?}"
        );
    }
    assert!(
        errors[3] < 0.2 * errors[0],
        "first-order endpoint ladder: {errors:?}"
    );
    eprintln!("independent continuous endpoint={reference:?}; coupled time errors={errors:?}");
}
