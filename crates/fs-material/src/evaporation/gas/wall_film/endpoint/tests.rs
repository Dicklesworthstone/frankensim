//! Independent original-NASA9 endpoints for zero-capable wall condensate.
use super::*;
use crate::evaporation::gas::tests::{Phase, binding, catalog, controls, gas, liquid, temperature};
use fs_thermochem::SpeciesId;

const R: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
const RA: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
const LATENT: f64 = 200_000.0;
const AREA: f64 = 0.001;

fn vapor() -> SpeciesId {
    SpeciesId::new("B").unwrap()
}
fn carrier() -> SpeciesId {
    SpeciesId::new("A").unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}
fn h(t: f64) -> f64 {
    liquid(t).specific_enthalpy_j_kg
}
fn vapor_h(t: f64) -> f64 {
    R * (4.0 * t + 0.0005 * t * t)
}
fn dh(a: f64, b: f64) -> f64 {
    R * (b - a) * (4.0 + 0.0005 * (a + b))
}
fn rate(tw: f64, tg: f64, pressure: f64, partial: f64) -> f64 {
    let saturation = 1500.0 * fs_math::det::exp(LATENT / R * (1.0 / 300.0 - 1.0 / tw));
    pressure / (R * tg) * ((pressure - partial) / (pressure - saturation)).ln()
}
fn coefficient(
    initial: &ClosedGasParcel,
    dm: f64,
    requested: f64,
    tg: f64,
    tw: f64,
    dt: f64,
) -> f64 {
    let vapor_mass = initial.species_mass_kg(&vapor()).unwrap() + dm;
    let pressure = (initial.species_mass_kg(&carrier()).unwrap() * RA + vapor_mass * R) * tg
        / initial.volume_m3();
    let partial = vapor_mass * R * tg / initial.volume_m3();
    let coefficient = requested / (AREA * dt * rate(tw, tg, pressure, partial));
    assert!(coefficient > 0.0);
    coefficient
}
// Integrate the original variable Cv, canceling each formation reference
// before subtracting the upstream vapor enthalpy carried with the mass stream.
fn required_gas_heat(initial: &ClosedGasParcel, dm: f64, tg: f64, tw: f64) -> f64 {
    let t0 = initial.state().state().temperature;
    let vapor_mass = initial.species_mass_kg(&vapor()).unwrap() + dm;
    initial.species_mass_kg(&carrier()).unwrap() * RA * 2.5 * (tg - t0)
        + vapor_mass * R * (tg - t0) * (3.0 + 0.0005 * (tg + t0))
        - dm * R * t0
        - dm * dh(t0, if dm < 0.0 { tg } else { tw })
}
fn surface<'a>(
    binding: &'a GasLiquidBinding,
    film: GasWallFilmInventory,
    tw: f64,
    coefficient: f64,
) -> GasWallFilmSurface<'a> {
    GasWallFilmSurface {
        law: &Phase,
        binding,
        film,
        area_m2: AREA,
        interface_liquid_enthalpy_j_kg: h(tw),
        temperature_from_enthalpy: &temperature,
        mass_transfer_coefficient_m_s: coefficient,
    }
}

#[test]
fn dry_condensation_uses_endpoint_gas_affinity_and_upstream_enthalpy_in_one_balance() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial = gas(&source, 320.0, 5000.0, 0.001);
    let snapshot = initial.clone();
    let dm = -0.1 * initial.species_mass_kg(&vapor()).unwrap();
    let (tw, tg, dt) = (300.0, 314.0, 0.5);
    let k = coefficient(&initial, dm, dm, tg, tw, dt);
    let heat = required_gas_heat(&initial, dm, tg, tw);
    let input = surface(&binding, GasWallFilmInventory::Dry, tw, k);
    let endpoint = GasWallFilmEndpoint::new(input, &initial, &source, dt, controls()).unwrap();
    let trial = endpoint.evaluate(h(tw), dm, tg, heat, || false).unwrap();
    endpoint.check_balance(&trial).unwrap();
    near(trial.gas_mass_residual_kg(), 0.0, 1e-14);
    near(trial.gas_energy_residual_j(), 0.0, 1e-8);
    assert_eq!(trial.heat_into_gas_j(), heat);
    let candidate = endpoint.finalize(&trial, || false).unwrap();
    let step = candidate.step();
    near(candidate.gas().state().state().temperature, tg, 2e-8);
    near(step.transferred_mass_kg, dm, 1e-14);
    near(candidate.film().mass_kg(), -dm, 1e-14);
    assert!(step.film_formed);
    assert_eq!(step.film_sensible_energy_j, 0.0);
    assert_eq!(step.dry_out_time_s, None);
    near(
        step.wall_heat_into_system_j,
        dm * (LATENT + dh(tw, tg)),
        1e-8,
    );
    near(step.energy_into_gas_j, dm * vapor_h(tg), 1e-8);
    assert!((step.wall_heat_into_system_j - dm * (LATENT + dh(tw, 320.0))).abs() > 0.03);
    let frozen = AREA * dt * k * rate(tw, 320.0, 100_000.0, 5000.0);
    assert!((step.transferred_mass_kg - frozen).abs() > 1e-7);
    near(
        step.film_enthalpy_change_j + step.gas_internal_energy_change_j
            - step.wall_heat_into_system_j
            - heat,
        0.0,
        1e-8,
    );
    near(
        candidate.gas().mass_kg() + candidate.film().mass_kg(),
        initial.mass_kg(),
        1e-14,
    );
    near(
        candidate.gas().species_mass_kg(&carrier()).unwrap(),
        initial.species_mass_kg(&carrier()).unwrap(),
        1e-14,
    );
    assert_eq!(candidate.gas().point().axes().get("local-lot"), Some(&11.0));
    assert_eq!(initial, snapshot);
}

#[test]
fn complete_dryout_keeps_outgoing_vapor_heat_and_then_has_exactly_zero_film_storage() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial = gas(&source, 290.0, 500.0, 0.001);
    let mass = 2e-7;
    let (tw, tg, dt) = (300.0, 294.0, 1.0);
    let film = GasWallFilmInventory::Wet {
        mass_kg: mass,
        specific_enthalpy_j_kg: h(295.0),
    };
    let k = coefficient(&initial, mass, 3.0 * mass, tg, tw, dt);
    let input = surface(&binding, film, 295.0, k);
    let endpoint = GasWallFilmEndpoint::new(input, &initial, &source, dt, controls()).unwrap();
    let heat = required_gas_heat(&initial, mass, tg, tw);
    let trial = endpoint.evaluate(h(tw), mass, tg, heat, || false).unwrap();
    endpoint.check_balance(&trial).unwrap();
    let candidate = endpoint.finalize(&trial, || false).unwrap();
    let step = candidate.step();
    assert_eq!(candidate.film(), GasWallFilmInventory::Dry);
    assert_eq!(candidate.film().specific_enthalpy_j_kg(), None);
    assert_eq!(step.liquid_mass_kg, [mass, 0.0]);
    assert_eq!(step.transferred_mass_kg, mass);
    assert_eq!(step.duration_s, dt);
    near(step.requested_transferred_mass_kg, 3.0 * mass, 1e-14);
    // This is the event for the endpoint-frozen rate, not a claim about the
    // continuous trajectory's exhaustion time inside the numerical interval.
    near(step.dry_out_time_s.unwrap(), dt / 3.0, 1e-8);
    near(step.phase_active_duration_s, dt / 3.0, 1e-8);
    near(
        step.wall_heat_into_system_j,
        mass * (vapor_h(tw) - h(295.0)),
        1e-8,
    );
    near(step.film_sensible_energy_j, mass * dh(295.0, tw), 1e-10);
    near(step.film_enthalpy_change_j, -mass * h(295.0), 1e-10);
    near(
        step.film_enthalpy_change_j + step.gas_internal_energy_change_j
            - step.wall_heat_into_system_j
            - heat,
        0.0,
        1e-8,
    );
    let dry = GasWallFilmEndpoint::new(
        surface(&binding, GasWallFilmInventory::Dry, tw, k),
        candidate.gas(),
        &source,
        dt,
        controls(),
    )
    .unwrap();
    let idle = dry
        .evaluate(
            h(tw),
            0.0,
            candidate.gas().state().state().temperature,
            0.0,
            || false,
        )
        .unwrap();
    dry.check_balance(&idle).unwrap();
    assert!(idle.step().requested_transferred_mass_kg > 0.0);
    assert_eq!(idle.step().transferred_mass_kg, 0.0);
    assert_eq!(idle.step().phase_active_duration_s, 0.0);
    assert_eq!(idle.step().wall_heat_into_system_j, 0.0);
    let idle = dry.finalize(&idle, || false).unwrap();
    assert_eq!(idle.film(), GasWallFilmInventory::Dry);
    assert_eq!(idle.gas(), candidate.gas());
}

#[test]
fn independent_gas_and_liquid_datums_leave_wet_film_capacity_and_actual_phase_heat_invariant() {
    let mut baseline: Option<[f64; 4]> = None;
    for (gas_shift, liquid_shift) in [(0.0, 0.0), (1e12, 1_048_576.0)] {
        let source = catalog(gas_shift);
        let binding = binding(gas_shift - liquid_shift);
        let initial = gas(&source, 320.0, 5000.0, 0.001);
        let old_mass = 1e-5;
        let dm = -0.05 * initial.species_mass_kg(&vapor()).unwrap();
        let (tw, tg, dt) = (304.0, 316.0, 0.5);
        let k = coefficient(&initial, dm, dm, tg, tw, dt);
        let heat = required_gas_heat(&initial, dm, tg, tw);
        let inverse = |value| temperature(value - liquid_shift);
        let input = GasWallFilmSurface {
            film: GasWallFilmInventory::Wet {
                mass_kg: old_mass,
                specific_enthalpy_j_kg: h(300.0) + liquid_shift,
            },
            interface_liquid_enthalpy_j_kg: h(300.0) + liquid_shift,
            temperature_from_enthalpy: &inverse,
            ..surface(&binding, GasWallFilmInventory::Dry, 300.0, k)
        };
        let mut config = controls();
        config.latent_compatibility_tolerance_j_kg = 0.01;
        let endpoint = GasWallFilmEndpoint::new(input, &initial, &source, dt, config).unwrap();
        let trial = endpoint
            .evaluate(h(tw) + liquid_shift, dm, tg, heat, || false)
            .unwrap();
        endpoint.check_balance(&trial).unwrap();
        let candidate = endpoint.finalize(&trial, || false).unwrap();
        let step = candidate.step();
        near(candidate.film().mass_kg(), old_mass - dm, 1e-14);
        near(
            inverse(candidate.film().specific_enthalpy_j_kg().unwrap()).unwrap(),
            tw,
            1e-10,
        );
        near(step.film_sensible_energy_j, old_mass * dh(300.0, tw), 1e-9);
        near(
            step.wall_heat_into_system_j,
            old_mass * dh(300.0, tw) + dm * (LATENT + dh(tw, tg)),
            1e-8,
        );
        near(step.latent_energy_j, dm * LATENT, 1e-8);
        near(step.energy_into_gas_j, dm * (vapor_h(tg) + gas_shift), 1e-6);
        assert!(step.sum_absolute_energy_residual_j <= config.gas.energy_tolerance_j);
        let actual = [
            candidate.gas().state().state().temperature,
            candidate.film().mass_kg(),
            step.wall_heat_into_system_j,
            step.latent_energy_j,
        ];
        if let Some(baseline) = &baseline {
            for i in 0..4 {
                near(actual[i], baseline[i], if i == 1 { 1e-14 } else { 2e-8 });
            }
        } else {
            baseline = Some(actual);
        }
    }
}

#[test]
fn unbalanced_trials_original_source_refusal_and_final_cancellation_publish_nothing() {
    let source = catalog(0.0);
    let binding = binding(0.0);
    let initial = gas(&source, 320.0, 5000.0, 0.001);
    let snapshot = initial.clone();
    let dm = -0.1 * initial.species_mass_kg(&vapor()).unwrap();
    let (tw, tg, dt) = (300.0, 314.0, 0.5);
    let k = coefficient(&initial, dm, dm, tg, tw, dt);
    let heat = required_gas_heat(&initial, dm, tg, tw);
    let input = surface(&binding, GasWallFilmInventory::Dry, tw, k);
    let endpoint = GasWallFilmEndpoint::new(input, &initial, &source, dt, controls()).unwrap();
    let wrong_mass = endpoint
        .evaluate(h(tw), dm + 1e-7, tg, heat, || false)
        .unwrap();
    assert!(wrong_mass.maximum_mass_residual_kg() > 1e-8);
    assert!(endpoint.check_balance(&wrong_mass).is_err());
    assert!(endpoint.finalize(&wrong_mass, || false).is_err());
    let wrong_heat = endpoint
        .evaluate(h(tw), dm, tg, heat + 1.0, || false)
        .unwrap();
    near(wrong_heat.gas_energy_residual_j().abs(), 1.0, 1e-8);
    assert!(endpoint.check_balance(&wrong_heat).is_err());
    assert!(endpoint.finalize(&wrong_heat, || false).is_err());
    assert!(endpoint.evaluate(h(950.0), dm, tg, heat, || false).is_err());
    let mut incompatible = binding.clone();
    incompatible.liquid_enthalpy_offset_j_kg += 1000.0;
    assert!(
        GasWallFilmEndpoint::new(
            GasWallFilmSurface {
                binding: &incompatible,
                ..input
            },
            &initial,
            &source,
            dt,
            controls()
        )
        .is_err()
    );
    let trial = endpoint.evaluate(h(tw), dm, tg, heat, || false).unwrap();
    let mut polls = 0usize;
    endpoint
        .finalize(&trial, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut count = 0usize;
    let cancelled = endpoint.finalize(&trial, || {
        count += 1;
        count == polls
    });
    assert!(matches!(cancelled, Err(GasEvaporationError::Cancelled)));
    assert_eq!(count, polls);
    assert_eq!(initial, snapshot);
    assert_eq!(input.film, GasWallFilmInventory::Dry);
}
