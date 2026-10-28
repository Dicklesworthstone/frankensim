use super::*;
use crate::thermochemical_gas::tests::{model, point};

fn energy() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn controls() -> GasNozzleStepSolve {
    GasNozzleStepSolve { nozzle: energy(), maximum_donor_fraction: 0.1,
        transfer: GasMassTransferSolve { donor: energy(), receiver: energy(),
            mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 } }
}
fn nozzle() -> GasNozzle { GasNozzle { area_m2: 1e-6, discharge_coefficient: 0.8 } }
fn gas(x: f64, slope: f64, t: f64, p: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x, 1.0-x], slope), point(t, p), 0.001).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a-b).abs() <= tol*b.abs().max(1e-12), "{a} != {b}");
}

#[test]
fn source_model_matches_calorically_perfect_subsonic_and_choked_references() {
    let upstream = gas(0.75, 0.0, 500.0, 200000.0);
    let source = upstream.state().state();
    let gamma = source.gamma;
    let critical = (2.0/(gamma+1.0)).powf(gamma/(gamma-1.0));
    for ratio in [0.95_f64, 0.8, 0.55, 0.5, 0.2] {
        let flow = upstream.nozzle_flow_to(200000.0*ratio, nozzle(), energy(), || false).unwrap();
        let exit_ratio = ratio.max(critical);
        let t = 500.0*exit_ratio.powf((gamma-1.0)/gamma);
        let rho = 200000.0*exit_ratio/(source.specific_gas_constant*t);
        let v = (2.0*source.specific_heat_cp*(500.0-t)).sqrt();
        near(flow.throat().state().temperature, t, 1e-9);
        near(flow.throat().state().pressure, 200000.0*exit_ratio, 1e-9);
        near(flow.velocity_m_s, v, 1e-9);
        near(flow.mass_flow_kg_s, 0.8e-6*rho*v, 1e-9);
        assert_eq!(flow.choked, ratio <= critical);
        if flow.choked { near(flow.mach_number, 1.0, 1e-9); }
        assert!(flow.entropy_energy_residual_j_kg.abs() <= 2.0*energy().energy_tolerance_j_kg);
    }
    let a = upstream.nozzle_flow_to(40000.0, nozzle(), energy(), || false).unwrap();
    let b = upstream.nozzle_flow_to(20000.0, nozzle(), energy(), || false).unwrap();
    near(a.mass_flow_kg_s, b.mass_flow_kg_s, 1e-9);
}

#[test]
fn variable_heat_capacity_changes_choking_without_a_frozen_gamma_substitution() {
    let upstream = gas(0.75, 0.002, 800.0, 200000.0);
    let flow = upstream.nozzle_flow_to(40000.0, nozzle(), energy(), || false).unwrap();
    assert!(flow.choked);
    // Independent integral of the synthetic linear cp/R curve, not a nozzle
    // lookup or rescaling of the old acoustic speed.
    let (a, b) = (3.25, 0.0015);
    let (mut lo, mut hi) = (200.0, 800.0);
    for _ in 0..90 {
        let t = 0.5*(lo+hi);
        let cp = a+b*t;
        let residual = a*(t-800.0)+0.5*b*(t*t-800.0*800.0)+0.5*cp/(cp-1.0)*t;
        if residual < 0.0 { lo = t; } else { hi = t; }
    }
    let t: f64 = 0.5*(lo+hi);
    let p = 200000.0*(a*(t/800.0).ln()+b*(t-800.0)).exp();
    near(flow.throat().state().temperature, t, 1e-9);
    near(flow.throat().state().pressure, p, 1e-9);
    let frozen_t = 1600.0/(upstream.state().state().gamma+1.0);
    assert!((t-frozen_t).abs() > 1.0);
    assert_ne!(flow.throat().state().gamma, upstream.state().state().gamma);
    let other = gas(0.25, 0.002, 800.0, 200000.0)
        .nozzle_flow_to(40000.0, nozzle(), energy(), || false).unwrap();
    assert_ne!(flow.mass_flow_kg_s, other.mass_flow_kg_s);
    assert_ne!(flow.throat().state().sound_speed, other.throat().state().sound_speed);
}

#[test]
fn area_scaling_and_closed_paths_do_not_mutate_sources_or_invent_flow() {
    let upstream = gas(0.75, 0.001, 500.0, 200000.0);
    let a = upstream.nozzle_flow_to(100000.0, nozzle(), energy(), || false).unwrap();
    let b = upstream.nozzle_flow_to(100000.0, GasNozzle { area_m2: 2e-6, ..nozzle() }, energy(), || false).unwrap();
    assert_eq!(a.throat(), b.throat());
    near(b.mass_flow_kg_s, 2.0*a.mass_flow_kg_s, 1e-12);
    for (p, area, coefficient) in [(200000.0, 1e-6, 0.8), (100000.0, 0.0, 0.8), (100000.0, 1e-6, 0.0)] {
        let flow = upstream.nozzle_flow_to(p, GasNozzle { area_m2: area, discharge_coefficient: coefficient }, energy(), || false).unwrap();
        assert_eq!(flow.throat(), upstream.state());
        assert_eq!(flow.mass_flow_kg_s, 0.0);
        assert_eq!(flow.iterations, 0);
    }
    // A small pressure drop needs only the actual subsonic source interval.
    let narrow = GasEnergySolve { minimum_temperature_k: 490.0, ..energy() };
    assert!(upstream.nozzle_flow_to(198000.0, nozzle(), narrow, || false).is_ok());
    assert!(upstream.nozzle_flow_to(20000.0, nozzle(), narrow, || false).is_err());
}

#[test]
fn pressure_driven_pair_transfers_actual_rate_species_and_total_enthalpy_once() {
    let donor = gas(0.75, 0.001, 500.0, 300000.0);
    let recipient = gas(0.25, 0.001, 400.0, 100000.0);
    let (trial, flow) = donor.through_nozzle_to(&recipient, nozzle(), 0.0002, controls(), || false).unwrap();
    assert_eq!(trial.step().transferred_mass_kg, flow.mass_flow_kg_s*0.0002);
    assert_eq!(trial.step().source_enthalpy_j_kg, flow.upstream_enthalpy_j_kg);
    assert!(trial.donor().state().state().pressure < donor.state().state().pressure);
    assert!(trial.receiver().state().state().pressure > recipient.state().state().pressure);
    assert!(trial.step().energy_residual_j.abs() <= controls().transfer.energy_tolerance_j);
    assert!(trial.step().maximum_species_mass_residual_kg <= controls().transfer.mass_tolerance_kg);
    assert_ne!(trial.receiver().state().state().specific_gas_constant, recipient.state().state().specific_gas_constant);
}

#[test]
fn finite_tank_flow_refines_and_retains_every_accepted_inventory() {
    let initial_a = gas(0.75, 0.001, 500.0, 300000.0);
    let initial_b = gas(0.25, 0.001, 400.0, 100000.0);
    let mut pressures = Vec::new();
    for n in [16, 32, 64] {
        let (mut a, mut b) = (initial_a.clone(), initial_b.clone());
        for _ in 0..n {
            let (trial, _) = a.through_nozzle_to(&b, nozzle(), 0.04/f64::from(n), controls(), || false).unwrap();
            let (da, db, step) = trial.into_parts();
            assert!(step.energy_residual_j.abs() < 1e-8);
            a = da; b = db;
        }
        pressures.push(a.state().state().pressure);
        assert!(((a.mass_kg()-initial_a.mass_kg())+(b.mass_kg()-initial_b.mass_kg())).abs() < 1e-14);
        assert!((super::super::stored_change(&initial_a, &a)+super::super::stored_change(&initial_b, &b)).abs() < 1e-7);
    }
    assert!((pressures[0]-pressures[1]).abs() > 1.9*(pressures[1]-pressures[2]).abs());
}

#[test]
fn source_limits_oversized_steps_and_late_cancellation_refuse_without_half_updates() {
    let a = gas(0.75, 0.001, 500.0, 300000.0);
    let b = gas(0.25, 0.001, 400.0, 100000.0);
    let old = (a.clone(), b.clone());
    for p in [-1.0, f64::NAN, f64::INFINITY, 400000.0] {
        assert!(a.nozzle_flow_to(p, nozzle(), energy(), || false).is_err());
    }
    assert!(a.through_nozzle_to(&b, nozzle(), 1000.0, controls(), || false).is_err());
    // Less than 0.2% of the donor mass still overshoots near equilibrium:
    // inventory caps alone are NOT a pressure-direction stability check.
    let almost_equal = gas(0.25, 0.001, 400.0, 299000.0);
    assert!(matches!(a.through_nozzle_to(&almost_equal, nozzle(), 0.04, controls(), || false),
        Err(GasEnergyError::Invalid("nozzle timestep crosses pressure equilibrium; reduce the physical step"))));
    for bad in [GasNozzle { area_m2: -1.0, ..nozzle() },
        GasNozzle { discharge_coefficient: f64::NAN, ..nozzle() },
        GasNozzle { discharge_coefficient: 1.1, ..nozzle() }] {
        assert!(a.nozzle_flow_to(100000.0, bad, energy(), || false).is_err());
    }
    let limited = GasEnergySolve { maximum_iterations: 1, ..energy() };
    assert!(a.nozzle_flow_to(100000.0, nozzle(), limited, || false).is_err());
    let mut total = 0;
    a.through_nozzle_to(&b, nozzle(), 0.0002, controls(), || { total += 1; false }).unwrap();
    for stop in [1, total/2, total] {
        let mut count = 0;
        assert!(a.through_nozzle_to(&b, nozzle(), 0.0002, controls(), || { count += 1; count == stop }).is_err());
        assert_eq!((&a, &b), (&old.0, &old.1));
    }
}

#[test]
fn vacuum_back_pressure_chokes_without_querying_a_zero_pressure_material_state() {
    let upstream = gas(0.75, 0.002, 800.0, 300000.0);
    let vacuum = upstream.nozzle_flow_to(0.0, nozzle(), energy(), || false).unwrap();
    let low = upstream.nozzle_flow_to(30000.0, nozzle(), energy(), || false).unwrap();
    assert!(vacuum.choked);
    assert!(vacuum.throat().state().pressure > 0.0);
    near(vacuum.mach_number, 1.0, 1e-9);
    assert_eq!(vacuum.throat(), low.throat());
    assert_eq!(vacuum.mass_flow_kg_s, low.mass_flow_kg_s);
}

#[test]
fn formation_reference_does_not_change_nozzle_speed_or_choking() {
    use crate::thermochemical_gas::{ThermochemicalGasModel,
        tests::{species_with_energy_reference, transport}};
    use fs_thermochem::{mixture::FrozenIdealGasMixtureModelV1, Composition, CompositionBasis};
    let parcel = |reference| {
        let thermodynamics = FrozenIdealGasMixtureModelV1::new(
            vec![species_with_energy_reference("A", 0.028, 3.5, 0.002, reference),
                species_with_energy_reference("B", 0.004, 2.5, 0.0, reference)],
            Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap()).unwrap();
        let model = ThermochemicalGasModel::new(thermodynamics,
            vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.15)], None).unwrap();
        ClosedGasParcel::new(model, point(800.0, 300000.0), 0.001).unwrap()
    };
    let ordinary = parcel(0.0);
    let shifted = parcel(1e24);
    for back in [0.0, 40000.0, 240000.0, 299999.0] {
        let a = ordinary.nozzle_flow_to(back, nozzle(), energy(), || false).unwrap();
        let b = shifted.nozzle_flow_to(back, nozzle(), energy(), || false).unwrap();
        assert_eq!(a.throat().state(), b.throat().state());
        assert_eq!(a.velocity_m_s, b.velocity_m_s);
        assert_eq!(a.mass_flow_kg_s, b.mass_flow_kg_s);
        assert_eq!(a.choked, b.choked);
        assert_ne!(a.upstream_enthalpy_j_kg, b.upstream_enthalpy_j_kg);
    }
}
