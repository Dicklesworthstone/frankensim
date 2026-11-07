//! G0/G1 conservative-state inversion; continuing-vessel tests live in parcel.
use super::*;
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::nonideal::PengRobinsonSpeciesBinding;
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_qty::Pressure;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis, SpeciesId};

fn bound() -> ThermochemicalPengRobinsonModel {
    bound_with_reference(0.0)
}

pub(super) fn bound_with_reference(reference_j_kg: f64) -> ThermochemicalPengRobinsonModel {
    let source = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.028, 3.5, 0.001, reference_j_kg),
                species_with_energy_reference("B", 0.004, 2.5, 0.0, reference_j_kg),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        ).unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.15)],
        None,
    ).unwrap();
    ThermochemicalPengRobinsonModel::new(source, vec![
        PengRobinsonSpeciesBinding { species: SpeciesId::new("A").unwrap(),
            critical_temperature: Temperature::new(190.564),
            critical_pressure: Pressure::new(4.5992e6), acentric_factor: 0.01142 },
        PengRobinsonSpeciesBinding { species: SpeciesId::new("B").unwrap(),
            critical_temperature: Temperature::new(150.687),
            critical_pressure: Pressure::new(4.863e6), acentric_factor: -0.00219 },
    ], vec![0.0, 0.03, 0.03, 0.0]).unwrap()
}
fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 280.0, maximum_temperature_k: 450.0,
        energy_tolerance_j_kg: 1.0e-7, maximum_iterations: 128 }
}

#[test]
fn g1_absolute_and_reference_energy_recover_one_conserved_density_state() {
    let model = bound();
    let density = 47.889_003_167_633_395;
    let context = point(300.0, 5.0e6).with("lot", 7.0).unwrap();
    let final_point = point(360.0, 5.0e6).with("lot", 7.0).unwrap();
    let expected = model.evaluate_density_temperature(&final_point, density).unwrap();
    let initial = model.evaluate_density_temperature(&context, density).unwrap();
    let absolute = model.reconstruct_density_energy(&context, density,
        expected.specific_internal_energy_j_kg(), controls(), || false).unwrap();
    let relative = model.reconstruct_density_energy_change(&context, density,
        initial.specific_internal_energy_j_kg(),
        expected.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg(),
        controls(), || false).unwrap();
    for solved in [absolute, relative] {
        assert_eq!(solved.state().density_kg_m3().to_bits(), density.to_bits());
        assert!((solved.state().temperature_k() - 360.0).abs() < 2.0e-9);
        assert!((solved.state().pressure_pa() - 6_226_134.941_447_334).abs() < 1.0e-4);
        assert!(solved.residual_j_kg().abs() <= controls().energy_tolerance_j_kg);
        assert_eq!(solved.state().source().components()[0].query().axes().get("lot"), Some(&7.0));
    }
}

#[test]
fn g0_invalid_density_unbracketed_energy_and_cancellation_do_not_return_states() {
    let model = bound();
    let context = point(300.0, 5.0e6);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(model.reconstruct_density_energy(&context, bad, 100_000.0, controls(), || false).is_err());
    }
    assert!(matches!(model.reconstruct_density_energy(&context, 50.0, 1.0e12,
        controls(), || false), Err(GasEnergyError::OutsideBracket { .. })));
    assert!(matches!(model.reconstruct_density_energy(&context, 50.0, 100_000.0,
        controls(), || true), Err(GasEnergyError::Cancelled)));
    assert!(model.reconstruct_density_internal_energy_increment(&point(200.0, 5.0e6),
        50.0, 1.0, controls(), || false).is_err());
    let zero = model.reconstruct_density_internal_energy_increment(&context, 50.0,
        0.0, controls(), || false).unwrap();
    assert_eq!(zero.state().temperature_k(), 300.0);
    assert_eq!(zero.state().density_kg_m3(), 50.0);
    assert_eq!(zero.iterations(), 0);
}

#[test]
fn g1_joule_expansion_retains_original_density_departure_and_matches_oracle() {
    let model = bound();
    let rho0 = 47.889_003_167_633_395;
    let rho1 = rho0 / 1.25;
    let context = point(300.0, 5.0e6);
    // Independent 75-digit Decimal PR energy balance, retaining the original
    // residual energy at rho0: Delta u_ideal + uR(rho1,T1)-uR(rho0,T0)=0.
    // This is a constitutive endpoint; the parcel owns free-expansion admission.
    let solved = model.reconstruct_density_internal_energy_increment_between_densities(
        &context, rho0, rho1, 0.0, controls(), || false).unwrap();
    let state = solved.state();
    assert_eq!(state.density_kg_m3().to_bits(), rho1.to_bits());
    assert!((state.temperature_k() - 294.872_058_857_231_76).abs() < 2.0e-9);
    assert!((state.pressure_pa() - 3_975_006.042_948_880_2).abs() < 1.0e-4);
    assert!(solved.residual_j_kg().abs() <= controls().energy_tolerance_j_kg);
    let entropy = model.specific_entropy_change_between_densities(rho0, 300.0,
        rho1, state.temperature_k()).unwrap();
    assert!((entropy - 78.087_722_432_685_2).abs() < 1.0e-8);
    // Fixed-u expansion is not fixed-T for this nonideal material.
    assert!((state.temperature_k() - 300.0).abs() > 5.0);
    let returned = model.reconstruct_density_internal_energy_increment_between_densities(
        &point(state.temperature_k(), state.pressure_pa()), rho1, rho0,
        0.0, controls(), || false).unwrap();
    assert!((returned.state().temperature_k() - 300.0).abs() < 2.0e-9);
}

#[test]
fn g1_general_energy_increment_keeps_both_density_endpoints_and_heat_direction() {
    let model = bound();
    let (rho0, rho1) = (40.0, 60.0);
    let context = point(300.0, 5.0e6);
    let isothermal_change = model.specific_internal_energy_change_between_densities(
        rho0, 300.0, rho1, 300.0).unwrap();
    assert!(isothermal_change.abs() > 1000.0);
    let isothermal = model.reconstruct_density_internal_energy_increment_between_densities(
        &context, rho0, rho1, isothermal_change, controls(), || false).unwrap();
    assert!((isothermal.state().temperature_k() - 300.0).abs() < 2.0e-9);
    let target = model.specific_internal_energy_change_between_densities(
        rho0, 300.0, rho1, 350.0).unwrap();
    let returned = model.specific_internal_energy_change_between_densities(
        rho1, 350.0, rho0, 300.0).unwrap();
    assert_eq!(returned, -target);
    let solved = model.reconstruct_density_internal_energy_increment_between_densities(
        &context, rho0, rho1, target, controls(), || false).unwrap();
    assert!((solved.state().temperature_k() - 350.0).abs() < 2.0e-9);
    // Existing fixed-density interface delegates to the same conserved law.
    assert_eq!(model.specific_internal_energy_change_at_density(rho0, 300.0, 350.0).unwrap(),
        model.specific_internal_energy_change_between_densities(rho0, 300.0, rho0, 350.0).unwrap());
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(model.reconstruct_density_internal_energy_increment_between_densities(
            &context, bad, rho1, target, controls(), || false).is_err());
        assert!(model.reconstruct_density_internal_energy_increment_between_densities(
            &context, rho0, bad, target, controls(), || false).is_err());
    }
}

#[test]
fn paired_energy_increment_retains_original_endpoint_cv_and_formation_invariance() {
    let model = bound();
    for ((rho0, t0), (rho1, t1)) in [
        ((40.0, 300.0), (60.0, 350.0)),
        ((60.0, 350.0), (40.0, 300.0)),
        ((40.0, 300.0), (60.0, 300.0)),
        ((40.0, 300.0), (40.0, 300.0)),
        ((47.889_003_167_633_395, 300.0), (47.889_003_172_422_3, 300.000_000_01)),
    ] {
        let paired = model.specific_internal_energy_change_and_cv_between_densities(
            rho0, t0, rho1, t1,
        ).unwrap();
        let full = model.evaluate_density_temperature(&point(t1, 5e6), rho1).unwrap();
        let cv = paired.final_cv_j_kg_k();
        assert!(cv > 0.0);
        assert!((cv / full.specific_heat_cv_j_kg_k() - 1.0).abs() < 128.0 * f64::EPSILON);
        let dt = 0.001;
        let plus = model.specific_internal_energy_change_between_densities(
            rho0, t0, rho1, t1 + dt,
        ).unwrap();
        let minus = model.specific_internal_energy_change_between_densities(
            rho0, t0, rho1, t1 - dt,
        ).unwrap();
        assert!(((plus - minus) / (2.0 * dt) / cv - 1.0).abs() < 1e-9);
        assert_eq!(model.specific_internal_energy_change_between_densities(
            rho1, t1, rho0, t0,
        ).unwrap(), -paired.internal_energy_change_j_kg());
        for datum in [-1e24, 1e24] {
            let shifted = bound_with_reference(datum);
            assert_eq!(shifted.specific_internal_energy_change_and_cv_between_densities(
                rho0, t0, rho1, t1,
            ).unwrap(), paired);
        }
    }
}

#[test]
fn paired_energy_increment_keeps_both_endpoint_domains() {
    let model = bound();
    for density in [0.0, -1.0, 1e9, f64::NAN, f64::INFINITY] {
        assert!(model.specific_internal_energy_change_and_cv_between_densities(
            density, 300.0, 40.0, 350.0,
        ).is_err());
        assert!(model.specific_internal_energy_change_and_cv_between_densities(
            40.0, 300.0, density, 350.0,
        ).is_err());
    }
    for temperature in [199.0, 2001.0, f64::NAN, f64::INFINITY] {
        assert!(model.specific_internal_energy_change_and_cv_between_densities(
            40.0, temperature, 60.0, 350.0,
        ).is_err());
        assert!(model.specific_internal_energy_change_and_cv_between_densities(
            40.0, 300.0, 60.0, temperature,
        ).is_err());
    }
}

#[test]
fn original_law_enthalpy_closes_energy_and_flow_work_and_reverses() {
    let model = bound();
    let rho = 47.889_003_167_633_395;
    for ((rho0, t0), (rho1, t1)) in [
        ((40.0, 300.0), (60.0, 350.0)),
        ((rho, 300.0), (rho, 360.0)),
        ((rho, 300.0), (0.5 * rho, 220.692_821_099_187_96)),
        ((rho, 300.0), (rho * (1.0 + 1e-8), 300.000_001)),
    ] {
        let initial = model.evaluate_density_temperature(&point(t0, 5e6), rho0).unwrap();
        let final_state = model.evaluate_density_temperature(&point(t1, 5e6), rho1).unwrap();
        let change = model.specific_enthalpy_change_between_densities(rho0, t0, rho1, t1).unwrap();
        let du = model.specific_internal_energy_change_between_densities(rho0, t0, rho1, t1).unwrap();
        let flow_work = final_state.pressure_pa() / rho1 - initial.pressure_pa() / rho0;
        assert!((change - du - flow_work).abs() < 2e-8);
        // With the unshifted datum, absolute subtraction is independently
        // resolved well enough to check h = u + p/rho at the material boundary.
        assert!((change - (final_state.specific_enthalpy_j_kg()
            - initial.specific_enthalpy_j_kg())).abs() < 2e-8);
        assert_eq!(model.specific_enthalpy_change_between_densities(rho1, t1, rho0, t0).unwrap(), -change);
        assert_eq!(model.specific_enthalpy_change_between_densities(rho0, t0, rho0, t0).unwrap(), 0.0);
        for reference in [-1e24, 1e24] {
            let shifted = bound_with_reference(reference);
            assert_eq!(shifted.specific_enthalpy_change_between_densities(rho0, t0, rho1, t1).unwrap(), change);
        }
    }
}

#[test]
fn weak_real_gas_energy_and_enthalpy_match_independent_finite_changes() {
    let rho = 47.889_003_167_633_395;
    // 85-digit Decimal PR/NASA differences with the exact f64 coordinates.
    // These changes are thousands of coordinate ULPs above resolution, but
    // subtraction of rounded residual energy/pressure loses about 0.02%.
    for reference in [0.0, -1e24, 1e24] {
        let model = bound_with_reference(reference);
        for (t1, rho1, expected_u, expected_h) in [
            (300.000_000_01, rho, 9.591_597_381_493_554_6e-6, 1.387_949_738_165_587_2e-5),
            (300.0, 47.889_003_167_681_29, -2.422_959_241_806_747_6e-8, -3.137_111_222_645_165_1e-8),
            (300.000_000_01, 47.889_003_172_422_3, 7.168_717_227_411_026_5e-6, 1.074_248_855_739_969_5e-5),
        ] {
            let du = model.specific_internal_energy_change_between_densities(rho, 300.0, rho1, t1).unwrap();
            let dh = model.specific_enthalpy_change_between_densities(rho, 300.0, rho1, t1).unwrap();
            assert!((du / expected_u - 1.0).abs() < 1e-12);
            assert!((dh / expected_h - 1.0).abs() < 1e-12);
            assert_eq!(model.specific_internal_energy_change_between_densities(rho1, t1, rho, 300.0).unwrap(), -du);
            assert_eq!(model.specific_enthalpy_change_between_densities(rho1, t1, rho, 300.0).unwrap(), -dh);
        }
    }
}

#[test]
fn weak_heating_reconstructs_without_relaxing_the_requested_energy_tolerance() {
    let rho = 47.889_003_167_633_395;
    let target = 9.591_597_381_493_554_6e-6; // Independent finite-energy oracle.
    let controls = GasEnergySolve { energy_tolerance_j_kg: 1e-15, ..controls() };
    for reference in [0.0, -1e24, 1e24] {
        let model = bound_with_reference(reference);
        let solved = model.reconstruct_density_internal_energy_increment(
            &point(300.0, 5e6), rho, target, controls, || false).unwrap();
        assert_eq!(solved.state().temperature_k().to_bits(), 300.000_000_01_f64.to_bits());
        assert!(solved.residual_j_kg().abs() <= controls.energy_tolerance_j_kg);
        assert_eq!(solved.state().density_kg_m3().to_bits(), rho.to_bits());
    }
}
