//! Native G1/G4 paths with retained material density and resolved gas state.
use super::*;
use super::super::tests::body;
use crate::cavity_phs::piston::tests::gas as parcel;
use fs_gas::dynamics::{Boundaries, Boundary, CaloricDomain, FrozenGas};
use fs_gas::multicomponent::MixtureGas;
use fs_qty::{Length, Temperature};
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;

const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting, right: Boundary::Reflecting,
};
const DT: Time = Time::new(1e-6);
fn domain() -> CaloricDomain {
    CaloricDomain::new(Temperature::new(200.0), Temperature::new(1800.0))
}
fn frozen() -> Euler1d {
    let source = parcel(0.75, 0.001).model().thermodynamics().clone();
    let gas = FrozenGas::new(source, domain(), || false).unwrap();
    let q = gas.conserved(100000.0 / (400.0 * gas.gas_constant()), 0.0, 100000.0).unwrap();
    Euler1d::new(gas, Length::new(0.01), 0.01, vec![q; 8], CLOSED, || false).unwrap()
}
fn mixture() -> MixtureEuler1d {
    let models = parcel(0.75, 0.001).model().thermodynamics().components().to_vec();
    let ids: Vec<_> = models.iter().map(|model| model.species().clone()).collect();
    let source = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(vec![ElementId::new("H").unwrap()], ids.clone(),
            vec![vec![1; ids.len()]]).unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(), models,
    ).unwrap();
    let gas = MixtureGas::new(source, domain(), || false).unwrap();
    let states = (0..8).map(|i| {
        let f = if i < 4 { 0.8 } else { 0.2 };
        let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * (f / 0.028 + (1.0-f) / 0.004);
        let rho = 100000.0 / (400.0 * r);
        gas.conserved(&[rho*f, rho*(1.0-f)], 0.0, 100000.0).unwrap()
    }).collect();
    MixtureEuler1d::new(gas, Length::new(0.01), 0.01, states, CLOSED, || false).unwrap()
}
fn controls() -> PistonControl {
    PistonControl { midpoint_velocity_bracket_m_s: [-1.0, 1.0],
        momentum_tolerance_n_s: 1e-10, energy_tolerance_j: 1e-8, maximum_iterations: 96 }
}
fn load() -> PistonLoad {
    PistonLoad { exterior_pressure_pa: 105000.0, external_force_n: -1.0 }
}
fn suspension() -> PistonSuspension {
    PistonSuspension::new(8000.0, 0.8, 0.079).unwrap()
}
fn equilibrium() -> ContactEquilibriumControl {
    ContactEquilibriumControl { relative_tolerance: 1e-10, maximum_energy_defect_j: 1e-3 }
}
fn scratch() -> ContactWorkspace {
    ContactWorkspace::new(8, ContactWorkspace::required_workspace_bytes(8).unwrap(), || false).unwrap()
}

#[test]
fn sourced_mass_changes_resolved_motion_and_matches_both_native_frozen_paths() {
    for suspended in [false, true] {
        let mut speeds = Vec::new();
        for density in [1000.0, 8000.0] {
            let specimen = body(density);
            let mut wall = MaterialResolvedPiston::new(specimen.clone(), 0.0).unwrap();
            let mut native = Piston::new(specimen.inertia().mass(), 0.0).unwrap();
            let (mut gas, mut reference) = (frozen(), frozen());
            let g = StepControl::default();
            let report = if suspended {
                let actual = wall.step_suspended(&mut gas, DT, g, load(), suspension(), controls(), || false).unwrap();
                let expected = reference.step_with_suspended_piston(&mut native, DT, g, load(), suspension(), controls(), || false).unwrap();
                assert_eq!(actual, expected);
                assert!(actual.damping_energy_j >= 0.0);
                actual.piston
            } else {
                let actual = wall.step(&mut gas, DT, g, load(), controls(), || false).unwrap();
                let expected = reference.step_with_piston(&mut native, DT, g, load(), controls(), || false).unwrap();
                assert_eq!(actual, expected);
                actual
            };
            assert_eq!(wall.body(), &specimen);
            assert_eq!(wall.piston(), &native);
            assert_eq!(wall.piston().mass_kg(), wall.body().inertia().mass());
            assert_eq!(gas.states(), reference.states());
            assert_eq!(gas.mesh(), reference.mesh());
            assert_eq!(gas.time(), reference.time());
            assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
            assert!(report.combined_momentum_residual_n_s.abs() <= controls().momentum_tolerance_n_s);
            assert!(report.gas.right_boundary_change.energy > 0.0);
            speeds.push(wall.piston().velocity_m_s());
        }
        assert!(speeds[0] < speeds[1] && speeds[1] < 0.0);
    }
}

#[test]
fn every_mixture_route_retains_native_species_fields_and_complete_reports() {
    for mode in 0..4 {
        let specimen = body(4000.0);
        let mut wall = MaterialResolvedPiston::new(specimen.clone(), 0.08).unwrap();
        let mut native = Piston::new(specimen.inertia().mass(), 0.08).unwrap();
        let (mut gas, mut reference) = (mixture(), mixture());
        let before = gas.totals();
        let (mut a, mut b) = (scratch(), scratch());
        let g = StepControl::default();
        match mode {
            0 => assert_eq!(
                wall.step_mixture(&mut gas, DT, g, load(), controls(), || false).unwrap(),
                reference.step_with_piston(&mut native, DT, g, load(), controls(), || false).unwrap()),
            1 => assert_eq!(
                wall.step_mixture_suspended(&mut gas, DT, g, load(), suspension(), controls(), || false).unwrap(),
                reference.step_with_suspended_piston(&mut native, DT, g, load(), suspension(), controls(), || false).unwrap()),
            2 => assert_eq!(
                wall.step_mixture_with_contact_equilibrium(&mut gas, DT, g, load(), controls(), equilibrium(), &mut a, || false).unwrap(),
                reference.step_with_piston_and_contact_equilibrium(&mut native, DT, g, load(), controls(), equilibrium(), &mut b, || false).unwrap()),
            _ => assert_eq!(
                wall.step_mixture_suspended_with_contact_equilibrium(&mut gas, DT, g, load(), suspension(), controls(), equilibrium(), &mut a, || false).unwrap(),
                reference.step_with_suspended_piston_and_contact_equilibrium(&mut native, DT, g, load(), suspension(), controls(), equilibrium(), &mut b, || false).unwrap()),
        }
        assert_eq!(wall.body(), &specimen);
        assert_eq!(wall.piston(), &native);
        assert_eq!(gas.mesh(), reference.mesh());
        assert_eq!(gas.time(), reference.time());
        for i in 0..gas.cell_count() { assert_eq!(gas.cell(i), reference.cell(i)); }
        for (after, initial) in gas.totals().species_masses.iter().zip(&before.species_masses) {
            assert!((after-initial).abs() < 1e-14);
        }
    }
}

#[test]
fn final_cancellation_and_native_input_refusal_preserve_all_accepted_owners() {
    let specimen = body(4000.0);
    assert!(MaterialResolvedPiston::new(specimen.clone(), f64::NAN).is_err());
    let mut wall = MaterialResolvedPiston::new(specimen, 0.08).unwrap();
    let (mut gas, mut reference) = (frozen(), frozen());
    let mut baseline = wall.clone();
    let g = StepControl::default();
    wall.step(&mut gas, DT, g, load(), controls(), || false).unwrap();
    baseline.step(&mut reference, DT, g, load(), controls(), || false).unwrap();
    let original = wall.clone();
    let states = gas.states().to_vec();
    let mesh = gas.mesh().clone();
    let time = gas.time();
    let mut polls = 0;
    let expected = baseline.step_suspended(&mut reference, DT, g, load(), suspension(), controls(),
        || { polls += 1; false }).unwrap();
    let mut visits = 0;
    let cancelled = wall.step_suspended(&mut gas, DT, g, load(), suspension(), controls(),
        || { visits += 1; visits == polls });
    assert_eq!(cancelled, Err(PistonError::Cancelled));
    assert_eq!(wall, original);
    assert_eq!(gas.states(), states);
    assert_eq!(gas.mesh(), &mesh);
    assert_eq!(gas.time(), time);
    let invalid = PistonLoad { exterior_pressure_pa: f64::NAN, ..load() };
    assert!(wall.step(&mut gas, DT, g, invalid, controls(), || false).is_err());
    assert_eq!(wall, original);
    assert_eq!(gas.states(), states);
    assert_eq!(gas.mesh(), &mesh);
    assert_eq!(gas.time(), time);
    let actual = wall.step_suspended(&mut gas, DT, g, load(), suspension(), controls(), || false).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(wall, baseline);
    assert_eq!(gas.states(), reference.states());
    assert_eq!(gas.mesh(), reference.mesh());
    assert_eq!(gas.time(), reference.time());
}
