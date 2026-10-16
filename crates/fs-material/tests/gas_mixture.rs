//! Numerical regressions of the declared ideal-gas/Wilke/WMS models, not
//! experimental qualification of synthetic species or arbitrary atmospheres.
use fs_material::gas::{ConductivityModel, GasSpec, GasState, saturation_pressure_water_pa};
use fs_material::gas_mixture::{MAX_GAS_MIXTURE_COMPONENTS, mix_gas_states};
use fs_qty::semantic::{Composition, CompositionBasis};

fn composition(basis: CompositionBasis, fractions: Vec<f64>) -> Composition {
    Composition::new(basis, fractions).expect("valid test composition")
}

fn specs() -> [GasSpec; 3] {
    // Deliberately synthetic component transport data; no named-source claim.
    [
        GasSpec::dry_air_ussa1976(),
        GasSpec { molar_mass: 0.004, gamma: 5.0 / 3.0, sutherland_beta: 1.0e-6,
            sutherland_s: 80.0, conductivity: ConductivityModel::Eucken },
        GasSpec { molar_mass: 0.044, gamma: 1.3, sutherland_beta: 1.3e-6,
            sutherland_s: 220.0, conductivity: ConductivityModel::Eucken },
    ]
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 2.0e-12 * a.abs().max(b.abs()), "{a} != {b}");
}

fn same_properties(a: GasState, b: GasState) {
    for (a, b) in [
        (a.density, b.density), (a.sound_speed, b.sound_speed), (a.gamma, b.gamma),
        (a.specific_heat_cp, b.specific_heat_cp), (a.dynamic_viscosity, b.dynamic_viscosity),
        (a.thermal_conductivity, b.thermal_conductivity), (a.prandtl, b.prandtl),
        (a.characteristic_impedance, b.characteristic_impedance),
    ] { near(a, b); }
}

#[test]
fn pure_limit_uses_the_existing_evaluator_bitwise() {
    for spec in specs() {
        let x = composition(CompositionBasis::MoleFraction, vec![1.0]);
        let mixture = GasState::try_new_mixture(&[spec], &x, 300.0, 101325.0, None).unwrap();
        assert_eq!(mixture, GasState::try_new(&spec, 300.0, 101325.0).unwrap());
    }
}

#[test]
fn moist_air_matches_existing_binary_thermodynamics_and_transport() {
    for (t, p, rh) in [(293.15, 101325.0, 0.5), (310.15, 90000.0, 0.8)] {
        let xw = rh * saturation_pressure_water_pa(t).unwrap() / p;
        let x = composition(CompositionBasis::MoleFraction, vec![1.0 - xw, xw]);
        let mixed = GasState::try_new_mixture(
            &[GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()], &x, t, p, Some(1),
        ).unwrap();
        let binary = GasState::try_new_moist_air(t, p, rh).unwrap();
        same_properties(mixed, binary);
        assert_eq!(mixed.water_mole_fraction, binary.water_mole_fraction);
    }
}

#[test]
fn mass_and_mole_bases_describe_the_same_physical_mixture() {
    let specs = specs();
    let mole = composition(CompositionBasis::MoleFraction, vec![0.5, 0.25, 0.25]);
    let mass_sum: f64 = specs.iter().zip(mole.fractions()).map(|(s, x)| s.molar_mass * x).sum();
    let mass = composition(CompositionBasis::MassFraction,
        specs.iter().zip(mole.fractions()).map(|(s, x)| s.molar_mass * x / mass_sum).collect());
    same_properties(
        GasState::try_new_mixture(&specs, &mole, 320.0, 101325.0, None).unwrap(),
        GasState::try_new_mixture(&specs, &mass, 320.0, 101325.0, None).unwrap(),
    );
}

#[test]
fn pressure_changes_loading_and_absorption_not_ideal_gas_pitch() {
    let x = composition(CompositionBasis::MoleFraction, vec![0.5, 0.25, 0.25]);
    let a = GasState::try_new_mixture(&specs(), &x, 300.0, 50000.0, None).unwrap();
    let b = GasState::try_new_mixture(&specs(), &x, 300.0, 100000.0, None).unwrap();
    near(b.density, 2.0 * a.density);
    near(b.characteristic_impedance, 2.0 * a.characteristic_impedance);
    near(b.sound_speed, a.sound_speed);
    near(b.dynamic_viscosity, a.dynamic_viscosity);
    near(b.thermal_conductivity, a.thermal_conductivity);
    near(b.stokes_kirchhoff_absorption(6000.0), 0.5 * a.stokes_kirchhoff_absorption(6000.0));
}

#[test]
fn composition_changes_the_complete_medium_not_just_sound_speed() {
    let x = composition(CompositionBasis::MoleFraction, vec![0.5, 0.5]);
    let mixed = GasState::try_new_mixture(&specs()[..2], &x, 300.0, 101325.0, None).unwrap();
    let air = GasState::try_new(&specs()[0], 300.0, 101325.0).unwrap();
    assert!(mixed.sound_speed > air.sound_speed);
    assert!(mixed.density < air.density);
    assert_ne!(mixed.specific_heat_cp, air.specific_heat_cp);
    assert_ne!(mixed.thermal_conductivity, air.thermal_conductivity);
    assert_ne!(mixed.dynamic_viscosity, air.dynamic_viscosity);
    near(mixed.prandtl, mixed.dynamic_viscosity * mixed.specific_heat_cp / mixed.thermal_conductivity);
    near(mixed.sound_speed * mixed.sound_speed, mixed.gamma * mixed.pressure / mixed.density);
}

#[test]
fn input_permutation_preserves_physics_and_replay_is_exact() {
    let s = specs();
    let a = composition(CompositionBasis::MoleFraction, vec![0.5, 0.25, 0.25]);
    let b = composition(CompositionBasis::MoleFraction, vec![0.25, 0.5, 0.25]);
    let first = GasState::try_new_mixture(&s, &a, 300.0, 101325.0, None).unwrap();
    same_properties(first, GasState::try_new_mixture(&[s[2], s[0], s[1]], &b, 300.0, 101325.0, None).unwrap());
    assert_eq!(first, GasState::try_new_mixture(&s, &a, 300.0, 101325.0, None).unwrap());
}

#[test]
fn inconsistent_and_mismatched_states_fail_closed() {
    let a = GasState::try_new(&specs()[0], 300.0, 101325.0).unwrap();
    let x = composition(CompositionBasis::MoleFraction, vec![0.5, 0.5]);
    let mut corrupt = a;
    corrupt.density *= 2.0;
    assert!(mix_gas_states(&[a, corrupt], &x, None).is_err());
    let b = GasState::try_new(&specs()[1], 301.0, 101325.0).unwrap();
    assert!(mix_gas_states(&[a, b], &x, None).is_err());
    assert!(mix_gas_states(&[a, a], &x, Some(2)).is_err());
    assert!(mix_gas_states(&[a], &x, None).is_err());
    let oversized = vec![a; MAX_GAS_MIXTURE_COMPONENTS + 1];
    assert!(mix_gas_states(&oversized, &x, None).is_err());
}
