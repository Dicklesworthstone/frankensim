//! Independent work identities for the native tensor-Voigt return map.
use fs_material::{SmallStrainLaw, Voigt};
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::{J2Plasticity, J2State};
use fs_material::plastic::energy::{EnergyError, EnergyTolerances};

const TOL: EnergyTolerances = EnergyTolerances { energy_j_m3: 1e-11, yield_pa: 1e-9 };

fn law(yield_pa: f64, hardening: f64) -> J2Plasticity {
    J2Plasticity::new(IsotropicElastic::new(1200.0, 0.2, 0.04).unwrap(),
        yield_pa, hardening).unwrap()
}
fn shear(xy: f64) -> Voigt { [0.0, 0.0, 0.0, xy, 0.0, 0.0] }
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-12 * b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn pure_shear_work_has_tensor_factor_two_and_stores_energetic_hardening() {
    let material = law(3.0, 100.0);
    let step = material.advance_with_energy(&[0.0; 6], &material.initial_state(),
        &shear(0.01), TOL).unwrap();
    // E=1200, nu=.2 -> mu=500. A pure-shear yield state has q=sqrt(3)*tau.
    let alpha = (3.0_f64.sqrt() * 10.0 - 3.0) / 1600.0;
    let tau = (3.0 + 100.0 * alpha) / 3.0_f64.sqrt();
    near(step.state.alpha, alpha);
    near(step.stress_pa[3], tau);
    near(step.energy.mechanical_work_j_m3, 2.0 * tau * 0.01);
    near(step.energy.final_state.elastic_j_m3, tau * tau / 1000.0);
    near(step.energy.final_state.hardening_j_m3, 50.0 * alpha * alpha);
    near(step.energy.physical_dissipation_j_m3, 3.0 * alpha);
    let elastic_shear = tau / 1000.0;
    near(step.energy.numerical_dissipation_j_m3,
        1000.0 * elastic_shear * elastic_shear + 50.0 * alpha * alpha);
    near(step.energy.residual_j_m3, 0.0);
    assert!(step.energy.final_state.hardening_j_m3 > 0.0);
    assert_eq!(step.state, material.update_state(&shear(0.01), &material.initial_state()));
}

#[test]
fn elastic_loading_and_unloading_have_no_physical_heat() {
    let material = law(3.0, 100.0);
    let elastic = material.advance_with_energy(&[0.0; 6], &material.initial_state(),
        &shear(0.0005), TOL).unwrap();
    near(elastic.stress_pa[3], 0.5);
    near(elastic.energy.mechanical_work_j_m3, 0.0005);
    near(elastic.energy.final_state.elastic_j_m3, 0.00025);
    near(elastic.energy.numerical_dissipation_j_m3, 0.00025);
    assert_eq!(elastic.energy.physical_dissipation_j_m3, 0.0);
    let loaded = material.advance_with_energy(&[0.0; 6], &material.initial_state(),
        &shear(0.01), TOL).unwrap();
    let unloaded = material.advance_with_energy(&loaded.strain, &loaded.state,
        &shear(0.009), TOL).unwrap();
    assert_eq!(unloaded.state, loaded.state);
    assert_eq!(unloaded.energy.physical_dissipation_j_m3, 0.0);
    assert!(unloaded.energy.mechanical_work_j_m3 < 0.0);
    near(unloaded.energy.initial.hardening_j_m3, unloaded.energy.final_state.hardening_j_m3);
    near(unloaded.energy.residual_j_m3, 0.0);
}

#[test]
fn reversal_preserves_accumulated_plastic_history_and_the_full_work_identity() {
    let material = law(3.0, 100.0);
    let loaded = material.advance_with_energy(&[0.0; 6], &material.initial_state(),
        &shear(0.01), TOL).unwrap();
    let reversed = material.advance_with_energy(&loaded.strain, &loaded.state,
        &shear(-0.012), TOL).unwrap();
    assert!(reversed.state.alpha > loaded.state.alpha);
    assert!(reversed.stress_pa[3] < 0.0);
    assert!(reversed.energy.physical_dissipation_j_m3 > 0.0);
    near(reversed.energy.initial.total_j_m3(), loaded.energy.final_state.total_j_m3());
    near(reversed.energy.mechanical_work_j_m3,
        reversed.energy.stored_change_j_m3()
            + reversed.energy.physical_dissipation_j_m3
            + reversed.energy.numerical_dissipation_j_m3);
}

#[test]
fn lowering_only_yield_at_fixed_total_strain_releases_stored_energy() {
    let cold = law(3.0, 100.0);
    let old = cold.advance_with_energy(&[0.0; 6], &cold.initial_state(),
        &shear(0.01), TOL).unwrap();
    let hot = law(1.0, 100.0);
    assert!(hot.validate_committed_history(&old.strain, &old.state, TOL.yield_pa).is_err());
    let step = hot.advance_with_energy(&old.strain, &old.state, &old.strain, TOL).unwrap();
    assert_eq!(step.energy.mechanical_work_j_m3, 0.0);
    assert!(step.energy.stored_change_j_m3() < 0.0);
    assert!(step.energy.physical_dissipation_j_m3 > 0.0);
    near(-step.energy.stored_change_j_m3(),
        step.energy.physical_dissipation_j_m3 + step.energy.numerical_dissipation_j_m3);
}

#[test]
fn nonfinite_traceful_unreachable_or_large_strain_history_refuses() {
    let material = law(3.0, 100.0);
    for state in [
        J2State { alpha: f64::NAN, plastic_strain: [0.0; 6] },
        J2State { alpha: -0.1, plastic_strain: [0.0; 6] },
        J2State { alpha: 0.01, plastic_strain: [0.001, 0.0, 0.0, 0.0, 0.0, 0.0] },
        J2State { alpha: 0.0, plastic_strain: shear(0.001) },
    ] {
        assert!(matches!(material.advance_with_energy(&[0.0; 6], &state,
            &shear(0.01), TOL), Err(EnergyError::Invalid(_))));
    }
    for strain in [shear(f64::NAN), shear(0.04)] {
        assert!(material.advance_with_energy(&[0.0; 6], &material.initial_state(),
            &strain, TOL).is_err());
    }
    assert!(material.validate_committed_history(&shear(0.01),
        &material.initial_state(), TOL.yield_pa).is_err());
}
