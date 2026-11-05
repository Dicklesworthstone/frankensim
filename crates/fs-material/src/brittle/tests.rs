use super::*;
use crate::tensor::{rotate, rotation};

fn parameters() -> BrittleCrackBandParameters {
    BrittleCrackBandParameters {
        tensile_strength_pa: 3.0e6,
        shear_strength_pa: 6.0e6,
        mode_i_fracture_energy_j_m2: 120.0,
        mode_ii_fracture_energy_j_m2: 360.0,
        bk_exponent: 1.6,
    }
}
fn settings(width: f64, normal: [f64; 3]) -> BrittleCrackBandSettings {
    BrittleCrackBandSettings {
        band_width_m: width, normal, maximum_strain: 0.05, energy_tolerance_j_m3: 1e-6,
    }
}
fn law(width: f64, normal: [f64; 3]) -> BrittleCrackBand {
    BrittleCrackBand::new(
        IsotropicElastic::new(3.0e9, 0.25, 0.05).unwrap(), parameters(), settings(width, normal),
    ).unwrap()
}
fn close(a: f64, b: f64, relative: f64, floor: f64) {
    assert!((a-b).abs() <= relative*a.abs().max(b.abs()).max(floor), "{a:e} != {b:e}");
}
fn close_tensor(a: &Voigt, b: &Voigt, relative: f64, floor: f64) {
    for i in 0..6 { close(a[i], b[i], relative, floor); }
}

#[test]
fn g0_intact_projection_recovers_isotropic_stress_energy_and_tangent_including_auxetics() {
    let strain = [2e-6, -1e-6, 1.5e-6, 0.3e-6, -0.2e-6, 0.4e-6];
    for poisson in [-0.8, -0.2, 0.0, 0.25, 0.49] {
        for normal in [[1.0, 0.0, 0.0], [1.0, 2.0, -3.0]] {
            let elastic = IsotropicElastic::new(3e9, poisson, 0.05).unwrap();
            let band = BrittleCrackBand::new(elastic.clone(), parameters(), settings(0.004, normal)).unwrap();
            let step = band.trial(&strain, &band.initial_state()).unwrap();
            let stress = elastic.stress(&strain, &());
            close_tensor(&step.stress_pa, &stress, 1e-12, 1e4);
            close(step.state.stored_energy_j_m3(), 0.5*contract(&strain, &stress), 1e-12, 1.0);
            assert_eq!(step.state.damage(), 0.0);
            assert_eq!(step.energy.dissipated, 0.0);
            for (actual, expected) in step.branch_tangent_pa.iter().zip(elastic.stiffness()) {
                close_tensor(actual, &expected, 1e-12, 3e9);
            }
        }
    }
}

#[test]
fn g1_free_sided_compression_does_not_initiate_damage_through_poisson_expansion() {
    for poisson in [-0.4, 0.25, 0.49] {
        let elastic = IsotropicElastic::new(3e9, poisson, 0.05).unwrap();
        let band = BrittleCrackBand::new(elastic.clone(), parameters(), settings(0.004, [1.0, 0.0, 0.0])).unwrap();
        let strain = [-0.02, poisson*0.02, poisson*0.02, 0.0, 0.0, 0.0];
        let step = band.trial(&strain, &band.initial_state()).unwrap();
        close_tensor(&step.stress_pa, &[-6e7, 0.0, 0.0, 0.0, 0.0, 0.0], 2e-12, 6e7);
        assert_eq!(step.state.damage(), 0.0);
        assert_eq!(step.energy.dissipated, 0.0);
        assert!(!step.state.failed());
    }
}

#[test]
fn g1_pure_mode_fracture_work_is_original_toughness_per_band_width() {
    for width in [0.004, 0.008, 0.012] {
        let band = law(width, [1.0, 0.0, 0.0]);
        for shear in [false, true] {
            let (strength, toughness) = if shear { (6e6, 360.0) } else { (3e6, 120.0) };
            let failure_mode_strain = 2.0*toughness/(strength*width);
            let mut state = band.initial_state();
            let mut work = 0.0;
            let mut dissipated = 0.0;
            for i in 1..=64 {
                let q = 1.05*failure_mode_strain*f64::from(i)/64.0;
                let mut strain = [0.0; 6];
                if shear { strain[3] = 0.5*q; } else { strain[0] = q; }
                let step = band.trial(&strain, &state).unwrap();
                work += step.energy.work;
                dissipated += step.energy.dissipated;
                assert!(step.state.damage() >= state.damage());
                state = step.state;
            }
            assert!(state.failed());
            assert_eq!(state.damage(), 1.0);
            close(work*width, toughness, 2e-12, 1.0);
            close(dissipated*width, toughness, 2e-12, 1.0);
            close(state.dissipated_j_m3()*width, toughness, 2e-12, 1.0);
            close(state.stored_energy_j_m3(), 0.0, 1e-10, 1.0);
        }
    }
}

#[test]
fn g1_declared_plane_normal_and_shear_strengths_are_actual_peak_tractions() {
    let band = law(0.004, [1.0, 0.0, 0.0]);
    let (lambda, mu) = band.elastic().lame();
    for (index, onset, expected) in [(0, 3e6/(lambda+2.0*mu), 3e6), (3, 6e6/(2.0*mu), 6e6)] {
        let mut strain = [0.0; 6];
        strain[index] = onset;
        let step = band.trial(&strain, &band.initial_state()).unwrap();
        close(step.stress_pa[index], expected, 1e-12, 1.0);
    }
}

#[test]
fn g3_unloading_keeps_damage_and_failed_band_recloses_without_healing() {
    let band = law(0.004, [1.0, 0.0, 0.0]);
    let initial = band.initial_state();
    let damaged = band.trial(&[0.012, 0.0, 0.0, 0.0, 0.0, 0.0], &initial).unwrap();
    assert!(damaged.state.damage() > 0.0 && damaged.state.damage() < 1.0);
    let unloaded = band.trial(&[0.004, 0.0, 0.0, 0.0, 0.0, 0.0], &damaged.state).unwrap();
    assert_eq!(unloaded.state.damage(), damaged.state.damage());
    assert_eq!(unloaded.state.dissipated_j_m3(), damaged.state.dissipated_j_m3());
    assert!(unloaded.energy.work < 0.0);
    close(unloaded.stress_pa[0], damaged.stress_pa[0]/3.0, 1e-12, 1.0);
    let reloaded = band.trial(&[0.01, 0.0, 0.0, 0.0, 0.0, 0.0], &unloaded.state).unwrap();
    assert_eq!(reloaded.state.damage(), damaged.state.damage());
    let failed = band.trial(&[0.022, 0.0, 0.0, 0.0, 0.0, 0.0], &reloaded.state).unwrap();
    assert!(failed.state.failed());
    let strain = [-0.01, 0.0025, 0.0025, 0.0, 0.0, 0.0];
    let compressed = band.trial(&strain, &failed.state).unwrap();
    assert!(compressed.state.failed());
    assert_eq!(compressed.state.damage(), 1.0);
    assert_eq!(compressed.energy.dissipated, 0.0);
    close_tensor(&compressed.stress_pa, &band.elastic().stress(&strain, &()), 1e-12, 3e7);
    let reopened = band.trial(&[0.03, 0.0, 0.0, 0.005, 0.0, 0.0], &compressed.state).unwrap();
    close_tensor(&reopened.stress_pa, &[0.0; 6], 1e-12, 1.0);
    assert!(reopened.state.failed());
    assert_eq!(reopened.energy.dissipated, 0.0);
    assert_eq!(initial.damage(), 0.0);
}

#[test]
fn g3_failed_plane_keeps_only_in_plane_elasticity_under_opening_and_shear() {
    let band = law(0.004, [1.0, 0.0, 0.0]);
    let failed = band.trial(&[0.022, 0.0, 0.0, 0.0, 0.0, 0.0], &band.initial_state()).unwrap();
    let strain = [0.025, 0.003, -0.001, 0.002, 0.001, -0.003];
    let step = band.trial(&strain, &failed.state).unwrap();
    let e = band.elastic().youngs;
    let nu = band.elastic().poisson;
    let d = e/(1.0-nu*nu);
    let expected = [0.0, d*(strain[1]+nu*strain[2]), d*(strain[2]+nu*strain[1]),
        0.0, e/(1.0+nu)*strain[4], 0.0];
    close_tensor(&step.stress_pa, &expected, 1e-12, 1.0);
    assert!(step.state.stored_energy_j_m3() > 0.0);
}

#[test]
fn g3_rotating_plane_and_tensor_history_rotates_stress_and_preserves_energy() {
    let axis = normalize([0.3, -0.4, 0.5]).unwrap();
    let rotation = rotation(axis, 0.73);
    let normal = [rotation[0][0], rotation[1][0], rotation[2][0]];
    let a = law(0.004, [1.0, 0.0, 0.0]);
    let b = law(0.004, normal);
    let mut state_a = a.initial_state();
    let mut state_b = b.initial_state();
    for strain in [
        [0.0001, 0.0, 0.0, 0.0002, 0.0001, 0.0003],
        [0.006, -0.001, 0.0005, 0.002, 0.0008, -0.001],
        [0.003, 0.001, 0.0005, -0.003, 0.0005, 0.002],
        [-0.005, 0.001, 0.001, 0.001, 0.0, 0.0],
        [0.028, 0.0, 0.0, -0.008, 0.0, 0.005],
    ] {
        let step_a = a.trial(&strain, &state_a).unwrap();
        let rotated_strain = rotate(&strain, &rotation);
        let step_b = b.trial(&rotated_strain, &state_b).unwrap();
        close_tensor(&step_b.stress_pa, &rotate(&step_a.stress_pa, &rotation), 2e-11, 1e6);
        close(step_a.state.damage(), step_b.state.damage(), 2e-11, 1.0);
        close(step_a.state.dissipated_j_m3(), step_b.state.dissipated_j_m3(), 2e-11, 1.0);
        close(step_a.energy.work, step_b.energy.work, 2e-11, 1.0);
        close(step_a.state.stored_energy_j_m3(), step_b.state.stored_energy_j_m3(), 2e-11, 1.0);
        state_a = step_a.state;
        state_b = step_b.state;
    }
}

fn check_tangent(band: &BrittleCrackBand, old: &BrittleState<'_>, strain: Voigt) {
    let step = band.trial(&strain, old).unwrap();
    let tangent = step.tangent_pa.expect("smooth test point");
    let perturbation = 1e-8;
    for j in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        plus[j] += perturbation;
        minus[j] -= perturbation;
        let sp = band.trial(&plus, old).unwrap().stress_pa;
        let sm = band.trial(&minus, old).unwrap().stress_pa;
        for i in 0..6 {
            close(tangent[i][j], (sp[i]-sm[i])/(2.0*perturbation), 2e-6, 3e3);
        }
    }
}

#[test]
fn g1_consistent_tensor_shear_tangent_covers_mixed_onset_history_loading_and_unloading() {
    let band = law(0.004, [1.0, 2.0, 3.0]);
    let initial = band.initial_state();
    let strain = [0.006, -0.001, 0.0005, 0.002, 0.0008, -0.001];
    check_tangent(&band, &initial, strain);
    let accepted = band.trial(&strain, &initial).unwrap();
    assert!(accepted.state.damage() > 0.0 && accepted.state.damage() < 1.0);
    check_tangent(&band, &accepted.state, strain.map(|x| 1.25*x));
    check_tangent(&band, &accepted.state, strain.map(|x| 0.65*x));
    assert!(!band.small_strain_law().admissibility().tangent_symmetric);
}

#[test]
fn g0_invalid_width_domain_plane_and_history_refuse_without_mutating_state() {
    let band = law(0.004, [1.0, 0.0, 0.0]);
    let initial = band.initial_state();
    let identical = band.clone();
    assert_eq!(identical.trial(&[0.0; 6], &initial), Err(BrittleError::StateMismatch));
    assert!(matches!(band.trial(&[0.051, 0.0, 0.0, 0.0, 0.0, 0.0], &initial),
        Err(BrittleError::StrainLimit { .. })));
    assert!(band.trial(&[f64::NAN; 6], &initial).is_err());
    for width in [0.0, -1.0, f64::INFINITY, 0.2] {
        assert!(BrittleCrackBand::new(band.elastic().clone(), parameters(), settings(width, [1.0, 0.0, 0.0])).is_err());
    }
    for normal in [[0.0; 3], [f64::NAN, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0]] {
        assert!(BrittleCrackBand::new(band.elastic().clone(), parameters(), settings(0.004, normal)).is_err());
    }
    let too_wide = BrittleCrackBandSettings { maximum_strain: 0.06, ..settings(0.004, [1.0, 0.0, 0.0]) };
    assert!(BrittleCrackBand::new(band.elastic().clone(), parameters(), too_wide).is_err());
    let borrowed = band.small_strain_law();
    assert_eq!(borrowed.card().discrepancy_rel, f64::INFINITY);
    assert!(borrowed.card().calibration.is_none());
    let invalid = [0.1, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert!(borrowed.stress(&invalid, &initial).iter().all(|x| x.is_nan()));
    assert!(borrowed.tangent(&invalid, &initial).iter().flatten().all(|x| x.is_nan()));
    let refused = borrowed.update_state(&invalid, &initial);
    assert!(!refused.is_valid());
    assert!(refused.damage().is_nan());
    assert!(band.trial(&[0.0; 6], &refused).is_err());
    assert_eq!(initial, band.initial_state());
}
