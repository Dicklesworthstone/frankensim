//! Original calibrated elastic strain bounds in checked constitutive consumers.

use fs_material::plastic::reduction::{PlaneStressReduction, ReductionConfig};
use fs_material::{IsotropicElastic, OrthotropicElastic, SmallStrainLaw, Voigt};

fn check_bound(law: &impl SmallStrainLaw<State = ()>, limit: f64) {
    let inside = [
        0.4 * limit,
        -0.2 * limit,
        0.1 * limit,
        0.2 * limit,
        0.0,
        0.0,
    ];
    let outside = [0.0, 0.0, 0.0, 0.8 * limit, 0.0, 0.0];
    assert!(law.admits_strain(&inside));
    assert!(law.admits_strain(&[limit, 0.0, 0.0, 0.0, 0.0, 0.0]));
    assert!(!law.admits_strain(&[limit.next_up(), 0.0, 0.0, 0.0, 0.0, 0.0]));
    assert!(!law.admits_strain(&outside)); // BOTH tensor shear entries matter.
    // Force the borrowed-law implementation rather than relying on autoderef.
    assert!(SmallStrainLaw::admits_strain(&law, &inside));
    assert!(!SmallStrainLaw::admits_strain(&law, &outside));
    let rotation = [[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
    assert!(law.admits_strain(&fs_material::tensor::rotate(&inside, &rotation)));
    assert!(!law.admits_strain(&fs_material::tensor::rotate(&outside, &rotation)));
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for component in 0..6 {
            let mut strain = [0.0; 6];
            strain[component] = bad;
            assert!(!law.admits_strain(&strain));
        }
    }
}

#[test]
fn calibrated_domains_retain_tensor_norm_rotation_and_extreme_scales() {
    for limit in [1e-280, 0.001, 1e280] {
        check_bound(&IsotropicElastic::new(1000.0, 0.3, limit).unwrap(), limit);
        check_bound(
            &OrthotropicElastic::new(
                [1000.0, 1300.0, 1700.0],
                [0.2, 0.1, 0.15],
                [350.0, 450.0, 550.0],
                limit,
            )
            .unwrap(),
            limit,
        );
    }
}

#[test]
fn checked_admission_preserves_raw_extrapolation_and_model_evidence() {
    let law = IsotropicElastic::new(1000.0, 0.3, 0.001).unwrap();
    let strain: Voigt = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert!(!law.admits_strain(&strain));
    let raw = law.stress(&strain, &());
    assert!(raw.iter().all(|v| v.is_finite()));
    assert!(
        law.tangent(&strain, &())
            .iter()
            .flatten()
            .all(|v| v.is_finite())
    );
    let evidence = fs_material::evidence_stress(&law, &strain, &());
    assert_eq!(evidence.value, raw);
    assert!(!evidence.model.in_domain);
}

#[test]
fn transverse_relaxation_cannot_extend_the_calibrated_elastic_domain() {
    let law = IsotropicElastic::new(1000.0, 0.3, 0.001).unwrap();
    let reduction = PlaneStressReduction::new(
        &law,
        ReductionConfig {
            max_iterations: 40,
            stress_tolerance_pa: 1e-12,
            strain_limit: 0.05,
        },
    )
    .unwrap();
    let old = reduction.initial_state();
    let accepted = reduction.evaluate([0.0005, 0.0005, 0.0], &old).unwrap();
    let ezz = -2.0 * 0.3 / 0.7 * 0.0005;
    assert!((accepted.state.strain[2] - ezz).abs() < 1e-15);
    assert!(law.admits_strain(&accepted.state.strain));
    for stress in &accepted.full_stress_pa[2..] {
        assert!(stress.abs() < 1e-12);
    }
    // The prescribed in-plane norm is admissible, but its actual relaxed
    // thickness strain would make the full tensor norm exceed 0.001.
    assert!(law.admits_strain(&[0.00065, 0.00065, 0.0, 0.0, 0.0, 0.0]));
    assert!(reduction.evaluate([0.00065, 0.00065, 0.0], &old).is_err());
    assert_eq!(old, reduction.initial_state());
    // A previous source law's strain is only a numerical transverse warm
    // start. It need not belong to the NEW law's calibrated strain domain.
    let mut warm_start = old;
    warm_start.strain = [0.002, 0.0, 0.002, 0.0, 0.0, 0.0];
    let before = warm_start.clone();
    let current = reduction
        .evaluate([0.0003, 0.0003, 0.0], &warm_start)
        .unwrap();
    assert!(law.admits_strain(&current.state.strain));
    assert_eq!(warm_start, before);
    warm_start.strain[0] = 0.06;
    assert!(reduction.evaluate([0.0; 3], &warm_start).is_err());
}
