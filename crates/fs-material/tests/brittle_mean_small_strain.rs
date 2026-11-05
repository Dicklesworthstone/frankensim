//! Native mean-force assembly retains the original brittle endpoint history.
use fs_material::brittle::{
    BrittleCrackBand, BrittleCrackBandParameters, BrittleCrackBandSettings, BrittleError,
};
use fs_material::elastic::IsotropicElastic;
use fs_material::tensor::contract;
use fs_material::{SmallStrainLaw, Voigt};

const E: f64 = 3e9;

fn law(normal: [f64; 3]) -> BrittleCrackBand {
    BrittleCrackBand::new(
        IsotropicElastic::new(E, 0.25, 0.05).unwrap(),
        BrittleCrackBandParameters {
            tensile_strength_pa: 3e6,
            shear_strength_pa: 6e6,
            mode_i_fracture_energy_j_m2: 120.0,
            mode_ii_fracture_energy_j_m2: 360.0,
            bk_exponent: 1.6,
        },
        BrittleCrackBandSettings {
            band_width_m: 0.004,
            normal,
            maximum_strain: 0.05,
            energy_tolerance_j_m3: 1e-6,
        },
    )
    .unwrap()
}

fn close(a: f64, b: f64, relative: f64, floor: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(floor),
        "{a:e} != {b:e}"
    );
}

#[test]
fn intact_force_is_the_mean_and_tangent_is_its_exact_half_stiffness() {
    let band = law([1.0, 2.0, 3.0]);
    let adapter = band.mean_small_strain_law();
    assert_eq!(adapter.initial_state(), band.initial_state());
    let old_strain = [0.0001, 0.00008, 0.00004, 0.00002, 0.0, -0.00001];
    let end = [0.0002, 0.00005, 0.00003, -0.00001, 0.00003, 0.00002];
    let old = band
        .trial(&old_strain, &band.initial_state())
        .unwrap()
        .state;
    let saved = old;
    let (path, tangent) = adapter.trial(&end, &old).unwrap();
    let mean_strain = core::array::from_fn(|i| 0.5 * (old_strain[i] + end[i]));
    let expected = band.elastic().stress(&mean_strain, &());
    let stress = adapter.stress(&end, &old);
    assert_eq!(stress, path.mean_stress_pa);
    assert_ne!(stress, path.endpoint.stress_pa);
    for i in 0..6 {
        close(stress[i], expected[i], 3e-12, 1.0);
    }
    let tangent = tangent.unwrap();
    assert_eq!(adapter.tangent(&end, &old), tangent);
    let elastic = band.elastic().stiffness();
    for i in 0..6 {
        for j in 0..6 {
            close(tangent[i][j], 0.5 * elastic[i][j], 3e-12, E);
        }
    }
    assert_eq!(adapter.update_state(&end, &old), path.endpoint.state);
    assert_eq!(path.endpoint, band.trial(&end, &old).unwrap());
    assert_eq!(old, saved);
}

#[test]
fn force_work_and_accepted_history_survive_loading_unloading_rupture_and_reclosure() {
    let band = law([1.0, 0.0, 0.0]);
    let adapter = band.mean_small_strain_law();
    let mut old = adapter.initial_state();
    let mut work = 0.0;
    for end in [
        [0.012, 0.0, 0.0, 0.002, 0.0, 0.0],
        [0.004, 0.0, 0.0, 0.0005, 0.0, 0.0],
        [0.01, 0.0002, -0.0001, -0.002, 0.0, 0.001],
        [0.028, 0.0, 0.0, 0.002, 0.0, 0.0],
        [-0.006, 0.0015, 0.0015, 0.0, 0.0, 0.0],
        [0.03, 0.0, 0.0, 0.003, 0.0, 0.0],
    ] {
        let saved = old;
        let (path, _) = adapter.trial(&end, &old).unwrap();
        let force = adapter.stress(&end, &old);
        let increment: Voigt = core::array::from_fn(|i| end[i] - old.strain()[i]);
        assert_eq!(force, path.mean_stress_pa);
        close(
            contract(&force, &increment),
            path.endpoint.energy.work,
            2e-11,
            1.0,
        );
        let accepted = adapter.update_state(&end, &old);
        assert_eq!(accepted, path.endpoint.state);
        assert_eq!(accepted, band.trial(&end, &old).unwrap().state);
        assert!(accepted.damage() >= old.damage());
        assert_eq!(old, saved);
        work += contract(&force, &increment);
        old = accepted;
    }
    assert!(old.failed());
    close(
        work,
        old.stored_energy_j_m3() + old.dissipated_j_m3(),
        2e-11,
        1.0,
    );
}

#[test]
fn mixed_mode_jacobian_differentiates_mean_force_at_fixed_original_history() {
    let band = law([1.0, 2.0, 3.0]);
    let adapter = band.mean_small_strain_law();
    let end = [0.006, -0.001, 0.0005, 0.002, 0.0008, -0.001];
    let virgin = band
        .trial(
            &[0.0002, 0.0, 0.0, 0.00004, 0.0, 0.0],
            &band.initial_state(),
        )
        .unwrap()
        .state;
    let retained = adapter.update_state(&end, &virgin);
    assert!(retained.damage() > 0.0 && retained.damage() < 1.0);
    for (old, target) in [
        (virgin, end),
        (retained, end.map(|x| 1.25 * x)),
        (retained, end.map(|x| 0.65 * x)),
    ] {
        let saved = old;
        let tangent = adapter.tangent(&target, &old);
        let h = 1e-8;
        for j in 0..6 {
            let mut plus = target;
            let mut minus = target;
            plus[j] += h;
            minus[j] -= h;
            let plus = adapter.stress(&plus, &old);
            let minus = adapter.stress(&minus, &old);
            for i in 0..6 {
                close(tangent[i][j], (plus[i] - minus[i]) / (2.0 * h), 3e-6, 3e3);
            }
        }
        assert_eq!(old, saved);
    }
}

#[test]
fn nonsmooth_newton_starts_and_source_refusals_never_fabricate_a_tangent_or_state() {
    let band = law([1.0, 0.0, 0.0]);
    let adapter = band.mean_small_strain_law();
    let virgin = adapter.initial_state();
    let damaged = adapter.update_state(&[0.006, 0.0, 0.0, 0.0, 0.0, 0.0], &virgin);
    assert!(damaged.damage() > 0.0 && damaged.damage() < 1.0);
    for old in [virgin, damaged] {
        let (path, tangent) = adapter.trial(old.strain(), &old).unwrap();
        assert!(tangent.is_none());
        assert!(
            adapter
                .tangent(old.strain(), &old)
                .iter()
                .flatten()
                .all(|x| x.is_nan())
        );
        assert!(
            adapter
                .stress(old.strain(), &old)
                .iter()
                .all(|x| x.is_finite())
        );
        assert_eq!(path.path_work_j_m3, 0.0);
        assert_eq!(adapter.update_state(old.strain(), &old), old);
    }
    let clone = band.clone();
    let foreign = clone.initial_state();
    assert!(matches!(
        adapter.trial(&[0.0; 6], &foreign),
        Err(BrittleError::StateMismatch)
    ));
    for (strain, old) in [
        ([0.0; 6], foreign),
        ([0.036, 0.036, 0.0, 0.0, 0.0, 0.0], virgin),
    ] {
        assert!(adapter.trial(&strain, &old).is_err());
        assert!(adapter.stress(&strain, &old).iter().all(|x| x.is_nan()));
        assert!(
            adapter
                .tangent(&strain, &old)
                .iter()
                .flatten()
                .all(|x| x.is_nan())
        );
        let marker = adapter.update_state(&strain, &old);
        assert!(!marker.is_valid());
        assert!(adapter.trial(&[0.0; 6], &marker).is_err());
    }
    assert_eq!(virgin, band.initial_state());
}
