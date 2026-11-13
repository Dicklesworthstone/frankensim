//! Measured-curve shape, permanent axial strain and native energy survive reduction.
use fs_material::fiber::Uniaxial;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::reduction::{UniaxialJ2, UniaxialTabulated};
use fs_material::plastic::tabulated::{
    HardeningPoint, TabulatedError, TabulatedHardening, TabulatedJ2,
};
use fs_material::IsotropicElastic;

fn curve(points: &[(f64, f64)]) -> TabulatedHardening {
    TabulatedHardening::new(
        points
            .iter()
            .map(|&(plastic_strain, yield_stress_pa)| HardeningPoint {
                plastic_strain,
                yield_stress_pa,
            })
            .collect(),
    )
    .unwrap()
}
fn metal(points: &[(f64, f64)], kinematic: f64) -> UniaxialTabulated {
    let base = J2Plasticity::new_mixed(
        IsotropicElastic::new(1000.0, 0.3, 0.05).unwrap(),
        points[0].1,
        0.0,
        kinematic,
    )
    .unwrap();
    UniaxialTabulated::new(TabulatedJ2::new(base, curve(points)).unwrap()).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, allowance {tolerance:e}"
    );
}
fn gates() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-10,
        yield_pa: 1e-9,
    }
}

#[test]
fn each_measured_segment_drives_exact_axial_stress_tangent_and_native_energy() {
    let law = metal(&[(0.0, 1.0), (0.002, 1.2), (0.006, 3.2), (0.02, 4.6)], 0.0);
    let old = law.initial_state();
    for (alpha, yield_stress, h) in [(0.001, 1.1, 100.0), (0.004, 2.2, 500.0), (0.01, 3.6, 100.0)] {
        let strain = alpha + yield_stress / 1000.0;
        for sign in [-1.0, 1.0] {
            let response = law.response(sign * strain, &old).unwrap();
            near(response.state.equivalent_plastic_strain(), alpha, 1e-16);
            near(response.state.plastic_strain(), sign * alpha, 1e-16);
            near(response.stress_pa, sign * yield_stress, 1e-12);
            near(response.tangent_pa, 1000.0 * h / (1000.0 + h), 1e-11);
            let native = law
                .backbone()
                .response(&response.strain, &old.tensor_history())
                .unwrap();
            for i in 1..6 {
                near(native.stress_pa[i], 0.0, 1e-12);
            }
            near(native.state.alpha, alpha, 1e-16);
            let energy = law
                .advance_with_energy(0.0, &old, sign * strain, gates())
                .unwrap();
            near(energy.energy.physical_dissipation_j_m3, alpha, 1e-14);
            near(energy.energy.residual_j_m3, 0.0, 1e-13);
            assert!(energy.energy.final_state.hardening_j_m3 > 0.0);
            let unloaded = law
                .response(response.state.plastic_strain(), &response.state)
                .unwrap();
            near(unloaded.stress_pa, 0.0, 1e-14);
            assert_eq!(unloaded.state, response.state);
        }
    }
    assert_eq!(law.initial_state(), old);
}

#[test]
fn linear_table_matches_mixed_j2_through_loading_unloading_and_reversal() {
    let tabulated = metal(&[(0.0, 1.0), (0.03, 7.0)], 80.0);
    let original = UniaxialJ2::new(
        J2Plasticity::new_mixed(
            IsotropicElastic::new(1000.0, 0.3, 0.05).unwrap(),
            1.0,
            200.0,
            80.0,
        )
        .unwrap(),
    )
    .unwrap();
    let mut a = tabulated.initial_state();
    let mut b = original.initial_state();
    for strain in [0.0005, 0.003, 0.008, 0.005, 0.0, -0.006, -0.002, 0.005] {
        let x = tabulated.response(strain, &a).unwrap();
        let y = original.response(strain, &b);
        near(x.stress_pa, y.stress_pa, 2e-12);
        near(x.tangent_pa, y.tangent_pa, 2e-10);
        near(x.state.plastic_strain(), y.state.plastic_strain(), 2e-15);
        near(
            x.state.equivalent_plastic_strain(),
            y.state.equivalent_plastic_strain(),
            2e-15,
        );
        a = x.state;
        b = y.state;
    }
}

#[test]
fn measured_limit_is_retained_during_elastic_unloading_and_refuses_further_flow() {
    let law = metal(&[(0.0, 1.0), (0.002, 1.2)], 0.0);
    let virgin = law.initial_state();
    let loaded = law.response(0.0032, &virgin).unwrap();
    near(loaded.state.equivalent_plastic_strain(), 0.002, 1e-17);
    assert!(matches!(
        law.response(0.0033, &loaded.state),
        Err(TabulatedError::Domain { .. })
    ));
    assert!(matches!(
        law.response(0.0007, &loaded.state),
        Err(TabulatedError::Domain { .. })
    ));
    let unloaded = law.response(0.002, &loaded.state).unwrap();
    near(unloaded.stress_pa, 0.0, 1e-14);
    assert_eq!(unloaded.state, loaded.state);
    assert!(
        law.advance_with_energy(0.0032, &loaded.state, 0.0033, gates())
            .is_err()
    );
    assert!(law.stress(0.0033, &loaded.state).is_nan());
    assert!(law.response(f64::NAN, &virgin).is_err());
    assert_eq!(law.initial_state(), virgin);
}

#[test]
fn an_accepted_state_re_evaluated_at_its_own_strain_is_accepted_again() {
    // Steel-scale coefficients (those of fs-solid's frame_ductile fibers): E
    // 200 GPa, Prager 18 GPa, a table rising 2 GPa from 200 MPa. At the same
    // strain the re-evaluated trial stress is assembled apart from the yield
    // value and can land an ulp above it; where its exact return is below the
    // resolution of the accumulated plastic strain, the accepted history is
    // the nearest representable one.
    let (e, kinematic) = (200e9, 18e9);
    let base = J2Plasticity::new_mixed(
        IsotropicElastic::new(e, 0.3, 0.04).unwrap(),
        200e6,
        0.0,
        kinematic,
    )
    .unwrap();
    let law = UniaxialTabulated::new(
        TabulatedJ2::new(base, curve(&[(0.0, 200e6), (0.05, 300e6)])).unwrap(),
    )
    .unwrap();
    let virgin = law.initial_state();
    let mut kept = 0;
    for k in 1..=400 {
        let strain = 0.001 + 0.005 * f64::from(k) / 400.0;
        for sign in [-1.0, 1.0] {
            let loaded = law.response(sign * strain, &virgin).unwrap();
            let alpha = loaded.state.equivalent_plastic_strain();
            let plastic = loaded.state.plastic_strain();
            let again = law.response(sign * strain, &loaded.state).unwrap();
            near(
                again.stress_pa,
                loaded.stress_pa,
                32.0 * f64::EPSILON * loaded.stress_pa.abs(),
            );
            let advance = again.state.equivalent_plastic_strain() - alpha;
            assert!((0.0..=1e-18).contains(&advance), "{advance:e}");
            let trial = e * (sign * strain - plastic) - kinematic * plastic;
            let radius = law.backbone().hardening().yield_stress_pa(alpha).unwrap();
            if trial.abs() > radius && advance == 0.0 {
                assert_eq!(again.state, loaded.state);
                kept += 1;
            }
        }
    }
    assert!(
        kept > 0,
        "no re-evaluated trial landed above its yield value"
    );
}
