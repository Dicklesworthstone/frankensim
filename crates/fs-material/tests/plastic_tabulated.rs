//! G0/G1/G3: independent measured-curve returns, cyclic histories and work.
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::tabulated::{
    HardeningPoint, TabulatedError, TabulatedHardening, TabulatedJ2,
};
use fs_material::plastic::{J2Plasticity, J2State};
use fs_material::tensor::{contract, trace};
use fs_material::{SmallStrainLaw, Voigt};

const TOL: EnergyTolerances = EnergyTolerances {
    energy_j_m3: 1e-10,
    yield_pa: 1e-9,
};

fn points() -> Vec<HardeningPoint> {
    [
        (0.0, 3.0),
        (0.002, 4.0),
        (0.006, 4.0),
        (0.012, 7.0),
        (0.04, 8.4),
    ]
    .map(|(plastic_strain, yield_stress_pa)| HardeningPoint {
        plastic_strain,
        yield_stress_pa,
    })
    .to_vec()
}

fn curve() -> TabulatedHardening {
    TabulatedHardening::new(points()).unwrap()
}

fn law(kinematic: f64) -> TabulatedJ2 {
    TabulatedJ2::new(
        J2Plasticity::new_mixed(
            IsotropicElastic::new(1200.0, 0.2, 0.05).unwrap(),
            3.0,
            0.0,
            kinematic,
        )
        .unwrap(),
        curve(),
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; allowance {tolerance:.3e}"
    );
}

fn dev(value: f64) -> Voigt {
    [value, -0.5 * value, -0.5 * value, 0.0, 0.0, 0.0]
}

#[test]
fn g1_table_interpolation_slopes_and_exact_piecewise_integrals() {
    let c = curve();
    for (a, y, slope) in [
        (0.0, 3.0, 500.0),
        (0.001, 3.5, 500.0),
        (0.002, 4.0, 0.0),
        (0.005, 4.0, 0.0),
        (0.006, 4.0, 500.0),
        (0.01, 6.0, 500.0),
        (0.012, 7.0, 50.0),
        (0.04, 8.4, 50.0),
    ] {
        near(c.yield_stress_pa(a).unwrap(), y, 2e-15);
        near(c.tangent_pa(a).unwrap(), slope, 1e-13);
    }
    // Trapezoids of the HARDENING force Y-Y0: triangle, plateau, ramp.
    near(
        c.hardening_energy_j_m3(0.012).unwrap(),
        0.001 + 0.004 + 0.015,
        1e-17,
    );
    near(
        c.hardening_energy_j_m3(0.04).unwrap(),
        0.02 + 0.028 * 4.7,
        3e-17,
    );
    let da = 0.009;
    let integral = c.hardening_energy_increment_j_m3(0.001, da).unwrap();
    near(integral, 0.00075 + 0.004 + 0.008, 1e-17);
    // R(end) times the width, less the integral, over the interval the
    // history moves: 0.001 + 0.009 rounds to 0.009999999999999998, not 0.01.
    let end = 0.001 + da;
    near(
        c.numerical_loss_j_m3(0.001, da).unwrap(),
        (c.yield_stress_pa(end).unwrap() - 3.0) * (end - 0.001) - integral,
        1e-17,
    );
    assert!(c.yield_stress_pa(0.040001).is_err());
    assert!(c.tangent_pa(-1e-20).is_err());
    assert!(c.hardening_energy_j_m3(f64::NAN).is_err());
}

#[test]
fn g1_scalar_consistency_crosses_every_segment_and_retains_a_yield_plateau() {
    let c = curve();
    for old in [0.0, 0.001, 0.002, 0.004, 0.006, 0.011, 0.012] {
        for target in [0.0015, 0.002, 0.003, 0.006, 0.009, 0.012, 0.019, 0.04] {
            if target <= old {
                continue;
            }
            let trial = c.yield_stress_pa(target).unwrap() + 1660.0 * (target - old);
            let delta = c.solve_increment(trial, old, 1660.0).unwrap();
            near(delta, target - old, 2e-17);
            near(
                trial - 1660.0 * delta,
                c.yield_stress_pa(old + delta).unwrap(),
                2e-14,
            );
        }
    }
    assert_eq!(c.solve_increment(8.4, 0.04, 1660.0).unwrap(), 0.0);
    assert!(matches!(
        c.solve_increment(8.5, 0.04, 1660.0),
        Err(TabulatedError::Domain { .. })
    ));
}

#[test]
fn g1_last_measured_knot_survives_scalar_and_native_return_rounding() {
    for kinematic in [0.0, 160.0] {
        for measured in [
            curve(),
            TabulatedHardening::new(vec![
                HardeningPoint {
                    plastic_strain: 0.0,
                    yield_stress_pa: 3.0,
                },
                HardeningPoint {
                    plastic_strain: 0.04,
                    yield_stress_pa: 3.0,
                },
            ])
            .unwrap(),
        ] {
            let m = TabulatedJ2::new(law(kinematic).backbone().clone(), measured).unwrap();
            let old = J2State {
                plastic_strain: [0.0; 6],
                alpha: 0.012,
            };
            m.validate_committed_history(&[0.0; 6], &old, TOL.yield_pa)
                .unwrap();
            let maximum = m.hardening().maximum_plastic_strain();
            let delta = maximum - old.alpha;
            let radius = m.hardening().yield_stress_pa(maximum).unwrap();
            // A value assembled from the physical equation can round just
            // above the separately assembled last-knot consistency value.
            let total = ((radius + (1500.0 + kinematic) * delta) / 1500.0).next_up();
            let strain = dev(total);
            let response = m.response(&strain, &old).unwrap();
            let step = m
                .advance_with_energy(&[0.0; 6], &old, &strain, TOL)
                .unwrap();
            assert!(response.state.alpha <= maximum && step.state.alpha <= maximum);
            near(response.state.alpha, maximum, 3e-16);
            near(step.state.alpha, response.state.alpha, 3e-16);
            near(step.energy.residual_j_m3, 0.0, 1e-13);
            assert_eq!(m.update_state(&strain, &old), response.state);
            assert!(m.response(&dev(total + 1e-6), &old).is_err());
            assert_eq!(old.alpha, 0.012);
        }
    }
}

#[test]
fn g1_native_tensor_return_follows_tabulated_flow_stress_without_a_fit() {
    let m = law(160.0);
    let virgin = m.initial_state();
    // mu=500 Pa. Prescribe alpha and independently solve the consistency
    // equation for a deviatoric total strain. eps_p_xx = alpha in this loading.
    for alpha in [1e-8, 0.001, 0.003, 0.005, 0.008, 0.014, 0.022] {
        let y = m.hardening().yield_stress_pa(alpha).unwrap();
        let strain = dev((y + 1660.0 * alpha) / 1500.0);
        let step = m
            .advance_with_energy(&[0.0; 6], &virgin, &strain, TOL)
            .unwrap();
        near(step.state.alpha, alpha, 2e-17);
        near(step.state.plastic_strain[0], alpha, 2e-17);
        near(step.state.plastic_strain[1], -0.5 * alpha, 2e-17);
        near(step.stress_pa[0], (2.0 / 3.0) * (y + 160.0 * alpha), 2e-13);
        near(
            m.yield_function_with_state(&step.stress_pa, &step.state),
            0.0,
            2e-13,
        );
        near(
            step.energy.final_state.hardening_j_m3,
            m.hardening().hardening_energy_j_m3(alpha).unwrap(),
            2e-16,
        );
        near(step.energy.physical_dissipation_j_m3, 3.0 * alpha, 2e-16);
        near(step.energy.residual_j_m3, 0.0, 1e-13);
    }
}

fn check_tangent(m: &TabulatedJ2, strain: &Voigt, state: &J2State) {
    let response = m.response(strain, state).unwrap();
    let h = 1e-7;
    for j in 0..6 {
        let mut plus = *strain;
        let mut minus = *strain;
        plus[j] += h;
        minus[j] -= h;
        let sp = m.stress(&plus, state);
        let sm = m.stress(&minus, state);
        for i in 0..6 {
            near(response.tangent_pa[i][j], (sp[i] - sm[i]) / (2.0 * h), 1e-5);
        }
    }
}

#[test]
fn g1_nonproportional_cyclic_tangent_and_energy_keep_complete_history() {
    let m = law(160.0);
    let mut state = m.initial_state();
    let mut old = [0.0; 6];
    for strain in [
        [0.0002, -0.0001, -0.0001, 0.0001, 0.0, 0.0],
        [0.009, -0.003, -0.006, 0.002, 0.001, 0.0],
        [0.0085, -0.0028, -0.0057, 0.0019, 0.0009, 0.0],
        [-0.004, 0.001, 0.003, -0.002, 0.002, 0.001],
        [0.010, -0.006, -0.004, 0.005, -0.002, 0.001],
    ] {
        check_tangent(&m, &strain, &state);
        let accepted = state.clone();
        let step = m.advance_with_energy(&old, &state, &strain, TOL).unwrap();
        assert_eq!(state, accepted);
        assert!(step.state.alpha >= state.alpha);
        near(trace(&step.state.plastic_strain), 0.0, 3e-17);
        let da = step.state.alpha - state.alpha;
        let dp: Voigt =
            core::array::from_fn(|i| step.state.plastic_strain[i] - state.plastic_strain[i]);
        let de: Voigt = core::array::from_fn(|i| strain[i] - old[i] - dp[i]);
        let de_stress = m.backbone().elastic.stress(&de, &());
        let expected_numerical = 0.5 * contract(&de, &de_stress)
            + 160.0 / 3.0 * contract(&dp, &dp)
            + m.hardening().numerical_loss_j_m3(state.alpha, da).unwrap();
        near(
            step.energy.numerical_dissipation_j_m3,
            expected_numerical,
            1e-14,
        );
        near(step.energy.physical_dissipation_j_m3, 3.0 * da, 1e-15);
        assert!(step.energy.numerical_dissipation_j_m3 >= 0.0);
        near(step.energy.residual_j_m3, 0.0, 1e-13);
        near(
            m.stored_energy(&strain, &step.state).unwrap().total_j_m3(),
            step.energy.final_state.total_j_m3(),
            1e-14,
        );
        old = strain;
        state = step.state;
    }
}

#[test]
fn g3_linear_table_matches_native_linear_isotropic_and_prager_response() {
    let linear = J2Plasticity::new_mixed(
        IsotropicElastic::new(1200.0, 0.2, 0.05).unwrap(),
        3.0,
        80.0,
        160.0,
    )
    .unwrap();
    let mut base = linear.clone();
    base.hardening = 0.0;
    let table = TabulatedHardening::new(vec![
        HardeningPoint {
            plastic_strain: 0.0,
            yield_stress_pa: 3.0,
        },
        HardeningPoint {
            plastic_strain: 0.05,
            yield_stress_pa: 7.0,
        },
    ])
    .unwrap();
    let m = TabulatedJ2::new(base, table).unwrap();
    let mut state = m.initial_state();
    for strain in [
        dev(0.012),
        [-0.006, 0.002, 0.004, -0.003, 0.002, 0.001],
        dev(0.016),
    ] {
        let actual = m.response(&strain, &state).unwrap();
        let expected_stress = linear.stress(&strain, &state);
        let expected_tangent = linear.tangent(&strain, &state);
        let expected_state = linear.update_state(&strain, &state);
        for i in 0..6 {
            near(actual.stress_pa[i], expected_stress[i], 2e-13);
            near(
                actual.state.plastic_strain[i],
                expected_state.plastic_strain[i],
                2e-17,
            );
            for j in 0..6 {
                near(actual.tangent_pa[i][j], expected_tangent[i][j], 2e-12);
            }
        }
        near(actual.state.alpha, expected_state.alpha, 3e-17);
        state = actual.state;
    }
}

#[test]
fn g0_plateau_unloading_and_out_of_domain_trials_do_not_erase_plastic_state() {
    let m = law(0.0);
    let first = dev((4.0 + 1500.0 * 0.004) / 1500.0);
    let initial = m.initial_state();
    let loaded = m
        .advance_with_energy(&[0.0; 6], &initial, &first, TOL)
        .unwrap();
    let zero_stress = loaded.state.plastic_strain;
    let unloaded = m
        .advance_with_energy(&first, &loaded.state, &zero_stress, TOL)
        .unwrap();
    assert_eq!(unloaded.state, loaded.state);
    assert!(unloaded.stress_pa.iter().all(|s| s.abs() < 1e-12));
    assert_eq!(unloaded.energy.physical_dissipation_j_m3, 0.0);

    let mut shorter = points();
    shorter.truncate(3);
    let short = TabulatedJ2::new(
        m.backbone().clone(),
        TabulatedHardening::new(shorter).unwrap(),
    )
    .unwrap();
    let before = unloaded.state.clone();
    assert!(matches!(
        short.response(&dev(0.02), &before),
        Err(TabulatedError::Domain { .. })
    ));
    assert!(
        short
            .advance_with_energy(&zero_stress, &before, &dev(0.02), TOL)
            .is_err()
    );
    assert_eq!(before, unloaded.state);
    assert!(short.stress(&dev(0.02), &before).iter().all(|x| x.is_nan()));
}

#[test]
fn g0_curve_admission_refuses_softening_bad_knots_and_ambiguous_backbones() {
    for bad in [
        vec![],
        vec![points()[0]],
        vec![points()[1], points()[2]],
        vec![points()[0], points()[0]],
        vec![
            points()[0],
            HardeningPoint {
                plastic_strain: 0.1,
                yield_stress_pa: 2.0,
            },
        ],
        vec![
            points()[0],
            HardeningPoint {
                plastic_strain: f64::NAN,
                yield_stress_pa: 4.0,
            },
        ],
        vec![
            points()[0],
            HardeningPoint {
                plastic_strain: 0.1,
                yield_stress_pa: f64::INFINITY,
            },
        ],
    ] {
        assert!(TabulatedHardening::new(bad).is_err());
    }
    let mut base = law(0.0).backbone().clone();
    base.hardening = 1.0;
    assert!(TabulatedJ2::new(base.clone(), curve()).is_err());
    base.hardening = 0.0;
    base.yield_stress = 2.0;
    assert!(TabulatedJ2::new(base, curve()).is_err());
    let m = law(0.0);
    let invalid = J2State {
        plastic_strain: [0.001; 6],
        alpha: 0.002,
    };
    assert!(m.response(&[0.0; 6], &invalid).is_err());
}

#[test]
fn g3_scalar_stress_units_do_not_change_the_plastic_trajectory() {
    let original = curve();
    for scale in [1e-12, 1e-3, 1e6, 1e100] {
        let scaled = TabulatedHardening::new(
            points()
                .into_iter()
                .map(|mut p| {
                    p.yield_stress_pa *= scale;
                    p
                })
                .collect(),
        )
        .unwrap();
        let delta = original.solve_increment(31.0, 0.0015, 1660.0).unwrap();
        let actual = scaled
            .solve_increment(31.0 * scale, 0.0015, 1660.0 * scale)
            .unwrap();
        near(actual, delta, 1e-17);
        near(
            scaled
                .hardening_energy_increment_j_m3(0.0015, actual)
                .unwrap()
                / scale,
            original
                .hardening_energy_increment_j_m3(0.0015, delta)
                .unwrap(),
            1e-16,
        );
    }
}

#[test]
fn g0_small_intervals_preserve_hardening_storage_and_quadratic_numerical_loss() {
    let c = TabulatedHardening::new(vec![
        HardeningPoint {
            plastic_strain: 0.0,
            yield_stress_pa: 100.0,
        },
        HardeningPoint {
            plastic_strain: 1.0,
            yield_stress_pa: 300.0,
        },
    ])
    .unwrap();
    let delta = 1e-12;
    near(
        c.strength_increment_pa(0.5, delta).unwrap(),
        200.0 * delta,
        1e-25,
    );
    near(
        c.hardening_energy_increment_j_m3(0.5, delta).unwrap(),
        (100.0 + 100.0 * delta) * delta,
        1e-25,
    );
    near(
        c.numerical_loss_j_m3(0.5, delta).unwrap(),
        100.0 * delta * delta,
        1e-37,
    );
}
