use super::*;
use fs_qty::AmountConcentration;
use gradient::BinarySquareGradientV1;

fn pair_law(a: f64, b: f64) -> BinaryRegularSolutionV1 {
    BinaryRegularSolutionV1::new(
        [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
        a,
        b,
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g0_zero_excess_cp_retains_original_pair_owner_and_arithmetic() {
    for (a, b) in [(0.4, 600.0), (3.2, -450.0), (0.0, 0.0)] {
        let original = pair_law(a, b);
        for coefficient in [0.0, -0.0] {
            let extended = original
                .clone()
                .with_excess_heat_capacity(coefficient, Temperature::new(450.0))
                .unwrap();
            assert_eq!(extended, original);
            assert_eq!(extended.excess_heat_capacity_coefficient_j_per_mol_k(), 0.0);
            for t in [100.0, 300.0, 450.0, 1000.0] {
                for x in [0.0, 1e-20, 0.25, 0.5, 1.0] {
                    assert_eq!(
                        extended.evaluate(Temperature::new(t), x).unwrap(),
                        original.evaluate(Temperature::new(t), x).unwrap()
                    );
                }
            }
        }
    }
}

#[test]
fn g1_excess_cp_state_matches_independent_high_precision_thermodynamic_anchor() {
    let law = pair_law(0.4, 600.0)
        .with_excess_heat_capacity(40.0, Temperature::new(300.0))
        .unwrap();
    let state = law.evaluate(Temperature::new(600.0), 0.25).unwrap();
    // Independent 70-digit decimal evaluation of H0+C*(T-T0),
    // S0+C*ln(T/T0), and G=H-T*S, including ideal mixing entropy.
    close(state.chi(), 0.470_789_348_967_709_8, 0.0, 2e-14);
    close(
        state.mixing_enthalpy_j_per_mol(),
        3185.377_044_542_239_5,
        0.0,
        2e-15,
    );
    close(
        state.mixing_entropy_j_per_mol_k(),
        9.250_533_696_644_975,
        0.0,
        2e-14,
    );
    close(
        state.mixing_gibbs_j_per_mol(),
        -2364.943_173_444_745_5,
        0.0,
        2e-14,
    );
}

#[test]
fn g1_signed_excess_cp_is_both_enthalpy_and_entropy_temperature_response() {
    for coefficient in [-80.0, 40.0] {
        let original = pair_law(0.4, 600.0);
        let law = original
            .clone()
            .with_excess_heat_capacity(coefficient, Temperature::new(300.0))
            .unwrap();
        assert_eq!(
            law.excess_heat_capacity_coefficient_j_per_mol_k(),
            coefficient
        );
        for t in [200.0, 300.0, 700.0] {
            for x in [0.07, 0.25, 0.5, 0.8] {
                let state = law.evaluate(Temperature::new(t), x).unwrap();
                let step = 0.01;
                let lower = law.evaluate(Temperature::new(t - step), x).unwrap();
                let upper = law.evaluate(Temperature::new(t + step), x).unwrap();
                let expected_cp = coefficient * x * (1.0 - x);
                close(
                    (upper.mixing_enthalpy_j_per_mol() - lower.mixing_enthalpy_j_per_mol())
                        / (2.0 * step),
                    expected_cp,
                    2e-9,
                    2e-10,
                );
                close(
                    t * (upper.mixing_entropy_j_per_mol_k() - lower.mixing_entropy_j_per_mol_k())
                        / (2.0 * step),
                    expected_cp,
                    2e-9,
                    2e-9,
                );
                close(
                    (upper.mixing_gibbs_j_per_mol() - lower.mixing_gibbs_j_per_mol())
                        / (2.0 * step),
                    -state.mixing_entropy_j_per_mol_k(),
                    2e-8,
                    2e-9,
                );
                close(
                    state.mixing_gibbs_j_per_mol(),
                    state.mixing_enthalpy_j_per_mol() - t * state.mixing_entropy_j_per_mol_k(),
                    1e-10,
                    3e-14,
                );
                if t == 300.0 {
                    let reference = original.evaluate(Temperature::new(t), x).unwrap();
                    assert_eq!(state.chi(), reference.chi());
                    close(
                        state.mixing_enthalpy_j_per_mol(),
                        reference.mixing_enthalpy_j_per_mol(),
                        1e-12,
                        2e-15,
                    );
                    close(
                        state.mixing_entropy_j_per_mol_k(),
                        reference.mixing_entropy_j_per_mol_k(),
                        2e-15,
                        2e-15,
                    );
                }
            }
        }
    }
}

#[test]
fn g1_excess_cp_supports_closed_and_reentrant_miscibility_temperature_windows() {
    for (a, coefficient, phase_counts) in [(2.2, 100.0, [1, 2, 1]), (1.8, -100.0, [2, 1, 2])] {
        let law = pair_law(a, 0.0)
            .with_excess_heat_capacity(coefficient, Temperature::new(400.0))
            .unwrap();
        for (t, expected_phases) in [200.0, 400.0, 800.0].into_iter().zip(phase_counts) {
            let state = law
                .equilibrate(
                    Temperature::new(t),
                    0.5,
                    BinaryCoexistenceOptionsV1::default(),
                    || false,
                )
                .unwrap();
            assert_eq!(state.phases().len(), expected_phases);
            close(
                state.mixing_gibbs_j_per_mol(),
                state.mixing_enthalpy_j_per_mol() - t * state.mixing_entropy_j_per_mol_k(),
                2e-10,
                4e-14,
            );
        }
    }
}

#[test]
fn g0_near_datum_chi_retains_quadratic_cp_response_and_entropy_increment() {
    let reference = 300.0_f64;
    let coefficient = 1e30;
    let law = pair_law(0.0, 0.0)
        .with_excess_heat_capacity(coefficient, Temperature::new(reference))
        .unwrap();
    assert_eq!(law.chi(Temperature::new(reference)).unwrap(), 0.0);
    for t in [reference.next_down(), reference.next_up()] {
        let relative = (t - reference) / reference;
        // The independent leading Taylor limit is accurate far beyond double
        // precision here. Subtracting H/(RT)-S/R loses this nonzero response.
        let expected = -(coefficient / (2.0 * R)) * relative * relative;
        let actual = law.chi(Temperature::new(t)).unwrap();
        assert!(actual < 0.0);
        close(actual, expected, 0.0, 2e-15);
    }
    let law = pair_law(1e20, 600.0)
        .with_excess_heat_capacity(37.0, Temperature::new(300.0))
        .unwrap();
    for (from, to) in [(600.0_f64, 600.0_f64.next_up()), (100.0, 1000.0)] {
        let forward = law
            .interaction_entropy_change_j_per_mol_k(Temperature::new(from), Temperature::new(to))
            .unwrap();
        let reverse = law
            .interaction_entropy_change_j_per_mol_k(Temperature::new(to), Temperature::new(from))
            .unwrap();
        assert!(forward > 0.0);
        assert_eq!(reverse, -forward);
        let expected = if to == from.next_up() {
            37.0 * (to - from) / from
        } else {
            37.0 * (to / from).ln()
        };
        close(forward, expected, 0.0, 2e-14);
    }
}

#[test]
fn g1_isothermal_diffusion_heat_and_entropy_use_the_same_extended_interaction() {
    for coefficient in [-40.0, 80.0] {
        for t in [200.0, 600.0] {
            let thermodynamics = pair_law(0.4, 600.0)
                .with_excess_heat_capacity(coefficient, Temperature::new(300.0))
                .unwrap();
            let law = BinarySquareGradientV1::new(
                thermodynamics,
                Temperature::new(t),
                AmountConcentration::new(1000.0),
                1e-6,
            )
            .unwrap();
            for (old, new) in [
                ([0.1, 0.9], [0.4, 0.6]),
                ([0.5, 0.5], [0.50001, 1.0 - 0.50001]),
                ([0.0, 1.0], [0.2, 0.8]),
                ([1e-20_f64, 1.0], [1e-20_f64.next_up(), 1.0]),
            ] {
                let h = law.mixing_enthalpy_change_j_per_mol(old, new).unwrap();
                let s = law.mixing_entropy_change_j_per_mol_k(old, new).unwrap();
                let g = law.mixing_free_energy_change_j_per_mol(old, new).unwrap();
                assert_eq!(law.mixing_enthalpy_change_j_per_mol(new, old).unwrap(), -h);
                assert_eq!(law.mixing_entropy_change_j_per_mol_k(new, old).unwrap(), -s);
                let rounding = 40.0 * f64::EPSILON * (h.abs() + (t * s).abs());
                close(g, h - t * s, rounding, 3e-13);
                // An independent composition product difference preserves the
                // real dilute increment even when each majority rounds to one.
                let delta_q = (new[0] - old[0]) * (0.5 * new[1] + 0.5 * old[1])
                    + (new[1] - old[1]) * (0.5 * new[0] + 0.5 * old[0]);
                let expected = (R * 600.0 + coefficient * (t - 300.0)) * delta_q;
                close(h, expected, 0.0, 2e-10);
            }
        }
    }
}

#[test]
fn g0_excess_cp_source_domain_and_pure_limits_are_explicit() {
    let original = pair_law(0.4, 600.0);
    for coefficient in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            original
                .clone()
                .with_excess_heat_capacity(coefficient, Temperature::new(300.0))
                .is_err()
        );
    }
    for t in [f64::NAN, f64::INFINITY, 0.0, 99.0, 1001.0] {
        assert!(
            original
                .clone()
                .with_excess_heat_capacity(40.0, Temperature::new(t))
                .is_err()
        );
    }
    for reference in [100.0, 1000.0] {
        let law = original
            .clone()
            .with_excess_heat_capacity(-40.0, Temperature::new(reference))
            .unwrap();
        for t in [100.0, 300.0, 1000.0] {
            for x in [0.0, 1.0] {
                let state = law.evaluate(Temperature::new(t), x).unwrap();
                assert_eq!(state.mixing_enthalpy_j_per_mol(), 0.0);
                assert_eq!(state.mixing_entropy_j_per_mol_k(), 0.0);
                assert_eq!(state.mixing_gibbs_j_per_mol(), 0.0);
            }
        }
        for t in [99.0, 1001.0, f64::NAN] {
            assert!(law.chi(Temperature::new(t)).is_err());
            assert!(
                law.interaction_enthalpy_j_per_mol(Temperature::new(t))
                    .is_err()
            );
            assert!(
                law.interaction_entropy_j_per_mol_k(Temperature::new(t))
                    .is_err()
            );
            assert!(
                law.interaction_entropy_change_j_per_mol_k(
                    Temperature::new(t),
                    Temperature::new(300.0)
                )
                .is_err()
            );
        }
    }
}
