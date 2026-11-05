use super::*;
use crate::SpeciesId;

fn model(chi: f64) -> BinarySquareGradientV1 {
    BinarySquareGradientV1::new(
        BinaryRegularSolutionV1::new(
            [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
            chi,
            0.0,
            Temperature::new(250.0),
            Temperature::new(800.0),
        )
        .unwrap(),
        Temperature::new(300.0),
        AmountConcentration::new(1000.0),
        1e-6,
    )
    .unwrap()
}
fn pair(x: f64) -> [f64; 2] {
    [x, 1.0 - x]
}
fn close(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_exchange_potential_is_the_original_homogeneous_free_energy_derivative() {
    for chi in [-3.0, 0.0, 2.0, 3.0, 10.0] {
        let owner = model(chi);
        for x in [0.01, 0.2, 0.49, 0.5, 0.7, 0.99] {
            let state = owner.evaluate_composition(pair(x)).unwrap();
            assert_eq!(
                state,
                owner
                    .thermodynamics()
                    .evaluate(owner.temperature(), x)
                    .unwrap()
            );
            let mu = owner.exchange_potential_j_per_mol(pair(x)).unwrap();
            let [a, b] = state.reduced_chemical_potentials().map(Option::unwrap);
            close(mu, owner.rt_j_per_mol() * (a - b), 2e-11, 2e-14);
            let step = 1e-6;
            let plus = owner
                .evaluate_composition(pair(x + step))
                .unwrap()
                .mixing_gibbs_j_per_mol();
            let minus = owner
                .evaluate_composition(pair(x - step))
                .unwrap()
                .mixing_gibbs_j_per_mol();
            close(mu, (plus - minus) / (2.0 * step), 1e-5, 1e-8);
        }
    }
    let owner = model(3.0);
    assert!(
        owner
            .evaluate_composition(pair(0.5))
            .unwrap()
            .thermodynamic_factor()
            < 0.0
    );
    assert_eq!(owner.exchange_potential_j_per_mol(pair(0.5)).unwrap(), 0.0);
}

#[test]
fn g0_direct_free_energy_increments_reverse_and_retain_trace_changes() {
    for chi in [-3.0, 0.0, 2.0, 3.0] {
        let owner = model(chi);
        for (a, b) in [
            (0.0, 0.2),
            (0.1, 0.9),
            (0.499, 0.501),
            (0.8, 0.80001),
            (0.9, 1.0),
        ] {
            let change = owner
                .mixing_free_energy_change_j_per_mol(pair(a), pair(b))
                .unwrap();
            let initial = owner
                .evaluate_composition(pair(a))
                .unwrap()
                .mixing_gibbs_j_per_mol();
            let final_energy = owner
                .evaluate_composition(pair(b))
                .unwrap()
                .mixing_gibbs_j_per_mol();
            close(change, final_energy - initial, 2e-12, 2e-12);
            assert_eq!(
                owner
                    .mixing_free_energy_change_j_per_mol(pair(b), pair(a))
                    .unwrap(),
                -change
            );
            close(
                change,
                owner
                    .mixing_free_energy_change_j_per_mol(pair(1.0 - a), pair(1.0 - b))
                    .unwrap(),
                3e-12,
                2e-11,
            );
        }
        let a = [1e-20_f64, 1.0];
        let b = [a[0].next_up(), 1.0];
        let change = owner.mixing_free_energy_change_j_per_mol(a, b).unwrap();
        assert!(change < 0.0);
        let expected = owner.exchange_potential_j_per_mol(a).unwrap() * (b[0] - a[0]);
        close(change, expected, 0.0, 1e-14);
        assert_eq!(
            owner.mixing_free_energy_change_j_per_mol(b, a).unwrap(),
            -change
        );
    }
}

#[test]
fn g1_critical_energy_and_potential_retain_their_higher_order_terms() {
    let owner = model(2.0);
    let x = 0.50001_f64;
    let d = x - 0.5;
    // Independent Taylor limit: g(1/2+d)-g(1/2)=RT*(4d^4/3+32d^6/15+...).
    let expected = owner.rt_j_per_mol() * (4.0 / 3.0 * d.powi(4) + 32.0 / 15.0 * d.powi(6));
    let increment = owner
        .mixing_free_energy_change_j_per_mol(pair(0.5), pair(x))
        .unwrap();
    close(increment, expected, 0.0, 2e-14);
    assert!(increment > 0.0);
    let x = 0.5000001_f64;
    let z = 2.0 * (x - 0.5);
    let expected = owner.rt_j_per_mol() * (2.0 / 3.0 * z.powi(3) + 2.0 / 5.0 * z.powi(5));
    close(
        owner.exchange_potential_j_per_mol(pair(x)).unwrap(),
        expected,
        0.0,
        2e-14,
    );
}

#[test]
fn g1_convex_split_bounds_the_original_free_energy_for_both_chi_signs() {
    for chi in [-5.0, 0.0, 3.0, 20.0] {
        let owner = model(chi);
        for old in [0.0, 0.01, 0.2, 0.5, 0.9, 1.0] {
            for new in [0.001, 0.1, 0.4, 0.7, 0.999] {
                let mu = owner
                    .convex_split_exchange_potential_j_per_mol(pair(new), pair(old))
                    .unwrap();
                let energy = owner
                    .mixing_free_energy_change_j_per_mol(pair(old), pair(new))
                    .unwrap();
                assert!(energy <= mu * (new - old) + 2e-10);
                let curvature = owner.convex_curvature_j_per_mol(pair(new)).unwrap();
                assert!(curvature > 0.0);
                let h = 1e-7;
                let plus = owner
                    .convex_split_exchange_potential_j_per_mol(pair(new + h), pair(old))
                    .unwrap();
                let minus = owner
                    .convex_split_exchange_potential_j_per_mol(pair(new - h), pair(old))
                    .unwrap();
                close(curvature, (plus - minus) / (2.0 * h), 1e-3, 1e-7);
            }
        }
    }
}

#[test]
fn g0_source_domain_and_continuum_regularization_are_explicit() {
    let owner = model(3.0);
    for (c, kappa, t) in [
        (0.0, 1e-6, 300.0),
        (1000.0, -1.0, 300.0),
        (1000.0, 0.0, 300.0),
        (1000.0, 1e-6, 900.0),
    ] {
        assert!(
            BinarySquareGradientV1::new(
                owner.thermodynamics().clone(),
                Temperature::new(t),
                AmountConcentration::new(c),
                kappa
            )
            .is_err()
        );
    }
    let ideal = model(0.0);
    assert!(
        BinarySquareGradientV1::new(
            ideal.thermodynamics().clone(),
            ideal.temperature(),
            AmountConcentration::new(1000.0),
            0.0
        )
        .is_ok()
    );
    for invalid in [[0.0, 1.0], [1.0, 0.0], [0.4, 0.4], [f64::NAN, 0.5]] {
        assert!(owner.exchange_potential_j_per_mol(invalid).is_err());
    }
    assert!(
        owner
            .exchange_potential_j_per_mol([1e-300, 1.0])
            .unwrap()
            .is_finite()
    );
    assert!(owner.evaluate_composition([0.0, 1.0]).is_ok());
}

fn caloric_model(a: f64, b: f64, temperature: f64) -> BinarySquareGradientV1 {
    BinarySquareGradientV1::new(
        BinaryRegularSolutionV1::new(
            [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
            a,
            b,
            Temperature::new(250.0),
            Temperature::new(800.0),
        )
        .unwrap(),
        Temperature::new(temperature),
        AmountConcentration::new(1000.0),
        1e-6,
    )
    .unwrap()
}

#[test]
fn g0_direct_caloric_increments_close_h_minus_t_s_and_reverse_exactly() {
    for (a, b, temperature) in [(0.3, 900.0, 300.0), (-0.2, -300.0, 450.0)] {
        let owner = caloric_model(a, b, temperature);
        for (old, new) in [
            (0.0, 0.2),
            (0.1, 0.9),
            (0.3, 0.30001),
            (0.8, 0.8001),
            (0.9, 1.0),
        ] {
            let old = pair(old);
            let new = pair(new);
            let h = owner.mixing_enthalpy_change_j_per_mol(old, new).unwrap();
            let s = owner.mixing_entropy_change_j_per_mol_k(old, new).unwrap();
            let g = owner.mixing_free_energy_change_j_per_mol(old, new).unwrap();
            // The identity may subtract two larger increments. Its rounding
            // allowance follows those increments, not only their difference.
            let rounding = 32.0 * f64::EPSILON * (h.abs() + (temperature * s).abs());
            close(g, h - temperature * s, rounding, 2e-12);
            let before = owner.evaluate_composition(old).unwrap();
            let after = owner.evaluate_composition(new).unwrap();
            close(
                h,
                after.mixing_enthalpy_j_per_mol() - before.mixing_enthalpy_j_per_mol(),
                2e-12,
                2e-11,
            );
            close(
                s,
                after.mixing_entropy_j_per_mol_k() - before.mixing_entropy_j_per_mol_k(),
                2e-14,
                2e-11,
            );
            assert_eq!(
                owner.mixing_enthalpy_change_j_per_mol(new, old).unwrap(),
                -h
            );
            assert_eq!(
                owner.mixing_entropy_change_j_per_mol_k(new, old).unwrap(),
                -s
            );
        }
    }
}

#[test]
fn g0_caloric_increments_admit_pure_limits_and_preserve_each_trace_coordinate() {
    let owner = caloric_model(0.7, 1200.0, 300.0);
    assert_eq!(
        owner
            .mixing_enthalpy_change_j_per_mol(pair(0.0), pair(1.0))
            .unwrap(),
        0.0
    );
    assert_eq!(
        owner
            .mixing_entropy_change_j_per_mol_k(pair(0.0), pair(1.0))
            .unwrap(),
        0.0
    );
    let trace = 1e-30_f64;
    for (old, new) in [([0.0, 1.0], [trace, 1.0]), ([1.0, 0.0], [1.0, trace])] {
        let h = owner.mixing_enthalpy_change_j_per_mol(old, new).unwrap();
        let s = owner.mixing_entropy_change_j_per_mol_k(old, new).unwrap();
        close(h, (1200.0 * trace) * R, 0.0, 2e-14);
        close(
            s,
            -R * ((-trace).ln_1p() + trace * trace.ln() + 0.7 * trace),
            0.0,
            2e-14,
        );
        assert_eq!(
            owner.mixing_enthalpy_change_j_per_mol(new, old).unwrap(),
            -h
        );
        assert_eq!(
            owner.mixing_entropy_change_j_per_mol_k(new, old).unwrap(),
            -s
        );
    }
    for (old, new) in [
        ([trace, 1.0], [trace.next_up(), 1.0]),
        ([1.0, trace], [1.0, trace.next_up()]),
    ] {
        let d = trace.next_up() - trace;
        close(
            owner.mixing_enthalpy_change_j_per_mol(old, new).unwrap(),
            (1200.0 * d) * R,
            0.0,
            2e-14,
        );
        close(
            owner.mixing_entropy_change_j_per_mol_k(old, new).unwrap(),
            -R * (trace.ln() + 0.7) * d,
            0.0,
            2e-14,
        );
    }
    for invalid in [[0.4, 0.4], [-0.1, 1.1], [f64::NAN, 0.5]] {
        assert!(
            owner
                .mixing_enthalpy_change_j_per_mol(invalid, pair(0.5))
                .is_err()
        );
        assert!(
            owner
                .mixing_entropy_change_j_per_mol_k(pair(0.5), invalid)
                .is_err()
        );
    }
}

#[test]
fn g1_caloric_increments_retain_quadratic_enthalpy_and_critical_quartic_entropy() {
    let entropy = caloric_model(2.0, 0.0, 300.0);
    let shifted = 0.50000001_f64;
    let d = shifted - 0.5;
    let expected = -R * (4.0 / 3.0 * d.powi(4) + 32.0 / 15.0 * d.powi(6));
    let change = entropy
        .mixing_entropy_change_j_per_mol_k(pair(0.5), pair(shifted))
        .unwrap();
    assert!(change < 0.0);
    close(change, expected, 0.0, 2e-14);
    assert!(
        change.abs()
            < f64::EPSILON
                * entropy
                    .evaluate_composition(pair(0.5))
                    .unwrap()
                    .mixing_entropy_j_per_mol_k()
                    .abs()
    );
    assert_eq!(
        entropy
            .mixing_entropy_change_j_per_mol_k(pair(shifted), pair(0.5))
            .unwrap(),
        -change
    );

    let enthalpy = caloric_model(0.0, 600.0, 300.0);
    let shifted = 0.5_f64.next_up();
    let d = shifted - 0.5;
    let change = enthalpy
        .mixing_enthalpy_change_j_per_mol(pair(0.5), pair(shifted))
        .unwrap();
    close(change, -(600.0 * d * d) * R, 0.0, 2e-15);
    assert!(change < 0.0);
    assert_eq!(
        enthalpy
            .evaluate_composition(pair(0.5))
            .unwrap()
            .mixing_enthalpy_j_per_mol(),
        enthalpy
            .evaluate_composition(pair(shifted))
            .unwrap()
            .mixing_enthalpy_j_per_mol()
    );
    assert_eq!(
        enthalpy
            .mixing_enthalpy_change_j_per_mol(pair(shifted), pair(0.5))
            .unwrap(),
        -change
    );
}

#[test]
fn g0_direct_enthalpy_need_not_form_an_overflowing_reference_energy() {
    let owner = caloric_model(0.0, 1e308, 300.0);
    let old = pair(0.5);
    let new = pair(0.50000001);
    let d = new[0] - old[0];
    assert!(!(R * 1e308).is_finite());
    assert!(owner.evaluate_composition(old).is_err());
    let change = owner.mixing_enthalpy_change_j_per_mol(old, new).unwrap();
    close(change, -((1e308 * d) * d) * R, 0.0, 3e-15);
    assert_eq!(
        owner.mixing_enthalpy_change_j_per_mol(new, old).unwrap(),
        -change
    );
}
