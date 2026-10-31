use super::*;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn model(masses: &[f64]) -> MaxwellStefanV1 {
    MaxwellStefanV1::new(
        masses
            .iter()
            .enumerate()
            .map(|(i, &mass)| (id(&format!("S{i:03}")), MolarMass::new(mass)))
            .collect(),
    )
    .unwrap()
}
fn matrix(n: usize, d: f64) -> Vec<f64> {
    let mut values = vec![d; n * n];
    for i in 0..n {
        values[i * n + i] = 0.0;
    }
    values
}
fn run(
    model: &MaxwellStefanV1,
    x: &[f64],
    gradient: &[f64],
    pressure_gradient: f64,
    d: &[f64],
) -> (Vec<f64>, Vec<f64>, f64, Option<f64>) {
    let mut work = model.workspace().unwrap();
    let value = model
        .evaluate(
            x,
            gradient,
            Temperature::new(500.0),
            Pressure::new(R * 500.0),
            pressure_gradient,
            d,
            MaxwellStefanSolveV1::default(),
            &mut work,
            || false,
        )
        .unwrap();
    assert!(value.maximum_gradient_residual_per_m() < 2e-12);
    assert!(value.mass_frame_residual_kg_m2_s().abs() < 2e-13);
    (
        value.mass_fluxes().to_vec(),
        value.molar_fluxes().to_vec(),
        value.maximum_gradient_residual_per_m(),
        value.drag_entropy_production_w_m3_k(),
    )
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-20),
        "actual={actual:e},expected={expected:e}"
    );
}

#[test]
fn g1_binary_unequal_masses_match_analytic_barycentric_flux() {
    let masses = [0.002, 0.028];
    let owner = model(&masses);
    let x = [0.25, 0.75];
    let d = matrix(2, 0.5);
    let (j, z, _, entropy) = run(&owner, &x, &[0.2, -0.2], 0.0, &d);
    let mean = x[0] * masses[0] + x[1] * masses[1];
    let expected = -0.5 * 0.2 * masses[0] * masses[1] / mean;
    near(j[0], expected, 2e-14);
    near(j[1], -expected, 2e-14);
    // Equal and opposite molar fluxes would violate the actual gas mass frame.
    assert!((z[0] + z[1]).abs() > 0.1);
    let expected_entropy = -R * (z[0] * 0.2 / x[0] + z[1] * (-0.2) / x[1]);
    near(entropy.unwrap(), expected_entropy, 2e-14);
}

#[test]
fn g1_molar_binary_matches_fick_without_imposing_zero_mass_flux() {
    let owner = model(&[0.002, 0.028]);
    let mut workspace = owner.workspace().unwrap();
    let evaluated = owner
        .evaluate_molar(
            &[0.25, 0.75],
            &[0.2, -0.2],
            Temperature::new(500.0),
            Pressure::new(R * 500.0),
            0.0,
            &matrix(2, 0.5),
            MaxwellStefanMolarSolveV1::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    near(evaluated.molar_fluxes()[0], -0.1, 2e-14);
    near(evaluated.molar_fluxes()[1], 0.1, 2e-14);
    near(evaluated.mass_fluxes()[0], -0.0002, 2e-14);
    near(evaluated.mass_fluxes()[1], 0.0028, 2e-14);
    assert!(evaluated.molar_frame_residual_mol_m2_s().abs() < 1e-16);
    assert!(sum(evaluated.mass_fluxes().iter().copied()).abs() > 0.002);
    near(
        evaluated.drag_entropy_production_w_m3_k().unwrap(),
        -R * (-0.1 * 0.2 / 0.25 + 0.1 * -0.2 / 0.75),
        3e-14,
    );
}

#[test]
fn g3_reference_frame_change_preserves_pair_friction_and_entropy() {
    let owner = model(&[0.002, 0.028, 0.044]);
    let x = [0.2, 0.3, 0.5];
    let gradient = [0.1, -0.3, 0.2];
    let pairs = [0.0, 2e-5, 3e-5, 2e-5, 0.0, 7e-6, 3e-5, 7e-6, 0.0];
    let (_, barycentric, _, entropy) = run(&owner, &x, &gradient, 0.7, &pairs);
    let drift = sum(barycentric.iter().copied());
    let mut workspace = owner.workspace().unwrap();
    let address = workspace.molar_flux.as_ptr();
    let evaluated = owner
        .evaluate_molar(
            &x,
            &gradient,
            Temperature::new(500.0),
            Pressure::new(R * 500.0),
            0.7,
            &pairs,
            MaxwellStefanMolarSolveV1::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    for i in 0..3 {
        near(
            evaluated.molar_fluxes()[i],
            barycentric[i] - x[i] * drift,
            1e-12,
        );
    }
    near(
        evaluated.drag_entropy_production_w_m3_k().unwrap(),
        entropy.unwrap(),
        1e-12,
    );
    assert_eq!(evaluated.molar_fluxes().as_ptr(), address);
    assert!(evaluated.maximum_gradient_residual_per_m() < 1e-12);
    assert!(evaluated.molar_frame_residual_mol_m2_s().abs() < 1e-19);
}

#[test]
fn g0_molar_primary_flux_survives_mass_underflow_and_cancelled_workspace_reuse() {
    let owner = model(&[0.002, 0.028]);
    let mut workspace = owner.workspace().unwrap();
    let gradient = f64::from_bits(2);
    let pair = matrix(2, 1.0);
    let mut polls = 0;
    let evaluated = owner
        .evaluate_molar(
            &[0.5, 0.5],
            &[gradient, -gradient],
            Temperature::new(500.0),
            Pressure::new(R * 500.0),
            0.0,
            &pair,
            MaxwellStefanMolarSolveV1::default(),
            &mut workspace,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(evaluated.molar_fluxes(), &[-gradient, gradient]);
    assert!(evaluated.mass_fluxes().iter().all(|&value| value == 0.0));
    assert_eq!(evaluated.maximum_gradient_residual_per_m(), 0.0);
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert_eq!(
            owner
                .evaluate_molar(
                    &[0.5, 0.5],
                    &[gradient, -gradient],
                    Temperature::new(500.0),
                    Pressure::new(R * 500.0),
                    0.0,
                    &pair,
                    MaxwellStefanMolarSolveV1::default(),
                    &mut workspace,
                    || {
                        seen += 1;
                        seen == stop
                    }
                )
                .unwrap_err(),
            MaxwellStefanErrorV1::Cancelled
        );
    }
    let evaluated = owner
        .evaluate_molar(
            &[1.0, 0.0],
            &[0.1, -0.1],
            Temperature::new(500.0),
            Pressure::new(R * 500.0),
            0.0,
            &pair,
            MaxwellStefanMolarSolveV1::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    near(evaluated.molar_fluxes()[0], -0.1, 2e-14);
    near(evaluated.molar_fluxes()[1], 0.1, 2e-14);
    assert_eq!(evaluated.drag_entropy_production_w_m3_k(), None);
}

#[test]
fn g1_ternary_crossdiffusion_moves_species_without_its_own_gradient() {
    let owner = model(&[1.0, 1.0, 1.0]);
    let d = [0.0, 1.0, 2.0, 1.0, 0.0, 3.0, 2.0, 3.0, 0.0];
    let (j, _, _, entropy) = run(&owner, &[1.0 / 3.0; 3], &[1.0, -1.0, 0.0], 0.0, &d);
    // Independent solution after eliminating z3=-z1-z2:
    // 4z1-z2=-6, 5z2-2z1=9, at total concentration1 mol/m³.
    for (actual, expected) in j.iter().zip([-7.0 / 6.0, 4.0 / 3.0, -1.0 / 6.0]) {
        near(*actual, expected, 4e-14);
    }
    assert!(entropy.unwrap() > 0.0);
}

#[test]
fn g1_original_rows_recover_a_manufactured_unequal_mass_ternary_flux() {
    let masses = [0.002, 0.028, 0.044];
    let owner = model(&masses);
    let x = [0.2, 0.3, 0.5];
    let d = [0.0, 2e-5, 3e-5, 2e-5, 0.0, 7e-6, 3e-5, 7e-6, 0.0];
    let original = [1e-5, -4e-6, -6e-6];
    let molar = [
        original[0] / masses[0],
        original[1] / masses[1],
        original[2] / masses[2],
    ];
    let mut gradient = [0.0; 3];
    for i in 0..3 {
        for j in 0..3 {
            if i != j {
                gradient[i] -= (x[j] * molar[i] - x[i] * molar[j]) / d[i * 3 + j];
            }
        }
    }
    let (actual, _, residual, _) = run(&owner, &x, &gradient, 0.0, &d);
    for (a, b) in actual.iter().zip(original) {
        near(*a, b, 3e-13);
    }
    assert!(residual < 2e-12);
}

#[test]
fn g0_exact_zero_species_have_finite_inward_flux_and_structural_zeros_remain_zero() {
    let masses = [0.002, 0.028, 0.044];
    let owner = model(&masses);
    let d = [0.0, 2e-5, 3e-5, 2e-5, 0.0, 7e-6, 3e-5, 7e-6, 0.0];
    let (j, z, _, entropy) = run(&owner, &[0.0, 1.0, 0.0], &[0.1, -0.1, 0.0], 0.0, &d);
    near(z[0], -2e-6, 2e-14);
    near(j[0], -4e-9, 2e-14);
    near(j[1], 4e-9, 2e-14);
    assert_eq!(j[2], 0.0);
    assert_eq!(entropy, None); // No artificial finite chemical potential at x=0.
    let (j, _, _, entropy) = run(&owner, &[0.0, 1.0, 0.0], &[0.0; 3], 1000.0, &d);
    assert!(j.iter().all(|x| *x == 0.0));
    assert_eq!(entropy, Some(0.0));
}

#[test]
fn g1_barodiffusion_uses_mass_fraction_and_vanishes_for_equal_masses() {
    let masses = [0.002, 0.028];
    let owner = model(&masses);
    let x = [0.25, 0.75];
    let mean = x[0] * masses[0] + x[1] * masses[1];
    let driving = (x[0] - x[0] * masses[0] / mean) * 0.3;
    let (j, _, _, entropy) = run(&owner, &x, &[0.0, 0.0], 0.3, &matrix(2, 0.5));
    let expected = -0.5 * driving * masses[0] * masses[1] / mean;
    near(j[0], expected, 3e-14);
    near(j[1], -expected, 3e-14);
    assert!(entropy.unwrap() > 0.0);
    let (j, _, _, _) = run(
        &model(&[0.028, 0.028]),
        &x,
        &[0.0, 0.0],
        0.3,
        &matrix(2, 0.5),
    );
    assert!(j.iter().all(|x| x.abs() < 1e-18));
}

#[test]
fn g3_canonical_permutation_and_pressure_diffusivity_scaling_preserve_flux() {
    let owner = MaxwellStefanV1::new(vec![
        (id("B"), MolarMass::new(0.028)),
        (id("A"), MolarMass::new(0.002)),
    ])
    .unwrap();
    let canonical = MaxwellStefanV1::new(vec![
        (id("A"), MolarMass::new(0.002)),
        (id("B"), MolarMass::new(0.028)),
    ])
    .unwrap();
    assert_eq!(owner, canonical);
    let x = [0.25, 0.75];
    let gradient = [0.2, -0.2];
    let d = matrix(2, 0.5);
    let (expected, _, _, _) = run(&owner, &x, &gradient, 0.0, &d);
    let mut work = canonical.workspace().unwrap();
    let pointer = work.mass_flux.as_ptr();
    for (p_scale, d_scale) in [(1.0, 1.0), (2.0, 0.5), (0.125, 8.0)] {
        let pair = matrix(2, 0.5 * d_scale);
        let result = canonical
            .evaluate(
                &x,
                &gradient,
                Temperature::new(500.0),
                Pressure::new(R * 500.0 * p_scale),
                0.0,
                &pair,
                MaxwellStefanSolveV1::default(),
                &mut work,
                || false,
            )
            .unwrap();
        for (a, b) in result.mass_fluxes().iter().zip(&expected) {
            near(*a, *b, 2e-14);
        }
        assert_eq!(result.mass_fluxes().as_ptr(), pointer);
    }
}

#[test]
fn g0_bounds_original_equation_refusals_and_reuse_after_cancellation() {
    assert!(MaxwellStefanV1::new(vec![]).is_err());
    assert!(MaxwellStefanV1::new(vec![(id("A"), MolarMass::new(0.0))]).is_err());
    assert!(
        MaxwellStefanV1::new(vec![
            (id("A"), MolarMass::new(1.0)),
            (id("A"), MolarMass::new(1.0))
        ])
        .is_err()
    );
    assert!(
        MaxwellStefanV1::new(
            (0..129)
                .map(|i| (id(&format!("S{i:03}")), MolarMass::new(1.0)))
                .collect()
        )
        .is_err()
    );
    let owner = model(&[0.002, 0.028]);
    let mut work = owner.workspace().unwrap();
    let t = Temperature::new(500.0);
    let p = Pressure::new(R * 500.0);
    for (x, g, d) in [
        (vec![0.2, 0.2], vec![0.1, -0.1], matrix(2, 0.5)),
        (vec![0.5, 0.5], vec![0.1, 0.1], matrix(2, 0.5)),
        (vec![-0.1, 1.1], vec![0.1, -0.1], matrix(2, 0.5)),
        (vec![0.5, 0.5], vec![0.1, -0.1], vec![0.0, 0.5, 0.6, 0.0]),
        (vec![0.5, 0.5], vec![0.1, -0.1], vec![0.0, 0.0, 0.0, 0.0]),
    ] {
        assert!(
            owner
                .evaluate(
                    &x,
                    &g,
                    t,
                    p,
                    0.0,
                    &d,
                    MaxwellStefanSolveV1::default(),
                    &mut work,
                    || false
                )
                .is_err()
        );
    }
    let n = MAX_MAXWELL_STEFAN_SPECIES_V1;
    let owner = model(&vec![1.0; n]);
    let mut work = owner.workspace().unwrap();
    let x = vec![1.0 / n as f64; n];
    let mut g = vec![0.0; n];
    g[0] = 0.1;
    g[n - 1] = -0.1;
    let d = matrix(n, 0.5);
    let mut polls = 0;
    let error = owner
        .evaluate(
            &x,
            &g,
            t,
            p,
            0.0,
            &d,
            MaxwellStefanSolveV1::default(),
            &mut work,
            || {
                polls += 1;
                polls > 2 * n + 10
            },
        )
        .unwrap_err();
    assert_eq!(error, MaxwellStefanErrorV1::Cancelled);
    let result = owner
        .evaluate(
            &x,
            &g,
            t,
            p,
            0.0,
            &d,
            MaxwellStefanSolveV1::default(),
            &mut work,
            || false,
        )
        .unwrap();
    near(result.mass_fluxes()[0], -0.05, 1e-12);
    near(result.mass_fluxes()[n - 1], 0.05, 1e-12);
    assert!(result.mass_frame_residual_kg_m2_s().abs() < 1e-15);
}

#[test]
fn g1_small_driving_forces_retain_normal_and_subnormal_binary_fick_flux() {
    let owner = model(&[1.0, 1.0]);
    for (gradient, diffusivity) in [(1e-20, 1e5), (f64::from_bits(32), 1.0)] {
        let (mass, _, residual, _) = run(
            &owner,
            &[0.5, 0.5],
            &[gradient, -gradient],
            0.0,
            &matrix(2, diffusivity),
        );
        // Equal molecular masses and concentration 1 mol/m³ give j=-D grad x.
        // Both driving forces are far below the absolute residual tolerance;
        // that tolerance must not become a cutoff for representable transport.
        let expected = diffusivity * gradient;
        assert!(mass[0] < 0.0 && mass[1] > 0.0);
        if expected.is_subnormal() {
            assert!(mass[0].abs().to_bits().abs_diff(expected.to_bits()) <= 1);
            assert!(mass[1].to_bits().abs_diff(expected.to_bits()) <= 1);
        } else {
            assert!((mass[0] / -expected - 1.0).abs() < 2e-14);
            assert!((mass[1] / expected - 1.0).abs() < 2e-14);
        }
        assert!(residual <= gradient * 1e-13 + f64::from_bits(1));
    }
}

#[test]
fn g0_rounded_zero_flux_is_admitted_only_by_original_equation_budgets() {
    let owner = model(&[0.028, 0.028]);
    let mut work = owner.workspace().unwrap();
    let gradient = f64::from_bits(4);
    let pair = matrix(2, 1e-5);
    let controls = MaxwellStefanSolveV1::default();
    let result = owner
        .evaluate(
            &[0.5, 0.5],
            &[gradient, -gradient],
            Temperature::new(500.0),
            Pressure::new(R * 500.0),
            0.0,
            &pair,
            controls,
            &mut work,
            || false,
        )
        .unwrap();
    // The exact nonzero mass flux lies below half the smallest binary64 value.
    // Its rounded zero leaves the actual original-row residual equal to grad x.
    assert!(result.mass_fluxes().iter().all(|&j| j == 0.0));
    assert_eq!(result.maximum_gradient_residual_per_m(), gradient);
    let strict = MaxwellStefanSolveV1 {
        absolute_gradient_tolerance_per_m: f64::from_bits(1),
        ..controls
    };
    let error = owner
        .evaluate(
            &[0.5, 0.5],
            &[gradient, -gradient],
            Temperature::new(500.0),
            Pressure::new(R * 500.0),
            0.0,
            &pair,
            strict,
            &mut work,
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, MaxwellStefanErrorV1::Residual {
        maximum_gradient_residual_per_m, mass_frame_residual_kg_m2_s: 0.0,
    } if maximum_gradient_residual_per_m == gradient));
}

#[test]
fn g1_inactive_tiny_diffusivity_does_not_erase_a_representable_response() {
    let owner = model(&[1.0, 1.0, 1.0]);
    let gradient = 1e-30;
    let pair = [0.0, 1.0, 1.0, 1.0, 0.0, 1e-300, 1.0, 1e-300, 0.0];
    let (mass, _, residual, _) = run(
        &owner,
        &[1.0, 0.0, 0.0],
        &[gradient, -gradient, 0.0],
        0.0,
        &pair,
    );
    // Species 2 has no driving force and remains absent. D12 therefore cannot
    // change the independent binary result j0=-D01*grad x0, j1=-j0. Yet it sets
    // the matrix normalization: c*D_min*drive=1e-330 rounds to zero while the
    // dimensionless solved amplitudes are +/-1e300. Their final flux is normal.
    assert!((mass[0] / -gradient - 1.0).abs() < 3e-14);
    assert!((mass[1] / gradient - 1.0).abs() < 3e-14);
    assert_eq!(mass[2], 0.0);
    assert!(residual < 1e-43);
}
