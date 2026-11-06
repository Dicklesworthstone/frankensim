use super::*;
use crate::solution::BinaryRegularSolutionV1;
use fs_qty::MolarMass;

fn ids(n: usize) -> Vec<SpeciesId> {
    (0..n)
        .map(|i| SpeciesId::new(format!("S{i:03}")).unwrap())
        .collect()
}
fn matrix(n: usize, pairs: &[(usize, usize, f64)]) -> Vec<f64> {
    let mut values = vec![0.0; n * n];
    for &(i, j, value) in pairs {
        values[i * n + j] = value;
        values[j * n + i] = value;
    }
    values
}
fn model(
    n: usize,
    a: &[(usize, usize, f64)],
    b: &[(usize, usize, f64)],
) -> MulticomponentRegularSolutionV1 {
    MulticomponentRegularSolutionV1::new(
        ids(n),
        matrix(n, a),
        matrix(n, b),
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap()
}
fn evaluate(
    law: &MulticomponentRegularSolutionV1,
    t: f64,
    x: &[f64],
) -> MulticomponentSolutionStateV1 {
    law.evaluate(Temperature::new(t), x, || false).unwrap()
}
fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g0_binary_reduction_retains_original_thermodynamics_and_trace_limits() {
    for chi in [-4.0, 0.0, 1.5, 2.0, 3.0] {
        let law = model(2, &[(0, 1, chi - 3.0)], &[(0, 1, 900.0)]);
        let binary = BinaryRegularSolutionV1::new(
            [ids(2)[0].clone(), ids(2)[1].clone()],
            chi - 3.0,
            900.0,
            Temperature::new(100.0),
            Temperature::new(1000.0),
        )
        .unwrap();
        for x in [
            [0.0, 1.0],
            [1e-100, 1.0],
            [0.1, 0.9],
            [0.5, 0.5],
            [0.9, 0.1],
            [1.0, 1e-100],
            [1.0, 0.0],
        ] {
            let state = evaluate(&law, 300.0, &x);
            let reference = binary
                .evaluate_composition(Temperature::new(300.0), x)
                .unwrap();
            near(
                state.mixing_gibbs_j_per_mol(),
                reference.mixing_gibbs_j_per_mol(),
                3e-14,
                0.0,
            );
            near(
                state.mixing_enthalpy_j_per_mol(),
                reference.mixing_enthalpy_j_per_mol(),
                3e-14,
                0.0,
            );
            near(
                state.mixing_entropy_j_per_mol_k(),
                reference.mixing_entropy_j_per_mol_k(),
                3e-14,
                0.0,
            );
            for i in 0..2 {
                near(
                    state.log_activity_coefficients()[i],
                    reference.log_activity_coefficients()[i],
                    3e-14,
                    0.0,
                );
                match (
                    state.reduced_chemical_potentials()[i],
                    reference.reduced_chemical_potentials()[i],
                ) {
                    (Some(a), Some(b)) => near(a, b, 3e-14, 1e-15),
                    (None, None) => {}
                    _ => panic!("pure/trace species identity changed"),
                }
            }
            let force = state.driving_forces(&[0.2, -0.2], || false).unwrap();
            for (a, b) in force
                .values()
                .iter()
                .zip(reference.driving_forces(0.2).unwrap())
            {
                near(*a, b, 2e-14, 1e-16);
            }
        }
    }
}

#[test]
fn g1_ternary_chemical_potentials_are_derivatives_of_extensive_gibbs_energy() {
    let law = model(
        3,
        &[(0, 1, 0.7), (0, 2, -0.4), (1, 2, 1.2)],
        &[(0, 1, 450.0), (0, 2, -180.0), (1, 2, 120.0)],
    );
    let amounts = [0.2, 0.3, 0.5];
    let t = 350.0;
    let state = evaluate(&law, t, &amounts);
    let total_g = |amounts: [f64; 3]| {
        let n: f64 = amounts.iter().sum();
        let x = amounts.map(|v| v / n);
        n * evaluate(&law, t, &x).mixing_gibbs_j_per_mol()
    };
    let total_h = |amounts: [f64; 3]| {
        let n: f64 = amounts.iter().sum();
        let x = amounts.map(|v| v / n);
        n * evaluate(&law, t, &x).mixing_enthalpy_j_per_mol()
    };
    for i in 0..3 {
        let step = 1e-5;
        let derivative = |f: &dyn Fn([f64; 3]) -> f64| {
            let sample = |offset: f64| {
                let mut v = amounts;
                v[i] += offset * step;
                f(v)
            };
            (sample(-2.0) - 8.0 * sample(-1.0) + 8.0 * sample(1.0) - sample(2.0)) / (12.0 * step)
        };
        near(
            derivative(&total_g),
            R * t * state.reduced_chemical_potentials()[i].unwrap(),
            3e-10,
            1e-7,
        );
        near(
            derivative(&total_h),
            state.partial_molar_mixing_enthalpies_j_per_mol()[i],
            3e-10,
            1e-8,
        );
    }
    near(
        amounts
            .iter()
            .zip(state.reduced_chemical_potentials())
            .map(|(x, mu)| x * mu.unwrap())
            .sum::<f64>()
            * R
            * t,
        state.mixing_gibbs_j_per_mol(),
        3e-15,
        1e-12,
    );
    near(
        amounts
            .iter()
            .zip(state.partial_molar_mixing_enthalpies_j_per_mol())
            .map(|(x, h)| x * h)
            .sum(),
        state.mixing_enthalpy_j_per_mol(),
        3e-15,
        1e-12,
    );
}

#[test]
fn g1_temperature_derivatives_close_mixing_enthalpy_and_entropy() {
    let law = model(
        4,
        &[(0, 1, 0.8), (0, 3, -0.9), (2, 3, 1.0)],
        &[(0, 1, 200.0), (1, 2, 800.0), (1, 3, -700.0)],
    );
    let x = [0.1, 0.2, 0.3, 0.4];
    let (t, dt) = (400.0, 0.01);
    let state = evaluate(&law, t, &x);
    let lo = evaluate(&law, t - dt, &x);
    let hi = evaluate(&law, t + dt, &x);
    near(
        -(hi.mixing_gibbs_j_per_mol() - lo.mixing_gibbs_j_per_mol()) / (2.0 * dt),
        state.mixing_entropy_j_per_mol_k(),
        2e-10,
        1e-10,
    );
    near(
        state.mixing_enthalpy_j_per_mol() - t * state.mixing_entropy_j_per_mol_k(),
        state.mixing_gibbs_j_per_mol(),
        3e-15,
        1e-12,
    );
    for i in 0..4 {
        let h =
            -R * t * t * (hi.log_activity_coefficients()[i] - lo.log_activity_coefficients()[i])
                / (2.0 * dt);
        near(
            h,
            state.partial_molar_mixing_enthalpies_j_per_mol()[i],
            2e-9,
            1e-8,
        );
    }
}

#[test]
fn g1_complete_cross_component_forces_obey_gibbs_duhem_and_hessian() {
    let law = model(3, &[(0, 1, 3.0), (0, 2, -0.5), (1, 2, 0.6)], &[]);
    let x = [0.2, 0.3, 0.5];
    let gradient = [0.4, -0.1, -0.3];
    let h = 1e-5;
    let low = evaluate(
        &law,
        300.0,
        &std::array::from_fn::<_, 3, _>(|i| x[i] - h * gradient[i]),
    );
    let high = evaluate(
        &law,
        300.0,
        &std::array::from_fn::<_, 3, _>(|i| x[i] + h * gradient[i]),
    );
    let state = evaluate(&law, 300.0, &x);
    let forces = state.driving_forces(&gradient, || false).unwrap();
    let mut curvature = 0.0;
    for i in 0..3 {
        let derivative = (high.reduced_chemical_potentials()[i].unwrap()
            - low.reduced_chemical_potentials()[i].unwrap())
            / (2.0 * h);
        near(forces.values()[i], x[i] * derivative, 2e-9, 1e-10);
        curvature += gradient[i] * derivative;
    }
    near(forces.values().iter().sum(), 0.0, 0.0, 2e-16);
    near(
        state.reduced_curvature(&gradient).unwrap(),
        curvature,
        2e-9,
        1e-10,
    );
    // The third component has zero concentration gradient, yet responds to
    // changes in the other two activities: independent Fick laws would miss it.
    let cross = state.driving_forces(&[0.1, -0.1, 0.0], || false).unwrap();
    assert!(cross.values()[2].abs() > 0.01);
}

#[test]
fn g1_cosolvent_and_temperature_change_the_actual_demixing_curvature() {
    let law = model(3, &[], &[(0, 1, 900.0)]);
    let binary_feed = evaluate(&law, 300.0, &[0.5, 0.5, 0.0]);
    near(
        binary_feed.reduced_curvature(&[1.0, -1.0, 0.0]).unwrap(),
        -2.0,
        0.0,
        1e-14,
    );
    let with_solvent = evaluate(&law, 300.0, &[0.1, 0.1, 0.8]);
    near(
        with_solvent.reduced_curvature(&[1.0, -1.0, 0.0]).unwrap(),
        14.0,
        0.0,
        1e-14,
    );
    let warmed = evaluate(&law, 600.0, &[0.5, 0.5, 0.0]);
    near(
        warmed.reduced_curvature(&[1.0, -1.0, 0.0]).unwrap(),
        1.0,
        0.0,
        1e-14,
    );
    // One negative direction is sufficient to demonstrate homogeneous
    // instability; this does not assign any equilibrium phase amounts.
    let ternary = model(3, &[(0, 1, 3.5), (0, 2, 3.5), (1, 2, 3.5)], &[]);
    let symmetric = evaluate(&ternary, 300.0, &[1.0 / 3.0; 3]);
    near(
        symmetric.reduced_curvature(&[1.0, -1.0, 0.0]).unwrap(),
        -1.0,
        0.0,
        2e-14,
    );
}

#[test]
fn g1_both_maxwell_stefan_frames_consume_the_same_ternary_liquid_forces() {
    let law = model(3, &[(0, 1, 1.3), (0, 2, -0.4), (1, 2, 0.8)], &[]);
    let x = [0.2, 0.3, 0.5];
    let gradient = [0.1, -0.3, 0.2];
    let state = evaluate(&law, 340.0, &x);
    let masses = [0.018, 0.046, 0.092];
    let transport =
        MaxwellStefanV1::new(ids(3).into_iter().zip(masses.map(MolarMass::new)).collect()).unwrap();
    let diffusivity = 8e-10;
    let d = matrix(
        3,
        &[
            (0, 1, diffusivity),
            (0, 2, diffusivity),
            (1, 2, diffusivity),
        ],
    );
    let concentration = AmountConcentration::new(18000.0);
    let forces = state.driving_forces(&gradient, || false).unwrap();
    let mut workspace = transport.workspace().unwrap();
    let molar = state
        .molar_diffusion(
            &transport,
            &gradient,
            concentration,
            &d,
            MaxwellStefanMolarSolveV1::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    let molar_flux: [f64; 3] = molar.molar_fluxes().try_into().unwrap();
    for i in 0..3 {
        near(
            molar_flux[i],
            -concentration.value() * diffusivity * forces.values()[i],
            3e-14,
            1e-20,
        );
    }
    let physical_entropy: f64 = -R
        * (0..3)
            .map(|i| molar_flux[i] * forces.values()[i] / x[i])
            .sum::<f64>();
    near(
        molar.drag_entropy_production_w_m3_k().unwrap(),
        physical_entropy,
        3e-14,
        1e-20,
    );
    assert!(physical_entropy > 0.0);
    let mass = state
        .mass_diffusion(
            &transport,
            &gradient,
            concentration,
            &d,
            MaxwellStefanSolveV1::default(),
            &mut workspace,
            || false,
        )
        .unwrap();
    let mean_mass: f64 = (0..3).map(|i| x[i] * masses[i]).sum();
    let frame_shift: f64 = (0..3).map(|i| molar_flux[i] * masses[i]).sum::<f64>() / mean_mass;
    for i in 0..3 {
        near(
            mass.molar_fluxes()[i],
            molar_flux[i] - x[i] * frame_shift,
            3e-14,
            1e-20,
        );
    }
    near(mass.mass_fluxes().iter().sum(), 0.0, 0.0, 1e-20);
    near(
        mass.drag_entropy_production_w_m3_k().unwrap(),
        physical_entropy,
        3e-14,
        1e-20,
    );
}

#[test]
fn g0_pure_absent_and_multiple_trace_components_are_not_floored() {
    let law = model(
        3,
        &[(0, 1, 3.0), (0, 2, 2.0), (1, 2, -0.4)],
        &[(0, 1, 300.0)],
    );
    let pure = evaluate(&law, 300.0, &[1.0, 0.0, 0.0]);
    assert_eq!(pure.mixing_gibbs_j_per_mol(), 0.0);
    assert_eq!(pure.reduced_chemical_potentials(), &[Some(0.0), None, None]);
    assert_eq!(pure.log_activity_coefficients(), &[0.0, 4.0, 2.0]);
    let trace = evaluate(&law, 300.0, &[1.0, 1e-100, 2e-100]);
    assert_eq!(trace.mole_fractions(), &[1.0, 1e-100, 2e-100]);
    near(
        trace.reduced_chemical_potentials()[0].unwrap(),
        -3e-100,
        1e-14,
        0.0,
    );
    // Carrier gamma = (x_B+x_C)*(chi_AB*x_B+chi_AC*x_C) - chi_BC*x_B*x_C.
    near(trace.log_activity_coefficients()[0], 24.8e-200, 2e-14, 0.0);
    assert!(trace.reduced_chemical_potentials()[1].unwrap().is_finite());
    assert!(trace.reduced_curvature(&[0.0, 1.0, -1.0]).unwrap() > 1e99);
    assert!(pure.reduced_curvature(&[-1.0, 1.0, 0.0]).is_err());
    assert_eq!(pure.reduced_curvature(&[0.0, 0.0, 0.0]).unwrap(), 0.0);
    // A finite one-sided diffusion force does not assert a finite absent mu.
    assert_eq!(
        pure.driving_forces(&[-0.2, 0.1, 0.1], || false)
            .unwrap()
            .values(),
        &[-0.2, 0.1, 0.1]
    );
}

#[test]
fn g3_species_permutation_reorders_every_original_pair_without_changing_physics() {
    let original = model(
        3,
        &[(0, 1, 1.1), (0, 2, -0.7), (1, 2, 2.4)],
        &[(0, 2, 200.0), (1, 2, -450.0)],
    );
    let permutation = [2, 0, 1];
    let reindex = |v: &[f64]| {
        permutation
            .iter()
            .flat_map(|&i| permutation.iter().map(move |&j| v[3 * i + j]))
            .collect()
    };
    let offered = MulticomponentRegularSolutionV1::new(
        permutation.iter().map(|&i| ids(3)[i].clone()).collect(),
        reindex(original.chi_constants()),
        reindex(original.chi_inverse_temperature_k()),
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap();
    assert_eq!(original, offered);
    assert_eq!(
        evaluate(&original, 350.0, &[0.2, 0.3, 0.5]),
        evaluate(&offered, 350.0, &[0.2, 0.3, 0.5])
    );
}

#[test]
fn g0_single_species_and_maximum_species_are_real_supported_limits() {
    let single = model(1, &[], &[]);
    let state = evaluate(&single, 300.0, &[1.0]);
    assert_eq!(state.mixing_gibbs_j_per_mol(), 0.0);
    assert_eq!(
        state.driving_forces(&[0.0], || false).unwrap().values(),
        &[0.0]
    );
    let n = MAX_REGULAR_SOLUTION_SPECIES_V1;
    let maximum = model(n, &[], &[]);
    let x = vec![1.0 / n as f64; n];
    let state = evaluate(&maximum, 300.0, &x);
    near(
        state.mixing_gibbs_j_per_mol(),
        -R * 300.0 * (n as f64).ln(),
        3e-15,
        1e-12,
    );
    let mut gradient = vec![0.0; n];
    gradient[0] = 0.1;
    gradient[n - 1] = -0.1;
    assert_eq!(
        state.driving_forces(&gradient, || false).unwrap().values(),
        gradient
    );
    assert!(
        MulticomponentRegularSolutionV1::new(
            ids(n + 1),
            vec![0.0; (n + 1) * (n + 1)],
            vec![0.0; (n + 1) * (n + 1)],
            Temperature::new(100.0),
            Temperature::new(1000.0)
        )
        .is_err()
    );
}

#[test]
fn g0_malformed_pair_sources_compositions_and_tangents_refuse() {
    let new = |species, a, b, lo, hi| {
        MulticomponentRegularSolutionV1::new(
            species,
            a,
            b,
            Temperature::new(lo),
            Temperature::new(hi),
        )
    };
    assert!(new(vec![], vec![], vec![], 100.0, 1000.0).is_err());
    assert!(
        new(
            vec![ids(1)[0].clone(); 2],
            vec![0.0; 4],
            vec![0.0; 4],
            100.0,
            1000.0
        )
        .is_err()
    );
    assert!(new(ids(2), vec![0.0; 3], vec![0.0; 4], 100.0, 1000.0).is_err());
    for a in [
        vec![1.0, 0.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0, 0.0],
        vec![0.0, f64::NAN, f64::NAN, 0.0],
    ] {
        assert!(new(ids(2), a, vec![0.0; 4], 100.0, 1000.0).is_err());
    }
    assert!(new(ids(1), vec![0.0], vec![0.0], 0.0, 1000.0).is_err());
    let law = model(3, &[], &[]);
    for x in [
        vec![0.5, 0.5],
        vec![0.2, 0.3, 0.51],
        vec![-0.1, 0.5, 0.6],
        vec![0.0, f64::NAN, 1.0],
    ] {
        assert!(law.evaluate(Temperature::new(300.0), &x, || false).is_err());
    }
    for t in [99.0, 1001.0, f64::NAN] {
        assert_eq!(
            law.evaluate(Temperature::new(t), &[0.2, 0.3, 0.5], || false),
            Err(MulticomponentSolutionErrorV1::OutsideTemperatureDomain)
        );
    }
    let state = evaluate(&law, 300.0, &[0.2, 0.3, 0.5]);
    for tangent in [
        vec![1.0, -1.0],
        vec![0.1, 0.2, 0.3],
        vec![f64::INFINITY, -1.0, 1.0],
    ] {
        assert!(state.driving_forces(&tangent, || false).is_err());
        assert!(state.reduced_curvature(&tangent).is_err());
    }
}

#[test]
fn g0_original_species_are_required_by_diffusion_even_with_same_array_length() {
    let state = evaluate(&model(2, &[], &[]), 300.0, &[0.5, 0.5]);
    let transport = MaxwellStefanV1::new(vec![
        (
            SpeciesId::new("unrelated-a").unwrap(),
            MolarMass::new(0.018),
        ),
        (
            SpeciesId::new("unrelated-b").unwrap(),
            MolarMass::new(0.046),
        ),
    ])
    .unwrap();
    let mut workspace = transport.workspace().unwrap();
    assert!(matches!(
        state.molar_diffusion(
            &transport,
            &[0.1, -0.1],
            AmountConcentration::new(1000.0),
            &[0.0, 1e-9, 1e-9, 0.0],
            MaxwellStefanMolarSolveV1::default(),
            &mut workspace,
            || false
        ),
        Err(MulticomponentSolutionErrorV1::Invalid(_))
    ));
}

#[test]
fn g4_cancelled_evaluation_or_forces_never_publish_a_partial_response() {
    let law = model(4, &[(0, 1, 1.3), (1, 3, 0.4)], &[]);
    for trigger in [1, 2, 5, 9] {
        let mut calls = 0;
        assert_eq!(
            law.evaluate(Temperature::new(300.0), &[0.1, 0.2, 0.3, 0.4], || {
                calls += 1;
                calls == trigger
            }),
            Err(MulticomponentSolutionErrorV1::Cancelled)
        );
    }
    let state = evaluate(&law, 300.0, &[0.1, 0.2, 0.3, 0.4]);
    for trigger in [1, 3, 7, 10] {
        let mut calls = 0;
        assert_eq!(
            state.driving_forces(&[0.1, -0.1, 0.2, -0.2], || {
                calls += 1;
                calls == trigger
            }),
            Err(MulticomponentSolutionErrorV1::Cancelled)
        );
    }
    assert_eq!(state, evaluate(&law, 300.0, &[0.1, 0.2, 0.3, 0.4]));
}
