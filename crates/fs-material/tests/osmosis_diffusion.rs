//! Independent liquid friction and original-free-energy checks through public owners.
use fs_material::osmosis::{
    AqueousSolution,
    diffusion::{SolutionDiffusionError, SolutionDiffusionSolve, SolutionMaxwellStefan},
    nonideal::VirialSolutionModel,
};
use fs_qty::MolarMass;
use fs_thermochem::{
    SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    transport::maxwell_stefan::MaxwellStefanMolarSolveV1,
};
use std::sync::Arc;

const T: f64 = 298.15;
const VBAR: f64 = 18e-6;
fn model(reference: f64, b: Vec<f64>) -> Arc<VirialSolutionModel> {
    let n = (b.len() as f64).sqrt() as usize;
    Arc::new(VirialSolutionModel::new(T, VBAR, reference, vec![200.0; n], b).unwrap())
}
fn species(name: &str, mass: f64) -> (SpeciesId, MolarMass) {
    (SpeciesId::new(name).unwrap(), MolarMass::new(mass))
}
fn solution(model: &Arc<VirialSolutionModel>, c: &[f64]) -> AqueousSolution {
    AqueousSolution::nonideal(model.clone(), 1e-3, c.iter().map(|v| 1e-3 * v).collect()).unwrap()
}
fn controls() -> SolutionDiffusionSolve {
    SolutionDiffusionSolve {
        friction: MaxwellStefanMolarSolveV1::default(),
        entropy_tolerance_w_m3_k: 1e-10,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_binary_uses_dilute_solution_mu_and_exact_solvent_frame() {
    let source = model(1.0, vec![0.0]);
    // Solvent sorts AFTER the original solute: slot zero cannot be assumed.
    let law = SolutionMaxwellStefan::new(
        source.clone(),
        species("water", 0.018),
        vec![species("A", 0.06)],
        vec![0.0, 1e-9, 1e-9, 0.0],
    )
    .unwrap();
    let (left, right) = (solution(&source, &[100.0]), solution(&source, &[20.0]));
    let mut workspace = law.workspace().unwrap();
    let mut flux = [0.0];
    let report = law
        .flux_into(
            &left,
            &right,
            0.005,
            controls(),
            &mut workspace,
            &mut flux,
            || false,
        )
        .unwrap();
    let average = 80.0 / 5.0_f64.ln();
    let x0 = (1.0 / VBAR) / (1.0 / VBAR + average);
    let expected = 1e-9 * 80.0 / x0 / 0.005;
    close(flux[0], expected, 2e-18);
    assert_eq!(report.solvent_molar_flux_mol_m2_s, 0.0);
    close(
        report.chemical_entropy_production_w_m3_k.unwrap(),
        flux[0] * R * 5.0_f64.ln() / 0.005,
        1e-13,
    );
    assert!(report.chemical_entropy_production_w_m3_k.unwrap() > 0.0);
    law.flux_into(
        &right,
        &left,
        0.005,
        controls(),
        &mut workspace,
        &mut flux,
        || false,
    )
    .unwrap();
    close(flux[0], -expected, 2e-18);
}

#[test]
fn g1_g3_cross_interactions_drive_uphill_diffusion_with_positive_friction() {
    let b = vec![0.001, 0.008, 0.008, 0.1];
    let source = model(1.0, b);
    let d01 = 1e-9;
    let d02 = 3e-9;
    let d12 = 5e-10;
    let law = SolutionMaxwellStefan::new(
        source.clone(),
        species("water", 0.018),
        vec![species("Z", 0.08), species("A", 0.04)],
        vec![0.0, d01, d02, d01, 0.0, d12, d02, d12, 0.0],
    )
    .unwrap();
    let (left, right) = (
        solution(&source, &[20.0, 100.0]),
        solution(&source, &[21.0, 1.0]),
    );
    let mut workspace = law.workspace().unwrap();
    let mut flux = [0.0; 2];
    let report = law
        .flux_into(
            &left,
            &right,
            0.01,
            controls(),
            &mut workspace,
            &mut flux,
            || false,
        )
        .unwrap();
    // Independent elimination of the two solvent-fixed friction rows.
    let mean1 = 1.0 / (21.0_f64 / 20.0).ln();
    let mean2 = 99.0 / 100.0_f64.ln();
    let c = 1.0 / VBAR + mean1 + mean2;
    let (x0, x1, x2) = ((1.0 / VBAR) / c, mean1 / c, mean2 / c);
    let force = [
        (21.0_f64 / 20.0).ln() + 0.001 - 0.008 * 99.0,
        (1.0_f64 / 100.0).ln() + 0.008 - 0.1 * 99.0,
    ];
    let (a11, a12, a21, a22) = (
        x0 / d01 + x2 / d12,
        -x1 / d12,
        -x2 / d12,
        x0 / d02 + x1 / d12,
    );
    let rhs = [-mean1 * force[0] / 0.01, -mean2 * force[1] / 0.01];
    let determinant = a11 * a22 - a12 * a21;
    let expected = [
        (rhs[0] * a22 - a12 * rhs[1]) / determinant,
        (a11 * rhs[1] - a21 * rhs[0]) / determinant,
    ];
    for i in 0..2 {
        close(flux[i], expected[i], 3e-17);
    }
    assert!(
        flux[0] > 0.0,
        "A solute must move into its higher-concentration side"
    );
    let entropy = -R * (flux[0] * force[0] + flux[1] * force[1]) / 0.01;
    close(
        report.chemical_entropy_production_w_m3_k.unwrap(),
        entropy,
        1e-11,
    );
    close(
        report.friction_entropy_production_w_m3_k.unwrap(),
        entropy,
        1e-11,
    );
    assert!(entropy > 0.0);
}

#[test]
fn g3_original_species_order_and_concentration_reference_are_preserved() {
    let original_b = vec![0.001, 0.008, 0.008, 0.1];
    let run = |reference, permuted: bool| {
        let b = if permuted {
            vec![0.1, 0.008, 0.008, 0.001]
        } else {
            original_b.clone()
        };
        let source = model(reference, b);
        let (ids, binary, cl, cr) = if permuted {
            (
                vec![species("A", 0.04), species("Z", 0.08)],
                vec![0.0, 3e-9, 1e-9, 3e-9, 0.0, 5e-10, 1e-9, 5e-10, 0.0],
                vec![100.0, 20.0],
                vec![1.0, 21.0],
            )
        } else {
            (
                vec![species("Z", 0.08), species("A", 0.04)],
                vec![0.0, 1e-9, 3e-9, 1e-9, 0.0, 5e-10, 3e-9, 5e-10, 0.0],
                vec![20.0, 100.0],
                vec![21.0, 1.0],
            )
        };
        let law = SolutionMaxwellStefan::new(source.clone(), species("water", 0.018), ids, binary)
            .unwrap();
        let mut flux = [0.0; 2];
        law.flux_into(
            &solution(&source, &cl),
            &solution(&source, &cr),
            0.01,
            controls(),
            &mut law.workspace().unwrap(),
            &mut flux,
            || false,
        )
        .unwrap();
        flux
    };
    let expected = run(1.0, false);
    assert_eq!(expected, run(1e200, false));
    let permuted = run(1.0, true);
    close(permuted[0], expected[1], 2e-18);
    close(permuted[1], expected[0], 2e-18);
}

#[test]
fn g1_reference_cancelled_changes_follow_original_free_energy_and_resolve_adjacent_states() {
    let source = model(1.0, vec![0.001, 0.008, 0.008, 0.1]);
    let shifted = model(1e200, vec![0.001, 0.008, 0.008, 0.1]);
    let before = [20.0, 100.0];
    let after = [20.5, 99.5];
    let expected = source.mixing_free_energy_density_j_m3(&after).unwrap()
        - source.mixing_free_energy_density_j_m3(&before).unwrap();
    let actual = source
        .reference_cancelled_free_energy_change_j_m3(&before, &after)
        .unwrap();
    close(actual, expected, 1e-8);
    assert_eq!(
        actual,
        shifted
            .reference_cancelled_free_energy_change_j_m3(&before, &after)
            .unwrap()
    );
    let ideal = model(1e200, vec![0.0]);
    let a = 100.0_f64;
    let b = a.next_up();
    let mut difference = [0.0];
    ideal
        .chemical_potential_differences_into(&[a], &[b], &mut difference)
        .unwrap();
    let expected = R * T * ((b - a) / a).ln_1p();
    assert!(difference[0] > 0.0);
    close(difference[0], expected, 1e-27);
}

#[test]
fn g1_positive_dilute_chemical_forces_do_not_lose_subnormal_concentration_ratios() {
    let source = model(1.0, vec![0.0]);
    let left = f64::from_bits(2);
    let right = f64::from_bits(3);
    let mut difference = [0.0];
    source
        .chemical_potential_differences_into(&[left], &[right], &mut difference)
        .unwrap();
    close(difference[0], R * T * 1.5_f64.ln(), 1e-12);
    source
        .chemical_potential_differences_into(&[right], &[left], &mut difference)
        .unwrap();
    close(difference[0], -R * T * 1.5_f64.ln(), 1e-12);
}

#[test]
fn g1_ideal_fronts_keep_the_fick_limit_from_equal_to_absent_endpoints() {
    let source = model(1.0, vec![0.0]);
    let law = SolutionMaxwellStefan::new(
        source.clone(),
        species("water", 0.018),
        vec![species("A", 0.06)],
        vec![0.0, 1e-9, 1e-9, 0.0],
    )
    .unwrap();
    let mut workspace = law.workspace().unwrap();
    for (a, b) in [
        (100.0_f64, 100.0_f64),
        (100.0, 100.0_f64.next_down()),
        (100.0, 20.0),
        (100.0, 1e-12),
        (100.0, 1e-200),
        (100.0, 0.0),
        (1e-10, 1e-200),
        (0.0, 0.0),
    ] {
        // Unit volume retains adjacent concentration inputs exactly.
        let left = AqueousSolution::nonideal(source.clone(), 1.0, vec![a]).unwrap();
        let right = AqueousSolution::nonideal(source.clone(), 1.0, vec![b]).unwrap();
        let mut flux = [f64::NAN];
        let report = law
            .flux_into(
                &left,
                &right,
                0.005,
                controls(),
                &mut workspace,
                &mut flux,
                || false,
            )
            .unwrap();
        let mean = if b == 0.0 {
            0.0
        } else if a == b {
            a
        } else if a - b < b {
            (a - b) / ((a - b) / b).ln_1p()
        } else {
            (a - b) / (a.ln() - b.ln())
        };
        // Independent binary solvent-fixed reduction. The dilute limit is
        // D Delta(c)/L; the retained finite-solvent factor is explicit.
        let fick = 1e-9 * (a - b) / 0.005;
        let expected = fick * (1.0 + VBAR * mean);
        close(flux[0], expected, expected.abs() * 2e-12 + 1e-34);
        if a != b && b == 0.0 {
            assert_eq!(report.chemical_entropy_production_w_m3_k, None);
            assert_eq!(report.friction_entropy_production_w_m3_k, None);
        } else {
            assert!(report.chemical_entropy_production_w_m3_k.unwrap() >= 0.0);
            assert!(report.friction_entropy_production_w_m3_k.unwrap() >= 0.0);
        }
        if a <= 1e-10 {
            close(flux[0], fick, fick.abs() * 2e-12 + 1e-34);
        }
        law.flux_into(
            &right,
            &left,
            0.005,
            controls(),
            &mut workspace,
            &mut flux,
            || false,
        )
        .unwrap();
        close(flux[0], -expected, expected.abs() * 2e-12 + 1e-34);
    }
}

#[test]
fn g1_absent_spectators_and_counterdiffusing_fronts_never_require_seeded_solute() {
    let source = model(1.0, vec![0.001, 0.002, 0.002, 0.01]);
    let law = SolutionMaxwellStefan::new(
        source.clone(),
        species("water", 0.018),
        vec![species("A", 0.06), species("B", 0.04)],
        vec![0.0, 1e-9, 3e-9, 1e-9, 0.0, 5e-10, 3e-9, 5e-10, 0.0],
    )
    .unwrap();
    let mut workspace = law.workspace().unwrap();
    let mut flux = [f64::NAN; 2];
    let report = law
        .flux_into(
            &solution(&source, &[100.0, 0.0]),
            &solution(&source, &[20.0, 0.0]),
            0.005,
            controls(),
            &mut workspace,
            &mut flux,
            || false,
        )
        .unwrap();
    let mean = 80.0 / 5.0_f64.ln();
    let expected = 1e-9 * 80.0 * (1.0 + 0.001 * mean) * (1.0 + VBAR * mean) / 0.005;
    close(flux[0], expected, 2e-18);
    assert_eq!(flux[1], 0.0);
    assert!(report.chemical_entropy_production_w_m3_k.unwrap() > 0.0);
    assert!(report.friction_entropy_production_w_m3_k.unwrap() > 0.0);

    let report = law
        .flux_into(
            &solution(&source, &[100.0, 0.0]),
            &solution(&source, &[0.0, 20.0]),
            0.005,
            controls(),
            &mut workspace,
            &mut flux,
            || false,
        )
        .unwrap();
    // At exact disjoint fronts the log-mean solute mobilities vanish while
    // their ideal force products retain Delta(c). Solvent drag remains finite.
    close(flux[0], 1e-9 * 100.0 / 0.005, 2e-18);
    close(flux[1], -3e-9 * 20.0 / 0.005, 2e-18);
    assert_eq!(report.chemical_entropy_production_w_m3_k, None);
    assert_eq!(report.friction_entropy_production_w_m3_k, None);
    assert_eq!(report.solvent_molar_flux_mol_m2_s, 0.0);
}

#[test]
fn g0_missing_pair_source_mismatch_and_cancellation_refuse_atomically() {
    let source = model(1.0, vec![0.0]);
    for pairs in [vec![], vec![0.0, 1e-9, 2e-9, 0.0], vec![0.0, 0.0, 0.0, 0.0]] {
        assert!(
            SolutionMaxwellStefan::new(
                source.clone(),
                species("water", 0.018),
                vec![species("A", 0.06)],
                pairs
            )
            .is_err()
        );
    }
    let law = SolutionMaxwellStefan::new(
        source.clone(),
        species("water", 0.018),
        vec![species("A", 0.06)],
        vec![0.0, 1e-9, 1e-9, 0.0],
    )
    .unwrap();
    let left = solution(&source, &[100.0]);
    let mut workspace = law.workspace().unwrap();
    let mut output = [123.0];
    for right in [
        solution(&model(2.0, vec![0.0]), &[20.0]),
        AqueousSolution::nonideal(source.clone(), 2.0, vec![f64::from_bits(1)]).unwrap(),
    ] {
        assert!(
            law.flux_into(
                &left,
                &right,
                0.01,
                controls(),
                &mut workspace,
                &mut output,
                || false
            )
            .is_err()
        );
        assert_eq!(output, [123.0]);
    }
    let right = solution(&source, &[20.0]);
    let mut polls = 0;
    law.flux_into(
        &left,
        &right,
        0.01,
        controls(),
        &mut workspace,
        &mut output,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    output = [123.0];
    let mut called = 0;
    assert_eq!(
        law.flux_into(
            &left,
            &right,
            0.01,
            controls(),
            &mut workspace,
            &mut output,
            || {
                called += 1;
                called == polls
            }
        ),
        Err(SolutionDiffusionError::Cancelled)
    );
    assert_eq!(output, [123.0]);
}
