//! Prescribed thin-solid reference changes preserve memory and physical work.

use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

const G0: [f64; 9] = [1.08, 0.04, 0.16, 0.0, 0.95, -0.09, 0.025, -0.015, 1.03];
const G1: [f64; 9] = [1.13, 0.07, 0.19, 0.02, 0.98, -0.06, 0.04, -0.02, 1.09];
const A: [f64; 4] = [1.27, 0.13, -0.04, 0.94];

fn law(g: [f64; 9], kappa: f64) -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6e4,
                c01: 4e4,
                kappa,
            },
            4.0,
        )
        .unwrap(),
        vec![
            FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 1.6e5 },
                tau_s: 0.12,
            },
            FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 0.9e5 },
                tau_s: 0.3,
            },
        ],
    )
    .unwrap()
    .with_reference_deformation(ReferenceDeformation::try_new(g).unwrap())
}

fn close(a: f64, b: f64, scale: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= 5e-10 * scale,
        "{a:e} versus {b:e}; scale={scale:e}"
    );
}

fn settings() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 1e-12,
        ..PlaneStressSettings::default()
    }
}

fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|i| (0..3).map(|k| a[3 * (i / 3) + k] * b[3 * k + i % 3]).sum())
}

#[test]
fn finite_event_reuses_native_transactions_without_erasing_memory_or_creating_heat() {
    let old = law(G0, 8e5);
    let mut target = law(G1, 8e5);
    target.branches[0].tau_s *= 0.4;
    let accepted = old
        .step_plane_stress_response(&old.initial_state(), &A, 0.09, None, settings())
        .unwrap()
        .material
        .state;
    let snapshot = accepted.clone();
    assert!(accepted.dissipated > 0.0);
    let event = target
        .change_plane_stress_reference_from(&old, &accepted, None, settings())
        .unwrap();
    let fixed = target.change_reference_from(&old, &accepted).unwrap();
    let direct = target
        .step_plane_stress_response(&fixed.state, &A, 0.0, None, settings())
        .unwrap();
    assert_eq!(event.material, direct.material);
    assert_eq!(event.reference_work_pa, fixed.reference_work_pa);
    assert_eq!(event.constraint_reference_work_pa, 0.0);
    assert_eq!(accepted, snapshot);
    assert_eq!(
        event.material.state.inverse_viscous_cauchy_green,
        accepted.inverse_viscous_cauchy_green
    );
    assert_eq!(
        event.material.state.dissipated.to_bits(),
        accepted.dissipated.to_bits()
    );
    assert_eq!(event.material.dissipated_pa, 0.0);
    for i in [0, 1, 3, 4] {
        assert_eq!(
            event.material.state.deformation_gradient[i],
            accepted.deformation_gradient[i]
        );
    }
    for i in [2, 5, 8] {
        close(event.material.first_piola_pa[i], 0.0, 1e6);
    }
    close(
        event.material.state.stored - accepted.stored,
        event.reference_work_pa + event.material.loading_work_pa,
        1e6,
    );
    close(
        event.material.state.stored - accepted.stored,
        event.reference_work_pa + event.transverse_endpoint_work_pa - event.endpoint_work_defect_pa,
        1e6,
    );
    close(
        event.material.state.work + event.material.state.reference_work,
        event.material.state.stored + event.material.state.dissipated,
        1e6,
    );
    assert!(event.endpoint_work_defect_pa.abs() > 1.0);
    let resumed = target
        .step_plane_stress_response(&event.material.state, &A, 0.04, None, settings())
        .unwrap();
    let erased = target
        .step_plane_stress_response(&target.initial_state(), &A, 0.04, None, settings())
        .unwrap();
    assert!(resumed.material.dissipated_pa > 0.0);
    assert_ne!(
        resumed.material.state.inverse_viscous_cauchy_green,
        event.material.state.inverse_viscous_cauchy_green
    );
    assert!((resumed.plane_stress.piola[0] - erased.plane_stress.piola[0]).abs() > 100.0);
}

#[test]
fn exact_event_keeps_pressure_reference_actuation_separate_and_never_evaluates_bulk_penalty() {
    let mut results = Vec::new();
    for kappa in [8e5, 1e300] {
        let old = law(G0, kappa);
        let target = law(G1, kappa);
        let accepted = old
            .step_incompressible_plane_stress_response(
                &old.initial_remainder_state(),
                &A,
                0.09,
                None,
                settings(),
            )
            .unwrap()
            .material
            .state;
        let snapshot = accepted.clone();
        let event = target
            .change_incompressible_plane_stress_reference_from(&old, &accepted, None, settings())
            .unwrap();
        let fixed = target
            .change_remainder_reference_from(&old, &accepted)
            .unwrap();
        let direct = target
            .step_incompressible_plane_stress_response(&fixed.state, &A, 0.0, None, settings())
            .unwrap();
        assert_eq!(event.material, direct.material);
        assert_eq!(event.reference_work_pa, fixed.reference_work_pa);
        assert_eq!(accepted, snapshot);
        assert_eq!(
            event.material.state.inverse_viscous_cauchy_green,
            accepted.inverse_viscous_cauchy_green
        );
        assert_eq!(
            event.material.state.dissipated.to_bits(),
            accepted.dissipated.to_bits()
        );
        assert_eq!(event.material.dissipated_pa, 0.0);
        let f = event.plane_stress.deformation_gradient;
        close(determinant(&f), determinant(&G1), 1.0);
        assert!(event.plane_stress.pressure_pa.abs() > 100.0);
        let delta_j = determinant(&G1) - determinant(&G0);
        close(
            event.constraint_reference_work_pa,
            event.plane_stress.pressure_pa * delta_j,
            1e6,
        );
        // Independently recover pressure actuation from the traction split and
        // the actual thickness motion, rather than from the response account.
        let transverse_remainder_work: f64 = [2, 5, 8]
            .iter()
            .map(|&i| event.material.first_piola_pa[i] * (f[i] - accepted.deformation_gradient[i]))
            .sum();
        close(
            transverse_remainder_work - event.transverse_endpoint_work_pa,
            event.constraint_reference_work_pa,
            1e6,
        );
        for i in [2, 5, 8] {
            close(event.plane_stress.full_piola[i], 0.0, 1e6);
        }
        close(
            event.material.state.stored - accepted.stored,
            event.reference_work_pa
                + event.constraint_reference_work_pa
                + event.transverse_endpoint_work_pa
                - event.endpoint_work_defect_pa,
            1e6,
        );
        close(
            event.material.state.work + event.material.state.reference_work,
            event.material.state.stored + event.material.state.dissipated,
            1e6,
        );
        results.push(event);
    }
    assert_eq!(results[0].material, results[1].material);
    assert_eq!(
        results[0].plane_stress.full_piola,
        results[1].plane_stress.full_piola
    );
    assert_eq!(
        results[0].constraint_reference_work_pa,
        results[1].constraint_reference_work_pa
    );
    assert_eq!(
        results[0].endpoint_work_defect_pa,
        results[1].endpoint_work_defect_pa
    );
}

#[test]
fn exact_pressure_reference_term_is_first_order_but_corrected_endpoint_defect_is_second_order() {
    let old = law(G0, 8e5);
    let accepted = old
        .step_incompressible_plane_stress_response(
            &old.initial_remainder_state(),
            &A,
            0.09,
            None,
            settings(),
        )
        .unwrap()
        .material
        .state;
    assert!(accepted.dissipated > 0.0);
    let mut corrected = Vec::new();
    let mut uncorrected = Vec::new();
    for h in [0.002, 0.001, 0.0005] {
        let target = law(G0.map(|x| (1.0 + h) * x), 8e5);
        let event = target
            .change_incompressible_plane_stress_reference_from(&old, &accepted, None, settings())
            .unwrap();
        assert!(event.plane_stress.pressure_pa.abs() > 100.0);
        corrected.push(event.endpoint_work_defect_pa.abs());
        uncorrected
            .push((event.transverse_endpoint_work_pa - event.material.loading_work_pa).abs());
        assert!(
            event.constraint_reference_work_pa.abs() > 10.0 * event.endpoint_work_defect_pa.abs()
        );
        assert_eq!(event.material.dissipated_pa, 0.0);
    }
    for pair in corrected.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.7..4.3).contains(&ratio),
            "corrected defect {corrected:?}, ratio {ratio}"
        );
    }
    for pair in uncorrected.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (1.8..2.2).contains(&ratio),
            "uncorrected work {uncorrected:?}, ratio {ratio}"
        );
    }
}

#[test]
fn event_stress_and_work_are_objective_for_rotated_membranes_with_nonsymmetric_references() {
    let rotation = [0.8, -0.6, 0.0, 0.6, 0.8, 0.0, 0.0, 0.0, 1.0];
    let rotated_a = [
        0.8 * A[0] - 0.6 * A[2],
        0.8 * A[1] - 0.6 * A[3],
        0.6 * A[0] + 0.8 * A[2],
        0.6 * A[1] + 0.8 * A[3],
    ];
    let old = law(G0, 8e5);
    let target = law(G1, 8e5);
    let full = old
        .step_plane_stress_response(&old.initial_state(), &A, 0.09, None, settings())
        .unwrap()
        .material
        .state;
    let rotated_full = old
        .step_plane_stress_response(&old.initial_state(), &rotated_a, 0.09, None, settings())
        .unwrap()
        .material
        .state;
    let exact = old
        .step_incompressible_plane_stress_response(
            &old.initial_remainder_state(),
            &A,
            0.09,
            None,
            settings(),
        )
        .unwrap()
        .material
        .state;
    let rotated_exact = old
        .step_incompressible_plane_stress_response(
            &old.initial_remainder_state(),
            &rotated_a,
            0.09,
            None,
            settings(),
        )
        .unwrap()
        .material
        .state;
    let x = target
        .change_plane_stress_reference_from(&old, &full, None, settings())
        .unwrap();
    let rx = target
        .change_plane_stress_reference_from(&old, &rotated_full, None, settings())
        .unwrap();
    let y = target
        .change_incompressible_plane_stress_reference_from(&old, &exact, None, settings())
        .unwrap();
    let ry = target
        .change_incompressible_plane_stress_reference_from(&old, &rotated_exact, None, settings())
        .unwrap();
    for (stress, rotated_stress, f, rotated_f) in [
        (
            x.material.first_piola_pa,
            rx.material.first_piola_pa,
            x.material.state.deformation_gradient,
            rx.material.state.deformation_gradient,
        ),
        (
            y.plane_stress.full_piola,
            ry.plane_stress.full_piola,
            y.material.state.deformation_gradient,
            ry.material.state.deformation_gradient,
        ),
    ] {
        let expected_stress = multiply(&rotation, &stress);
        let expected_f = multiply(&rotation, &f);
        for i in 0..9 {
            close(rotated_stress[i], expected_stress[i], 1e6);
            close(rotated_f[i], expected_f[i], 1.0);
        }
    }
    for (value, rotated) in [
        (x.reference_work_pa, rx.reference_work_pa),
        (x.endpoint_work_defect_pa, rx.endpoint_work_defect_pa),
        (y.reference_work_pa, ry.reference_work_pa),
        (
            y.constraint_reference_work_pa,
            ry.constraint_reference_work_pa,
        ),
        (y.endpoint_work_defect_pa, ry.endpoint_work_defect_pa),
        (y.plane_stress.pressure_pa, ry.plane_stress.pressure_pa),
    ] {
        close(value, rotated, 1e6);
    }
    for (plain, rotated) in [
        (
            &x.material.state.inverse_viscous_cauchy_green,
            &rx.material.state.inverse_viscous_cauchy_green,
        ),
        (
            &y.material.state.inverse_viscous_cauchy_green,
            &ry.material.state.inverse_viscous_cauchy_green,
        ),
    ] {
        for (a, b) in plain.iter().flatten().zip(rotated.iter().flatten()) {
            close(*a, *b, 1.0);
        }
    }
}

#[test]
fn reverse_and_kinetic_only_events_retain_memory_loss_and_total_energy_balance() {
    let old = law(G0, 8e5);
    let target = law(G1, 8e5);
    let accepted = old
        .step_plane_stress_response(&old.initial_state(), &A, 0.09, None, settings())
        .unwrap()
        .material
        .state;
    let remainder = old
        .step_incompressible_plane_stress_response(
            &old.initial_remainder_state(),
            &A,
            0.09,
            None,
            settings(),
        )
        .unwrap()
        .material
        .state;
    let forward = target
        .change_plane_stress_reference_from(&old, &accepted, None, settings())
        .unwrap();
    let reverse = old
        .change_plane_stress_reference_from(&target, &forward.material.state, None, settings())
        .unwrap();
    let rf = target
        .change_incompressible_plane_stress_reference_from(&old, &remainder, None, settings())
        .unwrap();
    let rr = old
        .change_incompressible_plane_stress_reference_from(
            &target,
            &rf.material.state,
            None,
            settings(),
        )
        .unwrap();
    for (before_f, after_f, before_stored, after_stored) in [
        (
            accepted.deformation_gradient,
            reverse.material.state.deformation_gradient,
            accepted.stored,
            reverse.material.state.stored,
        ),
        (
            remainder.deformation_gradient,
            rr.material.state.deformation_gradient,
            remainder.stored,
            rr.material.state.stored,
        ),
    ] {
        for i in 0..9 {
            close(before_f[i], after_f[i], 1.0);
        }
        close(before_stored, after_stored, 1e6);
    }
    assert_eq!(
        reverse.material.state.inverse_viscous_cauchy_green,
        accepted.inverse_viscous_cauchy_green
    );
    assert_eq!(
        rr.material.state.inverse_viscous_cauchy_green,
        remainder.inverse_viscous_cauchy_green
    );
    assert_eq!(reverse.material.state.dissipated, accepted.dissipated);
    assert_eq!(rr.material.state.dissipated, remainder.dissipated);
    close(
        forward.reference_work_pa
            + forward.material.loading_work_pa
            + reverse.reference_work_pa
            + reverse.material.loading_work_pa,
        0.0,
        1e6,
    );
    close(
        rf.reference_work_pa
            + rf.material.loading_work_pa
            + rr.reference_work_pa
            + rr.material.loading_work_pa,
        0.0,
        1e6,
    );
    let mut kinetic = old.clone();
    for branch in &mut kinetic.branches {
        branch.tau_s *= 2.0;
    }
    let same = kinetic
        .change_plane_stress_reference_from(
            &old,
            &accepted,
            Some([
                accepted.deformation_gradient[2],
                accepted.deformation_gradient[5],
                accepted.deformation_gradient[8],
            ]),
            settings(),
        )
        .unwrap();
    let same_r = kinetic
        .change_incompressible_plane_stress_reference_from(&old, &remainder, None, settings())
        .unwrap();
    assert_eq!(same.reference_work_pa, 0.0);
    assert_eq!(same_r.reference_work_pa, 0.0);
    assert_eq!(same_r.constraint_reference_work_pa, 0.0);
    assert_eq!(
        same.material.state.inverse_viscous_cauchy_green,
        accepted.inverse_viscous_cauchy_green
    );
    assert_eq!(
        same_r.material.state.inverse_viscous_cauchy_green,
        remainder.inverse_viscous_cauchy_green
    );
    assert_eq!(same.material.dissipated_pa, 0.0);
    assert_eq!(same_r.material.dissipated_pa, 0.0);
    close(same.material.loading_work_pa, 0.0, 1e6);
    close(same_r.material.loading_work_pa, 0.0, 1e6);
}

#[test]
fn incompatible_laws_non_thin_histories_and_failed_transverse_proposals_are_atomic() {
    let old = law(G0, 8e5);
    let target = law(G1, 8e5);
    let accepted = old
        .step_plane_stress_response(&old.initial_state(), &A, 0.09, None, settings())
        .unwrap()
        .material
        .state;
    let remainder = old
        .step_incompressible_plane_stress_response(
            &old.initial_remainder_state(),
            &A,
            0.09,
            None,
            settings(),
        )
        .unwrap()
        .material
        .state;
    let snapshot = (accepted.clone(), remainder.clone());
    let mut reordered = target.clone();
    reordered.branches.reverse();
    let mut domain = target.clone();
    domain.equilibrium.stretch_limit = 5.0;
    let mut energy = target.clone();
    energy.branches[0].energy = BranchEnergy::NeoHookean { mu_pa: 2e5 };
    let mut kinetics = target.clone();
    kinetics.branches[0].tau_s = f64::NAN;
    for invalid in [reordered, domain, energy, kinetics] {
        assert!(
            invalid
                .change_plane_stress_reference_from(&old, &accepted, None, settings())
                .is_err()
        );
        assert!(
            invalid
                .change_incompressible_plane_stress_reference_from(
                    &old,
                    &remainder,
                    None,
                    settings()
                )
                .is_err()
        );
    }
    // The fixed-F transaction itself succeeds; failures below occur afterward,
    // when the transverse seed or solve cannot provide the complete endpoint.
    assert!(target.change_reference_from(&old, &accepted).is_ok());
    assert!(
        target
            .change_remainder_reference_from(&old, &remainder)
            .is_ok()
    );
    assert!(
        target
            .change_plane_stress_reference_from(
                &old,
                &accepted,
                Some([f64::NAN, 0.0, 1.0]),
                settings(),
            )
            .is_err()
    );
    assert!(
        target
            .change_incompressible_plane_stress_reference_from(
                &old,
                &remainder,
                Some([0.0, f64::NAN]),
                settings(),
            )
            .is_err()
    );
    let no_updates = PlaneStressSettings {
        max_iterations: 0,
        ..settings()
    };
    assert!(
        target
            .change_plane_stress_reference_from(&old, &accepted, None, no_updates)
            .is_err()
    );
    assert!(
        target
            .change_incompressible_plane_stress_reference_from(&old, &remainder, None, no_updates)
            .is_err()
    );
    // An admitted general 3-D history is not silently projected to a membrane.
    let full_3d = old
        .step_response(&old.initial_state(), &G0, 0.0)
        .unwrap()
        .state;
    assert!(
        target
            .change_plane_stress_reference_from(&old, &full_3d, None, settings())
            .is_err()
    );
    let mut off_volume_f = remainder.deformation_gradient;
    off_volume_f[8] *= 1.02;
    let off_volume = old
        .step_remainder_response(&remainder, &off_volume_f, 0.0)
        .unwrap()
        .state;
    assert!(
        target
            .change_incompressible_plane_stress_reference_from(&old, &off_volume, None, settings())
            .is_err()
    );
    assert_eq!((accepted, remainder), snapshot);
}
