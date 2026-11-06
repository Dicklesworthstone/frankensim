//! G0/G3: native material histories evolve through retained plane-stress states.

use fs_material::hyper::incompressible_plane_stress::{
    incompressible_plane_stress, incompressible_plane_stress_with_seed,
};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

fn close(a: f64, b: f64, scale: f64, tolerance: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= tolerance * scale,
        "{a:e} versus {b:e}; scale {scale:e}"
    );
}

fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|i| (0..3).map(|k| a[3 * (i / 3) + k] * b[3 * k + i % 3]).sum())
}

fn cofactor(f: &[f64; 9]) -> [f64; 9] {
    [
        f[4] * f[8] - f[5] * f[7],
        f[5] * f[6] - f[3] * f[8],
        f[3] * f[7] - f[4] * f[6],
        f[2] * f[7] - f[1] * f[8],
        f[0] * f[8] - f[2] * f[6],
        f[1] * f[6] - f[0] * f[7],
        f[1] * f[5] - f[2] * f[4],
        f[2] * f[3] - f[0] * f[5],
        f[0] * f[4] - f[1] * f[3],
    ]
}

fn determinant(f: &[f64; 9]) -> f64 {
    let c = cofactor(f);
    f[0] * c[0] + f[1] * c[1] + f[2] * c[2]
}

fn settings() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 2e-11,
        ..PlaneStressSettings::default()
    }
}

fn reference() -> ReferenceDeformation {
    ReferenceDeformation::try_new([1.08, 0.04, 0.16, 0.0, 0.95, -0.09, 0.0, 0.0, 1.03]).unwrap()
}

fn law(kappa: f64) -> ViscoHyperelastic {
    let elastic = Hyperelastic::new(
        HyperelasticModel::MooneyRivlin {
            c10: 6e4,
            c01: 4e4,
            kappa,
        },
        4.0,
    )
    .unwrap();
    ViscoHyperelastic::new(
        elastic,
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 1.6e5 },
            tau_s: 0.12,
        }],
    )
    .unwrap()
    .with_reference_deformation(reference())
}

#[test]
fn g0_finite_loading_hold_and_unload_use_one_native_material_interval_per_step() {
    let law = law(8e5);
    let mut accepted = law.initial_state();
    let unloaded = [1.08, 0.04, 0.0, 0.95];
    let loaded = [1.24, 0.08, 0.0, 0.95];
    let mut held_stress = Vec::new();
    let mut held_thickness = Vec::new();
    for (a, dt) in [
        (unloaded, 0.0),
        (loaded, 0.03),
        (loaded, 0.06),
        (loaded, 0.06),
        (unloaded, 0.03),
        (unloaded, 0.06),
    ] {
        let before = accepted.clone();
        let response = law
            .step_plane_stress_response(&accepted, &a, dt, None, settings())
            .unwrap();
        assert_eq!(accepted, before);
        let original = law
            .step_response(&accepted, &response.plane_stress.deformation_gradient, dt)
            .unwrap();
        assert_eq!(response.material, original);
        for (i, j) in [0, 1, 3, 4].into_iter().enumerate() {
            close(
                response.plane_stress.piola[i],
                original.first_piola_pa[j],
                1e6,
                1e-12,
            );
        }
        for i in [2, 5, 8] {
            close(original.first_piola_pa[i], 0.0, 1e6, 5e-10);
        }
        close(
            response.plane_stress.energy_density_pa,
            original.incremental_potential_pa,
            1e6,
            1e-12,
        );
        close(
            original.loading_work_pa,
            original.state.stored - accepted.stored + original.dissipated_pa,
            1e6,
            2e-12,
        );
        assert!(original.dissipated_pa >= 0.0);
        assert!(original.state.dissipated >= accepted.dissipated);
        assert_eq!(original.state.reference_work, accepted.reference_work);
        if a == loaded {
            held_stress.push(response.plane_stress.piola[0]);
            held_thickness.push(response.plane_stress.deformation_gradient[8]);
        }
        accepted = response.material.state;
    }
    assert!(held_stress[2] < held_stress[1] && held_stress[1] < held_stress[0]);
    assert!((held_thickness[2] - held_thickness[0]).abs() > 1e-6);
    assert!(accepted.dissipated > 0.0);
    close(
        accepted.work + accepted.reference_work,
        accepted.stored + accepted.dissipated,
        1e6,
        2e-12,
    );
}

#[test]
fn g0_exact_cycles_retain_pressure_stress_and_remainder_energy_without_bulk_penalty_heat() {
    let law = law(8e5);
    let extreme = law.clone().equilibrium;
    let extreme = ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6e4,
                c01: 4e4,
                kappa: 1e300,
            },
            extreme.stretch_limit,
        )
        .unwrap(),
        law.branches.clone(),
    )
    .unwrap()
    .with_reference_deformation(reference());
    let mut accepted = law.initial_remainder_state();
    let mut extreme_state = extreme.initial_remainder_state();
    let mut nonzero_pressure = false;
    for (a, dt) in [
        ([1.08, 0.04, 0.0, 0.95], 0.0),
        ([1.2, 0.12, -0.03, 0.97], 0.03),
        ([1.2, 0.12, -0.03, 0.97], 0.08),
        ([1.08, 0.04, 0.0, 0.95], 0.04),
    ] {
        let before = accepted.clone();
        let response = law
            .step_incompressible_plane_stress_response(&accepted, &a, dt, None, settings())
            .unwrap();
        let other = extreme
            .step_incompressible_plane_stress_response(&extreme_state, &a, dt, None, settings())
            .unwrap();
        assert_eq!(accepted, before);
        assert_eq!(response.material, other.material);
        assert_eq!(
            response.plane_stress.full_piola,
            other.plane_stress.full_piola
        );
        let f = response.plane_stress.deformation_gradient;
        let cof = cofactor(&f);
        close(determinant(&f), reference().volume_ratio(), 1.0, 2e-14);
        let original = law.step_remainder_response(&accepted, &f, dt).unwrap();
        assert_eq!(response.material, original);
        for i in 0..9 {
            close(
                response.plane_stress.full_piola[i],
                original.first_piola_pa[i] - response.plane_stress.pressure_pa * cof[i],
                1e6,
                2e-12,
            );
        }
        for i in [2, 5, 8] {
            close(response.plane_stress.full_piola[i], 0.0, 1e6, 5e-10);
        }
        close(
            original.loading_work_pa,
            original.state.stored - accepted.stored + original.dissipated_pa,
            1e6,
            2e-12,
        );
        close(
            response.plane_stress.energy_density_pa,
            original.incremental_potential_pa,
            1e6,
            2e-12,
        );
        assert!(original.dissipated_pa >= 0.0);
        nonzero_pressure |= response.plane_stress.pressure_pa.abs() > 1.0;
        accepted = response.material.state;
        extreme_state = other.material.state;
    }
    assert!(nonzero_pressure);
    assert!(accepted.dissipated > 0.0);
    close(
        accepted.work + accepted.reference_work,
        accepted.stored + accepted.dissipated,
        1e6,
        2e-12,
    );
}

#[test]
fn g3_rotated_trajectories_retain_the_same_material_memory_and_rotate_physical_stress() {
    let law = law(8e5);
    let mut full = law.initial_state();
    let mut rotated_full = law.initial_state();
    let mut rem = law.initial_remainder_state();
    let mut rotated_rem = law.initial_remainder_state();
    let r = [0.8, -0.6, 0.0, 0.6, 0.8, 0.0, 0.0, 0.0, 1.0];
    for (a, dt) in [
        ([1.08, 0.04, 0.0, 0.95], 0.0),
        ([1.19, 0.13, -0.02, 0.96], 0.04),
        ([1.19, 0.13, -0.02, 0.96], 0.07),
    ] {
        let rotated = [
            0.8 * a[0] - 0.6 * a[2],
            0.8 * a[1] - 0.6 * a[3],
            0.6 * a[0] + 0.8 * a[2],
            0.6 * a[1] + 0.8 * a[3],
        ];
        let a_full = law
            .step_plane_stress_response(&full, &a, dt, None, settings())
            .unwrap();
        let b_full = law
            .step_plane_stress_response(&rotated_full, &rotated, dt, None, settings())
            .unwrap();
        let a_rem = law
            .step_incompressible_plane_stress_response(&rem, &a, dt, None, settings())
            .unwrap();
        let b_rem = law
            .step_incompressible_plane_stress_response(&rotated_rem, &rotated, dt, None, settings())
            .unwrap();
        for (f, g, p, q, ea, eb, da, db) in [
            (
                a_full.plane_stress.deformation_gradient,
                b_full.plane_stress.deformation_gradient,
                a_full.material.first_piola_pa,
                b_full.material.first_piola_pa,
                a_full.material.state.stored,
                b_full.material.state.stored,
                a_full.material.dissipated_pa,
                b_full.material.dissipated_pa,
            ),
            (
                a_rem.plane_stress.deformation_gradient,
                b_rem.plane_stress.deformation_gradient,
                a_rem.plane_stress.full_piola,
                b_rem.plane_stress.full_piola,
                a_rem.material.state.stored,
                b_rem.material.state.stored,
                a_rem.material.dissipated_pa,
                b_rem.material.dissipated_pa,
            ),
        ] {
            let rf = multiply(&r, &f);
            let rp = multiply(&r, &p);
            for i in 0..9 {
                close(g[i], rf[i], 1.0, 3e-9);
                close(q[i], rp[i], 1e6, 3e-9);
            }
            close(ea, eb, 1e6, 2e-10);
            close(da, db, 1e6, 2e-10);
        }
        for (a, b) in a_full
            .material
            .state
            .inverse_viscous_cauchy_green
            .iter()
            .flatten()
            .zip(
                b_full
                    .material
                    .state
                    .inverse_viscous_cauchy_green
                    .iter()
                    .flatten(),
            )
        {
            close(*a, *b, 1.0, 3e-9);
        }
        for (a, b) in a_rem
            .material
            .state
            .inverse_viscous_cauchy_green
            .iter()
            .flatten()
            .zip(
                b_rem
                    .material
                    .state
                    .inverse_viscous_cauchy_green
                    .iter()
                    .flatten(),
            )
        {
            close(*a, *b, 1.0, 3e-9);
        }
        close(
            a_rem.plane_stress.pressure_pa,
            b_rem.plane_stress.pressure_pa,
            1e6,
            3e-9,
        );
        full = a_full.material.state;
        rotated_full = b_full.material.state;
        rem = a_rem.material.state;
        rotated_rem = b_rem.material.state;
    }
}

#[test]
fn g0_converged_transverse_warm_starts_do_not_add_newton_or_material_time_steps() {
    let law = law(8e5);
    let full = law.initial_state();
    let rem = law.initial_remainder_state();
    let a = [1.2, 0.09, -0.02, 0.96];
    let dt = 0.05;
    let first = law
        .step_plane_stress_response(&full, &a, dt, None, settings())
        .unwrap();
    let f = first.plane_stress.deformation_gradient;
    let repeated = law
        .step_plane_stress_response(&full, &a, dt, Some([f[2], f[5], f[8]]), settings())
        .unwrap();
    assert_eq!(repeated.material, first.material);
    assert_eq!(repeated.plane_stress.iterations, 0);
    let exact = law
        .step_incompressible_plane_stress_response(&rem, &a, dt, None, settings())
        .unwrap();
    let f = exact.plane_stress.deformation_gradient;
    let repeated = law
        .step_incompressible_plane_stress_response(&rem, &a, dt, Some([f[2], f[5]]), settings())
        .unwrap();
    assert_eq!(repeated.material, exact.material);
    assert_eq!(repeated.plane_stress.iterations, 0);
    let frozen = law.freeze_remainder(&rem, dt).unwrap();
    let baseline = incompressible_plane_stress(&frozen, &a, settings()).unwrap();
    let seeded =
        incompressible_plane_stress_with_seed(&frozen, &a, [f[2], f[5]], settings()).unwrap();
    assert_eq!(seeded.iterations, 0);
    for i in 0..4 {
        for j in 0..4 {
            close(seeded.tangent[i][j], baseline.tangent[i][j], 1e6, 2e-10);
        }
    }
    let fallback =
        incompressible_plane_stress_with_seed(&frozen, &a, [1e5, -1e5], settings()).unwrap();
    for i in 0..9 {
        close(fallback.full_piola[i], baseline.full_piola[i], 1e6, 2e-10);
    }
}

#[test]
fn g0_zero_duration_jumps_change_strain_and_stress_without_relaxation_loss_or_memory_reset() {
    let law = law(8e5);
    let loaded = [1.21, 0.1, -0.02, 0.98];
    let next = [1.1, 0.02, 0.03, 0.96];
    let full = law
        .step_plane_stress_response(&law.initial_state(), &loaded, 0.06, None, settings())
        .unwrap()
        .material
        .state;
    let rem = law
        .step_incompressible_plane_stress_response(
            &law.initial_remainder_state(),
            &loaded,
            0.06,
            None,
            settings(),
        )
        .unwrap()
        .material
        .state;
    let response = law
        .step_plane_stress_response(&full, &next, 0.0, None, settings())
        .unwrap();
    let exact = law
        .step_incompressible_plane_stress_response(&rem, &next, 0.0, None, settings())
        .unwrap();
    assert_eq!(response.material.dissipated_pa, 0.0);
    assert_eq!(exact.material.dissipated_pa, 0.0);
    assert_eq!(response.material.state.dissipated, full.dissipated);
    assert_eq!(exact.material.state.dissipated, rem.dissipated);
    assert_eq!(
        response.material.state.inverse_viscous_cauchy_green,
        full.inverse_viscous_cauchy_green
    );
    assert_eq!(
        exact.material.state.inverse_viscous_cauchy_green,
        rem.inverse_viscous_cauchy_green
    );
    assert_ne!(
        response.material.state.deformation_gradient,
        full.deformation_gradient
    );
    assert_ne!(
        exact.material.state.deformation_gradient,
        rem.deformation_gradient
    );
}

#[test]
fn g0_refused_intervals_and_final_ledger_overflow_leave_both_accepted_histories_unchanged() {
    let law = law(8e5);
    let full = law.initial_state();
    let rem = law.initial_remainder_state();
    let before = full.clone();
    let before_rem = rem.clone();
    let a = [1.2, 0.09, -0.02, 0.96];
    for dt in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            law.step_plane_stress_response(&full, &a, dt, None, settings())
                .is_err()
        );
        assert!(
            law.step_incompressible_plane_stress_response(&rem, &a, dt, None, settings())
                .is_err()
        );
    }
    for invalid in [
        [1.0, 0.0, 0.0, -1.0],
        [10.0, 0.0, 0.0, 1.0],
        [f64::NAN, 0.0, 0.0, 1.0],
    ] {
        assert!(
            law.step_plane_stress_response(&full, &invalid, 0.03, None, settings())
                .is_err()
        );
        assert!(
            law.step_incompressible_plane_stress_response(&rem, &invalid, 0.03, None, settings())
                .is_err()
        );
    }
    let exhausted = PlaneStressSettings {
        max_iterations: 0,
        ..settings()
    };
    assert!(
        law.step_plane_stress_response(&full, &a, 0.03, None, exhausted)
            .is_err()
    );
    assert!(
        law.step_incompressible_plane_stress_response(&rem, &a, 0.03, None, exhausted)
            .is_err()
    );
    assert!(
        law.step_plane_stress_response(&full, &a, 0.03, Some([0.0, 0.0, -1.0]), settings())
            .is_err()
    );
    assert!(
        law.step_incompressible_plane_stress_response(
            &rem,
            &a,
            0.03,
            Some([f64::NAN, 0.0]),
            settings()
        )
        .is_err()
    );
    let rebound = law
        .clone()
        .with_reference_deformation(ReferenceDeformation::identity());
    assert!(
        rebound
            .step_plane_stress_response(&full, &a, 0.03, None, settings())
            .is_err()
    );
    assert!(
        rebound
            .step_incompressible_plane_stress_response(&rem, &a, 0.03, None, settings())
            .is_err()
    );
    assert_eq!(full, before);
    assert_eq!(rem, before_rem);

    // Public finite ledger values can still overflow only after the converged
    // reduction. The native final accounting gate must remain transactional.
    let mut overflow = full.clone();
    overflow.work = f64::MAX;
    overflow.stored = -f64::MAX;
    let mut overflow_rem = rem.clone();
    overflow_rem.work = f64::MAX;
    overflow_rem.stored = -f64::MAX;
    let snapshot = overflow.clone();
    let snapshot_rem = overflow_rem.clone();
    assert!(law.freeze(&overflow, 0.03).is_ok());
    assert!(law.freeze_remainder(&overflow_rem, 0.03).is_ok());
    assert!(
        law.step_plane_stress_response(&overflow, &a, 0.03, None, settings())
            .is_err()
    );
    assert!(
        law.step_incompressible_plane_stress_response(&overflow_rem, &a, 0.03, None, settings())
            .is_err()
    );
    assert_eq!(overflow, snapshot);
    assert_eq!(overflow_rem, snapshot_rem);
}

#[test]
fn g0_default_seed_uses_reference_volume_even_when_virgin_total_thickness_is_negative() {
    let g = ReferenceDeformation::try_new([-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0]).unwrap();
    let law = law(8e5).with_reference_deformation(g);
    let response = law
        .step_plane_stress_response(
            &law.initial_state(),
            &[1.0, 0.0, 0.0, 1.0],
            0.0,
            None,
            settings(),
        )
        .unwrap();
    let exact = law
        .step_incompressible_plane_stress_response(
            &law.initial_remainder_state(),
            &[1.0, 0.0, 0.0, 1.0],
            0.0,
            None,
            settings(),
        )
        .unwrap();
    assert_eq!(
        response.plane_stress.deformation_gradient,
        [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    );
    assert_eq!(
        exact.plane_stress.deformation_gradient,
        response.plane_stress.deformation_gradient
    );
    close(response.material.state.stored, 0.0, 1e6, 1e-14);
    close(exact.material.state.stored, 0.0, 1e6, 1e-14);
    assert_eq!(response.material.dissipated_pa, 0.0);
    assert_eq!(exact.material.dissipated_pa, 0.0);
}
