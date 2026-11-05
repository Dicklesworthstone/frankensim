//! Prescribed eigenstrain work with the original retained Maxwell histories.

use fs_material::MaterialError;
use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::hyper::{Hyperelastic, HyperelasticModel, OgdenTerm};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

const I: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const G0: [f64; 9] = [1.1, 0.08, 0.0, -0.03, 0.94, 0.02, 0.0, 0.0, 1.05];
const G1: [f64; 9] = [1.2, 0.16, -0.03, 0.02, 0.88, 0.09, 0.0, -0.04, 1.12];
const F: [f64; 9] = [1.24, 0.12, 0.02, 0.01, 0.98, 0.04, 0.0, -0.02, 1.03];
const KAPPA: f64 = 8e5;

fn law(g: [f64; 9]) -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 3e4,
                c01: 2e4,
                kappa: KAPPA,
            },
            4.0,
        )
        .unwrap(),
        vec![
            FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: 3e4 },
                tau_s: 0.2,
            },
            FiniteStrainBranch {
                energy: BranchEnergy::Ogden {
                    terms: vec![OgdenTerm {
                        mu_pa: 2e4,
                        alpha: 1.6,
                    }],
                },
                tau_s: 0.07,
            },
        ],
    )
    .unwrap()
    .with_reference_deformation(ReferenceDeformation::try_new(g).unwrap())
}

fn near(a: f64, b: f64, scale: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= 2e-11 * scale,
        "{a:e} vs {b:e}, scale {scale:e}"
    );
}

fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}

#[test]
fn held_deformation_reference_change_preserves_memory_and_separates_all_work() {
    let old = law(G0);
    let target = law(G1);
    let accepted = old
        .step_response(&old.initial_state(), &F, 0.06)
        .unwrap()
        .state;
    let rem = old
        .step_remainder_response(&old.initial_remainder_state(), &F, 0.06)
        .unwrap()
        .state;
    let snapshot = (accepted.clone(), rem.clone());
    assert!(accepted.dissipated > 0.0);
    assert!(target.freeze(&accepted, 0.0).is_err());
    let change = target.change_reference_from(&old, &accepted).unwrap();
    let rem_change = target.change_remainder_reference_from(&old, &rem).unwrap();
    assert_eq!(change.state.deformation_gradient, F);
    assert_eq!(
        change.state.inverse_viscous_cauchy_green,
        accepted.inverse_viscous_cauchy_green
    );
    assert_eq!(change.state.work.to_bits(), accepted.work.to_bits());
    assert_eq!(
        change.state.dissipated.to_bits(),
        accepted.dissipated.to_bits()
    );
    assert_eq!(rem_change.state.deformation_gradient, F);
    assert_eq!(
        rem_change.state.inverse_viscous_cauchy_green,
        rem.inverse_viscous_cauchy_green
    );
    assert_eq!(rem_change.state.work.to_bits(), rem.work.to_bits());
    assert_eq!(
        rem_change.state.dissipated.to_bits(),
        rem.dissipated.to_bits()
    );
    near(
        change.reference_work_pa,
        target.free_energy(&change.state, &F).unwrap() - accepted.stored,
        1e6,
    );
    near(
        rem_change.reference_work_pa,
        target.remainder_free_energy(&rem_change.state, &F).unwrap() - rem.stored,
        1e6,
    );
    assert_eq!(change.state.reference_work, change.reference_work_pa);
    assert_eq!(
        rem_change.state.reference_work,
        rem_change.reference_work_pa
    );
    assert_eq!(
        target
            .freeze(&change.state, 0.0)
            .unwrap()
            .first_piola(&F)
            .unwrap(),
        change.first_piola_pa
    );
    assert_eq!(
        target
            .freeze_remainder(&rem_change.state, 0.0)
            .unwrap()
            .remainder_stress_and_potential(&F)
            .unwrap()
            .0,
        rem_change.first_piola_pa
    );
    assert!(old.freeze(&change.state, 0.0).is_err());
    // Actual physical relaxation follows the instantaneous reference event.
    let next_f = [1.2, 0.08, 0.01, 0.0, 1.02, 0.04, 0.0, 0.01, 1.0];
    let next = target.step_response(&change.state, &next_f, 0.04).unwrap();
    let next_rem = target
        .step_remainder_response(&rem_change.state, &next_f, 0.04)
        .unwrap();
    assert_eq!(next.state.reference_work, change.state.reference_work);
    assert_eq!(
        next_rem.state.reference_work,
        rem_change.state.reference_work
    );
    assert!(next.dissipated_pa > 0.0 && next_rem.dissipated_pa > 0.0);
    for (work, reference, stored, dissipated) in [
        (
            change.state.work,
            change.state.reference_work,
            change.state.stored,
            change.state.dissipated,
        ),
        (
            rem_change.state.work,
            rem_change.state.reference_work,
            rem_change.state.stored,
            rem_change.state.dissipated,
        ),
        (
            next.state.work,
            next.state.reference_work,
            next.state.stored,
            next.state.dissipated,
        ),
        (
            next_rem.state.work,
            next_rem.state.reference_work,
            next_rem.state.stored,
            next_rem.state.dissipated,
        ),
    ] {
        near(work + reference - stored - dissipated, 0.0, 1e6);
    }
    assert_eq!((accepted, rem), snapshot);
}

#[test]
fn reference_cycles_are_reversible_and_kinetic_only_changes_are_exact_history_noops() {
    let old = law(G0);
    let target = law(G1);
    let accepted = old
        .step_response(&old.initial_state(), &F, 0.06)
        .unwrap()
        .state;
    let rem = old
        .step_remainder_response(&old.initial_remainder_state(), &F, 0.06)
        .unwrap()
        .state;
    let forward = target.change_reference_from(&old, &accepted).unwrap();
    let back = old.change_reference_from(&target, &forward.state).unwrap();
    let forward_rem = target.change_remainder_reference_from(&old, &rem).unwrap();
    let back_rem = old
        .change_remainder_reference_from(&target, &forward_rem.state)
        .unwrap();
    near(forward.reference_work_pa + back.reference_work_pa, 0.0, 1e6);
    near(
        forward_rem.reference_work_pa + back_rem.reference_work_pa,
        0.0,
        1e6,
    );
    near(back.state.stored, accepted.stored, 1e6);
    near(back_rem.state.stored, rem.stored, 1e6);
    assert_eq!(
        back.state.inverse_viscous_cauchy_green,
        accepted.inverse_viscous_cauchy_green
    );
    assert_eq!(back.state.work.to_bits(), accepted.work.to_bits());
    assert_eq!(
        back.state.dissipated.to_bits(),
        accepted.dissipated.to_bits()
    );
    let mut kinetic = old.clone().with_local_newton(2e-13, 96).unwrap();
    for branch in &mut kinetic.branches {
        branch.tau_s *= 3.0;
    }
    for source in [&old, &kinetic] {
        let change = kinetic.change_reference_from(source, &accepted).unwrap();
        let rem_change = kinetic
            .change_remainder_reference_from(source, &rem)
            .unwrap();
        assert_eq!(change.reference_work_pa.to_bits(), 0.0_f64.to_bits());
        assert_eq!(rem_change.reference_work_pa.to_bits(), 0.0_f64.to_bits());
        assert_eq!(change.state, accepted);
        assert_eq!(rem_change.state, rem);
        for (a, b) in [
            (change.state.work, accepted.work),
            (change.state.reference_work, accepted.reference_work),
            (change.state.stored, accepted.stored),
            (change.state.dissipated, accepted.dissipated),
        ] {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
}

#[test]
fn reference_events_refuse_incompatible_parameters_bindings_and_invalid_endpoints_atomically() {
    let old = law(G0);
    let target = law(G1);
    let accepted = old
        .step_response(&old.initial_state(), &F, 0.04)
        .unwrap()
        .state;
    let rem = old
        .step_remainder_response(&old.initial_remainder_state(), &F, 0.04)
        .unwrap()
        .state;
    let snapshot = (accepted.clone(), rem.clone());
    let mut altered_energy = target.clone();
    altered_energy.branches[0].energy = BranchEnergy::NeoHookean { mu_pa: 6e4 };
    let mut reordered = target.clone();
    reordered.branches.reverse();
    let mut altered_domain = target.clone();
    altered_domain.equilibrium.stretch_limit = 5.0;
    let mut altered_equilibrium = target.clone();
    altered_equilibrium.equilibrium.model = HyperelasticModel::MooneyRivlin {
        c10: 4e4,
        c01: 2e4,
        kappa: KAPPA,
    };
    let mut invalid_kinetics = target.clone();
    invalid_kinetics.branches[0].tau_s = f64::NAN;
    for invalid in [
        altered_energy,
        reordered,
        altered_domain,
        altered_equilibrium,
        invalid_kinetics,
    ] {
        assert!(matches!(
            invalid.change_reference_from(&old, &accepted),
            Err(MaterialError::Parameters { .. })
        ));
        assert!(matches!(
            invalid.change_remainder_reference_from(&old, &rem),
            Err(MaterialError::Parameters { .. })
        ));
    }
    assert!(target.change_reference_from(&law(I), &accepted).is_err());
    assert!(
        target
            .change_remainder_reference_from(&law(I), &rem)
            .is_err()
    );
    let inadmissible = target.clone().with_reference_deformation(
        ReferenceDeformation::try_new([0.01, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap(),
    );
    assert!(inadmissible.change_reference_from(&old, &accepted).is_err());
    assert!(
        inadmissible
            .change_remainder_reference_from(&old, &rem)
            .is_err()
    );
    let mut bad = accepted.clone();
    bad.reference_work = f64::NAN;
    assert!(target.change_reference_from(&old, &bad).is_err());
    assert!(old.step_response(&bad, &F, 0.01).is_err());
    let mut bad_rem = rem.clone();
    bad_rem.reference_work = f64::INFINITY;
    assert!(
        target
            .change_remainder_reference_from(&old, &bad_rem)
            .is_err()
    );
    assert!(old.step_remainder_response(&bad_rem, &F, 0.01).is_err());
    assert_eq!((accepted, rem), snapshot);
}

#[test]
fn prescribed_reference_work_is_objective_under_spatial_rotation() {
    let old = law(G0);
    let target = law(G1);
    let rotation = [0.8, -0.6, 0.0, 0.6, 0.8, 0.0, 0.0, 0.0, 1.0];
    let rf = multiply(&rotation, &F);
    let accepted = old
        .step_response(&old.initial_state(), &F, 0.06)
        .unwrap()
        .state;
    let rotated = old
        .step_response(&old.initial_state(), &rf, 0.06)
        .unwrap()
        .state;
    let rem = old
        .step_remainder_response(&old.initial_remainder_state(), &F, 0.06)
        .unwrap()
        .state;
    let rotated_rem = old
        .step_remainder_response(&old.initial_remainder_state(), &rf, 0.06)
        .unwrap()
        .state;
    let change = target.change_reference_from(&old, &accepted).unwrap();
    let rotated_change = target.change_reference_from(&old, &rotated).unwrap();
    let rem_change = target.change_remainder_reference_from(&old, &rem).unwrap();
    let rotated_rem_change = target
        .change_remainder_reference_from(&old, &rotated_rem)
        .unwrap();
    near(
        change.reference_work_pa,
        rotated_change.reference_work_pa,
        1e6,
    );
    near(
        rem_change.reference_work_pa,
        rotated_rem_change.reference_work_pa,
        1e6,
    );
    near(change.state.stored, rotated_change.state.stored, 1e6);
    for (actual, expected) in rotated_change
        .first_piola_pa
        .iter()
        .zip(multiply(&rotation, &change.first_piola_pa))
    {
        near(*actual, expected, 1e6);
    }
    for (actual, expected) in rotated_rem_change
        .first_piola_pa
        .iter()
        .zip(multiply(&rotation, &rem_change.first_piola_pa))
    {
        near(*actual, expected, 1e6);
    }
}

#[test]
fn held_isotropic_expansion_has_the_reference_volume_work_and_no_viscous_heat() {
    let old = law(I);
    let stretch: f64 = 1.25;
    let jg = stretch.powi(3);
    let target = law([stretch, 0.0, 0.0, 0.0, stretch, 0.0, 0.0, 0.0, stretch]);
    let accepted = old.initial_state();
    let change = target.change_reference_from(&old, &accepted).unwrap();
    let expected_work = jg * 0.5 * KAPPA * (1.0 / jg - 1.0).powi(2);
    near(change.reference_work_pa, expected_work, 1e6);
    near(change.state.stored, expected_work, 1e6);
    assert_eq!(change.state.work, 0.0);
    assert_eq!(change.state.dissipated, 0.0);
    assert_eq!(
        change.state.inverse_viscous_cauchy_green,
        accepted.inverse_viscous_cauchy_green
    );
    for i in 0..9 {
        near(
            change.first_piola_pa[i],
            if i % 4 == 0 {
                KAPPA * (1.0 / jg - 1.0)
            } else {
                0.0
            },
            1e6,
        );
    }
    let back = old.change_reference_from(&target, &change.state).unwrap();
    near(back.reference_work_pa, -expected_work, 1e6);
    near(back.state.reference_work, 0.0, 1e6);
    assert_eq!(back.state.dissipated, 0.0);
    // A remainder proposal deliberately does not assert mixed admissibility.
    let rem = target
        .change_remainder_reference_from(&old, &old.initial_remainder_state())
        .unwrap();
    let mixed = target
        .freeze_remainder(&rem.state, 0.0)
        .unwrap()
        .mixed_state(&I, 0.0, VolumetricTreatment::ExactIncompressible)
        .unwrap();
    near(mixed.pressure_residual, 1.0 - jg, 1.0);
    near(rem.reference_work_pa, 0.0, 1e6);
    assert_eq!(rem.state.dissipated, 0.0);
}

mod configurational_work {
    //! Independent configurational-work derivative checks for the original-volume potential.
    use fs_material::hyper::reference::ReferenceDeformation;
    use fs_material::hyper::{Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};
    use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

    fn determinant(g: &[f64; 9]) -> f64 {
        g[0] * (g[4] * g[8] - g[5] * g[7]) - g[1] * (g[3] * g[8] - g[5] * g[6])
            + g[2] * (g[3] * g[7] - g[4] * g[6])
    }
    fn inverse(g: &[f64; 9]) -> [f64; 9] {
        let j = determinant(g);
        [
            g[4] * g[8] - g[5] * g[7],
            g[2] * g[7] - g[1] * g[8],
            g[1] * g[5] - g[2] * g[4],
            g[5] * g[6] - g[3] * g[8],
            g[0] * g[8] - g[2] * g[6],
            g[2] * g[3] - g[0] * g[5],
            g[3] * g[7] - g[4] * g[6],
            g[1] * g[6] - g[0] * g[7],
            g[0] * g[4] - g[1] * g[3],
        ]
        .map(|x| x / j)
    }
    fn law() -> ViscoHyperelastic {
        ViscoHyperelastic::new(
            Hyperelastic::new(
                HyperelasticModel::Ogden {
                    terms: vec![OgdenTerm {
                        mu_pa: 1.1e5,
                        alpha: 1.7,
                    }],
                    volumetric: Volumetric::LogSquared { kappa: 8e5 },
                },
                3.0,
            )
            .unwrap(),
            vec![
                FiniteStrainBranch {
                    energy: BranchEnergy::NeoHookean { mu_pa: 3e4 },
                    tau_s: 0.2,
                },
                FiniteStrainBranch {
                    energy: BranchEnergy::Ogden {
                        terms: vec![OgdenTerm {
                            mu_pa: 2e4,
                            alpha: 1.6,
                        }],
                    },
                    tau_s: 0.07,
                },
            ],
        )
        .unwrap()
    }
    // d_G Psi0 = G^-T (Psi0 I - F^T P0), at fixed total F and fixed natural metrics.
    // This derivative follows from Psi0=det(G)*Psi(F G^-1), independently of the
    // implementation of the finite reference-work increment.
    fn conjugate(g: &[f64; 9], f: &[f64; 9], p: &[f64; 9], w: f64) -> [f64; 9] {
        let gi = inverse(g);
        core::array::from_fn(|ij| {
            let (i, j) = (ij / 3, ij % 3);
            (0..3)
                .map(|a| {
                    gi[3 * a + i]
                        * ((if a == j { w } else { 0.0 })
                            - (0..3).map(|k| f[3 * k + a] * p[3 * k + j]).sum::<f64>())
                })
                .sum()
        })
    }

    #[test]
    fn evolving_reference_work_gradient_matches_configurational_stress_all_nine_components() {
        let g = [1.15, 0.12, -0.03, 0.02, 0.94, 0.06, 0.01, -0.025, 1.07];
        let f = [1.26, 0.19, 0.04, -0.01, 0.89, 0.07, 0.025, -0.02, 1.12];
        let owner = law().with_reference_deformation(ReferenceDeformation::try_new(g).unwrap());
        let accepted = owner
            .step_response(&owner.initial_state(), &f, 0.11)
            .unwrap()
            .state;
        let remainder = owner
            .step_remainder_response(&owner.initial_remainder_state(), &f, 0.11)
            .unwrap()
            .state;
        assert_ne!(
            accepted.inverse_viscous_cauchy_green[0],
            [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
        );
        let stress = owner
            .step_response(&accepted, &f, 0.0)
            .unwrap()
            .first_piola_pa;
        let remainder_stress = owner
            .step_remainder_response(&remainder, &f, 0.0)
            .unwrap()
            .first_piola_pa;
        let expected = conjugate(&g, &f, &stress, accepted.stored);
        let expected_remainder = conjugate(&g, &f, &remainder_stress, remainder.stored);
        let h = 2e-6;
        for ij in 0..9 {
            let mut gp = g;
            let mut gm = g;
            gp[ij] += h;
            gm[ij] -= h;
            let plus = owner
                .clone()
                .with_reference_deformation(ReferenceDeformation::try_new(gp).unwrap());
            let minus = owner
                .clone()
                .with_reference_deformation(ReferenceDeformation::try_new(gm).unwrap());
            let wp = plus.change_reference_from(&owner, &accepted).unwrap();
            let wm = minus.change_reference_from(&owner, &accepted).unwrap();
            let rp = plus
                .change_remainder_reference_from(&owner, &remainder)
                .unwrap();
            let rm = minus
                .change_remainder_reference_from(&owner, &remainder)
                .unwrap();
            let finite_difference = (wp.reference_work_pa - wm.reference_work_pa) / (2.0 * h);
            let remainder_difference = (rp.reference_work_pa - rm.reference_work_pa) / (2.0 * h);
            for (actual, expect, label) in [
                (finite_difference, expected[ij], "full"),
                (remainder_difference, expected_remainder[ij], "remainder"),
            ] {
                assert!(
                    (actual - expect).abs() <= 2e-7 * 1e5_f64.max(expect.abs()),
                    "{label} G[{ij}] derivative {actual:e} versus {expect:e}"
                );
            }
        }
    }
}
