//! Actual momentum, constitutive-memory and time-accuracy checks, using the
//! original finite-strain solver and an independent small-strain Maxwell oracle.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
use fs_solid::hyper3::HyperTetSettings;
use fs_solid::hyper3::visco::dynamics::{FiniteViscoMotionControls, FiniteViscoTetMotion};
use fs_solid::hyper3::visco::{
    FiniteViscoStepControls, FiniteViscoTetError, FiniteViscoTetSolid, ViscoelasticLawField,
};
use fs_solid::linear3::TetAssemblyBudget;

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const VOLUME: f64 = 1.0 / 6.0;
const DENSITY: f64 = 60.0;
const MU: f64 = 1200.0;
const BRANCH_MU: f64 = 3600.0;
const TAU: f64 = 0.07;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 42,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn material(viscous: bool) -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: MU,
                lambda: 2.0 * MU,
            },
            4.0,
        )
        .unwrap(),
        if viscous {
            vec![FiniteStrainBranch {
                energy: BranchEnergy::NeoHookean { mu_pa: BRANCH_MU },
                tau_s: TAU,
            }]
        } else {
            vec![]
        },
    )
    .unwrap()
}
fn settings() -> HyperTetSettings {
    HyperTetSettings {
        load_steps: 1,
        relative_tolerance: 1e-11,
        linear_tolerance: 1e-12,
        ..HyperTetSettings::default()
    }
}
fn controls() -> FiniteViscoMotionControls {
    FiniteViscoMotionControls {
        material: FiniteViscoStepControls {
            equilibrium_work_tolerance_j: 2e-9,
            maximum_endpoint_work_defect_j: 1.0,
        },
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}
fn supports() -> Vec<(usize, f64)> {
    (0..12)
        .filter(|&dof| dof != 3)
        .map(|dof| (dof, 0.0))
        .collect()
}
fn initial<'a>(
    law: &'a ViscoHyperelastic,
    extension: f64,
    velocity: [[f64; 3]; 4],
    cx: &Cx<'_>,
) -> FiniteViscoTetMotion<'a> {
    let mut solid = FiniteViscoTetSolid::try_new(
        &NODES,
        &TETS,
        ViscoelasticLawField::Uniform(law),
        TetAssemblyBudget::default(),
        cx,
    )
    .unwrap();
    if extension != 0.0 {
        let fixed: Vec<_> = (0..12)
            .map(|dof| (dof, if dof == 3 { extension } else { 0.0 }))
            .collect();
        solid = solid
            .advance(
                0.0,
                &fixed,
                &[],
                [0.0; 3],
                settings(),
                controls().material,
                cx,
            )
            .unwrap()
            .solid;
    }
    FiniteViscoTetMotion::try_new(solid, &[DENSITY], &velocity, cx).unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} versus {expected:e}, allowance={tolerance:e}"
    );
}

#[test]
fn free_rigid_acceleration_uses_original_mass_and_no_artificial_supports() {
    with_cx(|cx| {
        let law = material(true);
        let v0 = [0.2, -0.1, 0.05];
        let acceleration = [0.3, -0.2, 0.1];
        let dt = 0.02;
        let mut motion = initial(&law, 0.0, [v0; 4], cx);
        for interval in 1..=3 {
            let count = f64::from(interval);
            let expected_v: [f64; 3] =
                core::array::from_fn(|c| v0[c] + count * dt * acceleration[c]);
            let expected_u: [f64; 3] = core::array::from_fn(|c| {
                count * dt * v0[c] + 0.5 * count * (count + 1.0) * dt * dt * acceleration[c]
            });
            let step = motion
                .advance(
                    dt,
                    &[],
                    &[],
                    acceleration.map(|a| DENSITY * a),
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            assert!(step.equilibrium.reactions_n.is_empty());
            for node in 0..4 {
                for c in 0..3 {
                    close(
                        step.motion.solid().displacement_m()[node][c],
                        expected_u[c],
                        2e-11,
                    );
                    close(step.motion.velocity_mps()[node][c], expected_v[c], 2e-9);
                }
            }
            close(step.motion.reference_mass_kg(), DENSITY * VOLUME, 1e-13);
            assert_eq!(step.motion.reference_density_kg_m3(), &[DENSITY]);
            close(step.material_energy.dissipated_j, 0.0, 1e-15);
            close(step.material_energy.stored_after_j, 0.0, 1e-12);
            close(
                step.inertia.numerical_loss_j,
                0.5 * DENSITY * VOLUME * dt * dt * acceleration.iter().map(|a| a * a).sum::<f64>(),
                1e-10,
            );
            close(step.balance_error_j, 0.0, 2e-9);
            motion = step.motion;
        }
    });
}

#[test]
fn deforming_viscoelastic_oscillator_matches_independent_maxwell_be_recurrence() {
    with_cx(|cx| {
        let law = material(true);
        let amplitude = 1e-4;
        let dt = 0.0025;
        let fixed = supports();
        let mut motion = initial(&law, amplitude, [[0.0; 3]; 4], cx);
        // Independent one-free-coordinate limit: original consistent diagonal
        // mass rho*V/10, equilibrium constrained modulus lambda+2mu, and each
        // ISOCHORIC Maxwell branch constrained modulus 4mu_b/3.
        let mass = DENSITY * VOLUME / 10.0;
        let stiffness = 4.0 * MU * VOLUME;
        let branch_stiffness = 4.0 * BRANCH_MU * VOLUME / 3.0;
        let (mut u, mut v, mut q) = (amplitude, 0.0, branch_stiffness * amplitude);
        let (mut material_loss, mut kinetic_numerical_loss) = (0.0, 0.0);
        let mut crossed_reference = false;
        for _ in 0..40 {
            let ratio = 1.0 + dt / TAU;
            v = (mass * v - dt * (stiffness * u + q / ratio))
                / (mass + dt * dt * (stiffness + branch_stiffness / ratio));
            u += dt * v;
            q = (q + branch_stiffness * dt * v) / ratio;
            let old_state = motion.solid().element_states()[0].clone();
            let step = motion
                .advance(dt, &fixed, &[], [0.0; 3], settings(), controls(), cx)
                .unwrap();
            // The nonlinear finite-strain correction is bounded explicitly at
            // this small amplitude, instead of modifying the production law.
            close(
                step.motion.solid().displacement_m()[1][0],
                u,
                5e-4 * amplitude,
            );
            close(step.motion.velocity_mps()[1][0], v, 2e-6);
            let actual = &step.motion.solid().element_states()[0];
            let native = law
                .step_response(&old_state, &actual.deformation_gradient, dt)
                .unwrap();
            assert_eq!(*actual, native.state);
            close(
                step.material_energy.dissipated_j,
                VOLUME * native.dissipated_pa,
                1e-15,
            );
            let energy = &step.material_energy;
            close(
                energy.dead_load_endpoint_work_j + energy.support_endpoint_work_j,
                energy.stored_after_j - energy.stored_before_j
                    + energy.dissipated_j
                    + energy.endpoint_work_defect_j
                    + step.inertia.kinetic_change_j
                    + step.inertia.numerical_loss_j
                    + step.balance_error_j,
                2e-9,
            );
            close(
                energy.inertial_endpoint_work_j,
                step.inertia.inertial_work_j,
                0.0,
            );
            assert!(energy.dissipated_j >= 0.0 && step.inertia.numerical_loss_j >= 0.0);
            material_loss += energy.dissipated_j;
            kinetic_numerical_loss += step.inertia.numerical_loss_j;
            crossed_reference |= step.motion.solid().displacement_m()[1][0] < 0.0;
            motion = step.motion;
        }
        assert!(
            crossed_reference,
            "actual deforming motion must cross its reference configuration"
        );
        assert!(material_loss > 1e-8 && kinetic_numerical_loss > 1e-8);
        close(
            motion.solid().element_states()[0].dissipated * VOLUME,
            material_loss,
            1e-14,
        );
    });
}

#[test]
fn elastic_motion_time_refinement_is_first_order_against_continuous_oscillator() {
    with_cx(|cx| {
        let law = material(false);
        let amplitude = 1e-4;
        let duration = 0.1;
        let omega = (40.0 * MU / DENSITY).sqrt();
        let expected_u = amplitude * (omega * duration).cos();
        let expected_v = -amplitude * omega * (omega * duration).sin();
        let fixed = supports();
        let mut errors = Vec::new();
        for intervals in [20, 40, 80] {
            let dt = duration / f64::from(intervals);
            let mut motion = initial(&law, amplitude, [[0.0; 3]; 4], cx);
            for _ in 0..intervals {
                let step = motion
                    .advance(dt, &fixed, &[], [0.0; 3], settings(), controls(), cx)
                    .unwrap();
                assert_eq!(step.material_energy.dissipated_j, 0.0);
                motion = step.motion;
            }
            let displacement_error = motion.solid().displacement_m()[1][0] - expected_u;
            let velocity_error = (motion.velocity_mps()[1][0] - expected_v) / omega;
            let error = displacement_error.hypot(velocity_error) / amplitude;
            println!(
                "visco-elastic-time-rung intervals={intervals} dt={dt:e} normalized_error={error:e}"
            );
            errors.push(error);
        }
        for pair in errors.windows(2) {
            let order = (pair[0] / pair[1]).log2();
            assert!(
                (0.85..1.15).contains(&order),
                "observed order={order}, errors={errors:?}"
            );
        }
        assert!(errors[2] < 0.06);
    });
}

#[test]
fn failed_numerical_loss_budget_preserves_entire_accepted_motion_and_retries() {
    with_cx(|cx| {
        let law = material(true);
        let motion = initial(&law, 1e-4, [[0.0; 3]; 4], cx);
        let old_u = motion.solid().displacement_m().to_vec();
        let old_v = motion.velocity_mps().to_vec();
        let old_states = motion.solid().element_states().to_vec();
        let fixed = supports();
        let strict = FiniteViscoMotionControls {
            maximum_numerical_kinetic_loss_j: 0.0,
            ..controls()
        };
        assert!(matches!(
            motion.advance(0.0025, &fixed, &[], [0.0; 3], settings(), strict, cx),
            Err(FiniteViscoTetError::EnergyGate {
                quantity: "backward-Euler numerical kinetic loss",
                ..
            })
        ));
        assert_eq!(motion.solid().displacement_m(), old_u);
        assert_eq!(motion.velocity_mps(), old_v);
        assert_eq!(motion.solid().element_states(), old_states);
        assert_eq!(motion.solid().time_s(), 0.0);
        let a = motion
            .advance(0.0025, &fixed, &[], [0.0; 3], settings(), controls(), cx)
            .unwrap();
        let b = motion
            .advance(0.0025, &fixed, &[], [0.0; 3], settings(), controls(), cx)
            .unwrap();
        assert_eq!(
            a.motion.solid().element_states(),
            b.motion.solid().element_states()
        );
        assert_eq!(
            a.motion.solid().displacement_m(),
            b.motion.solid().displacement_m()
        );
        assert_eq!(a.motion.velocity_mps(), b.motion.velocity_mps());
        assert_eq!(a.material_energy, b.material_energy);
        assert_eq!(a.inertia, b.inertia);
        assert!(
            motion
                .advance(0.0, &fixed, &[], [0.0; 3], settings(), controls(), cx)
                .is_err()
        );
    });
}
