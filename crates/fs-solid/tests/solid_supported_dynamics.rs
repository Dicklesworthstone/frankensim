//! Prescribed moving supports with the original consistent mass cross terms.
//! Continuous relative-coordinate solutions supply independent force and
//! actuator-work integrals; no production energy difference defines an oracle.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::spectral::CompleteSpectralSolidControls;
use fs_solid::dynamics::supported::{SupportedSolidControls, SupportedSpectralSolid};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::Coo;

#[path = "solid_supported_dynamics/frictional_impacts.rs"]
mod frictional_impacts;
#[path = "solid_supported_dynamics/normal_impacts.rs"]
mod normal_impacts;

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 62,
                kernel_id: 41,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn controls() -> SupportedSolidControls {
    SupportedSolidControls {
        spectral: CompleteSpectralSolidControls {
            maximum_free_dofs: 8,
            maximum_preparation_bytes: 2 * 1024 * 1024,
            maximum_preparation_multiply_terms: 2_000_000,
            maximum_step_multiply_terms: 2_000_000,
            maximum_relative_operator_error: 1e-12,
            maximum_absolute_eigen_residual: 0.0,
            maximum_relative_state_error: 1e-12,
            maximum_work_error_j: 1e-10,
            maximum_steps: 1000,
        },
        maximum_total_dofs: 8,
        maximum_preparation_bytes: 4 * 1024 * 1024,
        maximum_preparation_multiply_terms: 10_000_000,
        maximum_step_multiply_terms: 2_000_000,
        maximum_free_force_residual_n: 1e-10,
        maximum_work_error_j: 1e-10,
    }
}

fn pencil() -> TetElasticAssembly {
    let (mut mass, mut stiffness) = (Coo::new(2, 2), Coo::new(2, 2));
    for (row, entries) in [[2.0, 0.35], [0.35, 1.5]].iter().enumerate() {
        for (column, &entry) in entries.iter().enumerate() {
            mass.push(row, column, entry);
        }
    }
    for (row, entries) in [[18.0, -18.0], [-18.0, 18.0]].iter().enumerate() {
        for (column, &entry) in entries.iter().enumerate() {
            stiffness.push(row, column, entry);
        }
    }
    TetElasticAssembly {
        stiffness: stiffness.assemble(),
        mass: mass.assemble(),
        // These are original global labels, not positions in the two-row pencil.
        free_dofs: vec![3, 8],
        total_mass_kg: 4.2,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0],
    }
}

fn prepare(
    dt: f64,
    u: [f64; 2],
    v: [f64; 2],
    acceleration: f64,
    controls: SupportedSolidControls,
    cx: &Cx<'_>,
) -> SupportedSpectralSolid {
    SupportedSpectralSolid::prepare_with_basis(
        &pencil(),
        &[8],
        &[9.0],
        &[1.0 / 2.0_f64.sqrt()],
        dt,
        &u,
        &v,
        &[acceleration],
        controls,
        cx,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "actual={actual:e}, expected={expected:e}, allowance={allowance:e}"
    );
}

fn energy(u: [f64; 2], v: [f64; 2]) -> f64 {
    9.0 * (u[0] - u[1]).powi(2) + v[0] * v[0] + 0.35 * v[0] * v[1] + 0.75 * v[1] * v[1]
}

fn integral(start: f64, end: f64, f: impl Fn(f64) -> f64) -> f64 {
    // Independent quadrature of the closed-form continuous trajectory.
    let panels = 1024;
    let h = (end - start) / panels as f64;
    let mut sum = f(start) + f(end);
    for i in 1..panels {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(start + i as f64 * h);
    }
    h * sum / 3.0
}

#[derive(Clone, Copy)]
struct Reference {
    support_u: f64,
    support_v: f64,
    support_a: f64,
    jerk: f64,
    relative_u: f64,
    relative_v: f64,
    f0: f64,
    f_rate: f64,
    support_f0: f64,
    support_f_rate: f64,
}

impl Reference {
    fn fixture() -> Self {
        Self {
            support_u: 0.03,
            support_v: -0.2,
            support_a: 0.4,
            jerk: -0.6,
            relative_u: 0.12,
            relative_v: 0.08,
            f0: 0.7,
            f_rate: 0.25,
            support_f0: -0.5,
            support_f_rate: 0.11,
        }
    }

    fn state(self, t: f64) -> ([f64; 2], [f64; 2], [f64; 2]) {
        let y = self.support_u
            + self.support_v * t
            + 0.5 * self.support_a * t * t
            + self.jerk * t * t * t / 6.0;
        let vy = self.support_v + self.support_a * t + 0.5 * self.jerk * t * t;
        let ay = self.support_a + self.jerk * t;
        // r=x-y solves 2*r''+18*r=f_x-(2+0.35)*y''.
        let offset = (self.f0 - 2.35 * self.support_a) / 18.0;
        let slope = (self.f_rate - 2.35 * self.jerk) / 18.0;
        let a = self.relative_u - offset;
        let b = (self.relative_v - slope) / 3.0;
        let (sin, cos) = (3.0 * t).sin_cos();
        let r = offset + slope * t + a * cos + b * sin;
        let vr = slope - 3.0 * a * sin + 3.0 * b * cos;
        let ar = -9.0 * (a * cos + b * sin);
        ([y + r, y], [vy + vr, vy], [ay + ar, ay])
    }

    fn force(self, t: f64) -> [f64; 2] {
        [
            self.f0 + self.f_rate * t,
            self.support_f0 + self.support_f_rate * t,
        ]
    }

    fn reaction(self, t: f64) -> f64 {
        let (u, _, a) = self.state(t);
        0.35 * a[0] + 1.5 * a[1] + 18.0 * (u[1] - u[0]) - self.force(t)[1]
    }

    fn applied_power(self, t: f64) -> f64 {
        let f = self.force(t);
        let (_, v, _) = self.state(t);
        f[0] * v[0] + f[1] * v[1]
    }

    fn support_power(self, t: f64) -> f64 {
        self.reaction(t) * self.state(t).1[1]
    }
}

fn bits(s: &SupportedSpectralSolid) -> Vec<u64> {
    s.displacement_m()
        .iter()
        .chain(s.velocity_mps())
        .chain(s.support_acceleration_mps2())
        .copied()
        .chain([
            s.time_s(),
            s.stored_energy_j(),
            s.kinetic_energy_j(),
            s.cumulative_applied_work_j(),
            s.cumulative_support_work_j(),
        ])
        .map(f64::to_bits)
        .chain([s.accepted_steps()])
        .collect()
}

#[test]
fn g1_accelerating_supports_retain_cross_mass_reaction_and_independent_work() {
    with_cx(false, |cx| {
        let reference = Reference::fixture();
        let (u0, v0, a0) = reference.state(0.0);
        let dt = 0.037;
        let mut solid = prepare(dt, u0, v0, a0[1], controls(), cx);
        assert_eq!(solid.support_dofs(), &[8]);
        assert_eq!(solid.free_dofs(), &[3]);
        let (mut applied, mut actuator, mut maximum_balance) = (0.0, 0.0, 0.0_f64);
        for interval in 0..40 {
            let (start, end) = (interval as f64 * dt, (interval + 1) as f64 * dt);
            let (u, v, a) = reference.state(end);
            let report = solid
                .advance_affine_supports(
                    &[a[1]],
                    &reference.force(start),
                    &reference.force(end),
                    cx,
                )
                .unwrap();
            let applied_work = integral(start, end, |t| reference.applied_power(t));
            let support_work = integral(start, end, |t| reference.support_power(t));
            applied += applied_work;
            actuator += support_work;
            for i in 0..2 {
                near(solid.displacement_m()[i], u[i], 2e-12);
                near(solid.velocity_mps()[i], v[i], 3e-12);
            }
            near(report.support_reaction_n[0], reference.reaction(end), 2e-11);
            near(report.energy.applied_work_j, applied_work, 3e-12);
            near(report.energy.support_work_j, support_work, 3e-12);
            near(report.energy.cumulative_applied_work_j, applied, 5e-12);
            near(report.energy.cumulative_support_work_j, actuator, 5e-12);
            near(
                report.energy.stored_after_j + report.energy.kinetic_after_j,
                energy(u, v),
                5e-12,
            );
            near(energy(u, v) - energy(u0, v0), applied + actuator, 5e-12);
            maximum_balance = maximum_balance.max(report.energy.cumulative_balance_error_j.abs());
        }
        assert!(
            applied.abs() > 0.01 && actuator.abs() > 0.01,
            "both physical sources must perform live work"
        );
        println!(
            "moving support: applied_work_j={applied:e}, support_work_j={actuator:e}, maximum_balance_j={maximum_balance:e}"
        );
    });
}

#[test]
fn g1_consistent_cross_mass_preserves_uniform_accelerating_translation() {
    with_cx(false, |cx| {
        let (u0, v0, a0, jerk, dt) = (0.2, 0.3, -2.0, 4.0, 0.025);
        let mut solid = prepare(dt, [u0; 2], [v0; 2], a0, controls(), cx);
        let a = |t: f64| a0 + jerk * t;
        // The full inertial load uses both consistent off-diagonal entries.
        let force = |t: f64| [2.35 * a(t), 1.85 * a(t)];
        for interval in 0..40 {
            let (start, end) = (interval as f64 * dt, (interval + 1) as f64 * dt);
            let expected_u = u0 + v0 * end + a0 * end * end / 2.0 + jerk * end * end * end / 6.0;
            let expected_v = v0 + a0 * end + jerk * end * end / 2.0;
            let report = solid
                .advance_affine_supports(&[a(end)], &force(start), &force(end), cx)
                .unwrap();
            for i in 0..2 {
                near(solid.displacement_m()[i], expected_u, 2e-12);
                near(solid.velocity_mps()[i], expected_v, 2e-12);
            }
            near(report.support_reaction_n[0], 0.0, 2e-11);
            near(report.energy.stored_after_j, 0.0, 2e-14);
            near(report.energy.support_work_j, 0.0, 2e-12);
            near(
                report.energy.kinetic_after_j,
                0.5 * 4.2 * expected_v.powi(2),
                4e-12,
            );
        }
    });
}

#[test]
fn g3_partitioning_and_constrained_forces_preserve_the_physical_trajectory() {
    with_cx(false, |cx| {
        let reference = Reference::fixture();
        let changed = Reference {
            support_f0: 3.1,
            support_f_rate: -1.7,
            ..reference
        };
        let (u0, v0, a0) = reference.state(0.0);
        let dt = 0.24;
        let mut coarse = prepare(dt, u0, v0, a0[1], controls(), cx);
        let mut fine = prepare(dt / 8.0, u0, v0, a0[1], controls(), cx);
        let mut alternate = prepare(dt, u0, v0, a0[1], controls(), cx);
        let acceleration = reference.state(dt).2[1];
        let base = coarse
            .advance_affine_supports(
                &[acceleration],
                &reference.force(0.0),
                &reference.force(dt),
                cx,
            )
            .unwrap();
        let other = alternate
            .advance_affine_supports(&[acceleration], &changed.force(0.0), &changed.force(dt), cx)
            .unwrap();
        for i in 0..8 {
            let (start, end) = (i as f64 * dt / 8.0, (i + 1) as f64 * dt / 8.0);
            fine.advance_affine_supports(
                &[reference.state(end).2[1]],
                &reference.force(start),
                &reference.force(end),
                cx,
            )
            .unwrap();
        }
        for i in 0..2 {
            near(coarse.displacement_m()[i], fine.displacement_m()[i], 2e-13);
            near(coarse.velocity_mps()[i], fine.velocity_mps()[i], 2e-13);
            assert_eq!(
                coarse.displacement_m()[i].to_bits(),
                alternate.displacement_m()[i].to_bits()
            );
            assert_eq!(
                coarse.velocity_mps()[i].to_bits(),
                alternate.velocity_mps()[i].to_bits()
            );
        }
        near(
            coarse.cumulative_applied_work_j(),
            fine.cumulative_applied_work_j(),
            2e-13,
        );
        near(
            coarse.cumulative_support_work_j(),
            fine.cumulative_support_work_j(),
            2e-13,
        );
        let added_work = integral(0.0, dt, |t| {
            (changed.force(t)[1] - reference.force(t)[1]) * reference.state(t).1[1]
        });
        near(
            other.energy.applied_work_j - base.energy.applied_work_j,
            added_work,
            2e-13,
        );
        near(
            other.energy.support_work_j - base.energy.support_work_j,
            -added_work,
            2e-13,
        );
        near(
            other.support_reaction_n[0] - base.support_reaction_n[0],
            -(changed.force(dt)[1] - reference.force(dt)[1]),
            2e-12,
        );
    });
}

#[test]
fn g4_support_admission_and_failed_steps_publish_no_partial_motion_or_work() {
    with_cx(false, |cx| {
        let original = pencil();
        let basis = [1.0 / 2.0_f64.sqrt()];
        for supports in [vec![], vec![1], vec![8, 8], vec![8, 3], vec![3, 8]] {
            assert!(
                SupportedSpectralSolid::prepare_with_basis(
                    &original,
                    &supports,
                    &[9.0],
                    &basis,
                    0.02,
                    &[0.1, 0.0],
                    &[0.2, 0.0],
                    &vec![0.0; supports.len()],
                    controls(),
                    cx
                )
                .is_err()
            );
        }
        // The free subproblem is admissible in both cases. Its basis cannot
        // certify the original inertia or stiffness in prescribed directions.
        for invalid_mass in [true, false] {
            let mut invalid = original.clone();
            let mut matrix = Coo::new(2, 2);
            let entries = if invalid_mass {
                [[2.0, 2.0], [2.0, 1.5]]
            } else {
                // Zero diagonals do not make these nonzero cross terms PSD.
                [[0.0, 1.0], [1.0, 0.0]]
            };
            for (row, entries) in entries.iter().enumerate() {
                for (column, &entry) in entries.iter().enumerate() {
                    matrix.push(row, column, entry);
                }
            }
            if invalid_mass {
                invalid.mass = matrix.assemble();
            } else {
                invalid.stiffness = matrix.assemble();
            }
            assert!(
                SupportedSpectralSolid::prepare_with_basis(
                    &invalid,
                    &[8],
                    &[if invalid_mass { 9.0 } else { 0.0 }],
                    &basis,
                    0.02,
                    &[0.0; 2],
                    &[0.0; 2],
                    &[0.0],
                    controls(),
                    cx,
                )
                .is_err()
            );
        }
        for kind in 0..4 {
            let mut limited = controls();
            match kind {
                0 => limited.maximum_total_dofs = 1,
                1 => limited.maximum_preparation_bytes = 1,
                2 => limited.maximum_preparation_multiply_terms = 1,
                _ => limited.maximum_step_multiply_terms = 1,
            }
            assert!(
                SupportedSpectralSolid::prepare_with_basis(
                    &original,
                    &[8],
                    &[9.0],
                    &basis,
                    0.02,
                    &[0.1, 0.0],
                    &[0.2, 0.0],
                    &[0.0],
                    limited,
                    cx
                )
                .is_err()
            );
        }
        let make = || prepare(0.02, [0.1, 0.0], [0.2, 0.0], 0.0, controls(), cx);
        let (mut solid, mut clean) = (make(), make());
        let before = bits(&solid);
        for (acceleration, start, end) in [
            (vec![f64::NAN], vec![0.0, 0.0], vec![0.0, 0.0]),
            (vec![0.0], vec![0.0], vec![0.0, 0.0]),
            (vec![0.0], vec![0.0, 0.0], vec![f64::INFINITY, 0.0]),
            (vec![0.0], vec![f64::MAX, 0.0], vec![-f64::MAX, 0.0]),
        ] {
            assert!(
                solid
                    .advance_affine_supports(&acceleration, &start, &end, cx)
                    .is_err()
            );
            assert_eq!(bits(&solid), before);
        }
        with_cx(true, |cancelled| {
            assert!(
                solid
                    .advance_affine_supports(&[0.1], &[0.1, -0.1], &[0.2, -0.2], cancelled)
                    .is_err()
            );
        });
        assert_eq!(bits(&solid), before);
        solid
            .advance_affine_supports(&[0.1], &[0.1, -0.1], &[0.2, -0.2], cx)
            .unwrap();
        clean
            .advance_affine_supports(&[0.1], &[0.1, -0.1], &[0.2, -0.2], cx)
            .unwrap();
        assert_eq!(bits(&solid), bits(&clean));

        // An explicitly loose basis admission permits this inaccurate mode,
        // but the original full-system work balance must still reject motion.
        let mut strict = controls();
        strict.spectral.maximum_relative_operator_error = 0.02;
        strict.spectral.maximum_work_error_j = 1.0;
        strict.maximum_free_force_residual_n = 1.0;
        strict.maximum_work_error_j = 1e-16;
        let mut inaccurate = SupportedSpectralSolid::prepare_with_basis(
            &original,
            &[8],
            &[9.09],
            &basis,
            0.1,
            &[0.1, 0.0],
            &[0.2, 0.0],
            &[0.0],
            strict,
            cx,
        )
        .unwrap();
        let unchanged = bits(&inaccurate);
        assert!(
            inaccurate
                .advance_affine_supports(&[0.0], &[0.0; 2], &[0.0; 2], cx)
                .is_err()
        );
        assert_eq!(bits(&inaccurate), unchanged);
        let mut one_step = controls();
        one_step.spectral.maximum_steps = 1;
        let mut bounded = prepare(0.02, [0.1, 0.0], [0.2, 0.0], 0.0, one_step, cx);
        bounded
            .advance_affine_supports(&[0.0], &[0.0; 2], &[0.0; 2], cx)
            .unwrap();
        let final_state = bits(&bounded);
        assert!(
            bounded
                .advance_affine_supports(&[0.0], &[0.0; 2], &[0.0; 2], cx)
                .is_err()
        );
        assert_eq!(bits(&bounded), final_state);
    });
}

#[test]
fn g4_endpoint_energy_cannot_become_negative_within_an_operator_allowance() {
    with_cx(false, |cx| {
        let mut original = pencil();
        let mut stiffness = Coo::new(2, 2);
        stiffness.push(0, 0, 1e12);
        stiffness.push(1, 1, -0.1);
        original.stiffness = stiffness.assemble();
        let mut mass = Coo::new(2, 2);
        mass.push(0, 0, 1.0);
        mass.push(1, 1, 1.0);
        original.mass = mass.assemble();

        // The global operator allowance is 1 N/m, so this small negative
        // support stiffness lies inside that numerical admission tolerance.
        // Initially H_elastic=0; prescribed velocity would make it -0.05 J.
        // Its exactly matching negative reaction work must not legitimize
        // publishing a negative stored energy outside the energy allowance.
        let mut solid = SupportedSpectralSolid::prepare_with_basis(
            &original,
            &[8],
            &[1e12],
            &[1.0],
            1.0,
            &[0.0; 2],
            &[0.0, 1.0],
            &[0.0],
            controls(),
            cx,
        )
        .unwrap();
        let original_state = bits(&solid);
        assert!(
            solid
                .advance_affine_supports(&[0.0], &[0.0; 2], &[0.0; 2], cx)
                .is_err()
        );
        assert_eq!(bits(&solid), original_state);
    });
}

#[test]
fn g1_support_impulse_retains_cross_mass_reaction_work_and_continuation() {
    with_cx(false, |cx| {
        let dt = 0.02;
        let acceleration = 0.4;
        let force = [0.7, -0.3];
        let mut solid = prepare(
            dt,
            [0.12, -0.03],
            [0.22, -0.11],
            acceleration,
            controls(),
            cx,
        );
        solid
            .advance_affine_supports(&[acceleration], &force, &force, cx)
            .unwrap();
        let u: [f64; 2] = solid.displacement_m().try_into().unwrap();
        let before: [f64; 2] = solid.velocity_mps().try_into().unwrap();
        let applied_before = solid.cumulative_applied_work_j();
        let support_before = solid.cumulative_support_work_j();
        let j = [0.13, -0.29];
        let after_support = 0.45;
        let dc = after_support - before[1];
        let df = (j[0] - 0.35 * dc) / 2.0;
        let after = [before[0] + df, after_support];
        let reaction = 0.35 * df + 1.5 * dc - j[1];
        let mean = [0.5 * (before[0] + after[0]), 0.5 * (before[1] + after[1])];
        let applied_work = j[0] * mean[0] + j[1] * mean[1];
        let support_work = reaction * mean[1];
        let report = solid
            .apply_momentum_impulse(&j, &[after_support], 1e-12, cx)
            .unwrap();
        assert_eq!(solid.displacement_m(), &u);
        assert_eq!(solid.support_acceleration_mps2(), &[acceleration]);
        assert_eq!(solid.accepted_steps(), 1);
        assert_eq!(solid.time_s(), dt);
        assert_eq!(report.time_s, dt);
        for i in 0..2 {
            near(solid.velocity_mps()[i], after[i], 2e-14);
            near(report.impulse_residual_kg_m_s[i], 0.0, 1e-13);
        }
        near(report.momentum_change_kg_m_s[0], j[0], 2e-14);
        near(
            report.momentum_change_kg_m_s[1],
            0.35 * df + 1.5 * dc,
            2e-14,
        );
        near(report.support_reaction_kg_m_s[0], reaction, 2e-14);
        near(report.energy.applied_work_j, applied_work, 2e-14);
        near(report.energy.support_work_j, support_work, 2e-14);
        near(report.kinetic_change_j, applied_work + support_work, 2e-14);
        near(
            solid.cumulative_applied_work_j(),
            applied_before + applied_work,
            2e-14,
        );
        near(
            solid.cumulative_support_work_j(),
            support_before + support_work,
            2e-14,
        );
        near(report.energy.cumulative_balance_error_j, 0.0, 2e-12);

        // The native continuation must start at the event's new velocities,
        // retaining the original step clock, support acceleration and all modes.
        let reference = Reference {
            support_u: u[1],
            support_v: after[1],
            support_a: acceleration,
            jerk: 0.0,
            relative_u: u[0] - u[1],
            relative_v: after[0] - after[1],
            f0: force[0],
            f_rate: 0.0,
            support_f0: force[1],
            support_f_rate: 0.0,
        };
        let step = solid
            .advance_affine_supports(&[acceleration], &force, &force, cx)
            .unwrap();
        let (expected_u, expected_v, _) = reference.state(dt);
        for i in 0..2 {
            near(solid.displacement_m()[i], expected_u[i], 2e-12);
            near(solid.velocity_mps()[i], expected_v[i], 2e-12);
        }
        near(step.energy.cumulative_balance_error_j, 0.0, 2e-12);
        assert_eq!(solid.accepted_steps(), 2);
    });
}

#[test]
fn g3_support_impulse_preserves_uniform_boost_and_source_partition() {
    with_cx(false, |cx| {
        let mut boosted = prepare(0.04, [0.2; 2], [0.1; 2], 0.0, controls(), cx);
        let dv = 0.5;
        let report = boosted
            .apply_momentum_impulse(&[2.35 * dv, 1.85 * dv], &[0.6], 1e-12, cx)
            .unwrap();
        for &v in boosted.velocity_mps() {
            near(v, 0.6, 2e-14);
        }
        near(report.support_reaction_kg_m_s[0], 0.0, 2e-14);
        near(report.energy.support_work_j, 0.0, 2e-14);
        near(
            report.kinetic_change_j,
            0.5 * 4.2 * (0.6 * 0.6 - 0.1 * 0.1),
            2e-14,
        );
        boosted
            .advance_affine_supports(&[0.0], &[0.0; 2], &[0.0; 2], cx)
            .unwrap();
        for &u in boosted.displacement_m() {
            near(u, 0.224, 2e-12);
        }

        let mut a = prepare(0.02, [0.1, -0.02], [0.2, 0.4], 0.0, controls(), cx);
        let mut b = prepare(0.02, [0.1, -0.02], [0.2, 0.4], 0.0, controls(), cx);
        let first = a
            .apply_momentum_impulse(&[0.2, 0.0], &[0.4], 1e-12, cx)
            .unwrap();
        let second = b
            .apply_momentum_impulse(&[0.2, 5.0], &[0.4], 1e-12, cx)
            .unwrap();
        assert_eq!(a.velocity_mps(), b.velocity_mps());
        near(
            second.support_reaction_kg_m_s[0] - first.support_reaction_kg_m_s[0],
            -5.0,
            2e-14,
        );
        near(
            second.energy.applied_work_j - first.energy.applied_work_j,
            2.0,
            2e-14,
        );
        near(
            second.energy.support_work_j - first.energy.support_work_j,
            -2.0,
            2e-14,
        );
        near(second.kinetic_change_j, first.kinetic_change_j, 0.0);
        // A no-op event must not round the accepted basis or cumulative accounts.
        let before = bits(&b);
        b.apply_momentum_impulse(&[0.0; 2], &[0.4], 0.0, cx)
            .unwrap();
        assert_eq!(bits(&b), before);
        // An exact no-op also retains the original energy reference when one
        // store is below the ULP of the other. It cannot invent a 0.75 J debt.
        let mut separated = prepare(
            0.02,
            [2.0_f64.powi(50), 0.0],
            [0.0, 1.0],
            0.0,
            controls(),
            cx,
        );
        let before = bits(&separated);
        separated
            .apply_momentum_impulse(&[0.0; 2], &[1.0], 0.0, cx)
            .unwrap();
        assert_eq!(bits(&separated), before);
    });
}

#[test]
fn g4_support_impulse_refuses_before_either_owner_publishes() {
    with_cx(false, |cx| {
        let mut solid = prepare(0.02, [0.1, 0.0], [0.0; 2], 0.0, controls(), cx);
        let mut witness = prepare(0.02, [0.1, 0.0], [0.0; 2], 0.0, controls(), cx);
        let before = bits(&solid);
        for (j, support, allowance) in [
            (vec![0.0], vec![0.0], 1e-12),
            (vec![0.0; 2], vec![], 1e-12),
            (vec![f64::NAN, 0.0], vec![0.0], 1e-12),
            (vec![0.0; 2], vec![f64::INFINITY], 1e-12),
            (vec![0.0; 2], vec![0.0], -1.0),
        ] {
            assert!(
                solid
                    .apply_momentum_impulse(&j, &support, allowance, cx)
                    .is_err()
            );
            assert_eq!(bits(&solid), before);
        }
        // The free trial is admissible, but this support's large applied and
        // opposite reaction impulses cannot retain its 0.035 kg m/s cross-mass
        // increment. The outer original-row audit must refuse AFTER that trial.
        let error = solid
            .apply_momentum_impulse(&[0.2, 2.0_f64.powi(80)], &[0.0], 1e-12, cx)
            .unwrap_err();
        assert!(matches!(
            error,
            fs_solid::dynamics::supported::SupportedSolidError::Impulse(
                fs_solid::dynamics::LinearSolidImpulseError::Momentum { .. }
            )
        ));
        assert_eq!(bits(&solid), before);
        with_cx(true, |cancelled| {
            assert!(
                solid
                    .apply_momentum_impulse(&[0.2, 0.0], &[0.3], 1e-12, cancelled)
                    .is_err()
            );
        });
        assert_eq!(bits(&solid), before);
        solid
            .advance_affine_supports(&[0.0], &[0.0; 2], &[0.0; 2], cx)
            .unwrap();
        witness
            .advance_affine_supports(&[0.0], &[0.0; 2], &[0.0; 2], cx)
            .unwrap();
        assert_eq!(
            bits(&solid),
            bits(&witness),
            "hidden free candidate must also roll back"
        );
    });
}
