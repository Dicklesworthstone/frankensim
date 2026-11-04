//! Finite contact forces on prescribed native supports. Independent reaction
//! integrals and complete original momentum check the physical work partition.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::sustained::NaturalMapScales;
use fs_contact::sustained::solve::SolveError;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::spectral::CompleteSpectralSolidControls;
use fs_solid::dynamics::supported::{
    SupportedSolidControls, SupportedSolidError, SupportedSolidStep, SupportedSpectralSolid,
};
use fs_solid::dynamics::sustained::{
    SolidSustainedContact, SolidSustainedControls, SolidSustainedError, SolidSustainedPath,
    SolidSustainedStep,
};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::Coo;

const N: [f64; 3] = [0.0, 1.0, 0.0];
const X: [f64; 3] = [1.0, 0.0, 0.0];
const Z: [f64; 3] = [0.0, 0.0, 1.0];
const I: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
const BASIS: [f64; 4] = [1.0, 0.0, 0.0, 1.0];

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 62,
                    kernel_id: 92,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}
fn native_controls() -> SupportedSolidControls {
    SupportedSolidControls {
        spectral: CompleteSpectralSolidControls {
            maximum_free_dofs: 24,
            maximum_preparation_bytes: 64 * 1024 * 1024,
            maximum_preparation_multiply_terms: 1_000_000_000,
            maximum_step_multiply_terms: 10_000_000,
            maximum_relative_operator_error: 1e-12,
            maximum_absolute_eigen_residual: 0.0,
            maximum_relative_state_error: 1e-12,
            maximum_work_error_j: 1e-10,
            maximum_steps: 1000,
        },
        maximum_total_dofs: 24,
        maximum_preparation_bytes: 64 * 1024 * 1024,
        maximum_preparation_multiply_terms: 1_000_000_000,
        maximum_step_multiply_terms: 10_000_000,
        maximum_free_force_residual_n: 1e-10,
        maximum_work_error_j: 1e-9,
    }
}
fn controls() -> SolidSustainedControls {
    SolidSustainedControls {
        solve: Default::default(),
        maximum_workspace_bytes: 1024 * 1024,
        maximum_response_multiply_terms: 10_000_000,
        displacement_tolerance_m: 1e-10,
        velocity_tolerance_mps: 1e-9,
        force_tolerance_n: 1e-8,
        complementarity_tolerance_j: 1e-10,
        maximum_work_error_j: 1e-9,
    }
}
fn assembly(k: &[f64], m: &[f64]) -> TetElasticAssembly {
    let mut stiffness = Coo::new(3, 3);
    let mut mass = Coo::new(3, 3);
    for i in 0..3 {
        for j in 0..3 {
            stiffness.push(i, j, k[i * 3 + j]);
            mass.push(i, j, m[i * 3 + j]);
        }
    }
    TetElasticAssembly {
        stiffness: stiffness.assemble(),
        mass: mass.assemble(),
        free_dofs: vec![0, 1, 2],
        total_mass_kg: 1.0,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0],
    }
}
fn point(h: f64, boost: [f64; 2], cx: &Cx<'_>) -> SupportedSpectralSolid {
    SupportedSpectralSolid::prepare_with_basis(
        &assembly(&[0.0; 9], &I),
        &[2],
        &[0.0; 2],
        &BASIS,
        h,
        &[0.0; 3],
        &[3.0 + boost[0], 0.0, 4.0 + boost[1]],
        &[0.0],
        native_controls(),
        cx,
    )
    .unwrap()
}
fn contact(h: f64, mu: f64, surface: [f64; 3]) -> SolidSustainedContact<'static> {
    SolidSustainedContact {
        id: 31,
        velocity_weights: [&N, &X, &Z],
        gap_m: 0.0,
        prescribed_surface_velocity_mps: surface,
        friction_coefficient: mu,
        scales: NaturalMapScales {
            normal_n_per_m: 2.0 / (h * h),
            tangent_n_per_m: 2.0 / (h * h),
        },
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}
fn bits(s: &SupportedSpectralSolid) -> Vec<u64> {
    s.displacement_m()
        .iter()
        .chain(s.velocity_mps())
        .chain(s.support_acceleration_mps2())
        .copied()
        .chain([
            s.time_s(),
            s.cumulative_applied_work_j(),
            s.cumulative_support_work_j(),
            s.stored_energy_j(),
            s.kinetic_energy_j(),
        ])
        .map(f64::to_bits)
        .collect()
}
fn energy<P>(r: &SolidSustainedStep<P, SupportedSolidStep>) {
    near(r.energy_balance_error_j, 0.0, 1e-9);
    near(r.contact_work_balance_error_j, 0.0, 1e-9);
    near(r.native.energy.balance_error_j, 0.0, 1e-9);
    near(r.native.energy.cumulative_balance_error_j, 0.0, 1e-9);
    assert_eq!(
        r.support_work_j.to_bits(),
        r.native.energy.support_work_j.to_bits()
    );
    near(
        r.native.energy.applied_work_j,
        r.external_work_j
            + r.prescribed_surface_work_j
            + r.normal_contact_work_j
            + r.frictional_work_j,
        1e-9,
    );
}
fn admit_steady(
    path: &SolidSustainedPath<'_, '_, SupportedSpectralSolid>,
) -> Result<(), &'static str> {
    // Original M=I, K=0, free applied forces zero, and support acceleration
    // identically zero: velocity is constant over the WHOLE physical interval.
    near(path.total_force_n[0], 0.0, 1e-10);
    near(path.total_force_n[1], 0.0, 1e-10);
    assert_eq!(path.before.support_acceleration_mps2(), [0.0]);
    for i in 0..3 {
        near(
            path.velocity_after_mps[i],
            path.before.velocity_mps()[i],
            1e-10,
        );
        near(
            path.displacement_change_m[i],
            path.duration_s * path.before.velocity_mps()[i],
            1e-10,
        );
    }
    Ok(())
}

#[test]
fn g1_contact_loads_on_prescribed_coordinates_retain_actuator_work_and_boosts() {
    with_cx(|cx, _| {
        let h = 0.125;
        for boost in [[0.0, 0.0], [7.0, -2.0]] {
            let mut s = point(h, boost, cx);
            let mut direct = point(h, boost, cx);
            let r = s
                .advance_sustained_contacts(
                    &[0.0],
                    &[contact(h, 0.5, [0.0, boost[0], boost[1]])],
                    &[[10.0, -3.0, -4.0]],
                    &[3.0, -10.0, 0.0],
                    controls(),
                    cx,
                    admit_steady,
                )
                .unwrap();
            assert_eq!(s.represented_dofs(), [0, 1, 2]);
            assert_eq!(s.free_dofs(), [0, 1]);
            near(r.native.support_reaction_n[0], 4.0, 1e-10);
            near(r.external_work_j, 3.0 * (3.0 + boost[0]) * h, 1e-10);
            near(r.support_work_j, 4.0 * (4.0 + boost[1]) * h, 1e-10);
            near(
                r.prescribed_surface_work_j,
                (-3.0 * boost[0] - 4.0 * boost[1]) * h,
                1e-10,
            );
            near(r.frictional_dissipation_j, 25.0 * h, 1e-10);
            for change in &r.momentum_change_kg_m_s {
                near(*change, 0.0, 1e-10);
            }
            energy(&r);
            direct
                .advance_affine_supports(&[0.0], &[0.0, 0.0, -4.0], &[0.0, 0.0, -4.0], cx)
                .unwrap();
            assert_eq!(bits(&s), bits(&direct));
            s.advance_affine_supports(&[0.5], &[1.0, 0.0, 0.3], &[-2.0, 0.0, 1.5], cx)
                .unwrap();
            direct
                .advance_affine_supports(&[0.5], &[1.0, 0.0, 0.3], &[-2.0, 0.0, 1.5], cx)
                .unwrap();
            assert_eq!(bits(&s), bits(&direct));
        }
    });
}

fn accelerating(h: f64, cx: &Cx<'_>) -> SupportedSpectralSolid {
    // Coordinates: first node's X/Y, and a prescribed second node's X.
    let mut source = assembly(
        &[4.0, 0.0, -4.0, 0.0, 0.0, 0.0, -4.0, 0.0, 4.0],
        &[1.0, 0.0, 0.5, 0.0, 1.0, 0.0, 0.5, 0.0, 2.0],
    );
    source.free_dofs = vec![0, 1, 3];
    SupportedSpectralSolid::prepare_with_basis(
        &source,
        &[3],
        &[4.0, 0.0],
        &BASIS,
        h,
        &[0.0; 3],
        &[2.0, 0.0, 0.5],
        &[1.0],
        native_controls(),
        cx,
    )
    .unwrap()
}
#[test]
fn g1_accelerating_support_work_uses_original_cross_mass_and_reaction_integral() {
    with_cx(|cx, _| {
        let h = 0.125_f64;
        let (sin, cos) = (2.0 * h).sin_cos();
        let z = 0.5 * h + 0.5 * h * h;
        let x = z + 0.125 * (1.0 - cos) + 0.75 * sin;
        let vx = 0.5 + h + 0.25 * sin + 1.5 * cos;
        let vz = 0.5 + h;
        let cosine_moment = 0.5 * sin / 2.0 + h * sin / 2.0 + (cos - 1.0) / 4.0;
        let sine_moment = 0.5 * (1.0 - cos) / 2.0 - h * cos / 2.0 + sin / 4.0;
        let zero = [0.0; 3];
        let c = SolidSustainedContact {
            velocity_weights: [&N, &X, &zero],
            ..contact(h, 0.2, [0.0; 3])
        };
        for applied_support_load in [0.0, 3.0] {
            let mut s = accelerating(h, cx);
            let mut direct = accelerating(h, cx);
            let r = s
                .advance_sustained_contacts(
                    &[1.0],
                    &[c],
                    &[[10.0, -2.0, 0.0]],
                    &[4.0, -10.0, applied_support_load],
                    controls(),
                    cx,
                    |p| {
                        // x'=.5+t+.25*sin(2t)+1.5*cos(2t). The original
                        // |x''|<=1+.5+3 bounds x'>2-4.5*h>0 throughout.
                        assert!(2.0 - 4.5 * h > 0.0);
                        near(p.total_force_n[0], 2.0, 1e-12);
                        near(p.total_force_n[1], 0.0, 1e-12);
                        assert_eq!(p.before.support_acceleration_mps2(), [1.0]);
                        Ok::<_, &'static str>(())
                    },
                )
                .unwrap();
            near(s.displacement_m()[0], x, 1e-12);
            near(s.displacement_m()[2], z, 1e-13);
            near(s.velocity_mps()[0], vx, 1e-12);
            near(s.velocity_mps()[2], vz, 1e-13);
            // Rz(t)=2-Fz+.75*cos(2t)-4.5*sin(2t), integrated against
            // z'=.5+t. This witness does not use a body-energy difference.
            let support_work =
                (2.0 - applied_support_load) * z + 0.75 * cosine_moment - 4.5 * sine_moment;
            near(
                r.native.support_reaction_n[0],
                2.0 - applied_support_load + 0.75 * cos - 4.5 * sin,
                1e-12,
            );
            near(r.support_work_j, support_work, 1e-12);
            near(r.external_work_j, 4.0 * x + applied_support_load * z, 1e-12);
            near(
                r.native.energy.applied_work_j,
                2.0 * x + applied_support_load * z,
                1e-12,
            );
            near(r.frictional_dissipation_j, 2.0 * x, 1e-12);
            near(
                r.momentum_change_kg_m_s[0],
                vx - 2.0 + 0.5 * (vz - 0.5),
                1e-12,
            );
            near(
                r.momentum_change_kg_m_s[2],
                0.5 * (vx - 2.0) + 2.0 * (vz - 0.5),
                1e-12,
            );
            energy(&r);
            let load = [2.0, 0.0, applied_support_load];
            direct
                .advance_affine_supports(&[1.0], &load, &load, cx)
                .unwrap();
            assert_eq!(bits(&s), bits(&direct));
        }
    });
}

#[test]
fn g4_supported_trials_preserve_both_owners_on_late_rejection_and_cancelled_retry() {
    with_cx(|cx, gate| {
        let h = 0.125;
        let mut s = point(h, [0.0; 2], cx);
        let before = bits(&s);
        let contacts = [contact(h, 0.5, [0.0; 3])];
        let initial = [[10.0, -3.0, -4.0]];
        for config in [
            SolidSustainedControls {
                maximum_workspace_bytes: 0,
                ..controls()
            },
            SolidSustainedControls {
                maximum_response_multiply_terms: 0,
                ..controls()
            },
        ] {
            assert!(
                s.advance_sustained_contacts(
                    &[0.0],
                    &contacts,
                    &initial,
                    &[3.0, -10.0, 0.0],
                    config,
                    cx,
                    admit_steady
                )
                .is_err()
            );
            assert_eq!(bits(&s), before);
        }
        assert!(
            s.advance_sustained_contacts(
                &[],
                &contacts,
                &initial,
                &[3.0, -10.0, 0.0],
                controls(),
                cx,
                admit_steady
            )
            .is_err()
        );
        assert_eq!(bits(&s), before);
        let result = s.advance_sustained_contacts(
            &[0.0],
            &contacts,
            &initial,
            &[3.0, -10.0, 0.0],
            controls(),
            cx,
            |_| Err::<(), _>("interior support geometry refused"),
        );
        assert!(matches!(
            result,
            Err(SolidSustainedError::Path(
                "interior support geometry refused"
            ))
        ));
        assert_eq!(bits(&s), before);
        let result = s.advance_sustained_contacts(
            &[0.0],
            &contacts,
            &initial,
            &[3.0, -10.0, 0.0],
            controls(),
            cx,
            |p| {
                admit_steady(p)?;
                gate.request();
                Ok::<_, &'static str>(())
            },
        );
        assert!(matches!(result, Err(SolidSustainedError::Cancelled)));
        assert_eq!(bits(&s), before);
        with_cx(|fresh, _| {
            let r = s
                .advance_sustained_contacts(
                    &[0.0],
                    &contacts,
                    &initial,
                    &[3.0, -10.0, 0.0],
                    controls(),
                    fresh,
                    admit_steady,
                )
                .unwrap();
            let mut direct = point(h, [0.0; 2], fresh);
            direct
                .advance_affine_supports(&[0.0], &[0.0, 0.0, -4.0], &[0.0, 0.0, -4.0], fresh)
                .unwrap();
            assert_eq!(bits(&s), bits(&direct));
            energy(&r);
            s.advance_affine_supports(&[-0.4], &[1.0, 0.0, 0.3], &[0.0, 0.0, -0.7], fresh)
                .unwrap();
            direct
                .advance_affine_supports(&[-0.4], &[1.0, 0.0, 0.3], &[0.0, 0.0, -0.7], fresh)
                .unwrap();
            assert_eq!(bits(&s), bits(&direct));
        });
    });
}

#[test]
fn g4_original_full_supported_work_gate_remains_authoritative() {
    with_cx(|cx, _| {
        let h = 0.125;
        let mut cfg = native_controls();
        cfg.spectral.maximum_relative_operator_error = 1e-2;
        cfg.spectral.maximum_relative_state_error = 1e-2;
        cfg.spectral.maximum_work_error_j = 1e-2;
        cfg.maximum_work_error_j = 1e-14;
        let mut basis = BASIS;
        basis[0] = 1.0001;
        let mut s = SupportedSpectralSolid::prepare_with_basis(
            &assembly(&[0.0; 9], &I),
            &[2],
            &[0.0; 2],
            &basis,
            h,
            &[0.0; 3],
            &[3.0, 0.0, 4.0],
            &[0.0],
            cfg,
            cx,
        )
        .unwrap();
        let before = bits(&s);
        let result = s.advance_sustained_contacts(
            &[0.0],
            &[contact(h, 0.5, [0.0; 3])],
            &[[10.0, -3.0, -4.0]],
            &[3.0, -10.0, 0.0],
            SolidSustainedControls {
                maximum_work_error_j: 0.1,
                ..controls()
            },
            cx,
            |_| Ok::<_, &'static str>(()),
        );
        match result {
            Err(SolidSustainedError::Solve(error)) => assert!(matches!(
                error.as_ref(),
                SolveError::Response(SupportedSolidError::Energy { .. })
            )),
            other => panic!("expected original whole-supported energy refusal, got {other:?}"),
        }
        assert_eq!(bits(&s), before);
    });
}
