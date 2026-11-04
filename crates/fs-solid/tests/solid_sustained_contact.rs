//! Complete native trajectories under finite contact forces. Independent point
//! solutions and an elastic/damped closed form check work and actual motion.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::sustained::NaturalMapScales;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::sustained::{
    SolidSustainedContact, SolidSustainedControls, SolidSustainedError, SolidSustainedPath,
    SolidSustainedStep,
};
use fs_solid::dynamics::{ExactLinearSolid, ExactLinearSolidControls};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::Coo;

const N: [f64; 3] = [0.0, 1.0, 0.0];
const T1: [f64; 3] = [1.0, 0.0, 0.0];
const T2: [f64; 3] = [0.0, 0.0, 1.0];

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 62,
                    kernel_id: 73,
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
fn native_controls() -> ExactLinearSolidControls {
    ExactLinearSolidControls {
        maximum_free_dofs: 24,
        maximum_preparation_bytes: 64 * 1024 * 1024,
        maximum_exponential_multiply_terms: 1_000_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_work_error_j: 1e-9,
        maximum_steps: 1000,
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
fn assembly(k: &[f64], m: &[f64], n: usize) -> TetElasticAssembly {
    let mut stiffness = Coo::new(n, n);
    let mut mass = Coo::new(n, n);
    for i in 0..n {
        for j in 0..n {
            stiffness.push(i, j, k[i * n + j]);
            mass.push(i, j, m[i * n + j]);
        }
    }
    TetElasticAssembly {
        stiffness: stiffness.assemble(),
        mass: mass.assemble(),
        free_dofs: (0..n).collect(),
        total_mass_kg: 1.0,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0],
    }
}
fn point(h: f64, velocity: [f64; 3], cx: &Cx<'_>) -> ExactLinearSolid {
    ExactLinearSolid::prepare(
        &assembly(&[0.0; 9], &[2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0], 3),
        h,
        &[0.0; 3],
        &velocity,
        native_controls(),
        cx,
    )
    .unwrap()
}
fn contact(h: f64, mu: f64, surface: [f64; 3]) -> SolidSustainedContact<'static> {
    SolidSustainedContact {
        id: 7,
        velocity_weights: [&N, &T1, &T2],
        gap_m: 0.0,
        prescribed_surface_velocity_mps: surface,
        friction_coefficient: mu,
        scales: NaturalMapScales {
            normal_n_per_m: 4.0 / (h * h),
            tangent_n_per_m: 4.0 / (h * h),
        },
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e} (tol={tol:e})");
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn bits(s: &ExactLinearSolid) -> Vec<u64> {
    s.displacement_m()
        .iter()
        .chain(s.velocity_mps())
        .copied()
        .chain([
            s.time_s(),
            s.cumulative_work_j(),
            s.cumulative_dissipated_energy_j(),
            s.stored_energy_j(),
            s.kinetic_energy_j(),
        ])
        .map(f64::to_bits)
        .collect()
}
// Independent whole-path owner for a mass-two free point. The normal gap is
// quadratic; slip velocity is affine. Endpoint-only numerical force closure
// is insufficient, so the interior normal minimum is evaluated explicitly.
fn admit_point(path: &SolidSustainedPath<'_, '_>) -> Result<&'static str, &'static str> {
    let h = path.duration_s;
    let acceleration: Vec<_> = path.total_force_n.iter().map(|f| f / 2.0).collect();
    for i in 0..3 {
        near(
            path.velocity_after_mps[i],
            path.before.velocity_mps()[i] + h * acceleration[i],
            1e-10,
        );
        near(
            path.displacement_change_m[i],
            h * path.before.velocity_mps()[i] + 0.5 * h * h * acceleration[i],
            1e-11,
        );
    }
    for (c, f) in path.contacts.iter().zip(path.forces_n) {
        let vn = dot(c.velocity_weights[0], path.before.velocity_mps())
            - c.prescribed_surface_velocity_mps[0];
        let an = dot(c.velocity_weights[0], &acceleration);
        let gap = |t: f64| c.gap_m + vn * t + 0.5 * an * t * t;
        let mut minimum = gap(0.0).min(gap(h));
        if an > 0.0 {
            let t = -vn / an;
            if t > 0.0 && t < h {
                minimum = minimum.min(gap(t));
            }
        }
        if minimum < -1e-11 {
            return Err("interior penetration");
        }
        if f[0] > 1e-7 && (vn.abs() > 1e-9 || an.abs() > 1e-8 || c.gap_m > 1e-10) {
            return Err("loaded point does not dwell");
        }
        let cap = c.friction_coefficient * f[0];
        if cap > 1e-7 {
            let ft = f[1].hypot(f[2]);
            for t in [0.0, h] {
                let vt: [f64; 2] = core::array::from_fn(|a| {
                    dot(c.velocity_weights[a + 1], path.before.velocity_mps())
                        - c.prescribed_surface_velocity_mps[a + 1]
                        + t * dot(c.velocity_weights[a + 1], &acceleration)
                });
                if ft < cap - 1e-7 && vt[0].hypot(vt[1]) > 1e-9 {
                    return Err("not sticking");
                }
                if f[1] * vt[0] + f[2] * vt[1] > 1e-9 {
                    return Err("sliding reverses");
                }
            }
        }
    }
    Ok("analytic continuous point path")
}
fn energy<P>(s: &SolidSustainedStep<P>) {
    near(s.native.energy.balance_error_j, 0.0, 1e-9);
    near(s.native.energy.cumulative_balance_error_j, 0.0, 1e-9);
    near(s.contact_work_balance_error_j, 0.0, 1e-9);
    near(s.energy_balance_error_j, 0.0, 1e-9);
    near(s.normal_contact_work_j, 0.0, 1e-9);
    near(s.passivity_residual_j, 0.0, 1e-9);
    near(
        s.native.energy.supplied_work_j,
        s.external_work_j
            + s.prescribed_surface_work_j
            + s.normal_contact_work_j
            + s.frictional_work_j,
        1e-9,
    );
}

#[test]
fn g1_resting_sticking_and_breakaway_use_continuous_forces_without_impact_loss() {
    with_cx(|cx, _| {
        for load in [0.0, 5.0, 10.0, 12.0] {
            let h = 0.125;
            let mut s = point(h, [0.0; 3], cx);
            let r = s
                .advance_sustained_contacts(
                    &[contact(h, 0.5, [0.0; 3])],
                    &[[0.0; 3]],
                    &[load, -20.0, 0.0],
                    controls(),
                    cx,
                    admit_point,
                )
                .unwrap();
            let ax = (load - 10.0_f64).max(0.0) / 2.0;
            near(r.solve.forces_n[0][0], 20.0, 1e-8);
            near(r.solve.forces_n[0][1], -load.min(10.0), 1e-8);
            near(s.velocity_mps()[0], ax * h, 1e-9);
            near(s.displacement_m()[0], 0.5 * ax * h * h, 1e-10);
            near(r.frictional_dissipation_j, 10.0 * 0.5 * ax * h * h, 1e-9);
            assert_eq!(r.native.energy.dissipated_j, 0.0);
            assert_eq!(s.accepted_steps(), 1);
            assert_eq!(r.contact_ids, [7]);
            assert!(r.solve.response_evaluations >= 2);
            energy(&r);
        }
    });
}

#[test]
fn g1_two_tangent_sliding_has_circular_coulomb_work_and_no_axis_clipping() {
    with_cx(|cx, _| {
        let h = 0.125;
        let mut s = point(h, [1.5, 0.0, 2.0], cx);
        let r = s
            .advance_sustained_contacts(
                &[contact(h, 0.5, [0.0; 3])],
                &[[20.0, 0.0, 0.0]],
                &[0.0, -20.0, 0.0],
                controls(),
                cx,
                admit_point,
            )
            .unwrap();
        for (a, b) in r.solve.forces_n[0].iter().zip([20.0, -6.0, -8.0]) {
            near(*a, b, 1e-8);
        }
        near(s.velocity_mps()[0], 1.125, 1e-9);
        near(s.velocity_mps()[2], 1.5, 1e-9);
        near(s.displacement_m()[0], 0.1640625, 1e-10);
        near(s.displacement_m()[2], 0.21875, 1e-10);
        near(r.frictional_dissipation_j, 2.734375, 1e-9);
        near(
            r.native.energy.kinetic_before_j - r.native.energy.kinetic_after_j,
            2.734375,
            1e-9,
        );
        energy(&r);
    });
}

#[test]
fn g1_moving_surfaces_supply_signed_work_separate_from_relative_friction() {
    with_cx(|cx, _| {
        let h = 0.125;
        for boost in [0.0, 7.0] {
            let mut s = point(h, [boost, 1.0, 0.0], cx);
            let r = s
                .advance_sustained_contacts(
                    &[contact(h, 0.5, [1.0, 2.0 + boost, 0.0])],
                    &[[20.0, 10.0, 0.0]],
                    &[0.0, -20.0, 0.0],
                    controls(),
                    cx,
                    admit_point,
                )
                .unwrap();
            near(s.velocity_mps()[0], boost + 0.625, 1e-10);
            near(s.velocity_mps()[1], 1.0, 1e-10);
            near(r.external_work_j, -2.5, 1e-10);
            near(r.prescribed_surface_work_j, 5.0 + 1.25 * boost, 1e-10);
            near(r.frictional_dissipation_j, 2.109375, 1e-9);
            near(
                r.native.energy.supplied_work_j,
                0.390625 + 1.25 * boost,
                1e-9,
            );
            assert_eq!(s.cumulative_dissipated_energy_j(), 0.0);
            energy(&r);
        }
    });
}

#[test]
fn g1_consistent_cross_mass_changes_normal_load_and_keeps_full_momentum() {
    with_cx(|cx, _| {
        let h = 0.02;
        let m = [2.0, 0.5, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 3.0];
        let mut s = ExactLinearSolid::prepare(
            &assembly(&[0.0; 9], &m, 3),
            h,
            &[0.0; 3],
            &[1.0, 0.0, 0.0],
            native_controls(),
            cx,
        )
        .unwrap();
        let r = s
            .advance_sustained_contacts(
                &[contact(h, 0.2, [0.0; 3])],
                &[[20.0, -4.0, 0.0]],
                &[12.0, -20.0, 0.0],
                controls(),
                cx,
                |p| {
                    // 2 ax=12-mu*(20+ax/2), ay=az=0 from the complete
                    // original [[2,.5],[.5,1]] inertia, not a diagonal proxy.
                    let a = 80.0 / 21.0;
                    near(p.velocity_after_mps[0], 1.0 + h * a, 1e-9);
                    near(p.velocity_after_mps[1], 0.0, 1e-9);
                    near(p.displacement_change_m[0], h + 0.5 * h * h * a, 1e-11);
                    Ok::<_, &'static str>("constant positive slip and zero normal acceleration")
                },
            )
            .unwrap();
        near(r.solve.forces_n[0][0], 460.0 / 21.0, 1e-8);
        near(r.solve.forces_n[0][1], -92.0 / 21.0, 1e-8);
        near(r.momentum_change_kg_m_s[0], 2.0 * h * 80.0 / 21.0, 1e-9);
        near(r.momentum_change_kg_m_s[1], 0.5 * h * 80.0 / 21.0, 1e-9);
        energy(&r);
    });
}

#[test]
fn g3_simultaneous_contacts_share_original_inertia_and_are_order_independent() {
    with_cx(|cx, _| {
        let h = 0.125;
        let mut mass = [0.0; 36];
        for a in 0..2 {
            for b in 0..2 {
                for axis in 0..3 {
                    mass[(3 * a + axis) * 6 + 3 * b + axis] = if a == b { 2.0 } else { 1.0 };
                }
            }
        }
        let pencil = assembly(&[0.0; 36], &mass, 6);
        let rows = [
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        ];
        let c = [
            SolidSustainedContact {
                id: 3,
                velocity_weights: [&rows[0], &rows[1], &rows[2]],
                friction_coefficient: 0.1,
                ..contact(h, 0.0, [0.0; 3])
            },
            SolidSustainedContact {
                id: 9,
                velocity_weights: [&rows[3], &rows[4], &rows[5]],
                friction_coefficient: 0.2,
                ..contact(h, 0.0, [0.0; 3])
            },
        ];
        let mut baseline = None;
        for reverse in [false, true] {
            let mut s = ExactLinearSolid::prepare(
                &pencil,
                h,
                &[0.0; 6],
                &[1.0, 0.0, 0.0, 2.0, 0.0, 0.0],
                native_controls(),
                cx,
            )
            .unwrap();
            let contacts = if reverse { [c[1], c[0]] } else { c };
            let r = s
                .advance_sustained_contacts(
                    &contacts,
                    &[[0.0; 3]; 2],
                    &[0.0, -10.0, 0.0, 0.0, -20.0, 0.0],
                    controls(),
                    cx,
                    |p| {
                        for (j, v, a) in [(0, 1.0, 2.0 / 3.0), (3, 2.0, -7.0 / 3.0)] {
                            near(p.velocity_after_mps[j], v + h * a, 1e-9);
                            near(p.displacement_change_m[j], h * v + 0.5 * h * h * a, 1e-10);
                            assert!(v + h * a > 0.0);
                        }
                        Ok::<_, &'static str>("two original coupled point paths")
                    },
                )
                .unwrap();
            assert_eq!(r.contact_ids, [3, 9]);
            near(r.solve.forces_n[0][1], -1.0, 1e-8);
            near(r.solve.forces_n[1][1], -4.0, 1e-8);
            energy(&r);
            if let Some(ref previous) = baseline {
                assert_eq!(bits(&s), *previous);
            } else {
                baseline = Some(bits(&s));
            }
        }
    });
}

#[test]
fn g1_elastic_viscous_motion_retains_exact_native_flow_and_separate_losses() {
    with_cx(|cx, _| {
        let h = 0.01;
        let k = [8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let m = [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0];
        let c = [0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let pencil = assembly(&k, &m, 3);
        let viscosity = assembly(&[0.0; 9], &c, 3).mass;
        let make = || {
            ExactLinearSolid::prepare_viscous(
                &pencil,
                &viscosity,
                h,
                &[0.0; 3],
                &[1.0, 0.0, 0.0],
                native_controls(),
                cx,
            )
            .unwrap()
        };
        let mut s = make();
        let mut original = make();
        let r = s
            .advance_sustained_contacts(
                &[contact(h, 0.1, [0.0; 3])],
                &[[20.0, -2.0, 0.0]],
                &[4.0, -20.0, 0.0],
                controls(),
                cx,
                |p| {
                    // A continuous bootstrap: |v|<2 implies |u|<2h and
                    // |a| <= (|F|+c*2+k*2h)/m. Its integrated bound stays
                    // strictly inside |v|<2 and above zero throughout [0,h].
                    let excursion = (2.0 + 0.4 * 2.0 + 8.0 * 2.0 * h) / 2.0 * h;
                    assert!(1.0 + excursion < 2.0 && 1.0 - excursion > 0.0);
                    near(p.total_force_n[1], 0.0, 1e-12);
                    Ok::<_, &'static str>("continuous no-reversal bootstrap")
                },
            )
            .unwrap();
        original.advance(&[2.0, 0.0, 0.0], cx).unwrap();
        assert_eq!(bits(&s), bits(&original));
        let beta = 3.99_f64.sqrt();
        let exact = |t: f64| {
            let (sin, cos) = (beta * t).sin_cos();
            let a = -0.25;
            let b = 0.975 / beta;
            let decay = (-0.1 * t).exp();
            (
                0.25 + decay * (a * cos + b * sin),
                decay * (-0.1 * (a * cos + b * sin) + beta * (-a * sin + b * cos)),
            )
        };
        near(s.displacement_m()[0], exact(h).0, 1e-13);
        near(s.velocity_mps()[0], exact(h).1, 1e-13);
        let mut quadrature = 0.0;
        for i in 0..=128 {
            let weight = if i == 0 || i == 128 {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            quadrature += weight * 0.4 * exact(h * i as f64 / 128.0).1.powi(2);
        }
        quadrature *= h / (3.0 * 128.0);
        near(r.native.energy.dissipated_j, quadrature, 1e-13);
        near(r.frictional_dissipation_j, 2.0 * exact(h).0, 1e-12);
        near(r.external_work_j, 4.0 * exact(h).0, 1e-12);
        assert!(r.native.energy.dissipated_j > 0.0 && r.frictional_dissipation_j > 0.0);
        energy(&r);
        s.advance(&[-1.0, 0.0, 0.0], cx).unwrap();
        original.advance(&[-1.0, 0.0, 0.0], cx).unwrap();
        assert_eq!(bits(&s), bits(&original));
    });
}

#[test]
fn g0_impacts_and_unresolved_sliding_stops_refuse_without_creating_work() {
    with_cx(|cx, _| {
        let h = 0.125;
        for (v, impact) in [([0.0, -1.0, 0.0], true), ([0.1, 0.0, 0.0], false)] {
            let mut s = point(h, v, cx);
            let before = bits(&s);
            let result = s.advance_sustained_contacts(
                &[contact(h, 0.5, [0.0; 3])],
                &[[20.0, -3.2, 0.0]],
                &[0.0, -20.0, 0.0],
                controls(),
                cx,
                admit_point,
            );
            if impact {
                assert!(matches!(
                    result,
                    Err(SolidSustainedError::ImpactRequired { .. })
                ));
            } else {
                assert!(
                    result.is_err(),
                    "a hidden stop/reversal must not become sustained sticking"
                );
            }
            assert_eq!(bits(&s), before);
        }
        let mut s = point(h, [0.0; 3], cx);
        let r = s
            .advance_sustained_contacts(
                &[contact(h, 0.5, [0.0; 3])],
                &[[20.0, 0.0, 0.0]],
                &[0.0, 6.0, 0.0],
                controls(),
                cx,
                admit_point,
            )
            .unwrap();
        near(r.solve.forces_n[0][0], 0.0, 1e-10);
        near(s.velocity_mps()[1], 3.0 * h, 1e-10);
        near(r.frictional_dissipation_j, 0.0, 1e-10);
        energy(&r);
    });
}

#[test]
fn g4_path_refusal_cancellation_and_native_budgets_leave_exact_retry() {
    with_cx(|cx, gate| {
        let h = 0.125;
        let mut s = point(h, [1.0, 0.0, 0.0], cx);
        let before = bits(&s);
        let contacts = [contact(h, 0.5, [0.0; 3])];
        for config in [
            SolidSustainedControls {
                maximum_workspace_bytes: 0,
                ..controls()
            },
            SolidSustainedControls {
                maximum_response_multiply_terms: 0,
                ..controls()
            },
            SolidSustainedControls {
                solve: fs_contact::sustained::solve::SolveControls {
                    maximum_response_evaluations: 1,
                    ..Default::default()
                },
                ..controls()
            },
        ] {
            assert!(
                s.advance_sustained_contacts(
                    &contacts,
                    &[[0.0; 3]],
                    &[12.0, -20.0, 0.0],
                    config,
                    cx,
                    admit_point
                )
                .is_err()
            );
            assert_eq!(bits(&s), before);
        }
        let rejected = s.advance_sustained_contacts(
            &contacts,
            &[[20.0, -10.0, 0.0]],
            &[12.0, -20.0, 0.0],
            controls(),
            cx,
            |_| Err::<(), _>("whole-path geometry refused"),
        );
        assert!(matches!(
            rejected,
            Err(SolidSustainedError::Path("whole-path geometry refused"))
        ));
        assert_eq!(bits(&s), before);
        let cancelled = s.advance_sustained_contacts(
            &contacts,
            &[[20.0, -10.0, 0.0]],
            &[12.0, -20.0, 0.0],
            controls(),
            cx,
            |p| {
                admit_point(p)?;
                gate.request();
                Ok::<_, &'static str>(())
            },
        );
        assert!(matches!(cancelled, Err(SolidSustainedError::Cancelled)));
        assert_eq!(bits(&s), before);
        with_cx(|fresh, _| {
            let r = s
                .advance_sustained_contacts(
                    &contacts,
                    &[[20.0, -10.0, 0.0]],
                    &[12.0, -20.0, 0.0],
                    controls(),
                    fresh,
                    admit_point,
                )
                .unwrap();
            let mut reference = point(h, [1.0, 0.0, 0.0], fresh);
            reference.advance(&[2.0, 0.0, 0.0], fresh).unwrap();
            assert_eq!(bits(&s), bits(&reference));
            energy(&r);
        });
    });
}
