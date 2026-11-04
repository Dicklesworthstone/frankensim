//! Sustained contact on every original spectral coordinate. Closed-form rigid
//! and elastic trajectories independently check signed physical work.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::sustained::NaturalMapScales;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::spectral::{
    CompleteSpectralSolid, CompleteSpectralSolidControls, CompleteSpectralSolidError,
};
use fs_solid::dynamics::sustained::{
    SolidSustainedContact, SolidSustainedControls, SolidSustainedError, SolidSustainedPath,
    SolidSustainedStep,
};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::Coo;

const N: [f64; 3] = [0.0, 1.0, 0.0];
const T1: [f64; 3] = [1.0, 0.0, 0.0];
const T2: [f64; 3] = [0.0, 0.0, 1.0];
const I: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 62,
                    kernel_id: 91,
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
fn native_controls() -> CompleteSpectralSolidControls {
    CompleteSpectralSolidControls {
        maximum_free_dofs: 24,
        maximum_preparation_bytes: 64 * 1024 * 1024,
        maximum_preparation_multiply_terms: 1_000_000_000,
        maximum_step_multiply_terms: 10_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_absolute_eigen_residual: 0.0,
        maximum_relative_state_error: 1e-12,
        maximum_work_error_j: 1e-10,
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
fn point(h: f64, velocity: [f64; 3], cx: &Cx<'_>) -> CompleteSpectralSolid {
    CompleteSpectralSolid::prepare_with_basis(
        &assembly(&[0.0; 9], &I, 3),
        &[0.0; 3],
        &I,
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
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn bits(s: &CompleteSpectralSolid) -> Vec<u64> {
    s.displacement_m()
        .iter()
        .chain(s.velocity_mps())
        .copied()
        .chain([
            s.time_s(),
            s.cumulative_work_j(),
            s.stored_energy_j(),
            s.kinetic_energy_j(),
        ])
        .map(f64::to_bits)
        .collect()
}
fn energy<P>(s: &SolidSustainedStep<P>) {
    near(s.energy_balance_error_j, 0.0, 1e-9);
    near(s.contact_work_balance_error_j, 0.0, 1e-9);
    near(s.native.energy.balance_error_j, 0.0, 1e-9);
    near(s.native.energy.cumulative_balance_error_j, 0.0, 1e-9);
    near(
        s.native.energy.supplied_work_j,
        s.external_work_j
            + s.prescribed_surface_work_j
            + s.normal_contact_work_j
            + s.frictional_work_j,
        1e-9,
    );
    assert_eq!(s.native.energy.dissipated_j, 0.0);
    assert_eq!(s.native.energy.cumulative_dissipated_j, 0.0);
}

// For a complete zero-K pencil, acceleration is constant. Checking the vertex
// of each quadratic gap and both ends of each affine tangential velocity
// certifies the entire interval, including the absence of hidden slip stops.
fn admit_rigid(
    path: &SolidSustainedPath<'_, '_, CompleteSpectralSolid>,
    inverse_mass: &[f64],
) -> Result<(), &'static str> {
    let n = path.total_force_n.len();
    let h = path.duration_s;
    let a: Vec<_> = (0..n)
        .map(|i| dot(&inverse_mass[i * n..(i + 1) * n], path.total_force_n))
        .collect();
    for i in 0..n {
        near(
            path.velocity_after_mps[i],
            path.before.velocity_mps()[i] + h * a[i],
            1e-9,
        );
        near(
            path.displacement_change_m[i],
            h * path.before.velocity_mps()[i] + 0.5 * h * h * a[i],
            1e-10,
        );
    }
    for (c, f) in path.contacts.iter().zip(path.forces_n) {
        let vn = dot(c.velocity_weights[0], path.before.velocity_mps())
            - c.prescribed_surface_velocity_mps[0];
        let an = dot(c.velocity_weights[0], &a);
        let gap = |t: f64| c.gap_m + vn * t + 0.5 * an * t * t;
        let mut minimum = gap(0.0).min(gap(h));
        if an > 0.0 {
            let t = -vn / an;
            if t > 0.0 && t < h {
                minimum = minimum.min(gap(t));
            }
        }
        if minimum < -1e-10 {
            return Err("interior penetration");
        }
        if f[0] > 1e-7 && (c.gap_m > 1e-10 || vn.abs() > 1e-9 || an.abs() > 1e-8) {
            return Err("loaded trace does not dwell");
        }
        let cap = c.friction_coefficient * f[0];
        if cap > 1e-7 {
            for t in [0.0, h] {
                let v: [f64; 2] = core::array::from_fn(|j| {
                    dot(c.velocity_weights[j + 1], path.before.velocity_mps())
                        - c.prescribed_surface_velocity_mps[j + 1]
                        + t * dot(c.velocity_weights[j + 1], &a)
                });
                if f[1].hypot(f[2]) < cap - 1e-7 && v[0].hypot(v[1]) > 1e-9 {
                    return Err("not sticking");
                }
                if (f[1] * v[0] + f[2] * v[1] + cap * v[0].hypot(v[1])).abs() > 1e-8 {
                    return Err("sliding direction reverses");
                }
            }
        }
    }
    Ok(())
}

#[test]
fn g1_spectral_sticking_breakaway_and_circular_sliding_retain_original_work() {
    with_cx(|cx, _| {
        let h = 0.125;
        for load in [0.0_f64, 3.0, 4.0, 6.0] {
            let mut s = point(h, [0.0; 3], cx);
            let r = s
                .advance_sustained_contacts(
                    &[contact(h, 0.5, [0.0; 3])],
                    &[[8.0, -load.min(4.0), 0.0]],
                    &[load, -8.0, 0.0],
                    controls(),
                    cx,
                    |p| admit_rigid(p, &I),
                )
                .unwrap();
            let acceleration = (load - 4.0_f64).max(0.0);
            near(r.solve.forces_n[0][0], 8.0, 1e-8);
            near(r.solve.forces_n[0][1], -load.min(4.0), 1e-8);
            near(s.velocity_mps()[0], acceleration * h, 1e-9);
            near(s.displacement_m()[0], 0.5 * acceleration * h * h, 1e-10);
            near(
                r.frictional_dissipation_j,
                4.0 * s.displacement_m()[0],
                1e-9,
            );
            energy(&r);
        }
        let mut s = point(h, [1.5, 0.0, 2.0], cx);
        let r = s
            .advance_sustained_contacts(
                &[contact(h, 0.5, [0.0; 3])],
                &[[8.0, -2.4, -3.2]],
                &[0.0, -8.0, 0.0],
                controls(),
                cx,
                |p| admit_rigid(p, &I),
            )
            .unwrap();
        near(s.velocity_mps()[0], 1.2, 1e-9);
        near(s.velocity_mps()[2], 1.6, 1e-9);
        near(r.frictional_dissipation_j, 1.125, 1e-9);
        energy(&r);
    });
}

#[test]
fn g1_moving_surface_and_common_boost_keep_separate_signed_work() {
    with_cx(|cx, _| {
        let h = 0.125;
        for boost in [0.0, 5.0, -3.0] {
            let mut s = point(h, [2.0 + boost, 0.75, 0.0], cx);
            let r = s
                .advance_sustained_contacts(
                    &[contact(h, 0.5, [0.75, boost, 0.0])],
                    &[[8.0, -4.0, 0.0]],
                    &[0.0, -8.0, 0.0],
                    controls(),
                    cx,
                    |p| admit_rigid(p, &I),
                )
                .unwrap();
            near(s.velocity_mps()[0], 1.5 + boost, 1e-10);
            near(s.velocity_mps()[1], 0.75, 1e-10);
            near(r.external_work_j, -0.75, 1e-10);
            near(r.prescribed_surface_work_j, 0.75 - 0.5 * boost, 1e-10);
            near(r.frictional_dissipation_j, 0.875, 1e-10);
            near(r.momentum_change_kg_m_s[0], -0.5, 1e-10);
            near(r.native.energy.supplied_work_j, -0.875 - 0.5 * boost, 1e-10);
            energy(&r);
        }
    });
}

#[test]
fn g1_full_consistent_cross_mass_controls_normal_force_and_momentum() {
    with_cx(|cx, _| {
        // Phi is complete and M-orthonormal with exact integer coefficients.
        // M_xy=-1: constrained ay=0 requires Fn=4 and Ft=-2 under Fy=-3.
        let m = [2.0, -1.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let inverse = [1.0, 1.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 1.0];
        let basis = [1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let h = 0.25;
        let mut s = CompleteSpectralSolid::prepare_with_basis(
            &assembly(&[0.0; 9], &m, 3),
            &[0.0; 3],
            &basis,
            h,
            &[0.0; 3],
            &[2.0, 0.0, 0.0],
            native_controls(),
            cx,
        )
        .unwrap();
        let r = s
            .advance_sustained_contacts(
                &[contact(h, 0.5, [0.0; 3])],
                &[[4.0, -2.0, 0.0]],
                &[0.0, -3.0, 0.0],
                controls(),
                cx,
                |p| admit_rigid(p, &inverse),
            )
            .unwrap();
        near(r.solve.forces_n[0][0], 4.0, 1e-9);
        near(r.solve.forces_n[0][1], -2.0, 1e-9);
        near(s.velocity_mps()[0], 1.75, 1e-10);
        near(s.velocity_mps()[1], 0.0, 1e-10);
        near(r.momentum_change_kg_m_s[0], -0.5, 1e-10);
        near(r.momentum_change_kg_m_s[1], 0.25, 1e-10);
        near(r.frictional_dissipation_j, 0.9375, 1e-10);
        energy(&r);
    });
}

#[test]
fn g1_every_elastic_and_rigid_mode_survives_contact_and_native_continuation() {
    with_cx(|cx, _| {
        let k = [4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0];
        let h = 0.125;
        let make = |basis: &[f64], values: &[f64]| {
            CompleteSpectralSolid::prepare_with_basis(
                &assembly(&k, &I, 3),
                values,
                basis,
                h,
                &[0.0, 0.0, 0.4],
                &[2.0, 0.0, 0.0],
                native_controls(),
                cx,
            )
            .unwrap()
        };
        let mut a = make(&I, &[4.0, 0.0, 9.0]);
        let mut reference = make(&I, &[4.0, 0.0, 9.0]);
        let mut reordered = make(
            &[0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
            &[9.0, 4.0, 0.0],
        );
        // The third coordinate is an independent internal elastic displacement,
        // not this contact point's second tangent; its complete mode must remain.
        let zero = [0.0; 3];
        let c = SolidSustainedContact {
            velocity_weights: [&N, &T1, &zero],
            ..contact(h, 0.2, [0.0; 3])
        };
        for s in [&mut a, &mut reordered] {
            let r = s
                .advance_sustained_contacts(
                    &[c],
                    &[[10.0, -2.0, 0.0]],
                    &[0.0, -10.0, 0.0],
                    controls(),
                    cx,
                    |p| {
                        // x'=2*cos(2t)-sin(2t)>0 for all 0<=t<=h; y=0.
                        assert!(2.0 * (2.0 * h).cos() - (2.0 * h).sin() > 0.0);
                        near(p.total_force_n[1], 0.0, 1e-12);
                        Ok::<_, &'static str>(())
                    },
                )
                .unwrap();
            near(
                s.displacement_m()[0],
                (2.0 * h).sin() - 0.5 * (1.0 - (2.0 * h).cos()),
                1e-13,
            );
            near(
                s.velocity_mps()[0],
                2.0 * (2.0 * h).cos() - (2.0 * h).sin(),
                1e-13,
            );
            near(s.displacement_m()[2], 0.4 * (3.0 * h).cos(), 1e-13);
            near(s.velocity_mps()[2], -1.2 * (3.0 * h).sin(), 1e-13);
            near(
                r.frictional_dissipation_j,
                2.0 * s.displacement_m()[0],
                1e-12,
            );
            energy(&r);
        }
        reference.advance(&[-2.0, 0.0, 0.0], cx).unwrap();
        assert_eq!(bits(&a), bits(&reference));
        assert_eq!(bits(&a), bits(&reordered));
        for s in [&mut a, &mut reference, &mut reordered] {
            s.advance_affine(&[1.0, 0.0, -0.2], &[-0.5, 0.0, 0.3], cx)
                .unwrap();
        }
        assert_eq!(bits(&a), bits(&reference));
        assert_eq!(bits(&a), bits(&reordered));
    });
}

#[test]
fn g4_geometry_cancellation_budgets_and_impacts_preserve_exact_native_retry() {
    with_cx(|cx, gate| {
        let h = 0.125;
        let mut s = point(h, [2.0, 0.0, 0.0], cx);
        let before = bits(&s);
        for controls in [
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
                    &[contact(h, 0.5, [0.0; 3])],
                    &[[8.0, -4.0, 0.0]],
                    &[0.0, -8.0, 0.0],
                    controls,
                    cx,
                    |p| admit_rigid(p, &I)
                )
                .is_err()
            );
            assert_eq!(bits(&s), before);
        }
        let result = s.advance_sustained_contacts(
            &[contact(h, 0.5, [0.0; 3])],
            &[[8.0, -4.0, 0.0]],
            &[0.0, -8.0, 0.0],
            controls(),
            cx,
            |_| Err::<(), _>("interior geometry refused"),
        );
        assert!(matches!(
            result,
            Err(SolidSustainedError::Path("interior geometry refused"))
        ));
        assert_eq!(bits(&s), before);
        let result = s.advance_sustained_contacts(
            &[contact(h, 0.5, [0.0; 3])],
            &[[8.0, -4.0, 0.0]],
            &[0.0, -8.0, 0.0],
            controls(),
            cx,
            |p| {
                admit_rigid(p, &I)?;
                gate.request();
                Ok::<_, &'static str>(())
            },
        );
        assert!(matches!(result, Err(SolidSustainedError::Cancelled)));
        assert_eq!(bits(&s), before);
        with_cx(|fresh, _| {
            let r = s
                .advance_sustained_contacts(
                    &[contact(h, 0.5, [0.0; 3])],
                    &[[8.0, -4.0, 0.0]],
                    &[0.0, -8.0, 0.0],
                    controls(),
                    fresh,
                    |p| admit_rigid(p, &I),
                )
                .unwrap();
            let mut direct = point(h, [2.0, 0.0, 0.0], fresh);
            direct.advance(&[-4.0, 0.0, 0.0], fresh).unwrap();
            assert_eq!(bits(&s), bits(&direct));
            energy(&r);
            let mut impact = point(h, [0.0, -1.0, 0.0], fresh);
            let before = bits(&impact);
            assert!(matches!(
                impact.advance_sustained_contacts(
                    &[contact(h, 0.5, [0.0; 3])],
                    &[[0.0; 3]],
                    &[0.0, -8.0, 0.0],
                    controls(),
                    fresh,
                    |p| admit_rigid(p, &I)
                ),
                Err(SolidSustainedError::ImpactRequired { .. })
            ));
            assert_eq!(bits(&impact), before);
        });
    });
}

#[test]
fn g4_approximate_basis_cannot_bypass_the_original_work_gate() {
    with_cx(|cx, _| {
        let h = 0.125;
        let config = CompleteSpectralSolidControls {
            maximum_relative_operator_error: 1e-2,
            maximum_relative_state_error: 1e-2,
            maximum_work_error_j: 1e-14,
            ..native_controls()
        };
        let mut basis = I;
        basis[0] = 1.0001;
        let mut s = CompleteSpectralSolid::prepare_with_basis(
            &assembly(&[0.0; 9], &I, 3),
            &[0.0; 3],
            &basis,
            h,
            &[0.0; 3],
            &[2.0, 0.0, 0.0],
            config,
            cx,
        )
        .unwrap();
        let before = bits(&s);
        // Deliberately loose outer allowance isolates the ORIGINAL native gate.
        let result = s.advance_sustained_contacts(
            &[contact(h, 0.5, [0.0; 3])],
            &[[8.0, -4.0, 0.0]],
            &[0.0, -8.0, 0.0],
            SolidSustainedControls {
                maximum_work_error_j: 0.1,
                ..controls()
            },
            cx,
            |_| Ok::<_, &'static str>(()),
        );
        assert!(matches!(
            result,
            Err(SolidSustainedError::Native(
                CompleteSpectralSolidError::Energy { .. }
            ))
        ));
        assert_eq!(bits(&s), before);
    });
}
