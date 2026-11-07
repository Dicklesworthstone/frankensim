//! Actual tensor-Prony dispersion through the original structural harmonic factor.
use fs_exec::Cx;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicParameterControls, HarmonicPortControls, HarmonicPortParameterDirection,
    HarmonicResponseControls, HarmonicResponseError,
};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::visco::*;

fn with_cx(run: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 43,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 2e-12,
    }
}

fn tangent() -> HarmonicParameterControls {
    HarmonicParameterControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 1_000_000,
    }
}

fn diagonal(values: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}

fn scalar_law() -> TensorProny {
    TensorProny::new_mandel(
        diagonal([6.0; 6]),
        vec![PronyBranch::new(diagonal([18.0, 0.0, 0.0, 0.0, 0.0, 0.0]), 0.1).unwrap()],
        0.02,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(diagonal([1.2, 0.0, 0.0, 0.0, 0.0, 0.0]))
    .unwrap()
}

fn solid(law: TensorProny, free: &[usize], cx: &Cx<'_>) -> ViscoElements {
    let material = TetViscoMaterial {
        law,
        density_kg_m3: 60.0,
        material_state_identity: [43; 32],
    };
    let fixed: Vec<_> = (0..12).filter(|dof| !free.contains(dof)).collect();
    ViscoElements::new(
        SolidMeshRef::Tet4 {
            nodes_m: &[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            elements: &[[0, 1, 2, 3]],
        },
        TetViscoField::Uniform(&material),
        &fixed,
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 1000,
        },
        cx,
    )
    .unwrap()
}

fn near(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1e-15),
        "{actual:?} != {expected:?}"
    );
}

fn relative(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-11 * expected.abs().max(1e-15),
        "{actual:e} != {expected:e}"
    );
}

fn inner(a: C64, b: C64) -> f64 {
    (a.conj() * b).re
}

// Unit tetra: M_33=rho*V/10=1, K_eq=1, branch K=3, tau=.1, Kelvin C=.2.
fn scalar_pencil(w: f64) -> (C64, C64, f64, f64) {
    let x = 0.1 * w;
    let denominator = 1.0 + x * x;
    let storage = 1.0 + 3.0 * x * x / denominator;
    let storage_slope = 0.6 * x / denominator.powi(2);
    let loss = 0.2 * w + 3.0 * x / denominator;
    let loss_slope = 0.2 + 0.3 * (1.0 - x * x) / denominator.powi(2);
    (
        C64::new(storage - w * w, -loss),
        C64::new(storage_slope - 2.0 * w, -loss_slope),
        loss,
        loss_slope,
    )
}

#[test]
fn actual_material_dispersion_and_source_phase_give_analytic_motion_and_heat_slopes() {
    with_cx(|cx, _| {
        let source = solid(scalar_law(), &[3], cx);
        let df = C64::new(3e-5, -1e-5);
        for w in [0.4, 1.0, 1.1, 2.0, 10.0] {
            let f = C64::new(0.001, -0.0004) + df.scale(w);
            let (a, da, loss, dloss) = scalar_pencil(w);
            let q = f / a;
            let dq = (df - da * q) / a;
            let v = C64::new(0.0, -w) * q;
            let dv = C64::new(0.0, -1.0) * q + C64::new(0.0, -w) * dq;
            let expected_heat = 0.5 * w * loss * q.norm_sq();
            let expected_heat_slope =
                0.5 * (loss + w * dloss) * q.norm_sq() + w * loss * inner(q, dq);
            let prepared = source
                .prepare_harmonic_response(w, controls(), cx)
                .unwrap()
                .with_frequency_derivative(24, cx)
                .unwrap();
            assert_eq!(prepared.free_dofs(), &[3]);
            assert_eq!(prepared.omega_rad_s(), w);
            let result = prepared
                .solve(&[f], Some(&[df]), None, tangent(), cx)
                .unwrap();
            near(result.response.displacement_m[0], q);
            near(result.displacement_derivative[0], dq);
            near(result.velocity_derivative[0], dv);
            relative(result.response.average_material_heat_w, expected_heat);
            relative(result.average_material_heat_derivative, expected_heat_slope);
            relative(
                result.average_input_power_derivative,
                0.5 * (inner(df, v) + inner(f, dv)),
            );
            relative(result.average_input_power_derivative, expected_heat_slope);
            assert!(result.power_balance_residual_derivative.abs() < 1e-11 * expected_heat.abs());
            assert!(result.relative_derivative_residual < 2e-12);
            assert!(result.port_load.is_none());
            assert_eq!(
                prepared.base().solve(&[f], cx).unwrap().displacement_m,
                result.response.displacement_m
            );
            let replay = prepared
                .solve(&[f], Some(&[df]), None, tangent(), cx)
                .unwrap();
            assert_eq!(
                replay.displacement_derivative,
                result.displacement_derivative
            );
        }
    });
}

#[test]
fn anisotropic_cross_terms_and_consistent_mass_retain_the_complete_frequency_direction() {
    with_cx(|cx, _| {
        let mut eq = diagonal([6.0, 9.0, 12.0, 15.0, 18.0, 21.0]);
        eq[0][3] = 1.2;
        eq[3][0] = 1.2;
        let outer = |v: [f64; 6], scale: f64| -> Mandel6 {
            core::array::from_fn(|i| core::array::from_fn(|j| scale * v[i] * v[j]))
        };
        let law = TensorProny::new_mandel(
            eq,
            vec![PronyBranch::new(outer([1.0, 0.0, 0.0, 2.0, 0.0, 0.0], 12.0), 0.3).unwrap()],
            0.02,
        )
        .unwrap()
        .with_kelvin_viscosity_mandel(outer([1.0, 0.0, 0.0, -1.0, 0.0, 0.0], 0.6))
        .unwrap();
        let source = solid(law, &[3, 6], cx);
        let s = core::f64::consts::FRAC_1_SQRT_2;
        // Analytic V*B^T*C*B, including the Mandel sqrt(2) shear factors.
        let mass = [1.0, 0.5, 0.5, 1.0];
        let equilibrium = [1.0, 0.2 * s, 0.2 * s, 1.25];
        let branch = [2.0, 4.0 * s, 4.0 * s, 4.0];
        let kelvin = [0.1, -0.1 * s, -0.1 * s, 0.05];
        let q0 = [C64::new(4e-4, -2e-4), C64::new(-1e-4, 3e-4)];
        let dq = [C64::new(2e-5, -1e-5), C64::new(-3e-5, 1e-5)];
        for w in [0.4, 1.5, 4.0, 10.0] {
            let x = 0.3 * w;
            let den = 1.0 + x * x;
            let q: [C64; 2] = core::array::from_fn(|i| q0[i] + dq[i].scale(w));
            let loss: [f64; 4] = core::array::from_fn(|i| kelvin[i] * w + branch[i] * x / den);
            let dloss: [f64; 4] =
                core::array::from_fn(|i| kelvin[i] + branch[i] * 0.3 * (1.0 - x * x) / den.powi(2));
            let a: [C64; 4] = core::array::from_fn(|i| {
                C64::new(
                    equilibrium[i] + branch[i] * x * x / den - w * w * mass[i],
                    -loss[i],
                )
            });
            let da: [C64; 4] = core::array::from_fn(|i| {
                C64::new(
                    branch[i] * 0.6 * x / den.powi(2) - 2.0 * w * mass[i],
                    -dloss[i],
                )
            });
            let force: [C64; 2] = core::array::from_fn(|i| a[2 * i] * q[0] + a[2 * i + 1] * q[1]);
            let df: [C64; 2] = core::array::from_fn(|i| {
                da[2 * i] * q[0] + da[2 * i + 1] * q[1] + a[2 * i] * dq[0] + a[2 * i + 1] * dq[1]
            });
            let prepared = source
                .prepare_harmonic_response(w, controls(), cx)
                .unwrap()
                .with_frequency_derivative(96, cx)
                .unwrap();
            let result = prepared
                .solve(&force, Some(&df), None, tangent(), cx)
                .unwrap();
            let mut heat_slope = 0.0;
            for i in 0..2 {
                near(result.response.displacement_m[i], q[i]);
                near(result.displacement_derivative[i], dq[i]);
                for j in 0..2 {
                    let at = 2 * i + j;
                    heat_slope += 0.5 * (loss[at] + w * dloss[at]) * inner(q[i], q[j])
                        + 0.5 * w * loss[at] * (inner(dq[i], q[j]) + inner(q[i], dq[j]));
                }
            }
            relative(result.average_material_heat_derivative, heat_slope);
            relative(result.average_input_power_derivative, heat_slope);
        }
    });
}

#[test]
fn dispersive_boundary_and_attachment_slopes_keep_material_heat_separate() {
    with_cx(|cx, _| {
        let source = solid(scalar_law(), &[3], cx);
        let df = C64::new(3e-5, -1e-5);
        let force = |w| C64::new(0.001, -0.0004) + df.scale(w);
        let map = |w| 1.0 + 0.1 * w;
        let impedance = |w| C64::new(0.3, 0.2 / w - 0.05 * w);
        let prepare = |w| {
            source
                .prepare_harmonic_response_with_port_impedance(
                    w,
                    &[map(w)],
                    1,
                    &[impedance(w)],
                    controls(),
                    HarmonicPortControls {
                        maximum_ports: 1,
                        maximum_projection_terms: 10_000,
                        maximum_payload_bytes: 1 << 20,
                    },
                    cx,
                )
                .unwrap()
        };
        let w = 2.0;
        let prepared = prepare(w).with_frequency_derivative(24, cx).unwrap();
        let result = prepared
            .solve(
                &[force(w)],
                Some(&[df]),
                Some(HarmonicPortParameterDirection {
                    force_map: &[0.1],
                    impedance: &[C64::new(0.0, -0.2 / (w * w) - 0.05)],
                }),
                tangent(),
                cx,
            )
            .unwrap();
        let h = 1e-5;
        let plus = prepare(w + h).solve(&[force(w + h)], cx).unwrap();
        let minus = prepare(w - h).solve(&[force(w - h)], cx).unwrap();
        let difference = |a: C64, b: C64| (a - b).scale(0.5 / h);
        let near_difference = |actual: C64, expected: C64| {
            assert!((actual - expected).abs() < 3e-8 * expected.abs().max(1e-15))
        };
        near_difference(
            result.displacement_derivative[0],
            difference(plus.displacement_m[0], minus.displacement_m[0]),
        );
        near_difference(
            result.velocity_derivative[0],
            difference(plus.velocity_m_s[0], minus.velocity_m_s[0]),
        );
        let p = result.port_load.unwrap();
        let (a, b) = (plus.port_load.unwrap(), minus.port_load.unwrap());
        near_difference(p.flow_derivative[0], difference(a.flow[0], b.flow[0]));
        near_difference(p.effort_derivative[0], difference(a.effort[0], b.effort[0]));
        for (actual, expected) in [
            (
                result.average_material_heat_derivative,
                (plus.average_material_heat_w - minus.average_material_heat_w) / (2.0 * h),
            ),
            (
                p.average_absorbed_power_derivative,
                (a.average_absorbed_power_w - b.average_absorbed_power_w) / (2.0 * h),
            ),
            (
                result.average_input_power_derivative,
                (plus.average_input_power_w - minus.average_input_power_w) / (2.0 * h),
            ),
        ] {
            near_difference(C64::from_re(actual), C64::from_re(expected));
        }
        assert!(result.power_balance_residual_derivative.abs() < 1e-17);
    });
}

#[test]
fn derivative_budgets_cancellation_and_full_cycle_material_admission_are_preserved() {
    with_cx(|cx, gate| {
        let source = solid(scalar_law().with_band((0.5, 10.0)).unwrap(), &[3], cx);
        assert!(
            source
                .prepare_harmonic_response(11.0, controls(), cx)
                .is_err()
        );
        assert!(matches!(
            source
                .prepare_harmonic_response(2.0, controls(), cx)
                .unwrap()
                .with_frequency_derivative(23, cx),
            Err(ViscoHarmonicError::Harmonic(
                HarmonicResponseError::Budget { requested: 24, .. }
            ))
        ));
        let prepared = source
            .prepare_harmonic_response(2.0, controls(), cx)
            .unwrap()
            .with_frequency_derivative(24, cx)
            .unwrap();
        let q = C64::new(0.015, 0.015);
        let force = scalar_pencil(2.0).0 * q;
        assert!(matches!(
            prepared.solve(&[force], None, None, tangent(), cx),
            Err(ViscoHarmonicError::Element(ViscoElementError::Material(_)))
        ));
        let safe = force.scale(0.5);
        let result = prepared.solve(&[safe], None, None, tangent(), cx).unwrap();
        near(result.response.displacement_m[0], q.scale(0.5));
        assert!(matches!(
            prepared.solve(
                &[safe],
                None,
                None,
                HarmonicParameterControls {
                    maximum_work: 0,
                    ..tangent()
                },
                cx
            ),
            Err(ViscoHarmonicError::Harmonic(
                HarmonicResponseError::WorkBudget { .. }
            ))
        ));
        assert_eq!(
            prepared
                .solve(&[safe], None, None, tangent(), cx)
                .unwrap()
                .displacement_derivative,
            result.displacement_derivative
        );
        let pending = source
            .prepare_harmonic_response(0.5, controls(), cx)
            .unwrap();
        gate.request();
        assert!(pending.with_frequency_derivative(24, cx).is_err());
        assert!(matches!(
            prepared.solve(&[safe], None, None, tangent(), cx),
            Err(ViscoHarmonicError::Harmonic(
                HarmonicResponseError::Cancelled
            ))
        ));
    });
}
