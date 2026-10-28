use super::*;
use fs_material::gas::GasSpec;
use fs_material::visco::GeneralizedMaxwell;
use fs_phs::WallPin;

const WORK: usize = 100_000_000;
fn duct() -> DgDuct {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.003,
                length: 0.08,
            }],
        },
        &gas,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::Matched,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()) + 1e-23,
        "{actual:e} != {expected:e}"
    );
}

/// A lossless source must export the very same skew generator as its existing
/// spectral propagator, including at maximum admitted polynomial order and in
/// variable-area geometry. A tiny indefinite symmetric remainder is not heat.
#[test]
fn continuous_lossless_generator_matches_the_original_skew_time_owner() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    for segment in [
        Segment::Cylinder {
            radius: 0.0073,
            length: 0.57,
        },
        Segment::Cone {
            inlet_radius: 0.004,
            outlet_radius: 0.020,
            length: 0.60,
        },
    ] {
        for order in [1, 3, MAX_DG_ORDER] {
            for outlet in [
                DgBoundary::Closed,
                DgBoundary::IdealOpen,
                DgBoundary::Port(PortId(9)),
            ] {
                let dg = DgDuct::new(
                    &Duct {
                        segments: vec![segment.clone()],
                    },
                    &gas,
                    DgMesh { elements: 2, order },
                    DgBoundary::Port(PortId(7)),
                    outlet,
                )
                .unwrap();
                let view = dg.continuous_operator();
                let l = view.generator_per_s();
                let n = view.dimension();
                let m = dg.field_dofs();
                let g = dg.g_block.as_ref().unwrap();
                for i in 0..n {
                    for j in 0..n {
                        assert_eq!(l[i * n + j] + l[j * n + i], 0.0, "order={order}, ({i},{j})");
                    }
                }
                for i in 0..m {
                    for j in 0..m {
                        assert_eq!(l[i * n + m + j], -g[i * m + j]);
                        assert_eq!(l[(m + i) * n + j], g[j * m + i]);
                    }
                }
                let propagator = dg.propagator(1.0 / 48000.0, 1e-12).unwrap();
                let PropagatorKernel::SkewModal {
                    sigmas,
                    u_cols,
                    v_cols,
                    ..
                } = &propagator.kernel
                else {
                    panic!("lossless time owner must retain its original skew kernel");
                };
                let scale = g.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
                let mut defect = 0.0_f64;
                for i in 0..m {
                    for j in 0..m {
                        let spectral: f64 = (0..m)
                            .map(|k| u_cols[i * m + k] * sigmas[k] * v_cols[j * m + k])
                            .sum();
                        defect = defect.max((spectral - g[i * m + j]).abs() / scale);
                    }
                }
                assert!(
                    defect < 2e-11,
                    "existing time generator defect={defect:e}, order={order}"
                );
            }
        }
    }
}

/// The exact transpose choice is restricted to independently verified SBP
/// roundoff; it may not turn an actually nonreciprocal source into a lossless one.
#[test]
fn continuous_lossless_source_refuses_a_nonroundoff_transpose_defect() {
    let g = [1.0, 2.0, 3.0, 4.0];
    let mut h = [-1.0, -3.0, -2.0, -4.0];
    DgDuct::admit_lossless_transpose_pair(&g, &h, 2, 2).unwrap();
    h[1] += 4.0 * f64::EPSILON;
    DgDuct::admit_lossless_transpose_pair(&g, &h, 2, 2).unwrap();
    h[1] += 1e-8;
    assert!(DgDuct::admit_lossless_transpose_pair(&g, &h, 2, 2).is_err());
    h[1] = f64::NAN;
    assert!(DgDuct::admit_lossless_transpose_pair(&g, &h, 2, 2).is_err());
}

#[test]
fn continuous_original_operator_retains_all_memory_and_independent_power() {
    let dg = duct()
        .with_viscothermal(16, [50.0, 2000.0], 0.01)
        .unwrap()
        .with_compliant_wall_relaxation(
            WallPin {
                surface_density: 0.2,
                stiffness_per_area: 8e7,
                resistance: 30.0,
            },
            &GeneralizedMaxwell::new(2e9, vec![(8e8, 0.002)]).unwrap(),
        )
        .unwrap();
    let view = dg.continuous_operator();
    assert_eq!(view.operator_identity(), dg.operator_identity());
    assert_eq!(view.dimension(), dg.total_dofs());
    assert!(view.dimension() > 2 * dg.field_dofs());
    let flow = C64::new(1e-7, 4e-8);
    let response = dg
        .harmonic_response(3000.0, &[(PortId(7), flow)], WORK, 1e-11, &mut || false)
        .unwrap();
    assert!(response.mean_wall_energy_j > 0.0);
    assert!(response.wall_dissipation_w > 0.0);
    assert!(response.matched_absorption_w > 0.0);
    assert!(response.input_power_w > 0.0);
    near(
        response.input_power_w,
        response.wall_dissipation_w + response.matched_absorption_w,
        2e-10,
    );
    assert_eq!(response.ports[0].inflow_m3_s, flow);
    // The port-column dot and original mass-map reconstruction may round
    // independently; this fixture has three nonzero projection coefficients.
    near(
        response.ports[0].pressure_pa.re,
        response.pressure_pa[0].re,
        8.0 * f64::EPSILON,
    );
    near(
        response.ports[0].pressure_pa.im,
        response.pressure_pa[0].im,
        8.0 * f64::EPSILON,
    );
    // Check the original generator against independent physical loss squares.
    let n = view.dimension();
    let mut generator_loss = 0.0;
    for i in 0..n {
        let rate = (0..n).fold(C64::ZERO, |sum, j| {
            sum + response.energy_coordinates[j].scale(view.generator_per_s()[i * n + j])
        });
        generator_loss -= 0.5 * (response.energy_coordinates[i].conj() * rate).re;
    }
    near(
        generator_loss,
        response.wall_dissipation_w + response.matched_absorption_w,
        2e-10,
    );
    let observed = view
        .observe_harmonic(
            3000.0,
            &response.energy_coordinates,
            &[(PortId(7), flow)],
            WORK,
            1e-11,
            &mut || false,
        )
        .unwrap();
    assert_eq!(observed, response);
    let mut wrong = response.energy_coordinates.clone();
    wrong[0] = wrong[0].scale(1.001);
    assert!(
        view.observe_harmonic(
            3000.0,
            &wrong,
            &[(PortId(7), flow)],
            WORK,
            1e-11,
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn harmonic_phasor_is_the_refined_limit_of_the_original_duct_time_owner() {
    let dg = duct();
    let omega = core::f64::consts::TAU * 1000.0;
    let force = 1e-7;
    let response = dg
        .harmonic_response(
            omega,
            &[(PortId(7), C64::from_re(force))],
            WORK,
            1e-11,
            &mut || false,
        )
        .unwrap();
    let mut errors = Vec::new();
    for rate in [16000, 32000, 64000] {
        let dt = 1.0 / f64::from(rate);
        let propagator = dg.propagator(dt, 1e-12).unwrap();
        let mut state = dg.zero_state();
        for (x, z) in state.tilde.iter_mut().zip(&response.energy_coordinates) {
            *x = z.re;
        }
        state.initial_energy_j = state.stored_energy_j();
        state.nodal_dirty = true;
        for tick in 0..rate / 250 {
            let flow = force * det::cos(omega * (f64::from(tick) + 0.5) * dt);
            propagator
                .step(
                    &mut state,
                    &[PortInput {
                        id: PortId(7),
                        effort: flow,
                    }],
                )
                .unwrap();
        }
        let error = det::sqrt(
            state
                .tilde
                .iter()
                .zip(&response.energy_coordinates)
                .map(|(x, z)| (x - z.re).powi(2))
                .sum::<f64>(),
        );
        let scale = det::sqrt(
            response
                .energy_coordinates
                .iter()
                .map(|z| z.norm_sq())
                .sum::<f64>(),
        );
        errors.push(error / scale);
    }
    assert!(
        errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "clock errors: {errors:?}"
    );
    assert!(errors[2] < 0.002, "clock errors: {errors:?}");
}

#[test]
fn harmonic_band_budget_cancellation_and_input_refusals_preserve_the_source() {
    let dg = duct().with_viscothermal(16, [50.0, 2000.0], 0.01).unwrap();
    let identity = dg.operator_identity();
    let columns = dg.continuous_operator().ports().to_vec();
    let flow = [(PortId(7), C64::from_re(1e-7))];
    for omega in [0.0, f64::NAN, 1.0, 20000.0] {
        assert!(
            dg.harmonic_response(omega, &flow, WORK, 1e-11, &mut || false)
                .is_err()
        );
    }
    assert!(
        dg.harmonic_response(1000.0, &flow, 0, 1e-11, &mut || false)
            .is_err()
    );
    assert!(
        dg.harmonic_response(1000.0, &flow, WORK, 0.0, &mut || false)
            .is_err()
    );
    assert!(
        dg.harmonic_response(1000.0, &[(PortId(9), C64::ONE)], WORK, 1e-11, &mut || false)
            .is_err()
    );
    assert!(
        dg.harmonic_response(
            1000.0,
            &[(PortId(7), C64::new(f64::NAN, 0.0))],
            WORK,
            1e-11,
            &mut || false
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        dg.harmonic_response(1000.0, &flow, WORK, 1e-11, &mut || {
            polls += 1;
            polls == 10
        }),
        Err(DuctError::Cancelled { completed_steps: 0 })
    ));
    let zero = dg
        .harmonic_response(1000.0, &[], WORK, 1e-11, &mut || false)
        .unwrap();
    assert_eq!(zero.mean_stored_energy_j, 0.0);
    assert_eq!(zero.input_power_w, 0.0);
    assert_eq!(zero.wall_dissipation_w, 0.0);
    assert_eq!(zero.matched_absorption_w, 0.0);
    assert_eq!(dg.operator_identity(), identity);
    assert_eq!(dg.continuous_operator().ports(), columns);
}
