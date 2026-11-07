//! Original rational coordinates, physical boundary outputs and signed slopes.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::phs_orbit::{
    PeriodicPortInput,
    rational::{
        PhsMimoError, PhsMimoLimits, PreparedPhsMimo,
        sideband::{PhsMimoSidebandError, PhsMimoSidebandSource},
    },
    sideband::{PhsPortSidebandControls, PhsPortSidebandError},
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_orbit::{
    HbOrbit,
    sideband::{
        SidebandControls, SidebandError, SidebandScales,
        sensitivity::SidebandCarrierDerivativeControls,
    },
};
use fs_phs::{PortHamiltonian, Storage};
use fs_vfit::{
    mimo::{
        MimoModel,
        load::{MimoLoadLimits, PreparedMimoLoad},
    },
    model::PoleTerm,
};
use std::{cell::Cell, rc::Rc};

const LIMITS: PhsMimoLimits = PhsMimoLimits {
    maximum_states: 32,
    maximum_ports: 8,
    maximum_harmonics: 31,
    maximum_hb_unknowns: 1024,
    maximum_products: 200_000_000,
    maximum_workspace_bytes: 4 << 20,
};
fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate, arena,
            StreamKey { seed: 419, kernel_id: 103, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        ))
    });
}
fn source(model: &MimoModel) -> PreparedMimoLoad {
    PreparedMimoLoad::new(
        model, 0.05,
        MimoLoadLimits {
            band_rad_s: [0.1, 100.0],
            sampled_passivity_points: 33,
            max_states: 32,
            max_memory_bytes: 1 << 20,
        },
    ).unwrap()
}
fn scalar(direct: f64, dark: bool) -> MimoModel {
    let mut model = MimoModel {
        ports: 1,
        poles: vec![PoleTerm::Real { pole: -2.0, residue: 0.0 }],
        residues: vec![vec![C64::from_re(2.0)]],
        d: vec![direct], e: vec![0.0], reciprocal: true,
    };
    if dark {
        model.poles.push(PoleTerm::Real { pole: -3.0, residue: 0.0 });
        model.residues.push(vec![C64::ZERO]);
    }
    model
}
fn controls(ports: usize, harmonics: usize) -> PhsPortSidebandControls {
    PhsPortSidebandControls {
        sideband: SidebandControls {
            harmonics,
            samples_per_period: 128,
            maximum_unknowns: 512,
            maximum_payload_bytes: 64 << 20,
            maximum_work: 1_000_000_000,
            maximum_base_residual: 1e-9,
            maximum_tangent_sampling_defect: 1e-9,
            maximum_relative_residual: 1e-10,
            maximum_omitted_residual: 1e-8,
        },
        output_harmonics: harmonics,
        initial_samples: 32,
        max_samples: 512,
        max_hessian_actions: 50_000,
        max_workspace_bytes: 16 << 20,
        max_observation_work: 1_000_000_000,
        envelope_tolerances: vec![1e-9; ports],
    }
}
fn derivative_controls() -> SidebandCarrierDerivativeControls {
    SidebandCarrierDerivativeControls {
        maximum_payload_bytes: 16 << 20,
        maximum_work: 500_000_000,
        maximum_relative_residual: 1e-10,
        maximum_omitted_residual: 1e-8,
    }
}
fn scales(dimension: usize) -> SidebandScales {
    SidebandScales { states: vec![1.0; dimension], equations: vec![1.0; dimension] }
}
fn zero_orbit(dimension: usize, omega: f64) -> HbOrbit {
    HbOrbit {
        omega, coeffs: vec![vec![C64::ZERO]; dimension],
        residual_trace: vec![0.0], residual: 0.0, aft_samples: 0, sampling_defect: 0.0,
    }
}
fn complex(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn original_dark_states_singular_d_signed_slopes_and_dc_survive_the_same_factor() {
    with_cx(|cx| {
        let (mass, stiffness, resistance, map, omega) = (2.0, 8.0, 0.3, 0.7, 2.0);
        let host = fs_phs::mass_spring_damper(mass, stiffness, resistance).unwrap();
        let model = scalar(0.0, true);
        let mut load = source(&model);
        load.step_mean_into(&[0.2], &mut [0.0]).unwrap();
        let history = load.pole_states().to_vec();
        let input = PeriodicPortInput::new(omega, vec![vec![C64::ZERO]]).unwrap();
        let mut loaded = PreparedPhsMimo::new(
            &host, &load, &input, &[0.0, map], LIMITS, cx,
        ).unwrap();
        let base = zero_orbit(4, omega);
        let mut source = PhsMimoSidebandSource::new(
            &mut loaded, &base, &controls(1, 1), &[1e-9], &[1e-9],
        ).unwrap();
        let mut tangent = source.prepare(&scales(4), || false).unwrap();
        let forcing = [C64::new(0.4, 0.1), C64::new(0.7, -0.2), C64::new(-0.2, 0.3)];
        let slope = [C64::new(0.05, 0.02), C64::new(-0.1, 0.03), C64::new(0.04, -0.01)];
        let i = C64::new(0.0, 1.0);
        for carrier in [-2.0, 0.0, 1.3, 2.0] {
            let mut factor = tangent.factor(carrier, || false).unwrap();
            let pair = factor.solve_with_carrier_derivative(
                &forcing, Some(&slope), derivative_controls(),
                &[1e-9], &[1e-9], &[1e-9], || false,
            ).unwrap();
            let r = &pair.response;
            let dr = &pair.angular_derivative;
            assert_eq!(r.response.state_dimension, 2);
            assert_eq!(r.source_dimension, 2);
            assert_eq!(r.response.state_response.coefficients.len(), 6);
            assert_eq!(r.source_coefficients.len(), 6);
            assert_eq!(r.source_band_rad_s, [0.1, 100.0]);
            assert_eq!(r.probe_includes_model_dc, carrier != 1.3);
            for band in 0..3 {
                let w = (band as f64 - 1.0).mul_add(omega, carrier);
                let s = C64::new(0.0, w);
                let pole = s + C64::from_re(2.0);
                let z = C64::from_re(2.0) / pole;
                let dz = i.scale(-2.0) / (pole * pole);
                let denominator = C64::from_re(stiffness) + s * s.scale(mass)
                    + s.scale(resistance) + s * z.scale(map * map);
                let d_denominator = i * (s.scale(2.0 * mass)
                    + C64::from_re(resistance) + z.scale(map * map))
                    + s * dz.scale(map * map);
                let q = forcing[band] / denominator;
                let dq = slope[band] / denominator - q * d_denominator / denominator;
                let velocity = s * q;
                let d_velocity = i * q + s * dq;
                complex(r.response.state_response.coefficients[2 * band], q, 2e-10);
                complex(dr.response.state_response.coefficients[2 * band], dq, 2e-10);
                complex(r.response.state_response.coefficients[2 * band + 1], velocity.scale(mass), 2e-10);
                complex(dr.response.state_response.coefficients[2 * band + 1], d_velocity.scale(mass), 2e-10);
                complex(r.response.port_coefficients[band], velocity, 2e-10);
                complex(dr.response.port_coefficients[band], d_velocity, 2e-10);
                let flow = velocity.scale(map);
                let d_flow = d_velocity.scale(map);
                complex(r.boundary_flow_coefficients[band], flow, 2e-10);
                complex(dr.boundary_flow_coefficients[band], d_flow, 2e-10);
                for (coordinate, pole) in [2.0, 3.0].into_iter().enumerate() {
                    let denominator = s + C64::from_re(pole);
                    let state = flow / denominator;
                    let derivative = d_flow / denominator - i * flow / (denominator * denominator);
                    complex(r.source_coefficients[2 * band + coordinate], state, 2e-10);
                    complex(dr.source_coefficients[2 * band + coordinate], derivative, 2e-10);
                }
                complex(r.boundary_effort_coefficients[band], r.source_coefficients[2 * band].scale(2.0), 2e-10);
                complex(dr.boundary_effort_coefficients[band], dr.source_coefficients[2 * band].scale(2.0), 2e-10);
            }
            assert!(r.response.state_response.relative_residual < 1e-10);
            assert!(dr.response.state_response.relative_residual < 1e-10);
            assert_eq!(pair.total_hessian_actions, 28 * 32);
            let again = factor.solve(&forcing, || false).unwrap();
            for (&a, &b) in again.source_coefficients.iter().zip(&r.source_coefficients) {
                complex(a, b, 1e-12);
            }
        }
        assert_eq!(load.pole_states(), history);
    });
}

#[test]
fn full_mutual_pair_realization_observes_boundary_ports_distinct_from_drive_ports() {
    with_cx(|cx| {
        let real = [0.8, -0.25, -0.25, 0.6];
        let pair = [0.4, 0.12, 0.12, 0.3];
        let direct = [0.9, 0.15, 0.15, 0.7];
        let model = MimoModel {
            ports: 2,
            poles: vec![
                PoleTerm::Real { pole: -2.0, residue: 0.0 },
                PoleTerm::Pair { pole: C64::new(-0.4, 1.8), residue: C64::ZERO },
            ],
            residues: vec![real.map(C64::from_re).to_vec(), pair.map(C64::from_re).to_vec()],
            d: direct.to_vec(), e: vec![0.0; 4], reciprocal: true,
        };
        let load = source(&model);
        let host = fs_phs::mass_spring_damper(2.0, 8.0, 0.3).unwrap();
        let input = PeriodicPortInput::new(2.0, vec![vec![C64::ZERO]]).unwrap();
        let map = [1.0, -0.4];
        let mut loaded = PreparedPhsMimo::new(
            &host, &load, &input, &[0.0, 0.0, map[0], map[1]], LIMITS, cx,
        ).unwrap();
        let base = zero_orbit(8, 2.0);
        let mut source = PhsMimoSidebandSource::new(
            &mut loaded, &base, &controls(1, 0), &[1e-9; 2], &[1e-9; 2],
        ).unwrap();
        let mut tangent = source.prepare(&scales(8), || false).unwrap();
        let i = C64::new(0.0, 1.0);
        let forcing = C64::new(0.7, -0.2);
        for carrier in [-1.3, 0.0, 1.3] {
            let mut factor = tangent.factor(carrier, || false).unwrap();
            let response = factor.solve_with_carrier_derivative(
                &[forcing], None, derivative_controls(), &[1e-9],
                &[1e-9; 2], &[1e-9; 2], || false,
            ).unwrap();
            let (r, dr) = (&response.response, &response.angular_derivative);
            assert_eq!(r.response.port_count, 1);
            assert_eq!(r.boundary_port_count, 2);
            assert_eq!(r.source_dimension, 6);
            let s = C64::new(0.0, carrier);
            let a = s + C64::from_re(0.4);
            let b = s + C64::from_re(2.0);
            let pair_den = a * a + C64::from_re(1.8 * 1.8);
            let mut z = [C64::ZERO; 4];
            let mut dz = z;
            for j in 0..4 {
                z[j] = C64::from_re(direct[j]) + C64::from_re(real[j]) / b
                    + a.scale(2.0 * pair[j]) / pair_den;
                dz[j] = i.scale(-real[j]) / (b * b)
                    + i * (C64::from_re(1.8 * 1.8) - a * a).scale(2.0 * pair[j])
                        / (pair_den * pair_den);
            }
            let (mut connection, mut d_connection) = (C64::ZERO, C64::ZERO);
            for row in 0..2 {
                for col in 0..2 {
                    connection = connection + z[2 * row + col].scale(map[row] * map[col]);
                    d_connection = d_connection + dz[2 * row + col].scale(map[row] * map[col]);
                }
            }
            let denominator = C64::from_re(8.0) + s * s.scale(2.0)
                + s.scale(0.3) + s * connection;
            let d_denominator = i * (s.scale(4.0) + C64::from_re(0.3) + connection)
                + s * d_connection;
            let q = forcing / denominator;
            let dq = -q * d_denominator / denominator;
            let velocity = s * q;
            let d_velocity = i * q + s * dq;
            complex(r.response.port_coefficients[0], velocity, 3e-10);
            complex(dr.response.port_coefficients[0], d_velocity, 3e-10);
            for port in 0..2 {
                let flow = velocity.scale(map[port]);
                let d_flow = d_velocity.scale(map[port]);
                complex(r.boundary_flow_coefficients[port], flow, 3e-10);
                complex(dr.boundary_flow_coefficients[port], d_flow, 3e-10);
                let states = [flow / b, flow * a / pair_den, flow.scale(-1.8) / pair_den];
                let slopes = [
                    d_flow / b - i * flow / (b * b),
                    d_flow * a / pair_den + i * flow
                        * (C64::from_re(1.8 * 1.8) - a * a) / (pair_den * pair_den),
                    d_flow.scale(-1.8) / pair_den
                        + i * flow * a.scale(3.6) / (pair_den * pair_den),
                ];
                for coordinate in 0..3 {
                    complex(r.source_coefficients[2 * coordinate + port], states[coordinate], 3e-10);
                    complex(dr.source_coefficients[2 * coordinate + port], slopes[coordinate], 3e-10);
                }
                let (mut effort, mut d_effort) = (C64::ZERO, C64::ZERO);
                for col in 0..2 {
                    effort = effort + z[2 * port + col] * velocity.scale(map[col]);
                    d_effort = d_effort + z[2 * port + col] * d_velocity.scale(map[col])
                        + dz[2 * port + col] * velocity.scale(map[col]);
                }
                complex(r.boundary_effort_coefficients[port], effort, 3e-10);
                complex(dr.boundary_effort_coefficients[port], d_effort, 3e-10);
            }
        }
    });
}

struct Quartic;
impl Storage for Quartic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.6 * x[0] * x[0] + 0.2 * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = 1.2 * x[0] + 0.8 * x[0].powi(3);
    }
}

#[test]
fn nonlinear_hessian_conversion_and_source_derivative_use_one_shared_observer() {
    with_cx(|cx| {
        let host = PortHamiltonian::new(
            1, 1, vec![0.0], vec![0.4], vec![1.0], Box::new(Quartic),
        ).unwrap();
        let model = scalar(0.2, false);
        let load = source(&model);
        let (omega, amplitude, map) = (2.0, 0.2, 0.5);
        let x1 = C64::from_re(amplitude / 2.0);
        let g1 = x1.scale(1.2 + 0.75 * 0.8 * amplitude * amplitude);
        let g3 = C64::from_re(0.8 * amplitude * amplitude * amplitude / 8.0);
        let z1 = g1.scale(map) / C64::new(2.0, omega);
        let z3 = g3.scale(map) / C64::new(2.0, 3.0 * omega);
        let u1 = C64::new(0.0, omega) * x1 + g1.scale(0.4 + 0.2 * map * map)
            + z1.scale(2.0 * map);
        let u3 = g3.scale(0.4 + 0.2 * map * map) + z3.scale(2.0 * map);
        let input = PeriodicPortInput::new(
            omega, vec![vec![C64::ZERO, u1, C64::ZERO, u3]],
        ).unwrap();
        let mut base = zero_orbit(2, omega);
        base.coeffs = vec![
            vec![C64::ZERO, x1, C64::ZERO, C64::ZERO],
            vec![C64::ZERO, z1, C64::ZERO, z3],
        ];
        let hessian_calls = Rc::new(Cell::new(0));
        let hessian = |x: &[f64], direction: &[f64], out: &mut [f64]| {
            hessian_calls.set(hessian_calls.get() + 1);
            out[0] = (1.2 + 2.4 * x[0] * x[0]) * direction[0];
            true
        };
        let mut loaded = PreparedPhsMimo::new(
            &host, &load, &input, &[map], LIMITS, cx,
        ).unwrap().with_storage_hessian(&hessian);
        let mut c = controls(1, 6);
        c.output_harmonics = 8;
        let mut source = PhsMimoSidebandSource::new(
            &mut loaded, &base, &c, &[1e-9], &[1e-9],
        ).unwrap();
        let mut tangent = source.prepare(&scales(2), || false).unwrap();
        assert!(tangent.base_residual() < 1e-10);
        let before = hessian_calls.get();
        let mut factor = tangent.factor(1.3, || false).unwrap();
        let mut forcing = vec![C64::ZERO; 13];
        forcing[6] = C64::new(0.3, -0.1);
        let pair = factor.solve_with_carrier_derivative(
            &forcing, None, derivative_controls(), &[1e-9], &[1e-9], &[1e-9], || false,
        ).unwrap();
        assert_eq!(hessian_calls.get() - before, pair.total_hessian_actions);
        assert_eq!(pair.total_hessian_actions, 28 * c.initial_samples);
        let h0 = 1.2 + 1.5 * 0.8 * amplitude * amplitude;
        let h2 = 0.75 * 0.8 * amplitude * amplitude;
        for response in [&pair.response, &pair.angular_derivative] {
            let at = |k: isize| {
                if k.unsigned_abs() <= 6 {
                    response.response.state_response.coefficients[(k + 6) as usize]
                } else {
                    C64::ZERO
                }
            };
            for k in -8_isize..=8 {
                let output = at(k).scale(h0) + (at(k - 2) + at(k + 2)).scale(h2);
                let band = (k + 8) as usize;
                let state = if k.unsigned_abs() <= 6 {
                    response.source_coefficients[(k + 6) as usize]
                } else {
                    C64::ZERO
                };
                complex(response.response.port_coefficients[band], output, 3e-10);
                complex(response.boundary_flow_coefficients[band], output.scale(map), 3e-10);
                complex(response.boundary_effort_coefficients[band], state.scale(2.0) + output.scale(0.2 * map), 3e-10);
            }
        }
        assert!(pair.response.response.state_response.coefficients[8].abs() > 1e-7);
        let i = C64::new(0.0, 1.0);
        for k in -6_isize..=6 {
            let (state_band, output_band) = ((k + 6) as usize, (k + 8) as usize);
            let s = C64::new(0.0, (k as f64).mul_add(omega, 1.3));
            let (r, dr) = (&pair.response, &pair.angular_derivative);
            let x = r.response.state_response.coefficients[state_band];
            let dx = dr.response.state_response.coefficients[state_band];
            let z = r.source_coefficients[state_band];
            let dz = dr.source_coefficients[state_band];
            complex((s + C64::from_re(2.0)) * z, r.boundary_flow_coefficients[output_band], 3e-10);
            complex((s + C64::from_re(2.0)) * dz + i * z, dr.boundary_flow_coefficients[output_band], 3e-10);
            complex(s * x, forcing[state_band] - r.response.port_coefficients[output_band].scale(0.4)
                - r.boundary_effort_coefficients[output_band].scale(map), 3e-10);
            complex(s * dx + i * x, -dr.response.port_coefficients[output_band].scale(0.4)
                - dr.boundary_effort_coefficients[output_band].scale(map), 3e-10);
        }
    });
}

#[test]
fn admission_rechecks_joint_base_and_keeps_original_band_cancellation_and_pair_caps() {
    with_cx(|cx| {
        let host = fs_phs::mass_spring_damper(2.0, 8.0, 0.3).unwrap();
        let load = source(&scalar(0.0, false));
        let input = PeriodicPortInput::new(2.0, vec![vec![C64::ZERO]]).unwrap();
        let mut loaded = PreparedPhsMimo::new(
            &host, &load, &input, &[0.0, 1.0], LIMITS, cx,
        ).unwrap();
        let base = zero_orbit(3, 2.0);
        let mut tiny = controls(1, 0);
        tiny.max_workspace_bytes = 1;
        assert!(matches!(
            PhsMimoSidebandSource::new(&mut loaded, &base, &tiny, &[1e-9], &[1e-9]),
            Err(PhsMimoSidebandError::Physical(PhsPortSidebandError::Budget(_)))
        ));
        let mut bad = base.clone();
        bad.coeffs[0][0] = C64::from_re(0.1);
        {
            let mut source = PhsMimoSidebandSource::new(
                &mut loaded, &bad, &controls(1, 0), &[1e-9], &[1e-9],
            ).unwrap();
            assert!(matches!(
                source.prepare(&scales(3), || false),
                Err(PhsMimoSidebandError::Physical(PhsPortSidebandError::Sideband(
                    SidebandError::BaseUnresolved { .. }
                )))
            ));
        }
        let mut c = controls(1, 0);
        c.max_hessian_actions = 14 * c.initial_samples;
        let mut source = PhsMimoSidebandSource::new(
            &mut loaded, &base, &c, &[1e-9], &[1e-9],
        ).unwrap();
        let mut tangent = source.prepare(&scales(3), || false).unwrap();
        for carrier in [0.05, 99.0] {
            assert!(matches!(
                tangent.factor(carrier, || false),
                Err(PhsMimoSidebandError::Causal(PhsMimoError::Source(_)))
            ));
        }
        assert!(matches!(
            tangent.factor(1.3, || true),
            Err(PhsMimoSidebandError::Physical(PhsPortSidebandError::Cancelled))
        ));
        let mut factor = tangent.factor(0.0, || false).unwrap();
        assert!(matches!(
            factor.solve_with_carrier_derivative(
                &[C64::from_re(1.0)], None, derivative_controls(),
                &[1e-9], &[1e-9], &[1e-9], || false,
            ),
            Err(PhsMimoSidebandError::Physical(PhsPortSidebandError::Budget(_)))
        ));
        assert!(matches!(
            factor.solve(&[C64::from_re(1.0)], || true),
            Err(PhsMimoSidebandError::Physical(PhsPortSidebandError::Cancelled))
        ));
        assert!(factor.solve(&[C64::from_re(1.0)], || false).unwrap().probe_includes_model_dc);
    });
}

#[test]
fn static_source_has_an_empty_source_chart_and_missing_hessian_is_not_approximated() {
    with_cx(|cx| {
        let static_model = MimoModel {
            ports: 1, poles: vec![], residues: vec![],
            d: vec![0.6], e: vec![0.0], reciprocal: true,
        };
        let load = source(&static_model);
        let input = PeriodicPortInput::new(2.0, vec![vec![C64::ZERO]]).unwrap();
        let host = PortHamiltonian::new(
            1, 1, vec![0.0], vec![1.0], vec![1.0],
            Box::new(fs_phs::QuadraticStorage::new(vec![2.0], 1).unwrap()),
        ).unwrap();
        let mut loaded = PreparedPhsMimo::new(&host, &load, &input, &[1.0], LIMITS, cx).unwrap();
        let base = zero_orbit(1, 2.0);
        {
            let mut source = PhsMimoSidebandSource::new(
                &mut loaded, &base, &controls(1, 0), &[1e-9], &[1e-9],
            ).unwrap();
            let mut tangent = source.prepare(&scales(1), || false).unwrap();
            let mut factor = tangent.factor(-1.3, || false).unwrap();
            let r = factor.solve_with_carrier_derivative(
                &[C64::from_re(1.0)], None, derivative_controls(),
                &[1e-9], &[1e-9], &[1e-9], || false,
            ).unwrap();
            let denominator = C64::new(3.2, -1.3);
            let x = C64::from_re(1.0) / denominator;
            let dx = C64::new(0.0, -1.0) / (denominator * denominator);
            assert_eq!(r.response.source_dimension, 0);
            assert!(r.response.source_coefficients.is_empty());
            assert!(r.angular_derivative.source_coefficients.is_empty());
            complex(r.response.response.state_response.coefficients[0], x, 2e-10);
            complex(r.angular_derivative.response.state_response.coefficients[0], dx, 2e-10);
            complex(r.response.boundary_effort_coefficients[0], x.scale(1.2), 2e-10);
            complex(r.angular_derivative.boundary_effort_coefficients[0], dx.scale(1.2), 2e-10);
        }
        let nonlinear = PortHamiltonian::new(
            1, 1, vec![0.0], vec![1.0], vec![1.0], Box::new(Quartic),
        ).unwrap();
        let mut loaded = PreparedPhsMimo::new(&nonlinear, &load, &input, &[1.0], LIMITS, cx).unwrap();
        let mut source = PhsMimoSidebandSource::new(
            &mut loaded, &base, &controls(1, 0), &[1e-9], &[1e-9],
        ).unwrap();
        assert!(matches!(
            source.prepare(&scales(1), || false),
            Err(PhsMimoSidebandError::Physical(PhsPortSidebandError::Sideband(
                SidebandError::TangentUnavailable
            )))
        ));
    });
}
