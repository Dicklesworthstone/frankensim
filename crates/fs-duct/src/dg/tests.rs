use super::*;
use crate::profile::{RadiusProfile, RiccatiSteps};
use crate::{Duct, LossModel, Segment, Termination, input_impedance};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_math::det;
use fs_mms::{LadderSide, OrderGate, RefinementLadder};
use std::sync::Arc;
use std::time::Instant;

/// A scalar decay has independent closed forms for the homogeneous advance,
/// held-input response and integrated response. The latter two must retain
/// every Taylor factorial before the exact scaling/doubling recurrences.
#[test]
fn dg_dense_phi_functions_match_scalar_decay_and_held_input_integrals() {
    for z in [-0.25_f64, -0.5, -4.0, -32.0] {
        let exp = det::exp(z);
        let expected = [exp, (exp - 1.0) / z, (exp - 1.0 - z) / (z * z)];
        for extra_halving in [0, 1] {
            let shift = series_shift(&[z], 1, 1.0) + extra_halving;
            let (phi0, phi1, phi2, _) = series_phis(&[z], 1, 1.0, shift).unwrap();
            for (order, (actual, reference)) in [phi0[0], phi1[0], phi2[0]]
                .into_iter().zip(expected).enumerate()
            {
                assert!((actual - reference).abs() <= 2e-13 * reference.abs() + 1e-30,
                    "phi_{order}({z}), shift={shift}: {actual:e} != {reference:e}");
            }
        }
    }
}

/// The same three functions on a damped rotation exercise oscillatory forcing
/// and off-diagonal coupling against the independent scalar complex solution.
#[test]
fn dg_dense_phi_functions_match_damped_rotation() {
    for z in [C64::new(-0.25, 0.5), C64::new(-4.0, 3.0), C64::new(-32.0, 12.0)] {
        let exp = C64::new(det::cos(z.im), det::sin(z.im)).scale(det::exp(z.re));
        let expected = [exp, (exp - C64::ONE) / z, (exp - C64::ONE - z) / (z * z)];
        let matrix = [z.re, -z.im, z.im, z.re];
        for extra_halving in [0, 1] {
            let shift = series_shift(&matrix, 2, 1.0) + extra_halving;
            let (phi0, phi1, phi2, _) = series_phis(&matrix, 2, 1.0, shift).unwrap();
            for (order, (actual, reference)) in [phi0, phi1, phi2]
                .into_iter().zip(expected).enumerate()
            {
                for (value, exact) in actual.into_iter()
                    .zip([reference.re, -reference.im, reference.im, reference.re])
                {
                    assert!((value - exact).abs() <= 2e-13 * exact.abs() + 1e-30,
                        "rotation phi_{order}({z:?}), shift={shift}: {value:e} != {exact:e}");
                }
            }
        }
    }
}

fn dry_air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("dry air at 293.15 K")
}

fn cylinder_fixture() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0073,
            length: 0.57,
        }],
    }
}

fn cone_fixture() -> Duct {
    Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.004,
            outlet_radius: 0.020,
            length: 0.60,
        }],
    }
}

/// G5: a uniform source retains the original assembly and propagation bits,
/// including the dissipative characteristic boundary and a supplied profile.
#[test]
fn dg_gas_field_uniform_is_bitwise_single_state() {
    let field = crate::profile::gas_field::tests::field([300.0; 2], [0.8; 2], 101325.0);
    let gas = *field.uniform_state().unwrap().state();
    let geometries = [
        cylinder_fixture(),
        Duct {
            segments: vec![Segment::Profile {
                profile: Arc::new(RadiusProfile::exponential(0.0073, 0.8, 0.57).unwrap()),
                steps: RiccatiSteps(8),
            }],
        },
    ];
    let mesh = DgMesh {
        elements: 4,
        order: 3,
    };
    for geometry in geometries {
        for outlet in [DgBoundary::IdealOpen, DgBoundary::Matched] {
            let original = DgDuct::new(&geometry, &gas, mesh, DgBoundary::Closed, outlet).unwrap();
            let bound =
                DgDuct::new_with_gas_field(&geometry, &field, mesh, DgBoundary::Closed, outlet)
                    .unwrap();
            assert_eq!(original.l_tilde, bound.l_tilde);
            assert_eq!(original.g_block, bound.g_block);
            assert_eq!(original.node_positions(), bound.node_positions());
            assert_eq!(bound.gas_field(), Some(&field));
            assert_eq!(bound.inlet_gas(), &gas);
            assert_eq!(bound.outlet_gas(), &gas);
            let initial = |duct: &DgDuct| {
                duct.project_initial_state(|x| det::sin(core::f64::consts::PI * x / 0.57), |_| 0.0)
            };
            let mut before = initial(&original);
            let mut after = initial(&bound);
            let a = original.propagator(1.0 / 48000.0, 1e-13).unwrap();
            let b = bound.propagator(1.0 / 48000.0, 1e-13).unwrap();
            for _ in 0..64 {
                assert_eq!(
                    a.step(&mut before, &[]).unwrap(),
                    b.step(&mut after, &[]).unwrap()
                );
            }
            a.sync_nodal(&mut before).unwrap();
            b.sync_nodal(&mut after).unwrap();
            assert_eq!(before.pressure(), after.pressure());
            assert_eq!(before.volume_flow(), after.volume_flow());
        }
    }
}

/// G1/G3: for constant composition and linear temperature, K=gamma*p is
/// constant while rho is proportional to 1/T. Integrating rho's exact log
/// tests material reevaluation, and would reject interpolation of rho itself.
#[test]
fn dg_gas_field_masses_match_analytic_thermodynamics() {
    let field = crate::profile::gas_field::tests::field([300.0, 600.0], [0.8; 2], 101325.0);
    let dg = DgDuct::new_with_gas_field(
        &cylinder_fixture(),
        &field,
        DgMesh {
            elements: 3,
            order: 3,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    let ones = vec![1.0; dg.field_dofs()];
    let zeros = vec![0.0; dg.field_dofs()];
    let norm_squared = |v: Vec<f64>| v.iter().map(|x| x * x).sum::<f64>();
    let mass_p = norm_squared(dg.to_energy_coords(&ones, &zeros));
    let mass_q = norm_squared(dg.to_energy_coords(&zeros, &ones));
    let inlet = field.state_at(0.0).unwrap();
    let area = core::f64::consts::PI * 0.0073_f64.powi(2);
    let bulk = inlet.state().density * inlet.state().sound_speed.powi(2);
    let expected_p = area * 0.57 / bulk;
    let expected_q = inlet.state().density * 0.57 * det::ln(2.0) / area;
    assert!((mass_p / expected_p - 1.0).abs() < 1e-12);
    assert!((mass_q / expected_q - 1.0).abs() < 1e-12);
    let interpolated_density_mass = 0.75 * inlet.state().density * 0.57 / area;
    assert!((mass_q / interpolated_density_mass - 1.0).abs() > 0.05);
}

/// G0: the physical boundary power uses local endpoint impedances, even when
/// their gases differ. The quadratic form is an independent reconstruction of
/// the matched-boundary dissipation before any time-step ledger is consulted.
#[test]
fn dg_gas_field_matched_boundaries_use_local_impedance() {
    let field = crate::profile::gas_field::tests::field([300.0, 600.0], [0.95, 0.1], 101325.0);
    let dg = DgDuct::new_with_gas_field(
        &cylinder_fixture(),
        &field,
        DgMesh {
            elements: 4,
            order: 3,
        },
        DgBoundary::Matched,
        DgBoundary::Matched,
    )
    .unwrap();
    let area = core::f64::consts::PI * 0.0073_f64.powi(2);
    let zi = dg.inlet_gas().density * dg.inlet_gas().sound_speed / area;
    let zo = dg.outlet_gas().density * dg.outlet_gas().sound_speed / area;
    assert!((zi / zo - 1.0).abs() > 0.1);
    let positions = dg.node_positions();
    let p: Vec<_> = positions.iter().map(|x| 1.0 + x / 0.57).collect();
    let q: Vec<_> = positions
        .iter()
        .map(|x| (3.0 + 4.0 * x / 0.57) * 1e-7)
        .collect();
    let state = dg.to_energy_coords(&p, &q);
    let mut rate = 0.0;
    for (i, &xi) in state.iter().enumerate() {
        for (j, &xj) in state.iter().enumerate() {
            rate += xi * dg.l_tilde[i * state.len() + j] * xj;
        }
    }
    let expected = -0.5 * (1.0 / zi + zi * (3e-7_f64).powi(2) + 4.0 / zo + zo * (7e-7_f64).powi(2));
    assert!(
        (rate / expected - 1.0).abs() < 1e-11,
        "rate={rate:e}, expected={expected:e}"
    );
}

/// G0/G3: material gradients change physical storage, while the original
/// reciprocal time operator still conserves energy and accounts for port work.
#[test]
fn dg_gas_field_energy_and_port_work_close() {
    let field = crate::profile::gas_field::tests::field([300.0, 450.0], [0.9, 0.2], 101325.0);
    let mesh = DgMesh {
        elements: 5,
        order: 3,
    };
    for inlet in [DgBoundary::Closed, DgBoundary::Port(PortId(0))] {
        let dg = DgDuct::new_with_gas_field(
            &cylinder_fixture(),
            &field,
            mesh,
            inlet,
            DgBoundary::IdealOpen,
        )
        .unwrap();
        let propagator = dg.propagator(1.0 / 48000.0, 1e-13).unwrap();
        let mut state =
            dg.project_initial_state(|x| det::cos(core::f64::consts::PI * x / 1.14), |_| 0.0);
        let initial = state.initial_energy_j();
        let mut work = 0.0;
        for step in 0..1000 {
            let inputs = if matches!(inlet, DgBoundary::Port(_)) {
                vec![PortInput {
                    id: PortId(0),
                    effort: 1e-6 * det::sin(0.01 * f64::from(step)),
                }]
            } else {
                vec![]
            };
            let output = propagator.step(&mut state, &inputs).unwrap();
            assert!(output.ledger.balance_checked);
            for (input, output) in inputs.iter().zip(&output.port_outputs) {
                work += input.effort * output.delta_x;
            }
        }
        propagator.sync_nodal(&mut state).unwrap();
        let physical = dg.to_energy_coords(state.pressure(), state.volume_flow());
        let energy = 0.5 * physical.iter().map(|x| x * x).sum::<f64>();
        assert!((energy - initial - work).abs() / (initial + work.abs()) < 1e-12);
    }
}

/// G1: the first pole converges to the independent Frobenius solution of the
/// linear-temperature Webster equation; the gas gradient remains audible and
/// cannot be replaced by the inlet state at every quadrature point.
#[test]
fn dg_gas_field_resonance_converges_to_analytic_gradient() {
    use crate::profile::gas_field::tests::{field, linear_temperature_exact};
    let source = field([300.0, 330.0], [0.8; 2], 101325.0);
    let gas = *source.state_at(0.0).unwrap().state();
    let uniform_omega = core::f64::consts::PI * gas.sound_speed / (2.0 * 0.57);
    let (mut low, mut high) = (uniform_omega * 0.9, uniform_omega * 1.2);
    for _ in 0..48 {
        let middle = 0.5 * (low + high);
        if linear_temperature_exact(&gas, 30.0, middle)[3].re > 0.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    let exact = 0.5 * (low + high);
    assert!((exact / uniform_omega - 1.0).abs() > 0.005);
    let mut errors = Vec::new();
    for elements in [2, 4, 8] {
        let dg = DgDuct::new_with_gas_field(
            &cylinder_fixture(),
            &source,
            DgMesh { elements, order: 3 },
            DgBoundary::Closed,
            DgBoundary::IdealOpen,
        )
        .unwrap();
        let propagator = dg.propagator(1.0 / 48000.0, 1e-13).unwrap();
        let PropagatorKernel::SkewModal { sigmas, .. } = &propagator.kernel else {
            panic!("lossless boundaries require the existing skew propagator");
        };
        let first = sigmas
            .iter()
            .copied()
            .filter(|x| *x > 0.0)
            .min_by(f64::total_cmp)
            .unwrap();
        errors.push((first / exact - 1.0).abs());
    }
    assert!(errors[0] > 20.0 * errors[1], "{errors:?}");
    assert!(errors[1] > 20.0 * errors[2], "{errors:?}");
    assert!(errors[2] < 1e-8, "{errors:?}");
}

/// G0/G3: source knots remain actual element boundaries, even inside a cone.
/// A budget cannot silently erase a knot or resize the supplied field.
#[test]
fn dg_gas_field_coverage_and_material_stations() {
    use fs_material::gas_species::input::ThermochemicalPhase;
    let base = crate::profile::gas_field::tests::field([300.0, 450.0], [0.8; 2], 101325.0);
    let middle = base.primitive_states()[0].clone();
    let field = AxialGasField::new(
        vec![0.0, 0.13, 0.57],
        101325.0,
        vec![
            base.primitive_states()[0].clone(),
            middle,
            base.primitive_states()[1].clone(),
        ],
        base.model().clone(),
        fs_matdb::QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap();
    let cone = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.005,
            outlet_radius: 0.015,
            length: 0.57,
        }],
    };
    assert!(
        DgDuct::new_with_gas_field(
            &cone,
            &field,
            DgMesh {
                elements: 1,
                order: 3
            },
            DgBoundary::Closed,
            DgBoundary::Closed
        )
        .is_err()
    );
    assert!(
        DgDuct::new_with_gas_field(
            &cone_fixture(),
            &field,
            DgMesh {
                elements: 4,
                order: 3
            },
            DgBoundary::Closed,
            DgBoundary::Closed
        )
        .is_err()
    );
    let dg = DgDuct::new_with_gas_field(
        &cone,
        &field,
        DgMesh {
            elements: 4,
            order: 3,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    assert!(dg.elements.iter().any(|element| element.x_left == 0.13));
    for element in &dg.elements {
        for (&x, &area) in element.nodes_x.iter().zip(&element.areas) {
            let radius = 0.005 + 0.01 * x / 0.57;
            assert!((area / (core::f64::consts::PI * radius * radius) - 1.0).abs() < 1e-14);
        }
    }
}

/// G0: LGL weights integrate polynomials up to degree 2p - 1 exactly (1e-14),
/// endpoints are included, and discrete integration by parts holds to 1e-14.
#[test]
fn dg_lgl_nodes_exact() {
    let start = Instant::now();
    for p in 1..=MAX_DG_ORDER {
        let rule = lgl_rule(p).expect("admitted LGL order");
        assert_eq!(rule.nodes.len(), p + 1);
        assert_eq!(rule.weights.len(), p + 1);
        assert!((rule.nodes[0] - (-1.0)).abs() <= 1e-15);
        assert!((rule.nodes[p] - 1.0).abs() <= 1e-15);
        for w in rule.nodes.windows(2) {
            assert!(w[1] > w[0], "LGL nodes must be strictly ascending");
        }
        for &w in &rule.weights {
            assert!(w > 0.0, "LGL weights must be strictly positive");
        }

        let max_deg = 2 * p - 1;
        let mut max_quad_err = 0.0f64;
        for deg in 0..=max_deg {
            let exact = if deg % 2 == 0 {
                2.0 / ((deg + 1) as f64)
            } else {
                0.0
            };
            let approx: f64 = rule
                .nodes
                .iter()
                .zip(rule.weights.iter())
                .map(|(&xi, &w)| w * det::powi(xi, deg as i32))
                .sum();
            let err = (approx - exact).abs();
            max_quad_err = max_quad_err.max(err);
            assert!(
                err <= 1e-14,
                "p={p}, deg={deg}: LGL quadrature error {err:e} > 1e-14"
            );
        }

        // Discrete integration by parts: S + S^T = diag(-1, 0, ..., 0, +1) where S_{ij} = w_i D_{ij}.
        let n = p + 1;
        let mut max_ibp_err = 0.0f64;
        for i in 0..n {
            for j in 0..n {
                let s_ij = rule.weights[i] * rule.diff_matrix[i * n + j];
                let s_ji = rule.weights[j] * rule.diff_matrix[j * n + i];
                let expected = if i == 0 && j == 0 {
                    -1.0
                } else if i == p && j == p {
                    1.0
                } else {
                    0.0
                };
                max_ibp_err = max_ibp_err.max((s_ij + s_ji - expected).abs());
            }
        }
        assert!(
            max_ibp_err <= 1e-13,
            "p={p}: discrete IBP error {max_ibp_err:e} > 1e-13"
        );

        eprintln!(
            "{{\"suite\":\"fs-duct\",\"case\":\"dg_lgl_nodes_exact\",\"inputs\":\"analytic_legendre_lobatto\",\"axis\":\"order\",\"rung\":{p},\"order\":{p},\"observed\":{max_quad_err:.3e},\"ibp_err\":{max_ibp_err:.3e},\"oracle\":\"2/(deg+1)\",\"tol\":1e-14,\"verdict\":\"pass\",\"elapsed_us\":{}}}",
            start.elapsed().as_micros()
        );
    }
}

/// G0: the element mass matrices integrate the exact area. With `Σ_i ℓ_i = 1`, a unit field
/// maps through the energy coordinates to `‖M^{1/2} 1‖² = 1ᵀ M 1`, which is `∫ A/(ρc²) dx` for
/// `p` and `∫ ρ/A dx` for `q`, in closed form for a cone `r = r_0 + α x`.
#[test]
fn dg_mass_matrices_integrate_the_exact_area() {
    let start = Instant::now();
    let gas = dry_air();
    let dg = DgDuct::new(
        &cone_fixture(),
        &gas,
        DgMesh {
            elements: 8,
            order: 4,
        },
        DgBoundary::Closed,
        DgBoundary::IdealOpen,
    )
    .expect("cone DG");
    let (r0, alpha) = (0.004, 0.016 / 0.60);
    let (rho, bulk) = (gas.density, gas.density * gas.sound_speed * gas.sound_speed);
    let pi = core::f64::consts::PI;
    let m = dg.field_dofs;
    let (ones, zeros) = (vec![1.0; m], vec![0.0; m]);
    let unit_p = dg.to_energy_coords(&ones, &zeros);
    let unit_q = dg.to_energy_coords(&zeros, &ones);
    assert!(
        unit_p[m..].iter().chain(&unit_q[..m]).all(|&t| t == 0.0),
        "fields do not mix in the energy map"
    );
    let n_local = 5;
    let mut worst = 0.0_f64;
    for (e, elem) in dg.elements.iter().enumerate() {
        let (ra, rb) = (
            r0 + alpha * elem.x_left,
            r0 + alpha * (elem.x_left + elem.length),
        );
        let exact_q = rho / (pi * alpha) * (1.0 / ra - 1.0 / rb);
        let exact_p = pi / bulk * (rb * rb * rb - ra * ra * ra) / (3.0 * alpha);
        let slice = e * n_local..(e + 1) * n_local;
        let sum_p: f64 = unit_p[slice.clone()].iter().map(|t| t * t).sum();
        let sum_q: f64 = unit_q[m + slice.start..m + slice.end]
            .iter()
            .map(|t| t * t)
            .sum();
        worst = worst
            .max((sum_q - exact_q).abs() / exact_q)
            .max((sum_p - exact_p).abs() / exact_p);
    }
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dg_mass_exact_area\",\"inputs\":\"cone_4to20mm_0.6m_closed_form_integrals\",\"elements\":8,\"order\":4,\"observed\":{worst:.3e},\"oracle\":\"closed form\",\"tol\":1e-13,\"verdict\":\"{}\",\"elapsed_us\":{}}}",
        if worst <= 1e-13 { "pass" } else { "fail" },
        start.elapsed().as_micros()
    );
    assert!(
        worst <= 1e-13,
        "mass matrices vs exact area integrals {worst:e}"
    );
}

/// The extension of `f` (given on `[0, L]`) across both ends, `F(-y) = s0 F(y)` and
/// `F(2L - y) = sl F(y)`: even at a rigid end, odd at a pressure-release end.
fn image(f: &dyn Fn(f64) -> f64, mut y: f64, length: f64, s0: f64, sl: f64) -> f64 {
    let mut sign = 1.0;
    loop {
        if y < 0.0 {
            y = -y;
            sign *= s0;
        } else if y > length {
            y = 2.0 * length - y;
            sign *= sl;
        } else {
            return sign * f(y);
        }
    }
}

/// Exact steps from a projected pulse at rest to `t_final`; returns the L2 pressure error, the
/// DOF count and the propagator certificate.
#[allow(clippy::too_many_arguments)]
fn ir_error(
    duct: &Duct,
    mesh: DgMesh,
    inlet: DgBoundary,
    outlet: DgBoundary,
    p0: &dyn Fn(f64) -> f64,
    p_exact: &dyn Fn(f64) -> f64,
    t_final: f64,
    steps: usize,
) -> (f64, usize, f64) {
    let dg = DgDuct::new(duct, &dry_air(), mesh, inlet, outlet).expect("DG duct");
    let prop = dg
        .propagator(t_final / steps as f64, 1e-13)
        .expect("certified propagator");
    let mut state = dg.project_initial_state(p0, |_| 0.0);
    for _ in 0..steps {
        prop.step(&mut state, &[]).expect("step");
    }
    prop.sync_nodal(&mut state).unwrap();
    (
        dg.l2_pressure_error(&state, p_exact).unwrap(),
        dg.total_dofs(),
        prop.certificate_residual(),
    )
}

/// h-ladders gated at order p + 1 (fs-mms, ±0.2) for p in {2, 3, 4} with rung-1 liveness, then
/// the p-ladder at K = 16, which must fall by at least 2x per order. The h-ladder starts at
/// K = 24: at K = 16 (element = σ) p = 2 measured a local order of 4.0 (pre-asymptotic), so
/// rungs were added at the fine end rather than the gate widened.
#[allow(clippy::too_many_arguments)]
fn h_and_p_ladders(
    case: &str,
    inputs: &str,
    duct: &Duct,
    inlet: DgBoundary,
    outlet: DgBoundary,
    p0: &dyn Fn(f64) -> f64,
    p_exact: &dyn Fn(f64) -> f64,
    t_final: f64,
    length: f64,
) {
    const H_RUNGS: [usize; 5] = [24, 32, 48, 64, 96];
    const IR_STEPS: usize = 6;
    for p in [2usize, 3, 4] {
        let (mut hs, mut errs) = (Vec::new(), Vec::new());
        for &k in &H_RUNGS {
            let t0 = Instant::now();
            let mesh = DgMesh {
                elements: k,
                order: p,
            };
            let (err, dofs, cert) =
                ir_error(duct, mesh, inlet, outlet, p0, p_exact, t_final, IR_STEPS);
            hs.push(length / k as f64);
            errs.push(err);
            eprintln!(
                "{{\"suite\":\"fs-duct\",\"case\":\"{case}\",\"inputs\":\"{inputs}\",\"axis\":\"elements\",\"rung\":{k},\"elements\":{k},\"order\":{p},\"dofs\":{dofs},\"dt\":{:.4e},\"l2_err\":{err:.6e},\"certificate\":{cert:.2e},\"elapsed_us\":{}}}",
                t_final / IR_STEPS as f64,
                t0.elapsed().as_micros()
            );
        }
        assert!(
            errs[0] > 1e-6,
            "{case} p={p}: rung-1 error {} must exceed 1e-6 for ladder liveness",
            errs[0]
        );
        let ladder = RefinementLadder::new(hs, errs).expect("valid refinement ladder");
        let verdict = OrderGate {
            theoretical: (p + 1) as f64,
        }
        .check(&format!("{case}-p{p}"), LadderSide::Primal, &ladder)
        .unwrap_or_else(|e| panic!("{case}: h-convergence gate failed for p={p}: {e}"));
        eprintln!("{}", verdict.json_line(true));
    }
    let p_errs: Vec<f64> = (1..=6)
        .map(|p| {
            let mesh = DgMesh {
                elements: 16,
                order: p,
            };
            ir_error(duct, mesh, inlet, outlet, p0, p_exact, t_final, IR_STEPS).0
        })
        .collect();
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"{case}\",\"inputs\":\"{inputs}\",\"axis\":\"order\",\"elements\":16,\"rungs\":[1,2,3,4,5,6],\"l2_err\":{p_errs:?}}}"
    );
    for w in p_errs.windows(2) {
        assert!(
            w[1] < 0.5 * w[0],
            "{case}: p-convergence must at least halve the error per order: {p_errs:?}"
        );
    }
}

/// G1/G2: a Gaussian pressure pulse at rest (centre 0.4 L, σ = L/16) in the closed/open
/// 0.57 m x 7.3 mm cylinder, after 0.75 L/c: each half has reflected once, evenly at the
/// closed inlet and oddly at the open outlet. Exact: `(F(x - ct) + F(x + ct))/2`.
#[test]
fn dg_cylinder_ir_matches_dalembert() {
    let gas = dry_air();
    let length = 0.57;
    let c = gas.sound_speed;
    let (centre, sigma) = (0.4 * length, length / 16.0);
    let p0 = move |x: f64| 100.0 * det::exp(-0.5 * ((x - centre) / sigma).powi(2));
    let t_final = 0.75 * length / c;
    let p_exact = |x: f64| {
        0.5 * (image(&p0, x - c * t_final, length, 1.0, -1.0)
            + image(&p0, x + c * t_final, length, 1.0, -1.0))
    };
    h_and_p_ladders(
        "dg_cylinder_ir_matches_dalembert",
        "cylinder_0.57m_r7.3mm_closed_open_gaussian_analytic_dalembert_images",
        &cylinder_fixture(),
        DgBoundary::Closed,
        DgBoundary::IdealOpen,
        &p0,
        &p_exact,
        t_final,
        length,
    );
}

/// G1/G2: in the open/open conical frustum (4 -> 20 mm over 0.6 m), `ψ = (x_apex + x) p` obeys
/// the 1-D wave equation exactly, with odd images at both pressure-release ends. A Gaussian ψ
/// at rest (centre 0.4 L, σ = L/16) after 0.75 L/c. `Segment::Profile` (a linear power profile)
/// reproduces `Segment::Cone` to roundoff.
#[test]
fn dg_cone_ir_matches_spherical_exact() {
    let gas = dry_air();
    let (r_in, r_out, length) = (0.004, 0.020, 0.60);
    let x_apex = r_in * length / (r_out - r_in);
    let c = gas.sound_speed;
    let (centre, sigma) = (0.4 * length, length / 16.0);
    let psi0 = move |x: f64| 25.0 * det::exp(-0.5 * ((x - centre) / sigma).powi(2));
    let p0 = move |x: f64| psi0(x) / (x_apex + x);
    let t_final = 0.75 * length / c;
    let p_exact = |x: f64| {
        0.5 * (image(&psi0, x - c * t_final, length, -1.0, -1.0)
            + image(&psi0, x + c * t_final, length, -1.0, -1.0))
            / (x_apex + x)
    };
    h_and_p_ladders(
        "dg_cone_ir_matches_spherical_exact",
        "cone_4to20mm_0.6m_open_open_gaussian_analytic_spherical_rp_images",
        &cone_fixture(),
        DgBoundary::IdealOpen,
        DgBoundary::IdealOpen,
        &p0,
        &p_exact,
        t_final,
        length,
    );

    let mesh = DgMesh {
        elements: 8,
        order: 4,
    };
    let profile_duct = Duct {
        segments: vec![Segment::Profile {
            profile: Arc::new(
                RadiusProfile::power(r_in, x_apex, 1.0, length).expect("linear power profile"),
            ),
            steps: RiccatiSteps(64),
        }],
    };
    let run = |duct: &Duct| {
        let dg = DgDuct::new(
            duct,
            &gas,
            mesh,
            DgBoundary::IdealOpen,
            DgBoundary::IdealOpen,
        )
        .expect("DG");
        let prop = dg.propagator(t_final / 6.0, 1e-13).expect("propagator");
        let mut state = dg.project_initial_state(p0, |_| 0.0);
        for _ in 0..6 {
            prop.step(&mut state, &[]).expect("step");
        }
        prop.sync_nodal(&mut state).unwrap();
        state.pressure().to_vec()
    };
    let (cone, profile) = (run(&cone_fixture()), run(&profile_duct));
    let scale = cone.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
    let diff = cone
        .iter()
        .zip(&profile)
        .fold(0.0_f64, |a, (x, y)| a.max((x - y).abs()));
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dg_profile_equals_cone\",\"inputs\":\"power_profile_exponent_1_vs_cone\",\"elements\":8,\"order\":4,\"observed\":{:.3e},\"tol\":1e-12,\"verdict\":\"{}\"}}",
        diff / scale,
        if diff <= 1e-12 * scale {
            "pass"
        } else {
            "fail"
        }
    );
    assert!(
        diff <= 1e-12 * scale,
        "Segment::Profile linear taper vs Segment::Cone: {:e}",
        diff / scale
    );
}

/// `½‖M^{1/2}[p; q]‖²` of the synced nodal state: the energy computed independently of the
/// propagator's internal coordinates.
fn nodal_energy(dg: &DgDuct, prop: &DgPropagator, state: &mut DgState) -> f64 {
    prop.sync_nodal(state).unwrap();
    let tilde = dg.to_energy_coords(state.pressure(), state.volume_flow());
    0.5 * tilde.iter().map(|v| v * v).sum::<f64>()
}

/// G0: 10^5 exact steps keep the physical energy to 1e-12 relative on closed/closed and
/// closed/open ducts; the ledger's energy is the reconstructed state's (checked against the
/// nodal state). Two ports exchange exactly `e · ΔX`, and the balance closes to 1e-12. A matched
/// end is passive, and its ledger is marked as not a conservation check.
#[test]
fn dg_energy_ledger_lossless() {
    let start = Instant::now();
    let gas = dry_air();
    let duct = cone_fixture();
    let dt = 1.0 / 48_000.0;
    let mesh = DgMesh {
        elements: 6,
        order: 3,
    };
    for (label, outlet) in [
        ("closed_closed", DgBoundary::Closed),
        ("closed_open", DgBoundary::IdealOpen),
    ] {
        let t0 = Instant::now();
        let dg = DgDuct::new(&duct, &gas, mesh, DgBoundary::Closed, outlet).expect("DG duct");
        let prop = dg.propagator(dt, 1e-13).expect("propagator");
        let mut state = dg.project_initial_state(
            |x| 50.0 * det::sin(core::f64::consts::PI * x / 0.60),
            |x| 1e-4 * det::cos(core::f64::consts::PI * x / 0.60),
        );
        let e0 = state.initial_energy_j();
        assert!(e0 > 0.0);
        assert!((nodal_energy(&dg, &prop, &mut state) - e0).abs() <= 1e-15 * e0);
        let mut max_drift = 0.0f64;
        let mut last = 0.0;
        for _ in 0..100_000 {
            let out = prop.step(&mut state, &[]).expect("step");
            assert!(out.ledger.balance_checked);
            max_drift = max_drift.max(out.ledger.relative_energy_drift);
            last = out.ledger.stored_energy_j;
        }
        let e_nodal = nodal_energy(&dg, &prop, &mut state);
        let ledger_vs_nodal = (e_nodal - last).abs() / e0;
        max_drift = max_drift.max((e_nodal - e0).abs() / e0);
        eprintln!(
            "{{\"suite\":\"fs-duct\",\"case\":\"dg_energy_ledger_lossless\",\"inputs\":\"{label}_cone_4to20mm_100k_steps\",\"elements\":{},\"order\":{},\"dofs\":{},\"dt\":{dt:.4e},\"observed\":{max_drift:.3e},\"ledger_vs_nodal\":{ledger_vs_nodal:.3e},\"certificate\":{:.2e},\"oracle\":0.0,\"tol\":1e-12,\"verdict\":\"{}\",\"elapsed_us\":{}}}",
            mesh.elements,
            mesh.order,
            dg.total_dofs(),
            prop.certificate_residual(),
            if max_drift < 1e-12 { "pass" } else { "fail" },
            t0.elapsed().as_micros()
        );
        assert!(
            ledger_vs_nodal <= 1e-13,
            "{label}: ledger vs nodal energy {ledger_vs_nodal:e}"
        );
        assert!(
            max_drift < 1e-12,
            "{label}: 10^5-step energy drift {max_drift:e} >= 1e-12"
        );
    }

    // Two ports exchanging work over 500 steps.
    let (p_in, p_out) = (PortId(1), PortId(2));
    let dg = DgDuct::new(
        &duct,
        &gas,
        mesh,
        DgBoundary::Port(p_in),
        DgBoundary::Port(p_out),
    )
    .expect("two-port DG duct");
    let prop = dg.propagator(dt, 1e-13).expect("two-port propagator");
    let mut state = dg.zero_state();
    let (mut work, mut max_drift) = (0.0f64, 0.0f64);
    for step in 0..500usize {
        let t = (step as f64) * dt;
        let u1 = 2.0e-5 * det::sin(core::f64::consts::TAU * 440.0 * t);
        let u2 = -1.0e-5 * det::cos(core::f64::consts::TAU * 660.0 * t);
        let inputs = [
            PortInput {
                id: p_in,
                effort: u1,
            },
            PortInput {
                id: p_out,
                effort: u2,
            },
        ];
        let out = prop.step(&mut state, &inputs).expect("port step");
        let w = u1 * out.port_outputs[0].delta_x + u2 * out.port_outputs[1].delta_x;
        assert!(
            (w - out.ledger.port_work_step_j).abs() <= 1e-15 * w.abs().max(f64::MIN_POSITIVE),
            "step port work must equal sum(e * delta_x)"
        );
        work += w;
        max_drift = max_drift.max(out.ledger.relative_energy_drift);
    }
    let e_nodal = nodal_energy(&dg, &prop, &mut state);
    let nodal_balance = (e_nodal - work).abs() / work.abs();
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dg_energy_ledger_ports\",\"inputs\":\"two_port_cone_500_steps\",\"elements\":{},\"order\":{},\"dt\":{dt:.4e},\"port_work_j\":{work:.6e},\"observed\":{max_drift:.3e},\"nodal_balance\":{nodal_balance:.3e},\"tol\":1e-12,\"verdict\":\"{}\"}}",
        mesh.elements,
        mesh.order,
        if max_drift < 1e-12 && nodal_balance < 1e-12 {
            "pass"
        } else {
            "fail"
        }
    );
    assert!(work.abs() > 1e-12, "liveness: the ports must exchange work");
    assert!(
        max_drift < 1e-12,
        "port energy balance drift {max_drift:e} >= 1e-12"
    );
    assert!(
        nodal_balance < 1e-12,
        "nodal energy vs port work {nodal_balance:e}"
    );

    // A matched end absorbs: stored energy never rises, and the ledger says it is not a check.
    let dg = DgDuct::new(&duct, &gas, mesh, DgBoundary::Closed, DgBoundary::Matched)
        .expect("matched DG duct");
    let prop = dg.propagator(dt, 1e-13).expect("matched propagator");
    let mut state = dg.project_initial_state(
        |x| 50.0 * det::exp(-0.5 * ((x - 0.3) / 0.04).powi(2)),
        |_| 0.0,
    );
    let e0 = state.initial_energy_j();
    let mut previous = e0;
    let mut worst_rise = 0.0f64;
    for _ in 0..2000 {
        let out = prop.step(&mut state, &[]).expect("matched step");
        assert!(!out.ledger.balance_checked);
        worst_rise = worst_rise.max((out.ledger.stored_energy_j - previous) / e0);
        previous = out.ledger.stored_energy_j;
    }
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dg_energy_ledger_matched\",\"inputs\":\"closed_matched_cone_2000_steps\",\"observed\":{{\"worst_rise\":{worst_rise:.3e},\"final_over_initial\":{:.4e}}},\"tol\":1e-13,\"elapsed_us\":{}}}",
        previous / e0,
        start.elapsed().as_micros()
    );
    assert!(
        worst_rise <= 1e-13,
        "a matched end must be passive: rise {worst_rise:e}"
    );
    assert!(
        previous < 0.5 * e0,
        "liveness: the matched end must absorb the pulse"
    );
}

/// In-place radix-2 Cooley-Tukey FFT (forward: sum x[n] exp(+i 2 pi k n / N) matching e^{-i omega t}).
fn fft_radix2_minus_i_omega(buf: &mut [C64]) {
    let n = buf.len();
    assert!(n.is_power_of_two());
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            buf.swap(i, j);
        }
    }
    let mut len = 2usize;
    while len <= n {
        let half = len / 2;
        let ang = core::f64::consts::TAU / (len as f64);
        for start in (0..n).step_by(len) {
            for k in 0..half {
                let theta = ang * (k as f64);
                let w = C64::new(det::cos(theta), det::sin(theta));
                let u = buf[start + k];
                let v = buf[start + k + half] * w;
                buf[start + k] = u + v;
                buf[start + k + half] = u - v;
            }
        }
        len <<= 1;
    }
}

/// G2: FFT of a converged DG impulse response vs `fs_duct::input_impedance(.., LossModel::Lossless, ..)`
/// matches within 0.1 cent on the first 6 resonance peaks.
#[test]
fn dg_impedance_matches_tmm_lossless() {
    let t0 = Instant::now();
    let gas = dry_air();
    let duct = cylinder_fixture();
    let length = 0.57;
    let port_in = PortId(0);
    let mesh = DgMesh {
        elements: 20,
        order: 4,
    };
    let dg = DgDuct::new(
        &duct,
        &gas,
        mesh,
        DgBoundary::Port(port_in),
        DgBoundary::IdealOpen,
    )
    .expect("DG cylinder for impedance");

    let fs = 65_536.0;
    let dt = 1.0 / fs;
    let prop = dg.propagator(dt, 1e-13).expect("certified propagator");
    let f1_analytic = gas.sound_speed / (4.0 * length);
    let guesses: Vec<f64> = (0..6).map(|n| (2 * n + 1) as f64 * f1_analytic).collect();

    // The Lossless TMM pole near each guess. The TMM refuses at the pole itself, so it is
    // located from two flanking samples: near an isolated pole 1/|Z| = k |f - f_p|, hence
    // f_p = (f_1/|Z_2| + f_2/|Z_1|) / (1/|Z_1| + 1/|Z_2|).
    let tmm_pole = |f_guess: f64| -> f64 {
        let inv = |f_hz: f64| {
            let omega = core::f64::consts::TAU * f_hz;
            let z = input_impedance(
                &duct,
                &gas,
                omega,
                LossModel::Lossless,
                Termination::IdealOpen,
            )
            .expect("TMM impedance")
            .impedance;
            1.0 / z.abs()
        };
        let (f_1, f_2) = (f_guess - 0.05, f_guess + 0.07);
        let (s_1, s_2) = (inv(f_1), inv(f_2));
        (f_1 * s_2 + f_2 * s_1) / (s_1 + s_2)
    };
    let cents = |a: f64, b: f64| 1200.0 * det::ln(a / b) / core::f64::consts::LN_2;

    // The semi-discrete impedance's own poles: the modal frequencies sigma_k / 2 pi of the
    // port-coupled modes, lowest first.
    let PropagatorKernel::SkewModal {
        sigmas,
        modal_ports,
        ..
    } = &prop.kernel
    else {
        panic!("a port/open cylinder propagates on the skew kernel");
    };
    let coupling = &modal_ports[0].1;
    let strongest = coupling.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
    let mut modal: Vec<f64> = sigmas
        .iter()
        .zip(coupling)
        .filter(|(s, b)| **s > 0.0 && b.abs() > 1e-6 * strongest)
        .map(|(s, _)| s / core::f64::consts::TAU)
        .collect();
    modal.sort_by(f64::total_cmp);

    // FFT of the windowed impulse response: unit volume-flow impulse 1/dt held over the first
    // step, the port pressure sampled at every step end, window exp(-sigma t) with
    // sigma = 18 / T_ir. 1/|Z|^2 of an isolated pole is a parabola in f around the peak.
    let fft_peaks = |n_ir: usize, n_fft: usize| -> Vec<f64> {
        let damp_sigma = 18.0 / ((n_ir as f64) * dt);
        let mut state = dg.zero_state();
        let mut spec = vec![C64::ZERO; n_fft];
        for n in 0..n_ir {
            let effort = if n == 0 { 1.0 / dt } else { 0.0 };
            let out = prop
                .step(
                    &mut state,
                    &[PortInput {
                        id: port_in,
                        effort,
                    }],
                )
                .expect("step");
            let t = ((n + 1) as f64) * dt;
            spec[n] = C64::new(
                out.port_outputs[0].end_conjugate * det::exp(-damp_sigma * t) * dt,
                0.0,
            );
        }
        fft_radix2_minus_i_omega(&mut spec);
        let df = fs / (n_fft as f64);
        guesses
            .iter()
            .map(|&f_guess| {
                let k_lo = ((f_guess - 25.0) / df).floor().max(2.0) as usize;
                let k_hi = ((f_guess + 25.0) / df).ceil() as usize;
                let m2 = |k: usize| spec[k].re * spec[k].re + spec[k].im * spec[k].im;
                let best = (k_lo..=k_hi).fold(k_lo, |b, k| if m2(k) > m2(b) { k } else { b });
                let (y_m, y_0, y_p) = (1.0 / m2(best - 1), 1.0 / m2(best), 1.0 / m2(best + 1));
                let denom = y_m - 2.0 * y_0 + y_p;
                let shift = if denom.abs() > 1e-30 {
                    0.5 * (y_m - y_p) / denom
                } else {
                    0.0
                };
                ((best as f64) + shift) * df
            })
            .collect()
    };
    // 1 s of response (sigma = 18 /s), and 0.25 s (sigma = 72 /s) to expose the window bias.
    let long = fft_peaks(65_536, 131_072);
    let short = fft_peaks(16_384, 65_536);

    let poles: Vec<f64> = guesses.iter().map(|&f| tmm_pole(f)).collect();
    let rows: Vec<(f64, f64, f64, f64)> = (0..6)
        .map(|n| {
            (
                poles[n],
                cents(modal[n], poles[n]),
                cents(long[n], poles[n]),
                cents(short[n], poles[n]),
            )
        })
        .collect();
    let worst = |pick: fn(&(f64, f64, f64, f64)) -> f64| {
        rows.iter().fold(0.0_f64, |a, r| a.max(pick(r).abs()))
    };
    let (max_modal, max_fft, max_short) = (worst(|r| r.1), worst(|r| r.2), worst(|r| r.3));
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dg_impedance_matches_tmm_lossless\",\"inputs\":\"cylinder_0.57m_r7.3mm_port_open_vs_tmm\",\"elements\":{},\"order\":{},\"dofs\":{},\"dt\":{dt:.4e},\"observed\":{{\"max_cents_modal\":{max_modal:.4e},\"max_cents_fft_1s\":{max_fft:.4e},\"max_cents_fft_0.25s\":{max_short:.4e},\"rows_tmm_hz_modal_fft1s_fft025s_cents\":{rows:?}}},\"oracle\":\"fs_duct::input_impedance_lossless\",\"tol\":0.1,\"verdict\":\"{}\",\"elapsed_us\":{}}}",
        mesh.elements,
        mesh.order,
        dg.total_dofs(),
        if max_modal <= 0.1 && max_fft <= 0.1 {
            "pass"
        } else {
            "fail"
        },
        t0.elapsed().as_micros()
    );
    assert!(
        max_modal <= 0.1,
        "DG poles vs Lossless TMM {rows:?} exceed 0.1 c"
    );
    assert!(
        max_fft <= 0.1,
        "DG FFT peaks vs Lossless TMM {rows:?} exceed 0.1 c"
    );
}

/// G1: the Δt axis with a coupled port. The port holds the step average of a smooth inflow
/// (a Gaussian-windowed 900 Hz burst), the operator itself is exact, so the port pressure at
/// t = 2 ms converges at order 2 in Δt (fs-mms, ±0.2) against 16x the finest rung.
#[test]
fn dg_port_coupling_dt_ladder() {
    let gas = dry_air();
    let port = PortId(3);
    let dg = DgDuct::new(
        &cylinder_fixture(),
        &gas,
        DgMesh {
            elements: 12,
            order: 3,
        },
        DgBoundary::Port(port),
        DgBoundary::IdealOpen,
    )
    .expect("port DG");
    let inflow = |t: f64| {
        1e-5 * det::exp(-((t - 1e-3) / 3e-4).powi(2)) * det::sin(core::f64::consts::TAU * 900.0 * t)
    };
    // Eight-point Gauss–Legendre step average.
    let (gx, gw) = gauss_legendre_rule(8);
    let average = |t0: f64, dt: f64| -> f64 {
        gx.iter()
            .zip(&gw)
            .map(|(x, w)| 0.5 * w * inflow(t0 + 0.5 * dt * (x + 1.0)))
            .sum()
    };
    let t_end = 2e-3;
    let run = |steps: usize| -> f64 {
        let dt = t_end / steps as f64;
        let prop = dg.propagator(dt, 1e-13).expect("propagator");
        let mut state = dg.zero_state();
        let mut last = 0.0;
        for n in 0..steps {
            let effort = average(n as f64 * dt, dt);
            last = prop
                .step(&mut state, &[PortInput { id: port, effort }])
                .expect("step")
                .port_outputs[0]
                .end_conjugate;
        }
        last
    };
    let rungs = [40usize, 80, 160, 320];
    let reference = run(16 * rungs[3]);
    let errs: Vec<f64> = rungs.iter().map(|&n| (run(n) - reference).abs()).collect();
    for (n, e) in rungs.iter().zip(&errs) {
        eprintln!(
            "{{\"suite\":\"fs-duct\",\"case\":\"dg_port_coupling_dt_ladder\",\"inputs\":\"cylinder_port_open_step_averaged_burst\",\"axis\":\"dt\",\"rung\":{n},\"dt\":{:.4e},\"observed\":{e:.6e},\"oracle\":\"16x finest\",\"reference_pa\":{reference:.6e}}}",
            t_end / *n as f64
        );
    }
    assert!(
        errs[0] > 1e-9 * reference.abs(),
        "liveness: the coarse Δt must show error"
    );
    let ladder = RefinementLadder::new(rungs.iter().map(|&n| t_end / n as f64).collect(), errs)
        .expect("dt ladder");
    let verdict = OrderGate { theoretical: 2.0 }
        .check("dg-port-dt-ladder", LadderSide::Primal, &ladder)
        .unwrap_or_else(|e| panic!("port Δt ladder: {e}"));
    eprintln!("{}", verdict.json_line(true));
}

/// G0: the lossless (skew) and matched (series) propagators pass a 1e-13 certificate; an
/// unachievable tolerance refuses by name on both paths.
#[test]
fn dg_propagator_certificate() {
    let t0 = Instant::now();
    let gas = dry_air();
    let duct = cone_fixture();
    let mesh = DgMesh {
        elements: 8,
        order: 4,
    };
    let mut observed = Vec::new();
    for (label, outlet) in [
        ("skew_closed_open", DgBoundary::IdealOpen),
        ("skew_closed_closed", DgBoundary::Closed),
        ("series_closed_matched", DgBoundary::Matched),
    ] {
        let dg = DgDuct::new(&duct, &gas, mesh, DgBoundary::Closed, outlet).expect("DG");
        let prop = dg
            .propagator(1.0e-5, 1.0e-13)
            .unwrap_or_else(|e| panic!("{label}: propagator must pass a 1e-13 certificate: {e}"));
        observed.push((label, prop.certificate_residual()));
        assert!(prop.certificate_residual() <= 1e-13);
        assert!(matches!(
            dg.propagator(1.0e-5, 1.0e-30),
            Err(DuctError::PropagatorCertificate { .. })
        ));
    }
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dg_propagator_certificate\",\"inputs\":\"cone_skew_and_matched\",\"elements\":{},\"order\":{},\"dt\":1e-5,\"observed\":{observed:?},\"tol\":1e-13,\"verdict\":\"pass\",\"elapsed_us\":{}}}",
        mesh.elements,
        mesh.order,
        t0.elapsed().as_micros()
    );
}

/// G0: caps admit at the cap and refuse at cap + 1 (order, elements, DOFs); A <= 0, a tone
/// hole, an empty duct, shared port ids, dt <= 0, tol <= 0, unknown or non-finite port input.
#[test]
fn dg_refusals() {
    let t0 = Instant::now();
    let gas = dry_air();
    let duct = cylinder_fixture();
    let new = |duct: &Duct, elements: usize, order: usize, inlet, outlet| {
        DgDuct::new(duct, &gas, DgMesh { elements, order }, inlet, outlet)
    };
    let (closed, open) = (DgBoundary::Closed, DgBoundary::IdealOpen);
    assert!(new(&duct, 2, MAX_DG_ORDER, closed, closed).is_ok());
    assert!(new(&duct, 2, MAX_DG_ORDER + 1, closed, closed).is_err());
    assert!(new(&duct, 2, 0, closed, closed).is_err());
    // K = MAX with p = 1 sits exactly at the DOF cap; p = 2 is above it.
    assert_eq!(2 * MAX_DG_ELEMENTS * 2, MAX_DG_PROPAGATOR_DOFS);
    assert!(new(&duct, MAX_DG_ELEMENTS, 1, closed, closed).is_ok());
    assert!(new(&duct, MAX_DG_ELEMENTS, 2, closed, closed).is_err());
    assert!(new(&duct, MAX_DG_ELEMENTS + 1, 1, closed, closed).is_err());
    assert!(new(&duct, 0, 2, closed, closed).is_err());
    for radius in [0.0, -0.01, f64::NAN] {
        let bad = Duct {
            segments: vec![Segment::Cylinder {
                radius,
                length: 0.5,
            }],
        };
        assert!(new(&bad, 4, 2, closed, closed).is_err());
    }
    let hole = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.0073,
                length: 0.3,
            },
            Segment::ToneHole {
                hole_radius: 0.003,
                chimney_height: 0.004,
                bore_radius: 0.0073,
                state: crate::HoleState::Open,
            },
        ],
    };
    assert!(new(&hole, 4, 2, closed, open).is_err());
    assert!(matches!(
        new(&Duct { segments: vec![] }, 4, 2, closed, closed),
        Err(DuctError::EmptyDuct)
    ));
    assert!(
        new(
            &duct,
            4,
            2,
            DgBoundary::Port(PortId(5)),
            DgBoundary::Port(PortId(5))
        )
        .is_err()
    );

    let dg = new(&duct, 4, 2, DgBoundary::Port(PortId(9)), open).expect("valid DG");
    assert!(dg.propagator(0.0, 1e-13).is_err());
    assert!(dg.propagator(-1e-5, 1e-13).is_err());
    assert!(dg.propagator(f64::NAN, 1e-13).is_err());
    assert!(dg.propagator(1e-5, 0.0).is_err());
    let prop = dg.propagator(1e-5, 1e-13).expect("propagator");
    let mut state = dg.zero_state();
    let unknown = PortInput {
        id: PortId(10),
        effort: 1.0,
    };
    let non_finite = PortInput {
        id: PortId(9),
        effort: f64::INFINITY,
    };
    assert!(prop.step(&mut state, &[unknown]).is_err());
    assert!(prop.step(&mut state, &[non_finite]).is_err());
    let foreign = new(&duct, 5, 2, closed, open)
        .expect("other DG")
        .zero_state();
    assert!(prop.step(&mut foreign.clone(), &[]).is_err());

    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dg_refusals\",\"inputs\":\"cap_and_cap_plus_1_boundaries\",\"observed\":\"refused\",\"oracle\":\"Err(DuctError)\",\"tol\":0.0,\"verdict\":\"refused\",\"elapsed_us\":{}}}",
        t0.elapsed().as_micros()
    );
}

/// G5: Bitwise-identical repeatability across independent builds and stepping runs.
#[test]
fn dg_bitwise_repeat() {
    let t0 = Instant::now();
    let run_once = || -> (Vec<u64>, Vec<u64>, u64, u64) {
        let gas = dry_air();
        let duct = cone_fixture();
        let port = PortId(7);
        let dg = DgDuct::new(
            &duct,
            &gas,
            DgMesh {
                elements: 6,
                order: 3,
            },
            DgBoundary::Port(port),
            DgBoundary::IdealOpen,
        )
        .expect("DG");
        let prop = dg.propagator(1.0 / 48_000.0, 1e-13).expect("prop");
        let mut state = dg.project_initial_state(|x| 10.0 * x, |_| 0.0);
        let mut last_out = None;
        for step in 0..64usize {
            let u = 1.0e-5 * det::sin(0.2 * (step as f64));
            last_out = Some(
                prop.step(
                    &mut state,
                    &[PortInput {
                        id: port,
                        effort: u,
                    }],
                )
                .expect("step"),
            );
        }
        prop.sync_nodal(&mut state).unwrap();
        let out = last_out.expect("ran 64 steps");
        let p_bits = state.pressure().iter().map(|v| v.to_bits()).collect();
        let q_bits = state.volume_flow().iter().map(|v| v.to_bits()).collect();
        (
            p_bits,
            q_bits,
            out.port_outputs[0].delta_x.to_bits(),
            out.ledger.stored_energy_j.to_bits(),
        )
    };

    let a = run_once();
    let b = run_once();
    assert_eq!(a, b, "DG propagation must be bitwise repeatable");
    eprintln!(
        "{{\"suite\":\"fs-duct\",\"case\":\"dg_bitwise_repeat\",\"inputs\":\"cone_port_64_steps\",\"observed\":0,\"oracle\":0,\"tol\":0,\"verdict\":\"pass\",\"elapsed_us\":{}}}",
        t0.elapsed().as_micros()
    );
}

/// Equal dimensions do not make energy coordinates portable across materials,
/// geometry or boundaries. Every refusing consumer must leave the state intact.
#[test]
fn dg_state_binding_refuses_foreign_physics_transactionally() {
    let geometry = cylinder_fixture();
    let gas = dry_air();
    let mesh = DgMesh {
        elements: 3,
        order: 2,
    };
    let make = |geometry: &Duct, gas: &GasState, inlet, outlet| {
        DgDuct::new(geometry, gas, mesh, inlet, outlet).unwrap()
    };
    let source = make(&geometry, &gas, DgBoundary::Closed, DgBoundary::Closed);
    let source_prop = source.propagator(1.0 / 48000.0, 1e-13).unwrap();
    let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 330.0, 101325.0).unwrap();
    let dense = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 202650.0).unwrap();
    let longer = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0073,
            length: 0.7,
        }],
    };
    let wider = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.009,
            length: 0.57,
        }],
    };
    let others = [
        make(&geometry, &hot, DgBoundary::Closed, DgBoundary::Closed),
        make(&geometry, &dense, DgBoundary::Closed, DgBoundary::Closed),
        make(&longer, &gas, DgBoundary::Closed, DgBoundary::Closed),
        make(&wider, &gas, DgBoundary::Closed, DgBoundary::Closed),
        make(&geometry, &gas, DgBoundary::IdealOpen, DgBoundary::Closed),
        make(&geometry, &gas, DgBoundary::Closed, DgBoundary::Matched),
    ];
    for other in others {
        assert_eq!(source.total_dofs(), other.total_dofs());
        assert_ne!(source.operator_identity(), other.operator_identity());
        let propagator = other.propagator(1.0 / 48000.0, 1e-13).unwrap();
        for dirty in [false, true] {
            let mut state = source.project_initial_state(|x| 1.0 + x, |_| 0.0);
            if dirty {
                source_prop.step(&mut state, &[]).unwrap();
            }
            let before = state.clone();
            assert!(propagator.step(&mut state, &[]).is_err());
            assert_eq!(state, before);
            assert!(propagator.sync_nodal(&mut state).is_err());
            assert_eq!(state, before);
            assert!(other.l2_pressure_error(&state, |_| 0.0).is_err());
            assert_eq!(state, before);
        }
    }

    // Renaming a port leaves the spatial matrix unchanged but changes its
    // coupling contract, so its identity must participate independently.
    let a = make(
        &geometry,
        &gas,
        DgBoundary::Port(PortId(7)),
        DgBoundary::Closed,
    );
    let b = make(
        &geometry,
        &gas,
        DgBoundary::Port(PortId(8)),
        DgBoundary::Closed,
    );
    assert_eq!(a.l_tilde, b.l_tilde);
    assert_ne!(a.operator_identity(), b.operator_identity());
    let mut state = a.zero_state();
    let before = state.clone();
    assert!(
        b.propagator(1.0 / 48000.0, 1e-13)
            .unwrap()
            .step(&mut state, &[])
            .is_err()
    );
    assert_eq!(state, before);
}

/// The binding is reproducible content identity, so separately assembled
/// copies remain compatible without relying on addresses or instance counters.
#[test]
fn dg_state_binding_accepts_independent_identical_assemblies() {
    let build = || {
        DgDuct::new(
            &cylinder_fixture(),
            &dry_air(),
            DgMesh {
                elements: 4,
                order: 3,
            },
            DgBoundary::Closed,
            DgBoundary::IdealOpen,
        )
        .unwrap()
    };
    let a = build();
    let b = build();
    assert_eq!(a.operator_identity(), b.operator_identity());
    let mut state = a.project_initial_state(|x| 1.0 - x, |_| 0.0);
    let mut expected = state.clone();
    let pa = a.propagator(1.0 / 48000.0, 1e-13).unwrap();
    let pb = b.propagator(1.0 / 48000.0, 1e-13).unwrap();
    for _ in 0..17 {
        assert_eq!(
            pa.step(&mut expected, &[]).unwrap(),
            pb.step(&mut state, &[]).unwrap()
        );
    }
    pb.sync_nodal(&mut state).unwrap();
    pa.sync_nodal(&mut expected).unwrap();
    assert_eq!(state, expected);
    assert_eq!(
        a.l2_pressure_error(&state, |_| 0.0).unwrap(),
        b.l2_pressure_error(&state, |_| 0.0).unwrap(),
    );
}

/// Changing dt must preserve the already elapsed phase, physical state and
/// cumulative port work instead of reinterpreting all old steps with the new dt.
#[test]
fn dg_state_binding_step_size_changes_preserve_evolution_and_work() {
    let field = crate::profile::gas_field::tests::field([300.0, 330.0], [0.8; 2], 101325.0);
    for with_port in [false, true] {
        let inlet = if with_port {
            DgBoundary::Port(PortId(7))
        } else {
            DgBoundary::Closed
        };
        let dg = DgDuct::new_with_gas_field(
            &cylinder_fixture(),
            &field,
            DgMesh {
                elements: 3,
                order: 3,
            },
            inlet,
            DgBoundary::Closed,
        )
        .unwrap();
        let dt = 1.0 / 192000.0;
        let fine = dg.propagator(dt, 1e-13).unwrap();
        let coarse = dg.propagator(2.0 * dt, 1e-13).unwrap();
        let mut state =
            dg.project_initial_state(|x| det::cos(core::f64::consts::PI * x / 0.57), |_| 0.0);
        let mut expected = state.clone();
        let input = if with_port {
            vec![PortInput {
                id: PortId(7),
                effort: 1e-8,
            }]
        } else {
            vec![]
        };
        for (multiple, steps) in [(2, 6), (1, 9), (2, 4)] {
            let propagator = if multiple == 1 { &fine } else { &coarse };
            // Invalid inputs must not reset the cache when switching dt.
            let before = state.clone();
            assert!(
                propagator
                    .step(
                        &mut state,
                        &[PortInput {
                            id: PortId(999),
                            effort: 0.0
                        }]
                    )
                    .is_err()
            );
            assert_eq!(state, before);
            for _ in 0..steps {
                propagator.step(&mut state, &input).unwrap();
                for _ in 0..multiple {
                    fine.step(&mut expected, &input).unwrap();
                }
            }
            let error = det::sqrt(
                state
                    .tilde
                    .iter()
                    .zip(&expected.tilde)
                    .map(|(a, b)| (a - b) * (a - b))
                    .sum::<f64>(),
            );
            let scale = det::sqrt(expected.tilde.iter().map(|v| v * v).sum::<f64>());
            assert!(
                error / scale < 2e-11,
                "relative state error = {}",
                error / scale
            );
            let work_scale = expected.initial_energy_j + expected.cumulative_port_work_j.abs();
            assert!(
                (state.cumulative_port_work_j - expected.cumulative_port_work_j).abs() / work_scale
                    < 2e-11
            );
            assert!(
                (state.stored_energy_j() - expected.stored_energy_j()).abs() / work_scale < 2e-11
            );
        }
        assert_eq!(state.step_count(), 19);
        assert_eq!(expected.step_count(), 29);
    }
}
