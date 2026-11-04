//! Original 1D field owner with fixed molecular composition/thermal baths.
use super::*;
use crate::multicomponent::Cell;

fn ends() -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

fn balance(report: &DissipationReport) {
    for k in 0..report.before.species_count + 2 {
        close(report.conservation_residual.get(k), 0.0, 8e-12);
        close(
            report.boundary_change.get(k),
            report.left_boundary_change.get(k) + report.right_boundary_change.get(k),
            8e-12,
        );
    }
}

fn advance(
    stabilized: bool,
    op: &mut DissipationOperator<ConstantTransport>,
    flow: &mut MixtureEuler1d,
    duration: f64,
    controls: DissipationControl,
    baths: [Option<Cell<'_>>; 2],
    cancelled: impl FnMut() -> bool,
) -> Result<DissipationReport, DissipationError> {
    if stabilized {
        op.advance_stabilized_with_reservoirs_for(
            flow,
            Time::new(duration),
            StabilizedDissipationControl {
                transport: controls,
                ..Default::default()
            },
            baths,
            cancelled,
        )
        .map(|report| report.transport)
    } else {
        op.advance_with_reservoirs_for(flow, Time::new(duration), controls, baths, cancelled)
    }
}

#[test]
fn g1_binary_baths_follow_dirichlet_decay_to_linear_steady_diffusion() {
    let count = 20;
    let diffusivity = 0.1;
    let duration = 8.0;
    let pi = core::f64::consts::PI;
    // Exact eigenvalue of the independently written cell-centered Dirichlet
    // Laplacian: the boundary row is (q1 - 3*q0 + 2*qb)/dx².
    let eigenvalue =
        -4.0 * diffusivity * (count as f64 * (pi / (2.0 * count as f64)).sin()).powi(2);
    let amplitude = 0.12 * (eigenvalue * duration).exp();
    for stabilized in [false, true] {
        let g = gas(2, 0.0, false);
        let left = g.conserved(&[0.8, 0.2], 0.0, 2.0).unwrap();
        let right = g.conserved(&[0.2, 0.8], 0.0, 2.0).unwrap();
        let mut flow = field(g, count, ends(), |x| {
            let a = 0.8 - 0.6 * x + 0.12 * (pi * x).sin();
            (vec![a, 1.0 - a], 0.0, 2.0)
        });
        let mut op = operator(&flow, 0.0, 0.0, 0.0, diffusivity);
        let report = advance(
            stabilized,
            &mut op,
            &mut flow,
            duration,
            DissipationControl {
                maximum_step_s: 0.05,
                ..Default::default()
            },
            [Some(left.as_cell()), Some(right.as_cell())],
            || false,
        )
        .unwrap();
        balance(&report);
        assert_eq!(flow.time().value(), 0.0);
        assert_eq!(report.elapsed.value(), duration);
        assert!(report.left_boundary_change.species_masses[0] > 0.0);
        assert!(report.right_boundary_change.species_masses[0] < 0.0);
        close(report.left_boundary_change.mass(), 0.0, 8e-13);
        close(report.right_boundary_change.mass(), 0.0, 8e-13);
        close(report.after.mass(), report.before.mass(), 8e-13);
        close(report.after.energy, report.before.energy, 8e-13);
        for i in 0..count {
            let x = (i as f64 + 0.5) / count as f64;
            let linear = 0.8 - 0.6 * x;
            let q = flow.cell(i).unwrap();
            close(
                q.partial_densities[0],
                linear + amplitude * (pi * x).sin(),
                4e-7,
            );
            close(q.partial_densities[0], linear, 5e-5);
            close(q.density(), 1.0, 8e-13);
            close(
                flow.gas().primitive(q, || false).unwrap().temperature,
                2.0,
                2e-11,
            );
        }
    }
}

#[test]
fn g0_nonuniform_thermal_baths_use_actual_half_cell_distances_and_signed_heat() {
    for stabilized in [false, true] {
        let g = gas(1, 0.0, false);
        let widths = [0.2, 0.3, 0.5];
        let temperatures = [2.1, 2.35, 2.75];
        let states = temperatures
            .iter()
            .map(|t| g.conserved(&[1.0], 0.0, *t).unwrap())
            .collect();
        let left = g.conserved(&[1.0], 0.0, 2.0).unwrap();
        let right = g.conserved(&[1.0], 0.0, 3.0).unwrap();
        let mesh = Mesh1d::new(
            vec![0.0, 0.2, 0.5, 1.0],
            vec![0.7; 4],
            widths.iter().map(|w| 0.7 * w).collect(),
        )
        .unwrap();
        let mut flow = MixtureEuler1d::new_on_mesh(g, mesh, states, ends(), || false).unwrap();
        let initial = flow.states.clone();
        let mut op = operator(&flow, 0.0, 0.0, 0.05, 0.0);
        let report = advance(
            stabilized,
            &mut op,
            &mut flow,
            0.04,
            DissipationControl::default(),
            [Some(left.as_cell()), Some(right.as_cell())],
            || false,
        )
        .unwrap();
        balance(&report);
        // Independently fixed Fourier flux -lambda*(3-2)/length at both ends.
        close(
            report.left_boundary_change.energy,
            -0.05 * 0.7 * 0.04,
            2e-13,
        );
        close(
            report.right_boundary_change.energy,
            0.05 * 0.7 * 0.04,
            2e-13,
        );
        close(report.boundary_change.energy, 0.0, 2e-13);
        assert_eq!(report.boundary_change.mass(), 0.0);
        for (&actual, expected) in flow.states.iter().zip(initial) {
            close(actual, expected, 2e-13);
        }
    }
}

#[test]
fn g0_bath_enthalpy_and_zero_stationary_traction_work_ignore_formation_gauge() {
    for stabilized in [false, true] {
        let mut baseline = None;
        for formation in [0.0, 1e18] {
            let g = gas(2, formation, true);
            // A has R_specific=1, Cp=3.5; B has R_specific=.5, Cp=2.
            // The source energy gauge subtracts Cv*Tmin for each species.
            let bath_temperature = 3.0;
            let bath = g
                .conserved(&[0.2, 0.8], 0.0, 0.6 * bath_temperature)
                .unwrap();
            let mut flow = field(g, 8, ends(), |_| (vec![0.7, 0.3], 0.4, 1.7));
            let mut op = operator(&flow, 0.02, 0.01, 0.0, 0.03);
            let report = advance(
                stabilized,
                &mut op,
                &mut flow,
                0.03,
                DissipationControl {
                    maximum_step_s: 0.01,
                    ..Default::default()
                },
                [Some(bath.as_cell()), None],
                || false,
            )
            .unwrap();
            balance(&report);
            let transfer = report.left_boundary_change;
            assert!(transfer.species_masses[0].abs() > 1e-6);
            assert!(transfer.momentum < 0.0);
            close(transfer.mass(), 0.0, 2e-13);
            let enthalpy = [
                3.5 * bath_temperature - 2.5 * 0.001,
                2.0 * bath_temperature - 1.5 * 0.001,
            ];
            close(
                transfer.energy,
                transfer.species_masses[0] * enthalpy[0] + transfer.species_masses[1] * enthalpy[1],
                2e-12,
            );
            // No right bath and a transmissive end: all exchange is left.
            assert_eq!(report.right_boundary_change.energy, 0.0);
            assert_eq!(report.right_boundary_change.momentum, 0.0);
            if let Some((states, boundary)) = &baseline {
                assert_eq!(&flow.states, states);
                assert_eq!(&report.boundary_change, boundary);
            } else {
                baseline = Some((flow.states.clone(), report.boundary_change));
                close(report.absolute_energy_residual_j, 0.0, 2e-12);
            }
        }
    }
}

#[test]
fn g4_reservoir_refusal_late_cancellation_and_work_limits_preserve_field_and_reuse() {
    for stabilized in [false, true] {
        let g = gas(2, 0.0, false);
        let bath = g.conserved(&[0.8, 0.2], 0.0, 2.0).unwrap();
        let controls = DissipationControl {
            maximum_step_s: 0.01,
            ..Default::default()
        };
        let make = || field(g.clone(), 8, ends(), |_| (vec![0.5, 0.5], 0.0, 2.0));
        let mut expected = make();
        let mut reference = operator(&expected, 0.0, 0.0, 0.0, 0.03);
        let mut polls = 0;
        let report = advance(
            stabilized,
            &mut reference,
            &mut expected,
            0.05,
            controls,
            [Some(bath.as_cell()), None],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
        assert!(report.accepted_steps > 1);
        let mut flow = make();
        let initial = flow.states.clone();
        let mut op = operator(&flow, 0.0, 0.0, 0.0, 0.03);
        let mut seen = 0;
        assert_eq!(
            advance(
                stabilized,
                &mut op,
                &mut flow,
                0.05,
                controls,
                [Some(bath.as_cell()), None],
                || {
                    seen += 1;
                    seen == polls
                }
            ),
            Err(DissipationError::Cancelled)
        );
        assert_eq!(flow.states, initial);
        assert_eq!(flow.time().value(), 0.0);
        let exhausted = advance(
            stabilized,
            &mut op,
            &mut flow,
            0.05,
            DissipationControl {
                maximum_steps: 1,
                ..controls
            },
            [Some(bath.as_cell()), None],
            || false,
        );
        assert!(matches!(exhausted, Err(DissipationError::WorkLimit { .. })));
        assert_eq!(flow.states, initial);
        let invalid_densities = [1.0];
        for invalid in [
            Cell {
                momentum: 0.1,
                ..bath.as_cell()
            },
            Cell {
                partial_densities: &invalid_densities,
                ..bath.as_cell()
            },
            Cell {
                energy: -1.0,
                ..bath.as_cell()
            },
        ] {
            assert!(
                advance(
                    stabilized,
                    &mut op,
                    &mut flow,
                    0.05,
                    controls,
                    [Some(invalid), None],
                    || false
                )
                .is_err()
            );
            assert_eq!(flow.states, initial);
        }
        let repeated = advance(
            stabilized,
            &mut op,
            &mut flow,
            0.05,
            controls,
            [Some(bath.as_cell()), None],
            || false,
        )
        .unwrap();
        assert_eq!(flow.states, expected.states);
        assert_eq!(report, repeated);
        let mut periodic = field(g.clone(), 8, Boundaries::Periodic, |_| {
            (vec![0.5, 0.5], 0.0, 2.0)
        });
        let initial = periodic.states.clone();
        let mut op = operator(&periodic, 0.0, 0.0, 0.0, 0.03);
        assert!(matches!(
            advance(
                stabilized,
                &mut op,
                &mut periodic,
                0.05,
                controls,
                [Some(bath.as_cell()), None],
                || false
            ),
            Err(DissipationError::Invalid(_))
        ));
        assert_eq!(periodic.states, initial);
    }
}

#[test]
fn g0_absent_molecular_baths_preserve_original_reports_and_poll_sequence() {
    for stabilized in [false, true] {
        for boundaries in [ends(), Boundaries::Periodic] {
            let g = gas(1, 0.0, false);
            let make = || {
                field(g.clone(), 8, boundaries, |x| {
                    (vec![1.0], 0.0, 2.0 + 0.2 * x)
                })
            };
            let mut a = make();
            let mut b = make();
            let mut oa = operator(&a, 0.0, 0.0, 0.05, 0.0);
            let mut ob = operator(&b, 0.0, 0.0, 0.05, 0.0);
            let (mut pa, mut pb) = (0, 0);
            let controls = DissipationControl::default();
            let ra = if stabilized {
                oa.advance_stabilized_for(
                    &mut a,
                    Time::new(0.01),
                    StabilizedDissipationControl::default(),
                    || {
                        pa += 1;
                        false
                    },
                )
                .unwrap()
                .transport
            } else {
                oa.advance_for(&mut a, Time::new(0.01), controls, || {
                    pa += 1;
                    false
                })
                .unwrap()
            };
            let rb = advance(
                stabilized,
                &mut ob,
                &mut b,
                0.01,
                controls,
                [None; 2],
                || {
                    pb += 1;
                    false
                },
            )
            .unwrap();
            assert_eq!(a.states, b.states);
            assert_eq!(ra, rb);
            assert_eq!(pa, pb);
        }
    }
}
