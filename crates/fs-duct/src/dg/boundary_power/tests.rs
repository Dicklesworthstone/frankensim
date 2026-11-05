use super::*;
use fs_material::gas::GasSpec;

const PORT: PortId = PortId(37);
const RADIUS: f64 = 0.01;
const LENGTH: f64 = 0.3;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn duct(inlet: DgBoundary, outlet: DgBoundary) -> DgDuct {
    DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH,
            }],
        },
        &gas(),
        DgMesh {
            elements: 1,
            order: 2,
        },
        inlet,
        outlet,
    )
    .unwrap()
}

fn impedance() -> f64 {
    let gas = gas();
    gas.density * gas.sound_speed / (core::f64::consts::PI * RADIUS * RADIUS)
}

fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 << 20,
        max_multiply_terms: 40_000_000_000,
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn source_coefficients_partition_the_original_boundary_power() {
    let zc = 1.3e6;
    for ratio in [1e-12, 0.02, 1.0, 30.0, 1e12] {
        let resistance = ratio * zc;
        let boundary = DgBoundary::resistive(resistance).unwrap();
        for normal in [-1.0, 1.0] {
            let (load, mismatch) = source_coefficients(boundary, zc, normal).unwrap().unwrap();
            for (p, q) in [(23.0, 2e-6), (-7.0, 4e-6), (10.0, 0.0)] {
                let physical = (load[0] * p + load[1] * q).powi(2);
                let numerical = (mismatch[0] * p + mismatch[1] * q).powi(2);
                let flux = (p + normal * zc * q) / (resistance + zc);
                close(physical, resistance * flux * flux, 2e-14);
                let native =
                    p * p / (resistance + zc) + resistance * zc / (resistance + zc) * q * q;
                close(physical + numerical, native, 2e-14);
            }
        }
    }
    for boundary in [
        DgBoundary::Closed,
        DgBoundary::IdealOpen,
        DgBoundary::Port(PORT),
    ] {
        assert_eq!(source_coefficients(boundary, zc, 1.0).unwrap(), None);
    }
}

/// The two-coordinate boundary quadratic has independent exponential integrals.
/// This checks the new source rows against actual native source integration and
/// DenseSeries phi functions, including weak physical load and interference.
#[test]
fn source_integrals_agree_with_analytic_decay_and_native_dense_exponential() {
    let zc = 1.3e6;
    let duration = 0.7;
    let initial = [0.6, -0.2];
    for ratio in [1e-12, 0.02, 1.0, 30.0, 1e12] {
        let boundary = DgBoundary::resistive(ratio * zc).unwrap();
        let (alpha, beta) = boundary.characteristic_weights(zc).unwrap();
        let generator = [-beta, 0.0, 0.0, -alpha];
        let (phi0, _, _, _) = series_phis(
            &generator,
            2,
            duration,
            series_shift(&generator, 2, duration),
        )
        .unwrap();
        let endpoint = [
            phi0[0] * initial[0] + phi0[1] * initial[1],
            phi0[2] * initial[0] + phi0[3] * initial[1],
        ];
        close(endpoint[0], initial[0] * det::exp(-beta * duration), 1e-13);
        close(endpoint[1], initial[1] * det::exp(-alpha * duration), 1e-13);
        for normal in [-1.0, 1.0] {
            let (load, mismatch) = source_coefficients(boundary, zc, normal).unwrap().unwrap();
            // Physical p=sqrt(z)*x0 and axial q=x1/sqrt(z).
            let rows = [load, mismatch].map(|row| [row[0] * det::sqrt(zc), row[1] / det::sqrt(zc)]);
            let integrals = rows.map(|row| {
                PreparedSourceLoss::new_controlled(
                    2,
                    0,
                    &generator,
                    &[],
                    &row,
                    duration,
                    budget(),
                    || false,
                )
                .unwrap()
                .evaluate_controlled(&initial, &[], || false)
                .unwrap()
            });
            let integral = |rate: f64| -det::expm1(-rate * duration) / rate;
            let cross = 2.0 * normal * initial[0] * initial[1] * integral(alpha + beta);
            let physical = alpha
                * beta
                * (initial[0] * initial[0] * integral(2.0 * beta)
                    + initial[1] * initial[1] * integral(2.0 * alpha)
                    + cross);
            let numerical = beta * beta * initial[0] * initial[0] * integral(2.0 * beta)
                + alpha * alpha * initial[1] * initial[1] * integral(2.0 * alpha)
                - alpha * beta * cross;
            close(integrals[0], physical, 3e-12);
            close(integrals[1], numerical, 3e-12);
            let stored_change = 0.5
                * (initial[0] * initial[0] + initial[1] * initial[1]
                    - endpoint[0] * endpoint[0]
                    - endpoint[1] * endpoint[1]);
            close(integrals[0] + integrals[1], stored_change, 3e-13);
        }
    }
}

/// G1: uniform steady flow into a physical resistance has load power R*Q^2
/// and no boundary mismatch. Check both normals and strongly unequal loads.
#[test]
fn steady_load_work_is_physical_pressure_times_outward_flow() {
    let dt = 3e-4;
    let flow = 2e-6;
    for ratio in [0.02, 1.0, 30.0] {
        let resistance = ratio * impedance();
        for load_at_inlet in [false, true] {
            let load = if ratio == 1.0 {
                DgBoundary::Matched
            } else {
                DgBoundary::resistive(resistance).unwrap()
            };
            let (inlet, outlet, axial) = if load_at_inlet {
                (load, DgBoundary::Port(PORT), -flow)
            } else {
                (DgBoundary::Port(PORT), load, flow)
            };
            let duct = duct(inlet, outlet);
            let propagator = duct.propagator(dt, 1e-10).unwrap();
            let energy = duct
                .prepare_boundary_energy(&propagator, budget(), &mut || false)
                .unwrap();
            let mut state = duct.project_initial_state(|_| resistance * flow, |_| axial);
            let inputs = [
                PortInput {
                    id: PORT,
                    effort: 0.75 * flow,
                },
                PortInput {
                    id: PORT,
                    effort: 0.25 * flow,
                },
            ];
            let observed = duct
                .step_with_boundary_energy(
                    &propagator,
                    &energy,
                    &mut state,
                    &inputs,
                    1e-20,
                    1e-9,
                    &mut || false,
                )
                .unwrap();
            close(
                observed.boundary.load_absorbed_j,
                resistance * flow * flow * dt,
                3e-10,
            );
            assert!(observed.boundary.numerical_loss_j < 1e-18 * observed.boundary.load_absorbed_j);
            close(
                observed.step.ledger.port_work_step_j,
                resistance * flow * flow * dt,
                3e-10,
            );
            assert!(observed.step.ledger.balance_checked);
            assert_eq!(state.step_count(), 1);
            for &pressure in state.pressure() {
                close(pressure, resistance * flow, 3e-9);
            }
        }
    }
}

/// Independent physical flux reconstruction, deliberately not the source rows.
fn boundary_powers(duct: &DgDuct, state: &DgState) -> [f64; 2] {
    let (pressure, flow) = duct.from_energy_coords(&state.tilde);
    let z = impedance();
    let mut result = [0.0; 2];
    for (boundary, index, normal) in [
        (duct.inlet, 0, -1.0),
        (duct.outlet, duct.field_dofs - 1, 1.0),
    ] {
        let resistance = match boundary {
            DgBoundary::Matched => z,
            DgBoundary::Resistive(r) => r.resistance_pa_s_m3(),
            _ => continue,
        };
        let p = pressure[index];
        let q = normal * flow[index];
        let load_flow = (p + z * q) / (resistance + z);
        let load_pressure = resistance * load_flow;
        result[0] += load_pressure * load_flow;
        result[1] += ((p - load_pressure).powi(2) + z * z * (q - load_flow).powi(2)) / (2.0 * z);
    }
    result
}

/// G1/G3: independent time quadrature over actual boundary fluxes separates
/// physical absorption from a finite upwind penalty during a transient.
#[test]
fn transient_boundary_integrals_match_flux_quadrature_and_close_native_step() {
    let duct = duct(
        DgBoundary::resistive(0.2 * impedance()).unwrap(),
        DgBoundary::Matched,
    );
    let dt = 4e-4;
    let propagator = duct.propagator(dt, 1e-10).unwrap();
    let prepared = duct
        .prepare_boundary_energy(&propagator, budget(), &mut || false)
        .unwrap();
    let initial = duct.project_initial_state(|x| 40.0 * (1.0 + x / LENGTH), |x| 1e-6 * x / LENGTH);
    let observed = prepared
        .integrated_energy_j(&initial, &[], &mut || false)
        .unwrap();
    let intervals = 256;
    let small = duct.propagator(dt / intervals as f64, 1e-10).unwrap();
    let mut state = initial.clone();
    let mut quadrature = [0.0; 2];
    for index in 0..=intervals {
        let weight = if index == 0 || index == intervals {
            1.0
        } else if index % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let power = boundary_powers(&duct, &state);
        for channel in 0..2 {
            quadrature[channel] += weight * power[channel] * small.dt() / 3.0;
        }
        if index < intervals {
            small.step(&mut state, &[]).unwrap();
        }
    }
    close(observed.load_absorbed_j, quadrature[0], 1e-8);
    close(observed.numerical_loss_j, quadrature[1], 1e-8);
    assert!(observed.numerical_loss_j > 0.01 * observed.total_j());
    close(
        observed.total_j(),
        initial.stored_energy_j() - state.stored_energy_j(),
        2e-10,
    );

    let mut original = initial.clone();
    let unchecked = propagator.step(&mut original, &[]).unwrap();
    assert!(!unchecked.ledger.balance_checked);
    let mut checked = initial.clone();
    let checked_output = duct
        .step_with_boundary_energy(
            &propagator,
            &prepared,
            &mut checked,
            &[],
            1e-20,
            1e-9,
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        checked.tilde, original.tilde,
        "original trajectory is unchanged"
    );
    assert_eq!(checked_output.step.port_outputs, unchecked.port_outputs);
    assert!(checked_output.step.ledger.balance_checked);
    close(checked.cumulative_dissipated_j(), observed.total_j(), 1e-14);
    assert!(checked_output.step.ledger.energy_balance_residual_j.abs() < 1e-9 * observed.total_j());
}

#[test]
fn weak_physical_load_is_not_replaced_by_numerical_loss() {
    let duct = duct(
        DgBoundary::Closed,
        DgBoundary::resistive(1e-12 * impedance()).unwrap(),
    );
    let propagator = duct.propagator(1e-7, 1e-10).unwrap();
    let prepared = duct
        .prepare_boundary_energy(&propagator, budget(), &mut || false)
        .unwrap();
    let state = duct.project_initial_state(|_| 30.0, |_| 0.0);
    let energy = prepared
        .integrated_energy_j(&state, &[], &mut || false)
        .unwrap();
    assert!(energy.load_absorbed_j > 0.0);
    assert!(energy.load_absorbed_j < 1e-10 * energy.numerical_loss_j);
    close(
        energy.load_absorbed_j,
        30.0_f64.powi(2) / impedance() * 1e-12 * propagator.dt(),
        0.01,
    );
}

#[test]
fn no_characteristic_load_has_exactly_zero_boundary_energy() {
    let duct = duct(DgBoundary::Port(PORT), DgBoundary::Closed);
    let propagator = duct.propagator(1e-4, 1e-9).unwrap();
    let prepared = duct
        .prepare_boundary_energy(&propagator, budget(), &mut || false)
        .unwrap();
    let mut state = duct.zero_state();
    let output = duct
        .step_with_boundary_energy(
            &propagator,
            &prepared,
            &mut state,
            &[PortInput {
                id: PORT,
                effort: 1e-6,
            }],
            1e-20,
            1e-8,
            &mut || false,
        )
        .unwrap();
    assert_eq!(output.boundary, DgBoundaryEnergy::default());
    assert!(output.step.ledger.balance_checked);
}

#[test]
fn checked_lossless_steps_preserve_the_native_modal_phase_anchor() {
    let duct = duct(DgBoundary::Closed, DgBoundary::IdealOpen);
    let propagator = duct.propagator(2e-5, 1e-9).unwrap();
    let prepared = duct
        .prepare_boundary_energy(&propagator, budget(), &mut || false)
        .unwrap();
    let mut original = duct.project_initial_state(|x| 10.0 * (1.0 - x / LENGTH), |_| 0.0);
    let mut checked = original.clone();
    for _ in 0..16 {
        propagator.step(&mut original, &[]).unwrap();
        duct.step_with_boundary_energy(
            &propagator,
            &prepared,
            &mut checked,
            &[],
            1e-20,
            1e-8,
            &mut || false,
        )
        .unwrap();
        assert_eq!(checked.tilde, original.tilde);
        assert_eq!(checked.modal_anchor_p, original.modal_anchor_p);
        assert_eq!(checked.modal_anchor_q, original.modal_anchor_q);
        assert_eq!(checked.modal_anchor_dt_bits, original.modal_anchor_dt_bits);
        assert_eq!(checked.unforced_steps, original.unforced_steps);
    }
}

#[test]
fn checked_step_refuses_wrong_clock_and_corrupt_advance_without_publication() {
    let duct = duct(DgBoundary::Port(PORT), DgBoundary::Matched);
    let propagator = duct.propagator(1e-4, 1e-10).unwrap();
    let prepared = duct
        .prepare_boundary_energy(&propagator, budget(), &mut || false)
        .unwrap();
    let original = duct.project_initial_state(|x| 10.0 + x, |_| 0.0);
    let mut state = original.clone();
    let wrong_clock = duct.propagator(2e-4, 1e-10).unwrap();
    assert!(
        duct.step_with_boundary_energy(
            &wrong_clock,
            &prepared,
            &mut state,
            &[],
            0.0,
            1e-9,
            &mut || false,
        )
        .is_err()
    );
    assert_eq!(state, original);

    // A test-only faulty exponential would be hidden by deficit-inferred heat.
    // Independent original boundary-source integration must reject it.
    let mut corrupted = propagator.clone();
    let PropagatorKernel::DenseSeries { phi0, .. } = &mut corrupted.kernel else {
        panic!()
    };
    phi0[0] += 0.1;
    assert!(
        duct.step_with_boundary_energy(
            &corrupted,
            &prepared,
            &mut state,
            &[],
            0.0,
            1e-9,
            &mut || false,
        )
        .is_err()
    );
    assert_eq!(state, original);

    let mut visits = 0;
    let mut candidate = original.clone();
    duct.step_with_boundary_energy(
        &propagator,
        &prepared,
        &mut candidate,
        &[],
        1e-20,
        1e-9,
        &mut || {
            visits += 1;
            false
        },
    )
    .unwrap();
    let mut current = 0;
    assert!(matches!(
        duct.step_with_boundary_energy(
            &propagator,
            &prepared,
            &mut state,
            &[],
            1e-20,
            1e-9,
            &mut || {
                current += 1;
                current == visits
            },
        ),
        Err(DuctError::Cancelled { .. })
    ));
    assert_eq!(state, original);
}

#[test]
fn source_owner_budget_input_and_full_loss_scope_are_enforced() {
    let duct = duct(DgBoundary::Port(PORT), DgBoundary::Matched);
    let propagator = duct.propagator(1e-4, 1e-10).unwrap();
    let prepared = duct
        .prepare_boundary_energy(&propagator, budget(), &mut || false)
        .unwrap();
    let state = duct.zero_state();
    let foreign = super::tests::duct(DgBoundary::Closed, DgBoundary::Matched);
    assert!(
        prepared
            .integrated_energy_j(&foreign.zero_state(), &[], &mut || false)
            .is_err()
    );
    assert!(
        duct.prepare_boundary_energy(
            &foreign.propagator(1e-4, 1e-10).unwrap(),
            budget(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        duct.prepare_boundary_energy(
            &propagator,
            DenseExponentialBudget {
                max_dense_bytes: 1,
                ..budget()
            },
            &mut || false
        )
        .is_err()
    );
    for inputs in [
        vec![PortInput {
            id: PortId(999),
            effort: 0.0,
        }],
        vec![PortInput {
            id: PORT,
            effort: f64::NAN,
        }],
        vec![
            PortInput {
                id: PORT,
                effort: f64::MAX
            };
            2
        ],
    ] {
        assert!(
            prepared
                .integrated_energy_j(&state, &inputs, &mut || false)
                .is_err()
        );
    }
    let wall = duct.with_viscothermal(8, [100.0, 4000.0], 1e-3).unwrap();
    let propagator = wall.propagator(1e-6, 1e-9).unwrap();
    let prepared = wall
        .prepare_boundary_energy(&propagator, budget(), &mut || false)
        .unwrap();
    let mut state = wall.zero_state();
    assert!(
        prepared
            .integrated_energy_j(&state, &[], &mut || false)
            .is_ok()
    );
    assert!(
        wall.step_with_boundary_energy(
            &propagator,
            &prepared,
            &mut state,
            &[],
            1e-20,
            1e-9,
            &mut || false
        )
        .is_err()
    );
    assert_eq!(state.step_count(), 0);
}
