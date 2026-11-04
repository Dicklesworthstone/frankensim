//! Original NSF/load coupling checks; authored numerical ports, not measured walls.
use super::*;
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfSpec, WallThermal};
use fs_material::gas::GasSpec;
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::{ComplianceRelaxationTerm, RelaxationTerm};

pub(super) fn limits(cells: usize) -> NsfDgLimits {
    NsfDgLimits {
        substeps: 1,
        maximum_cells: cells,
        maximum_cell_steps: 64 * cells,
        maximum_evaluations: 64,
        maximum_steps: 256,
        incoming_bracket_m_s: [-0.05, 0.05],
        pressure_absolute_tolerance_pa: 1e-7,
        pressure_relative_tolerance: 1e-10,
        // These checks separately inspect the covariance. This permissive
        // allowance is not an accuracy or transparent-boundary certificate.
        projection_absolute_tolerance_j: 0.0,
        projection_relative_tolerance: 1.0,
    }
}

pub(super) fn grid(dims: [usize; 3], axis: usize, pulse_to_high: Option<bool>) -> NsfGrid3<D3q27> {
    let gas = GasSpec::dry_air_ussa1976();
    let g = GasState::try_new(&gas, 300.0, 100000.0).unwrap();
    let mut periodic = [true; 3];
    periodic[axis] = false;
    let face = NsfCharacteristicBoundary {
        reference_pressure_pa: g.pressure,
        reference_temperature_k: g.temperature,
        reference_velocity_m_s: [0.0; 3],
        incoming_characteristic_m_s: 0.0,
    };
    NsfGrid3::<D3q27>::new(
        NsfSpec {
            dims,
            periodic,
            gas,
            temperature_k: g.temperature,
            pressure_pa: g.pressure,
            dx_m: 5e-5,
            temperature_window_k: [250.0, 600.0],
            maximum_mach: 0.2,
            maximum_temperature_ratio: 1.0,
            bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "authored causal-port numerical fixture; not gas validation",
            maximum_cells: dims.iter().product(),
            wall: WallThermal::Adiabatic,
            workers: 1,
            transport_points: 17,
        },
        |_| false,
        |x| {
            let mut velocity = [0.0; 3];
            let amplitude = pulse_to_high.map_or(0.0, |high| {
                let distance = if high {
                    x[axis]
                } else {
                    dims[axis] - 1 - x[axis]
                } as f64;
                let z = (distance - 0.3 * dims[axis] as f64) / (0.1 * dims[axis] as f64);
                let a = 1e-5 * fs_math::det::exp(-z * z);
                velocity[axis] =
                    if high { 1.0 } else { -1.0 } * (g.gamma / 3.0).sqrt() * a / g.gamma;
                a
            });
            (
                1.0 + amplitude / g.gamma,
                velocity,
                g.temperature * (1.0 + (g.gamma - 1.0) / g.gamma * amplitude),
            )
        },
    )
    .unwrap()
    .with_characteristic_axis(NsfCharacteristicAxis {
        axis,
        faces: [face; 2],
        maximum_relative_perturbation: 0.05,
    })
    .unwrap()
}

fn resistance(value: f64) -> RelaxationImpedanceSpec {
    RelaxationImpedanceSpec::new(
        SeriesImpedanceSpec {
            resistance_pa_s_m3: value,
            inertance_pa_s2_m3: 0.0,
            compliance_m3_pa: None,
        },
        &[],
    )
    .unwrap()
}

fn reactive(port: f64, dt: f64) -> RelaxationImpedanceSpec {
    RelaxationImpedanceSpec::new(
        SeriesImpedanceSpec {
            resistance_pa_s_m3: 0.2 * port,
            inertance_pa_s2_m3: 4.0 * port * dt,
            compliance_m3_pa: Some(4.0 * dt / port),
        },
        &[RelaxationTerm {
            resistance_pa_s_m3: 0.3 * port,
            rate_per_s: 1.0 / (3.0 * dt),
        }],
    )
    .unwrap()
    .with_compliance_terms(&[ComplianceRelaxationTerm {
        resistance_pa_s_m3: 0.4 * port,
        rate_per_s: 1.0 / (5.0 * dt),
    }])
    .unwrap()
}

fn port(grid: &NsfGrid3<D3q27>, high: bool) -> f64 {
    let g = grid.characteristic_reference_gas(high).unwrap();
    let traces = grid.open_boundary_traces(high).unwrap();
    g.density * g.sound_speed / (traces.len() as f64 * traces[0].area_m2)
}

fn snapshot(owner: &NsfImpedance<D3q27>) -> (u64, usize, u64, Vec<f64>) {
    let load = owner.load();
    let state = load.base_state();
    let mut values = vec![
        state.inertive_flow_m3_s,
        state.compliance_pressure_pa,
        load.stored_energy_j(),
    ];
    values.extend_from_slice(load.branch_flows_m3_s());
    values.extend_from_slice(load.branch_pressures_pa());
    (
        owner.zone().state_hash(),
        owner.zone().steps(),
        owner.accepted_steps(),
        values,
    )
}

#[test]
fn g1_matched_load_replays_native_outgoing_waves_on_both_sides_and_every_axis() {
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        for high in [false, true] {
            let mut dims = [1; 3];
            dims[axis] = 16;
            let original = grid(dims, axis, Some(high));
            let z = port(&original, high);
            let mut baseline = grid(dims, axis, Some(high));
            let mut owner =
                NsfImpedance::new(original, resistance(z), high, limits(16), &gate).unwrap();
            for _ in 0..16 {
                let initial_mass = owner.zone().total_mass_kg();
                let initial_energy = owner.zone().total_energy_j();
                let frame = owner
                    .step(NsfBoundaryDrive::IncomingCharacteristic(0.0), 0.0, &gate)
                    .unwrap();
                baseline.set_incoming_characteristics([0.0; 2]).unwrap();
                baseline.step().unwrap();
                assert!(frame.incoming_characteristic_m_s.abs() < 1e-8);
                assert!(frame.pressure_residual_pa.abs() < 2e-7);
                assert!(
                    (frame.flow_mass_change_kg
                        + frame.faces.iter().map(|f| f.mass_out_kg).sum::<f64>())
                    .abs()
                        < 2e-10 * initial_mass
                );
                assert!(
                    (frame.flow_energy_change_j
                        + frame
                            .faces
                            .iter()
                            .map(|f| f.total_energy_out_j)
                            .sum::<f64>())
                    .abs()
                        < 2e-10 * initial_energy
                );
                assert_eq!(frame.load.port.stored_energy_j, 0.0);
                assert!(frame.load.port.dissipated_energy_j >= 0.0);
            }
            for index in 0..16 {
                let mut x = [0; 3];
                x[axis] = index;
                let (rho, u, t) = owner.zone().cell(x);
                let (other_rho, other_u, other_t) = baseline.cell(x);
                assert!((rho - other_rho).abs() < 1e-10);
                assert!((u[axis] - other_u[axis]).abs() < 1e-10);
                assert!((t - other_t).abs() < 1e-7);
            }
        }
    }
}

#[test]
fn g3_reactive_storage_returns_work_and_both_relaxation_histories_survive() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        let original = grid([16, 1, 1], 0, None);
        let z = port(&original, high);
        let dt = original.law().dt_s;
        let mut owner =
            NsfImpedance::new(original, reactive(z, dt), high, limits(16), &gate).unwrap();
        let (mut work, mut dissipated, mut returned, mut peak_storage) = (0.0, 0.0, false, 0.0_f64);
        for step in 0..96 {
            let rise = if step < 4 { 2.0 } else { 0.0 };
            let frame = owner
                .step(NsfBoundaryDrive::IncomingCharacteristic(0.0), rise, &gate)
                .unwrap();
            let acoustic = frame.load.port;
            work += acoustic.supplied_work_j;
            dissipated += acoustic.dissipated_energy_j;
            peak_storage = peak_storage.max(acoustic.stored_energy_j);
            returned |= step >= 4 && acoustic.supplied_work_j < 0.0;
            let scale = acoustic.stored_energy_j
                + acoustic.storage_change_j.abs()
                + acoustic.dissipated_energy_j
                + acoustic.supplied_work_j.abs();
            assert!(frame.load_energy_residual_j.abs() < 1e-10 * scale.max(1e-300));
            assert!(acoustic.dissipated_energy_j >= 0.0);
            let expected_source = -rise * frame.load_exchange.normal_volume_flow_m3_s * dt;
            assert!(
                (frame.source_work_j - expected_source).abs()
                    < 1e-14 * (frame.source_work_j.abs() + expected_source.abs()).max(1e-300)
            );
            let decomposition = frame.projection_covariance_j
                + frame.pressure_residual_pa * frame.load_exchange.normal_volume_flow_m3_s * dt
                + acoustic.pressure_pa * frame.flow_residual_m3_s * dt;
            assert!(
                (frame.acoustic_work_residual_j - decomposition).abs()
                    < 1e-10
                        * (frame.load_exchange.acoustic_pressure_work_out_j.abs()
                            + frame.source_work_j.abs()
                            + acoustic.supplied_work_j.abs())
                        .max(1e-300)
            );
            if step == 3 {
                assert!(owner.load().base_state().inertive_flow_m3_s != 0.0);
                assert!(owner.load().base_state().compliance_pressure_pa != 0.0);
                assert!(owner.load().branch_flows_m3_s()[0] != 0.0);
                assert!(owner.load().branch_pressures_pa()[0] != 0.0);
            }
        }
        assert!(peak_storage > 0.0);
        assert!(
            returned,
            "a reactive load must be able to return stored energy"
        );
        assert!(
            (owner.load().stored_energy_j() + dissipated - work).abs()
                < 1e-9 * (peak_storage + dissipated + work.abs())
        );
        let accepted = snapshot(&owner);
        let (fluid, load) = owner.into_parts();
        assert_eq!(
            (fluid.state_hash(), fluid.steps()),
            (accepted.0, accepted.1)
        );
        assert_eq!(load.stored_energy_j(), accepted.3[2]);
    }
}

#[test]
fn g1_independent_apertures_keep_original_area_backgrounds_and_load_memory() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        let background = [
            NsfBoundaryDrive::PressureOffset(-0.25),
            NsfBoundaryDrive::PressureOffset(0.125),
        ];
        let patch = NsfCharacteristicPatch {
            start: [0, 0],
            end: [1, 1],
        };
        let mut original = grid([12, 2, 2], 0, None)
            .with_characteristic_patch(false, patch, NsfBoundaryDrive::IncomingCharacteristic(0.0))
            .unwrap()
            .with_characteristic_patch(true, patch, NsfBoundaryDrive::IncomingCharacteristic(0.0))
            .unwrap();
        original.set_boundary_drives(background).unwrap();
        let g = original.characteristic_reference_gas(high).unwrap();
        let area = original.law().dx_m.powi(2);
        let z = g.density * g.sound_speed / area;
        let dt = original.law().dt_s;
        let mut owner =
            NsfImpedance::new_on_patches(original, reactive(z, dt), high, limits(48), &gate)
                .unwrap();
        assert_eq!(owner.area_m2(), area);
        assert_eq!(owner.port_impedance_pa_s_m3(), z);
        for _ in 0..4 {
            let frame = owner
                .step(NsfBoundaryDrive::PressureOffset(0.5), 0.0, &gate)
                .unwrap();
            assert_eq!(frame.coupled_patch, Some(frame.load_exchange));
            assert!(frame.source_patch.is_some());
            for side in 0..2 {
                assert_eq!(
                    owner.patch_background_drive_on(side != 0),
                    Some(background[side])
                );
                for trace in owner.zone().open_boundary_traces(side != 0).unwrap() {
                    if trace.cell[1] != 0 || trace.cell[2] != 0 {
                        let p = if side == 0 { 99999.75 } else { 100000.125 };
                        assert!((trace.pressure_pa - p).abs() < 1e-8);
                    }
                }
            }
        }
    }
}

#[test]
fn g4_failed_roots_cancellation_and_exhausted_budgets_publish_no_memory() {
    let gate = CancelGate::new_clock_free();
    let original = grid([12, 1, 1], 0, None);
    let z = port(&original, true);
    let dt = original.law().dt_s;
    let mut budget = limits(12);
    budget.incoming_bracket_m_s = [-1e-12, 1e-12];
    let mut owner = NsfImpedance::new(original, reactive(z, dt), true, budget, &gate).unwrap();
    let before = snapshot(&owner);
    assert!(
        owner
            .step(NsfBoundaryDrive::IncomingCharacteristic(0.0), 2.0, &gate)
            .is_err()
    );
    assert_eq!(snapshot(&owner), before);
    assert!(
        owner
            .step(
                NsfBoundaryDrive::IncomingCharacteristic(0.0),
                f64::NAN,
                &gate
            )
            .is_err()
    );
    assert_eq!(snapshot(&owner), before);
    gate.request();
    assert!(matches!(
        owner.step(NsfBoundaryDrive::IncomingCharacteristic(0.0), 0.0, &gate),
        Err(NsfDgError::Cancelled)
    ));
    assert_eq!(snapshot(&owner), before);

    let gate = CancelGate::new_clock_free();
    let original = grid([12, 1, 1], 0, None);
    let mut budget = limits(12);
    budget.maximum_steps = 1;
    let mut owner = NsfImpedance::new(original, reactive(z, dt), true, budget, &gate).unwrap();
    let frame = owner
        .step(NsfBoundaryDrive::IncomingCharacteristic(0.0), 1.0, &gate)
        .unwrap();
    assert!(frame.evaluations <= budget.maximum_evaluations);
    assert_eq!(owner.zone().steps(), budget.substeps);
    let accepted = snapshot(&owner);
    assert!(
        owner
            .step(NsfBoundaryDrive::IncomingCharacteristic(0.0), 0.0, &gate)
            .is_err()
    );
    assert_eq!(snapshot(&owner), accepted);
}

#[test]
fn g0_wave_transform_reuses_original_load_memory_and_counts_the_source_once() {
    for z in [1.0, 417.0, 2e8] {
        let load = RelaxationImpedance::new(reactive(z, 1e-4), z, 1e-4).unwrap();
        for incident in [-0.5, 0.0, 2.0] {
            for source in [-1.0, 0.0, 0.25] {
                let native = load.preview_step(incident).unwrap();
                let pressure = native.port.pressure_pa + source;
                let flow = native.port.flow_m3_s;
                let (candidate, actual) = candidate_load(load, z, pressure, flow, source).unwrap();
                let scale = incident.abs().max(1.0);
                assert!((actual.port.pressure_pa - native.port.pressure_pa).abs() < 1e-12 * scale);
                assert!((actual.port.flow_m3_s - native.port.flow_m3_s).abs() < 1e-12 * scale / z);
                assert_eq!(load.stored_energy_j(), 0.0);
                let stored = candidate.stored_energy_j();
                assert!(
                    (stored - actual.port.stored_energy_j).abs()
                        < 8.0
                            * f64::EPSILON
                            * (stored.abs() + actual.port.stored_energy_j.abs()).max(1e-300)
                );
                assert!(
                    (pressure * flow * 1e-4 - source * flow * 1e-4 - actual.port.supplied_work_j)
                        .abs()
                        < 1e-11
                            * (actual.port.supplied_work_j.abs()
                                + source.abs() * flow.abs() * 1e-4)
                                .max(1e-300)
                );
            }
        }
    }
}
