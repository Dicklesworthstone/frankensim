//! Transaction and quadrature checks using the unchanged original NSF step.
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{
    NsfCharacteristicAxis, NsfCharacteristicBoundary, NsfGrid3, NsfSpec, WallThermal,
};
use fs_material::gas::{GasSpec, GasState};
use std::cell::Cell;

fn grid(reference_speed_m_s: f64) -> NsfGrid3<D3q27> {
    let gas = GasSpec::dry_air_ussa1976();
    let state = GasState::try_new(&gas, 300.0, 100000.0).unwrap();
    let velocity_unit = (3.0 * state.specific_gas_constant * state.temperature).sqrt();
    let face = NsfCharacteristicBoundary {
        reference_pressure_pa: state.pressure, reference_temperature_k: state.temperature,
        reference_velocity_m_s: [reference_speed_m_s, 0.0, 0.0],
        incoming_characteristic_m_s: 0.0,
    };
    NsfGrid3::<D3q27>::new(NsfSpec {
        dims: [12, 2, 2], periodic: [false, true, true], gas,
        temperature_k: 300.0, pressure_pa: 100000.0, dx_m: 5e-5,
        temperature_window_k: [250.0, 600.0], maximum_mach: 0.2,
        maximum_temperature_ratio: 1.0, bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012 N2 ratio; numerical transaction fixture",
        maximum_cells: 48, wall: WallThermal::Adiabatic, workers: 1, transport_points: 17,
    }, |_| false, |_| (1.0, [reference_speed_m_s / velocity_unit, 0.0, 0.0], 300.0))
        .unwrap().with_characteristic_axis(NsfCharacteristicAxis {
            axis: 0, faces: [face; 2], maximum_relative_perturbation: 0.05,
        }).unwrap()
}

#[test]
fn interval_accept_retains_exact_original_populations_charts_and_accounts() {
    let mut previewed = grid(0.0);
    let mut direct = grid(0.0);
    direct.set_incoming_characteristics([-0.1, 0.03]).unwrap();
    let mut pressure = [0.0; 2];
    let mut flow = [0.0; 2];
    let mut work = [0.0; 2];
    let mut mass = [0.0; 2];
    for _ in 0..3 {
        direct.step().unwrap();
        for (i, face) in direct.last_boundary_exchange().unwrap().iter().enumerate() {
            pressure[i] += face.mean_acoustic_pressure_pa;
            flow[i] += face.normal_volume_flow_m3_s;
            work[i] += face.acoustic_pressure_work_out_j;
            mass[i] += face.mass_out_kg;
        }
    }
    let trial = previewed.preview_characteristic_interval([-0.1, 0.03], 3, 144, || false).unwrap();
    assert_eq!(trial.grid().state_hash(), direct.state_hash());
    for (i, face) in trial.boundary_exchange().iter().enumerate() {
        assert_eq!(face.mean_acoustic_pressure_pa, pressure[i] / 3.0);
        assert_eq!(face.normal_volume_flow_m3_s, flow[i] / 3.0);
        assert_eq!(face.acoustic_pressure_work_out_j, work[i]);
        assert_eq!(face.mass_out_kg, mass[i]);
    }
    trial.accept();
    assert_eq!(previewed.steps(), 3);
    assert_eq!(previewed.cumulative_ledger(), direct.cumulative_ledger());
    assert_eq!(previewed.last_boundary_exchange(), direct.last_boundary_exchange());
    for x in 0..12 { assert_eq!(previewed.cell([x, 0, 0]), direct.cell([x, 0, 0])); }
}

#[test]
fn discarded_cancelled_and_refused_intervals_leave_every_accepted_history_untouched() {
    let mut zone = grid(0.0);
    zone.set_incoming_characteristics([-0.03, 0.0]).unwrap();
    zone.step().unwrap();
    let before = (zone.state_hash(), zone.steps(), zone.cumulative_ledger(),
        *zone.characteristic_axis().unwrap(), *zone.last_boundary_exchange().unwrap(),
        zone.open_boundary_traces(false).unwrap(), zone.cell([0, 0, 0]));
    drop(zone.preview_characteristic_interval([-0.1, 0.02], 3, 144, || false).unwrap());
    let calls = Cell::new(0);
    assert!(zone.preview_characteristic_interval([-0.1, 0.02], 3, 144, || {
        let n = calls.get(); calls.set(n + 1); n >= 2
    }).is_err()); // One private step has executed before cancellation.
    assert!(calls.get() >= 3);
    assert!(zone.preview_characteristic_interval([0.0; 2], 3, 143, || false).is_err());
    assert!(zone.preview_characteristic_interval([f64::NAN, 0.0], 3, 144, || false).is_err());
    assert_eq!(before, (zone.state_hash(), zone.steps(), zone.cumulative_ledger(),
        *zone.characteristic_axis().unwrap(), *zone.last_boundary_exchange().unwrap(),
        zone.open_boundary_traces(false).unwrap(), zone.cell([0, 0, 0])));
}

#[test]
fn relative_face_work_is_accumulated_before_large_static_pressure_work() {
    let mut zone = grid(5.0);
    let reference = zone.characteristic_reference_gas(true).unwrap();
    let trial = zone.preview_characteristic_interval([0.0; 2], 4, 192, || false).unwrap();
    for face in trial.boundary_exchange() {
        assert!(face.pressure_work_out_j.abs() > 0.0);
        assert!(face.mean_acoustic_pressure_pa.abs() < 1e-9 * reference.pressure);
        assert!(face.acoustic_pressure_work_out_j.abs()
            < 1e-9 * face.pressure_work_out_j.abs());
    }
}
