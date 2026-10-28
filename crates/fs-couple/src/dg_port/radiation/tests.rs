use super::*;
use fs_duct::dg::{DgMesh, DgSidePort};
use fs_duct::{Duct, Segment};
use fs_material::gas_species::input::read_gas_record;

const PORTS: [PortId; 2] = [PortId(7), PortId(1)];

fn exterior(text: &str) -> GasRecord {
    read_gas_record(&mut text.split_ascii_whitespace(), None, None).unwrap()
}

fn duct(gas: &GasRecord) -> DgDuct {
    DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.002, length: 0.05 }] },
        gas.state(), DgMesh { elements: 2, order: 2 },
        DgBoundary::Port(PortId(0)), DgBoundary::Port(PortId(1)),
    ).unwrap().with_side_ports(&[DgSidePort { id: PORTS[0], downstream_element: 1 }]).unwrap()
}

fn openings() -> [(PortId, CircularOpening); 2] {
    // Deliberately reverse the original DG [inlet, outlet, side] order. These
    // are the physical source owner's established two-piston bake geometry.
    [
        (PORTS[0], CircularOpening { center_m: [0.015, 0.0], radius_m: 0.0015 }),
        (PORTS[1], CircularOpening { center_m: [0.0, 0.0], radius_m: 0.002 }),
    ]
}

#[test]
fn source_binding_requires_all_original_ports_actual_axial_area_and_static_pressure() {
    let gas = exterior("293.15 101325 0");
    let duct = duct(&gas);
    let openings = openings();
    admit_openings(&duct, &[], &gas, &openings).unwrap();
    assert!(admit_openings(&duct, &[], &gas, &openings[..1]).is_err());
    for refused in [
        [(PortId(0), openings[0].1), openings[1]],
        [(PortId(99), openings[0].1), openings[1]],
        [openings[0], openings[0]],
        [openings[0], (PORTS[1], CircularOpening { radius_m: 0.001, ..openings[1].1 })],
        [openings[0], (PORTS[1], CircularOpening { radius_m: f64::INFINITY, ..openings[1].1 })],
    ] {
        assert!(admit_openings(&duct, &[], &gas, &refused).is_err());
    }
    assert!(admit_openings(&duct, &[], &exterior("293.15 90000 0"), &openings).is_err());
    let changed = exterior("mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2");
    admit_openings(&duct, &[], &changed, &openings).unwrap();
    assert_ne!(gas.state().sound_speed, changed.state().sound_speed);
}

#[test]
fn computed_mutual_load_and_microphone_keep_the_same_complete_ordered_flow_basis() {
    let gas = exterior("293.15 101325 0");
    let openings = openings();
    let controls = ComputedBaffledInput::parse(
        "20 65 4 8 12 20 0.01 67108864 2 100 1e-8 64 1048576 257",
    ).unwrap();
    let mut load = DgPortLoad::new_with_computed_baffled_ports(
        duct(&gas), 8000, 1e-11, 64, &gas, controls, &openings,
        &CancelGate::new_clock_free(),
    ).unwrap();
    assert_eq!(load.terminal_ports(), &PORTS);
    let source = load.radiation_record().unwrap();
    assert_eq!(source.gas_record(), &gas);
    assert_eq!(source.apertures(), &[openings[0].1, openings[1].1]);
    let matrix = load.terminal().unwrap().source_model().eval_iw(core::f64::consts::TAU * 1100.0);
    assert!(matrix[1].abs() > 0.005 * matrix[0].abs(), "mutual reaction must remain in the source");
    let config = CircularOutletReceiver {
        position_m: [0.1, 0.05, 1.0], radial_rings: 8, angular_points: 32,
        maximum_frequency_hz: 1000.0,
    };
    let mut combined = load.baffled_receiver(8000, config).unwrap();
    let mut individual: Vec<_> = openings.iter().map(|(_, aperture)|
        BaffledPressure::circular_openings_multirate_differentiating(
            &[*aperture], 8000, 8000, config, gas.state().into(),
        ).unwrap()
    ).collect();
    let mut energy = [0.0; 2];
    for sample in 0..512 {
        load.preview_step(if sample < 8 { 1e-9 } else { 0.0 }).unwrap().accept();
        let frames = load.last_terminal_frames();
        assert_eq!([frames[0].port, frames[1].port], PORTS);
        assert_eq!(load.last_terminal_frame(), Some(&frames[1]));
        let flows = [frames[0].outward_flow_m3_s, frames[1].outward_flow_m3_s];
        let pressure = combined.step_trace(&flows).unwrap();
        let a = individual[0].step_trace(&flows[..1]).unwrap();
        let b = individual[1].step_trace(&flows[1..]).unwrap();
        assert!((pressure - (a + b)).abs() <= 1e-12 * (a.abs() + b.abs()).max(1e-300));
        energy[0] += a * a;
        energy[1] += b * b;
    }
    assert!(energy.iter().all(|value| *value > 0.0), "both actual openings must contribute pressure");
}

#[path = "tests/chimneys.rs"]
mod chimneys;

#[test]
fn computed_radiation_preserves_fractional_physical_steps_and_integer_clock_bits() {
    let gas = exterior("293.15 101325 0");
    let openings = openings();
    let controls = ComputedBaffledInput::parse(
        "20 65 4 8 12 20 0.01 67108864 2 100 1e-8 64 1048576 257",
    ).unwrap();
    let gate = CancelGate::new_clock_free();
    let prepare = |dt| DgPortLoad::new_with_computed_baffled_branches_at_dt(
        duct(&gas), dt, 1e-11, 64, &gas, controls, &openings, Vec::new(), &gate,
    );
    let mut original = DgPortLoad::new_with_computed_baffled_branches(
        duct(&gas), 8001, 1e-11, 64, &gas, controls, &openings, Vec::new(), &gate,
    ).unwrap();
    // Reciprocal roundoff would incorrectly ceil this legacy clock to 8002 Hz.
    let mut same = prepare(1.0 / 8001.0).unwrap();
    assert_eq!(original.radiation_record().unwrap().declaration_identity(),
        same.radiation_record().unwrap().declaration_identity());
    assert_eq!(original.radiation_record().unwrap().realization_identity(),
        same.radiation_record().unwrap().realization_identity());
    for step in 0..16 {
        let flow = if step < 3 { 1e-9 } else { 0.0 };
        assert_eq!(original.preview_step(flow).unwrap().accept(),
            same.preview_step(flow).unwrap().accept());
        assert_eq!(original.last_terminal_frames(), same.last_terminal_frames());
        for (a, b) in original.state().pressure().iter().zip(same.state().pressure()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
    // A physical material clock need not be the reciprocal of a whole hertz.
    // The extra source-band margin must never enter either live propagator.
    let dt = 1.0 / 8000.25;
    let mut fractional = prepare(dt).unwrap();
    assert_eq!(fractional.propagator().dt().to_bits(), dt.to_bits());
    let terminal = fractional.terminal().unwrap();
    assert_eq!(terminal.time_step_s().to_bits(), dt.to_bits());
    assert_eq!(terminal.limits().band_rad_s[1].to_bits(),
        (core::f64::consts::PI / dt).to_bits());
    assert!(terminal.reports().1.passive);
    let source = fractional.radiation_record().unwrap();
    assert_eq!(source.band_hz(), [20.0, 4000.5]);
    assert_eq!(source.gas_record(), &gas);
    assert_eq!(source.apertures(), &[openings[0].1, openings[1].1]);
    assert_eq!(fractional.terminal_ports(), &PORTS);
    assert_eq!(fractional.duct().operator_identity(), same.duct().operator_identity());
    fractional.preview_step(1e-9).unwrap().accept();
    assert_eq!(fractional.state().step_count(), 1);
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, 1.0 / (f64::from(u32::MAX) + 1.0)] {
        assert!(prepare(dt).is_err());
    }
}
