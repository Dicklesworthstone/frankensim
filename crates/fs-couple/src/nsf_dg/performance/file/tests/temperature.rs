//! G0/G3: physical source temperature reaches the original coupled flow owner.
use super::*;

fn thermal_input(temperature_k: f64) -> String {
    fixture(reference(100000.0), "dry-air 300 100000").replace(
        &format!("boundary 1\nhold {SAMPLES} -0.1"),
        &format!("boundary pressure-temperature 1\nhold {SAMPLES} 20 {temperature_k}"),
    )
}

#[test]
fn g3_held_source_temperature_reaches_flow_and_retains_callback_identity() {
    let text = thermal_input(300.5);
    let file = read(&text, 8);
    let control = NsfDgDrive::Boundary(NsfBoundaryDrive::PressureTemperature {
        pressure_offset_pa: 20.0, inflow_temperature_k: 300.5,
    });
    assert!(file.renderer().drive_history().iter().all(|&value| value == control));
    let (whole, warm) = render(file, 8);
    let (partitioned, replay) = render(read(&text, 3), 3);
    assert_eq!(whole, partitioned);
    assert_eq!(warm.system().zone().state_hash(), replay.system().zone().state_hash());
    let (_, cooler) = render(read(&thermal_input(300.1), 8), 8);
    let warm_cell = warm.system().zone().cell([0, 1, 1]);
    let cool_cell = cooler.system().zone().cell([0, 1, 1]);
    assert!(warm_cell.2 > cool_cell.2, "source temperature must affect the actual fluid");
    assert_ne!(warm.system().zone().state_hash(), cooler.system().zone().state_hash());
    assert_eq!(warm.system().reference_gas(), cooler.system().reference_gas());
    assert_eq!(warm.system().load().duct().inlet_gas(), cooler.system().load().duct().inlet_gas());
    assert_eq!(warm.system().zone().steps(), SAMPLES);
}

#[test]
fn g0_temperature_history_preserves_physical_segment_endpoints_and_rejects_bad_units() {
    let text = thermal_input(300.5).replace(
        &format!("boundary pressure-temperature 1\nhold {SAMPLES} 20 300.5"),
        &format!("boundary pressure-temperature 2\nhold 12 20 300.5\nhold {SAMPLES} 10 300.2"),
    );
    let file = read(&text, 8);
    let history = file.renderer().drive_history();
    assert_eq!(history[11], NsfDgDrive::Boundary(NsfBoundaryDrive::PressureTemperature {
        pressure_offset_pa: 20.0, inflow_temperature_k: 300.5,
    }));
    assert_eq!(history[12], NsfDgDrive::Boundary(NsfBoundaryDrive::PressureTemperature {
        pressure_offset_pa: 10.0, inflow_temperature_k: 300.2,
    }));
    for record in ["20", "20 0", "20 249", "20 601", "20 NaN", "20 300 extra"] {
        let invalid = thermal_input(300.5).replace("20 300.5", record);
        assert!(ResolvedDgFilePerformance::from_bytes(invalid.as_bytes(), 8,
            Arc::new(CancelGate::new_clock_free())).is_err(), "accepted {record}");
    }
    let ambiguous = thermal_input(300.5).replace("boundary pressure-temperature", "boundary pressure");
    assert!(ResolvedDgFilePerformance::from_bytes(ambiguous.as_bytes(), 8,
        Arc::new(CancelGate::new_clock_free())).is_err());
}

#[test]
fn g0_canonical_warm_source_retains_explicit_pressure_and_temperature() {
    let file = read(include_str!("../../../../../examples/resolved-dg-warm-source.performance"), 8);
    assert_eq!(file.info().samples, 2048);
    assert_eq!(file.renderer().drive_history()[0],
        NsfDgDrive::Boundary(NsfBoundaryDrive::PressureTemperature {
            pressure_offset_pa: 20.0, inflow_temperature_k: 303.0,
        }));
    assert!(file.renderer().observer_gas().is_some());
}
