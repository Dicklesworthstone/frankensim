use super::*;

fn pressure_input() -> String {
    fixture(reference(100000.0), "dry-air 300 100000")
        .replace("boundary 1", "boundary pressure 1")
        .replace("hold 32 -0.1", "hold 32 20")
}

#[test]
fn physical_pressure_records_reach_original_faces_and_keep_callback_identity() {
    let text = pressure_input();
    let file = read(&text, 8);
    assert!(
        file.renderer()
            .drive_history()
            .iter()
            .all(|&x| x == NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(20.0)))
    );
    let (whole, solved) = render(file, 8);
    let (parts, partitioned) = render(read(&text, 3), 3);
    assert_eq!(whole, parts);
    assert_eq!(
        solved.system().zone().state_hash(),
        partitioned.system().zone().state_hash()
    );
    assert_eq!(
        solved.system().load().state(),
        partitioned.system().load().state()
    );
    let source = solved.last_frame().unwrap().faces[0];
    assert!((source.mean_acoustic_pressure_pa - 20.0).abs() < 1e-9);
    assert!(source.normal_volume_flow_m3_s < 0.0);
    let work = 20.0 * source.normal_volume_flow_m3_s * solved.system().load().propagator().dt();
    assert!((source.acoustic_pressure_work_out_j - work).abs() < 1e-12 * work.abs());

    let original = fixture(reference(100000.0), "dry-air 300 100000");
    let explicit = original.replace("boundary 1", "boundary characteristic 1");
    let (a, a_state) = render(read(&original, 8), 8);
    let (b, b_state) = render(read(&explicit, 3), 3);
    assert_eq!(a, b);
    assert_eq!(
        a_state.system().zone().state_hash(),
        b_state.system().zone().state_hash()
    );
    assert_ne!(
        whole, a,
        "a pressure is not the same physical input as a characteristic velocity"
    );
    for invalid in [
        text.replace("boundary pressure 1", "boundary velocity 1"),
        text.replace("boundary pressure 1", "boundary pressure 1 extra"),
        text.replace("hold 32 20", "hold 32 NaN"),
    ] {
        assert!(
            ResolvedDgFilePerformance::from_bytes(
                invalid.as_bytes(),
                8,
                Arc::new(CancelGate::new_clock_free())
            )
            .is_err()
        );
    }
}

fn sourced_pressure(
    cp: f64,
    beta: f64,
    molecular: Option<&str>,
    block: usize,
) -> Result<ResolvedDgFilePerformance, ResolvedDgInputError> {
    let a = support::packs("synthetic-alpha", 0.028, cp, 200.0);
    let b = support::packs("synthetic-beta", 0.044, 4.0, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let declaration = format!(
        "thermochemical 300 100000 mole 2 synthetic-alpha {} synthetic-beta {beta}",
        1.0 - beta
    );
    let mut text = fixture(reference(100000.0), &declaration)
        .replace("boundary 1", "boundary pressure 1")
        .replace("hold 32 -0.1", "hold 32 20");
    if let Some(law) = molecular {
        text = text.replace(
            "observation inlet-pressure",
            &format!("dg_molecular {law}\nobservation inlet-pressure"),
        );
    }
    ResolvedDgFilePerformance::from_bytes_with_sources(
        text.as_bytes(),
        block,
        Arc::new(CancelGate::new_clock_free()),
        Some(&sources),
        |_| Err("this fixture has no material wall".into()),
    )
}

#[test]
fn data_only_composition_and_calorics_change_flow_at_the_same_blowing_pressure() {
    let baseline = sourced_pressure(3.5, 0.2, None, 8).unwrap();
    let cp = sourced_pressure(4.5, 0.2, None, 8).unwrap();
    let composition = sourced_pressure(3.5, 0.8, None, 8).unwrap();
    let dx = baseline.renderer().system().zone().law().dx_m;
    for file in [&baseline, &cp, &composition] {
        assert_eq!(file.renderer().system().zone().law().dx_m, dx);
        assert!(
            file.renderer()
                .drive_history()
                .iter()
                .all(|&x| x == NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(20.0)))
        );
    }
    let (a, solved) = render(baseline, 8);
    let (b, cp) = render(cp, 8);
    let (c, composition) = render(composition, 3);
    assert_ne!(a, b);
    assert_ne!(a, c);
    for result in [&solved, &cp, &composition] {
        assert!(
            (result.last_frame().unwrap().faces[0].mean_acoustic_pressure_pa - 20.0).abs() < 1e-9
        );
        assert_eq!(
            result
                .system()
                .zone()
                .thermochemical_state([3, 1, 1])
                .unwrap()
                .unwrap()
                .source_packs()
                .len(),
            2
        );
    }
    assert_ne!(
        solved.last_frame().unwrap().faces[0].normal_volume_flow_m3_s,
        composition.last_frame().unwrap().faces[0].normal_volume_flow_m3_s
    );
}

#[test]
fn native_resolved_bore_retains_original_species_molecular_memory_and_refusals() {
    // Authored constitutive data only. The NSF source keeps its original
    // caloric/transport law; this selects the separate DG region's causal loss.
    const LAW: &str =
        "declared 2 synthetic-alpha 1 vibration 1 1000 1 power 10 300 0 250 600 synthetic-beta 0";
    let plain = sourced_pressure(3.5, 0.2, None, 8).unwrap();
    let lossy = sourced_pressure(3.5, 0.2, Some(LAW), 8).unwrap();
    let repeated = sourced_pressure(3.5, 0.2, Some(LAW), 3).unwrap();
    let plain_dofs = plain.renderer().system().load().duct().total_dofs();
    let duct = lossy.renderer().system().load().duct();
    assert_eq!(duct.molecular_relaxation().unwrap().gas(), duct.inlet_gas());
    assert!(duct.total_dofs() > plain_dofs);
    let (a, _) = render(plain, 8);
    let (b, solved) = render(lossy, 8);
    let (c, partitioned) = render(repeated, 3);
    assert_ne!(a, b);
    assert_eq!(b, c);
    assert_eq!(
        solved.system().load().state(),
        partitioned.system().load().state()
    );
    assert_eq!(
        solved.system().zone().state_hash(),
        partitioned.system().zone().state_hash()
    );
    let missing_species = LAW
        .replace("declared 2", "declared 1")
        .replace(" synthetic-beta 0", "");
    assert!(sourced_pressure(3.5, 0.2, Some(&missing_species), 8).is_err());
    let unsupported_temperature = LAW.replace("250 600", "350 600");
    assert!(sourced_pressure(3.5, 0.2, Some(&unsupported_temperature), 8).is_err());
    let snapshot = pressure_input().replace(
        "observation inlet-pressure",
        &format!("dg_molecular {LAW}\nobservation inlet-pressure"),
    );
    assert!(
        ResolvedDgFilePerformance::from_bytes(
            snapshot.as_bytes(),
            8,
            Arc::new(CancelGate::new_clock_free())
        )
        .is_err()
    );
}

#[test]
fn canonical_pressure_source_keeps_the_original_exterior_observer() {
    let file = read(
        include_str!("../../../../../examples/resolved-dg-pressure.performance"),
        8,
    );
    assert_eq!(file.info().dimensions, [16, 3, 3]);
    assert_eq!(file.info().solid_cells, 4);
    assert!(matches!(
        file.renderer().observation(),
        ApertureObservation::TubeBaffled(_)
    ));
    assert!(file.renderer().baffled_receiver().is_some());
    assert!(file.renderer().observer_gas().is_some());
    let history = file.renderer().drive_history();
    assert_eq!(history[0], NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(20.0)));
    assert_eq!(history[128], NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(-20.0)));
    assert_eq!(history[256], NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(0.0)));
}

const OBSERVER_LAW: &str = "declared 2 N2 0 CO2 1 vibration 1 1000 1 power 10 300 0 250 600";

fn molecular_observer_input() -> String {
    // The fluid zone remains dry air. Only the independent exterior contains
    // this deliberately synthetic molecular process; it is not a measured
    // carbon-dioxide relaxation law or a flute validation fixture.
    include_str!("../../../../../examples/resolved-dg-pressure.performance")
        .replace("output 2048 20", "output 32 20")
        .replace("hold 128 20", "hold 8 20")
        .replace("hold 256 -20", "hold 16 -20")
        .replace("hold 2048 0", "hold 32 0")
        .replace(
            "dg_exterior 300 100000 0",
            "dg_exterior mixture 300 100000 mole 2 N2 0.9 CO2 0.1",
        )
        .replace(
            "observer_gas stokes-kirchhoff 4096 1000000 50000000",
            &format!("observer_gas molecular 4096 1000000 50000000 {OBSERVER_LAW}"),
        )
}

#[test]
fn resolved_molecular_observer_keeps_downstream_species_with_a_distinct_upstream_exterior() {
    let text = molecular_observer_input()
        .replace("coupling 64 9216", "coupling 128 18432")
        .replace("boundary pressure", "boundary upstream-pressure")
        .replace("observation baffled", concat!(
            "upstream_dg dg_tube 0.03 0.00896143253443043 baffled 1 2 32 50 2000 0.01 1e-11 256\n",
            "dg_radiation computed-baffled 20 65 4 8 12 20 0.1 67108864 2 100 1e-8 64 1048576 257\n",
            "dg_exterior mixture 300 100000 mole 2 N2 0.8 O2 0.2\n",
            "observation baffled"));
    let file = read(&text, 8);
    let system = file.renderer().system();
    let downstream = system.load().radiation_record().unwrap().gas_record();
    let upstream = system
        .upstream_load()
        .unwrap()
        .radiation_record()
        .unwrap()
        .gas_record();
    let report = file.renderer().observer_gas().unwrap();
    assert_eq!(&report.gas, downstream.state());
    assert_ne!(&report.gas, upstream.state());
    assert_eq!(file.gas_records().last().unwrap().state(), upstream.state());
}
