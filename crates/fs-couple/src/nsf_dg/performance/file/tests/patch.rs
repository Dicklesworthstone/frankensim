//! Native declarations reach the original split-aperture flow/duct owners.
//! These are numerical integration fixtures, not measured flute geometry.
use super::*;

fn input(background_pa: f64) -> String {
    let gas = reference(100000.0);
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    let whole_radius = (9.0 * dx * dx / core::f64::consts::PI).sqrt();
    let patch_radius = (dx * dx / core::f64::consts::PI).sqrt();
    fixture(gas, "dry-air 300 100000")
        .replace(
            "characteristic 0 high 0.05",
            &format!("characteristic 0 high 0.05\ncoupled_patch 1 2 1 2 pressure {background_pa}"),
        )
        .replace(
            &format!("dg_tube 0.1 {whole_radius:.17e}"),
            &format!("dg_tube 0.1 {patch_radius:.17e}"),
        )
        .replace("boundary 1", "boundary pressure 1")
        .replace("hold 32 -0.1", "hold 32 2")
}

#[test]
fn native_aperture_keeps_its_geometry_background_and_physical_playback() {
    let text = input(-0.25);
    let file = read(&text, 8);
    let system = file.renderer().system();
    assert_eq!(
        system.coupled_patch(),
        Some(NsfCharacteristicPatch {
            start: [1, 1],
            end: [2, 2],
        })
    );
    assert_eq!(
        system.patch_background_drive(),
        Some(NsfBoundaryDrive::PressureOffset(-0.25))
    );
    assert_eq!(system.area_m2(), system.zone().law().dx_m.powi(2));
    let (whole, complete) = render(file, 8);
    let (parts, partitioned) = render(read(&text, 3), 3);
    assert_eq!(whole, parts);
    assert!(whole.iter().any(|p| p.abs() > 1e-8));
    assert_eq!(
        complete.system().zone().state_hash(),
        partitioned.system().zone().state_hash()
    );
    let frame = complete.last_frame().unwrap();
    let patch = frame.coupled_patch.unwrap();
    assert_eq!(frame.inlet_flow_m3_s, patch.normal_volume_flow_m3_s);
    assert_ne!(
        frame.inlet_flow_m3_s,
        frame.faces[1].normal_volume_flow_m3_s
    );
    assert!((frame.inlet_pressure_pa - patch.mean_acoustic_pressure_pa).abs() < 1e-6);
    for trace in complete.system().zone().open_boundary_traces(true).unwrap() {
        if trace.cell[1] != 1 || trace.cell[2] != 1 {
            assert!((trace.pressure_pa - 99999.75).abs() < 1e-7);
        }
    }
    let (ambient, altered) = render(read(&input(0.0), 8), 8);
    assert_ne!(
        whole, ambient,
        "surrounding pressure must alter the original coupled motion"
    );
    assert_ne!(
        complete.system().zone().state_hash(),
        altered.system().zone().state_hash()
    );
}

#[test]
fn native_aperture_refuses_invalid_bounds_controls_and_inlet_area() {
    let text = input(0.0);
    let gas = reference(100000.0);
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    let radius = (dx * dx / core::f64::consts::PI).sqrt();
    for invalid in [
        text.replace(
            "coupled_patch 1 2 1 2 pressure 0",
            "coupled_patch 2 1 1 2 pressure 0",
        ),
        text.replace(
            "coupled_patch 1 2 1 2 pressure 0",
            "coupled_patch 1 4 1 2 pressure 0",
        ),
        text.replace(
            "coupled_patch 1 2 1 2 pressure 0",
            "coupled_patch 1 2 1 2 pressure NaN",
        ),
        text.replace(
            "coupled_patch 1 2 1 2 pressure 0",
            "coupled_patch 1 2 1 2 inferred-ambient",
        ),
        text.replace(
            &format!("dg_tube 0.1 {radius:.17e}"),
            &format!("dg_tube 0.1 {:.17e}", 2.0 * radius),
        ),
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
    let characteristic = text.replace(
        "coupled_patch 1 2 1 2 pressure 0",
        "coupled_patch 1 2 1 2 characteristic 0.001",
    );
    assert_eq!(
        read(&characteristic, 8)
            .renderer()
            .system()
            .patch_background_drive(),
        Some(NsfBoundaryDrive::IncomingCharacteristic(0.001))
    );
}

#[test]
fn canonical_split_aperture_reaches_the_original_exterior_observer() {
    let file = read(
        include_str!("../../../../../examples/resolved-dg-split-aperture.performance"),
        8,
    );
    assert!(file.renderer().system().coupled_patch().is_some());
    assert!(file.renderer().baffled_receiver().is_some());
    assert!(file.renderer().observer_gas().is_some());
    let mut renderer = file.into_renderer();
    renderer.block(&mut [0.0; 8]).unwrap();
    assert!(renderer.last_frame().unwrap().coupled_patch.is_some());
    assert_eq!(renderer.system().accepted_steps(), 8);
}

#[test]
fn both_native_apertures_share_one_flow_owner_and_keep_the_selected_bore_flux() {
    let text = input(-0.25).replace(
        "initial rest",
        &format!(
            "reservoir_aperture 1 high {:.17e} 0 0 0 300 1\nreservoir_cell 5 2 1\ninitial rest",
            reference(100000.0).density
        ),
    );
    let mut renderer = read(&text, 8).into_renderer();
    let initial_mass = renderer.system().zone().total_mass_kg();
    for _ in 0..4 {
        renderer.block(&mut [0.0]).unwrap();
        let frame = renderer.last_frame().unwrap();
        assert_eq!(
            frame.inlet_flow_m3_s,
            frame.coupled_patch.unwrap().normal_volume_flow_m3_s
        );
        assert!(renderer.system().zone().reservoir_aperture().is_some());
        let balance = frame.flow_mass_change_kg
            + frame.reservoir_mass_out_kg
            + frame.faces.iter().map(|face| face.mass_out_kg).sum::<f64>();
        assert!(balance.abs() / initial_mass < 2e-11);
    }
}

#[test]
fn independent_native_source_and_bore_patches_replay_si_geometry_at_any_callback_size() {
    let gas = reference(100000.0);
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    for high in [false, true] {
        let text = input(-0.25)
            .replace("characteristic 0 high", if high { "characteristic 0 high" } else { "characteristic 0 low" })
            .replace("initial rest", "source_patch 0 2 0 1 pressure 0.375\ninitial rest");
        let physical = text
            .replace("coupled_patch 1 2 1 2 pressure -0.25",
                &format!("coupled_patch_m {dx:.17e} {:.17e} {dx:.17e} {:.17e} pressure -0.25", 2.0*dx, 2.0*dx))
            .replace("source_patch 0 2 0 1 pressure 0.375",
                &format!("source_patch_m 0 {:.17e} 0 {dx:.17e} pressure 0.375", 2.0*dx));
        let file = read(&text, 8);
        assert_eq!(file.renderer().system().coupled_patch(), Some(NsfCharacteristicPatch { start: [1, 1], end: [2, 2] }));
        assert_eq!(file.renderer().system().source_patch(), Some(NsfCharacteristicPatch { start: [0, 0], end: [2, 1] }));
        let (whole, a) = render(file, 8);
        let (parts, b) = render(read(&physical, 3), 3);
        assert_eq!(whole, parts);
        assert_eq!(a.system().zone().state_hash(), b.system().zone().state_hash());
        assert_eq!(a.system().load().state(), b.system().load().state());
        let frame = a.last_frame().unwrap();
        assert_eq!(frame.inlet_flow_m3_s, frame.coupled_patch.unwrap().normal_volume_flow_m3_s);
        assert_ne!(frame.inlet_flow_m3_s, frame.faces[usize::from(high)].normal_volume_flow_m3_s);
        assert_ne!(frame.source_patch.unwrap().normal_volume_flow_m3_s,
            frame.faces[usize::from(!high)].normal_volume_flow_m3_s);
        assert_eq!(a.system().patch_background_drive_on(high), Some(NsfBoundaryDrive::PressureOffset(-0.25)));
        assert_eq!(a.system().patch_background_drive_on(!high), Some(NsfBoundaryDrive::PressureOffset(0.375)));
        let reversed = text.replace(
            "coupled_patch 1 2 1 2 pressure -0.25\nsource_patch 0 2 0 1 pressure 0.375",
            "source_patch 0 2 0 1 pressure 0.375\ncoupled_patch 1 2 1 2 pressure -0.25");
        let (reordered, c) = render(read(&reversed, 8), 8);
        assert_eq!(whole, reordered);
        assert_eq!(a.system().zone().state_hash(), c.system().zone().state_hash());
    }
}

#[test]
fn native_dual_apertures_refuse_repeated_roles_including_mixed_units() {
    for duplicate in ["coupled_patch 0 1 0 1 pressure 0",
        "source_patch 0 1 0 1 pressure 0", "source_patch_m 0 0.001 0 0.001 pressure 0"]
    {
        let text = input(-0.25).replace("initial rest",
            &format!("source_patch 0 2 0 1 pressure 0.375\n{duplicate}\ninitial rest"));
        assert!(ResolvedDgFilePerformance::from_bytes(text.as_bytes(), 8,
            Arc::new(CancelGate::new_clock_free())).is_err());
    }
}
