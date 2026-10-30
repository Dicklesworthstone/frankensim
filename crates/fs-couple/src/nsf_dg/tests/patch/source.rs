//! A localized original source is independent of its surrounding opening.
use super::*;

fn source_zone(high: bool, substeps: usize) -> NsfGrid3<D3q27> {
    bare_zone(high, substeps)
        .with_characteristic_patch(!high, PATCH, NsfBoundaryDrive::PressureOffset(BACKGROUND_PA))
        .unwrap()
}

fn owner(high: bool, substeps: usize, upstream: bool) -> NsfDg<D3q27> {
    let gate = CancelGate::new_clock_free();
    let zone = source_zone(high, substeps);
    let main = load(&zone, reference(), FACE_CELLS as f64, 1.0, DgBoundary::Matched);
    let owner = NsfDg::new_with_source_patch(zone, main, high, budget(substeps), &gate).unwrap();
    if upstream {
        let source = load(owner.zone(), reference(), PATCH_CELLS as f64, 1.0, DgBoundary::Closed);
        owner.with_upstream_load(source, &gate).unwrap()
    } else { owner }
}

// Recompute original boundary quadrature independently of both coupling roots
// and the patch-exchange accumulator. No whole-face flow stands in for a patch.
fn advance_original(
    fluid: &mut NsfGrid3<D3q27>, high: bool, substeps: usize,
    incoming_main: f64, source_drive: NsfBoundaryDrive,
) -> (f64, f64, f64) {
    let mut drives = [NsfBoundaryDrive::PressureOffset(BACKGROUND_PA); 2];
    drives[usize::from(high)] = NsfBoundaryDrive::IncomingCharacteristic(incoming_main);
    fluid.set_boundary_drives_with_patch(drives, source_drive).unwrap();
    let (mut pressure, mut flow, mut work) = (0.0, 0.0, 0.0);
    let dt = fluid.law().dt_s;
    for _ in 0..substeps {
        let before = fluid.open_boundary_traces(!high).unwrap();
        fluid.step().unwrap();
        let after = fluid.open_boundary_traces(!high).unwrap();
        for (old, new) in before.iter().zip(&after) {
            assert_eq!(old.cell, new.cell);
            if (PATCH.start[0]..PATCH.end[0]).contains(&old.cell[1])
                && (PATCH.start[1]..PATCH.end[1]).contains(&old.cell[2])
            {
                pressure += 0.5 * (old.pressure_pa + new.pressure_pa);
                flow += 0.5 * (old.outward_normal_velocity_m_s + new.outward_normal_velocity_m_s)
                    * old.area_m2;
                work += 0.5 * ((old.pressure_pa - reference().pressure) * old.outward_normal_velocity_m_s
                    + (new.pressure_pa - reference().pressure) * new.outward_normal_velocity_m_s)
                    * old.area_m2 * dt;
            } else {
                assert_eq!(new.pressure_pa, reference().pressure + BACKGROUND_PA);
            }
        }
    }
    (pressure / (PATCH_CELLS * substeps) as f64, flow / substeps as f64, work)
}

fn check_patch(frame: &NsfDgFrame, expected: (f64, f64, f64), high: bool, dt: f64) -> f64 {
    assert!(frame.coupled_patch.is_none());
    let patch = frame.source_patch.unwrap();
    assert!((patch.mean_pressure_pa - expected.0).abs() < 1e-8);
    assert!((patch.normal_volume_flow_m3_s - expected.1).abs()
        < 1e-23 + 1e-11 * expected.1.abs());
    assert!((patch.acoustic_pressure_work_out_j - expected.2).abs()
        < 1e-25 + 1e-11 * expected.2.abs());
    assert_eq!(frame.inlet_flow_m3_s, frame.faces[usize::from(high)].normal_volume_flow_m3_s);
    let face = frame.faces[usize::from(!high)];
    let escaping = face.normal_volume_flow_m3_s - patch.normal_volume_flow_m3_s;
    let remainder_work = face.acoustic_pressure_work_out_j - patch.acoustic_pressure_work_out_j;
    let expected_work = BACKGROUND_PA * escaping * dt;
    let scale = face.acoustic_pressure_work_out_j.abs() + patch.acoustic_pressure_work_out_j.abs()
        + expected_work.abs();
    assert!((remainder_work - expected_work).abs() < 256.0 * f64::EPSILON * scale.max(1e-30));
    escaping
}

#[test]
fn source_aperture_controls_only_its_original_cells_and_preserves_main_face_work() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        for substeps in [1, 3] {
            let mut owner = owner(high, substeps, false);
            assert_eq!(owner.source_patch(), Some(PATCH));
            assert_eq!(owner.coupled_patch(), None);
            assert_eq!(owner.area_m2(), FACE_CELLS as f64 * owner.zone().law().dx_m.powi(2));
            let mut fluid = source_zone(high, substeps);
            let mut main = load(&fluid, reference(), FACE_CELLS as f64, 1.0, DgBoundary::Matched);
            let mut largest_escape: f64 = 0.0;
            for drive in [NsfBoundaryDrive::PressureOffset(2.0),
                NsfBoundaryDrive::IncomingCharacteristic(0.002), NsfBoundaryDrive::PressureOffset(-0.1)]
            {
                let mass = owner.zone().total_mass_kg();
                let energy = owner.zone().total_energy_j();
                let frame = owner.step_with_boundary(drive, &gate).unwrap();
                let expected = advance_original(&mut fluid, high, substeps,
                    frame.incoming_characteristic_m_s, drive);
                largest_escape = largest_escape.max(check_patch(&frame, expected, high, main.propagator().dt()).abs());
                let acoustic = main.preview_step(frame.inlet_flow_m3_s).unwrap();
                let inlet_work = acoustic.inlet_work_j();
                acoustic.accept();
                assert_eq!(owner.zone().state_hash(), fluid.state_hash());
                assert_eq!(owner.load().state(), main.state());
                assert_eq!(owner.zone().characteristic_patch_drive(), Some(drive));
                assert_eq!(owner.zone().boundary_drives().unwrap()[usize::from(!high)],
                    NsfBoundaryDrive::PressureOffset(BACKGROUND_PA));
                assert_eq!(frame.acoustic_work_residual_j,
                    frame.faces[usize::from(high)].acoustic_pressure_work_out_j - inlet_work);
                let mass_out: f64 = frame.faces.iter().map(|f| f.mass_out_kg).sum();
                let energy_out: f64 = frame.faces.iter().map(|f| f.total_energy_out_j).sum();
                assert!((frame.flow_mass_change_kg + mass_out).abs() < 3e-12 * mass);
                assert!((frame.flow_energy_change_j + energy_out).abs() < 3e-12 * energy);
            }
            assert!(largest_escape > 1e-13, "source remainder must carry its own flow");
        }
    }
}

#[test]
fn finite_upstream_uses_source_aperture_flow_and_work_without_absorbing_its_remainder() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        for substeps in [1, 3] {
            let mut owner = owner(high, substeps, true);
            let mut fluid = source_zone(high, substeps);
            let mut main = load(&fluid, reference(), FACE_CELLS as f64, 1.0, DgBoundary::Matched);
            let mut source = load(&fluid, reference(), PATCH_CELLS as f64, 1.0, DgBoundary::Closed);
            let mut largest_escape: f64 = 0.0;
            for rise in [2.0, 1.0, -1.0] {
                let frame = owner.step_with_upstream_pressure_rise(rise, &gate).unwrap();
                let upstream = frame.upstream.as_ref().unwrap();
                let expected = advance_original(&mut fluid, high, substeps, frame.incoming_characteristic_m_s,
                    NsfBoundaryDrive::IncomingCharacteristic(upstream.incoming_characteristic_m_s));
                largest_escape = largest_escape.max(check_patch(&frame, expected, high, main.propagator().dt()).abs());
                let patch = frame.source_patch.unwrap();
                assert_eq!(upstream.inlet_flow_m3_s, patch.normal_volume_flow_m3_s);
                assert!(upstream.pressure_residual_pa.abs() < 1e-6);
                assert!(frame.pressure_residual_pa.abs() < 1e-6);
                main.preview_step(frame.inlet_flow_m3_s).unwrap().accept();
                let acoustic = source.preview_step(patch.normal_volume_flow_m3_s).unwrap();
                let inlet_work = acoustic.inlet_work_j();
                acoustic.accept();
                let work = -rise * (patch.normal_volume_flow_m3_s * source.propagator().dt());
                assert_eq!(upstream.source_work_j, work);
                assert_eq!(upstream.acoustic_work_residual_j, patch.acoustic_pressure_work_out_j + work - inlet_work);
                assert_eq!(upstream.mass_projection_difference_kg, patch.mass_out_kg
                    - reference().density * (patch.normal_volume_flow_m3_s * source.propagator().dt()));
                assert_eq!(owner.zone().state_hash(), fluid.state_hash());
                assert_eq!(owner.load().state(), main.state());
                assert_eq!(owner.upstream_load().unwrap().state(), source.state());
            }
            assert!(largest_escape > 1e-13);
        }
    }
}

#[test]
fn source_aperture_refusals_leave_both_ducts_and_original_boundary_history_unchanged() {
    let gate = CancelGate::new_clock_free();
    let mut model = owner(true, 3, true);
    model.step_with_upstream_pressure_rise(2.0, &gate).unwrap();
    let snapshot = |model: &NsfDg<D3q27>| (
        model.zone().state_hash(), model.zone().steps(), model.zone().cumulative_ledger(),
        model.zone().boundary_drives(), model.zone().characteristic_patch_drive(),
        model.zone().last_boundary_exchange().copied(), model.zone().last_patch_exchange().copied(),
        model.load().state().clone(), model.upstream_load().unwrap().state().clone(), model.accepted_steps(),
    );
    let before = snapshot(&model);
    model.limits.projection_absolute_tolerance_j = 0.0;
    model.limits.projection_relative_tolerance = 0.0;
    assert!(matches!(model.step_with_upstream_pressure_rise(1.5, &gate), Err(NsfDgError::Projection { .. })));
    assert_eq!(snapshot(&model), before);
    model.limits.projection_relative_tolerance = 1.0;
    assert!(model.step_with_upstream_pressure_rise(f64::NAN, &gate).is_err());
    assert_eq!(snapshot(&model), before);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(model.step_with_upstream_pressure_rise(1.5, &cancelled), Err(NsfDgError::Cancelled)));
    assert_eq!(snapshot(&model), before);
    let mut uninterrupted = owner(true, 3, true);
    uninterrupted.step_with_upstream_pressure_rise(2.0, &gate).unwrap();
    let resumed = model.step_with_upstream_pressure_rise(1.5, &gate).unwrap();
    let expected = uninterrupted.step_with_upstream_pressure_rise(1.5, &gate).unwrap();
    assert_eq!(resumed.source_patch, expected.source_patch);
    assert_eq!(resumed.inlet_pressure_pa, expected.inlet_pressure_pa);
    assert_eq!(snapshot(&model), snapshot(&uninterrupted));
}

#[test]
fn source_aperture_refuses_wrong_face_and_upstream_whole_face_or_foreign_gas() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        let zone = patched_zone(high, 1);
        let main = load(&zone, reference(), FACE_CELLS as f64, 1.0, DgBoundary::Matched);
        assert!(matches!(NsfDg::new_with_source_patch(zone, main, high, budget(1), &gate), Err(NsfDgError::Input(_))));
        for (gas, cells) in [(reference(), FACE_CELLS as f64),
            (GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 120000.0).unwrap(), PATCH_CELLS as f64)]
        {
            let model = owner(high, 1, false);
            let upstream = load(model.zone(), gas, cells, 1.0, DgBoundary::Closed);
            assert!(matches!(model.with_upstream_load(upstream, &gate), Err(NsfDgError::Input(_))));
        }
    }
}
