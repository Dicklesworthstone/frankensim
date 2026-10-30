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

const SOURCE_PATCH: NsfCharacteristicPatch = NsfCharacteristicPatch {
    start: [0, 1], end: [1, 4],
};
const SOURCE_CELLS: usize = 3;
const SOURCE_REMAINDER_PA: f64 = 0.375;

fn dual_zone(high: bool, substeps: usize) -> NsfGrid3<D3q27> {
    let mut zone = bare_zone(high, substeps);
    let mut backgrounds = [NsfBoundaryDrive::PressureOffset(BACKGROUND_PA); 2];
    backgrounds[usize::from(!high)] = NsfBoundaryDrive::PressureOffset(SOURCE_REMAINDER_PA);
    zone.set_boundary_drives(backgrounds).unwrap();
    zone.with_characteristic_patch(high, PATCH, backgrounds[usize::from(high)]).unwrap()
        .with_characteristic_patch(!high, SOURCE_PATCH, backgrounds[usize::from(!high)]).unwrap()
}

fn dual_owner(high: bool, substeps: usize, upstream: bool) -> NsfDg<D3q27> {
    let gate = CancelGate::new_clock_free();
    let zone = dual_zone(high, substeps);
    let main = load(&zone, reference(), PATCH_CELLS as f64, 1.0, DgBoundary::Matched);
    let owner = NsfDg::new_on_patches(zone, main, high, budget(substeps), &gate).unwrap();
    if upstream {
        let source = load(owner.zone(), reference(), SOURCE_CELLS as f64, 1.0, DgBoundary::Closed);
        owner.with_upstream_load(source, &gate).unwrap()
    } else { owner }
}

#[test]
fn independent_source_and_bore_apertures_replay_original_owners_and_separate_both_remainders() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        for substeps in [1, 3] {
            for upstream in [false, true] {
                let mut owner = dual_owner(high, substeps, upstream);
                let mut fluid = dual_zone(high, substeps);
                let mut main = load(&fluid, reference(), PATCH_CELLS as f64, 1.0, DgBoundary::Matched);
                let mut source = load(&fluid, reference(), SOURCE_CELLS as f64, 1.0, DgBoundary::Closed);
                let side = usize::from(high);
                let other = 1 - side;
                let mut shapes = [PATCH; 2]; shapes[other] = SOURCE_PATCH;
                let mut counts = [PATCH_CELLS; 2]; counts[other] = SOURCE_CELLS;
                let backgrounds = fluid.boundary_drives().unwrap();
                assert_eq!(owner.source_patch(), Some(SOURCE_PATCH));
                assert_eq!(owner.coupled_patch(), Some(PATCH));
                assert_eq!(owner.patch_background_drive(), None);
                assert_eq!(owner.area_m2(), PATCH_CELLS as f64 * fluid.law().dx_m.powi(2));
                let mut escaped = [0.0_f64; 2];
                for pressure in [2.0, 1.0, -0.1] {
                    let mass = fluid.total_mass_kg();
                    let energy = fluid.total_energy_j();
                    let frame = if upstream { owner.step_with_upstream_pressure_rise(pressure, &gate) }
                        else { owner.step_with_boundary(NsfBoundaryDrive::PressureOffset(pressure), &gate) }.unwrap();
                    let mut drives = [Some(NsfBoundaryDrive::IncomingCharacteristic(frame.incoming_characteristic_m_s)); 2];
                    drives[other] = Some(if let Some(upstream) = &frame.upstream {
                        NsfBoundaryDrive::IncomingCharacteristic(upstream.incoming_characteristic_m_s)
                    } else { NsfBoundaryDrive::PressureOffset(pressure) });
                    fluid.set_boundary_drives_with_patches(backgrounds, drives).unwrap();
                    // Independent trace quadrature; neither patch accumulator
                    // nor root output supplies the expected pressure/flow/work.
                    let mut expected = [[0.0; 3]; 2];
                    for _ in 0..substeps {
                        let before = [fluid.open_boundary_traces(false).unwrap(), fluid.open_boundary_traces(true).unwrap()];
                        fluid.step().unwrap();
                        for face in 0..2 {
                            let after = fluid.open_boundary_traces(face != 0).unwrap();
                            for (old, new) in before[face].iter().zip(after) {
                                let patch = shapes[face];
                                if (patch.start[0]..patch.end[0]).contains(&old.cell[1])
                                    && (patch.start[1]..patch.end[1]).contains(&old.cell[2])
                                {
                                    expected[face][0] += 0.5 * (old.pressure_pa + new.pressure_pa);
                                    expected[face][1] += 0.5 * (old.outward_normal_velocity_m_s + new.outward_normal_velocity_m_s) * old.area_m2;
                                    expected[face][2] += 0.5 * ((old.pressure_pa - reference().pressure) * old.outward_normal_velocity_m_s
                                        + (new.pressure_pa - reference().pressure) * new.outward_normal_velocity_m_s) * old.area_m2 * fluid.law().dt_s;
                                } else {
                                    let remainder = if face == side { BACKGROUND_PA } else { SOURCE_REMAINDER_PA };
                                    assert_eq!(new.pressure_pa, reference().pressure + remainder);
                                }
                            }
                        }
                    }
                    let mut patches = [frame.coupled_patch.unwrap(); 2]; patches[other] = frame.source_patch.unwrap();
                    for face in 0..2 {
                        let patch = patches[face];
                        assert!((patch.mean_pressure_pa - expected[face][0] / (counts[face] * substeps) as f64).abs() < 1e-8);
                        assert!((patch.normal_volume_flow_m3_s - expected[face][1] / substeps as f64).abs()
                            < 1e-23 + 1e-11 * patch.normal_volume_flow_m3_s.abs());
                        assert!((patch.acoustic_pressure_work_out_j - expected[face][2]).abs()
                            < 1e-25 + 1e-11 * expected[face][2].abs());
                        let remainder = if face == side { BACKGROUND_PA } else { SOURCE_REMAINDER_PA };
                        let escape = frame.faces[face].normal_volume_flow_m3_s - patch.normal_volume_flow_m3_s;
                        escaped[face] = escaped[face].max(escape.abs());
                        let work = frame.faces[face].acoustic_pressure_work_out_j - patch.acoustic_pressure_work_out_j;
                        let expected_work = remainder * escape * main.propagator().dt();
                        let scale = frame.faces[face].acoustic_pressure_work_out_j.abs() + patch.acoustic_pressure_work_out_j.abs() + expected_work.abs();
                        assert!((work - expected_work).abs() < 512.0 * f64::EPSILON * scale.max(1e-30));
                        assert_eq!(owner.patch_background_drive_on(face != 0), Some(backgrounds[face]));
                    }
                    assert_eq!(frame.inlet_flow_m3_s, patches[side].normal_volume_flow_m3_s);
                    main.preview_step(frame.inlet_flow_m3_s).unwrap().accept();
                    if let Some(upstream) = &frame.upstream {
                        assert_eq!(upstream.inlet_flow_m3_s, patches[other].normal_volume_flow_m3_s);
                        let candidate = source.preview_step(upstream.inlet_flow_m3_s).unwrap();
                        let work = candidate.inlet_work_j(); candidate.accept();
                        assert_eq!(upstream.source_work_j, -pressure * (upstream.inlet_flow_m3_s * source.propagator().dt()));
                        assert_eq!(upstream.acoustic_work_residual_j,
                            patches[other].acoustic_pressure_work_out_j + upstream.source_work_j - work);
                        assert_eq!(owner.upstream_load().unwrap().state(), source.state());
                    }
                    assert_eq!(owner.zone().state_hash(), fluid.state_hash());
                    assert_eq!(owner.load().state(), main.state());
                    assert!((frame.flow_mass_change_kg + frame.faces.iter().map(|f| f.mass_out_kg).sum::<f64>()).abs() < 3e-12 * mass);
                    assert!((frame.flow_energy_change_j + frame.faces.iter().map(|f| f.total_energy_out_j).sum::<f64>()).abs() < 3e-12 * energy);
                }
                assert!(escaped.into_iter().all(|flow| flow > 1e-13));
            }
        }
    }
}

#[test]
fn dual_aperture_refusal_preserves_both_controls_fluxes_and_acoustic_histories() {
    let gate = CancelGate::new_clock_free();
    let mut model = dual_owner(true, 3, true);
    model.step_with_upstream_pressure_rise(2.0, &gate).unwrap();
    let snapshot = |m: &NsfDg<D3q27>| (
        m.zone().state_hash(), m.zone().cumulative_ledger(), m.zone().boundary_drives(),
        [m.zone().characteristic_patch_drive_on(false), m.zone().characteristic_patch_drive_on(true)],
        [m.zone().last_patch_exchange_on(false).copied(), m.zone().last_patch_exchange_on(true).copied()],
        m.load().state().clone(), m.upstream_load().unwrap().state().clone(), m.accepted_steps(),
    );
    let before = snapshot(&model);
    model.limits.projection_absolute_tolerance_j = 0.0;
    model.limits.projection_relative_tolerance = 0.0;
    assert!(matches!(model.step_with_upstream_pressure_rise(1.5, &gate), Err(NsfDgError::Projection { .. })));
    assert_eq!(snapshot(&model), before);
    model.limits.projection_relative_tolerance = 1.0;
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(model.step_with_upstream_pressure_rise(1.5, &cancelled), Err(NsfDgError::Cancelled)));
    assert_eq!(snapshot(&model), before);
    let mut uninterrupted = dual_owner(true, 3, true);
    uninterrupted.step_with_upstream_pressure_rise(2.0, &gate).unwrap();
    model.step_with_upstream_pressure_rise(1.5, &gate).unwrap();
    uninterrupted.step_with_upstream_pressure_rise(1.5, &gate).unwrap();
    assert_eq!(snapshot(&model), snapshot(&uninterrupted));
}

#[test]
fn dual_apertures_require_explicit_roles_and_each_loads_actual_area() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        for role in 0..3 {
            let zone = dual_zone(high, 1);
            let main = load(&zone, reference(), PATCH_CELLS as f64, 1.0, DgBoundary::Matched);
            let result = match role {
                0 => NsfDg::new(zone, main, high, budget(1), &gate),
                1 => NsfDg::new_on_patch(zone, main, high, budget(1), &gate),
                _ => NsfDg::new_with_source_patch(zone, main, high, budget(1), &gate),
            };
            assert!(matches!(result, Err(NsfDgError::Input(_))));
        }
        let zone = dual_zone(high, 1);
        let main = load(&zone, reference(), FACE_CELLS as f64, 1.0, DgBoundary::Matched);
        assert!(matches!(NsfDg::new_on_patches(zone, main, high, budget(1), &gate), Err(NsfDgError::Input(_))));
        let model = dual_owner(high, 1, false);
        let source = load(model.zone(), reference(), PATCH_CELLS as f64, 1.0, DgBoundary::Closed);
        assert!(matches!(model.with_upstream_load(source, &gate), Err(NsfDgError::Input(_))));
    }
}

fn reservoir_pair(high: bool, substeps: usize, upstream: bool) -> NsfDg<D3q27> {
    let gas = reference();
    let state = NsfOpenState { density_kg_m3: gas.density,
        velocity_m_s: [0.0; 3], temperature_k: gas.temperature };
    let zone = dual_zone(high, substeps).with_reservoir_apertures(vec![
        fs_lbm::nsf3::NsfReservoirAperture { axis: 1, high: false,
            cells: vec![[3, 0, 1], [4, 0, 1]], state },
        fs_lbm::nsf3::NsfReservoirAperture { axis: 1, high: true,
            cells: vec![[4, 3, 2]], state },
    ]).unwrap();
    let main = load(&zone, gas, PATCH_CELLS as f64, 1.0, DgBoundary::Matched);
    let gate = CancelGate::new_clock_free();
    let owner = NsfDg::new_on_patches(zone, main, high, budget(substeps), &gate).unwrap();
    if upstream {
        let source = load(owner.zone(), gas, SOURCE_CELLS as f64, 1.0, DgBoundary::Matched);
        owner.with_upstream_load(source, &gate).unwrap()
    } else { owner }
}

#[test]
fn indexed_breath_replays_original_controls_and_preserves_independent_escape_accounts() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        for upstream in [false, true] {
            for substeps in [1, 3] {
                let mut owner = reservoir_pair(high, substeps, upstream);
                let original = reservoir_pair(high, substeps, upstream);
                let mut fluid = original.zone;
                let mut main = original.load;
                let mut source = original.upstream;
                let mut control_values = [fluid.reservoir_aperture_at(0).unwrap().state,
                    fluid.reservoir_aperture_at(1).unwrap().state];
                for (index, pressure, temperature) in [(0, 100002.0, 300.01), (1, 99999.0, 300.0), (0, 100000.0, 300.0)] {
                    let frame = owner.step_with_reservoir_aperture(index, pressure, temperature, &gate).unwrap();
                    assert_eq!(frame.drive, NsfDgDrive::ReservoirAperture { aperture: index,
                        pressure_pa: pressure, temperature_k: temperature });
                    control_values[index] = fluid.reservoir_aperture_state_for(index, pressure, temperature).unwrap();
                    let mut incoming = [NsfBoundaryDrive::IncomingCharacteristic(frame.incoming_characteristic_m_s); 2];
                    incoming[usize::from(!high)] = NsfBoundaryDrive::IncomingCharacteristic(
                        frame.upstream.as_ref().map_or(0.0, |source| source.incoming_characteristic_m_s));
                    let backgrounds = fluid.boundary_drives().unwrap();
                    let trial = fluid.preview_boundary_interval_with_patches_and_reservoirs(
                        backgrounds, incoming.map(Some), &[(index, control_values[index])],
                        substeps, GRID_CELLS * substeps, || false).unwrap();
                    assert_eq!(trial.boundary_exchange(), &frame.faces);
                    assert_eq!(trial.reservoir_flows(), frame.reservoir_flows);
                    assert_eq!(trial.reservoir_exchanges(), frame.reservoir_exchanges);
                    assert_eq!(trial.reservoir_mass_out_kg(), frame.reservoir_mass_out_kg);
                    let balance = frame.flow_mass_change_kg + frame.reservoir_mass_out_kg
                        + frame.faces.iter().map(|face| face.mass_out_kg).sum::<f64>();
                    assert!(balance.abs() <= 3e-12 * trial.previous_mass_kg());
                    for opening in 0..2 {
                        assert_eq!(trial.grid().reservoir_aperture_at(opening).unwrap().state, control_values[opening]);
                        assert_eq!(frame.reservoir_flows.iter().filter(|tile| tile.aperture == opening).count(),
                            if opening == 0 { 2 } else { 1 });
                        assert_eq!(frame.reservoir_exchanges[opening].aperture, opening);
                    }
                    assert!(!frame.reservoir_flows.iter().any(|tile| tile.aperture == 2));
                    assert!(frame.reservoir_exchanges.get(2).is_none());
                    trial.accept();
                    main.preview_step(frame.inlet_flow_m3_s).unwrap().accept();
                    if let Some(source) = &mut source {
                        let upstream = frame.upstream.as_ref().unwrap();
                        assert_eq!(upstream.source_work_j, 0.0);
                        source.load.preview_step(upstream.inlet_flow_m3_s).unwrap().accept();
                        assert_eq!(owner.upstream_load().unwrap().state(), source.load.state());
                    }
                    assert_eq!(owner.zone().state_hash(), fluid.state_hash());
                    assert_eq!(owner.load().state(), main.state());
                }
            }
        }
    }
}

#[test]
fn indexed_breath_refusals_preserve_every_opening_and_both_loads() {
    let gate = CancelGate::new_clock_free();
    let mut owner = reservoir_pair(true, 1, true);
    owner.step_with_reservoir_aperture(0, 100002.0, 300.01, &gate).unwrap();
    let snapshot = |owner: &NsfDg<D3q27>| (
        owner.zone().state_hash(), owner.zone().cumulative_ledger(),
        owner.zone().last_reservoir_exchanges().to_vec(), owner.zone().last_reservoir_flows().to_vec(),
        [owner.zone().reservoir_aperture_at(0).unwrap().state, owner.zone().reservoir_aperture_at(1).unwrap().state],
        owner.load().state().clone(), owner.upstream_load().unwrap().state().clone(), owner.accepted_steps(),
    );
    let before = snapshot(&owner);
    assert!(owner.step_with_reservoir(100002.0, 300.01, &gate).is_err());
    for (index, pressure, temperature) in [(2, 100002.0, 300.01), (0, f64::NAN, 300.0), (1, 100000.0, 1000.0)] {
        assert!(owner.step_with_reservoir_aperture(index, pressure, temperature, &gate).is_err());
        assert_eq!(snapshot(&owner), before);
    }
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(owner.step_with_reservoir_aperture(1, 100000.0, 300.0, &cancelled), Err(NsfDgError::Cancelled)));
    assert_eq!(snapshot(&owner), before);
    let mut replay = reservoir_pair(true, 1, true);
    replay.step_with_reservoir_aperture(0, 100002.0, 300.01, &gate).unwrap();
    owner.step_with_reservoir_aperture(1, 99999.0, 300.0, &gate).unwrap();
    replay.step_with_reservoir_aperture(1, 99999.0, 300.0, &gate).unwrap();
    assert_eq!(snapshot(&owner), snapshot(&replay));
}
