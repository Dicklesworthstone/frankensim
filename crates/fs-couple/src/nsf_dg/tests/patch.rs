//! A resolved aperture and its escape opening retain distinct physical owners.
use super::*;
mod source;

const DIMS: [usize; 3] = [8, 4, 4];
const GRID_CELLS: usize = 128;
const FACE_CELLS: usize = 16;
const PATCH_CELLS: usize = 8;
const BACKGROUND_PA: f64 = -0.25;
const PATCH: NsfCharacteristicPatch = NsfCharacteristicPatch {
    start: [1, 0],
    end: [3, 4],
};

fn reference() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap()
}

fn bare_zone(high: bool, substeps: usize) -> NsfGrid3<D3q27> {
    let reference = reference();
    let velocity = (3.0 * reference.specific_gas_constant * reference.temperature).sqrt();
    let c = (reference.gamma / 3.0).sqrt();
    let sign = if high { 1.0 } else { -1.0 };
    let face = NsfCharacteristicBoundary {
        reference_pressure_pa: reference.pressure,
        reference_temperature_k: reference.temperature,
        reference_velocity_m_s: [0.0; 3],
        incoming_characteristic_m_s: 0.0,
    };
    let mut zone = NsfGrid3::<D3q27>::new(
        NsfSpec {
            dims: DIMS,
            periodic: [false; 3],
            gas: GasSpec::dry_air_ussa1976(),
            temperature_k: reference.temperature,
            pressure_pa: reference.pressure,
            dx_m: velocity / f64::from(RATE) / substeps as f64,
            temperature_window_k: [250.0, 600.0],
            maximum_mach: 0.2,
            maximum_temperature_ratio: 1.0,
            bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "authored aperture composition fixture; not jet qualification",
            maximum_cells: GRID_CELLS,
            wall: WallThermal::Adiabatic,
            workers: 1,
            transport_points: 17,
        },
        |_| false,
        |x| {
            let distance = if high { DIMS[0] - 1 - x[0] } else { x[0] };
            let z = (distance as f64 - 3.0) / 4.0;
            let transverse = 0.5 + 0.1 * x[1] as f64 + 0.2 * x[2] as f64;
            let pulse = 1e-5 * transverse * fs_math::det::exp(-z * z);
            (
                1.0 + pulse / reference.gamma,
                [sign * c * pulse / reference.gamma, 0.0, 0.0],
                reference.temperature * (1.0 + (reference.gamma - 1.0) / reference.gamma * pulse),
            )
        },
    )
    .unwrap()
    .with_characteristic_axis(NsfCharacteristicAxis {
        axis: 0,
        faces: [face; 2],
        maximum_relative_perturbation: 0.05,
    })
    .unwrap();
    zone.set_boundary_drives([NsfBoundaryDrive::PressureOffset(BACKGROUND_PA); 2])
        .unwrap();
    zone
}

fn patched_zone(high: bool, substeps: usize) -> NsfGrid3<D3q27> {
    bare_zone(high, substeps)
        .with_characteristic_patch(high, PATCH, NsfBoundaryDrive::PressureOffset(BACKGROUND_PA))
        .unwrap()
}

fn budget(substeps: usize) -> NsfDgLimits {
    let mut budget = limits(substeps);
    budget.maximum_cells = GRID_CELLS;
    budget.maximum_evaluations = 128;
    budget.maximum_cell_steps = GRID_CELLS * substeps * budget.maximum_evaluations;
    budget
}

fn coupled(high: bool, substeps: usize, upstream: bool) -> NsfDg<D3q27> {
    let gate = CancelGate::new_clock_free();
    let zone = patched_zone(high, substeps);
    let main = load(
        &zone,
        reference(),
        PATCH_CELLS as f64,
        1.0,
        DgBoundary::Matched,
    );
    let owner = NsfDg::new_on_patch(zone, main, high, budget(substeps), &gate).unwrap();
    if upstream {
        let source = load(
            owner.zone(),
            reference(),
            FACE_CELLS as f64,
            1.0,
            DgBoundary::Closed,
        );
        owner.with_upstream_load(source, &gate).unwrap()
    } else {
        owner
    }
}

fn check_complement(owner: &NsfDg<D3q27>, high: bool, frame: &NsfDgFrame) -> f64 {
    let background = NsfBoundaryDrive::PressureOffset(BACKGROUND_PA);
    assert_eq!(owner.patch_background_drive(), Some(background));
    assert_eq!(
        owner.zone().boundary_drives().unwrap()[usize::from(high)],
        background
    );
    for trace in owner.zone().open_boundary_traces(high).unwrap() {
        if !(PATCH.start[0]..PATCH.end[0]).contains(&trace.cell[1]) {
            assert_eq!(trace.pressure_pa, reference().pressure + BACKGROUND_PA);
        }
    }
    let face = frame.faces[usize::from(high)];
    let patch = frame.coupled_patch.unwrap();
    let escaping_flow = face.normal_volume_flow_m3_s - patch.normal_volume_flow_m3_s;
    let escaping_work = face.acoustic_pressure_work_out_j - patch.acoustic_pressure_work_out_j;
    let expected = BACKGROUND_PA * escaping_flow * owner.load().propagator().dt();
    let scale = face.acoustic_pressure_work_out_j.abs()
        + patch.acoustic_pressure_work_out_j.abs()
        + expected.abs();
    assert!((escaping_work - expected).abs() <= 256.0 * f64::EPSILON * scale.max(1e-30));
    escaping_flow
}

#[test]
fn oriented_aperture_uses_its_area_flow_and_work_without_absorbing_the_escape_opening() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        for substeps in [1, 3] {
            let mut owner = coupled(high, substeps, false);
            let mut fluid = patched_zone(high, substeps);
            let mut main = load(
                &fluid,
                reference(),
                PATCH_CELLS as f64,
                1.0,
                DgBoundary::Matched,
            );
            assert_eq!(owner.coupled_patch(), Some(PATCH));
            assert_eq!(
                owner.area_m2(),
                PATCH_CELLS as f64 * fluid.law().dx_m.powi(2)
            );
            let mut largest_escape: f64 = 0.0;
            for pressure in [2.0, 1.0, -0.25] {
                let boundary = NsfBoundaryDrive::PressureOffset(pressure);
                let mass = owner.zone().total_mass_kg();
                let energy = owner.zone().total_energy_j();
                let frame = owner.step_with_boundary(boundary, &gate).unwrap();
                let patch = frame.coupled_patch.unwrap();
                assert_eq!(frame.inlet_flow_m3_s, patch.normal_volume_flow_m3_s);
                assert!(frame.pressure_residual_pa.abs() < 1e-6);

                // Replay the accepted controls through the original owners,
                // independently of the coupled solver's candidate lifetime.
                let mut drives = [boundary; 2];
                drives[usize::from(high)] = NsfBoundaryDrive::PressureOffset(BACKGROUND_PA);
                let raw = fluid
                    .preview_boundary_interval_with_patch(
                        drives,
                        NsfBoundaryDrive::IncomingCharacteristic(frame.incoming_characteristic_m_s),
                        substeps,
                        GRID_CELLS * substeps,
                        || false,
                    )
                    .unwrap();
                assert_eq!(raw.patch_exchange(), Some(&patch));
                assert_eq!(raw.boundary_exchange(), &frame.faces);
                raw.accept();
                let acoustic = main.preview_step(patch.normal_volume_flow_m3_s).unwrap();
                let inlet_work = acoustic.inlet_work_j();
                acoustic.accept();
                assert_eq!(fluid.state_hash(), owner.zone().state_hash());
                assert_eq!(main.state(), owner.load().state());
                assert_eq!(
                    frame.acoustic_work_residual_j,
                    patch.acoustic_pressure_work_out_j - inlet_work
                );
                assert_eq!(
                    frame.mass_projection_difference_kg,
                    patch.mass_out_kg
                        - reference().density * (frame.inlet_flow_m3_s * main.propagator().dt())
                );

                // The complete face already contains the aperture; total NSF
                // conservation must count it once, including escaped energy.
                let mass_out: f64 = frame.faces.iter().map(|face| face.mass_out_kg).sum();
                let energy_out: f64 = frame.faces.iter().map(|face| face.total_energy_out_j).sum();
                assert!((frame.flow_mass_change_kg + mass_out).abs() < 3e-12 * mass);
                assert!((frame.flow_energy_change_j + energy_out).abs() < 3e-12 * energy);
                assert!((frame.flow.open_boundary_energy_j + energy_out).abs() < 3e-12 * energy);
                largest_escape = largest_escape.max(check_complement(&owner, high, &frame).abs());
            }
            assert!(
                largest_escape > 1e-13,
                "the fixed complement must carry its own flow"
            );
        }
    }
}

#[test]
fn finite_upstream_retains_the_complete_face_while_main_bore_receives_only_the_patch() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        let mut owner = coupled(high, 1, true);
        let mut fluid = patched_zone(high, 1);
        let mut main = load(
            &fluid,
            reference(),
            PATCH_CELLS as f64,
            1.0,
            DgBoundary::Matched,
        );
        let mut source = load(
            &fluid,
            reference(),
            FACE_CELLS as f64,
            1.0,
            DgBoundary::Closed,
        );
        let expected_area = FACE_CELLS as f64 * fluid.law().dx_m.powi(2);
        assert!(
            (owner.upstream_load().unwrap().duct().inlet_area_m2() - expected_area).abs()
                < 32.0 * f64::EPSILON * expected_area
        );
        let mut largest_escape: f64 = 0.0;
        for rise in [2.0, 2.0, -1.0, 0.0] {
            let frame = owner.step_with_upstream_pressure_rise(rise, &gate).unwrap();
            let upstream = frame.upstream.as_ref().unwrap();
            let patch = frame.coupled_patch.unwrap();
            assert_eq!(upstream.face, usize::from(!high));
            assert_eq!(frame.inlet_flow_m3_s, patch.normal_volume_flow_m3_s);
            assert_eq!(
                upstream.inlet_flow_m3_s,
                frame.faces[upstream.face].normal_volume_flow_m3_s
            );
            assert!(frame.pressure_residual_pa.abs() < 1e-6);
            assert!(upstream.pressure_residual_pa.abs() < 1e-6);

            let mut drives =
                [NsfBoundaryDrive::IncomingCharacteristic(upstream.incoming_characteristic_m_s); 2];
            drives[usize::from(high)] = NsfBoundaryDrive::PressureOffset(BACKGROUND_PA);
            let raw = fluid
                .preview_boundary_interval_with_patch(
                    drives,
                    NsfBoundaryDrive::IncomingCharacteristic(frame.incoming_characteristic_m_s),
                    1,
                    GRID_CELLS,
                    || false,
                )
                .unwrap();
            assert_eq!(raw.patch_exchange(), Some(&patch));
            assert_eq!(raw.boundary_exchange(), &frame.faces);
            raw.accept();
            main.preview_step(frame.inlet_flow_m3_s).unwrap().accept();
            source
                .preview_step(upstream.inlet_flow_m3_s)
                .unwrap()
                .accept();
            assert_eq!(fluid.state_hash(), owner.zone().state_hash());
            assert_eq!(main.state(), owner.load().state());
            assert_eq!(source.state(), owner.upstream_load().unwrap().state());
            let expected_work = -rise * (upstream.inlet_flow_m3_s * main.propagator().dt());
            assert_eq!(upstream.source_work_j, expected_work);
            largest_escape = largest_escape.max(check_complement(&owner, high, &frame).abs());
        }
        assert!(largest_escape > 1e-13);
    }
}

#[test]
fn rejected_live_aperture_projection_preserves_both_duct_histories_and_all_boundary_accounts() {
    let gate = CancelGate::new_clock_free();
    let mut owner = coupled(true, 3, true);
    owner.step_with_upstream_pressure_rise(2.0, &gate).unwrap();
    let snapshot = |owner: &NsfDg<D3q27>| {
        (
            owner.zone().state_hash(),
            owner.zone().steps(),
            owner.zone().cumulative_ledger(),
            (
                owner.zone().boundary_drives(),
                owner.zone().characteristic_patch_drive(),
                owner.patch_background_drive(),
                owner.zone().characteristic_axis().copied(),
            ),
            (
                owner.zone().last_boundary_exchange().copied(),
                owner.zone().last_patch_exchange().copied(),
            ),
            (
                owner.zone().open_boundary_traces(false),
                owner.zone().open_boundary_traces(true),
            ),
            owner.load().state().clone(),
            owner.upstream_load().unwrap().state().clone(),
            owner.accepted_steps(),
        )
    };
    let before = snapshot(&owner);
    owner.limits.projection_absolute_tolerance_j = 0.0;
    owner.limits.projection_relative_tolerance = 0.0;
    assert!(matches!(
        owner.step_with_upstream_pressure_rise(1.5, &gate),
        Err(NsfDgError::Projection { .. })
    ));
    assert_eq!(before, snapshot(&owner));
    owner.limits.projection_relative_tolerance = 1.0;

    // Continuing after refusal must match a run with no rejected interval,
    // including both original acoustic histories and accepted patch controls.
    let mut uninterrupted = coupled(true, 3, true);
    uninterrupted
        .step_with_upstream_pressure_rise(2.0, &gate)
        .unwrap();
    let resumed = owner.step_with_upstream_pressure_rise(1.5, &gate).unwrap();
    let expected = uninterrupted
        .step_with_upstream_pressure_rise(1.5, &gate)
        .unwrap();
    assert_eq!(resumed.coupled_patch, expected.coupled_patch);
    assert_eq!(resumed.inlet_pressure_pa, expected.inlet_pressure_pa);
    assert_eq!(
        resumed.upstream.unwrap().inlet_pressure_pa,
        expected.upstream.unwrap().inlet_pressure_pa
    );
    assert_eq!(snapshot(&owner), snapshot(&uninterrupted));
}

#[test]
fn aperture_admission_rejects_hidden_wrong_side_and_complete_face_area_mismatches() {
    let gate = CancelGate::new_clock_free();
    for high in [false, true] {
        let zone = patched_zone(high, 1);
        let main = load(
            &zone,
            reference(),
            FACE_CELLS as f64,
            1.0,
            DgBoundary::Matched,
        );
        assert!(matches!(
            NsfDg::new_on_patch(zone, main, high, budget(1), &gate),
            Err(NsfDgError::Input(_))
        ));

        let zone = patched_zone(high, 1);
        let main = load(
            &zone,
            reference(),
            PATCH_CELLS as f64,
            1.0,
            DgBoundary::Matched,
        );
        assert!(matches!(
            NsfDg::new(zone, main, high, budget(1), &gate),
            Err(NsfDgError::Input(_))
        ));

        let zone = patched_zone(high, 1);
        let main = load(
            &zone,
            reference(),
            PATCH_CELLS as f64,
            1.0,
            DgBoundary::Matched,
        );
        assert!(matches!(
            NsfDg::new_on_patch(zone, main, !high, budget(1), &gate),
            Err(NsfDgError::Input(_))
        ));

        let zone = bare_zone(high, 1);
        let main = load(
            &zone,
            reference(),
            FACE_CELLS as f64,
            1.0,
            DgBoundary::Matched,
        );
        assert!(matches!(
            NsfDg::new_on_patch(zone, main, high, budget(1), &gate),
            Err(NsfDgError::Input(_))
        ));

        let owner = coupled(high, 1, false);
        let source = load(
            owner.zone(),
            reference(),
            PATCH_CELLS as f64,
            1.0,
            DgBoundary::Closed,
        );
        assert!(matches!(
            owner.with_upstream_load(source, &gate),
            Err(NsfDgError::Input(_))
        ));
    }
}
