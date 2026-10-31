//! G0/G3 original-owner composition and rejection checks, not jet qualification.
use super::*;
use super::super::thermoviscous::NsfThermoviscousDg;
use super::super::thermoviscous::performance::NsfThermoviscousDgPerformance;
use crate::cavity_modal::interior::thermoviscous::{ThermoviscousInterior, source_medium};
use crate::cavity_modal::interior::thermoviscous::dg::{ThermoviscousDg, ThermoviscousDgLimits};
use crate::cavity_modal::interior::thermoviscous::dg::load::ThermoviscousDgLoad;
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_feec::acoustics::helmholtz::multiport::BoundaryMode;
use fs_feec::acoustics::thermoviscous::{ThermoviscousBoundary, ThermoviscousFe, ThermoviscousLimits};
use fs_feec::acoustics::thermoviscous::time::ThermoviscousTimeLimits;
use fs_feec::fixtures::kuhn_cube;
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};
use fs_rep_mesh::TetComplex;

mod radiation;
mod upstream;
mod radiation_observation;

const BAND: [f64; 2] = [20.0, 2000.0];

fn reference() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap()
}

// A manufactured cube has two opposed duct ports and one distinct source port.
// All three areas are side², so its source can match an original NSF tile set.
fn field_source<'c>(mesh: &'c TetComplex, unit: &[[f64; 3]], side: f64,
    gas: GasState) -> ThermoviscousInterior<'c>
{
    let xyz: Vec<_> = unit.iter().map(|p| p.map(|x| x * side)).collect();
    let geometry = HelmholtzFe::new(mesh, &xyz, 2, Medium {
        density_kg_m3: gas.density, sound_speed_m_s: gas.sound_speed,
    }).unwrap();
    let faces = geometry.boundary_faces();
    let opening = |axis: usize, coordinate: f64| faces.iter().copied().filter(|face| {
        face.iter().all(|&v| (xyz[v as usize][axis] - coordinate).abs() < 1e-14)
    }).collect::<Vec<_>>();
    let left = opening(0, 0.0);
    let right = opening(0, side);
    let input = opening(1, side);
    let walls: Vec<_> = faces.into_iter().filter(|face| {
        !left.contains(face) && !right.contains(face) && !input.contains(face)
    }).collect();
    let field = ThermoviscousFe::new(mesh, &xyz, source_medium(&gas, 0.0).unwrap(),
        &ThermoviscousBoundary {
            no_slip_faces: walls.clone(), isothermal_faces: walls,
            ports: vec![BoundaryMode::uniform(&left), BoundaryMode::uniform(&right),
                BoundaryMode::uniform(&input)],
        }, ThermoviscousLimits::default(), &mut || false,
    ).unwrap();
    // This original synthetic reservoir isolates the retained caloric channel;
    // it is not a measured gas-absorption law or an instrument parameter.
    let law = read_relaxation(&mut
        "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500"
            .split_ascii_whitespace()).unwrap();
    let medium = RelaxingMedium::new(gas,
        vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    ThermoviscousInterior::from_relaxation(field, law.evaluate(&medium).unwrap(), BAND).unwrap()
}

fn acoustic<'a, 'c>(source: &'a ThermoviscousInterior<'c>, area: f64,
    dt: f64, cork: DgBoundary) -> ThermoviscousDg<'a, 'c>
{
    let gate = CancelGate::new_clock_free();
    let duct = |length, end| {
        let duct = DgDuct::new(&Duct { segments: vec![Segment::Cylinder {
            radius: (area / std::f64::consts::PI).sqrt(), length,
        }] }, source.gas(), DgMesh { elements: 2, order: 2 },
            DgBoundary::Port(PortId(7)), end).unwrap()
            .with_molecular_relaxation(source.relaxation().unwrap().clone()).unwrap();
        DgPortLoad::new(duct, dt, 1e-11, MAX_DG_PROPAGATOR_DOFS).unwrap()
    };
    let prepared = source.prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false).unwrap();
    let mut initial = prepared.zero_state();
    initial.molecular_temperature_k[0][0] = 0.0002;
    ThermoviscousDg::new(source, prepared, duct(0.03, DgBoundary::Matched), 1,
        ThermoviscousDgLimits::default(), &gate).unwrap()
        .with_additional_duct(0, duct(0.012, cork), &gate).unwrap()
        .with_initial_state(initial).unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    let scale = actual.abs().max(expected.abs()).max(1e-20);
    assert!((actual - expected).abs() <= tolerance * scale,
        "{actual:e} != {expected:e}, scale {scale:e}");
}

macro_rules! snapshot {
    ($owner:expr) => {{
        let owner = $owner;
        let field = owner.load().system();
        (owner.zone().state_hash(), owner.zone().cumulative_ledger(),
            owner.zone().boundary_drives(), owner.accepted_steps(), field.state().clone(),
            field.load().state().clone(), field.load_at_port(0).unwrap().state().clone(),
            field.accepted_steps(), owner.load().port_response(&CancelGate::new_clock_free()).unwrap())
    }};
}

#[test]
fn external_flow_port_keeps_the_complete_schur_response_and_dropped_trials() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let side: f64 = 0.004;
    let area = side * side;
    let dt = 1.0 / f64::from(RATE);
    let source = field_source(&mesh, &unit, side, reference());
    let make = || acoustic(&source, area, dt, DgBoundary::Closed);
    assert!(ThermoviscousDgLoad::new(make(), 1, &gate).is_err());
    assert!(ThermoviscousDgLoad::new(make(), 3, &gate).is_err());
    let mut advanced = make();
    advanced.step(&[0.0, 0.0, 0.1], &[], &gate).unwrap();
    assert!(ThermoviscousDgLoad::new(advanced, 2, &gate).is_err());
    let mut load = ThermoviscousDgLoad::new(make(), 2, &gate).unwrap();
    let mut pressure_driven = make();
    assert_eq!(load.port(), 2);
    near(load.area_m2(), area, 64.0 * f64::EPSILON);
    assert_eq!(load.duration_s(), dt);
    let before = (load.system().state().clone(), load.system().load().state().clone(),
        load.system().load_at_port(0).unwrap().state().clone(), load.port_response(&gate).unwrap());
    {
        let trial = load.preview_step(2e-8, &gate).unwrap();
        assert_eq!(trial.frame().step, 1);
        assert!(trial.frame().stored_energy_j > 0.0);
    }
    assert_eq!(before, (load.system().state().clone(), load.system().load().state().clone(),
        load.system().load_at_port(0).unwrap().state().clone(), load.port_response(&gate).unwrap()));
    assert_eq!(load.system().accepted_steps(), 0);
    assert!(load.system().last_frame().is_none());
    for flow in [2e-8, -1e-8, 0.4e-8, 1e-8] {
        let response = load.port_response(&gate).unwrap();
        assert!(response.impedance_pa_s_m3 > 0.0);
        let pressure = response.history_pressure_pa + response.impedance_pa_s_m3 * flow;
        let trial = load.preview_step(flow, &gate).unwrap();
        near(trial.mean_pressure_pa(), pressure, 1e-12);
        near(trial.frame().field.mean_inward_flow_m3_s[2], flow, 1e-9);
        near(trial.inlet_work_j(), pressure * (flow * dt), 1e-9);
        // The preexisting pressure-driven composition independently eliminates
        // the internal ports at this same external pressure.
        let old = pressure_driven.step(&[0.0, 0.0, pressure], &[], &gate).unwrap();
        near(old.field.mean_inward_flow_m3_s[2], flow, 1e-8);
        let scale = pressure.abs() + old.mean_pressure_pa.abs()
            + old.additional_ducts[0].mean_pressure_pa.abs();
        assert!((old.mean_pressure_pa - trial.frame().mean_pressure_pa).abs()
            <= 1e-8 * scale.max(1e-20));
        assert!((old.additional_ducts[0].mean_pressure_pa
            - trial.frame().additional_ducts[0].mean_pressure_pa).abs()
            <= 1e-8 * scale.max(1e-20));
        trial.accept();
    }
    let before = (load.system().state().clone(), load.system().load().state().clone(),
        load.system().load_at_port(0).unwrap().state().clone());
    assert!(load.preview_step(f64::NAN, &gate).is_err());
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(load.preview_step(1e-8, &cancelled).is_err());
    assert_eq!(before, (load.system().state().clone(), load.system().load().state().clone(),
        load.system().load_at_port(0).unwrap().state().clone()));
    assert_eq!(load.system().accepted_steps(), 4);
}

#[test]
fn resolved_interval_replays_original_populations_complete_fields_and_both_ducts() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let gas = reference();
    let dt = 1.0 / f64::from(RATE);
    for high in [false, true] {
        for substeps in [1, 3] {
            let fluid = zone(GasSpec::dry_air_ussa1976(), gas, high, substeps);
            let side = fluid.law().dx_m;
            let source = field_source(&mesh, &unit, side, gas);
            let system = acoustic(&source, side * side, dt, DgBoundary::Closed);
            let mut owner = NsfThermoviscousDg::new(fluid,
                ThermoviscousDgLoad::new(system, 2, &gate).unwrap(), high,
                limits(substeps), &gate).unwrap();
            let mut original_fluid = zone(GasSpec::dry_air_ussa1976(), gas, high, substeps);
            let prepared = source.prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false).unwrap();
            let mut molecular_heat = 0.0;
            for pressure in [2.0, 0.0, -1.0] {
                let initial = owner.load().system().state().clone();
                let mut ducts = [owner.load().system().load().clone(),
                    owner.load().system().load_at_port(0).unwrap().clone()];
                let boundary = NsfBoundaryDrive::PressureOffset(pressure);
                let frame = owner.step_with_boundary(boundary, &gate).unwrap();
                let mut drives = [boundary; 2];
                drives[usize::from(high)] =
                    NsfBoundaryDrive::IncomingCharacteristic(frame.incoming_characteristic_m_s);
                let fluid = original_fluid.preview_boundary_interval(drives, substeps,
                    CELLS * substeps, || false).unwrap();
                assert_eq!(fluid.boundary_exchange(), &frame.faces);
                fluid.accept();
                assert_eq!(original_fluid.state_hash(), owner.zone().state_hash());
                let complete = owner.load().system().last_frame().unwrap();
                let cork = &complete.additional_ducts[0];
                let replay = prepared.step(&initial,
                    &[cork.mean_pressure_pa, complete.mean_pressure_pa, frame.inlet_pressure_pa],
                    &[], &mut || false).unwrap();
                assert_eq!(complete.field.state, replay.state);
                assert_eq!(complete.field.mean_inward_flow_m3_s, replay.mean_inward_flow_m3_s);
                near(replay.mean_inward_flow_m3_s[2], frame.inlet_flow_m3_s, 1e-9);
                for (i, port) in [1, 0].into_iter().enumerate() {
                    let trial = ducts[i].preview_step(-replay.mean_inward_flow_m3_s[port]).unwrap();
                    let output = trial.accept();
                    assert_eq!(if i == 0 { &complete.duct } else { &cork.duct }, &output);
                    assert_eq!(owner.load().system().load_at_port(port).unwrap().state(), ducts[i].state());
                }
                let exchange = frame.faces[usize::from(high)];
                assert_eq!(frame.acoustic_work_residual_j,
                    exchange.acoustic_pressure_work_out_j - complete.authored_external_work_j);
                assert!(frame.pressure_residual_pa.abs() < 1e-6);
                assert!(complete.balance_residual_j.abs() <= complete.balance_tolerance_j);
                molecular_heat += complete.field.work.molecular_relaxation_heat_j;
            }
            assert!(molecular_heat > 0.0);
            assert_eq!(owner.accepted_steps(), 3);
            assert_eq!(owner.load().system().accepted_steps(), 3);
            assert_eq!(owner.zone().steps(), 3 * substeps);
        }
    }
}

#[test]
fn finite_cork_return_changes_the_same_resolved_flow_and_main_bore() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let gas = reference();
    let fluid = zone(GasSpec::dry_air_ussa1976(), gas, true, 1);
    let side = fluid.law().dx_m;
    let source = field_source(&mesh, &unit, side, gas);
    let make = |cork| NsfThermoviscousDg::new(
        zone(GasSpec::dry_air_ussa1976(), gas, true, 1),
        ThermoviscousDgLoad::new(acoustic(&source, side * side,
            1.0 / f64::from(RATE), cork), 2, &gate).unwrap(), true, limits(1), &gate).unwrap();
    let mut closed = make(DgBoundary::Closed);
    let mut matched = make(DgBoundary::Matched);
    let mut difference: f64 = 0.0;
    for _ in 0..32 {
        let a = closed.step(0.0, &gate).unwrap();
        let b = matched.step(0.0, &gate).unwrap();
        difference = difference.max((a.inlet_pressure_pa - b.inlet_pressure_pa).abs());
    }
    assert!(difference > 1e-8, "finite cork reaction must reach the nonlinear source: {difference:e}");
    assert_ne!(closed.zone().state_hash(), matched.zone().state_hash());
    assert_ne!(closed.load().system().state(), matched.load().system().state());
    assert_ne!(closed.load().system().load().state(), matched.load().system().load().state());
}

#[test]
fn late_projection_refusal_and_cancellation_preserve_every_original_owner() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let gas = reference();
    let fluid = zone(GasSpec::dry_air_ussa1976(), gas, true, 3);
    let side = fluid.law().dx_m;
    let source = field_source(&mesh, &unit, side, gas);
    let make = |budget| NsfThermoviscousDg::new(
        zone(GasSpec::dry_air_ussa1976(), gas, true, 3),
        ThermoviscousDgLoad::new(acoustic(&source, side * side,
            1.0 / f64::from(RATE), DgBoundary::Closed), 2, &gate).unwrap(),
        true, budget, &gate).unwrap();
    let mut exact = limits(3); exact.projection_relative_tolerance = 0.0;
    let mut refused = make(exact);
    let before = snapshot!(&refused);
    assert!(matches!(refused.step(0.0, &gate), Err(NsfDgError::Projection { .. })));
    assert_eq!(before, snapshot!(&refused));
    assert!(refused.load().system().last_frame().is_none());
    let mut owner = make(limits(3));
    owner.step(0.0, &gate).unwrap();
    let before = snapshot!(&owner);
    assert!(owner.step_with_boundary(NsfBoundaryDrive::PressureOffset(f64::NAN), &gate).is_err());
    assert_eq!(before, snapshot!(&owner));
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(owner.step(0.0, &cancelled), Err(NsfDgError::Cancelled)));
    assert_eq!(before, snapshot!(&owner));
    owner.step(0.0, &gate).unwrap();
    assert_eq!(owner.accepted_steps(), 2);
    assert_eq!(owner.load().system().load().state().step_count(), 2);
    assert_eq!(owner.load().system().load_at_port(0).unwrap().state().step_count(), 2);
}

#[test]
fn resolved_binding_requires_the_original_gas_area_clock_and_cold_histories() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let gas = reference();
    let fluid = zone(GasSpec::dry_air_ussa1976(), gas, true, 1);
    let dx = fluid.law().dx_m;
    for (state, side, dt) in [
        (GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 120000.0).unwrap(), dx, 1.0 / f64::from(RATE)),
        (gas, dx * 1.01, 1.0 / f64::from(RATE)),
        (gas, dx, 2.0 / f64::from(RATE)),
    ] {
        let source = field_source(&mesh, &unit, side, state);
        let load = ThermoviscousDgLoad::new(acoustic(&source, side * side, dt,
            DgBoundary::Closed), 2, &gate).unwrap();
        assert!(NsfThermoviscousDg::new(zone(GasSpec::dry_air_ussa1976(), gas, true, 1),
            load, true, limits(1), &gate).is_err());
    }
    let source = field_source(&mesh, &unit, dx, gas);
    let mut load = ThermoviscousDgLoad::new(acoustic(&source, dx * dx,
        1.0 / f64::from(RATE), DgBoundary::Closed), 2, &gate).unwrap();
    load.preview_step(1e-8, &gate).unwrap().accept();
    assert!(NsfThermoviscousDg::new(zone(GasSpec::dry_air_ussa1976(), gas, true, 1),
        load, true, limits(1), &gate).is_err());
}

#[test]
fn selected_resolved_tiles_feed_the_field_while_the_escape_remainder_keeps_its_work() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let gas = reference();
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    let patch = NsfCharacteristicPatch { start: [0, 0], end: [1, 2] };
    let background = NsfBoundaryDrive::PressureOffset(-0.25);
    let make_zone = |high: bool| {
        let sign = if high { 1.0 } else { -1.0 };
        let face = NsfCharacteristicBoundary {
            reference_pressure_pa: gas.pressure, reference_temperature_k: gas.temperature,
            reference_velocity_m_s: [0.0; 3], incoming_characteristic_m_s: 0.0,
        };
        let mut fluid = NsfGrid3::<D3q27>::new(NsfSpec {
            dims: [8, 2, 2], periodic: [false; 3], gas: GasSpec::dry_air_ussa1976(),
            temperature_k: gas.temperature, pressure_pa: gas.pressure, dx_m: dx,
            temperature_window_k: [250.0, 600.0], maximum_mach: 0.2,
            maximum_temperature_ratio: 1.0, bulk_to_shear_ratio: 0.73,
            bulk_ratio_source: "authored complete-field patch transaction fixture",
            maximum_cells: CELLS, wall: WallThermal::Adiabatic, workers: 1, transport_points: 17,
        }, |_| false, |x| {
            let distance = if high { 7 - x[0] } else { x[0] };
            let z = (distance as f64 - 3.0) / 4.0;
            let pulse = 1e-5 * (0.7 + 0.1 * x[1] as f64) * fs_math::det::exp(-z * z);
            (1.0 + pulse / gas.gamma,
                [sign * (gas.gamma / 3.0).sqrt() * pulse / gas.gamma, 0.0, 0.0],
                gas.temperature * (1.0 + (gas.gamma - 1.0) / gas.gamma * pulse))
        }).unwrap().with_characteristic_axis(NsfCharacteristicAxis {
            axis: 0, faces: [face; 2], maximum_relative_perturbation: 0.05,
        }).unwrap();
        fluid.set_boundary_drives([background; 2]).unwrap();
        fluid.with_characteristic_patch(high, patch, background).unwrap()
    };
    let area = 2.0 * dx * dx;
    let source = field_source(&mesh, &unit, area.sqrt(), gas);
    let make_load = || ThermoviscousDgLoad::new(acoustic(&source, area,
        1.0 / f64::from(RATE), DgBoundary::Closed), 2, &gate).unwrap();
    for high in [false, true] {
        assert!(NsfThermoviscousDg::new(make_zone(high), make_load(), high,
            limits(1), &gate).is_err());
        let mut owner = NsfThermoviscousDg::new_on_patch(make_zone(high), make_load(),
            high, limits(1), &gate).unwrap();
        let mut original = make_zone(high);
        let mut escaping_flow: f64 = 0.0;
        for pressure in [2.0, 0.0, -1.0] {
            let mass = owner.zone().total_mass_kg();
            let energy = owner.zone().total_energy_j();
            let boundary = NsfBoundaryDrive::PressureOffset(pressure);
            let frame = owner.step_with_boundary(boundary, &gate).unwrap();
            let selected = frame.coupled_patch.unwrap();
            assert_eq!(frame.field_exchange, selected);
            assert_eq!(frame.inlet_flow_m3_s, selected.normal_volume_flow_m3_s);
            assert_eq!(owner.zone().boundary_drives().unwrap()[usize::from(high)], background);
            let mut drives = [boundary; 2]; drives[usize::from(high)] = background;
            let replay = original.preview_boundary_interval_with_patch(drives,
                NsfBoundaryDrive::IncomingCharacteristic(frame.incoming_characteristic_m_s),
                1, CELLS, || false).unwrap();
            assert_eq!(replay.patch_exchange(), Some(&selected));
            assert_eq!(replay.boundary_exchange(), &frame.faces);
            replay.accept();
            assert_eq!(original.state_hash(), owner.zone().state_hash());
            let complete = frame.faces[usize::from(high)];
            let complement_flow = complete.normal_volume_flow_m3_s - selected.normal_volume_flow_m3_s;
            let complement_work = complete.acoustic_pressure_work_out_j - selected.acoustic_pressure_work_out_j;
            let expected = -0.25 * complement_flow / f64::from(RATE);
            let scale = complete.acoustic_pressure_work_out_j.abs()
                + selected.acoustic_pressure_work_out_j.abs() + expected.abs();
            assert!((complement_work - expected).abs() <= 256.0 * f64::EPSILON * scale.max(1e-30));
            assert_eq!(frame.acoustic_work_residual_j,
                selected.acoustic_pressure_work_out_j - frame.inlet_work_j);
            assert!((frame.flow_mass_change_kg + frame.faces.iter().map(|f| f.mass_out_kg).sum::<f64>())
                .abs() < 3e-12 * mass);
            assert!((frame.flow_energy_change_j + frame.faces.iter().map(|f| f.total_energy_out_j).sum::<f64>())
                .abs() < 3e-12 * energy);
            escaping_flow = escaping_flow.max(complement_flow.abs());
        }
        assert!(escaping_flow > 1e-15, "the surrounding opening must carry separate flow");
    }
}

fn playback_field<'c>(mesh: &'c TetComplex, unit: &[[f64; 3]]) -> ThermoviscousInterior<'c> {
    let gas = reference();
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    field_source(mesh, unit, dx, gas)
}

fn playback_system<'a, 'c>(source: &'a ThermoviscousInterior<'c>, gate: &CancelGate)
    -> NsfThermoviscousDg<'a, 'c, D3q27>
{
    let fluid = zone(GasSpec::dry_air_ussa1976(), *source.gas(), true, 1);
    let area = fluid.law().dx_m.powi(2);
    let load = ThermoviscousDgLoad::new(acoustic(source, area, 1.0 / f64::from(RATE),
        DgBoundary::Closed), 2, gate).unwrap();
    NsfThermoviscousDg::new(fluid, load, true, limits(1), gate).unwrap()
}

#[test]
fn selected_duct_playback_matches_original_intervals_and_callback_bits() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let source = playback_field(&mesh, &unit);
    let controls: Vec<_> = (0..12).map(|sample| NsfDgDrive::Boundary(
        NsfBoundaryDrive::PressureOffset(if sample < 6 { 2.0 } else { 0.0 }))).collect();
    let make = |block| NsfThermoviscousDgPerformance::new_with_drives(
        playback_system(&source, &gate), 0, ApertureObservation::Inlet, RATE,
        controls.clone(), block, Arc::clone(&gate)).unwrap();
    let mut whole = make(12);
    let mut parts = make(5);
    let mut direct = playback_system(&source, &gate);
    let mut expected = [0.0; 12];
    for (output, &drive) in expected.iter_mut().zip(&controls) {
        direct.step_with_drive(drive, &gate).unwrap();
        let frame = direct.load().system().last_frame().unwrap();
        *output = frame.additional_ducts.iter().find(|duct| duct.port == 0)
            .unwrap().mean_pressure_pa;
    }
    let mut a = [0.0; 12]; let mut b = [0.0; 12];
    whole.block(&mut a).unwrap();
    for chunk in b.chunks_mut(5) { parts.block(chunk).unwrap(); }
    assert_eq!(a.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(b.map(f64::to_bits), expected.map(f64::to_bits));
    assert_eq!(snapshot!(whole.system()), snapshot!(parts.system()));
    assert_eq!(snapshot!(whole.system()), snapshot!(&direct));
    assert_eq!(whole.observation_port(), 0);
    assert_eq!(whole.drive_history(), controls.as_slice());
    assert_eq!(whole.samples_rendered(), 12);
    assert_eq!(whole.last_frame().unwrap().step, 12);
    assert_eq!(whole.remaining_samples(), 0);
    assert!(whole.validate_sample_count(1).is_err());
}

#[test]
fn selected_duct_output_clock_averages_original_intervals_without_retiming_physics() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let source = playback_field(&mesh, &unit);
    let controls = vec![NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(1.0)); 8];
    let mut direct = playback_system(&source, &gate);
    let original: Vec<_> = controls.iter().map(|&drive| {
        direct.step_with_drive(drive, &gate).unwrap();
        direct.load().system().last_frame().unwrap().additional_ducts[0].mean_pressure_pa
    }).collect();
    let mut renderer = NsfThermoviscousDgPerformance::new_resampled_with_drives(
        playback_system(&source, &gate), 0, ApertureObservation::Inlet, RATE / 2,
        controls, 4, 3, Arc::clone(&gate)).unwrap();
    let mut output = [0.0; 4];
    renderer.block(&mut output[..1]).unwrap();
    renderer.block(&mut output[1..]).unwrap();
    for (&actual, pair) in output.iter().zip(original.chunks_exact(2)) {
        let expected = 0.5 * pair[0] + 0.5 * pair[1];
        let scale = pair[0].abs() + pair[1].abs() + actual.abs();
        assert!((actual - expected).abs() <= 64.0 * f64::EPSILON * scale.max(1e-20));
    }
    assert_eq!(snapshot!(renderer.system()), snapshot!(&direct));
    assert_eq!(renderer.samples_rendered(), 4);
    assert_eq!(renderer.system().accepted_steps(), 8);
    assert_eq!(renderer.remaining_samples(), 0);
}

#[test]
fn cancellation_refuses_cached_output_and_preserves_the_completed_prefix_and_all_owners() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let source = playback_field(&mesh, &unit);
    let mut renderer = NsfThermoviscousDgPerformance::new_resampled_with_drives(
        playback_system(&source, &gate), 0, ApertureObservation::Inlet, RATE * 3,
        vec![NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(1.0)); 2],
        6, 3, Arc::clone(&gate)).unwrap();
    let mut prefix = [0.0];
    renderer.block(&mut prefix).unwrap();
    let prefix_bits = prefix[0].to_bits();
    let before = snapshot!(renderer.system());
    assert_eq!(renderer.system().accepted_steps(), 1);
    gate.request();
    // Both requested outputs would reuse that already accepted interval.
    let mut suffix = [123.0, 456.0];
    assert!(renderer.block(&mut suffix).is_err());
    assert_eq!(suffix, [123.0, 456.0]);
    assert_eq!(prefix[0].to_bits(), prefix_bits);
    assert_eq!(snapshot!(renderer.system()), before);
    assert_eq!(renderer.samples_rendered(), 1);
    assert_eq!(renderer.remaining_samples(), 5);
    assert!(matches!(renderer.block(&mut suffix), Err(crate::render::RenderError::Poisoned)));
}

#[test]
fn playback_admission_requires_an_attached_observation_and_complete_physical_controls() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let source = playback_field(&mesh, &unit);
    let valid = NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(1.0));
    for (port, rate, controls) in [
        (2, RATE, vec![valid]),
        (0, RATE, Vec::new()),
        (0, RATE, vec![NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(f64::NAN))]),
        (0, RATE, vec![NsfDgDrive::UpstreamPressureRise(1.0)]),
        (0, RATE / 2, vec![valid]),
    ] {
        assert!(NsfThermoviscousDgPerformance::new_with_drives(
            playback_system(&source, &gate), port, ApertureObservation::Inlet, rate,
            controls, 4, Arc::clone(&gate)).is_err());
    }
    assert!(NsfThermoviscousDgPerformance::new_resampled_with_drives(
        playback_system(&source, &gate), 0, ApertureObservation::Inlet, RATE / 2,
        vec![valid; 2], 4, 4, Arc::clone(&gate)).is_err());
    let mut advanced = playback_system(&source, &gate);
    advanced.step_with_drive(valid, &gate).unwrap();
    assert!(NsfThermoviscousDgPerformance::new_with_drives(advanced, 0,
        ApertureObservation::Inlet, RATE, vec![valid], 4, gate).is_err());
}

// G0/G3/G4: finite material sources retain their original storage while the
// complete acoustic field and both finite passages feed back into the fluid.
// The small chambers and Cartesian channel are manufactured coupling fixtures.
fn finite_medium(thermochemical: bool) -> (NsfGasMedium, GasState) {
    if thermochemical {
        let source = ThermochemicalNsfGas::new(
            super::source::model(-4000.0, 100500.0, 1.0),
            super::source::point(300.0, 100000.0),
            GasEnergySolve { minimum_temperature_k: 290.0, maximum_temperature_k: 320.0,
                energy_tolerance_j_kg: 1e-6, maximum_iterations: 64 },
        ).unwrap();
        let gas = *source.model().evaluate(source.context()).unwrap().state();
        (source.into(), gas)
    } else {
        use fs_material::gas_species::input::{GasRecord, read_gas_record};
        use fs_qty::semantic::{Composition, CompositionBasis};
        let GasRecord::Mixture(gas, parts) = read_gas_record(&mut
            "mixture 300 100000 mole 3 N2 0.78 O2 0.21 Ar 0.01"
                .split_ascii_whitespace(), None, None).unwrap()
            else { panic!("fixture requires its original named mixture"); };
        let source = fs_lbm::nsf3::NamedMixtureNsfGas::new(
            parts.iter().map(|part| part.0).collect(),
            Composition::new(CompositionBasis::MoleFraction,
                parts.iter().map(|part| part.1).collect()).unwrap(),
            gas.temperature, gas.pressure,
        ).unwrap();
        (source.into(), gas)
    }
}

fn finite_zone(medium: NsfGasMedium, gas: GasState) -> NsfGrid3<D3q27> {
    let face = NsfCharacteristicBoundary {
        reference_pressure_pa: gas.pressure, reference_temperature_k: gas.temperature,
        reference_velocity_m_s: [0.0; 3], incoming_characteristic_m_s: 0.0,
    };
    NsfGrid3::new(NsfSpec {
        dims: [12, 3, 3], periodic: [false; 3], gas: medium,
        temperature_k: gas.temperature, pressure_pa: gas.pressure,
        dx_m: (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE),
        temperature_window_k: [290.0, 320.0], maximum_mach: 0.2,
        maximum_temperature_ratio: 1.0, bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "manufactured finite-source composition fixture",
        maximum_cells: 108, wall: WallThermal::Adiabatic, workers: 1, transport_points: 17,
    }, |_| false, |_| (1.0, [0.0; 3], gas.temperature)).unwrap()
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis: 0, faces: [face; 2], maximum_relative_perturbation: 0.05,
        }).unwrap()
        .with_reservoir_aperture(fs_lbm::nsf3::NsfReservoirAperture {
            axis: 1, high: false, cells: vec![[5, 0, 1]],
            state: fs_lbm::nsf3::NsfOpenState {
                density_kg_m3: (gas.pressure + 2.0)
                    / (gas.specific_gas_constant * gas.temperature),
                velocity_m_s: [0.0; 3], temperature_k: gas.temperature,
            },
        }).unwrap()
}

fn finite_limits() -> NsfFiniteReservoirLimits {
    NsfFiniteReservoirLimits {
        maximum_relative_mass_exchange: 0.1, maximum_relative_energy_exchange: 0.1,
        mass_absolute_tolerance_kg: 1e-16, energy_absolute_tolerance_j: 1e-9,
        relative_tolerance: 1e-9,
    }
}

fn finite_field_system<'a, 'c>(field: &'a ThermoviscousInterior<'c>,
    medium: NsfGasMedium, high: bool, gate: &CancelGate)
    -> NsfThermoviscousDg<'a, 'c, D3q27>
{
    let fluid = finite_zone(medium, *field.gas());
    let area = 9.0 * fluid.law().dx_m.powi(2);
    let mut budget = limits(1);
    budget.maximum_cells = 108;
    budget.maximum_cell_steps = 108 * budget.maximum_evaluations;
    NsfThermoviscousDg::new(fluid, ThermoviscousDgLoad::new(
        acoustic(field, area, 1.0 / f64::from(RATE), DgBoundary::Closed), 2, gate,
    ).unwrap(), high, budget, gate).unwrap()
}

#[test]
fn finite_material_source_replays_the_full_field_and_conserves_original_heat_and_work() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    for thermochemical in [false, true] {
        let (medium, gas) = finite_medium(thermochemical);
        let side = 3.0 * finite_zone(medium.clone(), gas).law().dx_m;
        let field = field_source(&mesh, &unit, side, gas);
        for high in [false, true] {
            let mut owner = finite_field_system(&field, medium.clone(), high, &gate)
                .with_finite_reservoir(0, 1e-5, finite_limits(), &gate).unwrap();
            let mut replay = finite_field_system(&field, medium.clone(), high, &gate);
            for (rate, heat) in [(0.0, 0.0), (-2e-7, 0.01), (2e-7, -0.005), (0.0, 0.0)] {
                let before = owner.finite_reservoir_inventory().unwrap();
                let drive = if rate == 0.0 && heat == 0.0 {
                    NsfDgDrive::FiniteReservoir { aperture: 0 }
                } else {
                    NsfDgDrive::FiniteReservoirActuation {
                        aperture: 0, volume_rate_m3_s: rate, heat_power_w: heat,
                        maximum_relative_volume_change: 0.01,
                    }
                };
                let frame = owner.step_with_drive(drive, &gate).unwrap();
                let actual = frame.finite_reservoir.as_ref().unwrap();
                let after = owner.finite_reservoir_inventory().unwrap();
                // The source holds its accepted material state while the same
                // NSF + complete-field solve runs. Replaying that state through
                // the prescribed primitive must reproduce every physical bit.
                let original = replay.step_with_reservoir_aperture(0,
                    before.gas.pressure, before.gas.temperature, &gate).unwrap();
                assert_eq!(snapshot!(&owner), snapshot!(&replay));
                assert_eq!(frame.field_exchange, original.field_exchange);
                assert_eq!(frame.reservoir_exchanges, original.reservoir_exchanges);
                assert_eq!(actual.source_exchange, frame.reservoir_exchanges[0]);
                let (heat_j, work_j) = actual.actuation.map_or((0.0, 0.0), |receipt|
                    (receipt.step.supplied_heat_j, receipt.step.work_by_gas_j));
                assert!((after.mass_kg - before.mass_kg
                    - actual.source_exchange.mass_out_kg).abs() <= 1e-16);
                assert!((after.sensible_energy_j - before.sensible_energy_j
                    - actual.source_exchange.total_energy_out_j - heat_j + work_j).abs() <= 1e-9);
                assert_eq!(actual.initial_pressure_pa, before.gas.pressure);
                assert_eq!(actual.final_pressure_pa, after.gas.pressure);
                assert_eq!(actual.final_mass_kg, after.mass_kg);
                assert!(actual.mass_residual_kg.abs() <= 1e-16);
                assert!(actual.energy_residual_j.abs() <= 1e-9);
                let dt = owner.load().duration_s();
                assert!((after.volume_m3 - before.volume_m3 - rate * dt).abs()
                    <= 4.0 * f64::EPSILON * before.volume_m3);
                assert_eq!(heat_j, heat * dt);
                if rate != 0.0 { assert_eq!(work_j.is_sign_positive(), rate.is_sign_positive()); }
                match owner.zone().medium() {
                    NsfGasMedium::Thermochemical(source) => {
                        assert_eq!(after.source_kind, NsfFiniteReservoirSource::Thermochemical);
                        assert_eq!(owner.finite_reservoir().unwrap().model(), source.model());
                    }
                    NsfGasMedium::NamedMixture(source) => {
                        assert_eq!(after.source_kind, NsfFiniteReservoirSource::NamedMixture);
                        assert!(owner.finite_reservoir().is_none());
                        assert_eq!(after.gas, source.at(after.gas.temperature, after.gas.pressure).unwrap());
                    }
                    _ => panic!("finite source lost its original material owner"),
                }
            }
        }
    }
}

#[test]
fn finite_source_failure_and_cancellation_preserve_chamber_fluid_field_and_every_duct() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let side = 3.0 * finite_zone(medium.clone(), gas).law().dx_m;
    let field = field_source(&mesh, &unit, side, gas);
    let mut owner = finite_field_system(&field, medium, true, &gate)
        .with_finite_reservoir(0, 1e-5, finite_limits(), &gate).unwrap();
    owner.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).unwrap();
    let before = (snapshot!(&owner), owner.finite_reservoir_inventory(),
        owner.zone().reservoir_aperture_at(0).unwrap().clone(),
        owner.zone().last_reservoir_exchanges().to_vec());
    // The rate is finite and belongs to the right owner, so its oversized
    // physical volume increment is refused only after fluid/acoustic trials.
    let invalid = NsfDgDrive::FiniteReservoirActuation {
        aperture: 0, volume_rate_m3_s: -1.0, heat_power_w: 0.0,
        maximum_relative_volume_change: 0.01,
    };
    assert!(owner.validate_drive(invalid).is_ok());
    assert!(matches!(owner.step_with_drive(invalid, &gate), Err(NsfDgError::Input(
        "finite reservoir actuation loses an increment or exceeds its positive bounded volume movement; refine the interval"))));
    assert_eq!(before, (snapshot!(&owner), owner.finite_reservoir_inventory(),
        owner.zone().reservoir_aperture_at(0).unwrap().clone(),
        owner.zone().last_reservoir_exchanges().to_vec()));
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(owner.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 },
        &cancelled), Err(NsfDgError::Cancelled)));
    assert_eq!(before, (snapshot!(&owner), owner.finite_reservoir_inventory(),
        owner.zone().reservoir_aperture_at(0).unwrap().clone(),
        owner.zone().last_reservoir_exchanges().to_vec()));
    for invalid in [
        NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(1.0)),
        NsfDgDrive::ReservoirAperture { aperture: 0, pressure_pa: 100001.0, temperature_k: 300.0 },
        NsfDgDrive::FiniteReservoir { aperture: 1 },
    ] { assert!(owner.validate_drive(invalid).is_err()); }
    assert_eq!(owner.accepted_steps(), 1);
}

#[test]
fn finite_source_feedback_reaches_playback_without_callback_dependence() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let side = 3.0 * finite_zone(medium.clone(), gas).law().dx_m;
    let field = field_source(&mesh, &unit, side, gas);
    let make = || finite_field_system(&field, medium.clone(), true, &gate)
        .with_finite_reservoir(0, 1e-5, finite_limits(), &gate).unwrap();
    let initial = make().finite_reservoir_inventory().unwrap();
    let controls = vec![NsfDgDrive::FiniteReservoirActuation {
        aperture: 0, volume_rate_m3_s: -2e-7, heat_power_w: 0.01,
        maximum_relative_volume_change: 0.01,
    }; 24];
    let render = |system, controls| NsfThermoviscousDgPerformance::new_with_drives(
        system, 1, ApertureObservation::Inlet, RATE, controls, 24, Arc::clone(&gate),
    ).unwrap();
    let mut whole = render(make(), controls.clone());
    let mut parts = render(make(), controls);
    let mut held = render(finite_field_system(&field, medium.clone(), true, &gate),
        vec![NsfDgDrive::ReservoirAperture { aperture: 0,
            pressure_pa: initial.gas.pressure, temperature_k: initial.gas.temperature }; 24]);
    let mut a = [0.0; 24]; let mut b = [0.0; 24]; let mut fixed = [0.0; 24];
    whole.block(&mut a).unwrap(); held.block(&mut fixed).unwrap();
    for chunk in b.chunks_mut(5) { parts.block(chunk).unwrap(); }
    assert_eq!(a, b);
    assert_eq!(a[0], fixed[0]);
    assert_ne!(a, fixed, "finite source work/heat must reach the original pressure observation");
    assert_eq!(snapshot!(whole.system()), snapshot!(parts.system()));
    assert_eq!(whole.system().finite_reservoir_inventory(), parts.system().finite_reservoir_inventory());
    assert_ne!(whole.system().finite_reservoir_inventory().unwrap().gas.pressure, initial.gas.pressure);
    assert_eq!(whole.samples_rendered(), 24);
    assert_eq!(whole.system().accepted_steps(), 24);
    assert!(whole.last_frame().unwrap().finite_reservoir.as_ref().unwrap().actuation.is_some());
}

#[test]
fn complete_field_finite_source_requires_a_cold_owned_material_aperture() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let side = 3.0 * finite_zone(medium.clone(), gas).law().dx_m;
    let field = field_source(&mesh, &unit, side, gas);
    let make = || finite_field_system(&field, medium.clone(), true, &gate);
    assert!(make().validate_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }).is_err());
    for (aperture, volume) in [(1, 1e-5), (0, 0.0), (0, f64::NAN)] {
        assert!(make().with_finite_reservoir(aperture, volume, finite_limits(), &gate).is_err());
    }
    let owner = make().with_finite_reservoir(0, 1e-5, finite_limits(), &gate).unwrap();
    assert_eq!(owner.finite_reservoir_aperture(), Some(0));
    assert_eq!(owner.finite_reservoir_limits(), Some(&finite_limits()));
    assert!(owner.with_finite_reservoir(0, 1e-5, finite_limits(), &gate).is_err());
    let mut advanced = make();
    advanced.step(0.0, &gate).unwrap();
    assert!(advanced.with_finite_reservoir(0, 1e-5, finite_limits(), &gate).is_err());
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(make().with_finite_reservoir(0, 1e-5, finite_limits(), &cancelled),
        Err(NsfDgError::Cancelled)));
    let legacy = reference();
    let medium = NsfGasMedium::from(GasSpec::dry_air_ussa1976());
    let side = 3.0 * finite_zone(medium.clone(), legacy).law().dx_m;
    let legacy_field = field_source(&mesh, &unit, side, legacy);
    assert!(finite_field_system(&legacy_field, medium, true, &gate)
        .with_finite_reservoir(0, 1e-5, finite_limits(), &gate).is_err());
}
