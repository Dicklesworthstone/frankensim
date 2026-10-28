//! Native cards and mechanical key programs retain the original whole source.
use super::*;
use crate::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture};
use crate::bernoulli_aperture::network::boundaries::{BoundaryAperture, BoundaryDrive};
use crate::bernoulli_aperture::two_sided::{TwoSidedAperture, TwoSidedDrive};
use crate::pcm_wav::baffled::openings::duct::{DuctPressure, DuctView};
use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue};

const KEYED: &str = include_str!("../../../../../../examples/material_modal_valve/paired-keyed-sampled.performance");
const SIDES: [ValveSide; 2] = [ValveSide::Upstream, ValveSide::Downstream];

fn source() -> String {
    KEYED.replace("audio 768000 48000 2400 0.005", "audio 768000 48000 256 0.005")
}
fn load(text: &str, block: usize) -> Result<MaterialModalValvePerformance, PlateValveInputError> {
    MaterialModalValvePerformance::from_bytes_with_sources(text.as_bytes(), block,
        &CancelGate::new_clock_free(), None,
        |path| Ok(pack(if path == "key.matpack" { 1800.0 } else { 900.0 }, false)))
}

#[test]
fn g0_native_keys_share_original_card_cache_and_actual_lowered_chimney_medium() {
    let mut calls = 0;
    let performance = MaterialModalValvePerformance::from_bytes_with_sources(source().as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |path| {
            assert_eq!(path, "blade.matpack"); calls += 1; Ok(pack(900.0, false))
        }).unwrap();
    assert_eq!(calls, 1, "all reed regions and both keys share one immutable source cache");
    assert_eq!(performance.material_records().len(), 4);
    assert_eq!(performance.boundary_records().len(), 2);
    assert_eq!(valve(&performance).modal_states().len(), 6);
    assert_eq!(valve(&performance).viscosity_matrix_n_s_m().unwrap().len(), 36);
    let system = performance.renderer().source().system();
    assert!(system.supply_channel().is_some());
    assert!(system.graph(ValveSide::Downstream).sections.len()
        > system.source(ValveSide::Downstream).graph.sections.len());
    for (side, record) in performance.boundary_records() {
        assert!(std::ptr::eq(record.pack(), performance.material_records()[0].pack()));
        let keys = system.boundary_apertures(*side).unwrap();
        let key = keys.aperture(record.node()).unwrap();
        let section = system.graph(*side).sections.iter()
            .position(|section| section.nodes.contains(&record.node())).unwrap();
        assert_eq!(key.spec().density_kg_m3, system.section_gases(*side)[section].density);
        assert_eq!(key.spec().impedance_pa_s_m3,
            system.represented_sections(*side)[section].impedance_pa_s_m3);
        assert_eq!(key.viscous_slit().unwrap().gas(), Some(&system.section_gases(*side)[section]));
        assert_eq!(key.spec().time_step_s.to_bits(), valve(&performance).spec().time_step_s.to_bits());
        assert_eq!(key.spec().max_steps, 256 * 16);
        assert_eq!(key.accepted_steps(), 0);
        assert!(key.plate_reduction().unwrap().material_chart().is_some());
        assert_eq!(performance.renderer().pending_boundary_controls(*side).len(), 2);
        assert_eq!(performance.renderer().held_boundary_drives(*side), &[BoundaryDrive {
            node: 3, exterior_pressure_pa: 0.0, opening_force_n: -0.0,
        }]);
    }
    assert_eq!(system.boundary_apertures(ValveSide::Downstream).unwrap().gas(3).unwrap().temperature, 313.15);
}

#[test]
fn g3_only_key_material_or_local_gas_changes_reach_reed_pressure_and_pcm() {
    let text = source();
    let base = load(&text, 37).unwrap();
    let specimens: Vec<_> = valve(&base).material_plates().iter()
        .map(|plate| plate.source().specimen_identity()).collect();
    let (pressure, state) = render(base);
    assert_eq!(pressure, render(load(&text, 1).unwrap()).0);
    assert_eq!(pressure, render(load(&text, 128).unwrap()).0);
    let pcm = encode_pcm16_wav(&pressure, 48_000, 0.005).unwrap().0;
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pcm[44..].iter().any(|byte| *byte != 0));
    for side in SIDES {
        assert!(state.pending_boundary_controls(side).is_empty());
        assert_eq!(state.applied_boundary_controls(side).len(), 2);
        assert_eq!(state.source().system().boundary_apertures(side).unwrap()
            .aperture(3).unwrap().accepted_steps(), 256 * 16);
    }
    assert!(state.source().last_frame().unwrap().balance_residual_j().abs() < 1e-12);
    for changed in [
        text.replace("boundary_plate blade.matpack", "boundary_plate key.matpack"),
        text.replace("gas 313.15 0 viscothermal", "gas 293.15 0 viscothermal"),
        text.replace("boundary_control 321 3 0.003 0.3", "boundary_control 321 3 0 0.3"),
        text.replace("boundary_control 321 3 0.003 0.3", "boundary_control 321 3 0.003 0"),
    ] {
        let changed = load(&changed, 37).unwrap();
        assert_eq!(valve(&changed).material_plates().iter()
            .map(|plate| plate.source().specimen_identity()).collect::<Vec<_>>(), specimens);
        let (other, other_state) = render(changed);
        assert_ne!(pressure, other);
        assert_ne!(pcm, encode_pcm16_wav(&other, 48_000, 0.005).unwrap().0);
        assert_ne!(state.source().last_frame().unwrap().aperture.state,
            other_state.source().last_frame().unwrap().aperture.state);
    }
}

#[test]
fn key_maxwell_memory_and_supplied_mutual_radiation_survive_the_same_native_path() {
    let text = source()
        .replace("boundary_loss elastic", concat!("boundary_loss maxwell 1 0 1000 0.1 0.1\n",
            "boundary_term boundary_modulus boundary_time"))
        .replace("downstream\nnetwork 4 3 1048576", "downstream\nnetwork 5 4 1048576")
        .replace("duct_node reflection 0\nduct_node driven", "duct_node reflection 0\nduct_node driven\nduct_node reflection 0")
        .replace("propagation fractional", "duct_section 1 4 0.12 0.007 0.002\npropagation fractional")
        .replace("openings 1\nopening downstream 2 0 0", concat!("openings 2\nopening downstream 2 0 0\n",
            "opening downstream 4 0.03 0\nradiation mimo 1\n",
            "radiation_group downstream mutual.radiation 256 8388608 65"));
    let artifact = include_str!("../../../../../../examples/material_modal_valve/mutual.radiation");
    let build = |block| MaterialModalValvePerformance::from_bytes_with_all_sources(text.as_bytes(), block,
        &CancelGate::new_clock_free(), None, |_| Ok(pack(900.0, false)),
        |_| Ok(artifact.as_bytes().to_vec())).unwrap();
    let original = build(37);
    assert_eq!(original.mimo_radiation_records().len(), 1);
    assert_eq!(original.boundary_records().len(), 2);
    let (pressure, rendered) = render(original);
    assert_eq!(pressure, render(build(1)).0);
    let system = rendered.source().system();
    assert_eq!(system.mimo_environment(ValveSide::Downstream).unwrap().radiation().openings
        .iter().map(|opening| opening.node).collect::<Vec<_>>(), vec![2, 4]);
    assert!(system.mimo_pole_states(ValveSide::Downstream).unwrap().iter()
        .any(|state| state.re != 0.0 || state.im != 0.0));
    assert!(system.supply_channel().is_some());
    assert!(system.fractional_sections(ValveSide::Downstream).is_some());
    assert_eq!(system.viscothermal_sections(ValveSide::Downstream).len(), 2);
    for side in SIDES {
        let key = system.boundary_apertures(side).unwrap().aperture(3).unwrap();
        assert_eq!(key.relaxation().unwrap().branches().len(), 1);
        assert!(key.relaxation().unwrap().dissipated_energy_j() > 0.0);
    }
    assert!(rendered.source().last_frame().unwrap().balance_residual_j().abs() < 1e-12);
}

pub(super) fn direct_system(source: &MaterialModalValvePerformance) -> TwoSidedAperture {
    let gate = CancelGate::new_clock_free();
    let original = source.renderer().source().system();
    let reed = valve(source);
    let spec = reed.spec();
    let plates = reed.material_plates().iter().map(|plate| (
        plate.source().clone(), plate.options().clone(), plate.viscosity().unwrap().clone(),
        crate::plate_viscosity::PlateViscosityLimits { max_modes: 32,
            max_projection_terms: 1_000_000, reference_resistance_n_s_m: 1.0 },
    )).collect();
    let contact = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0],
        1e8, 2.0, "direct manufactured reed".into()).unwrap().with_internal_loss(0.5).unwrap();
    let mut valve = ModalAperture::from_viscous_material_plates_with_propagator(
        plates, spec.rest_opening_m, *reed.gas().unwrap(), spec.impedance_pa_s_m3,
        spec.time_step_s, spec.max_steps, spec.limits, reed.modal_states().to_vec(), contact,
        reed.blade_propagator(), &gate).unwrap();
    valve.set_material_viscous_slit(0.001).unwrap();
    let channel = original.supply_channel().unwrap().slit().unwrap();
    let channel = SlitChannel::from_gas(channel.gas(), *channel.spec()).unwrap();
    let policies = SIDES.map(|side| original.fractional_sections(side).unwrap().0);
    let mut system = TwoSidedAperture::with_material_walls_and_supply_channel(
        valve.into(), SIDES.map(|side| original.source(side).clone()), [None, None], channel,
        TwoSidedLimits { maximum_energy_j: 1.0, energy_absolute_tolerance_j: 1e-16,
            energy_relative_tolerance: 1e-7, flow_absolute_tolerance_m3_s: 1e-16,
            flow_relative_tolerance: 1e-8 }, Some(policies), &gate).unwrap();
    for side in SIDES {
        let key = original.boundary_apertures(side).unwrap().aperture(3).unwrap();
        let mut aperture = DynamicAperture::from_plate_with_propagator(
            key.plate_reduction().unwrap().clone(), key.spec().density_kg_m3,
            key.spec().impedance_pa_s_m3, spec.time_step_s, spec.max_steps,
            ApertureState { opening_m: 0.0002, opening_velocity_m_s: 0.0 },
            key.contact_law().clone(), BladePropagator::ExactZoh).unwrap();
        aperture.set_viscous_slit(key.viscous_slit().unwrap()).unwrap();
        system = system.with_boundary_apertures(side,
            vec![BoundaryAperture { node: 3, aperture }], 1, &gate).unwrap();
    }
    system
}

#[test]
fn g1_native_subframe_key_signs_and_events_match_direct_physics_and_receiver() {
    let phrase = GestureSchedule::try_new(1, vec![GestureTrack { id: "mouth".into(),
        target: GestureTarget::BlowingPressure, initial: GestureValue::PressurePa(5.0),
        events: vec![] }]).unwrap();
    let text = format!("{}\nschedule\n{}", source().split_once("\nschedule\n").unwrap().0,
        String::from_utf8(phrase.to_canonical_bytes()).unwrap());
    let native = load(&text, 37).unwrap();
    let mut direct = direct_system(&native);
    let mut receiver = DuctPressure::new(&SIDES.map(|side| DuctView {
        graph: direct.graph(side), gases: direct.section_gases(side),
    }), native.renderer().source().radiation().clone()).unwrap();
    let mut expected = Vec::new();
    for sample in 0..4096 {
        let up = [BoundaryDrive { node: 3,
            exterior_pressure_pa: if (193..769).contains(&sample) { 0.1 } else { 0.0 },
            opening_force_n: if (193..769).contains(&sample) { -0.002 } else { -0.0 } }];
        let down = [BoundaryDrive { node: 3,
            exterior_pressure_pa: if (321..1505).contains(&sample) { 0.3 } else { 0.0 },
            opening_force_n: if (321..1505).contains(&sample) { -0.003 } else { -0.0 } }];
        direct.step_with_boundaries(TwoSidedDrive { source_pressure_pa: 5.0,
            blade_forces_n: [0.0, 0.0] }, [&up, &down]).unwrap();
        if let Some(pressure) = receiver.push_frame(|network, node|
            direct.node_frame(SIDES[network], node).map(|frame| frame.net_flow_into_node_m3_s)).unwrap() {
            expected.push(pressure);
        }
    }
    let (actual, rendered) = render(native);
    assert_eq!(actual, expected);
    assert_eq!(encode_pcm16_wav(&actual, 48_000, 0.005).unwrap().0,
        encode_pcm16_wav(&expected, 48_000, 0.005).unwrap().0);
    for side in SIDES {
        assert_eq!(direct.boundary_apertures(side).unwrap().aperture(3).unwrap().state(),
            rendered.source().system().boundary_apertures(side).unwrap().aperture(3).unwrap().state());
    }
}

#[test]
fn incomplete_keys_wrong_controls_and_missing_material_refuse_native_admission() {
    let text = source();
    for invalid in [
        text.replace("boundary_side downstream", "boundary_side upstream"),
        text.replace("boundary 3 instrument bore-key", "boundary 2 instrument bore-key"),
        text.replace("boundary 3 instrument bore-key", "boundary 3 instrument tract-key"),
        text.replace("boundary_control 321 3 0.003 0.3", "boundary_control 321 2 0.003 0.3"),
        text.replace("boundary_control 321 3 0.003 0.3", "boundary_control 4096 3 0.003 0.3"),
        text.replace("boundary_control 321 3 0.003 0.3", "boundary_control 321 3 -0.003 0.3"),
        text.replace("boundary_control 1505 3 0 0", "boundary_control 321 3 0 0"),
        text.replace("boundary_plate blade.matpack", "boundary_plate ../blade.matpack"),
        text.replace("compile_limits 100000 8192", "compile_limits 1 8192"),
        text.replace("compile_limits 100000 8192", "compile_limits 100000 1"),
        text.replace("compile_limits 100000 8192", "control_source instrument reed\ncompile_limits 100000 8192"),
        text.replace("compile_limits 100000 8192", "upstream_boundaries\nboundaries none\ncompile_limits 100000 8192"),
    ] { assert!(load(&invalid, 37).is_err()); }
    let (prefix, suffix) = text.split_once("boundary_groups 2\n").unwrap();
    let (_, controls) = suffix.split_once("compile_limits ").unwrap();
    assert!(load(&format!("{prefix}compile_limits {controls}"), 37).is_err());
    let missing = text.replace("boundary_plate blade.matpack", "boundary_plate absent.matpack");
    assert!(matches!(MaterialModalValvePerformance::from_bytes_with_sources(missing.as_bytes(), 37,
        &CancelGate::new_clock_free(), None, |path| if path == "absent.matpack" {
            Err("missing physical key material".into())
        } else { Ok(pack(900.0, false)) }), Err(PlateValveInputError::MaterialSource { .. })));
}
