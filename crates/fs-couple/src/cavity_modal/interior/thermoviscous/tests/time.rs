use super::*;
use super::super::dg::{ThermoviscousDg, ThermoviscousDgLimits};
use crate::dg_port::DgPortLoad;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, PortId, MAX_DG_PROPAGATOR_DOFS};
use fs_exec::CancelGate;

fn time_dg_load(gas: &GasState, duration_s: f64, outlet: DgBoundary) -> DgPortLoad {
    let duct = DgDuct::new(
        &fs_duct::Duct {
            segments: vec![Segment::Cylinder {
                radius: (HEIGHT * WIDTH / std::f64::consts::PI).sqrt(),
                length: 0.05,
            }],
        },
        gas,
        DgMesh { elements: 2, order: 2 },
        DgBoundary::Port(PortId(7)),
        outlet,
    )
    .unwrap();
    DgPortLoad::new(duct, duration_s, 1e-11, MAX_DG_PROPAGATOR_DOFS).unwrap()
}

#[test]
fn original_gas_source_reaches_causal_field_work_without_caloric_shortcuts() {
    let (mesh, points0) = kuhn_cube(1);
    let xyz = points(&points0);
    let gas = gas();
    let original =
        ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, true), gas, BAND).unwrap();
    let prepared = original
        .prepare_time(
            1.0 / 48_000.0,
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        prepared.inner().medium(),
        source_medium(original.gas(), 0.0).unwrap()
    );
    let step = prepared
        .step(&prepared.zero_state(), &[1.0, 0.0], &[], &mut || false)
        .unwrap();
    assert!(step.work.pressure_work_j > 0.0 && step.work.final_energy_j > 0.0);
    assert!(step.work.balance_residual_j().abs() < 1e-10 * step.work.pressure_work_j);
    assert!(
        original
            .prepare_time(
                1.0 / BAND[1],
                ThermoviscousTimeLimits::default(),
                &mut || false
            )
            .is_err()
    );
    let molecular = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, true),
        source(gas, true, 1.0),
        BAND,
    )
    .unwrap();
    let prepared = molecular
        .prepare_time(
            1.0 / 48_000.0,
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let source = molecular.relaxation().unwrap();
    assert_eq!(
        prepared.frozen_cp_j_kg_k(),
        source.frozen_cv_j_kg_k() + gas.specific_gas_constant
    );
    assert_eq!(
        prepared.caloric_reservoirs().len(),
        source.processes().len()
    );
    for (state, original) in prepared.caloric_reservoirs().iter().zip(source.processes()) {
        assert_eq!(
            state.heat_capacity_j_kg_k.to_bits(),
            original.heat_capacity_j_kg_k.to_bits()
        );
        assert_eq!(
            state.relaxation_time_s.to_bits(),
            original.relaxation_time_s.to_bits()
        );
    }
    // The same material law, including the explicit phasor convention change.
    let frequency = 400.0;
    let mut cp = C64::from_re(prepared.frozen_cp_j_kg_k());
    for state in prepared.caloric_reservoirs() {
        cp = cp
            + C64::new(
                1.0,
                std::f64::consts::TAU * frequency * state.relaxation_time_s,
            )
            .recip()
            .scale(state.heat_capacity_j_kg_k);
    }
    close(
        cp,
        molecular.heat_capacity(frequency, &mut || false).unwrap(),
        cp.abs(),
        4e-15,
    );
    let mut initial = prepared.zero_state();
    // Mesh vertex zero is an isothermal wall coefficient. Its molecular state
    // must be admitted independently of the constrained translational field.
    assert_eq!(initial.temperature_k[0], 0.0);
    initial.molecular_temperature_k[0][0] = 1.0;
    let step = prepared
        .step(&initial, &[0.0, 0.0], &[], &mut || false)
        .unwrap();
    assert!(step.work.molecular_relaxation_heat_j > 0.0);
    assert!(step.work.final_energy_j < step.work.initial_energy_j);
    assert!(step.work.balance_residual_j().abs() < 1e-10 * step.work.initial_energy_j);
}

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::{PreparedThermalTransport, ThermalTransportSourceQuery};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn wall_transport(rho: f64, cp: f64, k: f64) -> PreparedThermalTransport {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let keys = [
        PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0])),
        PropertyKey::new("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
        PropertyKey::new("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0])),
    ];
    let mut claims = ClaimSet::new();
    for (key, value) in keys.iter().zip([rho, cp, k]) {
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Scalar {
                    value,
                    dims: key.quantity().dims(),
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature,
                    280.0,
                    360.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic fixed-geometry radial thermal source".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic wall".into(),
            phase: "solid".into(),
            process: "authored numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let lower = QueryPoint::new()
        .with_quantity("T", temperature, 280.0)
        .unwrap();
    let upper = QueryPoint::new()
        .with_quantity("T", temperature, 360.0)
        .unwrap();
    PreparedThermalTransport::try_new(
        &card,
        ThermalTransportSourceQuery {
            density: &keys[0],
            specific_heat: &keys[1],
            conductivity: &keys[2],
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 3],
        },
    )
    .unwrap()
}

#[test]
fn source_wall_material_changes_actual_shared_time_fields_and_physical_heat() {
    use super::super::wall::MaterialWallRegion;
    use fs_feec::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
    use fs_feec::acoustics::helmholtz::boundary_layer::thermal_wall::{
        ThermalBacking, ThermalWallDepth,
    };
    let (mesh, points0) = kuhn_cube(1);
    let xyz = points(&points0);
    let gas = gas();
    let geometry = HelmholtzFe::new(
        &mesh,
        &xyz,
        2,
        Medium {
            density_kg_m3: gas.density,
            sound_speed_m_s: gas.sound_speed,
        },
    )
    .unwrap();
    let face = geometry
        .boundary_faces()
        .into_iter()
        .find(|f| f.iter().all(|&v| xyz[v as usize][1] == 0.0))
        .unwrap();
    let mut heats = Vec::new();
    let mut identities = Vec::new();
    for k in [0.003, 0.012] {
        let wall = wall_transport(1000.0, 2.0, k);
        let sample = wall.at_temperature(gas.temperature).unwrap();
        let source = ThermoviscousInterior::from_relaxation(
            build(&mesh, &xyz, &gas, false),
            source(gas, true, 1.0),
            BAND,
        )
        .unwrap()
        .with_conductive_walls(
            &[MaterialWallRegion {
                faces: &[face],
                sample: &sample,
                depth: ThermalWallDepth::Slab {
                    thickness_m: 0.00002,
                    backing: ThermalBacking::Adiabatic,
                },
                minimum_tangential_scale_m: 0.004,
                geometry: BoundaryLayerGeometry {
                    minimum_curvature_radius_m: f64::INFINITY,
                    minimum_half_gap_m: 0.0001,
                    maximum_layer_ratio: 0.05,
                },
            }],
            &mut || false,
        )
        .unwrap();
        assert_eq!(source.wall_sources(), &[sample]);
        assert!(
            source
                .prepare_time(1e-7, ThermoviscousTimeLimits::default(), &mut || false)
                .is_err()
        );
        let prepared = source
            .prepare_time_with_wall_diffusion(
                1e-7,
                &[2],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(prepared.maximum_wall_frequency_hz(), Some(BAND[1]));
        assert_eq!(
            prepared.inner().conductive_walls()[0]
                .wall
                .transport
                .conductivity_w_m_k,
            k
        );
        assert_eq!(
            prepared.inner().conductive_walls()[0]
                .wall
                .transport
                .volumetric_heat_capacity_j_m3_k,
            sample.volumetric_heat_capacity_j_m3_k()
        );
        assert_eq!(
            prepared.caloric_reservoirs().len(),
            source.relaxation().unwrap().processes().len()
        );
        let mut initial = prepared.zero_state();
        initial.conductive_wall_temperature_k[0][0] = 1.0;
        let step = prepared
            .step(&initial, &[0.0, 0.0], &[], &mut || false)
            .unwrap();
        assert!(step.work.solid_conductive_heat_j > 0.0);
        assert!(step.state.temperature_k.iter().any(|t| *t != 0.0));
        assert!(step.work.balance_residual_j().abs() < 1e-10 * step.work.initial_energy_j);
        // The original finite wall and molecular fields also survive actual
        // reciprocal bore loading; their heat must move pressure and bore state.
        let gate = CancelGate::new_clock_free();
        let mut coupled = ThermoviscousDg::new(
            &source,
            prepared,
            time_dg_load(&gas, 1e-7, DgBoundary::Closed),
            1,
            ThermoviscousDgLimits::default(),
            &gate,
        )
        .unwrap()
        .with_initial_state(initial)
        .unwrap();
        let frame = coupled.step(&[0.0, 0.0], &[], &gate).unwrap().clone();
        assert!(frame.field.work.solid_conductive_heat_j > 0.0);
        assert!(frame.field.state.temperature_k.iter().any(|t| *t != 0.0));
        assert!(frame.field.state.molecular_temperature_k[0].iter().any(|t| *t != 0.0));
        assert!(coupled.load().state().pressure().iter().any(|p| *p != 0.0));
        assert_eq!(coupled.prepared().caloric_reservoirs().len(), source.relaxation().unwrap().processes().len());
        assert!(frame.balance_residual_j.abs() <= frame.balance_tolerance_j);
        heats.push(step.work.solid_conductive_heat_j);
        identities.push(sample.path_identity());
    }
    assert_ne!(identities[0], identities[1]);
    assert!(heats[1] > 3.8 * heats[0]);
}

#[test]
fn resolved_molecular_bore_coupling_replays_both_original_owners_and_closes_work() {
    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let source = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, true), source(gas, true, 1.0), BAND,
    ).unwrap();
    let dt = 1.0 / 96_000.0;
    let prepared = source.prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false).unwrap();
    let mut initial = prepared.zero_state();
    initial.molecular_temperature_k[0][0] = 0.01;
    let gate = CancelGate::new_clock_free();
    let mut coupled = ThermoviscousDg::new(
        &source, prepared, time_dg_load(&gas, dt, DgBoundary::Closed), 1,
        ThermoviscousDgLimits::default(), &gate,
    ).unwrap().with_initial_state(initial).unwrap();
    let mut molecular_heat = 0.0;
    let mut inlet_work = 0.0;
    for drive in [1.0, 1.0, 0.0, -0.5, 0.0] {
        let original_field = coupled.state().clone();
        let mut original_bore = coupled.load().clone();
        let frame = coupled.step(&[drive, 0.0], &[], &gate).unwrap().clone();
        let replay = coupled.prepared().step(
            &original_field, &[drive, frame.mean_pressure_pa], &[], &mut || false,
        ).unwrap();
        let bore = original_bore.preview_step(frame.mean_flow_m3_s).unwrap().accept();
        assert_eq!(coupled.state(), &replay.state);
        assert_eq!(coupled.load().state(), original_bore.state());
        assert_eq!(frame.duct, bore);
        assert_eq!(frame.field.mean_inward_flow_m3_s, replay.mean_inward_flow_m3_s);
        let scale = frame.mean_flow_m3_s.abs().max(1e-30);
        assert!((frame.mean_flow_m3_s + replay.mean_inward_flow_m3_s[1]).abs() < 2e-12 * scale);
        let work = -frame.mean_pressure_pa * replay.mean_inward_flow_m3_s[1] * dt;
        let bore_work = frame.duct.ledger.port_work_step_j;
        assert!((work - bore_work).abs() <= 2e-12 * work.abs().max(1e-30));
        assert!(frame.interface_work_residual_j.abs() <= frame.balance_tolerance_j);
        assert!(frame.balance_residual_j.abs() <= frame.balance_tolerance_j);
        assert!((frame.storage_change_j + frame.material_dissipation_j + frame.terminal_work_j
            - frame.authored_external_work_j).abs() <= frame.balance_tolerance_j);
        molecular_heat += frame.field.work.molecular_relaxation_heat_j;
        inlet_work += work.abs();
    }
    assert!(molecular_heat > 0.0 && inlet_work > 0.0);
    assert_eq!(coupled.accepted_steps(), 5);
}

#[test]
fn resolved_field_responds_to_bore_termination_and_failed_steps_publish_nothing() {
    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let source = ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, true), gas, BAND).unwrap();
    let dt = 1.0 / 96_000.0;
    let gate = CancelGate::new_clock_free();
    let make = |outlet, limits| ThermoviscousDg::new(
        &source,
        source.prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false).unwrap(),
        time_dg_load(&gas, dt, outlet), 1, limits, &gate,
    ).unwrap();
    let mut closed = make(DgBoundary::Closed, ThermoviscousDgLimits::default());
    let mut matched = make(DgBoundary::Matched, ThermoviscousDgLimits::default());
    for sample in 0..64 {
        let drive = if sample < 16 { 1.0 } else { 0.0 };
        closed.step(&[drive, 0.0], &[], &gate).unwrap();
        matched.step(&[drive, 0.0], &[], &gate).unwrap();
    }
    let reaction = closed.state().pressure_pa.iter().zip(&matched.state().pressure_pa)
        .map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    assert!(reaction > 1e-7, "downstream termination must alter the original FE pressure: {reaction:e}");

    let accepted_field = closed.state().clone();
    let accepted_bore = closed.load().state().clone();
    let accepted_response = closed.load().port_response().unwrap();
    for traction in [[f64::NAN, 0.0], [0.0, 1.0]] {
        assert!(closed.step(&traction, &[], &gate).is_err());
        assert_eq!(closed.state(), &accepted_field);
        assert_eq!(closed.load().state(), &accepted_bore);
        assert_eq!(closed.load().port_response().unwrap(), accepted_response);
        assert_eq!(closed.accepted_steps(), 64);
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(closed.step(&[0.0, 0.0], &[], &cancelled).is_err());
    assert_eq!(closed.state(), &accepted_field);
    assert_eq!(closed.load().state(), &accepted_bore);
    assert_eq!(closed.last_frame().unwrap().step, 64);
    closed.step(&[0.0, 0.0], &[], &gate).unwrap();

    // Force a late physical-energy refusal after private FE and bore proposals.
    let mut bounded = make(DgBoundary::Closed, ThermoviscousDgLimits {
        maximum_total_energy_j: 1e-100, ..ThermoviscousDgLimits::default()
    });
    let initial_field = bounded.state().clone();
    let initial_bore = bounded.load().state().clone();
    assert!(bounded.step(&[1.0, 0.0], &[], &gate).is_err());
    assert_eq!(bounded.state(), &initial_field);
    assert_eq!(bounded.load().state(), &initial_bore);
    assert_eq!(bounded.accepted_steps(), 0);
    assert!(bounded.last_frame().is_none());
    bounded.step(&[0.0, 0.0], &[], &gate).unwrap();
    assert_eq!(bounded.accepted_steps(), 1);
}

#[test]
fn resolved_bore_binding_requires_original_source_clock_and_owned_port() {
    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let source = ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, true), gas, BAND).unwrap();
    let other = ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, true), gas, BAND).unwrap();
    let gate = CancelGate::new_clock_free();
    let dt = 1.0 / 96_000.0;
    let prepare = || source.prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false).unwrap();
    assert!(ThermoviscousDg::new(&other, prepare(), time_dg_load(&gas, dt, DgBoundary::Closed),
        1, ThermoviscousDgLimits::default(), &gate).is_err());
    assert!(ThermoviscousDg::new(&source, prepare(), time_dg_load(&gas, 0.5 * dt, DgBoundary::Closed),
        1, ThermoviscousDgLimits::default(), &gate).is_err());
    assert!(ThermoviscousDg::new(&source, prepare(), time_dg_load(&gas, dt, DgBoundary::Closed),
        2, ThermoviscousDgLimits::default(), &gate).is_err());
    let different = GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, gas.pressure).unwrap();
    assert!(ThermoviscousDg::new(&source, prepare(), time_dg_load(&different, dt, DgBoundary::Closed),
        1, ThermoviscousDgLimits::default(), &gate).is_err());
    // The same FE pointer and gas are insufficient when bare preparation
    // drops the source's original molecular caloric fields.
    let molecular = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, true), super::source(gas, true, 1.0), BAND,
    ).unwrap();
    let bare = molecular.inner().prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false).unwrap();
    assert!(ThermoviscousDg::new(&molecular, bare, time_dg_load(&gas, dt, DgBoundary::Closed),
        1, ThermoviscousDgLimits::default(), &gate).is_err());
}

#[test]
fn resolved_pressure_playback_retains_full_states_across_callback_partitions() {
    use super::super::dg::performance::{
        ThermoviscousDgControl, ThermoviscousDgPerformance, ThermoviscousDgPerformanceConfig,
    };
    use crate::bernoulli_aperture::performance::ApertureObservation;
    use crate::pcm_wav::observation::PressureRenderer;
    use std::sync::Arc;

    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let source = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, true), source(gas, true, 1.0), BAND,
    ).unwrap();
    let dt = 1.0 / 96_000.0;
    let gate = Arc::new(CancelGate::new_clock_free());
    let controls: Vec<_> = (0..24).map(|i| ThermoviscousDgControl {
        traction_pa: vec![if i < 8 { 0.1 } else { 0.0 }, 0.0],
        wall_forces_n_per_sqrt_kg: vec![],
    }).collect();
    let make = |maximum_block| {
        let prepared = source.prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false).unwrap();
        let mut initial = prepared.zero_state();
        initial.molecular_temperature_k[0][0] = 0.005;
        let system = ThermoviscousDg::new(
            &source, prepared, time_dg_load(&gas, dt, DgBoundary::Closed), 1,
            ThermoviscousDgLimits::default(), &gate,
        ).unwrap().with_initial_state(initial).unwrap();
        ThermoviscousDgPerformance::new(system, controls.clone(), ApertureObservation::Inlet,
            ThermoviscousDgPerformanceConfig {
                sample_rate_hz: 96_000, maximum_block, maximum_control_values: 48,
            }, Arc::clone(&gate)).unwrap()
    };
    let mut whole = make(24);
    let mut divided = make(5);
    let mut a = [0.0; 24];
    let mut b = [0.0; 24];
    whole.block(&mut a).unwrap();
    for chunk in b.chunks_mut(5) { divided.block(chunk).unwrap(); }
    assert_eq!(a, b);
    assert!(a.iter().any(|p| *p != 0.0));
    assert_eq!(whole.system().state(), divided.system().state());
    assert_eq!(whole.system().load().state(), divided.system().load().state());
    assert_eq!(whole.controls(), divided.controls());
    assert_eq!(whole.samples_rendered(), 24);
    assert_eq!(whole.remaining_samples(), 0);
    assert_eq!(a[23].to_bits(), whole.last_frame().unwrap().mean_pressure_pa.to_bits());
    assert!(whole.validate_sample_count(1).is_err());
}

#[test]
fn resolved_pressure_admits_complete_controls_and_preserves_cancelled_prefix() {
    use super::super::dg::performance::{
        ThermoviscousDgControl, ThermoviscousDgPerformance, ThermoviscousDgPerformanceConfig,
    };
    use crate::bernoulli_aperture::performance::ApertureObservation;
    use crate::pcm_wav::observation::PressureRenderer;
    use crate::render::RenderError;
    use std::sync::Arc;

    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let source = ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, true), gas, BAND).unwrap();
    let dt = 1.0 / 96_000.0;
    let gate = Arc::new(CancelGate::new_clock_free());
    let make = |maximum_steps| ThermoviscousDg::new(
        &source,
        source.prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false).unwrap(),
        time_dg_load(&gas, dt, DgBoundary::Closed), 1,
        ThermoviscousDgLimits { maximum_steps, ..ThermoviscousDgLimits::default() }, &gate,
    ).unwrap();
    let config = ThermoviscousDgPerformanceConfig {
        sample_rate_hz: 96_000, maximum_block: 4, maximum_control_values: 16,
    };
    let valid = ThermoviscousDgControl {
        traction_pa: vec![1.0, 0.0], wall_forces_n_per_sqrt_kg: vec![],
    };
    for bad in [
        ThermoviscousDgControl { traction_pa: vec![1.0], ..valid.clone() },
        ThermoviscousDgControl { traction_pa: vec![0.0, 1.0], ..valid.clone() },
        ThermoviscousDgControl { traction_pa: vec![f64::NAN, 0.0], ..valid.clone() },
        ThermoviscousDgControl { wall_forces_n_per_sqrt_kg: vec![1.0], ..valid.clone() },
    ] {
        assert!(ThermoviscousDgPerformance::new(make(8), vec![bad], ApertureObservation::Inlet,
            config, Arc::clone(&gate)).is_err());
    }
    assert!(ThermoviscousDgPerformance::new(make(1), vec![valid.clone(); 2], ApertureObservation::Inlet,
        config, Arc::clone(&gate)).is_err());
    for bad in [
        ThermoviscousDgPerformanceConfig { sample_rate_hz: 48_000, ..config },
        ThermoviscousDgPerformanceConfig { maximum_control_values: 1, ..config },
    ] {
        assert!(ThermoviscousDgPerformance::new(make(8), vec![valid.clone()], ApertureObservation::Inlet,
            bad, Arc::clone(&gate)).is_err());
    }
    let mut stream = ThermoviscousDgPerformance::new(make(8), vec![valid; 8], ApertureObservation::Inlet,
        config, Arc::clone(&gate)).unwrap();
    let mut prefix = [0.0; 2];
    stream.block(&mut prefix).unwrap();
    let accepted_field = stream.system().state().clone();
    let accepted_bore = stream.system().load().state().clone();
    gate.request();
    let mut tail = [f64::MAX; 2];
    assert!(stream.block(&mut tail).is_err());
    assert_eq!(tail, [f64::MAX; 2]);
    assert_eq!(stream.system().state(), &accepted_field);
    assert_eq!(stream.system().load().state(), &accepted_bore);
    assert_eq!(stream.samples_rendered(), 2);
    assert_eq!(stream.system().accepted_steps(), 2);
    assert_eq!(stream.remaining_samples(), 6);
    assert!(matches!(stream.block(&mut tail), Err(RenderError::Poisoned)));
    assert_eq!(tail, [f64::MAX; 2]);
}
