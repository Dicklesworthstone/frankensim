use super::*;
use super::super::dg::{ThermoviscousDg, ThermoviscousDgLimits};
use crate::dg_port::DgPortLoad;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, PortId, MAX_DG_PROPAGATOR_DOFS};
use fs_exec::CancelGate;

pub(super) fn time_dg_load(gas: &GasState, duration_s: f64, outlet: DgBoundary) -> DgPortLoad {
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

// Two opposite openings carry finite ducts; a separate side opening supplies
// authored work without overwriting either of the solved interface pressures.
fn three_port_field<'c>(
    mesh: &'c TetComplex,
    xyz: &[[f64; 3]],
    gas: &GasState,
) -> ThermoviscousFe<'c> {
    let geometry = HelmholtzFe::new(mesh, xyz, 2, Medium {
        density_kg_m3: gas.density, sound_speed_m_s: gas.sound_speed,
    }).unwrap();
    let faces = geometry.boundary_faces();
    let opening = |axis: usize, coordinate: f64| faces.iter().copied().filter(|face| {
        face.iter().all(|&v| (xyz[v as usize][axis] - coordinate).abs() < 1e-14)
    }).collect::<Vec<_>>();
    let left = opening(0, 0.0);
    let right = opening(0, LENGTH);
    let drive = opening(1, HEIGHT);
    let walls: Vec<_> = faces.into_iter().filter(|face| {
        !left.contains(face) && !right.contains(face) && !drive.contains(face)
    }).collect();
    ThermoviscousFe::new(mesh, xyz, source_medium(gas, 0.0).unwrap(),
        &ThermoviscousBoundary {
            no_slip_faces: walls.clone(), isothermal_faces: walls,
            ports: vec![BoundaryMode::uniform(&left), BoundaryMode::uniform(&right),
                BoundaryMode::uniform(&drive)],
        }, ThermoviscousLimits::default(), &mut || false,
    ).unwrap()
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

pub(super) fn wall_transport(rho: f64, cp: f64, k: f64) -> PreparedThermalTransport {
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
        .with_additional_duct(0, time_dg_load(&gas, 1e-7, DgBoundary::Closed), &gate)
        .unwrap()
        .with_initial_state(initial)
        .unwrap();
        let frame = coupled.step(&[0.0, 0.0], &[], &gate).unwrap().clone();
        assert!(frame.field.work.solid_conductive_heat_j > 0.0);
        assert!(frame.field.state.temperature_k.iter().any(|t| *t != 0.0));
        assert!(frame.field.state.molecular_temperature_k[0].iter().any(|t| *t != 0.0));
        assert!(coupled.load().state().pressure().iter().any(|p| *p != 0.0));
        assert!(coupled.load_at_port(0).unwrap().state().pressure().iter().any(|p| *p != 0.0));
        assert_eq!(frame.additional_ducts.len(), 1);
        assert_eq!(frame.additional_ducts[0].port, 0);
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
fn multiple_bores_use_the_full_field_response_and_replay_every_original_owner() {
    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let source = ThermoviscousInterior::from_relaxation(
        three_port_field(&mesh, &xyz, &gas), source(gas, true, 1.0), BAND,
    ).unwrap();
    let dt = 1.0 / 96_000.0;
    let gate = CancelGate::new_clock_free();
    let prepared = source.prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false).unwrap();
    let mut initial = prepared.zero_state();
    initial.molecular_temperature_k[0][0] = 0.002;
    let mut coupled = ThermoviscousDg::new(
        &source, prepared, time_dg_load(&gas, dt, DgBoundary::Closed), 1,
        ThermoviscousDgLimits::default(), &gate,
    ).unwrap().with_additional_duct(0, time_dg_load(&gas, dt, DgBoundary::Matched), &gate)
        .unwrap().with_initial_state(initial).unwrap();
    assert_eq!(coupled.coupled_ports(), &[1, 0]);
    assert_eq!(coupled.additional_loads().map(|(port, _)| port).collect::<Vec<_>>(), vec![0]);
    assert!(coupled.load_at_port(2).is_none());
    let mut molecular_heat = 0.0;
    let mut work_at_each_inlet = [0.0; 2];
    for drive in [0.4, 0.4, 0.0, -0.2, 0.0] {
        let initial = coupled.state().clone();
        let mut original = [coupled.load().clone(), coupled.load_at_port(0).unwrap().clone()];
        let (free, mobility) = coupled.prepared().average_relations(
            &initial, &[0.0, 0.0, drive], &[], &[1, 0], &mut || false,
        ).unwrap();
        assert!(mobility[1] != 0.0 && mobility[2] != 0.0,
            "the fixture must retain cross-port field response");
        let relation = [original[0].port_response().unwrap(), original[1].port_response().unwrap()];
        let z = [relation[0].impedance_pa_s_m3, relation[1].impedance_pa_s_m3];
        // An explicit two-by-two solve detects an incorrectly diagonalized FE
        // port map, independently of the coupler's general multi-port factor.
        let m = [1.0 + z[0] * mobility[0], z[0] * mobility[1],
            z[1] * mobility[2], 1.0 + z[1] * mobility[3]];
        let rhs = [relation[0].history_pressure_pa - z[0] * free[0],
            relation[1].history_pressure_pa - z[1] * free[1]];
        let determinant = m[0] * m[3] - m[1] * m[2];
        let expected_pressure = [(rhs[0] * m[3] - m[1] * rhs[1]) / determinant,
            (m[0] * rhs[1] - m[2] * rhs[0]) / determinant];
        let frame = coupled.step(&[0.0, 0.0, drive], &[], &gate).unwrap().clone();
        assert_eq!(frame.additional_ducts.len(), 1);
        let additional = &frame.additional_ducts[0];
        assert_eq!(additional.port, 0);
        let pressures = [frame.mean_pressure_pa, additional.mean_pressure_pa];
        let flows = [frame.mean_flow_m3_s, additional.mean_flow_m3_s];
        let pressure_scale = pressures.iter().chain(&expected_pressure)
            .map(|p| p.abs()).fold(1e-30, f64::max);
        for (actual, expected) in pressures.iter().zip(expected_pressure) {
            assert!((actual - expected).abs() < 1e-10 * pressure_scale);
        }
        let replay = coupled.prepared().step(&initial,
            &[pressures[1], pressures[0], drive], &[], &mut || false,
        ).unwrap();
        assert_eq!(frame.field.state, replay.state);
        assert_eq!(frame.field.mean_inward_flow_m3_s, replay.mean_inward_flow_m3_s);
        let mut interface_residual = 0.0;
        let mut bore_energy = 0.0;
        for (index, port) in [1, 0].into_iter().enumerate() {
            assert_eq!(flows[index].to_bits(), (-replay.mean_inward_flow_m3_s[port]).to_bits());
            let trial = original[index].preview_step(flows[index]).unwrap();
            let work = trial.inlet_work_j();
            let residual = work + pressures[index] * (replay.mean_inward_flow_m3_s[port] * dt);
            assert!(residual.abs() <= frame.balance_tolerance_j);
            work_at_each_inlet[index] += work.abs();
            interface_residual += residual;
            let output = trial.accept();
            bore_energy += output.ledger.stored_energy_j;
            assert_eq!(coupled.load_at_port(port).unwrap().state(), original[index].state());
            assert_eq!(if index == 0 { &frame.duct } else { &additional.duct }, &output);
            if index == 1 {
                assert_eq!(additional.interface_work_residual_j.to_bits(), residual.to_bits());
            }
        }
        assert!((frame.interface_work_residual_j - interface_residual).abs() <= frame.balance_tolerance_j);
        let energy = replay.work.final_energy_j + bore_energy;
        assert!((frame.stored_energy_j - energy).abs() <= 1e-12 * energy.abs().max(1e-30));
        assert!((frame.storage_change_j + frame.material_dissipation_j + frame.terminal_work_j
            - frame.authored_external_work_j).abs() <= frame.balance_tolerance_j);
        molecular_heat += frame.field.work.molecular_relaxation_heat_j;
    }
    assert!(molecular_heat > 0.0);
    assert!(work_at_each_inlet.iter().all(|work| *work > 0.0));
    assert_eq!(coupled.accepted_steps(), 5);
}

#[test]
fn finite_upstream_bore_reacts_on_the_field_and_downstream_bore_atomically() {
    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let source = ThermoviscousInterior::from_gas(three_port_field(&mesh, &xyz, &gas), gas, BAND).unwrap();
    let dt = 1.0 / 96_000.0;
    let gate = CancelGate::new_clock_free();
    let make = |limits| ThermoviscousDg::new(
        &source,
        source.prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false).unwrap(),
        time_dg_load(&gas, dt, DgBoundary::Closed), 1, limits, &gate,
    ).unwrap();
    let mut prescribed = make(ThermoviscousDgLimits::default());
    let mut finite = make(ThermoviscousDgLimits::default())
        .with_additional_duct(0, time_dg_load(&gas, dt, DgBoundary::Closed), &gate).unwrap();
    for sample in 0..64 {
        let drive = if sample < 16 { 1.0 } else { 0.0 };
        prescribed.step(&[0.0, 0.0, drive], &[], &gate).unwrap();
        finite.step(&[0.0, 0.0, drive], &[], &gate).unwrap();
    }
    let field_reaction = finite.state().pressure_pa.iter().zip(&prescribed.state().pressure_pa)
        .map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    let bore_reaction = finite.load().state().pressure().iter().zip(prescribed.load().state().pressure())
        .map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    assert!(field_reaction > 1e-7, "finite upstream history must alter field pressure: {field_reaction:e}");
    assert!(bore_reaction > 1e-7, "finite upstream history must reach the other bore: {bore_reaction:e}");
    let accepted_field = finite.state().clone();
    let accepted_bores = [finite.load().state().clone(), finite.load_at_port(0).unwrap().state().clone()];
    let accepted_responses = [finite.load().port_response().unwrap(), finite.load_at_port(0).unwrap().port_response().unwrap()];
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    for (traction, request) in [([1.0, 0.0, 0.0], &gate), ([0.0, 1.0, 0.0], &gate),
        ([0.0, 0.0, f64::NAN], &gate), ([0.0, 0.0, 0.0], &cancelled)] {
        assert!(finite.step(&traction, &[], request).is_err());
        assert_eq!(finite.state(), &accepted_field);
        for (index, port) in [1, 0].into_iter().enumerate() {
            let load = finite.load_at_port(port).unwrap();
            assert_eq!(load.state(), &accepted_bores[index]);
            assert_eq!(load.port_response().unwrap(), accepted_responses[index]);
        }
        assert_eq!(finite.accepted_steps(), 64);
        assert_eq!(finite.last_frame().unwrap().step, 64);
    }
    finite.step(&[0.0, 0.0, 0.0], &[], &gate).unwrap();
    let mut bounded = make(ThermoviscousDgLimits {
        maximum_total_energy_j: 1e-100, ..ThermoviscousDgLimits::default()
    }).with_additional_duct(0, time_dg_load(&gas, dt, DgBoundary::Closed), &gate).unwrap();
    let initial_field = bounded.state().clone();
    let initial_bores = [bounded.load().state().clone(), bounded.load_at_port(0).unwrap().state().clone()];
    // This drive reaches all private candidates before the total-energy check.
    assert!(bounded.step(&[0.0, 0.0, 1.0], &[], &gate).is_err());
    assert_eq!(bounded.state(), &initial_field);
    for (index, port) in [1, 0].into_iter().enumerate() {
        assert_eq!(bounded.load_at_port(port).unwrap().state(), &initial_bores[index]);
    }
    assert_eq!(bounded.accepted_steps(), 0);
    assert!(bounded.last_frame().is_none());
    bounded.step(&[0.0, 0.0, 0.0], &[], &gate).unwrap();
    assert_eq!(bounded.accepted_steps(), 1);
    assert!(make(ThermoviscousDgLimits::default())
        .with_additional_duct(1, time_dg_load(&gas, dt, DgBoundary::Closed), &gate).is_err());
    assert!(make(ThermoviscousDgLimits::default())
        .with_additional_duct(0, time_dg_load(&gas, 0.5 * dt, DgBoundary::Closed), &gate).is_err());
    let other = GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, gas.pressure).unwrap();
    assert!(make(ThermoviscousDgLimits::default())
        .with_additional_duct(0, time_dg_load(&other, dt, DgBoundary::Closed), &gate).is_err());
    let mut advanced = time_dg_load(&gas, dt, DgBoundary::Closed);
    advanced.preview_step(0.0).unwrap().accept();
    assert!(make(ThermoviscousDgLimits::default()).with_additional_duct(0, advanced, &gate).is_err());
    assert!(prescribed.with_additional_duct(0, time_dg_load(&gas, dt, DgBoundary::Closed), &gate).is_err());
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
fn additional_bore_pressure_playback_uses_its_own_frame_and_retains_all_histories() {
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
    let source = ThermoviscousInterior::from_relaxation(
        three_port_field(&mesh, &xyz, &gas), source(gas, true, 1.0), BAND,
    ).unwrap();
    let dt = 1.0 / 96_000.0;
    let gate = Arc::new(CancelGate::new_clock_free());
    let controls: Vec<_> = (0..24).map(|i| ThermoviscousDgControl {
        traction_pa: vec![0.0, 0.0, if i < 8 { 0.1 } else { 0.0 }],
        wall_forces_n_per_sqrt_kg: vec![],
    }).collect();
    let system = || {
        let prepared = source.prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false).unwrap();
        let mut initial = prepared.zero_state();
        initial.molecular_temperature_k[0][0] = 0.002;
        ThermoviscousDg::new(
            &source, prepared, time_dg_load(&gas, dt, DgBoundary::Closed), 1,
            ThermoviscousDgLimits::default(), &gate,
        ).unwrap().with_additional_duct(0, time_dg_load(&gas, dt, DgBoundary::Matched), &gate)
            .unwrap().with_initial_state(initial).unwrap()
    };
    let config = ThermoviscousDgPerformanceConfig {
        sample_rate_hz: 96_000, maximum_block: 24, maximum_control_values: 72,
    };
    let mut whole = ThermoviscousDgPerformance::new_at_port(
        system(), 0, controls.clone(), ApertureObservation::Inlet, config, Arc::clone(&gate),
    ).unwrap();
    let mut divided = ThermoviscousDgPerformance::new_at_port(
        system(), 0, controls.clone(), ApertureObservation::Inlet,
        ThermoviscousDgPerformanceConfig { maximum_block: 5, ..config }, Arc::clone(&gate),
    ).unwrap();
    let mut a = [0.0; 24];
    let mut b = [0.0; 24];
    whole.block(&mut a).unwrap();
    for chunk in b.chunks_mut(5) { divided.block(chunk).unwrap(); }
    assert_eq!(a, b);
    assert!(a.iter().any(|pressure| *pressure != 0.0));
    assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
    assert_eq!(whole.observation_port(), 0);
    assert_eq!(whole.system().state(), divided.system().state());
    for port in [1, 0] {
        assert_eq!(whole.system().load_at_port(port).unwrap().state(),
            divided.system().load_at_port(port).unwrap().state());
    }
    assert_eq!(a[23].to_bits(), whole.last_frame().unwrap().additional_ducts[0].mean_pressure_pa.to_bits());
    assert_eq!(whole.samples_rendered(), 24);
    assert_eq!(whole.remaining_samples(), 0);
    assert!(ThermoviscousDgPerformance::new_at_port(system(), 2, controls.clone(),
        ApertureObservation::Inlet, config, Arc::clone(&gate)).is_err());
    for port in [1, 0] {
        let mut bad = controls[0].clone();
        bad.traction_pa[port] = 1.0;
        assert!(ThermoviscousDgPerformance::new_at_port(system(), 0, vec![bad],
            ApertureObservation::Inlet, config, Arc::clone(&gate)).is_err());
    }
    let mut stream = ThermoviscousDgPerformance::new_at_port(system(), 0, controls,
        ApertureObservation::Inlet, config, Arc::clone(&gate)).unwrap();
    let mut prefix = [0.0; 2];
    stream.block(&mut prefix).unwrap();
    assert_eq!(prefix, a[..2]);
    let accepted_field = stream.system().state().clone();
    let accepted_bores = [stream.system().load().state().clone(),
        stream.system().load_at_port(0).unwrap().state().clone()];
    gate.request();
    let mut tail = [f64::MAX; 2];
    assert!(stream.block(&mut tail).is_err());
    assert_eq!(tail, [f64::MAX; 2]);
    assert_eq!(stream.system().state(), &accepted_field);
    for (index, port) in [1, 0].into_iter().enumerate() {
        assert_eq!(stream.system().load_at_port(port).unwrap().state(), &accepted_bores[index]);
    }
    assert_eq!(stream.samples_rendered(), 2);
    assert_eq!(stream.system().accepted_steps(), 2);
    assert_eq!(stream.remaining_samples(), 22);
    assert!(matches!(stream.block(&mut tail), Err(RenderError::Poisoned)));
    assert_eq!(tail, [f64::MAX; 2]);
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
