//! Reed, material key and every chimney wave advance as one physical system.
use super::*;
use fs_couple::bernoulli_aperture::dynamic::relaxation::{
    InitialApertureMemory, PlateRelaxationRegion, PlateRelaxationSpec,
};
use fs_couple::bernoulli_aperture::network::boundaries::{BoundaryAperture, BoundaryDrive};
use fs_couple::bernoulli_aperture::paired::BladePropagator;
use fs_couple::bernoulli_aperture::plate::{PlateApertureOptions, PlateApertureReduction};
use fs_material::gas::{GasSpec, GasState};
use fs_material::visco::GeneralizedMaxwell;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};

pub(super) const RELEASE: u64 = 2048;
pub(super) const CYCLE_STEPS: u64 = 6144;

pub(super) fn geometry() -> (TubeNetworkSpec, Vec<GasState>) {
    let gases: Vec<_> = [293.15, 313.15, 283.15].into_iter().map(|temperature|
        GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101325.0).unwrap()).collect();
    let sections = [(0, 1, 8, 0.007), (1, 2, 12, 0.009), (1, 3, 5, 0.003)]
        .into_iter().zip(&gases).map(|((a, b, delay, radius), gas)| TubeSection {
            nodes: [a, b], length_m: delay as f64 * (gas.sound_speed * DT),
            radius_m: radius, max_length_error_m: 1e-15,
        }).collect();
    (TubeNetworkSpec { nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
        NetworkNode::Termination { reflection: -0.8 }, NetworkNode::Driven],
        sections, sound_speed_m_s: gases[0].sound_speed, max_wave_memory_bytes: 1 << 20 }, gases)
}

pub(super) fn key_plate() -> PlateApertureReduction {
    let mesh = PlateMesh::rectangle(0.025, 0.01, 6, 2);
    let chart = PlateChart::with_boundary_and_regions(mesh,
        PlateSection::isotropic(4e9, 0.3, 0.0003, 900.0).unwrap(),
        vec![0, 7, 14], vec![]).unwrap();
    PlateApertureReduction::from_chart(chart, PlateApertureOptions {
        assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::Clamped },
        eigenvalue_window: (1.0, (core::f64::consts::TAU * 450.0).powi(2)),
        mode_index: 0, eigensolver: SliceOptions::default(), max_nodes: 100, max_triangles: 100,
        slit_edges: vec![[6, 13], [13, 20]], rest_opening_m: 0.0002,
        damping_ratio: 0.0, max_slit_mode_variation: 0.25, max_slope: 0.1,
    }, &CancelGate::new_clock_free()).unwrap()
}

pub(super) fn key(plate: &PlateApertureReduction, gas: GasState, z: f64, steps: u64)
    -> BoundaryAperture
{
    let relaxation = PlateRelaxationSpec { regions: vec![PlateRelaxationRegion {
        triangles: (0..plate.chart().mesh.tris.len()).collect(),
        material: GeneralizedMaxwell::new(4e9, vec![(2e9, 0.001), (1e9, 0.01)]).unwrap(),
        poisson_ratio: 0.3, band_hz: (0.0, 1000.0),
        provenance: "manufactured key material; not measured reed-instrument data".into(),
    }], max_branches: 2, max_dt_over_tau: 0.1, max_angular_step: 0.1 };
    let lay = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e10, 2.0,
        "manufactured moving-key contact".into()).unwrap().with_internal_loss(0.5).unwrap();
    let aperture = DynamicAperture::from_plate_with_propagator(plate.clone(), gas.density,
        z, DT, steps, ApertureState { opening_m: plate.options().rest_opening_m,
            opening_velocity_m_s: 0.0 }, lay, BladePropagator::ExactZoh).unwrap()
        .with_plate_relaxation(relaxation, InitialApertureMemory::Relaxed).unwrap();
    BoundaryAperture { node: 3, aperture }
}

fn unkeyed(steps: u64) -> (ApertureNetwork, Vec<GasState>) {
    let (graph, gases) = geometry();
    let spec = DynamicApertureSpec {
        aperture: BernoulliAperture { rest_opening_m: 4e-4, width_m: 0.013, closing_pressure_pa: 6000.0 },
        mass_kg: 1e-5, stiffness_n_m: 500.0, damping_ratio: 0.35,
        density_kg_m3: gases[0].density, impedance_pa_s_m3: graph.inlet_impedance_with_gases(&gases).unwrap(),
        time_step_s: DT, max_steps: steps,
    };
    let reed = DynamicAperture::new_with_propagator(spec, initial(),
        Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 1e8, 2.0,
            "manufactured inlet reed contact".into()).unwrap().with_internal_loss(5.0).unwrap(),
        BladePropagator::ExactZoh).unwrap();
    (ApertureNetwork::with_section_gases(reed, graph, gases.clone()).unwrap(), gases)
}

pub(super) fn make(plate: &PlateApertureReduction, steps: u64) -> ApertureNetwork {
    let (network, gases) = unkeyed(steps);
    let z = network.represented_sections()[2].impedance_pa_s_m3;
    network.with_boundary_apertures(&gases, vec![key(plate, gases[2], z, steps)], 1).unwrap()
}

pub(super) fn drives(step: u64) -> (TubeDrive, f64, [BoundaryDrive; 1]) {
    let blowing = step < 1024;
    (TubeDrive { upstream_pressure_pa: if blowing { 100.0 } else { 0.0 },
        body_flow_m3_s: if blowing { 2e-7 * (step as f64 * 0.023).sin() } else { 0.0 } },
        if blowing { 0.001 * (step as f64 * 0.019).sin() } else { 0.0 },
        [BoundaryDrive { node: 3, exterior_pressure_pa: if blowing { 1.0 } else { 0.0 },
            opening_force_n: if (512..RELEASE).contains(&step) { -0.03 } else { 0.0 } }])
}

pub(super) fn check_balance(frame: ApertureNetworkFrame, before: f64) {
    let keys = frame.boundaries.expect("the material key participates in the accepted frame");
    assert_eq!(frame.network, keys.network);
    let storage = frame.aperture.stored_energy_j + keys.stored_energy_j;
    assert!((frame.stored_energy_j - storage).abs() <= 8.0 * f64::EPSILON * storage.max(1e-20));
    let scale = before + frame.stored_energy_j + frame.dissipated_energy_j.abs()
        + frame.upstream_work_j.abs() + frame.body_work_j.abs()
        + frame.aperture.mechanical_work_j.abs() + keys.mechanical_work_j.abs()
        + keys.exterior_work_j.abs() + frame.network.coupled_terminal_work_j.abs();
    assert!(frame.balance_residual_j().abs() <= 3e-8 * scale.max(1e-20));
    assert!((frame.storage_change_j - (frame.stored_energy_j - before)).abs()
        <= 3e-8 * scale.max(1e-20));
}

#[test]
fn single_reed_key_closes_and_releases_stored_motion_after_blowing_stops() {
    let plate = key_plate();
    let mut model = make(&plate, CYCLE_STEPS);
    let (_, gases) = geometry();
    assert_eq!(model.boundary_apertures().unwrap().gas(3), Some(&gases[2]));
    assert_ne!(gases[0], gases[2]);
    let (mut closed, mut reopened) = (false, false);
    let (mut key_work, mut reed_work, mut body_work, mut exterior_work) = (0.0, 0.0, 0.0, 0.0);
    let mut closure_pressure = 0.0_f64;
    for step in 0..CYCLE_STEPS {
        let before = model.stored_energy_j();
        if step == RELEASE {
            assert!(before > 0.0);
            let memory = model.boundary_apertures().unwrap().aperture(3).unwrap().relaxation().unwrap();
            assert!(memory.memory_sqrt_j().iter().any(|x| *x != 0.0));
        }
        let (drive, force, keys) = drives(step);
        let frame = model.step_with_boundaries(drive, force, &keys).unwrap();
        check_balance(frame, before);
        let key = model.boundary_apertures().unwrap().aperture_frame(3).unwrap();
        assert_eq!(key.step, frame.aperture.step);
        assert_eq!(model.node_frame(3).unwrap().pressure_pa, key.bore_pressure_pa);
        assert!((model.stored_energy_j() - frame.stored_energy_j).abs()
            <= 8.0 * f64::EPSILON * frame.stored_energy_j.max(1e-20));
        key_work += frame.boundaries.unwrap().mechanical_work_j.abs();
        reed_work += frame.aperture.mechanical_work_j.abs();
        body_work += frame.body_work_j.abs();
        exterior_work += frame.boundaries.unwrap().exterior_work_j.abs();
        if (512..RELEASE).contains(&step) && key.midpoint_opening_m == 0.0 {
            closed = true;
            assert_eq!(key.jet_flow_m3_s, 0.0);
            closure_pressure = closure_pressure.max(key.bore_pressure_pa.abs());
        }
        if step >= RELEASE {
            assert!(frame.stored_energy_j <= before + 3e-8 * before.max(1e-20));
            reopened |= key.midpoint_opening_m > 1e-4;
        }
    }
    assert!(closed && reopened);
    assert!(closure_pressure > 0.01, "closure preserves the real chimney pressure, not a reset load");
    assert!(key_work > 0.0 && reed_work > 0.0 && body_work > 0.0 && exterior_work > 0.0);
    assert!(model.boundary_apertures().unwrap().aperture(3).unwrap()
        .relaxation().unwrap().dissipated_energy_j() > 0.0);
}

#[test]
fn single_reed_refusal_and_budget_resume_preserve_key_memory_and_all_waves() {
    let plate = key_plate();
    // Equal density/sound speed cannot stand in for the complete transport gas.
    // A different viscosity must not silently enter a retained chimney's key.
    let (cold, gases) = unkeyed(64);
    let mut wrong_transport = gases[2];
    wrong_transport.dynamic_viscosity *= 2.0;
    wrong_transport.prandtl *= 2.0;
    wrong_transport.validate_coherence().unwrap();
    let mut wrong_key = key(&plate, gases[2], cold.represented_sections()[2].impedance_pa_s_m3, 64);
    wrong_key.aperture.set_viscous_slit(
        fs_couple::bernoulli_aperture::viscous::ViscousSlit::from_gas(&wrong_transport, 0.0003).unwrap()
    ).unwrap();
    assert!(cold.with_boundary_apertures(&gases, vec![wrong_key], 1).is_err());
    let mut actual = make(&plate, 64);
    let mut reference = make(&plate, 96);
    let (drive, force, keys) = drives(600);
    for _ in 0..32 {
        assert_eq!(actual.step_with_boundaries(drive, force, &keys).unwrap(),
            reference.step_with_boundaries(drive, force, &keys).unwrap());
    }
    let state = actual.aperture().state();
    let energy = actual.stored_energy_j();
    let nodes: Vec<_> = (0..4).map(|node| *actual.node_frame(node).unwrap()).collect();
    let key_state = actual.boundary_apertures().unwrap().aperture(3).unwrap().state();
    let memory = actual.boundary_apertures().unwrap().aperture(3).unwrap()
        .relaxation().unwrap().memory_sqrt_j().to_vec();
    for invalid in [vec![], vec![BoundaryDrive { node: 2, ..keys[0] }],
        vec![BoundaryDrive { opening_force_n: f64::NAN, ..keys[0] }]] {
        assert!(actual.step_with_boundaries(drive, force, &invalid).is_err());
    }
    assert_eq!(actual.aperture().state(), state);
    assert_eq!(actual.aperture().accepted_steps(), 32);
    assert_eq!(actual.stored_energy_j(), energy);
    let boundary = actual.boundary_apertures().unwrap().aperture(3).unwrap();
    assert_eq!(boundary.state(), key_state);
    assert_eq!(boundary.accepted_steps(), 32);
    assert_eq!(boundary.relaxation().unwrap().memory_sqrt_j(), memory);
    for (node, frame) in nodes.iter().enumerate() { assert_eq!(actual.node_frame(node), Some(frame)); }
    for _ in 32..64 {
        assert_eq!(actual.step_with_boundaries(drive, force, &keys).unwrap(),
            reference.step_with_boundaries(drive, force, &keys).unwrap());
    }
    assert!(actual.step_with_boundaries(drive, force, &keys).is_err());
    assert!(actual.extend_step_budget(32).is_err());
    assert!(actual.extend_step_budget(64).is_err());
    assert_eq!(actual.aperture().accepted_steps(), 64);
    assert_eq!(actual.boundary_apertures().unwrap().aperture(3).unwrap().accepted_steps(), 64);
    actual.extend_step_budget(96).unwrap();
    assert_eq!(actual.boundary_apertures().unwrap().aperture(3).unwrap().spec().max_steps, 96);
    for _ in 64..96 {
        assert_eq!(actual.step_with_boundaries(drive, force, &keys).unwrap(),
            reference.step_with_boundaries(drive, force, &keys).unwrap());
    }
}
