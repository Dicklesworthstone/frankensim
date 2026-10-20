//! The original plate/contact/Maxwell state at a driven chimney, with one clock.
use super::*;
use fs_couple::bernoulli_aperture::dynamic::relaxation::{
    InitialApertureMemory, PlateRelaxationRegion, PlateRelaxationSpec,
};
use fs_couple::bernoulli_aperture::dynamic::{ApertureState, DynamicAperture};
use fs_couple::bernoulli_aperture::network::boundaries::{
    BoundaryAperture, BoundaryApertures, BoundaryDrive,
};
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::paired::BladePropagator;
use fs_dcontact::Obstacle;
use fs_material::gas::{GasSpec, GasState};
use fs_material::visco::GeneralizedMaxwell;
use fs_vfit::mimo::MimoModel;
use fs_vfit::waveguide::network::{DrivenWave, WaveguideNetwork, multiport::MimoTerminalSpec};

const DT: f64 = 1e-5;
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn graph() -> TubeNetworkSpec {
    let gas = gas();
    TubeNetworkSpec {
        nodes: vec![
            NetworkNode::Inlet,
            NetworkNode::Junction,
            NetworkNode::Termination { reflection: -0.8 },
            NetworkNode::Driven,
        ],
        sections: [(0, 1, 8, 0.007), (1, 2, 12, 0.007), (1, 3, 5, 0.003)]
            .into_iter()
            .map(|(a, b, delay, radius)| TubeSection {
                nodes: [a, b],
                length_m: delay as f64 * (gas.sound_speed * DT),
                radius_m: radius,
                max_length_error_m: 1e-15,
            })
            .collect(),
        sound_speed_m_s: gas.sound_speed,
        max_wave_memory_bytes: 1 << 20,
    }
}
fn aperture(z: f64) -> DynamicAperture {
    let (chart, mut options) = fixture(4e9, 900.0);
    options.damping_ratio = 0.0;
    let plate =
        PlateApertureReduction::from_chart(chart, options, &CancelGate::new_clock_free()).unwrap();
    let relaxation = PlateRelaxationSpec {
        regions: vec![PlateRelaxationRegion {
            triangles: (0..plate.chart().mesh.tris.len()).collect(),
            material: GeneralizedMaxwell::new(4e9, vec![(2e9, 0.001), (1e9, 0.01)]).unwrap(),
            poisson_ratio: 0.3,
            band_hz: (0.0, 1000.0),
            provenance: "manufactured boundary material".into(),
        }],
        max_branches: 2,
        max_dt_over_tau: 0.1,
        max_angular_step: 0.1,
    };
    let initial = ApertureState {
        opening_m: plate.options().rest_opening_m,
        opening_velocity_m_s: 0.0,
    };
    let contact = Obstacle::new(
        vec![-1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        1e10,
        2.0,
        "manufactured boundary lay".into(),
    )
    .unwrap()
    .with_internal_loss(0.5)
    .unwrap();
    DynamicAperture::from_plate_with_propagator(
        plate,
        gas().density,
        z,
        DT,
        8192,
        initial,
        contact,
        BladePropagator::ExactZoh,
    )
    .unwrap()
    .with_plate_relaxation(relaxation, InitialApertureMemory::Relaxed)
    .unwrap()
}
fn model() -> (WaveguideNetwork, BoundaryApertures) {
    let graph = graph();
    let gases = vec![gas(); 3];
    let (network, sections) = graph
        .realize_waveguide(&gases, DT, &CancelGate::new_clock_free())
        .unwrap();
    let network = network
        .with_mimo_terminals(MimoTerminalSpec {
            nodes: vec![2],
            model: MimoModel {
                ports: 1,
                poles: vec![],
                residues: vec![],
                d: vec![2e5],
                e: vec![0.0],
                reciprocal: true,
            },
            band_rad_s: [1.0, core::f64::consts::PI / DT],
            sampled_passivity_points: 65,
            max_states: 4,
            max_memory_bytes: 1 << 20,
        })
        .unwrap();
    let boundaries = BoundaryApertures::new(
        &graph,
        &gases,
        &network,
        vec![BoundaryAperture {
            node: 3,
            aperture: aperture(sections[2].impedance_pa_s_m3),
        }],
        1,
    )
    .unwrap();
    (network, boundaries)
}
fn drive(force: f64) -> [BoundaryDrive; 1] {
    [BoundaryDrive {
        node: 3,
        exterior_pressure_pa: 0.0,
        opening_force_n: force,
    }]
}

#[test]
fn force_closure_retains_chimney_pressure_contact_and_material_history_with_radiation() {
    let (mut network, mut boundaries) = model();
    let (mut closed, mut reopened, mut radiation_work, mut mechanical_work) =
        (false, false, 0.0, 0.0);
    let mut trapped_pressure: f64 = 0.0;
    for step in 0..8192 {
        let closing = (512..4096).contains(&step);
        let before = network.stored_energy_j() + boundaries.stored_energy_j();
        let frame = boundaries
            .step(
                &mut network,
                3.0 * (step as f64 * 0.017).sin(),
                &drive(if closing { -0.03 } else { 0.0 }),
            )
            .unwrap();
        let aperture = boundaries.aperture_frame(3).unwrap();
        assert_eq!(aperture.step, step + 1);
        assert_eq!(network.accepted_steps(), aperture.step);
        let scale = before
            + frame.stored_energy_j
            + frame.dissipated_energy_j
            + frame.network.inlet_work_j.abs()
            + frame.network.coupled_terminal_work_j.abs()
            + frame.exterior_work_j.abs()
            + frame.mechanical_work_j.abs();
        assert!(frame.balance_residual_j().abs() <= 3e-8 * scale.max(1e-20));
        assert!((frame.storage_change_j - (frame.stored_energy_j - before)).abs() < 1e-16);
        assert_eq!(
            network.node_frame(3).unwrap().pressure_pa,
            aperture.bore_pressure_pa
        );
        radiation_work += frame.network.coupled_terminal_work_j;
        mechanical_work += frame.mechanical_work_j.abs();
        if closing && aperture.midpoint_opening_m == 0.0 {
            closed = true;
            assert_eq!(aperture.jet_flow_m3_s, 0.0);
            trapped_pressure = trapped_pressure.max(aperture.bore_pressure_pa.abs());
        }
        if step > 7000 && aperture.midpoint_opening_m > 1e-4 {
            reopened = true;
        }
    }
    assert!(
        closed && reopened,
        "force must physically close and release the original aperture"
    );
    assert!(
        trapped_pressure > 0.01,
        "closed chimney must retain its real acoustic pressure"
    );
    assert!(radiation_work > 0.0 && mechanical_work > 0.0);
    let material = boundaries.aperture(3).unwrap().relaxation().unwrap();
    assert!(material.dissipated_energy_j() > 0.0);
    assert_eq!(material.memory_sqrt_j().len(), 2);
}

#[test]
fn invalid_drives_and_an_independent_graph_advance_do_not_publish_mechanical_history() {
    let (mut network, mut boundaries) = model();
    for _ in 0..32 {
        boundaries.step(&mut network, 2.0, &drive(-0.01)).unwrap();
    }
    let old_frame = *boundaries.last_frame().unwrap();
    let old_wave = network.driven_boundary(3).unwrap();
    let old_state = boundaries.aperture(3).unwrap().state();
    let old_memory = boundaries
        .aperture(3)
        .unwrap()
        .relaxation()
        .unwrap()
        .memory_sqrt_j()
        .to_vec();
    for drives in [
        vec![],
        drive(f64::NAN).to_vec(),
        vec![BoundaryDrive {
            node: 2,
            exterior_pressure_pa: 0.0,
            opening_force_n: 0.0,
        }],
    ] {
        assert!(boundaries.step(&mut network, 2.0, &drives).is_err());
        assert_eq!(boundaries.last_frame(), Some(&old_frame));
        assert_eq!(network.driven_boundary(3).unwrap(), old_wave);
        assert_eq!(boundaries.aperture(3).unwrap().state(), old_state);
        assert_eq!(
            boundaries
                .aperture(3)
                .unwrap()
                .relaxation()
                .unwrap()
                .memory_sqrt_j(),
            old_memory
        );
    }
    network
        .step_with_driven(
            0.0,
            &[DrivenWave {
                node: 3,
                outgoing_pressure_pa: 0.0,
            }],
        )
        .unwrap();
    assert!(boundaries.step(&mut network, 2.0, &drive(0.0)).is_err());
    assert_eq!(boundaries.aperture(3).unwrap().accepted_steps(), 32);
    assert_eq!(boundaries.aperture(3).unwrap().state(), old_state);
    assert_eq!(
        boundaries
            .aperture(3)
            .unwrap()
            .relaxation()
            .unwrap()
            .memory_sqrt_j(),
        old_memory
    );
}

#[test]
fn declared_local_geometry_cannot_be_substituted_under_a_matching_aperture_impedance() {
    let mut graph = graph();
    let gases = vec![gas(); 3];
    let (network, represented) = graph
        .realize_waveguide(&gases, DT, &CancelGate::new_clock_free())
        .unwrap();
    graph.sections[2].radius_m *= 2.0;
    assert!(
        BoundaryApertures::new(
            &graph,
            &gases,
            &network,
            vec![BoundaryAperture {
                node: 3,
                aperture: aperture(represented[2].impedance_pa_s_m3)
            }],
            1
        )
        .is_err()
    );
}
