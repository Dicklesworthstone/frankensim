//! The retained material-modal reed and a material-card key share one performer.
use super::*;
use fs_couple::bernoulli_aperture::{
    dynamic::{ApertureState, DynamicAperture},
    network::{NetworkNode, TubeSection, boundaries::{BoundaryAperture, BoundaryDrive}},
    paired::BladePropagator,
    performance::{AperturePerformanceConfig, boundaries::{BoundaryControl, BoundaryProgram}},
    plate::{PlateApertureOptions, PlateApertureReduction},
    two_sided::{LoadedAperture, TwoSidedAperture, TwoSidedLimits, ValveSide,
        combined_inlet_state, observer::{TwoSidedPressure, schedule::ScheduledTwoSidedPressure}},
    viscous::ViscousSlit,
};
use fs_couple::pcm_wav::{baffled::CircularOutletReceiver, observation::PressureRenderer};
use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureTrack, GestureValue};

#[test]
fn material_card_key_change_reaches_original_modal_reed_and_pressure_without_source_reduction_loss() {
    let render = |density, block| {
        let dt = 1.0 / 48_000.0;
        let gate = CancelGate::new_clock_free();
        let mut ducts = sides();
        ducts[1].graph.nodes[1] = NetworkNode::Junction;
        for (node, kind, delay, radius) in [(2, NetworkNode::Termination { reflection: -0.8 }, 12.0, 0.006),
            (3, NetworkNode::Driven, 5.0, 0.003)] {
            ducts[1].graph.nodes.push(kind);
            ducts[1].graph.sections.push(TubeSection { nodes: [1, node],
                length_m: delay * gas().sound_speed * dt, radius_m: radius, max_length_error_m: 1e-12 });
            ducts[1].gases.push(gas());
        }
        let (gas, z) = combined_inlet_state(&ducts).unwrap();
        let mut plates = vec![input(1.0, 900.0), input(0.0, 1200.0)];
        for plate in &mut plates { plate.1.maximum_angular_step = 2.0; }
        let identities: Vec<_> = plates.iter().map(|p| p.0.specimen_identity()).collect();
        let valve = ModalAperture::from_viscous_material_plates_with_propagator(
            plates, 4e-4, gas, z, dt, 256, limits(), vec![BladeState::default(); 6],
            contact(), BladePropagator::ExactZoh, &gate,
        ).unwrap();
        let original_spec = valve.spec().clone();
        let original_viscosity = valve.viscosity_matrix_n_s_m().unwrap().to_vec();
        let key_source = input(0.0, density);
        let key_identity = key_source.0.specimen_identity();
        let plate = PlateApertureReduction::from_material_chart(key_source.0, PlateApertureOptions {
            assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::Clamped },
            eigenvalue_window: (1.0, (std::f64::consts::TAU * 450.0).powi(2)),
            mode_index: 0, eigensolver: fs_plate::SliceOptions::default(),
            max_nodes: 100, max_triangles: 100, slit_edges: key_source.1.slit_edges,
            rest_opening_m: 0.0002, damping_ratio: 0.01,
            max_slit_mode_variation: 0.25, max_slope: 0.1,
        }, &gate).unwrap();
        let key_z = ducts[1].graph.sections[2].realize(gas.sound_speed, gas.density, dt)
            .unwrap().impedance_pa_s_m3;
        let mut key = DynamicAperture::from_plate_with_propagator(plate, gas.density, key_z,
            dt, 256, ApertureState { opening_m: 0.0002, opening_velocity_m_s: 0.0 },
            contact(), BladePropagator::ExactZoh).unwrap();
        key.set_viscous_slit(ViscousSlit::from_gas(&gas, 0.0003).unwrap()).unwrap();
        let system = TwoSidedAperture::new(valve.into(), ducts, TwoSidedLimits {
            maximum_energy_j: 1.0, energy_absolute_tolerance_j: 1e-16,
            energy_relative_tolerance: 1e-7, flow_absolute_tolerance_m3_s: 1e-16,
            flow_relative_tolerance: 1e-8,
        }, &gate).unwrap().with_boundary_apertures(ValveSide::Downstream,
            vec![BoundaryAperture { node: 3, aperture: key }], 1, &gate).unwrap();
        let source = TwoSidedPressure::new(system, ValveSide::Downstream, 2,
            CircularOutletReceiver { position_m: [0.0, 0.0, 0.25], radial_rings: 4,
                angular_points: 16, maximum_frequency_hz: 1000.0 }, 48_000, 64).unwrap();
        let phrase = GestureSchedule::try_new(1, vec![GestureTrack { id: "mouth".into(),
            target: GestureTarget::BlowingPressure, initial: GestureValue::PressurePa(5.0),
            events: vec![] }]).unwrap();
        let initial = BoundaryDrive { node: 3, exterior_pressure_pa: 0.0, opening_force_n: 0.0 };
        let mut performer = ScheduledTwoSidedPressure::new(source, phrase, AperturePerformanceConfig {
            sample_rate_hz: 48_000, samples: 256, max_block: 64, max_compile_work: 10000, max_controls: 32,
        }).unwrap().with_boundary_controls(ValveSide::Downstream, BoundaryProgram {
            initial: vec![initial], controls: vec![
                BoundaryControl { sample: 17, drive: BoundaryDrive { opening_force_n: -0.005, ..initial } },
                BoundaryControl { sample: 97, drive: initial },
            ],
        }).unwrap();
        let mut pressure = vec![0.0; 256];
        for chunk in pressure.chunks_mut(block) { performer.block(chunk).unwrap(); }
        let system = performer.source().system();
        let LoadedAperture::Modal(valve) = system.valve() else { panic!("original material-modal owner") };
        assert_eq!(valve.spec(), &original_spec);
        assert_eq!(valve.viscosity_matrix_n_s_m().unwrap(), original_viscosity);
        for (plate, identity) in valve.material_plates().iter().zip(identities) {
            assert_eq!(plate.source().specimen_identity(), identity);
            assert_eq!(plate.viscosity().unwrap().chart_identity(), identity);
        }
        let key = system.boundary_apertures(ValveSide::Downstream).unwrap().aperture(3).unwrap();
        assert_eq!(key.plate_reduction().unwrap().material_chart().unwrap().specimen_identity(), key_identity);
        assert_eq!(key.accepted_steps(), 256);
        assert!(performer.source().last_frame().unwrap().balance_residual_j().abs() < 1e-12);
        (pressure, performer.source().last_frame().unwrap().aperture.state)
    };
    let base = render(900.0, 1);
    assert_eq!(base, render(900.0, 37));
    let changed = render(1800.0, 37);
    assert_ne!(base.0, changed.0);
    assert_ne!(base.1, changed.1);
}
