use super::*;
use fs_couple::bernoulli_aperture::dynamic::{ApertureDrive, ApertureState, DynamicAperture, DynamicApertureSpec};
use fs_couple::bernoulli_aperture::BernoulliAperture;
use fs_couple::bernoulli_aperture::dynamic::force::PlateForceFootprint;
use fs_couple::bernoulli_aperture::tube::{ApertureTube, TubeDrive, UniformTubeSpec};
use fs_couple::bernoulli_aperture::network::{ApertureNetwork, NetworkNode, TubeNetworkSpec, TubeSection};
use fs_dcontact::Obstacle;

fn free_lay() -> Obstacle {
    Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 0.0, 2.0,
        "explicit zero-contact analytical force test".into()).unwrap()
}
#[test]
fn mechanical_force_matches_independent_midpoint_momentum_and_work() {
    let spec = DynamicApertureSpec {
        aperture: BernoulliAperture { rest_opening_m: 0.0002, width_m: 0.01, closing_pressure_pa: 1000.0 },
        mass_kg: 1e-5, stiffness_n_m: 500.0, damping_ratio: 0.02,
        density_kg_m3: 1.2, impedance_pa_s_m3: 1e6, time_step_s: 1e-5, max_steps: 10,
    };
    let old = ApertureState { opening_m: 0.00021, opening_velocity_m_s: 0.003 };
    let h = spec.time_step_s * 0.5;
    let damping = 2.0 * spec.damping_ratio * (spec.stiffness_n_m * spec.mass_kg).sqrt();
    let area = spec.stiffness_n_m * spec.aperture.rest_opening_m / spec.aperture.closing_pressure_pa;
    for force in [-0.002, 0.0, 0.002] {
        // Cancel the independently predicted swept flow. This leaves zero dp,
        // making the exact linear midpoint mechanical solution an independent oracle.
        let vm = (spec.mass_kg*old.opening_velocity_m_s + h*(force-spec.stiffness_n_m*(old.opening_m-spec.aperture.rest_opening_m)))
            / (spec.mass_kg + h*damping + h*h*spec.stiffness_n_m);
        let mut valve = DynamicAperture::new(spec, old, free_lay()).unwrap();
        let f = valve.step_with_force(ApertureDrive { body_flow_m3_s: area*vm, ..ApertureDrive::default() }, force).unwrap();
        near(f.state.opening_m, old.opening_m+spec.time_step_s*vm, 1e-11);
        near(f.state.opening_velocity_m_s, 2.0*vm-old.opening_velocity_m_s, 1e-9);
        if force != 0.0 { near(f.mechanical_work_j, force*vm*spec.time_step_s, 1e-9); }
        assert!(f.bore_pressure_pa.abs() < 1e-9);
        let scale=f.stored_energy_j.abs()+f.storage_change_j.abs()+f.mechanical_work_j.abs()+f.dissipated_energy_j;
        assert!(f.balance_residual_j().abs() < 1e-9*scale);
    }
}

#[test]
fn point_and_patch_ports_use_actual_mesh_shapes_and_the_same_virtual_work() {
    let plate = reduction(4e9, 900.0);
    let port = plate.force_port(PlateForceFootprint::Patch((0..plate.chart().mesh.tris.len()).collect())).unwrap();
    near(port.area_m2().unwrap(), 0.025*0.01, 1e-12);
    near(port.coefficient()*port.area_m2().unwrap(), plate.pressure_area_m2(), 1e-12);
    for (force, velocity) in [(-0.01, 0.1), (0.02, -0.3)] {
        near(port.generalized_force_n(force).unwrap()*velocity, force*port.velocity_m_s(velocity).unwrap(), 1e-12);
    }
    let root = plate.force_port(PlateForceFootprint::Node(0)).unwrap();
    assert_eq!(root.generalized_force_n(0.02).unwrap(), 0.0);
    let tip = plate.options().slit_edges[0][0];
    let point = plate.force_port(PlateForceFootprint::Node(tip)).unwrap();
    near(point.coefficient(), plate.shape_per_opening()[tip][0], 1e-12);
    for bad in [PlateForceFootprint::Node(usize::MAX), PlateForceFootprint::Patch(vec![]),
        PlateForceFootprint::Patch(vec![0,0]), PlateForceFootprint::Patch(vec![usize::MAX])] {
        assert!(plate.force_port(bad).is_err());
    }
    assert!(root.generalized_force_n(f64::INFINITY).is_err());
    assert!(port.velocity_m_s(f64::NAN).is_err());
}

fn coupled_plate() -> (PlateApertureReduction, UniformTubeSpec) {
    (reduction(4e9, 900.0), UniformTubeSpec { length_m: 0.25, radius_m: 0.007,
        sound_speed_m_s: 343.0, terminal_reflection: -0.8, max_length_error_m: 0.002, max_wave_memory_bytes: 1<<20 })
}
fn valve(plate: PlateApertureReduction, tube: UniformTubeSpec) -> DynamicAperture {
    let state=ApertureState { opening_m:plate.options().rest_opening_m, opening_velocity_m_s:0.0 };
    let lay=Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],1e8,2.0,
        "explicit mechanical-force contact fixture".into()).unwrap().with_internal_loss(5.0).unwrap();
    DynamicAperture::from_plate(plate,1.2,tube.characteristic_impedance(1.2).unwrap(),1.0/96000.0,2048,state,lay).unwrap()
}

#[test]
fn mechanical_press_and_release_drive_flow_contact_and_reciprocal_tube_work() {
    let (plate,spec)=coupled_plate();
    let force=-4.0*plate.stiffness_n_m()*plate.options().rest_opening_m;
    let mut model=ApertureTube::new(valve(plate,spec),spec).unwrap();
    let (mut contact,mut pressure,mut work)=(false,0.0_f64,0.0);
    for n in 0..2048 {
        let f=model.step_with_force(TubeDrive::default(),if n<1024 {force} else {0.0}).unwrap();
        work+=f.aperture.mechanical_work_j;
        contact |= f.aperture.state.opening_m < 0.0;
        pressure=pressure.max(f.aperture.bore_pressure_pa.abs());
        let scale=f.stored_energy_j+f.storage_change_j.abs()+f.dissipated_energy_j+f.aperture.mechanical_work_j.abs();
        assert!(f.balance_residual_j().abs() < 1e-8*scale.max(1e-30));
        assert_eq!(f.upstream_work_j,0.0);
        assert_eq!(f.body_work_j,0.0);
    }
    assert!(contact && pressure>1e-4 && work>0.0);
    assert!(model.stored_energy_j()>0.0,"force release must retain vibration and propagation");
}

#[test]
fn forced_network_matches_tube_and_rejected_force_preserves_wave_and_material_state() {
    let (plate,spec)=coupled_plate();
    let mut tube=ApertureTube::new(valve(plate.clone(),spec),spec).unwrap();
    let graph=TubeNetworkSpec {nodes:vec![NetworkNode::Inlet, NetworkNode::Termination {reflection:spec.terminal_reflection}],
        sections:vec![TubeSection {nodes:[0,1],length_m:spec.length_m,radius_m:spec.radius_m,max_length_error_m:spec.max_length_error_m}],
        sound_speed_m_s:spec.sound_speed_m_s,max_wave_memory_bytes:spec.max_wave_memory_bytes};
    let mut network=ApertureNetwork::new(valve(plate,spec),graph).unwrap();
    for n in 0..512 {
        if n==200 {
            let before=(network.aperture().state(),network.stored_energy_j().to_bits());
            assert!(network.step_with_force(TubeDrive::default(),f64::NAN).is_err());
            assert_eq!((network.aperture().state(),network.stored_energy_j().to_bits()),before);
            assert_eq!(network.aperture().accepted_steps(),200);
        }
        let force=if n<256 {-0.001} else {0.0};
        let a=tube.step_with_force(TubeDrive::default(),force).unwrap();
        let b=network.step_with_force(TubeDrive::default(),force).unwrap();
        near(a.aperture.state.opening_m,b.aperture.state.opening_m,1e-9);
        assert!((a.aperture.bore_pressure_pa-b.aperture.bore_pressure_pa).abs()<1e-8);
        let scale=b.stored_energy_j+b.storage_change_j.abs()+b.dissipated_energy_j+b.aperture.mechanical_work_j.abs();
        assert!(b.balance_residual_j().abs()<1e-8*scale.max(1e-30));
    }
}

#[test]
fn zero_force_preserves_original_sample_arithmetic() {
    let (plate,spec)=coupled_plate();
    let mut old=ApertureTube::new(valve(plate.clone(),spec),spec).unwrap();
    let mut forced=ApertureTube::new(valve(plate,spec),spec).unwrap();
    for _ in 0..256 {
        let drive=TubeDrive {upstream_pressure_pa:5.0,body_flow_m3_s:0.0};
        let a=old.step(drive).unwrap();let b=forced.step_with_force(drive,0.0).unwrap();
        assert_eq!(a,b);
        assert_eq!(b.aperture.mechanical_work_j.to_bits(),0.0_f64.to_bits());
    }
}
