//! Independent double-reed blades use the same material-key network transaction.
use super::*;
use super::keyed::{CYCLE_STEPS, RELEASE, check_balance, drives, geometry, key, key_plate};
use fs_couple::bernoulli_aperture::network::boundaries::BoundaryDrive;
use fs_couple::bernoulli_aperture::paired::{
    ApertureBlade, BladePropagator, BladeState, PairedAperture, PairedApertureLimits, PairedApertureSpec,
};
use fs_couple::bernoulli_aperture::paired_network::PairedApertureNetwork;
use fs_couple::bernoulli_aperture::plate::PlateApertureReduction;

fn unkeyed(steps: u64) -> PairedApertureNetwork {
    let (graph, gases) = geometry();
    let z = graph.inlet_impedance_with_gases(&gases).unwrap();
    let spec = PairedApertureSpec { rest_opening_m: 4e-4, width_m: 0.013,
        blades: [ApertureBlade { mass_kg: 1e-5, stiffness_n_m: 350.0,
            damping_n_s_m: 0.03, pressure_area_m2: 1e-4 },
            ApertureBlade { mass_kg: 2e-5, stiffness_n_m: 700.0,
                damping_n_s_m: 0.04, pressure_area_m2: 1.2e-4 }],
        density_kg_m3: gases[0].density, impedance_pa_s_m3: z, time_step_s: DT, max_steps: steps,
        limits: PairedApertureLimits { maximum_abs_displacement_m: 0.01,
            maximum_abs_velocity_m_s: 100.0, maximum_energy_j: 1.0,
            max_evaluations: 256, velocity_tolerance_m_s: 1e-15, relative_tolerance: 1e-13 } };
    let reed = PairedAperture::new_with_propagator(spec, [BladeState::default(); 2],
        Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 2e9, 2.0,
            "manufactured double-reed shared contact".into()).unwrap(), BladePropagator::ExactZoh).unwrap();
    PairedApertureNetwork::new(reed, graph, gases, &CancelGate::new_clock_free()).unwrap()
}

fn make(plate: &PlateApertureReduction, steps: u64) -> PairedApertureNetwork {
    let network = unkeyed(steps);
    let boundary = key(plate, network.section_gases()[2],
        network.represented_sections()[2].impedance_pa_s_m3, steps);
    network.with_boundary_apertures(vec![boundary], 1).unwrap()
}

#[test]
fn asymmetric_double_reed_and_material_key_share_the_same_pressure_work_and_ringdown() {
    let plate = key_plate();
    let mut model = make(&plate, CYCLE_STEPS);
    assert_eq!(model.boundary_apertures().unwrap().gas(3), Some(&model.section_gases()[2]));
    let (mut closed, mut reopened, mut asymmetric) = (false, false, false);
    let (mut key_work, mut blade_work) = (0.0, 0.0);
    for step in 0..CYCLE_STEPS {
        let before = model.stored_energy_j();
        if step == RELEASE {
            assert!(before > 0.0);
            assert!(model.boundary_apertures().unwrap().aperture(3).unwrap().relaxation().unwrap()
                .memory_sqrt_j().iter().any(|x| *x != 0.0));
        }
        let (drive, force, keys) = drives(step);
        let frame = model.step_with_boundaries(drive, [force, -0.7 * force], &keys).unwrap();
        check_balance(frame.coupled, before);
        let boundary = model.boundary_apertures().unwrap().aperture_frame(3).unwrap();
        assert_eq!(frame.mechanics.aperture, frame.coupled.aperture);
        assert_eq!(boundary.step, model.aperture().accepted_steps());
        assert_eq!(model.node_frame(3).unwrap().pressure_pa, boundary.bore_pressure_pa);
        assert!((model.stored_energy_j() - frame.coupled.stored_energy_j).abs()
            <= 8.0 * f64::EPSILON * frame.coupled.stored_energy_j.max(1e-20));
        asymmetric |= frame.mechanics.blades[0] != frame.mechanics.blades[1];
        key_work += frame.coupled.boundaries.unwrap().mechanical_work_j.abs();
        blade_work += frame.coupled.aperture.mechanical_work_j.abs();
        if (512..RELEASE).contains(&step) && boundary.midpoint_opening_m == 0.0 {
            closed = true;
            assert_eq!(boundary.jet_flow_m3_s, 0.0);
        }
        if step >= RELEASE {
            assert!(frame.coupled.stored_energy_j <= before + 3e-8 * before.max(1e-20));
            reopened |= boundary.midpoint_opening_m > 1e-4;
        }
    }
    assert!(closed && reopened && asymmetric);
    assert!(key_work > 0.0 && blade_work > 0.0);
    assert!(model.boundary_apertures().unwrap().aperture(3).unwrap()
        .relaxation().unwrap().dissipated_energy_j() > 0.0);
}

#[test]
fn paired_key_binding_refusal_and_budget_resume_preserve_all_three_mechanical_states() {
    let plate = key_plate();
    let cold = unkeyed(64);
    let wrong_gas = key(&plate, cold.section_gases()[0],
        cold.represented_sections()[2].impedance_pa_s_m3, 64);
    assert!(cold.with_boundary_apertures(vec![wrong_gas], 1).is_err());
    let cold = unkeyed(64);
    let wrong_budget = key(&plate, cold.section_gases()[2],
        cold.represented_sections()[2].impedance_pa_s_m3, 63);
    assert!(cold.with_boundary_apertures(vec![wrong_budget], 1).is_err());
    let mut actual = make(&plate, 64);
    let mut reference = make(&plate, 96);
    let (drive, force, keys) = drives(600);
    let forces = [force, -0.7 * force];
    for _ in 0..32 {
        assert_eq!(actual.step_with_boundaries(drive, forces, &keys).unwrap(),
            reference.step_with_boundaries(drive, forces, &keys).unwrap());
    }
    let blades = *actual.aperture().blade_states();
    let energy = actual.stored_energy_j();
    let nodes: Vec<_> = (0..4).map(|node| *actual.node_frame(node).unwrap()).collect();
    let boundary = actual.boundary_apertures().unwrap().aperture(3).unwrap();
    let key_state = boundary.state();
    let memory = boundary.relaxation().unwrap().memory_sqrt_j().to_vec();
    assert!(actual.step_with_boundaries(drive, [0.0, f64::NAN], &keys).is_err());
    assert!(actual.step_with_boundaries(drive, forces, &[]).is_err());
    assert!(actual.step_with_boundaries(drive, forces,
        &[BoundaryDrive { exterior_pressure_pa: f64::NAN, ..keys[0] }]).is_err());
    assert_eq!(actual.aperture().blade_states(), &blades);
    assert_eq!(actual.aperture().accepted_steps(), 32);
    assert_eq!(actual.stored_energy_j(), energy);
    let boundary = actual.boundary_apertures().unwrap().aperture(3).unwrap();
    assert_eq!(boundary.state(), key_state);
    assert_eq!(boundary.accepted_steps(), 32);
    assert_eq!(boundary.relaxation().unwrap().memory_sqrt_j(), memory);
    for (node, frame) in nodes.iter().enumerate() { assert_eq!(actual.node_frame(node), Some(frame)); }
    for _ in 32..64 {
        assert_eq!(actual.step_with_boundaries(drive, forces, &keys).unwrap(),
            reference.step_with_boundaries(drive, forces, &keys).unwrap());
    }
    assert!(actual.step_with_boundaries(drive, forces, &keys).is_err());
    assert!(actual.extend_step_budget(64).is_err());
    assert_eq!(actual.aperture().accepted_steps(), 64);
    assert_eq!(actual.boundary_apertures().unwrap().aperture(3).unwrap().accepted_steps(), 64);
    actual.extend_step_budget(96).unwrap();
    assert_eq!(actual.boundary_apertures().unwrap().aperture(3).unwrap().spec().max_steps, 96);
    for _ in 64..96 {
        assert_eq!(actual.step_with_boundaries(drive, forces, &keys).unwrap(),
            reference.step_with_boundaries(drive, forces, &keys).unwrap());
    }
}
