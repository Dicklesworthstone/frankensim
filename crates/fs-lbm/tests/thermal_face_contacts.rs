//! Synthetic checks of the actual spatial energy solver, not material measurements.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, CompactComponent, ContactResistance, EnergyConfig, FaceContactResistance,
    FlowFace, FlowField, FluidProperties, SolidMaterial, ThermalFace, ThermalSetup,
    Voxel, VoxelDomain, solve_energy,
};

fn fluid() -> FluidProperties {
    FluidProperties { density_kg_m3: 1.0, specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0, kinematic_viscosity_m2_s: 1.0 }
}
fn flux_boundary() -> ThermalSetup {
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::HeatFlux(4.0);
    faces[1] = ThermalFace::Temperature(300.0);
    ThermalSetup::new(faces)
}
fn joint(a: usize, b: usize, r: f64) -> FaceContactResistance {
    FaceContactResistance { cells: [a, b], resistance_m2_k_w: r }
}

#[test]
fn g1_g3_same_material_joints_obey_graded_series_resistance_and_reverse_identically() {
    let widths = [0.5, 1.0, 2.0, 0.5];
    let domain = VoxelDomain::graded_from_fn(
        [widths.to_vec(), vec![2.0], vec![3.0]], |_| Voxel::Solid(0),
    ).unwrap();
    let before = domain.clone();
    let flow = FlowField::quiescent(&domain);
    let k = [1.0, 2.0, 3.0, 4.0];
    let field = k.map(|x| [x, 2.0*x, 3.0*x]);
    let mut setup = flux_boundary();
    setup.face_contacts = vec![joint(0, 1, 0.5), joint(2, 3, 1.0)];
    let config = EnergyConfig::default();
    let gate = CancelGate::new();
    let solve = |setup: &ThermalSetup| config.solve_with_conductivity(
        &domain, &fluid(), &[SolidMaterial::new("same material", 1.0)],
        &flow, setup, &field, &gate,
    ).unwrap();
    let result = solve(&setup);
    for c in 0..4 {
        let bulk = 0.5*widths[c]/k[c] + ((c+1)..4).map(|n| widths[n]/k[n]).sum::<f64>();
        let joints = if c == 0 { 1.5 } else if c <= 2 { 1.0 } else { 0.0 };
        assert!((result.temperature[c] - (300.0 + 4.0*(bulk + joints))).abs() < 1e-8);
    }
    assert!(result.report.balance.relative_residual < 1e-9);
    setup.face_contacts = vec![joint(3, 2, 1.0), joint(1, 0, 0.5)];
    assert_eq!(result, solve(&setup));
    assert_eq!(domain, before);
}

#[test]
fn g1_solid_to_fluid_observer_uses_the_same_resistive_face_as_assembly() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |p| {
        if p[0] < 1.0 { Voxel::Solid(0) } else { Voxel::Fluid }
    }).unwrap();
    let mut setup = flux_boundary();
    setup.faces[0] = ThermalFace::Temperature(310.0);
    setup.face_contacts.push(joint(0, 1, 2.0));
    let solution = solve_energy(&domain, &fluid(), &[SolidMaterial::new("solid", 1.0)],
        &FlowField::quiescent(&domain), &setup, &EnergyConfig::default(), &CancelGate::new()).unwrap();
    assert!((solution.temperature[0] - 308.75).abs() < 1e-9);
    assert!((solution.temperature[1] - 301.25).abs() < 1e-9);
    assert!((solution.solid_to_fluid_heat_w(&domain) - 2.5).abs() < 1e-9);
    assert!(solution.report.balance.relative_residual < 1e-9);
}

#[test]
fn g1_g3_transient_joint_exchanges_energy_without_a_reservoir() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.face_contacts.push(joint(0, 1, 3.0));
    let original = setup.clone();
    let dt = 0.2;
    let (solution, record) = EnergyConfig::default().step_with_conductivity(
        &domain, &fluid(), &[SolidMaterial::new("solid", 1.0)],
        &FlowField::quiescent(&domain), &setup, &[[1.0; 3]; 2],
        &[2.0, 5.0], &[320.0, 300.0], dt, &CancelGate::new(),
    ).unwrap();
    let heat = dt * 0.25 * 20.0 / (1.0 + dt * 0.25 * (0.5 + 0.2));
    assert!((solution.temperature[0] - (320.0 - heat/2.0)).abs() < 1e-9);
    assert!((solution.temperature[1] - (300.0 + heat/5.0)).abs() < 1e-9);
    assert!(record.closure_j.abs() < 1e-9);
    assert_eq!(record.source_j, 0.0);
    assert_eq!(record.boundary_outflow_j, 0.0);
    assert_eq!(setup, original);
}

#[test]
fn g5_pair_and_face_contacts_share_identical_flux_arithmetic() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |p| Voxel::Solid(u16::from(p[0] > 1.0))).unwrap();
    let solids = [SolidMaterial::new("a", 2.0), SolidMaterial::new("b", 8.0)];
    let mut setup = flux_boundary();
    setup.contacts.push(ContactResistance { materials: (0, 1), resistance_m2_k_w: 0.7 });
    let solve = |setup: &ThermalSetup| solve_energy(&domain, &fluid(), &solids,
        &FlowField::quiescent(&domain), setup, &EnergyConfig::default(), &CancelGate::new()).unwrap();
    let pair = solve(&setup);
    setup.contacts.clear();
    setup.face_contacts.push(joint(0, 1, 0.7));
    let face = solve(&setup);
    assert_eq!(pair.temperature, face.temperature);
    assert_eq!(pair.report, face.report);
    assert_eq!(pair.conductivity, face.conductivity);
}

#[test]
fn g0_ambiguous_and_nonphysical_faces_refuse() {
    let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |p| Voxel::Solid(u16::from(p[0] > 1.0))).unwrap();
    let solids = [SolidMaterial::new("a", 1.0), SolidMaterial::new("b", 1.0)];
    let check = |setup: &ThermalSetup| solve_energy(&domain, &fluid(), &solids,
        &FlowField::quiescent(&domain), setup, &EnergyConfig::default(), &CancelGate::new());
    for contacts in [vec![joint(0, 2, 1.0)], vec![joint(0, 3, 1.0)],
        vec![joint(0, 0, 1.0)], vec![joint(0, 1, 1.0), joint(1, 0, 1.0)],
        vec![joint(0, 1, -1.0)], vec![joint(0, 1, f64::NAN)],
        vec![joint(0, 1, f64::INFINITY)]]
    {
        let mut setup = flux_boundary();
        setup.face_contacts = contacts;
        assert!(matches!(check(&setup), Err(ChtError::InvalidInput { field: "thermal.face_contacts", .. })));
    }
    let mut setup = flux_boundary();
    setup.face_contacts.push(joint(0, 1, 1.0));
    setup.contacts.push(ContactResistance { materials: (0, 1), resistance_m2_k_w: 1.0 });
    assert!(check(&setup).is_err());
    setup.contacts.clear();
    setup.compact_components.push(CompactComponent { lo: [1,0,0], hi: [2,1,1],
        board_face: fs_lbm::Face3::ALL[0], power_w: 0.0,
        junction_to_case_k_w: 1.0, junction_to_board_k_w: 1.0, heat_capacity_j_k: None });
    assert!(check(&setup).is_err());
}

#[test]
fn g0_g4_through_flow_and_requested_cancellation_are_not_ignored() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let mut boundaries = [FlowFace::Wall; 6];
    boundaries[0] = FlowFace::Fixed;
    boundaries[1] = FlowFace::Fixed;
    let flow = FlowField::from_face_velocity(&domain, boundaries,
        |axis, _| if axis == 0 { 1.0 } else { 0.0 }).unwrap();
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.faces[0] = ThermalFace::Inflow { temperature: 300.0 };
    setup.faces[1] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    setup.face_contacts.push(joint(0, 1, 1.0));
    assert!(matches!(solve_energy(&domain, &fluid(), &[], &flow, &setup,
        &EnergyConfig::default(), &CancelGate::new()),
        Err(ChtError::InvalidInput { field: "thermal.face_contacts", .. })));
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(solve_energy(&domain, &fluid(), &[], &FlowField::quiescent(&domain),
        &flux_boundary(), &EnergyConfig::default(), &gate), Err(ChtError::Cancelled));
}
