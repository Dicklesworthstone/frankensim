//! Synthetic checks of the actual spatial energy solver, not material measurements.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ChtError, CompactComponent, ContactResistance, EddyWallTreatment, EnergyConfig, FaceContactResistance,
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
fn g1_sgs_wall_heat_is_molecular_while_fluid_mixing_and_open_faces_keep_eddy_transport() {
    let domain = VoxelDomain::graded([vec![0.5, 1.0, 1.5], vec![1.0], vec![1.0]]).unwrap();
    let flow = FlowField::quiescent(&domain);
    for (face, wall_resistance) in [
        (ThermalFace::Temperature(300.0), 0.75),
        (ThermalFace::Convective { h: 2.0, ambient: 300.0 }, 1.25),
        (ThermalFace::Inflow { temperature: 300.0 }, 0.75 / 11.0),
    ] {
        let mut setup = flux_boundary();
        setup.faces[1] = face;
        setup.eddy_conductivity_w_m_k = vec![10.0; 3];
        setup.eddy_wall_treatment = EddyWallTreatment::Molecular;
        let solution = EnergyConfig::default().solve_with_conductivity(
            &domain, &fluid(), &[], &flow, &setup, &[[1.0, 2.0, 3.0]; 3],
            &CancelGate::new(),
        ).unwrap();
        let resistances = [wall_resistance + 2.0 / 11.0,
            wall_resistance + 1.25 / 11.0, wall_resistance];
        for (actual, resistance) in solution.temperature.iter().zip(resistances) {
            near(*actual, 300.0 + 4.0 * resistance, 1e-8);
        }
        assert_eq!(solution.conductivity, vec![[11.0, 12.0, 13.0]; 3]);
        assert_eq!(solution.wall_conductivity(), &[[1.0, 2.0, 3.0]; 3]);
        assert!(solution.report.balance.relative_residual < 1e-9);
    }
}

#[test]
fn g1_sgs_solid_interfaces_and_compact_surfaces_keep_the_molecular_half_cell() {
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |p| {
        if p[0] < 1.0 { Voxel::Solid(0) } else { Voxel::Fluid }
    }).unwrap();
    let flow = FlowField::quiescent(&domain);
    let mut setup = flux_boundary();
    setup.eddy_conductivity_w_m_k = vec![0.0, 10.0];
    setup.eddy_wall_treatment = EddyWallTreatment::Molecular;
    let solids = [SolidMaterial::new("synthetic solid", 2.0)];
    let solve = |setup: &ThermalSetup| solve_energy(&domain, &fluid(), &solids,
        &flow, setup, &EnergyConfig::default(), &CancelGate::new()).unwrap();
    let result = solve(&setup);
    near(result.temperature[0], 305.0, 1e-8);
    near(result.temperature[1], 302.0, 1e-8);
    near(result.solid_to_fluid_heat_w(&domain), 4.0, 1e-8);
    setup.eddy_wall_treatment = EddyWallTreatment::EffectiveConductivity;
    let legacy = solve(&setup);
    near(legacy.temperature[0], 301.0 + 4.0 / 11.0, 1e-8);
    near(legacy.temperature[1], 300.0 + 2.0 / 11.0, 1e-8);

    let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |p| {
        if p[0] > 1.0 && p[0] < 2.0 { Voxel::Solid(0) } else { Voxel::Fluid }
    }).unwrap();
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature = vec![(0, 300.0), (2, 300.0)];
    setup.eddy_conductivity_w_m_k = vec![10.0, 0.0, 10.0];
    setup.eddy_wall_treatment = EddyWallTreatment::Molecular;
    setup.compact_components.push(CompactComponent {
        lo: [1, 0, 0], hi: [2, 1, 1], board_face: fs_lbm::d3q19::Face3::XMin,
        power_w: 20.0, junction_to_case_k_w: 1.0, junction_to_board_k_w: 1.0,
        heat_capacity_j_k: None,
    });
    let result = solve_energy(&domain, &fluid(), &solids,
        &FlowField::quiescent(&domain), &setup, &EnergyConfig::default(),
        &CancelGate::new()).unwrap();
    // Two parallel 1 + 0.5 K/W paths, without an invented package bulk leg.
    near(result.junctions[0].temperature_k, 315.0, 1e-8);
    near(result.junctions[0].case_w, 10.0, 1e-8);
    near(result.junctions[0].board_w, 10.0, 1e-8);
    near(result.solid_to_fluid_heat_w(&domain), 20.0, 1e-8);
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
        board_face: fs_lbm::Face3::ALL[2], power_w: 0.0,
        junction_to_case_k_w: 1.0, junction_to_board_k_w: 1.0, heat_capacity_j_k: None });
    // A physical case/board contact is supported; this joint instead touches
    // the package's adiabatic x side and must not create a new heat path.
    assert!(matches!(check(&setup), Err(ChtError::InvalidInput {
        field: "thermal.face_contacts", ref reason,
    }) if reason.contains("compact sides")));
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

fn source_resistance(a: f64, b: f64) -> f64 {
    0.5 + (a - 320.0) / 1024.0 + 3.0 * (b - 320.0) / 1024.0
}

fn source_point(a: f64, b: f64) -> fs_matdb::QueryPoint {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    let kelvin = fs_qty::QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ));
    fs_matdb::QueryPoint::new().with_quantity("TA", kelvin, a).unwrap()
        .with_quantity("TB", kelvin, b).unwrap()
        .with_quantity("normal_pressure", fs_qty::QuantitySpec::dimensional(fs_qty::Pressure::DIMS), 1e6).unwrap()
}

fn source_card(a: [f64; 2], b: [f64; 2], reverse: bool) -> std::sync::Arc<fs_matdb::InterfaceSystemCard> {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, SurfaceSpec, SystemContext, UncertaintyModel};
    use fs_material::thermal::interface::{AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as DIMS,
        AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as PROPERTY};
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    let kelvin = fs_qty::QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let pascal = fs_qty::QuantitySpec::dimensional(fs_qty::Pressure::DIMS);
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(PROPERTY, DIMS),
        value: PropertyValue::Surface2D {
            axes: ["TA".into(), "TB".into()], axis_dims: [fs_qty::Temperature::DIMS; 2],
            samples: [[a[0],b[0]], [a[1],b[0]], [a[0],b[1]], [a[1],b[1]]]
                .into_iter().map(|p| (p, source_resistance(p[0],p[1]))).collect(),
            triangles: vec![[0,1,2],[1,3,2]], dims: DIMS,
        },
        // Each validity axis carries the quantity that source_point() gives it: a query
        // axis typed differently from its claim's axis is outside the claim's domain.
        validity: ValidityDomain::unconstrained().with_quantity("TA",kelvin,a[0],a[1])
            .with_quantity("TB",kelvin,b[0],b[1]).with_quantity("normal_pressure",pascal,1e6,1e6),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance {
            source: "manufactured ordered surface-resistance law; not measured material data".into(),
            license: "internal-test-use".into(), artifact: None,
        },
    }).unwrap();
    let surface = |name: &str| SurfaceSpec { material: MaterialStateId {
        chemistry: name.into(), phase: "fixture".into(), process: "test".into(), revision: 0,
    }, texture_frame: "fixed reference normal".into() };
    let [a,b] = if reverse { ["B","A"] } else { ["A","B"] };
    std::sync::Arc::new(InterfaceSystemCard::assemble(surface(a), surface(b), SystemContext {
        medium: "dry".into(), third_body: None, environment: "test".into(), history: "unaged".into(),
    }, claims, vec![]).unwrap())
}

fn source_contact(cells: [usize; 2], card: std::sync::Arc<fs_matdb::InterfaceSystemCard>,
    point: &fs_matdb::QueryPoint, reverse: bool,
) -> fs_lbm::conjugate::SourceFaceContact {
    fs_lbm::conjugate::SourceFaceContact::new(cells, card, point,
        if reverse { ["TB","TA"] } else { ["TA","TB"] },
        fs_material::state_point::MaterialPropertySelection::SingleClaimOnly).unwrap()
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a:.16e} != {b:.16e}, tolerance={tolerance:.3e}");
}

#[test]
fn g1_sgs_source_contact_queries_molecular_surfaces_not_eddy_shifted_ones() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let surfaces = [350.0, 320.0];
    let q = 30.0 / source_resistance(surfaces[0], surfaces[1]);
    let cells = [surfaces[0] + 0.5 * q, surfaces[1] - 0.5 * q];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature = vec![(0, cells[0]), (1, cells[1])];
    setup.eddy_conductivity_w_m_k = vec![10.0; 2];
    setup.eddy_wall_treatment = EddyWallTreatment::Molecular;
    setup.source_face_contacts.push(source_contact([0, 1],
        source_card([348.0, 352.0], [318.0, 322.0], false),
        &source_point(350.0, 320.0), false));
    let result = solve_energy(&domain, &fluid(), &[], &FlowField::quiescent(&domain),
        &setup, &EnergyConfig::default(), &CancelGate::new()).unwrap();
    let state = &result.source_contact_states[0];
    for (actual, expected) in state.surface_temperature_k.iter().zip(surfaces) {
        near(*actual, expected, 2e-7);
    }
    near(state.heat_rate_first_to_second_w, q, 2e-8);
    near(state.source.state_coordinate("TA").unwrap(), 350.0, 2e-7);
    near(state.source.state_coordinate("TB").unwrap(), 320.0, 2e-7);
    assert_eq!(result.conductivity, vec![[11.0; 3]; 2]);
}

#[test]
fn g1_source_contact_transient_matches_independent_surface_flux_and_closed_energy() {
    let domain = VoxelDomain::from_fn(2,1,1,1.0,|p| Voxel::Solid(u16::from(p[0]>1.0))).unwrap();
    let solids = [SolidMaterial::new("A",2.0),SolidMaterial::new("B",4.0)];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic;6]);
    setup.source_face_contacts.push(source_contact([0,1],
        source_card([256.,384.],[256.,384.],false), &source_point(320.,320.), false));
    let before = setup.clone();
    let previous = [360.,300.]; let capacity = [2.,5.]; let dt = 0.35;
    // Eliminate the two stored-energy equations. The actual surfaces include
    // each bulk half-cell drop; the source is never evaluated at cell centres.
    let surface = |q: f64| [previous[0]-dt*q/capacity[0]-0.25*q,
        previous[1]+dt*q/capacity[1]+0.125*q];
    let (mut lo,mut hi) = (0.0,60.0/(0.375+dt*(0.5+0.2)));
    for _ in 0..80 {
        let q = 0.5*(lo+hi); let s = surface(q);
        if s[0]-s[1]-source_resistance(s[0],s[1])*q > 0.0 {lo=q;} else {hi=q;}
    }
    let q = 0.5*(lo+hi); let expected = [previous[0]-dt*q/2.,previous[1]+dt*q/5.];
    let (solution,record) = EnergyConfig::default().step_with_conductivity(
        &domain,&fluid(),&solids,&FlowField::quiescent(&domain),&setup,
        &[[2.;3],[4.;3]],&capacity,&previous,dt,&CancelGate::new()).unwrap();
    for i in 0..2 {near(solution.temperature[i],expected[i],2e-8);}
    assert_eq!(solution.source_contact_states.len(),1);
    let state = &solution.source_contact_states[0];
    assert_eq!(state.cells,[0,1]);
    for i in 0..2 {near(state.surface_temperature_k[i],surface(q)[i],2e-8);}
    near(state.heat_rate_first_to_second_w,q,2e-8);
    near(state.resistance_m2_k_w,source_resistance(surface(q)[0],surface(q)[1]),2e-11);
    near(state.applied_resistance_m2_k_w,
        solution.face_contact_resistance_m2_k_w[&(0,1)],0.0);
    for (axis,temperature) in [("TA",state.surface_temperature_k[0]),("TB",state.surface_temperature_k[1])] {
        near(state.source.state_coordinate(axis).unwrap(),temperature,0.0);
    }
    near(state.source.state_coordinate("normal_pressure").unwrap(),1e6,0.0);
    near(state.contact_temperature_residual_k,0.0,
        1.1*setup.source_contact_controls.temperature_tolerance_k);
    assert!(state.resistance_relative_residual <= setup.source_contact_controls.resistance_relative_tolerance);
    assert!(solution.source_contact_iterations > 1);
    near(record.stored_energy_change_j,0.,1e-8);
    near(record.closure_j,0.,1e-8);
    near(record.source_j,0.,0.0);near(record.boundary_outflow_j,0.,0.0);
    assert_eq!(previous,[360.,300.]);assert_eq!(setup,before);
}

#[test]
fn g1_g3_source_contact_uses_nonlinear_solid_surfaces_and_consistent_fluid_heat() {
    let domain = VoxelDomain::graded_from_fn([vec![0.5,0.8],vec![0.2],vec![0.3]],
        |p| if p[0]<0.5 {Voxel::Solid(0)} else {Voxel::Fluid}).unwrap();
    let solid = SolidMaterial::new("A",1.0)
        .with_conductivity_table(&[(250.,0.5),(330.,1.3),(370.,1.7),(450.,2.5)]);
    let fluid = FluidProperties {conductivity_w_m_k:1.25,..fluid()};
    let potential = |t: f64| {let x=t-300.;x+0.005*x*x};
    let inverse = |u: f64| 300.+2.*u/(1.+(1.+0.02*u).sqrt());
    let surfaces = [350.,320.];
    let q = 30./source_resistance(surfaces[0],surfaces[1]);
    let solid_cell = inverse(potential(surfaces[0])+q*0.25);
    let fluid_cell = surfaces[1]-q*0.4/2.;
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic;6]);
    setup.fixed_temperature = vec![(0,solid_cell)];
    setup.faces[1] = ThermalFace::Temperature(fluid_cell-q*0.4/2.);
    setup.eddy_conductivity_w_m_k = vec![0.,0.75];
    let mut physical = None;
    for reverse in [false,true] {
        setup.source_face_contacts = vec![source_contact(if reverse {[1,0]}else{[0,1]},
            source_card([348.,352.],[318.,322.],reverse),&source_point(350.,320.),reverse)];
        let before = setup.clone();
        let solution = solve_energy(&domain,&fluid,std::slice::from_ref(&solid),
            &FlowField::quiescent(&domain),&setup,&EnergyConfig::default(),&CancelGate::new()).unwrap();
        near(solution.temperature[0],solid_cell,1e-8);
        near(solution.temperature[1],fluid_cell,2e-7);
        assert!(solid_cell>352. && fluid_cell<318.,"bulk centres are outside the interface source box");
        let state=&solution.source_contact_states[0];
        let expected=if reverse {[320.,350.]}else{[350.,320.]};
        for i in 0..2 {near(state.surface_temperature_k[i],expected[i],2e-7);}
        near(state.source.state_coordinate("TA").unwrap(),350.,2e-7);
        near(state.source.state_coordinate("TB").unwrap(),320.,2e-7);
        near(state.heat_rate_first_to_second_w,if reverse {-0.06*q}else{0.06*q},2e-8);
        near(solution.solid_to_fluid_heat_w(&domain),0.06*q,2e-8);
        assert!(solution.report.balance.relative_residual<1e-9);
        if let Some(ref previous) = physical {
            let previous: &Vec<f64> = previous;
            for (&a,&b) in previous.iter().zip(&solution.temperature) {near(a,b,1e-9);}
        }
        physical=Some(solution.temperature.clone());
        assert_eq!(setup,before);
    }
}

#[test]
fn g0_g4_source_contact_refusals_preserve_setup_and_transient_history() {
    use fs_lbm::conjugate::SourceFaceContact;
    use fs_material::state_point::MaterialPropertySelection;
    let domain=VoxelDomain::from_fn(2,1,1,1.0,|_|Voxel::Solid(0)).unwrap();
    let solids=[SolidMaterial::new("A",1.0)];
    let point=source_point(320.,320.);
    let card=source_card([256.,384.],[256.,384.],false);
    for axes in [["TA","TA"],["TA","missing"]] {
        assert!(SourceFaceContact::new([0,1],card.clone(),&point,axes,
            MaterialPropertySelection::SingleClaimOnly).is_err());
    }
    let mut setup=ThermalSetup::new([ThermalFace::Adiabatic;6]);
    setup.source_face_contacts.push(source_contact([0,1],card,&point,false));
    let initial=[360.,300.];
    let solve=|setup:&ThermalSetup,gate:&CancelGate| EnergyConfig::default().step_with_conductivity(
        &domain,&fluid(),&solids,&FlowField::quiescent(&domain),setup,
        &[[1.;3];2],&[2.,5.],&initial,0.35,gate);
    let mut duplicate=setup.clone();duplicate.face_contacts.push(joint(1,0,0.2));
    assert!(solve(&duplicate,&CancelGate::new()).is_err());
    let mut unsupported=setup.clone();
    unsupported.source_face_contacts[0]=source_contact([0,1],
        source_card([319.,321.],[319.,321.],false),&point,false);
    assert!(solve(&unsupported,&CancelGate::new()).is_err());
    let mut exhausted=setup.clone();exhausted.source_contact_controls.max_iterations=1;
    assert!(matches!(solve(&exhausted,&CancelGate::new()),Err(ChtError::SolverNotConverged{..})));
    let before=setup.clone();let gate=CancelGate::new();gate.request();
    assert_eq!(solve(&setup,&gate),Err(ChtError::Cancelled));
    assert_eq!(setup,before);assert_eq!(initial,[360.,300.]);
}
