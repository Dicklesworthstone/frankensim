//! Explicit compact mouths on original source faces. These synthetic soft
//! specimens exercise the shared owners; they are not full instrument models.
use super::*;
use fs_couple::air_path::cylinder::CylinderFlowResolution;
use fs_couple::bowed_string::runtime::reciprocal::bank::cavity::{
    ShellCavityOpening, ShellCavityOptions,
};
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::shear_shell::pressure::ShellPressureSurfaces;
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace};

fn with_volume<R>(atmosphere: GasState,
    f: impl FnOnce(&HelmholtzFe<'_>, [u32;3]) -> R) -> R
{
    let surface = SpherePanels::icosphere(0.1,1).unwrap();
    let mut indices = BTreeMap::new();
    let mut positions = Vec::new();
    let faces: Vec<[u32;3]> = surface.triangles().unwrap().iter().map(|triangle|
        triangle.map(|point| {
            let key = point.map(|x| if x == 0.0 { 0 } else { x.to_bits() });
            *indices.entry(key).or_insert_with(|| {
                positions.push(point); (positions.len()-1) as u32
            })
        })).collect();
    // The opposite pole is the sole moving nodal translation. This mouth's
    // actual entire rim is constrained, as required by the fixed-neck chart.
    let face = *faces.iter().find(|face| face.iter()
        .all(|&v| positions[v as usize][0] < 0.01)).unwrap();
    let center = positions.len() as u32;
    positions.push([0.0;3]);
    let mut tets: Vec<_> = faces.into_iter().map(|[a,b,c]|[center,a,b,c]).collect();
    for tet in &mut tets {
        let [a,b,c] = [tet[1],tet[2],tet[3]].map(|v|positions[v as usize]);
        let sign = a[0]*(b[1]*c[2]-b[2]*c[1])
            - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]);
        if sign < 0.0 { tet.swap(2,3); }
    }
    let cells = fs_rep_mesh::TetComplex::from_tets(positions.len(),tets);
    let fe = HelmholtzFe::new(&cells,&positions,1,Medium {
        density_kg_m3:atmosphere.density,sound_speed_m_s:atmosphere.sound_speed,
    }).unwrap();
    f(&fe,face)
}

fn body(fe: &HelmholtzFe<'_>, omitted: &[[u32;3]], rho: f64,
    gate: &CancelGate) -> SharedMaterialShellBody
{
    let mut missing: Vec<_> = omitted.iter().map(|face| {
        let mut key = *face; key.sort_unstable(); key
    }).collect();
    missing.sort_unstable();
    let mut triangles = Vec::new();
    for face in fe.boundary_faces() {
        let mut key = face; key.sort_unstable();
        if missing.binary_search(&key).is_ok() { continue; }
        let (p,normal,_) = fe.boundary_geometry(face).unwrap();
        let a = [p[1][0]-p[0][0],p[1][1]-p[0][1],p[1][2]-p[0][2]];
        let b = [p[2][0]-p[0][0],p[2][1]-p[0][1],p[2][2]-p[0][2]];
        let cross = [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
        let mut triangle = face.map(|v|v as usize);
        if cross.iter().zip(normal).map(|(a,b)|a*b).sum::<f64>() < 0.0 { triangle.swap(1,2); }
        triangles.push(triangle);
    }
    // All original surface vertices remain; only the original center is gas.
    let nodes = fe.positions()[..fe.positions().len()-1].to_vec();
    let moving = (0..nodes.len()).max_by(|&a,&b|nodes[a][0].total_cmp(&nodes[b][0])).unwrap();
    let vectors = nodes.iter().map(|p| {
        let r = p[0].hypot(p[1]).hypot(p[2]); p.map(|x|x/r)
    }).collect();
    let mesh = ShellMesh::new(nodes,triangles).unwrap();
    let directors = NodalDirectors::new(&mesh,vectors,0.05).unwrap();
    let constraints = (0..mesh.nodes.len()).map(|node| {
        let mut mask = [true;6]; if node == moving { mask[0] = false; }
        ShellConstraint {node,mask,frame:ConstraintFrame::Global}
    }).collect();
    let budget = MitcBudget {max_nodes:mesh.nodes.len(),max_triangles:mesh.tris.len()};
    let region = ShellRegionMaterial {region:PlateRegionMaterial {
        region:PlateRegion {name:"skin".into(),triangle_indices:(0..mesh.tris.len()).collect()},
        material:soft_material_scaled(rho,0.01),model:PlateMaterialModel::Orthotropic12,
        thickness_constraint:PlateThicknessConstraint::FixedThickness(0.003),material_angle_rad:0.17,
    },shear_correction:5.0/6.0,axis1_global:[1.0,0.173,0.519]};
    let facets: Vec<_> = (0..mesh.tris.len()).map(|triangle|
        OrientedShellFacet {triangle,side:ShellFace::Positive}).collect();
    let specimen = compile_material_shear_shell(mesh,vec![region],budget,gate).unwrap()
        .with_pressure_surfaces(ShellPressureSurfaces {exterior:facets.clone(),cavity:facets},gate).unwrap();
    let mut port = vec![0.0;6*specimen.mesh().nodes.len()]; port[6*moving] = 1.0;
    SharedMaterialShellBody {
        specimen,directors,constraints,mass:ShellMass::Consistent,drilling_alpha:1e-3,
        eigenvalue_window:(0.0,1e9),n_modes:1,force_ports:vec![port],
        bridges:vec![SharedPlateBridge {station_fraction:0.8,stiffness_n_m:3000.0,
            damping_n_s_m:0.1,rest_extension_m:0.0}],
        assembly_budget:budget,
        projection_budget:ShellProjectionBudget {max_modes:1,max_projection_terms:128_000_000,
            reference_resistance_n_s_m:1.0},
        modal_budget:ModalAcousticTimeBudget::audible_reference(),listener_m:1.0,
    }
}
fn bank(cx: &Cx<'_>, gate: &CancelGate, fe: &HelmholtzFe<'_>,
    omitted: &[[u32;3]], atmosphere: GasState, rho: f64) -> ReciprocalStringBank
{
    let mut settings = options(0);
    settings.network.max_modes = 64; settings.network.max_connections = 64;
    settings.network.max_setup_terms = 1_048_576;
    ReciprocalStringBank::from_material_shear_shell(cx,vec![member(7800.0).source],
        atmosphere,body(fe,omitted,rho,gate),None,law(),settings,gate).unwrap()
}
fn interior() -> ShellCavityOptions {
    ShellCavityOptions {
        maximum_hz:1.0,eigensolver:fs_modal::SliceOptions::default(),
        maximum_pressure_modes:1,maximum_dofs:64,maximum_boundary_faces:80,
        maximum_projection_terms:10_000_000,
        acoustic_budget:ModalAcousticTimeBudget::audible_reference(),
    }
}
fn opening(face:[u32;3]) -> ShellCavityOpening {
    ShellCavityOpening {boundary_faces:vec![face],internal_length_m:0.04,
        internal_resistance_pa_s_m3:1000.0}
}
fn exterior() -> ShellBemOptions {
    ShellBemOptions {minimum_hz:1.0,maximum_hz:300.0,frequency_samples:65,
        maximum_surface_panels:80,maximum_sweep_terms:100_000_000,maximum_poles:64}
}
fn load() -> Model {
    // Explicit synthetic passive matrix for an independent work oracle; only
    // the BEM test below claims a load derived from the actual source geometry.
    Model {ports:2,poles:vec![
        Pole {omega:500.0,zeta:0.2,coupling:vec![120.0,-80.0]},
        Pole {omega:1400.0,zeta:0.15,coupling:vec![-30.0,70.0]},
    ]}
}

#[test]
fn original_face_integrals_set_neck_mass_and_uniform_compression_after_gas_changes() {
    with_cx(|cx,gate| {
        for pressure in [101_325.0,202_650.0] {
            let atmosphere = gas(pressure);
            with_volume(atmosphere,|fe,face| {
                let opened = bank(cx,gate,fe,&[face],atmosphere,700.0)
                    .with_material_shell_cavity_openings(fe,interior(),vec![opening(face)],gate).unwrap();
                let cavity = opened.source().spatial_cavity().unwrap();
                assert_eq!(cavity.opening(0),Some(&opening(face)));
                assert_eq!(cavity.coupling().cavity_modes(),1);
                assert_eq!(cavity.coupling().total_modes(),2);
                assert!(cavity.coupling().neck_accepts_exterior_load(0).unwrap());
                let port = cavity.coupling().neck_radiation_port(0).unwrap();
                let (_,area) = fe.boundary_load(&[face]).unwrap();
                assert!((port.area_m2-area).abs() < 1e-14*area);
                let weight = (area/(atmosphere.density*0.04)).sqrt();
                assert!((port.volume_weight_m2_per_sqrt_kg/weight-1.0).abs() < 1e-14);
                let reference = atmosphere.sound_speed*(area/(cavity.volume_m3()*0.04)).sqrt();
                let actual = cavity.coupling().neck_fixed_wall_angular_frequency_rad_s(0).unwrap();
                assert!((actual/reference-1.0).abs() < 1e-12);
                assert_eq!(opened.spatial_cavity_neck_exterior_pressure(0).unwrap(),0.0);
                assert_eq!(opened.source().material_shell_image().unwrap().source().mesh().tris.len(),79);
                assert_eq!(fe.boundary_face_count(),80);
            });
        }
    });
}

#[test]
fn mouth_exterior_pressure_reacts_in_the_same_string_body_air_energy_balance() {
    with_cx(|cx,gate| with_volume(gas(101_325.0),|fe,face| {
        let build = || bank(cx,gate,fe,&[face],gas(101_325.0),700.0)
            .with_material_shell_cavity_openings(fe,interior(),vec![opening(face)],gate).unwrap();
        let mut bare = build();
        let mut loaded = build().with_body_radiation_load(load(),2,gate).unwrap();
        let cavity = loaded.source().spatial_cavity().unwrap();
        let port = cavity.coupling().neck_radiation_port(0).unwrap();
        assert_eq!(loaded.source().body_radiation().unwrap().source_body_coordinates(),&[0,port.coordinate]);
        let mut scratch = vec![0.0;cavity.observation_state_len()];
        let mut force = vec![0.0;loaded.system().mode_count()]; force[4] = 0.02;
        let released = FrictionDrive {speed_m_s:0.0,normal_force_n:0.0};
        let mut peak_pressure = 0.0_f64;
        let mut peak_flow = 0.0_f64;
        for _ in 0..96 {
            let frame = *loaded.step(&force,released).unwrap();
            bare.step(&force,released).unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            let neck = loaded.spatial_cavity_neck_observation(0,&mut scratch).unwrap();
            let pressure = loaded.spatial_cavity_neck_exterior_pressure(0).unwrap();
            let state = loaded.system().radiation_state();
            let reaction = -load().poles.iter().zip(state.chunks_exact(2))
                .map(|(pole,state)|pole.coupling[1]*state[1]).sum::<f64>();
            let power = pressure*neck.volume_flow_m3_s;
            let work = reaction*scratch[2*port.coordinate+1];
            assert!((power+work).abs() <= 1e-12*(1.0+power.abs()+work.abs()));
            peak_pressure = peak_pressure.max(pressure.abs());
            peak_flow = peak_flow.max(neck.volume_flow_m3_s.abs());
        }
        assert!(peak_pressure > 0.0 && peak_flow > 0.0);
        assert_ne!(loaded.system().components()[0].states(),bare.system().components()[0].states());
        assert_ne!(loaded.system().components()[2].states(),bare.system().components()[2].states());
        assert!(loaded.system().radiation_energy_j().unwrap() > 0.0);
        assert!(loaded.system().coupled_frame().unwrap().radiation_dissipation_j > 0.0);
        let mut zero = load();
        for pole in &mut zero.poles { pole.coupling.fill(0.0); }
        let mut silent = build().with_body_radiation_load(zero,2,gate).unwrap();
        let mut unloaded = build();
        for _ in 0..8 {
            assert_eq!(silent.step(&force,released).unwrap(),unloaded.step(&force,released).unwrap());
            assert_eq!(silent.spatial_cavity_neck_exterior_pressure(0).unwrap(),0.0);
        }
    }));
}

#[test]
fn original_bem_boundary_closes_with_a_fluid_mouth_and_keeps_mutual_reaction() {
    // Pending native integration gate: independent localized columns expose
    // collocation reciprocity error. A source refusal is a failure of this
    // qualification, never permission to symmetrize or weaken the 1e-12 gate.
    with_cx(|cx,gate| with_volume(gas(101_325.0),|fe,face| {
        let (loaded,report) = bank(cx,gate,fe,&[face],gas(101_325.0),700.0)
            .with_material_shell_cavity_and_bem(fe,interior(),vec![opening(face)],exterior(),gate).unwrap();
        assert_eq!(report.body_modes,1); assert_eq!(report.opening_ports,1);
        assert_eq!(report.surface_panels,80);
        assert!(report.reciprocity_relative_residual <= 1e-12);
        let binding = loaded.source().body_radiation().unwrap();
        assert_eq!(binding.source_body_coordinates(),&[0,1]);
        assert_eq!(binding.model().ports,2);
        assert_eq!(binding.material_shell_bem().unwrap().loading_identity,report.loading_identity);
        let z = binding.model().impedance(2.0*core::f64::consts::PI*150.0).unwrap();
        assert!((z[1]-z[2]).abs() <= 1e-12*(1.0+z[1].abs()));
        assert!(z[1].abs() > 0.0,"body-mouth pressure reaction must not be diagonalized away");
        assert!(z[3].abs() > 0.0,"physical neck receives exterior inertia and radiation");
    }));
}

#[test]
fn rebuilt_string_gas_states_preserve_the_original_body_and_neck_source_addresses() {
    with_cx(|cx,gate| with_volume(gas(101_325.0),|fe,face| {
        let flow = CylinderFlowResolution {minimum_wavenumber_radius:0.01,
            maximum_wavenumber_radius:4.0,intervals:2};
        let build = || bank(cx,gate,fe,&[face],gas(101_325.0),700.0);
        let mut first = build().with_material_shell_cavity_openings(fe,interior(),vec![opening(face)],gate)
            .unwrap().with_body_radiation_load(load(),2,gate).unwrap()
            .with_cylinder_gas_loading(flow,gate).unwrap();
        let mut second = build().with_cylinder_gas_loading(flow,gate).unwrap()
            .with_material_shell_cavity_openings(fe,interior(),vec![opening(face)],gate).unwrap()
            .with_body_radiation_load(load(),2,gate).unwrap();
        assert_eq!(first.system().mode_count(),second.system().mode_count());
        assert_eq!(first.source().body_radiation().unwrap().source_body_coordinates(),&[0,1]);
        assert_eq!(first.source().body_connections(),second.source().body_connections());
        let mut force = vec![0.0;first.system().mode_count()]; force[4] = 0.02;
        let released = FrictionDrive {speed_m_s:0.0,normal_force_n:0.0};
        for _ in 0..24 {
            assert_eq!(first.step(&force,released).unwrap(),second.step(&force,released).unwrap());
            assert_eq!(first.spatial_cavity_neck_exterior_pressure(0).unwrap(),
                second.spatial_cavity_neck_exterior_pressure(0).unwrap());
        }
    }));
}

#[test]
fn solid_caps_undeclared_faces_moving_rims_duplicate_mouths_and_partial_loads_refuse() {
    with_cx(|cx,gate| with_volume(gas(101_325.0),|fe,face| {
        let build = |removed:&[[u32;3]]|bank(cx,gate,fe,removed,gas(101_325.0),700.0);
        assert!(build(&[]).with_material_shell_cavity_openings(fe,interior(),vec![opening(face)],gate).is_err());
        assert!(build(&[face]).with_material_shell_cavity(fe,interior(),gate).is_err());
        assert!(build(&[face]).with_material_shell_cavity_openings(
            fe,interior(),vec![opening(face),opening(face)],gate).is_err());
        let moving = (0..fe.positions().len()-1)
            .max_by(|&a,&b|fe.positions()[a][0].total_cmp(&fe.positions()[b][0])).unwrap() as u32;
        let moving_face = fe.boundary_faces().into_iter().find(|face|face.contains(&moving)).unwrap();
        assert!(build(&[moving_face]).with_material_shell_cavity_openings(
            fe,interior(),vec![opening(moving_face)],gate).is_err());
        let opened = ||build(&[face]).with_material_shell_cavity_openings(
            fe,interior(),vec![opening(face)],gate).unwrap();
        let mut partial = load(); partial.ports = 1;
        for pole in &mut partial.poles { pole.coupling.pop(); }
        assert!(opened().with_body_radiation_load(partial,2,gate).is_err());
        let mut too_wide = exterior(); too_wide.maximum_hz = 1000.0;
        assert!(opened().prepare_material_shell_bem(too_wide,gate).is_err());
        let mut too_small = exterior(); too_small.maximum_surface_panels = 79;
        assert!(opened().prepare_material_shell_bem(too_small,gate).is_err());
    }));
}

