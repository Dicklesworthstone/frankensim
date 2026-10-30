//! Actual source shell/pressure FE coupling; this synthetic fixture is not
//! a measured instrument or a claim of converged full-band acoustics.
use super::*;
use fs_couple::air_path::cylinder::CylinderFlowResolution;
use fs_couple::bowed_string::runtime::reciprocal::bank::cavity::ShellCavityOptions;
use fs_couple::shear_shell::pressure::ShellPressureSurfaces;
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace};

#[path = "cavity/thermal.rs"]
mod thermal;

fn with_volume<R>(n: usize, atmosphere: GasState, f: impl FnOnce(&HelmholtzFe<'_>) -> R) -> R {
    assert!((1..=2).contains(&n));
    // A closed icosahedral shell has strictly acute incident-facet angles.
    // Mild unequal axis scales separate symmetric eigenvalues without changing
    // its conforming topology or weakening the original director admission.
    let golden = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut positions: Vec<[f64; 3]> = [
        [-1.0,golden,0.0], [1.0,golden,0.0], [-1.0,-golden,0.0], [1.0,-golden,0.0],
        [0.0,-1.0,golden], [0.0,1.0,golden], [0.0,-1.0,-golden], [0.0,1.0,-golden],
        [golden,0.0,-1.0], [golden,0.0,1.0], [-golden,0.0,-1.0], [-golden,0.0,1.0],
    ].into_iter().map(|p| [0.1*p[0],0.11*p[1],0.09*p[2]]).collect();
    let center = positions.len() as u32;
    positions.push([0.0;3]);
    let faces = [
        [0,11,5], [0,5,1], [0,1,7], [0,7,10], [0,10,11],
        [1,5,9], [5,11,4], [11,10,2], [10,7,6], [7,1,8],
        [3,9,4], [3,4,2], [3,2,6], [3,6,8], [3,8,9],
        [4,9,5], [2,4,11], [6,2,10], [8,6,7], [9,8,1],
    ];
    let mut tetrahedra: Vec<[u32;4]> = faces.into_iter()
        .map(|[a,b,c]| [center,a,b,c]).collect();
    if n == 2 {
        // Split only cell interiors: adjacent tets retain the exact same face
        // triangulation and the shell keeps its original 12 boundary vertices.
        for tet in std::mem::take(&mut tetrahedra) {
            let centroid: [f64;3] = std::array::from_fn(|axis|
                tet.iter().map(|&v| positions[v as usize][axis]).sum::<f64>() / 4.0);
            let interior = positions.len() as u32;
            positions.push(centroid);
            for omitted in 0..4 {
                let mut child = [interior;4];
                let mut next = 1;
                for (local, vertex) in tet.into_iter().enumerate() {
                    if local != omitted { child[next] = vertex; next += 1; }
                }
                tetrahedra.push(child);
            }
        }
    }
    for tet in &mut tetrahedra {
        let origin = positions[tet[0] as usize];
        let [a,b,c] = [tet[1],tet[2],tet[3]].map(|v| {
            let point = positions[v as usize];
            [point[0]-origin[0],point[1]-origin[1],point[2]-origin[2]]
        });
        let signed_volume = a[0]*(b[1]*c[2]-b[2]*c[1])
            - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]);
        if signed_volume < 0.0 { tet.swap(2,3); }
    }
    let complex = fs_rep_mesh::TetComplex::from_tets(positions.len(),tetrahedra);
    let fe = HelmholtzFe::new(&complex, &positions, 1, Medium {
        density_kg_m3: atmosphere.density, sound_speed_m_s: atmosphere.sound_speed,
    }).unwrap();
    f(&fe)
}
fn resolution() -> ShellCavityOptions {
    ShellCavityOptions {
        maximum_hz: 6000.0, eigensolver: fs_modal::SliceOptions::default(),
        maximum_pressure_modes: 64, maximum_dofs: 64, maximum_boundary_faces: 64,
        maximum_projection_terms: 10_000_000,
        acoustic_budget: ModalAcousticTimeBudget::audible_reference(),
    }
}
fn closed_body(fe: &HelmholtzFe<'_>, rho: f64, gate: &CancelGate) -> SharedMaterialShellBody {
    let mut triangles = Vec::new();
    for face in fe.boundary_faces() {
        let (p, normal, _) = fe.boundary_geometry(face).unwrap();
        let a = [p[1][0]-p[0][0],p[1][1]-p[0][1],p[1][2]-p[0][2]];
        let b = [p[2][0]-p[0][0],p[2][1]-p[0][1],p[2][2]-p[0][2]];
        let cross = [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
        let mut triangle = face.map(|v| v as usize);
        if cross.iter().zip(normal).map(|(a,b)| a*b).sum::<f64>() < 0.0 { triangle.swap(1,2); }
        triangles.push(triangle);
    }
    // Only boundary nodes belong to the shell; volume-interior FE vertices
    // remain gas degrees of freedom, not fictitious free structural nodes.
    let mut map = std::collections::BTreeMap::new();
    let mut nodes = Vec::new();
    for triangle in &mut triangles {
        for vertex in triangle {
            *vertex = *map.entry(*vertex).or_insert_with(|| {
                nodes.push(fe.positions()[*vertex]); nodes.len()-1
            });
        }
    }
    let vectors = nodes.iter().map(|p| {
        let normal = *p;
        let length = normal.iter().map(|v| v*v).sum::<f64>().sqrt();
        normal.map(|v| v/length)
    }).collect();
    let mesh = ShellMesh::new(nodes, triangles).unwrap();
    let directors = NodalDirectors::new(&mesh, vectors, 0.05).unwrap();
    let moving = mesh.nodes.iter().position(|p| p[0] > 0.0 && p[2] == 0.0).unwrap();
    let constraints = (0..mesh.nodes.len()).map(|node| {
        let mut mask = [true;6]; if node == moving { mask[0] = false; }
        ShellConstraint { node, mask, frame: ConstraintFrame::Global }
    }).collect();
    let budget = MitcBudget { max_nodes: mesh.nodes.len(), max_triangles: mesh.tris.len() };
    let facets: Vec<_> = (0..mesh.tris.len()).map(|triangle|
        OrientedShellFacet { triangle, side: ShellFace::Positive }).collect();
    let input = ShellRegionMaterial {
        region: PlateRegionMaterial {
            region: PlateRegion { name: "skin".into(), triangle_indices: (0..mesh.tris.len()).collect() },
            material: shell_material(rho,1000.0), model: PlateMaterialModel::Orthotropic12,
            thickness_constraint: PlateThicknessConstraint::FixedThickness(0.003), material_angle_rad: 0.17,
        }, shear_correction: 5.0/6.0, axis1_global: [1.0,0.173,0.519],
    };
    let specimen = compile_material_shear_shell(mesh,vec![input],budget,gate).unwrap()
        .with_pressure_surfaces(ShellPressureSurfaces { exterior: facets.clone(), cavity: facets }, gate).unwrap();
    let mut port = vec![0.0;6*specimen.mesh().nodes.len()]; port[6*moving] = 1.0;
    // One source loss law, one force footprint, one scalar modal factor and
    // both complete pressure maps. The 20-element fixture needs 20,007,612
    // terms with the native MITC condensation allowance, not the former 1e6.
    let projection_work = fs_plate::shell::mitc::viscosity::projection_terms(
        budget.max_nodes,budget.max_triangles,1).unwrap()
        + 6*budget.max_nodes + 4
        + 2*fs_plate::shell::pressure::projection_terms(
            budget.max_nodes,budget.max_triangles,1).unwrap();
    SharedMaterialShellBody {
        specimen, directors, constraints, mass: ShellMass::Consistent, drilling_alpha: 1e-3,
        eigenvalue_window: (0.0,1e10), n_modes: 1, force_ports: vec![port],
        bridges: vec![SharedPlateBridge { station_fraction: 0.8, stiffness_n_m: 3000.0,
            damping_n_s_m: 0.1, rest_extension_m: 0.0 }],
        assembly_budget: budget,
        projection_budget: ShellProjectionBudget { max_modes: 1, max_projection_terms: projection_work,
            reference_resistance_n_s_m: 1.0 },
        modal_budget: ModalAcousticTimeBudget::audible_reference(), listener_m: 1.0,
    }
}
fn bank(cx: &Cx<'_>, gate: &CancelGate, fe: &HelmholtzFe<'_>, atmosphere: GasState,
    rho: f64) -> ReciprocalStringBank
{
    from_body(cx,gate,atmosphere,closed_body(fe,rho,gate))
}
fn from_body(cx: &Cx<'_>, gate: &CancelGate, atmosphere: GasState,
    body: SharedMaterialShellBody) -> ReciprocalStringBank
{
    let mut o = options(0); o.network.max_modes = 64; o.network.max_connections = 64;
    o.network.max_setup_terms = 1_048_576;
    ReciprocalStringBank::from_material_shear_shell(cx,vec![member(7800.0).source],
        atmosphere,body,None,law(),o,gate).unwrap()
}
fn uniform_stiffness(bank: &ReciprocalStringBank) -> f64 {
    let count = bank.source().spatial_cavity().unwrap().coupling().cavity_modes();
    let connections = bank.source().body_connections();
    let spring = &connections[connections.len()-count];
    spring.stiffness_n_m * spring.left.shapes[0] * spring.left.shapes[0]
}

#[test]
fn spatial_pressure_reacts_on_the_original_material_body_and_string() {
    with_cx(|cx,gate| with_volume(2,gas(101_325.0),|fe| {
        let mut bare = bank(cx,gate,fe,gas(101_325.0),700.0);
        let mut loaded = bank(cx,gate,fe,gas(101_325.0),700.0)
            .with_material_shell_cavity(fe,resolution(),gate).unwrap();
        let cavity = loaded.source().spatial_cavity().unwrap();
        assert!(cavity.coupling().cavity_modes() > 8);
        assert_eq!(cavity.coupling().total_modes(),
            cavity.coupling().structural_modes()+cavity.coupling().cavity_modes()-1);
        let volume = (10.0/3.0)*(3.0+5.0_f64.sqrt())*0.1*0.11*0.09;
        assert!((cavity.volume_m3()-volume).abs() < 1e-14);
        let expected = loaded.source().material_shell_image().unwrap()
            .cavity_volume_shapes_m2_per_sqrt_kg()[0] / cavity.volume_m3().sqrt();
        assert!((cavity.overlap_rows()[0]-expected).abs() < 1e-10*expected.abs());
        assert_eq!(bare.system().components()[1].modes(),loaded.system().components()[1].modes());
        assert_eq!(bare.source().body_connections()[0],loaded.source().body_connections()[0]);
        let mut scratch = vec![0.0;cavity.observation_state_len()];
        let mut pressure = vec![0.0;cavity.coupling().cavity_modes()];
        let mut force = vec![0.0;loaded.system().mode_count()]; force[4] = 0.02;
        let mut bare_force = vec![0.0;bare.system().mode_count()]; bare_force[4] = 0.02;
        let released = FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };
        let mut sound_difference = 0.0;
        for _ in 0..64 {
            let a = *bare.step(&bare_force,released).unwrap();
            let b = *loaded.step(&force,released).unwrap();
            assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
            sound_difference += (a.observer_pressure_pa-b.observer_pressure_pa).abs();
        }
        loaded.spatial_cavity_pressures_into(&mut scratch,&mut pressure).unwrap();
        assert!(pressure.iter().any(|p| p.abs() > 1e-9));
        assert!(loaded.system().components()[2].states().iter().any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
        assert!(sound_difference > 1e-12);
        assert_ne!(bare.system().components()[0].states(),loaded.system().components()[0].states());
    }));
}

#[test]
fn atmosphere_and_material_substitutions_recompute_the_actual_cavity_load() {
    with_cx(|cx,gate| {
        let air = gas(101_325.0);
        let ordinary = with_volume(1,air,|fe| bank(cx,gate,fe,air,700.0)
            .with_material_shell_cavity(fe,resolution(),gate).unwrap());
        let pressurized = with_volume(1,gas(202_650.0),|fe| bank(cx,gate,fe,gas(202_650.0),700.0)
            .with_material_shell_cavity(fe,resolution(),gate).unwrap());
        let dense = with_volume(1,air,|fe| bank(cx,gate,fe,air,1400.0)
            .with_material_shell_cavity(fe,resolution(),gate).unwrap());
        let a = ordinary.source().spatial_cavity().unwrap();
        let b = pressurized.source().spatial_cavity().unwrap();
        let d = dense.source().spatial_cavity().unwrap();
        assert_eq!(a.basis_identity(),b.basis_identity());
        assert_ne!(a.basis_identity(),d.basis_identity());
        assert!((uniform_stiffness(&pressurized)/uniform_stiffness(&ordinary)-2.0).abs() < 1e-10);
        assert!((uniform_stiffness(&dense)/uniform_stiffness(&ordinary)-0.5).abs() < 1e-8);
        use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
        use fs_qty::semantic::{Composition, CompositionBasis};
        let composition = Composition::new(CompositionBasis::MoleFraction,vec![0.4,0.1,0.5]).unwrap();
        let mixture = resolve_gas_mixture(
            &[GasSpecies::Nitrogen,GasSpecies::Oxygen,GasSpecies::CarbonDioxide],
            &composition,293.15,101_325.0).unwrap();
        let mixed = with_volume(1,mixture,|fe| bank(cx,gate,fe,mixture,700.0)
            .with_material_shell_cavity(fe,resolution(),gate).unwrap());
        let m = mixed.source().spatial_cavity().unwrap();
        assert_eq!(a.basis_identity(),m.basis_identity());
        assert_eq!(*m.gas(),mixture);
        assert_ne!(mixture.sound_speed,air.sound_speed);
        let ratio = m.modal_report().frequency_hz[1]/a.modal_report().frequency_hz[1];
        assert!((ratio-mixture.sound_speed/air.sound_speed).abs() < 1e-7);
        let bulk_ratio = mixture.density*mixture.sound_speed.powi(2)/(air.density*air.sound_speed.powi(2));
        assert!((uniform_stiffness(&mixed)/uniform_stiffness(&ordinary)-bulk_ratio).abs() < 1e-8);
    });
}

#[test]
fn cavity_and_unsteady_string_gas_loading_commute_without_losing_components() {
    with_cx(|cx,gate| with_volume(1,gas(101_325.0),|fe| {
        let flow = CylinderFlowResolution { minimum_wavenumber_radius: 0.01,
            maximum_wavenumber_radius: 4.0, intervals: 2 };
        let mut first = bank(cx,gate,fe,gas(101_325.0),700.0)
            .with_material_shell_cavity(fe,resolution(),gate).unwrap()
            .with_cylinder_gas_loading(flow,gate).unwrap();
        let mut second = bank(cx,gate,fe,gas(101_325.0),700.0)
            .with_cylinder_gas_loading(flow,gate).unwrap()
            .with_material_shell_cavity(fe,resolution(),gate).unwrap();
        assert!(first.source().spatial_cavity().is_some());
        assert!(second.source().cylinder_gas_loading().is_some());
        assert_eq!(first.system().mode_count(),second.system().mode_count());
        assert_eq!(first.system().acceleration_pressure_weights(),second.system().acceleration_pressure_weights());
        assert_eq!(first.source().body_connections(),second.source().body_connections());
        let zero = vec![0.0;first.system().mode_count()];
        for _ in 0..32 { assert_eq!(first.step(&zero,drive()).unwrap(),second.step(&zero,drive()).unwrap()); }
        let cavity = first.source().spatial_cavity().unwrap();
        let mut scratch = vec![0.0;cavity.observation_state_len()];
        let mut a = vec![0.0;cavity.coupling().cavity_modes()]; let mut b = a.clone();
        first.spatial_cavity_pressures_into(&mut scratch,&mut a).unwrap();
        second.spatial_cavity_pressures_into(&mut scratch,&mut b).unwrap();
        assert_eq!(a,b);
    }));
}

#[test]
fn missing_faces_wrong_normals_stale_gas_and_incomplete_pressure_windows_refuse() {
    with_cx(|cx,gate| with_volume(1,gas(101_325.0),|fe| {
        for flip in [false,true] {
            let mut body = closed_body(fe,700.0,gate);
            let mut surfaces = body.specimen.pressure_surfaces().unwrap().clone();
            if flip { surfaces.cavity[0].side = ShellFace::Negative; } else { surfaces.cavity.pop(); }
            body.specimen = body.specimen.with_pressure_surfaces(surfaces,gate).unwrap();
            assert!(from_body(cx,gate,gas(101_325.0),body)
                .with_material_shell_cavity(fe,resolution(),gate).is_err());
        }
        assert!(bank(cx,gate,fe,gas(202_650.0),700.0)
            .with_material_shell_cavity(fe,resolution(),gate).is_err());
        let mut small = resolution(); small.maximum_pressure_modes = 1;
        assert!(bank(cx,gate,fe,gas(101_325.0),700.0)
            .with_material_shell_cavity(fe,small,gate).is_err());
        let loaded = bank(cx,gate,fe,gas(101_325.0),700.0)
            .with_material_shell_cavity(fe,resolution(),gate).unwrap();
        assert!(loaded.with_material_shell_cavity(fe,resolution(),gate).is_err());
        let cold = bank(cx,gate,fe,gas(101_325.0),700.0); gate.request();
        assert!(cold.with_material_shell_cavity(fe,resolution(),gate).is_err());
    }));
}



#[test]
fn prescribed_bridge_pulse_drives_the_same_body_cavity_and_strings_on_one_clock() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::{
        MultiBowPerformanceConfig,
        force::{MaterialForcePort, MaterialForceProgram, StringForceCoupling, StringForcePort},
    };
    use fs_couple::pcm_wav::observation::PressureRenderer;
    use fs_couple::render::schedule::force::{
        ForceInitialization, ForceRenderConfig, ModalForceEvent,
        coupled::contact::friction::multiple::MultiFrictionConfig,
    };
    with_cx(|cx,gate| with_volume(1,gas(101_325.0),|fe| {
        let make = || bank(cx,gate,fe,gas(101_325.0),700.0)
            .with_material_shell_cavity(fe,resolution(),gate).unwrap()
            .with_additional_bows(vec![],MultiFrictionConfig {
                max_ports:4,max_sweeps:32,max_setup_terms:1_048_576,
            },gate).unwrap();
        let mut direct = make();
        let string = StringForcePort {
            location:StringBowLocation {string:0,station_fraction:0.31},
            coupling:StringForceCoupling::Centreline,
        };
        // Separate actuators on the SAME bridge must add coherently; releasing
        // one does not release the other. A string actuator shares their clock.
        let ports = vec![MaterialForcePort::ShellBridge {bridge:0},
            MaterialForcePort::ShellBridge {bridge:0},MaterialForcePort::String(string)];
        let columns = direct.material_force_columns(&ports,1024,gate).unwrap();
        let source = direct.source().bank();
        let body = source.strings().len();
        let offset = direct.system().components()[..body].iter().map(|c|c.modes().len()).sum::<usize>();
        let row = &source.material_shell_image().unwrap().force_ports()[0];
        assert_eq!(&columns[0][offset..offset+row.len()],row.as_slice());
        assert!(columns[0][..offset].iter().all(|&f|f==0.0));
        // The real spatial-cavity constructor appends a separate acoustic
        // component. Its pressure receives no fictitious direct bridge force.
        assert!(direct.system().components().len()>body+1);
        assert!(columns[0][offset+row.len()..].iter().all(|&f|f==0.0));
        assert_eq!(columns[0],columns[1]);
        assert_eq!(columns[2],direct.string_force_columns(&[string],1024,gate).unwrap()[0]);
        assert!(direct.material_force_columns(&[MaterialForcePort::ShellBridge {bridge:1}],1024,gate).is_err());
        assert!(direct.material_force_columns(&ports,direct.system().mode_count()*ports.len()-1,gate).is_err());
        let event = |sample,port,force_n| ModalForceEvent {sample,voice:0,port,force_n};
        let program = MaterialForceProgram {
            ports,initial_force_n:vec![0.02,-0.01,0.0],
            // Unsorted events and equal-time overwrite exercise the original
            // grouped compiler. All three ports release at sample 20.
            events:vec![event(20,0,0.0),event(8,0,0.03),event(8,0,0.02),
                event(8,1,0.005),event(8,2,0.001),event(20,1,0.0),event(20,2,0.0)],
            initialization:ForceInitialization::RetainState,
            limits:ForceRenderConfig {sample_rate_hz:48_000,max_block:32,
                max_events:64,max_controls:1024,max_projection_terms:32_768},
        };
        let released = FrictionDrive {speed_m_s:0.0,normal_force_n:0.0};
        let mut expected = Vec::new();
        for sample in 0..96 {
            let held = if sample<8 {[0.02,-0.01,0.0]} else if sample<20 {[0.02,0.005,0.001]} else {[0.0;3]};
            let force:Vec<_> = (0..direct.system().mode_count()).map(|mode|
                columns.iter().zip(held).fold(0.0,|sum,(column,n)|sum+column[mode]*n)).collect();
            let frame = direct.step(&force,&[released]).unwrap();
            assert!(frame.energy_residual_j.abs()<=frame.energy_tolerance_j);
            if sample>=20 {assert_eq!(frame.external_work_j,0.0);}
            expected.push(frame.observer_pressure_pa.to_bits());
            // Before the string actuator turns on, its motion must already
            // come from reciprocal bridge/body coupling.
            if sample==7 {
                assert!(direct.system().components()[0].states().iter()
                    .any(|s|s.velocity_m_sqrt_kg_per_s!=0.0));
            }
        }
        assert!(expected[20..].iter().any(|&p|f64::from_bits(p).abs()>1e-12));
        assert!(direct.system().components()[0].states().iter().any(|s|s.velocity_m_sqrt_kg_per_s!=0.0));
        assert!(direct.system().components()[body+1].states().iter().any(|s|s.velocity_m_sqrt_kg_per_s!=0.0));
        let cavity = direct.source().bank().spatial_cavity().unwrap();
        let mut scratch = vec![0.0;cavity.observation_state_len()];
        let mut pressure = vec![0.0;cavity.coupling().cavity_modes()];
        direct.spatial_cavity_pressures_into(&mut scratch,&mut pressure).unwrap();
        assert!(pressure.iter().any(|p|p.abs()>1e-12));
        for chunk in [1,31] {
            let schedule = GestureSchedule::try_new(700,vec![GestureTrack {
                id:"released".into(),target:GestureTarget::BowStroke {string:0},
                initial:GestureValue::Bow {velocity_m_per_s:0.0,normal_force_n:0.0,station:0.2},
                events:vec![],
            }]).unwrap();
            let mut p = make().into_scheduled(schedule,vec!["released".into()],
                MultiBowPerformanceConfig {samples:96,max_block:32,max_work:4096,max_events:64}).unwrap()
                .with_material_forces(program.clone(),gate).unwrap();
            assert_eq!(p.material_force_program().unwrap().events,program.events);
            assert!(p.string_force_program().is_none());
            assert_eq!(p.projected_material_forces().unwrap().events(),p.projected_string_forces().unwrap().events());
            let mut output = vec![0.0;96];
            for block in output.chunks_mut(chunk) {p.block(block).unwrap();}
            assert_eq!(p.samples_rendered(),96);
            assert_eq!(p.bank().system().samples_rendered(),96);
            assert_eq!(output.iter().map(|p|p.to_bits()).collect::<Vec<_>>(),expected);
            for (a,b) in p.bank().system().components().iter().zip(direct.system().components()) {
                assert_eq!(a.states(),b.states());
            }
        }
    }));
}
