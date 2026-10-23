//! Actual source shell/pressure FE coupling; this synthetic fixture is not
//! a measured instrument or a claim of converged full-band acoustics.
use super::*;
use fs_couple::air_path::cylinder::CylinderFlowResolution;
use fs_couple::bowed_string::runtime::reciprocal::bank::cavity::ShellCavityOptions;
use fs_couple::shear_shell::pressure::ShellPressureSurfaces;
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace};

fn with_volume<R>(n: usize, atmosphere: GasState, f: impl FnOnce(&HelmholtzFe<'_>) -> R) -> R {
    let (complex, positions) = fs_feec::kuhn_cube(n);
    let positions: Vec<_> = positions.into_iter()
        .map(|p| [0.4*p[0], 0.3*p[1], 0.2*p[2]]).collect();
    let fe = HelmholtzFe::new(&complex, &positions, 1, Medium {
        density_kg_m3: atmosphere.density, sound_speed_m_s: atmosphere.sound_speed,
    }).unwrap();
    f(&fe)
}
fn resolution() -> ShellCavityOptions {
    ShellCavityOptions {
        maximum_hz: 6000.0, eigensolver: fs_modal::SliceOptions::default(),
        maximum_pressure_modes: 32, maximum_dofs: 64, maximum_boundary_faces: 64,
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
        let normal = [p[0]-0.2, p[1]-0.15, p[2]-0.1];
        let length = normal.iter().map(|v| v*v).sum::<f64>().sqrt();
        normal.map(|v| v/length)
    }).collect();
    let mesh = ShellMesh::new(nodes, triangles).unwrap();
    let directors = NodalDirectors::new(&mesh, vectors, 0.05).unwrap();
    let moving = mesh.nodes.iter().position(|p| *p == [0.4,0.3,0.2]).unwrap();
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
    SharedMaterialShellBody {
        specimen, directors, constraints, mass: ShellMass::Consistent, drilling_alpha: 1e-3,
        eigenvalue_window: (0.0,1e10), n_modes: 1, force_ports: vec![port],
        bridges: vec![SharedPlateBridge { station_fraction: 0.8, stiffness_n_m: 3000.0,
            damping_n_s_m: 0.1, rest_extension_m: 0.0 }],
        assembly_budget: budget,
        projection_budget: ShellProjectionBudget { max_modes: 1, max_projection_terms: 1_000_000,
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
        assert!((cavity.volume_m3()-0.4*0.3*0.2).abs() < 1e-14);
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
