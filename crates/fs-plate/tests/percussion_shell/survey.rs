use fs_plate::shell::{ShellMesh, ShellSupport};
use fs_plate::shell::profile::{ProfileBudget, ProfileStation, revolve};
use fs_plate::shell::survey::{IsotropicMaterial, MeshShell};
use fs_plate::shell::reduction::{ReductionBudget, ShellReduction};
use fs_plate::{ModePair, PlateSection};

fn budget() -> ProfileBudget {
    ProfileBudget { max_nodes: 1024, max_triangles: 2048, max_feature_evaluations: 0 }
}
fn bronze() -> IsotropicMaterial {
    IsotropicMaterial { young_pa: 112.6e9, poisson: 0.342, density_kg_m3: 8607.0 }
}
fn panel() -> ShellMesh {
    ShellMesh::new(vec![[0.,0.,0.], [0.1,0.,0.], [0.1,0.1,0.], [0.,0.1,0.]],
        vec![[0,1,2],[0,2,3]]).unwrap()
}
fn shell(mesh: ShellMesh, h: f64) -> MeshShell {
    let n = mesh.nodes.len(); let f = mesh.tris.len();
    MeshShell::new(mesh, vec![h;n], &vec![bronze();f], budget()).unwrap()
}

#[test]
fn nodal_thickness_and_facet_materials_reach_the_original_section_and_mass_laws() {
    let thickness = vec![0.001,0.002,0.004,0.003];
    let materials = [bronze(), IsotropicMaterial { young_pa: 90e9, density_kg_m3: 7800., ..bronze() }];
    let s = MeshShell::new(panel(), thickness.clone(), &materials, budget()).unwrap();
    assert_eq!(s.mesh, panel()); assert_eq!(s.nodal_thickness_m, thickness);
    let mut mass = 0.0;
    for (e,t) in s.mesh.tris.iter().enumerate() {
        let h = (thickness[t[0]]+thickness[t[1]]+thickness[t[2]])/3.0;
        let m = materials[e];
        let expected = PlateSection::isotropic(m.young_pa,m.poisson,h,m.density_kg_m3).unwrap();
        assert_eq!(s.sections[e].d,expected.d);
        mass += 0.005*h*m.density_kg_m3;
    }
    assert!((s.mass_kg-mass).abs() < 1e-15);
    assert!((s.max_edge_m-0.1*2.0_f64.sqrt()).abs() < 1e-15);
    let model = s.assemble(&[],ShellSupport::Free).unwrap();
    let mut vertical = vec![0.;model.free];
    for i in 0..s.mesh.nodes.len() { vertical[model.dof_map[6*i+2].unwrap()] = 1.; }
    let mut momentum = vertical.clone(); model.m.spmv(&vertical,&mut momentum);
    let modal_mass: f64 = vertical.iter().zip(&momentum).map(|(a,b)|a*b).sum();
    assert!((modal_mass-mass).abs() < 1e-15);
    let thin = shell(panel(),0.001); let thick = shell(panel(),0.002);
    let a = thin.assemble(&[],ShellSupport::Free).unwrap();
    let b = thick.assemble(&[],ShellSupport::Free).unwrap();
    assert!((thick.mass_kg/thin.mass_kg-2.).abs() < 1e-14);
    assert!((b.k.get(2,2)/a.k.get(2,2)-8.).abs() < 1e-12);
}

#[test]
fn nonaxisymmetric_reference_geometry_changes_stiffness_and_true_surface_mass() {
    let flat = shell(panel(),0.001);
    let mut mesh = panel(); mesh.nodes[2][2] = 0.009; mesh.nodes[3][0] = -0.003;
    let changed = shell(mesh.clone(),0.001);
    assert_eq!(changed.mesh,mesh);
    assert!(changed.mass_kg > flat.mass_kg);
    let a = flat.assemble(&[],ShellSupport::Free).unwrap();
    let b = changed.assemble(&[],ShellSupport::Free).unwrap();
    assert!((a.k.get(2,2)-b.k.get(2,2)).abs() > 1e-4);
    let mut translated = mesh;
    for p in &mut translated.nodes { p[0] += 0.25; p[1] -= 0.1; p[2] += 0.05; }
    let same = shell(translated,0.001);
    assert!((same.mass_kg-changed.mass_kg).abs() < 1e-14);
}

#[test]
fn a_sampled_profile_preserves_geometry_pencil_nonlinearity_and_two_sided_radiation() {
    let stations = [ProfileStation { radius_m: 0.006, height_m: 0.01, thickness_m: 0.0015 },
        ProfileStation { radius_m: 0.03, height_m: 0.006, thickness_m: 0.001 },
        ProfileStation { radius_m: 0.1, height_m: 0., thickness_m: 0.0005 }];
    let p = revolve(&stations,8,bronze().young_pa,bronze().poisson,bronze().density_kg_m3,
        &[],&[],budget()).unwrap();
    let sampled = MeshShell::new(p.mesh.clone(),p.nodal_thickness_m.clone(),
        &vec![bronze();p.mesh.tris.len()],budget()).unwrap();
    assert_eq!(sampled.mass_kg,p.mass_kg); assert_eq!(sampled.max_edge_m,p.max_edge_m);
    let old = p.assemble(&[],ShellSupport::Free).unwrap();
    let new = sampled.assemble(&[],ShellSupport::Free).unwrap();
    for i in 0..old.free { for j in 0..old.free {
        assert_eq!(old.k.get(i,j),new.k.get(i,j)); assert_eq!(old.m.get(i,j),new.m.get(i,j));
    }}
    // A physical translation observes the complete finite-thickness boundary,
    // including the mounting-hole wall. No pressure or contact oscillator is invented.
    let mut phi = vec![0.;new.free];
    for i in 0..sampled.mesh.nodes.len() { phi[new.dof_map[6*i+2].unwrap()] = 1./sampled.mass_kg.sqrt(); }
    let mut defect = phi.clone(); new.k.spmv(&phi,&mut defect);
    let residual = defect.iter().enumerate().map(|(i,x)|x*x/new.m.get(i,i)).sum::<f64>().sqrt();
    let modes = [ModePair { lambda: 0.,phi,residual,interval:(-residual,residual) }];
    let r = ShellReduction::new(&sampled.mesh,&sampled.sections,&new,&modes,
        ReductionBudget { max_modes: 1,max_facet_modes: 2048,relative_tolerance: 1e-5 }).unwrap();
    let surface = r.radiation_surface(&sampled.nodal_thickness_m,
        fs_plate::shell::reduction::radiation::RadiationSurfaceBudget { max_panels: 1024,max_panel_modes: 2048 }).unwrap();
    assert_eq!(surface.triangles().len(),2*sampled.mesh.tris.len()+32);
    assert_eq!(surface.normal_velocity_weights().len(),1);
    assert!(r.potential(&[0.0001]).abs() < 1e-12);
}

#[test]
fn disconnected_nonmanifold_duplicate_or_unreferenced_geometry_is_not_repaired() {
    let fails = |m: ShellMesh| {
        let (n,f) = (m.nodes.len(),m.tris.len());
        assert!(MeshShell::new(m,vec![0.001;n],&vec![bronze();f],budget()).is_err());
    };
    let mut m=panel(); m.tris[1]=[0,3,2]; fails(m); // shared edge has equal winding
    let mut m=panel(); m.tris.push(m.tris[0]); fails(m);
    let mut m=panel(); m.nodes.push([0.2,0.2,0.]); fails(m);
    let mut m=panel(); m.nodes.push(m.nodes[0]); m.tris[1][0]=4; fails(m);
    let mut m=panel(); m.tris[1][2]=999; fails(m);
    let mut m=panel(); m.nodes.extend([[0.2,0.,0.],[0.3,0.,0.],[0.2,0.1,0.]]); m.tris.push([4,5,6]); fails(m);
    let mut m=panel(); m.nodes.push([0.05,0.05,0.01]); m.tris.push([0,2,4]); fails(m); // three faces on an edge
    let mut m=panel(); m.nodes.extend([[-0.1,0.,0.],[0.,-0.1,0.]]); m.tris.push([0,4,5]); fails(m); // pinched vertex
    let mut m=panel(); m.nodes[0][0]=f64::NAN; fails(m);
}

#[test]
fn field_admission_and_original_mesh_budgets_remain_explicit() {
    for thickness in [vec![0.001;3],vec![0.;4],vec![f64::NAN;4]] {
        assert!(MeshShell::new(panel(),thickness,&[bronze();2],budget()).is_err());
    }
    assert!(MeshShell::new(panel(),vec![0.001;4],&[bronze();1],budget()).is_err());
    for material in [IsotropicMaterial { young_pa: -1.,..bronze() },
        IsotropicMaterial { poisson: 1.,..bronze() },
        IsotropicMaterial { density_kg_m3: 0.,..bronze() }] {
        assert!(MeshShell::new(panel(),vec![0.001;4],&[bronze(),material],budget()).is_err());
    }
    for b in [ProfileBudget { max_nodes: 3,..budget() },ProfileBudget { max_triangles: 1,..budget() }] {
        assert!(MeshShell::new(panel(),vec![0.001;4],&[bronze();2],b).is_err());
    }
}
