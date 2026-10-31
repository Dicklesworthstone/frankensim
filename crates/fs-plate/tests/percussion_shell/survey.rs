use fs_plate::shell::{CstDktMass, ShellMesh, ShellSupport};
use fs_plate::shell::profile::{ProfileBudget, ProfileStation, revolve};
use fs_plate::shell::survey::{IsotropicMaterial, SurveyMaterial, MeshShell};
use fs_plate::shell::reduction::{ReductionBudget, ShellReduction};
use fs_plate::{ModePair, PlateSection};

fn budget() -> ProfileBudget {
    ProfileBudget { max_nodes: 1024, max_triangles: 2048, max_feature_evaluations: 0 }
}
fn bronze() -> IsotropicMaterial {
    IsotropicMaterial { young_pa: 112.6e9, poisson: 0.342, density_kg_m3: 8607.0 }
}
/// Enum variants have no functional-update syntax: edit a copy of an orthotropic fixture.
fn edited(mut m: SurveyMaterial, edit: impl FnOnce(&mut f64, &mut f64, &mut [f64; 3])) -> SurveyMaterial {
    let SurveyMaterial::Orthotropic { e2_pa, g12_pa, axis1_global, .. } = &mut m else {
        panic!("orthotropic fixture expected")
    };
    edit(e2_pa, g12_pa, axis1_global);
    m
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
    // Profile assembly retains its original pencil; survey defaults now use
    // consistent inertia, so request the legacy comparison explicitly.
    let new = sampled.assemble_with_mass(&[],ShellSupport::Free,CstDktMass::Lumped).unwrap();
    for i in 0..old.free { for j in 0..old.free {
        assert_eq!(old.k.get(i,j),new.k.get(i,j)); assert_eq!(old.m.get(i,j),new.m.get(i,j));
    }}
    // A physical translation observes the complete finite-thickness boundary,
    // including the mounting-hole wall. No pressure or contact oscillator is invented.
    let mut phi = vec![0.;new.free];
    for i in 0..sampled.mesh.nodes.len() { phi[new.dof_map[6*i+2].unwrap()] = 1./sampled.mass_kg.sqrt(); }
    let mut defect = phi.clone(); new.k.spmv(&phi,&mut defect);
    let residual = defect.iter().enumerate().map(|(i,x)|x*x/new.m.get(i,i)).sum::<f64>().sqrt();
    // A reduction needs an elastic spectrum to scale its eigenpair check; the
    // translation rides beside the lowest computed elastic modes.
    let tau = core::f64::consts::TAU;
    let elastic = fs_plate::shell::modes_shell(&new,((tau*20.0).powi(2),(tau*20_000.0).powi(2)),
        &fs_plate::SliceOptions::default()).unwrap().modes;
    assert!(!elastic.is_empty());
    let mut modes = vec![ModePair { lambda: 0.,phi,residual,interval:(-residual,residual) }];
    modes.extend(elastic.into_iter().take(3));
    let r = ShellReduction::new(&sampled.mesh,&sampled.sections,&new,&modes,
        ReductionBudget { max_modes: modes.len(),max_facet_modes: 2048,relative_tolerance: 1e-5 }).unwrap();
    let surface = r.radiation_surface(&sampled.nodal_thickness_m,
        fs_plate::shell::reduction::radiation::RadiationSurfaceBudget { max_panels: 1024,max_panel_modes: 2048 }).unwrap();
    assert_eq!(surface.triangles().len(),2*sampled.mesh.tris.len()+32);
    assert_eq!(surface.normal_velocity_weights().len(),modes.len());
    let mut q = vec![0.;modes.len()]; q[0] = 0.0001;
    assert!(r.potential(&q).abs() < 1e-12);
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


#[test]
fn global_orthotropic_axis_reaches_each_facet_without_mesh_diagonal_grain_jumps() {
    let mesh=panel(); let h=vec![0.001;mesh.nodes.len()];
    let material=SurveyMaterial::Orthotropic {
        e1_pa:120e9,e2_pa:70e9,nu12:0.27,g12_pa:30e9,density_kg_m3:8200.,
        axis1_global:[1.,0.,0.],
    };
    let anisotropic=MeshShell::new_with_materials(mesh.clone(),h.clone(),&[material;2],budget()).unwrap();
    let isotropic=shell(mesh.clone(),0.001);
    assert!((anisotropic.mass_kg-0.01*0.001*8200.).abs()<1e-14);
    assert_ne!(anisotropic.sections[0].d,isotropic.sections[0].d);
    // The two triangles have different local first edges. Their stored local D
    // matrices must therefore differ so that the SAME global x material axis is
    // represented continuously rather than restarting at each triangle.
    assert_ne!(anisotropic.sections[0].d,anisotropic.sections[1].d);
    let a=anisotropic.assemble(&[],ShellSupport::Free).unwrap();
    let b=isotropic.assemble(&[],ShellSupport::Free).unwrap();
    assert_ne!(a.k.get(0,0),b.k.get(0,0));

    // Rotating both geometry and the declared physical axis is a rigid change
    // of frame: mass is identical and the material remains admissible.
    let mut rotated=mesh;
    for p in &mut rotated.nodes { let x=p[0];p[0]=-p[1];p[1]=x; }
    let rotated_material=edited(material, |_, _, axis| *axis=[0.,1.,0.]);
    let same=MeshShell::new_with_materials(rotated,h,&[rotated_material;2],budget()).unwrap();
    assert!((same.mass_kg-anisotropic.mass_kg).abs()<1e-14);
}

#[test]
fn orthotropic_survey_refuses_missing_or_normal_only_material_axes() {
    let base=SurveyMaterial::Orthotropic {
        e1_pa:120e9,e2_pa:70e9,nu12:0.27,g12_pa:30e9,density_kg_m3:8200.,
        axis1_global:[1.,0.,0.],
    };
    for material in [
        edited(base, |_, _, axis| *axis=[0.,0.,0.]),
        edited(base, |_, _, axis| *axis=[0.,0.,1.]),
        edited(base, |e2, _, _| *e2=0.0),
        edited(base, |_, g12, _| *g12=f64::NAN),
    ] {
        assert!(MeshShell::new_with_materials(panel(),vec![0.001;4],&[material;2],budget()).is_err());
    }
}

fn inertia_cross(a: [f64;3], b: [f64;3]) -> [f64;3] {
    [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
}
fn inertia_norm2(x: [f64;3]) -> f64 { x.into_iter().map(|v|v*v).sum() }
fn inertia_triangle() -> ShellMesh {
    ShellMesh::new(vec![[0.,0.,0.],[0.3,0.,0.],[0.07,0.19,0.]],vec![[0,1,2]]).unwrap()
}

#[test]
fn consistent_survey_mass_matches_exact_triangular_prism_rigid_inertia() {
    let h: f64 = 0.01;
    for turned in [false,true] {
        let mut mesh = inertia_triangle();
        if turned { for p in &mut mesh.nodes {
            *p=super::rotate(*p);
            for (x,shift) in p.iter_mut().zip([0.4,-0.2,0.7]) { *x+=shift; }
        }}
        let s = shell(mesh,h);
        let model = s.assemble(&[],ShellSupport::Free).unwrap();
        let center: [f64;3] = std::array::from_fn(|c|s.mesh.nodes.iter().map(|p|p[c]).sum::<f64>()/3.);
        let normal = s.mesh.facet(0).unwrap().frame[2];
        for (v,w) in [([0.2,-0.1,0.3],[0.4,0.5,-0.2]),
            ([0.;3],[1.,0.,0.]),([0.;3],[0.,1.,0.]),([0.;3],[0.,0.,1.]),
            ([0.2,-0.1,0.3],[0.;3])] {
            let mut q=vec![0.;18];
            for (node,p) in s.mesh.nodes.iter().enumerate() {
                let orbital=inertia_cross(w,*p);
                for c in 0..3 {q[6*node+c]=v[c]+orbital[c];q[6*node+3+c]=w[c];}
            }
            let wc=inertia_cross(w,center);
            let vc=std::array::from_fn(|c|v[c]+wc[c]);
            // Uniform triangular-prism central second moment, independently
            // integrated in volume: sum delta_i delta_i^T/12 + h² nn^T/12.
            let orbital:f64=s.mesh.nodes.iter().map(|p|
                inertia_norm2(inertia_cross(w,std::array::from_fn(|c|p[c]-center[c])))).sum();
            let exact=s.mass_kg*(inertia_norm2(vc)+orbital/12.
                +h*h/12.*inertia_norm2(inertia_cross(w,normal)));
            let got=super::energy(&model.m,&q);
            assert!((got/exact-1.).abs()<2e-13,"rigid kinetic energy {got} != {exact}");
        }
        // The old nodal quadrature has four times the centerline orbital
        // inertia, plus its fictitious section inertia about the facet normal.
        let mut spin=vec![0.;18];
        for (node,p) in s.mesh.nodes.iter().enumerate() {
            let v=inertia_cross(normal,std::array::from_fn(|c|p[c]-center[c]));
            spin[6*node..6*node+3].copy_from_slice(&v);
            spin[6*node+3..6*node+6].copy_from_slice(&normal);
        }
        let old=s.assemble_with_mass(&[],ShellSupport::Free,CstDktMass::Lumped).unwrap();
        let orbital=super::energy(&model.m,&spin);
        let old_orbital=super::energy(&old.m,&spin)-s.mass_kg*h*h/12.;
        assert!((old_orbital/orbital-4.).abs()<3e-13);
        assert_eq!(model.k,old.k);
    }
}

fn inertia_affine_state(mesh: &ShellMesh) -> Vec<f64> {
    mesh.nodes.iter().flat_map(|&[x,y,_]|[
        0.2+0.3*x-0.4*y,-0.1+0.5*x+0.2*y,0.4-0.3*x+0.6*y,
        0.3+0.7*x,-0.2-0.6*y,0.4+0.2*x-0.5*y,
    ]).collect()
}
fn inertia_refine(mesh: ShellMesh) -> ShellMesh {
    let mut nodes=mesh.nodes;let mut edges=std::collections::BTreeMap::new();let mut tris=Vec::new();
    for [a,b,c] in mesh.tris {
        let mut midpoint=|i:usize,j:usize| *edges.entry((i.min(j),i.max(j))).or_insert_with(|| {
            let p=std::array::from_fn(|c|0.5*(nodes[i][c]+nodes[j][c]));nodes.push(p);nodes.len()-1
        });
        let (ab,bc,ca)=(midpoint(a,b),midpoint(b,c),midpoint(c,a));
        tris.extend([[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]]);
    }
    ShellMesh::new(nodes,tris).unwrap()
}

#[test]
fn affine_kinetic_integral_is_frame_covariant_and_subdivision_invariant() {
    let h: f64=0.011;let mut mesh=inertia_triangle();
    let q=inertia_affine_state(&mesh);let area=0.3*0.19/2.;let rho=bronze().density_kg_m3;
    // Integral of an affine scalar squared, evaluated from its three vertex
    // values. The independent physical law is |v|² + h²|theta x n|²/12.
    let square=|values:[f64;3]|area/12.*(values.iter().map(|v|v*v).sum::<f64>()
        +values.iter().sum::<f64>().powi(2));
    let physical=rho*h*((0..3).map(|c|square([q[c],q[6+c],q[12+c]])).sum::<f64>()
        +h*h/12.*(3..5).map(|c|square([q[c],q[6+c],q[12+c]])).sum::<f64>());
    // Membrane spin is (v_x-u_y)/2 = 0.45. Only its difference from
    // drilling rotation enters the explicitly numerical kinetic gauge.
    let gauge=rho*h*h*h/12.*square([q[5]-0.45,q[11]-0.45,q[17]-0.45]);
    for _ in 0..3 {
        let source=shell(mesh.clone(),h);let model=source.assemble(&[],ShellSupport::Free).unwrap();
        let q=inertia_affine_state(&mesh);
        assert!((super::energy(&model.m,&q)/(physical+gauge)-1.).abs()<2e-13);
        let rotated=ShellMesh::new(mesh.nodes.iter().copied().map(super::rotate).collect(),mesh.tris.clone()).unwrap();
        let other=shell(rotated,h).assemble(&[],ShellSupport::Free).unwrap();
        let rotated_q:Vec<_>=q.chunks_exact(3).flat_map(|p|super::rotate([p[0],p[1],p[2]])).collect();
        assert!((super::energy(&other.m,&rotated_q)/(physical+gauge)-1.).abs()<2e-13);
        let supported=source.assemble(&[0],ShellSupport::Clamped).unwrap();
        let factor=fs_sparse::SymbolicLdlt::analyze(&supported.m,fs_sparse::DirectOrdering::Amd).unwrap()
            .factor(&supported.m,&fs_sparse::LdltOptions::default()).unwrap();
        assert_eq!(factor.inertia().positive,supported.free);
        assert_eq!(factor.inertia().negative,0);
        for i in 0..model.free {for j in 0..model.free {
            assert_eq!(model.m.get(i,j).to_bits(),model.m.get(j,i).to_bits());
        }}
        mesh=inertia_refine(mesh);
    }
}
