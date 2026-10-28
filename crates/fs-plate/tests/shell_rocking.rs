//! G0/G1 algebraic checks for the explicitly selected first-order rigid pose.
//! These are not finite-rotation or calibrated-cymbal accuracy claims.
use fs_plate::{ModePair,PlateSection};
use fs_plate::shell::{ShellMesh,ShellModel,ShellSupport,assemble_shell_sections};
use fs_plate::shell::reduction::{ShellReduction,ReductionBudget,RigidMotionLimit};
use fs_plate::shell::reduction::rigid::{mass_normalized_modes,mass_normalized_modes_with_options};
use fs_plate::shell::reduction::radiation::{ShellFace,RadiationSurfaceBudget};

const MOTIONS:[[f64;6];3]=[
    [0.,0.,1.,0.,0.,0.],[0.,0.,0.,1.,0.,0.],[0.,0.,0.,0.,1.,0.],
];
fn budget()->ReductionBudget {
    ReductionBudget{max_modes:32,max_facet_modes:100,relative_tolerance:1e-5}
}
fn fixture()->(ShellMesh,Vec<PlateSection>,Vec<f64>,ShellModel,Vec<ModePair>,ShellReduction) {
    // Curved, asymmetric and nonuniform: inertia is not that of an ideal disk.
    let mesh=ShellMesh::new(vec![[0.01,-0.02,0.003],[0.11,0.,0.013],
        [0.,0.10,0.022],[0.08,0.09,0.036]],vec![[0,1,3],[0,3,2]]).unwrap();
    let thickness=vec![0.001,0.0012,0.0008,0.0009];
    let sections:Vec<_>=mesh.tris.iter().enumerate().map(|(i,t)|
        PlateSection::isotropic(112.6e9,0.342,t.iter().map(|n|thickness[*n]/3.).sum(),
            if i==0{8607.}else{7800.}).unwrap()).collect();
    let model=assemble_shell_sections(&mesh,&sections,&[],ShellSupport::Free).unwrap();
    let n=model.free;let mut k=vec![0.;n*n];let mut m=k.clone();
    for i in 0..n {for j in 0..n {k[i*n+j]=model.k.get(i,j);m[i*n+j]=model.m.get(i,j);}}
    let elastic=fs_modal::eigh_gen_dense(&k,&m,n).unwrap().into_iter()
        .filter(|mode|mode.lambda>1e3).collect::<Vec<_>>();
    assert!(!elastic.is_empty());
    let mut modes=mass_normalized_modes(&mesh,&model,&MOTIONS).unwrap();modes.extend(elastic);
    let reduction=ShellReduction::new_with_linearized_rigid_motion(&mesh,&sections,&model,
        &modes,budget(),RigidMotionLimit{modes:3,maximum_rotation_rad:0.05}).unwrap();
    (mesh,sections,thickness,model,modes,reduction)
}
fn physical(motion:[f64;6],point:[f64;3])->[f64;3] {
    [motion[0]+motion[4]*point[2]-motion[5]*point[1],
        motion[1]+motion[5]*point[0]-motion[3]*point[2],
        motion[2]+motion[3]*point[1]-motion[4]*point[0]]
}
fn coordinates(mesh:&ShellMesh,model:&ShellModel,modes:&[ModePair],motion:[f64;6])->Vec<f64> {
    let mut nodal=vec![0.;model.free];
    for (node,p) in mesh.nodes.iter().enumerate() {
        nodal[6*node..6*node+3].copy_from_slice(&physical(motion,*p));
        nodal[6*node+3..6*node+6].copy_from_slice(&motion[3..]);
    }
    let mut mass_velocity=nodal.clone();model.m.spmv(&nodal,&mut mass_velocity);
    let mut q=vec![0.;modes.len()];
    for i in 0..3 {q[i]=modes[i].phi.iter().zip(&mass_velocity).map(|(a,b)|a*b).sum();}
    for i in 0..nodal.len() {
        let recovered=(0..3).map(|j|modes[j].phi[i]*q[j]).sum::<f64>();
        assert!((recovered-nodal[i]).abs()<2e-14*nodal[i].abs().max(1.));
    }
    q
}

#[test]
fn rigid_inertia_matches_independent_facet_mass_and_rotary_work() {
    let (mesh,sections,_,model,modes,r)=fixture();
    let motion=[0.,0.,0.3,-0.2,0.4,0.];
    let velocity=coordinates(&mesh,&model,&modes,motion);
    let mut expected=0.;let mut cross_inertia=0.;
    for (f,t) in mesh.tris.iter().enumerate() {
        let section=&sections[f];
        let mass=section.density*section.thickness*mesh.facet(f).unwrap().area_m2/3.;
        for node in t {
            let v=physical(motion,mesh.nodes[*node]);
            expected+=0.5*mass*(v.iter().map(|v|v*v).sum::<f64>()
                + section.thickness.powi(2)/12.*motion[3..].iter().map(|w|w*w).sum::<f64>());
            cross_inertia+=mass*mesh.nodes[*node][1];
        }
    }
    assert!(cross_inertia.abs()>1e-5,"heave/tilt inertia is coupled on this asymmetric mesh");
    let kinetic=0.5*velocity.iter().map(|v|v*v).sum::<f64>();
    assert!((kinetic-expected).abs()<2e-13*expected);
    let q=coordinates(&mesh,&model,&modes,[0.,0.,0.02,0.02,-0.03,0.]);
    let angle=r.rigid_rotation(&q).unwrap().unwrap();
    assert!((angle[0]-0.02).abs()<1e-14);assert!((angle[1]+0.03).abs()<1e-14);
    assert_eq!(angle[2],0.);
}

#[test]
fn rigid_pose_has_no_quartic_spring_and_does_not_change_elastic_energy_or_loss() {
    let (mesh,_,_,model,modes,r)=fixture();let n=r.mode_count();
    let pose=coordinates(&mesh,&model,&modes,[0.,0.,0.02,0.02,-0.03,0.]);
    assert_eq!(r.potential(&pose),0.);let mut gradient=vec![0.;n];
    r.gradient(&pose,&mut gradient);assert!(gradient.iter().all(|x|*x==0.));
    let mut q=vec![0.;n];q[3]=1e-6;let energy=r.potential(&q);
    r.gradient(&q,&mut gradient);assert!(energy>0.);
    let mut moved=q.clone();for i in 0..3 {moved[i]=pose[i];}
    assert_eq!(r.potential(&moved),energy);
    let mut actual=vec![0.;n];r.gradient(&moved,&mut actual);assert_eq!(actual,gradient);
    r.hessian_vector(&moved,&pose,&mut actual);assert!(actual.iter().all(|x|*x==0.));
    let bending=r.bending_stiffness(&[1.,1.],2*n*n).unwrap();
    for i in 0..3 {for j in 0..n {assert_eq!(bending[i*n+j],0.);assert_eq!(bending[j*n+i],0.);}}
    assert!((3..n).any(|i|bending[i*n+i]>0.));
}

#[test]
fn tilted_skin_contact_and_radiation_keep_actual_offset_levers_and_reciprocity() {
    let (mesh,_,thickness,model,modes,r)=fixture();let n=r.mode_count();
    let motion=[0.,0.,0.3,-0.2,0.4,0.];let velocity=coordinates(&mesh,&model,&modes,motion);
    let normal=[0.6,0.,0.8];let force=2.4;
    for face in [ShellFace::Positive,ShellFace::Negative] {
        let p=r.surface_point_port(&thickness,0,[0.2,0.3,0.5],face,normal).unwrap();
        let expected=physical(motion,p.position_m).iter().zip(normal).map(|(a,b)|a*b).sum::<f64>();
        let modal=p.weights.iter().zip(&velocity).map(|(a,b)|a*b).sum::<f64>();
        assert!((modal-expected).abs()<2e-14);
        let work=velocity.iter().zip(&p.weights).map(|(v,b)|v*b*force).sum::<f64>();
        assert!((work-force*expected).abs()<2e-14);
    }
    let skin=r.radiation_surface(&thickness,RadiationSurfaceBudget{
        max_panels:100,max_panel_modes:100*n}).unwrap();
    let mut vn=vec![0.;skin.triangles().len()];skin.normal_velocity(&velocity,&mut vn).unwrap();
    let pressure:Vec<_>=(0..vn.len()).map(|i|1.+i as f64).collect();let mut load=vec![0.;n];
    skin.pressure_force(&pressure,&mut load).unwrap();
    let power=load.iter().zip(&velocity).map(|(a,b)|a*b).sum::<f64>();
    let boundary=pressure.iter().zip(&vn).zip(skin.areas_m2()).map(|((p,v),a)|p*v*a).sum::<f64>();
    assert!((power+boundary).abs()<1e-13);
}

#[test]
fn supplied_angle_cap_and_geometric_nullity_are_real_refusals() {
    let (mesh,sections,_,model,modes,r)=fixture();
    let outside=coordinates(&mesh,&model,&modes,[0.,0.,0.,0.04,0.04,0.]);
    assert!(r.rigid_rotation(&outside).is_err(),"combined angle exceeds cap although each component fits");
    assert!(r.rigid_rotation(&[]).is_err());
    let mut nonfinite=vec![0.;modes.len()];nonfinite[3]=f64::NAN;
    assert!(r.rigid_rotation(&nonfinite).is_err());
    let limit=RigidMotionLimit{modes:3,maximum_rotation_rad:0.05};
    let mut forged=modes.clone();forged[1].phi[0]+=0.01;
    assert!(ShellReduction::new_with_linearized_rigid_motion(&mesh,&sections,&model,&forged,budget(),limit).is_err());
    let mut grounded=model.clone();let mut k=fs_sparse::Coo::new(model.free,model.free);
    for row in 0..model.free {let (columns,values)=model.k.row(row);
        for (&column,&value) in columns.iter().zip(values) {k.push(row,column,value);}}
    k.push(2,2,1e6);grounded.k=k.assemble();
    assert!(ShellReduction::new_with_linearized_rigid_motion(&mesh,&sections,&grounded,&modes,budget(),limit).is_err());
    for angle in [0.,-0.01,0.10001,f64::NAN] {
        assert!(ShellReduction::new_with_linearized_rigid_motion(&mesh,&sections,&model,&modes,budget(),
            RigidMotionLimit{modes:3,maximum_rotation_rad:angle}).is_err());
    }
    assert!(mass_normalized_modes(&mesh,&model,&[MOTIONS[0],MOTIONS[0]]).is_err());
}

const ALL_MOTIONS:[[f64;6];6]=[
    [1.,0.,0.,0.,0.,0.],[0.,1.,0.,0.,0.,0.],[0.,0.,1.,0.,0.,0.],
    [0.,0.,0.,1.,0.,0.],[0.,0.,0.,0.,1.,0.],[0.,0.,0.,0.,0.,1.],
];
const CONSISTENT_SECTIONS:[(f64,f64);2]=[(8607.,0.0011),(7800.,0.0009)];

fn consistent_fixture()->(ShellMesh,ShellModel) {
    let mesh=ShellMesh::new(vec![[0.01,-0.02,0.003],[0.11,0.,0.013],
        [0.,0.10,0.022],[0.08,0.09,0.036]],vec![[0,1,3],[0,3,2]]).unwrap();
    let n=6*mesh.nodes.len();let mut mass=fs_sparse::Coo::new(n,n);
    // Exact P1 consistent inertia, including physical through-thickness rotary
    // inertia. This manufactured free pencil isolates basis construction from
    // stiffness assembly; the oracle below integrates affine velocity directly.
    for (f,tri) in mesh.tris.iter().enumerate() {
        let (density,h)=CONSISTENT_SECTIONS[f];let area=mesh.facet(f).unwrap().area_m2;
        for (i,&a) in tri.iter().enumerate() {for (j,&b) in tri.iter().enumerate() {
            let integral=area/12.*if i==j {2.} else {1.};
            for c in 0..6 {
                let inertia=if c<3 {density*h} else {density*h.powi(3)/12.};
                mass.push(6*a+c,6*b+c,integral*inertia);
            }
        }}
    }
    (mesh,ShellModel{k:fs_sparse::Coo::new(n,n).assemble(),m:mass.assemble(),
        dof_map:(0..n).map(Some).collect(),free:n})
}
fn rigid_nodal(mesh:&ShellMesh,motion:[f64;6])->Vec<f64> {
    let mut nodal=vec![0.;6*mesh.nodes.len()];
    for (node,&point) in mesh.nodes.iter().enumerate() {
        nodal[6*node..6*node+3].copy_from_slice(&physical(motion,point));
        nodal[6*node+3..6*node+6].copy_from_slice(&motion[3..]);
    }
    nodal
}

#[test]
fn consistent_shell_inertia_preserves_rigid_continuum_work_and_cross_terms() {
    let (mesh,model)=consistent_fixture();
    let modes=mass_normalized_modes(&mesh,&model,&ALL_MOTIONS).unwrap();
    let mut work=vec![0.;model.free];
    for (i,mode) in modes.iter().enumerate() {
        model.m.spmv(&mode.phi,&mut work);
        for (j,other) in modes.iter().enumerate() {
            let product=other.phi.iter().zip(&work).map(|(a,b)|a*b).sum::<f64>();
            assert!((product-if i==j {1.} else {0.}).abs()<2e-13);
        }
        assert_eq!(mode.residual,0.);
    }
    let motion=[0.03,-0.07,0.3,-0.2,0.4,0.15];let nodal=rigid_nodal(&mesh,motion);
    model.m.spmv(&nodal,&mut work);
    let velocity:Vec<_>=modes.iter().map(|mode|
        mode.phi.iter().zip(&work).map(|(a,b)|a*b).sum::<f64>()).collect();
    for (i,&expected) in nodal.iter().enumerate() {
        let recovered=modes.iter().zip(&velocity).map(|(mode,q)|mode.phi[i]*q).sum::<f64>();
        assert!((recovered-expected).abs()<2e-14);
    }
    let mut expected=0.;
    for (f,tri) in mesh.tris.iter().enumerate() {
        let (density,h)=CONSISTENT_SECTIONS[f];let area=mesh.facet(f).unwrap().area_m2;
        let velocities=tri.map(|node|physical(motion,mesh.nodes[node]));
        let mut sum=[0.;3];let mut squared=0.;
        for v in velocities {for c in 0..3 {sum[c]+=v[c];squared+=v[c]*v[c];}}
        expected+=0.5*density*h*area/12.*(squared+sum.iter().map(|x|x*x).sum::<f64>())
            +0.5*density*h.powi(3)/12.*area*motion[3..].iter().map(|x|x*x).sum::<f64>();
    }
    let kinetic=0.5*velocity.iter().map(|v|v*v).sum::<f64>();
    assert!((kinetic-expected).abs()<2e-13*expected);
    let lumped=0.5*nodal.iter().enumerate().map(|(row,v)|
        model.m.row(row).1.iter().sum::<f64>()*v*v).sum::<f64>();
    assert!((lumped-expected).abs()>1e-4*expected,"off-diagonal inertia is load-bearing");
    let translation=rigid_nodal(&mesh,[0.,0.,1.,0.,0.,0.]);
    let rotation=rigid_nodal(&mesh,[0.,0.,0.,1.,0.,0.]);
    model.m.spmv(&rotation,&mut work);
    let cross=translation.iter().zip(&work).map(|(a,b)|a*b).sum::<f64>();
    assert!(cross.abs()>1e-5,"asymmetric geometry retains translation/rotation cross inertia");
}

#[test]
fn consistent_shell_rigid_residual_uses_full_inverse_mass_and_signed_stiffness() {
    let (mesh,mut model)=consistent_fixture();
    // K = s M has exact generalized eigenvalue s, including a negative s.
    // A geometric zero-frequency candidate must report |s| as its residual;
    // using diag(M)^-1 or silently dropping the stiffness cannot pass.
    for scale in [-7.,7.] {
        let mut k=fs_sparse::Coo::new(model.free,model.free);
        for row in 0..model.free {let (columns,values)=model.m.row(row);
            for (&column,&value) in columns.iter().zip(values) {k.push(row,column,scale*value);}}
        model.k=k.assemble();
        let modes=mass_normalized_modes(&mesh,&model,&ALL_MOTIONS).unwrap();
        for mode in &modes {
            assert_eq!(mode.lambda,0.);
            assert!((mode.residual-scale.abs()).abs()<2e-12);
            assert_eq!(mode.interval,(-mode.residual,mode.residual));
        }
        let mut residual=vec![0.;model.free];model.k.spmv(&modes[0].phi,&mut residual);
        let diagonal_norm=residual.iter().enumerate().map(|(i,r)|r*r/model.m.get(i,i)).sum::<f64>().sqrt();
        assert!((diagonal_norm-scale.abs()).abs()>0.1*scale.abs());
    }
}

#[test]
fn consistent_shell_rigid_admission_refuses_nonphysical_mass_and_respects_limits() {
    let (mesh,model)=consistent_fixture();
    let options=fs_sparse::LdltOptions{max_factor_bytes:Some(0),..Default::default()};
    assert!(matches!(mass_normalized_modes_with_options(&mesh,&model,&ALL_MOTIONS,&options,||false),
        Err(fs_plate::PlateError::Modal(fs_modal::ModalError::Factor{
            source:fs_sparse::LdltError::FactorBudget{..},..}))));
    let mut diagonal=model.clone();let mut mass=fs_sparse::Coo::new(model.free,model.free);
    for row in 0..model.free {mass.push(row,row,model.m.row(row).1.iter().sum());}
    diagonal.m=mass.assemble();
    assert!(mass_normalized_modes_with_options(&mesh,&diagonal,&ALL_MOTIONS,&options,||false).is_ok());

    let options=fs_sparse::LdltOptions::default();let mut polls=0;
    mass_normalized_modes_with_options(&mesh,&model,&ALL_MOTIONS,&options,||{polls+=1;false}).unwrap();
    for stop in 1..=polls {
        let mut current=0;
        assert!(matches!(mass_normalized_modes_with_options(&mesh,&model,&ALL_MOTIONS,&options,
            ||{current+=1;current==stop}),Err(fs_plate::PlateError::Cancelled)));
    }
    for (forward,reverse) in [(2.,2.),(1.,1.),(0.25,0.5),(f64::NAN,f64::NAN)] {
        let mut mass=fs_sparse::Coo::new(model.free,model.free);
        for i in 0..model.free {mass.push(i,i,1.);}
        mass.push(0,6,forward);mass.push(6,0,reverse);
        let invalid=ShellModel{m:mass.assemble(),..model.clone()};
        assert!(mass_normalized_modes(&mesh,&invalid,&ALL_MOTIONS).is_err(),
            "indefinite, singular, asymmetric and nonfinite inertia must refuse");
    }
}
