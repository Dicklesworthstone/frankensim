//! G0/G1 algebraic checks for the explicitly selected first-order rigid pose.
//! These are not finite-rotation or calibrated-cymbal accuracy claims.
use fs_plate::{ModePair,PlateSection};
use fs_plate::shell::{ShellMesh,ShellModel,ShellSupport,assemble_shell_sections};
use fs_plate::shell::reduction::{ShellReduction,ReductionBudget,RigidMotionLimit};
use fs_plate::shell::reduction::rigid::mass_normalized_modes;
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
