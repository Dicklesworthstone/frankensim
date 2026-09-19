//! Geometry and kinematic gates for real curved, tapered percussion bodies.
use fs_plate::{PlateSection, dkt_stiffness};
use fs_plate::shell::{ShellMesh, ShellSupport, assemble_shell, assemble_shell_sections};
use fs_plate::shell::profile::{ProfileStation, ProfileBudget, SurfaceIndentation, revolve};
fn energy(k:&fs_sparse::Csr,x:&[f64])->f64 {
    let mut y=vec![0.0;x.len()];k.spmv(x,&mut y);x.iter().zip(y).map(|(a,b)|a*b).sum()
}
fn section(h:f64)->PlateSection {PlateSection::isotropic(110e9,0.34,h,8800.0).unwrap()}
fn rotate(p:[f64;3])->[f64;3] { // fixed proper orthogonal frame
    [(p[0]+2.0*p[1]+2.0*p[2])/3.0,(2.0*p[0]+p[1]-2.0*p[2])/3.0,(-2.0*p[0]+2.0*p[1]-p[2])/3.0]
}
#[test]
fn shell_has_all_six_rigid_motions_on_curved_and_rotated_facets() {
    let nodes=vec![[0.0,0.0,0.0],[0.13,0.0,0.01],[0.02,0.11,0.03],[0.14,0.12,-0.01]];
    for nodes in [nodes.clone(),nodes.into_iter().map(rotate).collect()] {
        let mesh=ShellMesh::new(nodes,vec![[0,1,2],[1,3,2]]).unwrap();
        let model=assemble_shell(&mesh,&section(0.001),&[],ShellSupport::Free).unwrap();
        for c in 0..6 {
            let mut x=vec![0.0;model.free];
            for (i,&[px,py,pz]) in mesh.nodes.iter().enumerate() {
                if c<3 {x[6*i+c]=1.0;} else {
                    let displacement=match c {3=>[0.0,-pz,py],4=>[pz,0.0,-px],_=>[-py,px,0.0]};
                    x[6*i..6*i+3].copy_from_slice(&displacement);x[6*i+c]=1.0;
                }
            }
            let mut y=vec![0.0;x.len()];model.k.spmv(&x,&mut y);
            for row in 0..x.len() {
                let (cols,values)=model.k.row(row);
                let scale=cols.iter().zip(values).map(|(&j,a)|a.abs()*x[j].abs()).sum::<f64>();
                assert!(y[row].abs()<=2e-12*scale.max(1.0),"rigid motion {c}, row {row}: {}",y[row]);
            }
        }
    }
}
#[test]
fn physical_rotations_reproduce_the_existing_plate_curvature_energy() {
    let mesh=ShellMesh::new(vec![[0.0,0.0,0.0],[0.12,0.0,0.0],[0.03,0.09,0.0]],vec![[0,1,2]]).unwrap();
    let s=section(0.0015);let model=assemble_shell(&mesh,&s,&[],ShellSupport::Free).unwrap();
    let curvature=[0.7,-0.3,0.2];let mut physical=vec![0.0;18];let mut plate=[0.0;9];
    for (i,&[x,y,_]) in mesh.nodes.iter().enumerate() {
        let w=0.5*curvature[0]*x*x+0.5*curvature[1]*y*y+0.5*curvature[2]*x*y;
        let wx=curvature[0]*x+0.5*curvature[2]*y;
        let wy=curvature[1]*y+0.5*curvature[2]*x;
        physical[6*i+2]=w;physical[6*i+3]=wy;physical[6*i+4]=-wx;
        plate[3*i]=w;plate[3*i+1]=wx;plate[3*i+2]=wy;
    }
    let (k,_)=dkt_stiffness(&[0.0,0.12,0.03],&[0.0,0.0,0.09],&s.d,0).unwrap();
    let expected=(0..9).map(|i|(0..9).map(|j|plate[i]*k[9*i+j]*plate[j]).sum::<f64>()).sum::<f64>();
    assert!((energy(&model.k,&physical)-expected).abs()<1e-12*expected);
    let mut rotated=mesh.clone();rotated.nodes=rotated.nodes.into_iter().map(rotate).collect();
    let rotated_model=assemble_shell(&rotated,&s,&[],ShellSupport::Free).unwrap();
    let mut v=physical.clone();for i in 0..3 {for offset in [0,3] {
        let j=6*i+offset;v[j..j+3].copy_from_slice(&rotate([physical[j],physical[j+1],physical[j+2]]));
    }}
    assert!((energy(&rotated_model.k,&v)-expected).abs()<1e-8*expected);
}
#[test]
fn element_thickness_changes_actual_bending_and_inertia_without_a_mass_floor() {
    let mesh=ShellMesh::new(vec![[0.0,0.0,0.0],[0.1,0.0,0.0],[0.0,0.1,0.0]],vec![[0,1,2]]).unwrap();
    let a=assemble_shell(&mesh,&section(0.001),&[],ShellSupport::Free).unwrap();
    let b=assemble_shell_sections(&mesh,&[section(0.002)],&[],ShellSupport::Free).unwrap();
    assert!((b.m.get(0,0)/a.m.get(0,0)-2.0).abs()<1e-12);
    assert!((b.m.get(3,3)/a.m.get(3,3)-8.0).abs()<1e-12);
    assert!((b.k.get(2,2)/a.k.get(2,2)-8.0).abs()<1e-12);
    assert!(assemble_shell_sections(&mesh,&[],&[],ShellSupport::Free).is_err());
    let mut invalid=section(0.001);invalid.d[0]=f64::NAN;
    assert!(assemble_shell(&mesh,&invalid,&[],ShellSupport::Free).is_err());
}
fn budget()->ProfileBudget {ProfileBudget{max_nodes:10000,max_triangles:20000,max_feature_evaluations:100000}}
#[test]
fn revolved_annulus_and_unique_center_preserve_real_area_mass_and_taper() {
    for inner in [0.0,0.006] {
        let stations=[ProfileStation{radius_m:inner,height_m:0.0,thickness_m:0.001},
            ProfileStation{radius_m:0.254,height_m:0.0,thickness_m:0.001}];
        let mesh=revolve(&stations,64,110e9,0.34,8800.0,&[],&[],budget()).unwrap();
        let area=0.5*64.0*(2.0*std::f64::consts::PI/64.0).sin()*(0.254_f64.powi(2)-inner*inner);
        assert!((mesh.mass_kg-8800.0*0.001*area).abs()<1e-12);
        assert_eq!(mesh.inner_nodes.len(),if inner==0.0 {1}else{64});
        assert_eq!(mesh.mesh.nodes.len(),if inner==0.0 {65}else{128});
    }
    let stations=[ProfileStation{radius_m:0.006,height_m:0.03,thickness_m:0.002},
        ProfileStation{radius_m:0.06,height_m:0.02,thickness_m:0.0015},
        ProfileStation{radius_m:0.254,height_m:0.0,thickness_m:0.0008}];
    let clean=revolve(&stations,32,110e9,0.34,8800.0,&[],&[],budget()).unwrap();
    let dent=SurfaceIndentation{center_m:[0.06,0.0],radius_m:0.02,height_delta_m:-0.002,thickness_delta_m:-0.0001};
    let hammered=revolve(&stations,32,110e9,0.34,8800.0,&[],&[dent],budget()).unwrap();
    assert_ne!(clean.mesh.nodes,hammered.mesh.nodes);
    assert_ne!(clean.mass_kg.to_bits(),hammered.mass_kg.to_bits());
    let a=clean.assemble(&clean.inner_nodes,ShellSupport::Clamped).unwrap();
    let b=hammered.assemble(&hammered.inner_nodes,ShellSupport::Clamped).unwrap();
    assert_ne!(a.k,b.k);assert_ne!(a.m,b.m);
    assert!(hammered.underresolved_features>0);
}
#[test]
fn profile_budgets_and_nonpositive_feature_thickness_refuse() {
    let s=[ProfileStation{radius_m:0.0,height_m:0.0,thickness_m:0.001},
        ProfileStation{radius_m:0.1,height_m:0.0,thickness_m:0.001}];
    let mut small=budget();small.max_nodes=4;
    assert!(revolve(&s,32,110e9,0.34,8800.0,&[],&[],small).is_err());
    let dent=SurfaceIndentation{center_m:[0.0,0.0],radius_m:0.01,height_delta_m:0.0,thickness_delta_m:-0.002};
    assert!(revolve(&s,32,110e9,0.34,8800.0,&[],&[dent],budget()).is_err());
}
