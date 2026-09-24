use super::*;
use super::super::{TensionedDisk,TensionedDiskSpec};
use crate::{AssemblyOptions,EdgeSupport,PlateSection,assemble};
use crate::shell::profile::ProfileBudget;

fn spec() -> TensionedDiskSpec {
    TensionedDiskSpec {radius_m:0.1,thickness_m:0.0002,young_pa:4e9,poisson:0.38,
        density_kg_m3:1390.0,tension_n_m:1000.0,radial_intervals:3,azimuths:16}
}
fn budget() -> ProfileBudget {
    ProfileBudget {max_nodes:200,max_triangles:400,max_feature_evaluations:0}
}
fn field() -> TensionVariation {
    TensionVariation {constant_n_m:[250.0,-150.0,80.0],gradient_n_m2:[600.0,-200.0,100.0,400.0]}
}
fn close(a:f64,b:f64,relative:f64) {
    assert!((a-b).abs()<=relative*a.abs().max(b.abs()).max(1.0),"{a:e} != {b:e}");
}

#[test]
fn zero_field_preserves_uniform_pencil_and_nonzero_isotropic_field_matches_owner() {
    let original=TensionedDisk::new(spec(),budget()).unwrap();
    let zero=TensionedDisk::new_with_tension_variation(spec(),TensionVariation::default(),budget()).unwrap();
    assert_eq!(original.model.k,zero.model.k);assert_eq!(original.model.m,zero.model.m);
    let boundary:Vec<_>=(0..original.mesh.nodes.len()).filter(|&n|original.model.dof_map[3*n].is_none()).collect();
    let old=assemble(&original.mesh,&original.section,&boundary,&[],&AssemblyOptions {
        pretension:spec().tension_n_m,support:EdgeSupport::SimplySupported}).unwrap();
    assert_eq!(original.model.k,old.k);assert_eq!(original.model.m,old.m);
    let uniform=TensionVariation {constant_n_m:[200.0,200.0,0.0],..TensionVariation::default()};
    let varied=TensionedDisk::new_with_tension_variation(spec(),uniform,budget()).unwrap();
    let direct=TensionedDisk::new(TensionedDiskSpec {tension_n_m:1200.0,..spec()},budget()).unwrap();
    assert_eq!(varied.model.m,direct.model.m);assert_eq!(varied.mass_kg,direct.mass_kg);
    assert_eq!(varied.model.dof_map,direct.model.dof_map);
    for r in 0..varied.model.free {for c in 0..varied.model.free {
        close(varied.model.k.get(r,c),direct.model.k.get(r,c),2e-12);
    }}
}

#[test]
fn affine_tensor_integrates_exact_transverse_virtual_work_and_keeps_slope_dofs() {
    let mesh=PlateMesh::from_unstructured(vec![(0.0,0.0),(0.08,0.01),(-0.02,0.07)],vec![[0,1,2]]).unwrap();
    let section=PlateSection::isotropic(4e9,0.38,0.0002,1390.0).unwrap();
    let mut model=assemble(&mesh,&section,&[],&[],&AssemblyOptions {pretension:0.0,support:EdgeSupport::SimplySupported}).unwrap();
    let original=model.clone();let variation=field();variation.add_to(1000.0,&mesh,&mut model).unwrap();
    let [gx,gy]=[0.013,-0.009];
    let q:Vec<_>=mesh.nodes.iter().flat_map(|&(x,y)|[gx*x+gy*y,gx,gy]).collect();
    let mut work=0.0;
    for i in 0..9 {for j in 0..9 {
        let delta=model.k.get(i,j)-original.k.get(i,j);
        close(delta,model.k.get(j,i)-original.k.get(j,i),1e-12);
        if i%3!=0 || j%3!=0 {assert_eq!(delta,0.0);}
        work+=q[i]*delta*q[j];
    }}
    let area=0.5*(0.08*0.07+0.02*0.01);
    let [xx,yy,xy]=variation.resultant(1000.0,[0.06/3.0,0.08/3.0]);
    close(work,area*(xx*gx*gx+2.0*xy*gx*gy+yy*gy*gy),1e-12);
    assert_eq!(model.m,original.m);
}

#[test]
fn nonuniform_field_has_balanced_boundary_tractions_without_interior_forces() {
    let p=[[-0.08,-0.04],[0.09,-0.02],[0.01,0.095]];
    let variation=field();let mut force=[0.0;2];let mut moment=0.0;
    for i in 0..3 {
        let a=p[i];let b=p[(i+1)%3];let dx=b[0]-a[0];let dy=b[1]-a[1];
        for t in [0.5-0.5/3.0_f64.sqrt(),0.5+0.5/3.0_f64.sqrt()] {
            let x=[a[0]+t*dx,a[1]+t*dy];let [xx,yy,xy]=variation.resultant(1000.0,x);
            // Outward normal times edge length for the counterclockwise rim.
            let f=[0.5*(xx*dy-xy*dx),0.5*(xy*dy-yy*dx)];
            force[0]+=f[0];force[1]+=f[1];moment+=x[0]*f[1]-x[1]*f[0];
        }
    }
    assert!(force.iter().all(|f|f.abs()<1e-11));assert!(moment.abs()<1e-12);
    assert_ne!(variation.resultant(1000.0,p[0]),variation.resultant(1000.0,p[1]));
}

#[test]
fn anisotropic_tension_splits_a_real_head_mode_pair_without_changing_mass_or_supports() {
    let plain=TensionedDisk::new(spec(),budget()).unwrap();
    let tuned=TensionedDisk::new_with_tension_variation(spec(),
        TensionVariation {constant_n_m:[350.0,-350.0,100.0],..TensionVariation::default()},budget()).unwrap();
    let modes=|d:&TensionedDisk|fs_modal::slice_window(&d.model.k,&d.model.m,
        (0.0,(std::f64::consts::TAU*700.0).powi(2)),&crate::SliceOptions::default())
        .unwrap().modes.iter().map(|m|m.lambda).collect::<Vec<_>>();
    let a=modes(&plain);let b=modes(&tuned);
    assert!(a.len()>=3 && b.len()>=3);
    close(a[1],a[2],1e-5);
    assert!((b[2]-b[1]).abs()>0.01*b[1]);
    assert_eq!(plain.model.m,tuned.model.m);assert_eq!(plain.model.dof_map,tuned.model.dof_map);
    assert_eq!(plain.mesh.nodes,tuned.mesh.nodes);assert_eq!(plain.mass_kg,tuned.mass_kg);
}

#[test]
fn loss_of_tension_at_a_vertex_refuses_even_when_the_center_is_tensile() {
    let mut bad_fields=vec![TensionVariation {constant_n_m:[-1001.0,0.0,0.0],..TensionVariation::default()},
        TensionVariation {constant_n_m:[0.0,0.0,1000.0],..TensionVariation::default()},
        TensionVariation {gradient_n_m2:[20000.0,0.0,0.0,0.0],..TensionVariation::default()}];
    for i in 0..7 {
        let mut bad=TensionVariation::default();
        if i<3 {bad.constant_n_m[i]=f64::NAN;} else {bad.gradient_n_m2[i-3]=f64::INFINITY;}
        bad_fields.push(bad);
    }
    for bad in bad_fields {assert!(TensionedDisk::new_with_tension_variation(spec(),bad,budget()).is_err());}
    // Ordinary spatial variation remains legal; it is not replaced by its mean.
    let disk=TensionedDisk::new_with_tension_variation(spec(),field(),budget()).unwrap();
    assert_eq!(disk.tension_variation,field());
}
