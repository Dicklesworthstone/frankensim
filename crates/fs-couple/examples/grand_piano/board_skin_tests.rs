use super::*;
fn flat() -> ShellMesh {
    ShellMesh::new(vec![[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],[0.5,0.5,0.]],
        vec![[0,1,4],[1,2,4],[2,3,4],[3,0,4]]).unwrap()
}
fn normals(skin:&Skin)->Vec<[f64;3]> {
    skin.panel_triangles().iter().map(|&[a,b,c]| {
        let n=cross(sub(b,a),sub(c,a));let r=dot(n,n).sqrt();n.map(|v|v/r)
    }).collect()
}
fn rigid(mesh:ShellMesh)->MotionSurface {
    let theta=[0.2,-0.3,0.1];let shift=[0.7,-0.8,0.9];
    let mode=mesh.nodes.iter().map(|&p| {
        let u=cross(theta,p);[u[0]+shift[0],u[1]+shift[1],u[2]+shift[2],theta[0],theta[1],theta[2]]
    }).collect();MotionSurface::new(mesh,vec![mode]).unwrap()
}
fn check_rigid(skin:&Skin,motion:&MotionSurface) {
    let n=normals(skin);let w=skin.normal_weights(motion,&n).unwrap();
    let mut flux=0.;
    for (i,t) in skin.panel_triangles().iter().enumerate() {
        let center=std::array::from_fn(|c|t.iter().map(|p|p[c]).sum::<f64>()/3.);
        let rotation=cross([0.2,-0.3,0.1],center);
        let expected=dot(n[i],[0.7+rotation[0],-0.8+rotation[1],0.9+rotation[2]]);
        assert!((w[0][i]-expected).abs()<1e-12);
        let area=0.5*dot(cross(sub(t[1],t[0]),sub(t[2],t[0])),cross(sub(t[1],t[0]),sub(t[2],t[0]))).sqrt();
        flux+=area*w[0][i];
    }
    assert!(flux.abs()<1e-12,"closed rigid skin must not create volume: {flux:e}");
}
#[test]
fn finite_thickness_closes_both_faces_and_exactly_integrates_rigid_motion() {
    let mesh=flat();let skin=Skin::build(&mesh,&[0.008;4],0.02,1000).unwrap();
    assert!((skin.volume_m3-0.008).abs()<1e-14);
    assert!(skin.maximum_thickness_change_m<1e-15);assert_eq!(skin.maximum_offset_m,0.004);
    assert_eq!(skin.triangles.len(),24); // 8 top/bottom, 4 fan triangles per rim edge
    check_rigid(&skin,&rigid(mesh));
    assert_eq!(skin.obj(),skin.obj());assert!(skin.obj().contains("o soundboard_skin"));
}
#[test]
fn section_jumps_keep_their_walls_and_split_third_facet_levels_without_cracks() {
    let mesh=flat();let skin=Skin::build(&mesh,&[0.004,0.006,0.010,0.008],0.02,1000).unwrap();
    assert!((skin.volume_m3-0.007).abs()<1e-14);
    assert!(skin.triangles.len()>24);
    check_rigid(&skin,&rigid(mesh.clone()));
    let smooth=Skin::build(&mesh,&[0.007;4],0.02,1000).unwrap();
    assert_ne!(skin.vertices,smooth.vertices,"do not average supplied sections into a uniform slab");
    // Two diagonally high columns meet only along a vertical line: such a
    // nonmanifold step is refused rather than emitted as an acoustic boundary.
    assert!(Skin::build(&mesh,&[0.01,0.004,0.01,0.004],0.02,1000).is_err());
}
#[test]
fn supplied_crown_uses_normal_section_thickness_and_true_offset_rotation() {
    let mut mesh=flat();for p in &mut mesh.nodes {p[2]=0.04*p[0]+0.02*p[1];}
    let mesh=ShellMesh::new(mesh.nodes,mesh.tris).unwrap();
    let area=(1.+0.04_f64.powi(2)+0.02_f64.powi(2)).sqrt();
    let skin=Skin::build(&mesh,&[0.008;4],0.02,1000).unwrap();
    assert!((skin.section_volume_m3-area*0.008).abs()<1e-14);
    assert!((skin.volume_m3-skin.section_volume_m3).abs()<=area*HEIGHT_QUANTUM_M);
    check_rigid(&skin,&rigid(mesh));
}
#[test]
fn holes_stay_open_and_source_assignment_and_resource_failures_do_not_invent_geometry() {
    let mesh=ShellMesh::new(vec![[0.,0.,0.],[2.,0.,0.],[2.,2.,0.],[0.,2.,0.],
        [0.5,0.5,0.],[1.5,0.5,0.],[1.5,1.5,0.],[0.5,1.5,0.]],
        vec![[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]]).unwrap();
    let skin=Skin::build(&mesh,&[0.008;8],0.02,1000).unwrap();
    assert!((skin.volume_m3-3.*0.008).abs()<1e-14);
    assert_eq!(skin.triangles.len(),48); // eight outer+inner edges
    check_rigid(&skin,&rigid(mesh));
    let mesh=flat();
    for h in [vec![],vec![f64::NAN;4],vec![-1.;4],vec![0.;4],vec![0.1;4]] {
        assert!(Skin::build(&mesh,&h,0.02,1000).is_err());
    }
    assert!(Skin::build(&mesh,&[0.008;4],0.02,23).is_err());
    let motion=rigid(mesh);
    let mut source=String::from("frankensim-board-geometry-si-v1\n");
    for (e,t) in motion.mesh.tris.iter().enumerate() {
        writeln!(source,"triangle,{e},{},{},{},0.008,450,1e10,8e8,0.3,6e8,0",t[0],t[1],t[2]).unwrap();
    }
    assert!(Skin::from_source(&motion,&source,0.02,1000).is_ok());
    for bad in [source.replace("triangle,0,0,1,4","triangle,0,1,0,4"),
        source.replace("triangle,0,","triangle,99,"),source.replace("triangle,0,","triangle,1,"),
        source.replace("0.008","NaN"),String::from("bad header")] {
        assert!(Skin::from_source(&motion,&bad,0.02,1000).is_err());
    }
    let skin=Skin::from_source(&motion,&source,0.02,1000).unwrap();
    let mut moved=motion; moved.mesh.nodes[0][2]+=0.001;
    assert!(skin.normal_weights(&moved,&normals(&skin)).is_err());
}
