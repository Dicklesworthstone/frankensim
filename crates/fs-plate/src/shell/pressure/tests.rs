//! G0/G3 manufactured geometry/work checks, independent of an elastic eigensolver.
use super::*;

fn geometry() -> ShellMesh {
    // Four outward tetrahedron faces plus an independent internal triangle.
    // The last piece can carry mechanics without enclosing this fluid volume.
    ShellMesh::new(vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0],
        [0.1,0.1,0.1],[0.3,0.1,0.1],[0.1,0.3,0.1]],
        vec![[0,2,1],[0,1,3],[0,3,2],[1,2,3],[4,5,6]]).unwrap()
}
fn budget(mesh: &ShellMesh) -> ProfileBudget {
    ProfileBudget { max_nodes:mesh.nodes.len(), max_triangles:mesh.tris.len(),
        max_feature_evaluations:0 }
}
fn facets(count: usize, side: ShellFace) -> Vec<OrientedShellFacet> {
    (0..count).map(|triangle| OrientedShellFacet { triangle, side }).collect()
}
fn mode(mesh: &ShellMesh, displacement: impl Fn([f64;3]) -> [f64;3]) -> Vec<f64> {
    mesh.nodes.iter().flat_map(|&position| {
        let u = displacement(position);
        // Pure midsurface traction has no independent director moment.
        [u[0],u[1],u[2],0.3,-0.7,1.2]
    }).collect()
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-14*a.abs().max(b.abs()).max(1.0), "{a} != {b}");
}

#[test]
fn closed_pressure_work_annihilates_rigid_motion_and_recovers_affine_volume() {
    let mesh = geometry();
    let port = ShellPressurePort::new(&mesh, facets(4,ShellFace::Positive), budget(&mesh), ||false).unwrap();
    let motions = vec![mode(&mesh, |_|[0.7,-0.4,0.6]),
        mode(&mesh, |[x,y,z]|[-0.5*z-0.3*y,0.3*x-0.2*z,0.2*y+0.5*x]),
        mode(&mesh, |position|position)];
    let shapes = port.project_modes(&motions, 10_000, ||false).unwrap();
    close(shapes[0],0.0);
    close(shapes[1],0.0);
    // div(u)=3, tetrahedron V=1/6: delta V=3V=1/2 exactly.
    close(shapes[2],0.5);
    for (motion, volume) in motions.iter().zip(shapes) {
        let mechanical_work: f64 = port.physical_row_m2().iter().zip(motion)
            .map(|(pressure_load,velocity)|2.3*pressure_load*velocity).sum();
        close(mechanical_work,2.3*volume);
    }
    assert!(port.physical_row_m2().chunks_exact(6).all(|node|node[3..]==[0.0;3]));
}

#[test]
fn unselected_internal_structure_has_no_fictitious_fluid_volume() {
    let mesh = geometry();
    let selected = facets(4,ShellFace::Positive);
    let exterior = ShellPressurePort::new(&mesh,selected.clone(),budget(&mesh),||false).unwrap();
    let whole = ShellPressurePort::new(&mesh,facets(5,ShellFace::Positive),budget(&mesh),||false).unwrap();
    let mut internal_motion = vec![0.0;6*mesh.nodes.len()];
    for node in 4..7 { internal_motion[6*node+2]=1.0; }
    let motions = vec![internal_motion];
    assert_eq!(exterior.facets(),selected);
    assert_eq!(exterior.project_modes(&motions,10_000,||false).unwrap(),[0.0]);
    close(whole.project_modes(&motions,10_000,||false).unwrap()[0],0.02);
    assert!(exterior.physical_row_m2()[6*4..].iter().all(|v|*v==0.0));
    close(whole.area_m2()-exterior.area_m2(),0.02);
}

#[test]
fn orientation_reverses_both_conjugate_rows_and_explicit_empty_preserves_coordinates() {
    let mesh = geometry();
    let positive = ShellPressurePort::new(&mesh,facets(4,ShellFace::Positive),budget(&mesh),||false).unwrap();
    let negative = ShellPressurePort::new(&mesh,facets(4,ShellFace::Negative),budget(&mesh),||false).unwrap();
    let empty = ShellPressurePort::new(&mesh,vec![],budget(&mesh),||false).unwrap();
    let motions = vec![mode(&mesh,|p|p),mode(&mesh,|[x,y,z]|[y,z,x])];
    for (&a,&b) in positive.physical_row_m2().iter().zip(negative.physical_row_m2()) { close(a,-b); }
    let plus = positive.project_modes(&motions,10_000,||false).unwrap();
    let minus = negative.project_modes(&motions,10_000,||false).unwrap();
    for (a,b) in plus.into_iter().zip(minus) { close(a,-b); }
    assert_eq!(positive.area_m2(),negative.area_m2());
    assert_eq!(empty.area_m2(),0.0);
    assert_eq!(empty.physical_row_m2().len(),6*mesh.nodes.len());
    assert_eq!(empty.project_modes(&motions,10_000,||false).unwrap(),[0.0,0.0]);
}

#[test]
fn original_indices_unique_ownership_projection_extent_and_cancellation_are_admitted() {
    let mesh = geometry();
    let mut duplicate = facets(1,ShellFace::Positive);
    duplicate.push(OrientedShellFacet { triangle:0, side:ShellFace::Negative });
    assert!(ShellPressurePort::new(&mesh,duplicate,budget(&mesh),||false).is_err());
    assert!(ShellPressurePort::new(&mesh,vec![OrientedShellFacet {
        triangle:mesh.tris.len(),side:ShellFace::Positive }],budget(&mesh),||false).is_err());
    let mut short = budget(&mesh); short.max_nodes-=1;
    assert!(ShellPressurePort::new(&mesh,vec![],short,||false).is_err());
    assert!(matches!(ShellPressurePort::new(&mesh,vec![],budget(&mesh),||true),Err(PlateError::Cancelled)));
    let port = ShellPressurePort::new(&mesh,facets(4,ShellFace::Positive),budget(&mesh),||false).unwrap();
    let modes = vec![mode(&mesh,|p|p)];
    let work = projection_terms(mesh.nodes.len(),4,1).unwrap();
    assert!(port.project_modes(&modes,work-1,||false).is_err());
    close(port.project_modes(&modes,work,||false).unwrap()[0],0.5);
    assert!(port.project_modes(&[vec![0.0;6]],10_000,||false).is_err());
    let mut nonfinite = modes.clone(); nonfinite[0][0]=f64::NAN;
    assert!(port.project_modes(&nonfinite,10_000,||false).is_err());
    assert!(matches!(port.project_modes(&modes,work,||true),Err(PlateError::Cancelled)));
    assert!(projection_terms(usize::MAX,1,1).is_none());
}
