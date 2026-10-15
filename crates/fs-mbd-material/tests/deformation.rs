//! Real source inventories and dynamics under compatible, nonuniform P1 maps.
mod common;
use common::{close, density_card, point};
use fs_ga::{Se3, Vec3, Wrench};
use fs_mbd::articulated::{ArticulatedModel, BaseState, JointModel, Link, SpatialInertia, forward_dynamics};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd_material::{MaterialInertiaError, RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_mbd_material::deformation::TetrahedralMass;

const NODES: [[f64;3];5] = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.]];
const TETS: [[usize;4];2] = [[0,1,2,3],[0,2,1,4]];
fn vector([x,y,z]:[f64;3])->Vec3 { Vec3::new(x,y,z) }
fn reference(rho:[f64;2], tets:&[[usize;4]])->ResolvedBody {
    let cards=[density_card("a",rho[0],rho[0]),density_card("b",rho[1],rho[1])];
    let p=point(300.);
    let keys=["a","b"];
    let bindings:Vec<_>=tets.iter().enumerate().map(|(i,tet)|RegionBinding {
        key:keys[i],material:&cards[i],point:&p,
        selection:MaterialPropertySelection::SingleClaimOnly,
        geometry:VolumeMoments::tetrahedron(tet.map(|v|vector(NODES[v]))).unwrap(),
    }).collect();
    resolve_body(&bindings,2,||false).unwrap()
}
fn acceleration(inertia:SpatialInertia)->f64 {
    let model=ArticulatedModel::new(vec![
        Link::new("base",None,Se3::identity(),JointModel::FIXED,inertia),
        Link::new("body",Some(0),Se3::identity(),JointModel::revolute(Vec3::new(0.,0.,1.),None).unwrap(),inertia),
    ]).unwrap();
    forward_dynamics(&model,BaseState::stationary(Se3::identity()),&[0.],&[0.],&[1.],
        Vec3::default(),&[Wrench::default();2]).unwrap().generalized_acceleration[0]
}
// Independent degree-two tetrahedral quadrature: four equal-weight samples.
// This integrates the entire inertia at the body origin without using either
// VolumeMoments or the production parallel-axis reduction.
fn quadrature(nodes:&[[f64;3]],tets:&[[usize;4]],masses:&[f64])->(Vec3,[f64;9]) {
    let a=(5.0+3.0*5.0_f64.sqrt())/20.0;
    let b=(5.0-5.0_f64.sqrt())/20.0;
    let mut center=Vec3::default(); let mut about_origin=[0.;9];
    for (tet,&mass) in tets.iter().zip(masses) {
        for dominant in 0..4 {
            let mut q=Vec3::default();
            for j in 0..4 { q=q+vector(nodes[tet[j]]).scale(if j==dominant {a}else{b}); }
            center=center+q.scale(mass/4.);
            let [x,y,z]=[q.x,q.y,q.z];
            let row=[y*y+z*z,-x*y,-x*z,-x*y,x*x+z*z,-y*z,-x*z,-y*z,x*x+y*y];
            for (out,v) in about_origin.iter_mut().zip(row) { *out+=mass/4.*v; }
        }
    }
    let mass:f64 = masses.iter().sum();
    (center.scale(1./mass),about_origin)
}

#[test]
fn no_op_and_reset_retain_exact_mass_inventory() {
    let original=reference([6.,12.],&TETS);
    let mesh=TetrahedralMass::bind(&original,&NODES,&TETS,&["a","b"],5,||false).unwrap();
    let unchanged=mesh.deform(&NODES,||false).unwrap();
    assert_eq!(unchanged.inertia(),original.inertia());
    let mut hot=NODES;hot[3][2]=1.1;
    mesh.deform(&hot,||false).unwrap();
    assert_eq!(mesh.deform(&NODES,||false).unwrap(),unchanged);
    assert_eq!(unchanged.reference(),&original);
    assert_eq!(unchanged.tetrahedra(),TETS);
    assert_eq!(mesh.region_elements(),[0,1]);
}

#[test]
fn nonuniform_shear_matches_independent_volume_quadrature_and_real_dynamics() {
    let original=reference([6.,12.],&TETS);
    let mesh=TetrahedralMass::bind(&original,&NODES,&TETS,&["a","b"],5,||false).unwrap();
    let moved=NODES.map(|[x,y,z]|[1.2*x+0.3*y+0.1*z,0.8*y,if z>=0. {1.1*z}else{0.9*z}]);
    let body=mesh.deform(&moved,||false).unwrap();
    let masses:Vec<_>=original.regions().iter().map(|r|r.inertia().mass()).collect();
    let (center,origin)=quadrature(&moved,&TETS,&masses);
    close(center.x,body.inertia().center_of_mass().x);
    close(center.y,body.inertia().center_of_mass().y);
    close(center.z,body.inertia().center_of_mass().z);
    let [x,y,z]=[center.x,center.y,center.z]; let m=body.inertia().mass();
    let shift=[y*y+z*z,-x*y,-x*z,-x*y,x*x+z*z,-y*z,-x*z,-y*z,x*x+y*y];
    for i in 0..9 {close(body.inertia().inertia_com().m[i]+m*shift[i],origin[i]);}
    close(acceleration(body.inertia()),1./origin[8]);
    for (a,b) in original.regions().iter().zip(body.regions()) {
        assert_eq!(a.inertia().mass().to_bits(),b.inertia().mass().to_bits());
        close(b.density_kg_m3()*b.geometry().volume_m3(),a.inertia().mass());
    }
    // One shared face, but different physical volume changes in its two cells.
    close(body.regions()[0].geometry().volume_m3()/original.regions()[0].geometry().volume_m3(),1.2*0.8*1.1);
    close(body.regions()[1].geometry().volume_m3()/original.regions()[1].geometry().volume_m3(),1.2*0.8*0.9);
}

#[test]
fn data_only_density_change_and_cell_order_reach_the_same_solver() {
    let a=reference([6.,12.],&TETS);let b=reference([18.,36.],&TETS);
    let moved=NODES.map(|[x,y,z]|[x+0.4*y,y,z]); // unit-volume shear, not dilation
    let ma=TetrahedralMass::bind(&a,&NODES,&TETS,&["a","b"],5,||false).unwrap();
    let mb=TetrahedralMass::bind(&b,&NODES,&TETS,&["a","b"],5,||false).unwrap();
    let da=ma.deform(&moved,||false).unwrap();let db=mb.deform(&moved,||false).unwrap();
    close(db.inertia().mass()/da.inertia().mass(),3.);
    close(acceleration(db.inertia())/acceleration(da.inertia()),1./3.);
    let reverse=TetrahedralMass::bind(&a,&NODES,&[TETS[1],TETS[0]],&["b","a"],5,||false).unwrap();
    let dr=reverse.deform(&moved,||false).unwrap();
    assert_eq!(dr.inertia(),da.inertia());assert_eq!(dr.regions(),da.regions());
    assert_eq!(reverse.region_elements(),[1,0]);
}

#[test]
fn foreign_incomplete_and_invalid_reference_geometry_refuse() {
    let original=reference([6.,12.],&TETS);
    for keys in [vec!["a"],vec!["a","a"],vec!["a","foreign"]] {
        assert!(TetrahedralMass::bind(&original,&NODES,&TETS,&keys,5,||false).is_err());
    }
    assert!(TetrahedralMass::bind(&original,&NODES,&TETS,&["a","b"],4,||false).is_err());
    for tets in [[TETS[0],TETS[0]],[[0,1,2,8],TETS[1]],[[0,1,1,3],TETS[1]]] {
        assert!(TetrahedralMass::bind(&original,&NODES,&tets,&["a","b"],5,||false).is_err());
    }
    let mut foreign=NODES;foreign[3][2]=1.01;
    assert!(TetrahedralMass::bind(&original,&foreign,&TETS,&["a","b"],5,||false).is_err());
}

#[test]
fn inverted_collapsed_and_cancelled_updates_leave_reference_unchanged() {
    let original=reference([6.,12.],&TETS);let before=original.clone();
    let mesh=TetrahedralMass::bind(&original,&NODES,&TETS,&["a","b"],5,||false).unwrap();
    for z in [-1.,0.,f64::NAN,f64::INFINITY] {
        let mut bad=NODES;bad[3][2]=z;assert!(mesh.deform(&bad,||false).is_err());
    }
    assert!(mesh.deform(&NODES[..4],||false).is_err());
    let mut moved=NODES;moved[3][2]=1.01;
    let mut total=0;let success=mesh.deform(&moved,||{total+=1;false}).unwrap();
    for fail in [1,total/2,total] {
        let mut n=0;
        assert!(matches!(mesh.deform(&moved,||{n+=1;n==fail}),Err(MaterialInertiaError::Cancelled)));
    }
    assert_eq!(original,before);assert_eq!(mesh.deform(&moved,||false).unwrap(),success);
}

#[test]
fn negative_reference_orientation_does_not_create_negative_mass() {
    let tets=[[0,2,1,3],TETS[1]];
    let original=reference([6.,12.],&tets);
    let mesh=TetrahedralMass::bind(&original,&NODES,&tets,&["a","b"],5,||false).unwrap();
    let moved=NODES.map(|[x,y,z]|[-y,x,z]); // proper 90-degree rotation
    let body=mesh.deform(&moved,||false).unwrap();
    assert_eq!(body.inertia().mass(),original.inertia().mass());
    close(body.inertia().inertia_com().m[0],original.inertia().inertia_com().m[4]);
    close(body.inertia().inertia_com().m[8],original.inertia().inertia_com().m[8]);
}
