//! Source data -> clipped displacement -> forces, torque and actual mechanics.
mod common;
#[path = "hydrostatic/motion.rs"]
mod motion;
#[path = "hydrostatic/layered.rs"]
mod layered;
#[path = "hydrostatic/finite_tank.rs"]
mod finite_tank;
use common::{close, density_card, point};
use fs_ga::{Se3, So3, So3Tangent, Twist, Vec3};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::{Gravity, Pose, RigidBodyIntegrator, RigidBodyState};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_mbd_material::hydrostatic::{DisplacementMesh, HydrostaticError, HydrostaticMedium};

fn fluid(rho: f64) -> FluidTransportState {
    FluidTransportState::try_new(rho,4000.0,0.5,0.001).unwrap()
}
fn gravity() -> Gravity { Gravity::new(fs_mbd::Vec3::new(0.0,0.0,-10.0)).unwrap() }
fn origin() -> Vec3 { Vec3::new(0.0,0.0,0.0) }
fn pose(z: f64) -> Se3 { Se3::from_parts(So3::identity(),Vec3::new(0.0,0.0,z)).unwrap() }
fn body(density: f64, sides: Vec3) -> ResolvedBody {
    let card=density_card("synthetic-body",density,density);
    resolve_body(&[RegionBinding { key:"solid",geometry:VolumeMoments::cuboid(sides).unwrap(),
        material:&card,point:&point(300.0),selection:MaterialPropertySelection::SingleClaimOnly }],1,||false).unwrap()
}
fn cube() -> DisplacementMesh { DisplacementMesh::cuboid(Vec3::new(2.0,2.0,2.0),||false).unwrap() }
fn vec_close(a: Vec3,b: Vec3) { close(a.x,b.x); close(a.y,b.y); close(a.z,b.z); }

#[test]
fn dry_partial_full_and_waterline_contacts_have_analytic_volume_and_potential() {
    let water=fluid(1000.0); let medium=HydrostaticMedium::new(&water,gravity(),origin()).unwrap();
    let source=body(500.0,Vec3::new(2.0,2.0,2.0)); let body=source.principal_body(1e-12,||false).unwrap();
    let mesh=cube(); close(mesh.volume_m3(),8.0);
    for (z,volume,center_z,depth) in [(2.0,0.0,0.0,0.0),(1.0,0.0,0.0,0.0),
        (0.0,4.0,-0.5,0.5),(-1.0,8.0,-1.0,1.0),(-2.0,8.0,-2.0,2.0)] {
        let state=body.state_from_kinematics(pose(z),Twist::zero()).unwrap();
        let load=medium.load(&mesh,&body,state,||false).unwrap();
        close(load.displaced().volume_m3(),volume);
        close(load.wrench().force_world.z,10000.0*volume);
        close(load.pressure_potential_j(),10000.0*volume*depth);
        close(load.wrench().torque_body.norm_squared(),0.0);
        if volume==0.0 { assert_eq!(load.center_of_buoyancy_world(),None); }
        else { vec_close(load.center_of_buoyancy_world().unwrap(),Vec3::new(0.0,0.0,center_z)); }
    }
}

#[test]
fn all_tetrahedron_orders_recover_independent_sliced_simplex_moments() {
    let water=fluid(1000.0);
    let v=[origin(),Vec3::new(1.0,0.0,0.0),Vec3::new(0.0,1.0,0.0),Vec3::new(0.0,0.0,1.0)];
    for h in [0.25_f64,0.5,0.75] {
        let medium=HydrostaticMedium::new(&water,gravity(),Vec3::new(0.0,0.0,h)).unwrap();
        let volume=(1.0-(1.0-h).powi(3))/6.0;
        let x=(1.0-(1.0-h).powi(4))/(24.0*volume);
        let z=(h*h/4.0-h.powi(3)/3.0+h.powi(4)/8.0)/volume;
        for a in 0..4 { for b in 0..4 { for c in 0..4 { for d in 0..4 {
            if a==b || a==c || a==d || b==c || b==d || c==d { continue; }
            let mesh=DisplacementMesh::new(vec![[v[a],v[b],v[c],v[d]]],1,||false).unwrap();
            let wet=medium.displacement(&mesh,Se3::identity(),||false).unwrap();
            close(wet.volume_m3(),volume); vec_close(wet.centroid_body_m().unwrap(),Vec3::new(x,x,z));
        }}}}
    }
}

#[test]
fn two_wet_vertex_cut_and_oblique_gravity_are_not_wet_vertex_fractions() {
    let water=fluid(1000.0); let medium=HydrostaticMedium::new(&water,gravity(),origin()).unwrap();
    let mesh=DisplacementMesh::new(vec![[Vec3::new(-1.0,-1.0,-1.0),Vec3::new(1.0,1.0,-1.0),
        Vec3::new(1.0,-1.0,1.0),Vec3::new(-1.0,1.0,1.0)]],1,||false).unwrap();
    let wet=medium.displacement(&mesh,Se3::identity(),||false).unwrap();
    close(wet.volume_m3(),4.0/3.0); vec_close(wet.centroid_body_m().unwrap(),Vec3::new(0.0,0.0,-3.0/8.0));
    let oblique=Gravity::new(fs_mbd::Vec3::new(-10.0,-10.0,0.0)).unwrap();
    let medium=HydrostaticMedium::new(&water,oblique,origin()).unwrap();
    let wet=medium.displacement(&cube(),Se3::identity(),||false).unwrap();
    close(wet.volume_m3(),4.0); vec_close(wet.centroid_body_m().unwrap(),Vec3::new(-1.0/3.0,-1.0/3.0,0.0));
}

#[test]
fn body_density_and_fluid_density_change_actual_acceleration_without_solver_changes() {
    let mesh=cube();
    for (solid_rho,fluid_rho,acceleration) in [(500.0,1000.0,0.0),(1000.0,1000.0,-5.0),(1000.0,2000.0,0.0)] {
        let fluid=fluid(fluid_rho); let medium=HydrostaticMedium::new(&fluid,gravity(),origin()).unwrap();
        let source=body(solid_rho,Vec3::new(2.0,2.0,2.0)); let body=source.principal_body(1e-12,||false).unwrap();
        let state=body.state_from_kinematics(Se3::identity(),Twist::zero()).unwrap();
        let load=medium.load(&mesh,&body,state,||false).unwrap();
        let dt=0.001;
        let next=RigidBodyIntegrator::new(medium.gravity()).step(state,body.properties(),load.wrench(),dt).unwrap().state_after;
        close(next.center_of_mass_velocity_world(body.properties()).unwrap().z/dt,acceleration);
    }
}

#[test]
fn heterogeneous_mass_centroid_changes_torque_not_displaced_volume() {
    let left=density_card("light",500.0,500.0); let right=density_card("heavy",1500.0,1500.0);
    let half=VolumeMoments::cuboid(Vec3::new(1.0,2.0,2.0)).unwrap();
    let p=point(300.0);
    let source=resolve_body(&[
        RegionBinding { key:"left",material:&left,point:&p,selection:MaterialPropertySelection::SingleClaimOnly,
            geometry:half.transformed(Se3::from_parts(So3::identity(),Vec3::new(-0.5,0.0,0.0)).unwrap()).unwrap() },
        RegionBinding { key:"right",material:&right,point:&p,selection:MaterialPropertySelection::SingleClaimOnly,
            geometry:half.transformed(Se3::from_parts(So3::identity(),Vec3::new(0.5,0.0,0.0)).unwrap()).unwrap() },
    ],2,||false).unwrap();
    close(source.inertia().center_of_mass().x,0.25);
    let body=source.principal_body(1e-12,||false).unwrap();
    let state=body.state_from_kinematics(Se3::identity(),Twist::zero()).unwrap();
    let water=fluid(1000.0); let medium=HydrostaticMedium::new(&water,gravity(),origin()).unwrap();
    let load=medium.load(&cube(),&body,state,||false).unwrap();
    close(load.displaced().volume_m3(),4.0);
    vec_close(load.center_of_buoyancy_world().unwrap(),Vec3::new(0.0,0.0,-0.5));
    let torque=state.pose().orientation().rotate_body_to_world(load.wrench().torque_body);
    close(torque.x,0.0); close(torque.y,10000.0); close(torque.z,0.0);
}

#[test]
fn sealed_displacement_envelope_is_not_inferred_from_solid_material_volume() {
    let water=fluid(1000.0); let medium=HydrostaticMedium::new(&water,gravity(),origin()).unwrap();
    // Mass geometry is an explicitly idealized compact ballast; sealed hull
    // displacement is separately authored. No density rescaling fakes a hull.
    let source=body(1000.0,Vec3::new(1.0,1.0,1.0)); let body=source.principal_body(1e-12,||false).unwrap();
    let state=body.state_from_kinematics(pose(-2.0),Twist::zero()).unwrap();
    let load=medium.load(&cube(),&body,state,||false).unwrap();
    close(body.properties().mass(),1000.0); close(load.wrench().force_world.z,80000.0);
}

#[test]
fn pressure_potential_gradient_matches_force_and_principal_frame_torque() {
    let sides=Vec3::new(1.2,1.8,0.8); let source=body(800.0,sides);
    let body=source.principal_body(1e-12,||false).unwrap();
    let mesh=DisplacementMesh::cuboid(sides,||false).unwrap();
    let water=fluid(1000.0); let medium=HydrostaticMedium::new(&water,gravity(),origin()).unwrap();
    let rotation=So3::exp(So3Tangent::new(Vec3::new(0.15,0.3,0.1))).unwrap();
    let state=body.state_from_kinematics(Se3::from_parts(rotation,Vec3::new(0.1,-0.2,0.05)).unwrap(),Twist::zero()).unwrap();
    let load=medium.load(&mesh,&body,state,||false).unwrap();
    let h=1e-6;
    for axis in 0..3 {
        let unit=match axis {0=>fs_mbd::Vec3::new(1.0,0.0,0.0),1=>fs_mbd::Vec3::new(0.0,1.0,0.0),_=>fs_mbd::Vec3::new(0.0,0.0,1.0)};
        for rotational in [false,true] {
            let energy=|sign: f64| {
                let old=state.pose();
                let pos=if rotational {old.position_world()} else {old.position_world().add(unit.scale(sign*h))};
                let orient=if rotational {old.orientation().right_exp(unit.scale(sign*h)).unwrap()} else {old.orientation()};
                let changed=RigidBodyState::new(Pose::new(pos,orient).unwrap(),state.linear_momentum_world(),state.angular_momentum_body()).unwrap();
                medium.load(&mesh,&body,changed,||false).unwrap().pressure_potential_j()
            };
            let force=if rotational {load.wrench().torque_body} else {load.wrench().force_world};
            let expected=force.dot(unit); let derivative=-(energy(1.0)-energy(-1.0))/(2.0*h);
            assert!((derivative-expected).abs()<1e-5,"axis={axis}, rotation={rotational}: {derivative} != {expected}");
        }
    }
}

#[test]
fn shared_fluid_is_resolved_once_and_cancellation_returns_no_partial_result() {
    struct Counted(std::cell::Cell<usize>);
    impl FluidTransport for Counted {
        fn transport_state(&self)->Result<FluidTransportState,fs_material::MaterialError> {
            self.0.set(self.0.get()+1); Ok(fluid(1000.0))
        }
    }
    let source=Counted(std::cell::Cell::new(0));
    let medium=HydrostaticMedium::new(&source,gravity(),origin()).unwrap();
    assert!(std::ptr::eq(medium.source(),&source));
    let material=body(500.0,Vec3::new(2.0,2.0,2.0)); let body=material.principal_body(1e-12,||false).unwrap();
    let state=body.state_from_kinematics(Se3::identity(),Twist::zero()).unwrap(); let mesh=cube();
    let mut count=0;
    let expected=medium.load(&mesh,&body,state,||{count+=1;false}).unwrap();
    for stop in 1..=count {
        let mut n=0;
        assert!(matches!(medium.load(&mesh,&body,state,||{n+=1;n==stop}),Err(HydrostaticError::Cancelled)));
    }
    assert_eq!(medium.load(&mesh,&body,state,||false).unwrap(),expected);
    assert_eq!(source.0.get(),1);
    assert!(DisplacementMesh::new(vec![],1,||false).is_err());
    assert!(DisplacementMesh::new(mesh.cells().to_vec(),5,||false).is_err());
    assert!(HydrostaticMedium::new(&source,Gravity::ZERO,origin()).is_err());
    assert!(HydrostaticMedium::new(&source,gravity(),Vec3::new(f64::NAN,0.0,0.0)).is_err());
}
