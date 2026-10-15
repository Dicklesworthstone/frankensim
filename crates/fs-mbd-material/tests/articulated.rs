mod common;
use common::*;
use fs_ga::{Se3, So3, Vec3, Wrench};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::articulated::{ArticulatedError, BaseState, FreeFloatingBaseState, JointModel,
    forward_dynamics, forward_kinematics, free_floating_forward_dynamics, inverse_dynamics};
use fs_mbd_material::{MaterialInertiaError, RegionBinding, VolumeMoments};
use fs_mbd_material::articulated::{LinkMaterials, MaterialLink, ResolvedArticulated,
    rebind_articulated, resolve_articulated};

fn zero() -> Vec3 { Vec3::new(0.0,0.0,0.0) }
fn no_wrench() -> Wrench { Wrench::new(zero(),zero()) }
fn base() -> BaseState { BaseState::stationary(Se3::identity()) }
fn placed(x: f64) -> Se3 { Se3::from_parts(So3::identity(), Vec3::new(x,0.0,0.0)).unwrap() }
fn joint() -> JointModel { JointModel::revolute(Vec3::new(0.0,0.0,1.0), None).unwrap() }
fn bar<'a>(material: &'a MaterialCard, point: &'a QueryPoint, center: f64) -> RegionBinding<'a> {
    RegionBinding { key: "bar", material, point,
        geometry: VolumeMoments::cuboid(Vec3::new(1.0,0.2,0.2)).unwrap().transformed(placed(center)).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly }
}
fn lever(card: &MaterialCard) -> ResolvedArticulated {
    let p = point(300.0); let regions = [bar(card, &p, 1.0)];
    resolve_articulated(&[MaterialLink { name: "lever", parent: None,
        parent_from_child_zero: Se3::identity(), joint: joint(), regions: &regions }], 1, || false).unwrap()
}
fn chain(first: &MaterialCard, second: &MaterialCard) -> ResolvedArticulated {
    let p = point(300.0); let a = [bar(first,&p,0.5)]; let b = [bar(second,&p,0.5)];
    resolve_articulated(&[
        MaterialLink { name:"first", parent:None, parent_from_child_zero:Se3::identity(), joint:joint(), regions:&a },
        MaterialLink { name:"second", parent:Some(0), parent_from_child_zero:placed(1.0), joint:joint(), regions:&b },
    ], 2, || false).unwrap()
}

#[test]
fn data_only_density_swap_changes_actual_torque_acceleration_and_gravity_load() {
    let a = lever(&density_card("a",1000.0,1000.0));
    let b = lever(&density_card("b",3000.0,3000.0));
    let torque = 2.0;
    let accelerate = |m: &ResolvedArticulated| forward_dynamics(m.model(),base(),&[0.0],&[0.0],
        &[torque],zero(),&[no_wrench()]).unwrap().generalized_acceleration[0];
    let i_pivot = 40.0*((1.0+0.04)/12.0+1.0);
    close(accelerate(&a), torque/i_pivot);
    close(accelerate(&b), torque/(3.0*i_pivot));
    let angle: f64 = 0.37;
    let hold = |m: &ResolvedArticulated| inverse_dynamics(m.model(),base(),&[angle],&[0.0],&[0.0],
        Vec3::new(0.0,-9.81,0.0),&[no_wrench()]).unwrap().generalized_force[0];
    close(hold(&a),40.0*9.81*angle.cos());
    close(hold(&b),3.0*hold(&a));
    // Uniform gravity acceleration is mass-independent, unlike imposed torque.
    let fall = |m: &ResolvedArticulated| forward_dynamics(m.model(),base(),&[angle],&[0.0],&[0.0],
        Vec3::new(0.0,-9.81,0.0),&[no_wrench()]).unwrap().generalized_acceleration[0];
    close(fall(&a),fall(&b));
}

#[test]
fn changing_one_link_changes_coupled_motion_and_preserves_inverse_forward_closure() {
    let a = density_card("a",1000.0,1000.0); let b = density_card("b",3000.0,3000.0);
    let original = chain(&a,&a); let substituted = chain(&a,&b);
    let q = [0.3,-0.7]; let v = [0.4,-0.2]; let target = [0.8,-0.6];
    let gravity = Vec3::new(0.0,-9.81,0.0);
    let external = [no_wrench(), Wrench::new(Vec3::new(0.0,0.0,0.3),Vec3::new(0.4,-0.2,0.0))];
    for specimen in [&original,&substituted] {
        let inverse = inverse_dynamics(specimen.model(),base(),&q,&v,&target,gravity,&external).unwrap();
        let forward = forward_dynamics(specimen.model(),base(),&q,&v,&inverse.generalized_force,gravity,&external).unwrap();
        for (actual,expected) in forward.generalized_acceleration.iter().zip(target) { close(*actual,expected); }
    }
    let accelerate = |m: &ResolvedArticulated| forward_dynamics(m.model(),base(),&q,&v,
        &[1.0,-0.5],gravity,&external).unwrap().generalized_acceleration;
    let left = accelerate(&original); let right = accelerate(&substituted);
    // Independent planar two-link Lagrange equations, not another call to RNEA:
    // M11=I1+I2+m1*c²+m2*(L²+c²+2*L*c*cos(q2)),
    // M12=I2+m2*(c²+L*c*cos(q2)), M22=I2+m2*c²;
    // c=0.5, L=1 and I_i=m_i*(1²+0.2²)/12. Include gravity,
    // the Coriolis vector and the declared force at the second joint origin.
    for (actual,expected) in left.iter().zip([-9.93545972704422,7.940200164815278]) { close(*actual,expected); }
    for (actual,expected) in right.iter().zip([-8.648925572511027,5.202187401218664]) { close(*actual,expected); }
    assert!((left[0]-right[0]).abs() > 1e-4);
    assert!((left[1]-right[1]).abs() > 1e-4);
    assert_eq!(forward_kinematics(original.model(),base(),&q,&v).unwrap(),
        forward_kinematics(substituted.model(),base(),&q,&v).unwrap());
}

#[test]
fn complete_rebinding_keeps_tree_frames_and_joint_coordinates_without_stale_mass() {
    let a = density_card("a",1000.0,1000.0); let b = density_card("b",3000.0,3000.0);
    let original = chain(&a,&a); let saved = original.clone(); let p = point(300.0);
    let first = [bar(&a,&p,0.5)]; let second = [bar(&b,&p,0.5)];
    // Intentionally reversed; the model remains in its original parent-first order.
    let replaced = rebind_articulated(original.model(), &[
        LinkMaterials { link:"second", regions:&second }, LinkMaterials { link:"first", regions:&first },
    ], 2, || false).unwrap();
    assert_eq!(replaced,chain(&a,&b));
    assert_eq!(original,saved);
    assert_eq!(replaced.model().dof_count(),original.model().dof_count());
    for (new,old) in replaced.model().links().iter().zip(original.model().links()) {
        assert_eq!(new.name(),old.name()); assert_eq!(new.parent(),old.parent());
        assert_eq!(new.joint(),old.joint()); assert_eq!(new.parent_from_child_zero(),old.parent_from_child_zero());
    }
    assert_eq!(replaced.bodies()[1].inertia(),replaced.model().links()[1].inertia());
    assert_ne!(replaced.bodies()[1].regions()[0].material(),original.bodies()[1].regions()[0].material());
}

#[test]
fn free_floating_body_uses_the_same_resolved_material_inertia() {
    let p=point(300.0); let gravity=Vec3::new(0.0,0.0,-9.81);
    for density in [2.0,6.0] {
        let material = density_card("body",density,density);
        let regions=[RegionBinding { key:"solid", material:&material, point:&p,
            geometry:VolumeMoments::cuboid(Vec3::new(1.0,1.0,1.0)).unwrap(),
            selection:MaterialPropertySelection::SingleClaimOnly }];
        let body=resolve_articulated(&[MaterialLink { name:"free", parent:None,
            parent_from_child_zero:Se3::identity(),joint:JointModel::FIXED,regions:&regions }],1,||false).unwrap();
        let force=Wrench::new(Vec3::new(0.0,0.0,2.0),Vec3::new(4.0,0.0,0.0));
        let result=free_floating_forward_dynamics(body.model(),FreeFloatingBaseState::stationary(Se3::identity()),
            &[],&[],&[],gravity,&[force]).unwrap();
        close(result.base_spatial_acceleration_body.linear.x,4.0/density);
        close(result.base_spatial_acceleration_body.linear.z,-9.81);
        close(result.base_spatial_acceleration_body.angular.z,12.0/density);
    }
}

#[test]
fn invalid_full_bindings_budget_or_source_cannot_change_an_existing_mechanism() {
    let a=density_card("a",1.0,1.0); let original=chain(&a,&a); let saved=original.clone(); let p=point(300.0);
    let regions=[bar(&a,&p,0.5)];
    let binding=|name| LinkMaterials { link:name, regions:&regions };
    for bindings in [vec![binding("first")], vec![binding("first"),binding("first")],
        vec![binding("first"),binding("foreign")]] {
        assert!(rebind_articulated(original.model(),&bindings,2,||false).is_err());
    }
    let bindings=[binding("first"),binding("second")];
    assert!(rebind_articulated(original.model(),&bindings,1,||false).is_err());
    let outside=point(500.0); let bad=[bar(&a,&outside,0.5)];
    assert!(matches!(rebind_articulated(original.model(),&[binding("first"),
        LinkMaterials {link:"second",regions:&bad}],2,||false),Err(MaterialInertiaError::Material(_))));
    let mut polls=0;
    rebind_articulated(original.model(),&bindings,2,||{polls+=1;false}).unwrap();
    let final_poll=polls; let mut polls=0;
    assert_eq!(rebind_articulated(original.model(),&bindings,2,||{polls+=1;polls==final_poll}),
        Err(MaterialInertiaError::Cancelled));
    assert_eq!(original,saved);
}

#[test]
fn global_budget_counts_regions_not_just_links_and_original_topology_guard_remains() {
    let material=density_card("body",2.0,2.0); let p=point(300.0);
    let mut left=bar(&material,&p,-0.5); left.key="left";
    let mut right=bar(&material,&p,0.5); right.key="right";
    let pair=[left,right];
    let links=[MaterialLink { name:"root",parent:None,parent_from_child_zero:Se3::identity(),joint:joint(),regions:&pair },
        MaterialLink { name:"child",parent:Some(0),parent_from_child_zero:placed(2.0),joint:joint(),regions:&pair }];
    assert!(resolve_articulated(&links,3,||false).is_err());
    assert!(resolve_articulated(&links,4,||false).is_ok());
    let mut broken=links.clone(); broken[1].parent=Some(1);
    assert!(matches!(resolve_articulated(&broken,4,||false),
        Err(MaterialInertiaError::Dynamics(ArticulatedError::InvalidParent {..}))));
}
