//! Physical drive/reaction oracles on original constrained mechanism pencils.
use super::*;
use response::forced::{ConstrainedForcedControls,ConstrainedForcedError};
use fs_modal::quadratic_modes::PoleOptions;
use fs_modal::quadratic_modes::forced::QuadraticForcedOptions;

fn pole_options()->PoleOptions {
    PoleOptions {spectral_scale:1.0,step_tolerance:1e-11,residual_tolerance:1e-11,
        max_iterations:8,max_dimension:8,max_dense_bytes:4_000_000,max_work:100_000_000}
}
fn force_options()->QuadraticForcedOptions {
    QuadraticForcedOptions {maximum_relative_residual:1e-10,
        maximum_payload_bytes:1_000_000,maximum_work:100_000_000}
}
fn reaction_options()->ConstrainedForcedControls {
    ConstrainedForcedControls {maximum_reaction_power_w:1e-9}
}
fn cnear(a:C64,b:C64,tol:f64) {close(a.re,b.re,tol);close(a.im,b.im,tol);}
fn norm_squared(z:C64)->f64 {z.re*z.re+z.im*z.im}

fn geared_pencil(scale:f64)->response::ConstrainedQuadraticModel {
    let (model,connections)=two_sliders();
    let loads=SpatialMotionLoads {
        held:FreeFloatingMotionLoads {generalized_force:&[10.0,40.0],
            external_body_wrench:&[Wrench::default();3],gravity_world:Vec3::default()},
        world_points:&[],connections:&connections,
    };
    let mut rows=[row(vec![-2.0,1.0],0.1)];
    let mut root=static_owner::equilibrate(&model,&seed(vec![0.12,-0.04]),loads,
        EquilibriumBushings::Diagonal(&[]),&rows,controls(&[0.1;2],&[40.0;2]),
        joint_controls(&[20.0],&[0.1]),||false).unwrap();
    root.multipliers[0]/=scale;
    for a in &mut rows[0].coefficients {*a*=scale;}
    rows[0].position*=scale;
    rows[0].position_tolerance*=scale.abs();
    rows[0].velocity_tolerance*=scale.abs();
    response::prepare(&model,&root,loads,EquilibriumBushings::Diagonal(&[]),&rows,
        response_controls(&[0.1,10.0],&[0.2,0.4],&[40.0;2]),
        subspace_controls(),||false).unwrap().quadratic_model(1,64,||false).unwrap()
}

#[test]
fn geared_dc_and_harmonic_loads_recover_original_reactions_and_virtual_work() {
    let force=[C64::new(6.0,-1.0),C64::new(-2.0,3.0)];
    for scale in [1.0_f64,-3.0,1e-6,1e6] {
        let pencil=geared_pencil(scale);
        for omega in [0.0,2.0,8.0,30.0] {
            let factor=pencil.prepare_forced_response(omega,pole_options(),force_options(),
                reaction_options(),||false).unwrap();
            let response=factor.solve(&force,||false).unwrap();
            // q1=2*q0; native M=diag(2,3), C=diag(7,7), K=diag(100,200).
            let d0=C64::new(100.0-2.0*omega*omega,-7.0*omega);
            let d1=C64::new(200.0-3.0*omega*omega,-7.0*omega);
            let q0=(force[0]+force[1].scale(2.0))/(d0+d1.scale(4.0));
            let lambda=d1*q0.scale(2.0)-force[1];
            cnear(response.displacement[0],q0,2e-8);
            cnear(response.displacement[1],q0.scale(2.0),4e-8);
            cnear(response.reaction_effort[0],lambda.scale(-2.0),2e-5);
            cnear(response.reaction_effort[1],lambda,1e-5);
            cnear(response.multipliers[0].scale(scale),lambda,1e-5);
            close(response.balance.average_input_power_w,
                0.5*omega*omega*7.0*5.0*norm_squared(q0),3e-7);
            close(response.balance.average_reaction_power_w,0.0,1e-9);
            close(response.balance.average_stiffness_power_w,0.0,1e-9);
            close(response.balance.average_inertial_power_w,0.0,1e-9);
            assert!(response.balance.relative_force_residual<=1e-10);
            assert!(response.reduced_relative_force_residual<=1e-10);
            assert!(response.position_residuals[0].abs()<=1e-10*scale.abs());
        }
    }
}

#[test]
fn full_lock_and_no_constraint_limits_keep_physical_loads_and_empty_modes_distinct() {
    let (model,connection)=slider();
    let loads=SpatialMotionLoads {
        held:FreeFloatingMotionLoads {generalized_force:&[10.0],
            external_body_wrench:&[Wrench::default();2],gravity_world:Vec3::default()},
        world_points:&[],connections:&[connection],
    };
    let force=[C64::new(3.0,-4.0)];
    let rows=[row(vec![1.0],0.2)];
    let root=static_owner::equilibrate(&model,&seed(vec![0.0]),loads,
        EquilibriumBushings::Diagonal(&[]),&rows,controls(&[0.1],&[20.0]),
        joint_controls(&[20.0],&[0.1]),||false).unwrap();
    let locked=response::prepare(&model,&root,loads,EquilibriumBushings::Diagonal(&[]),
        &rows,response_controls(&[0.1],&[0.2],&[20.0]),subspace_controls(),||false)
        .unwrap().quadratic_model(0,0,||false).unwrap();
    assert_eq!(locked.source_evaluations(),0);
    for omega in [0.0,17.0] {
        let x=locked.prepare_forced_response(omega,pole_options(),force_options(),
            reaction_options(),||false).unwrap().solve(&force,||false).unwrap();
        assert!(x.reduced_displacement.is_empty());
        assert_eq!(x.displacement,vec![C64::ZERO]);
        assert_eq!(x.velocity,vec![C64::ZERO]);
        cnear(x.reaction_effort[0],force[0].scale(-1.0),2e-15);
        cnear(x.multipliers[0],force[0].scale(-1.0),2e-15);
        assert_eq!(x.balance.average_reaction_power_w,0.0);
    }
    let root=static_owner::equilibrate(&model,&seed(vec![0.0]),loads,
        EquilibriumBushings::Diagonal(&[]),&[],controls(&[0.1],&[20.0]),
        joint_controls(&[],&[]),||false).unwrap();
    let free=response::prepare(&model,&root,loads,EquilibriumBushings::Diagonal(&[]),
        &[],response_controls(&[0.1],&[0.2],&[20.0]),subspace_controls(),||false)
        .unwrap().quadratic_model(1,64,||false).unwrap();
    let x=free.prepare_forced_response(3.0,pole_options(),force_options(),
        reaction_options(),||false).unwrap().solve(&force,||false).unwrap();
    cnear(x.displacement[0],force[0]/C64::new(82.0,-21.0),2e-8);
    assert!(x.multipliers.is_empty());
    assert_eq!(x.reaction_effort,vec![C64::ZERO]);
}

fn rotating_pencil() -> response::ConstrainedQuadraticModel {

    let model = ArticulatedModel::new(vec![
        Link::new("support", None, Se3::identity(), JointModel::FIXED,
            inertia(3.0, Vec3::default())),
        Link::new("x carriage", Some(0), Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(1.0, Vec3::default())),
        Link::new("y carriage", Some(1), Se3::identity(),
            JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
            inertia(2.0, Vec3::default())),
        Link::new("locked z", Some(0), Se3::identity(),
            JointModel::prismatic(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(4.0, Vec3::default())),
    ]).unwrap();
    let bushings = [SpatialBushing {
        a: MaterialFrame { link: 0, link_from_attachment: Se3::identity() },
        b: MaterialFrame { link: 2, link_from_attachment: Se3::identity() },
        rest_translation_a: Vec3::default(),
        translation_stiffness: Vec3::new(100.0, 120.0, 0.0),
        translation_damping: Vec3::new(4.0, 5.0, 0.0),
        rotation_stiffness: Vec3::default(), rotation_damping: Vec3::default(),
    }];
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0, 0.0, 2.0],
            external_body_wrench: &[Wrench::default(); 4], gravity_world: Vec3::default(),
        }, world_points: &[], connections: &[],
    };
    let rows = [row(vec![0.0, 0.0, 1.0], 0.0)];
    let mut initial = seed(vec![0.0; 3]);
    initial.base.twist_body.angular = Vec3::new(0.0, 0.0, 3.0);
    let root = static_owner::equilibrate_rotating(&model, &initial, loads,
        EquilibriumBushings::Diagonal(&bushings), &rows,
        controls(&[1.0; 3], &[100.0; 3]), joint_controls(&[10.0], &[1.0]), || false).unwrap();
    let prepared = response::prepare(&model, &root, loads,
        EquilibriumBushings::Diagonal(&bushings), &rows,
        response_controls(&[1.0; 3], &[1.0; 3], &[100.0; 3]),
        subspace_controls(), || false).unwrap();
    let pencil = prepared.quadratic_model(2, 128, || false).unwrap();
    pencil
}

#[test]
fn spinning_two_axis_forcing_retains_gyroscopic_cross_terms_and_locked_reaction() {
    let pencil=rotating_pencil();
    let omega=2.0;
    let q=[C64::new(0.01,0.003),C64::new(-0.02,0.004)];
    // Original rotating tangent: M=diag(3,2), K=diag(73,102),
    // C=[[4,-12],[12,5]]. The separate z joint is ideally locked.
    let force=[
        C64::new(73.0-3.0*omega*omega,-4.0*omega)*q[0]+C64::new(0.0,12.0*omega)*q[1],
        C64::new(0.0,-12.0*omega)*q[0]+C64::new(102.0-2.0*omega*omega,-5.0*omega)*q[1],
        C64::new(3.0,-2.0),
    ];
    let x=pencil.prepare_forced_response(omega,pole_options(),force_options(),
        reaction_options(),||false).unwrap().solve(&force,||false).unwrap();
    cnear(x.displacement[0],q[0],2e-8);
    cnear(x.displacement[1],q[1],2e-8);
    assert_eq!(x.displacement[2],C64::ZERO);
    cnear(x.reaction_effort[0],C64::ZERO,1e-8);
    cnear(x.reaction_effort[1],C64::ZERO,1e-8);
    cnear(x.reaction_effort[2],force[2].scale(-1.0),1e-12);
    close(x.balance.average_input_power_w,
        0.5*omega*omega*(4.0*norm_squared(q[0])+5.0*norm_squared(q[1])),2e-8);
    close(x.balance.average_reaction_power_w,0.0,1e-10);
}

fn inverted_pencil() -> response::ConstrainedQuadraticModel {

    let model = ArticulatedModel::new(vec![
        Link::new("support", None, Se3::identity(), JointModel::FIXED,
            inertia(3.0, Vec3::default())),
        Link::new("inverted pendulum", Some(0), Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(2.0, Vec3::new(0.0, 1.0, 0.0))),
        Link::new("locked slider", Some(0), Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(4.0, Vec3::default())),
    ]).unwrap();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0, 3.0], external_body_wrench: &[Wrench::default(); 3],
            gravity_world: Vec3::new(0.0, -9.81, 0.0),
        }, world_points: &[], connections: &[],
    };
    let rows = [row(vec![0.0, 1.0], 0.2)];
    let root = static_owner::equilibrate(&model, &seed(vec![0.0, 0.2]), loads,
        EquilibriumBushings::Diagonal(&[]), &rows, controls(&[1.0; 2], &[20.0; 2]),
        joint_controls(&[10.0], &[1.0]), || false).unwrap();
    let prepared = response::prepare(&model, &root, loads, EquilibriumBushings::Diagonal(&[]),
        &rows, response_controls(&[1.0; 2], &[1.0; 2], &[20.0; 2]),
        subspace_controls(), || false).unwrap();
    let pencil = prepared.quadratic_model(1, 64, || false).unwrap();
    pencil
}

#[test]
fn unstable_pendulum_has_signed_compliance_without_claiming_transient_decay() {
    let pencil=inverted_pencil();
    let force=[C64::new(0.2,0.1),C64::new(3.0,-2.0)];
    for omega in [0.0,1.0,4.0] {
        let x=pencil.prepare_forced_response(omega,pole_options(),force_options(),
            reaction_options(),||false).unwrap().solve(&force,||false).unwrap();
        cnear(x.displacement[0],force[0].scale(1.0/(-19.62-3.0*omega*omega)),2e-8);
        assert_eq!(x.displacement[1],C64::ZERO);
        cnear(x.reaction_effort[1],force[1].scale(-1.0),2e-14);
        close(x.balance.average_input_power_w,0.0,1e-11);
        close(x.balance.average_reaction_power_w,0.0,1e-11);
    }
}

#[test]
fn invalid_physical_drives_budgets_and_final_cancellation_do_not_mutate_factors() {
    let pencil=geared_pencil(1.0);
    let before=(pencil.mass().to_vec(),pencil.stiffness().to_vec(),pencil.velocity_operator().to_vec());
    let mut small=force_options();small.maximum_payload_bytes=1;
    assert!(pencil.prepare_forced_response(2.0,pole_options(),small,reaction_options(),||false).is_err());
    let factor=pencil.prepare_forced_response(2.0,pole_options(),force_options(),
        reaction_options(),||false).unwrap();
    assert!(factor.solve(&[C64::ONE],||false).is_err());
    assert!(factor.solve(&[C64::new(f64::NAN,0.0),C64::ZERO],||false).is_err());
    let force=[C64::new(6.0,-1.0),C64::new(-2.0,3.0)];
    let mut checkpoints=0;
    let expected=factor.solve(&force,||{checkpoints+=1;false}).unwrap();
    let mut calls=0;
    assert!(matches!(factor.solve(&force,||{calls+=1;calls==checkpoints}),
        Err(ConstrainedForcedError::Cancelled)));
    assert_eq!(factor.solve(&force,||false).unwrap(),expected);
    assert_eq!((pencil.mass().to_vec(),pencil.stiffness().to_vec(),pencil.velocity_operator().to_vec()),before);
}
