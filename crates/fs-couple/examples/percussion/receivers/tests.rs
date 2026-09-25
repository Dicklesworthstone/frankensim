use super::*;
use super::super::{Boundary, Bake, zero_filter, OUTPUT_RATE, MECHANICAL_DT, SUBSTEPS, acceleration_fields};
use fs_bem::helmholtz::{Formulation, solve_radiation_batch};

fn options() -> Options {Options {relative_tolerance:1e-6, maximum_depth:14, maximum_kernel_evaluations:2_000_000}}
fn near(p:[f64;3], alpha:f64, front:[f64;3]) -> Receiver {
    Receiver::NearField(Microphone::new(p,0.001,options(),FirstOrder::new(alpha,front).unwrap()).unwrap())
}
// Outward, closed flattened body. This is a prescribed breathing surface, not
// a calibrated material or mode-frequency fixture. Its centre and enclosing
// sphere cannot decide whether a microphone is exterior to the actual faces.
fn box_boundary() -> Boundary {
    let points=[[-0.02,-0.01,-0.003],[0.02,-0.01,-0.003],[0.02,0.01,-0.003],[-0.02,0.01,-0.003],
        [-0.02,-0.01,0.003],[0.02,-0.01,0.003],[0.02,0.01,0.003],[-0.02,0.01,0.003]];
    let faces=[[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
        [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]];
    Boundary {triangles:faces.map(|t|t.map(|i|points[i])).to_vec(),weights:vec![vec![1.0;12]],state_modes:vec![0]}
}
fn radius(b:&Boundary) -> f64 {
    b.triangles.iter().flatten().map(|p|p.iter().map(|v|v*v).sum::<f64>().sqrt()).fold(0.0,f64::max)
}

#[test]
fn actual_clearance_admits_close_receivers_inside_sphere_and_rejects_interior() {
    let b=box_boundary();let surface=SpherePanels::from_triangles(b.triangles.clone()).unwrap();
    let medium=Medium::air();let dt=1.0/f64::from(OUTPUT_RATE);let gate=CancelGate::new_clock_free();
    let p=[0.005,0.0,0.008];
    assert!(Receiver::FinitePoint(p).propagation(radius(&b),medium,dt).is_err());
    let receivers=[near(p,1.0,[0.0,0.0,-1.0]),near([0.0,0.0,-0.025],1.0,[0.0,0.0,1.0])];
    let scene=Scene::new(&surface,&receivers,radius(&b),medium,dt,&gate).unwrap();
    assert_eq!(scene.gain(0),1.0);assert_eq!(scene.delay(0),0.0);
    assert!((scene.delay(1)-0.022/medium.sound_speed).abs()<1e-15);
    for p in [[0.0;3],[0.0,0.0,0.003],[0.0,0.0,0.0035]] {
        assert!(Scene::new(&surface,&[near(p,1.0,[0.0,0.0,1.0])],radius(&b),medium,dt,&gate).is_err());
    }
    gate.request();assert!(Scene::new(&surface,&receivers,radius(&b),medium,dt,&gate).is_err());
}

#[test]
fn prepared_receiver_matches_actual_integrated_fields_without_changing_source_load() {
    let b=box_boundary();let surface=SpherePanels::from_triangles(b.triangles.clone()).unwrap();
    let medium=Medium::air();let gate=CancelGate::new_clock_free();let dt=1.0/f64::from(OUTPUT_RATE);
    let point=[0.005,0.0,0.008];let receivers=[near(point,0.5,[0.0,0.0,-1.0]),near(point,0.5,[0.0,0.0,1.0])];
    let scene=Scene::new(&surface,&receivers,radius(&b),medium,dt,&gate).unwrap();
    for hz in [80.0,320.0,1000.0] {
        let w=core::f64::consts::TAU*hz;let k=w/medium.sound_speed;
        let fields=acceleration_fields(&b.weights,w);
        let solutions=solve_radiation_batch(&surface,k,medium,&[&fields[0]],Formulation::PlainCbie).unwrap();
        // Independently integrate the unchanged surface work; no receiver
        // projection or public visibility change to the loading owner is needed.
        let project=||solutions[0].pressure.iter().zip(surface.areas()).zip(&b.weights[0])
            .fold(C64::ZERO,|sum,((p,a),b)|sum+p.scale(a*b))*C64::new(0.0,-w);
        let before=project();
        let frequency=scene.prepare(k,&gate).unwrap();
        let direct=Geometry::new(&surface,&[point],0.001).unwrap()
            .prepare_velocity(k,medium,options()).unwrap().evaluate_velocity(&solutions[0]).unwrap();
        let front=frequency.response(0,&solutions[0]).unwrap();let rear=frequency.response(1,&solutions[0]).unwrap();
        let p=direct.scalar.pressure[0];
        assert!((front+rear-p).abs()<1e-11*p.abs().max(1e-20));
        let want=FirstOrder::new(0.5,[0.0,0.0,-1.0]).unwrap()
            .observe(p,direct.particle_velocity_m_s[0],medium).unwrap();
        assert!((front-want).abs()<1e-12*want.abs().max(1e-20));
        let after=project();
        assert_eq!(before,after);
    }
}

#[test]
fn unpeeled_short_flight_adds_no_fictitious_two_sample_delay() {
    let dt=1.0/f64::from(OUTPUT_RATE);let medium=Medium::air();let mut f=zero_filter(dt);f.d=3.0;
    let delay=peeled_delay(0.001,medium,dt).unwrap();assert_eq!(delay,0.0);
    let b=Bake {filters:vec![f],range_m:0.01,medium,propagation_delay_s:delay,pressure_gain:1.0};
    let mut observer=b.runtime().unwrap();
    assert_eq!(observer.step(&[2.0]).unwrap(),6.0);assert_eq!(observer.step(&[0.0]).unwrap(),0.0);
    for clearance in [0.001,0.02,0.1] {
        let peeled=peeled_delay(clearance,medium,dt).unwrap();
        let full=shift_to_enclosing_sphere(C64::new(0.2,-0.7),2000.0,clearance/medium.sound_speed);
        let residual=shift_to_enclosing_sphere(full,2000.0,-peeled);
        let restored=shift_to_enclosing_sphere(residual,2000.0,peeled);
        assert!((restored-full).abs()<1e-14);
    }
}

#[test]
fn real_head_geometry_and_source_addresses_are_retained_for_close_miking() {
    use fs_plate::{ModePair,shell::head::{TensionedDisk,TensionedDiskSpec},shell::profile::ProfileBudget};
    let film=||TensionedDisk::new(TensionedDiskSpec {radius_m:0.1,thickness_m:0.0002,
        young_pa:4e9,poisson:0.38,density_kg_m3:1390.0,tension_n_m:1000.0,radial_intervals:2,azimuths:8},
        ProfileBudget {max_nodes:100,max_triangles:200,max_feature_evaluations:0}).unwrap();
    let films=[film(),film()];let modes=films.iter().map(|f| {
        let mut phi=vec![0.0;f.model.free];
        for (i,&(x,y)) in f.mesh.nodes.iter().enumerate(){if let Some(k)=f.model.dof_map[3*i]{
            phi[k]=1.0-(x*x+y*y)/0.01;
        }}
        vec![ModePair {lambda:1.0,phi,residual:0.0,interval:(1.0,1.0)}]
    }).collect::<Vec<_>>();
    let b=Boundary::drum(&films,&modes,0.12,0.11).unwrap();
    let surface=SpherePanels::from_triangles(b.triangles.clone()).unwrap();let medium=Medium::air();
    let receivers=[near([0.0,0.0,0.08],1.0,[0.0,0.0,-1.0]),near([0.0,0.0,-0.08],1.0,[0.0,0.0,1.0])];
    let gate=CancelGate::new_clock_free();let scene=Scene::new(&surface,&receivers,radius(&b),medium,1.0/f64::from(OUTPUT_RATE),&gate).unwrap();
    let w=core::f64::consts::TAU*80.0;let fields=acceleration_fields(&b.weights,w);
    let solutions=solve_radiation_batch(&surface,w/medium.sound_speed,medium,
        &[&fields[0],&fields[1]],Formulation::PlainCbie).unwrap();
    let rows=scene.prepare(w/medium.sound_speed,&gate).unwrap();
    for solution in &solutions {for channel in 0..2 {assert!(rows.response(channel,solution).unwrap().abs()>0.0);}}
    assert_eq!(b.state_modes,[1,2]);
}

#[test]
fn geometry_quadrature_and_work_refusals_publish_no_mechanical_step() {
    let mut e=crate::drum(16,MECHANICAL_DT,true,true).unwrap();let before=e.system.state().to_vec();
    let gate=CancelGate::new_clock_free();let spec=super::super::stereo::radiation_spec::Spec::default();
    let bad=[near([0.0;3],1.0,[0.0,0.0,1.0])];
    assert!(super::super::stereo::render_receivers_with_spec(&mut e,1,20.0,&bad,spec,&gate).is_err());
    assert_eq!(before,e.system.state());
    let b=box_boundary();let surface=SpherePanels::from_triangles(b.triangles.clone()).unwrap();
    let receiver=Receiver::NearField(Microphone::new([0.005,0.0,0.008],0.001,
        Options {maximum_kernel_evaluations:80,..options()},FirstOrder::default()).unwrap());
    let rs=[receiver];let scene=Scene::new(&surface,&rs,radius(&b),Medium::air(),1.0/f64::from(OUTPUT_RATE),&gate).unwrap();
    assert!(scene.prepare(1.0,&gate).is_err());
    assert_eq!(before,e.system.state());
    assert_eq!(SUBSTEPS,16); // The receiver never selected a different mechanics clock.
}
