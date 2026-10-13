//! Source/observer integration on an actual closed thin polyhedron. Its imposed
//! vertical motion is a kinematic fixture, not a cymbal eigenmode or measurement.
use super::*;

fn thin_source() -> Boundary {
    let sphere=SpherePanels::icosphere(0.02,0).unwrap();
    let mut triangles=sphere.triangles().unwrap().to_vec();
    for p in triangles.iter_mut().flatten() { p[2]*=0.05; }
    let surface=SpherePanels::from_triangles(triangles.clone()).unwrap();
    let row:Vec<_>=surface.normals().iter().map(|n|n[2]).collect();
    Boundary {triangles,weights:vec![row],state_modes:vec![0]}
}
fn spec() -> radiation_spec::Spec {
    radiation_spec::Spec {image:radiation_spec::Image::Reciprocal,
        band_hz:[40.0,400.0],training_intervals:16,max_order:8,
        max_panels:80,max_dense_work:2_000_000_000,..Default::default()}
}

#[test]
fn image_selection_is_explicit_and_existing_files_retain_their_operator() {
    let text=include_str!("reciprocal-shell.radiation");
    let reciprocal=radiation_spec::Spec::parse(text).unwrap();
    assert_eq!(reciprocal.image,radiation_spec::Image::Reciprocal);
    let legacy=radiation_spec::Spec::parse(&text.replace("image,reciprocal\n","")).unwrap();
    assert_eq!(legacy.image,radiation_spec::Image::Collocation);
    assert_eq!(legacy.band_hz,reciprocal.band_hz);
    assert_eq!(legacy.frequencies().unwrap(),reciprocal.frequencies().unwrap());
    for replacement in ["image,unknown", "image,reciprocal\nimage,collocation", "image,reciprocal,extra"] {
        assert!(radiation_spec::Spec::parse(&text.replace("image,reciprocal",replacement)).is_err());
    }
    for name in ["splash-mic","splash-wav","hihat-mic","drum-modal-wav"] {
        reciprocal.admit_command(name).unwrap();
    }
    assert!(reciprocal.admit_command("splash").is_err());
}

#[test]
fn one_way_thin_source_uses_its_own_pressure_field_and_keeps_both_receivers() {
    let b=thin_source();let cfg=spec();let gate=CancelGate::new_clock_free();
    let receivers=[Receiver::FinitePoint([0.0,0.0,0.12]),Receiver::FinitePoint([0.1,0.0,0.08])];
    // Exercises explicit selection through the ordinary, non-feedback front
    // door. The reciprocal branch must not write into an absent impedance bank.
    let (baked,load)=bake_scene_with_spec(&b,&receivers,cfg,&gate,false).unwrap();
    assert!(load.is_none());assert_eq!(baked.len(),2);
    assert!(baked.iter().all(|b|b.filters.len()==1 && b.pressure_gain==1.0));
    let surface=SpherePanels::from_triangles(b.triangles.clone()).unwrap();
    let radius=b.triangles.iter().flatten().map(|p|p[0].hypot(p[1]).hypot(p[2]))
        .fold(0.0_f64,f64::max);
    let medium=Medium::air();let omega=cfg.frequencies().unwrap();
    let source=reciprocal_geometry(&surface,omega[omega.len()-1]/medium.sound_speed,
        cfg,1,omega.len(),&gate).unwrap();
    let scene=receivers::Scene::new(&surface,&receivers,radius,medium,1.0/f64::from(OUTPUT_RATE),&gate).unwrap();
    let observations=scene.reciprocal(&source,&gate).unwrap();
    for hz in [117.0,233.0,379.0] { // absent from training, selection and audit
        let w=core::f64::consts::TAU*hz;let fields=acceleration_fields(&b.weights,w);
        let solutions=source.solve_batch(w/medium.sound_speed,medium,&[&fields[0]],||false).unwrap();
        let receiver=observations.prepare(w/medium.sound_speed,&gate).unwrap();
        for channel in 0..2 {
            let truth=receiver.response(channel,&solutions[0]).unwrap();
            assert!(truth.abs()>1e-12);
            let actual=baked[channel].filters[0].eval(w).unwrap().conj();
            assert!((actual-truth).abs()<0.15*truth.abs(),"channel {channel}, {hz} Hz");
        }
    }
    assert_ne!(baked[0].propagation_delay_s,baked[1].propagation_delay_s);
}

#[test]
fn signed_velocity_and_acceleration_sources_have_the_same_reciprocal_work() {
    let mut b=thin_source();b.weights.push(b.weights[0].iter().map(|v|-0.4*v).collect());
    let surface=SpherePanels::from_triangles(b.triangles.clone()).unwrap();
    let gate=CancelGate::new_clock_free();let cfg=spec();let medium=Medium::air();
    let w=core::f64::consts::TAU*233.0;
    let source=reciprocal_geometry(&surface,w/medium.sound_speed,cfg,2,1,&gate).unwrap();
    let acc=acceleration_fields(&b.weights,w);
    let a=source.solve_batch(w/medium.sound_speed,medium,&[&acc[0],&acc[1]],||false).unwrap();
    let traces:Vec<&[C64]>=a.iter().map(|s|s.pressure.as_slice()).collect();
    let z=feedback::project_pressures(&surface,&b.weights,&traces,w).unwrap();
    let velocity:Vec<Vec<_>>=b.weights.iter().map(|row|row.iter().map(|v|C64::from_re(*v)).collect()).collect();
    let v=source.solve_batch(w/medium.sound_speed,medium,&[&velocity[0],&velocity[1]],||false).unwrap();
    for j in 0..2 {for i in 0..2 {
        let truth=v[j].pressure.iter().zip(surface.areas()).zip(&b.weights[i])
            .fold(C64::ZERO,|sum,((p,a),b)|sum+p.scale(a*b));
        assert!((z[2*i+j]-truth).abs()<1e-9*truth.abs());
    }}
    assert!(z[0].re>0.0 && z[1].re<0.0);
    assert!((z[1]-z[2]).abs()<1e-9*z[0].abs());
}

#[test]
fn refused_reciprocal_work_or_band_does_not_fall_back_to_collocation() {
    let b=thin_source();let gate=CancelGate::new_clock_free();
    let receiver=[Receiver::FinitePoint([0.0,0.0,0.12])];
    // Satisfies the collocation work preflight, but not the Galerkin owner.
    let low=radiation_spec::Spec {max_dense_work:1_000_000,..spec()};
    assert!(low.work(b.triangles.len(),1,1).is_ok());
    assert!(bake_scene_with_spec(&b,&receiver,low,&gate,false).is_err());
    gate.request();
    assert!(bake_scene_with_spec(&b,&receiver,spec(),&gate,false).is_err());
    let sphere=SpherePanels::icosphere(0.02,0).unwrap();
    let band=radiation_spec::Spec {band_hz:[40.0,10_000.0],max_dense_work:2_000_000_000,..spec()};
    let gate=CancelGate::new_clock_free();
    assert!(reciprocal_geometry(&sphere,core::f64::consts::TAU*band.band_hz[1]/Medium::air().sound_speed,
        band,1,65,&gate).is_err(),"beyond the owner's slab/resonance chart");
}
