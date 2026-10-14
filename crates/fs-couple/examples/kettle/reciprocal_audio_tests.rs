//! Focused image-selection and actual acoustic field regressions. The imposed
//! head shape below is a kinematic test fixture, not a measured instrument mode.
use super::*;
use fs_couple::render::plate::impact::cavity::meridian::Station;
use fs_plate::shell::{head::TensionedDisk,profile::ProfileBudget};

fn spec() -> Spec { Spec::read(include_str!("illustrative-parabolic.fsk")).unwrap() }
fn args(text:&str) -> Vec<String> {text.split_whitespace().map(String::from).collect()}

#[test]
fn reciprocal_audio_requires_explicit_work_and_preserves_tuning_force_indices() {
    let s=spec();let values="20 40 500 32 12 0.4 0.3 1.2";
    let original=Settings::read(&args(&format!("--wav {values}")),&s).unwrap();
    let reciprocal=Settings::read(&args(&format!("--wav-reciprocal {values} 2000000000000")),&s).unwrap();
    assert_eq!(original.image,Image::Collocation);assert_eq!(original.maximum_work,100_000_000_000);
    assert_eq!(reciprocal.image,Image::Reciprocal);assert_eq!(reciprocal.maximum_work,2_000_000_000_000);
    assert_eq!(original.band_hz,reciprocal.band_hz);assert_eq!(original.receiver_m,reciprocal.receiver_m);
    assert_eq!(original.order,reciprocal.order);assert_eq!(original.intervals,reciprocal.intervals);
    assert_eq!(original.full_scale_pa,reciprocal.full_scale_pa);
    for text in [String::new(),format!("--wav-reciprocal {values}"),format!("--wav {values} 2000000000000"),
        format!("--wav-reciprocal {values} 0"),format!("--wav-reciprocal {values} 100000000000001"),
        format!("--wav-reciprocal {values} NaN"),format!("--wav-reciprocal {values} -1")] {
        assert!(Settings::read(&args(&text),&s).is_err());
    }
    let mut command=args(&format!("spec.fsk 7200 --wav-reciprocal {values} 2000000000000 --edge-tensioner 0.8 4000 20 0.002 --edge-force 38400 80 --edge-force 76800 0"));
    let program=crate::tuning::option(&mut command,s.maximum_slope).unwrap().unwrap();
    assert_eq!(command.len(),12);
    assert_eq!(Settings::read(&command[2..],&s).unwrap().image,Image::Reciprocal);
    assert_eq!(program.force_n(38399),0.0);assert_eq!(program.force_n(38400),80.0);
    assert_eq!(program.force_n(76799),80.0);assert_eq!(program.force_n(76800),0.0);
}

fn fixture() -> Boundary {
    let mut s=spec();s.head.radial_intervals=2;s.head.azimuths=8;
    s.profile=vec![Station {z_m:0.0,radius_m:0.0},Station {z_m:0.1,radius_m:0.24},Station {z_m:0.28,radius_m:s.head.radius_m}];
    s.exterior=vec![Station {z_m:-0.002,radius_m:0.0},Station {z_m:0.1,radius_m:0.25},Station {z_m:0.28,radius_m:0.34}];
    let film=TensionedDisk::new(s.head,ProfileBudget {max_nodes:64,max_triangles:128,max_feature_evaluations:0}).unwrap();
    let mut phi=vec![0.0;film.model.free];
    for (i,&(x,y)) in film.mesh.nodes.iter().enumerate() {
        if let Some(k)=film.model.dof_map[3*i] {phi[k]=1.0-(x*x+y*y)/s.head.radius_m.powi(2);}
    }
    let modes=[ModePair {lambda:1.0,phi,residual:0.0,interval:(1.0,1.0)}];
    Boundary::new(&s,&film,&modes).unwrap()
}

#[test]
fn supplied_outer_wall_uses_original_reciprocal_pressure_and_signed_source_work() {
    let boundary=fixture();let surface=SpherePanels::from_triangles(boundary.triangles).unwrap();
    let gate=CancelGate::new_clock_free();let medium=Medium {density:1.2,sound_speed:343.0};
    let w=TAU*200.0;let k=w/medium.sound_speed;
    let source=field::prepare(&surface,k,3,1,2_000_000_000,&gate).unwrap();
    let one:Vec<_>=boundary.weights[0].iter().map(|b|C64::new(0.0,*b/w)).collect();
    let negative:Vec<_>=one.iter().map(|v|v.scale(-2.0)).collect();let zero=vec![C64::ZERO;one.len()];
    let solved=source.solve_batch(k,medium,&[&one,&negative,&zero],||false).unwrap();
    let point=[0.4,0.3,1.2];let truth=solved[0].pressure_at(point).unwrap();
    assert!(truth.abs()>1e-12);assert!(solved[0].outgoing_power_w>0.0);
    assert!((solved[1].outgoing_power_w-4.0*solved[0].outgoing_power_w).abs()<1e-9*solved[0].outgoing_power_w);
    let fields:Vec<_>=solved.into_iter().map(Field::Reciprocal).collect();
    for (field,scale) in fields.iter().zip([1.0,-2.0,0.0]) {
        field.admit(surface.areas().len()).unwrap();
        let actual=field.pressure_at(&surface,medium,point).unwrap();
        assert!((actual-truth.scale(scale)).abs()<1e-9*truth.abs());
    }
    // Full source/reaction work retains a signed off-diagonal entry. This
    // projects exact panel means but does not repurpose them as a P0 observer.
    let rows=[boundary.weights[0].clone(),boundary.weights[0].iter().map(|v|-2.0*v).collect()];
    let mut z=vec![C64::ZERO;4];
    for j in 0..2 {for i in 0..2 {
        let force=fields[j].trace().iter().zip(surface.areas()).zip(&rows[i])
            .fold(C64::ZERO,|sum,((p,a),b)|sum+p.scale(a*b));
        z[2*i+j]=force*C64::new(0.0,-w);
    }}
    admit_power(&z,2).unwrap();assert!(z[0].re>0.0 && z[1].re<0.0);
    assert!((z[1]-z[2]).abs()<1e-9*z[0].abs());
    let legacy=solve_radiation_batch(&surface,k,medium,&[&one],Formulation::PlainCbie).unwrap().remove(0);
    let expected=exterior_pressure_at_points_integrated(&surface,&legacy,medium,&[point]).unwrap()[0];
    let field=Field::Collocation(legacy);
    assert_eq!(field.pressure_at(&surface,medium,point).unwrap(),expected);
}

#[test]
fn reciprocal_budget_and_resonance_refuse_before_any_alternate_solve() {
    let b=fixture();let surface=SpherePanels::from_triangles(b.triangles).unwrap();let gate=CancelGate::new_clock_free();
    assert!(field::prepare(&surface,1.0,1,65,1,&gate).is_err());
    assert!(field::prepare(&surface,1.0,1,0,2_000_000_000,&gate).is_err());
    assert!(field::prepare(&surface,1.0,33,1,2_000_000_000,&gate).is_err());
    // The supplied depth is .282 m: this is beyond the slab exclusion.
    assert!(field::prepare(&surface,12.0,1,1,2_000_000_000,&gate).is_err());
    gate.request();assert!(field::prepare(&surface,1.0,1,1,2_000_000_000,&gate).is_err());
}
