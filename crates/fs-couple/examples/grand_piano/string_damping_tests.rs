//! Independent full-matrix and analytic checks of reciprocal intrinsic loss.
use super::*;
use fs_math::det;

fn fixture(duration: f64) -> (Vec<StringPort>, Vec<StringMode>, Vec<Range<usize>>, Prepared) {
    let strings=vec![
        StringPort {course:0,member:0,polarization:0,duplex:false,modes:0..2,
            bridge:vec![8.0,-3.0],hammer_lift:0.0,damper_lift:0.0,contact:None},
        StringPort {course:0,member:1,polarization:0,duplex:false,modes:2..3,
            bridge:vec![8.0,-3.0],hammer_lift:0.0,damper_lift:0.0,contact:None},
        StringPort {course:1,member:0,polarization:1,duplex:false,modes:3..5,
            bridge:vec![5.0,9.0],hammer_lift:0.0,damper_lift:0.0,contact:None},
    ];
    let modes=[(0.08,0),(-0.04,0),(0.07,1),(0.06,2),(-0.03,2)].into_iter()
        .map(|(beta,string)|StringMode {omega:1000.0,beta,a:1e6*beta,
            hammer_shape:0.0,damper_shape:0.0,string}).collect::<Vec<_>>();
    let groups=vec![0..2,2..3];
    let losses=[0.0,4.0,9.0,2.0,17.0].into_iter().enumerate()
        .map(|(k,strain_rate)|ModeLoss {mass_rate:0.6+0.1*k as f64,strain_rate}).collect();
    let p=Prepared::new(&strings,&modes,&groups,losses,vec![0.023,0.010,0.034],2,duration).unwrap();
    (strings,modes,groups,p)
}

// Assemble physical loss as outer products of displacement-port rows. This
// oracle contains no grouped Schur elimination or prepared response factors.
fn dense(strings: &[StringPort], modes: &[StringMode], p: &Prepared) -> Vec<f64> {
    let n=modes.len();let size=n+strings[0].bridge.len();let mut c=vec![0.0;size*size];
    for (si,s) in strings.iter().enumerate() {
        let mut endpoint=vec![0.0;size];endpoint[n..].copy_from_slice(&s.bridge);
        for i in 0..size {for j in 0..size {
            c[i*size+j]+=p.endpoint_drag[si]*endpoint[i]*endpoint[j];
        }}
        for k in s.modes.clone() {
            c[k*size+k]+=p.modes[k].mass_rate;
            let mut relative=vec![0.0;size];relative[k]=1.0;
            for (v,g) in relative[n..].iter_mut().zip(&s.bridge) {*v=-modes[k].beta*g;}
            for i in 0..size {for j in 0..size {
                c[i*size+j]+=p.modes[k].strain_rate*relative[i]*relative[j];
            }}
        }
    }
    c
}
fn quadratic(c: &[f64], v: &[f64]) -> f64 {
    let size=v.len();let mut value=0.0;
    for i in 0..size {for j in 0..size {value+=v[i]*c[i*size+j]*v[j];}}
    value
}

#[test]
fn reciprocal_flow_matches_full_midpoint_and_closes_physical_loss() {
    for dt in [1e-6,0.01,0.2,1.0] {
        let (strings,modes,groups,mut p)=fixture(dt);let size=modes.len()+2;
        let c=dense(&strings,&modes,&p);let mut matrix=identity(size);
        for (entry,loss) in matrix.iter_mut().zip(&c) {*entry+=0.5*dt*loss;}
        let inverse=inverse_spd(&matrix,size).unwrap();
        let initial:Vec<_>=(0..size).map(|k|0.1*det::sin(0.7*(k+1) as f64)).collect();
        let midpoint:Vec<f64>=(0..size).map(|i|(0..size).map(|j|
            inverse[i*size+j]*initial[j]).sum()).collect();
        let expected:Vec<_>=midpoint.iter().zip(&initial).map(|(m,v)|2.0*m-v).collect();
        let mut actual=initial.clone();let loss=p.apply(&strings,&modes,&groups,&mut actual);
        for (got,want) in actual.iter().zip(&expected) {assert!((got-want).abs()<2e-12);}
        let physical=dt*quadratic(&c,&midpoint);
        assert!(loss>0.0);
        assert!((loss-physical).abs()<1e-13+1e-10*physical);
        assert!((p.power(&strings,&modes,&initial)-quadratic(&c,&initial)).abs()<1e-13);
        let closed=loss+0.5*actual.iter().map(|v|v*v).sum::<f64>()
            -0.5*initial.iter().map(|v|v*v).sum::<f64>();
        assert!(closed.abs()<1e-15);
    }
}

#[test]
fn pure_chord_has_no_bending_loss_but_keeps_all_distributed_drag() {
    let mass=0.026;let r=2;
    let strings=vec![StringPort {course:0,member:0,polarization:0,duplex:false,modes:0..4,
        bridge:vec![5.0,-3.0],hammer_lift:0.0,damper_lift:0.0,contact:None}];
    let modes=(1..=4).map(|n| {
        let sign=if n%2==0 {-1.0}else{1.0};
        let beta=det::sqrt(2.0*mass)*sign/(n as f64*std::f64::consts::PI);
        StringMode {omega:1000.0*n as f64,beta,a:0.0,hammer_shape:0.0,damper_shape:0.0,string:0}
    }).collect::<Vec<_>>();
    let residual=mass/3.0-modes.iter().map(|m|m.beta*m.beta).sum::<f64>();
    assert!(residual>0.0);
    let groups=vec![0..1];let board=[0.04,-0.02];
    let bridge=strings[0].bridge.iter().zip(board).map(|(g,v)|g*v).sum::<f64>();
    let initial:Vec<_>=modes.iter().map(|m|m.beta*bridge).chain(board).collect();
    for a in [0.0,0.6] {
        let losses=modes.iter().map(|m|ModeLoss {mass_rate:a,strain_rate:1e-5*m.omega.powi(2)}).collect();
        let mut p=Prepared::new(&strings,&modes,&groups,losses,vec![a*residual],r,0.01).unwrap();
        let expected=a*mass/3.0*bridge*bridge;
        assert!((p.power(&strings,&modes,&initial)-expected).abs()<1e-16);
        if a==0.0 {
            let mut actual=initial.clone();let loss=p.apply(&strings,&modes,&groups,&mut actual);
            for (got,want) in actual.iter().zip(&initial) {assert!((got-want).abs()<2e-14);}
            assert!(loss.abs()<1e-15);
        }
    }
    // Even with no retained partials, the moving string's residual endpoint
    // mass still loses energy; it must not disappear with modal truncation.
    let empty=vec![StringPort {modes:0..0,..strings[0].clone()}];
    let drag=0.6*mass/3.0;let dt=0.03;
    let mut p=Prepared::new(&empty,&[],&groups,vec![],vec![drag],r,dt).unwrap();
    let mut actual=board;let loss=p.apply(&empty,&[],&groups,&mut actual);
    let norm2=empty[0].bridge.iter().map(|g|g*g).sum::<f64>();
    let scale=dt*drag*bridge/(1.0+0.5*dt*drag*norm2);
    for j in 0..r {assert!((actual[j]-(board[j]-scale*empty[0].bridge[j])).abs()<1e-14);}
    assert!(loss>0.0);
}

#[test]
fn reciprocal_flow_converges_quadratically_to_the_full_damping_exponential() {
    let (strings,modes,groups,p)=fixture(0.1);let size=modes.len()+2;
    let c=dense(&strings,&modes,&p);let eigen=fs_modal::eigh_gen_dense(&c,&identity(size),size).unwrap();
    let initial:Vec<_>=(0..size).map(|k|0.1*det::cos(0.3*k as f64)).collect();
    let time=0.2;let mut exact=vec![0.0;size];
    for e in eigen {
        assert!(e.lambda>0.0);
        let amplitude=e.phi.iter().zip(&initial).map(|(p,v)|p*v).sum::<f64>()*det::exp(-e.lambda*time);
        for (v,p) in exact.iter_mut().zip(e.phi) {*v+=amplitude*p;}
    }
    let mut previous=f64::INFINITY;
    for steps in [8,16,32] {
        let mut p=Prepared::new(&strings,&modes,&groups,p.modes.clone(),
            p.endpoint_drag.clone(),2,time/f64::from(steps)).unwrap();
        let mut actual=initial.clone();
        for _ in 0..steps {p.apply(&strings,&modes,&groups,&mut actual);}
        let error=det::sqrt(actual.iter().zip(&exact).map(|(a,b)|(a-b).powi(2)).sum());
        assert!(error>1e-10 && error<previous/3.5,"{steps} steps: {error:e}, previous {previous:e}");
        previous=error;
    }
}
