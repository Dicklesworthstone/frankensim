use super::*;

fn inputs(coincident:bool)->(Vec<Course>,Vec<BoardMode>) {
    let all=super::super::geometry::demonstration_scale().unwrap();
    let c=Course {unison:1,duplex_length_m:0.0,detune_cents:0.0,..all[48]};
    let courses=if coincident {vec![c,Course {midi:72,..c}]}else{vec![c]};
    let mut board=vec![BoardMode {frequency_hz:170.0,damping_ratio:0.012,bridge:[0.0;88],volume:0.1},
        BoardMode {frequency_hz:310.0,damping_ratio:0.018,bridge:[0.0;88],volume:-0.03}];
    board[0].bridge[48]=0.08;board[0].bridge[51]=0.06;
    board[1].bridge[48]=-0.04;board[1].bridge[51]=0.09;
    (courses,board)
}
fn model(coincident:bool,damping:bool)->BridgeResponse {
    let (c,b)=inputs(coincident);BridgeResponse::new(&c,&b,192_000,21_600.0,4,damping).unwrap()
}

// Polarize the original time bank's stored energy. This independently retains
// its endpoint inertia, loaded board, cross-potential and ALL string modes,
// rather than copying either Schur elimination or its new pole partition.
fn full(m:&BridgeResponse,hz:f64,key:u8,force:C64,z:Option<&[C64]>)->Vec<C64> {
    let n=m.bank.modes.len();let r=m.bank.board_count;let size=n+r;let w=TAU*hz;
    let zero=vec![0.0;size];let mut q=zero.clone();let mut energies=zero.clone();
    for i in 0..size {q[i]=1.0;energies[i]=m.bank.energy_at(&q,&zero);q[i]=0.0;}
    let mut a=vec![C64::ZERO;size*size];
    for i in 0..size {for j in i..size {
        let k=if i==j {2.0*energies[i]}else{
            q[i]=1.0;q[j]=1.0;let k=m.bank.energy_at(&q,&zero)-energies[i]-energies[j];q[i]=0.0;q[j]=0.0;k
        };
        a[i*size+j]=C64::new(k,0.0);a[j*size+i]=C64::new(k,0.0);
    }}
    for i in 0..size {a[i*size+i]=a[i*size+i]-C64::new(w*w,if i<n {w*m.string_c[i]}else{0.0});}
    for i in 0..r {for j in 0..r {
        let at=(n+i)*size+n+j;a[at]=a[at]+C64::new(0.0,-w*m.board_c[i*r+j]);
        if let Some(z)=z {a[at]=a[at]+C64::new(0.0,-w)*z[i*r+j];}
    }}
    let mut rhs=vec![C64::ZERO;size];
    for (v,g) in rhs[n..].iter_mut().zip(m.bridge_row(key).unwrap()) {*v=force.scale(*g);}
    lu_complex(&a,size).unwrap().solve(&mut rhs);rhs
}
fn compare(m:&BridgeResponse,hz:f64,key:u8,z:Option<&[C64]>)->Response {
    let force=C64::new(1.0,0.3);let expected=full(m,hz,key,force,z);
    let a=m.solve(hz,key,force,z).unwrap();
    let scale=expected.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
    for (got,want) in a.string_displacement.iter().chain(&a.board_displacement).zip(&expected) {
        assert!((*got-*want).abs()<1e-7*scale,"full physical pencil disagrees at {hz} Hz");
    }
    assert!(a.backward_error<1e-12);
    assert!(a.power_defect_w.abs()<1e-10+1e-7*a.input_w.abs());a
}

#[test]
fn exact_and_near_string_poles_solve_the_original_piano_without_artificial_loss() {
    let m=model(false,false);let hz=m.bank.modes[0].omega/TAU;
    let original_omegas=m.bank.modes.iter().map(|s|s.omega.to_bits()).collect::<Vec<_>>();
    let original_q=m.bank.q.clone();let original_v=m.bank.v.clone();
    let z=[C64::new(3.0,-TAU*hz*0.01),C64::new(0.5,-TAU*hz*0.002),
        C64::new(0.5,-TAU*hz*0.002),C64::new(2.0,-TAU*hz*0.015)];
    for radiation in [None,Some(z.as_slice())] {
        for relative in [0.0,-1e-14,1e-14,-1e-9,1e-9,-2e-7,2e-7] {
            let a=compare(&m,hz*(1.0+relative),69,radiation);
            if relative.abs()<1e-8 {assert_eq!(a.retained_string_poles,1);}
            assert_eq!(a.string_loss_w,0.0);assert_eq!(a.board_loss_w,0.0);
            assert!(a.string_displacement[0].abs()>0.0);
        }
    }
    let zero=m.solve(hz,69,C64::ZERO,Some(&z)).unwrap();
    assert_eq!(zero.retained_string_poles,1);assert_eq!(zero.input_w,0.0);assert_eq!(zero.backward_error,0.0);
    assert!(zero.string_displacement.iter().chain(&zero.board_displacement).all(|v|*v==C64::ZERO));
    assert_eq!(original_q,m.bank.q);assert_eq!(original_v,m.bank.v);
    assert_eq!(original_omegas,m.bank.modes.iter().map(|s|s.omega.to_bits()).collect::<Vec<_>>());
}

#[test]
fn coincident_partial_poles_keep_each_independent_bridge_coordinate_and_reciprocity() {
    let m=model(true,false);let hz=m.bank.modes[0].omega/TAU;
    let a=compare(&m,hz,69,None);let b=compare(&m,hz,72,None);
    assert_eq!(a.retained_string_poles,2);assert_eq!(b.retained_string_poles,2);
    for relative in [-1e-9,1e-9] {
        let nearby=compare(&m,hz*(1.0+relative),69,None);
        let scale=a.string_displacement.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
        for (x,y) in nearby.string_displacement.iter().zip(&a.string_displacement) {
            assert!((*x-*y).abs()<1e-5*scale);
        }
    }
    // Work-conjugate ports remain reciprocal. At the exact partial the bridge
    // may be an antiresonance, so avoid dividing by its near-zero velocity.
    assert!((a.bridge_velocity[1]-b.bridge_velocity[0]).abs()<1e-12);
}

#[test]
fn ordinary_damped_admittance_keeps_the_original_small_schur_image() {
    let m=model(true,true);
    for hz in [170.0,277.0,431.0,m.bank.modes[0].omega/TAU] {
        let a=compare(&m,hz,69,None);assert_eq!(a.retained_string_poles,0);
        assert!(a.input_w>0.0 && a.string_loss_w>0.0 && a.board_loss_w>0.0);
    }
}

#[test]
fn too_many_coincident_poles_refuse_without_deleting_voices_or_inventing_a_loss() {
    let (c,mut board)=inputs(false);board.truncate(1);board[0].bridge.fill(0.01);
    let courses=(21..86).map(|midi|Course {midi,unison:2,..c[0]}).collect::<Vec<_>>();
    let m=BridgeResponse::new(&courses,&board,192_000,21_600.0,1,false).unwrap();
    assert_eq!(m.bank.modes.len(),130);
    let hz=m.bank.modes[0].omega/TAU;
    let error=m.solve(hz,21,C64::ONE,None).unwrap_err();assert!(error.contains("128-coordinate"));
    assert_eq!(m.bank.modes.len(),130);assert!(m.string_c.iter().all(|v|*v==0.0));
    assert!(m.bank.q.iter().chain(&m.bank.v).all(|v|*v==0.0));
}
