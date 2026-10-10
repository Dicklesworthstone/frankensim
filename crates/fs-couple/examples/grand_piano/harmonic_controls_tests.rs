use super::*;

fn inputs()->(Vec<Course>,Vec<BoardMode>,Vec<Vec<f64>>) {
    let scale=super::super::geometry::demonstration_scale().unwrap();
    let courses=[48,51].map(|i|Course {unison:2,duplex_length_m:0.2,..scale[i]}).to_vec();
    let mut board=vec![BoardMode {frequency_hz:170.,damping_ratio:0.012,bridge:[0.;88],volume:0.1},
        BoardMode {frequency_hz:310.,damping_ratio:0.018,bridge:[0.;88],volume:-0.03}];
    board[0].bridge[48]=0.08;board[0].bridge[51]=0.06;
    board[1].bridge[48]=-0.04;board[1].bridge[51]=0.09;
    (courses,board,vec![vec![0.05,0.07],vec![-0.03,0.04]])
}
fn load(w:f64)->Vec<C64> {
    vec![C64::new(3.,-w*0.01),C64::new(0.5,-w*0.002),
        C64::new(0.5,-w*0.002),C64::new(2.,-w*0.015)]
}

// Independent complete mechanical equation: obtain stiffness by polarizing
// the ACTUAL time bank energy. This does not copy Schur elimination, endpoint
// stiffness assembly or string-pole bordering from the harmonic consumer.
fn full(m:&BridgeResponse,hz:f64,key:u8,z:Option<&[C64]>)->Vec<C64> {
    let n=m.bank.modes.len();let r=m.bank.board_count;let size=n+r;let w=TAU*hz;
    let zero=vec![0.;size];let mut q=zero.clone();let mut energy=zero.clone();
    for i in 0..size {q[i]=1.;energy[i]=m.bank.energy_at(&q,&zero);q[i]=0.;}
    let mut matrix=vec![C64::ZERO;size*size];
    for i in 0..size {for j in i..size {
        let k=if i==j {2.*energy[i]}else {
            q[i]=1.;q[j]=1.;let k=m.bank.energy_at(&q,&zero)-energy[i]-energy[j];q[i]=0.;q[j]=0.;k
        };
        matrix[i*size+j]=C64::new(k,0.);matrix[j*size+i]=C64::new(k,0.);
    }}
    for i in 0..size {
        matrix[i*size+i]=matrix[i*size+i]-C64::new(w*w,if i<n {w*m.string_c[i]}else{0.});
    }
    for i in 0..r {for j in 0..r {
        let at=(n+i)*size+n+j;matrix[at]=matrix[at]+C64::new(0.,-w*m.board_c[i*r+j]);
        if let Some(z)=z {matrix[at]=matrix[at]+C64::new(0.,-w)*z[i*r+j];}
    }}
    let mut rhs=vec![C64::ZERO;size];
    for (out,g) in rhs[n..].iter_mut().zip(m.bridge_row(key).unwrap()) {*out=C64::new(*g,0.);}
    lu_complex(&matrix,size).unwrap().solve(&mut rhs);rhs
}
fn compare(m:&BridgeResponse,hz:f64,key:u8,z:Option<&[C64]>)->Response {
    let expected=full(m,hz,key,z);let got=m.solve(hz,key,C64::ONE,z).unwrap();
    let scale=expected.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
    for (got,want) in got.string_displacement.iter().chain(&got.board_displacement).zip(&expected) {
        assert!((*got-*want).abs()<1e-7*scale,"selected physical bank disagrees with full equation at {hz} Hz");
    }
    assert!(got.backward_error<1e-12);
    assert!(got.power_defect_w.abs()<1e-10+1e-7*got.input_w.abs());got
}

#[test]
fn both_transverse_planes_unisons_and_duplexes_match_the_complete_energy_hessian() {
    let (courses,board,secondary)=inputs();
    for source in [false,true] {
        let m=BridgeResponse::new_with_string_damping(&courses,&board,192_000,21_600.,4,true,
            Some(&secondary),source).unwrap();
        assert_eq!(m.bank.strings.len(),16);assert_eq!(m.bank.contact_strings.len(),4);
        let z=load(TAU*277.);
        for radiation in [None,Some(z.as_slice())] {
            let a=compare(&m,277.,69,radiation);let b=compare(&m,277.,72,radiation);
            assert!((a.bridge_velocity[1]-b.bridge_velocity[0]).abs()<1e-10*a.bridge_velocity[1].abs());
            assert!(a.string_loss_w>0. && a.board_loss_w>0.);
            if radiation.is_some() {assert!(a.radiation_w>0.);}
            for port in &m.bank.strings {
                assert!(a.string_displacement[port.modes.clone()].iter().any(|q|q.abs()>0.),
                    "a retained unstruck direction/course/duplex lost its bridge response");
            }
        }
    }
}

#[test]
fn published_losses_match_the_actual_two_plane_time_bank_generator() {
    let (mut courses,board,mut secondary)=inputs();courses.truncate(1);secondary.truncate(1);
    courses[0].duplex_length_m=0.;let mut errors=Vec::new();
    for rate in [96_000,192_000] {
        let mut m=BridgeResponse::new_with_string_damping(&courses,&board,rate,21_600.,4,true,
            Some(&secondary),true).unwrap();
        let n=m.bank.modes.len();let r=m.bank.board_count;
        for i in 0..n {m.bank.v[i]=0.01*((i+1) as f64).sin();}
        m.bank.v[n]=0.006;m.bank.v[n+1]=-0.004;
        let mut expected:f64=(0..n).map(|i|m.string_c[i]*m.bank.v[i].powi(2)).sum();
        for i in 0..r {for j in 0..r {expected+=m.board_c[r*i+j]*m.bank.v[n+i]*m.bank.v[n+j];}}
        m.bank.predict();m.bank.finish(&vec![0.;m.bank.contact_strings.len()]);
        let measured=m.bank.last_modal_loss_j*f64::from(rate);
        errors.push((measured-expected).abs()/expected);
    }
    assert!(errors[1]<0.004,"published-loss continuous/time generator mismatch: {errors:?}");
    assert!(errors[1]<errors[0],"generator comparison must improve with rate: {errors:?}");
    let estimated=BridgeResponse::new_with_string_damping(&courses,&board,192_000,21_600.,4,true,
        Some(&secondary),false).unwrap();
    let source=BridgeResponse::new_with_string_damping(&courses,&board,192_000,21_600.,4,true,
        Some(&secondary),true).unwrap();
    let hz=source.bank.modes[0].omega/TAU;
    let a=estimated.solve(hz,69,C64::ONE,None).unwrap();let b=source.solve(hz,69,C64::ONE,None).unwrap();
    assert!((a.bridge_velocity[0]-b.bridge_velocity[0]).abs()>1e-6*b.bridge_velocity[0].abs());
    assert_eq!(estimated.bank.modes.iter().map(|m|m.omega).collect::<Vec<_>>(),
        source.bank.modes.iter().map(|m|m.omega).collect::<Vec<_>>());
}

#[test]
fn coincident_polarization_poles_keep_both_coordinates_and_lossless_power() {
    let (mut courses,board,mut secondary)=inputs();courses.truncate(1);secondary.truncate(1);
    courses[0].unison=1;courses[0].duplex_length_m=0.;courses[0].detune_cents=0.;
    let m=BridgeResponse::new_with_string_damping(&courses,&board,192_000,21_600.,4,false,
        Some(&secondary),true).unwrap();
    let hz=m.bank.modes[0].omega/TAU;let z=load(TAU*hz);
    assert!(m.string_c.iter().chain(&m.board_c).all(|c|*c==0.));
    for radiation in [None,Some(z.as_slice())] {
        for relative in [0.,-1e-9,1e-9] {
            let result=compare(&m,hz*(1.+relative),69,radiation);
            assert_eq!(result.retained_string_poles,2);
            assert_eq!(result.string_loss_w,0.);assert_eq!(result.board_loss_w,0.);
        }
    }
}

#[test]
fn zero_secondary_ports_preserve_the_original_primary_bridge_experiment() {
    let (mut courses,board,_)=inputs();courses.truncate(1);
    let original=BridgeResponse::new(&courses,&board,192_000,21_600.,4,true).unwrap();
    let selected=BridgeResponse::new_with_string_damping(&courses,&board,192_000,21_600.,4,true,
        Some(&[vec![0.;board.len()]]),false).unwrap();
    let z=load(TAU*277.);
    let a=original.solve(277.,69,C64::ONE,Some(&z)).unwrap();
    let b=selected.solve(277.,69,C64::ONE,Some(&z)).unwrap();
    assert_eq!(a.board_displacement,b.board_displacement);assert_eq!(a.bridge_velocity,b.bridge_velocity);
    assert_eq!(a.input_w,b.input_w);assert_eq!(a.string_loss_w,b.string_loss_w);
    assert_eq!(original.bridge_row(69).unwrap(),selected.bridge_row(69).unwrap());
    for (index,port) in original.bank.strings.iter().enumerate() {
        let added=&selected.bank.strings[index];
        assert_eq!(a.string_displacement[port.modes.clone()],b.string_displacement[added.modes.clone()]);
    }
    for port in selected.bank.strings.iter().filter(|s|s.polarization==1) {
        assert!(b.string_displacement[port.modes.clone()].iter().all(|q|*q==C64::ZERO));
    }
}
