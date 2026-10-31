//! Manufactured analytic and simultaneous-equation controls for mouth flow.
use super::*;
use fs_material::gas::GasSpec;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn cylinders(radii: &[f64]) -> Duct {
    Duct { segments: radii.iter().map(|&radius|
        Segment::Cylinder { radius, length: 0.03 }).collect() }
}
fn close(actual: C64, expected: C64, scale: f64) {
    assert!((actual-expected).abs() <= 5e-9 * scale.max(expected.abs()).max(1e-30),
        "actual {actual:?}, expected {expected:?}, scale {scale:e}");
}

#[test]
fn resistive_cylinder_has_analytic_phase_gain_and_decaying_higher_modes() {
    let gas = air();
    let radius = 0.006;
    let length = 0.03;
    let duct = cylinders(&[radius]);
    let zc = gas.density*gas.sound_speed/(core::f64::consts::PI*radius*radius);
    for modes in [1, 2, 4] {
        for frequency in [113.0, 471.0, 1307.0, 7003.0, 14999.0] {
            let omega = core::f64::consts::TAU*frequency;
            let zload = C64::from_re(0.7*zc);
            let response = mouth_response_load(&duct, &gas, omega, LossModel::Lossless,
                zload, modes, 1).unwrap();
            let original = mm_input_impedance_load(&duct, &gas, omega, LossModel::Lossless,
                zload, modes, 1).unwrap();
            assert_eq!(response.input.impedance_matrix, original.impedance_matrix);
            let outgoing = response.per_incident_plane_pressure(zc).unwrap();
            let phase = omega*length/gas.sound_speed;
            close(outgoing[0], C64::new(phase.cos(),phase.sin()).scale(2.0/(1.7*zc)), 1.0/zc);
            for &mode in &outgoing[1..] { assert_eq!(mode, C64::ZERO); }
            let wave = modal_wave(&gas,radius,omega,LossModel::Lossless,&j1_roots(modes-1).unwrap()).unwrap();
            for mode in 1..modes {
                close(response.flow_map[mode*modes+mode],crate::exp_i(wave.k[mode].scale(length)),1.0);
            }
        }
    }
    let long = Duct { segments: vec![Segment::Cylinder { radius: 0.002, length: 2.0 }] };
    let response = mouth_response_load(&long,&gas,100.0,LossModel::Lossless,C64::from_re(zc),4,1).unwrap();
    assert!(response.flow_map.iter().all(|&z| finite(z)));
    for mode in 1..4 { assert_eq!(response.flow_map[mode*4+mode],C64::ZERO); }
}

#[test]
fn one_mode_outlet_reproduces_the_scalar_chain_in_both_junction_directions() {
    let gas = air();
    for radii in [&[0.006,0.009,0.012][..], &[0.012,0.009,0.006][..]] {
        let duct = cylinders(radii);
        for frequency in [113.0,471.0,1307.0] {
            let omega = core::f64::consts::TAU*frequency;
            let scalar = crate::input_impedance(&duct,&gas,omega,LossModel::Bessel,Termination::IdealOpen).unwrap();
            let response = mouth_response(&duct,&gas,omega,LossModel::Bessel,Termination::IdealOpen,1,1).unwrap();
            close(response.flow_map[0],scalar.u_mouth_over_u_in,1.0);
            close(response.input.plane_impedance,scalar.impedance,scalar.impedance.abs());
        }
    }
}

// Independently solve ALL forward-left/backward-right waves at once. No
// recursive impedance conversion or production transfer helper appears here.
fn direct(duct: &Duct, gas: &GasState, omega: f64, loss: LossModel, n: usize)
    -> (Vec<C64>,Vec<C64>)
{
    let stations = stations_of(duct,1).unwrap();
    let gammas = j1_roots(n-1).unwrap();
    let waves: Vec<_> = stations.iter().map(|s|
        modal_wave(gas,s.radius,omega,loss,&gammas).unwrap()).collect();
    let dimension = 2*n*stations.len();
    let endpoint = |station: usize, left: bool| {
        let mut p = vec![C64::ZERO;n*dimension];
        let mut u = p.clone();
        for mode in 0..n {
            let e = crate::exp_i(waves[station].k[mode].scale(stations[station].length));
            let a = station*2*n+mode;
            let b = a+n;
            p[mode*dimension+a] = if left { C64::ONE } else { e };
            p[mode*dimension+b] = if left { e } else { C64::ONE };
            u[mode*dimension+a] = p[mode*dimension+a]/waves[station].zc[mode];
            u[mode*dimension+b] = -p[mode*dimension+b]/waves[station].zc[mode];
        }
        (p,u)
    };
    let (pin,uin) = endpoint(0,true);
    let (pout,uout) = endpoint(stations.len()-1,false);
    let mut a = vec![C64::ZERO;dimension*dimension];
    a[..n*dimension].copy_from_slice(&uin);
    let mut row = n;
    for station in 0..stations.len()-1 {
        let (p,u) = endpoint(station,false);
        let (q,v) = endpoint(station+1,true);
        let ru = stations[station].radius;
        let rd = stations[station+1].radius;
        let f = junction_projection(ru.min(rd),ru.max(rd),&gammas,n);
        for mode in 0..n {
            for col in 0..dimension {
                let mut pressure = if ru < rd { p[mode*dimension+col] } else { q[mode*dimension+col] };
                let mut flow = if ru < rd { v[mode*dimension+col] } else { u[mode*dimension+col] };
                for k in 0..n {
                    pressure = pressure - f[mode*n+k]*if ru < rd { q[k*dimension+col] } else { p[k*dimension+col] };
                    flow = flow - f[k*n+mode]*if ru < rd { u[k*dimension+col] } else { v[k*dimension+col] };
                }
                a[(row+mode)*dimension+col] = pressure;
                a[(row+n+mode)*dimension+col] = flow;
            }
        }
        row += 2*n;
    }
    let mut zload = waves.last().unwrap().zc.clone();
    // The scalar resistance is independent of modal loss; higher modes matched.
    let radius = stations.last().unwrap().radius;
    zload[0] = C64::from_re(0.7*gas.density*gas.sound_speed/(core::f64::consts::PI*radius*radius));
    for mode in 0..n {
        for col in 0..dimension {
            a[(row+mode)*dimension+col] = pout[mode*dimension+col]-zload[mode]*uout[mode*dimension+col];
        }
    }
    let scales: Vec<_> = a.chunks(dimension).map(|r|
        r.iter().map(|v|v.abs()).fold(0.0_f64,f64::max)).collect();
    for (r,&scale) in a.chunks_mut(dimension).zip(&scales) {
        for value in r { *value = value.scale(1.0/scale); }
    }
    let lu = lu_complex(&a,dimension).unwrap();
    let mut zin = vec![C64::ZERO;n*n];
    let mut transfer = zin.clone();
    for col in 0..n {
        let mut x = vec![C64::ZERO;dimension];
        x[col] = C64::from_re(1.0/scales[col]);
        lu.solve(&mut x);
        for mode in 0..n {
            for k in 0..dimension {
                zin[mode*n+col] = zin[mode*n+col]+pin[mode*dimension+k]*x[k];
                transfer[mode*n+col] = transfer[mode*n+col]+uout[mode*dimension+k]*x[k];
            }
        }
    }
    (zin,transfer)
}

#[test]
fn complete_flow_map_matches_simultaneous_modal_junction_equations() {
    let gas = air();
    for radii in [&[0.006,0.009][..],&[0.009,0.006][..],&[0.006,0.009,0.007,0.012][..]] {
        let duct = cylinders(radii);
        let r = *radii.last().unwrap();
        let load = C64::from_re(0.7*gas.density*gas.sound_speed/(core::f64::consts::PI*r*r));
        for n in [1,2,4] {
            for frequency in [113.0,1307.0,7003.0,14999.0] {
                for loss in [LossModel::Lossless,LossModel::Bessel] {
                    let omega = core::f64::consts::TAU*frequency;
                    let actual = mouth_response_load(&duct,&gas,omega,loss,load,n,1).unwrap();
                    let original = mm_input_impedance_load(&duct,&gas,omega,loss,load,n,1).unwrap();
                    assert_eq!(actual.input.impedance_matrix,original.impedance_matrix);
                    let (z,t) = direct(&duct,&gas,omega,loss,n);
                    let zscale = z.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
                    let tscale = t.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
                    for (&a,&b) in actual.input.impedance_matrix.iter().zip(&z) { close(a,b,zscale); }
                    for (&a,&b) in actual.flow_map.iter().zip(&t) { close(a,b,tscale); }
                }
            }
        }
    }
}

#[test]
fn coupled_flanged_load_keeps_input_response_and_all_modal_radiated_power() {
    let gas = air();
    let duct = cylinders(&[0.006,0.012]);
    let n = 3;
    for frequency in [471.0,7307.0] {
        let omega = core::f64::consts::TAU*frequency;
        let result = mouth_response_flanged(&duct,&gas,&gas,omega,LossModel::Lossless,n,1).unwrap();
        let original = mm_input_impedance_flanged(&duct,&gas,&gas,omega,LossModel::Lossless,n,1).unwrap();
        assert_eq!(result.input.impedance_matrix,original.impedance_matrix);
        let zload = crate::flanged::flanged_modal_impedance(&gas,0.012,omega,n).unwrap();
        let u: Vec<_> = (0..n).map(|i|result.flow_map[i*n]).collect();
        let mut power = 0.0;
        for i in 0..n {
            let mut p = C64::ZERO;
            for j in 0..n { p = p + zload[i*n+j]*u[j]; }
            power += p.re*u[i].re+p.im*u[i].im;
        }
        let input = result.input.plane_impedance.re;
        assert!((input-power).abs() < 5e-8*input.abs().max(1.0));
    }
}

#[test]
fn strict_table_scope_and_explicit_dc_remain_distinct() {
    let gas = air();
    let duct = cylinders(&[0.006]);
    let table = crate::TabulatedLoad::try_new(vec![10.0,1e5],vec![C64::from_re(2e6);2],"control",0.0)
        .unwrap().with_mouth(crate::TabulatedMouth { axial_position_m: 0.03,
            radius_m: 0.006, gas, dc_impedance: 2e6 }).unwrap();
    let result = mouth_response_tabulated(&duct,&gas,100.0,LossModel::Lossless,&table,2,1).unwrap();
    let original = mm_input_impedance_tabulated(&duct,&gas,100.0,LossModel::Lossless,&table,2,1).unwrap();
    assert_eq!(result.input.impedance_matrix,original.impedance_matrix);
    assert!(mouth_response_tabulated(&duct,&gas,1e-3,LossModel::Lossless,&table,2,1).is_err());
    assert!(mouth_response_load(&duct,&gas,1e-3,LossModel::Lossless,C64::from_re(2e6),2,1).is_ok());
    let other = cylinders(&[0.007]);
    assert!(mouth_response_tabulated(&other,&gas,100.0,LossModel::Lossless,&table,2,1).is_err());
    for reference in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(result.per_incident_plane_pressure(reference).is_err());
    }
    let closed = mouth_response(&duct,&gas,100.0,LossModel::Lossless,Termination::Closed,2,1).unwrap();
    assert!(closed.flow_map.iter().all(|v|v.abs()==0.0));
}

#[test]
fn rigid_cap_blocks_all_modal_mouth_flow_after_junctions() {
    // G1/G3: no mode may leak through a rigid cap, including propagating
    // higher modes and evanescent fields coupled by a change of radius.
    let gas = air();
    for radii in [&[0.05][..], &[0.035, 0.05][..], &[0.05, 0.035][..]] {
        let duct = cylinders(radii);
        for frequency in [300.0, 5000.0, 10000.0] {
            let closed = mouth_response(
                &duct, &gas, core::f64::consts::TAU * frequency,
                LossModel::Lossless, Termination::Closed, 3, 1,
            ).unwrap();
            assert!(
                closed.flow_map.iter().all(|&flow| flow == C64::ZERO),
                "rigid cap leaked at {frequency} Hz through {radii:?}: {:?}",
                closed.flow_map,
            );
        }
    }
}
