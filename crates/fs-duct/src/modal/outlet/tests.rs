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
fn direct(duct: &Duct, gas: &GasState, omega: f64, loss: LossModel, n: usize,
    ideal: Option<Termination>)
    -> (Vec<C64>,Vec<C64>)
{
    direct_source(duct, gas, omega, loss, n, ideal, None)
}

// `Some((references, plane_only))` prescribes real-reference incident waves
// instead of inlet flows. The global physical equations then remain finite
// when the inlet itself is a pressure antinode.
fn direct_source(duct: &Duct, gas: &GasState, omega: f64, loss: LossModel, n: usize,
    ideal: Option<Termination>, source: Option<(&[f64], bool)>)
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
    if let Some((references, plane_only)) = source {
        for mode in 0..n {
            if plane_only && mode > 0 { continue; }
            for col in 0..dimension {
                a[mode*dimension+col] = pin[mode*dimension+col]
                    + uin[mode*dimension+col].scale(references[mode]);
            }
        }
    }
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
            a[(row+mode)*dimension+col] = match ideal {
                Some(Termination::Closed) => uout[mode*dimension+col],
                Some(Termination::IdealOpen) => pout[mode*dimension+col],
                None => pout[mode*dimension+col]-zload[mode]*uout[mode*dimension+col],
                _ => unreachable!("only ideal boundaries are used by this oracle"),
            };
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
        if !source.is_some_and(|(_, plane_only)| plane_only && col > 0) {
            x[col] = C64::from_re(if source.is_some() { 2.0/scales[col] } else { 1.0/scales[col] });
        }
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
fn native_quarter_wave_mouth_flow_survives_an_exact_input_impedance_pole() {
    // Exact pressure-release quarter-wave cylinder with native Zc=2:
    // R_in=+1, H=2 exp(i*pi/2)/Zc=i. Input flow vanishes while the
    // mouth still moves; converting either Z_in or U_mouth/U_in must fail.
    let mut waves = MouthWaveResponse {
        input: ModalReflectionResponse { omega: 1.0, n_modes: 1,
            reflection_matrix: vec![C64::ONE], characteristic_impedances: vec![C64::from_re(2.0)],
            min_shear_number: f64::INFINITY, mouth_ka: 0.1,
            cutoffs_at_throat_hz: vec![], steps_used: 1 },
        inlet_radius_m: 0.01, mouth_radius_m: 0.01,
        native_flow_map: vec![C64::new(0.0, 1.0)],
    };
    assert!(matches!(impedance_from_reflection(&waves.input.reflection_matrix,
        &waves.input.characteristic_impedances, 1), Err(DuctError::Singular)));
    for reference in [0.5, 2.0, 7.0] {
        assert_eq!(waves.per_incident_plane_pressure(reference).unwrap(), vec![C64::new(0.0, 1.0)]);
        let incident = C64::new(0.4, -0.7);
        let result = waves.excite(&[reference], &[incident]).unwrap();
        close(result.reflected_pressure[0], incident, 1.0);
        close(result.mouth_flow[0], C64::new(0.0, 1.0)*incident, 1.0);
    }
    assert!(waves.excite(&[2.0], &[]).is_err());
    assert!(waves.excite(&[0.0], &[C64::ONE]).is_err());
    assert!(waves.excite(&[2.0], &[C64::new(f64::NAN, 0.0)]).is_err());
    waves.native_flow_map = vec![];
    assert!(waves.per_incident_plane_pressure(2.0).is_err());
}

#[test]
fn native_mouth_source_closures_match_independent_global_wave_equations() {
    let gas = air();
    let n = 3;
    let incident = [C64::new(1.0, 0.2), C64::new(-0.3, 0.4), C64::new(0.1, -0.5)];
    for radii in [[0.006, 0.009, 0.007], [0.012, 0.008, 0.010]] {
        let duct = cylinders(&radii);
        let geometry = PreparedModalDuct::new(&duct, n, 1).unwrap();
        let z0 = gas.density*gas.sound_speed/(core::f64::consts::PI*radii[0]*radii[0]);
        let load = C64::from_re(0.7*gas.density*gas.sound_speed/(core::f64::consts::PI*radii[2]*radii[2]));
        let references: Vec<_> = (0..n).map(|i| z0*(1.0+0.35*i as f64)).collect();
        for omega in [2400.0, 126_000.0] {
            for loss in [LossModel::Lossless, LossModel::Bessel] {
                let waves = geometry.mouth_waves_load(&gas, omega, loss, load).unwrap();
                let actual = waves.excite(&references, &incident).unwrap();
                let (pressure, flow) = direct_source(&duct, &gas, omega, loss, n, None,
                    Some((&references, false)));
                for row in 0..n {
                    let p = (0..n).fold(C64::ZERO, |sum, col| sum+pressure[row*n+col]*incident[col]);
                    let u = (0..n).fold(C64::ZERO, |sum, col| sum+flow[row*n+col]*incident[col]);
                    close(actual.reflected_pressure[row], p-incident[row], 1.0);
                    close(actual.mouth_flow[row], u, 1.0/z0);
                }
                let (_, flow) = direct_source(&duct, &gas, omega, loss, n, None,
                    Some((&references, true)));
                for (row, u) in waves.per_incident_plane_pressure(z0).unwrap().iter().enumerate() {
                    close(*u, flow[row*n], 1.0/z0);
                }
            }
        }
    }
}

#[test]
fn coherent_flanged_excitation_preserves_mouth_phase_and_radiated_work() {
    let gas = air();
    let exterior = GasState::try_new_moist_air(303.15, gas.pressure, 0.7).unwrap();
    let (radius, length, n) = (0.04, 0.03, 3);
    let duct = Duct { segments: vec![Segment::Cylinder { radius, length }] };
    let geometry = PreparedModalDuct::new(&duct, n, 1).unwrap();
    let z0 = gas.density*gas.sound_speed/(core::f64::consts::PI*radius*radius);
    let references = [z0, 0.7*z0, 1.4*z0];
    let incident = [C64::new(1.0, 0.2), C64::new(-0.3, 0.4), C64::new(0.1, -0.2)];
    for ka in [0.6, 4.6, 7.8] {
        let omega = ka*gas.sound_speed/radius;
        let waves = geometry.mouth_waves_flanged(&gas, &exterior, omega, LossModel::Lossless).unwrap();
        let result = waves.excite(&references, &incident).unwrap();
        let radiation = crate::flanged::flanged_modal_impedance(&exterior, radius, omega, n).unwrap();
        let load = ModalMouthImpedance::new(omega, radius, exterior, ModeCount(n), radiation.clone()).unwrap();
        let matrix_result = geometry.mouth_waves_matrix(&gas, omega, LossModel::Lossless, &load)
            .unwrap().excite(&references, &incident).unwrap();
        assert_eq!(matrix_result.reflected_pressure, result.reflected_pressure);
        assert_eq!(matrix_result.mouth_flow, result.mouth_flow);
        let wave = modal_wave(&gas, radius, omega, LossModel::Lossless, &j1_roots(n-1).unwrap()).unwrap();
        let mut radiation_power = 0.0;
        for row in 0..n {
            let u = result.mouth_flow[row];
            let p = (0..n).fold(C64::ZERO, |sum, col| sum+radiation[row*n+col]*result.mouth_flow[col]);
            radiation_power += 0.5*(p*u.conj()).re;
            // Independent physical pressure/flow propagation back to the inlet.
            let forward = crate::exp_i(wave.k[row].scale(length));
            let backward = crate::exp_i(wave.k[row].scale(-length));
            let cosine = (forward+backward).scale(0.5);
            let minus_i_sine = (backward-forward).scale(0.5);
            let pin = cosine*p + minus_i_sine*wave.zc[row]*u;
            let uin = minus_i_sine*p/wave.zc[row] + cosine*u;
            close(pin, incident[row]+result.reflected_pressure[row], 1.0);
            close(uin, (incident[row]-result.reflected_pressure[row]).scale(1.0/references[row]), 1.0/z0);
        }
        let input_power: f64 = incident.iter().zip(&result.reflected_pressure).zip(references)
            .map(|((&a,&b),z)| 0.5*(a.abs().powi(2)-b.abs().powi(2))/z).sum();
        assert!(radiation_power > 0.0);
        assert!((input_power-radiation_power).abs() < 5e-9*radiation_power);
    }
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
                    let (z,t) = direct(&duct,&gas,omega,loss,n,None);
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
fn modal_junctions_cross_local_pressure_and_flow_nodes() {
    // A quarter/half-wave downstream section puts a zero or pole of its
    // impedance exactly at the radius junction. The complete duct response
    // remains finite and must agree with a simultaneous pressure/flow solve.
    let gas = air();
    let downstream_length = 0.04;
    for radii in [[0.006,0.012], [0.012,0.006]] {
        let duct = Duct { segments: vec![
            Segment::Cylinder { radius: radii[0], length: 0.03 },
            Segment::Cylinder { radius: radii[1], length: downstream_length },
        ] };
        for n in [1,2,4] {
            for termination in [Termination::Closed, Termination::IdealOpen] {
                for phase in [core::f64::consts::FRAC_PI_2, core::f64::consts::PI] {
                    for detuning in [-1e-10, 0.0, 1e-10] {
                        let omega = gas.sound_speed*(phase+detuning)/downstream_length;
                        let actual = mouth_response(&duct,&gas,omega,LossModel::Lossless,
                            termination,n,1).unwrap();
                        let original = mm_input_impedance(&duct,&gas,omega,LossModel::Lossless,
                            termination,n,1).unwrap();
                        assert_eq!(actual.input.impedance_matrix,original.impedance_matrix);
                        let (z,t) = direct(&duct,&gas,omega,LossModel::Lossless,n,Some(termination));
                        let zscale = z.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
                        let tscale = t.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
                        for (&a,&b) in actual.input.impedance_matrix.iter().zip(&z) { close(a,b,zscale); }
                        for (&a,&b) in actual.flow_map.iter().zip(&t) { close(a,b,tscale.max(1.0)); }
                        if termination == Termination::Closed {
                            assert!(actual.flow_map.iter().all(|&v|v==C64::ZERO));
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn mixed_modal_nodes_preserve_junction_continuity_and_lossless_power() {
    // Exact +1/-1 reflection eigenvalues cannot pass through either an
    // impedance chart or an admittance chart. They are valid standing-wave
    // modes, and must coexist with mode conversion at a real radius change.
    let n = 3;
    let gammas = j1_roots(n-1).unwrap();
    let f = junction_projection(0.006,0.009,&gammas,n);
    let mut downstream_r = vec![C64::ZERO;n*n];
    downstream_r[0] = C64::ONE;
    downstream_r[n+1] = -C64::ONE;
    downstream_r[2*n+2] = C64::new(0.6,0.8);
    for expansion in [false,true] {
        let geometry = prepared::PreparedJunction {
            expansion, projection: f.clone(), transpose: transpose(&f,n),
        };
        for evanescent in [false,true] {
            let mut down = vec![C64::from_re(1.2e6),C64::from_re(2.1e6),C64::from_re(3.4e6)];
            let mut up = vec![C64::from_re(1.7e6),C64::from_re(2.8e6),C64::from_re(4.6e6)];
            if evanescent {
                down[1] = C64::new(0.0,-2.1e6);
                up[1] = C64::new(0.0,-2.8e6);
            }
            assert!(impedance_from_reflection(&downstream_r,&down,n).is_err());
            // Retaining identity as an observation extracts a_d = J a_u
            // from the SAME factorization used for the returned reflection.
            let mut transfer = Some(identity(n));
            let reflected = junction(&mut transfer,&downstream_r,&down,&up,&geometry,n).unwrap();
            let j = transfer.unwrap();
            let no_observation = junction(&mut None,&downstream_r,&down,&up,&geometry,n).unwrap();
            assert_eq!(reflected,no_observation);
            let mut p_down = downstream_r.clone();
            let mut p_up = reflected.clone();
            for mode in 0..n {
                p_down[mode*n+mode] = p_down[mode*n+mode]+C64::ONE;
                p_up[mode*n+mode] = p_up[mode*n+mode]+C64::ONE;
            }
            let p_down = mat_mul(&p_down,&j,n);
            let u_down = mat_mul(&terminal_flow(&downstream_r,&down,n),&j,n);
            let u_up = terminal_flow(&reflected,&up,n);
            let (p_left,p_right,u_left,u_right) = if expansion {
                (p_up.clone(),mat_mul(&f,&p_down,n),u_down.clone(),mat_mul(&geometry.transpose,&u_up,n))
            } else {
                (p_down,mat_mul(&f,&p_up,n),u_up.clone(),mat_mul(&geometry.transpose,&u_down,n))
            };
            for (&a,&b) in p_left.iter().zip(&p_right) { close(a,b,1.0); }
            for (&a,&b) in u_left.iter().zip(&u_right) { close(a,b,1e-6); }
            if !evanescent {
                // R* diag(1/Zc) R = diag(1/Zc) includes cross-mode energy,
                // so per-entry clipping could not satisfy this invariant.
                for row in 0..n {
                    for col in 0..n {
                        let mut power = C64::ZERO;
                        for mode in 0..n {
                            power = power+reflected[mode*n+row].conj()*reflected[mode*n+col]/up[mode];
                        }
                        let expected = if row==col { up[row].recip() } else { C64::ZERO };
                        close(power,expected,1e-6);
                    }
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

#[test]
fn modal_junction_crosses_simultaneous_pressure_and_flow_nodes() {
    // G1: mixed pressure/flow nodes cannot make a physical junction singular.
    let gas = air();
    let n = 2;
    let r_down = [C64::ONE, C64::ZERO, C64::ZERO, -C64::ONE];
    for radii in [&[0.035, 0.05][..], &[0.05, 0.035][..]] {
        let source = cylinders(radii);
        let prepared = PreparedModalDuct::new(&source, n, 1).unwrap();
        let data = prepared.junctions[0].as_ref().unwrap();
        let omega = core::f64::consts::TAU * 8000.0;
        let z_up = modal_wave(&gas, radii[0], omega, LossModel::Lossless,
            &prepared.gammas).unwrap().zc;
        let z_down = modal_wave(&gas, radii[1], omega, LossModel::Lossless,
            &prepared.gammas).unwrap().zc;
        let mut transmission = Some(vec![C64::ONE, C64::ZERO, C64::ZERO, C64::ONE]);
        let r_up = junction(&mut transmission, &r_down, &z_down, &z_up, data, n).unwrap();
        assert!(transmission.unwrap().iter().chain(&r_up).all(|&v| finite(v)));
    }
}
