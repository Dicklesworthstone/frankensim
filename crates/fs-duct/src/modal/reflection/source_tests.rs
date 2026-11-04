//! G1/G3: physical source/load oracles, coherent work and original modal closures.
use super::*;

fn response(reflection: Vec<C64>, impedances: Vec<C64>) -> ModalReflectionResponse {
    let n = impedances.len();
    ModalReflectionResponse {
        omega: 1000.0,
        n_modes: n,
        reflection_matrix: reflection,
        characteristic_impedances: impedances,
        min_shear_number: f64::INFINITY,
        mouth_ka: 0.1,
        cutoffs_at_throat_hz: (1..n).map(|i| 100.0*i as f64).collect(),
        steps_used: 1,
    }
}

fn pressure(pressure: C64, impedance: C64) -> PortExcitation {
    PortExcitation::Pressure {
        open_circuit_pressure_pa: pressure,
        impedance_pa_s_m3: impedance,
    }
}

fn close(actual: C64, expected: C64) {
    // These manufactured fields have unit-scale drives. In particular, a
    // closed-mode zero carries the cancellation error of those finite waves.
    assert!((actual-expected).abs() <= 2e-13*expected.abs().max(1.0),
        "actual={actual:?}, expected={expected:?}");
}

fn balanced(result: &ModalSourceResponse) {
    let scale = result.generator_powers_w.iter().map(|v| v.abs()).sum::<f64>()
        + result.source_dissipation_w.abs() + result.inlet_power_w.abs();
    assert!((result.generator_power_w-result.source_dissipation_w-result.inlet_power_w).abs()
        <= 2e-13*scale.max(1e-30));
    assert!(result.source_dissipations_w.iter().all(|&v| v >= 0.0));
}

fn transformer() -> ModalReflectionResponse {
    // Lossless transformer: P0=2 P1 and U1=-2 U0, with both ±1 reflection
    // eigenvalues. A full impedance or admittance matrix does not exist.
    response(vec![C64::ZERO, C64::from_re(2.0), C64::from_re(0.5), C64::ZERO],
        vec![C64::from_re(4.0), C64::ONE])
}

#[test]
fn complex_source_impedance_matches_the_independent_series_load_solution() {
    let load = C64::new(2.5, -1.2);
    let drive = C64::new(0.3, -1.1);
    // Native real, lossy and evanescent wave bases all describe the same
    // physical pressure/flow circuit; |R_native| is not a power criterion.
    for zc in [C64::from_re(4.0), C64::new(2.0, -0.2), C64::new(0.0, -0.4)] {
        let waves = response(vec![(load-zc)/(load+zc)], vec![zc]);
        for source in [C64::ZERO, C64::new(0.7, 1.5), C64::new(0.0, -2.0)] {
            let result = waves.drive_sources(&[pressure(drive, source)]).unwrap();
            let flow = drive/(load+source);
            close(result.inlet_flow[0], flow);
            close(result.inlet_pressure[0], load*flow);
            close(C64::from_re(result.inlet_power_w),
                C64::from_re(load.re*flow.norm_sq()*0.5));
            close(C64::from_re(result.source_dissipation_w),
                C64::from_re(source.re*flow.norm_sq()*0.5));
            balanced(&result);
        }
        let flow = C64::new(0.04, -0.03);
        let result = waves.drive_sources(&[PortExcitation::VolumeFlow(flow)]).unwrap();
        close(result.inlet_flow[0], flow);
        close(result.inlet_pressure[0], load*flow);
        assert_eq!(result.source_dissipation_w, 0.0);
        balanced(&result);
    }
}

#[test]
fn physical_sources_reproduce_plane_and_all_matched_wave_closures() {
    let waves = transformer();
    for references in [[4.0, 1.0], [0.5, 3.0], [20.0, 7.0]] {
        let incident = [C64::new(0.8, 0.3), C64::new(-0.2, 0.7)];
        let sources = core::array::from_fn::<_, 2, _>(|i|
            pressure(incident[i].scale(2.0), C64::from_re(references[i])));
        let (expected, native) = waves.excitation(&references, &incident).unwrap();
        let (actual, source_native) = waves.source_excitation(&sources).unwrap();
        for i in 0..2 {
            close(actual.inlet_pressure[i], incident[i]+expected[i]);
            close(actual.inlet_flow[i], (incident[i]-expected[i]).scale(1.0/references[i]));
            close(source_native[i], native[i]);
        }
        balanced(&actual);
        let (reflected, native) = waves.plane_excitation(references[0]).unwrap();
        let (actual, source_native) = waves.source_excitation(&[
            pressure(C64::from_re(2.0), C64::from_re(references[0])),
            PortExcitation::VolumeFlow(C64::ZERO),
        ]).unwrap();
        close(actual.inlet_pressure[0], C64::ONE+reflected);
        close(actual.inlet_flow[0], (C64::ONE-reflected).scale(1.0/references[0]));
        for i in 0..2 { close(source_native[i], native[i]); }
        close(actual.inlet_flow[1], C64::ZERO);
    }
}

#[test]
fn mixed_higher_mode_sources_retain_signed_generator_absorption() {
    let result = transformer().drive_sources(&[
        PortExcitation::VolumeFlow(C64::ONE),
        pressure(C64::ONE, C64::from_re(0.5)),
    ]).unwrap();
    // U1=-2 and P1=1-0.5*(-2)=2; P0=4. The flow generator supplies
    // 2 W, the pressure generator absorbs 1 W, and its resistance loses 1 W.
    close(result.inlet_flow[0], C64::ONE);
    close(result.inlet_flow[1], C64::from_re(-2.0));
    close(result.inlet_pressure[0], C64::from_re(4.0));
    close(result.inlet_pressure[1], C64::from_re(2.0));
    close(C64::from_re(result.generator_powers_w[0]), C64::from_re(2.0));
    close(C64::from_re(result.generator_powers_w[1]), C64::from_re(-1.0));
    close(C64::from_re(result.source_dissipation_w), C64::ONE);
    close(C64::from_re(result.inlet_power_w), C64::ZERO);
    balanced(&result);
}

#[test]
fn coherent_drives_keep_indirect_closed_modes_and_disconnected_free_modes() {
    let waves = response([
        0.1, 0.0, 0.2, 0.0,
        0.0, 0.1, 0.2, 0.0,
        0.2, 0.2, 0.1, 0.0,
        0.0, 0.0, 0.0, 1.0,
    ].map(C64::from_re).to_vec(), vec![C64::ONE; 4]);
    let impedances = [C64::new(0.6, 0.2), C64::new(0.3, -0.1)];
    let drives = [C64::new(0.8, 0.3), C64::new(-0.2, 0.7)];
    let conditions = |p: [C64; 2]| [
        pressure(p[0], impedances[0]), pressure(p[1], impedances[1]),
        PortExcitation::VolumeFlow(C64::ZERO), PortExcitation::VolumeFlow(C64::ZERO),
    ];
    let (actual, native) = waves.source_excitation(&conditions(drives)).unwrap();
    // Independent elimination: U2=0 gives a2=2(a0+a1)/9; the two driven
    // channels then have diagonal reflection 13/90 and off-diagonal 2/45.
    let r = 13.0/90.0;
    let a = impedances.map(|z| C64::from_re(1.0+r)+z.scale(1.0-r));
    let b = impedances.map(|z| (C64::ONE-z).scale(2.0/45.0));
    let det = a[0]*a[1]-b[0]*b[1];
    let a0 = (drives[0]*a[1]-b[0]*drives[1])/det;
    let a1 = (a[0]*drives[1]-drives[0]*b[1])/det;
    for (value, expected) in native.iter().zip([a0, a1, (a0+a1).scale(2.0/9.0), C64::ZERO]) {
        close(*value, expected);
    }
    close(actual.inlet_flow[2], C64::ZERO);
    assert_eq!(actual.inlet_pressure[3], C64::ZERO);
    assert_eq!(actual.inlet_flow[3], C64::ZERO);
    assert!(actual.inlet_power_w > 0.0);
    balanced(&actual);
    let left = waves.drive_sources(&conditions([drives[0], C64::ZERO])).unwrap();
    let right = waves.drive_sources(&conditions([C64::ZERO, drives[1]])).unwrap();
    for i in 0..4 {
        close(actual.inlet_pressure[i], left.inlet_pressure[i]+right.inlet_pressure[i]);
        close(actual.inlet_flow[i], left.inlet_flow[i]+right.inlet_flow[i]);
    }
    // Cross-source terms survive. Adding powers from isolated runs is wrong.
    assert!((actual.inlet_power_w-left.inlet_power_w-right.inlet_power_w).abs() > 1e-4);
}

#[test]
fn source_boundary_nodes_remain_finite_but_driven_resonances_refuse() {
    let mut waves = response(vec![C64::ONE, C64::ZERO, C64::ZERO, C64::ONE], vec![C64::ONE; 2]);
    let drive = C64::new(0.4, -0.7);
    let result = waves.drive_sources(&[
        pressure(drive, C64::new(2.0, -0.5)), PortExcitation::VolumeFlow(C64::ZERO),
    ]).unwrap();
    close(result.inlet_pressure[0], drive);
    assert_eq!(result.inlet_flow, vec![C64::ZERO; 2]);
    assert_eq!(result.inlet_pressure[1], C64::ZERO);
    let quiet = waves.drive_sources(&[PortExcitation::VolumeFlow(C64::ZERO); 2]).unwrap();
    assert_eq!(quiet.inlet_pressure, vec![C64::ZERO; 2]);
    assert!(matches!(waves.drive_sources(&[
        PortExcitation::VolumeFlow(C64::ONE), PortExcitation::VolumeFlow(C64::ZERO),
    ]), Err(DuctError::Singular)));
    for coupling in [1.0, 1e-100, f64::MIN_POSITIVE] {
        waves.reflection_matrix[2] = C64::from_re(coupling);
        assert!(matches!(waves.drive_sources(&[
            pressure(drive, C64::ONE), PortExcitation::VolumeFlow(C64::ZERO),
        ]), Err(DuctError::Singular)));
    }
    // A reactive source exactly cancels the reactive duct load: finite
    // passive resistance was not fabricated to damp the true resonance.
    let resonant = response(vec![C64::new(0.0, 1.0)], vec![C64::ONE]);
    assert!(matches!(resonant.drive_sources(&[pressure(C64::ONE, C64::new(0.0, -1.0))]),
        Err(DuctError::Singular)));
    let pressure_node = response(vec![-C64::ONE], vec![C64::ONE]);
    let result = pressure_node.drive_sources(&[pressure(drive, C64::from_re(2.0))]).unwrap();
    assert_eq!(result.inlet_pressure[0], C64::ZERO);
    close(result.inlet_flow[0], drive.scale(0.5));
}

#[test]
fn source_laws_refuse_invalid_inputs_and_keep_resistance_work_at_large_scale_ratios() {
    let waves = transformer();
    assert!(waves.drive_sources(&[]).is_err());
    for value in [C64::new(f64::NAN, 0.0), C64::new(0.0, f64::INFINITY)] {
        for source in [PortExcitation::VolumeFlow(value), pressure(value, C64::ONE), pressure(C64::ZERO, value)] {
            assert!(waves.drive_sources(&[source, PortExcitation::VolumeFlow(C64::ZERO)]).is_err());
        }
    }
    assert!(waves.drive_sources(&[
        pressure(C64::ZERO, C64::new(-1e-100, 1e100)), PortExcitation::VolumeFlow(C64::ZERO),
    ]).is_err());
    for (zc, zs) in [(1e-120, 1e120), (1e120, 1e-120)] {
        let waves = response(vec![C64::ZERO], vec![C64::from_re(zc)]);
        let result = waves.drive_sources(&[pressure(C64::ONE, C64::from_re(zs))]).unwrap();
        let flow = 1.0/(zc+zs);
        assert!((result.inlet_flow[0].re/flow-1.0).abs() < 5e-15);
        assert!((result.inlet_pressure[0].re/(zc*flow)-1.0).abs() < 5e-15);
        // This ordering remains representable even when U² underflows.
        let dissipation = (zs*flow)*flow*0.5;
        if dissipation != 0.0 {
            assert!((result.source_dissipation_w/dissipation-1.0).abs() < 5e-15);
        }
        balanced(&result);
    }
}
