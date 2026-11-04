use super::*;
use fs_material::gas::GasSpec;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn close(a: C64, b: C64, scale: f64) {
    assert!(
        (a - b).abs() <= 3.0e-10 * scale.max(b.abs()),
        "{a:?} != {b:?}"
    );
}

fn pressure(p: C64, z: C64) -> PortExcitation {
    PortExcitation::Pressure {
        open_circuit_pressure_pa: p,
        impedance_pa_s_m3: z,
    }
}

fn phase(k: C64, x: f64) -> C64 {
    let (s, c) = det::sin_cos(k.re * x);
    C64::new(c, s).scale(det::exp(-k.im * x))
}

fn power(p: &[C64], u: &[C64]) -> f64 {
    coherent_sum(p.iter().zip(u).map(|(&p, &u)| (p * u.conj()).scale(0.5))).re
}

#[test]
fn mixed_sources_reconstruct_exact_propagating_and_evanescent_cylinder_fields() {
    // G1: independent outgoing modal exponentials, including complex Zc. The
    // production constitutive law supplies k/Zc; no propagation map is reused.
    let gas = gas();
    let radius = 0.03;
    let profile = RadiusProfile::exponential(radius, 0.0, 0.035).unwrap();
    let omega = gas.sound_speed * 200.0;
    let roots = j1_roots_unchecked(2);
    let gammas = [0.0, roots[0], roots[1]];
    for loss in [LossModel::Lossless, LossModel::WideTube, LossModel::Bessel] {
        let prepared = PreparedModalProfile::new(
            &profile,
            &gas,
            omega,
            loss,
            MouthLoad::Matched,
            ModeCount(3),
            RiccatiSteps(48),
            &mut || false,
        )
        .unwrap();
        assert_eq!(prepared.gas_state(), &gas);
        assert_eq!(prepared.omega_rad_s(), omega);
        assert_eq!(prepared.steps_used(), 48);
        let zref = prepared.reference_impedance_pa_s_m3();
        let sources = [
            pressure(C64::new(2.0, -0.5), C64::new(0.7, -0.2).scale(zref)),
            PortExcitation::VolumeFlow(C64::new(1.0e-6, -2.0e-6)),
            pressure(C64::new(-0.4, 0.9), C64::new(0.3, 0.8).scale(zref)),
        ];
        let result = prepared.drive_sources(&sources, &mut || false).unwrap();
        let wave = segment_wave(&gas, radius, omega, loss).unwrap();
        for (mode, gamma) in gammas.into_iter().enumerate() {
            let kt = C64::from_re(gamma / radius);
            let k = crate::select_forward_decaying(
                ((wave.wavenumber - kt) * (wave.wavenumber + kt)).sqrt(),
            );
            let z = wave.characteristic_impedance * wave.wavenumber / k;
            let u0 = match sources[mode] {
                PortExcitation::VolumeFlow(u) => u,
                PortExcitation::Pressure {
                    open_circuit_pressure_pa: p,
                    impedance_pa_s_m3: zs,
                } => p / (z + zs),
            };
            for (station, &x) in result.field.x.iter().enumerate() {
                let u = u0 * phase(k, x);
                close(result.field.flow[station][mode], u, u0.abs());
                close(result.field.pressure[station][mode], z * u, (z * u0).abs());
            }
        }
        assert!(
            (result.inlet.generator_power_w
                - result.inlet.source_dissipation_w
                - result.inlet.inlet_power_w)
                .abs()
                < 2.0e-14
        );
    }
}

#[test]
fn physical_source_and_reflection_survive_an_exact_matched_transverse_cutoff() {
    // G1: k_1=0 has finite constant pressure and zero flow. Z is undefined,
    // while a finite-pressure source and its field remain perfectly regular.
    let gas = gas();
    let radius = j1_roots_unchecked(1)[0] / 128.0;
    let profile = RadiusProfile::exponential(radius, 0.0, 0.04).unwrap();
    let omega = gas.sound_speed * 128.0;
    assert_eq!(omega / gas.sound_speed, j1_roots_unchecked(1)[0] / radius);
    assert!(
        mm_profile_impedance(
            &profile,
            &gas,
            omega,
            LossModel::Lossless,
            MouthLoad::Matched,
            ModeCount(2),
            RiccatiSteps(32)
        )
        .is_err()
    );
    let prepared = PreparedModalProfile::new(
        &profile,
        &gas,
        omega,
        LossModel::Lossless,
        MouthLoad::Matched,
        ModeCount(2),
        RiccatiSteps(32),
        &mut || false,
    )
    .unwrap();
    let sources = [
        PortExcitation::VolumeFlow(C64::from_re(1.0e-6)),
        pressure(C64::new(2.0, -3.0), C64::from_re(1.0e6)),
    ];
    let result = prepared.drive_sources(&sources, &mut || false).unwrap();
    for (p, u) in result.field.pressure.iter().zip(&result.field.flow) {
        close(p[1], C64::new(2.0, -3.0), 1.0);
        assert_eq!(u[1], C64::ZERO);
    }
    let z = prepared.reference_impedance_pa_s_m3();
    let reflection = prepared.at_references(&[z, 2.0 * z]).unwrap();
    close(reflection[0], C64::ZERO, 1.0);
    assert_eq!(reflection[1], C64::ZERO);
    assert_eq!(reflection[2], C64::ZERO);
    close(reflection[3], C64::ONE, 1.0);
    // Unforced disconnected cutoff modes select zero forced amplitude. A
    // nonzero prescribed flow into the same zero-admittance channel refuses.
    let closed = [sources[0], PortExcitation::VolumeFlow(C64::ZERO)];
    let result = prepared.drive_sources(&closed, &mut || false).unwrap();
    assert!(result.field.pressure.iter().all(|p| p[1] == C64::ZERO));
    assert!(
        prepared
            .drive_sources(
                &[sources[0], PortExcitation::VolumeFlow(C64::ONE)],
                &mut || false
            )
            .is_err()
    );
}

#[test]
fn pressure_driven_half_wave_keeps_internal_motion_at_the_inlet_flow_node() {
    // G1: closed mouth, kL=pi: p(x)=p_open cos(kx),
    // U(x)=i p_open sin(kx)/Zc. Both endpoints have zero flow.
    let gas = gas();
    let length = 0.15;
    let profile = RadiusProfile::exponential(0.02, 0.0, length).unwrap();
    let k = core::f64::consts::PI / length;
    let prepared = PreparedModalProfile::new(
        &profile,
        &gas,
        k * gas.sound_speed,
        LossModel::Lossless,
        MouthLoad::Analytic(Termination::Closed),
        ModeCount(1),
        RiccatiSteps(64),
        &mut || false,
    )
    .unwrap();
    let z = prepared.reference_impedance_pa_s_m3();
    let p = C64::new(3.0, -0.7);
    let result = prepared
        .drive_sources(&[pressure(p, C64::from_re(z))], &mut || false)
        .unwrap();
    assert!(result.field.flow[0][0].abs() < 1.0e-12 * p.abs() / z);
    assert!(result.field.flow[32][0].abs() > 0.99 * p.abs() / z);
    for (i, &x) in result.field.x.iter().enumerate() {
        let (s, c) = det::sin_cos(k * x);
        close(result.field.pressure[i][0], p.scale(c), p.abs());
        close(
            result.field.flow[i][0],
            C64::new(0.0, s / z) * p,
            p.abs() / z,
        );
    }
}

#[test]
fn tapered_profile_sources_reuse_original_fields_and_full_reflection() {
    // G1/G3: nonzero taper couples the original modes. Compare the existing
    // prescribed-flow path and independent Z-to-reference port conversion.
    let gas = gas();
    let profile = RadiusProfile::exponential(0.025, 8.0, 0.06).unwrap();
    let omega = gas.sound_speed * 300.0;
    let sources = [
        C64::new(1.0e-6, 2.0e-6),
        C64::new(-0.8e-6, 0.4e-6),
        C64::ZERO,
    ];
    let prepared = PreparedModalProfile::new(
        &profile,
        &gas,
        omega,
        LossModel::Lossless,
        MouthLoad::Matched,
        ModeCount(3),
        RiccatiSteps(96),
        &mut || false,
    )
    .unwrap();
    let result = prepared
        .drive_sources(&sources.map(PortExcitation::VolumeFlow), &mut || false)
        .unwrap();
    let legacy = mm_profile_field(
        &profile,
        &gas,
        omega,
        LossModel::Lossless,
        MouthLoad::Matched,
        ModeCount(3),
        RiccatiSteps(96),
        &sources,
    )
    .unwrap();
    for (a, b) in result
        .field
        .pressure
        .iter()
        .flatten()
        .zip(legacy.pressure.iter().flatten())
    {
        close(*a, *b, 1.0);
    }
    for (a, b) in result
        .field
        .flow
        .iter()
        .flatten()
        .zip(legacy.flow.iter().flatten())
    {
        close(*a, *b, 1.0e-6);
    }
    let zref = prepared.reference_impedance_pa_s_m3();
    let references = [0.7 * zref, 1.3 * zref, 2.0 * zref];
    let reflection = prepared.at_references(&references).unwrap();
    let z = mm_profile_impedance(
        &profile,
        &gas,
        omega,
        LossModel::Lossless,
        MouthLoad::Matched,
        ModeCount(3),
        RiccatiSteps(96),
    )
    .unwrap()
    .impedance_matrix;
    // Independent boundary system (Z+D) U=2 a. Never use the new chart.
    let mut boundary = z.clone();
    for i in 0..3 {
        boundary[i * 3 + i] = boundary[i * 3 + i] + C64::from_re(references[i]);
    }
    let lu = fs_la::eigen_complex::lu_complex(&boundary, 3).unwrap();
    for j in 0..3 {
        let mut u = [C64::ZERO; 3];
        u[j] = C64::from_re(2.0);
        lu.solve(&mut u);
        for i in 0..3 {
            let b = (if i == j { C64::ONE } else { C64::ZERO }) - u[i].scale(references[i]);
            close(reflection[i * 3 + j], b, 1.0);
        }
    }
    let incoming = result.inlet.inlet_power_w;
    let outgoing = power(
        result.field.pressure.last().unwrap(),
        result.field.flow.last().unwrap(),
    );
    assert!(incoming > 0.0);
    assert!((incoming - outgoing).abs() < 3.0e-10 * incoming);
}

#[test]
fn mutual_mouth_load_and_coherent_sources_keep_actual_source_work() {
    // G1: the finite matrix is the actual terminal boundary, not diagonalized.
    let gas = gas();
    let profile = RadiusProfile::exponential(0.03, 4.0, 0.04).unwrap();
    let omega = gas.sound_speed * 220.0;
    let z =
        gas.density * gas.sound_speed / (core::f64::consts::PI * profile.outlet_radius().powi(2));
    let load = ModalMouthImpedance::new(
        omega,
        profile.outlet_radius(),
        gas,
        ModeCount(2),
        vec![
            C64::new(0.8, -0.2).scale(z),
            C64::from_re(0.3 * z),
            C64::from_re(0.3 * z),
            C64::new(1.1, 0.4).scale(z),
        ],
    )
    .unwrap();
    let prepared = PreparedModalProfile::with_mouth_matrix(
        &profile,
        &gas,
        omega,
        LossModel::Lossless,
        &load,
        ModeCount(2),
        RiccatiSteps(64),
        &mut || false,
    )
    .unwrap();
    let drive = [
        pressure(C64::new(2.0, 0.7), C64::new(0.4, -0.3).scale(z)),
        pressure(C64::new(-0.3, 1.0), C64::new(0.2, 0.1).scale(z)),
    ];
    let result = prepared.drive_sources(&drive, &mut || false).unwrap();
    let p = result.field.pressure.last().unwrap();
    let u = result.field.flow.last().unwrap();
    for i in 0..2 {
        close(
            p[i],
            load.impedance()[i * 2] * u[0] + load.impedance()[i * 2 + 1] * u[1],
            1.0,
        );
    }
    assert!((result.inlet.inlet_power_w - power(p, u)).abs() < 1.0e-14);
    assert!(
        (result.inlet.generator_power_w
            - result.inlet.source_dissipation_w
            - result.inlet.inlet_power_w)
            .abs()
            < 1.0e-14
    );
    let left = prepared
        .drive_sources(
            &[drive[0], pressure(C64::ZERO, C64::new(0.2, 0.1).scale(z))],
            &mut || false,
        )
        .unwrap();
    let right = prepared
        .drive_sources(
            &[pressure(C64::ZERO, C64::new(0.4, -0.3).scale(z)), drive[1]],
            &mut || false,
        )
        .unwrap();
    for i in 0..result.field.pressure.len() {
        for m in 0..2 {
            close(
                result.field.pressure[i][m],
                left.field.pressure[i][m] + right.field.pressure[i][m],
                1.0,
            );
        }
    }
    assert!(
        (result.inlet.inlet_power_w - left.inlet.inlet_power_w - right.inlet.inlet_power_w).abs()
            > 1.0e-9,
        "coherent sources must retain cross terms in acoustic power"
    );
}

#[test]
fn preparation_and_source_queries_refuse_bad_inputs_and_cancel_without_state_change() {
    let gas = gas();
    let profile = RadiusProfile::exponential(0.02, 5.0, 0.03).unwrap();
    let make = |cancelled: &mut dyn FnMut() -> bool| {
        PreparedModalProfile::new(
            &profile,
            &gas,
            1200.0,
            LossModel::Lossless,
            MouthLoad::Matched,
            ModeCount(2),
            RiccatiSteps(32),
            &mut || cancelled(),
        )
    };
    assert!(matches!(
        make(&mut || true),
        Err(DuctError::Cancelled { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        make(&mut || {
            polls += 1;
            polls > 8
        }),
        Err(DuctError::Cancelled { .. })
    ));
    let prepared = make(&mut || false).unwrap();
    let sources = [
        pressure(C64::ONE, C64::from_re(2.0e5)),
        PortExcitation::VolumeFlow(C64::ZERO),
    ];
    let before = prepared.drive_sources(&sources, &mut || false).unwrap();
    polls = 0;
    assert!(matches!(
        prepared.drive_sources(&sources, &mut || {
            polls += 1;
            polls > 8
        }),
        Err(DuctError::Cancelled { .. })
    ));
    let after = prepared.drive_sources(&sources, &mut || false).unwrap();
    assert_eq!(before.field.pressure, after.field.pressure);
    assert_eq!(before.field.flow, after.field.flow);
    assert!(
        prepared
            .drive_sources(&sources[..1], &mut || false)
            .is_err()
    );
    for bad_source in [
        pressure(C64::ONE, C64::from_re(-1.0)),
        pressure(C64::new(f64::NAN, 0.0), C64::ONE),
        PortExcitation::VolumeFlow(C64::new(f64::INFINITY, 0.0)),
    ] {
        assert!(
            prepared
                .drive_sources(&[bad_source, sources[1]], &mut || false)
                .is_err()
        );
    }
    assert!(prepared.at_references(&[1.0, 0.0]).is_err());
    assert!(prepared.at_references(&[1.0]).is_err());
}
