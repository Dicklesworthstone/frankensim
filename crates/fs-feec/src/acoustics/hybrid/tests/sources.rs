use super::*;

// Solve in physical pressure/volume-flow coordinates. No scattering,
// transformer waves or hybrid condensation is used by this oracle.
fn physical_source_closure(
    z: &[C64],
    y: &[C64],
    map: &[f64],
    source_flow: &[C64],
    h_ref: f64,
    l_ref: f64,
    incident: [C64; 2],
) -> Vec<C64> {
    let m = (map.len() as f64).sqrt() as usize;
    let n = m + 1;
    let hports = (1..n).collect::<Vec<_>>();
    let lports = (0..n).filter(|&p| p != 1).collect::<Vec<_>>();
    let mut a = vec![C64::ZERO; 4 * n * n];
    let mut rhs = vec![C64::ZERO; 2 * n];
    a[..n].copy_from_slice(&z[..n]);
    a[0] = a[0] + C64::from_re(h_ref);
    rhs[0] = incident[0].scale(2.0 * h_ref.sqrt());
    for i in 0..m {
        let row = 1 + i;
        a[row * 2 * n..row * 2 * n + n].copy_from_slice(&z[hports[i] * n..(hports[i] + 1) * n]);
        for j in 0..m {
            a[row * 2 * n + n + lports[j]] = C64::from_re(-map[j * m + i]);
        }
        let row = n + i;
        a[row * 2 * n + n..(row + 1) * 2 * n]
            .copy_from_slice(&y[lports[i] * n..(lports[i] + 1) * n]);
        rhs[row] = -source_flow[lports[i]];
        for j in 0..m {
            a[row * 2 * n + hports[j]] = C64::from_re(map[i * m + j]);
        }
    }
    for j in 0..n {
        a[(2 * n - 1) * 2 * n + n + j] = y[n + j].scale(l_ref);
    }
    a[(2 * n - 1) * 2 * n + n + 1] = a[(2 * n - 1) * 2 * n + n + 1] + C64::ONE;
    rhs[2 * n - 1] = incident[1].scale(2.0 * l_ref.sqrt()) - source_flow[1].scale(l_ref);
    lu_complex_with_cancel(&a, 2 * n, || false)
        .unwrap()
        .solve_with_cancel(&mut rhs, || false)
        .unwrap();
    rhs
}

fn source_work_closes(response: &HybridResponse) {
    let p = response.power;
    let scale =
        (p.external_input_w.abs() + p.source_input_w.abs() + p.dissipated_w().abs()).max(1e-30);
    assert!(p.balance_residual_w().abs() < 5e-8 * scale, "{p:?}");
    assert!(p.interface_residual_w().abs() < 2e-8 * scale, "{p:?}");
    let ps = response
        .helmholtz
        .pressure_pa
        .iter()
        .chain(&response.thermoviscous.pressure_drive_pa)
        .map(|z| z.abs())
        .fold(1e-30, f64::max);
    let qs = response
        .helmholtz
        .inward_flow_m3_s
        .iter()
        .chain(&response.thermoviscous.inward_flow_m3_s)
        .map(|z| z.abs())
        .fold(1e-30, f64::max);
    assert!(
        response
            .interface_pressure_residual_pa
            .iter()
            .all(|z| z.abs() < 2e-8 * ps)
    );
    assert!(
        response
            .interface_flow_residual_m3_s
            .iter()
            .all(|z| z.abs() < 2e-8 * qs)
    );
}

#[test]
fn g2_spatial_heat_drives_both_regions_with_independent_nonorthogonal_port_closure() {
    let size = [0.04, 0.02, 0.02];
    let lsize = [0.01, size[1], size[2]];
    let volume = lsize.iter().product::<f64>();
    let (hm, hp) = box_mesh([1, 2, 2], size, 0.0);
    let (lm, lp) = box_mesh([1, 2, 2], lsize, size[0]);
    let hf = HelmholtzFe::new(&hm, &hp, 2, medium()).unwrap();
    let lf = HelmholtzFe::new(&lm, &lp, 2, medium()).unwrap();
    let hs = shapes(&hf, 0.0, size[0], size[1], true, false);
    let ls = shapes(&lf, size[0], lsize[0], size[1], true, true);
    let h = InteriorMultiport::new(&hf, &complete(&hf), &hs, false, &|| false).unwrap();
    let l = ThermoviscousFe::new(
        &lm,
        &lp,
        gas(),
        &boundary(&lf, &ls, false),
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap();
    let source = l
        .volume_heat_source(
            4,
            &mut |_, p| Ok((1.0 + p[1] / size[1]) / volume),
            &mut || false,
        )
        .unwrap();
    let differential = l
        .volume_heat_source(
            4,
            &mut |_, p| Ok((2.0 * p[1] / size[1] - 1.0) / volume),
            &mut || false,
        )
        .unwrap();
    assert!(differential.integrated_shape().abs() < 1e-13);
    let sources = [
        (&source, C64::new(0.0012, -0.0007)),
        (&differential, C64::new(0.0009, 0.0003)),
    ];
    let hybrid = HybridInterior::new(&h, &l, &[1, 2], &[0, 2], &|| false).unwrap();
    let frequency = 1700.0;
    let cp = C64::new(970.0, -25.0);
    let hr = [1.1e6, 0.8e6, 2.5e6];
    let lr = [1.7e6, 0.9e6, 3.1e6];
    let prepared = hybrid.prepare(frequency, cp, &hr, &lr, &|| false).unwrap();
    let incident = [C64::new(0.0004, 0.0001), C64::new(-0.0003, 0.0007)];
    let response = prepared
        .response_with_sources(&incident, &[], &sources, &|| false)
        .unwrap();
    source_work_closes(&response);
    assert_eq!(
        response.power.source_input_w,
        response.thermoviscous.power.heat_source_work_w
    );
    let zero_traction = l
        .prepare_with_heat_capacity(frequency, cp, &mut || false)
        .unwrap()
        .response_with_heat_sources(&[C64::ZERO; 3], &[], &sources, &mut || false)
        .unwrap();
    let z = h.impedance(frequency, &|| false).unwrap();
    let y = l
        .admittance_with_heat_capacity(frequency, cp, &mut || false)
        .unwrap();
    let oracle = physical_source_closure(
        &z,
        &y,
        hybrid.interface().flow_map(),
        &zero_traction.inward_flow_m3_s,
        hr[0],
        lr[1],
        incident,
    );
    let qscale = oracle[..3].iter().map(|z| z.abs()).fold(1e-30, f64::max);
    let pscale = oracle[3..].iter().map(|z| z.abs()).fold(1e-30, f64::max);
    for i in 0..3 {
        assert!((response.helmholtz.inward_flow_m3_s[i] - oracle[i]).abs() < 2e-8 * qscale);
        assert!(
            (response.thermoviscous.pressure_drive_pa[i] - oracle[3 + i]).abs() < 2e-8 * pscale
        );
    }
    let emitted = prepared
        .response_with_sources(&[C64::ZERO; 2], &[], &sources, &|| false)
        .unwrap();
    source_work_closes(&emitted);
    assert!(emitted.power.external_input_w < 0.0);
    assert!(emitted.power.source_input_w > emitted.power.dissipated_w());
    assert!(emitted.helmholtz.amplitudes.iter().any(|z| z.abs() > 0.0));
    let unforced = prepared.response(&incident, &|| false).unwrap();
    for i in 0..2 {
        let expected = unforced.outgoing[i] + emitted.outgoing[i];
        assert!(
            (response.outgoing[i] - expected).abs()
                < 2e-8 * (unforced.outgoing[i].abs() + emitted.outgoing[i].abs())
        );
    }
    let wave_work = 0.5
        * (incident.iter().map(|z| z.norm_sq()).sum::<f64>()
            - response.outgoing.iter().map(|z| z.norm_sq()).sum::<f64>());
    assert!(
        (wave_work - response.power.external_input_w).abs()
            < 2e-9 * (wave_work.abs() + response.power.external_input_w.abs())
    );
    // An internal power-wave reference changes coordinates, never fields.
    let other = hybrid
        .prepare(
            frequency,
            cp,
            &[hr[0], 3.3e6, 0.4e6],
            &[0.6e6, lr[1], 1.9e6],
            &|| false,
        )
        .unwrap()
        .response_with_sources(&incident, &[], &sources, &|| false)
        .unwrap();
    for (&a, &b) in response.outgoing.iter().zip(&other.outgoing) {
        assert!((a - b).abs() < 2e-8 * (a.abs() + b.abs()));
    }
    source_work_closes(&other);
    let zero = prepared
        .response_with_sources(&incident, &[], &[], &|| false)
        .unwrap();
    assert_eq!(zero.helmholtz.amplitudes, unforced.helmholtz.amplitudes);
    assert_eq!(
        zero.thermoviscous.temperature_k,
        unforced.thermoviscous.temperature_k
    );
    assert!(
        prepared
            .response_with_sources(&incident[..1], &[], &sources, &|| false)
            .is_err()
    );
    assert!(matches!(
        prepared.response_with_sources(&incident, &[], &sources, &|| true),
        Err(HybridError::Cancelled)
    ));
}

#[test]
fn g1_heat_excites_a_sealed_hybrid_with_no_external_pressure_ports() {
    let size = [0.02, 0.02, 0.02];
    let (hm, hp) = box_mesh([1, 1, 1], size, 0.0);
    let (lm, lp) = box_mesh([1, 1, 1], size, size[0]);
    let hf = HelmholtzFe::new(&hm, &hp, 2, medium()).unwrap();
    let lf = HelmholtzFe::new(&lm, &lp, 2, medium()).unwrap();
    let hs = [BoundaryMode::uniform(&ends(&hf, size[0]))];
    let ls = [BoundaryMode::uniform(&ends(&lf, size[0]))];
    let h = InteriorMultiport::new(&hf, &complete(&hf), &hs, false, &|| false).unwrap();
    let l = ThermoviscousFe::new(
        &lm,
        &lp,
        gas(),
        &boundary(&lf, &ls, false),
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap();
    let source = l
        .volume_heat_source(
            4,
            &mut |_, _| Ok(1.0 / size.iter().product::<f64>()),
            &mut || false,
        )
        .unwrap();
    let sources = [(&source, C64::new(0.001, 0.0003))];
    let hybrid = HybridInterior::new(&h, &l, &[0], &[0], &|| false).unwrap();
    assert_eq!(hybrid.external_port_count(), 0);
    let frequency = 1300.0;
    let cp = C64::new(970.0, -45.0);
    let prepared = hybrid
        .prepare(frequency, cp, &[1.2e6], &[2.3e6], &|| false)
        .unwrap();
    assert!(prepared.scattering(&|| false).unwrap().is_empty());
    let response = prepared
        .response_with_sources(&[], &[], &sources, &|| false)
        .unwrap();
    assert!(response.outgoing.is_empty() && response.pressure_pa.is_empty());
    assert_eq!(response.power.external_input_w, 0.0);
    assert!(response.power.source_input_w > 0.0);
    assert!(response.power.thermoviscous.caloric_w > 0.0);
    source_work_closes(&response);
    // Scalar physical closure for the sealed interface: Q_H+Y P_L+Q_s=0,
    // P_L=Z Q_H. This independently fixes the source sign with no exterior.
    let z = h.impedance(frequency, &|| false).unwrap()[0];
    let y = l
        .admittance_with_heat_capacity(frequency, cp, &mut || false)
        .unwrap()[0];
    let q_source = l
        .prepare_with_heat_capacity(frequency, cp, &mut || false)
        .unwrap()
        .response_with_heat_sources(&[C64::ZERO], &[], &sources, &mut || false)
        .unwrap()
        .inward_flow_m3_s[0];
    let expected_flow = -q_source / (C64::ONE + y * z);
    let expected_pressure = z * expected_flow;
    assert!(relative(response.helmholtz.inward_flow_m3_s[0], expected_flow) < 2e-8);
    assert!(
        relative(
            response.thermoviscous.pressure_drive_pa[0],
            expected_pressure
        ) < 2e-8
    );
}
