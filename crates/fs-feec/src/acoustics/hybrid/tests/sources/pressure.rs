//! Physical pressure-region injection coupled to the original full LNS field.

use super::*;
use std::cell::Cell;

#[test]
fn g2_sources_in_both_regions_close_original_physical_equations_and_work() {
    let size = [0.04, 0.02, 0.02];
    let lsize = [0.01, size[1], size[2]];
    let hvolume = size.iter().product::<f64>();
    let lvolume = lsize.iter().product::<f64>();
    let (hm, hp) = box_mesh([1, 2, 2], size, 0.0);
    let (lm, lp) = box_mesh([1, 2, 2], lsize, size[0]);
    let hf = HelmholtzFe::new(&hm, &hp, 2, medium()).unwrap();
    let lf = HelmholtzFe::new(&lm, &lp, 2, medium()).unwrap();
    let hs = shapes(&hf, 0.0, size[0], size[1], true, false);
    let ls = shapes(&lf, size[0], lsize[0], size[1], true, true);
    let h = InteriorMultiport::new_full(&hf, &hs, hf.ndof(), &|| false).unwrap();
    let l = ThermoviscousFe::new(
        &lm, &lp, gas(), &boundary(&lf, &ls, false),
        ThermoviscousLimits::default(), &mut || false,
    ).unwrap();
    let calls = Cell::new(0);
    let injection = h.volume_velocity_source(
        &hf, 4, &mut |_, p| {
            calls.set(calls.get() + 1);
            Ok::<_, HeterogeneousError>((2.0 * p[1] / size[1] - 1.0) / hvolume)
        }, &mut || false,
    ).unwrap();
    assert!(injection.integrated_shape().abs() < 1e-13);
    let face = hf.boundary_faces().into_iter().find(|face| {
        face.iter().all(|&v| hp[v as usize][1] == 0.0)
    }).unwrap();
    let motion = h.quadratic_velocity_source(
        &hf, &[(face, [0.0, 0.0, 0.0, 1.0, 0.0, 0.0])], &mut || false,
    ).unwrap();
    let heat = l.volume_heat_source(
        4, &mut |_, p| Ok((1.0 + p[1] / size[1]) / lvolume), &mut || false,
    ).unwrap();
    let force = l.volume_force_source(
        4, &mut |_, _| Ok::<_, ThermoviscousError>([1.0 / lvolume, 0.0, 0.0]), &mut || false,
    ).unwrap();
    let hvolumes = [(&injection, C64::new(2e-7, -7e-8))];
    let hboundaries = [(&motion, C64::new(3e-4, 1e-4))];
    let heats = [(&heat, C64::new(0.0012, -0.0007))];
    let forces = [(&force, C64::new(2e-5, -7e-6))];
    let sampled = calls.get();
    let hybrid = HybridInterior::new(&h, &l, &[1, 2], &[0, 2], &|| false).unwrap();
    let frequency = 1700.0;
    let cp = C64::new(970.0, -25.0);
    let hr = [1.1e6, 0.8e6, 2.5e6];
    let lr = [1.7e6, 0.9e6, 3.1e6];
    let incident = [C64::new(0.0004, 0.0001), C64::new(-0.0003, 0.0007)];
    let prepared = hybrid.prepare(frequency, cp, &hr, &lr, &|| false).unwrap();
    let field = prepared.response_with_all_sources(
        &incident, &hboundaries, &hvolumes, &[], &forces, &heats, &|| false,
    ).unwrap();
    source_work_closes(&field);

    // P_H=Z Q_H+P_source and Q_L=Y P_L+Q_source, with physical external
    // Robin loads and reciprocal T. No loaded scattering or hybrid factor
    // participates in this independently assembled pressure/flow closure.
    let unloaded_h = h.prepare(frequency, &|| false).unwrap();
    let unloaded_l = l.prepare_with_heat_capacity(frequency, cp, &mut || false).unwrap();
    let hsource = unloaded_h.response_with_sources(
        &[C64::ZERO; 3], &hboundaries, &hvolumes, &|| false,
    ).unwrap();
    let lsource = unloaded_l.response_with_volume_sources(
        &[C64::ZERO; 3], &[], &forces, &heats, &mut || false,
    ).unwrap();
    let z = h.impedance(frequency, &|| false).unwrap();
    let y = l.admittance_with_heat_capacity(frequency, cp, &mut || false).unwrap();
    let oracle = physical_source_closure_both(
        &z, &y, hybrid.interface().flow_map(), &hsource.pressure_pa,
        &lsource.inward_flow_m3_s, hr[0], lr[1], incident,
    );
    let expected_h = unloaded_h.response_with_sources(
        &oracle[..3], &hboundaries, &hvolumes, &|| false,
    ).unwrap();
    let expected_l = unloaded_l.response_with_volume_sources(
        &oracle[3..], &[], &forces, &heats, &mut || false,
    ).unwrap();
    let hscale = expected_h.amplitudes.iter().map(|p| C64::new(p[0], p[1]).abs())
        .fold(1e-30_f64, f64::max);
    assert!(field.helmholtz.amplitudes.iter().zip(&expected_h.amplitudes)
        .all(|(p, q)| C64::new(p[0] - q[0], p[1] - q[1]).abs() < 2e-8 * hscale));
    for (actual, expected) in [
        (&field.thermoviscous.pressure_pa, &expected_l.pressure_pa),
        (&field.thermoviscous.temperature_k, &expected_l.temperature_k),
        (&field.thermoviscous.inward_flow_m3_s, &expected_l.inward_flow_m3_s),
    ] {
        let scale = expected.iter().map(|z| z.abs()).fold(1e-30_f64, f64::max);
        assert!(actual.iter().zip(expected).all(|(&a, &b)| (a - b).abs() < 2e-8 * scale));
    }
    let vscale = expected_l.velocity_m_s.iter().flatten().map(|v| v.abs())
        .fold(1e-30_f64, f64::max);
    assert!(field.thermoviscous.velocity_m_s.iter().flatten()
        .zip(expected_l.velocity_m_s.iter().flatten())
        .all(|(&a, &b)| (a - b).abs() < 2e-8 * vscale));
    assert_eq!(field.power.source_input_w,
        field.thermoviscous.moving_wall.external_input_w
            + field.thermoviscous.power.heat_source_work_w
            + field.thermoviscous.power.body_force_work_w
            + field.helmholtz.boundary_source_power_w
            + field.helmholtz.volume_source_power_w);

    let h_only = prepared.response_with_all_sources(
        &[C64::ZERO; 2], &hboundaries, &hvolumes, &[], &[], &[], &|| false,
    ).unwrap();
    source_work_closes(&h_only);
    assert!(h_only.power.external_input_w < 0.0);
    assert!(h_only.power.source_input_w > h_only.power.dissipated_w());
    assert!(h_only.thermoviscous.pressure_drive_pa.iter().any(|p| p.abs() > 1e-12));
    let l_only = prepared.response_with_volume_sources(
        &incident, &[], &forces, &heats, &|| false,
    ).unwrap();
    for i in 0..2 {
        let expected = h_only.outgoing[i] + l_only.outgoing[i];
        assert!((field.outgoing[i] - expected).abs()
            < 2e-8 * (h_only.outgoing[i].abs() + l_only.outgoing[i].abs()));
    }
    let other = hybrid.prepare(
        frequency, cp, &[hr[0], 3.3e6, 0.4e6], &[0.6e6, lr[1], 1.9e6], &|| false,
    ).unwrap().response_with_all_sources(
        &incident, &hboundaries, &hvolumes, &[], &forces, &heats, &|| false,
    ).unwrap();
    source_work_closes(&other);
    assert!(field.helmholtz.amplitudes.iter().zip(&other.helmholtz.amplitudes)
        .all(|(p, q)| C64::new(p[0] - q[0], p[1] - q[1]).abs() < 2e-8 * hscale));
    assert_eq!(calls.get(), sampled);
}
