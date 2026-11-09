//! Original distributed forces through the existing reciprocal hybrid closure.
use super::*;
use crate::acoustics::thermoviscous::body_force::MAX_BODY_FORCE_SOURCES;
use crate::highorder::simplex::SimplexSpace;
use std::cell::Cell;

fn same_fields(actual: &ThermoviscousField, expected: &ThermoviscousField) {
    for (a, b) in [
        (&actual.pressure_pa, &expected.pressure_pa),
        (&actual.temperature_k, &expected.temperature_k),
        (&actual.inward_flow_m3_s, &expected.inward_flow_m3_s),
    ] {
        assert_eq!(a.len(), b.len());
        let scale = b.iter().map(|z| z.abs()).fold(1e-30_f64, f64::max);
        assert!(a.iter().zip(b).all(|(a, b)| (*a - *b).abs() < 2e-8 * scale));
    }
    let scale = expected.velocity_m_s.iter().flatten()
        .map(|z| z.abs()).fold(1e-30_f64, f64::max);
    assert!(actual.velocity_m_s.iter().flatten()
        .zip(expected.velocity_m_s.iter().flatten())
        .all(|(a, b)| (*a - *b).abs() < 2e-8 * scale));
}

#[test]
fn g2_original_vector_forces_and_surface_heat_close_nonorthogonal_hybrid_fields() {
    let size = [0.04, 0.02, 0.02];
    let lsize = [0.01, size[1], size[2]];
    let volume = lsize.iter().product::<f64>();
    let (hm, hp) = box_mesh([1, 2, 2], size, 0.0);
    let (lm, lp) = box_mesh([1, 2, 2], lsize, size[0]);
    let hf = HelmholtzFe::new(&hm, &hp, 2, medium()).unwrap();
    let lf = HelmholtzFe::new(&lm, &lp, 2, medium()).unwrap();
    let hs = shapes(&hf, 0.0, size[0], size[1], true, false);
    let ls = shapes(&lf, size[0], lsize[0], size[1], true, true);
    let h = InteriorMultiport::new_full(&hf, &hs, hf.ndof(), &|| false).unwrap();
    let lb = boundary(&lf, &ls, false);
    let l = ThermoviscousFe::new(
        &lm, &lp, gas(), &lb, ThermoviscousLimits::default(), &mut || false,
    ).unwrap();
    let calls = Cell::new(0usize);
    let force = l.volume_force_source(
        4,
        &mut |_, p| {
            calls.set(calls.get() + 1);
            Ok::<_, ThermoviscousError>([
                (1.0 + p[1] / size[1]) / volume,
                0.3 * (2.0 * p[2] / size[2] - 1.0) / volume,
                0.0,
            ])
        },
        &mut || false,
    ).unwrap();
    let differential = l.volume_force_source(
        4,
        &mut |_, p| Ok::<_, ThermoviscousError>([
            (2.0 * p[1] / size[1] - 1.0) / volume, 0.0, 0.0,
        ]),
        &mut || false,
    ).unwrap();
    assert!(differential.integrated_shape().iter().all(|x| x.abs() < 1e-13));
    let heat_face = lb.no_slip_faces[0];
    let (_, _, area) = lf.boundary_geometry(heat_face).unwrap();
    let mut heat_shape = vec![0.0; SimplexSpace::new(&lm, 2).ndof];
    for v in heat_face {
        heat_shape[v as usize] = 1.0 / area;
    }
    let heat = l.surface_heat_source(&[heat_face], &heat_shape, &mut || false).unwrap();
    let forces = [
        (&force, C64::new(2e-5, -7e-6)),
        (&differential, C64::new(9e-6, 3e-6)),
    ];
    let heats = [(&heat, C64::new(0.0012, -0.0007))];
    let sampled = calls.get();
    let hybrid = HybridInterior::new(&h, &l, &[1, 2], &[0, 2], &|| false).unwrap();
    let frequency = 1700.0;
    let cp = C64::new(970.0, -25.0);
    let hr = [1.1e6, 0.8e6, 2.5e6];
    let lr = [1.7e6, 0.9e6, 3.1e6];
    let prepared = hybrid.prepare(frequency, cp, &hr, &lr, &|| false).unwrap();
    let incident = [C64::new(0.0004, 0.0001), C64::new(-0.0003, 0.0007)];
    let response = prepared.response_with_volume_sources(
        &incident, &[], &forces, &heats, &|| false,
    ).unwrap();
    source_work_closes(&response);
    assert!(response.thermoviscous.power.body_force_work_w.abs() > 0.0);
    assert_eq!(
        response.power.source_input_w,
        response.thermoviscous.moving_wall.external_input_w
            + response.thermoviscous.power.heat_source_work_w
            + response.thermoviscous.power.body_force_work_w,
    );

    // Original unloaded source flow enters the independent physical Z/Y
    // closure. Its equations contain neither power waves nor hybrid factors.
    let unloaded = l.prepare_with_heat_capacity(frequency, cp, &mut || false).unwrap();
    let source_only = unloaded.response_with_volume_sources(
        &[C64::ZERO; 3], &[], &forces, &heats, &mut || false,
    ).unwrap();
    let z = h.impedance(frequency, &|| false).unwrap();
    let y = l.admittance_with_heat_capacity(frequency, cp, &mut || false).unwrap();
    let oracle = physical_source_closure(
        &z, &y, hybrid.interface().flow_map(), &source_only.inward_flow_m3_s,
        hr[0], lr[1], incident,
    );
    let expected = unloaded.response_with_volume_sources(
        &oracle[3..], &[], &forces, &heats, &mut || false,
    ).unwrap();
    same_fields(&response.thermoviscous, &expected);
    let expected_h = h.prepare(frequency, &|| false).unwrap()
        .response(&oracle[..3], &|| false).unwrap();
    let scale = expected_h.amplitudes.iter().map(|z| z.abs()).fold(1e-30_f64, f64::max);
    assert!(response.helmholtz.amplitudes.iter().zip(&expected_h.amplitudes)
        .all(|(a, b)| (*a - *b).abs() < 2e-8 * scale));
    for i in 0..3 {
        assert!(relative(response.helmholtz.inward_flow_m3_s[i], oracle[i]) < 2e-8);
        assert!(relative(response.thermoviscous.pressure_drive_pa[i], oracle[3 + i]) < 2e-8);
    }

    let body_only = prepared.response_with_volume_sources(
        &[C64::ZERO; 2], &[], &forces, &[], &|| false,
    ).unwrap();
    source_work_closes(&body_only);
    assert!(body_only.outgoing[0].abs() > 0.0);
    assert!(body_only.power.external_input_w < 0.0);
    let emitted = 0.5 * body_only.outgoing.iter().map(|z| z.norm_sq()).sum::<f64>();
    let scale = body_only.power.source_input_w + emitted + body_only.power.dissipated_w();
    assert!((body_only.power.source_input_w - emitted - body_only.power.dissipated_w()).abs()
        < 5e-8 * scale);
    let heat_and_waves = prepared.response_with_sources(&incident, &[], &heats, &|| false).unwrap();
    for i in 0..2 {
        let expected = heat_and_waves.outgoing[i] + body_only.outgoing[i];
        let scale = heat_and_waves.outgoing[i].abs() + body_only.outgoing[i].abs();
        assert!((response.outgoing[i] - expected).abs() < 2e-8 * scale);
    }
    let changed_reference = hybrid.prepare(
        frequency, cp, &[hr[0], 3.3e6, 0.4e6], &[0.6e6, lr[1], 1.9e6], &|| false,
    ).unwrap().response_with_volume_sources(&incident, &[], &forces, &heats, &|| false).unwrap();
    same_fields(&changed_reference.thermoviscous, &response.thermoviscous);
    source_work_closes(&changed_reference);
    assert_eq!(calls.get(), sampled, "warm hybrid solves must not resample original force shapes");
}

#[test]
fn g1_sealed_body_force_response_preserves_original_work_binding_and_cancellation() {
    let size = [0.02, 0.02, 0.02];
    let volume = size.iter().product::<f64>();
    let (hm, hp) = box_mesh([1, 1, 1], size, 0.0);
    let (lm, lp) = box_mesh([1, 1, 1], size, size[0]);
    let hf = HelmholtzFe::new(&hm, &hp, 2, medium()).unwrap();
    let lf = HelmholtzFe::new(&lm, &lp, 2, medium()).unwrap();
    let hs = [BoundaryMode::uniform(&ends(&hf, size[0]))];
    let ls = [BoundaryMode::uniform(&ends(&lf, size[0]))];
    let h = InteriorMultiport::new_full(&hf, &hs, hf.ndof(), &|| false).unwrap();
    let lb = boundary(&lf, &ls, false);
    let make = || ThermoviscousFe::new(
        &lm, &lp, gas(), &lb, ThermoviscousLimits::default(), &mut || false,
    ).unwrap();
    let l = make();
    let force = l.volume_force_source(
        4, &mut |_, _| Ok::<_, ThermoviscousError>([1.0 / volume, 0.0, 0.0]), &mut || false,
    ).unwrap();
    let hybrid = HybridInterior::new(&h, &l, &[0], &[0], &|| false).unwrap();
    let frequency = 1300.0;
    let cp = C64::new(970.0, -45.0);
    let prepared = hybrid.prepare(frequency, cp, &[1.2e6], &[2.3e6], &|| false).unwrap();
    let amplitude = C64::new(2e-5, 3e-6);
    let forces = [(&force, amplitude)];
    let response = prepared.response_with_volume_sources(&[], &[], &forces, &[], &|| false).unwrap();
    assert!(response.outgoing.is_empty());
    assert_eq!(response.power.external_input_w, 0.0);
    assert!(response.power.source_input_w > 0.0);
    source_work_closes(&response);
    // Q_H+Y*P_L+Q_force=0 and P_L=Z*Q_H independently fix the sign.
    let z = h.impedance(frequency, &|| false).unwrap()[0];
    let y = l.admittance_with_heat_capacity(frequency, cp, &mut || false).unwrap()[0];
    let q_source = l.prepare_with_heat_capacity(frequency, cp, &mut || false).unwrap()
        .response_with_volume_sources(&[C64::ZERO], &[], &forces, &[], &mut || false)
        .unwrap().inward_flow_m3_s[0];
    let expected_flow = -q_source / (C64::ONE + y * z);
    assert!(relative(response.helmholtz.inward_flow_m3_s[0], expected_flow) < 2e-8);
    assert!(relative(response.thermoviscous.pressure_drive_pa[0], z * expected_flow) < 2e-8);
    let cancelled_sources = [(&force, amplitude), (&force, -amplitude)];
    let zero = prepared.response_with_volume_sources(
        &[], &[], &cancelled_sources, &[], &|| false,
    ).unwrap();
    assert!(zero.helmholtz.amplitudes.iter().all(|z| *z == C64::ZERO));
    assert!(zero.thermoviscous.velocity_m_s.iter().flatten().all(|z| *z == C64::ZERO));
    assert_eq!(zero.power.source_input_w, 0.0);

    let legacy = prepared.response_with_sources(&[], &[], &[], &|| false).unwrap();
    let empty = prepared.response_with_volume_sources(&[], &[], &[], &[], &|| false).unwrap();
    assert_eq!(legacy.helmholtz.amplitudes, empty.helmholtz.amplitudes);
    assert_eq!(legacy.thermoviscous.pressure_pa, empty.thermoviscous.pressure_pa);
    let other = make();
    let foreign = other.volume_force_source(
        4, &mut |_, _| Ok::<_, ThermoviscousError>([1.0 / volume, 0.0, 0.0]), &mut || false,
    ).unwrap();
    assert!(prepared.response_with_volume_sources(
        &[], &[], &[(&foreign, C64::ZERO)], &[], &|| false,
    ).is_err());
    assert!(prepared.response_with_volume_sources(
        &[], &[], &[(&force, C64::new(f64::NAN, 0.0))], &[], &|| false,
    ).is_err());
    let at_cap = vec![(&force, C64::ZERO); MAX_BODY_FORCE_SOURCES];
    assert!(prepared.response_with_volume_sources(&[], &[], &at_cap, &[], &|| false).is_ok());
    assert!(matches!(prepared.response_with_volume_sources(
        &[], &[], &vec![(&force, C64::ZERO); MAX_BODY_FORCE_SOURCES + 1], &[], &|| false,
    ), Err(HybridError::Thermoviscous(ThermoviscousError::Budget { .. }))));
    let polls = Cell::new(0usize);
    assert!(matches!(prepared.response_with_volume_sources(
        &[], &[], &forces, &[], &|| {
            polls.set(polls.get() + 1);
            polls.get() == 12
        },
    ), Err(HybridError::Cancelled | HybridError::Thermoviscous(ThermoviscousError::Cancelled))));
    let replay = prepared.response_with_volume_sources(&[], &[], &forces, &[], &|| false).unwrap();
    assert_eq!(replay.helmholtz.amplitudes, response.helmholtz.amplitudes);
    assert_eq!(replay.thermoviscous.pressure_pa, response.thermoviscous.pressure_pa);
    assert_eq!(replay.power.source_input_w, response.power.source_input_w);
}
