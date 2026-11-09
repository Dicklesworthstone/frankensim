//! G1/G3 complete FE/BEM closure with independently integrated source loads.

use super::*;
use fs_feec::SimplexSpace;
use fs_feec::acoustics::helmholtz::heterogeneous::HeterogeneousError;

const VOLUME: f64 = LENGTH * WIDTH * WIDTH;
const MIDSIDE: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];

// Kuhn's six affine tetrahedra have equal volume. The hierarchical P2 basis
// contains lambda_i and lambda_i lambda_j, whose normalized volume integrals
// are 1/4 and 1/20. No cached source column or source quadrature is read here.
fn uniform_volume_load(space: &SimplexSpace<'_>) -> Vec<f64> {
    assert_eq!(space.r, 2);
    assert_eq!(space.complex.tets.len(), 6);
    let mut load = vec![0.0; space.ndof];
    for tet in &space.complex.tets {
        for &vertex in tet {
            load[vertex as usize] += 1.0 / 24.0;
        }
        for i in 0..4 {
            for j in i + 1..4 {
                let mut edge = [tet[i], tet[j]];
                edge.sort_unstable();
                let index = space.complex.edges.binary_search(&edge).unwrap();
                load[space.edge_off + index] += 1.0 / 120.0;
            }
        }
    }
    load
}

// Integral lambda_0^a lambda_1^b lambda_2^c = 2 A a! b! c!/(a+b+c+2)!.
fn triangle_moment(area: f64, powers: [usize; 3]) -> f64 {
    const FACTORIAL: [f64; 7] = [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0];
    2.0 * area * powers.iter().map(|&p| FACTORIAL[p]).product::<f64>()
        / FACTORIAL[powers.iter().sum::<usize>() + 2]
}

fn midside_load(space: &SimplexSpace<'_>, face: [u32; 3], area: f64) -> Vec<f64> {
    let mut load = vec![0.0; space.ndof];
    for (i, vertex) in face.into_iter().enumerate() {
        let mut powers = [1, 1, 0];
        powers[i] += 1;
        load[vertex as usize] = 4.0 * triangle_moment(area, powers);
    }
    for (i, j) in [(0, 1), (1, 2), (2, 0)] {
        let mut edge = [face[i], face[j]];
        edge.sort_unstable();
        let index = space.complex.edges.binary_search(&edge).unwrap();
        let mut powers = [1, 1, 0];
        powers[i] += 1;
        powers[j] += 1;
        load[space.edge_off + index] = 4.0 * triangle_moment(area, powers);
    }
    load
}

fn work(load: &[f64], pressure: &[C64], amplitude: C64) -> f64 {
    let effort = load
        .iter()
        .zip(pressure)
        .fold(C64::ZERO, |sum, (&weight, &p)| sum + p.scale(weight));
    0.5 * (effort.conj() * amplitude).re
}

#[test]
fn distributed_drives_match_the_original_full_fe_bem_block_and_signed_work() {
    let gas = gas();
    let (mesh, points) = kuhn_cube(1);
    let points: Vec<_> = points
        .iter()
        .map(|p| [p[0] * LENGTH, p[1] * WIDTH, p[2] * WIDTH])
        .collect();
    let fe = HelmholtzFe::new(&mesh, &points, 2, medium(&gas)).unwrap();
    let space = SimplexSpace::new(&mesh, 2);
    assert_eq!(space.ndof, fe.ndof());
    let ports = ports(&fe);
    let exterior = surface(&fe);
    let selected = [2, 3, 4];
    let gate = CancelGate::new();
    let model = RadiatingInterior::new(
        &fe,
        &ports,
        &selected,
        &exterior,
        &gas,
        None,
        budget(),
        &gate,
    )
    .unwrap();
    // This is an interior feed on a non-opening face. The declared exterior
    // remains rigid there; no direct exterior radiation of an actuator is added.
    let face = patch(&fe, 2, 0.0)[0];
    let boundary = model
        .interior()
        .quadratic_velocity_source(&fe, &[(face, MIDSIDE)], &mut || false)
        .unwrap();
    let volume = model
        .interior()
        .volume_velocity_source(
            &fe,
            4,
            &mut |_, _| Ok::<_, HeterogeneousError>(1.0 / VOLUME),
            &mut || false,
        )
        .unwrap();
    let reference = gas.density * gas.sound_speed / (WIDTH * WIDTH);
    // Unequal power coordinates expose either missing reference-root factor.
    let references = [0.25, 3.0, 2.0, 0.5, 4.0].map(|x| x * reference);
    let prepared = model
        .prepare(FREQUENCY, &references, &route(), &gate)
        .unwrap();
    let feeds = [C64::new(1e-6, -0.3e-6), C64::new(0.4e-6, 0.7e-6)];
    let velocity = C64::new(0.002, -0.001);
    let flow = C64::new(-0.6e-6, 0.5e-6);
    let actual = prepared
        .response_with_sources_and_tolerance(
            &feeds,
            &[(&boundary, velocity)],
            &[(&volume, flow)],
            1e-8,
            &gate,
        )
        .unwrap();

    let native = aperture_impedance(
        &exterior,
        &source_apertures(&fe, &ports, &exterior, &selected),
        &[core::f64::consts::TAU * FREQUENCY],
        Medium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        },
        &route(),
        &gate,
    )
    .unwrap();
    let areas: Vec<_> = ports
        .iter()
        .map(|port| {
            let faces: Vec<_> = port.traces.iter().map(|&(face, _)| face).collect();
            fe.boundary_load(&faces).unwrap().1
        })
        .collect();
    let z: Vec<_> = native.rows[0]
        .impedance
        .iter()
        .enumerate()
        .map(|(i, z)| {
            z.conj()
                .scale(1.0 / (areas[selected[i / 3]] * areas[selected[i % 3]]))
        })
        .collect();
    let volume_load = uniform_volume_load(&space);
    let boundary_load = midside_load(&space, face, 0.5 * LENGTH * WIDTH);
    let load: Vec<_> = volume_load
        .iter()
        .zip(&boundary_load)
        .map(|(&v, &b)| flow.scale(v) + velocity.scale(b))
        .collect();
    let expected = direct_block_with_load(&fe, &ports, &selected, &[0, 1], &feeds, &z, Some(&load));
    for (&actual, &expected) in actual.interior.amplitudes.iter().zip(&expected) {
        near(C64::new(actual[0], actual[1]), expected, 1.0);
    }
    for (r, &mode) in selected.iter().enumerate() {
        near(
            actual.interior.inward_flow_m3_s[mode],
            expected[fe.ndof() + r],
            1e-6,
        );
        near(
            actual.exterior_pressure_pa[r],
            actual.interior.pressure_pa[mode],
            1.0,
        );
    }
    for (i, &flow) in feeds.iter().enumerate() {
        near(actual.interior.inward_flow_m3_s[i], flow, 1e-6);
    }
    let port_loads = source_loads(&fe, &ports);
    let port_work = feeds
        .iter()
        .enumerate()
        .map(|(i, &q)| work(&port_loads[i], &expected, q))
        .sum::<f64>();
    let boundary_work = work(&boundary_load, &expected, velocity);
    let volume_work = work(&volume_load, &expected, flow);
    let mut radiation_work = 0.0;
    for r in 0..3 {
        for s in 0..3 {
            radiation_work +=
                0.5 * (expected[fe.ndof() + r].conj() * z[r * 3 + s] * expected[fe.ndof() + s]).re;
        }
    }
    assert!(radiation_work > 0.0);
    for (actual, expected) in [
        (actual.source_power_w, port_work),
        (actual.interior.boundary_source_power_w, boundary_work),
        (actual.interior.volume_source_power_w, volume_work),
        (
            actual.total_source_power_w(),
            port_work + boundary_work + volume_work,
        ),
        (actual.radiation_power_w, radiation_work),
        (actual.power_balance_residual_w, 0.0),
    ] {
        near(C64::from_re(actual), C64::from_re(expected), radiation_work);
    }
    assert!(actual.port_relative_residual < 1e-10);
}

#[test]
fn zero_port_flow_still_radiates_distributed_sources_and_coherent_pairs_cancel() {
    let gas = gas();
    let (mesh, points) = kuhn_cube(1);
    let points: Vec<_> = points
        .iter()
        .map(|p| [p[0] * LENGTH, p[1] * WIDTH, p[2] * WIDTH])
        .collect();
    let fe = HelmholtzFe::new(&mesh, &points, 2, medium(&gas)).unwrap();
    let ports = ports(&fe);
    let exterior = surface(&fe);
    let gate = CancelGate::new();
    let model = RadiatingInterior::new(
        &fe,
        &ports,
        &[2, 3, 4],
        &exterior,
        &gas,
        None,
        budget(),
        &gate,
    )
    .unwrap();
    let face = patch(&fe, 2, 0.0)[0];
    let boundary = model
        .interior()
        .quadratic_velocity_source(&fe, &[(face, [0.0, 0.0, 0.0, 1.0, -1.0, 0.0])], &mut || {
            false
        })
        .unwrap();
    assert!(boundary.integrated_shape_m2().abs() < 1e-14 * LENGTH * WIDTH);
    let volume = model
        .interior()
        .volume_velocity_source(
            &fe,
            4,
            &mut |_, _| Ok::<_, HeterogeneousError>(1.0 / VOLUME),
            &mut || false,
        )
        .unwrap();
    assert!((volume.integrated_shape() - 1.0).abs() < 2e-13);
    let reference = gas.density * gas.sound_speed / (WIDTH * WIDTH);
    let prepared = model
        .prepare(FREQUENCY, &[reference; 5], &route(), &gate)
        .unwrap();
    let q = C64::new(1e-6, 0.3e-6);
    let v = C64::new(0.01, -0.004);
    let volume_only = prepared
        .response_with_sources(&[C64::ZERO; 2], &[], &[(&volume, q)], &gate)
        .unwrap();
    let boundary_only = prepared
        .response_with_sources(&[C64::ZERO; 2], &[(&boundary, v)], &[], &gate)
        .unwrap();
    assert_eq!(volume_only.interior.boundary_source_power_w, 0.0);
    assert_eq!(boundary_only.interior.volume_source_power_w, 0.0);
    for result in [&volume_only, &boundary_only] {
        assert_eq!(result.source_power_w, 0.0);
        assert!(result.radiation_power_w > 0.0);
        assert!(result.outward_flow_m3_s.iter().any(|q| q.abs() > 1e-14));
        near(
            C64::from_re(result.total_source_power_w()),
            C64::from_re(result.radiation_power_w),
            result.radiation_power_w,
        );
        near(
            C64::from_re(result.power_balance_residual_w),
            C64::ZERO,
            result.radiation_power_w,
        );
    }
    let dark = prepared
        .response_with_sources(
            &[C64::ZERO; 2],
            &[(&boundary, v), (&boundary, -v)],
            &[(&volume, q), (&volume, -q)],
            &gate,
        )
        .unwrap();
    assert!(dark.interior.amplitudes.iter().all(|p| *p == [0.0; 2]));
    assert!(dark.outward_flow_m3_s.iter().all(|q| *q == C64::ZERO));
    assert_eq!(dark.total_source_power_w(), 0.0);
    assert_eq!(dark.radiation_power_w, 0.0);
    assert_eq!(dark.power_balance_residual_w, 0.0);
}

#[test]
fn original_source_owner_support_and_cancellation_refusals_leave_factors_reusable() {
    let gas = gas();
    let (mesh, points) = kuhn_cube(1);
    let points: Vec<_> = points
        .iter()
        .map(|p| [p[0] * LENGTH, p[1] * WIDTH, p[2] * WIDTH])
        .collect();
    let fe = HelmholtzFe::new(&mesh, &points, 2, medium(&gas)).unwrap();
    let ports = ports(&fe);
    let exterior = surface(&fe);
    let gate = CancelGate::new();
    let model = RadiatingInterior::new(
        &fe,
        &ports,
        &[2, 3, 4],
        &exterior,
        &gas,
        None,
        budget(),
        &gate,
    )
    .unwrap();
    let other = RadiatingInterior::new(
        &fe,
        &ports,
        &[2, 3, 4],
        &exterior,
        &gas,
        None,
        budget(),
        &gate,
    )
    .unwrap();
    let foreign = other
        .interior()
        .volume_velocity_source(
            &fe,
            4,
            &mut |_, _| Ok::<_, HeterogeneousError>(1.0 / VOLUME),
            &mut || false,
        )
        .unwrap();
    let volume = model
        .interior()
        .volume_velocity_source(
            &fe,
            4,
            &mut |_, _| Ok::<_, HeterogeneousError>(1.0 / VOLUME),
            &mut || false,
        )
        .unwrap();
    // Geometrically identical FE substitutions still have a different owner.
    let other_fe = HelmholtzFe::new(&mesh, &points, 2, medium(&gas)).unwrap();
    let mut evaluations = 0;
    assert!(
        model
            .interior()
            .volume_velocity_source(
                &other_fe,
                4,
                &mut |_, _| {
                    evaluations += 1;
                    Ok::<_, HeterogeneousError>(1.0 / VOLUME)
                },
                &mut || false,
            )
            .is_err()
    );
    assert_eq!(evaluations, 0);
    for mode in [0, 2] {
        assert!(
            model
                .interior()
                .quadratic_velocity_source(&fe, &[(ports[mode].traces[0].0, MIDSIDE)], &mut || {
                    false
                },)
                .is_err()
        );
    }
    let reference = gas.density * gas.sound_speed / (WIDTH * WIDTH);
    let prepared = model
        .prepare(FREQUENCY, &[reference; 5], &route(), &gate)
        .unwrap();
    let sources = [(&volume, C64::new(1e-6, -0.2e-6))];
    let before = prepared
        .response_with_sources(&[C64::ZERO; 2], &[], &sources, &gate)
        .unwrap();
    // A zero requested amplitude cannot bypass original-owner admission.
    assert!(
        prepared
            .response_with_sources(&[C64::ZERO; 2], &[], &[(&foreign, C64::ZERO)], &gate,)
            .is_err()
    );
    assert!(
        prepared
            .response_with_sources(
                &[C64::ZERO; 2],
                &[],
                &[(&volume, C64::new(f64::NAN, 0.0))],
                &gate,
            )
            .is_err()
    );
    assert!(
        prepared
            .response_with_sources_and_tolerance(&[C64::ZERO; 2], &[], &sources, 0.0, &gate,)
            .is_err()
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        prepared.response_with_sources(&[C64::ZERO; 2], &[], &sources, &cancelled),
        Err(RadiationError::Cancelled)
    ));
    let after = prepared
        .response_with_sources(&[C64::ZERO; 2], &[], &sources, &gate)
        .unwrap();
    assert_eq!(after.interior.amplitudes, before.interior.amplitudes);
    assert_eq!(after.outward_flow_m3_s, before.outward_flow_m3_s);
    assert_eq!(after.total_source_power_w(), before.total_source_power_w());
    assert_eq!(after.radiation_power_w, before.radiation_power_w);
}
