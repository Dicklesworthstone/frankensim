use super::*;

fn options() -> OneWayShellObserverOptions {
    OneWayShellObserverOptions {
        receiver: ReceiverSpec {
            id: "front".into(),
            position_m: [0.0, 0.0, 1.0],
        },
        minimum_clearance_m: 0.05,
        pass_band_upper_hz: 4000.0,
        taper_stop_hz: 5000.0,
        filter_latency_samples: 512,
        fft_lengths: vec![512, 1024, 2048, 4096],
        tail_energy_tolerance: 1e-10,
        dtft_relative_tolerance: 1e-4,
        maximum_frequency_entries: 1 << 20,
        maximum_replay_terms: 1 << 22,
        panel_rule: PanelRule::Triangle { max_pairs: 1 << 24 },
        maximum_skin_panels: 8192,
        maximum_anchor_distance_m: 1e-6,
        maximum_boundary_terms: 1 << 46,
        workers: 16,
        rule: OneWayRule {
            maximum_loss_fraction: 0.1,
            maximum_shift_fraction: 0.1,
        },
    }
}

#[test]
fn the_declaration_admits_its_ceilings_and_refuses_one_past_each() {
    let rate = 144_000;
    options().validate(rate).unwrap();
    let edits: [(&str, fn(&mut OneWayShellObserverOptions)); 12] = [
        ("no workers", |o| o.workers = 0),
        ("257 workers", |o| o.workers = 257),
        ("past the dense panel cap", |o| {
            o.maximum_skin_panels = helmholtz::MAX_DENSE_PANELS + 1
        }),
        ("no boundary work", |o| o.maximum_boundary_terms = 0),
        ("past the boundary ceiling", |o| {
            o.maximum_boundary_terms = MAX_ONE_WAY_BOUNDARY_TERMS + 1
        }),
        ("negative loss fraction", |o| {
            o.rule.maximum_loss_fraction = -0.1
        }),
        ("nonfinite shift fraction", |o| {
            o.rule.maximum_shift_fraction = f64::NAN
        }),
        ("a tail past the shared 1e-10 cap", |o| {
            o.tail_energy_tolerance = 2e-10
        }),
        ("zero clearance", |o| o.minimum_clearance_m = 0.0),
        ("zero anchor distance", |o| {
            o.maximum_anchor_distance_m = 0.0
        }),
        ("an empty near-pair allowance", |o| {
            o.panel_rule = PanelRule::Triangle { max_pairs: 0 }
        }),
        ("a centroid transfer", |o| {
            o.panel_rule = PanelRule::Centroid
        }),
    ];
    for (what, edit) in edits {
        let mut changed = options();
        edit(&mut changed);
        assert!(changed.validate(rate).is_err(), "{what} was admitted");
    }
    let mut ceilings = options();
    ceilings.workers = 256;
    ceilings.maximum_skin_panels = helmholtz::MAX_DENSE_PANELS;
    ceilings.maximum_boundary_terms = MAX_ONE_WAY_BOUNDARY_TERMS;
    ceilings.validate(rate).unwrap();
}

/// An icosphere welded into one part whose single mode is a unit radial
/// breathing of every node, with the same vertices and triangles as its skin.
fn breathing(radius: f64, subdivisions: u32) -> (ShellMesh, Vec<Vec<f64>>, ShellAcousticSkin) {
    let sphere = SpherePanels::icosphere(radius, subdivisions).unwrap();
    let mut nodes: Vec<[f64; 3]> = Vec::new();
    let mut tris = Vec::new();
    for triangle in sphere.triangles().unwrap() {
        tris.push(triangle.map(|p| {
            nodes.iter().position(|q| *q == p).unwrap_or_else(|| {
                nodes.push(p);
                nodes.len() - 1
            })
        }));
    }
    let mut row = vec![0.0; 6 * nodes.len()];
    for (i, p) in nodes.iter().enumerate() {
        let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        for axis in 0..3 {
            row[6 * i + axis] = p[axis] / r;
        }
    }
    let mut anchors = vec![Vec::new(); nodes.len()];
    for (t, tri) in tris.iter().enumerate() {
        for (corner, &v) in tri.iter().enumerate() {
            if anchors[v].is_empty() {
                let mut barycentric = [0.0; 3];
                barycentric[corner] = 1.0;
                anchors[v].push(SkinAnchor {
                    weight: 1.0,
                    component: 0,
                    triangle: t,
                    barycentric,
                });
            }
        }
    }
    let skin = ShellAcousticSkin {
        vertices: nodes.clone(),
        anchors,
        triangles: tris.clone(),
        apertures: Vec::new(),
    };
    (ShellMesh { nodes, tris }, vec![row], skin)
}

#[test]
fn a_skin_on_its_own_part_carries_the_facet_trace_and_refuses_a_loose_binding() {
    let gate = CancelGate::new_clock_free();
    let (mesh, rows, skin) = breathing(0.1, 2);
    let parts = [(&mesh, rows.as_slice())];
    let traces = skin_traces(&parts, 1, &skin, 8192, 1e-12, &gate).unwrap();
    assert_eq!(traces.largest_anchor_distance_m, 0.0);
    for (t, tri) in mesh.tris.iter().enumerate() {
        let normal = traces.surface.normals()[t];
        let expected: f64 = (0..3)
            .map(|axis| {
                normal[axis] * (tri.iter().map(|&v| rows[0][6 * v + axis]).sum::<f64>() / 3.0)
            })
            .sum();
        assert_eq!(traces.modes[0][t], expected);
        assert!((0.97..=1.0).contains(&expected), "panel {t}: {expected}");
    }
    let refused = |skin: &ShellAcousticSkin, what: &str| {
        assert!(
            skin_traces(&parts, 1, skin, 8192, 1e-6, &gate).is_err(),
            "{what} was admitted"
        );
    };
    let mut moved = skin.clone();
    moved.vertices[0][0] += 2e-6;
    refused(&moved, "a vertex 2 um from its anchor");
    let mut loose = skin.clone();
    loose.anchors[0].clear();
    refused(&loose, "a structural vertex without an anchor");
    let mut heavy = skin.clone();
    heavy.anchors[0][0].weight = 1.0 + 1e-9;
    refused(&heavy, "anchor weights summing past one");
    let mut outside = skin.clone();
    outside.anchors[0][0].barycentric = [1.5, -0.5, 0.0];
    refused(&outside, "a barycentric point outside its triangle");
    let mut capped = skin.clone();
    capped.apertures = vec![SkinAperture {
        id: "a".into(),
        triangles: vec![0, 1],
        flow_share: 0.6,
    }];
    refused(&capped, "aperture shares summing to 0.6");
    capped.apertures.push(SkinAperture {
        id: "b".into(),
        triangles: vec![1],
        flow_share: 0.4,
    });
    refused(&capped, "a triangle in two apertures");
    capped.apertures[1].triangles = vec![2];
    let split = skin_traces(&parts, 1, &capped, 8192, 1e-6, &gate).unwrap();
    assert_eq!(
        split.modes[0][0], 0.0,
        "an aperture panel carries no structural motion"
    );
    let flow: f64 = split.apertures[0]
        .iter()
        .zip(split.surface.areas())
        .map(|(v, a)| v * a)
        .sum();
    assert!((flow - 1.0).abs() < 1e-12);
    assert_eq!(skin.triangles.len(), 320);
    assert!(
        skin_traces(&parts, 1, &skin, 319, 1e-6, &gate).is_err(),
        "320 panels past a 319-panel cap"
    );
}

/// A breathing sphere's radiation impedance under exp(-i w t), per unit normal
/// velocity: rho c A ((ka)^2 - i ka) / (1 + (ka)^2). Its resistance and its
/// added mass rho A a / (1 + (ka)^2) fix the receipt's signs.
#[test]
fn a_breathing_sphere_radiates_its_analytic_resistance_and_added_mass() {
    let gate = CancelGate::new_clock_free();
    let medium = Medium {
        density: 1.2,
        sound_speed: 343.0,
    };
    let radius = 0.1;
    let (mesh, rows, skin) = breathing(radius, 3);
    let parts = [(&mesh, rows.as_slice())];
    let traces = skin_traces(&parts, 1, &skin, 8192, 1e-12, &gate).unwrap();
    let areas = traces.surface.areas();
    let area: f64 = areas.iter().sum();
    let mean = traces.modes[0]
        .iter()
        .zip(areas)
        .map(|(t, a)| t * a)
        .sum::<f64>()
        / area;
    // ka 0.5 lies below the sphere's plain-CBIE limit; ka 2.5 lies above it, on ruled Burton-Miller.
    for ka in [0.5, 2.5] {
        let omega = ka * medium.sound_speed / radius;
        let (rows_z, limit) = impedance_rows(
            &traces.surface,
            &traces.modes,
            &[omega],
            medium,
            PanelRule::Triangle { max_pairs: 1 << 24 },
            4,
            &gate,
        )
        .unwrap();
        assert_eq!(ka < limit * radius, ka == 0.5);
        let z = rows_z[0][0];
        let scale = medium.density * medium.sound_speed * area * mean * mean / (1.0 + ka * ka);
        let (resistance, reactance) = (scale * ka * ka, -scale * ka);
        assert!(
            (z.re / resistance - 1.0).abs() < 0.03,
            "ka {ka}: R {} against {resistance}",
            z.re
        );
        assert!(
            (z.im / reactance - 1.0).abs() < 0.03,
            "ka {ka}: X {} against {reactance}",
            z.im
        );
        let bound = mode_bounds(
            &[omega],
            &[0.01],
            &rows_z,
            OneWayRule {
                maximum_loss_fraction: 0.1,
                maximum_shift_fraction: 0.1,
            },
        )
        .unwrap();
        let added_mass = -z.im / omega;
        assert!(
            added_mass > 0.0 && bound[0].added_mass_shift < 0.0,
            "air mass lowers the frequency"
        );
        assert_eq!(bound[0].added_mass_shift, -0.5 * added_mass);
        assert_eq!(bound[0].radiation_zeta, z.re / (2.0 * omega));
        assert_eq!(bound[0].nearest_spacing_hz, f64::INFINITY);
    }
}

#[test]
fn the_rule_keeps_a_lightly_radiating_mode_one_way_and_names_the_others() {
    let tau = core::f64::consts::TAU;
    let omega = [tau * 100.0, tau * 200.0, tau * 210.0];
    let wood = [0.01, 0.01, 0.01];
    let rule = OneWayRule {
        maximum_loss_fraction: 0.1,
        maximum_shift_fraction: 0.1,
    };
    // Mode 0: radiation 5 % of its wood loss and no reactance. Mode 1: 50 %. Mode 2: no
    // resistance, but a shift of 5e-3 f (1.05 Hz), past 10 % of its 4.2 Hz bandwidth.
    let z = vec![
        vec![
            C64::new(2.0 * omega[0] * 0.0005, 0.0),
            C64::new(1.0, 0.0),
            C64::ZERO,
        ],
        vec![C64::ZERO, C64::new(2.0 * omega[1] * 0.005, 0.0), C64::ZERO],
        vec![
            C64::ZERO,
            C64::new(0.0, 3.0),
            C64::new(0.0, 2.0 * 5e-3 * omega[2]),
        ],
    ];
    let bound = mode_bounds(&omega, &wood, &z, rule).unwrap();
    assert!(bound[0].one_way && !bound[1].one_way && !bound[2].one_way);
    assert!((bound[0].radiation_zeta - 0.0005).abs() < 1e-15);
    assert!(
        (bound[2].added_mass_shift - 5e-3).abs() < 1e-15,
        "{}",
        bound[2].added_mass_shift
    );
    assert!((bound[2].half_power_bandwidth_hz - 4.2).abs() < 1e-12);
    assert!((bound[1].nearest_spacing_hz - 10.0).abs() < 1e-9);
    assert!((bound[0].nearest_spacing_hz - 100.0).abs() < 1e-9);
    assert_eq!(bound[0].largest_mixing_mode, Some(1));
    let expected = omega[0] * 1.0 / (omega[1] * omega[1] - omega[0] * omega[0]);
    assert!((bound[0].largest_mixing - expected).abs() <= 1e-15 * expected);
    assert_eq!(bound[2].largest_mixing_mode, Some(1));
    let mut lossless = wood;
    lossless[0] = 0.0;
    assert!(
        !mode_bounds(&omega, &lossless, &z, rule).unwrap()[0].one_way,
        "a lossless mode cannot hide radiation"
    );
    assert!(mode_bounds(&omega, &wood, &z[..2], rule).is_err());
    assert!(
        mode_bounds(&[omega[0], omega[0], omega[2]], &wood, &z, rule).is_err(),
        "an exact degeneracy mixes without bound"
    );
}
