use super::*;

#[test]
fn pressure_projection_preserves_virtual_work_and_mode_sign() {
    // Nonuniform z motion at two free nodes; node 2 is clamped.
    let phi = [2.0, 4.0];
    let free = [2, 5];
    let loads = [[0.0, 0.0, 0.3], [0.0, 0.0, 0.6], [0.0, 0.0, 99.0]];
    let mean = mean_face_translation(&phi, &free, &[0, 1], [0.0, 0.0, 1.0]);
    let area = pressure_area(&phi, &free, &loads, mean);
    assert!((area - 1.0).abs() < 1e-14);
    let (pressure, velocity) = (7.0, 0.2);
    let nodal_power = pressure * (0.3 * (2.0 / 3.0) + 0.6 * (4.0 / 3.0)) * velocity;
    assert!((pressure * area * velocity - nodal_power).abs() < 1e-14);
    assert_eq!(pressure_area(&[-2.0, -4.0], &free, &loads, -mean), area);
    let opposite = loads.map(|v| v.map(|x| -x));
    assert_eq!(pressure_area(&phi, &free, &opposite, mean), -area);
}

#[test]
fn face_mean_preserves_partial_constraints_and_coordinate_permutations() {
    // Two face nodes move 2 and 4 along z; both x components are fixed.
    // G1: the nodal mean is 3, not the sum 6. The old denominator
    // understated effective mass (and stiffness) by a factor of four.
    let z_mean = mean_face_translation(&[2.0, 4.0], &[2, 5], &[0, 1], [0.0, 0.0, 1.0]);
    assert_eq!(z_mean, 3.0);
    assert_eq!((1.0 / z_mean).powi(2), 1.0 / 9.0);
    // G3: relabel z as x, preserving the physical displacements.
    let x_mean = mean_face_translation(&[2.0, 4.0], &[0, 3], &[0, 1], [1.0, 0.0, 0.0]);
    assert_eq!(z_mean, x_mean);
    // Releasing a transverse zero component cannot change the mean.
    assert_eq!(
        mean_face_translation(&[0.0, 2.0, 4.0], &[0, 2, 5], &[0, 1], [0.0, 0.0, 1.0]),
        z_mean
    );
    // A fully clamped face node contributes zero, but still belongs to
    // the declared face. Non-face motion does not enter its mean.
    assert_eq!(
        mean_face_translation(&[2.0, 4.0, 99.0], &[2, 5, 11], &[0, 1, 2], [0.0, 0.0, 1.0]),
        2.0
    );
}
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{MaterialPropertySelection, resolve_isotropic_elastic_state_point};
use fs_qty::{Dims, Pressure};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 0,
                kernel_id: 77,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"fs-solid\",\"case\":\"{case}\",\"verdict\":\"{}\",\"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

/// Synthetic tissue-ish card (AUTHORED, Estimate: no tissue pack
/// exists in data/matdb — the recorded upgrade path).
fn tissue_card(e_pa: f64, rho: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", fs_qty::Density::DIMS, rho),
        ("young_modulus", Pressure::DIMS, e_pa),
        ("poisson_ratio", Dims::NONE, 0.30),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                provenance: Provenance {
                    source: format!("authored tissue-like {name} (Estimate)"),
                    license: "CC0-1.0".to_owned(),
                    artifact: None,
                },
                observations: Vec::new(),
            })
            .expect("claim");
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored-soft-tissue".to_owned(),
            phase: "solid".to_owned(),
            process: "synthetic".to_owned(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("card")
}

/// Freudenthal 6-tet-per-hex slab: `nx x ny x nz` hexes over
/// `l x w x t` metres. Returns (nodes, tets).
fn slab_mesh(
    l: f64,
    w: f64,
    t: f64,
    nx: usize,
    ny: usize,
    nz: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let node = |i: usize, j: usize, k: usize| -> usize { (i * (ny + 1) + j) * (nz + 1) + k };
    let mut nodes = Vec::new();
    for i in 0..=nx {
        for j in 0..=ny {
            for k in 0..=nz {
                nodes.push([
                    l * i as f64 / nx as f64,
                    w * j as f64 / ny as f64,
                    t * k as f64 / nz as f64,
                ]);
            }
        }
    }
    let mut tets = Vec::new();
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let v = [
                    node(i, j, k),
                    node(i + 1, j, k),
                    node(i, j + 1, k),
                    node(i + 1, j + 1, k),
                    node(i, j, k + 1),
                    node(i + 1, j, k + 1),
                    node(i, j + 1, k + 1),
                    node(i + 1, j + 1, k + 1),
                ];
                // Freudenthal decomposition around the 0-7 diagonal.
                for tet in [
                    [v[0], v[1], v[3], v[7]],
                    [v[0], v[3], v[2], v[7]],
                    [v[0], v[2], v[6], v[7]],
                    [v[0], v[6], v[4], v[7]],
                    [v[0], v[4], v[5], v[7]],
                    [v[0], v[5], v[1], v[7]],
                ] {
                    tets.push(tet);
                }
            }
        }
    }
    (nodes, tets)
}

/// Clamp the x = 0 face; the orifice face is the BOTTOM-surface
/// strip (z = 0) within one element column of the free tip — the
/// underside that faces the opposing plane.
fn clamp_and_tip(nodes: &[[f64; 3]], l: f64, strip_m: f64) -> (Vec<usize>, Vec<usize>) {
    let mut fixed = Vec::new();
    let mut tip = Vec::new();
    for (n, p) in nodes.iter().enumerate() {
        if p[0].abs() < 1e-12 {
            fixed.extend_from_slice(&[3 * n, 3 * n + 1, 3 * n + 2]);
        }
        if p[2].abs() < 1e-12 && p[0] > l - strip_m - 1e-12 {
            tip.push(n);
        }
    }
    (fixed, tip)
}

struct Fixture {
    nodes: Vec<[f64; 3]>,
    tets: Vec<[usize; 4]>,
    fixed: Vec<usize>,
    tip: Vec<usize>,
    l: f64,
    w: f64,
    t: f64,
}

fn cantilever_fixture() -> Fixture {
    let (l, w, t) = (0.010f64, 0.012f64, 0.004f64);
    let (nodes, tets) = slab_mesh(l, w, t, 6, 6, 3);
    let (fixed, tip) = clamp_and_tip(&nodes, l, l / 6.0);
    Fixture {
        nodes,
        tets,
        fixed,
        tip,
        l,
        w,
        t,
    }
}

fn resolve_tissue(e_pa: f64, rho: f64) -> fs_material::state_point::IsotropicElasticStatePoint {
    let card = tissue_card(e_pa, rho);
    let point = QueryPoint::new().with("T", 293.15).expect("point");
    resolve_isotropic_elastic_state_point(&card, &point, MaterialPropertySelection::SingleClaimOnly)
        .expect("resolve")
}

#[test]
fn rl_001_cantilever_oracle() {
    with_cx(|cx| {
        let fx = cantilever_fixture();
        let e_pa = 1.0e6;
        let rho = 1050.0;
        let state = resolve_tissue(e_pa, rho);
        let request = ValveCardRequest {
            element_order: ValveElementOrder::Tet4,
            nodes_m: &fx.nodes,
            tetrahedra: &fx.tets,
            fixed_dofs: &fx.fixed,
            material: &state,
            orifice: OrificeSpec {
                face_nodes: &fx.tip,
                opening_axis: [0.0, 0.0, 1.0],
                width_axis: [0.0, 1.0, 0.0],
                opposing_plane_offset_m: -0.002,
            },
            window_hz: (20.0, 2000.0),
            retained_compliance_floor: 0.5,
            loss_factor: 0.1,
            source_id: "lab/cantilever-slab/24x12x4mm/v1",
        };
        let card = reduce_valve(&request, cx).expect("reduce");
        for line in card.debug_lines() {
            println!("{line}");
        }
        // Euler-Bernoulli fundamental (thin-beam limit; the thick P1
        // slab deviates and the band is authored from measurement).
        let i_area = fx.w * fx.t * fx.t * fx.t / 12.0;
        let a_area = fx.w * fx.t;
        let f_eb = (1.875f64.powi(2) / core::f64::consts::TAU)
            * (e_pa * i_area / (rho * a_area * fx.l.powi(4))).sqrt();
        let f1 = card.modes[0].frequency_hz;
        let f_ratio = f1 / f_eb;
        // Tip-normalized cantilever fundamental: m_eff ~ m_total/4.
        let mass_ratio = card.modes[0].mass_kg / card.total_mass_kg;
        // k = lambda m identity (construction).
        let lam = (core::f64::consts::TAU * f1).powi(2);
        let k_identity = (card.modes[0].stiffness_n_m - lam * card.modes[0].mass_kg).abs()
            / card.modes[0].stiffness_n_m;
        // Static tip stiffness 3EI/L^3 vs the reduced stiffness.
        let k_static = 3.0 * e_pa * i_area / fx.l.powi(3);
        let k_ratio = card.modes[0].stiffness_n_m / k_static;
        // Orifice geometry measured from the mesh.
        // Bottom strip: gap = 0 - (-0.002); width = full w; channel
        // thickness = the one-element strip length along the flow.
        let geom_ok = (card.width_m - fx.w).abs() < 1e-12
            && (card.effective_thickness_m - fx.l / 6.0).abs() < 1e-9
            && (card.rest_gap_m - 0.002).abs() < 1e-9;
        println!(
            "{{\"suite\":\"fs-solid\",\"case\":\"rl-001-numbers\",\"f1_hz\":{f1:.2},\
                 \"f_eb_hz\":{f_eb:.2},\"f_ratio\":{f_ratio:.3},\"mass_ratio\":{mass_ratio:.4},\
                 \"k_ratio_vs_static\":{k_ratio:.3},\
                 \"retained\":{:.4}}}",
            card.retained_compliance_fraction
        );
        // AUTHORED bands from the executed measure run: the thick
        // slab (L/t = 2.5) sits ABOVE the slender Euler-Bernoulli
        // line (plate-strip stiffening E/(1-nu^2) is +4.8% on f, P1
        // discretization adds the rest; measured +11.4%), the
        // BOTTOM-STRIP-normalized effective mass sits above the
        // slender-beam tip-normalized 0.25 (the strip's mean
        // deflection is below the tip's; measured 0.355), the
        // reduced stiffness exceeds the slender 3EI/L^3 consistently
        // (measured 1.814 = f_ratio^2 * mass_ratio/0.25, the
        // self-consistency identity), and the retained low modes
        // carry essentially the whole opening compliance (0.9805).
        let pass = k_identity < 1e-12
            && geom_ok
            && f_ratio > 1.0
            && f_ratio < 1.35
            && mass_ratio > 0.30
            && mass_ratio < 0.42
            && k_ratio > 1.5
            && k_ratio < 2.2
            && card.retained_compliance_fraction > 0.95;
        verdict(
            "rl-001-cantilever-oracle",
            pass,
            &format!(
                "f1 {f1:.1} Hz vs EB {f_eb:.1} (ratio {f_ratio:.3}); m_eff/m \
                     {mass_ratio:.3}; k vs 3EI/L^3 ratio {k_ratio:.3}; k=lam*m residual \
                     {k_identity:.2e}; geometry {geom_ok}"
            ),
        );
    });
}

#[test]
fn rl_002_card_drives_the_valve_island() {
    with_cx(|cx| {
        // Lip-ish leaflet: the slab with a parabolic thickness taper
        // (an AUTHORED lip-shaped fixture, honestly labeled; a real
        // lip scan slots through the same seam).
        let fx = cantilever_fixture();
        let mut nodes = fx.nodes.clone();
        for p in &mut nodes {
            let s = p[0] / fx.l;
            p[2] *= 1.0 - 0.4 * s * s;
        }
        let state = resolve_tissue(0.8e6, 1050.0);
        let request = ValveCardRequest {
            nodes_m: &nodes,
            tetrahedra: &fx.tets,
            element_order: ValveElementOrder::Tet4,
            fixed_dofs: &fx.fixed,
            material: &state,
            orifice: OrificeSpec {
                face_nodes: &fx.tip,
                opening_axis: [0.0, 0.0, 1.0],
                width_axis: [0.0, 1.0, 0.0],
                opposing_plane_offset_m: -0.0015,
            },
            window_hz: (20.0, 2000.0),
            retained_compliance_floor: 0.5,
            loss_factor: 0.15,
            source_id: "lab/lip-tapered-slab/v1",
        };
        let card = reduce_valve(&request, cx).expect("reduce");
        // Smoke loop: the card's 1-DOF island under a steady blowing
        // pressure through the existing fs-phs machinery.
        let m0 = &card.modes[0];
        let msd = fs_phs::mass_spring_damper(m0.mass_kg, m0.stiffness_n_m, m0.damping_n_s_m)
            .expect("msd admits");
        let face_area = card.width_m * card.effective_thickness_m;
        let p_blow = 500.0f64;
        let rho_air = 1.204;
        let dt = 1.0 / 48_000.0;
        let mut x = vec![0.0f64, 0.0];
        let mut worst_defect = 0.0f64;
        let mut opening = card.rest_gap_m;
        let mut flow = 0.0f64;
        for _ in 0..24_000 {
            let force = p_blow * face_area;
            let record = fs_phs::step(&msd, &x, &[force], dt).expect("step");
            worst_defect = worst_defect.max(record.supply_defect().abs());
            x = record.x;
            opening = (card.rest_gap_m + x[0]).max(0.0);
            flow = fs_phs::bernoulli_volume_flow(card.width_m, opening, p_blow, rho_air);
        }
        // Quasi-static equilibrium of the island: q* = F/k.
        let q_star = p_blow * face_area / m0.stiffness_n_m;
        let q_err = (x[0] - q_star).abs() / q_star.abs().max(1e-30);
        let flow_expected = fs_phs::bernoulli_volume_flow(
            card.width_m,
            (card.rest_gap_m + q_star).max(0.0),
            p_blow,
            rho_air,
        );
        let flow_err = (flow - flow_expected).abs() / flow_expected.max(1e-30);
        let pass = q_err < 1e-3 && flow_err < 1e-3 && worst_defect < 1e-8 && opening > 0.0;
        verdict(
            "rl-002-card-drives-the-island",
            pass,
            &format!(
                "settled q {:.4e} vs q* {q_star:.4e} (rel {q_err:.2e}); flow rel \
                     {flow_err:.2e}; worst supply defect {worst_defect:.2e}; opening {opening:.4e}",
                x[0]
            ),
        );
    });
}

#[test]
fn rl_003_refusals() {
    with_cx(|cx| {
        let fx = cantilever_fixture();
        let state = resolve_tissue(1.0e6, 1050.0);
        let base = |face: &'static [usize]| ValveCardRequest {
            element_order: ValveElementOrder::Tet4,
            nodes_m: &fx.nodes,
            tetrahedra: &fx.tets,
            fixed_dofs: &fx.fixed,
            material: &state,
            orifice: OrificeSpec {
                face_nodes: face,
                opening_axis: [0.0, 0.0, 1.0],
                width_axis: [0.0, 1.0, 0.0],
                opposing_plane_offset_m: -0.002,
            },
            window_hz: (20.0, 2000.0),
            retained_compliance_floor: 0.5,
            loss_factor: 0.1,
            source_id: "lab/refusals/v1",
        };
        let empty = reduce_valve(&base(&[]), cx);
        let oob = reduce_valve(&base(&[usize::MAX]), cx);
        let duplicate = reduce_valve(&base(&[0, 0]), cx);
        assert!(matches!(duplicate, Err(ReduceError::Invalid { .. })));
        for width in [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [0.0, 1.0, 1.0]] {
            let mut oblique = base(&[0]);
            oblique.orifice.width_axis = width;
            assert!(matches!(
                reduce_valve(&oblique, cx),
                Err(ReduceError::Invalid {
                    what: "orifice opening and width axes must be perpendicular"
                })
            ));
        }
        for offset in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut nonfinite = base(&[0]);
            nonfinite.orifice.opposing_plane_offset_m = offset;
            assert!(matches!(
                reduce_valve(&nonfinite, cx),
                Err(ReduceError::Invalid {
                    what: "orifice opposing-plane offset must be finite"
                })
            ));
        }
        // Retained-compliance floor: a window that EXCLUDES the
        // compliant fundamental keeps only stiff high modes.
        let mut high = base(Box::leak(fx.tip.clone().into_boxed_slice()));
        high.window_hz = (2500.0, 6000.0);
        high.retained_compliance_floor = 0.5;
        let starved = reduce_valve(&high, cx);
        // A face behind the opposing plane has no positive gap.
        let mut behind = base(Box::leak(fx.tip.clone().into_boxed_slice()));
        behind.orifice.opposing_plane_offset_m = 1.0;
        let negative_gap = reduce_valve(&behind, cx);
        let pass = matches!(empty, Err(ReduceError::Invalid { .. }))
            && matches!(oob, Err(ReduceError::Invalid { .. }))
            && matches!(starved, Err(ReduceError::RetainedComplianceTooLow { .. }))
            && matches!(negative_gap, Err(ReduceError::Invalid { .. }));
        verdict(
            "rl-003-refusals",
            pass,
            &format!(
                "empty {} oob {} starved {} negative-gap {}",
                empty.is_err(),
                oob.is_err(),
                starved.is_err(),
                negative_gap.is_err()
            ),
        );
    });
}

#[test]
fn rl_004_card_round_trip_and_tamper() {
    with_cx(|cx| {
        let fx = cantilever_fixture();
        let state = resolve_tissue(1.0e6, 1050.0);
        let request = ValveCardRequest {
            element_order: ValveElementOrder::Tet4,
            nodes_m: &fx.nodes,
            tetrahedra: &fx.tets,
            fixed_dofs: &fx.fixed,
            material: &state,
            orifice: OrificeSpec {
                face_nodes: &fx.tip,
                opening_axis: [0.0, 0.0, 1.0],
                width_axis: [0.0, 1.0, 0.0],
                opposing_plane_offset_m: -0.002,
            },
            window_hz: (20.0, 2000.0),
            retained_compliance_floor: 0.5,
            loss_factor: 0.1,
            source_id: "lab/round-trip/v1",
        };
        let card = reduce_valve(&request, cx).expect("reduce");
        let bytes = card.to_canonical_bytes();
        let back = ValveCard::from_canonical_bytes(&bytes).expect("round trip");
        let identical = back == card;
        let mut empty = card.clone();
        empty.modes.clear();
        empty.identity = empty.recomputed_identity();
        assert!(matches!(
            ValveCard::from_canonical_bytes(&empty.to_canonical_bytes()),
            Err(ReduceError::Card {
                what: "mode count must be nonzero and fit the remaining input rows"
            })
        ));
        let mut suffixed = bytes.clone();
        suffixed.extend_from_slice(b"mode\t1\t1\t1\t1\t1\n");
        assert!(matches!(
            ValveCard::from_canonical_bytes(&suffixed),
            Err(ReduceError::Card {
                what: "trailing records after card identity"
            })
        ));
        // Tamper with one stiffness digit: the identity must refuse.
        let text = String::from_utf8(bytes).expect("utf8");
        for count in [card.modes.len() + 1, usize::MAX] {
            let oversized = text.replacen(
                &format!("modes\t{}\n", card.modes.len()),
                &format!("modes\t{count}\n"),
                1,
            );
            assert!(matches!(
                ValveCard::from_canonical_bytes(oversized.as_bytes()),
                Err(ReduceError::Card {
                    what: "mode count must be nonzero and fit the remaining input rows"
                })
            ));
        }
        let mut altered = text.clone();
        let pos = altered.find("rest_gap_m\t").expect("field");
        altered.replace_range(pos + 11..pos + 12, "9");
        let refused = ValveCard::from_canonical_bytes(altered.as_bytes());
        let area = 1.0e-4;
        let beating = card.beating_reed(area, 2000.0, 0.01).expect("reed");
        let wider = card
            .beating_reed(2.0 * area, 2000.0, 0.01)
            .expect("wider pressure face");
        assert_eq!(wider.closing_pressure_pa, 0.5 * beating.closing_pressure_pa);
        assert_eq!(wider.mass_kg, beating.mass_kg);
        assert_eq!(wider.stiffness_n_m, beating.stiffness_n_m);
        let recovered_area =
            beating.stiffness_n_m / beating.closing_pressure_pa * beating.rest_opening_m;
        assert!((recovered_area / area - 1.0).abs() < 1e-14);
        assert!(empty.beating_reed(area, 2000.0, 0.01).is_err());
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(card.beating_reed(bad, 2000.0, 0.01).is_err());
        }
        assert!(card.beating_reed(area, -1.0, 0.01).is_err());
        assert!(card.beating_reed(area, 2000.0, f64::NAN).is_err());
        let reed_ok = beating.mass_kg > 0.0
            && beating.stiffness_n_m > 0.0
            && beating.rest_opening_m == card.rest_gap_m
            && beating.damping_ratio == 0.5 * card.loss_factor;
        let pass = identical && matches!(refused, Err(ReduceError::Card { .. })) && reed_ok;
        verdict(
            "rl-004-round-trip-and-tamper",
            pass,
            &format!(
                "round-trip identical {identical}; tampered refused {}; BeatingReed \
                     mint ok {reed_ok}",
                refused.is_err()
            ),
        );
    });
}

#[test]
fn rl_005_tet10_cantilever_card() {
    use crate::linear3::tet10::{TET10_EDGES, promote_p1_to_p2};
    use std::collections::BTreeMap;
    let started = std::time::Instant::now();
    with_cx(|cx| {
        let state = resolve_tissue(1.0e6, 1050.0);
        let material = TetElasticMaterial::from_resolved_elastic_state(&state);
        let (l, w, t) = (0.010, 0.012, 0.004);
        let quadrature = TetQuadrature {
            stiffness_points_per_dir: 3,
            mass_points_per_dir: 4,
        };
        // The physical load is uniform unit total z force on x=L, with
        // exact P1/P2 trace integration. Nodal point-load averaging would
        // change the load under refinement and mask the element error.
        let load = |nodes: &[[f64; 3]], tets: &[[usize; 4]], quadratic: Option<&[[usize; 10]]>| {
            let mut faces = BTreeMap::<[usize; 3], (usize, [usize; 6])>::new();
            for (e, tet) in tets.iter().enumerate() {
                for local in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
                    let [a, b, c] = local;
                    let mut f = [tet[a], tet[b], tet[c], 0, 0, 0];
                    if let Some(q) = quadratic {
                        for (i, [u, v]) in [[a, b], [b, c], [a, c]].into_iter().enumerate() {
                            let edge = TET10_EDGES
                                .iter()
                                .position(|p| p.contains(&u) && p.contains(&v))
                                .unwrap();
                            f[3 + i] = q[e][4 + edge];
                        }
                    }
                    let mut key = [f[0], f[1], f[2]];
                    key.sort_unstable();
                    faces.entry(key).or_insert((0, f)).0 += 1;
                }
            }
            let mut full = vec![0.0; 3 * nodes.len()];
            for (_, (count, f)) in faces {
                if count != 1 || !f[..3].iter().all(|&i| (nodes[i][0] - l).abs() < 1e-12) {
                    continue;
                }
                let [a, b, c] = [nodes[f[0]], nodes[f[1]], nodes[f[2]]];
                let area =
                    ((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])).abs() / 2.0;
                let trace = if quadratic.is_some() {
                    &f[3..]
                } else {
                    &f[..3]
                };
                for &node in trace {
                    full[3 * node + 2] += area / (3.0 * w * t);
                }
            }
            assert!((full.iter().sum::<f64>() - 1.0).abs() < 1e-12);
            full
        };
        let compliance = |a: &crate::linear3::TetElasticAssembly, full: &[f64]| {
            let rhs: Vec<_> = a.free_dofs.iter().map(|&i| full[i]).collect();
            let symbolic =
                fs_sparse::SymbolicLdlt::analyze(&a.stiffness, fs_sparse::DirectOrdering::Amd)
                    .unwrap();
            let factor = symbolic
                .factor(&a.stiffness, &fs_sparse::LdltOptions::default())
                .unwrap();
            let u = factor.solve(&rhs);
            rhs.iter().zip(u).map(|(a, b)| a * b).sum::<f64>()
        };
        let mut ladder = Vec::new();
        let mut free_dofs = Vec::new();
        for r in [1, 2, 4] {
            let (p1, tets) = slab_mesh(l, w, t, 2 * r, 2 * r, r);
            let (nodes, elements) = promote_p1_to_p2(&p1, &tets, None, cx).unwrap();
            let corners: Vec<_> = elements.iter().map(|q| [q[0], q[1], q[2], q[3]]).collect();
            let (fixed, tip) = clamp_and_tip(&nodes, l, l / 6.0);
            let a = QuadraticTetProblem {
                nodes_m: &nodes,
                elements: &elements,
                materials: TetMaterialField::Uniform(&material),
                fixed_dofs: &fixed,
                budget: TetAssemblyBudget::standard(),
                quadrature,
            }
            .assemble(cx)
            .unwrap();
            ladder.push(compliance(&a, &load(&nodes, &corners, Some(&elements))));
            free_dofs.push(a.free_dofs.len());
            if r == 2 {
                let card = reduce_valve(
                    &ValveCardRequest {
                        nodes_m: &nodes,
                        tetrahedra: &corners,
                        element_order: ValveElementOrder::Tet10 {
                            elements: &elements,
                            quadrature,
                        },
                        fixed_dofs: &fixed,
                        material: &state,
                        orifice: OrificeSpec {
                            face_nodes: &tip,
                            opening_axis: [0.0, 0.0, 1.0],
                            width_axis: [0.0, 1.0, 0.0],
                            opposing_plane_offset_m: -0.002,
                        },
                        window_hz: (20.0, 450.0),
                        retained_compliance_floor: 0.1,
                        loss_factor: 0.1,
                        source_id: "analytic/rl001-slab/tet10",
                    },
                    cx,
                )
                .unwrap();
                assert_eq!(card.discretization, ValveDiscretization::Tet10(quadrature));
                assert_eq!(
                    ValveCard::from_canonical_bytes(&card.to_canonical_bytes()).unwrap(),
                    card
                );
                let mut changed = card.clone();
                changed.discretization = ValveDiscretization::Tet4;
                assert_ne!(changed.recomputed_identity(), card.identity);
                let mut bytes = card.to_canonical_bytes();
                bytes[25] = b'1';
                assert!(ValveCard::from_canonical_bytes(&bytes).is_err());
            }
        }
        // P1 grid [8,8,4] and TET10 grid [4,4,2] have identical DOF counts.
        let (nodes, tets) = slab_mesh(l, w, t, 8, 8, 4);
        let (fixed, _) = clamp_and_tip(&nodes, l, l / 6.0);
        let p1 = TetLinearElasticProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &fixed,
            budget: TetAssemblyBudget::standard(),
        }
        .assemble(cx)
        .unwrap();
        assert_eq!(p1.free_dofs.len(), free_dofs[1]);
        let linear = compliance(&p1, &load(&nodes, &tets, None));
        let e1 = (linear - ladder[2]).abs() / ladder[2];
        let e2 = (ladder[1] - ladder[2]).abs() / ladder[2];
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.3-W1a\",\"test\":\"rl_005_tet10_cantilever_card\",\"tier\":\"G2\",\"inputs\":{{\"fixture\":\"rl001-slab\",\"element\":\"tet10\",\"quadrature_n\":[3,4]}},\"provenance\":\"analytic slab; fine TET10 static reference\",\"rungs\":{{\"axis\":\"h\",\"values\":[1.0,0.5,0.25],\"dofs\":{free_dofs:?}}},\"observed\":{ladder:?},\"order\":{{\"comparison\":\"matched DOFs\"}},\"p1_error\":{e1},\"p2_error\":{e2},\"tolerance\":0.2,\"elapsed_s\":{}}}",
            started.elapsed().as_secs_f64()
        );
        assert!(
            e2 * 5.0 <= e1,
            "DEF-08 improvement must be >=5x: P1 {e1:e}, P2 {e2:e}"
        );
    });
}

#[path = "tests/thermal.rs"]
mod thermal;
