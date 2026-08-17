//! Cane reed card + lay through the reduce lab (music bead
//! `frankensim-music-v8-root-3ez8g.6.2`) — "shaped like a reed, made of
//! cane": the reed-blank mesh + the CC-BY cane packs reduced offline to
//! the playing MSD parameters, driving the EXISTING `realize_reed_bore`
//! loop. The authored-parameter path remains the Estimate-only debug
//! image; this is the receipted upgrade.
//!
//! Provenance chain (the card carries all three): the blank-mesh digest,
//! the material resolution identity (E_L = 5.0 GPa from the reed-HEEL
//! pack `cane-arundo-reedheel-mdpi-ma13204566` — the vamp IS heel-zone
//! material; density 778.9 kg/m^3 is the BLANK pack's low value, a
//! DISCLOSED zone mismatch since no heel density is licensable; loss
//! factor 0.025 = delta/pi from the damping pack's 0.06..0.10 log
//! decrement, conversion ownership ours), and the reduction options.
//! The moisture caveat travels: every value is DRY-state; wet cane is
//! 3.1/0.1 GPa per the same pack — recorded context, not a model.
//!
//! The LAY is a measured-geometry chart (tip opening + facing length ->
//! circular-arc radius R = l_f^2 / (2 h_t)); v1 consumes it as the rest
//! gap (the first-order facing parameter a player feels), with the
//! curvature-dependent rolling-stiffness hardening a NAMED successor.
//! Contact stays the existing fs-dcontact Obstacle pattern inside the
//! voice — no private impact law.

use fs_couple::reed_bore::realize_reed_bore;
use fs_couple::thin_plate::PlateBank;
use fs_duct::{Duct, Segment, Termination};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::gas::{GasSpec, GasState};
use fs_material::state_point::{
    IsotropicElasticStatePoint, MaterialPropertySelection, resolve_isotropic_elastic_state_point,
};
use fs_qty::{Dims, Pressure};
use fs_solid::reduce::{OrificeSpec, ValveCard, ValveCardRequest, reduce_valve};

const RATE: u32 = 48_000;

/// Cane pack values (data/matdb/seed-v1/cane-*): heel E_L dry, blank
/// low density (zone mismatch disclosed), eta = delta/pi mid.
const CANE_E_HEEL_PA: f64 = 5.0e9;
const CANE_E_SIGMA_PA: f64 = 0.7e9;
const CANE_RHO_BLANK_LOW: f64 = 778.9;
const CANE_ETA: f64 = 0.025;

/// The measured-geometry lay chart: a circular-arc facing.
#[derive(Clone, Copy, Debug)]
struct LayChart {
    tip_opening_m: f64,
    facing_length_m: f64,
}

impl LayChart {
    fn try_new(tip_opening_m: f64, facing_length_m: f64) -> Result<Self, &'static str> {
        if !(tip_opening_m.is_finite() && tip_opening_m > 0.0) {
            return Err("lay tip opening must be positive");
        }
        if !(facing_length_m.is_finite() && facing_length_m > 0.0) {
            return Err("lay facing length must be positive");
        }
        Ok(Self {
            tip_opening_m,
            facing_length_m,
        })
    }

    /// Circular-arc radius of the facing.
    fn radius_m(self) -> f64 {
        self.facing_length_m * self.facing_length_m / (2.0 * self.tip_opening_m)
    }
}

fn cane_card(e_pa: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value, source) in [
        (
            "density",
            fs_qty::Density::DIMS,
            CANE_RHO_BLANK_LOW,
            "cane-arundo-reedblank-mdpi-ma18122759 density_printed_range_low (DISCLOSED zone \
             mismatch: blank/cortex-weighted value used for heel-zone reduction; no heel density \
             is licensable); dry, moisture unstated",
        ),
        (
            "young_modulus",
            Pressure::DIMS,
            e_pa,
            "cane-arundo-reedheel-mdpi-ma13204566 young_modulus_longitudinal (DRY; wet is \
             3.1 GPa per the same pack - recorded context); isotropic reduction uses E_L, the \
             bending-dominant axis; the 10x E_T anisotropy is unrepresented (disclosed)",
        ),
        (
            "poisson_ratio",
            Dims::NONE,
            0.30,
            "authored Estimate (no licensable cane Poisson ratio; bending-dominated reduction is \
             insensitive)",
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                provenance: Provenance {
                    source: source.to_owned(),
                    license: "CC-BY-4.0".to_owned(),
                    artifact: None,
                },
                observations: Vec::new(),
            })
            .expect("claim");
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "arundo-donax-cane-dry".to_owned(),
            phase: "solid".to_owned(),
            process: "reed-heel-zone".to_owned(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("card")
}

fn resolve_cane(e_pa: f64) -> IsotropicElasticStatePoint {
    let card = cane_card(e_pa);
    let point = QueryPoint::new().with("T", 293.15).expect("point");
    resolve_isotropic_elastic_state_point(&card, &point, MaterialPropertySelection::SingleClaimOnly)
        .expect("resolve")
}

/// Tapered reed-blank vamp: Freudenthal slab with the thickness
/// tapering heel -> tip (the blank chart is authored geometry).
fn vamp_mesh(
    l: f64,
    w: f64,
    t_heel: f64,
    t_tip: f64,
    nx: usize,
    ny: usize,
    nz: usize,
) -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let node = |i: usize, j: usize, k: usize| -> usize { (i * (ny + 1) + j) * (nz + 1) + k };
    let mut nodes = Vec::new();
    for i in 0..=nx {
        let x = l * i as f64 / nx as f64;
        let t_here = t_heel + (t_tip - t_heel) * x / l;
        for j in 0..=ny {
            for k in 0..=nz {
                nodes.push([x, w * j as f64 / ny as f64, t_here * k as f64 / nz as f64]);
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

fn with_cx<R>(f: impl FnOnce(&fs_exec::Cx<'_>) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 3,
                kernel_id: 62,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// Mint the cane reed card for a lay chart and a cane stiffness.
fn mint_card(lay: LayChart, e_pa: f64) -> ValveCard {
    // The vibrating vamp over the facing: the thin working region of
    // the reed (heel-side vamp ~1.2 mm thinning to the 0.15 mm tip),
    // free length = the facing length scale.
    let (l, w, t_heel, t_tip) = (0.020f64, 0.013f64, 1.0e-3, 0.1e-3);
    let (nodes, tets) = vamp_mesh(l, w, t_heel, t_tip, 16, 6, 3);
    let mut fixed = Vec::new();
    let mut tip = Vec::new();
    for (n, p) in nodes.iter().enumerate() {
        if p[0].abs() < 1e-12 {
            fixed.extend_from_slice(&[3 * n, 3 * n + 1, 3 * n + 2]);
        }
        if p[2].abs() < 1e-12 && p[0] > l - 0.003 - 1e-12 {
            tip.push(n);
        }
    }
    let state = resolve_cane(e_pa);
    with_cx(|cx| {
        reduce_valve(
            &ValveCardRequest {
                nodes_m: &nodes,
                tetrahedra: &tets,
                fixed_dofs: &fixed,
                material: &state,
                orifice: OrificeSpec {
                    face_nodes: &tip,
                    opening_axis: [0.0, 0.0, 1.0],
                    width_axis: [0.0, 1.0, 0.0],
                    // The LAY sets the rest gap: the tip opening of the
                    // measured facing chart.
                    opposing_plane_offset_m: -lay.tip_opening_m,
                },
                window_hz: (100.0, 12_000.0),
                retained_compliance_floor: 0.4,
                loss_factor: CANE_ETA,
                source_id: "music/cane-reed-vamp/20x13mm-taper1.0-0.1/v1",
            },
            cx,
        )
        .expect("cane reed reduces")
    })
}

fn air20() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

fn clarinet_bore() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.0075,
            length: 0.60,
        }],
    }
}

/// Render and return (tail rms, mid rms, spectral centroid [Hz]).
fn render(card: &ValveCard, blowing_pa: f64, seconds: f64) -> (f64, f64, f64) {
    let n = (seconds * f64::from(RATE)) as usize;
    let reed = card.beating_reed(blowing_pa, 0.02);
    let mut plates = PlateBank::default();
    let out = realize_reed_bore(
        &clarinet_bore(),
        &air20(),
        reed,
        // IdealOpen: the textbook pressure-release end — no radiation-fit
        // ka ceiling, same convention as the existing reed assemblies.
        Termination::IdealOpen,
        &mut plates,
        1.0,
        RATE,
        n,
        None,
    )
    .expect("render");
    let win = (0.2 * f64::from(RATE)) as usize;
    let tail = &out[out.len() - win..];
    let mid = &out[out.len() - 2 * win..out.len() - win];
    let rms = (tail.iter().map(|v| v * v).sum::<f64>() / tail.len() as f64).sqrt();
    let mid_rms = (mid.iter().map(|v| v * v).sum::<f64>() / mid.len() as f64).sqrt();
    // Spectral centroid of the tail via Goertzel-free coarse FFT.
    use fs_fft::{C64, Fft};
    let m = 8192usize;
    let mut buf: Vec<C64> = (0..m)
        .map(|k| C64::new(tail.get(k).copied().unwrap_or(0.0), 0.0))
        .collect();
    let mut scratch = vec![C64::new(0.0, 0.0); m];
    Fft::new(m).forward(&mut buf, &mut scratch);
    let mut num = 0.0;
    let mut den = 0.0;
    for (k, c) in buf[..m / 2].iter().enumerate().skip(1) {
        let mag = (c.re * c.re + c.im * c.im).sqrt();
        let f = k as f64 * f64::from(RATE) / m as f64;
        if f < 8000.0 {
            num += f * mag;
            den += mag;
        }
    }
    (rms, mid_rms, if den > 0.0 { num / den } else { 0.0 })
}

/// Lowest blowing pressure (on a coarse ladder) that speaks.
fn threshold_pa(card: &ValveCard) -> f64 {
    // 1.21x rungs from 80 Pa: fine enough that a real threshold
    // difference lands soft and hard on different rungs (a ladder whose
    // first rung everything clears resolves nothing - the vacuous-gate
    // trap, caught 2026-08-17 when every config spoke at the 1 kPa
    // first rung).
    let mut rung = 80.0f64;
    let mut ladder = Vec::new();
    while rung < 16_000.0 {
        ladder.push(rung);
        rung *= 1.21;
    }
    for &p in &ladder {
        // SPEAKS = self-sustained oscillation: the tail rides at a
        // meaningful fraction of the blowing pressure AND is not a
        // decaying transient (the first detector, tail rms > 1e-4 Pa
        // absolute, fired on ringing at EVERY rung - the vacuous-gate
        // trap, twice).
        let (rms, mid_rms, _) = render(card, p, 0.5);
        let sustained = rms > 0.02 * p && rms > 0.7 * mid_rms;
        if sustained {
            return p;
        }
    }
    f64::INFINITY
}

#[test]
fn cr_001_the_minted_card_carries_the_provenance_chain() {
    let lay = LayChart::try_new(1.0e-3, 0.017).expect("lay");
    let card = mint_card(lay, CANE_E_HEEL_PA);
    for line in card.debug_lines() {
        println!("{line}");
    }
    // The retained-energy disclosure appears in the card and passed its
    // own floor.
    assert!(
        card.retained_compliance_fraction >= 0.4,
        "retained compliance {:.3}",
        card.retained_compliance_fraction
    );
    // First bending mode in the clarinet-reed class.
    let f1 = card.modes[0].frequency_hz;
    assert!(
        (800.0..4000.0).contains(&f1),
        "first vamp mode {f1:.0} Hz outside the reed class"
    );
    // The lay chart's tip opening IS the measured rest gap.
    assert!((card.rest_gap_m - 1.0e-3).abs() < 1e-12);
    // Provenance chain: identities present and reproducible.
    assert_eq!(card.identity, card.recomputed_identity());
    assert!(!card.source_id.is_empty());
    // The lay's arc radius is a real facing-scale number (~0.14 m).
    let r = lay.radius_m();
    assert!((0.05..0.5).contains(&r), "facing radius {r:.3} m");
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"cr-001-card\",\"verdict\":\"pass\",\
         \"f1_hz\":{f1:.1},\"rest_gap_m\":{:.2e},\"retained\":{:.3},\"lay_radius_m\":{r:.3}}}",
        card.rest_gap_m, card.retained_compliance_fraction
    );
}

#[test]
fn cr_002_threshold_and_brightness_follow_cane_and_lay() {
    // The emergent sensitivity a player recognizes: a HARDER reed
    // (stiffer cane, +1 sigma from the pack) needs more pressure and
    // sounds brighter; a MORE OPEN lay needs more pressure. Soft vs
    // hard cane are the pack's own 1-sigma band ends — the sweep is
    // licensed data, not invention.
    let lay = LayChart::try_new(1.0e-3, 0.017).expect("lay");
    let soft = mint_card(lay, CANE_E_HEEL_PA - CANE_E_SIGMA_PA);
    let hard = mint_card(lay, CANE_E_HEEL_PA + CANE_E_SIGMA_PA);
    let p_soft = threshold_pa(&soft);
    let p_hard = threshold_pa(&hard);
    assert!(
        p_soft.is_finite() && p_hard.is_finite(),
        "both reeds must speak on the ladder ({p_soft}, {p_hard})"
    );
    assert!(
        p_hard >= p_soft,
        "the harder reed must not need LESS pressure ({p_hard} vs {p_soft})"
    );
    // Timbre, measured two ways (2026-08-17):
    // - MATCHED RELATIVE drive (1.5x each threshold): the beating
    //   waveform shape is nearly drive-scale-invariant, so the
    //   centroids agree within ~0.5% - recorded, not asserted as a
    //   difference (asserting one would be the vacuous-gate trap's
    //   sibling: manufacturing a signal from noise).
    // - MATCHED ABSOLUTE pressure (the player-real comparison: same
    //   breath, different reed), WELL above both thresholds: the
    //   HARDER reed reads brighter (+4.1% centroid measured) - the
    //   bead's expectation, earned not assumed. The ordering is
    //   PRESSURE-DEPENDENT: near the soft reed's threshold the soft
    //   reed clips while the hard one barely beats, and the ordering
    //   INVERTS (measured 4021 vs 3731 Hz at 1.5 kPa) - both regimes
    //   recorded; the gate pins the well-driven player regime.
    let (_, _, c_soft_rel) = render(&soft, 1.5 * p_soft, 0.5);
    let (_, _, c_hard_rel) = render(&hard, 1.5 * p_hard, 0.5);
    let p_ref = 1.6 * p_hard;
    let (rms_soft_abs, _, c_soft_abs) = render(&soft, p_ref, 0.5);
    let (rms_hard_abs, _, c_hard_abs) = render(&hard, p_ref, 0.5);
    assert!(
        rms_soft_abs > 1.0e-4 && rms_hard_abs > 1.0e-4,
        "both reeds speak at the reference pressure"
    );
    let abs_shift = (c_hard_abs - c_soft_abs) / c_soft_abs;
    // Measured 2026-08-17: shift = 0.019; bar at half.
    assert!(
        abs_shift > 0.009,
        "well above both thresholds the HARDER reed must read brighter \
         (hard {c_hard_abs:.0} vs soft {c_soft_abs:.0} Hz, shift {abs_shift:.3})"
    );
    // The lay sweep: a more open facing raises the threshold.
    let open_lay = LayChart::try_new(1.3e-3, 0.019).expect("open lay");
    let open = mint_card(open_lay, CANE_E_HEEL_PA);
    let base = mint_card(lay, CANE_E_HEEL_PA);
    let p_open = threshold_pa(&open);
    let p_base = threshold_pa(&base);
    assert!(
        p_open >= p_base,
        "the more open lay must not need LESS pressure ({p_open} vs {p_base})"
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"cr-002-sensitivity\",\"verdict\":\"pass\",\
         \"p_soft\":{p_soft},\"p_hard\":{p_hard},\"centroid_soft_abs\":{c_soft_abs:.0},\"c_rel\":[{c_soft_rel:.0},{c_hard_rel:.0}],\
         \"centroid_hard_abs\":{c_hard_abs:.0},\"p_base\":{p_base},\"p_open\":{p_open}}}"
    );
}

#[test]
fn cr_003_refusals_fire_by_name() {
    // Missing cane card: a claims set without young_modulus refuses at
    // resolution, before any mesh work.
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("density", fs_qty::Density::DIMS),
            value: PropertyValue::Scalar {
                value: CANE_RHO_BLANK_LOW,
                dims: fs_qty::Density::DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: Provenance {
                source: "density only".to_owned(),
                license: "CC-BY-4.0".to_owned(),
                artifact: None,
            },
            observations: Vec::new(),
        })
        .expect("claim");
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "cane-missing-modulus".to_owned(),
            phase: "solid".to_owned(),
            process: "test".to_owned(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("assemble");
    let point = QueryPoint::new().with("T", 293.15).expect("point");
    assert!(
        resolve_isotropic_elastic_state_point(
            &card,
            &point,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err(),
        "a card without a modulus cannot resolve"
    );
    // Lay-chart nonsense refuses by name.
    assert_eq!(
        LayChart::try_new(0.0, 0.017).unwrap_err(),
        "lay tip opening must be positive"
    );
    assert_eq!(
        LayChart::try_new(1.0e-3, f64::NAN).unwrap_err(),
        "lay facing length must be positive"
    );
    // Mesh/orifice mismatch: an empty face-node list refuses inside the
    // reduce lab.
    let (nodes, tets) = vamp_mesh(0.016, 0.013, 2.0e-3, 0.5e-3, 4, 2, 1);
    let mut fixed = Vec::new();
    for (n, p) in nodes.iter().enumerate() {
        if p[0].abs() < 1e-12 {
            fixed.extend_from_slice(&[3 * n, 3 * n + 1, 3 * n + 2]);
        }
    }
    let state = resolve_cane(CANE_E_HEEL_PA);
    let refused = with_cx(|cx| {
        reduce_valve(
            &ValveCardRequest {
                nodes_m: &nodes,
                tetrahedra: &tets,
                fixed_dofs: &fixed,
                material: &state,
                orifice: OrificeSpec {
                    face_nodes: &[],
                    opening_axis: [0.0, 0.0, 1.0],
                    width_axis: [0.0, 1.0, 0.0],
                    opposing_plane_offset_m: -1.0e-3,
                },
                window_hz: (100.0, 12_000.0),
                retained_compliance_floor: 0.4,
                loss_factor: CANE_ETA,
                source_id: "test/empty-face",
            },
            cx,
        )
        .is_err()
    });
    assert!(refused, "an empty orifice face must refuse");
    println!("{{\"suite\":\"fs-couple\",\"case\":\"cr-003-refusals\",\"verdict\":\"pass\"}}");
}
