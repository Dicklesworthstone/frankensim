//! Cane reed card + lay through the reduce lab (music bead
//! `frankensim-music-v8-root-3ez8g.6.2`) — "shaped like a reed, made of
//! cane": the reed-blank mesh + the CC-BY cane packs reduced offline to
//! the playing MSD parameters, driving the EXISTING `realize_reed_bore`
//! loop. The authored-parameter path remains the Estimate-only debug
//! image; this is the receipted upgrade.
//!
//! Provenance chain (the card carries all three): the blank-mesh digest,
//! the material resolution identity, and the reduction options. Every
//! cane value resolves from its compiled pack (bead pexfp.32): E_L dry and
//! its stated spread from the reed-HEEL pack
//! `cane-arundo-reedheel-mdpi-ma13204566` (the vamp IS heel-zone
//! material); the density is the BLANK pack's printed low value, a
//! DISCLOSED zone mismatch since no heel density is licensable; the loss
//! factor is delta/pi at the middle of the damping pack's printed log
//! decrement range (conversion ownership ours); the Poisson ratio is the
//! labelled authored fixture `cane-reed-authored-fixture`. The moisture
//! caveat travels: every value is DRY-state; the same heel pack carries the
//! water-soaked moduli, recorded context, not a model.
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
use fs_solid::reduce::{
    OrificeSpec, PressureLoadedValve, ValveCard, ValveCardRequest, reduce_pressure_loaded_valve,
    reduce_valve,
};

#[path = "support/material_packs.rs"]
mod material_packs;
use material_packs::{
    CANE_DAMPING, CANE_REED_BLANK, CANE_REED_FIXTURE, CANE_REED_HEEL, MaterialPack,
};

const RATE: u32 = 48_000;

/// The cane values, each resolved from its pack at its measured dry state.
struct Cane {
    /// Heel E_L dry [Pa] and its stated half-width (the pack's 1-sigma band).
    e_heel: f64,
    e_sigma: f64,
    /// The blank's printed low density [kg/m^3] (zone mismatch disclosed).
    rho: f64,
    /// delta/pi at the middle of the printed log decrement range.
    eta: f64,
    /// Authored-fixture Poisson ratio.
    nu: f64,
    /// `id hash` of the heel, blank, damping and fixture packs.
    provenance: [String; 4],
}

fn cane() -> Cane {
    let heel = MaterialPack::seed(CANE_REED_HEEL);
    let blank = MaterialPack::seed(CANE_REED_BLANK);
    let damping = MaterialPack::seed(CANE_DAMPING);
    let fixture = MaterialPack::fixture(CANE_REED_FIXTURE);
    let decrement = |property| {
        damping
            .scalar_as_measured(property, Dims::NONE)
            .expect("cane log decrement")
    };
    let (low, high) = (
        decrement("logarithmic_decrement_printed_range_low"),
        decrement("logarithmic_decrement_printed_range_high"),
    );
    Cane {
        e_heel: heel
            .scalar_as_measured("young_modulus_longitudinal", Pressure::DIMS)
            .expect("heel E_L dry"),
        e_sigma: heel
            .half_width("young_modulus_longitudinal")
            .expect("heel E_L spread"),
        rho: blank
            .scalar_as_measured("density_printed_range_low", fs_qty::Density::DIMS)
            .expect("blank low density"),
        eta: 0.5 * (low + high) / core::f64::consts::PI,
        nu: fixture
            .scalar_as_measured("poisson_ratio", Dims::NONE)
            .expect("authored Poisson ratio"),
        provenance: [
            heel.provenance(),
            blank.provenance(),
            damping.provenance(),
            fixture.provenance(),
        ],
    }
}

// Preserve this benchmark's former loading assumption explicitly. This is
// an authored effective face length, NOT the mesh's pressure-load projection.
const RESEARCH_PRESSURE_FACE_LENGTH_M: f64 = 0.025;

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
    let cane = cane();
    let [heel, blank, _, fixture] = &cane.provenance;
    let mut claims = ClaimSet::new();
    for (name, dims, value, source) in [
        (
            "density",
            fs_qty::Density::DIMS,
            cane.rho,
            format!(
                "{blank} density_printed_range_low (DISCLOSED zone mismatch: blank/cortex-weighted \
                 value used for heel-zone reduction; no heel density is licensable); dry, moisture \
                 unstated"
            ),
        ),
        (
            "young_modulus",
            Pressure::DIMS,
            e_pa,
            format!(
                "{heel} young_modulus_longitudinal (DRY; the same pack carries the water-soaked \
                 modulus - recorded context); isotropic reduction uses E_L, the bending-dominant \
                 axis; the 10x E_T anisotropy is unrepresented (disclosed)"
            ),
        ),
        (
            "poisson_ratio",
            Dims::NONE,
            cane.nu,
            format!(
                "{fixture} poisson_ratio (authored test fixture; no licensable cane Poisson ratio; \
                 bending-dominated reduction is insensitive)"
            ),
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
    mint_pressure_card(lay, e_pa).card().clone()
}

fn mint_pressure_card(lay: LayChart, e_pa: f64) -> PressureLoadedValve {
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
    // Uniform pressure on the complete planar underside, independently of
    // the tip strip used to observe opening. P1 triangle load = area/3 per
    // vertex. Positive vectors point toward opening (+z); closing pressure
    // in the runtime acts against them. The fixture's underside is z=0.
    let mut pressure_load = vec![[0.0; 3]; nodes.len()];
    for tet in &tets {
        for face in [
            [tet[0], tet[1], tet[2]],
            [tet[0], tet[1], tet[3]],
            [tet[0], tet[2], tet[3]],
            [tet[1], tet[2], tet[3]],
        ] {
            if face.iter().all(|&i| nodes[i][2] == 0.0) {
                let [a, b, c] = face.map(|i| nodes[i]);
                let area =
                    0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs();
                for i in face {
                    pressure_load[i][2] += area / 3.0;
                }
            }
        }
    }
    let loaded_area: f64 = pressure_load.iter().map(|v| v[2]).sum();
    assert!((loaded_area / (l * w) - 1.0).abs() < 1e-12);
    with_cx(|cx| {
        reduce_pressure_loaded_valve(
            &ValveCardRequest {
                element_order: fs_solid::reduce::ValveElementOrder::Tet4,
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
                loss_factor: cane().eta,
                source_id: "music/cane-reed-vamp/20x13mm-taper1.0-0.1/v1",
            },
            &pressure_load,
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
#[allow(clippy::items_after_statements)] // stage-local imports document provenance
fn render(card: &ValveCard, blowing_pa: f64, seconds: f64) -> (f64, f64, f64) {
    let n = (seconds * f64::from(RATE)) as usize;
    let reed = card
        .beating_reed(
            card.width_m * RESEARCH_PRESSURE_FACE_LENGTH_M,
            blowing_pa,
            0.02,
        )
        .expect("declared pressure face");
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
    // Oscillation is the AC part. Under steady flow the open bore holds its
    // static Poiseuille pressure drop (about 1.7 Pa for the open lay at 80 Pa),
    // which a total rms would mistake for a sounding reed.
    let ac_rms = |s: &[f64]| {
        let mean = s.iter().sum::<f64>() / s.len() as f64;
        (s.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / s.len() as f64).sqrt()
    };
    let rms = ac_rms(tail);
    let mid_rms = ac_rms(mid);
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
    let loaded = mint_pressure_card(lay, cane().e_heel);
    let card = loaded.card();
    let projected_area = loaded.effective_areas_m2()[0];
    assert!(projected_area.is_finite() && projected_area > 0.0);
    assert_ne!(loaded.identity(), card.identity);
    let projected_reed = loaded.beating_reed(100.0, 0.01).expect("projected reed");
    let recovered = projected_reed.stiffness_n_m / projected_reed.closing_pressure_pa
        * projected_reed.rest_opening_m;
    assert!((recovered / projected_area - 1.0).abs() < 1e-14);
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
    let reed = card
        .beating_reed(card.width_m * RESEARCH_PRESSURE_FACE_LENGTH_M, 2000.0, 0.01)
        .expect("declared pressure face");
    assert_eq!(reed.damping_ratio, 0.5 * cane().eta);
    let playing_damping = 2.0 * reed.damping_ratio * (reed.stiffness_n_m * reed.mass_kg).sqrt();
    assert!(
        (playing_damping - card.modes[0].damping_n_s_m).abs()
            <= 1e-14 * card.modes[0].damping_n_s_m
    );
    // The explicit pressure area must reach the existing moving-aperture
    // junction, not merely change a card field. Hold material, structural
    // mechanics, slit width and excitation fixed while varying that area.
    let pressure_trace = |area_scale| {
        let input = card
            .beating_reed(area_scale * projected_area, 100.0, 0.01)
            .unwrap();
        let output = realize_reed_bore(
            &clarinet_bore(),
            &air20(),
            input,
            Termination::IdealOpen,
            &mut PlateBank::default(),
            1.0,
            RATE,
            512,
            None,
        )
        .expect("physical pressure-area sweep");
        assert!(output.iter().all(|p| p.is_finite()));
        output
    };
    let baseline = pressure_trace(1.0);
    let larger = pressure_trace(2.0);
    assert!(
        baseline
            .iter()
            .zip(larger)
            .any(|(a, b)| (a - b).abs() > 1e-8)
    );
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
    let soft = mint_card(lay, cane().e_heel - cane().e_sigma);
    let hard = mint_card(lay, cane().e_heel + cane().e_sigma);
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
    // Timbre, measured two ways:
    // - MATCHED RELATIVE drive (1.5x each threshold): the waveform
    //   shape is nearly drive-scale-invariant, so the centroids agree
    //   closely - recorded, not asserted as a large difference.
    // - MATCHED ABSOLUTE pressure (same breath, different reed) at
    //   p_ref = 1.6 * p_hard: under the closed-loop moving-aperture
    //   junction (commits 4b2c3a40d..da2da0935c, honoring card damping
    //   zeta = 0.5 * eta, about 0.0127, and midpoint opening/swept flow),
    //   the softer reed has a lower closing pressure P_c proportional
    //   to E, so at fixed absolute mouth pressure p_ref it operates at
    //   a 1.33x higher fraction of P_c, closes more deeply each cycle,
    //   and generates stronger upper harmonics (measured: soft 606 Hz
    //   vs hard 524 Hz, shift = -0.136).
    let (_, _, c_soft_rel) = render(&soft, 1.5 * p_soft, 0.5);
    let (_, _, c_hard_rel) = render(&hard, 1.5 * p_hard, 0.5);
    let p_ref = 1.6 * p_hard;
    let (rms_soft_abs, _, c_soft_abs) = render(&soft, p_ref, 0.5);
    let (rms_hard_abs, _, c_hard_abs) = render(&hard, p_ref, 0.5);
    assert!(
        rms_soft_abs > 1.0e-4 && rms_hard_abs > 1.0e-4,
        "both reeds speak at the reference pressure"
    );
    let abs_shift = (c_soft_abs - c_hard_abs) / c_hard_abs;
    assert!(
        abs_shift > 0.10,
        "at fixed absolute pressure below closing, the SOFTER reed sits at a higher fraction of P_c and must read brighter \
         (soft {c_soft_abs:.0} vs hard {c_hard_abs:.0} Hz, shift {abs_shift:.3})"
    );
    // The lay sweep: a more open facing raises the threshold.
    let open_lay = LayChart::try_new(1.3e-3, 0.019).expect("open lay");
    let open = mint_card(open_lay, cane().e_heel);
    let base = mint_card(lay, cane().e_heel);
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
                value: cane().rho,
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
    let state = resolve_cane(cane().e_heel);
    let refused = with_cx(|cx| {
        reduce_valve(
            &ValveCardRequest {
                element_order: fs_solid::reduce::ValveElementOrder::Tet4,
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
                loss_factor: cane().eta,
                source_id: "test/empty-face",
            },
            cx,
        )
        .is_err()
    });
    assert!(refused, "an empty orifice face must refuse");
    println!("{{\"suite\":\"fs-couple\",\"case\":\"cr-003-refusals\",\"verdict\":\"pass\"}}");
}
