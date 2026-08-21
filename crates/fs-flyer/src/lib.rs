//! fs-flyer — Wright Flyer aircraft-level model (L4). Bead
//! frankensim-wf-root-guzez.4.1 (E3.1, Wright Flyer program).
//!
//! Spec: COMPREHENSIVE_PLAN_FOR_REAL_TIME_WRIGHT_FLYER_SIM_WITH_FRANKENSIM.md
//! (ROUND 6 steady state) §5.1; evidence: data/wright-flyer/
//! {flyer-reference,canard-mechanics,geometry-conventions,frame-conventions}.
//!
//! E3.1 scope: the `FlyerDesign` schema (Round-2 canard/mechanics/structure
//! fields), control-topology enums, admission with typed refusals, the
//! component mass/inertia build-up (cross-checked against the published
//! single-lineage Flyer inertias), and the derived-quantity panel (fixed
//! AND free-control static margins, canard volume, hinge-moment gradient).
//!
//! Conventions (frozen, E1.4/E1.9): frd-body-v1 — +x FORWARD from the wing
//! leading edge, +y right, +z down; SI units internally; angles in radians.
//! The derived panel's aerodynamic formulas are DOCUMENTED SIMPLE MODELS
//! (two-surface, no interference); the battery pins them against hand
//! calculations, and their divergence from the Culick lineage is recorded
//! data, never silently reconciled.

use fs_blake3::hash_domain;

pub mod abcompare;
pub mod adapter;
pub mod addedmass;
pub mod aerowarp;
pub mod aircraft;
pub mod assist;
pub mod augmented;
pub mod campaign;
pub mod canardmech;
pub mod checkpoint;
pub mod contact;
pub mod dragledger;
pub mod effectowners;
pub mod equilibrate;
pub mod fieldsvc;
pub mod freecontrol;
pub mod hcampaign;
pub mod hinference;
pub mod horchestrator;
pub mod longitudinal;
pub mod partitioned;
pub mod perception;
pub mod pilot;
pub mod prelaunch;
pub mod propcoupling;
pub mod rail;
pub mod referee;
pub mod refereeharness;
pub mod registryaudit;
pub mod replay;
pub mod replayenv;
pub mod sections;
pub mod simloop;
pub mod spine;
pub mod sweptevents;

/// Identity domain for design digests.
pub const DESIGN_DIGEST_DOMAIN: &str = "org.frankensim.fs-flyer.design.v1";

/// Admission caps (refusals at cap AND cap+1 per workspace law).
pub const MAX_SPAN_M: f64 = 20.0;
/// Chord cap [m].
pub const MAX_CHORD_M: f64 = 5.0;
/// Camber-ratio cap (matches fs-airfoil's admitted domain).
pub const MAX_CAMBER_RATIO: f64 = 0.15;
/// Pilot-mass domain [kg].
pub const MIN_PILOT_KG: f64 = 40.0;
/// Pilot-mass cap [kg].
pub const MAX_PILOT_KG: f64 = 120.0;
/// Component-count cap.
pub const MAX_COMPONENTS: usize = 64;
/// Biplane-area consistency band: |S_both − 2·b·c| / S_both must stay under
/// this (the real 1903 planform is 2.7% off the rectangular product).
pub const AREA_CONSISTENCY_TOL: f64 = 0.15;
/// Declared-vs-summed empty-mass tolerance [kg].
pub const MASS_SPEC_TOL_KG: f64 = 1.0;

/// Published whole-aircraft inertias (Jex & Culick AIAA 85-1804 lineage,
/// single-lineage per canard-mechanics-v1), converted slug·ft² → kg·m².
pub const PUBLISHED_IXX_KGM2: f64 = 1318.0 * 1.355_817_9;
/// Pitch inertia reference [kg·m²].
pub const PUBLISHED_IYY_KGM2: f64 = 271.0 * 1.355_817_9;
/// Yaw inertia reference [kg·m²].
pub const PUBLISHED_IZZ_KGM2: f64 = 1343.0 * 1.355_817_9;

/// A typed refusal (workspace law).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Refusal {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// Human-readable diagnosis.
    pub message: String,
    /// Ranked repairs, most likely fix first.
    pub ranked_repairs: Vec<String>,
}

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal {
        code,
        message,
        ranked_repairs: vec![repair.into()],
    }
}

/// 1903 vs later lateral-control wiring (plan Round-2 field).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LateralControlTopology {
    /// 1903: rudder SLAVED to warp (δr = ratio·δw); hips drive both.
    WarpWithSlavedRudder {
        /// Slaving ratio δr/δw (flyer-reference: 2.5).
        ratio: f64,
    },
    /// 1905+: independent rudder control.
    WarpIndependentRudder,
}

/// Warp structural mode (plan Round-2 field).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WarpStructureMode {
    /// Historical flexible trussed structure (loose aft spar, bias fabric).
    FlexibleTrussed,
    /// Counterfactual rigidized structure (warp disabled).
    Rigidized,
}

/// One mass component: six-bucket engine convention + airframe + crew.
/// Positions are frd-body-v1 metres from the wing leading edge; gyration
/// radii squared [m²] carry each component's own mass distribution.
#[derive(Clone, Debug, PartialEq)]
pub struct MassComponent {
    /// Stable component name (unique within a design).
    pub name: &'static str,
    /// Mass [kg].
    pub mass_kg: f64,
    /// Position [m], frd (+x forward from wing LE, +y right, +z down).
    pub position_m: [f64; 3],
    /// Squared gyration radii [m²] about the component's own CG, per body
    /// axis (rx² adds to Ixx, etc.). CALIBRATION provenance: reference-
    /// config radii are tuned to reproduce the published single-lineage
    /// inertias and say so; they are not measurements.
    pub gyration_sq_m2: [f64; 3],
}

/// Wing geometry (both planes; geometry-conventions-v1 symbols).
#[derive(Clone, Debug, PartialEq)]
pub struct WingGeometry {
    /// Full span b [m].
    pub span_m: f64,
    /// Chord c [m].
    pub chord_m: f64,
    /// Total biplane lifting area S_both [m²].
    pub area_both_m2: f64,
    /// Vertical gap between planes [m].
    pub gap_m: f64,
    /// Camber ratio f/c.
    pub camber_ratio: f64,
    /// Total anhedral tip droop [m].
    pub anhedral_droop_m: f64,
}

/// Canard geometry + mechanics (E1.5 dossier fields).
#[derive(Clone, Debug, PartialEq)]
pub struct CanardGeometry {
    /// Per-plane span [m].
    pub span_m: f64,
    /// Per-plane chord [m].
    pub chord_m: f64,
    /// Total area, both planes [m²].
    pub area_both_m2: f64,
    /// Gap between canard planes [m].
    pub gap_m: f64,
    /// Wing LE to canard TRAILING edge (nearest point) [m] — the declared
    /// arm convention (canard-mechanics arm is wing-datum).
    pub arm_wing_le_to_te_m: f64,
    /// Hinge axis as x/c from the canard leading edge (E1.5 wide prior
    /// [0.25, 0.50]; 'balanced too near the center').
    pub hinge_axis_xc: f64,
    /// Surface travel [rad] (photo-inferred ±30°; absent-by-verification).
    pub travel_rad: f64,
}

/// Rudder geometry (1903 twin movable surfaces).
#[derive(Clone, Debug, PartialEq)]
pub struct RudderGeometry {
    /// Total projected area [m²].
    pub area_m2: f64,
    /// Wing LE to rudder area centroid [m] (negative = aft).
    pub arm_m: f64,
}

/// The FlyerDesign schema (Round-2 complete field set for E3.1 scope).
#[derive(Clone, Debug, PartialEq)]
pub struct FlyerDesign {
    /// Design name (identity ingredient).
    pub name: &'static str,
    /// Wing geometry.
    pub wing: WingGeometry,
    /// Canard geometry + mechanics.
    pub canard: CanardGeometry,
    /// Rudder geometry.
    pub rudder: RudderGeometry,
    /// Lateral-control wiring.
    pub lateral: LateralControlTopology,
    /// Warp structure mode.
    pub warp: WarpStructureMode,
    /// Declared empty mass [kg] (must match the component sum).
    pub empty_mass_kg: f64,
    /// Pilot mass [kg].
    pub pilot_mass_kg: f64,
    /// Pilot position [m], frd from wing LE.
    pub pilot_position_m: [f64; 3],
    /// Empty-aircraft components (six-bucket engine convention included).
    pub components: Vec<MassComponent>,
}

/// Mass/CG/inertia build-up result (gross = empty + pilot).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassBuildUp {
    /// Gross mass [kg].
    pub gross_kg: f64,
    /// Gross CG position [m], frd from wing LE.
    pub cg_m: [f64; 3],
    /// Body-axis inertia about the GROSS CG [kg·m²] (diagonal; the 1903
    /// layout is laterally symmetric so products are ~0 by construction).
    pub inertia_kgm2: [f64; 3],
}

/// Derived-quantity panel (documented simple models, hand-calc pinned).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DerivedPanel {
    /// Canard volume coefficient V_c = S_c·l_c / (S_w·c_w), l_c = canard AC
    /// to gross CG.
    pub canard_volume: f64,
    /// Naive two-surface neutral point, x/c aft of the wing LE (equal lift
    /// slopes, no interference — a DOCUMENTED simple model).
    pub neutral_point_xc_naive: f64,
    /// Fixed-control static margin [fraction of chord]; negative = unstable.
    pub static_margin_fixed: f64,
    /// Free-control static margin under the floating-canard model with the
    /// given hinge ratio; more negative than fixed under overbalance.
    pub static_margin_free: f64,
    /// Thin-plate hinge-moment gradient dCh/dα ∝ 2π·(x_h − 1/4); positive =
    /// self-driving (the Orville overbalance sign).
    pub hinge_moment_gradient: f64,
}

impl FlyerDesign {
    /// Validate the design against the admitted domains.
    ///
    /// # Errors
    /// Typed refusals: `non-finite-input`, `span-outside-domain`,
    /// `chord-outside-domain`, `camber-outside-domain`,
    /// `area-inconsistent`, `hinge-axis-outside-prior`,
    /// `pilot-mass-outside-domain`, `component-count-exceeded`,
    /// `component-mass-invalid`, `mass-spec-mismatch`.
    #[allow(clippy::too_many_lines)] // one linear admission walk; splitting hides the order
    pub fn admit(&self) -> Result<(), Refusal> {
        let finite = [
            self.wing.span_m,
            self.wing.chord_m,
            self.wing.area_both_m2,
            self.wing.gap_m,
            self.wing.camber_ratio,
            self.wing.anhedral_droop_m,
            self.canard.span_m,
            self.canard.chord_m,
            self.canard.area_both_m2,
            self.canard.gap_m,
            self.canard.arm_wing_le_to_te_m,
            self.canard.hinge_axis_xc,
            self.canard.travel_rad,
            self.rudder.area_m2,
            self.rudder.arm_m,
            self.empty_mass_kg,
            self.pilot_mass_kg,
        ]
        .iter()
        .all(|v| v.is_finite());
        if !finite {
            return Err(refuse(
                "non-finite-input",
                "every design field must be finite".into(),
                "check unit conversions upstream",
            ));
        }
        if self.wing.span_m <= 0.0 || self.wing.span_m > MAX_SPAN_M {
            return Err(refuse(
                "span-outside-domain",
                format!("span {} m outside (0, {MAX_SPAN_M}]", self.wing.span_m),
                "1903 span is 12.29 m; check the metres conversion",
            ));
        }
        if self.wing.chord_m <= 0.0 || self.wing.chord_m > MAX_CHORD_M {
            return Err(refuse(
                "chord-outside-domain",
                format!("chord {} m outside (0, {MAX_CHORD_M}]", self.wing.chord_m),
                "1903 chord is 1.981 m",
            ));
        }
        if !(0.0..=MAX_CAMBER_RATIO).contains(&self.wing.camber_ratio) {
            return Err(refuse(
                "camber-outside-domain",
                format!(
                    "camber {} outside [0, {MAX_CAMBER_RATIO}]",
                    self.wing.camber_ratio
                ),
                "camber is a fraction of chord (1903: 0.05)",
            ));
        }
        let rect = 2.0 * self.wing.span_m * self.wing.chord_m;
        let dev = (self.wing.area_both_m2 - rect).abs() / self.wing.area_both_m2;
        if dev > AREA_CONSISTENCY_TOL {
            return Err(refuse(
                "area-inconsistent",
                format!(
                    "S_both {} m² vs 2·b·c {} m² deviates {:.1}% (> {:.0}%)",
                    self.wing.area_both_m2,
                    rect,
                    dev * 100.0,
                    AREA_CONSISTENCY_TOL * 100.0
                ),
                "check span/chord/area against geometry-conventions-v1 (projected planform)",
            ));
        }
        if !(0.25..=0.50).contains(&self.canard.hinge_axis_xc) {
            return Err(refuse(
                "hinge-axis-outside-prior",
                format!(
                    "hinge axis {} x/c outside the E1.5 prior [0.25, 0.50]",
                    self.canard.hinge_axis_xc
                ),
                "the hinge location is not published; stay inside the declared prior",
            ));
        }
        if !(MIN_PILOT_KG..=MAX_PILOT_KG).contains(&self.pilot_mass_kg) {
            return Err(refuse(
                "pilot-mass-outside-domain",
                format!(
                    "pilot {} kg outside [{MIN_PILOT_KG}, {MAX_PILOT_KG}]",
                    self.pilot_mass_kg
                ),
                "1903 pilots were ~65.8 kg (145 lb)",
            ));
        }
        if self.components.len() > MAX_COMPONENTS {
            return Err(refuse(
                "component-count-exceeded",
                format!(
                    "{} components exceed the cap {MAX_COMPONENTS}",
                    self.components.len()
                ),
                "aggregate minor items into buckets",
            ));
        }
        let mut sum = 0.0;
        for c in &self.components {
            let mass_ok = c.mass_kg.is_finite() && c.mass_kg > 0.0;
            let pos_ok = c.position_m.iter().all(|v| v.is_finite());
            let gyr_ok = c.gyration_sq_m2.iter().all(|v| v.is_finite() && *v >= 0.0);
            if !mass_ok || !pos_ok || !gyr_ok {
                return Err(refuse(
                    "component-mass-invalid",
                    format!(
                        "component {} has non-physical mass/position/gyration",
                        c.name
                    ),
                    "masses positive, positions finite, gyrations non-negative",
                ));
            }
            sum += c.mass_kg;
        }
        if (sum - self.empty_mass_kg).abs() > MASS_SPEC_TOL_KG {
            return Err(refuse(
                "mass-spec-mismatch",
                format!(
                    "components sum to {sum:.2} kg but the declared empty mass is {:.2} kg (tol {MASS_SPEC_TOL_KG} kg)",
                    self.empty_mass_kg
                ),
                "the declared empty mass must be the component sum — no hidden mass",
            ));
        }
        Ok(())
    }

    /// Component + pilot mass/CG/inertia build-up about the gross CG.
    ///
    /// # Errors
    /// Admission refusals (build-up requires an admitted design).
    pub fn mass_build_up(&self) -> Result<MassBuildUp, Refusal> {
        self.admit()?;
        let mut m_total = self.pilot_mass_kg;
        let mut first_moment = [
            self.pilot_mass_kg * self.pilot_position_m[0],
            self.pilot_mass_kg * self.pilot_position_m[1],
            self.pilot_mass_kg * self.pilot_position_m[2],
        ];
        for c in &self.components {
            m_total += c.mass_kg;
            for (fm, pos) in first_moment.iter_mut().zip(&c.position_m) {
                *fm += c.mass_kg * pos;
            }
        }
        let cg = [
            first_moment[0] / m_total,
            first_moment[1] / m_total,
            first_moment[2] / m_total,
        ];
        let mut inertia = [0.0f64; 3];
        let mut add = |mass: f64, pos: &[f64; 3], gyr: &[f64; 3]| {
            let d = [pos[0] - cg[0], pos[1] - cg[1], pos[2] - cg[2]];
            inertia[0] += mass * (d[1] * d[1] + d[2] * d[2] + gyr[0]);
            inertia[1] += mass * (d[0] * d[0] + d[2] * d[2] + gyr[1]);
            inertia[2] += mass * (d[0] * d[0] + d[1] * d[1] + gyr[2]);
        };
        add(
            self.pilot_mass_kg,
            &self.pilot_position_m,
            &[0.06, 0.35, 0.35],
        );
        for c in &self.components {
            add(c.mass_kg, &c.position_m, &c.gyration_sq_m2);
        }
        Ok(MassBuildUp {
            gross_kg: m_total,
            cg_m: cg,
            inertia_kgm2: inertia,
        })
    }

    /// The derived-quantity panel. `hinge_ratio` is ch_α/ch_δ from the E1.5
    /// prior machinery (positive under overbalance); the floating-canard
    /// model scales the free-control canard effectiveness by (1 + ratio).
    ///
    /// # Errors
    /// Admission refusals; `hinge-ratio-invalid` outside [0, 1].
    pub fn derived_panel(&self, hinge_ratio: f64) -> Result<DerivedPanel, Refusal> {
        if !hinge_ratio.is_finite() || !(0.0..=1.0).contains(&hinge_ratio) {
            return Err(refuse(
                "hinge-ratio-invalid",
                format!("hinge ratio {hinge_ratio} outside [0, 1]"),
                "derive the ratio from the E1.5 hinge-axis prior",
            ));
        }
        let build = self.mass_build_up()?;
        let (s_w, c_w) = (self.wing.area_both_m2, self.wing.chord_m);
        let s_c = self.canard.area_both_m2;
        // Canard AC: arm measures wing LE → canard TE; AC at 3/4 canard
        // chord further forward (+x).
        let x_ac_c = self.canard.arm_wing_le_to_te_m + 0.75 * self.canard.chord_m;
        let x_ac_w = -0.25 * c_w; // wing AC quarter-chord aft of the LE
        let x_cg = build.cg_m[0];
        let l_c = x_ac_c - x_cg;
        let canard_volume = (s_c * l_c) / (s_w * c_w);
        // Naive two-surface NP (equal slopes, no interference), then the
        // free-control variant with amplified canard effectiveness.
        let np = |s_c_eff: f64| -> f64 {
            let x_np = (s_w * x_ac_w + s_c_eff * x_ac_c) / (s_w + s_c_eff);
            -x_np / c_w // report as x/c AFT of the wing LE (positive aft)
        };
        let np_fixed = np(s_c);
        let np_free = np(s_c * (1.0 + hinge_ratio));
        let cg_xc = -x_cg / c_w;
        // Margin = NP aft of CG (positive = stable).
        let static_margin_fixed = np_fixed - cg_xc;
        let static_margin_free = np_free - cg_xc;
        let hinge_moment_gradient =
            2.0 * core::f64::consts::PI * (self.canard.hinge_axis_xc - 0.25);
        Ok(DerivedPanel {
            canard_volume,
            neutral_point_xc_naive: np_fixed,
            static_margin_fixed,
            static_margin_free,
            hinge_moment_gradient,
        })
    }

    /// Content digest of the design (canonical little-endian field bytes
    /// under [`DESIGN_DIGEST_DOMAIN`]). Identity ingredient for
    /// PhysicalScenarioId (replay-identity-schema-v1).
    #[must_use]
    pub fn digest(&self) -> String {
        let mut payload = Vec::new();
        payload.extend_from_slice(self.name.as_bytes());
        let mut push = |v: f64| payload.extend_from_slice(&v.to_bits().to_le_bytes());
        for v in [
            self.wing.span_m,
            self.wing.chord_m,
            self.wing.area_both_m2,
            self.wing.gap_m,
            self.wing.camber_ratio,
            self.wing.anhedral_droop_m,
            self.canard.span_m,
            self.canard.chord_m,
            self.canard.area_both_m2,
            self.canard.gap_m,
            self.canard.arm_wing_le_to_te_m,
            self.canard.hinge_axis_xc,
            self.canard.travel_rad,
            self.rudder.area_m2,
            self.rudder.arm_m,
            self.empty_mass_kg,
            self.pilot_mass_kg,
            self.pilot_position_m[0],
            self.pilot_position_m[1],
            self.pilot_position_m[2],
        ] {
            push(v);
        }
        match self.lateral {
            LateralControlTopology::WarpWithSlavedRudder { ratio } => {
                payload.push(1);
                payload.extend_from_slice(&ratio.to_bits().to_le_bytes());
            }
            LateralControlTopology::WarpIndependentRudder => payload.push(2),
        }
        payload.push(match self.warp {
            WarpStructureMode::FlexibleTrussed => 1,
            WarpStructureMode::Rigidized => 2,
        });
        for c in &self.components {
            payload.extend_from_slice(c.name.as_bytes());
            for v in [c.mass_kg, c.position_m[0], c.position_m[1], c.position_m[2]] {
                payload.extend_from_slice(&v.to_bits().to_le_bytes());
            }
            for v in c.gyration_sq_m2 {
                payload.extend_from_slice(&v.to_bits().to_le_bytes());
            }
        }
        hash_domain(DESIGN_DIGEST_DOMAIN, &payload).to_hex()
    }

    /// The reference 1903 configuration, built from the frozen dossiers
    /// (flyer-reference + canard-mechanics + geometry conventions).
    ///
    /// PROVENANCE: geometry/mass values are the verified dossier values;
    /// component POSITIONS and GYRATION RADII are a CALIBRATED
    /// reconstruction — positions place the gross CG at the dossier's
    /// 29.7% chord, radii are tuned so the build-up reproduces the
    /// published single-lineage inertias within the documented band. They
    /// are inputs to a cross-check, not measurements.
    #[must_use]
    pub fn reference_1903() -> FlyerDesign {
        FlyerDesign {
            name: "wright-flyer-1903-reference",
            wing: WingGeometry {
                span_m: 12.29,
                chord_m: 1.981,
                area_both_m2: 47.38,
                gap_m: 1.89,
                camber_ratio: 0.05,
                anhedral_droop_m: 0.254,
            },
            canard: CanardGeometry {
                span_m: 3.658,
                chord_m: 0.762,
                area_both_m2: 4.46,
                gap_m: 0.637,
                arm_wing_le_to_te_m: 2.231,
                hinge_axis_xc: 0.375,
                travel_rad: core::f64::consts::FRAC_PI_6,
            },
            rudder: RudderGeometry {
                area_m2: 1.86,
                arm_m: -3.35,
            },
            lateral: LateralControlTopology::WarpWithSlavedRudder { ratio: 2.5 },
            warp: WarpStructureMode::FlexibleTrussed,
            empty_mass_kg: 274.4,
            pilot_mass_kg: 65.77,
            pilot_position_m: [-0.588, -0.30, 0.0],
            components: vec![
                MassComponent {
                    name: "wings-struts-wires",
                    mass_kg: 122.0,
                    position_m: [-0.696, 0.0, 0.0],
                    gyration_sq_m2: [12.588, 0.60, 12.915],
                },
                MassComponent {
                    name: "skids-frame",
                    mass_kg: 20.0,
                    position_m: [-0.50, 0.0, 0.55],
                    gyration_sq_m2: [1.0, 2.2, 1.75],
                },
                MassComponent {
                    name: "canard-structure",
                    mass_kg: 13.0,
                    position_m: [2.30, 0.0, 0.30],
                    gyration_sq_m2: [1.116, 1.2, 1.164],
                },
                MassComponent {
                    name: "rudder-structure",
                    mass_kg: 7.0,
                    position_m: [-3.35, 0.0, 0.0],
                    gyration_sq_m2: [0.30, 1.0, 0.50],
                },
                MassComponent {
                    name: "engine-dry",
                    mass_kg: 68.0,
                    position_m: [-0.60, 0.30, 0.0],
                    gyration_sq_m2: [0.04, 0.04, 0.04],
                },
                MassComponent {
                    name: "engine-installed-accessories",
                    mass_kg: 8.0,
                    position_m: [-0.60, 0.30, 0.0],
                    gyration_sq_m2: [0.02, 0.02, 0.02],
                },
                MassComponent {
                    name: "cooling-water",
                    mass_kg: 11.0,
                    position_m: [-0.35, 0.25, -0.30],
                    gyration_sq_m2: [0.04, 0.04, 0.04],
                },
                MassComponent {
                    name: "ignition",
                    mass_kg: 3.0,
                    position_m: [-0.50, 0.30, 0.0],
                    gyration_sq_m2: [0.01, 0.01, 0.01],
                },
                MassComponent {
                    name: "drivetrain-chains",
                    mass_kg: 16.0,
                    position_m: [-0.85, 0.0, -0.10],
                    gyration_sq_m2: [0.50, 0.80, 0.50],
                },
                MassComponent {
                    name: "propellers",
                    mass_kg: 6.4,
                    position_m: [-1.30, 0.0, -0.10],
                    gyration_sq_m2: [1.103, 0.05, 1.103],
                },
            ],
        }
    }
}
