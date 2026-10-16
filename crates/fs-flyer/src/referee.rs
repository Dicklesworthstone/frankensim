//! Early steady referee fixtures + the E10.1 receipt format (bead
//! wf-root-guzez.5.21, E4.9a). Plan §10: the referee lane surfaces
//! discrepancy EARLY — independent closed-form references checked
//! against the production solvers through the reference plane, with
//! every difference RECORDED as a receipt (differences reported, never
//! forced to vanish — the V-05 philosophy). This module defines the
//! harness receipt format E10.1 builds on and the first steady cases.
//!
//! INDEPENDENCE RULE: referee values in this module come from classical
//! closed forms implemented HERE (lifting-line, actuator-disk momentum)
//! sharing no computation path with fs-wing/fs-airscrew — only frozen
//! geometry constants cross the boundary, and each receipt carries an
//! independence note saying so.

use crate::Refusal;
use fs_blake3::hash_domain;

/// Receipt schema id (the E10.1 harness format, first minted here).
pub const RECEIPT_SCHEMA: &str = "org.frankensim.wf.referee-receipt.v1";

/// Relative-discrepancy magnitude cap (absurd-input guard: a referee
/// disagreeing by more than 100x means a broken fixture, not physics).
pub const MAX_REL_DISCREPANCY: f64 = 100.0;

/// One discrepancy receipt (the harness row).
#[derive(Clone, Debug, PartialEq)]
pub struct DiscrepancyReceipt {
    /// Stable case id.
    pub case_id: &'static str,
    /// The compared quantity.
    pub quantity: &'static str,
    /// Units.
    pub units: &'static str,
    /// Production-solver value.
    pub production: f64,
    /// Independent referee value.
    pub referee: f64,
    /// Signed relative discrepancy (production/referee − 1).
    pub rel_discrepancy: f64,
    /// Why the referee is independent.
    pub independence_note: &'static str,
    /// The declared comparison class (a BAND the pair is expected to
    /// share, or "reported-only" when no band is licensed).
    pub comparison_class: &'static str,
}

impl DiscrepancyReceipt {
    /// Build + validate a receipt.
    ///
    /// # Errors
    /// `referee-receipt-invalid` (non-finite values, |rel| beyond
    /// [`MAX_REL_DISCREPANCY`] — cap and cap+1, referee exactly zero,
    /// empty ids).
    pub fn new(
        case_id: &'static str,
        quantity: &'static str,
        units: &'static str,
        production: f64,
        referee: f64,
        independence_note: &'static str,
        comparison_class: &'static str,
    ) -> Result<Self, Refusal> {
        if case_id.is_empty()
            || quantity.is_empty()
            || !production.is_finite()
            || !referee.is_finite()
            || referee == 0.0
        {
            return Err(Refusal {
                code: "referee-receipt-invalid",
                message: format!(
                    "{case_id}/{quantity}: production {production:?}, referee {referee:?}"
                ),
                ranked_repairs: vec!["finite values; nonzero referee; named case".into()],
            });
        }
        let rel = production / referee - 1.0;
        if !rel.is_finite() || rel.abs() > MAX_REL_DISCREPANCY {
            return Err(Refusal {
                code: "referee-receipt-invalid",
                message: format!("{case_id}: relative discrepancy {rel:?} beyond the sanity cap"),
                ranked_repairs: vec!["a 100x disagreement is a broken fixture, not physics".into()],
            });
        }
        Ok(DiscrepancyReceipt {
            case_id,
            quantity,
            units,
            production,
            referee,
            rel_discrepancy: rel,
            independence_note,
            comparison_class,
        })
    }

    /// Render the E10.1 harness JSONL row (stable field order).
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        format!(
            "{{\"schema\":\"{RECEIPT_SCHEMA}\",\"case_id\":\"{}\",\"quantity\":\"{}\",\
             \"units\":\"{}\",\"production\":{},\"referee\":{},\"rel_discrepancy\":{},\
             \"independence\":\"{}\",\"comparison_class\":\"{}\"}}",
            self.case_id,
            self.quantity,
            self.units,
            self.production,
            self.referee,
            self.rel_discrepancy,
            self.independence_note,
            self.comparison_class
        )
    }

    /// Content digest over the canonical payload.
    #[must_use]
    pub fn digest(&self) -> String {
        let mut p = Vec::new();
        p.extend_from_slice(self.case_id.as_bytes());
        p.push(0);
        p.extend_from_slice(self.quantity.as_bytes());
        p.push(0);
        p.extend_from_slice(&self.production.to_bits().to_le_bytes());
        p.extend_from_slice(&self.referee.to_bits().to_le_bytes());
        hash_domain("org.frankensim.wf.referee-receipt.v1", &p).to_hex()
    }
}

/// Classical monoplane lifting-line CL (the INDEPENDENT reference):
/// CL = a₀·(α + 2h) / (1 + a₀/(π·AR·e)). Closed form, no panels.
#[must_use]
pub fn referee_liftingline_cl(alpha_rad: f64, camber_ratio: f64, ar: f64, oswald_e: f64) -> f64 {
    let a0 = core::f64::consts::TAU;
    a0 * (alpha_rad + 2.0 * camber_ratio) / (1.0 + a0 / (core::f64::consts::PI * ar * oswald_e))
}

/// Classical actuator-disk momentum thrust for a prop absorbing power P
/// at axial speed V (the INDEPENDENT propulsion reference): solves
/// P = T·(V + w), T = 2ρA(V + w)w by deterministic bisection on w.
///
/// # Errors
/// `referee-momentum-invalid` (non-physical inputs).
pub fn referee_momentum_thrust(
    power_w: f64,
    v_mps: f64,
    rho_kg_m3: f64,
    disk_area_m2: f64,
) -> Result<f64, Refusal> {
    if !(power_w > 0.0 && v_mps >= 0.0 && rho_kg_m3 > 0.0 && disk_area_m2 > 0.0)
        || !power_w.is_finite()
        || !v_mps.is_finite()
    {
        return Err(Refusal {
            code: "referee-momentum-invalid",
            message: format!("P {power_w:?}, V {v_mps:?}, rho {rho_kg_m3:?}, A {disk_area_m2:?}"),
            ranked_repairs: vec!["positive power/density/area; V >= 0".into()],
        });
    }
    // P(w) = 2*rho*A*(V+w)^2*w is monotone in w >= 0: bisect.
    let f = |w: f64| 2.0 * rho_kg_m3 * disk_area_m2 * (v_mps + w) * (v_mps + w) * w - power_w;
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    while f(hi) < 0.0 {
        hi *= 2.0;
        if hi > 1.0e6 {
            return Err(Refusal {
                code: "referee-momentum-invalid",
                message: "no momentum solution below 1e6 m/s".into(),
                ranked_repairs: vec!["check the power/area scales".into()],
            });
        }
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let w = 0.5 * (lo + hi);
    Ok(2.0 * rho_kg_m3 * disk_area_m2 * (v_mps + w) * w)
}

/// Prandtl's classical biplane correction: the effective span-efficiency
/// multiplier for two equal wings at gap/span ratio g/b (interpolation
/// of Munk's sigma; the INDEPENDENT biplane reference).
#[must_use]
pub fn referee_biplane_efficiency(gap_over_span: f64) -> f64 {
    // Munk interference factor sigma ~ (1 - 0.66*g/b)/(1.05 + 3.7*g/b)
    // (classical fit, valid 0.05 <= g/b <= 0.5).
    let sigma = (1.0 - 0.66 * gap_over_span) / (1.05 + 3.7 * gap_over_span);
    // Two equal wings: e_biplane = 1/(1+sigma) relative to the ideal.
    1.0 / (1.0 + sigma)
}
