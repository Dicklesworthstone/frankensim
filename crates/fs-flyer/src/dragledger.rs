//! Parasite-drag component ledger + required-power ledger (bead
//! wf-root-guzez.5.13.1, E4.6a-i; V-13a/V-13b). Plan §5.1.2: parasite
//! drag is a COMPONENT LEDGER — pilot, engine/radiators, skids, wires,
//! uprights, struts, chains, misc — each with area, drag coefficient,
//! coefficient source, and uncertainty; the ledger total carries an
//! EXPLICIT `unresolved_interference_drag` line (never hidden inside a
//! padded component), and the flat-plate aggregate remains a SEPARATELY
//! IDENTIFIED fallback mode that can never masquerade as the ledger.
//!
//! Ownership boundary (AeroEffectOwners discipline): this ledger owns
//! NON-LIFTING parasite drag only. Wing/canard profile drag belongs to
//! the section closures and induced drag to the 3-D induction owner —
//! summing those here would double-count.
//!
//! At the Flyer's power margin (~12 hp sustained), parasite-drag error
//! decides whether the aircraft flies at all — V-13a checks the assembled
//! power requirement against the independent power-balance band.

use crate::Refusal;
use fs_blake3::hash_domain;

/// Component-count cap (absurd-input guard).
pub const MAX_COMPONENTS: usize = 32;

/// One ledger line: a bluff/streamlined component with cited drag data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DragComponent {
    /// Stable component id.
    pub component_id: &'static str,
    /// Reference (frontal) area [m²].
    pub area_m2: f64,
    /// Drag coefficient on that area.
    pub cd: f64,
    /// Coefficient provenance (cited estimate class).
    pub cd_source: &'static str,
    /// One-sigma fractional uncertainty on the component's drag.
    pub uncertainty_frac: f64,
}

/// The ledger: components + the explicit unresolved-interference line.
#[derive(Clone, Debug, PartialEq)]
pub struct DragLedger {
    /// Component lines.
    pub components: Vec<DragComponent>,
    /// Unresolved interference as a fraction of the component sum —
    /// mutual-wake/junction drag no component can honestly claim.
    pub unresolved_interference_frac: f64,
}

/// One evaluated line (per-item oracle surface).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LedgerItem {
    /// Component id.
    pub component_id: &'static str,
    /// Drag at the evaluated state [N].
    pub drag_n: f64,
    /// Power absorbed at the evaluated state [W].
    pub power_w: f64,
    /// One-sigma drag uncertainty [N].
    pub sigma_n: f64,
}

/// The evaluated ledger.
#[derive(Clone, Debug, PartialEq)]
pub struct LedgerReport {
    /// Per-component lines, ledger order.
    pub items: Vec<LedgerItem>,
    /// Sum of component drags [N].
    pub component_sum_n: f64,
    /// The explicit unresolved-interference line [N].
    pub unresolved_interference_drag_n: f64,
    /// Total parasite drag [N] (= component sum + interference).
    pub total_parasite_n: f64,
    /// RSS one-sigma on the total [N] (components independent; the
    /// interference line carries its own full-magnitude sigma).
    pub sigma_total_n: f64,
    /// Total parasite power [W] (= total drag × V).
    pub power_w: f64,
    /// Ledger content digest (identity of WHAT was summed).
    pub ledger_digest: String,
}

/// The separately identified flat-plate fallback result — deliberately a
/// DIFFERENT type so it cannot be passed where a ledger report is due.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlatPlateAggregate {
    /// Equivalent flat-plate area [m²].
    pub flat_plate_area_m2: f64,
    /// Drag [N].
    pub drag_n: f64,
    /// Power [W].
    pub power_w: f64,
}

/// The registered Wright v1 ledger. Areas and coefficients are Estimated
/// reconstruction values (bluff-body handbook classes); the wide sigmas
/// say so. Flat-plate-equivalent sum ≈ 1.45 m² + 10% interference.
#[must_use]
pub fn wright_ledger_v1() -> DragLedger {
    DragLedger {
        components: vec![
            DragComponent {
                component_id: "pilot-prone",
                area_m2: 0.28,
                cd: 1.0,
                cd_source: "prone-man bluff estimate (Culick/Jex discussion class)",
                uncertainty_frac: 0.30,
            },
            DragComponent {
                component_id: "engine-radiators",
                area_m2: 0.18,
                cd: 1.2,
                cd_source: "bluff block + radiator handbook class",
                uncertainty_frac: 0.35,
            },
            DragComponent {
                component_id: "skids-frame",
                area_m2: 0.15,
                cd: 1.1,
                cd_source: "long slender bluff frame estimate",
                uncertainty_frac: 0.40,
            },
            DragComponent {
                component_id: "wires-rigging",
                area_m2: 0.27,
                cd: 1.1,
                cd_source: "subcritical circular cylinder, summed rigging length",
                uncertainty_frac: 0.40,
            },
            DragComponent {
                component_id: "uprights-struts",
                area_m2: 0.83,
                cd: 0.45,
                cd_source: "rounded-rectangle strut class, partially streamlined",
                uncertainty_frac: 0.45,
            },
            DragComponent {
                component_id: "chains-sprockets",
                area_m2: 0.05,
                cd: 1.0,
                cd_source: "bluff estimate",
                uncertainty_frac: 0.50,
            },
            DragComponent {
                component_id: "misc-fittings",
                area_m2: 0.08,
                cd: 1.0,
                cd_source: "bluff allowance",
                uncertainty_frac: 0.50,
            },
        ],
        unresolved_interference_frac: 0.10,
    }
}

impl DragLedger {
    /// Admit the ledger.
    ///
    /// # Errors
    /// `ledger-component-invalid` (empty id, non-finite/non-positive
    /// area or cd, uncertainty outside [0, 1]),
    /// `ledger-component-count-invalid` (0 or > [`MAX_COMPONENTS`]),
    /// `ledger-component-duplicate`,
    /// `ledger-interference-invalid` (outside [0, 0.5]).
    pub fn admit(&self) -> Result<(), Refusal> {
        if self.components.is_empty() || self.components.len() > MAX_COMPONENTS {
            return Err(Refusal {
                code: "ledger-component-count-invalid",
                message: format!(
                    "{} components (admitted 1..={MAX_COMPONENTS})",
                    self.components.len()
                ),
                ranked_repairs: vec!["merge sub-items into a named component".into()],
            });
        }
        for (i, c) in self.components.iter().enumerate() {
            let ok = !c.component_id.is_empty()
                && c.area_m2.is_finite()
                && c.area_m2 > 0.0
                && c.cd.is_finite()
                && c.cd > 0.0
                && c.uncertainty_frac.is_finite()
                && (0.0..=1.0).contains(&c.uncertainty_frac);
            if !ok {
                return Err(Refusal {
                    code: "ledger-component-invalid",
                    message: format!(
                        "component {i} ({}): area {:?} m², cd {:?}, sigma {:?}",
                        c.component_id, c.area_m2, c.cd, c.uncertainty_frac
                    ),
                    ranked_repairs: vec![
                        "positive finite area and cd; uncertainty_frac in [0, 1]".into(),
                    ],
                });
            }
            if self.components[..i]
                .iter()
                .any(|p| p.component_id == c.component_id)
            {
                return Err(Refusal {
                    code: "ledger-component-duplicate",
                    message: format!("duplicate component id {}", c.component_id),
                    ranked_repairs: vec!["one line per physical component".into()],
                });
            }
        }
        if !self.unresolved_interference_frac.is_finite()
            || !(0.0..=0.5).contains(&self.unresolved_interference_frac)
        {
            return Err(Refusal {
                code: "ledger-interference-invalid",
                message: format!(
                    "unresolved_interference_frac {:?} outside [0, 0.5]",
                    self.unresolved_interference_frac
                ),
                ranked_repairs: vec![
                    "interference beyond 50% of the sum means the ledger itself is wrong".into(),
                ],
            });
        }
        Ok(())
    }

    /// Ledger content digest (ModelId ingredient).
    #[must_use]
    pub fn digest(&self) -> String {
        let mut p = Vec::new();
        for c in &self.components {
            p.extend_from_slice(
                &u32::try_from(c.component_id.len())
                    .expect("short id")
                    .to_le_bytes(),
            );
            p.extend_from_slice(c.component_id.as_bytes());
            for v in [c.area_m2, c.cd, c.uncertainty_frac] {
                p.extend_from_slice(&v.to_bits().to_le_bytes());
            }
        }
        p.extend_from_slice(&self.unresolved_interference_frac.to_bits().to_le_bytes());
        hash_domain("org.frankensim.fs-flyer.drag-ledger.v1", &p).to_hex()
    }

    /// Evaluate the ledger at an air state.
    ///
    /// # Errors
    /// Admission refusals pass through; `air-state-invalid` (rho or V
    /// non-finite or non-positive).
    pub fn evaluate(&self, rho_kg_m3: f64, v_mps: f64) -> Result<LedgerReport, Refusal> {
        self.admit()?;
        if !(rho_kg_m3.is_finite() && rho_kg_m3 > 0.0 && v_mps.is_finite() && v_mps > 0.0) {
            return Err(Refusal {
                code: "air-state-invalid",
                message: format!("rho {rho_kg_m3:?}, V {v_mps:?}"),
                ranked_repairs: vec!["evaluate at a physical air state".into()],
            });
        }
        let q = 0.5 * rho_kg_m3 * v_mps * v_mps;
        let mut items = Vec::with_capacity(self.components.len());
        let mut sum = 0.0;
        let mut var = 0.0;
        for c in &self.components {
            let d = q * c.area_m2 * c.cd;
            let sigma = d * c.uncertainty_frac;
            items.push(LedgerItem {
                component_id: c.component_id,
                drag_n: d,
                power_w: d * v_mps,
                sigma_n: sigma,
            });
            sum += d;
            var += sigma * sigma;
        }
        let interference = sum * self.unresolved_interference_frac;
        // The interference line is uncertain at its own full magnitude.
        var += interference * interference;
        let total = sum + interference;
        Ok(LedgerReport {
            items,
            component_sum_n: sum,
            unresolved_interference_drag_n: interference,
            total_parasite_n: total,
            sigma_total_n: var.sqrt(),
            power_w: total * v_mps,
            ledger_digest: self.digest(),
        })
    }
}

/// The flat-plate aggregate FALLBACK mode (separately identified; plan
/// §5.1.2). Not a ledger: no per-component lines, no allocation claims.
///
/// # Errors
/// `flat-plate-area-invalid`, `air-state-invalid`.
pub fn flat_plate_aggregate(
    flat_plate_area_m2: f64,
    rho_kg_m3: f64,
    v_mps: f64,
) -> Result<FlatPlateAggregate, Refusal> {
    if !(flat_plate_area_m2.is_finite() && flat_plate_area_m2 > 0.0) {
        return Err(Refusal {
            code: "flat-plate-area-invalid",
            message: format!("f = {flat_plate_area_m2:?} m²"),
            ranked_repairs: vec!["positive finite equivalent area".into()],
        });
    }
    if !(rho_kg_m3.is_finite() && rho_kg_m3 > 0.0 && v_mps.is_finite() && v_mps > 0.0) {
        return Err(Refusal {
            code: "air-state-invalid",
            message: format!("rho {rho_kg_m3:?}, V {v_mps:?}"),
            ranked_repairs: vec!["evaluate at a physical air state".into()],
        });
    }
    let d = 0.5 * rho_kg_m3 * v_mps * v_mps * flat_plate_area_m2;
    Ok(FlatPlateAggregate {
        flat_plate_area_m2,
        drag_n: d,
        power_w: d * v_mps,
    })
}
