//! Conjugate coupling: the partitioned solid-conduction ↔ air-path fixed
//! point.
//!
//! Before this module the E05 thermal vertical was conduction-only against a
//! *declared* ambient: a caller picked one bulk air temperature, handed it to
//! every Robin row, and the solid never fed back. That is a one-way chain, not
//! a conjugate solve. Two things are missing from a one-way chain and both are
//! built here.
//!
//! **Air carries state.** [`OperatingPoint`](crate::OperatingPoint) and
//! [`BranchFlow`](crate::BranchFlow) have no temperature field — every
//! `Temperature` elsewhere in this crate is a *solid* temperature or a
//! caller-declared limit. Air that absorbs heat gets hotter along the channel,
//! and in an enclosure that rise is not negligible: it is the difference
//! between the upstream and downstream ends of the same heatsink seeing the
//! same coolant. [`AirPath`] adds that state as a 1-D stream-wise network.
//!
//! **The exchange is iterated, not declared.** The solid's surface temperature
//! sets how much heat the air takes; the air's temperature sets how much heat
//! the solid loses. [`solve_conjugate`] closes that loop as a partitioned
//! fixed point over the per-region reference temperature vector.
//!
//! # The segment law
//!
//! Per segment, with wall temperature `T_w`, coefficient `h`, wetted area `A`
//! and capacity rate `ṁ c_p`:
//!
//! ```text
//! NTU        = h A / (ṁ c_p)
//! ε          = 1 − e^(−NTU)
//! Q          = ṁ c_p (T_w − T_in) ε
//! T_ref,eff  = T_w − (T_w − T_in) · ε/NTU
//! T_out      = T_in + (T_w − T_in) ε
//! ```
//!
//! `T_ref,eff` is chosen so that `h A (T_w − T_ref,eff) ≡ Q` *exactly*: the
//! reference temperature handed to the Robin row reproduces the heat the air
//! actually absorbed, so the two sides cannot disagree about the same
//! interface by construction of the reference alone.
//!
//! This is the exponential law, not an arithmetic inlet/outlet mean. The mean
//! is not a cosmetic difference:
//!
//! * The exponential form is unconditionally bounded in `[T_in, T_w]`. The
//!   arithmetic-mean reference leaves the physical interval past `NTU ≈ 2` and
//!   destabilises the outer fixed point exactly where enclosures operate.
//! * It is *exact* under segment refinement. Since
//!   `e^(−NTU) = (e^(−NTU/N))^N`, splitting one uniform-wall channel into `N`
//!   equal segments must give the same outlet temperature for every `N`. That
//!   is a free exact falsifier with no tolerance to tune, and the
//!   arithmetic-mean model fails it with clean `N`-dependence. See
//!   `tests/conjugate.rs`.
//!
//! # What the flux balance proves
//!
//! [`FluxBalanceAudit`] reports two different things and they carry different
//! authority.
//!
//! The **per-region** balance (`solid_heat_rate_w` vs `air_heat_rate_w`) is an
//! *algebraic identity* at the fixed point: `T_ref,eff` was defined to make it
//! one, and for a uniform `h` the area-weighted mean wall temperature
//! reproduces `∫ h (T_h − T_ref) dA` exactly. It therefore falsifies WIRING —
//! a dropped face, a region bound to the wrong segment, a wrong `ṁ` — and NOT
//! physics. Away from the fixed point its residual is the coupling residual
//! itself, which is why it is also the loop's convergence witness.
//!
//! The **decomposition** cross-check ([`FluxBalanceAudit::
//! decomposition_residual_w`]) compares `Σ_j Q_solid,j` against
//! [`fs_conduction::EnergyBalance::robin_out_w`], which `fs-conduction`
//! accumulates in its own loop over every Robin face rather than re-deriving
//! from these parts. That catches a face owned by no declared region or
//! counted twice — a decomposition bug the per-region identity structurally
//! cannot see. It is still a wiring falsifier, not a conservation proof.
//!
//! The load-bearing *physical* checks are the closed-form heated-channel
//! comparison and the segment-refinement invariance above.
//!
//! # No-claim boundaries
//!
//! * **`h` is frozen across the loop.** The coupling variable is the reference
//!   temperature vector only. Air properties are not re-evaluated at the drifting
//!   film or bulk temperature, so a temperature-dependent `h` is outside this
//!   model. A caller wanting that must re-run the driver with a new [`AirPath`].
//! * **Relaxation is scalar.** [`fs_couple::AitkenRelaxation`] is a scalar Δ²
//!   relaxer; this driver projects the vector residual onto one area-weighted
//!   scalar and applies a single `ω` to the whole vector. This is NOT a vector
//!   interface accelerator: interface quasi-Newton (IQN-ILS) does not exist in
//!   `fs-couple` and is not called here.
//! * **One-dimensional air.** The path is a stream-wise chain of well-mixed
//!   segments. There is no lateral mixing, recirculation, buoyancy, or flow
//!   redistribution driven by heating, and no momentum coupling back to the
//!   operating point.
//! * **The seam is typed, not ledgered.** [`SEAM_PORT_KIND`] declares the
//!   `fs-couple` port kind and its effort dimension is checked, and
//!   [`fs_couple::EnergyAudit`] records every exchange's imbalance. The
//!   window-balance ledger path (`BoundaryTemperatureReference`,
//!   `WindowEvidenceRef`) is not wired.

use core::cmp::Ordering;

use fs_couple::{AitkenRelaxation, EnergyAudit, PortKind};
use fs_exec::Cx;
use fs_math::det;
use fs_qty::Dims;

use crate::AirflowError;

/// The `fs-couple` port kind this seam exchanges: temperature (effort) against
/// entropy flow.
///
/// The driver's coupling variable is a temperature, so this is a checked
/// declaration rather than a label — see [`seam_effort_dimensions`].
pub const SEAM_PORT_KIND: PortKind = PortKind::ThermalTemperatureEntropy;

/// Canonical effort dimensions of [`SEAM_PORT_KIND`].
///
/// Equal to [`Temperature::DIMS`] by the port algebra, which is what makes the
/// exchanged reference temperature a legal effort for this port rather than an
/// untyped scalar.
#[must_use]
pub fn seam_effort_dimensions() -> Dims {
    SEAM_PORT_KIND.canonical_effort_dimensions()
}

fn finite_positive(field: &'static str, value: f64) -> Result<f64, AirflowError> {
    if !(value.is_finite() && value > 0.0) {
        return Err(AirflowError::InvalidConjugateInput {
            field,
            value_bits: value.to_bits(),
        });
    }
    Ok(value)
}

fn finite(field: &'static str, value: f64) -> Result<f64, AirflowError> {
    if !value.is_finite() {
        return Err(AirflowError::InvalidConjugateInput {
            field,
            value_bits: value.to_bits(),
        });
    }
    Ok(value)
}

/// One ordered segment of the air path: the named solid Robin region it
/// exchanges heat with, that region's wetted area, and the coefficient the
/// correlation rung produced for it.
#[derive(Debug, Clone, PartialEq)]
pub struct AirSegment {
    region: String,
    area_m2: f64,
    htc_w_per_m2_k: f64,
}

impl AirSegment {
    /// Declare a segment.
    ///
    /// # Errors
    /// [`AirflowError::EmptyElementName`] for a blank region;
    /// [`AirflowError::InvalidConjugateInput`] for a non-finite or
    /// non-positive area or coefficient.
    pub fn new(
        region: &str,
        area_m2: f64,
        htc_w_per_m2_k: f64,
    ) -> Result<AirSegment, AirflowError> {
        if region.trim().is_empty() {
            return Err(AirflowError::EmptyElementName);
        }
        Ok(AirSegment {
            region: region.to_string(),
            area_m2: finite_positive("air segment area", area_m2)?,
            htc_w_per_m2_k: finite_positive(
                "air segment heat-transfer coefficient",
                htc_w_per_m2_k,
            )?,
        })
    }

    /// The solid Robin region this segment exchanges with.
    #[must_use]
    pub fn region(&self) -> &str {
        &self.region
    }

    /// Wetted area, m².
    #[must_use]
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    /// Heat-transfer coefficient, W/(m²·K).
    #[must_use]
    pub const fn htc_w_per_m2_k(&self) -> f64 {
        self.htc_w_per_m2_k
    }

    /// `h A / (ṁ c_p)` for this segment against a capacity rate.
    #[must_use]
    pub fn ntu(&self, capacity_rate_w_per_k: f64) -> f64 {
        self.htc_w_per_m2_k * self.area_m2 / capacity_rate_w_per_k
    }
}

/// `(ε, ε/NTU)` for `ε = 1 − e^(−NTU)`.
///
/// `ε/NTU` is the reference-temperature weight and tends to 1 as `NTU → 0`,
/// where the effective reference degenerates to the segment inlet — the
/// physically correct limit for a segment that exchanges no heat. `NTU = 0` is
/// returned as that limit directly rather than evaluated as `0/0`.
fn effectiveness(ntu: f64) -> (f64, f64) {
    if ntu == 0.0 {
        return (0.0, 1.0);
    }
    let eps = -det::expm1(-ntu);
    (eps, eps / ntu)
}

/// One marched segment's state.
#[derive(Debug, Clone, PartialEq)]
pub struct SegmentState {
    /// The solid Robin region name.
    pub region: String,
    /// Air temperature entering this segment, K.
    pub inlet_temperature_k: f64,
    /// Air temperature leaving this segment, K.
    pub outlet_temperature_k: f64,
    /// The effective Robin reference temperature, K: the value that makes
    /// `h A (T_w − T_ref) ≡ Q` for this segment.
    pub reference_temperature_k: f64,
    /// `h A / (ṁ c_p)`.
    pub ntu: f64,
    /// `1 − e^(−NTU)`.
    pub effectiveness: f64,
    /// Heat absorbed by the air over this segment, W (positive = into the
    /// air, matching the solid's positive-leaving convention).
    pub heat_rate_w: f64,
}

/// The result of marching the whole path.
#[derive(Debug, Clone, PartialEq)]
pub struct AirMarch {
    /// Per-segment states in path order.
    pub segments: Vec<SegmentState>,
    /// Air temperature leaving the last segment, K.
    pub outlet_temperature_k: f64,
    /// `Σ_j Q_j`, W.
    pub total_heat_rate_w: f64,
}

impl AirMarch {
    /// The effective reference temperatures in path order.
    #[must_use]
    pub fn reference_temperatures_k(&self) -> Vec<f64> {
        self.segments
            .iter()
            .map(|s| s.reference_temperature_k)
            .collect()
    }
}

/// The stream-wise 1-D air path over one solved flow branch.
#[derive(Debug, Clone, PartialEq)]
pub struct AirPath {
    inlet_temperature_k: f64,
    mass_flow_kg_per_s: f64,
    specific_heat_j_per_kg_k: f64,
    segments: Vec<AirSegment>,
}

impl AirPath {
    /// Declare a path. Segments are ordered from inlet to outlet; each one's
    /// inlet is the previous one's outlet.
    ///
    /// # Errors
    /// [`AirflowError::EmptyAirPath`] for no segments;
    /// [`AirflowError::DuplicateAirSegment`] when one region appears twice —
    /// two segments cannot both own the same solid trace, because the solve
    /// gives that trace exactly one Robin row;
    /// [`AirflowError::InvalidConjugateInput`] for non-finite or non-positive
    /// inlet temperature, mass flow, or specific heat.
    pub fn new(
        inlet_temperature_k: f64,
        mass_flow_kg_per_s: f64,
        specific_heat_j_per_kg_k: f64,
        segments: Vec<AirSegment>,
    ) -> Result<AirPath, AirflowError> {
        if segments.is_empty() {
            return Err(AirflowError::EmptyAirPath);
        }
        for (index, segment) in segments.iter().enumerate() {
            if let Some(prior) = segments[..index]
                .iter()
                .find(|other| other.region == segment.region)
            {
                return Err(AirflowError::DuplicateAirSegment {
                    region: prior.region.clone(),
                });
            }
        }
        Ok(AirPath {
            inlet_temperature_k: finite_positive("air inlet temperature", inlet_temperature_k)?,
            mass_flow_kg_per_s: finite_positive("air mass flow", mass_flow_kg_per_s)?,
            specific_heat_j_per_kg_k: finite_positive(
                "air specific heat",
                specific_heat_j_per_kg_k,
            )?,
            segments,
        })
    }

    /// Path segments in inlet-to-outlet order.
    #[must_use]
    pub fn segments(&self) -> &[AirSegment] {
        &self.segments
    }

    /// Inlet air temperature, K.
    #[must_use]
    pub const fn inlet_temperature_k(&self) -> f64 {
        self.inlet_temperature_k
    }

    /// `ṁ c_p`, W/K.
    #[must_use]
    pub fn capacity_rate_w_per_k(&self) -> f64 {
        self.mass_flow_kg_per_s * self.specific_heat_j_per_kg_k
    }

    /// The region names this path binds, in order.
    #[must_use]
    pub fn regions(&self) -> Vec<&str> {
        self.segments.iter().map(AirSegment::region).collect()
    }

    /// March the path against one wall temperature per segment.
    ///
    /// # Errors
    /// [`AirflowError::SolidResponseArity`] when the wall-temperature count
    /// does not match the segment count;
    /// [`AirflowError::InvalidConjugateInput`] for a non-finite wall
    /// temperature;
    /// [`AirflowError::NonFiniteCoupling`] if the march itself produces a
    /// non-finite temperature or heat rate.
    pub fn march(&self, wall_temperatures_k: &[f64]) -> Result<AirMarch, AirflowError> {
        if wall_temperatures_k.len() != self.segments.len() {
            return Err(AirflowError::SolidResponseArity {
                expected: self.segments.len(),
                found: wall_temperatures_k.len(),
            });
        }
        let capacity = self.capacity_rate_w_per_k();
        let mut inlet = self.inlet_temperature_k;
        let mut states = Vec::with_capacity(self.segments.len());
        let mut total = 0.0;
        for (segment, &wall) in self.segments.iter().zip(wall_temperatures_k) {
            finite("solid wall temperature", wall)?;
            let ntu = segment.ntu(capacity);
            let (eps, eps_over_ntu) = effectiveness(ntu);
            let driving = wall - inlet;
            let heat_rate_w = capacity * driving * eps;
            let outlet = inlet + driving * eps;
            let reference = wall - driving * eps_over_ntu;
            for (stage, value) in [
                ("segment outlet temperature", outlet),
                ("segment reference temperature", reference),
                ("segment heat rate", heat_rate_w),
            ] {
                if !value.is_finite() {
                    return Err(AirflowError::NonFiniteCoupling {
                        stage,
                        value_bits: value.to_bits(),
                    });
                }
            }
            states.push(SegmentState {
                region: segment.region.clone(),
                inlet_temperature_k: inlet,
                outlet_temperature_k: outlet,
                reference_temperature_k: reference,
                ntu,
                effectiveness: eps,
                heat_rate_w,
            });
            total += heat_rate_w;
            inlet = outlet;
        }
        Ok(AirMarch {
            segments: states,
            outlet_temperature_k: inlet,
            total_heat_rate_w: total,
        })
    }
}

/// One named Robin region's response from the solid side of the exchange.
///
/// This is the shape [`fs_conduction::RobinFlux`] already has; see
/// [`SolidRegionState::from_robin_flux`].
#[derive(Debug, Clone, PartialEq)]
pub struct SolidRegionState {
    /// The declared Robin region name.
    pub region: String,
    /// Region area, m².
    pub area_m2: f64,
    /// Area-weighted mean wall temperature, K.
    pub mean_wall_temperature_k: f64,
    /// `∫ h (T_h − T_ref) dA` over the region, W (positive = leaving the
    /// solid).
    pub heat_rate_w: f64,
}

impl SolidRegionState {
    /// Adopt one solved Robin region from a conduction report.
    #[must_use]
    pub fn from_robin_flux(flux: &fs_conduction::RobinFlux) -> SolidRegionState {
        SolidRegionState {
            region: flux.region.clone(),
            area_m2: flux.area_m2,
            mean_wall_temperature_k: flux.mean_wall_temperature_k,
            heat_rate_w: flux.heat_rate_w,
        }
    }
}

/// How the outer fixed point is relaxed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Relaxation {
    /// A constant `ω`. `ω = 1` is plain Gauss-Seidel staggering.
    Fixed {
        /// The relaxation factor.
        omega: f64,
    },
    /// Scalar Aitken Δ² dynamic relaxation from [`fs_couple::AitkenRelaxation`],
    /// driven by the area-weighted mean reference-temperature residual.
    Aitken {
        /// Starting `ω`.
        omega_init: f64,
        /// Magnitude cap on `ω`.
        omega_max: f64,
    },
}

impl Default for Relaxation {
    fn default() -> Relaxation {
        Relaxation::Fixed { omega: 1.0 }
    }
}

/// The outer-loop stop rule.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConjugateConfig {
    /// Converged when the largest per-region reference-temperature update is
    /// at or below this, K.
    pub temperature_tolerance_k: f64,
    /// Outer iterations allowed before a typed non-convergence refusal.
    pub max_iterations: usize,
    /// The largest per-region interface imbalance admitted at convergence, W.
    ///
    /// The temperature criterion alone CANNOT bound this: a small reference
    /// change in kelvin only bounds a heat rate once multiplied by a
    /// conductance, so a converged temperature is not evidence of a closed
    /// interface. A run that meets the temperature criterion but not this one
    /// refuses with [`AirflowError::ConjugateBalanceUnclosed`] rather than
    /// returning an answer whose two sides disagree about the heat crossing
    /// the same surface.
    pub balance_tolerance_w: f64,
    /// Relaxation scheme.
    pub relaxation: Relaxation,
}

impl Default for ConjugateConfig {
    fn default() -> ConjugateConfig {
        ConjugateConfig {
            temperature_tolerance_k: 1.0e-9,
            max_iterations: 100,
            balance_tolerance_w: 1.0e-6,
            relaxation: Relaxation::default(),
        }
    }
}

/// One outer iteration's record.
#[derive(Debug, Clone, PartialEq)]
pub struct ConjugateIteration {
    /// Zero-based outer iteration index.
    pub iteration: usize,
    /// `max_j |T*_j − T_j|` before relaxation, K — the fixed-point residual.
    pub max_reference_change_k: f64,
    /// The signed area-weighted mean residual driving scalar relaxation, K.
    pub scalar_residual_k: f64,
    /// The `ω` applied this iteration.
    pub relaxation_omega: f64,
    /// Air temperature leaving the path this iteration, K.
    pub air_outlet_temperature_k: f64,
    /// `Σ_j Q_solid,j`, W.
    pub solid_heat_rate_w: f64,
    /// `Σ_j Q_air,j`, W.
    pub air_heat_rate_w: f64,
    /// `Σ_j Q_solid,j − Σ_j Q_air,j`, W.
    pub interface_imbalance_w: f64,
    /// The reference temperatures this iteration was solved against, K.
    pub reference_temperatures_k: Vec<f64>,
}

/// One region's two independently computed heat rates.
#[derive(Debug, Clone, PartialEq)]
pub struct RegionBalance {
    /// The region name.
    pub region: String,
    /// `∫ h (T_h − T_ref) dA` from the conduction solve, W.
    pub solid_heat_rate_w: f64,
    /// `ṁ c_p ΔT` from the air march, W.
    pub air_heat_rate_w: f64,
    /// `solid − air`, W.
    pub imbalance_w: f64,
}

/// The interface conservation audit.
///
/// Read the module-level "What the flux balance proves" section before citing
/// any field here: the per-region rows are a wiring falsifier, not a
/// conservation proof.
#[derive(Debug, Clone, PartialEq)]
pub struct FluxBalanceAudit {
    /// Per-region rows in path order.
    pub regions: Vec<RegionBalance>,
    /// `max_j |imbalance_j|`, W.
    pub max_region_imbalance_w: f64,
    /// `Σ_j Q_solid,j`, W.
    pub solid_total_w: f64,
    /// `Σ_j Q_air,j`, W.
    pub air_total_w: f64,
    /// `solid_total_w − air_total_w`, W.
    pub interface_imbalance_w: f64,
    /// `Σ_j Q_solid,j − robin_out_w` when the caller supplied the conduction
    /// report's whole-domain Robin total, W.
    ///
    /// `fs-conduction` accumulates `robin_out_w` in its own loop over every
    /// Robin face, so a nonzero residual here means the declared regions do
    /// not partition the Robin boundary — a face is unowned or double-counted.
    /// `None` when the caller did not supply the total.
    pub decomposition_residual_w: Option<f64>,
}

/// A converged conjugate exchange.
#[derive(Debug, Clone, PartialEq)]
pub struct ConjugateSolution {
    /// The converged reference temperatures in path order, K.
    pub reference_temperatures_k: Vec<f64>,
    /// The air march at the converged wall temperatures.
    pub march: AirMarch,
    /// The solid's converged per-region response.
    pub solid: Vec<SolidRegionState>,
    /// The interface conservation audit.
    pub balance: FluxBalanceAudit,
    /// Every outer iteration in order.
    pub history: Vec<ConjugateIteration>,
    /// Outer iterations performed.
    pub iterations: usize,
    /// The worst interface-power imbalance the [`fs_couple::EnergyAudit`] saw
    /// across every exchange, W. NaN if any exchange produced a NaN balance —
    /// the audit poisons rather than silently dropping it.
    pub worst_recorded_imbalance_w: f64,
}

impl ConjugateSolution {
    /// Attach the whole-domain Robin total from a conduction report and
    /// recompute the decomposition cross-check.
    #[must_use]
    pub fn with_decomposition_cross_check(mut self, robin_out_w: f64) -> ConjugateSolution {
        self.balance.decomposition_residual_w = Some(self.balance.solid_total_w - robin_out_w);
        self
    }
}

/// Solve the partitioned conjugate exchange.
///
/// `solid` is the caller-supplied residual callback: it receives the current
/// reference temperature for every path segment, in path order, and must run
/// the solid side against those Robin rows and return one
/// [`SolidRegionState`] per segment in the same order. `fs-couple` has no
/// coupling-loop driver and its contract assigns the driver to the consumer,
/// so this is that driver; nothing here calls a pre-existing loop.
///
/// Each outer iteration begins with a cancellation checkpoint, so a cancelled
/// run stops at an iteration boundary with the iteration index that had not
/// yet run. The returned [`ConjugateSolution::history`] is a complete replay
/// record of the exchange.
///
/// # Errors
/// [`AirflowError::ConjugateNotConverged`] when `max_iterations` is exhausted
/// without meeting the criterion — a run that ran out of iterations is a
/// refusal carrying its diagnosis, never a returned answer;
/// [`AirflowError::Cancelled`] at a checkpoint;
/// [`AirflowError::SolidResponseArity`] or
/// [`AirflowError::SegmentRegionMismatch`] when the callback's response does
/// not line up with the declared path;
/// [`AirflowError::InvalidConjugateInput`] for a malformed stop rule or
/// relaxation;
/// [`AirflowError::NonFiniteCoupling`] if the exchange produces a non-finite
/// quantity; plus whatever the callback itself returns.
pub fn solve_conjugate<F>(
    cx: &Cx<'_>,
    path: &AirPath,
    config: &ConjugateConfig,
    solid: F,
) -> Result<ConjugateSolution, AirflowError>
where
    F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
{
    // The exchange starts from the inlet: no heat has been absorbed yet, so
    // every segment's reference is the inlet temperature. This is the NTU → 0
    // limit of the segment law, so iteration 0 is a member of the same family
    // as every later one rather than a separately invented seed.
    let seed = vec![path.inlet_temperature_k(); path.segments().len()];
    solve_conjugate_from(cx, path, config, &seed, solid)
}

/// Resume a conjugate exchange from a retained reference-temperature vector.
///
/// Every outer iteration is a checkpoint boundary, and
/// [`ConjugateIteration::reference_temperatures_k`] is the complete state that
/// iteration was solved against — the driver keeps no other iteration-carried
/// state, and the relaxer is re-seeded from the config. Resuming from record
/// `k` therefore reproduces the uninterrupted tail bitwise under
/// [`Relaxation::Fixed`]. Under [`Relaxation::Aitken`] the relaxer's own
/// `ω` history is NOT carried across the resume, so the tail is a valid
/// continuation but not a bitwise replay; that limitation is the scalar
/// relaxer's, not the driver's.
///
/// # Errors
/// [`AirflowError::SolidResponseArity`] when `initial_references_k` does not
/// have one entry per segment; otherwise exactly what [`solve_conjugate`]
/// returns.
pub fn solve_conjugate_from<F>(
    cx: &Cx<'_>,
    path: &AirPath,
    config: &ConjugateConfig,
    initial_references_k: &[f64],
    mut solid: F,
) -> Result<ConjugateSolution, AirflowError>
where
    F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
{
    finite_positive(
        "conjugate temperature tolerance",
        config.temperature_tolerance_k,
    )?;
    finite_positive("conjugate balance tolerance", config.balance_tolerance_w)?;
    if initial_references_k.len() != path.segments().len() {
        return Err(AirflowError::SolidResponseArity {
            expected: path.segments().len(),
            found: initial_references_k.len(),
        });
    }
    for &value in initial_references_k {
        finite("resumed reference temperature", value)?;
    }
    if config.max_iterations == 0 {
        return Err(AirflowError::InvalidConjugateInput {
            field: "conjugate max iterations",
            value_bits: 0,
        });
    }
    let mut aitken = admit_relaxation(config.relaxation)?;

    let mut audit = EnergyAudit::new();
    let mut history: Vec<ConjugateIteration> = Vec::new();
    let mut reference = initial_references_k.to_vec();
    let total_area: f64 = path.segments().iter().map(AirSegment::area_m2).sum();

    for iteration in 0..config.max_iterations {
        cx.checkpoint().map_err(|_| AirflowError::Cancelled {
            iteration,
            references_k: reference.clone(),
        })?;

        let response = solid(cx, &reference)?;
        check_response_matches_path(path, &response)?;

        let walls: Vec<f64> = response
            .iter()
            .map(|state| state.mean_wall_temperature_k)
            .collect();
        let march = path.march(&walls)?;

        let solid_total: f64 = response.iter().map(|state| state.heat_rate_w).sum();
        let air_total = march.total_heat_rate_w;
        let interface_imbalance = solid_total - air_total;
        audit.record(interface_imbalance);

        let updated = march.reference_temperatures_k();
        let (max_change, scalar_residual) =
            fixed_point_residual(path, &reference, &updated, total_area)?;

        let omega = match (&config.relaxation, aitken.as_mut()) {
            (Relaxation::Fixed { omega }, _) => *omega,
            (Relaxation::Aitken { .. }, Some(relaxer)) => relaxer.next_omega(scalar_residual),
            // Unreachable: the relaxer is constructed iff the scheme is Aitken.
            (Relaxation::Aitken { omega_init, .. }, None) => *omega_init,
        };

        history.push(ConjugateIteration {
            iteration,
            max_reference_change_k: max_change,
            scalar_residual_k: scalar_residual,
            relaxation_omega: omega,
            air_outlet_temperature_k: march.outlet_temperature_k,
            solid_heat_rate_w: solid_total,
            air_heat_rate_w: air_total,
            interface_imbalance_w: interface_imbalance,
            reference_temperatures_k: reference.clone(),
        });

        if max_change <= config.temperature_tolerance_k {
            // Poll AFTER the exchange too. Checking only at iteration entry
            // lets a cancellation requested during the final solid solve
            // publish a success, which is precisely the drain-then-finalize
            // contract this is supposed to honour.
            cx.checkpoint().map_err(|_| AirflowError::Cancelled {
                iteration,
                references_k: reference.clone(),
            })?;
            return admit_converged(
                converged_solution(
                    updated,
                    march,
                    response,
                    history,
                    iteration + 1,
                    audit.max_generation(),
                ),
                config.balance_tolerance_w,
            );
        }

        for (slot, &next) in reference.iter_mut().zip(&updated) {
            *slot += omega * (next - *slot);
        }
    }

    let last = history
        .last()
        .expect("max_iterations > 0 pushes at least one record");
    Err(AirflowError::ConjugateNotConverged {
        iterations: config.max_iterations,
        max_change_bits: last.max_reference_change_k.to_bits(),
        tolerance_bits: config.temperature_tolerance_k.to_bits(),
    })
}

/// The CHT ladder's correlation-rung transfer, defined over the real coupled
/// state instead of an anonymous vector.
///
/// `fs-ladder` ships [`fs_ladder::Refine1d`] as a demonstrator: a generic 1-D
/// coarsen/refine-by-2 whose `restrict ∘ prolongate = identity` holds because
/// injection undoes interpolation, for any numbers at all. It knows nothing
/// about what it is moving. This transfer moves the seam's actual state — one
/// wall temperature per stream-wise segment — and carries the refined air path
/// that the finer rung would discretize, so the rung boundary is a statement
/// about the model rather than about array indices.
///
/// `prolongate` refines each coarse segment into [`Self::factor`] equal-area
/// sub-segments carrying that segment's wall temperature: the coarse rung's
/// own piecewise-constant-wall assumption, made explicit. `restrict` takes the
/// mean back. The two compose to the identity on the coarse space, and the
/// composition is *meaningful*: it says the coarse rung is the
/// piecewise-constant-wall projection of the finer one.
///
/// The property [`fs_ladder::Refine1d`] cannot state: prolongating a wall
/// state and re-marching the air on [`Self::refined_path`] reproduces the
/// coarse march's outlet temperature and total heat rate EXACTLY, because
/// `e^(−NTU) = (e^(−NTU/N))^N`. The transfer preserves the air-side quantity
/// of interest under refinement, with no tolerance. See `tests/conjugate.rs`.
///
/// # No-claim boundary
///
/// The finer side is a REFINED CORRELATION STATE, not a RANS field. No RANS
/// rung exists in this workspace (it is bead `f85xj.5.8`, gated on wedge
/// ratification), so this operator defines the correlation rung's own
/// refinement semantics and the state shape a RANS rung would have to accept.
/// It is not evidence that a RANS rung is implemented, and the `RANS`/`LES`
/// rung cost hints remain unmeasured declarations because those rungs do not
/// exist to measure.
#[derive(Debug, Clone, PartialEq)]
pub struct SegmentRefinementTransfer {
    path: AirPath,
    factor: usize,
}

impl SegmentRefinementTransfer {
    /// Bind a refinement transfer to a coarse path.
    ///
    /// # Errors
    /// [`AirflowError::InvalidConjugateInput`] when `factor` is zero.
    pub fn new(path: AirPath, factor: usize) -> Result<SegmentRefinementTransfer, AirflowError> {
        if factor == 0 {
            return Err(AirflowError::InvalidConjugateInput {
                field: "segment refinement factor",
                value_bits: 0,
            });
        }
        Ok(SegmentRefinementTransfer { path, factor })
    }

    /// Sub-segments each coarse segment refines into.
    #[must_use]
    pub const fn factor(&self) -> usize {
        self.factor
    }

    /// The coarse path this transfer refines.
    #[must_use]
    pub const fn coarse_path(&self) -> &AirPath {
        &self.path
    }

    /// The refined path: every coarse segment split into [`Self::factor`]
    /// equal-area sub-segments carrying its coefficient, named
    /// `<region>#<sub-index>`.
    ///
    /// # Errors
    /// Whatever [`AirPath::new`] refuses; a refined path is structurally
    /// admissible whenever the coarse one is.
    pub fn refined_path(&self) -> Result<AirPath, AirflowError> {
        let mut segments = Vec::with_capacity(self.path.segments().len() * self.factor);
        for segment in self.path.segments() {
            let area = segment.area_m2() / self.factor as f64;
            for sub in 0..self.factor {
                segments.push(AirSegment::new(
                    &format!("{}#{sub}", segment.region()),
                    area,
                    segment.htc_w_per_m2_k(),
                )?);
            }
        }
        AirPath::new(
            self.path.inlet_temperature_k(),
            self.path.mass_flow_kg_per_s,
            self.path.specific_heat_j_per_kg_k,
            segments,
        )
    }
}

impl fs_ladder::Transfer for SegmentRefinementTransfer {
    fn prolongate(&self, coarse: &[f64]) -> Vec<f64> {
        // The infallible `Transfer` contract cannot refuse, so a state whose
        // arity disagrees with the bound path POISONS at the correct output
        // arity instead of scaling garbage into a plausible fine state; the
        // downstream response gates refuse non-finite walls.
        let segment_count = self.path.segments().len();
        if coarse.len() != segment_count {
            return vec![f64::NAN; segment_count * self.factor];
        }
        let mut fine = Vec::with_capacity(coarse.len() * self.factor);
        for &wall in coarse {
            fine.extend(std::iter::repeat_n(wall, self.factor));
        }
        fine
    }

    fn restrict(&self, fine: &[f64]) -> Vec<f64> {
        // Same poisoning rule as `prolongate`: `zip` + `chunks` would
        // otherwise silently drop surplus downstream entries or restrict a
        // partial block with weights built for a full one.
        let segment_count = self.path.segments().len();
        if fine.len() != segment_count * self.factor {
            return vec![f64::NAN; segment_count];
        }
        // NOT a plain mean. The coarse wall that reproduces a block's outlet
        // is the DOWNSTREAM-WEIGHTED mean
        //
        //   W = (1 - a) * sum_i a^(f-1-i) W_i / (1 - a^f),   a = e^(-NTU/f)
        //
        // because air entering sub-segment i has already been warmed by
        // everything upstream, so a hot sub-segment near the outlet moves the
        // block's exit temperature more than an equally hot one near the
        // inlet. The plain mean is the a -> 1 (NTU -> 0) limit of this, and
        // coincides with it exactly when the block's walls are uniform --
        // which is every state `prolongate` emits, and is why a round-trip
        // test alone cannot see the difference. On a nonuniform block it is
        // wrong by tens of kelvin (measured 15.65 K at NTU 2 over 4
        // sub-segments), so restricting a genuinely refined fine state would
        // silently corrupt the coarse rung.
        let capacity = self.path.capacity_rate_w_per_k();
        self.path
            .segments()
            .iter()
            .zip(fine.chunks(self.factor))
            .map(|(segment, block)| {
                if block.is_empty() {
                    return f64::NAN;
                }
                let sub_ntu = segment.ntu(capacity) / self.factor as f64;
                let a = det::exp(-sub_ntu);
                let span = block.len();
                let numerator: f64 = block
                    .iter()
                    .enumerate()
                    .map(|(i, wall)| det::pow(a, (span - 1 - i) as f64) * wall)
                    .sum();
                let denominator: f64 = (0..span).map(|i| det::pow(a, i as f64)).sum();
                if denominator == 0.0 {
                    return f64::NAN;
                }
                numerator / denominator
            })
            .collect()
    }
}

/// The fixed-point residual: `(max_j |dT_j|, area-weighted mean dT)`.
///
/// The max drives the convergence test; the signed area-weighted mean is the
/// scalar the (scalar-only) Aitken relaxer needs, since a max of absolute
/// values carries no sign for a delta-squared step to work with.
fn fixed_point_residual(
    path: &AirPath,
    current: &[f64],
    updated: &[f64],
    total_area: f64,
) -> Result<(f64, f64), AirflowError> {
    let mut max_change = 0.0_f64;
    let mut weighted = 0.0_f64;
    for ((segment, &next), &now) in path.segments().iter().zip(updated).zip(current) {
        let delta = next - now;
        if !delta.is_finite() {
            return Err(AirflowError::NonFiniteCoupling {
                stage: "reference temperature update",
                value_bits: delta.to_bits(),
            });
        }
        max_change = max_change.max(delta.abs());
        weighted += segment.area_m2() * delta;
    }
    Ok((max_change, weighted / total_area))
}

/// The temperature criterion does not bound watts, so a converged exchange
/// must still close its interface before it is an answer.
///
/// NaN-safe by construction: `max_region_imbalance_w` poisons to NaN if any
/// region did, and `NaN.partial_cmp(_)` is `None`, so an unorderable
/// imbalance takes the refusal branch.
fn admit_converged(
    solution: ConjugateSolution,
    balance_tolerance_w: f64,
) -> Result<ConjugateSolution, AirflowError> {
    if matches!(
        solution
            .balance
            .max_region_imbalance_w
            .partial_cmp(&balance_tolerance_w),
        Some(Ordering::Less | Ordering::Equal)
    ) {
        return Ok(solution);
    }
    Err(AirflowError::ConjugateBalanceUnclosed {
        iterations: solution.iterations,
        max_region_imbalance_bits: solution.balance.max_region_imbalance_w.to_bits(),
        tolerance_bits: balance_tolerance_w.to_bits(),
    })
}

/// The solid side must answer with one region per segment, in path order.
/// Right count with wrong association would attribute one surface's heat to
/// another surface's air, which no downstream number would reveal.
fn check_response_matches_path(
    path: &AirPath,
    response: &[SolidRegionState],
) -> Result<(), AirflowError> {
    if response.len() != path.segments().len() {
        return Err(AirflowError::SolidResponseArity {
            expected: path.segments().len(),
            found: response.len(),
        });
    }
    for (index, (segment, state)) in path.segments().iter().zip(response).enumerate() {
        if segment.region() != state.region {
            return Err(AirflowError::SegmentRegionMismatch {
                index,
                expected: segment.region().to_string(),
                found: state.region.clone(),
            });
        }
        // Fail closed on the NUMBERS, not just the labels. A NaN heat rate
        // would otherwise reach the audit, where a plain `f64::max` fold
        // SILENTLY DROPS it (`f64::max(0.0, NaN) == 0.0`) and reports a
        // perfectly balanced interface for a coupling that broke down. This
        // is the same trap `fs_couple::EnergyAudit` documents and poisons
        // against; the per-region fold poisons too, and this gate stops the
        // value earlier.
        finite("solid wall temperature", state.mean_wall_temperature_k)?;
        finite("solid region heat rate", state.heat_rate_w)?;
        finite_positive("solid region area", state.area_m2)?;
        // The area the solid integrated over must be the area the air path
        // was declared with, or the two sides are exchanging heat across
        // different surfaces and every downstream watt is mis-scaled.
        let declared = segment.area_m2();
        if (state.area_m2 - declared).abs() > 1.0e-9 * declared {
            return Err(AirflowError::SegmentAreaMismatch {
                index,
                region: state.region.clone(),
                declared_bits: declared.to_bits(),
                reported_bits: state.area_m2.to_bits(),
            });
        }
    }
    Ok(())
}

/// Validate a relaxation scheme and build its relaxer, if it has one.
fn admit_relaxation(relaxation: Relaxation) -> Result<Option<AitkenRelaxation>, AirflowError> {
    match relaxation {
        Relaxation::Fixed { omega } => {
            finite("conjugate relaxation omega", omega)?;
            Ok(None)
        }
        Relaxation::Aitken {
            omega_init,
            omega_max,
        } => {
            finite("conjugate relaxation omega", omega_init)?;
            finite_positive("conjugate relaxation omega cap", omega_max)?;
            Ok(Some(AitkenRelaxation::new(omega_init, omega_max)))
        }
    }
}

/// Assemble the converged exchange, pairing each region's two independently
/// computed heat rates.
fn converged_solution(
    reference_temperatures_k: Vec<f64>,
    march: AirMarch,
    solid: Vec<SolidRegionState>,
    history: Vec<ConjugateIteration>,
    iterations: usize,
    worst_recorded_imbalance_w: f64,
) -> ConjugateSolution {
    let regions: Vec<RegionBalance> = solid
        .iter()
        .zip(&march.segments)
        .map(|(state, segment)| RegionBalance {
            region: state.region.clone(),
            solid_heat_rate_w: state.heat_rate_w,
            air_heat_rate_w: segment.heat_rate_w,
            imbalance_w: state.heat_rate_w - segment.heat_rate_w,
        })
        .collect();
    // Poison rather than fold: `f64::max` drops NaN, so a plain fold would
    // report zero imbalance for a coupling that produced NaN. `NaN <= tol` is
    // false, so every downstream comparison fails closed.
    let max_region_imbalance_w = if regions.iter().any(|row| row.imbalance_w.is_nan()) {
        f64::NAN
    } else {
        regions
            .iter()
            .map(|row| row.imbalance_w.abs())
            .fold(0.0_f64, f64::max)
    };
    let solid_total_w: f64 = solid.iter().map(|state| state.heat_rate_w).sum();
    let air_total_w = march.total_heat_rate_w;
    ConjugateSolution {
        reference_temperatures_k,
        balance: FluxBalanceAudit {
            regions,
            max_region_imbalance_w,
            solid_total_w,
            air_total_w,
            interface_imbalance_w: solid_total_w - air_total_w,
            decomposition_residual_w: None,
        },
        march,
        solid,
        history,
        iterations,
        worst_recorded_imbalance_w,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_qty::Temperature;

    #[test]
    fn seam_effort_is_a_temperature() {
        assert_eq!(seam_effort_dimensions(), Temperature::DIMS);
    }

    #[test]
    fn zero_ntu_reference_degenerates_to_the_inlet() {
        // Exact limits, so compare bits rather than within a margin: an
        // approximate check here would hide the 0/0 this branch exists to
        // avoid.
        let (eps, ratio) = effectiveness(0.0);
        assert_eq!(eps.to_bits(), 0.0_f64.to_bits());
        assert_eq!(ratio.to_bits(), 1.0_f64.to_bits());
    }

    #[test]
    fn effectiveness_ratio_stays_near_one_for_tiny_ntu() {
        let (_, ratio) = effectiveness(1.0e-12);
        assert!((ratio - 1.0).abs() < 1.0e-9, "ratio = {ratio}");
    }
}
