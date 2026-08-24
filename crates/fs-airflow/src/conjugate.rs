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

fn finite_coupling(stage: &'static str, value: f64) -> Result<f64, AirflowError> {
    if !value.is_finite() {
        return Err(AirflowError::NonFiniteCoupling {
            stage,
            value_bits: value.to_bits(),
        });
    }
    Ok(value)
}

fn finite_sum(
    stage: &'static str,
    values: impl IntoIterator<Item = f64>,
) -> Result<f64, AirflowError> {
    let mut total = 0.0;
    for value in values {
        total = finite_coupling(stage, total + value)?;
    }
    Ok(total)
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
    /// non-positive area or coefficient, or when their positive product is not
    /// representable as a finite, positive conductance.
    pub fn new(
        region: &str,
        area_m2: f64,
        htc_w_per_m2_k: f64,
    ) -> Result<AirSegment, AirflowError> {
        if region.trim().is_empty() {
            return Err(AirflowError::EmptyElementName);
        }
        let area_m2 = finite_positive("air segment area", area_m2)?;
        let htc_w_per_m2_k =
            finite_positive("air segment heat-transfer coefficient", htc_w_per_m2_k)?;
        finite_positive("air segment conductance", area_m2 * htc_w_per_m2_k)?;
        Ok(AirSegment {
            region: region.to_string(),
            area_m2,
            htc_w_per_m2_k,
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
    ///
    /// This public compatibility helper is infallible raw arithmetic and does
    /// not admit an arbitrary caller-supplied capacity rate. It can therefore
    /// return zero or a non-finite value for an invalid or unrepresentable
    /// ratio. [`AirPath::new`] and [`AirPath::march`] use the checked exchange
    /// path instead.
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

fn admitted_exchange_terms(
    segment: &AirSegment,
    capacity_rate_w_per_k: f64,
) -> Result<(f64, f64, f64), AirflowError> {
    let conductance = finite_positive(
        "air segment conductance",
        segment.htc_w_per_m2_k * segment.area_m2,
    )?;
    let ntu = finite_positive("air segment NTU", conductance / capacity_rate_w_per_k)?;
    let (effectiveness, effectiveness_over_ntu) = effectiveness(ntu);
    Ok((
        ntu,
        finite_positive("air segment effectiveness", effectiveness)?,
        finite_positive(
            "air segment effectiveness-to-NTU ratio",
            effectiveness_over_ntu,
        )?,
    ))
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
    /// inlet temperature, mass flow, or specific heat, or when the derived
    /// capacity rate, total area, conductance/NTU, or effectiveness terms are
    /// not representable in their required finite, positive domain.
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
        let path = AirPath {
            inlet_temperature_k: finite_positive("air inlet temperature", inlet_temperature_k)?,
            mass_flow_kg_per_s: finite_positive("air mass flow", mass_flow_kg_per_s)?,
            specific_heat_j_per_kg_k: finite_positive(
                "air specific heat",
                specific_heat_j_per_kg_k,
            )?,
            segments,
        };
        let capacity = path.checked_capacity_rate_w_per_k()?;
        path.checked_total_area_m2()?;
        for segment in &path.segments {
            admitted_exchange_terms(segment, capacity)?;
        }
        Ok(path)
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
    ///
    /// Construction admits the product as finite and strictly positive, so
    /// every `AirPath` returned by [`AirPath::new`] satisfies this invariant.
    #[must_use]
    pub fn capacity_rate_w_per_k(&self) -> f64 {
        self.mass_flow_kg_per_s * self.specific_heat_j_per_kg_k
    }

    fn checked_capacity_rate_w_per_k(&self) -> Result<f64, AirflowError> {
        finite_positive("air capacity rate", self.capacity_rate_w_per_k())
    }

    fn checked_total_area_m2(&self) -> Result<f64, AirflowError> {
        let mut total = 0.0;
        for segment in &self.segments {
            total = finite_positive("total air segment area", total + segment.area_m2)?;
        }
        Ok(total)
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
        let capacity = self.checked_capacity_rate_w_per_k()?;
        let mut inlet = self.inlet_temperature_k;
        let mut states = Vec::with_capacity(self.segments.len());
        let mut total = 0.0;
        for (segment, &wall) in self.segments.iter().zip(wall_temperatures_k) {
            finite("solid wall temperature", wall)?;
            let (ntu, eps, eps_over_ntu) = admitted_exchange_terms(segment, capacity)?;
            let driving = finite_coupling("segment temperature difference", wall - inlet)?;
            let heat_rate_w = finite_coupling("segment heat rate", capacity * driving * eps)?;
            let outlet = finite_coupling("segment outlet temperature", inlet + driving * eps)?;
            let reference = finite_coupling(
                "segment reference temperature",
                wall - driving * eps_over_ntu,
            )?;
            states.push(SegmentState {
                region: segment.region.clone(),
                inlet_temperature_k: inlet,
                outlet_temperature_k: outlet,
                reference_temperature_k: reference,
                ntu,
                effectiveness: eps,
                heat_rate_w,
            });
            total = finite_coupling("total air heat rate", total + heat_rate_w)?;
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
    /// Area-weighted mean `T_ref` the solid actually solved against, K.
    ///
    /// `None` makes no claim and skips the wiring cross-check. A supplied
    /// value must agree with the reference this driver sent for the region,
    /// or the exchange refuses with
    /// [`AirflowError::ReferenceTemperatureMismatch`]: the watt gate alone
    /// cannot see a reference-wiring error smaller than `hA` times its
    /// tolerance, while this check is denominated directly in kelvin.
    /// [`fs_conduction::RobinFlux`] returns the value for free, so
    /// [`SolidRegionState::from_robin_flux`] always claims it.
    pub mean_reference_temperature_k: Option<f64>,
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
            mean_reference_temperature_k: Some(flux.mean_reference_temperature_k),
        }
    }
}

/// The largest fixed relaxation factor this driver admits.
///
/// `ω = 0` never moves the reference and stalls by construction — a config
/// fault the non-convergence refusal would otherwise misattribute — and
/// nothing above classical over-relaxation has convergence theory for this
/// staggered seam, so both refuse at admission instead of burning
/// `max_iterations` of solid solves first.
const MAX_FIXED_OMEGA: f64 = 2.0;

/// How the outer fixed point is relaxed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Relaxation {
    /// A constant `ω`. `ω = 1` is plain Gauss-Seidel staggering.
    Fixed {
        /// The relaxation factor, admitted in `(0, 2]`.
        omega: f64,
    },
    /// Scalar Aitken Δ² dynamic relaxation from [`fs_couple::AitkenRelaxation`],
    /// driven by the area-weighted mean reference-temperature residual.
    Aitken {
        /// Starting `ω`, admitted in `(0, omega_max]`.
        omega_init: f64,
        /// Magnitude cap on `ω`. The secant update may drive `ω` negative;
        /// the cap bounds its magnitude and is the caller's own risk budget.
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
    /// The absolute floor of the per-region interface-imbalance gate, W.
    ///
    /// The temperature criterion alone CANNOT bound this: a small reference
    /// change in kelvin only bounds a heat rate once multiplied by a
    /// conductance, so a converged temperature is not evidence of a closed
    /// interface. A run that meets the temperature criterion but not the
    /// balance gate refuses with [`AirflowError::ConjugateBalanceUnclosed`]
    /// rather than returning an answer whose two sides disagree about the
    /// heat crossing the same surface.
    ///
    /// This is a FLOOR, not the whole gate: the admitted imbalance is
    /// `max(balance_tolerance_w, balance_relative_tolerance * scale)` with
    /// `scale = max_j max(|Q_solid,j|, |Q_air,j|)`. The floor exists so an
    /// all-but-zero-power interface (scale ≈ 0) still admits floating-point
    /// dust instead of demanding exact zeros.
    pub balance_tolerance_w: f64,
    /// The relative half of the imbalance gate, dimensionless in `[0, 1)`.
    ///
    /// An absolute watt threshold alone is scale-dependent in both failure
    /// directions: at milliwatt scale a dropped-face wiring fault produces an
    /// imbalance far below any fixed constant chosen for watt-scale runs and
    /// PASSES, while at kilowatt scale legitimate floating-point residue can
    /// exceed the same constant and spuriously refuses. Scaling the gate by
    /// the interface's own heat-rate magnitude makes its sensitivity
    /// scale-free: wiring faults sit at `O(scale)` and floating-point residue
    /// at `O(1e-12 · scale)`, so a per-million threshold separates them at
    /// every wattage. `0.0` disables the relative term and gates on the
    /// absolute floor alone.
    pub balance_relative_tolerance: f64,
    /// Relaxation scheme.
    pub relaxation: Relaxation,
}

impl Default for ConjugateConfig {
    fn default() -> ConjugateConfig {
        ConjugateConfig {
            temperature_tolerance_k: 1.0e-9,
            max_iterations: 100,
            balance_tolerance_w: 1.0e-12,
            balance_relative_tolerance: 1.0e-6,
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
    /// The worst finite interface-power imbalance the
    /// [`fs_couple::EnergyAudit`] saw across every admitted exchange, W.
    /// Non-finite exchange arithmetic refuses before solution assembly;
    /// the audit's poisoning remains a defensive invariant for independent
    /// callers.
    pub worst_recorded_imbalance_w: f64,
}

impl ConjugateSolution {
    /// Attach the whole-domain Robin total from a conduction report and
    /// recompute the decomposition cross-check.
    ///
    /// This compatibility surface is infallible. A non-finite total or
    /// subtraction therefore poisons the optional residual with `NaN`; callers
    /// that need a structured refusal must validate both operands and the
    /// representability of their subtraction before calling.
    #[must_use]
    pub fn with_decomposition_cross_check(mut self, robin_out_w: f64) -> ConjugateSolution {
        let residual = self.balance.solid_total_w - robin_out_w;
        self.balance.decomposition_residual_w = Some(if residual.is_finite() {
            residual
        } else {
            f64::NAN
        });
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
/// [`AirflowError::ReferenceTemperatureMismatch`] when the callback claims a
/// solved-against reference that is not the one this driver sent;
/// [`AirflowError::InvalidConjugateInput`] for a malformed stop rule or
/// relaxation (an inert `ω ≤ 0` or `ω` beyond over-relaxation refuses here,
/// not as a non-convergence);
/// [`AirflowError::NonFiniteCoupling`] if the exchange produces a non-finite
/// quantity, including a reference vector the relaxation itself overflowed;
/// plus whatever the callback itself returns.
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
    admit_conjugate_inputs(path, config, initial_references_k)?;
    let mut aitken = admit_relaxation(config.relaxation)?;

    let mut audit = EnergyAudit::new();
    let mut history: Vec<ConjugateIteration> = Vec::new();
    let mut reference = initial_references_k.to_vec();
    let total_area = path.checked_total_area_m2()?;

    for iteration in 0..config.max_iterations {
        cx.checkpoint().map_err(|_| AirflowError::Cancelled {
            iteration,
            references_k: reference.clone(),
        })?;

        let response = solid(cx, &reference)?;
        check_response_matches_path(path, &reference, &response)?;

        let walls: Vec<f64> = response
            .iter()
            .map(|state| state.mean_wall_temperature_k)
            .collect();
        let march = path.march(&walls)?;

        let solid_total = finite_sum(
            "total solid heat rate",
            response.iter().map(|state| state.heat_rate_w),
        )?;
        let air_total = finite_coupling("total air heat rate", march.total_heat_rate_w)?;
        let interface_imbalance = finite_coupling(
            "air-solid interface heat-rate imbalance",
            solid_total - air_total,
        )?;
        audit.record(interface_imbalance);

        let updated = march.reference_temperatures_k();
        let (max_change, scalar_residual) =
            fixed_point_residual(path, &reference, &updated, total_area)?;

        let omega = finite_coupling(
            "relaxation omega",
            match (&config.relaxation, aitken.as_mut()) {
                (Relaxation::Fixed { omega }, _) => *omega,
                (Relaxation::Aitken { .. }, Some(relaxer)) => relaxer.next_omega(scalar_residual),
                // Unreachable: the relaxer is constructed iff the scheme is Aitken.
                (Relaxation::Aitken { omega_init, .. }, None) => *omega_init,
            },
        )?;

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
                )?,
                config,
            );
        }

        relax_references(&mut reference, &updated, omega)?;
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

/// Upfront refusals for the conjugate driver (stage of
/// [`solve_conjugate_from`]): tolerances, arity, resumed-reference
/// finiteness, and the zero-iteration configuration fault.
fn admit_conjugate_inputs(
    path: &AirPath,
    config: &ConjugateConfig,
    initial_references_k: &[f64],
) -> Result<(), AirflowError> {
    finite_positive(
        "conjugate temperature tolerance",
        config.temperature_tolerance_k,
    )?;
    finite_positive("conjugate balance tolerance", config.balance_tolerance_w)?;
    // NaN fails the range test, so this refusal is also the finiteness gate.
    // 1.0 would admit a fully open interface: refuse it as a config fault.
    if !(0.0..1.0).contains(&config.balance_relative_tolerance) {
        return Err(AirflowError::InvalidConjugateInput {
            field: "conjugate relative balance tolerance",
            value_bits: config.balance_relative_tolerance.to_bits(),
        });
    }
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
    Ok(())
}

/// Under-relax the reference vector toward `updated` (stage of
/// [`solve_conjugate_from`]).
///
/// Attribute an overflow to the stage that produced it. Without this, the
/// poisoned reference is first refused deep inside the caller's solid
/// callback (or by the response gate), blaming the solid for a
/// relaxation-configuration fault.
fn relax_references(
    reference: &mut [f64],
    updated: &[f64],
    omega: f64,
) -> Result<(), AirflowError> {
    for (slot, &next) in reference.iter_mut().zip(updated) {
        *slot += omega * (next - *slot);
    }
    for &slot in reference.iter() {
        if !slot.is_finite() {
            return Err(AirflowError::NonFiniteCoupling {
                stage: "relaxed reference temperature",
                value_bits: slot.to_bits(),
            });
        }
    }
    Ok(())
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

/// Cross-fidelity transfer between 1D correlation segment wall states and
/// resolved 2D/boundary-layer RANS temperature fields (bead `frankensim-wd35h`).
///
/// Prolongates a 1D vector of segment wall temperatures into a multi-node normal
/// temperature profile $T_i(y_j)$ across the viscous sublayer to the core stream,
/// and restricts a resolved near-wall RANS temperature field back to the segment wall
/// temperatures.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelationToRansFieldTransfer {
    segment_count: usize,
    normal_nodes: usize,
    inlet_temperature_k: f64,
}

impl CorrelationToRansFieldTransfer {
    /// Create a new cross-fidelity transfer between 1D correlation segments and RANS normal profiles.
    ///
    /// # Errors
    /// [`AirflowError::InvalidConjugateInput`] if `segment_count == 0` or `normal_nodes < 2`.
    pub fn new(
        segment_count: usize,
        normal_nodes: usize,
        inlet_temperature_k: f64,
    ) -> Result<Self, AirflowError> {
        if segment_count == 0 {
            return Err(AirflowError::InvalidConjugateInput {
                field: "segment_count",
                value_bits: 0,
            });
        }
        if normal_nodes < 2 {
            return Err(AirflowError::InvalidConjugateInput {
                field: "normal_nodes",
                value_bits: normal_nodes as u64,
            });
        }
        Ok(Self {
            segment_count,
            normal_nodes,
            inlet_temperature_k,
        })
    }

    /// Number of normal nodes per stream-wise segment in the fine RANS representation.
    #[must_use]
    pub const fn normal_nodes(&self) -> usize {
        self.normal_nodes
    }
}

impl fs_ladder::Transfer for CorrelationToRansFieldTransfer {
    fn prolongate(&self, coarse: &[f64]) -> Vec<f64> {
        if coarse.len() != self.segment_count {
            return vec![f64::NAN; self.segment_count * self.normal_nodes];
        }
        let mut fine = Vec::with_capacity(self.segment_count * self.normal_nodes);
        for &t_wall in coarse {
            fine.push(t_wall); // y=0 (wall node)
            for node in 1..self.normal_nodes {
                let frac = node as f64 / (self.normal_nodes - 1) as f64;
                // Linear/exponential boundary layer profile to core stream
                let t_field = t_wall + frac * (self.inlet_temperature_k - t_wall);
                fine.push(t_field);
            }
        }
        fine
    }

    fn restrict(&self, fine: &[f64]) -> Vec<f64> {
        if fine.len() != self.segment_count * self.normal_nodes {
            return vec![f64::NAN; self.segment_count];
        }
        let mut coarse = Vec::with_capacity(self.segment_count);
        for segment_idx in 0..self.segment_count {
            let wall_idx = segment_idx * self.normal_nodes;
            coarse.push(fine[wall_idx]);
        }
        coarse
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
    finite_positive("total air segment area", total_area)?;
    let mut max_change = 0.0_f64;
    let mut weighted = 0.0_f64;
    for ((segment, &next), &now) in path.segments().iter().zip(updated).zip(current) {
        let delta = finite_coupling("reference temperature update", next - now)?;
        max_change = max_change.max(delta.abs());
        let contribution = finite_coupling(
            "area-weighted reference residual contribution",
            segment.area_m2() * delta,
        )?;
        weighted = finite_coupling(
            "area-weighted reference residual accumulation",
            weighted + contribution,
        )?;
    }
    Ok((
        max_change,
        finite_coupling("area-weighted reference residual", weighted / total_area)?,
    ))
}

/// The temperature criterion does not bound watts, so a converged exchange
/// must still close its interface before it is an answer.
///
/// The admitted imbalance is the HYBRID threshold
/// `max(balance_tolerance_w, balance_relative_tolerance * scale)` with
/// `scale = max_j max(|Q_solid,j|, |Q_air,j|)` taken from the same response
/// being judged, so the gate's sensitivity does not depend on whether the
/// interface moves microwatts or kilowatts. The refusal reports the
/// EFFECTIVE threshold, not the configured floor.
///
/// Final assembly has already refused any non-finite region imbalance or
/// aggregate before this gate receives a solution. The explicit finite checks
/// on the derived scale and threshold remain defensive, so this comparison
/// only judges ordered finite values.
fn admit_converged(
    solution: ConjugateSolution,
    config: &ConjugateConfig,
) -> Result<ConjugateSolution, AirflowError> {
    let scale = finite_coupling(
        "interface heat-rate scale",
        solution
            .balance
            .regions
            .iter()
            .map(|row| row.solid_heat_rate_w.abs().max(row.air_heat_rate_w.abs()))
            .fold(0.0_f64, f64::max),
    )?;
    let threshold = finite_coupling(
        "interface imbalance threshold",
        config
            .balance_tolerance_w
            .max(config.balance_relative_tolerance * scale),
    )?;
    if matches!(
        solution
            .balance
            .max_region_imbalance_w
            .partial_cmp(&threshold),
        Some(Ordering::Less | Ordering::Equal)
    ) {
        return Ok(solution);
    }
    Err(AirflowError::ConjugateBalanceUnclosed {
        iterations: solution.iterations,
        max_region_imbalance_bits: solution.balance.max_region_imbalance_w.to_bits(),
        tolerance_bits: threshold.to_bits(),
    })
}

/// The solid side must answer with one region per segment, in path order.
/// Right count with wrong association would attribute one surface's heat to
/// another surface's air, which no downstream number would reveal.
fn check_response_matches_path(
    path: &AirPath,
    sent_references_k: &[f64],
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
        // Fail closed on the NUMBERS, not just the labels. A non-finite state
        // is attributed and refused here before it can enter aggregate
        // arithmetic, the energy audit, or final solution assembly.
        // `fs_couple::EnergyAudit` retains its own poisoning behavior as a
        // defensive invariant for independent callers.
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
        // When the solid claims which reference it solved against, hold it to
        // the reference this driver actually sent. The watt gate cannot see a
        // reference-wiring error smaller than hA times its tolerance; this
        // check is denominated directly in kelvin. The driver sends one
        // uniform reference per region, so the area-weighted mean must equal
        // it up to accumulation rounding. `None` makes no claim.
        if let Some(reported) = state.mean_reference_temperature_k {
            finite("solid reference temperature", reported)?;
            let sent = sent_references_k[index];
            if (reported - sent).abs() > 1.0e-9 * sent.abs().max(1.0) {
                return Err(AirflowError::ReferenceTemperatureMismatch {
                    index,
                    region: state.region.clone(),
                    sent_bits: sent.to_bits(),
                    reported_bits: reported.to_bits(),
                });
            }
        }
    }
    Ok(())
}

/// Validate a relaxation scheme and build its relaxer, if it has one.
///
/// Inert (`ω ≤ 0`) and beyond-over-relaxation factors refuse HERE, with a
/// config-attributed diagnostic, because admitting them produces a guaranteed
/// or near-guaranteed stall that [`AirflowError::ConjugateNotConverged`]
/// would then blame on convergence rather than on the configuration.
fn admit_relaxation(relaxation: Relaxation) -> Result<Option<AitkenRelaxation>, AirflowError> {
    match relaxation {
        Relaxation::Fixed { omega } => {
            finite("conjugate relaxation omega", omega)?;
            if !(omega > 0.0 && omega <= MAX_FIXED_OMEGA) {
                return Err(AirflowError::InvalidConjugateInput {
                    field: "conjugate relaxation omega",
                    value_bits: omega.to_bits(),
                });
            }
            Ok(None)
        }
        Relaxation::Aitken {
            omega_init,
            omega_max,
        } => {
            finite("conjugate relaxation omega", omega_init)?;
            finite_positive("conjugate relaxation omega cap", omega_max)?;
            if !(omega_init > 0.0 && omega_init <= omega_max) {
                return Err(AirflowError::InvalidConjugateInput {
                    field: "conjugate relaxation omega",
                    value_bits: omega_init.to_bits(),
                });
            }
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
) -> Result<ConjugateSolution, AirflowError> {
    let mut regions = Vec::with_capacity(solid.len());
    let mut max_region_imbalance_w = 0.0_f64;
    for (state, segment) in solid.iter().zip(&march.segments) {
        let imbalance_w = finite_coupling(
            "converged region interface heat-rate imbalance",
            state.heat_rate_w - segment.heat_rate_w,
        )?;
        max_region_imbalance_w = max_region_imbalance_w.max(imbalance_w.abs());
        regions.push(RegionBalance {
            region: state.region.clone(),
            solid_heat_rate_w: state.heat_rate_w,
            air_heat_rate_w: segment.heat_rate_w,
            imbalance_w,
        });
    }
    let solid_total_w = finite_sum(
        "converged solid heat-rate total",
        solid.iter().map(|state| state.heat_rate_w),
    )?;
    let air_total_w = finite_sum(
        "converged air heat-rate total",
        march.segments.iter().map(|state| state.heat_rate_w),
    )?;
    let interface_imbalance_w = finite_coupling(
        "converged interface heat-rate imbalance",
        solid_total_w - air_total_w,
    )?;
    let worst_recorded_imbalance_w = finite_coupling(
        "recorded interface heat-rate imbalance",
        worst_recorded_imbalance_w,
    )?;
    Ok(ConjugateSolution {
        reference_temperatures_k,
        balance: FluxBalanceAudit {
            regions,
            max_region_imbalance_w,
            solid_total_w,
            air_total_w,
            interface_imbalance_w,
            decomposition_residual_w: None,
        },
        march,
        solid,
        history,
        iterations,
        worst_recorded_imbalance_w,
    })
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

    #[test]
    fn residual_division_overflow_keeps_its_stage_attribution() {
        let path = AirPath::new(
            300.0,
            1.0,
            1.0,
            vec![AirSegment::new("wall", 1.0, 1.0).expect("segment")],
        )
        .expect("path");
        let error = fixed_point_residual(&path, &[0.0], &[1.0], f64::from_bits(1))
            .expect_err("an unrepresentable weighted mean must refuse");
        assert!(matches!(
            error,
            AirflowError::NonFiniteCoupling {
                stage: "area-weighted reference residual",
                ..
            }
        ));
    }

    fn synthetic_march(heat_rates: &[f64]) -> AirMarch {
        AirMarch {
            segments: heat_rates
                .iter()
                .enumerate()
                .map(|(index, &heat_rate_w)| SegmentState {
                    region: format!("row-{index}"),
                    inlet_temperature_k: 300.0,
                    outlet_temperature_k: 300.0,
                    reference_temperature_k: 300.0,
                    ntu: 1.0,
                    effectiveness: 0.5,
                    heat_rate_w,
                })
                .collect(),
            outlet_temperature_k: 300.0,
            total_heat_rate_w: 0.0,
        }
    }

    fn synthetic_solid(heat_rates: &[f64]) -> Vec<SolidRegionState> {
        heat_rates
            .iter()
            .enumerate()
            .map(|(index, &heat_rate_w)| SolidRegionState {
                region: format!("row-{index}"),
                area_m2: 1.0,
                mean_wall_temperature_k: 300.0,
                heat_rate_w,
                mean_reference_temperature_k: Some(300.0),
            })
            .collect()
    }

    #[test]
    fn converged_reductions_refuse_at_the_exact_recomputed_stage() {
        let solid_overflow = converged_solution(
            vec![300.0; 2],
            synthetic_march(&[f64::MAX, f64::MAX]),
            synthetic_solid(&[f64::MAX, f64::MAX]),
            Vec::new(),
            1,
            0.0,
        )
        .expect_err("the recomputed solid total must refuse");
        assert!(matches!(
            solid_overflow,
            AirflowError::NonFiniteCoupling {
                stage: "converged solid heat-rate total",
                ..
            }
        ));

        let air_overflow = converged_solution(
            vec![300.0; 2],
            synthetic_march(&[f64::MAX, f64::MAX]),
            synthetic_solid(&[0.0, 0.0]),
            Vec::new(),
            1,
            0.0,
        )
        .expect_err("the recomputed air total must refuse");
        assert!(matches!(
            air_overflow,
            AirflowError::NonFiniteCoupling {
                stage: "converged air heat-rate total",
                ..
            }
        ));

        let row = f64::MAX * 0.45;
        let interface_overflow = converged_solution(
            vec![300.0; 2],
            synthetic_march(&[-row, -row]),
            synthetic_solid(&[row, row]),
            Vec::new(),
            1,
            0.0,
        )
        .expect_err("the recomputed interface difference must refuse");
        assert!(matches!(
            interface_overflow,
            AirflowError::NonFiniteCoupling {
                stage: "converged interface heat-rate imbalance",
                ..
            }
        ));
    }
}
