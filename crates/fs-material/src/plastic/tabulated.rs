//! Bounded, piecewise-linear isotropic work hardening on the native J2 flow.
//!
//! The table contains absolute yield stress versus accumulated equivalent
//! plastic strain, not engineering stress versus total tensile strain. Its
//! nondecreasing strength may have plateaus and changes of slope. There is no
//! extrapolation beyond the final measurement. Native small-strain J2 uses
//! Cauchy stress; a finite-strain consumer must explicitly select its own
//! stress convention rather than silently relabel a small-strain calibration.
//!
//! Hardening is explicitly energetic: storage is the exact integral of
//! `Y(alpha)-Y(0)`, and physical yield dissipation is `Y(0) delta_alpha`.
//! This partition is a model choice, not heat-fraction data inferred from a
//! tensile curve. The backward-Euler hardening remainder is numerical loss.

/// Original material-card curves with explicit stress and plastic-strain meaning.
pub mod source;
/// Associative Perzyna flow over the same bounded measured hardening curve.
pub mod rate;
/// Shared affine-temperature table storage and its exact entropy derivative.
pub mod thermal;
/// Original strain/temperature observations on the affine thermal finite owner.
pub mod thermal_source;

use super::energy::{EnergyError, EnergyIncrement, EnergyTolerances, StoredEnergy};
use super::{J2Plasticity, J2State, Update};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use core::fmt;
use fs_evidence::{Ambition, ModelCard};
use std::sync::Arc;

/// One absolute flow-stress datum in the consumer's declared stress convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HardeningPoint {
    /// Accumulated equivalent plastic strain, dimensionless; the first is zero.
    pub plastic_strain: f64,
    /// Positive absolute yield stress [Pa], nondecreasing with plastic strain.
    pub yield_stress_pa: f64,
}

/// Fixed material calibration with slopes and exact primitive prepared once.
/// Point queries and consistency returns allocate no memory. Immutable arrays
/// are shared so preparing another physical step never copies a measured curve.
#[derive(Clone, Debug, PartialEq)]
pub struct TabulatedHardening {
    points: Arc<[HardeningPoint]>,
    slopes: Arc<[f64]>,
    stored: Arc<[f64]>,
    maximum_slope_pa: f64,
}

/// Invalid data, a measured-domain exit, or unrepresentable arithmetic.
#[derive(Clone, Debug, PartialEq)]
pub enum TabulatedError {
    /// Invalid calibration or return-map coefficients.
    Parameters(&'static str),
    /// Inadmissible native material-point history.
    State(&'static str),
    /// The requested history requires extrapolation of the measured curve.
    Domain {
        /// Requested accumulated plastic strain, or an upper return estimate.
        plastic_strain: f64,
        /// Last measured accumulated plastic strain.
        maximum: f64,
    },
    /// The operation cannot be represented by finite native arithmetic.
    Numerical(&'static str),
}

impl fmt::Display for TabulatedError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parameters(s) => write!(f, "tabulated hardening parameters: {s}"),
            Self::State(s) => write!(f, "tabulated hardening history: {s}"),
            Self::Domain {
                plastic_strain,
                maximum,
            } => write!(
                f,
                "plastic strain {plastic_strain} is outside the tabulated domain [0, {maximum}]"
            ),
            Self::Numerical(s) => write!(f, "tabulated hardening arithmetic: {s}"),
        }
    }
}
impl std::error::Error for TabulatedError {}

impl TabulatedHardening {
    /// Admit at least two finite knots starting at zero. Strain is strictly
    /// increasing and yield stress is positive and nondecreasing. Softening
    /// requires a regularized damage law and is deliberately refused here.
    pub fn new(points: Vec<HardeningPoint>) -> Result<Self, TabulatedError> {
        if points.len() < 2 || points[0].plastic_strain != 0.0 {
            return Err(TabulatedError::Parameters(
                "at least two knots and a zero plastic-strain origin are required",
            ));
        }
        if points.iter().any(|p| {
            !p.plastic_strain.is_finite()
                || p.plastic_strain < 0.0
                || !p.yield_stress_pa.is_finite()
                || p.yield_stress_pa <= 0.0
        }) {
            return Err(TabulatedError::Parameters(
                "strains must be finite and nonnegative; strengths finite and positive",
            ));
        }
        let initial = points[0].yield_stress_pa;
        let mut slopes = Vec::with_capacity(points.len() - 1);
        let mut stored = Vec::with_capacity(points.len());
        stored.push(0.0);
        for pair in points.windows(2) {
            let width = pair[1].plastic_strain - pair[0].plastic_strain;
            let rise = pair[1].yield_stress_pa - pair[0].yield_stress_pa;
            if width <= 0.0 || rise < 0.0 {
                return Err(TabulatedError::Parameters(
                    "strains must increase strictly and yield stresses must not decrease",
                ));
            }
            let slope = rise / width;
            let area = (0.5 * (pair[0].yield_stress_pa - initial)
                + 0.5 * (pair[1].yield_stress_pa - initial))
                * width;
            let cumulative = stored[stored.len() - 1] + area;
            if !slope.is_finite()
                || (rise > 0.0 && slope == 0.0)
                || !cumulative.is_finite()
                || (pair[1].yield_stress_pa > initial && area == 0.0)
            {
                return Err(TabulatedError::Parameters(
                    "table slopes and hardening integrals must be representable",
                ));
            }
            slopes.push(slope);
            stored.push(cumulative);
        }
        let maximum_slope_pa = slopes.iter().copied().fold(0.0, f64::max);
        Ok(Self {
            points: points.into(),
            slopes: slopes.into(),
            stored: stored.into(),
            maximum_slope_pa,
        })
    }

    /// Original immutable measurement knots, including the domain endpoint.
    pub fn points(&self) -> &[HardeningPoint] {
        &self.points
    }

    /// Initial yield stress, the dissipative offset of the energetic model.
    pub fn initial_yield_stress_pa(&self) -> f64 {
        self.points[0].yield_stress_pa
    }

    /// Last measured accumulated plastic strain; extrapolation is refused.
    pub fn maximum_plastic_strain(&self) -> f64 {
        self.points[self.points.len() - 1].plastic_strain
    }

    /// Largest prepared slope, for admission of a consumer's return stiffness.
    pub fn maximum_tangent_pa(&self) -> f64 {
        self.maximum_slope_pa
    }

    fn segment(&self, alpha: f64) -> Result<usize, TabulatedError> {
        let maximum = self.maximum_plastic_strain();
        if !alpha.is_finite() || alpha < 0.0 || alpha > maximum {
            return Err(TabulatedError::Domain {
                plastic_strain: alpha,
                maximum,
            });
        }
        Ok(self
            .points
            .partition_point(|p| p.plastic_strain <= alpha)
            .saturating_sub(1)
            .min(self.slopes.len() - 1))
    }

    fn value_on_segment(&self, alpha: f64, segment: usize) -> f64 {
        let lower = self.points[segment];
        let upper = self.points[segment + 1];
        // Use the nearer endpoint, retaining exact supplied values at knots.
        if alpha - lower.plastic_strain <= upper.plastic_strain - alpha {
            self.slopes[segment].mul_add(alpha - lower.plastic_strain, lower.yield_stress_pa)
        } else {
            self.slopes[segment].mul_add(alpha - upper.plastic_strain, upper.yield_stress_pa)
        }
    }

    /// Interpolated absolute yield stress, with no out-of-domain extrapolation.
    pub fn yield_stress_pa(&self, alpha: f64) -> Result<f64, TabulatedError> {
        Ok(self.value_on_segment(alpha, self.segment(alpha)?))
    }

    /// The right-hand derivative at interior knots; the final knot uses the
    /// last measured segment. At a kink no unique classical tangent exists.
    pub fn tangent_pa(&self, alpha: f64) -> Result<f64, TabulatedError> {
        Ok(self.slopes[self.segment(alpha)?])
    }

    fn interval(&self, alpha: f64, delta: f64) -> Result<(usize, usize, f64), TabulatedError> {
        let first = self.segment(alpha)?;
        if !delta.is_finite() || delta < 0.0 {
            return Err(TabulatedError::State(
                "plastic increment must be finite and nonnegative",
            ));
        }
        let end = alpha + delta;
        let last = self.segment(end)?;
        if delta > 0.0 && end == alpha {
            return Err(TabulatedError::Numerical(
                "plastic increment is below history resolution",
            ));
        }
        Ok((first, last, end))
    }

    /// Evaluate the strength increment directly; do not subtract two large
    /// absolute yield stresses for a small plastic increment.
    pub fn strength_increment_pa(&self, alpha: f64, delta: f64) -> Result<f64, TabulatedError> {
        let (first, last, end) = self.interval(alpha, delta)?;
        if first == last {
            return Ok(self.slopes[first] * delta);
        }
        let first_end = self.points[first + 1];
        let last_start = self.points[last];
        Ok(self.slopes[first] * (first_end.plastic_strain - alpha)
            + (last_start.yield_stress_pa - first_end.yield_stress_pa)
            + self.slopes[last] * (end - last_start.plastic_strain))
    }

    /// Exact integral of `Y(p)-Y(0)` from zero to alpha, per reference volume.
    pub fn hardening_energy_j_m3(&self, alpha: f64) -> Result<f64, TabulatedError> {
        let segment = self.segment(alpha)?;
        let point = self.points[segment];
        let width = alpha - point.plastic_strain;
        let value = self.stored[segment]
            + ((point.yield_stress_pa - self.initial_yield_stress_pa())
                + 0.5 * self.slopes[segment] * width)
                * width;
        if !value.is_finite() {
            return Err(TabulatedError::Numerical(
                "hardening energy is not representable",
            ));
        }
        Ok(value)
    }

    /// Exact interval integral, without subtracting two accumulated energies.
    /// This is hardening storage only; add `Y(0)*delta` for total flow work.
    pub fn hardening_energy_increment_j_m3(
        &self,
        alpha: f64,
        delta: f64,
    ) -> Result<f64, TabulatedError> {
        let (first, last, end) = self.interval(alpha, delta)?;
        let mut position = alpha;
        let mut energy = 0.0;
        for segment in first..=last {
            let next = end.min(self.points[segment + 1].plastic_strain);
            let width = if first == last {
                delta
            } else {
                next - position
            };
            let start_strength = (self.points[segment].yield_stress_pa
                - self.initial_yield_stress_pa())
                + self.slopes[segment] * (position - self.points[segment].plastic_strain);
            energy += (start_strength + 0.5 * self.slopes[segment] * width) * width;
            position = next;
        }
        if !energy.is_finite() {
            return Err(TabulatedError::Numerical(
                "hardening interval integral is not representable",
            ));
        }
        Ok(energy)
    }

    /// Stable nonnegative `R(alpha+delta)*delta - integral R`, where R is
    /// hardening stress. This is the backward-Euler loss, not plastic heat.
    pub fn numerical_loss_j_m3(&self, alpha: f64, delta: f64) -> Result<f64, TabulatedError> {
        let (first, last, end) = self.interval(alpha, delta)?;
        let mut position = alpha;
        let mut result = 0.0;
        // Integration by parts: integral (p-alpha) R'(p) dp. Each term is
        // nonnegative even across plateaus or slope increases/decreases.
        for segment in first..=last {
            let next = end.min(self.points[segment + 1].plastic_strain);
            let width = if first == last {
                delta
            } else {
                next - position
            };
            result += self.slopes[segment] * width * ((position - alpha) + 0.5 * width);
            position = next;
        }
        if !result.is_finite() {
            return Err(TabulatedError::Numerical(
                "hardening integration loss is not representable",
            ));
        }
        Ok(result)
    }

    /// Solve `trial - correction*delta = Y(alpha+delta)` exactly on each
    /// linear segment. Binary search selects the segment; no iterative tensor
    /// solve or material-point allocation is needed. `correction` is 3 mu
    /// for finite J2, or 3 mu plus the Prager modulus for small-strain J2.
    /// A trial within the yield value's roundoff whose return is below the
    /// resolution of `alpha` gets no increment, so an accepted state
    /// re-evaluated at its own strain keeps its history.
    pub fn solve_increment(
        &self,
        trial_pa: f64,
        alpha: f64,
        correction_pa: f64,
    ) -> Result<f64, TabulatedError> {
        let first = self.segment(alpha)?;
        if !trial_pa.is_finite()
            || trial_pa < 0.0
            || !correction_pa.is_finite()
            || correction_pa <= 0.0
        {
            return Err(TabulatedError::Parameters(
                "trial stress and return modulus must be finite and nonnegative/positive",
            ));
        }
        let yield_old = self.value_on_segment(alpha, first);
        let excess = trial_pa - yield_old;
        if excess <= 0.0 {
            return Ok(0.0);
        }
        let last = self.points.len() - 1;
        let required = |index: usize| {
            correction_pa * (self.points[index].plastic_strain - alpha)
                + (self.points[index].yield_stress_pa - yield_old)
        };
        let last_required = required(last);
        // Independently assembled endpoint stresses can differ by a few ulps.
        // Admit only that consistency roundoff; the solved increment below
        // remains clipped to measured support, never extrapolated beyond it.
        let endpoint_tolerance = 32.0 * f64::EPSILON * excess.max(last_required);
        if excess > last_required && excess - last_required > endpoint_tolerance {
            return Err(TabulatedError::Domain {
                plastic_strain: alpha + excess / correction_pa,
                maximum: self.maximum_plastic_strain(),
            });
        }
        let mut low = first + 1;
        let mut high = last;
        while low < high {
            let middle = low + (high - low) / 2;
            if required(middle) < excess {
                low = middle + 1;
            } else {
                high = middle;
            }
        }
        let segment = low - 1;
        let start = alpha.max(self.points[segment].plastic_strain);
        let crossed = start - alpha;
        let spent = if crossed == 0.0 {
            0.0
        } else {
            required(segment)
        };
        let modulus = correction_pa + self.slopes[segment];
        // The searched segment is a mathematical bracket. Enforce it against
        // final-addition roundoff so a final-knot return cannot leave the table.
        let upper = self.points[segment + 1].plastic_strain - alpha;
        let mut delta = (crossed + (excess - spent) / modulus)
            .max(crossed)
            .min(upper);
        if alpha + delta > self.points[segment + 1].plastic_strain {
            delta = delta.next_down();
        }
        // An accepted state re-evaluated at its own strain assembles its trial
        // stress independently of the yield value and can land an ulp above
        // it. That return is below the resolution of `alpha`: the nearest
        // representable history is the accepted one, with no increment. Admit
        // only yield-stress roundoff; any larger unrepresentable return is
        // refused below.
        let unrepresentable = alpha + delta == alpha;
        if modulus.is_finite()
            && delta.is_finite()
            && delta > 0.0
            && unrepresentable
            && excess <= 32.0 * f64::EPSILON * yield_old
        {
            return Ok(0.0);
        }
        if !modulus.is_finite() || !delta.is_finite() || delta <= 0.0 || unrepresentable {
            return Err(TabulatedError::Numerical(
                "plastic consistency return is not representable",
            ));
        }
        self.segment(alpha + delta)?;
        Ok(delta)
    }
}

/// Native small-strain associative J2 with a bounded isotropic yield table.
/// Its independent linear Prager backstress retains the original J2 history.
#[derive(Clone, Debug, PartialEq)]
pub struct TabulatedJ2 {
    backbone: J2Plasticity,
    hardening: TabulatedHardening,
}

/// One checked trial; the caller commits the returned history after acceptance.
#[derive(Clone, Debug, PartialEq)]
pub struct TabulatedJ2Response {
    /// Candidate native tensor-Voigt Cauchy stress [Pa].
    pub stress_pa: Voigt,
    /// Derivative with respect to native tensor-Voigt strain coordinates [Pa].
    pub tangent_pa: Tangent6,
    /// Original native J2 history, with no transverse components discarded.
    pub state: J2State,
    /// Nonnegative equivalent plastic-strain increment.
    pub plastic_increment: f64,
}

impl TabulatedJ2 {
    /// The table supplies the entire isotropic yield radius. To prevent double
    /// counting, the backbone must have H_iso=0 and yield stress equal to the
    /// first table value. Its elastic and linear Prager coefficients remain.
    pub fn new(
        backbone: J2Plasticity,
        hardening: TabulatedHardening,
    ) -> Result<Self, MaterialError> {
        backbone
            .validate_energy_law()
            .map_err(|e| MaterialError::Parameters {
                what: e.to_string(),
            })?;
        let (_, mu) = backbone.elastic.lame();
        if backbone.hardening != 0.0
            || backbone.yield_stress != hardening.initial_yield_stress_pa()
            || !(3.0 * mu + backbone.kinematic_hardening + hardening.maximum_tangent_pa())
                .is_finite()
        {
            return Err(MaterialError::Parameters {
                what: "tabulated J2 needs zero linear isotropic hardening, matching initial yield and finite return stiffness".into(),
            });
        }
        Ok(Self {
            backbone,
            hardening,
        })
    }

    /// Original elastic, initial-yield and independent Prager coefficients.
    pub const fn backbone(&self) -> &J2Plasticity {
        &self.backbone
    }
    /// Original bounded absolute yield curve.
    pub const fn hardening(&self) -> &TabulatedHardening {
        &self.hardening
    }

    /// Exact scalar restriction used by both tensor and axial consumers.
    /// Pass 3 mu or E respectively; the independent Prager modulus is added here.
    pub fn plastic_increment(
        &self,
        trial_pa: f64,
        alpha: f64,
        elastic_correction_pa: f64,
    ) -> Result<f64, TabulatedError> {
        self.hardening.solve_increment(
            trial_pa,
            alpha,
            elastic_correction_pa + self.backbone.kinematic_hardening,
        )
    }

    /// Endpoint isotropic hardening slope; excludes the Prager translation.
    pub fn isotropic_hardening_derivative_pa(&self, alpha: f64) -> Result<f64, TabulatedError> {
        self.hardening.tangent_pa(alpha)
    }

    /// Native effective-stress yield residual; invalid table history gives NaN.
    pub fn yield_function_with_state(&self, stress: &Voigt, state: &J2State) -> f64 {
        let Ok(radius) = self.hardening.yield_stress_pa(state.alpha) else {
            return f64::NAN;
        };
        self.backbone.yield_function_with_state(stress, state)
            - (radius - self.backbone.yield_stress)
    }

    /// Admit an accepted endpoint without projecting or advancing its history.
    pub fn validate_committed_history(
        &self,
        strain: &Voigt,
        state: &J2State,
        tolerance_pa: f64,
    ) -> Result<(), EnergyError> {
        self.backbone.validate_overstress_history(strain, state)?;
        self.hardening.yield_stress_pa(state.alpha).map_err(|_| {
            EnergyError::Invalid("history leaves the tabulated plastic-strain domain")
        })?;
        if !tolerance_pa.is_finite() || tolerance_pa <= 0.0 {
            return Err(EnergyError::Invalid(
                "yield allowance must be finite and positive",
            ));
        }
        let elastic = core::array::from_fn(|i| strain[i] - state.plastic_strain[i]);
        let stress = self.backbone.elastic.stress(&elastic, &());
        let residual = self.yield_function_with_state(&stress, state);
        if !residual.is_finite() {
            return Err(EnergyError::Invalid(
                "tabulated yield residual is not representable",
            ));
        }
        if residual > tolerance_pa {
            return Err(EnergyError::YieldResidual {
                residual_pa: residual,
                tolerance_pa,
            });
        }
        Ok(())
    }

    fn return_parameters(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<(J2Plasticity, f64), TabulatedError> {
        self.backbone
            .validate_overstress_history(strain, state)
            .map_err(|_| TabulatedError::State("invalid native J2 strain or history"))?;
        let mut shifted = self.backbone.clone();
        shifted.yield_stress = self.hardening.yield_stress_pa(state.alpha)?;
        let trial = shifted.radial_return(strain, state, 0.0);
        if !trial.plastic {
            return Ok((shifted, 0.0));
        }
        let (_, mu) = shifted.elastic.lame();
        let delta = self.plastic_increment(trial.q_trial, state.alpha, 3.0 * mu)?;
        if delta == 0.0 {
            // A roundoff excess with no representable hardening increment.
            return Ok((shifted, 0.0));
        }
        let secant = self.hardening.strength_increment_pa(state.alpha, delta)? / delta;
        if !secant.is_finite() {
            return Err(TabulatedError::Numerical(
                "hardening secant is not representable",
            ));
        }
        let secant =
            self.bounded_return_secant(trial.q_trial - shifted.yield_stress, state.alpha, secant)?;
        Ok((shifted, secant))
    }

    // The scalar root and the native secant-based return round independently.
    // A last-knot root must remain inside measured support in EVERY native
    // stress/state/energy evaluation, including a perfect-flow plateau. Correct
    // only machine-scale consistency error by raising the return denominator;
    // do not project history or evaluate an extrapolated material coefficient.
    fn bounded_return_secant(
        &self,
        excess: f64,
        alpha: f64,
        secant: f64,
    ) -> Result<f64, TabulatedError> {
        let (_, mu) = self.backbone.elastic.lame();
        let correction = 3.0 * mu + self.backbone.kinematic_hardening;
        let denominator = correction + secant;
        if !secant.is_finite() || secant < 0.0 || !denominator.is_finite() {
            return Err(TabulatedError::Numerical(
                "tabulated return modulus is not representable",
            ));
        }
        let maximum = self.hardening.maximum_plastic_strain();
        if excess <= 0.0 || alpha + excess / denominator <= maximum {
            return Ok(secant);
        }
        let mut remaining = maximum - alpha;
        if alpha + remaining > maximum {
            remaining = remaining.next_down();
        }
        let residual = excess - denominator * remaining;
        let tolerance = 32.0 * f64::EPSILON * excess.max(denominator * remaining);
        if remaining <= 0.0 || !residual.is_finite() || residual.abs() > tolerance {
            return Err(TabulatedError::Domain {
                plastic_strain: alpha + excess / denominator,
                maximum,
            });
        }
        let mut bounded_denominator = (excess / remaining).next_up();
        for _ in 0..3 {
            let bounded_secant = (bounded_denominator - correction).max(secant);
            let total = correction + bounded_secant;
            if total.is_finite() && alpha + excess / total <= maximum {
                return Ok(bounded_secant);
            }
            bounded_denominator = bounded_denominator.next_up();
        }
        Err(TabulatedError::Numerical(
            "bounded native table return is not representable",
        ))
    }

    fn returned(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<(J2Plasticity, f64, Update), TabulatedError> {
        let (shifted, secant) = self.return_parameters(strain, state)?;
        let update = shifted.radial_return(strain, state, secant);
        self.hardening
            .yield_stress_pa(state.alpha + update.dgamma)?;
        if update.stress.iter().any(|x| !x.is_finite()) {
            return Err(TabulatedError::Numerical(
                "returned stress is not representable",
            ));
        }
        Ok((shifted, secant, update))
    }

    /// Checked native radial return and consistent tangent. No extrapolation
    /// or mutation of the accepted plastic state occurs on a rejected trial.
    pub fn response(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<TabulatedJ2Response, TabulatedError> {
        let (shifted, secant, update) = self.returned(strain, state)?;
        let slope = self.hardening.tangent_pa(state.alpha + update.dgamma)?;
        let tangent = self.backbone.tangent_from_return(&update, slope);
        let next = shifted.update_state_with_overstress(strain, state, secant);
        self.backbone
            .validate_overstress_history(strain, &next)
            .map_err(|_| TabulatedError::State("returned native J2 history is inadmissible"))?;
        if tangent.iter().flatten().any(|x| !x.is_finite()) {
            return Err(TabulatedError::Numerical(
                "consistent tangent is not representable",
            ));
        }
        Ok(TabulatedJ2Response {
            stress_pa: update.stress,
            tangent_pa: tangent,
            state: next,
            plastic_increment: update.dgamma,
        })
    }

    /// Fixed-history inspection; this does not perform another plastic return.
    pub fn stored_energy(
        &self,
        strain: &Voigt,
        state: &J2State,
    ) -> Result<StoredEnergy, EnergyError> {
        self.backbone.validate_overstress_history(strain, state)?;
        let elastic = core::array::from_fn(|i| strain[i] - state.plastic_strain[i]);
        let mut stored = self
            .backbone
            .energy_of(&elastic, state.alpha, &state.plastic_strain)?;
        stored.hardening_j_m3 = self
            .hardening
            .hardening_energy_j_m3(state.alpha)
            .map_err(|_| EnergyError::Invalid("stored history leaves the tabulated domain"))?;
        if !stored.total_j_m3().is_finite() {
            return Err(EnergyError::Invalid(
                "tabulated stored energy is not representable",
            ));
        }
        Ok(stored)
    }

    /// Original native work account with exact table storage and integration
    /// loss. An empirical yield curve does not calibrate the heat partition.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<EnergyIncrement, EnergyError> {
        self.validate_committed_history(old_strain, old_state, tolerances.yield_pa)?;
        self.advance_with_energy_from_admitted_history(old_strain, old_state, strain, tolerances)
    }

    /// Audit a thermal endpoint after its owner admitted the old history at the
    /// old temperature. Thermal softening may place that history outside the
    /// new yield surface; return it once without erasing its permanent strain.
    /// The public isothermal path still requires an initially admissible state.
    pub(crate) fn advance_with_energy_from_admitted_history(
        &self,
        old_strain: &Voigt,
        old_state: &J2State,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<EnergyIncrement, EnergyError> {
        self.backbone.validate_overstress_history(old_strain, old_state)?;
        self.hardening.yield_stress_pa(old_state.alpha).map_err(|_| {
            EnergyError::Invalid("old history leaves the measured hardening domain")
        })?;
        let (shifted, secant) = self.return_parameters(strain, old_state).map_err(|_| {
            EnergyError::Invalid(
                "tabulated plastic return leaves its domain or is not representable",
            )
        })?;
        let mut next = shifted
            .advance_with_overstress_energy(old_strain, old_state, strain, tolerances, secant)?;
        self.validate_committed_history(strain, &next.state, tolerances.yield_pa)?;
        let delta = next.state.alpha - old_state.alpha;
        let residual = self.yield_function_with_state(&next.stress_pa, &next.state);
        if !residual.is_finite() {
            return Err(EnergyError::Invalid("nonfinite tabulated yield residual"));
        }
        if residual > tolerances.yield_pa || delta > 0.0 && residual.abs() > tolerances.yield_pa {
            return Err(EnergyError::YieldResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.yield_pa,
            });
        }
        next.energy.initial.hardening_j_m3 = self
            .hardening
            .hardening_energy_j_m3(old_state.alpha)
            .map_err(|_| EnergyError::Invalid("invalid initial hardening energy"))?;
        next.energy.final_state.hardening_j_m3 = self
            .hardening
            .hardening_energy_j_m3(next.state.alpha)
            .map_err(|_| EnergyError::Invalid("invalid final hardening energy"))?;
        next.energy.physical_dissipation_j_m3 = self.backbone.yield_stress * delta;
        next.energy.numerical_dissipation_j_m3 += self
            .hardening
            .numerical_loss_j_m3(old_state.alpha, delta)
            .map_err(|_| EnergyError::Invalid("invalid tabulated integration loss"))?;
        next.energy.residual_j_m3 = next.energy.mechanical_work_j_m3
            - next.energy.stored_change_j_m3()
            - next.energy.physical_dissipation_j_m3
            - next.energy.numerical_dissipation_j_m3;
        if !next.energy.initial.total_j_m3().is_finite()
            || !next.energy.final_state.total_j_m3().is_finite()
            || !next.energy.physical_dissipation_j_m3.is_finite()
            || !next.energy.numerical_dissipation_j_m3.is_finite()
            || !next.energy.residual_j_m3.is_finite()
        {
            return Err(EnergyError::Invalid(
                "tabulated energy account is not representable",
            ));
        }
        if next.energy.residual_j_m3.abs() > tolerances.energy_j_m3 {
            return Err(EnergyError::EnergyResidual {
                residual_j_m3: next.energy.residual_j_m3,
                tolerance_j_m3: tolerances.energy_j_m3,
            });
        }
        Ok(next)
    }
}

impl SmallStrainLaw for TabulatedJ2 {
    type State = J2State;
    fn initial_state(&self) -> J2State {
        self.backbone.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &J2State) -> Voigt {
        self.returned(strain, state)
            .map_or([f64::NAN; 6], |(_, _, update)| update.stress)
    }
    fn tangent(&self, strain: &Voigt, state: &J2State) -> Tangent6 {
        let Ok((_, _, update)) = self.returned(strain, state) else {
            return [[f64::NAN; 6]; 6];
        };
        let Ok(slope) = self.hardening.tangent_pa(state.alpha + update.dgamma) else {
            return [[f64::NAN; 6]; 6];
        };
        self.backbone.tangent_from_return(&update, slope)
    }
    fn update_state(&self, strain: &Voigt, state: &J2State) -> J2State {
        let Ok((shifted, secant, _)) = self.returned(strain, state) else {
            return J2State {
                plastic_strain: [f64::NAN; 6],
                alpha: f64::NAN,
            };
        };
        shifted.update_state_with_overstress(strain, state, secant)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            failure_envelope: "J2 with a bounded tabulated isotropic yield radius and independent Prager backstress",
            ..self.backbone.admissibility()
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.j2-tabulated-plasticity",
            "0.1.0",
            Ambition::Solid,
            vec![
                "small strain, associative rate-independent J2".into(),
                "bounded nondecreasing tabulated energetic isotropic hardening".into(),
                "independent linear Prager kinematic hardening".into(),
            ],
            self.backbone.card().validity.with(
                "accumulated-plastic-strain",
                0.0,
                self.hardening.maximum_plastic_strain(),
            ),
            vec![
                "no extrapolation, damage, softening or heat-fraction calibration".into(),
                "no finite-strain, thermal-coefficient evolution or nonlinear backstress recovery"
                    .into(),
            ],
            0.05,
        )
    }
}
