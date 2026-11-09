//! Reflection from authored correlated coating-thickness realizations.
//!
//! Each row is a mutually exclusive, independently illuminated column or a
//! quasistatic time realization, with the same incident phasor at the common
//! entrance plane. Original materials, frequency, angle and polarization stay
//! fixed; only the explicitly supplied finite-film thicknesses vary. Joint
//! rows preserve correlations between layers. No distribution or coherence
//! length is fitted, and rows are never replaced by a mean geometry.
//!
//! Mean reflected power is E[|r|²], not |E[r]|². The latter is the coherent
//! component at the common entrance plane; E[|r-E[r]|²] retains the remainder
//! without reclassifying it as heat. This independent-column/time model does
//! not solve lateral coupling, diffraction or generic rough-surface scattering.
//! Transmitted power is measured at each realization's own exit face: no
//! coherent transmitted field at a common downstream plane is inferred.
//!
//! Coherent films and their absorption use the original Maxwell operator.
//! The uniform round-trip-phase limit follows the intensity sums in
//! S. J. Byrnes, "Multilayer optical calculations", section 6,
//! <https://arxiv.org/abs/1603.02720>. Finite authored rows are a discrete
//! ensemble; continuous-distribution quadrature error is not certified here.

use super::{LayeredReflection, PlanarDielectricLayers, Polarization, ReflectionError, bad};
use fs_math::{c64::C64, dd::Dd};

/// One joint realization, in the original stack's physical layer order.
#[derive(Clone, Copy, Debug)]
pub struct LayerThicknessRealization<'a> {
    /// Finite nonnegative physical thickness of every original film [m].
    pub thicknesses_m: &'a [f64],
    /// Nonnegative probability or normalized independent-column area fraction.
    /// All rows must sum to one within binary64 rounding error.
    pub probability: f64,
}

/// Explicit preparation/evaluation bounds. No unbounded sampling is performed.
#[derive(Clone, Copy, Debug)]
pub struct ThicknessEnsembleBudget {
    /// Maximum number of authored rows, including zero-probability rows.
    pub maximum_realizations: usize,
    /// Maximum rows times max(original layer count, 1).
    pub maximum_layer_realizations: usize,
}

/// Reflected field statistics and independently accumulated native power.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThicknessEnsembleReflection {
    /// E[r], on the original common tangential-E axis and entrance plane.
    pub mean_reflected_tangential_e: C64,
    /// E[|r|²], the total mean reflected normal power fraction.
    pub mean_reflected_power_fraction: f64,
    /// |E[r]|², evaluated from the retained compensated field mean.
    pub coherent_reflected_power_fraction: f64,
    /// E[|r-E[r]|²]. This is reflected-field variance, never absorbed heat.
    pub fluctuating_reflected_power_fraction: f64,
    /// E[(R-E[R])²], for R=|r|². Units are squared power fraction.
    pub reflected_power_variance: f64,
    /// Mean outgoing normal flux at each realization's physical exit face.
    pub mean_transmitted_power_fraction: f64,
    /// Mean physical finite-film absorption from the original positive loss modes.
    pub mean_absorbed_power_fraction: f64,
    /// Largest absolute R+T+A-1 residual of any evaluated realization.
    pub maximum_absolute_power_balance_error: f64,
    /// Number of positive-probability rows actually solved.
    pub evaluated_realizations: usize,
}

/// Borrowed original laws and joint geometry samples. Preparation and queries
/// allocate nothing. Angle evaluation is O(rows * layers), using the native
/// fixed-size bounded slab workspace and one set of compensated moments.
pub struct PreparedLayerThicknessEnsemble<'r, 's, 'm> {
    source: &'s PlanarDielectricLayers<'m>,
    realizations: &'r [LayerThicknessRealization<'r>],
    layer_realization_work: usize,
}

impl<'m> PlanarDielectricLayers<'m> {
    /// Prepare an explicit probability/area-weighted thickness ensemble.
    ///
    /// Rows retain all cross-layer correlations and the original material
    /// source laws. Zero-thickness films remain valid identities; their source
    /// admission already belongs to the original prepared stack. Every row's
    /// geometry is checked, including zero-probability rows. Those dark rows
    /// do not require a scattering solve. Cancellation returns no partial result.
    pub fn prepare_thickness_ensemble<'r, 's>(
        &'s self,
        realizations: &'r [LayerThicknessRealization<'r>],
        budget: ThicknessEnsembleBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedLayerThicknessEnsemble<'r, 's, 'm>, ReflectionError> {
        poll(&mut cancelled)?;
        if realizations.is_empty() || realizations.len() > budget.maximum_realizations {
            return Err(bad(
                "thickness ensemble needs nonempty rows within its budget",
            ));
        }
        let layer_realization_work = realizations
            .len()
            .checked_mul(self.layer_count().max(1))
            .ok_or_else(|| bad("thickness ensemble work overflow"))?;
        if layer_realization_work > budget.maximum_layer_realizations {
            return Err(bad(
                "thickness ensemble exceeds its layer-realization budget",
            ));
        }
        let mut probability = Dd::ZERO;
        for row in realizations {
            poll(&mut cancelled)?;
            if row.thicknesses_m.len() != self.layer_count()
                || row.thicknesses_m.iter().any(|d| !d.is_finite() || *d < 0.0)
                || !row.probability.is_finite()
                || row.probability < 0.0
                || row.probability > 1.0
            {
                return Err(bad(
                    "thickness ensemble needs one finite nonnegative thickness per film and valid probabilities",
                ));
            }
            probability = probability + Dd::from_f64(row.probability);
        }
        if (probability - Dd::ONE).abs().to_f64() > 128.0 * f64::EPSILON {
            return Err(bad("thickness ensemble probabilities must sum to one"));
        }
        poll(&mut cancelled)?;
        Ok(PreparedLayerThicknessEnsemble {
            source: self,
            realizations,
            layer_realization_work,
        })
    }
}

impl<'r, 's, 'm> PreparedLayerThicknessEnsemble<'r, 's, 'm> {
    /// Original immutable material laws, frequency and nominal geometry.
    pub const fn source(&self) -> &'s PlanarDielectricLayers<'m> {
        self.source
    }
    /// Authored joint rows, retaining their order and correlations.
    pub const fn realizations(&self) -> &'r [LayerThicknessRealization<'r>] {
        self.realizations
    }
    /// Conservative per-query layer work, including zero-probability rows.
    pub const fn layer_realization_work(&self) -> usize {
        self.layer_realization_work
    }

    /// Solve each illuminated realization with the original coherent operator.
    ///
    /// No material query, material averaging, ray randomization or allocation
    /// occurs. Each realization independently passes the native power-balance
    /// check before contributing to the result. Weighted online central moments
    /// preserve weak fluctuations without subtracting nearly equal powers.
    /// Probabilities are normalized by their admitted roundoff-close sum.
    /// Cancellation is polled at every row and before publishing the result.
    pub fn evaluate(
        &self,
        incidence_cosine: f64,
        polarization: Polarization,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ThicknessEnsembleReflection, ReflectionError> {
        poll(&mut cancelled)?;
        let mut moments = Moments::new();
        for row in self.realizations {
            poll(&mut cancelled)?;
            if row.probability == 0.0 {
                continue;
            }
            let response = self.source.response_with_layer_thicknesses(
                incidence_cosine,
                polarization,
                Some(row.thicknesses_m),
            )?;
            moments.push(row.probability, response);
        }
        poll(&mut cancelled)?;
        moments.finish(self.source.layer_count())
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ReflectionError> {
    if cancelled() {
        Err(bad("thickness ensemble cancelled"))
    } else {
        Ok(())
    }
}

// Weighted Welford updates retain the field and intensity central moments in
// one pass. DD means preserve the small mean shift needed by later rows;
// fluctuating reflected power never comes from E[R] - |E[r]|².
struct Moments {
    weight: Dd,
    mean: [Dd; 5],
    field_variance_sum: Dd,
    power_variance_sum: Dd,
    maximum_balance_error: f64,
    count: usize,
}
impl Moments {
    fn new() -> Self {
        Self {
            weight: Dd::ZERO,
            mean: [Dd::ZERO; 5],
            field_variance_sum: Dd::ZERO,
            power_variance_sum: Dd::ZERO,
            maximum_balance_error: 0.0,
            count: 0,
        }
    }
    fn push(&mut self, probability: f64, response: LayeredReflection) {
        let weight = Dd::from_f64(probability);
        let total = self.weight + weight;
        let values = [
            response.reflected_tangential_e.re,
            response.reflected_tangential_e.im,
            response.reflected_power_fraction,
            response.transmitted_power_fraction,
            response.absorbed_power_fraction,
        ];
        let delta: [Dd; 5] = core::array::from_fn(|i| Dd::from_f64(values[i]) - self.mean[i]);
        let ratio = weight / total;
        let central_weight = weight * (self.weight / total);
        self.field_variance_sum =
            self.field_variance_sum + central_weight * (delta[0] * delta[0] + delta[1] * delta[1]);
        self.power_variance_sum = self.power_variance_sum + central_weight * delta[2] * delta[2];
        for (mean, difference) in self.mean.iter_mut().zip(delta) {
            *mean = *mean + ratio * difference;
        }
        let residual =
            (Dd::from_f64(values[2]) + Dd::from_f64(values[3]) + Dd::from_f64(values[4]) - Dd::ONE)
                .abs()
                .to_f64();
        self.maximum_balance_error = self.maximum_balance_error.max(residual);
        self.weight = total;
        self.count += 1;
    }
    fn finish(self, layers: usize) -> Result<ThicknessEnsembleReflection, ReflectionError> {
        let field = C64::new(self.mean[0].to_f64(), self.mean[1].to_f64());
        let coherent = (self.mean[0] * self.mean[0] + self.mean[1] * self.mean[1]).to_f64();
        let fluctuating = (self.field_variance_sum / self.weight).to_f64();
        let power_variance = (self.power_variance_sum / self.weight).to_f64();
        let [_, _, reflected, transmitted, absorbed] = self.mean.map(Dd::to_f64);
        let tolerance = 512.0 * f64::EPSILON * (layers + 1) as f64;
        if [
            coherent,
            fluctuating,
            power_variance,
            reflected,
            transmitted,
            absorbed,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
            || !field.re.is_finite()
            || !field.im.is_finite()
            || (reflected + transmitted + absorbed - 1.0).abs() > tolerance
            || (coherent + fluctuating - reflected).abs() > tolerance
        {
            return Err(bad(
                "thickness ensemble moments or power balance are unrepresentable",
            ));
        }
        Ok(ThicknessEnsembleReflection {
            mean_reflected_tangential_e: field,
            mean_reflected_power_fraction: reflected,
            coherent_reflected_power_fraction: coherent,
            fluctuating_reflected_power_fraction: fluctuating,
            reflected_power_variance: power_variance,
            mean_transmitted_power_fraction: transmitted,
            mean_absorbed_power_fraction: absorbed,
            maximum_absolute_power_balance_error: self.maximum_balance_error,
            evaluated_realizations: self.count,
        })
    }
}
