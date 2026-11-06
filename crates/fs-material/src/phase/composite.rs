//! Fixed-composition sensible-matrix / phase-change-material calorics.
//!
//! The constituents share temperature (local thermal equilibrium), retain their
//! original masses, and have additive volumes and enthalpy, with no reaction or
//! heat of mixing. Their original source conditions and pressure/work convention
//! must be compatible; this constructor does not infer those conditions. Density
//! changes imply free expansion at fixed mass, not creation of cell inventory.
//! Confinement, stress, transport, evolving geometry and pressure work remain with
//! their existing owners. This law can be registered directly on the thermal
//! timeline or its implicit residual through [`ThermalLaw`].
//!
//! Only the PCM can change phase. Its liquid fraction is a fraction of PCM mass;
//! the composite's liquid mass fraction is that value times the PCM mass fraction.
//! A fully melted PCM therefore leaves the matrix solid. No normalized
//! solid/liquid chart or fictitious all-liquid endpoint is constructed.

use super::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, PhaseIntervalDensityLaw, PhaseStateError,
};
use crate::thermal::solid::PureSolidThermalLaw;
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// Preparation ceiling on the combined original knot counts.
pub const MAX_PCM_COMPOSITE_KNOTS: usize = 8192;

/// One exact caloric breakpoint and its two original enthalpy coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PcmCompositeKnot {
    /// Composite specific enthalpy [J/kg of composite].
    pub specific_enthalpy_j_kg: f64,
    /// Common constituent temperature [K].
    pub temperature_k: f64,
    /// Original matrix enthalpy, in its unchanged datum [J/kg of matrix].
    pub matrix_specific_enthalpy_j_kg: f64,
    /// Original PCM enthalpy, in its unchanged datum [J/kg of PCM].
    pub pcm_specific_enthalpy_j_kg: f64,
}

/// One composite state, including the unchanged constituent evaluations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PcmCompositeState {
    /// Composite specific enthalpy [J/kg of composite].
    pub specific_enthalpy_j_kg: f64,
    /// Common constituent temperature [K].
    pub temperature_k: f64,
    /// State evaluated by the original all-solid matrix chart.
    pub matrix: EnthalpyPhaseKnot,
    /// State evaluated by the original PCM chart; its liquid fraction is local
    /// to PCM mass, not the total composite mass.
    pub pcm: EnthalpyPhaseKnot,
    /// Fixed PCM mass divided by fixed total mass.
    pub pcm_mass_fraction: f64,
    /// Current PCM volume divided by current total volume.
    pub pcm_volume_fraction: f64,
    /// Liquid mass divided by total mass, `pcm_mass_fraction * pcm.f_liquid`.
    pub liquid_mass_fraction: f64,
    /// Additive constituent specific volumes [m3/kg of composite].
    pub specific_volume_m3_kg: f64,
    /// Reciprocal of additive specific volume [kg/m3].
    pub bulk_density_kg_m3: f64,
}

/// Original state and its exact one-sided enthalpy derivatives.
/// Interior knots select the right interval; the upper endpoint selects left.
/// These generalized derivatives do not smooth a latent plateau or fit density.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PcmCompositeResponse {
    pub state: PcmCompositeState,
    /// dT/dH [kg K/J], exactly zero on a sharp PCM plateau.
    pub temperature_slope_kg_k_j: f64,
    /// d(h_matrix)/dH in the original constituent enthalpy coordinate.
    pub matrix_enthalpy_derivative: f64,
    /// d(h_PCM)/dH in the original constituent enthalpy coordinate.
    pub pcm_enthalpy_derivative: f64,
    /// Derivative of PCM volume fraction with respect to composite H [kg/J].
    pub pcm_volume_fraction_slope_kg_j: f64,
    /// Derivative of liquid mass / total composite mass [kg/J].
    pub liquid_mass_fraction_slope_kg_j: f64,
}

/// Immutable, prepared law for a sensible solid matrix containing one PCM.
///
/// Preparation merges the original caloric temperature breakpoints. Both source
/// charts already use piecewise-affine `T(h)`, so this preserves their admitted
/// calorics exactly in real arithmetic, including every latent plateau. Original
/// component enthalpy mappings are retained: queries evaluate the original density
/// interpolation laws and entropy integrals instead of fitting mixture density or
/// averaging heat capacities. Queries allocate nothing and take logarithmic time.
///
/// The matrix's all-solid prefix must have a unique sensible inverse and cover
/// the entire PCM temperature range. The PCM must run from zero to unit liquid
/// fraction. Two phase-changing constituents are intentionally outside this law:
/// coincident isothermal plateaus would require an additional phase-allocation
/// closure to determine density from total enthalpy.
#[derive(Clone, Debug)]
pub struct PcmComposite<'a> {
    matrix: PureSolidThermalLaw<'a>,
    pcm: &'a EquilibriumEnthalpyPhaseCurve,
    pcm_mass_fraction: f64,
    knots: Vec<PcmCompositeKnot>,
    maximum_slope: f64,
    identity: ContentHash,
}

impl<'a> PcmComposite<'a> {
    /// Prepare a fixed mass fraction `0 < pcm_mass_fraction < 1`.
    ///
    /// `lower_specific_enthalpy_j_kg` is the explicitly chosen COMPOSITE datum
    /// at the PCM's lowest enthalpy/temperature. With `w` the PCM mass fraction,
    /// the compiled coordinate is `H = H_lower + (1-w)*(h_m-h_m0) + w*(h_p-h_p0)`.
    /// Independently shifting either constituent's arbitrary enthalpy datum
    /// therefore leaves the physical law unchanged (subject to representability).
    /// The original curves, source receipts and source identities remain borrowed.
    ///
    /// No source extrapolation, fraction normalization or temperature clamping is
    /// performed. Unrepresentable ordered breakpoints or slopes refuse at prepare.
    pub fn try_new(
        matrix: &'a EquilibriumEnthalpyPhaseCurve,
        pcm: &'a EquilibriumEnthalpyPhaseCurve,
        pcm_mass_fraction: f64,
        lower_specific_enthalpy_j_kg: f64,
    ) -> Result<Self, PhaseStateError> {
        let invalid = |what| PhaseStateError::InvalidCurve { what };
        if !pcm_mass_fraction.is_finite() || pcm_mass_fraction <= 0.0 || pcm_mass_fraction >= 1.0 {
            return Err(invalid(
                "PCM composite mass fraction must be strictly between zero and one",
            ));
        }
        if !lower_specific_enthalpy_j_kg.is_finite() {
            return Err(PhaseStateError::NonFiniteSpecificEnthalpy);
        }
        let count = matrix.knots().len().checked_add(pcm.knots().len());
        if count.is_none_or(|count| count > MAX_PCM_COMPOSITE_KNOTS) {
            return Err(invalid(
                "PCM composite preparation exceeds its original-knot ceiling",
            ));
        }
        let matrix = PureSolidThermalLaw::new(matrix)?;
        let matrix_knots = matrix.curve().knots();
        let solid_count = matrix_knots.partition_point(|knot| knot.liquid_mass_fraction == 0.0);
        let matrix_knots = &matrix_knots[..solid_count];
        if matrix_knots
            .windows(2)
            .any(|pair| pair[1].temperature_k <= pair[0].temperature_k)
        {
            return Err(invalid(
                "PCM composite matrix must have a unique sensible solid inverse",
            ));
        }
        let pcm_knots = pcm.knots();
        let first = pcm_knots[0];
        let last = pcm_knots[pcm_knots.len() - 1];
        if first.liquid_mass_fraction != 0.0 || last.liquid_mass_fraction != 1.0 {
            return Err(invalid(
                "PCM constituent must cover its complete solid-to-liquid transition",
            ));
        }
        let matrix_h0 = matrix.specific_enthalpy_at_temperature(first.temperature_k)?;
        matrix.specific_enthalpy_at_temperature(last.temperature_k)?;

        let mut knots = Vec::with_capacity(count.unwrap_or(0));
        for knot in pcm_knots {
            knots.push(PcmCompositeKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: knot.temperature_k,
                matrix_specific_enthalpy_j_kg: matrix
                    .specific_enthalpy_at_temperature(knot.temperature_k)?,
                pcm_specific_enthalpy_j_kg: knot.specific_enthalpy_j_kg,
            });
        }
        // A matrix breakpoint inside a PCM sensible interval splits that interval.
        // At a PCM temperature knot, every plateau endpoint is already retained.
        for knot in matrix_knots {
            let temperature = knot.temperature_k;
            if temperature <= first.temperature_k || temperature >= last.temperature_k {
                continue;
            }
            let upper = pcm_knots.partition_point(|knot| knot.temperature_k < temperature);
            if pcm_knots[upper].temperature_k == temperature {
                continue;
            }
            let a = pcm_knots[upper - 1];
            let b = pcm_knots[upper];
            let alpha = (temperature - a.temperature_k) / (b.temperature_k - a.temperature_k);
            let h = lerp(a.specific_enthalpy_j_kg, b.specific_enthalpy_j_kg, alpha);
            if !h.is_finite() || h <= a.specific_enthalpy_j_kg || h >= b.specific_enthalpy_j_kg {
                return Err(invalid(
                    "matrix temperature breakpoint has no representable PCM enthalpy",
                ));
            }
            let state = pcm.values_at_specific_enthalpy(h)?;
            require_same_temperature(temperature, state.temperature_k)?;
            knots.push(PcmCompositeKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: temperature,
                matrix_specific_enthalpy_j_kg: knot.specific_enthalpy_j_kg,
                pcm_specific_enthalpy_j_kg: h,
            });
        }
        knots.sort_by(|a, b| {
            a.pcm_specific_enthalpy_j_kg
                .total_cmp(&b.pcm_specific_enthalpy_j_kg)
        });
        for knot in &mut knots {
            let matrix_increment = knot.matrix_specific_enthalpy_j_kg - matrix_h0;
            let pcm_increment = knot.pcm_specific_enthalpy_j_kg - first.specific_enthalpy_j_kg;
            let increment =
                (1.0 - pcm_mass_fraction) * matrix_increment + pcm_mass_fraction * pcm_increment;
            knot.specific_enthalpy_j_kg = lower_specific_enthalpy_j_kg + increment;
            if !matrix_increment.is_finite()
                || !pcm_increment.is_finite()
                || !knot.specific_enthalpy_j_kg.is_finite()
            {
                return Err(invalid(
                    "PCM composite enthalpy increments must be representable",
                ));
            }
        }
        let mut maximum_slope = 0.0_f64;
        for pair in knots.windows(2) {
            let dh = pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg;
            let dt = pair[1].temperature_k - pair[0].temperature_k;
            if !dh.is_finite()
                || dh <= 0.0
                || pair[1].pcm_specific_enthalpy_j_kg <= pair[0].pcm_specific_enthalpy_j_kg
            {
                return Err(invalid(
                    "PCM composite caloric breakpoints are not resolvably ordered",
                ));
            }
            let slope = dt / dh;
            if !slope.is_finite() || (dt > 0.0 && slope == 0.0) {
                return Err(PhaseStateError::UnrepresentableTemperatureDerivative);
            }
            maximum_slope = maximum_slope.max(slope);
        }
        let mut identity = DomainHasher::new("org.frankensim.fs-material.pcm-composite.v1");
        identity.update(matrix.curve().identity().as_bytes());
        identity.update(pcm.identity().as_bytes());
        identity.update(&pcm_mass_fraction.to_bits().to_le_bytes());
        identity.update(&lower_specific_enthalpy_j_kg.to_bits().to_le_bytes());
        Ok(Self {
            matrix,
            pcm,
            pcm_mass_fraction,
            knots,
            maximum_slope,
            identity: identity.finalize(),
        })
    }

    /// Original matrix source, including any unused higher-temperature branch.
    #[must_use]
    pub const fn matrix_curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve {
        self.matrix.curve()
    }

    /// Original PCM source, including its source receipts and density laws.
    #[must_use]
    pub const fn pcm_curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve {
        self.pcm
    }

    /// Fixed PCM fraction of total mass.
    #[must_use]
    pub const fn pcm_mass_fraction(&self) -> f64 {
        self.pcm_mass_fraction
    }

    /// Compiled caloric breakpoints; no density is refitted onto these knots.
    #[must_use]
    pub fn knots(&self) -> &[PcmCompositeKnot] {
        &self.knots
    }

    /// Inclusive composite enthalpy domain in the supplied composite datum.
    #[must_use]
    pub fn enthalpy_range_j_kg(&self) -> [f64; 2] {
        [
            self.knots[0].specific_enthalpy_j_kg,
            self.knots[self.knots.len() - 1].specific_enthalpy_j_kg,
        ]
    }

    /// Identity binds both original curves, mass fraction and composite datum.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Resolve the original constituent states and their additive-volume density.
    pub fn values_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<PcmCompositeState, PhaseStateError> {
        let mapped = self.map_enthalpy(h)?;
        let matrix = self
            .matrix
            .curve()
            .values_at_specific_enthalpy(mapped.matrix_specific_enthalpy_j_kg)?;
        let pcm = self
            .pcm
            .values_at_specific_enthalpy(mapped.pcm_specific_enthalpy_j_kg)?;
        require_same_temperature(mapped.temperature_k, matrix.temperature_k)?;
        require_same_temperature(mapped.temperature_k, pcm.temperature_k)?;
        let matrix_volume = (1.0 - self.pcm_mass_fraction) / matrix.bulk_density_kg_m3;
        let pcm_volume = self.pcm_mass_fraction / pcm.bulk_density_kg_m3;
        let specific_volume = matrix_volume + pcm_volume;
        let density = specific_volume.recip();
        if matrix.liquid_mass_fraction != 0.0
            || matrix_volume <= 0.0
            || pcm_volume <= 0.0
            || !specific_volume.is_finite()
            || !density.is_finite()
            || density <= 0.0
        {
            return Err(PhaseStateError::InvalidCurve {
                what: "PCM composite specific volume must be finite and positive",
            });
        }
        Ok(PcmCompositeState {
            specific_enthalpy_j_kg: h,
            temperature_k: mapped.temperature_k,
            matrix,
            pcm,
            pcm_mass_fraction: self.pcm_mass_fraction,
            pcm_volume_fraction: pcm_volume / specific_volume,
            liquid_mass_fraction: self.pcm_mass_fraction * pcm.liquid_mass_fraction,
            specific_volume_m3_kg: specific_volume,
            bulk_density_kg_m3: density,
        })
    }

    /// Common temperature on the original admitted constituent charts [K].
    pub fn temperature_at_specific_enthalpy(&self, h: f64) -> Result<f64, PhaseStateError> {
        self.values_at_specific_enthalpy(h)
            .map(|state| state.temperature_k)
    }

    /// Evaluate original state and analytic caloric/volume-fraction tangents.
    ///
    /// The matrix stores no sensible heat at fixed melting temperature:
    /// d(h_matrix)/dH is zero and d(h_PCM)/dH is the inverse PCM mass fraction.
    /// Volume fraction can still evolve because original phase densities differ.
    /// Queries allocate nothing and search only the original knot arrays.
    pub fn response_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<PcmCompositeResponse, PhaseStateError> {
        let state = self.values_at_specific_enthalpy(h)?;
        let upper = self
            .knots
            .partition_point(|knot| knot.specific_enthalpy_j_kg <= h)
            .min(self.knots.len() - 1);
        let (a, b) = (self.knots[upper - 1], self.knots[upper]);
        let dh = b.specific_enthalpy_j_kg - a.specific_enthalpy_j_kg;
        let slope = |difference: f64| -> Result<f64, PhaseStateError> {
            let value = difference / dh;
            if !value.is_finite() || (difference != 0.0 && value == 0.0) {
                return Err(composite_derivative_error());
            }
            Ok(value)
        };
        let temperature_slope = slope(b.temperature_k - a.temperature_k)?;
        let matrix_slope =
            slope(b.matrix_specific_enthalpy_j_kg - a.matrix_specific_enthalpy_j_kg)?;
        let pcm_slope = slope(b.pcm_specific_enthalpy_j_kg - a.pcm_specific_enthalpy_j_kg)?;
        let w = self.pcm_mass_fraction;
        let matrix_volume_slope = (1.0 - w)
            * constituent_volume_slope(
                self.matrix.curve(),
                a.matrix_specific_enthalpy_j_kg,
                state.matrix.bulk_density_kg_m3,
                matrix_slope,
            )?;
        let pcm_volume_slope = w * constituent_volume_slope(
            self.pcm,
            a.pcm_specific_enthalpy_j_kg,
            state.pcm.bulk_density_kg_m3,
            pcm_slope,
        )?;
        let phi = state.pcm_volume_fraction;
        let fraction_slope = ((1.0 - phi) * pcm_volume_slope - phi * matrix_volume_slope)
            / state.specific_volume_m3_kg;
        // Select the original interval containing this complete composite span;
        // this also selects left correctly at the composite's upper endpoint.
        let original = self.pcm.knots();
        let upper = original
            .partition_point(|knot| knot.specific_enthalpy_j_kg <= a.pcm_specific_enthalpy_j_kg)
            .min(original.len() - 1);
        let (pa, pb) = (original[upper - 1], original[upper]);
        let liquid_slope = w
            * pcm_slope
            * ((pb.liquid_mass_fraction - pa.liquid_mass_fraction)
                / (pb.specific_enthalpy_j_kg - pa.specific_enthalpy_j_kg));
        if ![
            matrix_volume_slope,
            pcm_volume_slope,
            fraction_slope,
            liquid_slope,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(composite_derivative_error());
        }
        Ok(PcmCompositeResponse {
            state,
            temperature_slope_kg_k_j: temperature_slope,
            matrix_enthalpy_derivative: matrix_slope,
            pcm_enthalpy_derivative: pcm_slope,
            pcm_volume_fraction_slope_kg_j: fraction_slope,
            liquid_mass_fraction_slope_kg_j: liquid_slope,
        })
    }

    /// Free volume of the same fixed total mass [m3]; no inventory is rescaled.
    /// The receiving geometry/work owner decides how to use this volume.
    pub fn volume_m3(&self, h: f64, mass_kg: f64) -> Result<f64, PhaseStateError> {
        if !mass_kg.is_finite() || mass_kg <= 0.0 {
            return Err(PhaseStateError::InvalidCurve {
                what: "PCM composite mass must be finite and positive",
            });
        }
        let volume = mass_kg * self.values_at_specific_enthalpy(h)?.specific_volume_m3_kg;
        if !volume.is_finite() || volume <= 0.0 {
            return Err(PhaseStateError::InvalidCurve {
                what: "PCM composite volume must be representable",
            });
        }
        Ok(volume)
    }

    /// Exact mass-weighted original entropy increments [J/(kg composite K)].
    /// Latent entropy is already in the original PCM integral and is counted once.
    pub fn specific_entropy_change_j_kg_k(
        &self,
        initial: f64,
        final_h: f64,
    ) -> Result<f64, PhaseStateError> {
        let a = self.values_at_specific_enthalpy(initial)?;
        let b = self.values_at_specific_enthalpy(final_h)?;
        let matrix_entropy = self.matrix.curve().specific_entropy_change_j_kg_k(
            a.matrix.specific_enthalpy_j_kg,
            b.matrix.specific_enthalpy_j_kg,
        )?;
        let pcm_entropy = self.pcm.specific_entropy_change_j_kg_k(
            a.pcm.specific_enthalpy_j_kg,
            b.pcm.specific_enthalpy_j_kg,
        )?;
        let entropy =
            (1.0 - self.pcm_mass_fraction) * matrix_entropy + self.pcm_mass_fraction * pcm_entropy;
        if !entropy.is_finite() {
            return Err(PhaseStateError::InvalidCurve {
                what: "PCM composite entropy increment must be representable",
            });
        }
        Ok(entropy)
    }

    fn map_enthalpy(&self, h: f64) -> Result<PcmCompositeKnot, PhaseStateError> {
        if !h.is_finite() {
            return Err(PhaseStateError::NonFiniteSpecificEnthalpy);
        }
        let [lower_j_kg, upper_j_kg] = self.enthalpy_range_j_kg();
        if h < lower_j_kg || h > upper_j_kg {
            return Err(PhaseStateError::OutsideEnthalpyDomain {
                specific_enthalpy_j_kg: h,
                lower_j_kg,
                upper_j_kg,
            });
        }
        let upper = self
            .knots
            .partition_point(|knot| knot.specific_enthalpy_j_kg < h);
        let b = self.knots[upper];
        if b.specific_enthalpy_j_kg == h {
            return Ok(b);
        }
        let a = self.knots[upper - 1];
        let alpha =
            (h - a.specific_enthalpy_j_kg) / (b.specific_enthalpy_j_kg - a.specific_enthalpy_j_kg);
        let mapped = PcmCompositeKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: lerp(a.temperature_k, b.temperature_k, alpha),
            matrix_specific_enthalpy_j_kg: lerp(
                a.matrix_specific_enthalpy_j_kg,
                b.matrix_specific_enthalpy_j_kg,
                alpha,
            ),
            pcm_specific_enthalpy_j_kg: lerp(
                a.pcm_specific_enthalpy_j_kg,
                b.pcm_specific_enthalpy_j_kg,
                alpha,
            ),
        };
        if mapped.pcm_specific_enthalpy_j_kg <= a.pcm_specific_enthalpy_j_kg
            || mapped.pcm_specific_enthalpy_j_kg >= b.pcm_specific_enthalpy_j_kg
            || (a.temperature_k < b.temperature_k
                && (mapped.matrix_specific_enthalpy_j_kg <= a.matrix_specific_enthalpy_j_kg
                    || mapped.matrix_specific_enthalpy_j_kg >= b.matrix_specific_enthalpy_j_kg))
        {
            return Err(PhaseStateError::InvalidCurve {
                what: "interior composite enthalpy has no resolvable constituent coordinates",
            });
        }
        Ok(mapped)
    }
}

fn composite_derivative_error() -> PhaseStateError {
    PhaseStateError::InvalidCurve {
        what: "original composite derivative is not representable",
    }
}

/// d(v_constituent)/dH, retaining the original interval's density law.
fn constituent_volume_slope(
    curve: &EquilibriumEnthalpyPhaseCurve,
    segment_start_h: f64,
    density: f64,
    enthalpy_slope: f64,
) -> Result<f64, PhaseStateError> {
    if enthalpy_slope == 0.0 {
        return Ok(0.0);
    }
    let knots = curve.knots();
    let upper = knots
        .partition_point(|k| k.specific_enthalpy_j_kg <= segment_start_h)
        .min(knots.len() - 1);
    let (a, b) = (knots[upper - 1], knots[upper]);
    let dh = b.specific_enthalpy_j_kg - a.specific_enthalpy_j_kg;
    let slope = match curve.density_law_on_interval(upper - 1)? {
        PhaseIntervalDensityLaw::AffineSpecificVolume => {
            ((a.bulk_density_kg_m3 - b.bulk_density_kg_m3)
                / a.bulk_density_kg_m3
                / b.bulk_density_kg_m3
                / dh)
                * enthalpy_slope
        }
        PhaseIntervalDensityLaw::AffineDensity => {
            ((a.bulk_density_kg_m3 - b.bulk_density_kg_m3) / dh / density / density)
                * enthalpy_slope
        }
    };
    if !slope.is_finite() || (a.bulk_density_kg_m3 != b.bulk_density_kg_m3 && slope == 0.0) {
        return Err(composite_derivative_error());
    }
    Ok(slope)
}

impl ThermalLaw for PcmComposite<'_> {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(*self.identity.as_bytes())
    }

    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.temperature_at_specific_enthalpy(h)
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_slope
    }

    fn specific_entropy_change_j_kg_k(
        &self,
        initial: f64,
        final_h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        PcmComposite::specific_entropy_change_j_kg_k(self, initial, final_h)
            .map(Some)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

fn require_same_temperature(expected: f64, actual: f64) -> Result<(), PhaseStateError> {
    if (expected - actual).abs() > 128.0 * f64::EPSILON * expected.max(actual) {
        return Err(PhaseStateError::InvalidCurve {
            what: "composite mapping cannot reproduce the original common temperature",
        });
    }
    Ok(())
}

fn lerp(a: f64, b: f64, alpha: f64) -> f64 {
    if alpha <= 0.5 {
        (b - a).mul_add(alpha, a)
    } else {
        (a - b).mul_add(1.0 - alpha, b)
    }
}

#[cfg(test)]
mod tests;
