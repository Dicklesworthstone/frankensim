//! Continuous authored two-fluid references and calorics on ONE isobar.
//!
//! A positive piecewise-linear Cp(T) supplies h(T) by integration and s(T)
//! by integrating Cp/T. The independent piecewise-linear expansion alpha(T)
//! supplies rho(T) = rho_0 exp(-integral alpha dT). Superfluid fraction,
//! thermodynamic adiabatic sound speed and normal viscosity are interpolated
//! between the same knots; [`TwoFluidState`] derives the remaining response.
//! Consequently neither entropy nor density is an unrelated interpolation.
//!
//! This is a constitutive family for caller-authored, single-composition data,
//! NOT an experimental helium fit. In particular, the reference points in
//! `helium` have different saturation pressures and must not be passed off as
//! one isobar. The caller owns source validity and the external pressure load.
//! Expansion changes the material's equilibrium volume; a thermal owner's fixed
//! heat-exchange geometry is reference geometry, not a rigid vessel constraint.
//! No phase boundary, zero-temperature limit, arbitrary-pressure EOS or update
//! of a channel's kinetic/line storage is inferred by a caloric query.
//!
//! Thermodynamic basis: dh=Cp dT, ds=Cp dT/T and d ln(rho)=-alpha dT at fixed
//! pressure. See Donnelly & Barenghi, J. Phys. Chem. Ref. Data 27, 1217 (1998),
//! sections 7, 8 and 18, DOI: 10.1063/1.556028, for the distinction between
//! saturation and isobaric calorics. Local two-fluid response follows Landau,
//! Phys. Rev. 60, 356 (1941), DOI: 10.1103/PhysRev.60.356, and the finite-expansion
//! sound relations already implemented by [`TwoFluidState`]. No source values
//! or interpolation-error estimate are inferred by this module.

use super::{TwoFluidProperties, TwoFluidState, invalid, nonnegative, positive};
use crate::MaterialError;
use fs_blake3::{ContentHash, DomainHasher};
use fs_math::det;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw, ThermalStorageBasis};

/// Admission/work limit; an immutable law has bounded storage and quadrature.
pub const MAX_SUPERFLUID_CALORIC_KNOTS: usize = 4096;
/// Hard limit for the safeguarded fallback after the analytic inverse estimate.
pub const MAX_SUPERFLUID_CALORIC_INVERSION_STEPS: usize = 64;

/// Absolute pressure and thermodynamic data at the FIRST temperature knot.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SuperfluidCaloricReference {
    /// Externally maintained absolute pressure [Pa], positive and finite.
    pub pressure_pa: f64,
    /// Density at the first knot [kg/m³], positive and finite.
    pub density_at_min_temperature_kg_m3: f64,
    /// Specific-enthalpy datum at the first knot [J/kg], signed and finite.
    pub specific_enthalpy_at_min_temperature_j_kg: f64,
    /// Specific-entropy datum at the first knot [J/(kg K)], positive and finite.
    /// Sound response depends on this physical datum; it is not an arbitrary
    /// additive offset that can be changed independently of the material.
    pub specific_entropy_at_min_temperature_j_kg_k: f64,
}

/// Authored equilibrium data on the declared isobar; all interpolation is linear.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SuperfluidCaloricKnot {
    /// Absolute temperature [K]; knots must be strictly increasing and positive.
    pub temperature_k: f64,
    /// Isobaric heat capacity per total mass [J/(kg K)], strictly positive.
    pub specific_heat_cp_j_kg_k: f64,
    /// Superfluid mass fraction, strictly between zero and one.
    pub superfluid_fraction: f64,
    /// Thermodynamic adiabatic sound speed [m/s], strictly positive. With finite
    /// expansion this is not either coupled two-fluid mode speed.
    pub adiabatic_sound_speed_m_s: f64,
    /// Volumetric isobaric expansion [1/K], signed and finite.
    pub thermal_expansion_per_k: f64,
    /// Normal-component dynamic shear viscosity [Pa s], nonnegative and finite.
    pub normal_viscosity_pa_s: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct PreparedKnot {
    authored: SuperfluidCaloricKnot,
    enthalpy: f64,
    entropy: f64,
    density: f64,
}

/// One checked point on an immutable isobaric caloric law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SuperfluidCaloricState {
    two_fluid: TwoFluidState,
    enthalpy: f64,
}

impl SuperfluidCaloricState {
    /// Original thermodynamically compatible finite-temperature response.
    #[must_use]
    pub const fn two_fluid(&self) -> &TwoFluidState {
        &self.two_fluid
    }
    /// Specific enthalpy in this law's declared datum [J/kg].
    #[must_use]
    pub const fn specific_enthalpy_j_kg(&self) -> f64 {
        self.enthalpy
    }
    /// Absolute temperature [K].
    #[must_use]
    pub fn temperature_k(&self) -> f64 {
        self.two_fluid.temperature_k()
    }
    /// Specific entropy per total mass [J/(kg K)].
    #[must_use]
    pub fn specific_entropy_j_kg_k(&self) -> f64 {
        self.two_fluid.entropy_j_kg_k()
    }
}

/// Immutable integrable isobaric reference family with native thermal ownership.
///
/// Preparation takes O(knots) storage/work; ordinary evaluations locate a segment
/// in O(log knots), allocate nothing, and integrate it analytically. Inversion
/// first tries a scaled quadratic root and then at most 64 bracket reductions.
/// Entropy differences integrate only crossed segments, avoiding subtraction of
/// large absolute entropy datums. All queries refuse domain extrapolation.
#[derive(Clone, Debug)]
pub struct IsobaricSuperfluidLaw {
    reference: SuperfluidCaloricReference,
    knots: Vec<PreparedKnot>,
    identity: ContentHash,
    maximum_temperature_slope: f64,
}

impl IsobaricSuperfluidLaw {
    /// Prepare an authored law. `material_identity` identifies the caller's
    /// material/source context; the returned identity also hashes all data,
    /// pressure, domain and caloric datums. It is not a calibration certificate.
    pub fn try_new(
        material_identity: ContentHash,
        reference: SuperfluidCaloricReference,
        knots: Vec<SuperfluidCaloricKnot>,
    ) -> Result<Self, MaterialError> {
        if !(2..=MAX_SUPERFLUID_CALORIC_KNOTS).contains(&knots.len())
            || !positive(reference.pressure_pa)
            || !positive(reference.density_at_min_temperature_kg_m3)
            || !reference
                .specific_enthalpy_at_min_temperature_j_kg
                .is_finite()
            || !positive(reference.specific_entropy_at_min_temperature_j_kg_k)
        {
            return Err(invalid(
                "isobaric superfluid law needs bounded knots, positive pressure/density/entropy, and finite enthalpy datum",
            ));
        }
        for (i, k) in knots.iter().enumerate() {
            if !positive(k.temperature_k)
                || !positive(k.specific_heat_cp_j_kg_k)
                || !positive(k.adiabatic_sound_speed_m_s)
                || !k.thermal_expansion_per_k.is_finite()
                || !nonnegative(k.normal_viscosity_pa_s)
                || !k.superfluid_fraction.is_finite()
                || !(0.0 < k.superfluid_fraction && k.superfluid_fraction < 1.0)
                || (i > 0 && k.temperature_k <= knots[i - 1].temperature_k)
            {
                return Err(invalid(
                    "isobaric superfluid knots require increasing positive T, Cp and sound speed, strict component fractions and finite expansion/viscosity",
                ));
            }
        }
        let min_cp = knots
            .iter()
            .map(|k| k.specific_heat_cp_j_kg_k)
            .fold(f64::INFINITY, f64::min);
        // Round outward: the thermal timeline requires an upper bound, not an
        // accidentally downward-rounded reciprocal at the capacity minimum.
        let maximum_temperature_slope = (1.0 / min_cp).next_up();
        if !positive(maximum_temperature_slope) {
            return Err(invalid(
                "isobaric superfluid global temperature slope is unrepresentable",
            ));
        }
        let identity = law_identity(material_identity, reference, &knots);
        let mut prepared = Vec::with_capacity(knots.len());
        prepared.push(PreparedKnot {
            authored: knots[0],
            enthalpy: reference.specific_enthalpy_at_min_temperature_j_kg,
            entropy: reference.specific_entropy_at_min_temperature_j_kg_k,
            density: reference.density_at_min_temperature_kg_m3,
        });
        for pair in knots.windows(2) {
            let previous = *prepared.last().expect("first knot exists");
            let width = pair[1].temperature_k - pair[0].temperature_k;
            let delta_h = width
                * (0.5 * pair[0].specific_heat_cp_j_kg_k + 0.5 * pair[1].specific_heat_cp_j_kg_k);
            let delta_s = segment_entropy(
                pair[0],
                pair[1],
                pair[0].temperature_k,
                pair[1].temperature_k,
            )?;
            let expansion = width
                * (0.5 * pair[0].thermal_expansion_per_k + 0.5 * pair[1].thermal_expansion_per_k);
            let density = expanded_density(previous.density, expansion)?;
            let enthalpy = previous.enthalpy + delta_h;
            let entropy = previous.entropy + delta_s;
            if !positive(delta_h)
                || !positive(delta_s)
                || !enthalpy.is_finite()
                || enthalpy <= previous.enthalpy
                || !positive(entropy)
                || entropy <= previous.entropy
            {
                return Err(invalid(
                    "isobaric superfluid caloric interval or datum is unrepresentable",
                ));
            }
            // Linear alpha can change sign inside an interval. Check the actual
            // stationary density as well as the endpoints, not just knot values.
            let alpha_a = pair[0].thermal_expansion_per_k;
            let alpha_b = pair[1].thermal_expansion_per_k;
            if (alpha_a < 0.0 && alpha_b > 0.0) || (alpha_a > 0.0 && alpha_b < 0.0) {
                let scale = alpha_a.abs().max(alpha_b.abs());
                let weight =
                    (alpha_a.abs() / scale) / (alpha_a.abs() / scale + alpha_b.abs() / scale);
                expanded_density(previous.density, (width * weight) * (0.5 * alpha_a))?;
            }
            prepared.push(PreparedKnot {
                authored: pair[1],
                enthalpy,
                entropy,
                density,
            });
        }
        let law = Self {
            reference,
            knots: prepared,
            identity,
            maximum_temperature_slope,
        };
        for knot in &law.knots {
            law.state_at_temperature_k(knot.authored.temperature_k)?;
        }
        Ok(law)
    }

    /// Pressure of every admitted state [Pa]. There is no other-pressure query.
    #[must_use]
    pub const fn pressure_pa(&self) -> f64 {
        self.reference.pressure_pa
    }
    /// Original reference data and physical entropy datum.
    #[must_use]
    pub const fn reference(&self) -> SuperfluidCaloricReference {
        self.reference
    }
    /// Identity including all constitutive inputs, domain, pressure and datums.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Closed admitted absolute-temperature interval [K].
    #[must_use]
    pub fn temperature_bounds_k(&self) -> [f64; 2] {
        [
            self.knots[0].authored.temperature_k,
            self.knots[self.knots.len() - 1].authored.temperature_k,
        ]
    }
    /// Closed admitted specific-enthalpy interval [J/kg].
    #[must_use]
    pub fn specific_enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        [
            self.knots[0].enthalpy,
            self.knots[self.knots.len() - 1].enthalpy,
        ]
    }
    /// Certified global upper bound on dT/dh over the authored law [kg K/J].
    #[must_use]
    pub const fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_temperature_slope
    }
    /// Number of authored temperature knots.
    #[must_use]
    pub fn knot_count(&self) -> usize {
        self.knots.len()
    }

    fn temperature_segment(&self, temperature_k: f64) -> Result<usize, MaterialError> {
        let [lower, upper] = self.temperature_bounds_k();
        if !temperature_k.is_finite() || temperature_k < lower || temperature_k > upper {
            return Err(invalid(
                "temperature is outside the isobaric superfluid law",
            ));
        }
        Ok(self
            .knots
            .partition_point(|k| k.authored.temperature_k <= temperature_k)
            .saturating_sub(1)
            .min(self.knots.len() - 2))
    }

    /// Specific enthalpy from the SAME integrated Cp used by thermal inversion.
    pub fn specific_enthalpy_at_temperature_k(
        &self,
        temperature_k: f64,
    ) -> Result<f64, MaterialError> {
        let i = self.temperature_segment(temperature_k)?;
        self.segment_enthalpy(i, temperature_k)
    }

    fn segment_enthalpy(&self, i: usize, temperature_k: f64) -> Result<f64, MaterialError> {
        let a = self.knots[i];
        let b = self.knots[i + 1];
        if temperature_k == a.authored.temperature_k {
            return Ok(a.enthalpy);
        }
        if temperature_k == b.authored.temperature_k {
            return Ok(b.enthalpy);
        }
        let x = temperature_k - a.authored.temperature_k;
        let w = x / (b.authored.temperature_k - a.authored.temperature_k);
        let cp = interpolate(
            a.authored.specific_heat_cp_j_kg_k,
            b.authored.specific_heat_cp_j_kg_k,
            w,
        );
        let h = a.enthalpy + x * (0.5 * a.authored.specific_heat_cp_j_kg_k + 0.5 * cp);
        if !h.is_finite() {
            return Err(invalid("isobaric superfluid enthalpy is unrepresentable"));
        }
        Ok(h)
    }

    /// Recover temperature without extrapolation. The analytic root is scaled
    /// before squaring; a bracketed fallback handles cancellation/roundoff.
    /// Bounded-work exhaustion or unrepresentable arithmetic is a refusal.
    pub fn temperature_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
    ) -> Result<f64, MaterialError> {
        let [minimum, maximum] = self.specific_enthalpy_bounds_j_kg();
        if !specific_h_j_kg.is_finite() || specific_h_j_kg < minimum || specific_h_j_kg > maximum {
            return Err(invalid("enthalpy is outside the isobaric superfluid law"));
        }
        let i = self
            .knots
            .partition_point(|k| k.enthalpy <= specific_h_j_kg)
            .saturating_sub(1)
            .min(self.knots.len() - 2);
        let a = self.knots[i];
        let b = self.knots[i + 1];
        if specific_h_j_kg == a.enthalpy {
            return Ok(a.authored.temperature_k);
        }
        if specific_h_j_kg == b.enthalpy {
            return Ok(b.authored.temperature_k);
        }
        let mut lower = a.authored.temperature_k;
        let mut upper = b.authored.temperature_k;
        let width = upper - lower;
        let cp_a = a.authored.specific_heat_cp_j_kg_k;
        let cp_b = b.authored.specific_heat_cp_j_kg_k;
        let scale = cp_a.max(cp_b);
        let normalized_cp = cp_a / scale;
        let q = ((specific_h_j_kg - a.enthalpy) / width) / scale;
        let discriminant = normalized_cp * normalized_cp + 2.0 * ((cp_b - cp_a) / scale) * q;
        let denominator = normalized_cp + det::sqrt(discriminant.max(0.0));
        let weight = 2.0 * q / denominator;
        let estimate = lower + width * weight;
        let mut candidate = if estimate.is_finite() && lower <= estimate && estimate <= upper {
            estimate
        } else {
            lower + 0.5 * (upper - lower)
        };
        for _ in 0..MAX_SUPERFLUID_CALORIC_INVERSION_STEPS {
            let h = self.segment_enthalpy(i, candidate)?;
            if h == specific_h_j_kg
                || (h - specific_h_j_kg).abs() * self.maximum_temperature_slope
                    <= 8.0 * f64::EPSILON * candidate.abs()
            {
                return Ok(candidate);
            }
            if h < specific_h_j_kg {
                lower = candidate;
            } else {
                upper = candidate;
            }
            let midpoint = lower + 0.5 * (upper - lower);
            if midpoint == lower || midpoint == upper {
                let h_lower = self.segment_enthalpy(i, lower)?;
                let h_upper = self.segment_enthalpy(i, upper)?;
                return Ok(
                    if (specific_h_j_kg - h_lower).abs() <= (h_upper - specific_h_j_kg).abs() {
                        lower
                    } else {
                        upper
                    },
                );
            }
            candidate = midpoint;
        }
        Err(invalid(
            "isobaric superfluid inverse exhausted its bounded bracket reductions",
        ))
    }

    /// Continuous temperature-dependent two-fluid properties and caloric storage.
    /// The local state remains subject to TwoFluidState's representability checks.
    pub fn state_at_temperature_k(
        &self,
        temperature_k: f64,
    ) -> Result<SuperfluidCaloricState, MaterialError> {
        let i = self.temperature_segment(temperature_k)?;
        let a = self.knots[i];
        let b = self.knots[i + 1];
        let x = temperature_k - a.authored.temperature_k;
        let weight = x / (b.authored.temperature_k - a.authored.temperature_k);
        let alpha = interpolate(
            a.authored.thermal_expansion_per_k,
            b.authored.thermal_expansion_per_k,
            weight,
        );
        let (density, entropy) = if temperature_k == b.authored.temperature_k {
            (b.density, b.entropy)
        } else if temperature_k == a.authored.temperature_k {
            (a.density, a.entropy)
        } else {
            (
                expanded_density(
                    a.density,
                    x * (0.5 * a.authored.thermal_expansion_per_k + 0.5 * alpha),
                )?,
                a.entropy
                    + segment_entropy(
                        a.authored,
                        b.authored,
                        a.authored.temperature_k,
                        temperature_k,
                    )?,
            )
        };
        let two_fluid = TwoFluidState::try_new(TwoFluidProperties {
            temperature_k,
            density_kg_m3: density,
            superfluid_fraction: interpolate(
                a.authored.superfluid_fraction,
                b.authored.superfluid_fraction,
                weight,
            ),
            entropy_j_kg_k: entropy,
            specific_heat_cp_j_kg_k: interpolate(
                a.authored.specific_heat_cp_j_kg_k,
                b.authored.specific_heat_cp_j_kg_k,
                weight,
            ),
            adiabatic_sound_speed_m_s: interpolate(
                a.authored.adiabatic_sound_speed_m_s,
                b.authored.adiabatic_sound_speed_m_s,
                weight,
            ),
            thermal_expansion_per_k: alpha,
            normal_viscosity_pa_s: interpolate(
                a.authored.normal_viscosity_pa_s,
                b.authored.normal_viscosity_pa_s,
                weight,
            ),
        })?;
        Ok(SuperfluidCaloricState {
            two_fluid,
            enthalpy: self.segment_enthalpy(i, temperature_k)?,
        })
    }

    /// Recover both thermodynamic response and temperature from stored enthalpy.
    /// The supplied storage is retained exactly; its inversion has floating-point
    /// resolution and does not silently change the thermal owner's joules.
    pub fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
    ) -> Result<SuperfluidCaloricState, MaterialError> {
        let temperature = self.temperature_at_specific_enthalpy(specific_h_j_kg)?;
        let mut state = self.state_at_temperature_k(temperature)?;
        state.enthalpy = specific_h_j_kg;
        Ok(state)
    }

    /// Integral of dh/T between two admitted storage values [J/(kg K)].
    /// Reversing the arguments reverses the sign; datums do not enter subtraction.
    pub fn specific_entropy_change_j_kg_k(
        &self,
        initial_h: f64,
        final_h: f64,
    ) -> Result<f64, MaterialError> {
        self.temperature_at_specific_enthalpy(initial_h)?;
        self.temperature_at_specific_enthalpy(final_h)?;
        if initial_h == final_h {
            return Ok(0.0);
        }
        let (lower, upper, sign) = if initial_h < final_h {
            (initial_h, final_h, 1.0)
        } else {
            (final_h, initial_h, -1.0)
        };
        let segment = |h| {
            self.knots
                .partition_point(|k| k.enthalpy <= h)
                .saturating_sub(1)
                .min(self.knots.len() - 2)
        };
        let first = segment(lower);
        let last = segment(upper);
        let mut total = 0.0;
        let mut compensation = 0.0;
        for i in first..=last {
            let a = self.knots[i];
            let b = self.knots[i + 1];
            let start = lower.max(a.enthalpy);
            let end = upper.min(b.enthalpy);
            if start < end {
                let t_start = self.temperature_at_specific_enthalpy(start)?;
                let t_end = self.temperature_at_specific_enthalpy(end)?;
                let entropy = if t_start == t_end {
                    // Distinct h may resolve more finely than absolute T.
                    // Preserve the actual exchanged storage even in this limit.
                    (end - start) / t_start
                } else {
                    let width = b.authored.temperature_k - a.authored.temperature_k;
                    let cp_start = interpolate(
                        a.authored.specific_heat_cp_j_kg_k,
                        b.authored.specific_heat_cp_j_kg_k,
                        (t_start - a.authored.temperature_k) / width,
                    );
                    let cp_end = interpolate(
                        a.authored.specific_heat_cp_j_kg_k,
                        b.authored.specific_heat_cp_j_kg_k,
                        (t_end - a.authored.temperature_k) / width,
                    );
                    let dh = (t_end - t_start) * (0.5 * cp_start + 0.5 * cp_end);
                    let ds = segment_entropy(a.authored, b.authored, t_start, t_end)?;
                    // ds/dh is the exact mean reciprocal temperature for a
                    // linear Cp. Rescale by ORIGINAL h, retaining tiny heat
                    // exchanges when inverse temperatures have rounded.
                    (ds / dh) * (end - start)
                };
                if !positive(entropy) {
                    return Err(invalid(
                        "isobaric superfluid entropy exchange is unrepresentable",
                    ));
                }
                let next = entropy - compensation;
                let accumulated = total + next;
                compensation = (accumulated - total) - next;
                total = accumulated;
            }
        }
        if !nonnegative(total) {
            return Err(invalid(
                "isobaric superfluid entropy difference is unrepresentable",
            ));
        }
        Ok(sign * total)
    }
}

impl ThermalLaw for IsobaricSuperfluidLaw {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::Enthalpy
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(*self.identity.as_bytes())
    }
    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        self.temperature_at_specific_enthalpy(specific_h_j_kg)
            .map_err(|_| DepositionError::MaterialDomain)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        IsobaricSuperfluidLaw::maximum_temperature_slope(self)
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h: f64,
        final_h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        IsobaricSuperfluidLaw::specific_entropy_change_j_kg_k(self, initial_h, final_h)
            .map(Some)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

fn interpolate(a: f64, b: f64, weight: f64) -> f64 {
    if weight == 0.0 {
        a
    } else if weight == 1.0 {
        b
    } else {
        let value = if a.is_sign_positive() == b.is_sign_positive() {
            a + (b - a) * weight
        } else {
            (1.0 - weight) * a + weight * b
        };
        value.clamp(a.min(b), a.max(b))
    }
}

fn expanded_density(reference: f64, integrated_expansion: f64) -> Result<f64, MaterialError> {
    let density = reference * det::exp(-integrated_expansion);
    if !integrated_expansion.is_finite() || !positive(density) {
        return Err(invalid(
            "integrated isobaric superfluid density is unrepresentable",
        ));
    }
    Ok(density)
}

// Positive endpoint weights for the exact integral of a linear Cp over [a,b].
// The series avoids losing both ln(1+u) and 1-ln(1+u)/u for tiny increments.
fn segment_entropy(
    a: SuperfluidCaloricKnot,
    b: SuperfluidCaloricKnot,
    start: f64,
    end: f64,
) -> Result<f64, MaterialError> {
    if start == end {
        return Ok(0.0);
    }
    let width = b.temperature_k - a.temperature_k;
    let cp_start = interpolate(
        a.specific_heat_cp_j_kg_k,
        b.specific_heat_cp_j_kg_k,
        (start - a.temperature_k) / width,
    );
    let cp_end = interpolate(
        a.specific_heat_cp_j_kg_k,
        b.specific_heat_cp_j_kg_k,
        (end - a.temperature_k) / width,
    );
    let u = (end - start) / start;
    if !positive(u) {
        return Err(invalid(
            "isobaric superfluid entropy interval is unrepresentable",
        ));
    }
    let (logarithm, upper_weight) = if u < 0.001 {
        // g(u)=1-ln(1+u)/u; retaining terms through u^7 gives <1e-24
        // absolute truncation error on this branch, before floating rounding.
        let g = u
            * (0.5
                + u * (-1.0 / 3.0
                    + u * (0.25 + u * (-0.2 + u * (1.0 / 6.0 + u * (-1.0 / 7.0 + u * 0.125))))));
        (u * (1.0 - g), g)
    } else {
        let logarithm = det::ln(1.0 + u);
        (logarithm, 1.0 - logarithm / u)
    };
    let lower_weight = logarithm - upper_weight;
    let entropy = cp_start * lower_weight + cp_end * upper_weight;
    if !nonnegative(lower_weight) || !nonnegative(upper_weight) || !positive(entropy) {
        return Err(invalid(
            "isobaric superfluid entropy quadrature is unrepresentable",
        ));
    }
    Ok(entropy)
}

fn law_identity(
    material_identity: ContentHash,
    reference: SuperfluidCaloricReference,
    knots: &[SuperfluidCaloricKnot],
) -> ContentHash {
    let mut hasher = DomainHasher::new("org.frankensim.fs-material.isobaric-superfluid-caloric.v1");
    hasher.update(material_identity.as_bytes());
    for value in [
        reference.pressure_pa,
        reference.density_at_min_temperature_kg_m3,
        reference.specific_enthalpy_at_min_temperature_j_kg,
        reference.specific_entropy_at_min_temperature_j_kg_k,
    ] {
        hasher.update(&value.to_bits().to_le_bytes());
    }
    hasher.update(&(knots.len() as u64).to_le_bytes());
    for knot in knots {
        for value in [
            knot.temperature_k,
            knot.specific_heat_cp_j_kg_k,
            knot.superfluid_fraction,
            knot.adiabatic_sound_speed_m_s,
            knot.thermal_expansion_per_k,
            knot.normal_viscosity_pa_s,
        ] {
            hasher.update(&value.to_bits().to_le_bytes());
        }
    }
    hasher.finalize()
}

#[cfg(test)]
#[path = "caloric_tests.rs"]
mod tests;
