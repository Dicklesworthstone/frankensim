//! Finite J2 plasticity with affine thermal coefficients and thermal stretch.
//!
//! All free energies refer to the ORIGINAL volume. The coefficients are
//! Kirchhoff/log-strain moduli and yield radii: at a freely expanded state,
//! current-volume Cauchy moduli differ by the thermal volume factor. This is
//! an explicitly authored constitutive model, not a conversion of ordinary
//! small-strain or current-volume material claims.
//!
//! `Fm=F/theta(T)`, `ln(theta)=integral(a(T)dT)`, where the explicitly
//! logarithmic expansion coefficient is `a=d ln(theta)/dT`. Plastic history
//! remains exactly the native `Cp^-1, alpha`. The existing exponential J2
//! return is evaluated at the endpoint temperature. At frozen temperature,
//! `P=P_native/theta` and `dP/dF=A_native/theta²` by the exact chain rule.
//!
//! For `g=ln(det(Fm))` and `e=dev(Ee)`, mechanical free energy is
//! `psi=mu(T)e:e+K(T)g²/2+H(T)alpha²/2`. Entropy is its negative temperature
//! derivative at FIXED total F/history; internal energy is `psi+T*s`.
//! Caloric storage stays with the supplied thermal owner. No latent heat,
//! pressure-flow conversion, configurational plastic entropy, damage, or
//! independent Taylor–Quinney coefficient is inferred here.
//!
//! A finite step retains physical yield dissipation, reversible caloric
//! exchange and a SIGNED mechanical endpoint-work defect separately. The
//! defect is never heated. This endpoint rule needs timestep refinement;
//! first-law closure alone does not establish temporal accuracy.
//!
//! Explicit Perzyna kinetics may be selected with [`FiniteJ2ThermalOverstress`].
//! Their Arrhenius viscosity changes flow and physical dissipation, not free
//! energy or energetic hardening. Every constitutive return then requires the
//! actual physical duration. The reduced kinetic potential and physical heat
//! are distinct: only the latter enters the thermal work account.
//!
//! [`AffineThermalFiniteJ2::new_norton`] explicitly selects threshold-free
//! Norton creep with the same finite kinematics, stress and heat owners. Its
//! rate is `alpha_dot=(s0/eta(T))*(q/s0)^n` for equivalent KIRCHHOFF stress q,
//! not an inferred conversion of an ordinary Cauchy-stress creep fit. This
//! selection has no yield or hardening energy. The native `Cp^-1, alpha` state,
//! plastic-named result fields and history bounds then describe creep: alpha
//! is accumulated equivalent creep, and physical loss is `q*delta_alpha`.
//! This is one secondary-creep mechanism, not simultaneous creep and yielding,
//! primary/tertiary creep, rupture or an empirical material calibration.
//! The [Norton/Arrhenius reference](https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_theory.06.039.html)
//! uses the rate prefactor A; the explicit normalization here is `A=s0/eta`.
//!
//! The explicit exponential/tangent expansion convention is described in
//! the [COMSOL finite-strain theory reference](https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_theory.06.031.html).
//! This reduced free-energy model is not the full independent-plastic-entropy
//! theory of Simo–Miehe, doi:10.1016/0045-7825(92)90170-O.

use super::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use super::{
    FiniteJ2, FiniteJ2Error, FiniteJ2Limits, FiniteJ2Parameters, FiniteJ2Response, FiniteJ2State,
};
use fs_math::det;

/// Exact model-card binding for the affine Kirchhoff/logarithmic thermal family.
pub mod source;

/// Explicit derivatives of the original-reference-volume constitutive coefficients.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteJ2TemperatureSlopes {
    pub shear_modulus_pa_k: f64,
    pub bulk_modulus_pa_k: f64,
    pub yield_stress_pa_k: f64,
    pub hardening_modulus_pa_k: f64,
}

/// One smooth, globally affine family over a declared positive-temperature span.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteJ2ThermalParameters {
    /// Stress-free thermal reference; theta=1 here.
    pub reference_temperature_k: f64,
    pub temperature_range_k: [f64; 2],
    pub reference: FiniteJ2Parameters,
    pub slopes: FiniteJ2TemperatureSlopes,
    /// Explicit logarithmic instantaneous expansion `d ln(theta)/dT` at T_ref.
    pub logarithmic_expansion_per_k: f64,
    /// Derivative of the logarithmic expansion coefficient; finite and constant.
    pub logarithmic_expansion_slope_per_k2: f64,
}

/// Affine elastic and thermal coefficients for explicit secondary Norton creep.
/// Yield and hardening parameters are absent because this mechanism has neither.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteNortonThermalParameters {
    pub reference_temperature_k: f64,
    pub temperature_range_k: [f64; 2],
    /// Original-reference-volume Hencky shear modulus at the reference temperature.
    pub shear_modulus_pa: f64,
    /// Original-reference-volume Hencky bulk modulus at the reference temperature.
    pub bulk_modulus_pa: f64,
    pub shear_modulus_slope_pa_k: f64,
    pub bulk_modulus_slope_pa_k: f64,
    /// Logarithmic instantaneous thermal expansion, `d ln(theta)/dT`.
    pub logarithmic_expansion_per_k: f64,
    pub logarithmic_expansion_slope_per_k2: f64,
}

/// Explicit power-law viscosity and thermal activation, in the native Kirchhoff
/// stress/original-volume convention. This is a dissipative law: its parameters
/// contribute no artificial stored energy or mechanical entropy. The driving
/// stress is Perzyna overstress for plasticity and full equivalent stress for
/// the explicitly selected Norton creep family.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteJ2ThermalOverstress {
    /// Positive viscosity at the family's stress-free reference temperature.
    pub reference_viscosity_pa_s: f64,
    /// Positive stress scale in `alpha_dot=(sigma_ref/eta)*(f/sigma_ref)^n`.
    pub stress_scale_pa: f64,
    /// Dimensionless Perzyna exponent, finite and at least one.
    pub exponent: f64,
    /// Explicit `Q/R` in kelvin, finite and nonnegative. Viscosity is
    /// `eta_ref*exp((Q/R)*(1/T-1/T_ref))`; zero disables thermal activation.
    pub activation_temperature_k: f64,
}

/// Immutable native family; source admission is separate from this constitutive law.
#[derive(Clone, Debug, PartialEq)]
pub struct AffineThermalFiniteJ2 {
    parameters: FiniteJ2ThermalParameters,
    limits: FiniteJ2Limits,
    minimum_mechanical_heat_capacity_j_m3_k: f64,
    overstress: Option<FiniteJ2ThermalOverstress>,
    norton_creep: bool,
}

/// Keep only one native backbone per frozen point, including kinetic points.
#[derive(Clone, Debug, PartialEq)]
enum FrozenNative {
    RateIndependent(FiniteJ2),
    /// Fixed-history thermodynamics are defined without a physical interval.
    /// A constitutive return from this variant always refuses.
    Inspection(FiniteJ2),
    Overstress(FiniteOverstressJ2),
}

/// One owned endpoint-temperature law, without a borrowed temporary source field.
#[derive(Clone, Debug, PartialEq)]
pub struct FrozenThermalFiniteJ2 {
    native: FrozenNative,
    slopes: FiniteJ2TemperatureSlopes,
    temperature_k: f64,
    log_thermal_stretch: f64,
    thermal_stretch: f64,
    inverse_thermal_stretch: f64,
    expansion_per_k: f64,
    expansion_slope_per_k2: f64,
}

/// Borrow the complete thermal family, preserving its total-F mapping.
pub trait FiniteThermalJ2Material {
    fn thermal_finite_j2(&self) -> &AffineThermalFiniteJ2;
}
impl FiniteThermalJ2Material for AffineThermalFiniteJ2 {
    fn thermal_finite_j2(&self) -> &AffineThermalFiniteJ2 {
        self
    }
}
impl<T: FiniteThermalJ2Material + ?Sized> FiniteThermalJ2Material for &T {
    fn thermal_finite_j2(&self) -> &AffineThermalFiniteJ2 {
        T::thermal_finite_j2(*self)
    }
}

/// Thermodynamic inspection at fixed history: this does not perform a return.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteJ2Thermodynamics {
    pub first_piola_pa: [f64; 9],
    /// Mixed derivative at fixed total F and Cp/alpha, not along a plastic return.
    pub first_piola_temperature_derivative_fixed_history_pa_k: [f64; 9],
    pub free_energy_j_m3: f64,
    pub entropy_j_m3_k: f64,
    pub internal_energy_j_m3: f64,
    /// Mechanical contribution at fixed total F/history; it can be negative.
    pub heat_capacity_j_m3_k: f64,
    /// Actual total volume ratio det(F).
    pub jacobian: f64,
    /// Elastic volume ratio det(Fm); plastic flow is isochoric.
    pub mechanical_jacobian: f64,
}

/// Complete native candidate and independently inspectable thermal work split.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteJ2ThermalIncrement {
    /// Native history and reduced potential, with stress/tangent mapped to total F.
    pub mechanical: FiniteJ2Response,
    pub old: FiniteJ2Thermodynamics,
    pub new: FiniteJ2Thermodynamics,
    /// Change of coefficients at old mechanical deformation and old history.
    pub coefficient_storage_change_j_m3: f64,
    /// Delta(T*s_mechanical), separate from free-energy storage.
    pub thermoelastic_storage_change_j_m3: f64,
    /// Actual endpoint total-Piola work P_1:(F_1-F_0).
    pub mechanical_work_j_m3: f64,
    /// Endpoint native-Piola work Pm_1:(Fm_1-Fm_0).
    pub mechanical_coordinate_work_j_m3: f64,
    /// Difference of total and mechanical-coordinate endpoint work.
    pub thermal_deformation_work_j_m3: f64,
    /// Signed mechanical-coordinate endpoint quadrature defect, never heat.
    pub endpoint_work_defect_j_m3: f64,
    /// Signed transfer into the separate caloric internal-energy chart.
    pub reversible_heat_j_m3: f64,
    /// Total work minus mechanical internal storage, caloric heat and defect.
    /// The spatial/coupling owner applies its declared dimensional work gate.
    pub balance_error_j_m3: f64,
}
impl FiniteJ2ThermalIncrement {
    /// Deposit once: physical inelastic heat plus signed reversible exchange.
    pub fn caloric_heat_j_m3(&self) -> f64 {
        self.mechanical.plastic_dissipation_j_m3 + self.reversible_heat_j_m3
    }
    /// Stable two-temperature mechanical internal-energy increment.
    pub fn internal_energy_change_j_m3(&self) -> f64 {
        self.new.free_energy_j_m3 - self.old.free_energy_j_m3
            + self.thermoelastic_storage_change_j_m3
    }
}

/// Refusals preserve all borrowed deformation, temperature and plastic history.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteJ2ThermalError {
    Material(FiniteJ2Error),
    Parameters(&'static str),
    Temperature {
        temperature_k: f64,
        minimum_k: f64,
        maximum_k: f64,
    },
    Numerical(&'static str),
    /// Selected physical kinetics cannot be advanced without an interval.
    DurationRequired,
}
impl From<FiniteJ2Error> for FiniteJ2ThermalError {
    fn from(error: FiniteJ2Error) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for FiniteJ2ThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => error.fmt(f),
            Self::Parameters(what) => write!(f, "finite thermoplastic parameters: {what}"),
            Self::Temperature {
                temperature_k,
                minimum_k,
                maximum_k,
            } => write!(
                f,
                "finite thermoplastic temperature {temperature_k} K outside [{minimum_k},{maximum_k}] K"
            ),
            Self::Numerical(what) => write!(f, "finite thermoplastic numerical refusal: {what}"),
            Self::DurationRequired => f.write_str(
                "finite thermoplastic overstress requires an explicit physical duration",
            ),
        }
    }
}
impl std::error::Error for FiniteJ2ThermalError {}

fn finite(value: f64) -> Result<f64, FiniteJ2ThermalError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FiniteJ2ThermalError::Numerical(
            "nonfinite thermal coefficient, stretch or account",
        ))
    }
}

impl AffineThermalFiniteJ2 {
    /// All native limits apply to MECHANICAL Fm, including minimum Jacobian.
    /// The spatial owner separately retains its total-F geometric Jacobian gate.
    pub fn new(
        parameters: FiniteJ2ThermalParameters,
        limits: FiniteJ2Limits,
    ) -> Result<Self, FiniteJ2ThermalError> {
        Self::build(parameters, limits, false)
    }

    /// Select secondary Norton creep with explicit positive kinetics and no
    /// yield threshold or hardening storage. Every return requires physical dt.
    /// The existing native metric/alpha state and plastic-named limits then
    /// describe creep; their full reachability and finite-strain checks remain.
    pub fn new_norton(
        parameters: FiniteNortonThermalParameters,
        limits: FiniteJ2Limits,
        kinetics: FiniteJ2ThermalOverstress,
    ) -> Result<Self, FiniteJ2ThermalError> {
        let p = parameters;
        let shared = FiniteJ2ThermalParameters {
            reference_temperature_k: p.reference_temperature_k,
            temperature_range_k: p.temperature_range_k,
            reference: FiniteJ2Parameters {
                shear_modulus_pa: p.shear_modulus_pa,
                bulk_modulus_pa: p.bulk_modulus_pa,
                yield_stress_pa: 0.0,
                hardening_modulus_pa: 0.0,
            },
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: p.shear_modulus_slope_pa_k,
                bulk_modulus_pa_k: p.bulk_modulus_slope_pa_k,
                yield_stress_pa_k: 0.0,
                hardening_modulus_pa_k: 0.0,
            },
            logarithmic_expansion_per_k: p.logarithmic_expansion_per_k,
            logarithmic_expansion_slope_per_k2: p.logarithmic_expansion_slope_per_k2,
        };
        Self::build(shared, limits, true)?.with_overstress(kinetics)
    }

    fn build(
        parameters: FiniteJ2ThermalParameters,
        limits: FiniteJ2Limits,
        norton_creep: bool,
    ) -> Result<Self, FiniteJ2ThermalError> {
        let p = &parameters;
        let [lo, hi] = p.temperature_range_k;
        if !lo.is_finite()
            || !hi.is_finite()
            || lo <= 0.0
            || hi <= lo
            || !p.reference_temperature_k.is_finite()
            || !(lo..=hi).contains(&p.reference_temperature_k)
            || [
                p.slopes.shear_modulus_pa_k,
                p.slopes.bulk_modulus_pa_k,
                p.slopes.yield_stress_pa_k,
                p.slopes.hardening_modulus_pa_k,
                p.logarithmic_expansion_per_k,
                p.logarithmic_expansion_slope_per_k2,
            ]
            .iter()
            .any(|x| !x.is_finite())
        {
            return Err(FiniteJ2ThermalError::Parameters(
                "finite affine derivatives and a positive finite T interval containing T_ref required",
            ));
        }
        let mut family = Self {
            parameters,
            limits,
            minimum_mechanical_heat_capacity_j_m3_k: 0.0,
            overstress: None,
            norton_creep,
        };
        family.native_backbone(p.reference)?;
        let low = family.at_temperature(lo)?;
        let high = family.at_temperature(hi)?;
        // The quadratic log stretch has at most one interior extremum.
        if p.logarithmic_expansion_slope_per_k2 != 0.0 {
            let vertex = p.reference_temperature_k
                - p.logarithmic_expansion_per_k / p.logarithmic_expansion_slope_per_k2;
            if (lo..=hi).contains(&vertex) {
                family.at_temperature(vertex)?;
            }
        }
        let max_bulk = low
            .backbone()
            .parameters()
            .bulk_modulus_pa
            .max(high.backbone().parameters().bulk_modulus_pa);
        let max_expansion = low.expansion_per_k.abs().max(high.expansion_per_k.abs());
        let max_log_volume = det::sqrt(3.0) * limits.maximum_elastic_log_strain;
        // mu''=H''=0. Bound the remaining volumetric/thermal contribution
        // throughout the full admitted elastic and temperature envelope.
        let capacity_magnitude = finite(
            hi * ((6.0 * p.slopes.bulk_modulus_pa_k.abs() * max_expansion
                + 3.0 * max_bulk * p.logarithmic_expansion_slope_per_k2.abs())
                * max_log_volume
                + 9.0 * max_bulk * max_expansion * max_expansion),
        )?;
        family.minimum_mechanical_heat_capacity_j_m3_k = -capacity_magnitude;
        Ok(family)
    }

    fn native_backbone(&self, parameters: FiniteJ2Parameters) -> Result<FiniteJ2, FiniteJ2Error> {
        if self.norton_creep {
            FiniteJ2::new_creep_backbone(parameters, self.limits)
        } else {
            FiniteJ2::new(parameters, self.limits)
        }
    }

    /// Whether Cp/alpha and native inelastic result fields describe the
    /// explicitly selected threshold-free Norton mechanism.
    pub const fn is_norton_creep(&self) -> bool {
        self.norton_creep
    }

    pub fn parameters(&self) -> &FiniteJ2ThermalParameters {
        &self.parameters
    }
    pub fn limits(&self) -> &FiniteJ2Limits {
        &self.limits
    }
    pub fn initial_state(&self) -> FiniteJ2State {
        FiniteJ2State {
            inverse_plastic_cauchy_green: super::IDENTITY6,
            alpha: 0.0,
        }
    }
    /// Add this lower bound to the ORIGINAL-volume caloric capacity. The
    /// source/coupling owner must require a strictly positive total capacity.
    pub fn minimum_mechanical_heat_capacity_j_m3_k(&self) -> f64 {
        self.minimum_mechanical_heat_capacity_j_m3_k
    }

    /// Select explicit power-law/Arrhenius kinetics without changing any elastic,
    /// energetic-hardening or thermal-expansion coefficient. The complete
    /// positive-temperature span must admit representable positive viscosity.
    pub fn with_overstress(
        mut self,
        parameters: FiniteJ2ThermalOverstress,
    ) -> Result<Self, FiniteJ2ThermalError> {
        FiniteOverstressParameters {
            viscosity_pa_s: parameters.reference_viscosity_pa_s,
            stress_scale_pa: parameters.stress_scale_pa,
            exponent: parameters.exponent,
        }
        .validate()?;
        if !parameters.activation_temperature_k.is_finite()
            || parameters.activation_temperature_k < 0.0
        {
            return Err(FiniteJ2ThermalError::Parameters(
                "finite thermal overstress requires a finite nonnegative activation temperature",
            ));
        }
        self.overstress = Some(parameters);
        // Arrhenius viscosity is monotone for nonnegative activation, so its
        // two endpoint values bound the entire admitted thermal interval.
        for temperature in self.parameters.temperature_range_k {
            self.overstress_viscosity_pa_s(temperature)?;
        }
        Ok(self)
    }

    /// The original explicitly selected kinetics; None is rate independence.
    /// Norton creep always has kinetics; its driving stress is full equivalent q.
    pub const fn overstress(&self) -> Option<&FiniteJ2ThermalOverstress> {
        self.overstress.as_ref()
    }

    /// Source viscosity at one admitted temperature, without inventing a step.
    /// The logarithmic evaluation admits a finite product even when forming
    /// its bare exponential multiplier would overflow or underflow.
    pub fn overstress_viscosity_pa_s(
        &self,
        temperature_k: f64,
    ) -> Result<Option<f64>, FiniteJ2ThermalError> {
        self.admit_temperature(temperature_k)?;
        let Some(parameters) = self.overstress else {
            return Ok(None);
        };
        let reference = self.parameters.reference_temperature_k;
        if temperature_k == reference || parameters.activation_temperature_k == 0.0 {
            return Ok(Some(parameters.reference_viscosity_pa_s));
        }
        // Equivalent to A*(1/T-1/T_ref), without subtracting close reciprocals
        // or forming a potentially overflowing T*T_ref product.
        let log_multiplier = finite(
            ((reference - temperature_k) / temperature_k)
                * (parameters.activation_temperature_k / reference),
        )?;
        let viscosity = det::exp(finite(
            det::ln(parameters.reference_viscosity_pa_s) + log_multiplier,
        )?);
        if !viscosity.is_finite() || viscosity <= 0.0 {
            return Err(FiniteJ2ThermalError::Numerical(
                "Arrhenius viscosity must remain representable and positive over its temperature domain",
            ));
        }
        Ok(Some(viscosity))
    }

    /// Duration-free constitutive law for the rate-independent selection.
    /// Explicit kinetics require [`Self::at_temperature_for_step`].
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<FrozenThermalFiniteJ2, FiniteJ2ThermalError> {
        if self.overstress.is_some() {
            return Err(FiniteJ2ThermalError::DurationRequired);
        }
        self.freeze(temperature_k, None)
    }

    /// Freeze one actual positive physical interval at the trial temperature.
    /// The rate-independent selection retains its original return arithmetic.
    pub fn at_temperature_for_step(
        &self,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<FrozenThermalFiniteJ2, FiniteJ2ThermalError> {
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(FiniteJ2ThermalError::Parameters(
                "finite thermoplastic step duration must be positive and finite",
            ));
        }
        self.freeze(temperature_k, Some(duration_s))
    }

    /// Inspect fixed-history thermodynamics or an accepted state without a
    /// physical interval. With selected kinetics, response and stress/potential
    /// calls on this value refuse instead of silently dropping rate dependence.
    pub fn at_temperature_for_inspection(
        &self,
        temperature_k: f64,
    ) -> Result<FrozenThermalFiniteJ2, FiniteJ2ThermalError> {
        self.freeze(temperature_k, None)
    }

    /// Read mechanical energy/entropy at fixed total deformation and history.
    pub fn thermodynamics(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
        temperature_k: f64,
    ) -> Result<FiniteJ2Thermodynamics, FiniteJ2ThermalError> {
        self.at_temperature_for_inspection(temperature_k)?
            .thermodynamics(f, state)
    }

    /// Check an accepted state at its actual temperature, without returning it.
    /// Selected kinetics admit finite positive overstress; rate independence
    /// retains its explicit yield-surface gate.
    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
        temperature_k: f64,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2ThermalError> {
        self.at_temperature_for_inspection(temperature_k)?
            .validate_committed(f, state, yield_tolerance_pa)
    }

    fn admit_temperature(&self, temperature_k: f64) -> Result<(), FiniteJ2ThermalError> {
        let p = &self.parameters;
        let [minimum_k, maximum_k] = p.temperature_range_k;
        if !temperature_k.is_finite() || !(minimum_k..=maximum_k).contains(&temperature_k) {
            return Err(FiniteJ2ThermalError::Temperature {
                temperature_k,
                minimum_k,
                maximum_k,
            });
        }
        Ok(())
    }

    fn freeze(
        &self,
        temperature_k: f64,
        duration_s: Option<f64>,
    ) -> Result<FrozenThermalFiniteJ2, FiniteJ2ThermalError> {
        self.admit_temperature(temperature_k)?;
        let p = &self.parameters;
        let delta = temperature_k - p.reference_temperature_k;
        let coefficient = |value: f64, slope: f64| finite(slope.mul_add(delta, value));
        let native = self.native_backbone(FiniteJ2Parameters {
            shear_modulus_pa: coefficient(
                p.reference.shear_modulus_pa,
                p.slopes.shear_modulus_pa_k,
            )?,
            bulk_modulus_pa: coefficient(p.reference.bulk_modulus_pa, p.slopes.bulk_modulus_pa_k)?,
            yield_stress_pa: coefficient(p.reference.yield_stress_pa, p.slopes.yield_stress_pa_k)?,
            hardening_modulus_pa: coefficient(
                p.reference.hardening_modulus_pa,
                p.slopes.hardening_modulus_pa_k,
            )?,
        })?;
        let expansion_per_k = coefficient(
            p.logarithmic_expansion_per_k,
            p.logarithmic_expansion_slope_per_k2,
        )?;
        let log_thermal_stretch = finite(
            delta
                * (0.5 * p.logarithmic_expansion_slope_per_k2)
                    .mul_add(delta, p.logarithmic_expansion_per_k),
        )?;
        let thermal_stretch = det::exp(log_thermal_stretch);
        let inverse_thermal_stretch = thermal_stretch.recip();
        if [
            thermal_stretch,
            inverse_thermal_stretch,
            thermal_stretch * thermal_stretch * thermal_stretch,
            inverse_thermal_stretch * inverse_thermal_stretch * inverse_thermal_stretch,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(FiniteJ2ThermalError::Numerical(
                "thermal stretch and volume ratio must be representable and positive",
            ));
        }
        let native = if let Some(kinetics) = self.overstress {
            if let Some(duration) = duration_s {
                FrozenNative::Overstress(FiniteOverstressJ2::new(
                    native,
                    FiniteOverstressParameters {
                        viscosity_pa_s: self
                            .overstress_viscosity_pa_s(temperature_k)?
                            .expect("selected kinetics have a viscosity"),
                        stress_scale_pa: kinetics.stress_scale_pa,
                        exponent: kinetics.exponent,
                    },
                    duration,
                )?)
            } else {
                FrozenNative::Inspection(native)
            }
        } else {
            FrozenNative::RateIndependent(native)
        };
        Ok(FrozenThermalFiniteJ2 {
            native,
            slopes: p.slopes,
            temperature_k,
            log_thermal_stretch,
            thermal_stretch,
            inverse_thermal_stretch,
            expansion_per_k,
            expansion_slope_per_k2: p.logarithmic_expansion_slope_per_k2,
        })
    }

    /// Endpoint-temperature return with the old accepted law checked at OLD T.
    /// Inspection at new coefficients never re-returns or resets old history.
    #[allow(clippy::too_many_arguments)]
    pub fn checked_increment(
        &self,
        old_f: &[f64; 9],
        old_state: &FiniteJ2State,
        old_temperature_k: f64,
        new_f: &[f64; 9],
        new_temperature_k: f64,
        yield_tolerance_pa: f64,
    ) -> Result<FiniteJ2ThermalIncrement, FiniteJ2ThermalError> {
        let old_law = self.at_temperature(old_temperature_k)?;
        let new_law = self.at_temperature(new_temperature_k)?;
        Self::increment(
            old_f,
            old_state,
            new_f,
            yield_tolerance_pa,
            old_law,
            new_law,
        )
    }

    /// Thermal and kinetic update over one explicitly supplied physical interval.
    /// Old history is inspected at its old temperature, without advancing it
    /// again or requiring it to lie inside the rate-independent yield surface.
    #[allow(clippy::too_many_arguments)]
    pub fn checked_increment_for_step(
        &self,
        old_f: &[f64; 9],
        old_state: &FiniteJ2State,
        old_temperature_k: f64,
        new_f: &[f64; 9],
        new_temperature_k: f64,
        yield_tolerance_pa: f64,
        duration_s: f64,
    ) -> Result<FiniteJ2ThermalIncrement, FiniteJ2ThermalError> {
        let old_law = self.at_temperature_for_inspection(old_temperature_k)?;
        let new_law = self.at_temperature_for_step(new_temperature_k, duration_s)?;
        Self::increment(
            old_f,
            old_state,
            new_f,
            yield_tolerance_pa,
            old_law,
            new_law,
        )
    }

    fn increment(
        old_f: &[f64; 9],
        old_state: &FiniteJ2State,
        new_f: &[f64; 9],
        yield_tolerance_pa: f64,
        old_law: FrozenThermalFiniteJ2,
        new_law: FrozenThermalFiniteJ2,
    ) -> Result<FiniteJ2ThermalIncrement, FiniteJ2ThermalError> {
        let old_temperature_k = old_law.temperature_k;
        let new_temperature_k = new_law.temperature_k;
        old_law.validate_committed(old_f, old_state, yield_tolerance_pa)?;
        let old_endpoint = old_law.inspect(old_f, old_state)?;
        let old_fm = old_law.mechanical_deformation(old_f)?.0;
        let (new_fm, new_j) = new_law.mechanical_deformation(new_f)?;
        let native = new_law.native_response(&new_fm, old_state)?;
        let native_piola = native.first_piola_pa;
        let mechanical = new_law.map_response(native, new_j)?;
        let new_endpoint = new_law.inspect(new_f, &mechanical.state)?;
        new_law.validate_committed(new_f, &mechanical.state, yield_tolerance_pa)?;
        let old = old_endpoint.thermodynamics;
        let new = new_endpoint.thermodynamics;
        let old_p = old_law.backbone().parameters();
        let new_p = new_law.backbone().parameters();
        let coefficient_storage_change_j_m3 = finite(
            (new_p.shear_modulus_pa - old_p.shear_modulus_pa) * old_endpoint.deviator_squared
                + 0.5
                    * (new_p.bulk_modulus_pa - old_p.bulk_modulus_pa)
                    * old_endpoint.log_volume.powi(2)
                + 0.5
                    * (new_p.hardening_modulus_pa - old_p.hardening_modulus_pa)
                    * old_state.alpha.powi(2),
        )?;
        let frozen_old_storage = finite(
            new_p.shear_modulus_pa * old_endpoint.deviator_squared
                + 0.5 * new_p.bulk_modulus_pa * old_endpoint.log_volume.powi(2)
                + 0.5 * new_p.hardening_modulus_pa * old_state.alpha.powi(2),
        )?;
        let delta_f = core::array::from_fn(|i| new_f[i] - old_f[i]);
        let delta_fm = core::array::from_fn(|i| new_fm[i] - old_fm[i]);
        let mechanical_work_j_m3 = finite(super::dot(&mechanical.first_piola_pa, &delta_f))?;
        let mechanical_coordinate_work_j_m3 = finite(super::dot(&native_piola, &delta_fm))?;
        // Integrate this interval directly, avoiding subtraction of large
        // log-stretch primitives relative to an old reference temperature.
        let delta_log_stretch = finite(
            (new_temperature_k - old_temperature_k)
                * (0.5 * old_law.expansion_per_k + 0.5 * new_law.expansion_per_k),
        )?;
        let thermal_deformation_work_j_m3 =
            finite(super::dot(&mechanical.first_piola_pa, old_f) * det::expm1(delta_log_stretch))?;
        let thermoelastic_storage_change_j_m3 =
            finite((new_temperature_k - old_temperature_k).mul_add(
                old.entropy_j_m3_k,
                new_temperature_k * (new.entropy_j_m3_k - old.entropy_j_m3_k),
            ))?;
        let endpoint_work_defect_j_m3 = finite(
            mechanical_coordinate_work_j_m3
                - (new.free_energy_j_m3 - frozen_old_storage)
                - mechanical.plastic_dissipation_j_m3,
        )?;
        let reversible_heat_j_m3 = finite(
            thermal_deformation_work_j_m3
                - coefficient_storage_change_j_m3
                - thermoelastic_storage_change_j_m3,
        )?;
        let balance_error_j_m3 = finite(
            mechanical_work_j_m3
                - (new.free_energy_j_m3 - old.free_energy_j_m3 + thermoelastic_storage_change_j_m3)
                - mechanical.plastic_dissipation_j_m3
                - reversible_heat_j_m3
                - endpoint_work_defect_j_m3,
        )?;
        finite(mechanical.plastic_dissipation_j_m3 + reversible_heat_j_m3)?;
        Ok(FiniteJ2ThermalIncrement {
            mechanical,
            old,
            new,
            coefficient_storage_change_j_m3,
            thermoelastic_storage_change_j_m3,
            mechanical_work_j_m3,
            mechanical_coordinate_work_j_m3,
            thermal_deformation_work_j_m3,
            endpoint_work_defect_j_m3,
            reversible_heat_j_m3,
            balance_error_j_m3,
        })
    }
}

struct Inspection {
    thermodynamics: FiniteJ2Thermodynamics,
    deviator_squared: f64,
    log_volume: f64,
}

impl FrozenThermalFiniteJ2 {
    fn backbone(&self) -> &FiniteJ2 {
        match &self.native {
            FrozenNative::RateIndependent(law) | FrozenNative::Inspection(law) => law,
            FrozenNative::Overstress(law) => law.backbone(),
        }
    }

    fn native_response(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
    ) -> Result<FiniteJ2Response, FiniteJ2ThermalError> {
        Ok(match &self.native {
            FrozenNative::RateIndependent(law) => law.response(f, old)?,
            FrozenNative::Overstress(law) => law.response(f, old)?,
            FrozenNative::Inspection(_) => return Err(FiniteJ2ThermalError::DurationRequired),
        })
    }

    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub fn thermal_stretch(&self) -> f64 {
        self.thermal_stretch
    }
    pub fn log_thermal_stretch(&self) -> f64 {
        self.log_thermal_stretch
    }

    fn mechanical_deformation(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 9], f64), FiniteJ2ThermalError> {
        let j = super::determinant(f);
        if f.iter().any(|x| !x.is_finite()) || !j.is_finite() || j <= 0.0 {
            return Err(FiniteJ2Error::Domain {
                quantity: "total det(F), positive orientation",
                value: j,
                limit: 0.0,
            }
            .into());
        }
        Ok((f.map(|x| x * self.inverse_thermal_stretch), j))
    }

    fn map_response(
        &self,
        mut response: FiniteJ2Response,
        j: f64,
    ) -> Result<FiniteJ2Response, FiniteJ2ThermalError> {
        response.first_piola_pa = response
            .first_piola_pa
            .map(|x| x * self.inverse_thermal_stretch);
        let scale = self.inverse_thermal_stretch * self.inverse_thermal_stretch;
        response.tangent_pa = response.tangent_pa.map(|row| row.map(|x| x * scale));
        response.jacobian = j;
        if response
            .first_piola_pa
            .iter()
            .chain(response.tangent_pa.iter().flatten())
            .any(|x| !x.is_finite())
        {
            return Err(FiniteJ2ThermalError::Numerical(
                "thermal Piola/tangent chain is not representable",
            ));
        }
        Ok(response)
    }

    pub fn response(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
    ) -> Result<FiniteJ2Response, FiniteJ2ThermalError> {
        let (fm, j) = self.mechanical_deformation(f)?;
        self.map_response(self.native_response(&fm, old)?, j)
    }

    pub fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &FiniteJ2State,
    ) -> Result<([f64; 9], f64), FiniteJ2ThermalError> {
        let (fm, _) = self.mechanical_deformation(f)?;
        let (piola, potential) = match &self.native {
            FrozenNative::RateIndependent(law) => law.stress_and_potential(&fm, old)?,
            FrozenNative::Overstress(law) => law.stress_and_potential(&fm, old)?,
            FrozenNative::Inspection(_) => return Err(FiniteJ2ThermalError::DurationRequired),
        };
        let piola = piola.map(|x| x * self.inverse_thermal_stretch);
        if piola.iter().any(|x| !x.is_finite()) {
            return Err(FiniteJ2ThermalError::Numerical(
                "thermal Piola chain is not representable",
            ));
        }
        Ok((piola, potential))
    }

    pub fn validate_history(&self, state: &FiniteJ2State) -> Result<(), FiniteJ2ThermalError> {
        Ok(self.backbone().validate_history(state)?)
    }

    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
        yield_tolerance_pa: f64,
    ) -> Result<(), FiniteJ2ThermalError> {
        let (fm, _) = self.mechanical_deformation(f)?;
        match &self.native {
            FrozenNative::RateIndependent(law) => {
                law.validate_committed(&fm, state, yield_tolerance_pa)?;
            }
            FrozenNative::Overstress(law) => {
                law.validate_committed(&fm, state, yield_tolerance_pa)?;
            }
            FrozenNative::Inspection(_) => {
                if !yield_tolerance_pa.is_finite() || yield_tolerance_pa < 0.0 {
                    return Err(FiniteJ2Error::Parameters(
                        "yield tolerance must be finite and nonnegative",
                    )
                    .into());
                }
                // Overstress is a valid committed kinetic state. Inspect the
                // original history, elastic domain, finite stress and thermal
                // storage without a return or a fictitious physical interval.
                self.inspect(f, state)?;
            }
        }
        Ok(())
    }

    pub fn thermodynamics(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
    ) -> Result<FiniteJ2Thermodynamics, FiniteJ2ThermalError> {
        Ok(self.inspect(f, state)?.thermodynamics)
    }

    fn inspect(
        &self,
        f: &[f64; 9],
        state: &FiniteJ2State,
    ) -> Result<Inspection, FiniteJ2ThermalError> {
        let (fm, jacobian) = self.mechanical_deformation(f)?;
        let endpoint = self.backbone().elastic_endpoint(&fm, state)?;
        self.backbone()
            .check_elastic_domain(endpoint.log_j, &endpoint.deviator)?;
        let p = self.backbone().parameters();
        let s = &self.slopes;
        let g = endpoint.log_j;
        let a = self.expansion_per_k;
        let shape = super::dot(&endpoint.deviator, &endpoint.deviator);
        let free_energy_j_m3 = self.backbone().energy(g, &endpoint.deviator, state.alpha)?;
        let entropy_j_m3_k = finite(
            3.0 * p.bulk_modulus_pa * a * g
                - s.shear_modulus_pa_k * shape
                - 0.5 * s.bulk_modulus_pa_k * g * g
                - 0.5 * s.hardening_modulus_pa_k * state.alpha * state.alpha,
        )?;
        let internal_energy_j_m3 =
            finite(self.temperature_k.mul_add(entropy_j_m3_k, free_energy_j_m3))?;
        let heat_capacity_j_m3_k = finite(
            self.temperature_k
                * ((6.0 * s.bulk_modulus_pa_k * a
                    + 3.0 * p.bulk_modulus_pa * self.expansion_slope_per_k2)
                    * g
                    - 9.0 * p.bulk_modulus_pa * a * a),
        )?;
        let mut tau = endpoint.deviator.map(|e| 2.0 * p.shear_modulus_pa * e);
        let mut tau_temperature = endpoint.deviator.map(|e| 2.0 * s.shear_modulus_pa_k * e);
        for i in [0, 4, 8] {
            tau[i] += p.bulk_modulus_pa * g;
            tau_temperature[i] += s.bulk_modulus_pa_k * g - 3.0 * p.bulk_modulus_pa * a;
        }
        let inverse_transpose = super::transpose(&super::inverse(f, jacobian)?);
        let first_piola_pa = super::multiply(&tau, &inverse_transpose);
        let derivative = super::multiply(&tau_temperature, &inverse_transpose);
        if first_piola_pa
            .iter()
            .chain(derivative.iter())
            .any(|x| !x.is_finite())
        {
            return Err(FiniteJ2ThermalError::Numerical(
                "fixed-history thermodynamic stress is not finite",
            ));
        }
        Ok(Inspection {
            thermodynamics: FiniteJ2Thermodynamics {
                first_piola_pa,
                first_piola_temperature_derivative_fixed_history_pa_k: derivative,
                free_energy_j_m3,
                entropy_j_m3_k,
                internal_energy_j_m3,
                heat_capacity_j_m3_k,
                jacobian,
                mechanical_jacobian: endpoint.jacobian,
            },
            deviator_squared: shape,
            log_volume: g,
        })
    }
}
