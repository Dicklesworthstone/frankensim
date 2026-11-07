//! Thermomechanical finite porous plasticity on the original compaction owner.
//!
//! The explicit decomposition is F=theta(T) Fe Fp, with isotropic logarithmic
//! expansion `d ln(theta)/dT=a(T)`. The native pressure-sensitive return receives
//! Fm=F/theta and retains its complete non-isochoric plastic metric. Its elastic
//! volume strain is g=ln(J/theta^3)-ln(Jp); plastic compaction is never mistaken
//! for thermal contraction. At fixed T the total-Piola tangent is the original
//! native tangent divided by theta^2, including all geometric/plastic terms.
//!
//! All energetic coefficients refer to the ORIGINAL volume and Kirchhoff stress,
//! as in the finite J2 thermal family. K, mu, Y0, H and optional Voce Q are
//! explicitly affine over a bounded temperature interval. Pressure sensitivity
//! and Voce b stay fixed. Free energy is mu*dev(Ee):dev(Ee)+K*g^2/2+H*alpha^2/2,
//! plus the original Voce primitive. Entropy is -partial(psi)/partial(T) at fixed
//! total F and plastic history; internal energy is psi+T*s. Caloric storage
//! belongs to the thermal owner and must have positive combined heat capacity.
//!
//! The finite-step account separates physical yield heat, reversible thermal
//! transfer and SIGNED endpoint mechanical work error. The last is never heat.
//! This is a temperature-dependent extension of the native dry, symmetric,
//! rate-independent porous model. No changing pore fraction, pore pressure,
//! densification law, phase change or Taylor--Quinney fraction is inferred.
//! The logarithmic expansion/free-energy convention follows the native finite
//! J2 owner; this is not the full independent-plastic-entropy theory of
//! Simo--Miehe, doi:10.1016/0045-7825(92)90170-O.

use super::{
    FinitePorous, FinitePorousError, FinitePorousLimits, FinitePorousParameters,
    FinitePorousResponse, FinitePorousState,
};
use crate::plastic::finite::{determinant, deviator, dot, expand, inverse, multiply, transpose};
use fs_math::det;

/// Shared coefficient units and affine derivatives, with no kinematic conversion.
pub use crate::plastic::finite::thermal::{
    FiniteJ2TemperatureSlopes as FinitePorousTemperatureSlopes,
    FiniteJ2ThermalVoceHardening as FinitePorousThermalVoceHardening,
};

/// Exact material-card binding for the evolving-temperature family.
pub mod source;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FinitePorousThermalParameters {
    pub reference_temperature_k: f64,
    pub temperature_range_k: [f64; 2],
    /// Original-volume/Kirchhoff coefficients; pressure sensitivity stays fixed.
    pub reference: FinitePorousParameters,
    pub slopes: FinitePorousTemperatureSlopes,
    /// Explicit instantaneous logarithmic expansion at T_ref, d ln(theta)/dT.
    pub logarithmic_expansion_per_k: f64,
    pub logarithmic_expansion_slope_per_k2: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AffineThermalFinitePorous {
    parameters: FinitePorousThermalParameters,
    limits: FinitePorousLimits,
    minimum_mechanical_heat_capacity_j_m3_k: f64,
    voce_hardening: Option<FinitePorousThermalVoceHardening>,
}

/// Borrow the complete family; an isothermal adapter cannot bypass its mapping.
pub trait FiniteThermalPorousMaterial {
    fn thermal_finite_porous(&self) -> &AffineThermalFinitePorous;
}
impl FiniteThermalPorousMaterial for AffineThermalFinitePorous {
    fn thermal_finite_porous(&self) -> &AffineThermalFinitePorous {
        self
    }
}
impl<T: FiniteThermalPorousMaterial + ?Sized> FiniteThermalPorousMaterial for &T {
    fn thermal_finite_porous(&self) -> &AffineThermalFinitePorous {
        T::thermal_finite_porous(*self)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct FrozenVoce {
    strength_pa: f64,
    rate_per_strain: f64,
    slope_pa_k: f64,
}

/// Complete endpoint-temperature law, sharing the original scalar/tensor return.
#[derive(Clone, Debug, PartialEq)]
pub struct FrozenThermalFinitePorous {
    native: FinitePorous,
    slopes: FinitePorousTemperatureSlopes,
    voce: Option<FrozenVoce>,
    temperature_k: f64,
    log_thermal_stretch: f64,
    thermal_stretch: f64,
    inverse_thermal_stretch: f64,
    expansion_per_k: f64,
    expansion_slope_per_k2: f64,
}

/// Fixed-history thermodynamics, without another constitutive return.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FinitePorousThermodynamics {
    pub first_piola_pa: [f64; 9],
    pub first_piola_temperature_derivative_fixed_history_pa_k: [f64; 9],
    pub free_energy_j_m3: f64,
    pub entropy_j_m3_k: f64,
    pub internal_energy_j_m3: f64,
    /// Mechanical contribution at fixed total F/history; it may be negative.
    pub heat_capacity_j_m3_k: f64,
    pub jacobian: f64,
    /// det(F/theta), containing both elastic and plastic volume changes.
    pub mechanical_jacobian: f64,
    pub plastic_jacobian: f64,
    /// det(Fe)=det(F/theta)/det(Fp), distinct from both volumes above.
    pub elastic_jacobian: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FinitePorousThermalIncrement {
    /// Native state and potential, with Piola/tangent mapped to total F.
    pub mechanical: FinitePorousResponse,
    pub old: FinitePorousThermodynamics,
    pub new: FinitePorousThermodynamics,
    pub coefficient_storage_change_j_m3: f64,
    /// Delta(T*s_mechanical), separate from the caloric chart.
    pub thermoelastic_storage_change_j_m3: f64,
    pub mechanical_work_j_m3: f64,
    pub mechanical_coordinate_work_j_m3: f64,
    pub thermal_deformation_work_j_m3: f64,
    /// Signed endpoint quadrature defect, never physical heating.
    pub endpoint_work_defect_j_m3: f64,
    pub reversible_heat_j_m3: f64,
    pub balance_error_j_m3: f64,
}
impl FinitePorousThermalIncrement {
    /// Deposit this signed transfer once in the separate caloric owner.
    pub fn caloric_heat_j_m3(&self) -> f64 {
        self.mechanical.plastic_dissipation_j_m3 + self.reversible_heat_j_m3
    }
    pub fn internal_energy_change_j_m3(&self) -> f64 {
        self.new.free_energy_j_m3 - self.old.free_energy_j_m3
            + self.thermoelastic_storage_change_j_m3
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum FinitePorousThermalError {
    Material(FinitePorousError),
    Parameters(&'static str),
    Temperature {
        temperature_k: f64,
        minimum_k: f64,
        maximum_k: f64,
    },
    Numerical(&'static str),
}
impl From<FinitePorousError> for FinitePorousThermalError {
    fn from(value: FinitePorousError) -> Self {
        Self::Material(value)
    }
}
impl core::fmt::Display for FinitePorousThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => error.fmt(f),
            Self::Parameters(what) => write!(f, "finite porous thermal parameters: {what}"),
            Self::Temperature {
                temperature_k,
                minimum_k,
                maximum_k,
            } => write!(
                f,
                "finite porous temperature {temperature_k} K outside [{minimum_k},{maximum_k}] K"
            ),
            Self::Numerical(what) => write!(f, "finite porous thermal numerical refusal: {what}"),
        }
    }
}
impl std::error::Error for FinitePorousThermalError {}

fn finite(value: f64) -> Result<f64, FinitePorousThermalError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FinitePorousThermalError::Numerical(
            "thermal coefficient, stretch or energy is unrepresentable",
        ))
    }
}
fn primitive(alpha: f64, rate: f64) -> f64 {
    crate::plastic::finite::voce::saturation_primitive(alpha, rate)
}

impl AffineThermalFinitePorous {
    /// Native Jacobian and strain bounds apply to mechanical F/theta. The
    /// spatial owner separately admits the actual total geometry.
    pub fn new(
        parameters: FinitePorousThermalParameters,
        limits: FinitePorousLimits,
    ) -> Result<Self, FinitePorousThermalError> {
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
            return Err(FinitePorousThermalError::Parameters(
                "finite affine derivatives and a positive T interval containing T_ref required",
            ));
        }
        FinitePorous::new(p.reference, limits)?;
        let mut family = Self {
            parameters,
            limits,
            minimum_mechanical_heat_capacity_j_m3_k: 0.0,
            voce_hardening: None,
        };
        let low = family.at_temperature(lo)?;
        let high = family.at_temperature(hi)?;
        if p.logarithmic_expansion_slope_per_k2 != 0.0 {
            let vertex = p.reference_temperature_k
                - p.logarithmic_expansion_per_k / p.logarithmic_expansion_slope_per_k2;
            if (lo..=hi).contains(&vertex) {
                family.at_temperature(vertex)?;
            }
        }
        let max_bulk = low
            .native
            .parameters()
            .bulk_modulus_pa
            .max(high.native.parameters().bulk_modulus_pa);
        let max_expansion = low.expansion_per_k.abs().max(high.expansion_per_k.abs());
        let max_elastic_volume = det::sqrt(3.0) * limits.maximum_elastic_log_strain;
        // The remaining coefficients, including optional Q, are affine at
        // fixed history. Only the thermal elastic-volume terms contribute Cv.
        let magnitude = finite(
            hi * ((6.0 * p.slopes.bulk_modulus_pa_k.abs() * max_expansion
                + 3.0 * max_bulk * p.logarithmic_expansion_slope_per_k2.abs())
                * max_elastic_volume
                + 9.0 * max_bulk * max_expansion * max_expansion),
        )?;
        family.minimum_mechanical_heat_capacity_j_m3_k = -magnitude;
        Ok(family)
    }

    pub const fn parameters(&self) -> &FinitePorousThermalParameters {
        &self.parameters
    }
    pub const fn limits(&self) -> &FinitePorousLimits {
        &self.limits
    }
    pub fn initial_state(&self) -> FinitePorousState {
        FinitePorousState {
            inverse_plastic_cauchy_green: super::IDENTITY6,
            alpha: 0.0,
        }
    }
    /// Lower bound at fixed total F/history over the full admitted envelope.
    /// Add to the original-volume caloric capacity and require a positive sum.
    /// This is not a uniqueness theorem for a coupled plastic/thermal return.
    pub const fn minimum_mechanical_heat_capacity_j_m3_k(&self) -> f64 {
        self.minimum_mechanical_heat_capacity_j_m3_k
    }
    pub const fn voce_hardening(&self) -> Option<&FinitePorousThermalVoceHardening> {
        self.voce_hardening.as_ref()
    }
    pub fn with_voce_hardening(
        mut self,
        reference_saturation_pa: f64,
        rate_per_strain: f64,
        saturation_slope_pa_k: f64,
    ) -> Result<Self, FinitePorousThermalError> {
        let hardening = FinitePorousThermalVoceHardening {
            reference_saturation_pa,
            rate_per_strain,
            saturation_slope_pa_k,
        };
        if !saturation_slope_pa_k.is_finite()
            || self.voce_hardening.is_some_and(|old| old != hardening)
        {
            return Err(FinitePorousThermalError::Parameters(
                "finite saturation derivative and no replacement of selected Voce coefficients required",
            ));
        }
        self.voce_hardening = Some(hardening);
        for temperature in self.parameters.temperature_range_k {
            self.at_temperature(temperature)?;
        }
        let entropy = finite(
            saturation_slope_pa_k
                * primitive(
                    self.limits.maximum_accumulated_plastic_strain,
                    rate_per_strain,
                ),
        )?;
        finite(self.parameters.temperature_range_k[1] * entropy)?;
        Ok(self)
    }

    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<FrozenThermalFinitePorous, FinitePorousThermalError> {
        let p = &self.parameters;
        let [minimum_k, maximum_k] = p.temperature_range_k;
        if !temperature_k.is_finite() || !(minimum_k..=maximum_k).contains(&temperature_k) {
            return Err(FinitePorousThermalError::Temperature {
                temperature_k,
                minimum_k,
                maximum_k,
            });
        }
        let delta = temperature_k - p.reference_temperature_k;
        let coefficient = |value: f64, slope: f64| finite(slope.mul_add(delta, value));
        let mut native = FinitePorous::new(
            FinitePorousParameters {
                shear_modulus_pa: coefficient(
                    p.reference.shear_modulus_pa,
                    p.slopes.shear_modulus_pa_k,
                )?,
                bulk_modulus_pa: coefficient(
                    p.reference.bulk_modulus_pa,
                    p.slopes.bulk_modulus_pa_k,
                )?,
                yield_stress_pa: coefficient(
                    p.reference.yield_stress_pa,
                    p.slopes.yield_stress_pa_k,
                )?,
                hardening_modulus_pa: coefficient(
                    p.reference.hardening_modulus_pa,
                    p.slopes.hardening_modulus_pa_k,
                )?,
                pressure_sensitivity: p.reference.pressure_sensitivity,
            },
            self.limits,
        )?;
        let voce = if let Some(v) = self.voce_hardening {
            let strength_pa = coefficient(v.reference_saturation_pa, v.saturation_slope_pa_k)?;
            native = native.with_voce_hardening(strength_pa, v.rate_per_strain)?;
            Some(FrozenVoce {
                strength_pa,
                rate_per_strain: v.rate_per_strain,
                slope_pa_k: v.saturation_slope_pa_k,
            })
        } else {
            None
        };
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
            return Err(FinitePorousThermalError::Numerical(
                "thermal stretch and volume must be positive and representable",
            ));
        }
        Ok(FrozenThermalFinitePorous {
            native,
            slopes: p.slopes,
            voce,
            temperature_k,
            log_thermal_stretch,
            thermal_stretch,
            inverse_thermal_stretch,
            expansion_per_k,
            expansion_slope_per_k2: p.logarithmic_expansion_slope_per_k2,
        })
    }

    pub fn thermodynamics(
        &self,
        f: &[f64; 9],
        state: &FinitePorousState,
        temperature_k: f64,
    ) -> Result<FinitePorousThermodynamics, FinitePorousThermalError> {
        self.at_temperature(temperature_k)?.thermodynamics(f, state)
    }

    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FinitePorousState,
        temperature_k: f64,
        yield_tolerance_pa: f64,
    ) -> Result<(), FinitePorousThermalError> {
        self.at_temperature(temperature_k)?
            .validate_committed(f, state, yield_tolerance_pa)
    }

    /// Old history is admitted at OLD T, then returned once at the proposed T.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn checked_increment(
        &self,
        old_f: &[f64; 9],
        old_state: &FinitePorousState,
        old_temperature_k: f64,
        new_f: &[f64; 9],
        new_temperature_k: f64,
        yield_tolerance_pa: f64,
    ) -> Result<FinitePorousThermalIncrement, FinitePorousThermalError> {
        let old_law = self.at_temperature(old_temperature_k)?;
        let new_law = self.at_temperature(new_temperature_k)?;
        old_law.validate_committed(old_f, old_state, yield_tolerance_pa)?;
        let old_endpoint = old_law.inspect(old_f, old_state)?;
        let old_fm = old_law.mechanical_deformation(old_f)?.0;
        let (new_fm, new_j) = new_law.mechanical_deformation(new_f)?;
        let native = new_law.native.response(&new_fm, old_state)?;
        let native_piola = native.first_piola_pa;
        let mechanical = new_law.map_response(native, new_j)?;
        let new_endpoint = new_law.inspect(new_f, &mechanical.state)?;
        new_law.validate_committed(new_f, &mechanical.state, yield_tolerance_pa)?;
        let old = old_endpoint.thermodynamics;
        let new = new_endpoint.thermodynamics;
        let p0 = old_law.native.parameters();
        let p1 = new_law.native.parameters();
        let mut coefficient_storage_change_j_m3 = finite(
            (p1.shear_modulus_pa - p0.shear_modulus_pa) * old_endpoint.deviator_squared
                + 0.5
                    * (p1.bulk_modulus_pa - p0.bulk_modulus_pa)
                    * old_endpoint.log_elastic_volume.powi(2)
                + 0.5
                    * (p1.hardening_modulus_pa - p0.hardening_modulus_pa)
                    * old_state.alpha.powi(2),
        )?;
        let mut frozen_old_storage = finite(
            p1.shear_modulus_pa * old_endpoint.deviator_squared
                + 0.5 * p1.bulk_modulus_pa * old_endpoint.log_elastic_volume.powi(2)
                + 0.5 * p1.hardening_modulus_pa * old_state.alpha.powi(2),
        )?;
        if let (Some(v0), Some(v1)) = (old_law.voce, new_law.voce) {
            let shape = primitive(old_state.alpha, v1.rate_per_strain);
            coefficient_storage_change_j_m3 = finite(
                coefficient_storage_change_j_m3 + (v1.strength_pa - v0.strength_pa) * shape,
            )?;
            frozen_old_storage = finite(frozen_old_storage + v1.strength_pa * shape)?;
        }
        let delta_f = core::array::from_fn(|i| new_f[i] - old_f[i]);
        let delta_fm = core::array::from_fn(|i| new_fm[i] - old_fm[i]);
        let mechanical_work_j_m3 = finite(dot(&mechanical.first_piola_pa, &delta_f))?;
        let mechanical_coordinate_work_j_m3 = finite(dot(&native_piola, &delta_fm))?;
        let delta_log_stretch = finite(
            (new_temperature_k - old_temperature_k)
                * (0.5 * old_law.expansion_per_k + 0.5 * new_law.expansion_per_k),
        )?;
        let thermal_deformation_work_j_m3 =
            finite(dot(&mechanical.first_piola_pa, old_f) * det::expm1(delta_log_stretch))?;
        let thermoelastic_storage_change_j_m3 =
            finite((new_temperature_k - old_temperature_k).mul_add(
                old.entropy_j_m3_k,
                new_temperature_k * (new.entropy_j_m3_k - old.entropy_j_m3_k),
            ))?;
        let heat = mechanical.plastic_dissipation_j_m3;
        let endpoint_work_defect_j_m3 = finite(
            mechanical_coordinate_work_j_m3 - (new.free_energy_j_m3 - frozen_old_storage) - heat,
        )?;
        let reversible_heat_j_m3 = finite(
            thermal_deformation_work_j_m3
                - coefficient_storage_change_j_m3
                - thermoelastic_storage_change_j_m3,
        )?;
        let balance_error_j_m3 = finite(
            mechanical_work_j_m3
                - (new.free_energy_j_m3 - old.free_energy_j_m3 + thermoelastic_storage_change_j_m3)
                - heat
                - reversible_heat_j_m3
                - endpoint_work_defect_j_m3,
        )?;
        finite(heat + reversible_heat_j_m3)?;
        Ok(FinitePorousThermalIncrement {
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
    thermodynamics: FinitePorousThermodynamics,
    deviator_squared: f64,
    log_elastic_volume: f64,
}

impl FrozenThermalFinitePorous {
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub const fn thermal_stretch(&self) -> f64 {
        self.thermal_stretch
    }
    pub const fn log_thermal_stretch(&self) -> f64 {
        self.log_thermal_stretch
    }

    fn mechanical_deformation(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 9], f64), FinitePorousThermalError> {
        let j = determinant(f);
        if f.iter().any(|x| !x.is_finite()) || !j.is_finite() || j <= 0.0 {
            return Err(FinitePorousError::Domain {
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
        mut response: FinitePorousResponse,
        j: f64,
    ) -> Result<FinitePorousResponse, FinitePorousThermalError> {
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
            return Err(FinitePorousThermalError::Numerical(
                "thermal Piola/tangent chain is unrepresentable",
            ));
        }
        Ok(response)
    }
    pub fn response(
        &self,
        f: &[f64; 9],
        old: &FinitePorousState,
    ) -> Result<FinitePorousResponse, FinitePorousThermalError> {
        let (fm, j) = self.mechanical_deformation(f)?;
        self.map_response(self.native.response(&fm, old)?, j)
    }
    pub fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &FinitePorousState,
    ) -> Result<([f64; 9], f64), FinitePorousThermalError> {
        let (fm, _) = self.mechanical_deformation(f)?;
        let (piola, potential) = self.native.stress_and_potential(&fm, old)?;
        let piola = piola.map(|x| x * self.inverse_thermal_stretch);
        if piola.iter().any(|x| !x.is_finite()) {
            return Err(FinitePorousThermalError::Numerical(
                "thermal Piola chain is unrepresentable",
            ));
        }
        Ok((piola, potential))
    }
    pub fn validate_history(
        &self,
        state: &FinitePorousState,
    ) -> Result<(), FinitePorousThermalError> {
        Ok(self.native.validate_history(state)?)
    }
    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FinitePorousState,
        yield_tolerance_pa: f64,
    ) -> Result<(), FinitePorousThermalError> {
        let (fm, _) = self.mechanical_deformation(f)?;
        Ok(self
            .native
            .validate_committed(&fm, state, yield_tolerance_pa)?)
    }
    pub fn thermodynamics(
        &self,
        f: &[f64; 9],
        state: &FinitePorousState,
    ) -> Result<FinitePorousThermodynamics, FinitePorousThermalError> {
        Ok(self.inspect(f, state)?.thermodynamics)
    }

    fn inspect(
        &self,
        f: &[f64; 9],
        state: &FinitePorousState,
    ) -> Result<Inspection, FinitePorousThermalError> {
        let (fm, jacobian) = self.mechanical_deformation(f)?;
        let endpoint = self.native.endpoint(&fm, state)?;
        self.native.elastic_bound(&endpoint.elastic)?;
        let elastic = expand(&endpoint.elastic);
        let dev = deviator(&elastic);
        let g = super::trace(&elastic);
        let shape = dot(&dev, &dev);
        let p = self.native.parameters();
        let s = &self.slopes;
        let a = self.expansion_per_k;
        let free_energy_j_m3 = self.native.energy(&endpoint.elastic, state.alpha)?;
        let mut entropy_j_m3_k = finite(
            3.0 * p.bulk_modulus_pa * a * g
                - s.shear_modulus_pa_k * shape
                - 0.5 * s.bulk_modulus_pa_k * g * g
                - 0.5 * s.hardening_modulus_pa_k * state.alpha * state.alpha,
        )?;
        if let Some(v) = self.voce {
            // Keep the derivative even when Q(T)=0 disables native saturation.
            entropy_j_m3_k =
                finite(entropy_j_m3_k - v.slope_pa_k * primitive(state.alpha, v.rate_per_strain))?;
        }
        let internal_energy_j_m3 =
            finite(self.temperature_k.mul_add(entropy_j_m3_k, free_energy_j_m3))?;
        let heat_capacity_j_m3_k = finite(
            self.temperature_k
                * ((6.0 * s.bulk_modulus_pa_k * a
                    + 3.0 * p.bulk_modulus_pa * self.expansion_slope_per_k2)
                    * g
                    - 9.0 * p.bulk_modulus_pa * a * a),
        )?;
        let mut tau = dev.map(|e| 2.0 * p.shear_modulus_pa * e);
        let mut tau_t = dev.map(|e| 2.0 * s.shear_modulus_pa_k * e);
        for i in [0, 4, 8] {
            tau[i] += p.bulk_modulus_pa * g;
            tau_t[i] += s.bulk_modulus_pa_k * g - 3.0 * p.bulk_modulus_pa * a;
        }
        let inverse_transpose = transpose(&inverse(f, jacobian).map_err(FinitePorousError::from)?);
        let first_piola_pa = multiply(&tau, &inverse_transpose);
        let derivative = multiply(&tau_t, &inverse_transpose);
        if first_piola_pa
            .iter()
            .chain(derivative.iter())
            .any(|x| !x.is_finite())
        {
            return Err(FinitePorousThermalError::Numerical(
                "fixed-history thermodynamic stress is unrepresentable",
            ));
        }
        let plastic_jacobian = super::volume(endpoint.log_plastic_jacobian)?;
        let elastic_jacobian = finite(endpoint.jacobian / plastic_jacobian)?;
        if elastic_jacobian <= 0.0 {
            return Err(FinitePorousThermalError::Numerical(
                "elastic volume ratio must remain represented and positive",
            ));
        }
        Ok(Inspection {
            thermodynamics: FinitePorousThermodynamics {
                first_piola_pa,
                first_piola_temperature_derivative_fixed_history_pa_k: derivative,
                free_energy_j_m3,
                entropy_j_m3_k,
                internal_energy_j_m3,
                heat_capacity_j_m3_k,
                jacobian,
                mechanical_jacobian: endpoint.jacobian,
                plastic_jacobian,
                elastic_jacobian,
            },
            deviator_squared: shape,
            log_elastic_volume: g,
        })
    }
}
