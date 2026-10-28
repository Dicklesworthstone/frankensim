//! Signed solid/vapor exchange shared by thermal and material consumers.
//!
//! Positive flux removes solid; negative flux deposits vapor. One accepted
//! transfer supplies latent heat, carried enthalpy and normal recession.
//! Pressure means volatile PARTIAL pressure. Finite gas reservoirs are ideal
//! and thermostated. Optional gas-film transport is a dilute, linear resistance
//! model, not Stefan flow, pore transport, nucleation or a resolved vapor plume.
//!
//! Scalar evaluation is allocation-free. Finite-reservoir exchange integrates
//! back-pressure exponentially at frozen surface temperature. The coupled
//! endpoint solves latent feedback implicitly against the ORIGINAL caloric
//! chart. This is a first-order reduced surface model, not a spatial Stefan
//! solver or moving mesh. Transcendentals use the shared deterministic math
//! owner; this module's executed repeatability evidence is same-ISA only.
use core::fmt;

/// Molar gas constant [J/(mol K)] from the exact SI defining constants.
pub const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618_153_24;

/// No candidate state is published on refusal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SublimationError {
    Invalid(&'static str),
    TemperatureOutOfRange {
        temperature_k: f64,
        minimum_k: f64,
        maximum_k: f64,
    },
    NonConverged,
    Cancelled,
}
impl fmt::Display for SublimationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "sublimation: {self:?}")
    }
}
impl std::error::Error for SublimationError {}

/// One constitutive evaluation; derivatives are with respect to surface kelvin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SublimationProperties {
    pub saturation_pressure_pa: f64,
    pub saturation_pressure_derivative_pa_k: f64,
    pub latent_heat_j_kg: f64,
    pub latent_heat_derivative_j_kg_k: f64,
    pub solid_density_kg_m3: f64,
    pub molar_mass_kg_mol: f64,
}
/// Immutable solid/vapor chart; material crates own the constitutive data.
pub trait SublimationLaw {
    fn temperature_range_k(&self) -> [f64; 2];
    fn properties(&self, temperature_k: f64) -> Result<SublimationProperties, SublimationError>;
}
/// Equal incoming/outgoing accommodation preserves detailed balance at equal T.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceKinetics {
    /// Probability in [0,1]; zero seals the interface.
    pub accommodation_coefficient: f64,
    /// None: free-molecular exchange with LOCAL gas. Some(h_m): dilute film
    /// between interface and FAR-FIELD bath [m/s], including zero insulation.
    /// Pressure conductance is h_m/(R_specific*T_g).
    pub mass_transfer_coefficient_m_s: Option<f64>,
}
/// Carrier gas does not contribute to volatile partial pressure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VaporReservoir {
    Prescribed {
        partial_pressure_pa: f64,
        temperature_k: f64,
    },
    Finite {
        vapor_mass_kg: f64,
        volume_m3: f64,
        temperature_k: f64,
    },
}
impl VaporReservoir {
    pub fn partial_pressure_pa(self, molar_mass_kg_mol: f64) -> Result<f64, SublimationError> {
        if !positive(molar_mass_kg_mol) {
            return Err(SublimationError::Invalid(
                "positive finite molar mass required",
            ));
        }
        match self {
            Self::Prescribed {
                partial_pressure_pa,
                temperature_k,
            } => {
                if !nonnegative(partial_pressure_pa) || !positive(temperature_k) {
                    return Err(SublimationError::Invalid(
                        "invalid prescribed vapor reservoir",
                    ));
                }
                Ok(partial_pressure_pa)
            }
            Self::Finite {
                vapor_mass_kg,
                volume_m3,
                temperature_k,
            } => {
                if !nonnegative(vapor_mass_kg) || !positive(volume_m3) || !positive(temperature_k) {
                    return Err(SublimationError::Invalid("invalid finite vapor reservoir"));
                }
                finite(
                    vapor_mass_kg / volume_m3
                        * (MOLAR_GAS_CONSTANT / molar_mass_kg_mol)
                        * temperature_k,
                )
            }
        }
    }
    #[must_use]
    pub fn temperature_k(self) -> f64 {
        match self {
            Self::Prescribed { temperature_k, .. } | Self::Finite { temperature_k, .. } => {
                temperature_k
            }
        }
    }
}
/// Fluxes per unit EXPOSED surface area.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceFlux {
    /// kg/(m² s), positive solid to vapor.
    pub mass_flux_kg_m2_s: f64,
    pub mass_flux_derivative_kg_m2_s_k: f64,
    /// W/m² removed from solid; negative releases latent heat.
    pub latent_heat_flux_w_m2: f64,
    pub latent_heat_flux_derivative_w_m2_k: f64,
    /// Outward-solid-normal velocity; negative during sublimation.
    pub outward_normal_velocity_m_s: f64,
}
/// Specific enthalpy retains the original caloric datum. Geometry is caller-owned.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolidSurfaceState {
    pub mass_kg: f64,
    pub specific_enthalpy_j_kg: f64,
    pub area_m2: f64,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolverConfig {
    pub relative_tolerance: f64,
    pub absolute_energy_tolerance_j: f64,
    pub max_iterations: usize,
}
impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-10,
            absolute_energy_tolerance_j: 1e-10,
            max_iterations: 96,
        }
    }
}
/// All accounts describe the SAME accepted mass. Vapor energy already includes
/// latent energy; adding the latter again would double-count it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceStep {
    pub solid: SolidSurfaceState,
    pub reservoir: VaporReservoir,
    pub transferred_mass_kg: f64,
    pub latent_energy_j: f64,
    /// Signed bath energy in the solid's caloric datum. The thermostated bath
    /// owns subsequent gas sensible-heat equilibration.
    pub vapor_energy_j: f64,
    pub energy_residual_j: f64,
    /// Positive thickness loss, negative frost growth, at endpoint density.
    pub recession_m: f64,
    pub surface_temperature_k: f64,
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn finite(x: f64) -> Result<f64, SublimationError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(SublimationError::Invalid(
            "nonfinite derived surface quantity",
        ))
    }
}
fn properties<L: SublimationLaw + ?Sized>(
    law: &L,
    t: f64,
) -> Result<SublimationProperties, SublimationError> {
    let [lo, hi] = law.temperature_range_k();
    if !positive(lo) || !positive(hi) || lo > hi {
        return Err(SublimationError::Invalid(
            "invalid constitutive temperature range",
        ));
    }
    if !t.is_finite() || t < lo || t > hi {
        return Err(SublimationError::TemperatureOutOfRange {
            temperature_k: t,
            minimum_k: lo,
            maximum_k: hi,
        });
    }
    let p = law.properties(t)?;
    if !nonnegative(p.saturation_pressure_pa)
        || !nonnegative(p.saturation_pressure_derivative_pa_k)
        || !positive(p.latent_heat_j_kg)
        || !p.latent_heat_derivative_j_kg_k.is_finite()
        || !positive(p.solid_density_kg_m3)
        || !positive(p.molar_mass_kg_mol)
    {
        return Err(SublimationError::Invalid(
            "inadmissible sublimation constitutive properties",
        ));
    }
    Ok(p)
}
// j=a*p_sat-b*p_v. Outgoing molecules sample T_s, incoming molecules T_g.
// Eliminate interface pressure: BOTH terms get 1/(1+b_kinetic/K_film).
fn coefficients(
    p: SublimationProperties,
    ts: f64,
    tg: f64,
    k: SurfaceKinetics,
) -> Result<(f64, f64), SublimationError> {
    if !k.accommodation_coefficient.is_finite()
        || !(0.0..=1.0).contains(&k.accommodation_coefficient)
        || k.mass_transfer_coefficient_m_s
            .is_some_and(|x| !nonnegative(x))
        || !positive(tg)
    {
        return Err(SublimationError::Invalid(
            "invalid surface kinetics or gas temperature",
        ));
    }
    let factor = k.accommodation_coefficient
        * (p.molar_mass_kg_mol / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT)).sqrt();
    let a = factor / ts.sqrt();
    let b = factor / tg.sqrt();
    let resistance = match k.mass_transfer_coefficient_m_s {
        None => 1.0,
        Some(0.0) => 0.0,
        Some(h) => {
            let film = h * p.molar_mass_kg_mol / (MOLAR_GAS_CONSTANT * tg);
            if film == 0.0 {
                0.0
            } else if b > film {
                let ratio = film / b;
                ratio / (1.0 + ratio)
            } else {
                1.0 / (1.0 + b / film)
            }
        }
    };
    Ok((finite(a * resistance)?, finite(b * resistance)?))
}
/// Signed Hertz–Knudsen exchange and analytic surface-temperature tangent.
/// The optional film is a series resistance, not vacuum kinetics at far field.
pub fn surface_flux<L: SublimationLaw + ?Sized>(
    law: &L,
    surface_temperature_k: f64,
    kinetics: SurfaceKinetics,
    reservoir: VaporReservoir,
) -> Result<SurfaceFlux, SublimationError> {
    let p = properties(law, surface_temperature_k)?;
    let pv = reservoir.partial_pressure_pa(p.molar_mass_kg_mol)?;
    let (a, b) = coefficients(
        p,
        surface_temperature_k,
        reservoir.temperature_k(),
        kinetics,
    )?;
    // Preserve exact equal-temperature detailed balance and small pressure gaps.
    let j = if surface_temperature_k == reservoir.temperature_k() {
        a * (p.saturation_pressure_pa - pv)
    } else {
        a * p.saturation_pressure_pa - b * pv
    };
    let dj = a
        * (p.saturation_pressure_derivative_pa_k
            - p.saturation_pressure_pa / (2.0 * surface_temperature_k));
    Ok(SurfaceFlux {
        mass_flux_kg_m2_s: finite(j)?,
        mass_flux_derivative_kg_m2_s_k: finite(dj)?,
        latent_heat_flux_w_m2: finite(j * p.latent_heat_j_kg)?,
        latent_heat_flux_derivative_w_m2_k: finite(
            dj * p.latent_heat_j_kg + j * p.latent_heat_derivative_j_kg_k,
        )?,
        outward_normal_velocity_m_s: finite(-j / p.solid_density_kg_m3)?,
    })
}
fn validate_surface(solid: SolidSurfaceState, duration_s: f64) -> Result<(), SublimationError> {
    if !positive(solid.mass_kg)
        || !solid.specific_enthalpy_j_kg.is_finite()
        || !nonnegative(solid.area_m2)
        || !nonnegative(duration_s)
    {
        return Err(SublimationError::Invalid(
            "invalid surface mass, enthalpy, area or duration",
        ));
    }
    Ok(())
}
// Frozen-T exact relaxation. expm1 preserves tiny steps. Limit mass once;
// all energy and geometry accounts are derived from that accepted transfer.
fn transfer(
    p: SublimationProperties,
    t: f64,
    k: SurfaceKinetics,
    solid: SolidSurfaceState,
    reservoir: VaporReservoir,
    dt: f64,
) -> Result<f64, SublimationError> {
    let pv = reservoir.partial_pressure_pa(p.molar_mass_kg_mol)?;
    let (a, b) = coefficients(p, t, reservoir.temperature_k(), k)?;
    if solid.area_m2 == 0.0 || dt == 0.0 || (a == 0.0 && b == 0.0) {
        return Ok(0.0);
    }
    let raw = match reservoir {
        VaporReservoir::Prescribed { temperature_k, .. } => {
            let j = if t == temperature_k {
                a * (p.saturation_pressure_pa - pv)
            } else {
                a * p.saturation_pressure_pa - b * pv
            };
            finite(j * solid.area_m2 * dt)?
        }
        VaporReservoir::Finite {
            vapor_mass_kg,
            volume_m3,
            temperature_k,
        } => {
            if b == 0.0 {
                return Err(SublimationError::Invalid(
                    "underflowed finite-reservoir conductance",
                ));
            }
            let pressure_per_mass =
                finite((MOLAR_GAS_CONSTANT / p.molar_mass_kg_mol) * temperature_k / volume_m3)?;
            if !positive(pressure_per_mass) {
                return Err(SublimationError::Invalid("invalid gas compressibility"));
            }
            let equilibrium_pressure =
                finite(p.saturation_pressure_pa * (temperature_k / t).sqrt())?;
            let exponent = solid.area_m2 * b * pressure_per_mass * dt;
            if exponent.is_nan() {
                return Err(SublimationError::Invalid("invalid gas relaxation time"));
            }
            let fraction = -fs_math::det::expm1(-exponent);
            finite((equilibrium_pressure - pv) / pressure_per_mass * fraction)?.max(-vapor_mass_kg)
        }
    };
    Ok(raw.min(solid.mass_kg))
}
/// Mass exchange at imposed surface T. The caller supplies the work needed to
/// HOLD that temperature. For freely evolving solid T use advance_surface.
pub fn isothermal_mass_transfer<L: SublimationLaw + ?Sized>(
    law: &L,
    surface_temperature_k: f64,
    kinetics: SurfaceKinetics,
    solid: SolidSurfaceState,
    reservoir: VaporReservoir,
    duration_s: f64,
) -> Result<f64, SublimationError> {
    validate_surface(solid, duration_s)?;
    transfer(
        properties(law, surface_temperature_k)?,
        surface_temperature_k,
        kinetics,
        solid,
        reservoir,
        duration_s,
    )
}
#[derive(Clone, Copy)]
struct Endpoint {
    h: f64,
    t: f64,
    p: SublimationProperties,
    dm: f64,
    residual: f64,
}

/// Implicit latent cooling/heating. The callback maps the SAME stored specific
/// enthalpy to kelvin; no surrogate heat capacity or hidden material is used.
/// Endpoint vapor carries h_end+L. Conservation reduces to
/// m0*(h_end-h0)+dm*L=0, avoiding subtraction of large energy datums.
/// A bounded bracketed solve refuses unsupported endpoints, never clamps T.
/// Each evaluation polls cancellation. A zero returned solid mass means
/// extinction: the consumer must handle its topology or refuse the proposal.
pub fn advance_surface<L: SublimationLaw + ?Sized>(
    law: &L,
    kinetics: SurfaceKinetics,
    solid: SolidSurfaceState,
    reservoir: VaporReservoir,
    duration_s: f64,
    config: SolverConfig,
    temperature: impl Fn(f64) -> Result<f64, SublimationError>,
    mut checkpoint: impl FnMut() -> Result<(), SublimationError>,
) -> Result<SurfaceStep, SublimationError> {
    validate_surface(solid, duration_s)?;
    if !positive(config.relative_tolerance)
        || config.relative_tolerance >= 1.0
        || !nonnegative(config.absolute_energy_tolerance_j)
        || config.max_iterations == 0
        || config.max_iterations > 4096
    {
        return Err(SublimationError::Invalid(
            "invalid sublimation nonlinear budget",
        ));
    }
    let mut evaluations = 0;
    let mut evaluate = |h: f64| -> Result<Endpoint, SublimationError> {
        checkpoint()?;
        if evaluations >= config.max_iterations {
            return Err(SublimationError::NonConverged);
        }
        evaluations += 1;
        finite(h)?;
        let t = temperature(h)?;
        let p = properties(law, t)?;
        let dm = transfer(p, t, kinetics, solid, reservoir, duration_s)?;
        let residual =
            finite(solid.mass_kg * (h - solid.specific_enthalpy_j_kg) + dm * p.latent_heat_j_kg)?;
        Ok(Endpoint {
            h,
            t,
            p,
            dm,
            residual,
        })
    };
    let initial = evaluate(solid.specific_enthalpy_j_kg)?;
    let tolerance = finite(
        config.absolute_energy_tolerance_j + config.relative_tolerance * initial.residual.abs(),
    )?;
    let endpoint = if initial.residual.abs() <= tolerance {
        initial
    } else {
        let mut near = initial;
        let mut trial_h = finite(initial.h - initial.residual / solid.mass_kg)?;
        let mut invalid_bound: Option<f64> = None;
        let far = loop {
            if trial_h == near.h {
                return Err(SublimationError::NonConverged);
            }
            match evaluate(trial_h) {
                Ok(value) => {
                    if value.residual.abs() <= tolerance
                        || value.residual.signum() != initial.residual.signum()
                    {
                        break value;
                    }
                    near = value;
                    trial_h = if let Some(bound) = invalid_bound {
                        near.h + (bound - near.h) * 0.5
                    } else {
                        finite(initial.h + 2.0 * (near.h - initial.h))?
                    };
                }
                Err(SublimationError::Cancelled) => return Err(SublimationError::Cancelled),
                Err(SublimationError::NonConverged) => return Err(SublimationError::NonConverged),
                Err(_) => {
                    invalid_bound = Some(trial_h);
                    trial_h = near.h + (trial_h - near.h) * 0.5;
                }
            }
        };
        if far.residual.abs() <= tolerance {
            far
        } else {
            let (mut left, mut right) = if near.h < far.h {
                (near, far)
            } else {
                (far, near)
            };
            loop {
                let h = left.h + (right.h - left.h) * 0.5;
                if h == left.h || h == right.h {
                    return Err(SublimationError::NonConverged);
                }
                let middle = evaluate(h)?;
                if middle.residual.abs() <= tolerance {
                    break middle;
                }
                if middle.residual.signum() == left.residual.signum() {
                    left = middle;
                } else {
                    right = middle;
                }
            }
        }
    };
    let remaining = finite(solid.mass_kg - endpoint.dm)?;
    let new_reservoir = match reservoir {
        VaporReservoir::Finite {
            vapor_mass_kg,
            volume_m3,
            temperature_k,
        } => VaporReservoir::Finite {
            vapor_mass_kg: finite(vapor_mass_kg + endpoint.dm)?,
            volume_m3,
            temperature_k,
        },
        bath => bath,
    };
    let latent = finite(endpoint.dm * endpoint.p.latent_heat_j_kg)?;
    let vapor_energy = finite(endpoint.dm * (endpoint.h + endpoint.p.latent_heat_j_kg))?;
    let recession = if solid.area_m2 == 0.0 {
        0.0
    } else {
        finite(endpoint.dm / solid.area_m2 / endpoint.p.solid_density_kg_m3)?
    };
    Ok(SurfaceStep {
        solid: SolidSurfaceState {
            mass_kg: remaining,
            specific_enthalpy_j_kg: endpoint.h,
            area_m2: solid.area_m2,
        },
        reservoir: new_reservoir,
        transferred_mass_kg: endpoint.dm,
        latent_energy_j: latent,
        vapor_energy_j: vapor_energy,
        energy_residual_j: endpoint.residual,
        recession_m: recession,
        surface_temperature_k: endpoint.t,
    })
}
#[cfg(test)]
mod tests;
