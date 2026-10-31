//! Finite-rate liquid evaporation and condensation through a carrier-gas film.
//!
//! The carrier has zero net molar flux through the liquid interface. Integrating
//! binary ideal-gas diffusion at constant total pressure and film temperature
//! gives `j = M p h_m/(R T_f) ln[(p-p_v)/(p-p_sat)]`, with `h_m = D/delta`.
//! This retains the Stefan flow correction instead of using a dilute pressure
//! difference at every concentration. See ASHRAE Fundamentals (2021), chapter 6,
//! "Diffusion of One Gas Through a Second Stagnant Gas", equations (11)--(17):
//! <https://handbook.ashrae.org/Handbooks/F21/SI/F21_Ch06/F21_Ch06_si.aspx>.
//!
//! Positive flux evaporates liquid; negative flux condenses onto an EXISTING
//! liquid surface. The caller supplies the gas-film temperature and coefficient;
//! neither is inferred from a material name or a heat-transfer coefficient.
//! The film approximation neglects thermal diffusion, gas nonideality, curvature,
//! interfacial kinetic resistance and nucleation. Pure vapor and boiling refuse.
//! A remote bath is prescribed, so its mass, momentum and sensible heat are not
//! evolved here. Finite-atmosphere owners must close those accounts themselves.
//!
//! The implicit surface step uses the original liquid enthalpy inverse. Area
//! stays fixed during a step; no mesh or liquid-flow state is changed. It is a
//! first-order, uniform-liquid model, not a spatial Stefan or droplet solver.
//! Scalar rates and the bounded nonlinear solve allocate no memory. Shared
//! deterministic elementary functions retain same-ISA arithmetic repeatability.

/// Molar gas constant from the exact SI defining constants [J/(mol K)].
pub const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618_153_24;

/// Every failure leaves the caller's accepted liquid and atmosphere unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EvaporationError {
    Invalid(&'static str),
    TemperatureOutOfRange {
        temperature_k: f64,
        minimum_k: f64,
        maximum_k: f64,
    },
    /// No stagnant carrier film exists at or above the local boiling pressure.
    Boiling {
        saturation_pressure_pa: f64,
        total_pressure_pa: f64,
    },
    /// The caller's pure-liquid caloric interval cannot contain this endpoint.
    Unbracketed {
        residual_j: [f64; 2],
    },
    /// A smaller interval is required before liquid extinction/topology change.
    InventoryLimit {
        /// Suggested retry below the initial frozen-rate inventory limit.
        /// Arbitrary caller caloric closures need not be monotone; this is not
        /// a guaranteed accepted step. All source and energy checks still apply.
        maximum_duration_s: f64,
    },
    Unrepresentable(&'static str),
    NonConverged,
    Cancelled,
}
impl core::fmt::Display for EvaporationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "liquid evaporation: {self:?}")
    }
}
impl std::error::Error for EvaporationError {}

/// One immutable liquid/vapor interface evaluation in coherent SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvaporationProperties {
    pub saturation_pressure_pa: f64,
    pub saturation_pressure_derivative_pa_k: f64,
    /// Vapor stream enthalpy minus liquid enthalpy at the SAME temperature.
    /// Pressure conventions and caloric references belong to the material law.
    pub latent_heat_j_kg: f64,
    pub liquid_density_kg_m3: f64,
    pub molar_mass_kg_mol: f64,
}

/// Constitutive owner of a pure-liquid interface and its bounded source support.
pub trait EvaporationLaw {
    fn temperature_range_k(&self) -> [f64; 2];
    fn properties(&self, temperature_k: f64) -> Result<EvaporationProperties, EvaporationError>;
}

/// Prescribed far-field volatile pressure and one admitted binary gas film.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasFilm {
    pub total_pressure_pa: f64,
    /// Volatile partial pressure, never the total pressure or relative humidity.
    pub vapor_partial_pressure_pa: f64,
    /// Constant representative film temperature; its approximation is caller-owned.
    pub temperature_k: f64,
    /// Binary `D/delta` in m/s, BEFORE the Stefan correction; zero seals the film.
    pub mass_transfer_coefficient_m_s: f64,
}
impl GasFilm {
    /// Explicit diffusivity/film-thickness ingress; no empirical correlation is inferred.
    pub fn from_diffusion(
        total_pressure_pa: f64,
        vapor_partial_pressure_pa: f64,
        temperature_k: f64,
        diffusivity_m2_s: f64,
        thickness_m: f64,
    ) -> Result<Self, EvaporationError> {
        if !nonnegative(diffusivity_m2_s) || !positive(thickness_m) {
            return Err(EvaporationError::Invalid(
                "finite nonnegative diffusivity and positive film thickness required",
            ));
        }
        let film = Self {
            total_pressure_pa,
            vapor_partial_pressure_pa,
            temperature_k,
            mass_transfer_coefficient_m_s: checked(diffusivity_m2_s / thickness_m)?,
        };
        if diffusivity_m2_s > 0.0 && film.mass_transfer_coefficient_m_s == 0.0 {
            return Err(EvaporationError::Unrepresentable(
                "positive diffusivity/thickness underflows to zero transport",
            ));
        }
        film.validate()?;
        Ok(film)
    }

    pub fn validate(self) -> Result<(), EvaporationError> {
        if !positive(self.total_pressure_pa)
            || !nonnegative(self.vapor_partial_pressure_pa)
            || self.vapor_partial_pressure_pa >= self.total_pressure_pa
            || !positive(self.temperature_k)
            || !nonnegative(self.mass_transfer_coefficient_m_s)
        {
            return Err(EvaporationError::Invalid(
                "positive total pressure/film temperature, nonnegative coefficient and volatile pressure below total pressure required",
            ));
        }
        Ok(())
    }
}

/// Signed fluxes per unit actually exposed liquid area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceFlux {
    pub mass_flux_kg_m2_s: f64,
    /// Analytic derivative with film temperature, coefficient and pressures HELD.
    pub mass_flux_derivative_kg_m2_s_k: f64,
    /// Positive heat removed from liquid; negative condensation releases heat.
    pub latent_heat_flux_w_m2: f64,
    /// Negative normal velocity recedes into liquid; positive grows condensate.
    pub outward_normal_velocity_m_s: f64,
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn checked(x: f64) -> Result<f64, EvaporationError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(EvaporationError::Unrepresentable(
            "nonfinite surface quantity",
        ))
    }
}

/// Validate the original source once, shared by rate and thermal consumers.
pub fn properties<L: EvaporationLaw + ?Sized>(
    law: &L,
    temperature_k: f64,
) -> Result<EvaporationProperties, EvaporationError> {
    let [minimum_k, maximum_k] = law.temperature_range_k();
    if !positive(minimum_k) || !positive(maximum_k) || minimum_k > maximum_k {
        return Err(EvaporationError::Invalid(
            "invalid liquid source temperature interval",
        ));
    }
    if !temperature_k.is_finite() || !(minimum_k..=maximum_k).contains(&temperature_k) {
        return Err(EvaporationError::TemperatureOutOfRange {
            temperature_k,
            minimum_k,
            maximum_k,
        });
    }
    let p = law.properties(temperature_k)?;
    if !nonnegative(p.saturation_pressure_pa)
        || !nonnegative(p.saturation_pressure_derivative_pa_k)
        || !positive(p.latent_heat_j_kg)
        || !positive(p.liquid_density_kg_m3)
        || !positive(p.molar_mass_kg_mol)
    {
        return Err(EvaporationError::Invalid(
            "invalid liquid saturation, latent heat, density or molar mass",
        ));
    }
    Ok(p)
}

// ln(a/b), retaining the ORIGINAL pressure difference when a and b nearly agree.
// The atanh series has |s|<=1/8; the omitted 21st-order term is below one ulp.
fn log_carrier_ratio(a: f64, b: f64, difference: f64) -> f64 {
    let maximum = a.max(b);
    let s = (difference / maximum) / (1.0 + a.min(b) / maximum);
    if s.abs() <= 0.125 {
        let s2 = s * s;
        let series = 1.0
            + s2 * (1.0 / 3.0
                + s2 * (1.0 / 5.0
                    + s2 * (1.0 / 7.0
                        + s2 * (1.0 / 9.0
                            + s2 * (1.0 / 11.0
                                + s2 * (1.0 / 13.0
                                    + s2 * (1.0 / 15.0 + s2 * (1.0 / 17.0 + s2 / 19.0))))))));
        2.0 * s * series
    } else {
        let ratio = a / b;
        if positive(ratio) {
            fs_math::det::ln(ratio)
        } else {
            fs_math::det::ln(a) - fs_math::det::ln(b)
        }
    }
}

fn flux(p: EvaporationProperties, film: GasFilm) -> Result<SurfaceFlux, EvaporationError> {
    film.validate()?;
    if p.saturation_pressure_pa >= film.total_pressure_pa {
        return Err(EvaporationError::Boiling {
            saturation_pressure_pa: p.saturation_pressure_pa,
            total_pressure_pa: film.total_pressure_pa,
        });
    }
    if film.mass_transfer_coefficient_m_s == 0.0 {
        return Ok(SurfaceFlux {
            mass_flux_kg_m2_s: 0.0,
            mass_flux_derivative_kg_m2_s_k: 0.0,
            latent_heat_flux_w_m2: 0.0,
            outward_normal_velocity_m_s: 0.0,
        });
    }
    let carrier_surface = film.total_pressure_pa - p.saturation_pressure_pa;
    let carrier_far = film.total_pressure_pa - film.vapor_partial_pressure_pa;
    let pressure_gap = p.saturation_pressure_pa - film.vapor_partial_pressure_pa;
    let scale = checked(
        film.mass_transfer_coefficient_m_s
            * (p.molar_mass_kg_mol / MOLAR_GAS_CONSTANT)
            * (film.total_pressure_pa / film.temperature_k),
    )?;
    if scale == 0.0 {
        return Err(EvaporationError::Unrepresentable(
            "positive gas-film conductance underflows",
        ));
    }
    let mass = checked(scale * log_carrier_ratio(carrier_far, carrier_surface, pressure_gap))?;
    if pressure_gap != 0.0 && mass == 0.0 {
        return Err(EvaporationError::Unrepresentable(
            "nonzero vapor-pressure drive underflows to zero mass flux",
        ));
    }
    Ok(SurfaceFlux {
        mass_flux_kg_m2_s: mass,
        mass_flux_derivative_kg_m2_s_k: checked(
            scale * (p.saturation_pressure_derivative_pa_k / carrier_surface),
        )?,
        latent_heat_flux_w_m2: checked(mass * p.latent_heat_j_kg)?,
        outward_normal_velocity_m_s: checked(-mass / p.liquid_density_kg_m3)?,
    })
}

/// Signed continuum transfer; gas-film resistance and Stefan flow share one rate.
/// Constant-work evaluation, with no state mutation or allocation.
pub fn surface_flux<L: EvaporationLaw + ?Sized>(
    law: &L,
    surface_temperature_k: f64,
    film: GasFilm,
) -> Result<SurfaceFlux, EvaporationError> {
    flux(properties(law, surface_temperature_k)?, film)
}

/// Pure-liquid mass and ORIGINAL specific enthalpy; exposed area is caller-owned.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidSurfaceState {
    pub mass_kg: f64,
    pub specific_enthalpy_j_kg: f64,
    pub area_m2: f64,
}

/// Explicit pure-liquid domain and bounded nonlinear work. No default physical
/// enthalpy interval is inferred from a temperature or a material name.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolverConfig {
    pub specific_enthalpy_bracket_j_kg: [f64; 2],
    pub relative_tolerance: f64,
    pub absolute_energy_tolerance_j: f64,
    /// Total source/caloric evaluations, including the initial state, in 1..=4096.
    pub max_iterations: usize,
}

/// All heat, mass and recession accounts use the SAME represented transfer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceStep {
    pub liquid: LiquidSurfaceState,
    /// Positive liquid-to-vapor mass, negative condensation.
    pub transferred_mass_kg: f64,
    pub latent_energy_j: f64,
    /// Signed stream enthalpy in the liquid datum, already INCLUDING latent heat.
    /// The prescribed bath owns subsequent vapor sensible-heat equilibration.
    pub vapor_energy_j: f64,
    pub energy_residual_j: f64,
    /// Represented transfer minus backward-Euler rate integral [kg].
    pub mass_rate_residual_kg: f64,
    /// Positive liquid thickness loss, negative condensate growth.
    pub recession_m: f64,
    pub surface_temperature_k: f64,
}

#[derive(Clone, Copy)]
struct Endpoint {
    h: f64,
    t: f64,
    mass: f64,
    dm: f64,
    raw_dm: f64,
    latent: f64,
    residual: f64,
    density: f64,
}

/// Backward-Euler latent cooling/heating against the original liquid `T(h)`.
///
/// Conservation reduces to `m0*(h1-h0) + dm*L(T1) = 0`: transported vapor
/// carries `h1+L`, and latent heat must not be debited a second time. Bracket
/// endpoints must belong to the SAME pure-liquid chart. A missing bracket,
/// lost mass increment, exhausted inventory or cancellation publishes no state.
/// Condensation draws from a prescribed reservoir of unlimited inventory; a
/// finite-gas owner must provide its own inventory and energy acceptance.
#[allow(clippy::too_many_arguments)]
pub fn advance_surface<L: EvaporationLaw + ?Sized>(
    law: &L,
    film: GasFilm,
    liquid: LiquidSurfaceState,
    duration_s: f64,
    config: SolverConfig,
    temperature: impl Fn(f64) -> Result<f64, EvaporationError>,
    mut checkpoint: impl FnMut() -> Result<(), EvaporationError>,
) -> Result<SurfaceStep, EvaporationError> {
    let [lo, hi] = config.specific_enthalpy_bracket_j_kg;
    if !positive(liquid.mass_kg)
        || !nonnegative(liquid.area_m2)
        || !nonnegative(duration_s)
        || !lo.is_finite()
        || !hi.is_finite()
        || lo > hi
        || !liquid.specific_enthalpy_j_kg.is_finite()
        || !(lo..=hi).contains(&liquid.specific_enthalpy_j_kg)
        || !positive(config.relative_tolerance)
        || config.relative_tolerance >= 1.0
        || !nonnegative(config.absolute_energy_tolerance_j)
        || config.max_iterations == 0
        || config.max_iterations > 4096
    {
        return Err(EvaporationError::Invalid(
            "invalid liquid state, duration, pure-liquid enthalpy bracket or solve budget",
        ));
    }
    film.validate()?;
    let mut evaluations = 0;
    let mut evaluate = |h: f64| -> Result<Endpoint, EvaporationError> {
        checkpoint()?;
        if evaluations >= config.max_iterations {
            return Err(EvaporationError::NonConverged);
        }
        evaluations += 1;
        let t = temperature(h)?;
        let p = properties(law, t)?;
        let rate = flux(p, film)?.mass_flux_kg_m2_s;
        let raw_dm = checked(rate * liquid.area_m2 * duration_s)?;
        let mass = checked(liquid.mass_kg - raw_dm)?;
        let dm = liquid.mass_kg - mass;
        let latent = checked(dm * p.latent_heat_j_kg)?;
        let residual = checked(liquid.mass_kg * (h - liquid.specific_enthalpy_j_kg) + latent)?;
        Ok(Endpoint {
            h,
            t,
            mass,
            dm,
            raw_dm,
            latent,
            residual,
            density: p.liquid_density_kg_m3,
        })
    };
    let balanced = |e: Endpoint| {
        let sensible = liquid.mass_kg * (e.h - liquid.specific_enthalpy_j_kg);
        e.residual.abs()
            <= config.absolute_energy_tolerance_j
                + config.relative_tolerance * sensible.abs().max(e.latent.abs())
    };
    let initial = evaluate(liquid.specific_enthalpy_j_kg)?;
    let endpoint = if balanced(initial) {
        initial
    } else {
        let bound = if initial.residual > 0.0 { lo } else { hi };
        if bound == initial.h {
            return Err(EvaporationError::Unbracketed {
                residual_j: [initial.residual; 2],
            });
        }
        let far = evaluate(bound)?;
        if balanced(far) {
            far
        } else {
            if far.residual.signum() == initial.residual.signum() {
                return Err(EvaporationError::Unbracketed {
                    residual_j: [initial.residual, far.residual],
                });
            }
            let (mut left, mut right) = if far.h < initial.h {
                (far, initial)
            } else {
                (initial, far)
            };
            loop {
                let h = left.h + 0.5 * (right.h - left.h);
                if !h.is_finite() || h == left.h || h == right.h {
                    return Err(EvaporationError::NonConverged);
                }
                let middle = evaluate(h)?;
                if balanced(middle) {
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
    if endpoint.mass <= 0.0 {
        return Err(EvaporationError::InventoryLimit {
            maximum_duration_s: (duration_s * (liquid.mass_kg / initial.raw_dm))
                .next_down()
                .max(0.0),
        });
    }
    if endpoint.raw_dm != 0.0 && endpoint.dm == 0.0 {
        return Err(EvaporationError::Unrepresentable(
            "nonzero transfer rounds away in liquid mass",
        ));
    }
    let step = SurfaceStep {
        liquid: LiquidSurfaceState {
            mass_kg: endpoint.mass,
            specific_enthalpy_j_kg: endpoint.h,
            area_m2: liquid.area_m2,
        },
        transferred_mass_kg: endpoint.dm,
        latent_energy_j: endpoint.latent,
        vapor_energy_j: checked(endpoint.dm * endpoint.h + endpoint.latent)?,
        energy_residual_j: endpoint.residual,
        mass_rate_residual_kg: checked(endpoint.dm - endpoint.raw_dm)?,
        recession_m: if liquid.area_m2 == 0.0 {
            0.0
        } else {
            checked(endpoint.dm / liquid.area_m2 / endpoint.density)?
        },
        surface_temperature_k: endpoint.t,
    };
    checkpoint()?;
    Ok(step)
}

#[cfg(test)]
#[path = "evaporation/tests.rs"]
mod tests;
