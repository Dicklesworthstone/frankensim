//! Finite interfacial kinetics in series with the original nonlinear Stefan film.
//!
//! The local vapor and liquid share the interface temperature `T_s`. With one
//! explicit, constant accommodation probability, the mass conductance is
//! `K = a sqrt(M/(2 pi R T_s))` (Hertz--Knudsen), or the same expression with
//! `a = 2 alpha/(2-alpha)` (Schrage). The unique interfacial pressure `p_i`
//! closes `j = K (p_sat-p_i) = C ln[(p-p_v)/(p-p_i)]`, where
//! `C = h_m M p/(R T_f)`. Positive `j` evaporates, negative `j` condenses.
//!
//! These are explicitly selected constitutive approximations, not measured
//! accommodation data or a resolved Knudsen layer. The local equal-temperature
//! assumption is distinct from the representative film temperature `T_f`.
//! Nonisothermal molecular heat, rarefied film transport, vapor-only flow,
//! boiling, curvature, and nucleation are outside this model.
//!
//! References: the interface formulas and their small-temperature-jump scope
//! are equations (7)--(8) in Lu et al., Nature Communications 10, 2368 (2019),
//! <https://doi.org/10.1038/s41467-019-10209-w>. The nonlinear stagnant-carrier
//! film is the original ASHRAE law cited by the parent module.

use super::{
    EvaporationError, EvaporationLaw, EvaporationProperties, GasFilm, MOLAR_GAS_CONSTANT,
    SurfaceFlux, checked, log_carrier_ratio, positive,
};
use crate::enthalpy::{
    DepositedCell, DepositionError, ThermalLaw, ThermalStorageBasis,
    advection::species::{CaloricComposition, CaloricEntropyAnchor},
};

/// A caller-selected immutable interfacial mass-transfer closure.
/// No finite accommodation probability is inferred or supplied by default.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum InterfaceKinetics {
    /// Local vapor/liquid equilibrium; preserves the original Stefan film.
    #[default]
    Equilibrium,
    /// Equal incoming/outgoing accommodation without the Schrage drift factor.
    HertzKnudsen {
        /// Dimensionless probability in [0,1]; zero seals the interface.
        accommodation_coefficient: f64,
    },
    /// Small-drift Schrage correction `2 alpha/(2-alpha)` to Hertz--Knudsen.
    Schrage {
        /// Dimensionless probability in [0,1]; zero seals the interface.
        accommodation_coefficient: f64,
    },
}

impl InterfaceKinetics {
    pub fn validate(self) -> Result<(), EvaporationError> {
        let alpha = match self {
            Self::Equilibrium => return Ok(()),
            Self::HertzKnudsen {
                accommodation_coefficient,
            }
            | Self::Schrage {
                accommodation_coefficient,
            } => accommodation_coefficient,
        };
        if !alpha.is_finite() || !(0.0..=1.0).contains(&alpha) {
            return Err(EvaporationError::Invalid(
                "interfacial accommodation must be finite and in [0,1]",
            ));
        }
        Ok(())
    }
}

/// Borrow the original evaporation/caloric owner and select interface kinetics.
/// This changes only mass-transfer resistance, never the material's saturation,
/// latent heat, density, caloric datum, source support, or storage interpretation.
pub struct KineticEvaporation<'a, L: EvaporationLaw + ?Sized> {
    source: &'a L,
    kinetics: InterfaceKinetics,
}
impl<L: EvaporationLaw + ?Sized> Copy for KineticEvaporation<'_, L> {}
impl<L: EvaporationLaw + ?Sized> Clone for KineticEvaporation<'_, L> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, L: EvaporationLaw + ?Sized> KineticEvaporation<'a, L> {
    pub fn try_new(source: &'a L, kinetics: InterfaceKinetics) -> Result<Self, EvaporationError> {
        kinetics.validate()?;
        Ok(Self { source, kinetics })
    }

    pub const fn source(&self) -> &'a L {
        self.source
    }
}

impl<L: EvaporationLaw + ?Sized> EvaporationLaw for KineticEvaporation<'_, L> {
    fn temperature_range_k(&self) -> [f64; 2] {
        self.source.temperature_range_k()
    }
    fn properties(&self, temperature_k: f64) -> Result<EvaporationProperties, EvaporationError> {
        self.source.properties(temperature_k)
    }
    fn interface_kinetics(&self) -> InterfaceKinetics {
        self.kinetics
    }
}

impl<L: EvaporationLaw + ThermalLaw + ?Sized> ThermalLaw for KineticEvaporation<'_, L> {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.source.storage_basis()
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        self.source.caloric_identity()
    }
    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        self.source.validate_cell(cell)
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.source.temperature_k(h)
    }
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        self.source.bulk_density_kg_m3(h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.source.maximum_temperature_slope()
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h_j_kg: f64,
        final_h_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.source
            .specific_entropy_change_j_kg_k(initial_h_j_kg, final_h_j_kg)
    }
    fn advected_composition(&self) -> Option<CaloricComposition<'_>> {
        self.source.advected_composition()
    }
    fn advected_entropy_anchor(&self) -> Option<CaloricEntropyAnchor<'_>> {
        self.source.advected_entropy_anchor()
    }
}

fn sealed() -> SurfaceFlux {
    SurfaceFlux {
        mass_flux_kg_m2_s: 0.0,
        mass_flux_derivative_kg_m2_s_k: 0.0,
        latent_heat_flux_w_m2: 0.0,
        outward_normal_velocity_m_s: 0.0,
    }
}

// Same positive product with a fallback ordering for intermediate overflow or
// underflow. Every accepted nonzero flux remains finite and representable.
fn product(a: f64, b: f64, c: f64) -> f64 {
    let direct = (a * b) * c;
    if positive(direct) {
        return direct;
    }
    let alternate = (a * c) * b;
    if positive(alternate) {
        alternate
    } else {
        (b * c) * a
    }
}

fn harmonic(a: f64, b: f64) -> f64 {
    let smaller = a.min(b);
    let larger = a.max(b);
    smaller / (1.0 + smaller / larger)
}

// p_sat-p_i is retained directly. Reconstructing it by subtracting two rounded
// absolute pressures would lose the kinetic resistance near equilibrium.
fn finish(
    p: EvaporationProperties,
    t: f64,
    scale: f64,
    conductance: f64,
    carrier_interface: f64,
    kinetic_pressure_drop: f64,
    mass: f64,
) -> Result<SurfaceFlux, EvaporationError> {
    let film_tangent = scale / carrier_interface;
    let effective = harmonic(conductance, film_tangent);
    if !positive(effective) {
        return Err(EvaporationError::Unrepresentable(
            "positive combined interface/film conductance underflows",
        ));
    }
    // Implicit differentiation of j_k=j_f includes K'(T_s)=-K/(2 T_s).
    let derivative = checked(
        effective * (p.saturation_pressure_derivative_pa_k - 0.5 * (kinetic_pressure_drop / t)),
    )?;
    Ok(SurfaceFlux {
        mass_flux_kg_m2_s: checked(mass)?,
        mass_flux_derivative_kg_m2_s_k: derivative,
        latent_heat_flux_w_m2: checked(mass * p.latent_heat_j_kg)?,
        outward_normal_velocity_m_s: checked(-mass / p.liquid_density_kg_m3)?,
    })
}

pub(super) fn flux(
    p: EvaporationProperties,
    t: f64,
    film: GasFilm,
    kinetics: InterfaceKinetics,
) -> Result<SurfaceFlux, EvaporationError> {
    kinetics.validate()?;
    film.validate()?;
    if !positive(t) {
        return Err(EvaporationError::Invalid(
            "positive finite interface temperature required",
        ));
    }
    if p.saturation_pressure_pa >= film.total_pressure_pa {
        return Err(EvaporationError::Boiling {
            saturation_pressure_pa: p.saturation_pressure_pa,
            total_pressure_pa: film.total_pressure_pa,
        });
    }
    let factor = match kinetics {
        InterfaceKinetics::Equilibrium => return super::flux(p, film),
        InterfaceKinetics::HertzKnudsen {
            accommodation_coefficient: alpha,
        } => alpha,
        InterfaceKinetics::Schrage {
            accommodation_coefficient: alpha,
        } => alpha * (2.0 / (2.0 - alpha)),
    };
    if factor == 0.0 || film.mass_transfer_coefficient_m_s == 0.0 {
        return Ok(sealed());
    }
    let conductance = checked(
        factor
            * (fs_math::det::sqrt(p.molar_mass_kg_mol)
                / fs_math::det::sqrt(2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT))
            / fs_math::det::sqrt(t),
    )?;
    let scale = checked(
        film.mass_transfer_coefficient_m_s
            * (p.molar_mass_kg_mol / MOLAR_GAS_CONSTANT)
            * (film.total_pressure_pa / film.temperature_k),
    )?;
    if conductance == 0.0 || scale == 0.0 {
        return Err(EvaporationError::Unrepresentable(
            "positive interface or gas-film conductance underflows",
        ));
    }
    let carrier_far = film.total_pressure_pa - film.vapor_partial_pressure_pa;
    let carrier_surface = film.total_pressure_pa - p.saturation_pressure_pa;
    let gap = p.saturation_pressure_pa - film.vapor_partial_pressure_pa;
    if gap == 0.0 {
        return finish(p, t, scale, conductance, carrier_far, 0.0, 0.0);
    }
    let equilibrium_mass = checked(scale * log_carrier_ratio(carrier_far, carrier_surface, gap))?;
    if equilibrium_mass == 0.0 {
        return Err(EvaporationError::Unrepresentable(
            "nonzero vapor-pressure drive underflows to zero film flux",
        ));
    }
    let film_secant = (equilibrium_mass / gap).abs();
    let kinetic_limited = conductance <= film_secant;
    let ratio = conductance.min(film_secant) / conductance.max(film_secant);
    let mut fraction = ratio / (1.0 + ratio);
    if fraction == 0.0 {
        // The smaller pressure drop is below the f64 range; its correction to
        // the represented limiting flux is consequently below roundoff.
        let (carrier, drop, mass) = if kinetic_limited {
            (carrier_far, gap, checked(conductance * gap)?)
        } else {
            (carrier_surface, 0.0, equilibrium_mass)
        };
        if mass == 0.0 {
            return Err(EvaporationError::Unrepresentable(
                "nonzero interfacial transfer underflows",
            ));
        }
        return finish(p, t, scale, conductance, carrier, drop, mass);
    }
    let (mut lower, mut upper) = (0.0, 1.0);
    let magnitude = gap.abs();
    // Solve the smaller resistance's pressure-drop fraction. This retains a
    // tiny film drop for slow kinetics and a tiny kinetic drop for a slow film.
    // The monotone scalar equation has a unique root; 128 safeguarded trials
    // bound work without allocations or an unbounded convergence loop.
    for _ in 0..128 {
        let (film_fraction, kinetic_fraction) = if kinetic_limited {
            (fraction, 1.0 - fraction)
        } else {
            (1.0 - fraction, fraction)
        };
        let carrier = kinetic_fraction * carrier_far + film_fraction * carrier_surface;
        let film_mass =
            checked(scale * log_carrier_ratio(carrier_far, carrier, gap * film_fraction))?.abs();
        let kinetic_mass = product(conductance, magnitude, kinetic_fraction);
        if kinetic_mass.is_finite() {
            let residual = if kinetic_limited {
                film_mass - kinetic_mass
            } else {
                kinetic_mass - film_mass
            };
            let tolerance = 32.0 * f64::EPSILON * film_mass.max(kinetic_mass);
            if residual.abs() <= tolerance {
                let mass = if kinetic_limited {
                    kinetic_mass
                } else {
                    film_mass
                };
                if mass == 0.0 {
                    return Err(EvaporationError::Unrepresentable(
                        "nonzero interfacial transfer underflows",
                    ));
                }
                return finish(
                    p,
                    t,
                    scale,
                    conductance,
                    carrier,
                    gap * kinetic_fraction,
                    mass.copysign(gap),
                );
            }
            if residual < 0.0 {
                lower = fraction;
            } else {
                upper = fraction;
            }
            let tangent = scale / carrier;
            let largest = conductance.max(tangent);
            let change =
                ((residual / largest) / magnitude) / (1.0 + conductance.min(tangent) / largest);
            let newton = fraction - change;
            if newton.is_finite() && newton > lower && newton < upper {
                fraction = newton;
                continue;
            }
        } else {
            // An overflowing kinetic stream brackets the finite joint flux.
            // Move toward a smaller kinetic drop in either parameterization.
            if kinetic_limited {
                lower = fraction;
            } else {
                upper = fraction;
            }
        }
        let middle = lower + 0.5 * (upper - lower);
        if middle == lower || middle == upper {
            return Err(EvaporationError::NonConverged);
        }
        fraction = middle;
    }
    Err(EvaporationError::NonConverged)
}

#[cfg(test)]
mod tests;
