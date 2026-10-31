//! Explicit instantaneous-equilibrium acoustic view of the original fluid state.
//!
//! Phases share pressure, temperature, chemical potential and velocity during
//! the perturbation. This homogeneous equilibrium limit allows latent mass/heat
//! transfer; it is distinct from the frozen-phase Wood response. It does not
//! determine the frequencies or spatial scales at which these assumptions hold.
//! Surface energy, slip, finite-rate relaxation and attenuation are excluded.

use super::{E, PureFluidEquilibriumPhase, PureFluidEquilibriumState};
use crate::MaterialError;
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};

/// Source-derived Jacobian for pressure and mass-specific entropy in (T,rho).
/// The response is equilibrium, not frozen quality, and supplies no mixed Cp.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumThermodynamicDerivatives {
    /// Complete equilibrium `(du/dT)_rho` [J/(kg K)].
    pub isochoric_heat_capacity_j_kg_k: f64,
    /// `(dp/dT)_rho` [Pa/K].
    pub pressure_temperature_at_density_pa_k: f64,
    /// `(dp/drho)_T` [m²/s²], zero inside the two-phase dome.
    pub pressure_density_at_temperature_m2_s2: f64,
    /// `(ds/dT)_rho` [J/(kg K²)].
    pub entropy_temperature_at_density_j_kg_k2: f64,
    /// `(ds/drho)_T` [J m³/(kg² K)].
    pub entropy_density_at_temperature_j_m3_kg2_k: f64,
    /// `(dT/drho)_s` [K m³/kg].
    pub temperature_density_at_entropy_k_m3_kg: f64,
    /// `(dp/drho)_s=c²` [m²/s²].
    pub pressure_density_at_entropy_m2_s2: f64,
    /// `(dp/du)_rho` for specific internal energy [kg/m³].
    pub pressure_energy_at_density_kg_m3: f64,
    /// `(dp/drho)_u` for specific internal energy [m²/s²].
    pub pressure_density_at_energy_m2_s2: f64,
}

/// Borrowed original state plus its explicitly chosen homogeneous-equilibrium
/// acoustic response. The base caloric state does not silently select an
/// acoustic relaxation policy; construct this view to request that limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HomogeneousEquilibriumAcoustics<'a> {
    source: &'a PureFluidEquilibriumState,
    derivatives: EquilibriumThermodynamicDerivatives,
    acoustic: FluidAcousticState,
}
impl<'a> HomogeneousEquilibriumAcoustics<'a> {
    /// Derive response from one immutable, already admitted phase equilibrium.
    /// Exact or unresolvable saturated endpoints refuse a bidirectional wave
    /// law; homogeneous states outside the dome use original root mechanics.
    pub fn try_new(source: &'a PureFluidEquilibriumState) -> Result<Self, E> {
        let properties = source.reference().properties();
        let mass = properties.molar_mass().value();
        let ideal_cv = properties.cv().value();
        let (cv, pt, pv, sound_speed) = match source.phase() {
            PureFluidEquilibriumPhase::Coexistence(mixture) => {
                let response = mixture.equilibrium_response(ideal_cv).map_err(E::Phase)?;
                (
                    response.molar_cv_j_mol_k,
                    response.dp_dt_v_pa_k,
                    response.dp_dv_t_pa_mol_m3,
                    response.sound_speed_m_s,
                )
            }
            PureFluidEquilibriumPhase::Liquid(root) | PureFluidEquilibriumPhase::Vapor(root) => {
                let response = root.caloric_mechanics(ideal_cv, mass).ok_or(E::Invalid(
                    "original homogeneous acoustic response is unstable or unrepresentable",
                ))?;
                (
                    response.molar_cv_j_mol_k,
                    response.dp_dt_v_pa_k,
                    response.dp_dv_t_pa_mol_m3,
                    response.sound_speed_m_s,
                )
            }
        };
        let density = source.density_kg_m3();
        let temperature = source.temperature_k();
        let specific_cv = cv / mass;
        let entropy_temperature = specific_cv / temperature;
        let entropy_density = -(pt / density) / density;
        let pressure_density = -(pv * (mass / density)) / density;
        let temperature_density = -entropy_density / entropy_temperature;
        let pressure_energy = pt / specific_cv;
        let pressure_density_energy = pressure_density
            + pressure_energy * ((temperature * pt - source.pressure_pa()) / density) / density;
        let acoustic = FluidAcousticState::try_new(density, sound_speed)?;
        let sound_squared = sound_speed * sound_speed;
        if [
            specific_cv,
            entropy_temperature,
            entropy_density,
            pressure_density,
            temperature_density,
            sound_squared,
            pressure_energy,
            pressure_density_energy,
        ]
        .into_iter()
        .any(|x| !x.is_finite())
            || specific_cv <= 0.0
            || entropy_temperature <= 0.0
            || pressure_density < 0.0
        {
            return Err(E::Invalid(
                "equilibrium pressure/entropy derivatives are unrepresentable",
            ));
        }
        Ok(Self {
            source,
            acoustic,
            derivatives: EquilibriumThermodynamicDerivatives {
                isochoric_heat_capacity_j_kg_k: specific_cv,
                pressure_temperature_at_density_pa_k: pt,
                pressure_density_at_temperature_m2_s2: pressure_density,
                entropy_temperature_at_density_j_kg_k2: entropy_temperature,
                entropy_density_at_temperature_j_m3_kg2_k: entropy_density,
                temperature_density_at_entropy_k_m3_kg: temperature_density,
                pressure_density_at_entropy_m2_s2: sound_squared,
                pressure_energy_at_density_kg_m3: pressure_energy,
                pressure_density_at_energy_m2_s2: pressure_density_energy,
            },
        })
    }
    /// The original state and complete thermochemical source receipt.
    #[must_use]
    pub const fn source(&self) -> &'a PureFluidEquilibriumState {
        self.source
    }
    /// Source-derived pressure/entropy derivatives in physical (T,rho) coordinates.
    #[must_use]
    pub const fn derivatives(&self) -> EquilibriumThermodynamicDerivatives {
        self.derivatives
    }
    /// Coherent equilibrium density, sound speed, bulk modulus and impedance.
    #[must_use]
    pub const fn state(&self) -> FluidAcousticState {
        self.acoustic
    }
}

impl FluidAcoustics for HomogeneousEquilibriumAcoustics<'_> {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.acoustic)
    }
}

impl PureFluidEquilibriumState {
    /// Explicitly request instantaneous phase/thermal equilibrium and common
    /// phase velocity for a small-signal acoustic perturbation.
    pub fn equilibrium_acoustics(&self) -> Result<HomogeneousEquilibriumAcoustics<'_>, E> {
        HomogeneousEquilibriumAcoustics::try_new(self)
    }
}
