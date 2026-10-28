//! Frozen-composition material closure of the existing NSF solver.
//!
//! The original NASA9 and transport owners evaluate every actual cell/face
//! state. No table or constant gamma replaces their caloric inversion.

use super::{MAX_VERIFIED_LATTICE_VISCOSITY, NsfError, NsfLaw, NsfSpec};
use fs_matdb::QueryPoint;
use fs_material::gas::{GasSpec, GasState};
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::{ResolvedThermochemicalGasState, ThermochemicalGasModel};

/// Original source laws and typed conditions of one frozen-composition gas.
///
/// Temperature and pressure evolve; all other query axes and species fractions
/// stay fixed. This does not provide chemistry, species transport or phase change.
#[derive(Clone, Debug)]
pub struct ThermochemicalNsfGas {
    model: ThermochemicalGasModel,
    context: QueryPoint,
    energy_solve: GasEnergySolve,
    reference: GasState,
    energy_reference_j_kg: f64,
}

impl ThermochemicalNsfGas {
    /// Bind the source owner and its bounded existing energy inversion.
    ///
    /// The conserved chart is u_s(T) = u_NASA9(T) - u_offset, with
    /// u_offset = u_NASA9(T_ref) - cv(T_ref) T_ref. Its value at the reference
    /// state equals the legacy cv*T energy, and arbitrary formation-energy
    /// offsets cannot enter the separate lattice mass / finite-volume energy
    /// fluxes. Negative NASA9 absolute energies are admitted.
    ///
    /// # Errors
    /// Original material support/quantity failures or invalid inversion controls.
    pub fn new(
        model: ThermochemicalGasModel,
        context: QueryPoint,
        energy_solve: GasEnergySolve,
    ) -> Result<Self, NsfError> {
        energy_solve.validate().map_err(gas_error)?;
        let resolved = model.evaluate(&context).map_err(gas_error)?;
        let reference = *resolved.state();
        let cv = reference.specific_heat_cp - reference.specific_gas_constant;
        let energy_reference_j_kg =
            resolved.specific_internal_energy_j_kg() - cv * reference.temperature;
        if !energy_reference_j_kg.is_finite() {
            return Err(NsfError::Input("unrepresentable gas energy reference"));
        }
        Ok(Self {
            model,
            context,
            energy_solve,
            reference,
            energy_reference_j_kg,
        })
    }

    /// Original NASA9 mixture and source-selected transport materials.
    #[must_use]
    pub const fn model(&self) -> &ThermochemicalGasModel {
        &self.model
    }

    /// Reference typed conditions; extra material coordinates remain fixed.
    #[must_use]
    pub const fn context(&self) -> &QueryPoint {
        &self.context
    }

    /// Explicit work, residual and temperature bounds of caloric inversion.
    #[must_use]
    pub const fn energy_solve(&self) -> GasEnergySolve {
        self.energy_solve
    }

    /// Fixed value subtracted from NASA9 specific internal energy [J/kg].
    ///
    /// Grid energies and manufactured energy sources use this sensible chart.
    /// A source carrying mass must also subtract this value times its mass
    /// supply when converting an absolute NASA9 energy supply to this chart.
    #[must_use]
    pub const fn energy_reference_j_kg(&self) -> f64 {
        self.energy_reference_j_kg
    }

    fn at(
        &self,
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<ResolvedThermochemicalGasState, NsfError> {
        let point = self
            .context
            .clone()
            .with("temperature", temperature_k)
            .and_then(|p| p.with("pressure", pressure_pa))
            .map_err(gas_error)?;
        self.model.evaluate(&point).map_err(gas_error)
    }

    fn cell(
        &self,
        law: &NsfLaw,
        density_lattice: f64,
        resolved: &ResolvedThermochemicalGasState,
    ) -> Result<NsfCellMaterial, NsfError> {
        let state = *resolved.state();
        let mu_lattice = state.dynamic_viscosity * law.transport_scale();
        // The relaxation rate depends on local kinematic viscosity, including
        // pressure/density changes, not only on the reference transport table.
        let nu = mu_lattice / density_lattice;
        if !nu.is_finite() || nu > MAX_VERIFIED_LATTICE_VISCOSITY {
            return Err(NsfError::Viscosity {
                observed: nu,
                maximum: MAX_VERIFIED_LATTICE_VISCOSITY,
            });
        }
        let sensible_energy_lattice = (resolved.specific_internal_energy_j_kg()
            - self.energy_reference_j_kg)
            / law.specific_energy_unit_j_kg();
        if !sensible_energy_lattice.is_finite() {
            return Err(NsfError::Input("unrepresentable sensible gas energy"));
        }
        Ok(NsfCellMaterial {
            state,
            mu_lattice,
            sensible_energy_lattice,
        })
    }

    pub(super) fn reconstruct(
        &self,
        law: &NsfLaw,
        density_lattice: f64,
        sensible_energy_lattice: f64,
    ) -> Result<(f64, NsfCellMaterial), NsfError> {
        let solved = self
            .model
            .reconstruct_density_energy_change(
                &self.context,
                density_lattice * law.rho_ref_kg_m3,
                self.energy_reference_j_kg,
                sensible_energy_lattice * law.specific_energy_unit_j_kg(),
                self.energy_solve,
                || false,
            )
            .map_err(gas_error)?;
        let material = self.cell(law, density_lattice, solved.state())?;
        Ok((law.theta_of_kelvin(material.state.temperature), material))
    }
}

/// Gas law accepted by the same NSF grid, without another solver type.
#[derive(Clone, Debug)]
pub enum NsfGasMedium {
    /// Existing constant-cv GasSpec path, with its original arithmetic.
    CaloricallyPerfect(GasSpec),
    /// Source-owned NASA9 caloric properties and local transport.
    Thermochemical(ThermochemicalNsfGas),
}

impl From<GasSpec> for NsfGasMedium {
    fn from(value: GasSpec) -> Self {
        Self::CaloricallyPerfect(value)
    }
}
impl From<ThermochemicalNsfGas> for NsfGasMedium {
    fn from(value: ThermochemicalNsfGas) -> Self {
        Self::Thermochemical(value)
    }
}

/// Small numerical cache; the original source cards stay in the medium owner.
#[derive(Clone, Copy, Debug)]
pub(super) struct NsfCellMaterial {
    pub(super) state: GasState,
    pub(super) mu_lattice: f64,
    pub(super) sensible_energy_lattice: f64,
}

fn gas_error(error: impl core::fmt::Display) -> NsfError {
    NsfError::Gas(error.to_string())
}

impl NsfGasMedium {
    pub(super) fn reference(&self, spec: &NsfSpec<Self, String>) -> Result<GasState, NsfError> {
        if let Self::Thermochemical(source) = self {
            let controls = source.energy_solve;
            if source.reference.temperature != spec.temperature_k
                || source.reference.pressure != spec.pressure_pa
                || controls.minimum_temperature_k != spec.temperature_window_k[0]
                || controls.maximum_temperature_k != spec.temperature_window_k[1]
                || spec.temperature_k < spec.temperature_window_k[0]
                || spec.temperature_k > spec.temperature_window_k[1]
            {
                return Err(NsfError::Input(
                    "source reference T/p and energy bracket must match the NSF reference and run window",
                ));
            }
        }
        self.at(spec.temperature_k, spec.pressure_pa)
    }

    pub(super) fn at(&self, temperature_k: f64, pressure_pa: f64) -> Result<GasState, NsfError> {
        match self {
            Self::CaloricallyPerfect(gas) => {
                GasState::try_new(gas, temperature_k, pressure_pa).map_err(gas_error)
            }
            Self::Thermochemical(source) => Ok(*source.at(temperature_k, pressure_pa)?.state()),
        }
    }

    pub(super) fn resolve(
        &self,
        law: &NsfLaw,
        density_lattice: f64,
        theta: f64,
    ) -> Result<Option<NsfCellMaterial>, NsfError> {
        match self {
            Self::CaloricallyPerfect(_) => {
                // A dilute cell has a larger relaxation time even though the
                // gas's dynamic viscosity is pressure-independent. The
                // reference-density table check cannot admit rarefied cells
                // or reservoir states on its own.
                let nu = law.mu_lattice(theta) / density_lattice;
                if !nu.is_finite() || nu <= 0.0 || nu > MAX_VERIFIED_LATTICE_VISCOSITY {
                    return Err(NsfError::Viscosity {
                        observed: nu,
                        maximum: MAX_VERIFIED_LATTICE_VISCOSITY,
                    });
                }
                Ok(None)
            }
            Self::Thermochemical(source) => {
                let temperature = law.kelvin_of_theta(theta);
                let pressure = law.pressure_pa(density_lattice, theta);
                let resolved = source.at(temperature, pressure)?;
                source.cell(law, density_lattice, &resolved).map(Some)
            }
        }
    }

    pub(super) fn source_state(
        &self,
        law: &NsfLaw,
        density_lattice: f64,
        theta: f64,
    ) -> Result<Option<ResolvedThermochemicalGasState>, NsfError> {
        match self {
            Self::CaloricallyPerfect(_) => Ok(None),
            Self::Thermochemical(source) => {
                let temperature = law.kelvin_of_theta(theta);
                let pressure = law.pressure_pa(density_lattice, theta);
                source.at(temperature, pressure).map(Some)
            }
        }
    }

    pub(super) fn face_conductivity(
        &self,
        law: &NsfLaw,
        theta: f64,
        pressure_pa: f64,
    ) -> Result<f64, NsfError> {
        match self {
            Self::CaloricallyPerfect(_) => Ok(law.kappa_lattice(theta)),
            Self::Thermochemical(source) => {
                let state = source.at(law.kelvin_of_theta(theta), pressure_pa)?;
                Ok(state.state().thermal_conductivity * law.transport_scale() / law.gas_constant)
            }
        }
    }
}

impl NsfLaw {
    pub(super) fn pressure_pa(&self, density_lattice: f64, theta: f64) -> f64 {
        density_lattice * self.rho_ref_kg_m3 * self.gas_constant * self.kelvin_of_theta(theta)
    }

    /// SI specific-energy unit [(m/s)^2 = J/kg] of the lattice chart.
    #[must_use]
    pub fn specific_energy_unit_j_kg(&self) -> f64 {
        (self.dx_m / self.dt_s).powi(2)
    }

    fn transport_scale(&self) -> f64 {
        self.dt_s / (self.rho_ref_kg_m3 * self.dx_m * self.dx_m)
    }
}
