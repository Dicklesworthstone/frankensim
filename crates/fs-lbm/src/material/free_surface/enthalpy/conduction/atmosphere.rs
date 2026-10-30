//! Finite atmospheric internal energy coupled to the original moving-liquid
//! thermal timeline. Gas composition, EOS and energy inversion keep their
//! existing source owner; film geometry remains an explicit physical input.

use super::*;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};

mod body;
pub(super) use body::BodyAtmosphere;
pub use body::{FreeSurfaceBodyGasExchange, FreeSurfaceBodyGasFace, FreeSurfaceBodyGasSample};

/// Finite gas account over one accepted native VOF/enthalpy interval.
///
/// Gas-supplied heat in `environment` or a solved `surface` appears in the liquid
/// subsystem's external heat, but is INTERNAL to liquid + walls + atmosphere.
/// Molecular heat retains its original account; it is not gas work or an
/// additional external heat source. No full VOF kinetic-energy claim is added.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceFiniteGasExchange {
    /// Actual mass times the change of the gas's sourced specific INTERNAL energy.
    pub internal_energy_change_j: f64,
    /// Actual gas storage change plus supplied liquid AND solid heat.
    pub gas_balance_defect_j: f64,
    /// Combined actual liquid/wall/gas storage minus authored, molecular,
    /// electrical and radiant heat, according to their original accounts.
    pub combined_thermal_balance_defect_j: f64,
    /// Accumulated absolute allowance from the original gas energy inversion.
    pub gas_energy_tolerance_j: f64,
    pub substeps: usize,
}

#[derive(Clone)]
pub(super) struct FiniteAtmosphere {
    pub(super) parcel: IsochoricGasParcel,
    pub(super) temperature_slope: f64,
    controls: GasEnergySolve,
    initial_specific_energy_j_kg: f64,
    substeps: usize,
}

impl FiniteAtmosphere {
    pub(super) fn try_new(
        parcel: &IsochoricGasParcel,
        controls: GasEnergySolve,
    ) -> Result<Self, Error> {
        let temperature_slope = parcel
            .model()
            .maximum_isochoric_temperature_slope(controls)?;
        let temperature = parcel.state().state().temperature;
        if temperature < controls.minimum_temperature_k
            || temperature > controls.maximum_temperature_k
        {
            return Err(Error::Invalid(
                "finite atmosphere starts outside its caloric bracket",
            ));
        }
        Ok(Self {
            parcel: parcel.clone(),
            temperature_slope,
            controls,
            initial_specific_energy_j_kg: parcel.state().specific_internal_energy_j_kg(),
            substeps: 0,
        })
    }

    pub(super) fn accept_heat(
        &mut self,
        heat_j: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        let old = self.parcel.state().clone();
        let mut checkpoint_error = None;
        let result = self.parcel.add_heat(heat_j, self.controls, || {
            if let Err(error) = check() {
                checkpoint_error = Some(error);
                true
            } else {
                false
            }
        });
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        result?;
        if !environment::same_sources(&old, self.parcel.state()) {
            return Err(Error::Invalid(
                "finite atmosphere left its original transport source claims",
            ));
        }
        self.substeps = self
            .substeps
            .checked_add(1)
            .ok_or(Error::Invalid("finite atmosphere substep count overflow"))?;
        Ok(())
    }

    pub(super) fn balance(
        &self,
        heat_to_condensed_j: f64,
        thermal: &EnergyAudit,
        thermal_tolerance_j: f64,
    ) -> Result<FreeSurfaceFiniteGasExchange, Error> {
        let stored = environment::finite(
            self.parcel.mass_kg()
                * (self.parcel.state().specific_internal_energy_j_kg()
                    - self.initial_specific_energy_j_kg),
        )?;
        let defect = environment::finite(stored + heat_to_condensed_j)?;
        let allowance = environment::finite(
            self.substeps as f64 * self.parcel.mass_kg() * self.controls.energy_tolerance_j_kg,
        )?;
        if defect.abs() > allowance {
            return Err(Error::Invalid(
                "finite atmosphere exceeded its original gas energy budget",
            ));
        }
        let combined = environment::finite(
            thermal.stored_energy_change_j + stored
                - thermal.incoming_energy_j
                - (thermal.external_heat_j - heat_to_condensed_j),
        )?;
        if combined.abs() > environment::finite(thermal_tolerance_j + allowance)? {
            return Err(Error::Invalid(
                "finite liquid/wall/atmosphere thermal exchange does not close",
            ));
        }
        Ok(FreeSurfaceFiniteGasExchange {
            internal_energy_change_j: stored,
            gas_balance_defect_j: defect,
            combined_thermal_balance_defect_j: combined,
            gas_energy_tolerance_j: allowance,
            substeps: self.substeps,
        })
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Exchange heat with a finite source-backed atmosphere, finite walls and
    /// the SAME native liquid conduction, rheology, momentum and mass transport.
    /// Walls with `with_gas_faces` exchange directly with this SAME gas parcel.
    /// Every solid has its own sourced conductivity and caloric capacity; one
    /// aggregate opposite gas debit includes solid and liquid films together.
    /// Attached wall radiation supplies external heat to its original finite
    /// solid; it is never included in that opposite gas debit.
    ///
    /// The well-mixed gas has fixed mass, composition and volume. Each actual
    /// thermal substep deposits aggregate opposite liquid/solid joules into its original NASA-9
    /// internal-energy owner. Updated gas T and EOS pressure re-evaluate the
    /// original film transport sources; source switches and unsupported states
    /// refuse. The native thermal CFL includes solid and liquid contact losses and
    /// the finite gas capacity, so small atmospheres cannot act as infinite baths.
    ///
    /// Provide currently exposed patches at every call. Geometry is held during
    /// the thermal stage, then the existing native VOF/enthalpy split advances
    /// once. Gas volume is an explicitly declared fixed-volume approximation;
    /// lattice compressibility never supplies a physical gas-volume change.
    /// This supplies no gas
    /// momentum, pressure traction, evaporation/species transfer, reaction,
    /// participating-gas radiation, moving solid or phase-density handoff. Film thickness is
    /// explicit; neither a gas name nor fill fraction invents convection.
    ///
    /// `interfaces` retain existing ordered liquid/wall contact laws. Optional
    /// molecular heating uses the original physical stress-power owner once.
    /// No finite gas, wall or liquid state publishes until all thermal/source
    /// and native flow/advection gates pass. Refine time and space for the
    /// existing first-order thermal/flow and frozen-film approximation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_finite_gas_environment(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        bodies: &mut [FreeSurfaceThermalBody<'_>],
        contacts: &[FreeSurfaceBodyContact],
        interfaces: &[FreeSurfaceBodyInterface<'_>],
        gas: &mut IsochoricGasParcel,
        exposed: &[FreeSurfaceGasFace],
        gas_controls: GasEnergySolve,
        with_viscous_heating: bool,
        gate: &CancelGate,
    ) -> Result<viscous::FreeSurfaceViscousConductionReport, Error> {
        poll(gate)?;
        if (exposed.is_empty() && bodies.iter().all(|body| body.gas_faces().is_empty()))
            || exposed.len() > controls.maximum_faces
        {
            return Err(Error::Invalid(
                "finite atmosphere requires bounded exposed liquid or solid patches",
            ));
        }
        let mut candidate = FiniteAtmosphere::try_new(gas, gas_controls)?;
        let environment: Vec<_> = exposed
            .iter()
            .map(|&face| FreeSurfaceGasContact {
                face,
                model: gas.model(),
                context: gas.point(),
            })
            .collect();
        let applied = if with_viscous_heating {
            self.viscous_sources(gate)?
        } else {
            Vec::new()
        };
        let conduction = self.advance_with_interfaces_and_environment(
            faces,
            heat,
            controls,
            &applied,
            bodies,
            contacts,
            interfaces,
            &environment,
            Some(&mut candidate),
            gate,
        )?;
        // The shared native transaction accepted every participant already.
        // Publication has no further query, solve or fallible operation.
        *gas = candidate.parcel;
        Ok(viscous::FreeSurfaceViscousConductionReport {
            conduction,
            applied,
        })
    }
}
