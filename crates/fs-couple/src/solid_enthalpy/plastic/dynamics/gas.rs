//! Finite source gas and the original inertial thermoplastic trajectory.
//!
//! Every gas/enthalpy proposal restarts the native consistent-mass solid solve
//! from the SAME accepted displacement, velocity and constitutive history.
//! Actual plastic/viscous and reversible heat enters the shared caloric
//! residual once. Kinetic storage and its distinct numerical loss enter the
//! enclosing work balance; neither numerical loss is deposited as heat.
//!
//! Gas mass, volume and composition, solid reference geometry and thermal
//! areas remain fixed. This thermal connection supplies no gas pressure
//! traction or moving-volume work. The original source/phase/small-strain
//! admission and physical-timestep accuracy limits still apply.

use super::*;
use crate::circuit_enthalpy::contacts::heat::{HeatNetworkSample, PreparedHeatNetwork};
use crate::circuit_enthalpy::material::{MaterialNetworkSample, PreparedMaterialNetwork};
use crate::solid_enthalpy::gas::{SolidGasContact, SolidGasError, SolidGasSolve};
use crate::solid_enthalpy::plastic::gas::{ThermoplasticGasCandidate, ThermoplasticGasReport};
use fs_material::thermochemical_gas::energy::IsochoricGasParcel;

/// Accept finite gas, thermal inventory, motion and native history together.
pub struct DynamicGasCandidate<'m, 't, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    gas: IsochoricGasParcel,
    dynamics: ThermoplasticDynamics<'m>,
    report: ThermoplasticDynamicsReport<ThermoplasticGasReport<S>>,
}

impl<'m, 't, S> DynamicGasCandidate<'m, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn gas(&self) -> &IsochoricGasParcel {
        &self.gas
    }
    pub const fn dynamics(&self) -> &ThermoplasticDynamics<'m> {
        &self.dynamics
    }
    pub const fn report(&self) -> &ThermoplasticDynamicsReport<ThermoplasticGasReport<S>> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        IsochoricGasParcel,
        ThermoplasticDynamics<'m>,
        ThermoplasticDynamicsReport<ThermoplasticGasReport<S>>,
    ) {
        (self.thermal, self.gas, self.dynamics, self.report)
    }
}

fn finish<'m, 't, S>(
    candidate: ThermoplasticGasCandidate<'m, 't, S>,
    motion: Option<Motion>,
    cx: &Cx<'_>,
) -> Result<DynamicGasCandidate<'m, 't, S>, SolidGasError> {
    let motion = motion.ok_or(SolidGasError::Invalid(
        "inertial gas/solid candidate must retain its native kinetic account",
    ))?;
    let (thermal, gas, solid, material) = candidate.into_parts();
    poll::<Infallible>(cx)?;
    Ok(DynamicGasCandidate {
        thermal,
        gas,
        dynamics: ThermoplasticDynamics {
            solid,
            velocity_mps: motion.velocity_mps,
        },
        report: ThermoplasticDynamicsReport {
            material,
            inertia: motion.energy,
        },
    })
}

impl<'m> ThermoplasticDynamics<'m> {
    /// Solve finite-gas heat, original thermal cells and actual nodal motion
    /// at one shared endpoint. Source density supplies consistent reference
    /// mass; source moduli, yield, hardening and selected overstress kinetics
    /// resolve at each proposed temperature and the actual physical interval.
    /// Endpoint forces and full support reactions supply mechanical work.
    ///
    /// Guesses are sorted CellId energy increments followed by total gas heat,
    /// in joules. The gas controls gate each original account and the complete
    /// balance, including kinetic storage and both numerical losses. The
    /// material report keeps constitutive loss separate from the inertia report.
    /// Every error or cancellation preserves all accepted inputs. Electrical
    /// bindings refuse because this entry point advances no circuit.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_isochoric_gas<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        transport: Option<&PreparedMaterialNetwork<'_, '_>>,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        yield_tolerance_pa: f64,
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<DynamicGasCandidate<'m, 't>, SolidGasError> {
        let (candidate, motion) = self.solid.advance_gas_transport(
            thermal,
            gas,
            target_s,
            external_gas_heat_j,
            contacts,
            transport,
            prescribed_m,
            nodal_forces_n,
            settings,
            yield_tolerance_pa,
            controls,
            initial_energy_increments_j,
            Some(&self.velocity_mps),
            cx,
        )?;
        finish(candidate, motion, cx)
    }

    /// Add the complete prepared heat assembly to the same gas/motion solve.
    /// Original contact resistance, integrated layers and closed radiation
    /// participate in the shared residual and retain their complete accepted
    /// source sample. Internal transfers are never repeated as external heat.
    /// All interval, source-admission and atomicity rules of
    /// the bulk gas entry point apply.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_isochoric_gas_heat_network<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        transport: &PreparedHeatNetwork<'_, '_>,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        yield_tolerance_pa: f64,
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<DynamicGasCandidate<'m, 't, HeatNetworkSample>, SolidGasError> {
        let (candidate, motion) = self.solid.advance_gas_transport(
            thermal,
            gas,
            target_s,
            external_gas_heat_j,
            contacts,
            transport,
            prescribed_m,
            nodal_forces_n,
            settings,
            yield_tolerance_pa,
            controls,
            initial_energy_increments_j,
            Some(&self.velocity_mps),
            cx,
        )?;
        finish(candidate, motion, cx)
    }
}
