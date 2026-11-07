//! Independently sealed gas cavities at supplied thermal porous endpoints.
//!
//! Material temperature/history, caloric heat and optional original solid inertia
//! remain with the thermal porous owner. The same native cavity batch supplies
//! current force and sparse/geometric plus acoustic rank-one tangents. Gas work
//! uses the shared dry/thermal interval account; no gas caloric law is copied.

use fs_exec::Cx;
use fs_material::plastic::finite::porous::thermal::{
    AffineThermalFinitePorous, FiniteThermalPorousMaterial,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

use super::{
    FiniteThermalPorousTetError, FiniteThermalPorousTetSolid, FiniteThermalPorousTetStep,
    FiniteThermoplasticControls,
};
use crate::hyper3::HyperTetSettings;
use crate::hyper3::plastic::porous::cavity::FinitePorousCavityEnergy;
use crate::hyper3::pressure::cavity::{AdiabaticCavity, AdiabaticCavityGas, TetCavityGasLoad};
use crate::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};

/// Complete accepted thermal/material endpoint and original-source gas states.
///
/// The whole-system work identity is dead + support + held-pressure work +
/// summed exterior work = material internal-energy change + ORIGINAL physical
/// caloric heat + summed gas storage changes + kinetic change + numerical kinetic
/// loss + `external_work_defect_j` + `step.energy.balance_error_j`.
/// Caloric heat already includes physical yield loss and reversible exchange;
/// gas storage, pressure work and numerical defects must not be added to it.
///
/// Continue through `step.solid.advance_with_cavities` with the same original
/// undeformed-volume gas references. Accepted deformed gas parcels are output
/// evidence and never replacement reference sources.
#[derive(Debug)]
pub struct FiniteThermalPorousCavitiesStep<
    'a,
    L: FiniteThermalPorousMaterial = AffineThermalFinitePorous,
    G = ClosedGasParcel,
> {
    /// Original thermal material/momentum step, including all cavity endpoint
    /// work in its equilibrium residual and complete balance error.
    pub step: FiniteThermalPorousTetStep<'a, L>,
    /// Independent accepted gas endpoints and tangents in original source order.
    pub cavities: Vec<TetCavityGasLoad<G>>,
    /// Direct gas storage, exterior work and endpoint correction for each source.
    pub cavity_energy: Vec<FinitePorousCavityEnergy>,
    /// Material defect plus held-pressure and individual cavity corrections.
    /// Physical caloric heat and numerical kinetic loss retain separate accounts.
    pub external_work_defect_j: f64,
}

impl<'a, L: FiniteThermalPorousMaterial> FiniteThermalPorousTetSolid<'a, L> {
    /// Advance thermal compaction against independently sealed original gases.
    ///
    /// Supplied endpoint temperatures enter the SAME original thermal family;
    /// accepted history is admitted at its old temperature. Current cavity
    /// pressures, follower geometry, acoustic rank-one factors, thermal porous
    /// tangents and optional full consistent M/dt² enter one native solve.
    /// Each ideal or Peng--Robinson source retains its own caloric, homogeneous
    /// phase and entropy admission. Gas volume follows total solid displacement,
    /// never the plastic volume or thermally reduced material deformation.
    ///
    /// Original gas references remain at their undeformed cavity volumes across
    /// increments. The shared interval owner uses direct caloric and compact
    /// geometric work, without subtracting absolute gas energies or stored
    /// reference-relative potentials. Cavity inventories and source order remain
    /// independent. Optional held pressure is additional to each cavity's fixed
    /// exterior pressure; the caller chooses the boundaries represented by each.
    ///
    /// Only original material yield heat plus signed reciprocal exchange enters
    /// an external caloric owner. This method prescribes solid temperatures and
    /// keeps gas evolution isentropic: it adds neither gas-wall heat exchange
    /// nor gas momentum/pressure waves. Numerical defects do not become heat.
    ///
    /// The material work budget bounds the sum of absolute material, pressure
    /// and individual cavity endpoint corrections; numerical kinetic loss has
    /// its independent budget. An empty cavity slice follows the original
    /// thermal pressure/dead-load path with identical accepted material/motion.
    ///
    /// # Errors
    /// Original source, temperature/history, cavity geometry, gas-domain, native
    /// solve or work-gate refusal, or cancellation. The borrowed thermal solid,
    /// original inertia and every original gas reference remain unchanged even
    /// after a converged coupled solve that fails a final acceptance gate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_cavities<G: AdiabaticCavityGas>(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        dt_s: f64,
        cavities: &[AdiabaticCavity<'_, G>],
        pressure: Option<(&ClosedTetPressureSurface, TetPressure<'_>)>,
        settings: HyperTetSettings,
        controls: FiniteThermoplasticControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermalPorousCavitiesStep<'a, L, G>, FiniteThermalPorousTetError> {
        self.advance_cavities_loaded(
            new_temperatures_k,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            dt_s,
            cavities,
            pressure,
            settings,
            controls,
            cx,
        )
    }
}
