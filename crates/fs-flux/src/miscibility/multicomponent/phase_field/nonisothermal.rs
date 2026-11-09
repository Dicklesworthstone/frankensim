//! Closed local-temperature demixing with original multicomponent calorics.
//!
//! The SAME homogeneous A+B/T regular-solution law and original pure-liquid Cp
//! curves supply chemical forces, transported partial enthalpies and local
//! temperature recovery. Constant energetic interfaces retain the graph law
//! `E_e = kappa A/(4 d) sum_s (x_js-x_is)^2` and carry no interface entropy.
//!
//! Each cell stores `U_i = H_i + sum_{e incident i} E_e/2`. At fixed local U,
//! the interface contribution to the entropy-conjugate species potential is
//! `sum_j K_ij (1/T_i+1/T_j) (x_is-x_js)/N_i`, with `K_ij=kappa A/(4 d)`.
//! Full Maxwell--Stefan transport uses this potential together with the original
//! composition effort and authored reciprocal Soret--Dufour contrasts.
//!
//! After simultaneous species and paired total-energy transfers, each changed
//! edge contributes half its energy increment to EACH endpoint's bulk caloric
//! inverse. A neighbor with unchanged composition can therefore change
//! temperature. Every candidate closes original independent inventories,
//! original homogeneous enthalpy PLUS interface energy, and original entropy.
//! The original carried-inventory caloric precision guard remains mandatory.
//!
//! This explicit method is not unconditionally entropy stable. A refused
//! interval must be refined. The positive-interface owner admits nonconvex
//! original laws without changing the ordinary homogeneous convexity rule.
//! Fixed graph geometry, common fixed partial molar volumes, static pressure,
//! common positive species support and constant energetic kappa are retained.
//! A prescribed-flow phase-field advection law is not part of this owner.

use super::{MulticomponentPhaseField, MulticomponentPhaseFieldNetwork};
use super::solve as interface_solve;
use super::super::{
    MulticomponentMiscibilityError, Sum, amount_sum, finite, positive,
};
use super::super::nonisothermal as thermal;
use super::super::solve as spatial_solve;
use thermal::{MulticomponentNonisothermalField, admit_controls, poll};
use thermal::solve as thermal_solve;
use fs_qty::Temperature;
use fs_thermochem::solution::multicomponent::caloric::MulticomponentLiquidCaloricModelV1;

pub use thermal::{
    MulticomponentNonisothermalError, MulticomponentNonisothermalFace,
    MulticomponentNonisothermalSolve,
};

type Result<T> = core::result::Result<T, MulticomponentNonisothermalError>;

pub mod adaptive;
mod solve;

/// Accepted local temperatures, independent inventories and energetic
/// interfaces. Neither the raw ordinary thermal field nor its nonconvex
/// zero-gradient network is exposed.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalPhaseField {
    network: MulticomponentPhaseFieldNetwork,
    bulk: MulticomponentNonisothermalField,
    workspace_bytes: usize,
}

/// Complete physical-interval inventories, energy and entropy receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalPhaseFieldStep {
    pub duration_s: f64,
    /// Canonical face then canonical species, positive left to right.
    pub transferred_mol: Vec<Vec<f64>>,
    /// Paired TOTAL local-energy transfer [J], including original partial
    /// enthalpy, authored reciprocal Dufour transport and Fourier heat.
    pub transferred_energy_j: Vec<f64>,
    pub maxwell_stefan_evaluations: usize,
    pub maximum_gradient_residual_per_m: f64,
    pub maximum_molar_frame_residual_mol_m2_s: f64,
    pub caloric_iterations: u64,
    pub component_residual_mol: Vec<f64>,
    pub maximum_capacity_residual_mol: f64,
    /// Original homogeneous enthalpy change plus energetic interface change.
    pub energy_residual_j: f64,
    pub entropy_change_j_per_k: f64,
    pub homogeneous_enthalpy_change_j: f64,
    pub interfacial_energy_change_j: f64,
    /// Initial constitutive rate from bulk-plus-interface chemical efforts.
    /// The complete finite-step entropy change is checked independently.
    pub entropy_production_rate_w_per_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MulticomponentNonisothermalPhaseFieldCandidate {
    field: MulticomponentNonisothermalPhaseField,
    step: MulticomponentNonisothermalPhaseFieldStep,
}

impl MulticomponentNonisothermalPhaseFieldCandidate {
    pub const fn field(&self) -> &MulticomponentNonisothermalPhaseField {
        &self.field
    }

    pub const fn step(&self) -> &MulticomponentNonisothermalPhaseFieldStep {
        &self.step
    }

    pub fn into_parts(
        self,
    ) -> (
        MulticomponentNonisothermalPhaseField,
        MulticomponentNonisothermalPhaseFieldStep,
    ) {
        (self.field, self.step)
    }
}

impl MulticomponentPhaseFieldNetwork {
    /// Initialize local temperature without losing accepted original amounts,
    /// their carried low words, interface geometry or elapsed physical time.
    /// This is an initialization, not a conservative heat-transfer event.
    ///
    /// The original caloric bracket, source identities, positive support and
    /// coefficient checks all remain active. Only homogeneous convexity is
    /// waived inside this privately retained positive-interface owner.
    #[allow(clippy::too_many_arguments)]
    pub fn bind_nonisothermal_phase_field(
        &self,
        initial: &MulticomponentPhaseField,
        caloric: MulticomponentLiquidCaloricModelV1,
        temperatures: Vec<Temperature>,
        faces: Vec<MulticomponentNonisothermalFace>,
        temperature_bracket: [Temperature; 2],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentNonisothermalPhaseField> {
        poll(&mut cancelled)?;
        if self != initial.network() || !positive(self.gradient_coefficient_j_m()) {
            return Err(MulticomponentNonisothermalError::Invalid(
                "same original positive-interface phase-field owner required",
            ));
        }
        // Retain the isothermal graph's interface storage allowance and add
        // compensated temperature-weighted gradients, reduced efforts, local
        // energy allocations and simultaneous wrapper candidates. The original
        // thermal bound already covers bulk candidates and caloric work.
        let interface_storage = self
            .workspace_bytes()
            .checked_sub(self.bulk.workspace_bytes())
            .ok_or(MulticomponentNonisothermalError::Invalid(
                "phase-field interface workspace bound",
            ))?;
        let extra = self
            .cells()
            .len()
            .checked_mul(self.species_count())
            .and_then(|v| v.checked_mul(64))
            .and_then(|v| {
                self.cells().len().checked_mul(64).and_then(|w| v.checked_add(w))
            })
            .and_then(|v| v.checked_add(interface_storage))
            .and_then(|v| {
                v.checked_add(3 * core::mem::size_of::<MulticomponentNonisothermalPhaseField>())
            })
            .ok_or(MulticomponentNonisothermalError::Invalid(
                "thermal phase-field workspace size overflow",
            ))?;
        let bulk_budget = maximum_workspace_bytes.checked_sub(extra).ok_or(
            MulticomponentNonisothermalError::Workspace {
                required: extra,
                maximum: maximum_workspace_bytes,
            },
        )?;
        let bulk = self
            .bulk
            .prepare_nonisothermal_field(
                &initial.bulk,
                caloric,
                temperatures,
                faces,
                temperature_bracket,
                bulk_budget,
                false,
                &mut cancelled,
            )
            .map_err(|error| match error {
                MulticomponentNonisothermalError::Workspace { required, .. } => {
                    MulticomponentNonisothermalError::Workspace {
                        required: required.saturating_add(extra),
                        maximum: maximum_workspace_bytes,
                    }
                }
                other => other,
            })?;
        let workspace_bytes = bulk.workspace_bytes().checked_add(extra).ok_or(
            MulticomponentNonisothermalError::Invalid(
                "thermal phase-field workspace size overflow",
            ),
        )?;
        let field = MulticomponentNonisothermalPhaseField {
            network: self.clone(),
            bulk,
            workspace_bytes,
        };
        field.interfacial_energy_j(&mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(field)
    }
}

impl MulticomponentNonisothermalPhaseField {
    /// The positive-interface graph; its construction temperature does not
    /// replace any of this field's actual local temperatures.
    pub const fn network(&self) -> &MulticomponentPhaseFieldNetwork {
        &self.network
    }

    pub fn caloric_model(&self) -> &MulticomponentLiquidCaloricModelV1 {
        self.bulk.caloric_model()
    }

    pub fn component_moles(&self) -> &[f64] {
        self.bulk.component_moles()
    }

    pub fn mole_fractions(&self) -> &[f64] {
        self.bulk.mole_fractions()
    }

    pub fn cell_component_moles(&self, cell: usize) -> Option<&[f64]> {
        self.bulk.cell_component_moles(cell)
    }

    pub fn cell_mole_fractions(&self, cell: usize) -> Option<&[f64]> {
        self.bulk.cell_mole_fractions(cell)
    }

    pub fn temperatures(&self) -> &[Temperature] {
        self.bulk.temperatures()
    }

    pub fn elapsed_s(&self) -> f64 {
        self.bulk.elapsed_s()
    }

    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    pub fn transport_faces(&self) -> &[MulticomponentNonisothermalFace] {
        self.bulk.transport_faces()
    }

    pub fn temperature_bracket(&self) -> [Temperature; 2] {
        self.bulk.temperature_bracket()
    }

    /// Actual nonnegative interfacial storage in canonical graph order.
    pub fn interfacial_energy_j(&self, mut cancelled: impl FnMut() -> bool) -> Result<f64> {
        let mut energy = Sum::default();
        for face in 0..self.network.faces().len() {
            poll(&mut cancelled)?;
            energy.add(solve::edge_energy_j(self, face)?)?;
        }
        poll(&mut cancelled)?;
        Ok(energy.value()?)
    }

    /// One synchronous full-Maxwell--Stefan and caloric interval, atomically
    /// admitted by total energy and original entropy. No physical subcycling
    /// or dependent-species reconstruction occurs inside this step.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        controls: MulticomponentNonisothermalSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MulticomponentNonisothermalPhaseFieldCandidate> {
        solve::step(self, duration_s, controls, cancelled)
    }

    fn balances(
        &self,
        next: &Self,
        controls: MulticomponentNonisothermalSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, f64, f64, f64, f64, f64)> {
        admit_controls(controls)?;
        // This original account retains the carried-inventory caloric
        // resolution guard; it does not require bulk enthalpy alone to close.
        let (species, capacity, homogeneous, entropy) =
            self.bulk.caloric_balances(&next.bulk, controls, cancelled)?;
        let interfacial = solve::interfacial_energy_change_j(self, next, cancelled)?;
        let energy =
            solve::admit_total_energy_entropy(homogeneous, interfacial, entropy, controls)?;
        Ok((species, capacity, energy, entropy, homogeneous, interfacial))
    }
}

#[cfg(test)]
mod tests;
