//! Original finite aqueous exchange, lateral diffusion and material motion.
//!
//! Two original soluble half steps surround the existing mechanical trial.
//! Each chemical stage uses its actual face areas; mechanics keeps each face's
//! intermediate amount fixed. The original finite bath supplies every chemical
//! potential and is published once with the complete mechanical endpoint.
//! Backward-Euler adsorption makes this composition FIRST order despite its
//! symmetric arrangement. Refine the macro step independently of the chemical
//! substep limit; no mechanical or coupled temporal-error estimate is implied.
//!
//! Temperature, well-mixed solvent volume, topology, coefficients and nodal
//! reference masses stay fixed. This represents a material surface exchanging
//! a dilute, inertia-free coating with an independently held-volume bath. It
//! does not infer enclosed bulk volume, solvent displacement, osmotic traction,
//! solute momentum, adsorption heat or gas/liquid phase-change work.

use super::*;
use crate::capillary::surface::surfactant::soluble::{
    SolubleSurfaceControls, SolubleSurfaceError, SolubleSurfaceKinetics, SolubleSurfaceReport,
    advance_soluble_surface, energy_change,
};
use fs_material::osmosis::AqueousSolution;

/// Original per-face Langmuir laws, finite-bath kinetics and held surface graph.
#[derive(Clone, Copy, Debug)]
pub struct SolubleSurfaceTransport<'a> {
    pub laws: &'a [SurfactantLaw],
    pub kinetics: SolubleSurfaceKinetics<'a>,
    pub edges: &'a [SurfacePartitionDiffusionEdge],
    /// Substep, root and pair limits cover BOTH chemical halves. Amount and
    /// energy allowances also apply once to the complete proposed endpoint.
    pub controls: SolubleSurfaceControls,
}

/// No error publishes a partial bath, face amount or mechanical trajectory.
#[derive(Debug)]
pub enum SolubleSurfaceDynamicsError {
    Mechanics(SurfaceDynamicsError),
    Transport {
        phase: SurfaceDiffusionPhase,
        error: SolubleSurfaceError,
    },
    /// Original-to-final inventory or chemical/mechanical account refused.
    Account(SolubleSurfaceError),
}
impl core::fmt::Display for SolubleSurfaceDynamicsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "moving soluble surface: {self:?}")
    }
}
impl std::error::Error for SolubleSurfaceDynamicsError {}
impl From<SurfaceDynamicsError> for SolubleSurfaceDynamicsError {
    fn from(error: SurfaceDynamicsError) -> Self {
        Self::Mechanics(error)
    }
}

/// Closed bulk/surface account with geometric work separated from chemistry.
#[derive(Clone, Debug)]
pub struct SolubleSurfaceDynamicsReport {
    pub mechanics: SurfaceDynamicsReport,
    pub before_motion: SolubleSurfaceReport,
    pub after_motion: SolubleSurfaceReport,
    pub total_substeps: usize,
    pub total_root_evaluations: usize,
    pub total_surface_pair_exchanges: usize,
    pub total_amount_before_mol: f64,
    pub total_amount_after_mol: f64,
    pub amount_defect_mol: f64,
    /// Independently reconstructed original-to-final bulk mixing + surface
    /// free energy + standard-potential offset times surface amount [J].
    pub free_energy_change_j: f64,
    /// Only the two held-geometry chemical stages, excluding geometric work.
    pub chemical_free_energy_change_j: f64,
    pub exchange_free_energy_dissipation_j: f64,
    pub surface_diffusion_free_energy_dissipation_j: f64,
    /// Endpoint free-energy change minus both chemical changes and the native
    /// fixed-amount surface change. This independently checks the split account.
    pub energy_account_residual_j: f64,
    /// Delta K + delta F_bulk+surface + surface viscous loss + adsorption and
    /// diffusion free-energy losses - external mechanical work [J]. Loss is an
    /// isothermal free-energy account, not an automatically deposited heat.
    pub energy_residual_j: f64,
}

impl CapillarySurfaceDynamics<'_> {
    /// Advance moving material faces and one original finite soluble bath.
    ///
    /// The first soluble half step changes the same amounts that supply the
    /// original mechanical force/tangent. The second uses the solved areas,
    /// including their changed packing capacities. The common standard-potential
    /// offset is constant at fixed amount and adds no capillary force. No face
    /// receives a duplicated bulk inventory or an independent reservoir debit.
    ///
    /// All original adsorption, diffusion, source-domain and mechanical gates
    /// remain active. Shared work limits, complete amount conservation, the
    /// independently reconstructed free-energy account and final cancellation
    /// precede publication of all three physical states. Trajectory history is
    /// retained through the original publisher and never cloned for the trial.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_soluble_surfactant(
        &self,
        state: &mut SecondOrderState,
        bulk: &mut AqueousSolution,
        amounts_mol: &mut [f64],
        transport: SolubleSurfaceTransport<'_>,
        viscosity: SurfaceViscosityField<'_>,
        external_forces_n: &[Vector],
        dt_s: f64,
        controls: SurfaceDynamicsControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SolubleSurfaceDynamicsReport, SolubleSurfaceDynamicsError> {
        if cancelled() {
            return Err(SurfaceDynamicsError::Cancelled.into());
        }
        if state.q.len() != 3 * self.masses_kg.len()
            || amounts_mol.len() != self.faces.len()
            || transport.laws.len() != self.faces.len()
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || 0.5 * dt_s == 0.0
        {
            return Err(SurfaceDynamicsError::Invalid(
                "complete material state, one original law/amount per face and representable positive halves required",
            )
            .into());
        }
        if self.faces.len() > transport.controls.maximum_faces
            || transport.edges.len() > transport.controls.maximum_edges
        {
            return Err(stage_error(
                SurfaceDiffusionPhase::BeforeMotion,
                SolubleSurfaceError::WorkBudget,
            ));
        }
        let before_areas = face_areas(self.faces, &state.q, &mut cancelled)?;
        let mut candidate_amounts = Vec::new();
        candidate_amounts
            .try_reserve_exact(amounts_mol.len())
            .map_err(|_| SurfaceDynamicsError::Allocation)?;
        candidate_amounts.extend_from_slice(amounts_mol);
        let mut candidate_bulk = bulk.clone();
        let before = advance_soluble_surface(
            &mut candidate_bulk,
            transport.laws,
            &before_areas,
            &mut candidate_amounts,
            transport.kinetics,
            transport.edges,
            0.5 * dt_s,
            transport.controls,
            &mut cancelled,
        )
        .map_err(|error| stage_error(SurfaceDiffusionPhase::BeforeMotion, error))?;
        let (candidate, mechanics) = {
            let mut coating = Vec::new();
            coating
                .try_reserve_exact(self.faces.len())
                .map_err(|_| SurfaceDynamicsError::Allocation)?;
            for (i, (&law, &amount_mol)) in
                transport.laws.iter().zip(&candidate_amounts).enumerate()
            {
                if i % 64 == 0 && cancelled() {
                    return Err(SurfaceDynamicsError::Cancelled.into());
                }
                coating.push(SurfactantFace { law, amount_mol });
            }
            self.surface_trial(
                state,
                SurfaceForces::Insoluble(SurfactantField::PerFace(&coating)),
                viscosity,
                external_forces_n,
                dt_s,
                controls,
                &mut cancelled,
            )?
        };
        let after_areas = face_areas(self.faces, &candidate.q, &mut cancelled)?;
        let remaining = remaining_controls(transport.controls, before)?;
        let after = advance_soluble_surface(
            &mut candidate_bulk,
            transport.laws,
            &after_areas,
            &mut candidate_amounts,
            transport.kinetics,
            transport.edges,
            0.5 * dt_s,
            remaining,
            &mut cancelled,
        )
        .map_err(|error| {
            let error = match error {
                SolubleSurfaceError::RootBudget {
                    completed_evaluations,
                } => match before.root_evaluations.checked_add(completed_evaluations) {
                    Some(completed_evaluations) => SolubleSurfaceError::RootBudget {
                        completed_evaluations,
                    },
                    None => SolubleSurfaceError::WorkBudget,
                },
                error => error,
            };
            stage_error(SurfaceDiffusionPhase::AfterMotion, error)
        })?;
        let mut defect = Sum::default();
        defect.add(
            candidate_bulk.solute_moles()[transport.kinetics.bulk_species]
                - bulk.solute_moles()[transport.kinetics.bulk_species],
        );
        for (i, (&old, &new)) in amounts_mol.iter().zip(&candidate_amounts).enumerate() {
            if i % 64 == 0 && cancelled() {
                return Err(SurfaceDynamicsError::Cancelled.into());
            }
            defect.add(new - old);
        }
        let defect = finite(defect.value())?;
        if defect.abs() > transport.controls.amount_tolerance_mol {
            return Err(SolubleSurfaceDynamicsError::Account(
                SolubleSurfaceError::AmountBalance { defect_mol: defect },
            ));
        }
        let change = energy_change(
            bulk,
            transport.kinetics.bulk_species,
            bulk.solute_moles(),
            candidate_bulk.solute_moles(),
            transport.laws,
            &before_areas,
            &after_areas,
            amounts_mol,
            &candidate_amounts,
            transport.kinetics.standard_potential_offset_j_mol,
            &mut cancelled,
        )
        .map_err(SolubleSurfaceDynamicsError::Account)?;
        let chemical = sum(&[before.free_energy_change_j, after.free_energy_change_j])?;
        if chemical > transport.controls.energy_tolerance_j {
            return Err(SolubleSurfaceDynamicsError::Account(
                SolubleSurfaceError::FreeEnergy { change_j: chemical },
            ));
        }
        let account = sum(&[change, -chemical, -mechanics.surface.energy_change_j])?;
        if account.abs() > transport.controls.energy_tolerance_j {
            return Err(SolubleSurfaceDynamicsError::Account(
                SolubleSurfaceError::EnergyAccount {
                    residual_j: account,
                },
            ));
        }
        let exchange_loss = sum(&[
            before.exchange_free_energy_dissipation_j,
            after.exchange_free_energy_dissipation_j,
        ])?;
        let diffusion_loss = sum(&[
            before.surface_diffusion_free_energy_dissipation_j,
            after.surface_diffusion_free_energy_dissipation_j,
        ])?;
        let energy = sum(&[
            mechanics.kinetic_change_j,
            change,
            mechanics.surface_viscous_dissipation_j,
            exchange_loss,
            diffusion_loss,
            -mechanics.external_work_j,
        ])?;
        if energy.abs() > controls.energy_tolerance_j {
            return Err(SurfaceDynamicsError::Balance {
                impulse_n_s: mechanics.impulse_residual_n_s,
                kinematic_m: mechanics.kinematic_residual_m,
                energy_j: energy,
            }
            .into());
        }
        let report = SolubleSurfaceDynamicsReport {
            total_substeps: add_work(before.substeps, after.substeps)?,
            total_root_evaluations: add_work(before.root_evaluations, after.root_evaluations)?,
            total_surface_pair_exchanges: add_work(
                before.surface_pair_exchanges,
                after.surface_pair_exchanges,
            )?,
            total_amount_before_mol: sum(&[
                before.bulk_amount_before_mol,
                before.surface_amount_before_mol,
            ])?,
            total_amount_after_mol: sum(&[
                after.bulk_amount_after_mol,
                after.surface_amount_after_mol,
            ])?,
            amount_defect_mol: defect,
            free_energy_change_j: change,
            chemical_free_energy_change_j: chemical,
            exchange_free_energy_dissipation_j: exchange_loss,
            surface_diffusion_free_energy_dissipation_j: diffusion_loss,
            energy_account_residual_j: account,
            energy_residual_j: energy,
            mechanics,
            before_motion: before,
            after_motion: after,
        };
        if cancelled() {
            return Err(SurfaceDynamicsError::Cancelled.into());
        }
        publish_state(state, candidate);
        *bulk = candidate_bulk;
        amounts_mol.copy_from_slice(&candidate_amounts);
        Ok(report)
    }
}

fn stage_error(
    phase: SurfaceDiffusionPhase,
    error: SolubleSurfaceError,
) -> SolubleSurfaceDynamicsError {
    SolubleSurfaceDynamicsError::Transport { phase, error }
}

fn remaining_controls(
    controls: SolubleSurfaceControls,
    first: SolubleSurfaceReport,
) -> Result<SolubleSurfaceControls, SolubleSurfaceDynamicsError> {
    let budget = || {
        stage_error(
            SurfaceDiffusionPhase::AfterMotion,
            SolubleSurfaceError::WorkBudget,
        )
    };
    let mut remaining = controls;
    remaining.maximum_substeps = controls
        .maximum_substeps
        .checked_sub(first.substeps)
        .ok_or_else(budget)?;
    remaining.maximum_surface_pair_exchanges = controls
        .maximum_surface_pair_exchanges
        .checked_sub(first.surface_pair_exchanges)
        .ok_or_else(budget)?;
    remaining.maximum_root_evaluations = controls
        .maximum_root_evaluations
        .checked_sub(first.root_evaluations)
        .ok_or_else(budget)?;
    if remaining.maximum_substeps == 0 {
        return Err(budget());
    }
    if remaining.maximum_root_evaluations == 0 {
        return Err(stage_error(
            SurfaceDiffusionPhase::AfterMotion,
            SolubleSurfaceError::RootBudget {
                completed_evaluations: first.root_evaluations,
            },
        ));
    }
    Ok(remaining)
}

fn add_work(a: usize, b: usize) -> Result<usize, SolubleSurfaceDynamicsError> {
    a.checked_add(b).ok_or(SolubleSurfaceDynamicsError::Account(
        SolubleSurfaceError::WorkBudget,
    ))
}
fn finite(value: f64) -> Result<f64, SolubleSurfaceDynamicsError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SolubleSurfaceDynamicsError::Account(
            SolubleSurfaceError::Unrepresentable,
        ))
    }
}
fn sum(values: &[f64]) -> Result<f64, SolubleSurfaceDynamicsError> {
    let mut total = Sum::default();
    for &value in values {
        total.add(value);
    }
    finite(total.value())
}

#[cfg(test)]
mod tests;
