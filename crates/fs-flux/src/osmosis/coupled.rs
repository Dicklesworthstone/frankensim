//! One finite nonideal water/species network and one isothermal potential.
//!
//! Water uses the original fully retaining membrane owner; independent
//! stationary-solvent MS connections use the original molecular owner. The
//! ordered water-then-species composition is first-order Lie splitting, with
//! the existing full/two-half-sweep temporal control. These are separate
//! physical pathways, not an inferred solvent-drag or cross-coefficient law.
//!
//! For the declared pressure response, Psi(V)=integral p(V) dV and
//! d(F_mix+Psi)/dV=p-Pi. The existing hydraulic flux is driven by that SAME
//! potential difference. Therefore the closed isothermal account is
//! Phi=F_mix+sum(Psi), not F_mix alone. Prescribed pressures contribute load
//! work through their exact primitive. No heat reservoir, caloric energy,
//! fluid inertia or pressure-dependent liquid kinetics is supplied here.

use super::{
    MAX_MEMBRANES, MembraneBoundary, OsmoticCompartment, OsmoticMembrane, OsmoticNetwork,
    OsmoticNetworkError, OsmoticStepControls, compensated_add,
    diffusion::{self, SolutionDiffusionControls, SolutionDiffusionLink},
    exchange, poll, positive,
};
use fs_material::osmosis::OsmosisError;

/// One set of physical whole-window budgets, separate from temporal accuracy.
#[derive(Clone, Copy, Debug)]
pub struct CoupledSolutionControls {
    /// Original molecular controls. In coupled mode its free-energy allowance
    /// bounds TOTAL mixing plus mechanical potential, including water stages;
    /// it is not a requirement that mixing free energy alone always decrease.
    pub diffusion: SolutionDiffusionControls,
    /// Absolute actual original-to-final closed solvent-volume residual [m³].
    pub water_balance_tolerance_m3: f64,
}

/// Reference-cancelled original-to-final closed-system potential changes [J].
/// Common concentration and pressure references cancel analytically; actual
/// species and water residuals are checked separately, without reference-sized
/// energy allowances. This report is not a caloric energy/heat balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledPotentialChange {
    pub mixing_free_energy_change_j: f64,
    /// Exact pressure primitive, including prescribed pressure-load work.
    /// Its negative is mechanical work delivered to the modeled solutions.
    pub mechanical_potential_change_j: f64,
    /// Compensated sum of all original source and mechanical changes.
    pub total_potential_change_j: f64,
}

impl OsmoticNetwork {
    /// Prepare one finite, closed water/species network from existing owners.
    /// Every compartment uses the same original virial model, temperature,
    /// solvent and species ordering. Water membranes must be bare and fully
    /// retaining; the original ideal leaky/film law is not made nonideal.
    /// Molecular links are independent stationary-solvent exchange pathways.
    ///
    /// Either path set may be empty for its physical limiting case. Within a
    /// path set, unordered endpoint pairs must be distinct. Water and molecular
    /// paths may join the same compartments. Their areas and coefficients are
    /// independent caller inputs, never inferred from each other.
    ///
    /// Uses the existing prepared buffers and complete-duration controller.
    /// Source-backed molecular links retain their documented transient-query
    /// allocation cost and full-card payload charge; constant data allocate no
    /// memory during successful advancement. Late failure publishes no partial
    /// volumes, amounts or clock.
    #[allow(clippy::too_many_arguments)]
    pub fn with_coupled_solution_transport(
        species: Vec<String>,
        cells: Vec<OsmoticCompartment>,
        mut membranes: Vec<OsmoticMembrane>,
        mut links: Vec<SolutionDiffusionLink>,
        controls: CoupledSolutionControls,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, OsmoticNetworkError> {
        if !positive(controls.water_balance_tolerance_m3)
            || membranes.len().saturating_add(links.len()) > MAX_MEMBRANES
        {
            return Err(OsmoticNetworkError::Invalid(
                "bounded combined graph and positive closed-water balance allowance required",
            ));
        }
        for membrane in &mut membranes {
            if !matches!(membrane.material, MembraneBoundary::WellMixed(_))
                || membrane
                    .membrane_law()
                    .solutes()
                    .iter()
                    .any(|(reflection, permeance)| {
                        reflection.value() != 1.0 || permeance.value() != 0.0
                    })
            {
                return Err(OsmoticNetworkError::Invalid(
                    "coupled nonideal water paths require bare fully retaining membranes",
                ));
            }
            if membrane.left > membrane.right {
                core::mem::swap(&mut membrane.left, &mut membrane.right);
            }
        }
        membranes.sort_by_key(|membrane| (membrane.left, membrane.right));
        if membranes
            .windows(2)
            .any(|pair| pair[0].left == pair[1].left && pair[0].right == pair[1].right)
        {
            return Err(OsmoticNetworkError::Invalid(
                "duplicate physical water face",
            ));
        }
        let bytes = diffusion::prepare_solution_diffusion(
            &species,
            &cells,
            &mut links,
            controls.diffusion,
            maximum_workspace_bytes,
        )?;
        let mut network = Self::new(species, cells, membranes, maximum_workspace_bytes - bytes)?;
        network.workspace_bytes += bytes;
        network.diffusion = links;
        network.diffusion_controls = Some(controls.diffusion);
        network.coupled_controls = Some(controls);
        Ok(network)
    }
}

pub(super) fn admit_water_balance(
    residual_m3: f64,
    maximum_m3: f64,
) -> Result<(), OsmoticNetworkError> {
    if !residual_m3.is_finite() || residual_m3.abs() > maximum_m3 {
        Err(OsmoticNetworkError::CoupledWaterBalance {
            residual_m3,
            maximum_m3,
        })
    } else {
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn sweep(
    cells: &mut [OsmoticCompartment],
    membranes: &[OsmoticMembrane],
    links: &mut [SolutionDiffusionLink],
    physical: CoupledSolutionControls,
    dt: f64,
    numerical: OsmoticStepControls,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), OsmoticNetworkError> {
    let pressure_reference = cells[0].pressure.reference_pressure_pa;
    for membrane in membranes {
        poll(cancelled)?;
        let (left, right) = (&cells[membrane.left], &cells[membrane.right]);
        // Original implicit water solve and original thermodynamic pressure.
        let endpoint = exchange(left, right, membrane, dt, numerical, evaluations, cancelled)?;
        let n = left.solution.solute_moles().len();
        let water_residual = (endpoint.volumes[0] - left.solution.solvent_volume_m3())
            + (endpoint.volumes[1] - right.solution.solvent_volume_m3());
        admit_water_balance(water_residual, physical.water_balance_tolerance_m3)?;
        let mut change = 0.0;
        let mut correction = 0.0;
        for (old, volume, amounts) in [
            (left, endpoint.volumes[0], &endpoint.left[..n]),
            (right, endpoint.volumes[1], &endpoint.right[..n]),
        ] {
            let (mixing, mechanical) = cell_change(old, volume, amounts, pressure_reference)?;
            compensated_add(&mut change, &mut correction, mixing);
            compensated_add(&mut change, &mut correction, mechanical);
        }
        finite(change)?;
        if change > physical.diffusion.free_energy_tolerance_j {
            return Err(OsmoticNetworkError::CoupledStepFreeEnergy {
                change_j: change,
                maximum_j: physical.diffusion.free_energy_tolerance_j,
            });
        }
        poll(cancelled)?;
        cells[membrane.left]
            .solution
            .replace_inventories(endpoint.volumes[0], &endpoint.left[..n])?;
        cells[membrane.right]
            .solution
            .replace_inventories(endpoint.volumes[1], &endpoint.right[..n])?;
    }
    // Dilution/concentration from water motion changes the actual coefficient
    // coordinates even when a solute amount was retained exactly.
    diffusion::qualify_sources(cells, links, cancelled)?;
    diffusion::sweep(cells, links, dt, physical.diffusion, evaluations, cancelled)
}

fn cell_change(
    before: &OsmoticCompartment,
    after_volume: f64,
    after_amounts: &[f64],
    pressure_reference: f64,
) -> Result<(f64, f64), OsmoticNetworkError> {
    let law = before
        .solution
        .nonideal_model()
        .ok_or(OsmoticNetworkError::Invalid(
            "coupled potential needs the original virial solution model",
        ))?;
    let mixing = law
        .reference_cancelled_inventory_free_energy_change_j(
            before.solution.solvent_volume_m3(),
            before.solution.solute_moles(),
            after_volume,
            after_amounts,
        )
        .map_err(OsmosisError::from)?;
    let mechanical = before.pressure.reference_cancelled_potential_change_j(
        before.solution.solvent_volume_m3(),
        after_volume,
        pressure_reference,
    )?;
    Ok((mixing, mechanical))
}

pub(super) fn potential_change(
    original: &[OsmoticCompartment],
    candidate: &[OsmoticCompartment],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<CoupledPotentialChange, OsmoticNetworkError> {
    let mut values = [0.0; 3];
    let mut corrections = [0.0; 3];
    let pressure_reference = original[0].pressure.reference_pressure_pa;
    for (before, after) in original.iter().zip(candidate) {
        poll(cancelled)?;
        if before.solution.nonideal_model() != after.solution.nonideal_model() {
            return Err(OsmoticNetworkError::Invalid(
                "coupled transport changed original thermodynamics",
            ));
        }
        let (mixing, mechanical) = cell_change(
            before,
            after.solution.solvent_volume_m3(),
            after.solution.solute_moles(),
            pressure_reference,
        )?;
        compensated_add(&mut values[0], &mut corrections[0], mixing);
        compensated_add(&mut values[1], &mut corrections[1], mechanical);
        compensated_add(&mut values[2], &mut corrections[2], mixing);
        compensated_add(&mut values[2], &mut corrections[2], mechanical);
    }
    Ok(CoupledPotentialChange {
        mixing_free_energy_change_j: finite(values[0])?,
        mechanical_potential_change_j: finite(values[1])?,
        total_potential_change_j: finite(values[2])?,
    })
}

fn finite(value: f64) -> Result<f64, OsmoticNetworkError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(OsmoticNetworkError::Invalid(
            "coupled potential change is not representable",
        ))
    }
}

#[cfg(test)]
mod tests;
