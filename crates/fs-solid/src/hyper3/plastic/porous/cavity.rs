//! Reversible sealed gases and irreversible porous compaction in one endpoint.
//!
//! The native cavity batch supplies current pressure, follower geometry and
//! each source's acoustic rank-one tangent to the original Newton solver.
//! Frozen material trials all borrow the accepted plastic history. Neither a
//! gas endpoint nor a candidate plastic metric is published until every source,
//! constitutive, equilibrium and work gate has passed.
//!
//! Gas references remain at the undeformed cavity volumes across increments.
//! Accepted deformed parcels are output evidence, not new reference sources.
//! Work uses direct original-law energy and compact geometric increments;
//! subtracting reference-relative cavity potentials would lose small steps.
//! The retained dynamics wrapper supplies original solid inertia and physical
//! time integration; this owner adds no gas momentum, transport, thermal
//! exchange or pore-fluid model. The gas occupies explicit resolved voids,
//! not an inferred evolving porosity inside the material law.

use fs_exec::Cx;
use fs_material::plastic::finite::porous::{FinitePorous, FinitePorousMaterial};
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_math::dd::Dd;

use super::{FinitePorousTetError, FinitePorousTetSolid, FinitePorousTetStep};
use crate::hyper3::plastic::{
    AcceptedMotion, FinitePlasticStepControls, checkpoint, invalid, validate_controls,
};
use crate::hyper3::pressure::cavity::equilibrium::CavityBatch;
use crate::hyper3::pressure::cavity::{AdiabaticCavity, AdiabaticCavityGas, TetCavityGasLoad};
use crate::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use crate::hyper3::{HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField};

/// One complete accepted porous/gas endpoint.
///
/// `step.solid` retains the full material history. The gas loads and their work
/// receipts remain in the caller's original cavity order and never mix mass or
/// species. Call `step.solid.advance_with_cavities` with the same original
/// reference bindings for the next increment.
#[derive(Debug)]
pub struct FinitePorousCavitiesStep<'a, L: FinitePorousMaterial = FinitePorous, G = ClosedGasParcel>
{
    /// Native porous endpoint; its equilibrium residual includes all cavities.
    pub step: FinitePorousTetStep<'a, L>,
    /// Native accepted gas loads, including unchanged inventories and tangents.
    pub cavities: Vec<TetCavityGasLoad<G>>,
    /// Direct work for each independent source over this accepted increment.
    pub cavity_energy: Vec<FinitePorousCavityEnergy>,
    /// Material endpoint defect plus held-pressure and individual cavity work
    /// corrections, excluding equilibrium error and any inertial endpoint work.
    /// This signed integration defect is not gas heat or plastic heat.
    pub external_work_defect_j: f64,
}

/// Work of one original gas source between accepted solid configurations.
///
/// Positive work acts on the solid. Gas compression increases stored gas energy
/// and gives negative gas work on the solid; fixed exterior pressure contributes
/// `-p_exterior * delta(V_gas)`. These contributions are kept separate from the
/// constitutive plastic loss and from endpoint-force quadrature error.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FinitePorousCavityEnergy {
    /// Original gas owner's direct stored-energy change, U_after - U_before [J].
    pub gas_stored_energy_change_j: f64,
    /// Gas work on the solid, exactly the negative of its energy change [J].
    pub gas_work_j: f64,
    /// Fixed exterior pressure work on this cavity's solid boundary [J].
    pub exterior_work_j: f64,
    /// Net gas-plus-exterior work on the solid over the increment [J].
    pub work_j: f64,
    /// Accepted endpoint cavity force dotted with the displacement increment [J].
    pub endpoint_work_j: f64,
    /// `work_j - endpoint_work_j`, computed in compensated precision [J].
    pub endpoint_work_correction_j: f64,
}

/// Only work already admitted by original gas and geometry owners reaches the
/// common porous finalizer. Its absolute correction preserves source separation.
pub(super) struct CavityEndpointWork {
    pub(super) endpoint_work_j: f64,
    pub(super) absolute_correction_j: f64,
}

impl<'a, L: FinitePorousMaterial> FinitePorousTetSolid<'a, L> {
    /// Equilibrate independently sealed cavities against retained compaction.
    ///
    /// Every cavity must use this exact original solid reference mesh and its
    /// original undeformed-volume gas parcel. The same source isentrope is
    /// reconstructed at all accepted and trial cavity volumes. The original
    /// ideal or Peng--Robinson gas supplies its own calorics, homogeneous phase
    /// admission and acoustic derivative. Reference parcels are never mutated.
    ///
    /// Optional held pressure retains its original current-area load, tangent
    /// and work rule. It is additional to each cavity's explicit exterior
    /// pressure, so the caller must choose which boundary each load represents.
    /// The original support and referential dead/body loads remain available.
    ///
    /// One native physical endpoint is solved from the accepted displacement
    /// and frozen porous histories. The work budget bounds the sum of absolute
    /// constitutive, held-pressure and individual cavity endpoint defects.
    /// Opposite chamber errors cannot conceal an under-resolved increment.
    /// Physical yield loss and gas energy remain separate, without synthetic
    /// damping, gas mixing or a second equilibrium solver.
    ///
    /// An empty cavity slice follows the original pressure/dead-load path.
    /// Source order is retained in the result; mechanically coupled chambers
    /// still retain independent masses, species, entropy and exterior pressure.
    ///
    /// # Errors
    /// Original material, mesh, gas-domain, entropy, native solve or work-gate
    /// refusal, repeated cavity faces, mismatched reference volumes, or
    /// cancellation. The accepted solid and every borrowed reference parcel
    /// remain unchanged, including refusals after a converged native solve.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_cavities<G: AdiabaticCavityGas>(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        cavities: &[AdiabaticCavity<'_, G>],
        pressure: Option<(&ClosedTetPressureSurface, TetPressure<'_>)>,
        settings: HyperTetSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<FinitePorousCavitiesStep<'a, L, G>, FinitePorousTetError> {
        self.advance_cavities_loaded(
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            cavities,
            pressure,
            settings,
            controls,
            None,
            cx,
        )
        .map(|(step, _)| step)
    }

    /// One shared source/geometry/work owner for static and inertial endpoints.
    #[allow(
        clippy::too_many_arguments,
        clippy::too_many_lines,
        clippy::type_complexity
    )]
    pub(super) fn advance_cavities_loaded<G: AdiabaticCavityGas>(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        cavities: &[AdiabaticCavity<'_, G>],
        pressure: Option<(&ClosedTetPressureSurface, TetPressure<'_>)>,
        settings: HyperTetSettings,
        controls: FinitePlasticStepControls,
        inertia: Option<&crate::inertia::BackwardEulerInertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<(FinitePorousCavitiesStep<'a, L, G>, Option<AcceptedMotion>), FinitePorousTetError>
    {
        checkpoint(cx)?;
        if cavities.is_empty() {
            let (step, motion) = self.advance_loaded(
                prescribed_m,
                nodal_forces_n,
                body_force_n_m3,
                settings,
                controls,
                inertia,
                pressure,
                cx,
            )?;
            let external_work_defect_j = step
                .pressure
                .as_ref()
                .map_or(step.energy.endpoint_work_defect_j, |work| {
                    work.external_work_defect_j
                });
            return Ok((
                FinitePorousCavitiesStep {
                    step,
                    cavities: Vec::new(),
                    cavity_energy: Vec::new(),
                    external_work_defect_j,
                },
                motion,
            ));
        }
        validate_controls(settings, controls)?;
        self.validate_accepted_history(controls, cx)?;
        let batch = CavityBatch::new(cavities, self.state.nodes_m, self.state.tetrahedra, cx)?;
        let pressure = if let Some((surface, field)) = pressure {
            if !surface.matches_reference(self.state.nodes_m, self.state.tetrahedra) {
                return Err(
                    invalid("porous pressure surface must bind the exact original mesh").into(),
                );
            }
            surface.validate_pressure(field, cx)?;
            Some((surface, field))
        } else {
            None
        };
        // Retain one original-source gas epoch at the already accepted solid.
        // The batch admits the complete solid only once for multiple cavities.
        let old_u = self.state.displacement_m.as_flattened();
        let before = batch.evaluate(old_u, cx)?;
        let frozen = Self::freeze(&self.material, &self.state.element_states, cx)?;
        let problem = HyperTetProblem {
            nodes_m: self.state.nodes_m,
            tetrahedra: self.state.tetrahedra,
            material: HyperelasticLawField::PerElement(&frozen),
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            budget: self.state.budget,
            settings,
        };
        let solution = problem.solve_loaded_from_with_cavities(
            Some(old_u),
            Some(&self.state.elements),
            pressure,
            None,
            inertia,
            batch,
            cx,
        )?;
        let new_u = solution.solid.displacement_m.as_flattened();
        let mut cavity_energy = Vec::with_capacity(cavities.len());
        let mut endpoint_work = Dd::ZERO;
        let mut correction = Dd::ZERO;
        let mut absolute_correction = Dd::ZERO;
        for ((source, before), after) in cavities.iter().zip(&before).zip(&solution.cavities) {
            checkpoint(cx)?;
            let gas_change = before
                .gas
                .stored_energy_change_to(&after.gas, || cx.checkpoint().is_err())
                .map_err(|error| match error {
                    GasEnergyError::Cancelled => {
                        FinitePorousTetError::from(HyperTetError::Cancelled)
                    }
                    other => invalid(&format!(
                        "porous cavity gas energy increment refused: {other}"
                    ))
                    .into(),
                })?;
            // This interval is the native source-caloric and compact geometric
            // work, not subtraction of two reference-relative potentials. Its
            // Dd value retains exterior work separately when the gas term is
            // large; neither operation repeats the entropy solve.
            let potential = source.potential_increment(before, after, cx)?;
            let gas_change = Dd::from_f64(gas_change);
            let exterior_work = gas_change - potential;
            let work = -potential;
            let mut endpoint = Dd::ZERO;
            for (dof, (&force, (&after, &before))) in after
                .force_n
                .iter()
                .zip(new_u.iter().zip(old_u))
                .enumerate()
            {
                if dof % 256 == 0 {
                    checkpoint(cx)?;
                }
                endpoint =
                    endpoint + Dd::from_f64(force) * (Dd::from_f64(after) - Dd::from_f64(before));
            }
            let source_correction = work - endpoint;
            let receipt = FinitePorousCavityEnergy {
                gas_stored_energy_change_j: gas_change.to_f64(),
                gas_work_j: (-gas_change).to_f64(),
                exterior_work_j: exterior_work.to_f64(),
                work_j: work.to_f64(),
                endpoint_work_j: endpoint.to_f64(),
                endpoint_work_correction_j: source_correction.to_f64(),
            };
            if [
                receipt.gas_stored_energy_change_j,
                receipt.exterior_work_j,
                receipt.work_j,
                receipt.endpoint_work_j,
                receipt.endpoint_work_correction_j,
            ]
            .iter()
            .any(|value| !value.is_finite())
            {
                return Err(invalid("finite porous cavity work required").into());
            }
            endpoint_work = endpoint_work + endpoint;
            correction = correction + source_correction;
            absolute_correction =
                absolute_correction + Dd::from_f64(receipt.endpoint_work_correction_j.abs());
            cavity_energy.push(receipt);
        }
        let force = problem.external(&self.state.elements, 1.0);
        let (step, motion) = self.accept_endpoint(
            solution.solid,
            &force,
            controls,
            inertia,
            pressure,
            Some(CavityEndpointWork {
                endpoint_work_j: endpoint_work.to_f64(),
                absolute_correction_j: absolute_correction.to_f64(),
            }),
            cx,
        )?;
        let held_pressure_correction = step
            .pressure
            .as_ref()
            .map_or(0.0, |work| work.endpoint_work_correction_j);
        let external_work_defect_j = (Dd::from_f64(step.energy.endpoint_work_defect_j)
            + Dd::from_f64(held_pressure_correction)
            + correction)
            .to_f64();
        if !external_work_defect_j.is_finite() {
            return Err(invalid("finite combined porous/cavity work defect required").into());
        }
        checkpoint(cx)?;
        Ok((
            FinitePorousCavitiesStep {
                step,
                cavities: solution.cavities,
                cavity_energy,
                external_work_defect_j,
            },
            motion,
        ))
    }
}
