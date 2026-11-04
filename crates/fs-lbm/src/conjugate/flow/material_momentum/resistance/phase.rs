//! Original phase/source resistance on authoritative current mass and momentum.

use super::*;
use crate::conjugate::flow::density_projection::admit_feedback_controls;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::phase::mushy::{MushyPhaseLaw, PhaseMobility};
use fs_material::phase::transport::ResolvedPhaseTransportState;

/// Original chart state and extensive inventory at one explicit state boundary.
#[derive(Clone, Copy, Debug)]
pub struct PhaseResistanceMomentumInput<'a> {
    /// The original M/P, body forcing and fresh authored volume constraint.
    pub mechanics: MaterialDensityMomentumInput<'a>,
    /// One original h-resolved state per fluid cell, None in fixed grid solids.
    pub phase_states: &'a [Option<EquilibriumPhaseState>],
}

/// Native mechanical gates and original chart/occupied-volume admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseResistanceMomentumControls {
    pub mechanics: MaterialDensityMomentumControls,
    /// Maximum absolute M/rho_chart(h) - original cell volume, m^3. An admitted
    /// representation difference never changes M or becomes a volume source.
    pub maximum_cell_volume_defect_m3: f64,
}

/// One fresh source-bound mechanical endpoint, without caloric state evolution.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedPhaseResistanceMomentumStep {
    /// Actual M/V inertia and conservative momentum feedback. Physical cell
    /// drag heat is reported for the enclosing caloric owner, not inserted here.
    pub evolution: ResistiveMaterialDensityMomentumStep,
    /// Frozen mu*K^-1, Pa s/m^2, from the exact retained original sources.
    /// Fully liquid and fixed-grid-solid entries are exactly zero.
    pub resistance_pa_s_m2: Vec<f64>,
    /// One complete original phase/material source receipt per fluid cell.
    /// Fixed grid solids contain None and never trigger a source query.
    pub original_sources: Vec<Option<ResolvedPhaseTransportState>>,
}

impl FlowField {
    /// Resolve declared mushy resistance from the ORIGINAL h/phase/source
    /// state, then apply one implicit native actual-density mechanical step.
    ///
    /// The supplied MushyPhaseLaw owns constituent densities, true liquid
    /// VOLUME fraction and inverse permeability. `resolve` is called exactly
    /// once for each admitted fluid cell with its original chart state. Bind
    /// that h through `resolve_phase_transport_state_for_volume_change` to
    /// retain current-density, viscosity and source receipts from one query.
    /// A kinematic-viscosity source is converted by that owner using the same
    /// chart-matching density. Frozen beta is then dynamic mu times K^-1.
    ///
    /// Both the complete returned mobility/law identity and source density must
    /// match this original state. Original M/rho_chart must fill its cell within
    /// the explicit volume budget. Native inertia remains the ACTUAL M/V; it
    /// is never replaced by the source rho or resistance-modified pressure metric.
    /// Unequal constituent densities and sensible chart expansion are admitted
    /// at this explicit original-state boundary without a new mixture rule.
    ///
    /// Exact constitutive solids refuse before any liquid property query:
    /// finite Darcy resistance is not an immobility constraint. Fixed grid
    /// solids retain the existing wall behavior. Fully liquid material has
    /// exact zero extra drag, and nonzero beta is never clipped or capped.
    ///
    /// This entry changes momentum only. Accepted physical resistance heat is
    /// retained in `evolution.evolution.report.cell_resistance_heat_j` for the
    /// enclosing caloric owner. It does not advance h/phase, prescribe thermal
    /// expansion, transport a phase inventory or reapply an existing pressure
    /// impulse. The volume target must belong to this fresh mechanical solve.
    /// Viscosity is frozen at the caller's declared source coordinates; no
    /// constitutive iteration or source shear-rate inference is introduced.
    ///
    /// # Errors
    /// Invalid source/state identity, unsupported phase support, source/chart
    /// density disagreement, occupied-volume defect, unrepresentable beta or
    /// any native balance/budget/cancellation refusal returns no candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn force_project_material_momentum_with_phase_resistance(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        law: &MushyPhaseLaw<'_>,
        input: PhaseResistanceMomentumInput<'_>,
        controls: PhaseResistanceMomentumControls,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ResolvedPhaseResistanceMomentumStep, ChtError> {
        self.force_project_material_momentum_with_phase_resistance_reservoirs(
            domain, faces, None, law, input, controls, gate, resolve,
        )
    }

    /// The same source-bound resistance mechanics with uniform plane pressures
    /// or a spatial reservoir profile. Pressure acts once through the original
    /// half-cell inertia/resistance and retains the native opening-work account.
    /// Mechanical reservoir pressure does not replace the source's fixed
    /// thermodynamic pressure or advance its enthalpy chart.
    ///
    /// # Errors
    /// Invalid pressure or any refusal of
    /// [`Self::force_project_material_momentum_with_phase_resistance`].
    #[allow(clippy::too_many_arguments)]
    pub fn force_project_material_momentum_with_phase_resistance_and_reservoir_pressure<'p>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        reservoir_pressure: impl Into<ReservoirPressure<'p>>,
        law: &MushyPhaseLaw<'_>,
        input: PhaseResistanceMomentumInput<'_>,
        controls: PhaseResistanceMomentumControls,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ResolvedPhaseResistanceMomentumStep, ChtError> {
        self.force_project_material_momentum_with_phase_resistance_reservoirs(
            domain,
            faces,
            Some(reservoir_pressure.into()),
            law,
            input,
            controls,
            gate,
            resolve,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn force_project_material_momentum_with_phase_resistance_reservoirs(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        reservoir_pressure: Option<ReservoirPressure<'_>>,
        law: &MushyPhaseLaw<'_>,
        input: PhaseResistanceMomentumInput<'_>,
        controls: PhaseResistanceMomentumControls,
        gate: &CancelGate,
        mut resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ResolvedPhaseResistanceMomentumStep, ChtError> {
        poll(gate)?;
        admit_feedback_controls(controls.mechanics)?;
        let n = domain.cell_count();
        if n > controls.mechanics.reconstruction.max_cells
            || n > controls.mechanics.forcing.max_cells
            || input.phase_states.len() != n
            || input.mechanics.mass_kg.len() != n
            || input.mechanics.momentum_kg_m_s.len() != n
            || input.mechanics.force_density_n_m3.len() != n
            || input.mechanics.target_outflow_m3_s.len() != n
        {
            return Err(refuse(
                "phase resistance needs bounded original cell inventories",
            ));
        }
        let volume_budget = controls.maximum_cell_volume_defect_m3;
        finite("flow.phase_resistance_volume_budget", volume_budget)?;
        if volume_budget < 0.0 {
            return Err(refuse("phase resistance volume budget must be nonnegative"));
        }
        let mut resistance = vec![0.0; n];
        let mut sources = Vec::with_capacity(n);
        for cell in 0..n {
            poll_every(gate, cell)?;
            if !domain.is_fluid(cell) {
                if input.phase_states[cell].is_some() || input.mechanics.mass_kg[cell] != 0.0 {
                    return Err(refuse("fixed grid solids cannot contain phase inventory"));
                }
                sources.push(None);
                continue;
            }
            let phase = input.phase_states[cell].ok_or_else(|| {
                refuse("phase resistance needs the original state in each fluid cell")
            })?;
            if phase.phase_curve_identity() != law.curve().identity()
                || phase.material_card_identity() != law.curve().material_card_identity()
            {
                return Err(refuse(
                    "phase resistance needs the original law's chart and material",
                ));
            }
            let expected = law
                .classify_at_specific_enthalpy(phase.specific_enthalpy_j_kg())
                .map_err(|error| refuse(&error.to_string()))?;
            if expected.phase() != phase {
                return Err(refuse(
                    "phase resistance needs the exact original h-resolved state",
                ));
            }
            // Phase topology is decided before querying any liquid-only source.
            let PhaseMobility::Mobile(mobile) = expected else {
                return Err(refuse(
                    "exact solid phase needs an explicit immobile momentum owner",
                ));
            };
            let mass = input.mechanics.mass_kg[cell];
            finite_positive("flow.phase_resistance_mass", mass)?;
            let occupied = mass / phase.bulk_density_kg_m3();
            finite_positive("flow.phase_resistance_occupied_volume", occupied)?;
            let defect = occupied - domain.volume(cell);
            finite("flow.phase_resistance_volume_defect", defect)?;
            if defect.abs() > volume_budget {
                return Err(refuse(
                    "original mass and phase density do not fill the cell",
                ));
            }
            let source = resolve(cell, &phase)?;
            poll(gate)?;
            if *source.mobility() != expected
                || *source.phase_state() != phase
                || source.source_density_kg_m3().to_bits() != phase.bulk_density_kg_m3().to_bits()
            {
                return Err(refuse(
                    "resistance source must match the exact original phase, law and density",
                ));
            }
            let transport = source.mobile_transport().ok_or_else(|| {
                refuse("mobile phase resistance needs its original fluid transport")
            })?;
            if transport.density_kg_m3().to_bits() != phase.bulk_density_kg_m3().to_bits() {
                return Err(refuse(
                    "resistance viscosity must retain the original current density",
                ));
            }
            let viscosity = transport.dynamic_viscosity_pa_s();
            finite_positive("flow.phase_resistance_dynamic_viscosity", viscosity)?;
            let inverse = mobile.inverse_permeability_per_m2;
            let beta = if inverse == 0.0 {
                0.0
            } else {
                viscosity * inverse
            };
            finite("flow.phase_resistance_beta", beta)?;
            if inverse != 0.0 && beta <= 0.0 {
                return Err(refuse("nonzero phase resistance is not representable"));
            }
            resistance[cell] = beta;
            sources.push(Some(source));
        }
        let evolution = match reservoir_pressure {
            Some(pressure) => self
                .force_project_material_momentum_with_resistance_and_reservoir_pressure(
                    domain,
                    faces,
                    pressure,
                    input.mechanics,
                    &resistance,
                    controls.mechanics,
                    gate,
                )?,
            None => self.force_project_material_momentum_with_resistance(
                domain,
                faces,
                input.mechanics,
                &resistance,
                controls.mechanics,
                gate,
            )?,
        };
        poll(gate)?;
        Ok(ResolvedPhaseResistanceMomentumStep {
            evolution,
            resistance_pa_s_m2: resistance,
            original_sources: sources,
        })
    }
}

#[cfg(test)]
mod tests;
