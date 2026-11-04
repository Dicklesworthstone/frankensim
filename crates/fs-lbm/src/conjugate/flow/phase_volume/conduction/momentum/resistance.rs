//! Original-source conduction and coupled phase/resistance heat-volume feedback.

use super::super::super::{
    ResistiveThermalVolumeMomentumControls, ResolvedResistiveThermalVolumeMomentumStep,
    ThermalVolumeMomentumInput,
};
use super::super::{ConductiveThermalVolumeControls, ConductiveThermalVolumeReport, refuse};
use crate::conjugate::flow::{FlowFace, FlowField, StaggeredBoundaryPressure, poll_every};
use crate::conjugate::{ChtError, VoxelDomain, finite, poll};
use fs_exec::CancelGate;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::phase::mushy::MushyPhaseLaw;
use fs_material::phase::transport::ResolvedPhaseTransportState;
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;
use fs_math::dd::Dd;

/// One set of original physical stage controls and the explicit conduction cap.
#[derive(Clone, Copy, Debug)]
pub struct ConductiveResistiveThermalVolumeMomentumControls {
    /// Shared pressure clock, material/volume, momentum and coupled heat gates.
    pub resistance: ResistiveThermalVolumeMomentumControls,
    pub max_conductive_faces: usize,
}

/// One accepted conductive, phase-dependent resistance and material endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductiveResistiveThermalVolumeMomentumStep {
    pub evolution: ResolvedResistiveThermalVolumeMomentumStep,
    /// Paired physical heat from the ORIGINAL conductivity states, in joules.
    /// These increments enter the caloric equation once, before any drag trials.
    pub conductive_enthalpy_increment_j: Vec<f64>,
    /// Original conductivity identity, explicit-step bound and the shared
    /// conduction/caloric/coupling defect account. Pressure work stays separate.
    pub conduction: ConductiveThermalVolumeReport,
}

impl FlowField {
    /// Conduct original-state heat once, then close phase-dependent resistance,
    /// physical drag heating and thermal-density volume flow on that same step.
    ///
    /// The supplied scalar conductivity is explicitly isotropic. Its original
    /// harmonic half-cell conductances and enthalpy-secant monotonicity bound
    /// are those of [`Self::advance_isotropic_conductive_thermal_volume`]. The
    /// paired heat joins authored net caloric increments once; bounded drag,
    /// phase and volume trials all restart from the original M/h/P inventory.
    /// The accepted native pressure field and mass routes carry momentum once.
    ///
    /// Conductivity and permeability must retain the same caloric chart and
    /// material card. Initial, incoming, accepted and accepted pre-transport
    /// caloric states INCLUDING drag heat must remain within that original
    /// conductivity source. Resolved transport preserves current phase density
    /// and material-source identity; mechanical inertia remains actual M/V.
    ///
    /// Conduction is an explicit original-state first-order split, evaluated
    /// once rather than at every nonlinear drag trial. External and fixed-solid
    /// thermal faces are insulated. Pressure/source work and numerical losses
    /// keep their existing accounts; they are not deposited as heat. Physical
    /// conduction assembly, native caloric transport and resistance coupling
    /// share both the original enthalpy and mechanical energy-defect budgets.
    ///
    /// # Errors
    /// Source/phase mismatch, unstable explicit conduction, source-bound exit,
    /// coupling exhaustion or any native budget/cancellation refusal returns
    /// no endpoint and changes none of the supplied physical inventories.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_conductive_phase_resistance_thermal_volume_momentum(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        conductivity: &PreparedPhaseConductivity<'_>,
        law: &MushyPhaseLaw<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ConductiveResistiveThermalVolumeMomentumControls,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ConductiveResistiveThermalVolumeMomentumStep, ChtError> {
        self.conduct_phase_resistance_thermal_volume_momentum(
            domain,
            faces,
            None,
            conductivity,
            law,
            input,
            controls,
            gate,
            resolve,
        )
    }

    /// The same coupled physical step with the existing staggered mechanical
    /// reservoir-pressure profile. Backflow material/velocity remain explicit;
    /// source thermodynamic pressure and enthalpy conventions are unchanged.
    ///
    /// # Errors
    /// An invalid pressure profile or any original conductive, material,
    /// resistance, thermal-volume, momentum or cancellation refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_conductive_phase_resistance_thermal_volume_momentum_with_pressure_profile(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: StaggeredBoundaryPressure<'_>,
        conductivity: &PreparedPhaseConductivity<'_>,
        law: &MushyPhaseLaw<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ConductiveResistiveThermalVolumeMomentumControls,
        gate: &CancelGate,
        resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ConductiveResistiveThermalVolumeMomentumStep, ChtError> {
        self.conduct_phase_resistance_thermal_volume_momentum(
            domain,
            faces,
            Some(profile),
            conductivity,
            law,
            input,
            controls,
            gate,
            resolve,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn conduct_phase_resistance_thermal_volume_momentum(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        profile: Option<StaggeredBoundaryPressure<'_>>,
        conductivity: &PreparedPhaseConductivity<'_>,
        law: &MushyPhaseLaw<'_>,
        input: ThermalVolumeMomentumInput<'_>,
        controls: ConductiveResistiveThermalVolumeMomentumControls,
        gate: &CancelGate,
        mut resolve: impl FnMut(
            usize,
            &EquilibriumPhaseState,
        ) -> Result<ResolvedPhaseTransportState, ChtError>,
    ) -> Result<ConductiveResistiveThermalVolumeMomentumStep, ChtError> {
        poll(gate)?;
        let curve = conductivity.caloric_curve();
        if curve.identity() != law.curve().identity()
            || curve.material_card_identity() != law.curve().material_card_identity()
        {
            return Err(refuse(
                "conductivity and phase resistance need the same original chart and material",
            ));
        }
        let accepted = self.conduct_then_advance(
            domain,
            conductivity,
            input.volume,
            ConductiveThermalVolumeControls {
                volume: controls.resistance.thermal.volume,
                max_conductive_faces: controls.max_conductive_faces,
            },
            gate,
            |volume| {
                let heated = ThermalVolumeMomentumInput { volume, ..input };
                let evolution = match profile {
                    Some(profile) => self
                        .advance_resolved_phase_resistance_thermal_volume_momentum_with_pressure_profile(
                            domain, faces, profile, law, heated, controls.resistance, gate, &mut resolve,
                        ),
                    None => self.advance_resolved_phase_resistance_thermal_volume_momentum(
                        domain, faces, law, heated, controls.resistance, gate, &mut resolve,
                    ),
                }?;
                // Transport can mix an out-of-source caloric intermediate back
                // into the endpoint range. Check the actual accepted split,
                // using the same combined increment and drag addition order.
                let bounds = conductivity.enthalpy_bounds_j_kg();
                for cell in 0..domain.cell_count() {
                    poll_every(gate, cell)?;
                    if !domain.is_fluid(cell) { continue; }
                    let old = volume.phase_states[cell]
                        .ok_or_else(|| refuse("conductive resistance needs its original phase"))?;
                    let base = volume.net_enthalpy_increment_j[cell];
                    let drag = evolution.evolution.caloric_resistance_heat_j[cell];
                    let combined = if drag == 0.0 { base } else { base + drag };
                    let advanced = curve
                        .advance_specific_energy(old, combined / volume.mass_kg[cell])
                        .map_err(|error| refuse(error.to_string()))?;
                    if !(bounds[0]..=bounds[1]).contains(&advanced.specific_enthalpy_j_kg()) {
                        return Err(refuse(
                            "drag-heated caloric state leaves the original conductivity source",
                        ));
                    }
                }
                poll(gate)?;
                Ok(evolution)
            },
            |step| {
                let native = &step.evolution;
                let coupling = &native.coupling;
                (
                    &native.thermal,
                    (Dd::from_f64(coupling.absolute_heat_coupling_defect_j)
                        + Dd::from_f64(coupling.absolute_constitutive_heat_defect_j))
                    .to_f64(),
                )
            },
        )?;
        let joint_defect = (Dd::from_f64(accepted.report.absolute_enthalpy_defect_j)
            + Dd::from_f64(
                accepted
                    .evolution
                    .evolution
                    .resistance
                    .absolute_energy_defect_j,
            ))
        .to_f64();
        finite("flow.conductive_resistance_energy_defect", joint_defect)?;
        if joint_defect
            > controls
                .resistance
                .thermal
                .mechanics
                .forcing
                .maximum_energy_defect_j
        {
            return Err(refuse(
                "conduction and resistance exceed the shared mechanical/caloric energy budget",
            ));
        }
        poll(gate)?;
        Ok(ConductiveResistiveThermalVolumeMomentumStep {
            evolution: accepted.evolution,
            conductive_enthalpy_increment_j: accepted.increments,
            conduction: accepted.report,
        })
    }
}

#[cfg(test)]
mod tests;
