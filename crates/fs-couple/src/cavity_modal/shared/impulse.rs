//! Finite momentum changes through the retained cavity's original physical ports.
//!
//! Surface and neck inertia receive the input. Compression, viscous diffusion,
//! thermal/molecular memory and proper radiation states retain their histories;
//! their finite forces act during subsequent evolution, not during a zero-time
//! jump. Distributed ducts and algebraic constraints need their own reset law.

use super::{ModalCouplingError, SharedCavity, network_invalid};
use crate::render::schedule::force::coupled::impulse::ModalImpulseReceipt;
use crate::render::schedule::force::coupled::impulse::impact::frictional::{
    ModalFrictionalImpactError, ModalFrictionalImpactReceipt, PreparedModalFrictionalImpacts,
};
use crate::render::schedule::force::coupled::impulse::impact::multiple::{
    ModalMultiImpactError, ModalMultiImpactReceipt, PreparedModalImpacts,
};
use crate::render::schedule::force::coupled::impulse::impact::{
    ModalImpactError, ModalImpactReceipt, PreparedModalImpact,
};
use fs_exec::CancelGate;

impl SharedCavity {
    /// Apply one simultaneous Coulomb event to the original physical surfaces.
    /// Normal and both tangent traces address the same surface pair at each
    /// contact. Zero tangents retain the original native locked-direction
    /// classification; they are not replaced by invented inertia or friction.
    ///
    /// Actual flow and every modal/contact gate precede one publication.
    /// Neck/viscous auxiliary velocities, compression, constitutive/radiation
    /// histories, q and clock remain continuous. Prepare against Self::network;
    /// acoustic auxiliaries cannot become contact partners. Total impact loss
    /// carries no automatic local heat/pressure or shock interpretation.
    pub fn apply_surface_frictional_impacts(
        &mut self,
        impacts: &mut PreparedModalFrictionalImpacts,
        gate: &CancelGate,
    ) -> Result<ModalFrictionalImpactReceipt, ModalFrictionalImpactError> {
        poll(gate)?;
        if impacts
            .component_pairs()
            .flatten()
            .any(|index| index >= self.neck_component)
        {
            return Err(network_invalid(
                "cavity Coulomb impacts require original physical surfaces",
            )
            .into());
        }
        let column = &self.column;
        let surface_modes = self.surface_modes;
        let neck_shape = self.source.neck_volume_shape_m2_per_sqrt_kg();
        let neck_area = self.source.neck_area_m2();
        self.system
            .apply_frictional_impacts_with_admission(impacts, gate, |components| {
                validate_volume_flows(
                    components
                        .iter()
                        .flat_map(|model| model.states())
                        .map(|state| state.velocity_m_sqrt_kg_per_s),
                    column,
                    surface_modes,
                    neck_shape,
                    neck_area,
                    gate,
                )
            })
    }

    /// Resolve one simultaneous frictionless contact set on the continuing
    /// physical surfaces. Every row shares the complete native Delassus solve;
    /// a free striker may have zero cavity-volume projection while retaining
    /// its finite inertia. Prepare against [`Self::network`].
    ///
    /// All surface states and the actual physical volume flows are admitted
    /// before one publication. Neck/viscous auxiliary velocities, compression,
    /// material/radiation histories, displacements and the clock stay unchanged.
    /// Contacts attached to acoustic auxiliaries refuse. Redundant rows retain
    /// the original solver's unresolved local load allocation; no local heat,
    /// shock, contact search or persistent-contact model is inferred.
    pub fn apply_surface_impacts(
        &mut self,
        impacts: &mut PreparedModalImpacts,
        gate: &CancelGate,
    ) -> Result<ModalMultiImpactReceipt, ModalMultiImpactError> {
        poll(gate)?;
        if impacts
            .component_pairs()
            .flatten()
            .any(|index| index >= self.neck_component)
        {
            return Err(network_invalid(
                "cavity impacts must attach only to original physical surfaces",
            )
            .into());
        }
        let column = &self.column;
        let surface_modes = self.surface_modes;
        let neck_shape = self.source.neck_volume_shape_m2_per_sqrt_kg();
        let neck_area = self.source.neck_area_m2();
        self.system
            .apply_normal_impacts_with_admission(impacts, gate, |components| {
                validate_volume_flows(
                    components
                        .iter()
                        .flat_map(|model| model.states())
                        .map(|state| state.velocity_m_sqrt_kg_per_s),
                    column,
                    surface_modes,
                    neck_shape,
                    neck_area,
                    gate,
                )
            })
    }

    /// Resolve one frictionless impact between the original moving surfaces.
    /// Prepare its physical attachments against [`Self::network`]. A free-mass
    /// striker can be an ordinary surface with zero cavity-volume projection;
    /// its finite inertia then participates in the same reciprocal impact.
    ///
    /// Surface motion changes immediately. The finite neck, viscous auxiliaries,
    /// compression, thermal/molecular and radiation histories, volume and clock
    /// stay continuous and respond through subsequent native cavity evolution.
    /// Actual physical flow and every original native jump gate succeed before
    /// either surface publishes. Acoustic auxiliary coordinates cannot be used
    /// as mechanical contact partners. No shock, contact locator or simultaneous
    /// impact law is implied; the original scalar Newton-law limits apply.
    pub fn apply_surface_impact(
        &mut self,
        impact: &PreparedModalImpact,
        gate: &CancelGate,
    ) -> Result<ModalImpactReceipt, ModalImpactError> {
        poll(gate)?;
        if impact
            .components()
            .iter()
            .any(|&index| index >= self.neck_component)
        {
            return Err(network_invalid(
                "cavity impact must attach only to original physical surfaces",
            )
            .into());
        }
        let column = &self.column;
        let surface_modes = self.surface_modes;
        let neck_shape = self.source.neck_volume_shape_m2_per_sqrt_kg();
        let neck_area = self.source.neck_area_m2();
        self.system
            .apply_normal_impact_with_admission(impact, gate, |components| {
                validate_volume_flows(
                    components
                        .iter()
                        .flat_map(|model| model.states())
                        .map(|state| state.velocity_m_sqrt_kg_per_s),
                    column,
                    surface_modes,
                    neck_shape,
                    neck_area,
                    gate,
                )
            })
    }

    /// Apply simultaneous surface momentum inputs and an exterior neck pressure
    /// impulse without restarting the coupled cavity or advancing its clock.
    ///
    /// `surface_impulses` contains signed generalized momenta [N s / sqrt(kg)]
    /// in the original surface/mode order, excluding neck and auxiliary modes.
    /// A physical point impulse uses the same supplied structural port column
    /// as a force. Inputs must already sum all loads on each coordinate.
    /// Positive exterior neck pressure impulse [Pa s] acts inward and gives
    /// the retained neck mass momentum `-neck_area * pressure_impulse` [N s].
    /// Internal cavity pressure is not an additional prescribed input.
    ///
    /// Only physical surface/neck velocities jump. Positions, cavity volume,
    /// compression pressure, viscous auxiliary velocities, thermal/molecular
    /// memory, radiation history and the preceding frame/port window remain
    /// unchanged. The receipt reports supplied work for this separate event;
    /// no damping heat or elapsed-time port window is invented.
    ///
    /// Original momentum, state, energy, pressure and connection-force limits
    /// apply, and resulting physical volume flows must remain representable.
    /// Cancellation or refusal publishes nothing. The retained force scratch
    /// is reused; no allocation or sample-dependent force pulse is introduced.
    pub fn apply_modal_impulse(
        &mut self,
        surface_impulses: &[f64],
        exterior_neck_pressure_impulse_pa_s: f64,
        gate: &CancelGate,
    ) -> Result<ModalImpulseReceipt, ModalCouplingError> {
        poll(gate)?;
        if surface_impulses.len() != self.surface_modes
            || surface_impulses.iter().any(|value| !value.is_finite())
        {
            return Err(network_invalid(
                "cavity impulse requires finite momentum for every original surface mode",
            ));
        }
        let neck = pressure_momentum(
            self.source.neck_volume_shape_m2_per_sqrt_kg(),
            exterior_neck_pressure_impulse_pa_s,
        )?;
        self.external.fill(0.0);
        self.external[..self.surface_modes].copy_from_slice(surface_impulses);
        self.external[self.surface_modes] = neck;
        self.apply_prepared_impulse(gate)
    }

    /// Apply one uniform exterior pressure impulse [Pa s] per original moving
    /// surface, plus the independently prescribed exterior neck impulse.
    ///
    /// The original signed outward volume shapes map `J_k = -A_k integral(p dt)`
    /// into modal momentum, exactly as pressure force maps into the same basis.
    /// Positive pressure acts inward; opposite shape signs remain opposite.
    /// Work is `-integral(p dt) * midpoint(outward volume rate)` [J], including
    /// signed energy extraction when the applied impulse opposes existing motion.
    /// All surfaces and the finite neck participate in one native transaction.
    ///
    /// This is a prescribed pressure integral on the admitted finite-inertia,
    /// small-signal cavity image. It does not impose an internal pressure jump,
    /// model a shock front, or reset an algebraic/distributed acoustic boundary.
    /// A zero volume shape has no response within the retained basis. Missing
    /// modal participation is not inferred from surface area or pressure alone.
    /// The state/history, budget and cancellation guarantees are those of
    /// [`Self::apply_modal_impulse`].
    pub fn apply_pressure_impulses(
        &mut self,
        exterior_surface_pressure_impulses_pa_s: &[f64],
        exterior_neck_pressure_impulse_pa_s: f64,
        gate: &CancelGate,
    ) -> Result<ModalImpulseReceipt, ModalCouplingError> {
        poll(gate)?;
        if exterior_surface_pressure_impulses_pa_s.len() != self.neck_component
            || exterior_surface_pressure_impulses_pa_s
                .iter()
                .any(|value| !value.is_finite())
        {
            return Err(network_invalid(
                "cavity pressure impulse requires one finite value per original surface",
            ));
        }
        let neck = pressure_momentum(
            self.source.neck_volume_shape_m2_per_sqrt_kg(),
            exterior_neck_pressure_impulse_pa_s,
        )?;
        self.external.fill(0.0);
        let mut offset = 0;
        for (surface, pressure) in self
            .system
            .components()
            .iter()
            .take(self.neck_component)
            .zip(exterior_surface_pressure_impulses_pa_s)
        {
            poll(gate)?;
            for _ in surface.modes() {
                if offset % 64 == 0 {
                    poll(gate)?;
                }
                self.external[offset] =
                    pressure_momentum(self.surface_volume_shapes[offset], *pressure)?;
                offset += 1;
            }
        }
        self.external[self.surface_modes] = neck;
        self.apply_prepared_impulse(gate)
    }

    fn apply_prepared_impulse(
        &mut self,
        gate: &CancelGate,
    ) -> Result<ModalImpulseReceipt, ModalCouplingError> {
        // Validate the actual observation's physical flow before native state
        // publication. Reusing its arithmetic also catches a large normalized
        // intermediate that cannot be represented even if a later area scale
        // would make the physical result small.
        validate_volume_flows(
            self.system
                .components()
                .iter()
                .flat_map(|component| component.states())
                .enumerate()
                .map(|(index, state)| state.velocity_m_sqrt_kg_per_s + self.external[index]),
            &self.column,
            self.surface_modes,
            self.source.neck_volume_shape_m2_per_sqrt_kg(),
            self.source.neck_area_m2(),
            gate,
        )?;
        self.system.apply_modal_impulse(&self.external, gate)
    }
}

fn validate_volume_flows(
    mut velocities: impl Iterator<Item = f64>,
    column: &[f64],
    surface_modes: usize,
    neck_shape: f64,
    neck_area: f64,
    gate: &CancelGate,
) -> Result<(), ModalCouplingError> {
    let mut surface_rate = 0.0;
    for (index, &shape) in column.iter().take(surface_modes).enumerate() {
        if index % 64 == 0 {
            poll(gate)?;
        }
        let velocity = velocities
            .next()
            .ok_or_else(|| network_invalid("cavity impulse surface layout changed"))?;
        surface_rate += shape * velocity;
        if !surface_rate.is_finite() {
            return Err(network_invalid("cavity impulse surface flow overflow"));
        }
    }
    let neck_velocity = velocities
        .next()
        .ok_or_else(|| network_invalid("cavity impulse neck layout changed"))?;
    if !(surface_rate * neck_area).is_finite() || !(neck_shape * neck_velocity).is_finite() {
        return Err(network_invalid(
            "cavity impulse physical volume flow overflow",
        ));
    }
    Ok(())
}

fn pressure_momentum(shape: f64, pressure_impulse: f64) -> Result<f64, ModalCouplingError> {
    let momentum = -shape * pressure_impulse;
    if !pressure_impulse.is_finite()
        || !momentum.is_finite()
        || (shape != 0.0 && pressure_impulse != 0.0 && momentum == 0.0)
    {
        return Err(network_invalid(
            "cavity pressure impulse projection is nonfinite or loses nonzero momentum",
        ));
    }
    Ok(momentum)
}

fn poll(gate: &CancelGate) -> Result<(), ModalCouplingError> {
    if gate.is_requested() {
        Err(ModalCouplingError::Cancelled)
    } else {
        Ok(())
    }
}
