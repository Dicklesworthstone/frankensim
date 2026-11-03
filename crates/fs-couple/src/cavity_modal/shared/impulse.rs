//! Finite momentum changes through the retained cavity's original physical ports.
//!
//! Surface and neck inertia receive the input. Compression, viscous diffusion,
//! thermal/molecular memory and proper radiation states retain their histories;
//! their finite forces act during subsequent evolution, not during a zero-time
//! jump. Distributed ducts and algebraic constraints need their own reset law.

use super::{ModalCouplingError, SharedCavity, network_invalid};
use crate::render::schedule::force::coupled::impulse::ModalImpulseReceipt;
use fs_exec::CancelGate;

impl SharedCavity {
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
        let mut surface_rate = 0.0;
        for (index, state) in self
            .system
            .components()
            .iter()
            .flat_map(|component| component.states())
            .take(self.surface_modes)
            .enumerate()
        {
            if index % 64 == 0 {
                poll(gate)?;
            }
            surface_rate +=
                self.column[index] * (state.velocity_m_sqrt_kg_per_s + self.external[index]);
            if !surface_rate.is_finite() {
                return Err(network_invalid("cavity impulse surface flow overflow"));
            }
        }
        let neck = self.system.components()[self.neck_component].states()[0];
        let neck_flow = self.source.neck_volume_shape_m2_per_sqrt_kg()
            * (neck.velocity_m_sqrt_kg_per_s + self.external[self.surface_modes]);
        if !(surface_rate * self.source.neck_area_m2()).is_finite() || !neck_flow.is_finite() {
            return Err(network_invalid(
                "cavity impulse physical volume flow overflow",
            ));
        }
        self.system.apply_modal_impulse(&self.external, gate)
    }
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
