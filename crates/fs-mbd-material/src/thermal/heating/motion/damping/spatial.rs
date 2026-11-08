//! Full spatial viscous recoil with finite, source-bound thermal feedback.
//!
//! The original scalar material coefficient c(T) multiplies an explicitly
//! declared passive geometric shape. The existing pair thermal residual closes
//! endpoint temperature and damping heat together. Original contact conduction,
//! free flight and expansion at conserved world momenta retain their owners.
//! This is the existing first-order source/free-flight/configuration split:
//! independent time refinement remains necessary.

use super::*;
use fs_mbd::spatial_loads::internal_impulse::viscous::{
    PreparedViscousConnector, ViscousConnectorControls, ViscousConnectorError,
    ViscousConnectorReceipt,
};

mod shape;
pub use shape::SpatialResistanceShape;

/// Existing thermal/source solve plus native full-wrench and twist admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialDampingSolve {
    pub thermal: DampingSolve,
    pub connector: ViscousConnectorControls,
}

/// Native six-coordinate source receipt with the same finite thermal accounts.
pub type SpatialDampingAudit = DampingAudit<ViscousConnectorReceipt>;

/// Two continuing original thermal bodies linked by a full passive resistance.
///
/// This retains the original MaterialDamper's source, ordered surfaces, actual
/// material attachments, temperature convention and heat partition. The shape
/// supplies declared geometric direction/cross-coupling and rotary length;
/// it is expressed in the original first body's material axes, independent of
/// any principal-frame changes caused by thermal expansion.
///
/// A fixed shape/area is an explicit lumped interface approximation. No spring,
/// evolving bonded-pad dimensions, contact active set, interfacial mass or
/// reciprocal thermoelastic free energy is inferred.
#[derive(Clone)]
pub struct SpatialDampedRigidPair<'a> {
    pair: DampedRigidPair<'a>,
    shape: SpatialResistanceShape,
}

impl<'a> DampedRigidPair<'a> {
    /// Replace the axial constitutive geometry while retaining both complete
    /// accepted histories, original source and any bound thermal contact.
    pub fn with_spatial_resistance(
        self,
        shape: SpatialResistanceShape,
    ) -> SpatialDampedRigidPair<'a> {
        SpatialDampedRigidPair { pair: self, shape }
    }
}

impl<'a> SpatialDampedRigidPair<'a> {
    pub const fn bodies(&self) -> &[HeatedRigidBody<'a>; 2] {
        self.pair.bodies()
    }
    pub const fn damper(&self) -> &MaterialDamper<'a> {
        self.pair.damper()
    }
    pub const fn shape(&self) -> &SpatialResistanceShape {
        &self.shape
    }
    pub const fn thermal_contact(&self) -> Option<ThermalInterfacePoint<'a>> {
        self.pair.thermal_contact()
    }

    /// Recover the complete thermal pair and original geometric law.
    pub fn into_parts(self) -> (DampedRigidPair<'a>, SpatialResistanceShape) {
        (self.pair, self.shape)
    }

    /// Continue one bounded physical interval through the original thermal solve.
    ///
    /// Each trial restarts from the same accepted states, resolves the ORIGINAL
    /// endpoint c(T), compiles its full spatial law in the actual principal
    /// frame, and integrates one native passive source exactly at fixed poses.
    /// Its independent positive heat enters the named finite cells exactly once.
    ///
    /// The existing native free flight then advances both poses. Original source
    /// expansion refreshes the complete mass/inertia at conserved world momenta.
    /// Its parameter work and fixed-inertia integration error stay separate from
    /// physical dissipation, under the existing independent energy gates.
    /// The source's DC and complete prepared decay-rate range must be admitted.
    /// Newton convergence alone does not establish time accuracy.
    ///
    /// Native, source, phase, energy, work, budget or cancellation refusal
    /// returns no changed body, thermal history, law or clock.
    pub fn step_to(
        &self,
        target_s: f64,
        controls: SpatialDampingSolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, SpatialDampingAudit), DampingError> {
        let (pair, audit) = self.pair.step_with_impulse(
            target_s,
            controls.thermal,
            initial_energy_increments_j,
            cancelled,
            |source, dt, thermal, checkpoint| {
                self.impulse(source, dt, thermal, controls.connector, checkpoint)
            },
        )?;
        Ok((
            Self {
                pair,
                shape: self.shape.clone(),
            },
            audit,
        ))
    }

    fn impulse(
        &self,
        source: &ResolvedLinearInterface,
        dt: f64,
        thermal: DampingSolve,
        mut controls: ViscousConnectorControls,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Dissipation<ViscousConnectorReceipt>, DampingError> {
        checkpoint()?;
        let bodies = self.bodies();
        let states = bodies.each_ref().map(|body| body.state());
        let properties = bodies
            .each_ref()
            .map(|body| body.principal_frame().properties());
        let arms = [
            attachment_arm(&bodies[0], self.damper().attachments[0])?,
            attachment_arm(&bodies[1], self.damper().attachments[1])?,
        ];
        let law = self
            .shape
            .at_coefficient(
                source.damping_n_s_m(),
                bodies[0].principal_frame().body_from_principal().rotation(),
            )
            .map_err(|error| DampingError::Spatial(error.into()))?;
        let prepared = native_checked(checkpoint, |cancelled| {
            PreparedViscousConnector::new(
                states.map(|state| state.pose()),
                properties,
                arms,
                law,
                dt,
                cancelled,
            )
        })?;
        let rate = prepared.maximum_decay_rate_per_s();
        source.require_angular_frequency_rad_s(0.0)?;
        source.require_angular_frequency_rad_s(rate)?;
        if number(rate * dt)? > thermal.maximum_relaxation_step {
            return Err(DampingError::Budget(
                "refine physical timestep for the full spatial relaxation rate",
            ));
        }
        if !controls.impulse.work_tolerance_j.is_finite() || controls.impulse.work_tolerance_j < 0.0
        {
            return Err(DampingError::Invalid(
                "native spatial work tolerance must be finite and nonnegative",
            ));
        }
        // A more permissive native allowance cannot bypass the existing
        // complete thermal/mechanical work gate.
        controls.impulse.work_tolerance_j = controls
            .impulse
            .work_tolerance_j
            .min(thermal.total_energy_tolerance_j);
        let receipt = native_checked(checkpoint, |cancelled| {
            prepared.relax(states, controls, cancelled)
        })?;
        Ok(Dissipation {
            states_after: receipt.impulse.bodies.map(|body| body.state_after),
            heat_j: receipt.dissipated_energy_j,
            residual_j: receipt.kinetic_heat_residual_j,
            receipt,
        })
    }
}

fn native_checked<T>(
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    operation: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, ViscousConnectorError>,
) -> Result<T, DampingError> {
    let mut failure = None;
    let result = operation(&mut || match checkpoint() {
        Ok(()) => false,
        Err(error) => {
            failure = Some(error);
            true
        }
    });
    match failure {
        Some(error) => Err(error.into()),
        None => result.map_err(Into::into),
    }
}
