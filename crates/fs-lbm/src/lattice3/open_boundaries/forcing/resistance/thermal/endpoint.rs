//! Original scalar heat with evolving resistance and skeleton velocity.

use super::*;

/// Independent native heat and evolving-law momentum/work accounts.
/// No temporal work remainder or drag/solid heat partition is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
#[must_use]
pub struct ResistanceOpenEndpointThermalStep3 {
    pub momentum: ResistanceOpenEndpointExchange3,
    pub heat: ThermalOpenBalance3,
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> ResistanceOpenGrid3<L, Q> {
    /// Advance original D3Q7 heat and next-endpoint porous laws atomically.
    ///
    /// Opening heat transport uses accepted physical velocity and force history;
    /// its closing source half-kick uses the velocity resolved with `next_laws`.
    /// The original law and work budget of `advance_serial` apply unchanged.
    /// Final laws, force/drag history, both population fields, finite thermal
    /// walls, flux history, receipts and time publish only after every gate.
    ///
    /// Laws are supplied at the next physical time, not at a guessed thermal
    /// state. This API does not infer K(T), drag heat or skeleton inertia. Heat
    /// remains the original fixed-reference-capacity scalar with explicit source.
    pub fn advance_serial_with_thermal(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        next_laws: Vec<DarcyForchheimer3>,
        maximum_absolute_work_remainder: f64,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
    ) -> Result<ResistanceOpenEndpointThermalStep3, OpenThermalError3<L::Error>> {
        self.advance_with_thermal_poll(
            thermal,
            collisions,
            next_laws,
            maximum_absolute_work_remainder,
            heat_per_step,
            thermal_relaxation,
            None,
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Same endpoint transaction with the existing uniform or per-cell SPD
    /// thermal diffusivity, resolved in grid axes independently of porous axes.
    pub fn advance_serial_with_thermal_tensor_diffusion(
        &mut self,
        thermal: &mut ThermalOpenField3,
        cx: &Cx<'_>,
        collisions: Option<&[L::Collision]>,
        next_laws: Vec<DarcyForchheimer3>,
        maximum_absolute_work_remainder: f64,
        heat_per_step: &[f64],
        diffusion: DiffusivityField<'_, 3>,
    ) -> Result<ResistanceOpenEndpointThermalStep3, OpenThermalError3<L::Error>> {
        self.advance_with_thermal_poll(
            thermal,
            collisions,
            next_laws,
            maximum_absolute_work_remainder,
            heat_per_step,
            None,
            Some(diffusion),
            &mut || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    pub(super) fn advance_with_thermal_poll(
        &mut self,
        thermal: &mut ThermalOpenField3,
        collisions: Option<&[L::Collision]>,
        next_laws: Vec<DarcyForchheimer3>,
        maximum_absolute_work_remainder: f64,
        heat_per_step: &[f64],
        thermal_relaxation: Option<&[f64]>,
        diffusion: Option<DiffusivityField<'_, 3>>,
        poll: &mut impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<ResistanceOpenEndpointThermalStep3, OpenThermalError3<L::Error>> {
        let (_, momentum, step, heat) = self.stage_with_thermal_poll(
            thermal,
            collisions,
            heat_per_step,
            thermal_relaxation,
            diffusion,
            Some((&next_laws, maximum_absolute_work_remainder)),
            poll,
        )?;
        poll()?;
        self.publish_endpoint_step(step, momentum, next_laws);
        thermal.publish(step, heat);
        Ok(ResistanceOpenEndpointThermalStep3 { momentum, heat })
    }
}
