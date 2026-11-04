//! Prescribed source-bearing molecular boundary states shared by both steppers.

use super::super::Cell;
use super::*;

#[derive(Clone, Copy)]
pub(in crate::multicomponent) struct ReservoirTrace<'a> {
    pub(in crate::multicomponent) cell: Cell<'a>,
    pub(in crate::multicomponent) primitive: Primitive,
    pub(in crate::multicomponent) heat_capacity_v: f64,
}
pub(in crate::multicomponent) type ReservoirEnds<'a> = [Option<ReservoirTrace<'a>>; 2];

pub(in crate::multicomponent) fn admit<'a>(
    gas: &MixtureGas,
    boundaries: Boundaries,
    reservoirs: [Option<Cell<'a>>; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ReservoirEnds<'a>, DissipationError> {
    poll(cancelled)?;
    if boundaries == Boundaries::Periodic {
        return Err(DissipationError::Invalid(
            "prescribed molecular reservoirs require physical ends",
        ));
    }
    let mut admitted = [None; 2];
    for (target, reservoir) in admitted.iter_mut().zip(reservoirs) {
        if let Some(cell) = reservoir {
            poll(cancelled)?;
            if cell.momentum != 0.0 {
                return Err(DissipationError::Invalid(
                    "stationary molecular reservoir has nonzero momentum",
                ));
            }
            let primitive = gas.primitive(cell, &mut *cancelled)?;
            let heat_capacity_v =
                evolve::heat_capacity_density(gas, cell.partial_densities, primitive)?;
            *target = Some(ReservoirTrace {
                cell,
                primitive,
                heat_capacity_v,
            });
        }
    }
    poll(cancelled)?;
    Ok(admitted)
}

impl<L: TransportLaw> DissipationOperator<L> {
    /// Advance molecular transport against prescribed stationary end baths.
    ///
    /// `Some` replaces that end's molecular condition for this interval;
    /// `None` retains it, and two absent baths execute [`Self::advance_for`]
    /// exactly. Reservoir cells use the same canonical species order and
    /// per-species energy gauges as this gas, just as for
    /// [`MixtureEuler1d::step_with_reservoirs`]. Borrowed cells carry no source
    /// identity token: the caller owns that correspondence. Their full caloric
    /// states are admitted through the original gas source before a trial.
    /// Periodic geometry and nonzero reservoir momentum are refused.
    ///
    /// Every SSP stage reevaluates Maxwell--Stefan, Fourier and stress fluxes
    /// between the evolving end cell and fixed bath over the physical half-cell
    /// distance. The shared spatial reservoir kernel uses bath enthalpies for
    /// species energy and bath velocity for traction work. No advective flux
    /// is added; a composing owner must apply its corresponding pressure
    /// reservoirs separately. Signed species, momentum and energy transfers
    /// enter the existing complete boundary balance.
    ///
    /// Reservoirs are not depleted, retained or allocated as new workspace.
    /// Publication remains atomic for the complete duration, and the flow
    /// clock is unchanged. Moving boundaries require a separate closure.
    pub fn advance_with_reservoirs_for(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: DissipationControl,
        reservoirs: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DissipationReport, DissipationError> {
        if reservoirs.iter().all(Option::is_none) {
            return self.advance_for(flow, duration, controls, cancelled);
        }
        let reservoirs = admit(&self.gas, self.boundaries, reservoirs, &mut cancelled)?;
        evolve::advance(self, flow, duration, controls, reservoirs, &mut cancelled)
    }

    /// The same prescribed molecular baths through the original RKL2 owner.
    ///
    /// Source admission, physical half-cell gradients, bath energy effort and
    /// transaction semantics match [`Self::advance_with_reservoirs_for`]. All
    /// intermediate and endpoint source evaluations use the shared face owner;
    /// boundary transfers follow the same signed recurrence as the field.
    /// Original macrostep accuracy and work limits remain in force.
    pub fn advance_stabilized_with_reservoirs_for(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        reservoirs: [Option<Cell<'_>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StabilizedDissipationReport, DissipationError> {
        self.advance_stabilized_with_reservoirs_budget(
            flow,
            duration,
            controls,
            reservoirs,
            &mut 0,
            &mut cancelled,
        )
    }

    // Coupled retries and full intervals retain every actual RHS evaluation,
    // including work discarded when another owner shortens the shared step.
    pub(in crate::multicomponent) fn advance_stabilized_with_reservoirs_budget(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: StabilizedDissipationControl,
        reservoirs: [Option<Cell<'_>>; 2],
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StabilizedDissipationReport, DissipationError> {
        if reservoirs.iter().all(Option::is_none) {
            return self.advance_stabilized_with_budget(
                flow,
                duration,
                controls,
                evaluations,
                cancelled,
            );
        }
        let reservoirs = admit(&self.gas, self.boundaries, reservoirs, cancelled)?;
        stabilized::advance(
            self,
            flow,
            duration,
            controls,
            reservoirs,
            evaluations,
            cancelled,
        )
    }
}
