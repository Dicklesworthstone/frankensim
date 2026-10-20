//! Original moving-boundary owners on either physical side of the valve.
use super::{TwoSidedAperture, TwoSidedDrive, TwoSidedFrame, ValveSide, invalid, poll};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::network::{
    boundaries::{BoundaryAperture, BoundaryApertures, BoundaryDrive},
    keyed,
};
use fs_exec::CancelGate;

impl TwoSidedAperture {
    /// Attach each declared driven node on one cold acoustic side to its
    /// original material/contact aperture. The actual lowered chimney supplies
    /// its local gas, impedance and integer/fractional transit, preserving bore
    /// walls, distributed losses, supply history and any MIMO radiation group.
    ///
    /// Every key shares the original valve clock and lifetime horizon. Keys on
    /// the other side remain independently attached; an unselected side retains
    /// all of its authored loads. Runtime drives are explicit for both sides.
    /// This is the existing moving-slit approximation, not a resolved circular
    /// pad, pivot/lever mechanism, or radiation model for the moving opening.
    pub fn with_boundary_apertures(
        mut self,
        side: ValveSide,
        inputs: Vec<BoundaryAperture>,
        max_apertures: usize,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        if self.valve.accepted_steps() != 0 || self.sides[side.index()].boundaries.is_some() {
            return Err(invalid("moving boundaries require one cold attachment per acoustic side"));
        }
        keyed::admit_horizon(&inputs, max_apertures, self.valve.max_steps())?;
        let duct = &mut self.sides[side.index()];
        duct.boundaries = Some(BoundaryApertures::new(
            &duct.graph, &duct.gases, &duct.network, inputs, max_apertures,
        )?);
        let energy = self.stored_energy_j();
        if !energy.is_finite() || energy > self.limits.maximum_energy_j {
            return Err(invalid("two-sided moving boundaries exceed the complete initial storage limit"));
        }
        poll(gate)?;
        Ok(self)
    }

    /// Original fixed nodes, local gas, source plates, contact/material state and
    /// last accepted observations for this side. No independent stepping access.
    #[must_use]
    pub fn boundary_apertures(&self, side: ValveSide) -> Option<&BoundaryApertures> {
        self.sides[side.index()].boundaries.as_ref()
    }

    /// Uncancelled one-sample convenience using the same paired transaction.
    pub fn step_with_boundaries(
        &mut self,
        drive: TwoSidedDrive,
        boundaries: [&[BoundaryDrive]; 2],
    ) -> Result<TwoSidedFrame, AcousticRealizeError> {
        self.step_with_boundaries_under_gate(drive, boundaries, &CancelGate::new_clock_free())
    }
}
