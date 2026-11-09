//! Live remote pressure source in the original nonlinear valve/tract/bore solve.
//!
//! The source owns one explicit Driven endpoint of the upstream duct. Its
//! initialized compliance/inertance can emit and recover traveling-wave energy;
//! feedback reaches the moving valve only after the declared physical transit.
//! Existing source-free and local series-pressure paths retain their equations.
//! Source, valve, optional supply passage, downstream keys and both networks
//! publish together behind the original whole-system work/flow/storage gates.
//!
//! The upstream side must have exactly this one driven endpoint and no moving
//! key owner. Passive branches, wall/viscothermal loads, fractional transit and
//! radiation remain the original graph's responsibility. This is fixed linear
//! source storage, not an anatomical or thermochemical reservoir model.

use super::{TwoSidedAperture, TwoSidedDrive, TwoSidedFrame, invalid, poll, wave};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::network::boundaries::BoundaryDrive;
use crate::upstream_airway::StoredPressureSource;
use fs_exec::CancelGate;

impl TwoSidedAperture {
    /// Attach original initialized source storage to the explicit remote
    /// upstream Driven node of a cold valve/duct assembly.
    ///
    /// The source wave reference must equal its actual adjacent duct impedance,
    /// and its clock must equal the valve/network clock. The node must already
    /// be authored as Driven: no passive or fitted termination is overwritten.
    /// Initial source energy joins the existing complete energy bound.
    pub fn with_stored_pressure_source(
        mut self,
        node: usize,
        source: StoredPressureSource,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        if self.valve.accepted_steps() != 0
            || self.stored_source.is_some()
            || self.sides[0].boundaries.is_some()
            || self
                .sides
                .iter()
                .any(|side| side.network.accepted_steps() != 0)
            || !self.sides[0]
                .network
                .driven_nodes()
                .eq(core::iter::once(node))
            || source.impedance().time_step_s().to_bits() != self.valve.time_step_s().to_bits()
        {
            return Err(invalid(
                "stored source requires one cold upstream driven endpoint on the shared clock",
            ));
        }
        let endpoint = self.sides[0].network.driven_boundary(node).map_err(wave)?;
        if endpoint.impedance_pa_s_m3.to_bits() != source.impedance().port_impedance().to_bits() {
            return Err(invalid(
                "stored source wave reference must match its actual adjacent upstream section",
            ));
        }
        self.stored_source = Some((node, source));
        let energy = self.stored_energy_j();
        if !energy.is_finite() || energy > self.limits.maximum_energy_j {
            return Err(invalid(
                "stored source exceeds the complete initial valve/duct storage bound",
            ));
        }
        poll(gate)?;
        Ok(self)
    }

    /// Original endpoint identity and accepted source state. No independent
    /// stepping access is exposed through the coupled owner.
    #[must_use]
    pub fn stored_pressure_source(&self) -> Option<(usize, &StoredPressureSource)> {
        self.stored_source
            .as_ref()
            .map(|(node, source)| (*node, source))
    }

    /// Advance the original nonlinear valve against both retained ducts and
    /// its remote stored source. Remote pressure is the ideal pressure source
    /// in series with that storage; the existing drive remains the distinct
    /// local series source and mechanical force input.
    pub fn step_with_stored_source_under_gate(
        &mut self,
        drive: TwoSidedDrive,
        remote_pressure_pa: f64,
        gate: &CancelGate,
    ) -> Result<TwoSidedFrame, AcousticRealizeError> {
        self.step_coupled_under_gate(drive, [&[], &[]], [0.0; 2], Some(remote_pressure_pa), gate)
    }

    /// The same source/valve transaction retaining downstream moving keys and
    /// original spatial blade-pressure ports. Upstream boundary drives must be
    /// empty because that side's sole driven endpoint belongs to the source.
    pub fn step_with_stored_source_and_boundaries_under_gate(
        &mut self,
        drive: TwoSidedDrive,
        remote_pressure_pa: f64,
        boundaries: [&[BoundaryDrive]; 2],
        blade_pressures_pa: [f64; 2],
        gate: &CancelGate,
    ) -> Result<TwoSidedFrame, AcousticRealizeError> {
        self.step_coupled_under_gate(
            drive,
            boundaries,
            blade_pressures_pa,
            Some(remote_pressure_pa),
            gate,
        )
    }
}
