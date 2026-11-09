//! Reusable terminal response on the original swept conductivity operator.
//!
//! Exactly N-1 native field solves prepare I=G V. Simultaneous voltage,
//! total-current and mixed drives reuse the original factors and recover
//! source-admitted local fields and Joule heat, including all cross terms.

use super::{
    AxisymmetricConduction, AxisymmetricConductionField, ConductionControls, ConductionError,
    PlanarOhmicResponse, ResolvedConductorState, Result, admit_terminal_traces, poll,
};
use crate::conduction::matrix::{PlanarConductanceMatrix, PlanarMixedConductance};

/// The same explicit terminal-count, solve, storage and dense-work bounds used
/// by the native scalar response. Conductance is in siemens for the supplied
/// angular span; the probe voltage must lie inside the retained material domain.
pub use crate::conduction::matrix::PlanarConductanceControls as AxisymmetricConductanceControls;
/// Bounds for the retained mixed-current factor on a fixed terminal partition.
pub use crate::conduction::matrix::PlanarMixedConductanceControls as AxisymmetricMixedConductanceControls;

/// Immutable source/geometry response for complete physical electrodes.
/// The matrix relates total inward terminal currents [A] to voltages [V].
/// Reprepare after geometry, material, contact or temperature changes.
pub struct AxisymmetricConductanceMatrix<'m, C = ResolvedConductorState> {
    problem: &'m AxisymmetricConduction<C>,
    response: PlanarConductanceMatrix<'m, C>,
}

/// A reusable arbitrary voltage/current partition on the same response.
/// Zero total current floats an equipotential electrode; its local surface
/// currents can remain nonzero, transferring current between its faces.
pub struct AxisymmetricMixedConductance<'s, 'm, C = ResolvedConductorState> {
    problem: &'m AxisymmetricConduction<C>,
    response: PlanarMixedConductance<'s, 'm, C>,
}

impl<C: PlanarOhmicResponse> AxisymmetricConduction<C> {
    /// Prepare a reciprocal conductance matrix with N-1 native basis solves.
    ///
    /// Terminals are disjoint unions of complete positive-area surfaces of
    /// revolution; the zero-area symmetry axis cannot be an electrode. A
    /// terminal can join disconnected mesh pieces as one ideal equipotential
    /// conductor. All terminals must share a conducting connection through
    /// the bulk, admitted finite contacts or a shared electrode.
    ///
    /// The native energy Gram matrix, reaction audit and equilibrated dense
    /// factors use the exact radial measure and original conductivity law.
    /// No lumped-resistor fit or uniform terminal current density is assumed.
    pub fn conductance_matrix(
        &self,
        terminals: &[&[usize]],
        controls: AxisymmetricConductanceControls,
        field_controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductanceMatrix<'_, C>> {
        poll(cancelled)?;
        if terminals.len() > controls.maximum_terminals
            || terminals.len().saturating_sub(1) > controls.maximum_field_solves
        {
            return Err(ConductionError::Invalid(
                "axisymmetric terminal count exceeds the response bounds",
            ));
        }
        admit_terminal_traces(self.mesh(), terminals, cancelled)?;
        let response =
            self.core
                .conductance_matrix(terminals, controls, field_controls, cancelled)?;
        Ok(AxisymmetricConductanceMatrix {
            problem: self,
            response,
        })
    }
}

impl<'m, C: PlanarOhmicResponse> AxisymmetricConductanceMatrix<'m, C> {
    pub fn problem(&self) -> &AxisymmetricConduction<C> {
        self.problem
    }

    /// Canonical vertex sets in the authored terminal order.
    pub fn terminal_vertices(&self) -> &[Vec<usize>] {
        self.response.terminal_vertices()
    }

    /// Full row-major conductance [S], including terminal zero: I=G V.
    /// Current and dissipation scale with the supplied angular span.
    pub fn conductance_s(&self) -> &[f64] {
        self.response.conductance_s()
    }

    /// Number of native field solves used during preparation, exactly N-1.
    pub fn field_solves(&self) -> usize {
        self.response.field_solves()
    }

    pub fn maximum_coefficient_residual_s(&self) -> f64 {
        self.response.maximum_coefficient_residual_s()
    }

    /// Recover actual simultaneous voltages without another field solve.
    /// Common voltage bias is removed before current/field/heat arithmetic;
    /// returned potentials use the original terminal-zero reporting gauge.
    /// Every actual field is rechecked against the retained source domain.
    pub fn field_at_voltages(
        &self,
        voltages_v: &[f64],
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionField<'m, C>> {
        let field = self.response.field_at_voltages(
            voltages_v,
            current_tolerance_a,
            power_tolerance_w,
            cancelled,
        )?;
        Ok(AxisymmetricConductionField {
            problem: self.problem,
            field,
        })
    }

    /// Prescribe conserved total currents into every equipotential electrode.
    /// A zero current is a floating electrode, not an insulating surface.
    /// The existing grounded factor is reused, followed by independent
    /// original-geometry current, source-work and local Joule-heat checks.
    pub fn field_at_currents(
        &self,
        currents_a: &[f64],
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionField<'m, C>> {
        let field = self.response.field_at_currents(
            currents_a,
            current_tolerance_a,
            power_tolerance_w,
            cancelled,
        )?;
        Ok(AxisymmetricConductionField {
            problem: self.problem,
            field,
        })
    }

    /// Factor an arbitrary nonempty proper subset of current electrodes once.
    /// Current input order follows this list; voltage input order follows its
    /// complement in original terminal order. A voltage electrode fixes gauge.
    pub fn prepare_current_control(
        &self,
        current_terminals: &[usize],
        controls: AxisymmetricMixedConductanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricMixedConductance<'_, 'm, C>> {
        let response =
            self.response
                .prepare_current_control(current_terminals, controls, cancelled)?;
        Ok(AxisymmetricMixedConductance {
            problem: self.problem,
            response,
        })
    }
}

impl<'m, C: PlanarOhmicResponse> AxisymmetricMixedConductance<'_, 'm, C> {
    pub fn voltage_terminals(&self) -> &[usize] {
        self.response.voltage_terminals()
    }

    pub fn current_terminals(&self) -> &[usize] {
        self.response.current_terminals()
    }

    /// Recover simultaneous voltage sources and current/floating electrodes
    /// with the retained factor. Voltage-source reactions supply the remaining
    /// current and work; returned potentials use original terminal zero.
    pub fn field_at_controls(
        &self,
        voltages_v: &[f64],
        currents_a: &[f64],
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionField<'m, C>> {
        let field = self.response.field_at_controls(
            voltages_v,
            currents_a,
            current_tolerance_a,
            power_tolerance_w,
            cancelled,
        )?;
        Ok(AxisymmetricConductionField {
            problem: self.problem,
            field,
        })
    }
}
