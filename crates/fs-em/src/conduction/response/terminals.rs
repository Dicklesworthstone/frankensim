//! Original multi-terminal controls inside the native scalar source closure.
//!
//! Each constitutive trial solves the complete grouped-electrode system once
//! at the actual drive. A floating terminal constrains its total inward current,
//! while allowing nonzero current transfer between its faces. No basis probe or
//! scaled accepted nonlinear field supplies the response.
//!
//! Geometry, bulk/contact source identities, relaxation, local constitutive
//! gates and cancellation stay with the original response owner. Terminal work
//! and reactions come from the same accepted native field.

use super::{
    ConductionContact, ConductionControls, ConductionResponseCandidate, ConductionResponseControls,
    ConductionResponseError, PlanarConductionResponse, ResolvedConductorState,
    WeightedAssemblyLimits, poll,
};

pub use crate::conduction::terminals::{
    ConductionTerminalControl, ConductionTerminalLimits, ConductionTerminalReport,
};

/// Held ideal voltage/current controls on the original ordered terminal sets.
/// Current controls specify total inward amperes; a zero current floats the
/// equipotential electrode. All quantities refer to the original physical
/// metric, including the authored angular span on an axisymmetric owner.
#[derive(Clone, Copy, Debug)]
pub struct ConductionTerminalDrive<'a> {
    pub terminals: &'a [&'a [usize]],
    pub control: ConductionTerminalControl<'a>,
    /// Bounds for the original grouped solve's terminal and workspace storage.
    pub limits: ConductionTerminalLimits,
}

impl PlanarConductionResponse {
    /// Close original scalar bulk and explicit contact laws at simultaneous
    /// terminal voltages, conserved total currents or mixed ideal controls.
    ///
    /// Every constitutive trial makes one actual-amplitude grouped field solve,
    /// without artificial probe states or a reusable nonlinear conductance
    /// matrix. Both the returned field and terminal report describe the same
    /// accepted state. Terminal voltages are relative to terminal zero; supplied
    /// power uses authored current controls and native voltage-source reactions.
    ///
    /// `maximum_field_solves` is an additional hard native-solve cap. The
    /// existing constitutive iteration cap and all source/field gates remain
    /// binding. Failed admission, budget exhaustion and cancellation leave
    /// this owner reusable.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_terminals_candidate(
        &self,
        drive: ConductionTerminalDrive<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        maximum_field_solves: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(ConductionResponseCandidate, ConductionTerminalReport), ConductionResponseError>
    {
        self.solve_terminals_with_materials_and_contacts_candidate(
            self.conductors().to_vec(),
            self.contacts().to_vec(),
            drive,
            assembly,
            field_controls,
            controls,
            maximum_field_solves,
            cancelled,
        )
    }

    /// The same terminal closure with refreshed original bulk/contact context.
    /// An outer thermal/phase owner supplies complete seeds; original source
    /// selections, fixed coordinates, geometry and electrical coordinate roles
    /// remain pinned. Only the driver report paired with the converged native
    /// field escapes the shared source iteration.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_terminals_with_materials_and_contacts_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        contacts: Vec<ConductionContact>,
        drive: ConductionTerminalDrive<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        maximum_field_solves: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(ConductionResponseCandidate, ConductionTerminalReport), ConductionResponseError>
    {
        poll(cancelled)?;
        controls.validate()?;
        if maximum_field_solves == 0 {
            return Err(ConductionResponseError::FieldBudget {
                completed: 0,
                required: 1,
                maximum: maximum_field_solves,
            });
        }
        let mut completed = 0;
        self.solve_coupled_with_materials_and_contacts_candidate(
            seeds,
            contacts,
            assembly,
            field_controls,
            controls,
            cancelled,
            |native, cancelled| {
                // The grouped native solve uses one actual operating state;
                // preflight before that solve, including the final admitted
                // iteration. No discarded constitutive work escapes this cap.
                if completed >= maximum_field_solves {
                    return Err(ConductionResponseError::FieldBudget {
                        completed,
                        required: 1,
                        maximum: maximum_field_solves,
                    });
                }
                let result = native
                    .solve_terminals(
                        drive.terminals,
                        drive.control,
                        drive.limits,
                        field_controls,
                        &mut || cancelled(),
                    )
                    .map_err(ConductionResponseError::from)?;
                completed += 1;
                Ok(result)
            },
        )
    }
}
