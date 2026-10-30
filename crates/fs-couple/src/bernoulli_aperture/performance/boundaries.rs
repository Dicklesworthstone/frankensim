//! Shared sample-addressed controls for retained physical graph boundaries.
//!
//! Values are actual reservoir pressures and work-conjugate opening forces.
//! This layer holds admitted values; the original aperture/network owners retain
//! the motion, closure, material memory, trapped pressure and wave propagation.
use super::{AperturePerformance, CoupledAperture, GestureCompileError, invalid};
use crate::bernoulli_aperture::network::boundaries::{
    BoundaryApertures, BoundaryDrive, BoundaryNetworkFrame,
};

/// One complete boundary drive applied before its named physical sample.
/// It remains held until another assignment names the same node. Zero force
/// releases the actuator; it does not erase any mechanical or acoustic state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryControl {
    /// Mechanical sample in the half-open admitted performance window.
    pub sample: u64,
    /// Original fixed node and physical pressure [Pa]/opening force [N].
    pub drive: BoundaryDrive,
}

/// Complete initial boundary conditions and finite changes on the source clock.
/// There is no inferred zero-pressure/zero-force program for an attached key.
#[derive(Clone, Debug)]
pub struct BoundaryProgram {
    /// Exactly one drive per attached boundary, in ascending graph-node order.
    pub initial: Vec<BoundaryDrive>,
    /// May be unsorted. Duplicate (sample, node) assignments refuse. A sample-zero
    /// assignment replaces the initial drive before the first physical solve.
    pub controls: Vec<BoundaryControl>,
}

/// One cold compiler and preallocated accepted/candidate state for both source
/// adapters. The cursor publishes only after physics AND observation succeed.
pub(in crate::bernoulli_aperture) struct BoundarySchedule {
    held: Vec<BoundaryDrive>,
    candidate: Vec<BoundaryDrive>,
    controls: Vec<BoundaryControl>,
    next: usize,
    candidate_next: usize,
    pub(in crate::bernoulli_aperture) work: u64,
}

impl BoundarySchedule {
    pub(in crate::bernoulli_aperture) fn compile(
        boundaries: &BoundaryApertures,
        mut program: BoundaryProgram,
        samples: u64,
        max_assignments: usize,
        max_work: u64,
    ) -> Result<Self, GestureCompileError> {
        let nodes = program.initial.len();
        let events = program.controls.len();
        if !(1..=64).contains(&nodes)
            || samples == 0
            || nodes
                .checked_add(events)
                .is_none_or(|n| n > max_assignments)
        {
            return Err(invalid(
                "boundary controls require 1..=64 initial drives and combined assignments within budget",
            ));
        }
        // Charge cold node validation, event lookup and sorting before any
        // additional storage or sorting work. A comparison-sort allowance grows
        // with the complete event count, not merely the number of changed nodes.
        let sort_levels = usize::BITS - events.saturating_sub(1).leading_zeros();
        let work = 2 * nodes as u128
            + events as u128 * (nodes as u128 + 4 * (u128::from(sort_levels) + 1));
        if work > u128::from(max_work) {
            return Err(invalid(
                "boundary controls exceed the shared auxiliary compilation work budget",
            ));
        }
        if !boundaries
            .nodes()
            .eq(program.initial.iter().map(|drive| drive.node))
            || program.initial.iter().any(|drive| !finite(drive))
        {
            return Err(invalid(
                "initial boundary drives must be finite and cover every original node in ascending order",
            ));
        }
        for control in &program.controls {
            if control.sample >= samples
                || !finite(&control.drive)
                || program
                    .initial
                    .binary_search_by_key(&control.drive.node, |drive| drive.node)
                    .is_err()
            {
                return Err(invalid(
                    "boundary controls require a fixed node, finite physical values and an in-window sample",
                ));
            }
        }
        program
            .controls
            .sort_unstable_by_key(|control| (control.sample, control.drive.node));
        if program.controls.windows(2).any(|pair| {
            (pair[0].sample, pair[0].drive.node) == (pair[1].sample, pair[1].drive.node)
        }) {
            return Err(invalid(
                "boundary controls must uniquely assign each (sample, node)",
            ));
        }
        let mut candidate = Vec::new();
        candidate
            .try_reserve_exact(nodes)
            .map_err(|_| invalid("cannot reserve complete boundary control candidates"))?;
        candidate.extend_from_slice(&program.initial);
        Ok(Self {
            held: program.initial,
            candidate,
            controls: program.controls,
            next: 0,
            candidate_next: 0,
            work: work as u64,
        })
    }

    pub(in crate::bernoulli_aperture) fn assignment_count(&self) -> usize {
        // The checked combined count was admitted before allocation.
        self.held.len() + self.controls.len()
    }
    pub(in crate::bernoulli_aperture) fn held(&self) -> &[BoundaryDrive] {
        &self.held
    }
    pub(in crate::bernoulli_aperture) fn applied(&self) -> &[BoundaryControl] {
        &self.controls[..self.next]
    }
    pub(in crate::bernoulli_aperture) fn pending(&self) -> &[BoundaryControl] {
        &self.controls[self.next..]
    }
    pub(in crate::bernoulli_aperture) fn prepare(&mut self, sample: u64) {
        self.candidate.copy_from_slice(&self.held);
        self.candidate_next = self.next;
        while let Some(control) = self.controls.get(self.candidate_next) {
            if control.sample != sample {
                break;
            }
            let slot = self
                .candidate
                .binary_search_by_key(&control.drive.node, |drive| drive.node)
                .expect("admitted controls preserve the complete fixed node set");
            self.candidate[slot] = control.drive;
            self.candidate_next += 1;
        }
    }
    pub(in crate::bernoulli_aperture) fn candidate(&self) -> &[BoundaryDrive] {
        &self.candidate
    }
    pub(in crate::bernoulli_aperture) fn accept(&mut self) {
        std::mem::swap(&mut self.held, &mut self.candidate);
        self.next = self.candidate_next;
    }
}

fn finite(drive: &BoundaryDrive) -> bool {
    drive.exterior_pressure_pa.is_finite() && drive.opening_force_n.is_finite()
}

impl AperturePerformance {
    /// Attach complete sample-addressed boundary pressures and opening forces
    /// before playback. Every physical key keeps its original material/contact
    /// owner. Events are held, with no additional interpolation or state reset.
    ///
    /// Initial drives and changes share max_controls with pressure and source
    /// force assignments. Boundary compilation and source-force projection share
    /// max_compile_work; canonical pressure compilation keeps its original
    /// independent work allowance. Attachment order cannot enlarge either bound.
    /// A keyed source refuses rendering until this explicit program is attached.
    ///
    /// # Errors
    /// Missing physical boundaries, repeated/late attachment, incomplete nodes,
    /// nonfinite values, duplicate/out-of-window events, or exhausted budgets.
    pub fn with_boundary_controls(
        mut self,
        program: BoundaryProgram,
    ) -> Result<Self, GestureCompileError> {
        if self.poisoned
            || self.system.aperture().accepted_steps() != 0
            || self.boundaries.is_some()
        {
            return Err(invalid(
                "boundary controls must bind once before performance playback",
            ));
        }
        let (assignments, work) = self.remaining_control_budgets()?;
        let CoupledAperture::Network(network) = &self.system else {
            return Err(invalid(
                "boundary controls require an explicitly keyed network",
            ));
        };
        let boundaries = network
            .boundary_apertures()
            .ok_or_else(|| invalid("boundary controls require attached physical apertures"))?;
        self.boundaries = Some(BoundarySchedule::compile(
            boundaries,
            program,
            self.config.samples,
            assignments,
            work,
        )?);
        Ok(self)
    }

    // Used by either cold attachment path, including force.rs. Pressure sampling
    // keeps its existing independent work budget; all auxiliary work is shared.
    pub(super) fn remaining_control_budgets(&self) -> Result<(usize, u64), GestureCompileError> {
        let mut assignments = self
            .config
            .max_controls
            .checked_sub(self.controls.len())
            .ok_or_else(|| invalid("pressure assignments exceed the combined storage budget"))?;
        let mut work = self.config.max_compile_work;
        if let Some(forces) = &self.forces {
            assignments = assignments
                .checked_sub(forces.events.len())
                .ok_or_else(|| {
                    invalid("pressure and source forces exceed the combined storage budget")
                })?;
            work = work.checked_sub(forces.work).ok_or_else(|| {
                invalid("source forces exceed the auxiliary compilation work budget")
            })?;
        }
        if let Some(boundaries) = &self.boundaries {
            assignments = assignments
                .checked_sub(boundaries.assignment_count())
                .ok_or_else(|| {
                    invalid(
                        "pressure, source forces and boundaries exceed the combined storage budget",
                    )
                })?;
            work = work.checked_sub(boundaries.work).ok_or_else(|| {
                invalid("source forces and boundaries exceed the auxiliary compilation work budget")
            })?;
        }
        if let Some(gas) = &self.finite_gas {
            assignments = assignments.checked_sub(gas.program.controls.len())
                .ok_or_else(|| invalid("heat and other source controls exceed the combined storage budget"))?;
            work = work.checked_sub(gas.program.controls.len() as u64)
                .ok_or_else(|| invalid("heat and other source controls exceed the auxiliary compilation work budget"))?;
        }
        if let Some(transport) = &self.finite_transport {
            assignments = assignments.checked_sub(transport.assignments)
                .ok_or_else(|| invalid("gas transport and other controls exceed combined storage"))?;
            work = work.checked_sub(transport.work)
                .ok_or_else(|| invalid("gas transport and other controls exceed auxiliary compilation work"))?;
        }
        if let Some(wall) = &self.finite_wall {
            assignments = assignments.checked_sub(wall.program.controls.len())
                .ok_or_else(|| invalid("wall heat and other source controls exceed the combined storage budget"))?;
            work = work.checked_sub(wall.program.controls.len() as u64)
                .ok_or_else(|| invalid("wall heat and other source controls exceed the auxiliary compilation work budget"))?;
        }
        Ok((assignments, work))
    }

    /// Accepted held physical drives, or explicit initial values before playback.
    #[must_use]
    pub fn held_boundary_drives(&self) -> &[BoundaryDrive] {
        self.boundaries.as_ref().map_or(&[], BoundarySchedule::held)
    }
    /// Changes published only with completely observed physical samples.
    #[must_use]
    pub fn applied_boundary_controls(&self) -> &[BoundaryControl] {
        self.boundaries
            .as_ref()
            .map_or(&[], BoundarySchedule::applied)
    }
    /// Remaining changes in admitted (sample, node) execution order.
    #[must_use]
    pub fn pending_boundary_controls(&self) -> &[BoundaryControl] {
        self.boundaries
            .as_ref()
            .map_or(&[], BoundarySchedule::pending)
    }
    /// Complete key/graph frame from the last successfully observed sample.
    #[must_use]
    pub const fn last_boundary_frame(&self) -> Option<&BoundaryNetworkFrame> {
        self.boundary_frame.as_ref()
    }
    /// Physical key-actuation work [J], separate from primary reed actuation.
    #[must_use]
    pub fn last_boundary_work_j(&self) -> f64 {
        self.boundary_frame
            .as_ref()
            .map_or(0.0, |frame| frame.mechanical_work_j)
    }
}
