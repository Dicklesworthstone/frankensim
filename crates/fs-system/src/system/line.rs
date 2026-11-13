//! The load line of a looped block.
//!
//! A component whose effort is not a function of the flow it receives (dry
//! friction, a closed contact) cannot be iterated on that flow. Its
//! neighbours, frozen at an iterate, answer an effort with a flow that is
//! affine in it,
//!
//! ```text
//! flow = released - compliance * effort
//! ```
//!
//! and the component is stepped on that line
//! ([`fs_ports::component::LoadLine`]). This module finds the line once per
//! step and keeps the block's unknowns consistent with it:
//!
//! - The compliance is how the block's other components, with every other
//!   unknown eliminated, turn an effort of the line component into flow at
//!   its own ports. It comes from one assembly of the block's Jacobian with
//!   the line component's own columns left out, and one back-substitution
//!   per effort.
//! - During the solve the line's flow inputs hold the released flows. Their
//!   drivers gain the term `compliance * own effort`, so the residual of a
//!   line input is still the residual of the actual flow, and Newton runs
//!   on a map that is a function: released flow in, effort out.
//! - When the block is accepted, the actual flows are written back, so the
//!   ports, the audit and the ledger see what crossed.
//! - The solve starts from the flows the last step ended with, which are
//!   not the released flows of this one. A component that locates an event
//!   (it asks for the step that ends at a stop) would locate it from the
//!   wrong flows there, and ask again in every shorter trial. The first
//!   time it asks before any update of a trial, the other components take
//!   one update with its efforts held, by the factors of the line's search
//!   ([`System::release`]), and it is asked again.

use fs_exec::Cx;

use super::{BlockWork, ROUNDING_UNITS, System, SystemError};
use crate::solve::Condensed;
use crate::wiring::Driver;

/// A symmetric part with an eigenvalue below this fraction of the matrix's
/// largest entry, negated, is refused. A neighbour without a declared
/// tangent is differenced, so the compliance is known to about the square
/// root of the rounding unit and no better.
const NEGATIVE_FRACTION: f64 = 1.0e-6;

/// What a looped block keeps for its component on a load line.
pub(super) struct Line {
    /// The component.
    pub(super) node: usize,
    /// Per coefficient of the line: the global index of its flow input,
    /// the global index of its effort output, and the flow's row among the
    /// block's unknowns (`usize::MAX` for a flow the block does not solve
    /// for, which has a zero row and column in the matrix).
    pub(super) places: Vec<(usize, usize, usize)>,
    /// Per coefficient of the line's end inputs, if the component names
    /// any: the global index of the input and its row among the block's
    /// unknowns (`usize::MAX` for one the block does not solve for).
    pub(super) ends: Vec<(usize, usize)>,
    /// The line's matrix of this step, row-major over the coefficients,
    /// followed by one row per end input: how the flow at the END of the
    /// step answers each effort.
    pub(super) compliance: Vec<f64>,
    /// Room for the factor of the symmetric-part test.
    factor: Vec<f64>,
}

impl Line {
    pub(super) fn new(
        node: usize,
        places: Vec<(usize, usize, usize)>,
        ends: Vec<(usize, usize)>,
    ) -> Self {
        let width = places.len();
        Self {
            node,
            compliance: vec![0.0; (width + ends.len()) * width],
            places,
            ends,
            factor: vec![0.0; width * width],
        }
    }

    /// Write the matrix into the drivers of the line's flows: the last
    /// terms of each are the component's own efforts.
    fn write_drivers(&self, drivers: &mut [Driver]) {
        let width = self.places.len();
        for (row, &(input, _, unknown)) in self.places.iter().enumerate() {
            if unknown == usize::MAX {
                continue;
            }
            if let Driver::OnLine { terms, own } = &mut drivers[input] {
                let first = terms.len() - *own;
                for (column, term) in terms[first..].iter_mut().enumerate() {
                    term.1 = self.compliance[row * width + column];
                }
            }
        }
        for (end, &(input, unknown)) in self.ends.iter().enumerate() {
            if unknown == usize::MAX {
                continue;
            }
            if let Driver::OnLine { terms, own } = &mut drivers[input] {
                let first = terms.len() - *own;
                for (column, term) in terms[first..].iter_mut().enumerate() {
                    term.1 = self.compliance[(width + end) * width + column];
                }
            }
        }
    }

    /// Column `effort` of the rows of the end inputs, from the
    /// back-substitution that gave that column of the line: how the END of
    /// each flow answers the effort. `answer` is over the block's unknowns.
    fn end_answers(&mut self, effort: usize, answer: &[f64]) {
        let width = self.places.len();
        for (end, &(_, unknown)) in self.ends.iter().enumerate() {
            self.compliance[(width + end) * width + effort] = if unknown == usize::MAX {
                0.0
            } else {
                -answer[unknown]
            };
        }
    }

    /// Whether the matrix is finite and its symmetric part is not negative
    /// beyond what the matrix is known to: `(C + C^T) / 2` shifted by
    /// [`NEGATIVE_FRACTION`] of the largest entry has a Cholesky factor.
    fn passive(&mut self) -> bool {
        let width = self.places.len();
        // Checked entry by entry: a maximum would pass over a NaN.
        if !self.compliance.iter().all(|value| value.is_finite()) {
            return false;
        }
        let largest = self.compliance.iter().fold(0.0_f64, |seen, value| seen.max(value.abs()));
        if largest == 0.0 {
            return true;
        }
        let shift = NEGATIVE_FRACTION * largest;
        for row in 0..width {
            for column in 0..=row {
                let mut sum = 0.5
                    * (self.compliance[row * width + column]
                        + self.compliance[column * width + row]);
                if row == column {
                    sum += shift;
                }
                for before in 0..column {
                    sum -= self.factor[row * width + before] * self.factor[column * width + before];
                }
                if row == column {
                    if !(sum > 0.0) {
                        return false;
                    }
                    self.factor[row * width + row] = sum.sqrt();
                } else {
                    self.factor[row * width + column] = sum / self.factor[column * width + column];
                }
            }
        }
        true
    }
}

impl System {
    /// Find the load line of `block` for this step and put the block's
    /// inputs on it. Called once the unknowns hold their starting values
    /// and before the first check. `Ok(Some(dt))` asks to shorten.
    pub(super) fn find_line(
        &mut self,
        block: usize,
        unknowns: &[usize],
        nodes: &[usize],
        work: &mut BlockWork,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Option<f64>, SystemError> {
        let mut line = self.lines[block]
            .take()
            .expect("a block asked for its line has one");
        let found = self.find_lent(&mut line, unknowns, nodes, work, dt_s, cx);
        self.lines[block] = Some(line);
        found
    }

    fn find_lent(
        &mut self,
        line: &mut Line,
        unknowns: &[usize],
        nodes: &[usize],
        work: &mut BlockWork,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Option<f64>, SystemError> {
        let width = line.places.len();
        // The neighbours are evaluated with the line component's outputs
        // held at what the last step accepted, and with no compliance yet
        // in the drivers: the flow inputs hold actual flows.
        line.compliance.fill(0.0);
        line.write_drivers(&mut self.drivers);
        let held = self.node_outputs(line.node);
        self.outputs[held.clone()].copy_from_slice(&self.accepted_outputs[held]);
        for &node in nodes {
            if node == line.node {
                continue;
            }
            if let Some(shorter) = self.propose(node, dt_s, cx)? {
                return Ok(Some(shorter));
            }
        }
        work.group_scale.fill(0.0);
        for &index in unknowns {
            let group = self.input_group[index];
            work.group_scale[group] = work.group_scale[group]
                .max(self.inputs[index].abs())
                .max(self.drivers[index].magnitude(&self.outputs, &self.boundary_values));
        }

        // I - T with the line component's own columns left out: its
        // outputs are then free, and a unit of one of them is carried
        // through everything else.
        work.condensed.clear();
        for (row, terms) in &work.varying {
            self.drivers[unknowns[*row]].gradient(&self.outputs, &mut work.gradient);
            for (&reader, &(_, slope)) in
                work.varying_terms[terms.clone()].iter().zip(work.gradient.iter())
            {
                if reader != usize::MAX {
                    work.readers[reader].2 = slope;
                }
            }
        }
        let mut own_column = usize::MAX;
        for (column, &index) in unknowns.iter().enumerate() {
            let node = self.input_owner[index];
            if node == line.node {
                own_column = column;
                continue;
            }
            let local = index - self.input_start[node];
            if let Some(shorter) =
                self.tangent_column(node, local, dt_s, &work.group_scale, true, cx)?
            {
                return Ok(Some(shorter));
            }
            let reading = &work.readers[work.column_readers[column].clone()];
            for (&(_, offset, slope), &slot) in
                reading.iter().zip(&work.slots[work.slot_start[column]..])
            {
                work.condensed.subtract(slot, slope * self.scratch_outputs[offset]);
            }
        }
        // Singular or not finite: the other components do not say what
        // flow answers an effort, so there is no line. The refusal names
        // the component that needed one.
        let component = line.node;
        let no_line = |system: &Self| SystemError::LoadLine {
            node: system.names[component].clone(),
            reason: "the other components of its block do not determine how the flows at its \
                     ports answer an effort: their Jacobian is singular or not finite",
        };
        if !work.condensed.factor() {
            return Err(no_line(self));
        }
        // Every row that reads an output of the line component, whichever
        // of its inputs the column is.
        let reading = work.column_readers[own_column].clone();
        let first_output = self.output_start[line.node];
        for effort in 0..width {
            let local = line.places[effort].1 - first_output;
            work.residual.fill(0.0);
            for &(row, offset, slope) in &work.readers[reading.clone()] {
                if offset == local {
                    work.residual[row] += slope;
                }
            }
            if !work.condensed.solve(&mut work.residual) {
                return Err(no_line(self));
            }
            for (flow, &(_, _, unknown)) in line.places.iter().enumerate() {
                line.compliance[flow * width + effort] = if unknown == usize::MAX {
                    0.0
                } else {
                    -work.residual[unknown]
                };
            }
            line.end_answers(effort, &work.residual);
        }
        if !line.passive() {
            return Err(SystemError::LoadLine {
                node: self.names[line.node].clone(),
                reason: "its neighbours do not answer an effort with less flow: the line's \
                         matrix is not finite or its symmetric part is negative",
            });
        }

        // On the line: the drivers carry the matrix, and the flow inputs
        // become the released flows of the efforts accepted last.
        line.write_drivers(&mut self.drivers);
        self.released_flows(line);
        Ok(None)
    }

    /// The flow inputs of `line` become the released flows of the efforts
    /// accepted last, and its end inputs the released ends: what the
    /// component is handed at its first proposal of the step goes together.
    /// The counterpart of [`Self::line_flows`].
    fn released_flows(&mut self, line: &Line) {
        let width = line.places.len();
        for (flow, &(input, _, unknown)) in line.places.iter().enumerate() {
            if unknown == usize::MAX {
                continue;
            }
            for (effort, &(_, output, _)) in line.places.iter().enumerate() {
                self.inputs[input] +=
                    line.compliance[flow * width + effort] * self.accepted_outputs[output];
            }
        }
        for (end, &(input, unknown)) in line.ends.iter().enumerate() {
            if unknown == usize::MAX {
                continue;
            }
            for (effort, &(_, output, _)) in line.places.iter().enumerate() {
                self.inputs[input] +=
                    line.compliance[(width + end) * width + effort] * self.accepted_outputs[output];
            }
        }
    }

    /// One update of `block`'s unknowns with the efforts of its line's
    /// component held at what the last step accepted, by the factors
    /// [`Self::find_line`] left in `condensed`: the block's Jacobian with
    /// the component's own columns left out. Where the other components
    /// are affine, the line's inputs then hold the released flows of this
    /// step, and every other unknown what goes with the held efforts.
    ///
    /// For the step's starting values only, before the block's own first
    /// assembly replaces those factors.
    pub(super) fn release(
        &mut self,
        block: usize,
        unknowns: &[usize],
        residual: &mut [f64],
        condensed: &mut Condensed,
    ) -> Result<(), SystemError> {
        let node = self.lines[block]
            .as_ref()
            .expect("a block asked to release its line has one")
            .node;
        // The component asked to shorten, which leaves its outputs
        // unspecified.
        let held = self.node_outputs(node);
        self.outputs[held.clone()].copy_from_slice(&self.accepted_outputs[held]);
        for (row, &index) in unknowns.iter().enumerate() {
            residual[row] =
                self.drivers[index].value(&self.outputs, &self.boundary_values) - self.inputs[index];
        }
        if !condensed.solve(residual) {
            return Err(SystemError::LoadLine {
                node: self.names[node].clone(),
                reason: "the other components of its block do not determine the flows its held \
                         efforts release: their Jacobian is singular or not finite",
            });
        }
        for (row, &index) in unknowns.iter().enumerate() {
            self.inputs[index] += residual[row];
        }
        Ok(())
    }

    /// For unknown `row` of `block`, if it is a flow on the block's load
    /// line: the flow that crosses at the current iterate, and the rounding
    /// that flow carries.
    ///
    /// During the solve the input holds the released flow; what crosses is
    /// that less what the efforts take away. It is a difference, and it
    /// carries the rounding of its terms, not of its own size: a held
    /// slider's flow is `r - c (r / c)`, zero or a unit in the last place
    /// of `r`. The rounding is given in the units the solve allows
    /// everywhere else, `ROUNDING_UNITS` of `eps` times the sum of the
    /// terms' magnitudes.
    ///
    /// The solve judges the row by the crossing flow, as it judges every
    /// other flow by its own size, and allows it that rounding; the
    /// interface audit is told the same rounding. (Judged by the released
    /// flow, a mass held under a push that changed from step to step was
    /// accepted with a flow inside the solve's relative tolerance of the
    /// released flow, far above the rounding the audit is told at a
    /// junction at rest, and the audit refused it. Told no rounding, the
    /// audit refused a mass held under a constant push for one push in
    /// seven.)
    pub(super) fn crossing(&self, block: usize, row: usize) -> Option<(f64, f64)> {
        let line = self.lines[block].as_ref()?;
        let width = line.places.len();
        // A flow of the line, or the end of one: the same difference.
        let (flow, input) = if let Some(flow) = line.places.iter().position(|place| place.2 == row)
        {
            (flow, line.places[flow].0)
        } else {
            let end = line.ends.iter().position(|end| end.1 == row)?;
            (width + end, line.ends[end].0)
        };
        let mut crossing = self.inputs[input];
        let mut terms = crossing.abs();
        for (effort, &(_, output, _)) in line.places.iter().enumerate() {
            let taken = line.compliance[flow * width + effort] * self.outputs[output];
            crossing -= taken;
            terms += taken.abs();
        }
        Some((crossing, ROUNDING_UNITS * f64::EPSILON * terms))
    }

    /// Replace the released flows of `block`'s line by the flows on the
    /// line at the efforts just accepted: the ports, the audit and the
    /// ledger see what crossed. The rounding of these flows is already in
    /// what the solve told the audit ([`Self::crossing`]).
    pub(super) fn line_flows(&mut self, block: usize) {
        let Some(line) = &self.lines[block] else {
            return;
        };
        let width = line.places.len();
        for (flow, &(input, _, unknown)) in line.places.iter().enumerate() {
            if unknown == usize::MAX {
                continue;
            }
            for (effort, &(_, output, _)) in line.places.iter().enumerate() {
                self.inputs[input] -= line.compliance[flow * width + effort] * self.outputs[output];
            }
        }
        for (end, &(input, unknown)) in line.ends.iter().enumerate() {
            if unknown == usize::MAX {
                continue;
            }
            for (effort, &(_, output, _)) in line.places.iter().enumerate() {
                self.inputs[input] -=
                    line.compliance[(width + end) * width + effort] * self.outputs[output];
            }
        }
    }
}
