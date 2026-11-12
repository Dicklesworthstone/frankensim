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

use fs_exec::Cx;

use super::{BlockWork, ROUNDING_UNITS, System, SystemError};
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
    /// The line's matrix of this step, row-major over the coefficients.
    pub(super) compliance: Vec<f64>,
    /// Room for the factor of the symmetric-part test.
    factor: Vec<f64>,
}

impl Line {
    pub(super) fn new(node: usize, places: Vec<(usize, usize, usize)>) -> Self {
        let width = places.len();
        Self {
            node,
            places,
            compliance: vec![0.0; width * width],
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
            if let Some(shorter) = self.tangent_column(node, local, dt_s, &work.group_scale, cx)? {
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
        for (flow, &(input, _, unknown)) in line.places.iter().enumerate() {
            if unknown == usize::MAX {
                continue;
            }
            for (effort, &(_, output, _)) in line.places.iter().enumerate() {
                self.inputs[input] +=
                    line.compliance[flow * width + effort] * self.accepted_outputs[output];
            }
        }
        Ok(None)
    }

    /// Replace the released flows of `block`'s line by the flows on the
    /// line at the efforts just accepted.
    ///
    /// A flow on the line is a difference of the released flow and what
    /// the efforts take away, and it carries the rounding of those terms,
    /// not of its own size: a held slider's flow is `r - c (r / c)`, zero
    /// or a unit in the last place of `r`. The interface audit is told
    /// that rounding, in the units the solve allows everywhere else; told
    /// nothing, it refused a mass held at rest for one push in seven.
    pub(super) fn line_flows(&mut self, block: usize) {
        let Some(line) = &self.lines[block] else {
            return;
        };
        let width = line.places.len();
        for (flow, &(input, _, unknown)) in line.places.iter().enumerate() {
            if unknown == usize::MAX {
                continue;
            }
            let mut terms = self.inputs[input].abs();
            for (effort, &(_, output, _)) in line.places.iter().enumerate() {
                let taken = line.compliance[flow * width + effort] * self.outputs[output];
                self.inputs[input] -= taken;
                terms += taken.abs();
            }
            self.input_allowance[input] += ROUNDING_UNITS * f64::EPSILON * terms;
        }
    }
}
