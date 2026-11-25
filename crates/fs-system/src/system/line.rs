//! The load lines of a looped block.
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
//! - The compliance is how the block's smooth components, with their
//!   unknowns eliminated and every load-line owner's effort held, turn an
//!   effort into flow at that owner's ports. One Jacobian assembly excludes
//!   all load-line owners' columns; one back-substitution per effort gives
//!   each owner's self-compliance. Other owners' efforts remain independent
//!   unknowns of the full interface solve, so their coupling is not lost.
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

use core::ops::Range;

use fs_exec::Cx;

use super::{BlockWork, ROUNDING_UNITS, System, SystemError};
use crate::solve::{Condensed, factor_dense, solve_factored};
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

    pub(super) fn workspace_bytes(&self) -> usize {
        use core::mem::size_of;
        self.places.len() * size_of::<(usize, usize, usize)>()
            + self.ends.len() * size_of::<(usize, usize)>()
            + self.compliance.len() * size_of::<f64>()
            + self.factor.len() * size_of::<f64>()
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
        let largest = self
            .compliance
            .iter()
            .fold(0.0_f64, |seen, value| seen.max(value.abs()));
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

/// Woodbury tear on the outputs of a looped block's load-line component.
///
/// The block's interface Jacobian is `J = M - B D_L P_L`, where `M = I - T_N`
/// leaves the line component's columns out (held in `BlockWork::line_condensed`),
/// `B` is the driver slopes reading the line component's active outputs, and
/// `D_L` is the tangent of those outputs with respect to the line component's
/// unknown inputs. With `Z = M^-1 B` and `W = I - D_L P_L Z`, solving
/// `J x = b` reduces to `M u_0 = b`, `W y = D_L P_L u_0`, and `x = u_0 + Z y`:
/// a dense solve of size `w` (the active outputs) instead of all line inputs.
pub(super) struct LineTear {
    active_outputs: Vec<usize>,
    line_columns: Vec<usize>,
    node_readers: Range<usize>,
    m: usize,
    d_line: Vec<f64>,
    z_cols: Vec<f64>,
    w_dense: Vec<f64>,
    w_pivots: Vec<usize>,
    rhs_y: Vec<f64>,
    b_saved: Vec<f64>,
}

impl LineTear {
    pub(super) fn plan(
        m: usize,
        line_node: usize,
        unknowns: &[usize],
        input_owner: &[usize],
        active_outputs: Vec<usize>,
        node_readers: Range<usize>,
    ) -> Option<Self> {
        let w = active_outputs.len();
        if w == 0 {
            return None;
        }
        let line_columns: Vec<usize> = unknowns
            .iter()
            .enumerate()
            .filter_map(|(col, &index)| (input_owner[index] == line_node).then_some(col))
            .collect();
        let lc = line_columns.len();
        Some(Self {
            active_outputs,
            line_columns,
            node_readers,
            m,
            d_line: vec![0.0; w * lc],
            z_cols: vec![0.0; w * m],
            w_dense: vec![0.0; w * w],
            w_pivots: vec![0; w],
            rhs_y: vec![0.0; w],
            b_saved: vec![0.0; m],
        })
    }

    pub(super) fn width(&self) -> usize {
        self.active_outputs.len()
    }

    pub(super) fn active_outputs(&self) -> &[usize] {
        &self.active_outputs
    }

    pub(super) fn line_column_count(&self) -> usize {
        self.line_columns.len()
    }

    pub(super) fn workspace_bytes(&self) -> usize {
        use core::mem::size_of;
        self.active_outputs.len() * size_of::<usize>()
            + self.line_columns.len() * size_of::<usize>()
            + self.d_line.len() * size_of::<f64>()
            + self.z_cols.len() * size_of::<f64>()
            + self.w_dense.len() * size_of::<f64>()
            + self.w_pivots.len() * size_of::<usize>()
            + self.rhs_y.len() * size_of::<f64>()
            + self.b_saved.len() * size_of::<f64>()
    }

    pub(super) fn set_column(&mut self, col_pos: usize, outputs: &[f64]) {
        let lc = self.line_columns.len();
        for (row, &offset) in self.active_outputs.iter().enumerate() {
            self.d_line[row * lc + col_pos] = outputs[offset];
        }
    }

    pub(super) fn compute_z(
        &mut self,
        line_condensed: &mut Condensed,
        readers: &[(usize, usize, f64)],
    ) -> bool {
        for (col, &target_offset) in self.active_outputs.iter().enumerate() {
            let z_col = &mut self.z_cols[col * self.m..(col + 1) * self.m];
            z_col.fill(0.0);
            for &(row, offset, slope) in &readers[self.node_readers.clone()] {
                if offset == target_offset {
                    z_col[row] += slope;
                }
            }
            if !line_condensed.solve(z_col) {
                return false;
            }
        }
        true
    }

    pub(super) fn factor_w(&mut self) -> bool {
        let w = self.active_outputs.len();
        let lc = self.line_columns.len();
        self.w_dense.fill(0.0);
        for row in 0..w {
            self.w_dense[row * w + row] = 1.0;
            for col in 0..w {
                let z_col = &self.z_cols[col * self.m..(col + 1) * self.m];
                let mut dot = 0.0;
                for (col_pos, &unknown_col) in self.line_columns.iter().enumerate() {
                    dot += self.d_line[row * lc + col_pos] * z_col[unknown_col];
                }
                self.w_dense[row * w + col] -= dot;
            }
        }
        factor_dense(&mut self.w_dense, &mut self.w_pivots, w)
    }

    pub(super) fn solve(
        &mut self,
        line_condensed: &mut Condensed,
        readers: &[(usize, usize, f64)],
        b: &mut [f64],
    ) -> bool {
        self.b_saved.copy_from_slice(b);
        if !line_condensed.solve(b) {
            return false;
        }
        let w = self.active_outputs.len();
        let lc = self.line_columns.len();
        for row in 0..w {
            let mut sum = 0.0;
            for (col_pos, &unknown_col) in self.line_columns.iter().enumerate() {
                sum += self.d_line[row * lc + col_pos] * b[unknown_col];
            }
            self.rhs_y[row] = sum;
        }
        if !solve_factored(&self.w_dense, &self.w_pivots, &mut self.rhs_y, w) {
            return false;
        }
        // Form `b_orig + B y` before applying `M^-1` so cancelling effort terms
        // subtract before propagating into quiet follower rows.
        b.copy_from_slice(&self.b_saved);
        for &(row, offset, slope) in &readers[self.node_readers.clone()] {
            if let Some(pos) = self.active_outputs.iter().position(|&o| o == offset) {
                b[row] += slope * self.rhs_y[pos];
            }
        }
        line_condensed.solve(b)
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
        let mut lines = core::mem::take(&mut self.lines[block]);
        let found = self.find_lent(block, &mut lines, unknowns, nodes, work, dt_s, cx);
        self.lines[block] = lines;
        found
    }

    #[allow(clippy::too_many_arguments)]
    fn find_lent(
        &mut self,
        block: usize,
        lines: &mut [Line],
        unknowns: &[usize],
        nodes: &[usize],
        work: &mut BlockWork,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Option<f64>, SystemError> {
        if lines.len() == 1 {
            return self.find_lent_single(&mut lines[0], unknowns, nodes, work, dt_s, cx);
        }
        // All load-line owners are held while the smooth neighbours are
        // condensed. Asking another non-smooth owner for a tangent here
        // would make a stuck constraint part of a supposedly released line.
        for line in lines.iter_mut() {
            line.compliance.fill(0.0);
            line.write_drivers(&mut self.drivers);
            let held = self.node_outputs(line.node);
            self.outputs[held.clone()].copy_from_slice(&self.accepted_outputs[held]);
        }
        // What is evaluated here and still holds at the solve's first pass
        // is kept for it (`Found`): the inputs each component with a
        // declared tangent was proposed at, and the columns its tangent
        // gives below.
        work.found.clear();
        for (position, &node) in nodes.iter().enumerate() {
            if self.line_of[node].0 == block {
                continue;
            }
            if let Some(shorter) = self.propose(node, dt_s, cx)? {
                return Ok(Some(shorter));
            }
            let inputs = self.node_inputs(node);
            work.found.keep_inputs(position, &self.inputs[inputs]);
        }
        work.group_scale.fill(0.0);
        for &index in unknowns {
            let group = self.input_group[index];
            work.group_scale[group] = work.group_scale[group]
                .max(self.inputs[index].abs())
                .max(self.drivers[index].magnitude(&self.outputs, &self.boundary_values));
        }

        // I - T with every load-line owner's columns left out: their
        // outputs are independent efforts carried through the smooth
        // neighbors. The same factors supply every owner's self response.
        work.condensed.clear();
        for (row, terms) in &work.varying {
            self.drivers[unknowns[*row]].gradient(&self.outputs, &mut work.gradient);
            for (&reader, &(_, slope)) in work.varying_terms[terms.clone()]
                .iter()
                .zip(work.gradient.iter())
            {
                if reader != usize::MAX {
                    work.readers[reader].2 = slope;
                }
            }
        }
        for (column, &index) in unknowns.iter().enumerate() {
            let node = self.input_owner[index];
            if self.line_of[node].0 == block {
                continue;
            }
            let local = index - self.input_start[node];
            if let Some(shorter) = self.tangent_column(node, local, dt_s, &work.group_scale, false, cx)? {
                return Ok(Some(shorter));
            }
            // A declared tangent's column depends on the inputs, the
            // accepted state and the step only.
            let values = &self.scratch_outputs[..self.node_outputs(node).len()];
            work.found.keep_column(column, values);
            let reading = &work.readers[work.column_readers[column].clone()];
            for (&(_, offset, slope), &slot) in
                reading.iter().zip(&work.slots[work.slot_start[column]..])
            {
                work.condensed
                    .subtract(slot, slope * self.scratch_outputs[offset]);
            }
        }
        // Singular or not finite: the other components do not say what
        // flow answers an effort, so there is no line. The refusal names
        // the component that needed one.
        let component = lines[0].node;
        let no_line = |system: &Self| SystemError::LoadLine {
            node: system.names[component].clone(),
            reason: "the other components of its block do not determine how the flows at its \
                     ports answer an effort: their Jacobian is singular or not finite",
        };
        if !work.condensed.factor() {
            return Err(no_line(self));
        }
        for line in lines.iter_mut() {
            let width = line.places.len();
            // Every row that reads this owner's output, whichever of its
            // inputs the representative column belongs to.
            let own_column = unknowns
                .iter()
                .rposition(|&input| self.input_owner[input] == line.node)
                .expect("a looped load-line owner has an unknown");
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
        }

        // Each owner's inputs release its own efforts only. The other
        // efforts are still carried by the original connection equations
        // and are solved together in the full Newton iteration.
        for line in lines {
            line.write_drivers(&mut self.drivers);
            self.released_flows(line);
        }
        Ok(None)
    }


    /// Single-load-line fast path using `work.line_condensed`, optional
    /// Woodbury `LineTear`, and `affine_line` factor/compliance caching at an
    /// unchanged `dt_s`. Multi-line blocks use `find_lent` above.
    #[allow(clippy::too_many_lines)]
    fn find_lent_single(
        &mut self,
        line: &mut Line,
        unknowns: &[usize],
        nodes: &[usize],
        work: &mut BlockWork,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Option<f64>, SystemError> {
        let width = line.places.len();
        let held = self.node_outputs(line.node);
        self.outputs[held.clone()].copy_from_slice(&self.accepted_outputs[held]);
        if work.affine_line && work.kept_line_dt_s.to_bits() == dt_s.to_bits() {
            work.found.clear_proposed();
            for (position, &node) in nodes.iter().enumerate() {
                if node == line.node {
                    continue;
                }
                if let Some(shorter) = self.propose(node, dt_s, cx)? {
                    return Ok(Some(shorter));
                }
                let inputs = self.node_inputs(node);
                work.found.keep_inputs(position, &self.inputs[inputs]);
            }
            line.write_drivers(&mut self.drivers);
            self.released_flows(line);
            return Ok(None);
        }
        work.kept_line_dt_s = f64::NAN;
        line.compliance.fill(0.0);
        line.write_drivers(&mut self.drivers);
        work.found.clear();
        for (position, &node) in nodes.iter().enumerate() {
            if node == line.node {
                continue;
            }
            if let Some(shorter) = self.propose(node, dt_s, cx)? {
                return Ok(Some(shorter));
            }
            let inputs = self.node_inputs(node);
            work.found.keep_inputs(position, &self.inputs[inputs]);
        }
        work.group_scale.fill(0.0);
        for &index in unknowns {
            let group = self.input_group[index];
            work.group_scale[group] = work.group_scale[group]
                .max(self.inputs[index].abs())
                .max(self.drivers[index].magnitude(&self.outputs, &self.boundary_values));
        }
        work.line_condensed.clear();
        for (row, terms) in &work.varying {
            self.drivers[unknowns[*row]].gradient(&self.outputs, &mut work.gradient);
            for (&reader, &(_, slope)) in work.varying_terms[terms.clone()]
                .iter()
                .zip(work.gradient.iter())
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
                self.tangent_column(node, local, dt_s, &work.group_scale, false, cx)?
            {
                return Ok(Some(shorter));
            }
            let values = &self.scratch_outputs[..self.node_outputs(node).len()];
            work.found.keep_column(column, values);
            let reading = &work.readers[work.column_readers[column].clone()];
            for (&(_, offset, slope), &slot) in
                reading.iter().zip(&work.line_slots[work.line_slot_start[column]..])
            {
                work.line_condensed
                    .subtract(slot, slope * self.scratch_outputs[offset]);
            }
        }
        let component = line.node;
        let no_line = |system: &Self| SystemError::LoadLine {
            node: system.names[component].clone(),
            reason: "the other components of its block do not determine how the flows at its \
                     ports answer an effort: their Jacobian is singular or not finite",
        };
        if !work.line_condensed.factor() {
            return Err(no_line(self));
        }
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
            if !work.line_condensed.solve(&mut work.residual) {
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
        line.write_drivers(&mut self.drivers);
        if work.affine_line {
            for (row, terms) in &work.varying {
                self.drivers[unknowns[*row]].gradient(&self.outputs, &mut work.gradient);
                for (&reader, &(_, slope)) in work.varying_terms[terms.clone()]
                    .iter()
                    .zip(work.gradient.iter())
                {
                    if reader != usize::MAX {
                        work.readers[reader].2 = slope;
                    }
                }
            }
            if let Some(tear) = &mut work.line_tear
                && !tear.compute_z(&mut work.line_condensed, &work.readers)
            {
                return Err(no_line(self));
            }
            work.kept_line_dt_s = dt_s;
        }
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

    /// Every input of the component on `block`'s line (the line's flows, the
    /// ends of its flows, and any other input, such as the span of a side
    /// it is told) takes one update of `block`'s unknowns with that
    /// component's efforts held at what the last step accepted, by the
    /// factors [`Self::find_line`] left in `condensed`: the block's Jacobian
    /// with the component's own columns left out. Where the other
    /// components are affine, its inputs then hold what they give this step
    /// at those efforts; where they are not, their first Newton iterate.
    /// Every other unknown keeps its starting value.
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
        // The efforts the last step accepted, whatever an earlier trial of
        // this step left in the outputs.
        for line in &self.lines[block] {
            let held = self.node_outputs(line.node);
            self.outputs[held.clone()].copy_from_slice(&self.accepted_outputs[held]);
        }
        for (row, &index) in unknowns.iter().enumerate() {
            residual[row] = self.drivers[index].value(&self.outputs, &self.boundary_values)
                - self.inputs[index];
        }
        if !condensed.solve(residual) {
            return Err(SystemError::LoadLine {
                node: self.names[self.lines[block][0].node].clone(),
                reason: "the other components of its block do not determine the flows its held \
                         efforts release: their Jacobian is singular or not finite",
            });
        }
        for (row, &index) in unknowns.iter().enumerate() {
            if self.line_of[self.input_owner[index]].0 == block {
                self.inputs[index] += residual[row];
            }
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
        for line in &self.lines[block] {
            let width = line.places.len();
            // A flow of this line, or its end: the same difference.
            let (flow, input) =
                if let Some(flow) = line.places.iter().position(|place| place.2 == row) {
                    (flow, line.places[flow].0)
                } else if let Some(end) = line.ends.iter().position(|end| end.1 == row) {
                    (width + end, line.ends[end].0)
                } else {
                    continue;
                };
            let mut crossing = self.inputs[input];
            let mut terms = crossing.abs();
            for (effort, &(_, output, _)) in line.places.iter().enumerate() {
                let taken = line.compliance[flow * width + effort] * self.outputs[output];
                crossing -= taken;
                terms += taken.abs();
            }
            return Some((crossing, ROUNDING_UNITS * f64::EPSILON * terms));
        }
        None
    }

    /// Whether input `index` is a flow of `block`'s load line, or the end of
    /// one, that no effort of the line component moves: its row of the
    /// line's matrix is zero. The block then gives it whatever the component
    /// does (a velocity source's flow), and the component reads it as a
    /// given end.
    ///
    /// Its column is not differenced. A perturbation of it is a jump of a
    /// given end that nothing in the block makes, and a component that
    /// locates the stop of the slip a jump starts asks, on both sides of
    /// the iterate, for the step that ends there: a stuck friction interface
    /// beside a velocity source had every step cut to the `2^-27` s the
    /// kernel's own perturbation of `2^-26` m/s took to stop (bead
    /// `.5.8.1`).
    pub(super) fn given_on_line(&self, block: usize, index: usize) -> bool {
        let (owner_block, position) = self.line_of[self.input_owner[index]];
        if owner_block != block {
            return false;
        }
        let line = &self.lines[block][position];
        let width = line.places.len();
        let row = if let Some(flow) = line.places.iter().position(|place| place.0 == index) {
            flow
        } else if let Some(end) = line.ends.iter().position(|end| end.0 == index) {
            width + end
        } else {
            return false;
        };
        line.compliance[row * width..(row + 1) * width]
            .iter()
            .all(|&value| value == 0.0)
    }

    /// Replace the released flows of `block`'s line by the flows on the
    /// line at the efforts just accepted: the ports, the audit and the
    /// ledger see what crossed. The rounding of these flows is already in
    /// what the solve told the audit ([`Self::crossing`]).
    pub(super) fn line_flows(&mut self, block: usize) {
        for line in &self.lines[block] {
            let width = line.places.len();
            for (flow, &(input, _, unknown)) in line.places.iter().enumerate() {
                if unknown == usize::MAX {
                    continue;
                }
                for (effort, &(_, output, _)) in line.places.iter().enumerate() {
                    self.inputs[input] -=
                        line.compliance[flow * width + effort] * self.outputs[output];
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
}
