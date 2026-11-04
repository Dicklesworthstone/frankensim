//! Shared original acoustic topology and polynomial entries for AC and poles.
use super::*;
use fs_duct::dg::continuous::DgContinuous;
use fs_vfit::mimo::MimoStateSpace;

/// Only the algebraic consumer differs: AC evaluates a+s*b immediately;
/// resonance retains the same real coefficients and their analytic derivative.
pub(super) trait Entries {
    fn dimension(&self) -> usize;
    fn set(&mut self, row: usize, column: usize, constant: f64, linear: f64);
    fn clear_row(&mut self, row: usize);
}

pub(super) struct AcEntries<'a> {
    pub values: &'a mut [C64],
    pub dimension: usize,
    pub s: C64,
}
impl Entries for AcEntries<'_> {
    fn dimension(&self) -> usize {
        self.dimension
    }
    fn set(&mut self, row: usize, column: usize, constant: f64, linear: f64) {
        self.values[row * self.dimension + column] = C64::from_re(constant) + self.s.scale(linear);
    }
    fn clear_row(&mut self, row: usize) {
        self.values[row * self.dimension..(row + 1) * self.dimension].fill(C64::ZERO);
    }
}

// D'= -i*linear, so -D'*x = i*linear*x. Each nonzero dynamic
// coefficient is assigned exactly once by append_source; rigid closures clear
// their original algebraic rows. No derivative matrix or terminal copy is kept.
struct FrequencyAction<'a> {
    rhs: &'a mut [C64],
    state: &'a [C64],
}
impl Entries for FrequencyAction<'_> {
    fn dimension(&self) -> usize {
        self.state.len()
    }
    fn set(&mut self, row: usize, column: usize, _: f64, linear: f64) {
        if linear != 0.0 {
            self.rhs[row] = self.rhs[row] + C64::new(0.0, linear) * self.state[column];
        }
    }
    fn clear_row(&mut self, row: usize) {
        self.rhs[row] = C64::ZERO;
    }
}

/// Source columns borrow the caller's views; no operator or physical basis is
/// copied, fitted, projected or inferred from a time-step transition.
pub(super) struct DgBinding<'a> {
    pub view: DgContinuous<'a>,
    pub branch_views: &'a [DgContinuous<'a>],
    pub pressure_column: &'a [f64],
    pub branch_offsets: Vec<usize>,
    pub branch_ports: Vec<(PortId, Option<PortId>)>,
    pub junction_columns: Vec<(&'a [f64], &'a [f64])>,
    pub terminal_columns: Vec<(usize, &'a [f64])>,
    pub terminal_ports: &'a [PortId],
    pub terminal_owner: Option<&'a PreparedMimoLoad>,
    pub channels: &'a [DgTerminalChannelLoad],
    pub acoustic_dimension: usize,
    pub terminal_dimension: usize,
    pub channel_dimension: usize,
}

pub(super) struct AppendedDg {
    pub terminal_source: Option<MimoStateSpace>,
    pub terminal_offset: usize,
    pub flow_offset: usize,
    pub junction_offset: usize,
    pub channel_offsets: Vec<usize>,
}

impl<'a> DgBinding<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        n: usize,
        view: &'a DgContinuous<'a>,
        inlet: PortId,
        flow_column: &[f64],
        terminal: Option<(&'a PreparedMimoLoad, &'a [PortId])>,
        branches: &[DgBranchLoad],
        branch_views: &'a [DgContinuous<'a>],
        channels: &'a [DgTerminalChannelLoad],
        gate: &CancelGate,
    ) -> Result<Self, DgHarmonicCouplingError> {
        poll(Some(gate))?;
        let terminal_ports = terminal.map_or(&[][..], |(_, ports)| ports);
        let sealed = branches
            .iter()
            .filter(|branch| branch.duct().outlet() == DgBoundary::Closed)
            .count();
        if flow_column.len() != n
            || flow_column.iter().any(|b| !b.is_finite())
            || !flow_column.iter().any(|b| *b != 0.0)
            || view.ports().len() != terminal_ports.len() + 1 + sealed
            || branches.len() != branch_views.len()
            || terminal.is_some_and(|(owner, ports)| owner.ports() != ports.len())
            || terminal_ports
                .iter()
                .enumerate()
                .any(|(i, id)| *id == inlet || terminal_ports[..i].contains(id))
        {
            return Err(invalid("continuous modal/DG coupling requires a complete original volume projection and one distinct terminal channel for every non-inlet port").into());
        }
        let mut channel_dimension = 0usize;
        for (index, channel) in channels.iter().enumerate() {
            poll(Some(gate))?;
            if !terminal_ports.contains(&channel.input().port)
                || channels[..index]
                    .iter()
                    .any(|other| other.input().port == channel.input().port)
            {
                return Err(invalid(
                    "continuous gaps require distinct original exterior terminal ports",
                )
                .into());
            }
            if let Some(original) = channel.channel() {
                let field = original.continuous_operator()?;
                channel_dimension = channel_dimension
                    .checked_add(field.dimension())
                    .and_then(|value| value.checked_add(1))
                    .ok_or_else(|| invalid("continuous terminal field dimension overflow"))?;
            }
        }
        let pressure_column = &view
            .ports()
            .iter()
            .find(|(id, _)| *id == inlet)
            .ok_or_else(|| invalid("continuous DG inlet is absent from its source basis"))?
            .1;
        let mut acoustic_dimension = view.dimension();
        let mut branch_offsets = Vec::with_capacity(branches.len());
        let mut branch_ports = Vec::with_capacity(branches.len());
        let mut junction_columns = Vec::with_capacity(branches.len());
        for (index, (branch, original)) in branches.iter().zip(branch_views).enumerate() {
            poll(Some(gate))?;
            let DgBoundary::Port(branch_inlet) = branch.duct().inlet() else {
                return Err(
                    invalid("continuous finite branch lost its original inlet flow port").into(),
                );
            };
            let branch_outlet = match branch.duct().outlet() {
                DgBoundary::Port(outlet) => Some(outlet),
                DgBoundary::Closed => None,
                _ => {
                    return Err(invalid(
                        "continuous finite branch lost its admitted distal boundary",
                    )
                    .into());
                }
            };
            if original.ports().len() != 1 + usize::from(branch_outlet.is_some())
                || branch_outlet == Some(branch_inlet)
                || terminal_ports.contains(&branch.port()) != branch_outlet.is_some()
                || branches[..index]
                    .iter()
                    .any(|other| other.port() == branch.port())
            {
                return Err(invalid("continuous finite branches require distinct owned junctions and complete original terminal channels").into());
            }
            let bore_column = &view
                .ports()
                .iter()
                .find(|(id, _)| *id == branch.port())
                .ok_or_else(|| invalid("continuous finite branch lost its main-bore junction"))?
                .1;
            let inlet_column = &original
                .ports()
                .iter()
                .find(|(id, _)| *id == branch_inlet)
                .ok_or_else(|| invalid("continuous finite branch lost its inlet pressure trace"))?
                .1;
            branch_offsets.push(acoustic_dimension);
            branch_ports.push((branch_inlet, branch_outlet));
            junction_columns.push((bore_column.as_slice(), inlet_column.as_slice()));
            acoustic_dimension = acoustic_dimension
                .checked_add(original.dimension())
                .ok_or_else(|| invalid("continuous finite branch dimension overflow"))?;
        }
        let terminal_columns = terminal_ports
            .iter()
            .map(|id| {
                let (original, offset, local_id) =
                    match branches.iter().position(|branch| branch.port() == *id) {
                        Some(index) => (
                            &branch_views[index],
                            branch_offsets[index],
                            branch_ports[index].1.ok_or_else(|| {
                                invalid("sealed branch cannot occupy an exterior harmonic channel")
                            })?,
                        ),
                        None => (view, 0, *id),
                    };
                original
                    .ports()
                    .iter()
                    .find(|(port, _)| *port == local_id)
                    .map(|(_, column)| (offset, column.as_slice()))
                    .ok_or_else(|| {
                        invalid("continuous terminal channel names an undeclared DG port")
                    })
            })
            .collect::<Result<_, _>>()?;
        let terminal_owner = terminal.map(|(owner, _)| owner);
        let terminal_dimension =
            terminal_owner.map_or(0, |owner| owner.continuous_operator().dimension());
        Ok(Self {
            view: *view,
            branch_views,
            pressure_column,
            branch_offsets,
            branch_ports,
            junction_columns,
            terminal_columns,
            terminal_ports,
            terminal_owner,
            channels,
            acoustic_dimension,
            terminal_dimension,
            channel_dimension,
        })
    }

    pub fn admit_frequency(
        &self,
        omega: f64,
        gate: &CancelGate,
    ) -> Result<(), DgHarmonicCouplingError> {
        poll(Some(gate))?;
        self.view
            .admit_frequency(omega)
            .map_err(|e| duct_error(e, gate))?;
        for view in self.branch_views {
            poll(Some(gate))?;
            view.admit_frequency(omega)
                .map_err(|e| duct_error(e, gate))?;
        }
        if let Some(owner) = self.terminal_owner {
            owner
                .continuous_operator()
                .admit_frequency(omega)
                .map_err(|e| invalid(e.0))?;
        }
        Ok(())
    }

    pub fn additional_dimension(&self) -> Result<usize, ModalCouplingError> {
        self.acoustic_dimension
            .checked_add(self.terminal_dimension)
            .and_then(|v| v.checked_add(self.terminal_ports.len()))
            .and_then(|v| v.checked_add(self.branch_views.len()))
            .and_then(|v| v.checked_add(self.channel_dimension))
            .ok_or_else(|| invalid("continuous terminal harmonic dimension overflow"))
    }

    pub fn append(
        &self,
        offset: usize,
        flow_column: &[f64],
        maximum_dense_bytes: usize,
        matrix: &mut impl Entries,
        gate: &CancelGate,
    ) -> Result<AppendedDg, DgHarmonicCouplingError> {
        let terminal_offset = offset + self.acoustic_dimension;
        let flow_offset = terminal_offset + self.terminal_dimension;
        let junction_offset = flow_offset + self.terminal_ports.len();
        if offset
            .checked_add(self.additional_dimension()?)
            .is_none_or(|end| end > matrix.dimension())
        {
            return Err(invalid("continuous acoustic entries exceed the complete pencil").into());
        }
        let terminal_source = self
            .terminal_owner
            .map(|owner| {
                owner
                    .continuous_operator()
                    .state_space(maximum_dense_bytes, &mut || gate.is_requested())
                    .map_err(|e| {
                        if gate.is_requested() {
                            ModalCouplingError::Cancelled
                        } else {
                            invalid(e.0)
                        }
                    })
            })
            .transpose()?;
        let channel_offsets =
            self.append_source(offset, flow_column, terminal_source.as_ref(), matrix, gate)?;
        Ok(AppendedDg {
            terminal_source,
            terminal_offset,
            flow_offset,
            junction_offset,
            channel_offsets,
        })
    }

    // All AC and tangent entries consume the SAME source coefficients and
    // original algebraic rows. The retained terminal realization is reused;
    // evaluating its derivative never re-realizes or inverts the terminal.
    fn append_source(
        &self,
        offset: usize,
        flow_column: &[f64],
        terminal_source: Option<&MimoStateSpace>,
        matrix: &mut impl Entries,
        gate: &CancelGate,
    ) -> Result<Vec<usize>, DgHarmonicCouplingError> {
        let terminal_offset = offset + self.acoustic_dimension;
        let flow_offset = terminal_offset + self.terminal_dimension;
        let junction_offset = flow_offset + self.terminal_ports.len();
        let channel_offset = junction_offset + self.branch_views.len();
        if offset
            .checked_add(self.additional_dimension()?)
            .is_none_or(|end| end > matrix.dimension())
        {
            return Err(invalid("continuous acoustic entries exceed the complete pencil").into());
        }
        let width = self.view.dimension();
        let source = self.view.generator_per_s();
        for i in 0..width {
            poll(Some(gate))?;
            for j in 0..width {
                matrix.set(
                    offset + i,
                    offset + j,
                    -source[i * width + j],
                    if i == j { 1.0 } else { 0.0 },
                );
            }
            for (j, &b) in flow_column.iter().enumerate() {
                let coupling = self.pressure_column[i] * b;
                matrix.set(j, offset + i, coupling, 0.0);
                matrix.set(offset + i, j, 0.0, -coupling);
            }
        }
        for (index, original) in self.branch_views.iter().enumerate() {
            let start = offset + self.branch_offsets[index];
            let width = original.dimension();
            let source = original.generator_per_s();
            for i in 0..width {
                poll(Some(gate))?;
                for j in 0..width {
                    matrix.set(
                        start + i,
                        start + j,
                        -source[i * width + j],
                        if i == j { 1.0 } else { 0.0 },
                    );
                }
            }
            let junction = junction_offset + index;
            for (i, &b) in self.junction_columns[index].0.iter().enumerate() {
                matrix.set(offset + i, junction, b, 0.0);
                matrix.set(junction, offset + i, b, 0.0);
            }
            for (i, &b) in self.junction_columns[index].1.iter().enumerate() {
                matrix.set(start + i, junction, -b, 0.0);
                matrix.set(junction, start + i, -b, 0.0);
            }
        }
        if let Some(load) = terminal_source {
            let ports = self.terminal_ports.len();
            for (j, &(relative, column)) in self.terminal_columns.iter().enumerate() {
                poll(Some(gate))?;
                for (i, &b) in column.iter().enumerate() {
                    matrix.set(offset + relative + i, flow_offset + j, b, 0.0);
                    matrix.set(flow_offset + j, offset + relative + i, b, 0.0);
                }
            }
            for i in 0..self.terminal_dimension {
                poll(Some(gate))?;
                for j in 0..self.terminal_dimension {
                    matrix.set(
                        terminal_offset + i,
                        terminal_offset + j,
                        -load.a[i * self.terminal_dimension + j],
                        if i == j { 1.0 } else { 0.0 },
                    );
                }
                for j in 0..ports {
                    matrix.set(
                        terminal_offset + i,
                        flow_offset + j,
                        -load.b[i * ports + j],
                        0.0,
                    );
                    matrix.set(
                        flow_offset + j,
                        terminal_offset + i,
                        -load.c[j * self.terminal_dimension + i],
                        0.0,
                    );
                }
            }
            for i in 0..ports {
                for j in 0..ports {
                    matrix.set(
                        flow_offset + i,
                        flow_offset + j,
                        -load.d[i * ports + j],
                        0.0,
                    );
                }
            }
        }
        let mut channel_offsets = Vec::with_capacity(self.channels.len());
        let mut next_channel = channel_offset;
        for channel in self.channels {
            poll(Some(gate))?;
            channel_offsets.push(next_channel);
            let port = self
                .terminal_ports
                .iter()
                .position(|id| *id == channel.input().port)
                .ok_or_else(|| invalid("continuous channel lost its original terminal row"))?;
            let flow_row = flow_offset + port;
            if let Some(original) = channel.channel() {
                let field = original.continuous_operator()?;
                let width = field.dimension();
                let drop = next_channel + width;
                let mass = field.cell_masses_kg();
                let area = field.cell_areas_m2();
                let shear = field.face_conductances_kg_s();
                matrix.set(flow_row, drop, -1.0, 0.0);
                matrix.set(drop, flow_row, 1.0, 0.0);
                for i in 0..width {
                    poll(Some(gate))?;
                    let row = next_channel + i;
                    matrix.set(row, row, shear[i] + shear[i + 1], mass[i]);
                    if i != 0 {
                        matrix.set(row, row - 1, -shear[i], 0.0);
                    }
                    if i + 1 != width {
                        matrix.set(row, row + 1, -shear[i + 1], 0.0);
                    }
                    matrix.set(row, drop, -area[i], 0.0);
                    matrix.set(drop, row, -area[i], 0.0);
                }
                next_channel = drop + 1;
            } else {
                matrix.clear_row(flow_row);
                matrix.set(flow_row, flow_row, 1.0, 0.0);
            }
        }
        Ok(channel_offsets)
    }

    pub fn apply_frequency_derivative(
        &self,
        offset: usize,
        flow_column: &[f64],
        appended: &AppendedDg,
        state: &[C64],
        rhs: &mut [C64],
        gate: &CancelGate,
    ) -> Result<(), DgHarmonicCouplingError> {
        poll(Some(gate))?;
        if state.len() != rhs.len() || state.iter().chain(rhs.iter()).any(|&v| !complex_finite(v)) {
            return Err(
                invalid("complete finite joint state and frequency forcing required").into(),
            );
        }
        self.append_source(
            offset,
            flow_column,
            appended.terminal_source.as_ref(),
            &mut FrequencyAction { rhs, state },
            gate,
        )?;
        if rhs.iter().any(|&v| !complex_finite(v)) {
            return Err(invalid("original acoustic frequency tangent forcing overflows").into());
        }
        poll(Some(gate))?;
        Ok(())
    }
}
