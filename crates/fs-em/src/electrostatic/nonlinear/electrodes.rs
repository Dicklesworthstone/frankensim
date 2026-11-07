//! Prepared physical electrode traces and their voltage/total-charge quotient.
use super::*;

/// Original terminal control: voltage [V], or total conductor free charge [C].
/// Charge control does not prescribe local normal displacement or zero voltage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ElectrodeControl {
    Voltage(f64),
    Charge(f64),
}

pub(super) struct Layout {
    pub terminals: Vec<Vec<usize>>,
    pub node_terminal: Vec<Option<usize>>,
    pub components: Vec<usize>,
}

/// Independent gauges, retained separately from the small physical voltages.
/// A conductor joining mesh pieces gives those pieces one electrical component.
#[derive(Clone)]
pub(super) struct VoltageReferences {
    pub first_terminal: usize,
    pub first_voltage_v: f64,
    pub terminal_by_vertex: Vec<usize>,
    pub voltage_by_vertex_v: Vec<f64>,
}

impl Layout {
    pub fn voltage_references(
        &self,
        controls: &[ElectrodeControl],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VoltageReferences, Error> {
        if controls.len() != self.terminals.len()
            || controls.iter().any(|c| match c {
                ElectrodeControl::Voltage(v) | ElectrodeControl::Charge(v) => !v.is_finite(),
            })
        {
            return Err(invalid("complete finite electrode controls required"));
        }
        let mut anchors = vec![None; self.components.len()];
        let mut first = None;
        for (t, control) in controls.iter().enumerate() {
            poll(cancelled)?;
            if let ElectrodeControl::Voltage(voltage) = *control {
                first.get_or_insert((t, voltage));
                anchors[self.components[self.terminals[t][0]]].get_or_insert((t, voltage));
            }
        }
        let missing_reference =
            || invalid("every dielectric component requires a voltage-controlled reference");
        let (first_terminal, first_voltage_v) = first.ok_or_else(missing_reference)?;
        let mut terminal_by_vertex = Vec::with_capacity(self.components.len());
        let mut voltage_by_vertex_v = Vec::with_capacity(self.components.len());
        for &component in &self.components {
            poll(cancelled)?;
            let (terminal, voltage) = anchors[component].ok_or_else(missing_reference)?;
            terminal_by_vertex.push(terminal);
            voltage_by_vertex_v.push(voltage);
        }
        Ok(VoltageReferences {
            first_terminal,
            first_voltage_v,
            terminal_by_vertex,
            voltage_by_vertex_v,
        })
    }
}

/// Original mesh/terminal preparation reused across bias and charge sweeps.
/// Returned fields retain this immutable layout and their original materials.
pub struct NonlinearElectrodes<'m> {
    pub(super) problem: &'m PlanarNonlinearElectrostatics,
    pub(super) layout: Arc<Layout>,
}

pub(super) struct Dofs {
    pub index: Vec<Option<usize>>,
    pub representatives: Vec<usize>,
    pub charge: Vec<f64>,
}
impl Dofs {
    pub fn len(&self) -> usize {
        self.representatives.len()
    }
    pub fn residual(
        &self,
        evaluation: &Evaluation,
        rhs: &[f64],
        out: &mut [f64],
        scratch: &mut [f64],
    ) -> Result<(), Error> {
        out.fill(0.0);
        scratch.fill(0.0);
        for (v, index) in self.index.iter().enumerate() {
            if let Some(i) = *index {
                add(&mut out[i], &mut scratch[i], evaluation.nodal[v] - rhs[v])?;
            }
        }
        for i in 0..out.len() {
            out[i] = finite(out[i] + scratch[i] - self.charge[i])?;
        }
        Ok(())
    }
}

impl PlanarNonlinearElectrostatics {
    /// Admit disjoint unions of complete exterior edges through the original
    /// geometric owner. A single conductor can join disconnected mesh pieces.
    /// The voltage references are checked when each control set is supplied.
    pub fn prepare_electrodes<'m>(
        &'m self,
        terminals: &[&[usize]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NonlinearElectrodes<'m>, Error> {
        poll(cancelled)?;
        if matches!(self.mesh.metric(), Metric2::Axisymmetric { .. }) {
            super::super::axisymmetric::admit_terminal_traces(&self.mesh, terminals, cancelled)?;
        } else {
            scalar_potential::terminal_boundary(
                &self.mesh,
                terminals,
                &vec![0.0; terminals.len()],
                cancelled,
            )?;
        }
        let components = scalar_potential::terminal_components(&self.mesh, terminals, cancelled)?;
        let mut node_terminal = vec![None; self.mesh.vertices().len()];
        for (t, terminal) in terminals.iter().enumerate() {
            for &v in *terminal {
                poll(cancelled)?;
                node_terminal[v] = Some(t);
            }
        }
        poll(cancelled)?;
        Ok(NonlinearElectrodes {
            problem: self,
            layout: Arc::new(Layout {
                terminals: terminals.iter().map(|v| v.to_vec()).collect(),
                node_terminal,
                components,
            }),
        })
    }
}

impl<'m> NonlinearElectrodes<'m> {
    pub fn problem(&self) -> &'m PlanarNonlinearElectrostatics {
        self.problem
    }
    pub fn terminals(&self) -> &[Vec<usize>] {
        &self.layout.terminals
    }

    pub(super) fn coordinates(
        &self,
        controls: &[ElectrodeControl],
        initial: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Dofs, Vec<f64>, VoltageReferences), Error> {
        let n = self.problem.mesh.vertices().len();
        if initial.len() != n || initial.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "complete finite relative initial potentials required",
            ));
        }
        let references = self.layout.voltage_references(controls, cancelled)?;
        for (t, control) in controls.iter().enumerate() {
            poll(cancelled)?;
            let trace = &self.layout.terminals[t];
            if matches!(control, ElectrodeControl::Charge(_))
                && trace.iter().any(|&v| initial[v] != initial[trace[0]])
            {
                return Err(invalid(
                    "a floating trace needs one shared initial potential",
                ));
            }
        }
        let mut index = vec![None; n];
        let mut representatives = Vec::new();
        let mut charge = Vec::new();
        let mut floating = vec![None; controls.len()];
        let mut potential = initial.to_vec();
        for v in 0..n {
            poll(cancelled)?;
            if let Some(t) = self.layout.node_terminal[v] {
                match controls[t] {
                    ElectrodeControl::Voltage(value) => {
                        let difference = finite(value - references.voltage_by_vertex_v[v])?;
                        potential[v] = difference;
                    }
                    ElectrodeControl::Charge(value) => {
                        let i = if let Some(i) = floating[t] {
                            i
                        } else {
                            let i = representatives.len();
                            representatives.push(v);
                            charge.push(value);
                            floating[t] = Some(i);
                            i
                        };
                        index[v] = Some(i);
                    }
                }
            } else {
                index[v] = Some(representatives.len());
                representatives.push(v);
                charge.push(0.0);
            }
        }
        Ok((
            Dofs {
                index,
                representatives,
                charge,
            },
            potential,
            references,
        ))
    }

    /// Solve the original nonlinear Gauss equations with mixed voltage/charge
    /// electrodes. Initial potentials are relative to the first voltage control
    /// in each electrical component; all absolute references are retained
    /// separately. Prescribed nodes take their component-relative boundary
    /// value, and each floating trace has one shared guess.
    /// Every trial must stay inside all material field supports. No accepted
    /// history is mutated and no unconverged/partially admitted field is returned.
    pub fn solve(
        &self,
        free_charge_density_c_m3: &[f64],
        controls: &[ElectrodeControl],
        initial_relative_voltage_v: &[f64],
        solve: NonlinearElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NonlinearElectricField<'m>, Error> {
        solve::solve(
            self,
            free_charge_density_c_m3,
            controls,
            initial_relative_voltage_v,
            solve,
            cancelled,
        )
    }
}
