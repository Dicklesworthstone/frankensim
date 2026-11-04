//! Lower field-derived self/mutual inductance into the original circuit DAE.
//!
//! Winding current orientation is the original graph branch's first-to-second
//! terminal direction. Total linkage is lambda=lambda0+L*i; the native circuit
//! retains incremental delta_lambda=L*i, including declared end/leakage
//! self-inductance. Fixed lambda0 creates no spurious voltage or current. All capacitor
//! storage, Kirchhoff stamps, algebraic constraints, ports, initial-condition
//! solves and time stepping remain owned by fs-phs. This module adds no circuit
//! solver or time integrator.

use super::{PlanarWindingBasis, WindingError};
use crate::{
    planar::poll,
    transfer::{CellSample, PlanarFrame},
};
use fs_phs::{
    DescriptorPortHamiltonian, QuadraticStorage, StepRecord,
    circuit::{Branch, CircuitDae, CircuitGraph, assemble_circuit},
};

/// The original circuit DAE with the exact field/material basis retained.
/// Differential magnetic coordinates are INCREMENTAL linkage relative to the
/// fixed zero-current magnet state. They map from original branch indices;
/// conjugate currents include every mutual term, not delta_lambda_i/L_ii.
pub struct WindingCircuit<'a> {
    basis: PlanarWindingBasis<'a>,
    dae: CircuitDae,
    graph: CircuitGraph,
    branch_indices: Vec<usize>,
    flux_coordinates: Vec<usize>,
}
impl<'a> PlanarWindingBasis<'a> {
    /// Bind EVERY inductor branch once, in winding order. Existing numeric
    /// branch inductances are replaced by field-derived diagonal values before
    /// Kirchhoff assembly; the FULL inverse magnetic matrix then replaces the
    /// original flux storage block. Capacitors and all algebraic rows retain
    /// their original ordering and ownership.
    ///
    /// This allocating dense circuit lane is bounded to 32 nodes, 64 branches,
    /// 16 transformers, 32 source ports and 128 extended states. A singular
    /// winding matrix was already refused during field extraction; ideal
    /// perfectly coupled transformers use the existing TransformerLink lane.
    pub fn into_circuit(
        self,
        mut graph: CircuitGraph,
        branch_indices: Vec<usize>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WindingCircuit<'a>, WindingError> {
        poll(cancelled)?;
        if graph.node_count > 32
            || graph.branches.len() > 64
            || graph.transformers.len() > 16
            || branch_indices.len() != self.winding_count()
        {
            return Err(WindingError::Invalid(
                "field circuit exceeds capacity or has incomplete bindings",
            ));
        }
        let inductors: Vec<_> = graph
            .branches
            .iter()
            .enumerate()
            .filter_map(|(i, (_, _, branch))| {
                matches!(branch, Branch::Inductor { .. }).then_some(i)
            })
            .collect();
        if inductors.len() != self.winding_count() {
            return Err(WindingError::Invalid(
                "every original inductor needs one field winding",
            ));
        }
        let mut bound = vec![false; graph.branches.len()];
        let mut slots = Vec::with_capacity(self.winding_count());
        let n = self.winding_count();
        for (winding, &branch) in branch_indices.iter().enumerate() {
            poll(cancelled)?;
            let slot = inductors.iter().position(|&index| index == branch).ok_or(
                WindingError::Invalid("winding binding does not name an inductor"),
            )?;
            if bound[branch] {
                return Err(WindingError::Invalid("duplicate winding branch binding"));
            }
            bound[branch] = true;
            slots.push(slot);
            graph.branches[branch].2 = Branch::Inductor {
                henries: self.inductance_h()[winding * n + winding],
            };
        }
        let mut capacitors = 0;
        let mut voltage_sources = 0;
        for &(_, _, branch) in &graph.branches {
            poll(cancelled)?;
            match branch {
                Branch::Capacitor { farads } => {
                    if !(farads.is_finite() && farads > 0.0 && farads.recip().is_finite()) {
                        return Err(WindingError::Invalid(
                            "capacitive storage must be representable",
                        ));
                    }
                    capacitors += 1;
                }
                Branch::Resistor { ohms } => {
                    if !(ohms.is_finite() && ohms > 0.0 && ohms.recip().is_finite()) {
                        return Err(WindingError::Invalid(
                            "resistive conductance must be representable",
                        ));
                    }
                }
                Branch::VoltageSource { port } | Branch::CurrentSource { port } => {
                    if port >= 32 {
                        return Err(WindingError::Invalid(
                            "source port exceeds bounded circuit capacity",
                        ));
                    }
                    if matches!(branch, Branch::VoltageSource { .. }) {
                        voltage_sources += 1;
                    }
                }
                Branch::Inductor { .. } => {}
            }
        }
        let expected_states = n
            + 2 * capacitors
            + graph.node_count.saturating_sub(1)
            + voltage_sources
            + graph.transformers.len();
        if expected_states > 128 {
            return Err(WindingError::Invalid(
                "extended circuit state exceeds bounded capacity",
            ));
        }
        let original = assemble_circuit(&graph)?;
        let state_count = original.system.state_dim();
        let differential_count = original.system.differential_dim();
        let flux_coordinates: Vec<_> = slots
            .iter()
            .map(|&slot| original.flux_index[slot])
            .collect();
        let original_q =
            original
                .system
                .storage()
                .exact_quadratic_matrix()
                .ok_or(WindingError::Invalid(
                    "original circuit did not declare complete quadratic storage",
                ))?;
        let mut q = vec![0.0; differential_count * differential_count];
        for i in 0..differential_count {
            q[i * differential_count..(i + 1) * differential_count].copy_from_slice(
                &original_q[i * state_count..i * state_count + differential_count],
            );
        }
        for i in 0..n {
            for j in 0..n {
                q[flux_coordinates[i] * differential_count + flux_coordinates[j]] =
                    self.inverse_inductance_per_h()[i * n + j];
            }
        }
        if q.iter().any(|v| !v.is_finite()) {
            return Err(WindingError::Invalid(
                "coupled differential storage is not representable",
            ));
        }
        let storage = QuadraticStorage::new(q, differential_count)?;
        let (j, r, g) = original.system.structure();
        let system = DescriptorPortHamiltonian::new(
            state_count,
            differential_count,
            original.ports,
            j.to_vec(),
            r.to_vec(),
            g.to_vec(),
            Box::new(storage),
        )?;
        let dae = CircuitDae {
            system,
            flux_index: original.flux_index,
            charge_index: original.charge_index,
            node_potential_index: original.node_potential_index,
            ports: original.ports,
        };
        poll(cancelled)?;
        Ok(WindingCircuit {
            basis: self,
            dae,
            graph,
            branch_indices,
            flux_coordinates,
        })
    }
}
impl WindingCircuit<'_> {
    /// Source-resolved field geometry, winding orientation and retained basis.
    pub fn basis(&self) -> &PlanarWindingBasis<'_> {
        &self.basis
    }
    /// Existing circuit owner in incremental linkage coordinates. Its energy
    /// is relative to the zero-current magnet state. For permanent magnets use
    /// this wrapper's checked initialization/step so physical recoil envelopes
    /// are admitted before publishing a state; the raw DAE has no field-domain
    /// policy of its own.
    pub fn dae(&self) -> &CircuitDae {
        &self.dae
    }
    /// Original graph topology and authored resistor/source/capacitor data.
    /// Inductor placeholders carry the extracted diagonal; full mutual storage
    /// remains in the original descriptor.
    pub fn graph(&self) -> &CircuitGraph {
        &self.graph
    }
    /// Original graph branch index corresponding to each winding.
    pub fn branch_indices(&self) -> &[usize] {
        &self.branch_indices
    }
    /// Original extended-state INCREMENTAL flux coordinate for each winding.
    /// Add `basis().bias_flux_linkages_wb()` only when reporting total linkage.
    pub fn flux_coordinates(&self) -> &[usize] {
        &self.flux_coordinates
    }

    /// Initialize from winding-order terminal currents and original capacitor-
    /// branch-order charges. The existing consistency solve changes only
    /// algebraic coordinates; incompatible supplied charge/source data refuse.
    pub fn consistent_initial_state(
        &self,
        currents_a: &[f64],
        charges_c: &[f64],
        inputs: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, WindingError> {
        poll(cancelled)?;
        let flux = self.basis.incremental_flux_linkages_wb(currents_a)?;
        if charges_c.len() != self.dae.charge_index.len()
            || inputs.len() != self.dae.ports
            || charges_c.iter().chain(inputs).any(|v| !v.is_finite())
        {
            return Err(WindingError::Invalid(
                "complete finite capacitor charges and source ports are required",
            ));
        }
        self.basis.admit_currents(currents_a, cancelled)?;
        let mut state = vec![0.0; self.dae.system.state_dim()];
        for (&index, &value) in self.flux_coordinates.iter().zip(&flux) {
            state[index] = value;
        }
        for (&index, &value) in self.dae.charge_index.iter().zip(charges_c) {
            state[index] = value;
        }
        let result = self.dae.consistent_initial_state(&state, inputs)?;
        self.admit_state(&result, cancelled)?;
        poll(cancelled)?;
        Ok(result)
    }
    /// Read physical terminal currents from a complete original circuit state.
    pub fn currents_a(&self, state: &[f64]) -> Result<Vec<f64>, WindingError> {
        self.basis
            .currents_from_incremental_flux_wb(&self.incremental_flux_linkages_wb(state)?)
    }
    /// Native delta_lambda=L*i [Wb], in winding order. Keeping this coordinate
    /// avoids subtracting large permanent-flux offsets during every step.
    pub fn incremental_flux_linkages_wb(&self, state: &[f64]) -> Result<Vec<f64>, WindingError> {
        if state.len() != self.dae.system.state_dim() || state.iter().any(|v| !v.is_finite()) {
            return Err(WindingError::Invalid(
                "complete finite original circuit state is required",
            ));
        }
        Ok(self.flux_coordinates.iter().map(|&i| state[i]).collect())
    }
    /// Total physical lambda=lambda0+delta_lambda [Wb], in winding order.
    pub fn total_flux_linkages_wb(&self, state: &[f64]) -> Result<Vec<f64>, WindingError> {
        let mut total = self.incremental_flux_linkages_wb(state)?;
        for (value, &bias) in total.iter_mut().zip(self.basis.bias_flux_linkages_wb()) {
            *value = super::finite(*value + bias)?;
        }
        Ok(total)
    }
    /// Native magnetic/capacitor storage plus the constant zero-current field
    /// reference W0 [J]. Stepping uses the incremental energy directly so a
    /// large W0 cannot obscure small electrical work through subtraction.
    pub fn stored_energy_with_bias_j(&self, state: &[f64]) -> Result<f64, WindingError> {
        self.incremental_flux_linkages_wb(state)?;
        super::finite(self.dae.system.hamiltonian(state) + self.basis.bias_energy_j())
    }
    /// Validate actual reconstructed B/H on every original material cell.
    /// This is source-domain admission, not a Kirchhoff consistency solve.
    pub fn admit_state(
        &self,
        state: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WindingError> {
        self.basis
            .admit_currents(&self.currents_a(state)?, cancelled)
    }
    /// Advance through the existing audited circuit DAE, then admit the whole
    /// endpoint field before returning it. Fixed recoil envelopes are convex
    /// in current; admitted endpoints also admit their linear interpolation.
    /// This does not bound the continuous physical trajectory between steps.
    /// Failure/cancellation returns no accepted replacement for the input.
    pub fn step_audited(
        &self,
        state: &[f64],
        inputs: &[f64],
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(StepRecord, f64), WindingError> {
        poll(cancelled)?;
        self.admit_state(state, cancelled)?;
        let candidate = self.dae.step_audited(state, inputs, dt)?;
        self.admit_state(&candidate.0.x, cancelled)?;
        poll(cancelled)?;
        Ok(candidate)
    }
    /// Reconstruct the field of an actual circuit state and use the original
    /// explicit geometry/frame transfer into moving-conductor/flow consumers.
    pub fn sample_field(
        &self,
        state: &[f64],
        frame: &PlanarFrame,
        samples: &[CellSample],
        position_tolerance_m: f64,
        max_samples: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, WindingError> {
        self.basis.sample_for_currents(
            &self.currents_a(state)?,
            frame,
            samples,
            position_tolerance_m,
            max_samples,
            cancelled,
        )
    }
}
