//! Full grounded field capacitance in the existing native Kirchhoff owner.
//!
//! Each nonreference electrode retains its TOTAL free-charge coordinate.
//! Q=C V and H=Q^T C^-1 Q/2 use the original field matrix, with no pairwise
//! fitting or sign clipping. A floating terminal remains a native node whose
//! KCL conserves charge; it is never eliminated with its charged energy lost.

use super::{CapacitanceResponse, PlanarCapacitanceMatrix, add, invalid};
use crate::electrostatic::{ElectrostaticError, finite, poll, positive};
use fs_phs::circuit::{Branch, capacitance::CapacitanceBlock};

/// Source-retaining field export for a full native capacitance block.
/// The returned branches and block jointly describe ONE dielectric store.
/// Branch diagonal values alone are not an independent capacitor network.
pub struct CoupledCapacitanceCircuit<'s, 'm> {
    source: &'s PlanarCapacitanceMatrix<'m>,
    reference_terminal: usize,
    coordinate_terminals: Vec<usize>,
    branches: Vec<(usize, usize, Branch)>,
    block: CapacitanceBlock,
    initial_charges_c: Vec<f64>,
    response: CapacitanceResponse,
    maximum_charge_residual_c: f64,
    energy_residual_j: f64,
}

impl<'m> PlanarCapacitanceMatrix<'m> {
    /// Transfer the complete grounded matrix to native electrode-charge states.
    /// Original terminal order is retained after removing the explicit reference
    /// row/column. Branches run from each retained terminal to that reference;
    /// their coordinates are total electrode charges, not fitted pair charges.
    /// Positive off-diagonal entries are admitted when the native block is SPD.
    ///
    /// The complete node map must be unique. All dielectric components, after
    /// joining pieces connected by the same equipotential electrode, must have
    /// this reference. Additional gauge nullspaces refuse before factorization.
    /// Charged floating terminals remain represented: connect no external
    /// current to their native node to conserve their original total charge.
    ///
    /// Native charge and H are checked against independent original weak
    /// reactions and superposed cell-field energy at the actual voltage. Source
    /// material and geometry stay frozen; parameter or mechanical work requires
    /// its own coupled owner. Append BOTH branches and the remapped block to the
    /// circuit, without another parallel copy of this dielectric storage.
    pub fn coupled_circuit_at_voltages<'s>(
        &'s self,
        reference_terminal: usize,
        terminal_nodes: &[usize],
        voltage_v: &[f64],
        charge_tolerance_c: f64,
        energy_tolerance_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CoupledCapacitanceCircuit<'s, 'm>, ElectrostaticError> {
        poll(cancelled)?;
        let count = self.terminals.len();
        if reference_terminal >= count
            || terminal_nodes.len() != count
            || voltage_v.len() != count
            || voltage_v.iter().any(|v| !v.is_finite())
            || !positive(charge_tolerance_c)
            || !positive(energy_tolerance_j)
            || terminal_nodes
                .iter()
                .enumerate()
                .any(|(i, &node)| node == usize::MAX || terminal_nodes[..i].contains(&node))
        {
            return Err(invalid(
                "coupled capacitance needs a complete unique node map, explicit reference, finite voltages and positive physical tolerances",
            ));
        }
        let terminals: Vec<_> = self.terminals.iter().map(Vec::as_slice).collect();
        let components = crate::scalar_potential::terminal_components(
            self.problem.mesh(),
            &terminals,
            cancelled,
        )?;
        let reference_component = components[self.terminals[reference_terminal][0]];
        for component in components {
            poll(cancelled)?;
            if component != reference_component {
                return Err(invalid(
                    "every coupled dielectric component needs the explicit circuit reference",
                ));
            }
        }
        let coordinate_terminals: Vec<_> = (0..count)
            .filter(|&terminal| terminal != reference_terminal)
            .collect();
        let size = coordinate_terminals.len();
        poll(cancelled)?;
        // One native preparation owns both SPD admission and the shared inverse.
        let block = CapacitanceBlock::try_new(
            (0..size).collect(), self.grounded_matrix_f(reference_terminal)?,
        ).map_err(|_| invalid("original grounded capacitance cannot form a representable positive-definite native block"))?;
        poll(cancelled)?;
        let relative_voltage: Vec<_> = coordinate_terminals
            .iter()
            .map(|&terminal| finite(voltage_v[terminal] - voltage_v[reference_terminal]))
            .collect::<Result<_, _>>()?;
        let initial_charges_c = block
            .charges_c(&relative_voltage)
            .map_err(|_| invalid("native coupled electrode charge is unrepresentable"))?;
        let native_energy = block
            .stored_energy_j(&initial_charges_c)
            .map_err(|_| invalid("native coupled electric energy is unrepresentable"))?;
        let response = self.response(voltage_v, cancelled)?;
        let mut terminal_charge = vec![0.0; count];
        let (mut reference_charge, mut correction) = (0.0, 0.0);
        let mut branches = Vec::with_capacity(size);
        for (coordinate, &terminal) in coordinate_terminals.iter().enumerate() {
            poll(cancelled)?;
            terminal_charge[terminal] = initial_charges_c[coordinate];
            add(
                &mut reference_charge,
                &mut correction,
                -initial_charges_c[coordinate],
            )?;
            branches.push((
                terminal_nodes[terminal],
                terminal_nodes[reference_terminal],
                Branch::Capacitor {
                    farads: block.matrix_f()[coordinate * size + coordinate],
                },
            ));
        }
        terminal_charge[reference_terminal] = reference_charge;
        let mut maximum_charge_residual_c = 0.0_f64;
        for (&native, &field) in terminal_charge.iter().zip(&response.terminal_charge_c) {
            maximum_charge_residual_c =
                maximum_charge_residual_c.max(finite(native - field)?.abs());
        }
        let energy_residual_j = finite(native_energy - response.energy_j)?;
        if maximum_charge_residual_c > charge_tolerance_c
            || energy_residual_j.abs() > energy_tolerance_j
        {
            return Err(invalid(
                "native coupled charges or stored energy disagree with the original dielectric field",
            ));
        }
        poll(cancelled)?;
        Ok(CoupledCapacitanceCircuit {
            source: self,
            reference_terminal,
            coordinate_terminals,
            branches,
            block,
            initial_charges_c,
            response,
            maximum_charge_residual_c,
            energy_residual_j,
        })
    }
}

impl<'s, 'm> CoupledCapacitanceCircuit<'s, 'm> {
    pub fn source(&self) -> &'s PlanarCapacitanceMatrix<'m> {
        self.source
    }
    pub fn reference_terminal(&self) -> usize {
        self.reference_terminal
    }
    /// Original terminal index for each native block charge, in order.
    pub fn coordinate_terminals(&self) -> &[usize] {
        &self.coordinate_terminals
    }
    /// Append these original node-map branches together and in this order.
    pub fn branches(&self) -> &[(usize, usize, Branch)] {
        &self.branches
    }
    /// Prepared block for returned branches starting at graph branch zero.
    /// Use `block_at_branch_offset` when other branches precede this export.
    pub fn capacitance_block(&self) -> &CapacitanceBlock {
        &self.block
    }
    /// Reuse the exact prepared C/inverse after appending at a graph offset.
    /// This changes graph indices only; it does not refactor or refit the field.
    pub fn block_at_branch_offset(
        &self,
        first_branch: usize,
    ) -> Result<CapacitanceBlock, ElectrostaticError> {
        let branches = (0..self.branches.len())
            .map(|i| {
                first_branch
                    .checked_add(i)
                    .filter(|&index| index != usize::MAX)
                    .ok_or_else(|| invalid("coupled capacitance graph branch offset overflows"))
            })
            .collect::<Result<_, _>>()?;
        self.block
            .with_branches(branches)
            .map_err(|_| invalid("coupled capacitance graph branch remapping is invalid"))
    }
    /// Physical electrode charges in returned branch order. After full graph
    /// assembly place them at the corresponding native `charge_index` entries.
    pub fn initial_charges_c(&self) -> &[f64] {
        &self.initial_charges_c
    }
    pub fn field_response(&self) -> &CapacitanceResponse {
        &self.response
    }
    pub fn maximum_charge_residual_c(&self) -> f64 {
        self.maximum_charge_residual_c
    }
    pub fn energy_residual_j(&self) -> f64 {
        self.energy_residual_j
    }
}
