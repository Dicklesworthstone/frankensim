//! Physical conductor charges on the original axisymmetric dielectric field.
//!
//! The shared scalar owner performs the native field solves, capacitance Gram
//! contraction, floating-block factorization and independent charge/energy
//! checks. All integrated quantities include the original angular span. These
//! wrappers preserve the radial geometry type through recovered fields and
//! native circuit storage; they expose no planar mechanical-load shortcut.

use super::{
    AxisymmetricElectricField, AxisymmetricElectrostatics, ElectrostaticError,
    admit_terminal_traces, poll,
};
use crate::electrostatic::network::{
    CoupledCapacitanceCircuit, PlanarCapacitanceMatrix,
    floating::{FloatingCapacitance, FloatingCapacitanceResponse},
};
use fs_phs::circuit::{Branch, capacitance::CapacitanceBlock};

pub use crate::electrostatic::network::{
    CapacitanceMatrixControls, CapacitanceResponse, floating::FloatingCapacitanceControls,
};

/// Complete Maxwell matrix with the original radial field and material owner.
/// No prescribed volume/sheet charge or impressed polarization is folded into
/// this homogeneous capacitance response. Such sources use the direct solve.
pub struct AxisymmetricCapacitanceMatrix<'m> {
    problem: &'m AxisymmetricElectrostatics,
    core: PlanarCapacitanceMatrix<'m>,
}

impl AxisymmetricElectrostatics {
    /// Extract all named conductors, including any eventual reference or
    /// floating conductor. Every conductor consists of disjoint complete
    /// exterior edges with positive revolved surface area. Unassigned physical
    /// boundaries are insulated; the symmetry axis remains a regularity edge.
    ///
    /// The original bounded unit-field preparation retains material sources,
    /// radial volume weights, component gauges and independent weak reactions.
    /// Geometry and dielectric state stay fixed throughout subsequent queries.
    pub fn capacitance_matrix(
        &self,
        terminals: &[&[usize]],
        controls: CapacitanceMatrixControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricCapacitanceMatrix<'_>, ElectrostaticError> {
        poll(cancelled)?;
        if terminals.len() > controls.maximum_terminals
            || terminals.len() > controls.maximum_field_solves
        {
            return Err(ElectrostaticError::Invalid(
                "axisymmetric electrode count exceeds field preparation budget",
            ));
        }
        admit_terminal_traces(self.mesh(), terminals, cancelled)?;
        Ok(AxisymmetricCapacitanceMatrix {
            problem: self,
            core: self
                .core
                .capacitance_matrix(terminals, controls, cancelled)?,
        })
    }
}

impl<'m> AxisymmetricCapacitanceMatrix<'m> {
    pub fn problem(&self) -> &'m AxisymmetricElectrostatics {
        self.problem
    }
    pub fn terminals(&self) -> &[Vec<usize>] {
        self.core.terminals()
    }
    /// Original terminal order, row-major square matrix [F], for the complete
    /// authored angular span. Each connected electrical component has its own
    /// common-potential nullspace; no implicit ground is removed here.
    pub fn maxwell_matrix_f(&self) -> &[f64] {
        self.core.maxwell_matrix_f()
    }
    pub fn maximum_coefficient_residual_f(&self) -> f64 {
        self.core.maximum_coefficient_residual_f()
    }
    pub fn maximum_gauge_residual_f(&self) -> f64 {
        self.core.maximum_gauge_residual_f()
    }
    /// Remove only the named reference row and column. A floating conductor
    /// instead requires `prepare_floating` and its total-charge constraint.
    pub fn grounded_matrix_f(&self, reference: usize) -> Result<Vec<f64>, ElectrostaticError> {
        self.core.grounded_matrix_f(reference)
    }
    /// Recover the original radial field without another FE solve. Potentials
    /// are relative to each electrical component's first terminal, so separate
    /// common biases never enter the field, charge or energy contractions.
    /// Independent actual-amplitude C/J tolerances apply after superposition.
    pub fn field_at_voltages(
        &self,
        voltage_v: &[f64],
        charge_tolerance_c: f64,
        energy_tolerance_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricElectricField<'m>, ElectrostaticError> {
        Ok(AxisymmetricElectricField {
            problem: self.problem,
            field: self.core.field_at_voltages(
                voltage_v,
                charge_tolerance_c,
                energy_tolerance_j,
                cancelled,
            )?,
        })
    }
    /// Terminal charges [C] and original cell-field storage [J], for the
    /// authored angular span, obtained from the retained unit fields.
    pub fn response(
        &self,
        voltage_v: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CapacitanceResponse, ElectrostaticError> {
        self.core.response(voltage_v, cancelled)
    }
    /// Factor the floating block once. Total floating charges use the supplied
    /// index order; retained electrodes keep their original order. Shared
    /// conductors may join separate dielectric pieces, but every resulting
    /// electrical component must retain a prescribed voltage reference.
    pub fn prepare_floating<'s>(
        &'s self,
        floating_terminals: &[usize],
        controls: FloatingCapacitanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricFloatingCapacitance<'s, 'm>, ElectrostaticError> {
        Ok(AxisymmetricFloatingCapacitance {
            source: self,
            core: self
                .core
                .prepare_floating(floating_terminals, controls, cancelled)?,
        })
    }
    /// Export the full grounded capacitance into the native circuit's
    /// electrode-charge coordinates. The named reference must connect every
    /// electrical component. A floating electrode remains a native node with
    /// conserved charge when no external current is connected to it.
    ///
    /// Append BOTH returned branches and their capacitance block. Together
    /// they represent this one dielectric store; do not add its field energy
    /// again. No pairwise fit or sign clipping replaces the original matrix.
    pub fn coupled_circuit_at_voltages<'s>(
        &'s self,
        reference_terminal: usize,
        terminal_nodes: &[usize],
        voltage_v: &[f64],
        charge_tolerance_c: f64,
        energy_tolerance_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricCoupledCapacitanceCircuit<'s, 'm>, ElectrostaticError> {
        Ok(AxisymmetricCoupledCapacitanceCircuit {
            source: self,
            core: self.core.coupled_circuit_at_voltages(
                reference_terminal,
                terminal_nodes,
                voltage_v,
                charge_tolerance_c,
                energy_tolerance_j,
                cancelled,
            )?,
        })
    }
}

/// Prepared total-charge constraint on the original axisymmetric conductors.
pub struct AxisymmetricFloatingCapacitance<'s, 'm> {
    source: &'s AxisymmetricCapacitanceMatrix<'m>,
    core: FloatingCapacitance<'s, 'm>,
}

impl<'s, 'm> AxisymmetricFloatingCapacitance<'s, 'm> {
    pub fn source(&self) -> &'s AxisymmetricCapacitanceMatrix<'m> {
        self.source
    }
    pub fn retained_terminals(&self) -> &[usize] {
        self.core.retained_terminals()
    }
    pub fn floating_terminals(&self) -> &[usize] {
        self.core.floating_terminals()
    }
    /// Differential Maxwell matrix at fixed floating charge [F]. Nonzero
    /// trapped charge also carries an affine charge offset and stored energy;
    /// use `solve` to retain them instead of treating this matrix as the state.
    pub fn reduced_matrix_f(&self) -> &[f64] {
        self.core.reduced_matrix_f()
    }
    pub fn maximum_coefficient_residual_f(&self) -> f64 {
        self.core.maximum_coefficient_residual_f()
    }
    /// Reuse the prepared factor and original fields at the requested total
    /// floating charges. Charge is for the authored angular span, not per
    /// radian. The native actual-field charge and energy gates precede return;
    /// cancellation or refusal leaves the preparation reusable.
    pub fn solve(
        &self,
        retained_voltage_v: &[f64],
        floating_charge_c: &[f64],
        charge_tolerance_c: f64,
        energy_tolerance_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricFloatingCapacitanceResponse<'m>, ElectrostaticError> {
        let FloatingCapacitanceResponse {
            component_relative_voltage_v,
            component_reference_voltage_v,
            terminal_charge_c,
            maximum_floating_charge_residual_c,
            terminal_energy_residual_j,
            field,
            ..
        } = self.core.solve(
            retained_voltage_v,
            floating_charge_c,
            charge_tolerance_c,
            energy_tolerance_j,
            cancelled,
        )?;
        Ok(AxisymmetricFloatingCapacitanceResponse {
            component_relative_voltage_v,
            component_reference_voltage_v,
            terminal_charge_c,
            maximum_floating_charge_residual_c,
            terminal_energy_residual_j,
            field: AxisymmetricElectricField {
                problem: self.source.problem,
                field,
            },
        })
    }
}

/// Accepted mixed-control state in complete original terminal order.
pub struct AxisymmetricFloatingCapacitanceResponse<'m> {
    /// Voltages relative to the first retained terminal of each electrical
    /// component. Keep these values separate from their common reference when
    /// supplying the native circuit or recovering small charge-induced fields.
    pub component_relative_voltage_v: Vec<f64>,
    /// Prescribed common reference for each original terminal [V].
    pub component_reference_voltage_v: Vec<f64>,
    /// Independent original weak conductor reactions for the angular span [C].
    pub terminal_charge_c: Vec<f64>,
    pub maximum_floating_charge_residual_c: f64,
    /// Field storage minus half of the original voltage-charge pairing [J].
    pub terminal_energy_residual_j: f64,
    field: AxisymmetricElectricField<'m>,
}

impl<'m> AxisymmetricFloatingCapacitanceResponse<'m> {
    /// The same checked radial field, moved into this typed result without a
    /// new solve, a new recovery, or a clone of its physical arrays.
    pub fn field(&self) -> &AxisymmetricElectricField<'m> {
        &self.field
    }
}

/// Native circuit export retaining the original radial field source.
pub struct AxisymmetricCoupledCapacitanceCircuit<'s, 'm> {
    source: &'s AxisymmetricCapacitanceMatrix<'m>,
    core: CoupledCapacitanceCircuit<'s, 'm>,
}

impl<'s, 'm> AxisymmetricCoupledCapacitanceCircuit<'s, 'm> {
    pub fn source(&self) -> &'s AxisymmetricCapacitanceMatrix<'m> {
        self.source
    }
    pub fn reference_terminal(&self) -> usize {
        self.core.reference_terminal()
    }
    pub fn coordinate_terminals(&self) -> &[usize] {
        self.core.coordinate_terminals()
    }
    /// Append in this order, together with the correspondingly indexed block.
    pub fn branches(&self) -> &[(usize, usize, Branch)] {
        self.core.branches()
    }
    pub fn capacitance_block(&self) -> &CapacitanceBlock {
        self.core.capacitance_block()
    }
    pub fn block_at_branch_offset(
        &self,
        first_branch: usize,
    ) -> Result<CapacitanceBlock, ElectrostaticError> {
        self.core.block_at_branch_offset(first_branch)
    }
    /// Total electrode charges in the returned branch order [C].
    pub fn initial_charges_c(&self) -> &[f64] {
        self.core.initial_charges_c()
    }
    pub fn field_response(&self) -> &CapacitanceResponse {
        self.core.field_response()
    }
    pub fn maximum_charge_residual_c(&self) -> f64 {
        self.core.maximum_charge_residual_c()
    }
    pub fn energy_residual_j(&self) -> f64 {
        self.core.energy_residual_j()
    }
}
