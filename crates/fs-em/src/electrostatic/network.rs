//! Source-bound multi-electrode electrostatics and native capacitor networks.
//!
//! Unit-electrode solutions come from the original FEEC/checked-CG owner.
//! Their electric-energy Gram matrix is independently checked against every
//! electrode charge reaction, reciprocity and the common-voltage nullspace.
//! A passive circuit reduction uses -C_ij for each pair; positive off-diagonal
//! entries refuse reduction rather than being clipped or fitted.

use super::{
    ElectricField, ElectrostaticControls, ElectrostaticError, PlanarElectrostatics, finite, poll,
    positive,
};
use fs_phs::circuit::Branch;

/// Independent original field and multi-electrode work/accuracy bounds.
#[derive(Clone, Copy, Debug)]
pub struct CapacitanceMatrixControls {
    pub field: ElectrostaticControls,
    pub maximum_terminals: usize,
    /// Exactly one native unit-voltage solve is needed per terminal.
    pub maximum_field_solves: usize,
    /// Total retained scalar field/charge samples across all unit solutions,
    /// counted as terminals*(2*vertices + 6*triangles). Matrix entries must fit
    /// this bound too. This bounds the repeated solves and dense Gram work.
    pub maximum_field_values: usize,
    /// Absolute [F] charge/energy coefficient, reciprocity and gauge allowance.
    /// Unit potentials are one volt; this is independent of each field solve.
    pub coefficient_tolerance_f: f64,
}

/// Complete terminal response at one finite voltage vector.
#[derive(Clone, Debug)]
pub struct CapacitanceResponse {
    /// Free electrode charges in the original terminal ordering [C].
    pub terminal_charge_c: Vec<f64>,
    /// Integral E.epsilon.E/2 of the superposed field on original cells [J].
    pub energy_j: f64,
}

/// All terminals, including the eventual circuit reference, remain represented.
/// Matrix entry (i,j) is charge on i for one volt on j and zero on the others.
/// Unassigned exterior edges remain the original natural insulated boundary.
pub struct PlanarCapacitanceMatrix<'m> {
    problem: &'m PlanarElectrostatics,
    terminals: Vec<Vec<usize>>,
    basis: Vec<ElectricField<'m>>,
    matrix_f: Vec<f64>,
    charge_matrix_f: Vec<f64>,
    maximum_coefficient_residual_f: f64,
    maximum_gauge_residual_f: f64,
}

/// A source-retaining reduction ready for the existing Kirchhoff circuit owner.
pub struct CapacitanceCircuit<'s, 'm> {
    source: &'s PlanarCapacitanceMatrix<'m>,
    branches: Vec<(usize, usize, Branch)>,
    initial_charges_c: Vec<f64>,
    response: CapacitanceResponse,
    maximum_charge_residual_c: f64,
    energy_residual_j: f64,
}

fn invalid(reason: &'static str) -> ElectrostaticError {
    ElectrostaticError::Invalid(reason)
}
fn add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<(), ElectrostaticError> {
    let corrected = finite(value - *correction)?;
    let next = finite(*sum + corrected)?;
    *correction = finite((next - *sum) - corrected)?;
    *sum = next;
    Ok(())
}

impl PlanarElectrostatics {
    /// Extract the full Maxwell capacitance matrix from the SAME material field.
    /// Every terminal is a disjoint nonempty union of complete exterior edges;
    /// all provided terminals participate in every solve, including a grounded
    /// reference if desired. No unlisted conductor is inferred or grounded.
    ///
    /// Each connected mesh component must contain a terminal trace. Remaining
    /// exterior edges are insulated. The matrix retains its common-voltage
    /// gauge; choose a reference only when reading a grounded submatrix or
    /// attaching a circuit. Material, geometry and source state remain frozen.
    pub fn capacitance_matrix(
        &self,
        terminals: &[&[usize]],
        controls: CapacitanceMatrixControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarCapacitanceMatrix<'_>, ElectrostaticError> {
        poll(cancelled)?;
        let count = terminals.len();
        let values_per_field = self
            .mesh
            .vertices()
            .len()
            .checked_mul(2)
            .and_then(|n| {
                self.mesh
                    .faces()
                    .len()
                    .checked_mul(6)
                    .and_then(|m| n.checked_add(m))
            })
            .ok_or_else(|| invalid("capacitance basis extent overflow"))?;
        let values = count
            .checked_mul(values_per_field)
            .ok_or_else(|| invalid("capacitance basis extent overflow"))?;
        let entries = count
            .checked_mul(count)
            .ok_or_else(|| invalid("capacitance matrix extent overflow"))?;
        if count < 2
            || count > controls.maximum_terminals
            || count > controls.maximum_field_solves
            || values > controls.maximum_field_values
            || entries > controls.maximum_field_values
            || !positive(controls.coefficient_tolerance_f)
        {
            return Err(invalid(
                "bounded electrode count, basis solves/storage and positive coefficient tolerance required",
            ));
        }
        let mut boundary = self.terminal_boundary(terminals, &vec![0.0; count], cancelled)?;
        let zero_charge = vec![0.0; self.mesh.faces().len()];
        let mut basis = Vec::with_capacity(count);
        let mut charge_matrix = vec![0.0; entries];
        for column in 0..count {
            poll(cancelled)?;
            let mut offset = 0;
            for (terminal, vertices) in terminals.iter().enumerate() {
                for entry in &mut boundary[offset..offset + vertices.len()] {
                    entry.1 = if terminal == column { 1.0 } else { 0.0 };
                }
                offset += vertices.len();
            }
            let field = self.solve(&zero_charge, &boundary, controls.field, cancelled)?;
            for (row, vertices) in terminals.iter().enumerate() {
                let (mut charge, mut correction) = (0.0, 0.0);
                for &vertex in vertices.iter() {
                    poll(cancelled)?;
                    let reaction = field.nodal_boundary_charge_c()[vertex].ok_or_else(|| {
                        invalid("electrode lost its original prescribed-node charge")
                    })?;
                    add(&mut charge, &mut correction, reaction)?;
                }
                charge_matrix[row * count + column] = charge;
            }
            basis.push(field);
        }
        let mut matrix_f = vec![0.0; entries];
        for row in 0..count {
            for column in row..count {
                let (mut coefficient, mut correction) = (0.0, 0.0);
                for cell in 0..self.mesh.faces().len() {
                    poll(cancelled)?;
                    let e = basis[row].electric_field_v_m()[cell];
                    let d = basis[column].displacement_c_m2()[cell];
                    let volume = self
                        .mesh
                        .face_measure(cell)
                        .ok_or_else(|| invalid("missing original dielectric volume"))?;
                    let pairing = finite(volume * e[0].mul_add(d[0], e[1] * d[1]))?;
                    add(&mut coefficient, &mut correction, pairing)?;
                }
                matrix_f[row * count + column] = coefficient;
                matrix_f[column * count + row] = coefficient;
            }
        }
        let mut maximum_coefficient_residual_f: f64 = 0.0;
        let mut maximum_gauge_residual_f: f64 = 0.0;
        for row in 0..count {
            let (mut gauge, mut correction) = (0.0, 0.0);
            for column in 0..count {
                poll(cancelled)?;
                let index = row * count + column;
                let coefficient = matrix_f[index];
                maximum_coefficient_residual_f = maximum_coefficient_residual_f
                    .max(finite(coefficient - charge_matrix[index])?.abs())
                    .max(finite(charge_matrix[index] - charge_matrix[column * count + row])?.abs());
                add(&mut gauge, &mut correction, coefficient)?;
            }
            maximum_gauge_residual_f = maximum_gauge_residual_f.max(gauge.abs());
            maximum_coefficient_residual_f = maximum_coefficient_residual_f
                .max(finite(matrix_f[row * count + row] - 2.0 * basis[row].energy_j())?.abs());
            if matrix_f[row * count + row] < 0.0 {
                return Err(invalid("negative electrode self energy"));
            }
        }
        if maximum_coefficient_residual_f > controls.coefficient_tolerance_f
            || maximum_gauge_residual_f > controls.coefficient_tolerance_f
        {
            return Err(invalid(
                "electrode charge, field energy, reciprocity or common-voltage gauge failed",
            ));
        }
        poll(cancelled)?;
        Ok(PlanarCapacitanceMatrix {
            problem: self,
            terminals: terminals.iter().map(|v| v.to_vec()).collect(),
            basis,
            matrix_f,
            charge_matrix_f: charge_matrix,
            maximum_coefficient_residual_f,
            maximum_gauge_residual_f,
        })
    }
}

impl<'m> PlanarCapacitanceMatrix<'m> {
    pub fn problem(&self) -> &'m PlanarElectrostatics {
        self.problem
    }
    pub fn terminals(&self) -> &[Vec<usize>] {
        &self.terminals
    }
    /// Original ordering, row-major square matrix, in farads.
    pub fn maxwell_matrix_f(&self) -> &[f64] {
        &self.matrix_f
    }
    /// Actual original unit-voltage fields, retaining material and charge data.
    pub fn unit_fields(&self) -> &[ElectricField<'m>] {
        &self.basis
    }
    pub fn maximum_coefficient_residual_f(&self) -> f64 {
        self.maximum_coefficient_residual_f
    }
    pub fn maximum_gauge_residual_f(&self) -> f64 {
        self.maximum_gauge_residual_f
    }

    /// Remove one explicit reference row/column, preserving all other terminal
    /// ordering. This does not eliminate a floating electrode or change its
    /// physical charge constraint. No matrix inverse or new solve is performed.
    pub fn grounded_matrix_f(&self, reference: usize) -> Result<Vec<f64>, ElectrostaticError> {
        let count = self.terminals.len();
        if reference >= count {
            return Err(invalid(
                "reference electrode is outside the original terminal ordering",
            ));
        }
        let mut reduced = Vec::with_capacity((count - 1) * (count - 1));
        for row in 0..count {
            if row == reference {
                continue;
            }
            for column in 0..count {
                if column != reference {
                    reduced.push(self.matrix_f[row * count + column]);
                }
            }
        }
        Ok(reduced)
    }

    fn voltage_differences(&self, voltage_v: &[f64]) -> Result<Vec<f64>, ElectrostaticError> {
        if voltage_v.len() != self.terminals.len() || voltage_v.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "one finite voltage per original electrode required",
            ));
        }
        voltage_v.iter().map(|v| finite(v - voltage_v[0])).collect()
    }

    /// Recover an actual spatial electric field from the original unit-voltage
    /// solutions, including native weak charge reactions and energy admission.
    /// This performs no new linear solve. Nodal potentials are reported relative
    /// to the FIRST electrode; a common voltage never enters field arithmetic.
    ///
    /// Supply instantaneous terminal voltages from the native circuit in this
    /// owner's original electrode order. The returned field can drive the same
    /// `ElectricField::maxwell_wrench` mechanical connection as a direct solve.
    /// Actual-voltage charge and energy tolerances are independent of the unit
    /// basis extraction: scaling an approximate basis cannot bypass acceptance.
    /// The field retains the frozen material, geometry and zero volume sources;
    /// motion or changing dielectric state still needs a coupled work owner.
    pub fn field_at_voltages(
        &self,
        voltage_v: &[f64],
        charge_tolerance_c: f64,
        energy_tolerance_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ElectricField<'m>, ElectrostaticError> {
        poll(cancelled)?;
        if !positive(charge_tolerance_c) || !positive(energy_tolerance_j) {
            return Err(invalid(
                "positive actual-field charge and energy tolerances required",
            ));
        }
        let differences = self.voltage_differences(voltage_v)?;
        let vertices = self.problem.mesh.vertices().len();
        let mut voltage = vec![0.0; vertices];
        for (v, potential) in voltage.iter_mut().enumerate() {
            poll(cancelled)?;
            for (field, &scale) in self.basis.iter().zip(&differences) {
                *potential = finite(scale.mul_add(field.voltage_v()[v], *potential))?;
            }
        }
        let mut prescribed = vec![None; vertices];
        for (terminal, &potential) in self.terminals.iter().zip(&differences) {
            poll(cancelled)?;
            for &v in terminal {
                prescribed[v] = Some(potential);
            }
        }
        self.problem.recover_field(
            voltage,
            &vec![0.0; self.problem.mesh.faces().len()],
            &vec![0.0; vertices],
            &prescribed,
            0.0,
            charge_tolerance_c,
            energy_tolerance_j,
            None,
            cancelled,
        )
    }

    /// Superpose actual unit fields and return their charge and physical energy.
    /// Voltage differences are formed before products; an arbitrary common
    /// voltage cannot manufacture electric field or stored energy.
    pub fn response(
        &self,
        voltage_v: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CapacitanceResponse, ElectrostaticError> {
        poll(cancelled)?;
        let count = self.terminals.len();
        let differences = self.voltage_differences(voltage_v)?;
        let mut terminal_charge_c = Vec::with_capacity(count);
        for row in 0..count {
            let (mut charge, mut correction) = (0.0, 0.0);
            for (column, &voltage) in differences.iter().enumerate() {
                poll(cancelled)?;
                // Retain the independent original weak charge reactions. A
                // circuit gate at large voltages must also see basis errors,
                // rather than compare the reduced matrix against itself.
                add(
                    &mut charge,
                    &mut correction,
                    finite(self.charge_matrix_f[row * count + column] * voltage)?,
                )?;
            }
            terminal_charge_c.push(charge);
        }
        let (mut energy_j, mut correction) = (0.0, 0.0);
        for cell in 0..self.problem.mesh.faces().len() {
            poll(cancelled)?;
            let mut e = [0.0; 2];
            for (field, &voltage) in self.basis.iter().zip(&differences) {
                for axis in 0..2 {
                    e[axis] =
                        finite(voltage.mul_add(field.electric_field_v_m()[cell][axis], e[axis]))?;
                }
            }
            let volume = self
                .problem
                .mesh
                .face_measure(cell)
                .ok_or_else(|| invalid("missing original dielectric volume"))?;
            let d = self.problem.permittivity[cell].displacement_c_m2(e);
            let d = [finite(d[0])?, finite(d[1])?];
            let energy = finite(0.5 * volume * e[0].mul_add(d[0], e[1] * d[1]))?;
            if energy < 0.0 || e != [0.0; 2] && energy == 0.0 {
                return Err(invalid("superposed electric energy is not representable"));
            }
            add(&mut energy_j, &mut correction, energy)?;
        }
        poll(cancelled)?;
        Ok(CapacitanceResponse {
            terminal_charge_c,
            energy_j,
        })
    }

    /// Lower to original Kirchhoff capacitor branches and consistent branch
    /// charges. Branch orientation is (terminal i, terminal j), i<j, and its
    /// initial charge is C_ij_pair*(V_i-V_j), matching the native circuit owner.
    /// The complete explicit node map must be one-to-one; combine electrically
    /// shorted traces into one electrode before extraction.
    ///
    /// Every nonzero mutual coefficient must be strictly negative. A positive
    /// off-diagonal from an unsuitable mesh/solve is NOT clipped. The Maxwell
    /// matrix remains available independently, but a passive capacitor network
    /// then refuses. Zero mutual coefficients create no branch. Every terminal
    /// must participate in at least one physical capacitor.
    ///
    /// Actual network charge and storage must match the retained field at these
    /// initial voltages within independent absolute C/J tolerances. The native
    /// circuit owns subsequent charge/current evolution and work accounting.
    /// Do not add this network beside another copy of the same dielectric store.
    pub fn circuit_at_voltages<'s>(
        &'s self,
        terminal_nodes: &[usize],
        voltage_v: &[f64],
        charge_tolerance_c: f64,
        energy_tolerance_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CapacitanceCircuit<'s, 'm>, ElectrostaticError> {
        poll(cancelled)?;
        let count = self.terminals.len();
        if terminal_nodes.len() != count
            || !positive(charge_tolerance_c)
            || !positive(energy_tolerance_j)
            || terminal_nodes
                .iter()
                .enumerate()
                .any(|(i, &node)| node == usize::MAX || terminal_nodes[..i].contains(&node))
        {
            return Err(invalid(
                "complete unique circuit-node map and positive charge/energy tolerances required",
            ));
        }
        let response = self.response(voltage_v, cancelled)?;
        let mut branches = Vec::new();
        let mut initial_charges_c = Vec::new();
        let mut charges = vec![0.0; count];
        let mut incident = vec![false; count];
        let (mut energy, mut correction) = (0.0, 0.0);
        for i in 0..count {
            for j in i + 1..count {
                poll(cancelled)?;
                let mutual = self.matrix_f[i * count + j];
                if mutual > 0.0 {
                    return Err(invalid(
                        "positive mutual capacitance coefficient cannot form a passive circuit; refine the field discretization",
                    ));
                }
                if mutual == 0.0 {
                    continue;
                }
                let capacitance = -mutual;
                if !positive(capacitance.recip()) {
                    return Err(invalid(
                        "pair capacitance is not representable by the native storage owner",
                    ));
                }
                let difference = finite(voltage_v[i] - voltage_v[j])?;
                let charge = finite(capacitance * difference)?;
                let stored = finite(0.5 * charge * difference)?;
                if difference != 0.0 && (charge == 0.0 || stored == 0.0) {
                    return Err(invalid("initial pair charge or energy disappeared"));
                }
                charges[i] = finite(charges[i] + charge)?;
                charges[j] = finite(charges[j] - charge)?;
                add(&mut energy, &mut correction, stored)?;
                incident[i] = true;
                incident[j] = true;
                branches.push((
                    terminal_nodes[i],
                    terminal_nodes[j],
                    Branch::Capacitor {
                        farads: capacitance,
                    },
                ));
                initial_charges_c.push(charge);
            }
        }
        if incident.iter().any(|&connected| !connected) {
            return Err(invalid(
                "an electrode has no capacitive connection to the other terminals",
            ));
        }
        let mut maximum_charge_residual_c: f64 = 0.0;
        for (&network, &field) in charges.iter().zip(&response.terminal_charge_c) {
            maximum_charge_residual_c =
                maximum_charge_residual_c.max(finite(network - field)?.abs());
        }
        let energy_residual_j = finite(energy - response.energy_j)?;
        if maximum_charge_residual_c > charge_tolerance_c
            || energy_residual_j.abs() > energy_tolerance_j
        {
            return Err(invalid(
                "passive circuit charge or stored energy disagrees with the original dielectric field",
            ));
        }
        poll(cancelled)?;
        Ok(CapacitanceCircuit {
            source: self,
            branches,
            initial_charges_c,
            response,
            maximum_charge_residual_c,
            energy_residual_j,
        })
    }
}
impl<'s, 'm> CapacitanceCircuit<'s, 'm> {
    pub fn source(&self) -> &'s PlanarCapacitanceMatrix<'m> {
        self.source
    }
    /// Native graph branch tuples; append these together in this exact order.
    pub fn branches(&self) -> &[(usize, usize, Branch)] {
        &self.branches
    }
    /// Native capacitor charge order for the returned branches only. Place each
    /// value at its corresponding CircuitDae charge_index after graph assembly.
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
