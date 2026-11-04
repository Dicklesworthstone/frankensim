//! Reusable multi-electrode response of the original tetrahedral dielectric.
//!
//! Q=C V and U=V.C.V/2 describe this source-free, fixed-geometry linear store.
//! Each electrode is an equipotential union of original exterior faces. A
//! shared electrode electrically joins otherwise disconnected dielectric
//! volumes. Separate electrical components retain separate voltage gauges.
//! Only native solved potentials are retained; every actual drive recovers
//! and checks the original cell field, weak reactions and dielectric energy.

use super::terminals::{electrode_faces, terminal_charge, voltage_boundary};
use super::*;
use std::collections::BTreeMap;

/// Explicit memory and numerical-work bounds for terminal response extraction.
#[derive(Clone, Copy, Debug)]
pub struct SpatialCapacitanceControls {
    pub maximum_terminals: usize,
    /// One native solve per nonreference terminal, N minus electrical components.
    pub maximum_field_solves: usize,
    /// Conservative scalar storage: (N+2)*vertices + 3*N*N. Native field
    /// recovery and the already admitted volume operator retain their own bounds.
    pub maximum_response_values: usize,
    /// Conservative Gram/projection work: N*N*(vertices+N+2).
    /// Sparse operator applications are additionally bounded by the native mesh.
    pub maximum_dense_work: usize,
    /// Positive numerical probe amplitude. It is not an extra physical source.
    pub reference_voltage_v: f64,
    /// Absolute [F] allowance between energy and native reaction coefficients.
    pub coefficient_tolerance_f: f64,
}

/// Frozen source-free dielectric and its full reciprocal Maxwell matrix.
/// Reprepare after geometry/material changes. Bulk free charge, impressed
/// polarization and charged sheets require the source-aware direct field API;
/// this homogeneous response does not silently discard or absorb such sources.
pub struct SpatialCapacitanceMatrix<'a> {
    problem: &'a SpatialElectrostatics,
    faces: Vec<Vec<usize>>,
    vertices: Vec<Vec<usize>>,
    references: Vec<usize>,
    coordinates: Vec<usize>,
    basis: Vec<Vec<f64>>,
    matrix_f: Vec<f64>,
    maximum_coefficient_residual_f: f64,
}

/// Actual-amplitude native field and independently recovered electrode charges.
/// The field owns the sole dielectric storage; a circuit representation of it
/// must not be counted as an additional physical energy contribution.
pub struct SpatialCapacitanceResponse<'a> {
    pub terminal_charge_c: Vec<f64>,
    pub maximum_charge_residual_c: f64,
    /// Native field energy minus one half of all electrode voltage-charge work.
    pub terminal_energy_residual_j: f64,
    field: SpatialElectricField<'a>,
}

impl<'a> SpatialCapacitanceResponse<'a> {
    pub fn field(&self) -> &SpatialElectricField<'a> {
        &self.field
    }
    pub fn into_field(self) -> SpatialElectricField<'a> {
        self.field
    }
}

impl SpatialElectrostatics {
    /// Extract the full multi-electrode Maxwell capacitance matrix on the same
    /// scalar/tensor dielectric operator. Unlisted exterior faces have D.n=0.
    /// All modeled end/fringing fields are retained; omitted exterior space is
    /// not supplied analytically. Electrode faces must be complete and disjoint,
    /// and distinct electrodes must not share a vertex.
    ///
    /// The first terminal of each electrical component is its numerical voltage
    /// reference. Only the remaining basis fields are solved. Energy Gram
    /// coefficients independently agree with original weak electrode reactions;
    /// no assumed pair-capacitor graph or fitted field is introduced.
    pub fn capacitance_matrix(
        &self,
        terminals: &[&[usize]],
        controls: SpatialCapacitanceControls,
        field_controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialCapacitanceMatrix<'_>> {
        check(cancelled)?;
        let n = terminals.len();
        let nv = self.mesh().vertex_count;
        let square = n.checked_mul(n);
        let values = n
            .checked_add(2)
            .and_then(|count| count.checked_mul(nv)?.checked_add(square?.checked_mul(3)?));
        let work = square.and_then(|square| square.checked_mul(nv.checked_add(n)?.checked_add(2)?));
        if n < 2
            || n > controls.maximum_terminals
            || values.is_none_or(|value| value > controls.maximum_response_values)
            || work.is_none_or(|value| value > controls.maximum_dense_work)
            || !positive(controls.reference_voltage_v)
            || !positive(controls.coefficient_tolerance_f)
        {
            return Err(invalid(
                "bounded terminal storage/work and positive probe/tolerance required",
            ));
        }
        let (faces, vertices) = electrode_faces(self, terminals, cancelled)?;
        let references = terminal_references(self, &vertices, cancelled)?;
        let coordinates: Vec<_> = (0..n).filter(|&t| references[t] != t).collect();
        let m = coordinates.len();
        if m > controls.maximum_field_solves {
            return Err(invalid(
                "terminal response exceeds the native field solve allowance",
            ));
        }
        let rho = vec![0.0; self.materials.len()];
        let mut basis = Vec::with_capacity(m);
        let mut reactions = vec![0.0; n * m];
        for (column, &terminal) in coordinates.iter().enumerate() {
            check(cancelled)?;
            let mut voltage = vec![0.0; n];
            voltage[terminal] = controls.reference_voltage_v;
            let boundary = voltage_boundary(&faces, &voltage);
            let field = self.solve(&rho, &boundary, field_controls, cancelled)?;
            for (row, vertices) in vertices.iter().enumerate() {
                reactions[row * m + column] = finite(
                    terminal_charge(&field, vertices, cancelled)? / controls.reference_voltage_v,
                )?;
            }
            basis.push(
                field
                    .relative_voltage_v()
                    .iter()
                    .map(|&value| finite(value / controls.reference_voltage_v))
                    .collect::<Result<Vec<_>>>()?,
            );
        }

        // The original stiffness contains the complete inhomogeneous tensor
        // Hodge. Its Gram form keeps every cross term and allows mutual signs
        // that a passive pair-capacitor graph cannot represent.
        let mut matrix_f = vec![0.0; n * n];
        for column in 0..m {
            let mut applied = vec![0.0; nv];
            for (v, out) in applied.iter_mut().enumerate() {
                check(cancelled)?;
                let (indices, values) = self.stiffness.row(v);
                let mut correction = 0.0;
                for (&w, &value) in indices.iter().zip(values) {
                    add(
                        out,
                        &mut correction,
                        finite(value * finite(basis[column][w] - basis[column][v])?)?,
                    )?;
                }
            }
            for row in 0..=column {
                let (mut value, mut correction) = (0.0, 0.0);
                for (&potential, &charge) in basis[row].iter().zip(&applied) {
                    check(cancelled)?;
                    add(&mut value, &mut correction, finite(potential * charge)?)?;
                }
                let (a, b) = (coordinates[row], coordinates[column]);
                matrix_f[a * n + b] = value;
                matrix_f[b * n + a] = value;
            }
        }
        let reference_terminals: Vec<_> = (0..n).filter(|&t| references[t] == t).collect();
        for &row in &coordinates {
            for &reference in &reference_terminals {
                check(cancelled)?;
                let (mut value, mut correction) = (0.0, 0.0);
                for &column in &coordinates {
                    if references[column] == reference {
                        add(&mut value, &mut correction, -matrix_f[row * n + column])?;
                    }
                }
                matrix_f[row * n + reference] = value;
                matrix_f[reference * n + row] = value;
            }
        }
        for (i, &row) in reference_terminals.iter().enumerate() {
            for &column in &reference_terminals[..=i] {
                check(cancelled)?;
                let (mut value, mut correction) = (0.0, 0.0);
                for &coordinate in &coordinates {
                    if references[coordinate] == row {
                        add(
                            &mut value,
                            &mut correction,
                            -matrix_f[coordinate * n + column],
                        )?;
                    }
                }
                matrix_f[row * n + column] = value;
                matrix_f[column * n + row] = value;
            }
        }
        let mut maximum_coefficient_residual_f = 0.0_f64;
        for row in 0..n {
            check(cancelled)?;
            let (mut gauge, mut correction) = (0.0, 0.0);
            if matrix_f[row * n + row] < 0.0
                || (references[row] != row && !positive(matrix_f[row * n + row]))
            {
                return Err(invalid(
                    "terminal has no representably positive electric storage",
                ));
            }
            for column in 0..n {
                add(&mut gauge, &mut correction, matrix_f[row * n + column])?;
                if references[row] != references[column] {
                    maximum_coefficient_residual_f =
                        maximum_coefficient_residual_f.max(matrix_f[row * n + column].abs());
                }
            }
            maximum_coefficient_residual_f = maximum_coefficient_residual_f.max(gauge.abs());
            for (column, &terminal) in coordinates.iter().enumerate() {
                maximum_coefficient_residual_f = maximum_coefficient_residual_f
                    .max(finite(matrix_f[row * n + terminal] - reactions[row * m + column])?.abs());
            }
        }
        if maximum_coefficient_residual_f > controls.coefficient_tolerance_f {
            return Err(invalid(
                "native terminal reactions, energy Gram matrix or component gauges disagree",
            ));
        }
        check(cancelled)?;
        Ok(SpatialCapacitanceMatrix {
            problem: self,
            faces,
            vertices,
            references,
            coordinates,
            basis,
            matrix_f,
            maximum_coefficient_residual_f,
        })
    }
}

impl<'a> SpatialCapacitanceMatrix<'a> {
    pub fn problem(&self) -> &'a SpatialElectrostatics {
        self.problem
    }
    pub fn terminal_faces(&self) -> &[Vec<usize>] {
        &self.faces
    }
    /// Full row-major Maxwell matrix [F], in original terminal order.
    pub fn maxwell_matrix_f(&self) -> &[f64] {
        &self.matrix_f
    }
    /// Canonical voltage-reference terminal for each original terminal.
    pub fn reference_terminals(&self) -> &[usize] {
        &self.references
    }
    /// Independent terminal voltage coordinates, excluding component references.
    pub fn coordinate_terminals(&self) -> &[usize] {
        &self.coordinates
    }
    pub fn field_solves(&self) -> usize {
        self.basis.len()
    }
    pub fn maximum_coefficient_residual_f(&self) -> f64 {
        self.maximum_coefficient_residual_f
    }

    fn voltage_differences(&self, voltage_v: &[f64]) -> Result<Vec<f64>> {
        if voltage_v.len() != self.faces.len() || voltage_v.iter().any(|v| !v.is_finite()) {
            return Err(invalid("one finite voltage per original terminal required"));
        }
        voltage_v
            .iter()
            .zip(&self.references)
            .map(|(&voltage, &reference)| finite(voltage - voltage_v[reference]))
            .collect()
    }

    /// Recover an actual-amplitude field with no new FE solve. Component-local
    /// voltage differences enter all arithmetic; the original supplied voltages
    /// remain available as observations. Actual-scale charge and energy checks
    /// cannot be bypassed by scaling an approximately extracted basis.
    pub fn response(
        &self,
        voltage_v: &[f64],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialCapacitanceResponse<'a>> {
        check(cancelled)?;
        let differences = self.voltage_differences(voltage_v)?;
        let boundary = voltage_boundary(&self.faces, voltage_v);
        let rho = vec![0.0; self.problem.materials.len()];
        let (sources, traces) = self
            .problem
            .prepare_boundary(&rho, &boundary, controls, cancelled)?;
        let mut relative = vec![0.0; self.problem.mesh().vertex_count];
        for (v, out) in relative.iter_mut().enumerate() {
            check(cancelled)?;
            let mut correction = 0.0;
            for (&coordinate, basis) in self.coordinates.iter().zip(&self.basis) {
                add(
                    out,
                    &mut correction,
                    finite(differences[coordinate] * basis[v])?,
                )?;
            }
        }
        let field = self.problem.recover_field(
            relative, sources, traces, &boundary, controls, None, cancelled,
        )?;
        let n = self.faces.len();
        let mut terminal_charge_c = Vec::with_capacity(n);
        let mut maximum_charge_residual_c = 0.0_f64;
        let (mut work, mut work_correction) = (0.0, 0.0);
        let (mut expected_work, mut expected_correction) = (0.0, 0.0);
        for (row, vertices) in self.vertices.iter().enumerate() {
            check(cancelled)?;
            let charge = terminal_charge(&field, vertices, cancelled)?;
            let (mut expected, mut correction) = (0.0, 0.0);
            for (column, &voltage) in differences.iter().enumerate() {
                add(
                    &mut expected,
                    &mut correction,
                    finite(self.matrix_f[row * n + column] * voltage)?,
                )?;
            }
            maximum_charge_residual_c =
                maximum_charge_residual_c.max(finite(charge - expected)?.abs());
            add(
                &mut work,
                &mut work_correction,
                finite(0.5 * differences[row] * charge)?,
            )?;
            add(
                &mut expected_work,
                &mut expected_correction,
                finite(0.5 * differences[row] * expected)?,
            )?;
            terminal_charge_c.push(charge);
        }
        let terminal_energy_residual_j = finite(field.energy_j() - work)?;
        if maximum_charge_residual_c > controls.charge_tolerance_c
            || terminal_energy_residual_j.abs() > controls.energy_tolerance_j
            || finite(field.energy_j() - expected_work)?.abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "actual native field disagrees with terminal charge or electric storage",
            ));
        }
        check(cancelled)?;
        Ok(SpatialCapacitanceResponse {
            field,
            terminal_charge_c,
            maximum_charge_residual_c,
            terminal_energy_residual_j,
        })
    }

    pub fn field_at_voltages(
        &self,
        voltage_v: &[f64],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialElectricField<'a>> {
        self.response(voltage_v, controls, cancelled)
            .map(SpatialCapacitanceResponse::into_field)
    }
}

// A physical component is a bulk-connected volume or a group of such volumes
// joined by the same equipotential terminal. No tiny numerical matrix entry
// decides whether a floating-charge constraint has an electrical reference.
fn terminal_references(
    problem: &SpatialElectrostatics,
    terminals: &[Vec<usize>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<usize>> {
    fn root(parent: &[usize], mut terminal: usize) -> usize {
        while parent[terminal] != terminal {
            terminal = parent[terminal];
        }
        terminal
    }
    let mut parent: Vec<_> = (0..terminals.len()).collect();
    let mut first = BTreeMap::new();
    for (terminal, vertices) in terminals.iter().enumerate() {
        for &v in vertices {
            check(cancelled)?;
            let component = problem.geometry.components[v];
            if let Some(&previous) = first.get(&component) {
                let a = root(&parent, previous);
                let b = root(&parent, terminal);
                parent[a.max(b)] = a.min(b);
            } else {
                first.insert(component, terminal);
            }
        }
    }
    (0..parent.len())
        .map(|terminal| {
            check(cancelled)?;
            Ok(root(&parent, terminal))
        })
        .collect()
}
