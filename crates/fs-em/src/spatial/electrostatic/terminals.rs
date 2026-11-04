//! Equipotential face electrodes and native circuit storage from the 3D field.
//!
//! A single unit-voltage solve determines the fixed two-terminal capacitance.
//! Subsequent voltage/charge drives reuse its relative potential, then recover
//! and check the actual-amplitude field, charge and energy with original laws.

use super::*;
use fs_phs::circuit::Branch;
use std::collections::BTreeSet;

/// Original geometry/material field and charge/energy-checked capacitance.
/// This object is the dielectric store; a derived native capacitor is its
/// circuit representation, not a second storage contribution to be added.
pub struct SpatialCapacitance<'a> {
    field: SpatialElectricField<'a>,
    faces: [Vec<usize>; 2],
    vertices: [Vec<usize>; 2],
    capacitance_f: f64,
}

pub(super) fn voltage_boundary(
    faces: &[Vec<usize>],
    voltages: &[f64],
) -> SpatialElectrostaticBoundary {
    SpatialElectrostaticBoundary {
        voltage: faces
            .iter()
            .zip(voltages)
            .flat_map(|(faces, voltage)| {
                faces.iter().map(move |&face| ElectrostaticVoltageFace {
                    face,
                    voltage_v: [*voltage; 3],
                })
            })
            .collect(),
        displacement: Vec::new(),
    }
}

pub(super) fn terminal_charge(
    field: &SpatialElectricField<'_>,
    vertices: &[usize],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let (mut charge, mut correction) = (0.0, 0.0);
    for &v in vertices {
        check(cancelled)?;
        let value = field.nodal_boundary_charge_c()[v]
            .ok_or_else(|| invalid("electrode lost its physical charge reaction"))?;
        add(&mut charge, &mut correction, value)?;
    }
    Ok(charge)
}

/// Common exact exterior-face admission for two- and multi-electrode stores.
pub(super) fn electrode_faces(
    problem: &SpatialElectrostatics,
    terminals: &[&[usize]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<Vec<usize>>, Vec<Vec<usize>>)> {
    let mut used = BTreeSet::new();
    let mut assigned = vec![None; problem.mesh().vertex_count];
    let mut vertices = Vec::with_capacity(terminals.len());
    let mut canonical_faces = Vec::with_capacity(terminals.len());
    for (terminal, faces) in terminals.iter().enumerate() {
        if faces.is_empty() || faces.len() > problem.geometry.faces.len() {
            return Err(invalid(
                "nonempty bounded exterior face electrodes required",
            ));
        }
        let mut selected = BTreeSet::new();
        let mut ordered_faces = BTreeSet::new();
        for &face in *faces {
            check(cancelled)?;
            if problem
                .geometry
                .faces
                .get(face)
                .is_none_or(|f| f.cells.len() != 1)
                || !used.insert(face)
            {
                return Err(invalid(
                    "electrodes require distinct original exterior faces",
                ));
            }
            ordered_faces.insert(face);
            for &v in &problem.mesh().faces[face] {
                let v = v as usize;
                if assigned[v].is_some_and(|other| other != terminal) {
                    return Err(invalid("different electrodes cannot share a vertex"));
                }
                assigned[v] = Some(terminal);
                selected.insert(v);
            }
        }
        vertices.push(selected.into_iter().collect());
        canonical_faces.push(ordered_faces.into_iter().collect());
    }
    Ok((canonical_faces, vertices))
}

impl SpatialElectrostatics {
    /// Extract capacitance between disjoint unions of complete exterior faces.
    /// Positive faces are at one volt, negative at zero; all other faces have
    /// zero normal displacement. The resolved domain includes all modeled
    /// fringing/end fields; omitted exterior space is not supplied analytically.
    pub fn capacitance(
        &self,
        positive_faces: &[usize],
        negative_faces: &[usize],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialCapacitance<'_>> {
        check(cancelled)?;
        let (_, vertices) = electrode_faces(self, &[positive_faces, negative_faces], cancelled)?;
        let faces = [positive_faces.to_vec(), negative_faces.to_vec()];
        let vertices: [Vec<usize>; 2] = vertices
            .try_into()
            .map_err(|_| invalid("two capacitance electrodes required"))?;
        let boundary = voltage_boundary(&faces, &[1.0, 0.0]);
        let field = self.solve(
            &vec![0.0; self.materials.len()],
            &boundary,
            controls,
            cancelled,
        )?;
        let capacitance_f = finite(2.0 * field.energy_j())?;
        let charge = terminal_charge(&field, &vertices[0], cancelled)?;
        if !positive(capacitance_f)
            || !positive(charge)
            || (charge - capacitance_f).abs() > controls.charge_tolerance_c
        {
            return Err(invalid(
                "positive field capacitance and independent electrode charge must agree",
            ));
        }
        check(cancelled)?;
        Ok(SpatialCapacitance {
            field,
            faces,
            vertices,
            capacitance_f,
        })
    }
}

impl<'a> SpatialCapacitance<'a> {
    pub fn unit_voltage_field(&self) -> &SpatialElectricField<'a> {
        &self.field
    }
    pub fn capacitance_f(&self) -> f64 {
        self.capacitance_f
    }
    pub fn terminal_faces(&self) -> &[Vec<usize>; 2] {
        &self.faces
    }
    /// Native Kirchhoff/pHS capacitor using the actual three-dimensional field.
    /// Geometry, material and external state are frozen. A later parameter
    /// change requires a new physical work account, not replacement of C alone.
    pub fn circuit_branch(&self) -> Branch {
        Branch::Capacitor {
            farads: self.capacitance_f,
        }
    }

    /// Ground the negative electrode and set the signed potential difference.
    /// No new linear solve or constitutive substitution is performed.
    pub fn field_at_voltage(
        &self,
        voltage_difference_v: f64,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialElectricField<'a>> {
        self.field_at_voltages([voltage_difference_v, 0.0], controls, cancelled)
    }

    /// Recover the actual stored field from two absolute electrode potentials.
    /// Component-relative coordinates keep a common bias out of E, D and U.
    pub fn field_at_voltages(
        &self,
        voltages_v: [f64; 2],
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialElectricField<'a>> {
        check(cancelled)?;
        for value in voltages_v {
            finite(value)?;
        }
        let difference = finite(voltages_v[0] - voltages_v[1])?;
        let boundary = voltage_boundary(&self.faces, &voltages_v);
        let problem = self.field.problem;
        let rho = vec![0.0; problem.materials.len()];
        let (sources, traces) = problem.prepare_boundary(&rho, &boundary, controls, cancelled)?;
        let relative = self
            .field
            .relative_voltage_v
            .iter()
            .map(|value| finite(*value * difference))
            .collect::<Result<Vec<_>>>()?;
        let field = problem.recover_field(
            relative, sources, traces, &boundary, controls, None, cancelled,
        )?;
        let expected_charge = finite(self.capacitance_f * difference)?;
        let expected_energy = finite(0.5 * expected_charge * difference)?;
        let actual_charge = terminal_charge(&field, &self.vertices[0], cancelled)?;
        if (actual_charge - expected_charge).abs() > controls.charge_tolerance_c
            || (field.energy_j() - expected_energy).abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "actual field charge/storage disagrees with original terminal capacitance",
            ));
        }
        check(cancelled)?;
        Ok(field)
    }

    /// Recover the spatial field from the native circuit's signed positive
    /// electrode charge. No independent fitted capacitor or charge is invented.
    pub fn field_at_charge(
        &self,
        positive_charge_c: f64,
        controls: ElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialElectricField<'a>> {
        self.field_at_voltage(
            finite(positive_charge_c / self.capacitance_f)?,
            controls,
            cancelled,
        )
    }
}
