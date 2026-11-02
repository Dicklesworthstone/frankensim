//! Equipotential face electrodes on the original three-dimensional Ohmic owner.
//!
//! One positive-reference solve determines a fixed linear material response.
//! Actual drives reuse relative potentials and recover native currents and
//! bulk/contact heat again, without assuming a uniform injection density.

use super::*;

/// A linear terminal response bound to the original source and geometry.
/// Material/contact refresh requires preparation on the replacement owner.
pub struct SpatialResistance<'a, C = ResolvedConductorState> {
    problem: &'a SpatialConduction<C>,
    field: SpatialConductionField,
    relative: Vec<f64>,
    faces: [Vec<usize>; 2],
    vertices: [Vec<usize>; 2],
    resistance_ohm: f64,
    reference_voltage_v: f64,
}

impl<C: SpatialConductorMaterial> SpatialConduction<C> {
    /// Extract resistance between disjoint unions of complete exterior faces.
    /// Positive faces are held at one volt, negative faces at zero, and other
    /// faces insulate. The pair must share a bulk/contact-connected component.
    pub fn resistance(
        &self,
        positive_faces: &[usize],
        negative_faces: &[usize],
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialResistance<'_, C>> {
        self.resistance_at_reference_voltage(
            positive_faces,
            negative_faces,
            1.0,
            controls,
            cancelled,
        )
    }

    /// Extraction at an explicitly chosen positive reference voltage. This
    /// changes the numerical excitation, never the fixed material/source law.
    pub fn resistance_at_reference_voltage(
        &self,
        positive_faces: &[usize],
        negative_faces: &[usize],
        reference_voltage_v: f64,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialResistance<'_, C>> {
        check(cancelled)?;
        if !positive(reference_voltage_v) {
            return Err(SpatialConductionError::Invalid(
                "positive finite terminal reference voltage required",
            ));
        }
        let (faces, vertices) = self.terminal_faces(positive_faces, negative_faces, cancelled)?;
        let boundary = voltage_boundary(&faces, &[reference_voltage_v, 0.0]);
        let traces = self.prepare_boundary(&boundary, controls, cancelled)?;
        let (relative, report) = solve_scalar_system(
            &self.stiffness,
            &traces.flux_load_wb,
            &traces.prescribed,
            controls.linear,
            cancelled,
        )?;
        let field = self.recover_field(
            relative.clone(),
            traces,
            controls.current_tolerance_a,
            controls.power_tolerance_w,
            report,
            cancelled,
        )?;
        let current = terminal_current(&field, &vertices[0], cancelled)?;
        let conductance = finite(current / reference_voltage_v)?;
        let resistance_ohm = finite(reference_voltage_v / current)?;
        if !positive(current) || !positive(conductance) || !positive(resistance_ohm) {
            return Err(SpatialConductionError::Invalid(
                "representably positive two-terminal conductance required",
            ));
        }
        check(cancelled)?;
        Ok(SpatialResistance {
            problem: self,
            field,
            relative,
            faces,
            vertices,
            resistance_ohm,
            reference_voltage_v,
        })
    }

    /// Signed TOTAL current [A] through ideal equipotential face electrodes.
    /// Positive current enters the positive faces and leaves the grounded
    /// negative faces. Local injection follows the native field. Zero current
    /// still requires a conducting pair, established by the reference solve.
    pub fn solve_current(
        &self,
        positive_faces: &[usize],
        negative_faces: &[usize],
        current_a: f64,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductionField> {
        self.solve_current_with_reference_voltage(
            positive_faces,
            negative_faces,
            current_a,
            1.0,
            controls,
            cancelled,
        )
    }

    /// The same total-current drive using a caller-selected reference voltage.
    pub fn solve_current_with_reference_voltage(
        &self,
        positive_faces: &[usize],
        negative_faces: &[usize],
        current_a: f64,
        reference_voltage_v: f64,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductionField> {
        check(cancelled)?;
        if !current_a.is_finite() {
            return Err(SpatialConductionError::Invalid(
                "finite total current required",
            ));
        }
        self.resistance_at_reference_voltage(
            positive_faces,
            negative_faces,
            reference_voltage_v,
            controls,
            cancelled,
        )?
        .field_at_current(
            current_a,
            controls.current_tolerance_a,
            controls.power_tolerance_w,
            cancelled,
        )
    }

    fn terminal_faces(
        &self,
        positive_faces: &[usize],
        negative_faces: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([Vec<usize>; 2], [Vec<usize>; 2])> {
        let (faces, vertices) =
            self.electrode_faces(&[positive_faces, negative_faces], cancelled)?;
        let faces = faces.try_into().map_err(|_| {
            SpatialConductionError::Invalid("two-terminal admission lost its original face pair")
        })?;
        let vertices = vertices.try_into().map_err(|_| {
            SpatialConductionError::Invalid("two-terminal admission lost its original vertex pair")
        })?;
        Ok((faces, vertices))
    }

    // Shared complete-face admission for two-terminal resistance and the
    // multi-electrode response. A common electrode may join separate bulk
    // components, while finite material contacts already join self.components.
    pub(super) fn electrode_faces(
        &self,
        terminals: &[&[usize]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<Vec<usize>>, Vec<Vec<usize>>)> {
        check(cancelled)?;
        if terminals.len() < 2
            || terminals.iter().any(|faces| faces.is_empty())
            || terminals
                .iter()
                .try_fold(0usize, |n, faces| n.checked_add(faces.len()))
                .is_none_or(|n| n > self.geometry.faces.len())
        {
            return Err(SpatialConductionError::Invalid(
                "at least two nonempty bounded terminal face sets required",
            ));
        }
        let mut faces = Vec::with_capacity(terminals.len());
        let mut vertices = Vec::with_capacity(terminals.len());
        let mut vertex_terminal = vec![None; self.mesh().vertex_count];
        let mut component_terminal: Vec<Option<usize>> = vec![None; self.mesh().vertex_count];
        let mut parent: Vec<_> = (0..terminals.len()).collect();
        fn root(parent: &[usize], mut terminal: usize) -> usize {
            while parent[terminal] != terminal {
                terminal = parent[terminal];
            }
            terminal
        }
        for (side, input) in terminals.iter().enumerate() {
            let mut side_faces = BTreeSet::new();
            let mut side_vertices = BTreeSet::new();
            for &face in *input {
                check(cancelled)?;
                if self
                    .geometry
                    .faces
                    .get(face)
                    .is_none_or(|f| f.cells.len() != 1)
                    || self.contact_faces.contains(&face)
                    || !side_faces.insert(face)
                {
                    return Err(SpatialConductionError::Invalid(
                        "terminal needs unique exterior faces, excluding internal contacts",
                    ));
                }
                for vertex in self.mesh().faces[face] {
                    let vertex = vertex as usize;
                    if vertex_terminal[vertex].is_some_and(|other| other != side) {
                        return Err(SpatialConductionError::Invalid(
                            "different terminal faces must not share vertices",
                        ));
                    }
                    vertex_terminal[vertex] = Some(side);
                    side_vertices.insert(vertex);
                    let component = self.components[vertex];
                    if let Some(other) = component_terminal[component] {
                        let a = root(&parent, side);
                        let b = root(&parent, other);
                        parent[a.max(b)] = a.min(b);
                    } else {
                        component_terminal[component] = Some(side);
                    }
                }
            }
            faces.push(side_faces.into_iter().collect());
            vertices.push(side_vertices.into_iter().collect());
        }
        if (1..terminals.len()).any(|t| root(&parent, t) != root(&parent, 0)) {
            return Err(SpatialConductionError::Invalid(
                "terminals have no common conducting bulk/contact/electrode connection",
            ));
        }
        check(cancelled)?;
        Ok((faces, vertices))
    }
}

impl<C: SpatialConductorMaterial> SpatialResistance<'_, C> {
    /// Accepted native field at the actual reference excitation.
    pub fn reference_field(&self) -> &SpatialConductionField {
        &self.field
    }
    pub fn reference_voltage_v(&self) -> f64 {
        self.reference_voltage_v
    }
    pub fn resistance_ohm(&self) -> f64 {
        self.resistance_ohm
    }

    /// Recover native E/J, contact work and residuals without another Krylov
    /// solve. Relative potentials avoid subtracting published absolute values.
    /// The returned field has no linear-solve report; the original report
    /// belongs only to the retained reference field at its actual amplitude.
    pub fn field_at_voltage(
        &self,
        voltage_v: f64,
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductionField> {
        check(cancelled)?;
        if !voltage_v.is_finite() || !positive(current_tolerance_a) || !positive(power_tolerance_w)
        {
            return Err(SpatialConductionError::Invalid(
                "finite terminal voltage and positive actual-field tolerances required",
            ));
        }
        let mut relative = Vec::with_capacity(self.relative.len());
        for &value in &self.relative {
            check(cancelled)?;
            relative.push(finite((value / self.reference_voltage_v) * voltage_v)?);
        }
        let boundary = voltage_boundary(&self.faces, &[voltage_v, 0.0]);
        let g = &self.problem.geometry;
        let traces = prepare_scalar_boundary_traces(
            &g.mesh,
            &g.faces,
            &self.problem.components,
            boundary.voltage.iter().map(|t| (t.face, t.voltage_v)),
            std::iter::empty(),
            current_tolerance_a,
            cancelled,
        )?;
        self.problem.recover_field(
            relative,
            traces,
            current_tolerance_a,
            power_tolerance_w,
            None,
            cancelled,
        )
    }

    /// Recover at signed current, independently checking both terminal
    /// reactions and authored I*V against actual Joule heat. A successful zero
    /// drive has zero field and heat; the reference resistance remains finite.
    pub fn field_at_current(
        &self,
        current_a: f64,
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialConductionField> {
        check(cancelled)?;
        if !current_a.is_finite() {
            return Err(SpatialConductionError::Invalid(
                "finite total current required",
            ));
        }
        let voltage = finite(current_a * self.resistance_ohm)?;
        let power = finite(current_a * voltage)?;
        if current_a != 0.0 && (voltage == 0.0 || !positive(power)) {
            return Err(SpatialConductionError::Invalid(
                "total-current drive is not representable",
            ));
        }
        let field =
            self.field_at_voltage(voltage, current_tolerance_a, power_tolerance_w, cancelled)?;
        let mut current_residual_a = 0.0_f64;
        for (vertices, expected) in self.vertices.iter().zip([current_a, -current_a]) {
            let defect = finite(terminal_current(&field, vertices, cancelled)? - expected)?;
            if defect.abs() > current_residual_a.abs() {
                current_residual_a = defect;
            }
        }
        let power_residual_w = finite(field.power_w - power)?;
        if current_residual_a.abs() > current_tolerance_a
            || power_residual_w.abs() > power_tolerance_w
        {
            return Err(SpatialConductionError::Balance {
                current_residual_a,
                maximum_free_current_residual_a: field.maximum_free_current_residual_a,
                power_residual_w,
                hodge_power_residual_w: field.hodge_power_residual_w,
            });
        }
        check(cancelled)?;
        Ok(field)
    }
}

pub(super) fn voltage_boundary(
    faces: &[Vec<usize>],
    voltages: &[f64],
) -> SpatialConductionBoundary {
    SpatialConductionBoundary {
        voltage: faces
            .iter()
            .zip(voltages.iter().copied())
            .flat_map(|(faces, voltage)| {
                faces.iter().map(move |&face| ConductionVoltageFace {
                    face,
                    voltage_v: [voltage; 3],
                })
            })
            .collect(),
        current: Vec::new(),
    }
}

pub(super) fn terminal_current(
    field: &SpatialConductionField,
    vertices: &[usize],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let mut current = 0.0;
    for &vertex in vertices {
        check(cancelled)?;
        current = finite(
            current
                + field.nodal_boundary_current_a[vertex].ok_or(SpatialConductionError::Invalid(
                    "lost original terminal reaction",
                ))?,
        )?;
    }
    Ok(current)
}

pub(super) fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
