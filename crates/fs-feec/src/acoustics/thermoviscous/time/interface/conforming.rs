//! Complete matching P2 traces of the original resolved equations.
//!
//! A Cartesian traction t [Pa] is applied to the left and -t to the right.
//! Conductive flux q [W/m2] is inward to the left and outward from the right.
//! The original surface mass makes their duals integral(phi v dA) and
//! integral(phi theta dA), theta=T'/T0. Thus the common-flow constraint is
//! equality of the COMPLETE active P2 velocity and temperature traces. The
//! heat term in the acoustic work ledger is theta*q, not the caloric flux q.
use super::*;

/// Original hierarchical P2 coordinates in a common physical ordering.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConformingThermoviscousInterface {
    /// Left and right global P2 DOF numbers: geometric vertices followed by
    /// geometric edges. Edge bubbles are matched by their endpoint positions,
    /// never by local node numbering or a midpoint-only coincidence test.
    pub trace_dofs: [Vec<usize>; 2],
    /// Indices into `trace_dofs[side]` carrying the active P2 test space.
    /// With DOFs constrained to zero on BOTH sides, traction is represented
    /// in that active space, modulo loads doing no work on the admissible
    /// trace. Omitted coefficients do NOT assert zero physical stress.
    /// A one-sided zero trace still constrains the free side.
    pub traction_indices: Vec<usize>,
    /// The analogous active P2 conductive-flux basis. Temperature is the
    /// translational temperature; molecular reservoirs retain their own fields.
    pub heat_flux_indices: Vec<usize>,
}

/// Same complete trajectory and original work owners as the normal connection.
pub struct PreparedConformingThermoviscousTimePair<'p, 'a, 'c> {
    inner: PreparedThermoviscousTimePair<'p, 'a, 'c>,
    layout: ConformingThermoviscousInterface,
    heat_flux_areas_m2: Vec<f64>,
}

/// Fully admitted fields and separately integrated physical interface transfers.
#[derive(Clone, Debug)]
pub struct ConformingThermoviscousTimePairStep {
    pub regions: [ThermoviscousTimeStep; 2],
    /// Integral(t.v dA dt), inward to each region [J].
    pub mechanical_interface_work_j: [f64; 2],
    /// Integral((T'/T0) q_in dA dt) [J], the signed transfer of acoustic
    /// available energy. This is distinct from actual caloric heat below.
    pub conductive_interface_work_j: [f64; 2],
    pub external_pressure_work_j: f64,
    /// P2 vector traction ON THE LEFT [Pa], in `traction_indices` order.
    /// The traction on the right is its negative in the same Cartesian frame.
    /// On a constrained trace this is the active-space Galerkin representative
    /// modulo work-null loads; omitted coefficients are not physical zero stress.
    pub mean_traction_pa: Vec<[f64; 3]>,
    pub endpoint_traction_pa: Vec<[f64; 3]>,
    /// P2 conductive flux INWARD TO THE LEFT [W/m2], `heat_flux_indices` order.
    /// On a constrained thermal trace this is only the active-space Galerkin
    /// representative. Unobservable flux through fixed-temperature DOFs is not
    /// inferred or silently assigned zero.
    pub mean_inward_heat_flux_w_m2: Vec<f64>,
    pub endpoint_inward_heat_flux_w_m2: Vec<f64>,
    /// Integral(q_in dA dt) [J] for each region, only when the COMPLETE thermal
    /// trace is active. Otherwise a work-null load can have a nonzero integral
    /// and the actual caloric transfer is unidentifiable (`None`). Identified
    /// transfers cancel and are NOT added to the quadratic acoustic ledger.
    pub caloric_heat_transfer_j: Option<[f64; 2]>,
}
impl ConformingThermoviscousTimePairStep {
    pub fn balance_residual_j(&self) -> f64 {
        self.regions
            .iter()
            .map(|r| r.work.delta_energy_j + r.work.dissipated_j() - r.work.structural_work_j)
            .sum::<f64>()
            - self.external_pressure_work_j
    }
}

impl<'a, 'c> PreparedThermoviscousTime<'a, 'c> {
    /// Join matching original face triangulations with complete P2 vector and
    /// temperature continuity. Every pressure mode touching the cut is removed
    /// from the external inputs; such a mode must lie entirely on the cut.
    /// Other pressure and structural-force coordinates retain their ordering.
    ///
    /// Original stationary constraints at the edge are retained, including a
    /// constraint present on only one side. A moving structural trace sharing
    /// a cut DOF is refused: its kinematic constraint must have one owner.
    /// Pressure remains in each region's original P1 volume space; it is not a
    /// new boundary constraint. Total traction, rather than pressure alone,
    /// balances across a viscous material interface. No mean-flow, geometry,
    /// species transport, mean thermal reservoir or contact resistance is added.
    pub fn connect_conforming<'p>(
        &'p self,
        right: &'p Self,
        left_faces: &[[u32; 3]],
        right_faces: &[[u32; 3]],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedConformingThermoviscousTimePair<'p, 'a, 'c>, ThermoviscousTimeError> {
        let regions = [self, right];
        let n = self.system.state_dim() + right.system.state_dim();
        if n > limits.max_states || n >= MAX_LINEAR_EXPONENTIAL_BLOCK_DIM / 2 {
            return Err(invalid(
                "complete conforming interface state capacity required",
            ));
        }
        let (layout, selected, columns, heat_flux_areas_m2) =
            binding(regions, [left_faces, right_faces], cancelled)?;
        let inner = PreparedThermoviscousTimePair::prepare(
            regions,
            selected,
            Vec::new(),
            Some(columns),
            limits,
            cancelled,
        )?;
        Ok(PreparedConformingThermoviscousTimePair {
            inner,
            layout,
            heat_flux_areas_m2,
        })
    }
}

impl<'p, 'a, 'c> PreparedConformingThermoviscousTimePair<'p, 'a, 'c> {
    pub fn interface(&self) -> &ConformingThermoviscousInterface {
        &self.layout
    }
    pub fn system(&self) -> &PortHamiltonian {
        self.inner.system()
    }
    pub fn duration_s(&self) -> f64 {
        self.inner.duration_s()
    }
    pub fn external_pressure_modes(&self) -> &[(usize, usize)] {
        self.inner.external_pressure_modes()
    }
    pub fn zero_states(&self) -> [ThermoviscousTimeState; 2] {
        self.inner.zero_states()
    }

    /// Both original states remain untouched on every refusal or cancellation.
    /// Initial trace mismatch is rejected; it is never projected away.
    pub fn step(
        &self,
        initial: [&ThermoviscousTimeState; 2],
        external_traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: [&[f64]; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConformingThermoviscousTimePairStep, ThermoviscousTimeError> {
        let frame = self.inner.evaluate(
            initial,
            external_traction_pa,
            wall_forces_n_per_sqrt_kg,
            cancelled,
        )?;
        let nv = self.layout.traction_indices.len();
        let vectors = |reaction: &[f64]| {
            (0..nv)
                .map(|i| std::array::from_fn(|c| reaction[c * nv + i]))
                .collect::<Vec<_>>()
        };
        let mean_traction_pa = vectors(&frame.mean_reaction);
        let endpoint_traction_pa = vectors(&frame.endpoint_reaction);
        let mean_inward_heat_flux_w_m2 = frame.mean_reaction[3 * nv..].to_vec();
        let endpoint_inward_heat_flux_w_m2 = frame.endpoint_reaction[3 * nv..].to_vec();
        let caloric_heat_transfer_j =
            if self.layout.heat_flux_indices.len() == self.layout.trace_dofs[0].len() {
                let transfer: f64 = mean_inward_heat_flux_w_m2
                    .iter()
                    .zip(&self.heat_flux_areas_m2)
                    .map(|(q, area)| q * area * self.duration_s())
                    .sum();
                if !transfer.is_finite() {
                    return Err(invalid(
                        "original conforming caloric heat transfer is nonfinite",
                    ));
                }
                Some([transfer, -transfer])
            } else {
                None
            };
        poll(cancelled)?;
        Ok(ConformingThermoviscousTimePairStep {
            regions: frame.regions,
            mechanical_interface_work_j: frame.mechanical_work_j,
            conductive_interface_work_j: frame.conductive_work_j,
            external_pressure_work_j: frame.external_pressure_work_j,
            mean_traction_pa,
            endpoint_traction_pa,
            mean_inward_heat_flux_w_m2,
            endpoint_inward_heat_flux_w_m2,
            caloric_heat_transfer_j,
        })
    }
}

type PointKey = [u64; 3];
type DofKey = (u8, PointKey, PointKey);
fn point_key(point: [f64; 3]) -> PointKey {
    point.map(|x| {
        if x == 0.0 {
            0.0_f64.to_bits()
        } else {
            x.to_bits()
        }
    })
}
fn face_key(model: &ThermoviscousFe<'_>, face: [u32; 3]) -> FaceKey {
    let mut key = face.map(|v| point_key(model.fe.positions()[v as usize]));
    key.sort_unstable();
    key
}
struct Surface {
    faces: BTreeMap<FaceKey, ([f64; 3], f64)>,
    dofs: BTreeMap<DofKey, usize>,
    selected: Vec<usize>,
}
fn surface(
    model: &ThermoviscousFe<'_>,
    selected_faces: &[[u32; 3]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Surface, ThermoviscousTimeError> {
    if selected_faces.is_empty() || selected_faces.len() > model.space.complex.faces.len() {
        return Err(invalid("nonempty original boundary interface required"));
    }
    let mut faces = BTreeMap::new();
    let mut dofs = BTreeMap::new();
    let mut insert = |key, global| -> Result<(), ThermoviscousTimeError> {
        if dofs
            .insert(key, global)
            .is_some_and(|prior| prior != global)
        {
            return Err(invalid(
                "coincident distinct original interface coordinates refused",
            ));
        }
        Ok(())
    };
    for &face in selected_faces {
        poll(cancelled)?;
        let (_, normal, area) = model
            .fe
            .boundary_geometry(face)
            .map_err(ThermoviscousError::from)?;
        if model.wall_faces.contains(&canonical(face))
            || faces
                .insert(face_key(model, face), (normal, area))
                .is_some()
        {
            return Err(invalid("unique original free interface faces required"));
        }
        for &v in &face {
            let p = point_key(model.fe.positions()[v as usize]);
            insert((0, p, p), v as usize)?;
        }
        for a in 0..3 {
            for b in a + 1..3 {
                let mut edge = [face[a], face[b]];
                edge.sort_unstable();
                let id = model
                    .space
                    .complex
                    .edges
                    .binary_search(&edge)
                    .map_err(|_| invalid("original P2 interface edge required"))?;
                let mut points = edge.map(|v| point_key(model.fe.positions()[v as usize]));
                points.sort_unstable();
                insert((1, points[0], points[1]), model.space.edge_off + id)?;
            }
        }
    }
    let mut selected = Vec::new();
    for (mode, traces) in model.physical_modes.iter().enumerate() {
        let overlaps = traces
            .iter()
            .filter(|trace| faces.contains_key(&trace.vertices))
            .count();
        if overlaps > 0 {
            if overlaps != traces.len() {
                return Err(invalid(
                    "one pressure mode cannot straddle internal and external faces",
                ));
            }
            selected.push(mode);
        }
    }
    if selected.is_empty() {
        return Err(invalid(
            "original interface pressure-face ownership required",
        ));
    }
    Ok(Surface {
        faces,
        dofs,
        selected,
    })
}

fn binding(
    regions: [&PreparedThermoviscousTime<'_, '_>; 2],
    faces: [&[[u32; 3]]; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<
    (
        ConformingThermoviscousInterface,
        [Vec<usize>; 2],
        InterfaceColumns,
        Vec<f64>,
    ),
    ThermoviscousTimeError,
> {
    let surfaces = [
        surface(regions[0].model, faces[0], cancelled)?,
        surface(regions[1].model, faces[1], cancelled)?,
    ];
    if surfaces[0].faces.keys().ne(surfaces[1].faces.keys())
        || surfaces[0].dofs.keys().ne(surfaces[1].dofs.keys())
    {
        return Err(invalid(
            "complete original P2 interface triangulations must match",
        ));
    }
    for (key, &(normal, area)) in &surfaces[0].faces {
        poll(cancelled)?;
        let (other, other_area) = surfaces[1].faces[key];
        let dot: f64 = normal.iter().zip(other).map(|(a, b)| a * b).sum();
        if (dot + 1.0).abs() > 1e-12 || (area - other_area).abs() > 1e-12 * area {
            return Err(invalid(
                "opposite interface normals and equal original areas required",
            ));
        }
    }
    let trace_dofs: [Vec<_>; 2] =
        std::array::from_fn(|s| surfaces[s].dofs.values().copied().collect());
    let m = trace_dofs[0].len();
    let sizes = regions.map(|region| region.system.state_dim());
    let starts = [0, sizes[0]];
    let n = sizes.iter().sum::<usize>();
    for side in 0..2 {
        if let Some(wall) = &regions[side].model.moving_wall {
            for &dof in &trace_dofs[side] {
                if wall
                    .velocity_shapes
                    .iter()
                    .any(|shape| shape[dof].iter().any(|x| *x != 0.0))
                {
                    return Err(invalid(
                        "moving structural and fluid interface traces require one kinematic owner",
                    ));
                }
            }
        }
    }
    let traction_indices: Vec<_> = (0..m)
        .filter(|&i| (0..2).any(|s| regions[s].model.velocity_map[trace_dofs[s][i]].is_some()))
        .collect();
    let heat_flux_indices: Vec<_> = (0..m)
        .filter(|&i| (0..2).any(|s| regions[s].model.temperature_map[trace_dofs[s][i]].is_some()))
        .collect();
    let mechanical = 3 * traction_indices.len();
    let count = mechanical + heat_flux_indices.len();
    if count == 0 || count > n {
        return Err(invalid(
            "nonzero complete original conforming trace rank required",
        ));
    }
    // Reuse the ORIGINAL P2 surface quadrature. Compare the complete physical
    // matrices under the geometry-derived permutation before sharing the left
    // matrix. No fitted area, diagonal/lumped mass or reduced trace is used.
    let mut mass = [vec![0.0; m * m], vec![0.0; m * m]];
    for side in 0..2 {
        let source = regions[side]
            .model
            .fe
            .boundary_mass(faces[side])
            .map_err(ThermoviscousError::from)?;
        let indices: BTreeMap<_, _> = trace_dofs[side]
            .iter()
            .enumerate()
            .map(|(i, &g)| (g, i))
            .collect();
        for (i, &global) in trace_dofs[side].iter().enumerate() {
            poll(cancelled)?;
            let (columns, values) = source.row(global);
            for (&column, &value) in columns.iter().zip(values) {
                if let Some(&j) = indices.get(&column) {
                    mass[side][i * m + j] = value;
                }
            }
        }
    }
    for (a, b) in mass[0].iter().zip(&mass[1]) {
        if !a.is_finite()
            || !b.is_finite()
            || (a - b).abs() > 512.0 * f64::EPSILON * (a.abs() + b.abs())
        {
            return Err(invalid(
                "original P2 physical surface mass must agree without truncation",
            ));
        }
    }
    let mut loads = vec![0.0; n * count];
    for side in 0..2 {
        let source = regions[side];
        let model = source.model;
        let size = sizes[side];
        let sign = if side == 0 { 1.0 } else { -1.0 };
        let mut physical = vec![0.0; size * count];
        for (i, &global) in trace_dofs[side].iter().enumerate() {
            if let Some(v) = model.velocity_map[global] {
                for c in 0..3 {
                    let row = c * model.velocity_count + v;
                    for (j, &trace) in traction_indices.iter().enumerate() {
                        physical[row * count + c * traction_indices.len() + j] =
                            sign * mass[0][i * m + trace] / source.scale[row];
                    }
                }
            }
            if let Some(t) = model.temperature_map[global] {
                let row = model.temperature_start + t;
                for (j, &trace) in heat_flux_indices.iter().enumerate() {
                    // The original thermal equation and storage use theta=T'/T0;
                    // its physical load is q, so NO extra T0 belongs here.
                    physical[row * count + mechanical + j] =
                        sign * mass[0][i * m + trace] / source.scale[row];
                }
            }
        }
        for row in 0..size {
            poll(cancelled)?;
            for a in 0..count {
                loads[(starts[side] + row) * count + a] = (0..size)
                    .map(|k| source.inverse_mass[row * size + k] * physical[k * count + a])
                    .sum();
            }
        }
    }
    let (areas, _) = regions[0]
        .model
        .fe
        .boundary_load(faces[0])
        .map_err(ThermoviscousError::from)?;
    let heat_flux_areas_m2 = heat_flux_indices
        .iter()
        .map(|&i| areas[trace_dofs[0][i]])
        .collect();
    let selected = surfaces.map(|surface| surface.selected);
    Ok((
        ConformingThermoviscousInterface {
            trace_dofs,
            traction_indices,
            heat_flux_indices,
        },
        selected,
        InterfaceColumns {
            count,
            mechanical,
            loads,
        },
        heat_flux_areas_m2,
    ))
}

#[cfg(test)]
mod tests;
