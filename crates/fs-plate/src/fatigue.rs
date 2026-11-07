//! Stress-life observation on the original, source-bound DKT bending field.
//!
//! Homogeneous isotropic sections are compiled from their retained material
//! states into the SAME plate pencil used by static and dynamic consumers.
//! A gauge follows one original midside quadrature point and signed depth.
//! For the native DKT convention, `epsilon(z) = z B q`; a positive Hessian of
//! transverse displacement compresses the positive-z face. Engineering twist
//! is converted by the original plane-stress stiffness, never treated as a
//! tensor shear component. Uniform membrane pretension contributes `T/h` to
//! both in-plane normal stresses, including at zero transverse displacement.
//!
//! The section relation is `D = Q h^3/12`, as in Roylance, *Laminated Composite
//! Plates*, MIT, equations (2), (17)-(19):
//! <https://web.mit.edu/course/3/3.11/www/modules/laminates.pdf>.
//! This adapter admits homogeneous elastic sections, not laminate-effective
//! properties, plastic stress reconstruction, material viscosity or thermal
//! cycling. Fixed-plane Miner damage remains an empirical estimate; it does
//! not create cracks, soften stiffness, deposit heat or certify service life.

use std::sync::Arc;

use fs_material::fatigue::{
    FatigueError, FatigueSummary, MeanStressCorrection, NormalStressDirection, RainflowLimits,
    RainflowMiner, SnCurve,
};
use fs_material::state_point::IsotropicElasticStatePoint;
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use crate::{
    AssemblyOptions, PlateError, PlateMesh, PlateModel, PlateSection, PlateSectionField,
    assemble_with_sections, dkt_curvature_quadrature, edge_coeffs, h_vectors,
};

/// One original homogeneous material through the local centered thickness.
#[derive(Debug, Clone)]
pub struct PlateFatigueSection {
    pub material: Arc<IsotropicElasticStatePoint>,
    pub thickness_m: f64,
}

/// A physical material point and a fixed in-plane normal-stress direction.
#[derive(Debug, Clone)]
pub struct PlateFatigueGauge {
    pub element: usize,
    /// Original DKT midside quadrature index, in `0..3`.
    pub surface_point: usize,
    /// Signed distance from the midsurface, within `[-h/2, h/2]` [m].
    pub depth_m: f64,
    /// Direction in the original mesh x/y frame; normalized once at admission.
    pub normal_in_mesh: [f64; 2],
    /// Exact retained absolute-temperature axis, including its quantity kind.
    pub temperature_axis: String,
    pub curve: Arc<SnCurve>,
    pub mean_stress: MeanStressCorrection,
}

/// Explicit geometry, small-deformation and bounded observation controls.
#[derive(Debug, Clone, Copy)]
pub struct PlateFatigueConfig {
    pub history: RainflowLimits,
    pub maximum_elements: usize,
    pub maximum_full_dofs: usize,
    pub maximum_gauges: usize,
    pub maximum_total_residue: usize,
    /// Bound on the complete quadratic DKT rotation field, in `(0, 0.05]`.
    pub maximum_slope: f64,
    /// Bound on each engineering strain component, including pretension.
    /// Checked throughout every triangle and its thickness; in `(0, 0.05]`.
    pub maximum_absolute_strain: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum PlateFatigueError {
    Invalid(&'static str),
    Budget,
    SourceMismatch { gauge: usize, what: &'static str },
    Gauge { gauge: usize, error: FatigueError },
    Plate(PlateError),
    Cancelled,
    Finished,
}
impl core::fmt::Display for PlateFatigueError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native DKT fatigue refused: {self:?}")
    }
}
impl std::error::Error for PlateFatigueError {}
impl From<PlateError> for PlateFatigueError {
    fn from(error: PlateError) -> Self {
        if error == PlateError::Cancelled {
            Self::Cancelled
        } else {
            Self::Plate(error)
        }
    }
}

#[derive(Debug)]
struct Element {
    dofs: [Option<usize>; 9],
    rotations: [[[f64; 9]; 2]; 6],
    corner_surface_strain: [[[f64; 9]; 3]; 3],
    membrane_strain: f64,
}

#[derive(Debug)]
struct Gauge {
    binding: PlateFatigueGauge,
    dofs: [Option<usize>; 9],
    stress: [f64; 9],
    membrane_stress_pa: f64,
}

/// Immutable physical section sources, native pencil and prepared stress rows.
/// No caller-supplied assembled matrix can silently replace these materials.
#[derive(Debug)]
pub struct FatiguePlate {
    model: PlateModel,
    sections: Vec<PlateFatigueSection>,
    elements: Vec<Element>,
    gauges: Vec<Gauge>,
    config: PlateFatigueConfig,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlateFatigueError> {
    if cancelled() {
        Err(PlateFatigueError::Cancelled)
    } else {
        Ok(())
    }
}
fn dot(row: &[f64; 9], values: &[f64; 9]) -> f64 {
    row.iter().zip(values).map(|(a, b)| a * b).sum()
}
fn local(dofs: &[Option<usize>; 9], q: &[f64]) -> [f64; 9] {
    dofs.map(|d| d.map_or(0.0, |i| q[i]))
}
fn gauge_error(gauge: usize, error: FatigueError) -> PlateFatigueError {
    if error == FatigueError::Cancelled {
        PlateFatigueError::Cancelled
    } else {
        PlateFatigueError::Gauge { gauge, error }
    }
}

impl FatiguePlate {
    /// Compile native DKT stiffness/inertia and signed stress recovery together.
    /// Sections cover every element. Supports are the original homogeneous
    /// plate supports; beams, laminates and imposed thermal strains are absent.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn assemble(
        mesh: &PlateMesh,
        sections: Vec<PlateFatigueSection>,
        boundary: &[usize],
        options: &AssemblyOptions,
        gauges: Vec<PlateFatigueGauge>,
        config: PlateFatigueConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PlateFatigueError> {
        poll(cancelled)?;
        let invalid = PlateFatigueError::Invalid;
        let full = mesh
            .nodes
            .len()
            .checked_mul(3)
            .ok_or(PlateFatigueError::Budget)?;
        if mesh.tris.len() > config.maximum_elements
            || full > config.maximum_full_dofs
            || gauges.is_empty()
            || gauges.len() > config.maximum_gauges
            || gauges
                .len()
                .checked_mul(config.history.maximum_residue)
                .is_none_or(|n| n > config.maximum_total_residue)
        {
            return Err(PlateFatigueError::Budget);
        }
        if sections.len() != mesh.tris.len()
            || !config.maximum_slope.is_finite()
            || !(0.0 < config.maximum_slope && config.maximum_slope <= 0.05)
            || !config.maximum_absolute_strain.is_finite()
            || !(0.0 < config.maximum_absolute_strain && config.maximum_absolute_strain <= 0.05)
            || !options.pretension.is_finite()
        {
            return Err(invalid(
                "complete sections and finite small-slope/strain bounds required",
            ));
        }
        mesh.validate()?;
        let mut native_sections = Vec::with_capacity(sections.len());
        for section in &sections {
            poll(cancelled)?;
            let m = &section.material;
            if m.resolved()
                .properties()
                .iter()
                .any(|p| !p.is_static_at_fixed_external_conditions())
            {
                return Err(invalid(
                    "plate elasticity must have fixed external source conditions",
                ));
            }
            native_sections.push(PlateSection::orthotropic_plane_stress(
                m.young_modulus_pa(),
                m.young_modulus_pa(),
                m.poisson_ratio(),
                m.young_modulus_pa() / (2.0 * (1.0 + m.poisson_ratio())),
                section.thickness_m,
                m.density_kg_m3(),
            )?);
        }
        // Validate all gauge/source bindings before the matrix assembly.
        for (index, gauge) in gauges.iter().enumerate() {
            poll(cancelled)?;
            let section = sections
                .get(gauge.element)
                .ok_or_else(|| invalid("fatigue element is outside the original plate"))?;
            if gauge.surface_point >= 3
                || !gauge.depth_m.is_finite()
                || gauge.depth_m.abs() > 0.5 * section.thickness_m
            {
                return Err(invalid(
                    "fatigue point must lie in its original DKT section",
                ));
            }
            let material = section.material.resolved();
            let scope = gauge.curve.scope();
            let mismatch = |what| PlateFatigueError::SourceMismatch { gauge: index, what };
            let absolute = material
                .axis_quantities()
                .get(&gauge.temperature_axis)
                .is_some_and(|q| {
                    *q == QuantitySpec::dimensional(Temperature::DIMS)
                        || q.semantic_type().is_some_and(|s| {
                            s.kind() == QuantityKind::AbsoluteTemperature
                                && s.form() == ValueForm::Static
                        })
                });
            if !absolute
                || material.state_coordinate(&gauge.temperature_axis) != Some(scope.temperature_k)
                || scope.material
                    != format!(
                        "{}; {}",
                        material.material().chemistry,
                        material.material().phase
                    )
                || scope.process != material.material().process
            {
                return Err(mismatch(
                    "fatigue material, process and absolute temperature must match the section",
                ));
            }
            if gauge.curve.source().is_some_and(|source| {
                source.request().temperature_axis != gauge.temperature_axis
                    || !source.matches_material_state(
                        material.card_identity(),
                        material.material(),
                        &material.query(),
                    )
            }) {
                return Err(mismatch(
                    "S-N source card, state or fixed physical context differs from the section",
                ));
            }
            NormalStressDirection::new([gauge.normal_in_mesh[0], gauge.normal_in_mesh[1], 0.0])
                .map_err(|e| gauge_error(index, e))?;
            RainflowMiner::new(gauge.curve.clone(), gauge.mean_stress, config.history)
                .map_err(|e| gauge_error(index, e))?;
        }
        poll(cancelled)?;
        let model = assemble_with_sections(
            mesh,
            PlateSectionField::PerElement(&native_sections),
            boundary,
            &[],
            options,
        )?;
        poll(cancelled)?;
        let mut elements = Vec::with_capacity(mesh.tris.len());
        let mut curvature = Vec::with_capacity(mesh.tris.len());
        for (e, tri) in mesh.tris.iter().enumerate() {
            poll(cancelled)?;
            let x = tri.map(|i| mesh.nodes[i].0);
            let y = tri.map(|i| mesh.nodes[i].1);
            let (b, twice_area) = dkt_curvature_quadrature(&x, &y, e)?;
            if !(twice_area > 0.0) {
                return Err(invalid("fatigue requires positively oriented triangles"));
            }
            let cf = edge_coeffs(&x, &y);
            let mut rotations = [[[0.0; 9]; 2]; 6];
            for (i, rotation) in rotations.iter_mut().enumerate() {
                let mut shape = [0.0; 6];
                shape[i] = 1.0;
                let (rx, ry) = h_vectors(&cf, &shape);
                *rotation = [rx, ry];
            }
            // Convert midside Lagrange values to quadratic Bernstein controls;
            // their convex hull bounds the full native rotation field.
            for (edge, (a, b)) in [(1, 2), (2, 0), (0, 1)].into_iter().enumerate() {
                for c in 0..2 {
                    for d in 0..9 {
                        rotations[3 + edge][c][d] = 2.0 * rotations[3 + edge][c][d]
                            - 0.5 * (rotations[a][c][d] + rotations[b][c][d]);
                    }
                }
            }
            let section = &sections[e];
            // Curvature is linear. Reconstruct its three vertex values from
            // midside quadrature, then bound both faces (and all interior z).
            let corner_surface_strain = std::array::from_fn(|corner| {
                let opposite = (corner + 1) % 3;
                std::array::from_fn(|a| {
                    std::array::from_fn(|j| {
                        0.5 * section.thickness_m
                            * (b[corner][a][j] + b[(corner + 2) % 3][a][j] - b[opposite][a][j])
                    })
                })
            });
            let membrane_strain = (1.0 - section.material.poisson_ratio())
                * (options.pretension / section.thickness_m)
                / section.material.young_modulus_pa();
            if rotations
                .iter()
                .flatten()
                .flatten()
                .chain(corner_surface_strain.iter().flatten().flatten())
                .any(|x| !x.is_finite())
                || !membrane_strain.is_finite()
            {
                return Err(invalid("unrepresentable plate recovery coefficients"));
            }
            elements.push(Element {
                dofs: std::array::from_fn(|i| model.dof_map[3 * tri[i / 3] + i % 3]),
                rotations,
                corner_surface_strain,
                membrane_strain,
            });
            curvature.push(b);
        }
        let mut prepared = Vec::with_capacity(gauges.len());
        for (i, binding) in gauges.into_iter().enumerate() {
            poll(cancelled)?;
            let section = &sections[binding.element];
            let material = &section.material;
            let normal = NormalStressDirection::new([
                binding.normal_in_mesh[0],
                binding.normal_in_mesh[1],
                0.0,
            ])
            .map_err(|e| gauge_error(i, e))?
            .unit_normal();
            let [nx, ny, _] = normal;
            let e = material.young_modulus_pa();
            let nu = material.poisson_ratio();
            let scale = e / (1.0 - nu * nu);
            let q = [
                scale * (nx * nx + nu * ny * ny),
                scale * (ny * ny + nu * nx * nx),
                e / (1.0 + nu) * nx * ny,
            ];
            let b = &curvature[binding.element][binding.surface_point];
            let stress = std::array::from_fn(|j| {
                binding.depth_m * (0..3).map(|a| q[a] * b[a][j]).sum::<f64>()
            });
            let membrane_stress_pa =
                (options.pretension / section.thickness_m) * (nx * nx + ny * ny);
            if stress.iter().any(|s| !s.is_finite()) || !membrane_stress_pa.is_finite() {
                return Err(invalid("unrepresentable plate stress recovery"));
            }
            prepared.push(Gauge {
                dofs: elements[binding.element].dofs,
                binding,
                stress,
                membrane_stress_pa,
            });
        }
        poll(cancelled)?;
        Ok(Self {
            model,
            sections,
            elements,
            gauges: prepared,
            config,
        })
    }

    pub fn model(&self) -> &PlateModel {
        &self.model
    }
    pub fn sections(&self) -> &[PlateFatigueSection] {
        &self.sections
    }
    pub fn gauge(&self, index: usize) -> Option<&PlateFatigueGauge> {
        self.gauges.get(index).map(|g| &g.binding)
    }
    pub const fn config(&self) -> PlateFatigueConfig {
        self.config
    }

    /// Recover signed stresses from the original reduced displacement/slopes.
    /// Recovery is read-only and contributes no sample or cycle by itself.
    pub fn stresses(
        &self,
        displacement: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, PlateFatigueError> {
        poll(cancelled)?;
        let invalid = PlateFatigueError::Invalid;
        if displacement.len() != self.model.free || displacement.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "fatigue requires finite original reduced displacement/slopes",
            ));
        }
        for element in &self.elements {
            poll(cancelled)?;
            let u = local(&element.dofs, displacement);
            if element.rotations.iter().any(|r| {
                let slope = dot(&r[0], &u).hypot(dot(&r[1], &u));
                !slope.is_finite() || slope > self.config.maximum_slope
            }) {
                return Err(invalid("DKT rotation field exceeds the small-slope domain"));
            }
            for corner in &element.corner_surface_strain {
                for (a, row) in corner.iter().enumerate() {
                    let bending = dot(row, &u);
                    let membrane = if a < 2 { element.membrane_strain } else { 0.0 };
                    if !bending.is_finite()
                        || !membrane.is_finite()
                        || membrane.abs() + bending.abs() > self.config.maximum_absolute_strain
                    {
                        return Err(invalid(
                            "plate section exceeds the declared small-strain domain",
                        ));
                    }
                }
            }
        }
        let mut stresses = Vec::with_capacity(self.gauges.len());
        for gauge in &self.gauges {
            poll(cancelled)?;
            let stress =
                gauge.membrane_stress_pa + dot(&gauge.stress, &local(&gauge.dofs, displacement));
            if !stress.is_finite() {
                return Err(invalid("plate fatigue stress overflow"));
            }
            stresses.push(stress);
        }
        Ok(stresses)
    }

    /// Start observation at an explicitly supplied accepted native state.
    pub fn observe<'a>(
        &'a self,
        displacement: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlateFatigueObservation<'a>, PlateFatigueError> {
        poll(cancelled)?;
        let counters = self
            .gauges
            .iter()
            .enumerate()
            .map(|(i, g)| {
                RainflowMiner::new(
                    g.binding.curve.clone(),
                    g.binding.mean_stress,
                    self.config.history,
                )
                .map_err(|e| gauge_error(i, e))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut observation = PlateFatigueObservation {
            plate: self,
            counters,
            stresses: Vec::new(),
            finished: false,
        };
        observation.push(displacement, cancelled)?;
        Ok(observation)
    }
}

/// Bounded histories on one immutable plate. Clone retains open reversals.
#[derive(Debug, Clone)]
pub struct PlateFatigueObservation<'a> {
    plate: &'a FatiguePlate,
    counters: Vec<RainflowMiner>,
    stresses: Vec<f64>,
    finished: bool,
}
impl<'a> PlateFatigueObservation<'a> {
    pub fn plate(&self) -> &'a FatiguePlate {
        self.plate
    }
    pub fn stresses_pa(&self) -> &[f64] {
        &self.stresses
    }
    pub const fn is_finished(&self) -> bool {
        self.finished
    }
    pub fn summaries(&self) -> Vec<FatigueSummary> {
        self.counters.iter().map(RainflowMiner::summary).collect()
    }

    /// Observe an accepted native state once; every gauge accepts atomically.
    /// Processing boundaries must not call `finish` between these samples.
    pub fn push(
        &mut self,
        displacement: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), PlateFatigueError> {
        if self.finished {
            return Err(PlateFatigueError::Finished);
        }
        let stresses = self.plate.stresses(displacement, cancelled)?;
        let mut next = self.counters.clone();
        for (i, (counter, &stress)) in next.iter_mut().zip(&stresses).enumerate() {
            counter
                .push_with_cancel(stress, cancelled)
                .map_err(|e| gauge_error(i, e))?;
        }
        poll(cancelled)?;
        self.counters = next;
        self.stresses = stresses;
        Ok(())
    }

    /// Read-only damage for the accepted prefix, including endpoint half cycles.
    pub fn preview_finish(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<FatigueSummary>, PlateFatigueError> {
        poll(cancelled)?;
        self.counters
            .iter()
            .enumerate()
            .map(|(i, c)| {
                c.preview_finish_with_cancel(cancelled)
                    .map_err(|e| gauge_error(i, e))
            })
            .collect()
    }

    /// Close the complete observation once; no partial gauge finalization.
    pub fn finish(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<FatigueSummary>, PlateFatigueError> {
        poll(cancelled)?;
        let mut next = self.counters.clone();
        let summaries = next
            .iter_mut()
            .enumerate()
            .map(|(i, c)| {
                c.finish_with_cancel(cancelled)
                    .map_err(|e| gauge_error(i, e))
            })
            .collect::<Result<Vec<_>, _>>()?;
        poll(cancelled)?;
        self.counters = next;
        self.finished = true;
        Ok(summaries)
    }
}
