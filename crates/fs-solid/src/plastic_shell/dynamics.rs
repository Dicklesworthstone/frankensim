//! Backward-Euler motion on the original retained CST/DKT plastic shell.
//!
//! The native shell inertia owns reference thickness/density integration and
//! signed translation/director coupling. The original solid inertia supplies
//! momentum, its exact derivative, residual scaling and the kinetic-work
//! identity. The shell's accepted histories and Newton owner remain unchanged.
//! Numerical relative-drilling inertia has its own energy account; it cannot
//! become material heat or physical normal-spin inertia.

use super::{Element, SmallStrainShellProblem, SmallStrainShellSolution, checkpoint, invalid};
use crate::inertia::{BackwardEulerInertia, InertialEnergy, ReferenceMass};
use crate::plastic3::SmallStrainTetError;
use fs_exec::Cx;
use fs_material::SmallStrainLaw;
use fs_plate::shell::inertia::{prepare_reference_mass, reference_kinetic_increment};
use fs_plate::shell::{ReferenceDensityMoments, ShellMesh};
use fs_sparse::Coo;

/// Original density in kg/m³, constant within each physical material layer.
/// Values describe the closed reference inventory, never the current volume.
#[derive(Debug, Clone, Copy)]
pub enum PlasticShellDensity<'a> {
    /// The same positive reference density in every element and layer.
    Uniform(f64),
    /// One density per layer, shared by every element's section stack.
    Layers(&'a [f64]),
    /// One complete layer-density vector per element, in connectivity order.
    PerElement(&'a [Vec<f64>]),
}

#[derive(Debug)]
struct FacetInventory {
    bounds: Vec<[f64; 2]>,
    densities_kg_m3: Vec<f64>,
    thickness_m: f64,
    moments: ReferenceDensityMoments,
}

/// Immutable original mass, bound to the actual reference mesh, layer bounds
/// and numerical coordinate length. Prepare once and retain across steps.
/// Only one global sparse mass is stored; physical/gauge energy partitions
/// use the original native kinetic quadrature and compact facet inventory.
#[derive(Debug)]
pub struct PlasticShellReferenceInertia<'m> {
    mesh: &'m ShellMesh,
    coordinate_length_m: f64,
    facets: Vec<FacetInventory>,
    mass: ReferenceMass,
}

impl<'m> PlasticShellReferenceInertia<'m> {
    /// Bind mass to a mesh independently of temporary source-derived section
    /// laws. This is the same original inventory preparation as the problem
    /// method; the separately supplied mesh must match the problem exactly.
    pub fn prepare<L: SmallStrainLaw>(
        mesh: &'m ShellMesh,
        problem: &SmallStrainShellProblem<'_, L>,
        density: PlasticShellDensity<'_>,
        cx: &Cx<'_>,
    ) -> Result<Self, SmallStrainTetError> {
        checkpoint(cx)?;
        if !std::ptr::eq(mesh, problem.mesh) && mesh != problem.mesh {
            return Err(invalid(
                "reference inertia mesh differs from its original shell",
            ));
        }
        let (facets, mass) = problem.reference_inertia_data(density, cx)?;
        Ok(Self {
            mesh,
            coordinate_length_m: problem.coordinate_length_m,
            facets,
            mass,
        })
    }

    /// Actual original material mass, excluding numerical drilling inertia.
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }

    /// Original physical layer densities [kg/m³], in section order. Coupled
    /// source owners can bind the exact inventory, including eccentric moments;
    /// equal total mass alone does not identify the reference inertia.
    pub fn reference_layer_densities(&self, element: usize) -> Option<&[f64]> {
        self.facets
            .get(element)
            .map(|facet| facet.densities_kg_m3.as_slice())
    }

    fn admit<L: SmallStrainLaw>(
        &self,
        problem: &SmallStrainShellProblem<'_, L>,
        cx: &Cx<'_>,
    ) -> Result<(), SmallStrainTetError> {
        checkpoint(cx)?;
        if (!std::ptr::eq(self.mesh, problem.mesh) && self.mesh != problem.mesh)
            || self.coordinate_length_m != problem.coordinate_length_m
            || self.facets.len() != problem.mesh.tris.len()
        {
            return Err(invalid(
                "retained shell inertia belongs to different reference geometry or coordinates",
            ));
        }
        for (e, facet) in self.facets.iter().enumerate() {
            checkpoint(cx)?;
            let layers = problem.sections.at(e).layers();
            if layers.len() != facet.bounds.len()
                || layers
                    .iter()
                    .zip(&facet.bounds)
                    .any(|(layer, bounds)| [layer.bottom_z_m, layer.top_z_m] != *bounds)
            {
                return Err(invalid(
                    "retained shell inertia belongs to different physical layer bounds",
                ));
            }
        }
        Ok(())
    }
}

/// One contribution to the reference kinetic integral, in joules. For the
/// drilling gauge all storage and losses are numerical. For material inertia
/// only `numerical_loss_j` is the backward-Euler discretization loss.
pub use fs_plate::shell::inertia::ReferenceKineticIncrement as ShellKineticEnergy;

/// The original total kinetic-work receipt and independently integrated
/// physical/gauge contributions. Their sums must agree with the actual mass
/// used in the momentum solve; no physical energy is obtained by subtraction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellInertialEnergy {
    pub total: InertialEnergy,
    pub physical: ShellKineticEnergy,
    pub drilling_gauge: ShellKineticEnergy,
}

/// Complete proposed motion. Failure never changes the old endpoint or mass.
#[derive(Debug)]
pub struct PlasticShellDynamicStep<S> {
    pub solution: SmallStrainShellSolution<S>,
    /// Translations [m/s] and physical axial angular velocities [rad/s].
    pub velocity: Vec<[f64; 6]>,
    pub inertia: ShellInertialEnergy,
}

impl<'m, L: SmallStrainLaw> SmallStrainShellProblem<'m, L> {
    /// Prepare the full original mass before support elimination. Density
    /// belongs to each section's actual layers about its declared reference
    /// surface. The signed first moment is retained for eccentric stacks.
    pub fn prepare_reference_inertia(
        &self,
        density: PlasticShellDensity<'_>,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellReferenceInertia<'m>, SmallStrainTetError> {
        PlasticShellReferenceInertia::prepare(self.mesh, self, density, cx)
    }

    fn reference_inertia_data(
        &self,
        density: PlasticShellDensity<'_>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<FacetInventory>, ReferenceMass), SmallStrainTetError> {
        let elements = self.prepare(cx)?;
        if let PlasticShellDensity::PerElement(values) = density
            && values.len() != elements.len()
        {
            return Err(invalid(
                "shell reference density needs one layer inventory per element",
            ));
        }
        let ndofs = 6 * self.mesh.nodes.len();
        let mut matrix = Coo::new(ndofs, ndofs);
        let mut facets = Vec::with_capacity(elements.len());
        let mut total_mass = 0.0;
        for (e, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            let layers = self.sections.at(e).layers();
            let values = match density {
                PlasticShellDensity::Uniform(_) => None,
                PlasticShellDensity::Layers(values) => Some(values),
                PlasticShellDensity::PerElement(values) => Some(values[e].as_slice()),
            };
            if values.is_some_and(|values| values.len() != layers.len()) {
                return Err(invalid(
                    "shell reference density count differs from physical layer count",
                ));
            }
            let thickness = layers.last().expect("admitted section").top_z_m - layers[0].bottom_z_m;
            let mut q = [0.0; 3];
            let mut bounds = Vec::with_capacity(layers.len());
            let mut densities_kg_m3 = Vec::with_capacity(layers.len());
            for (i, layer) in layers.iter().enumerate() {
                let rho = match density {
                    PlasticShellDensity::Uniform(rho) => rho,
                    _ => values.expect("layer density field")[i],
                };
                if !rho.is_finite() || rho <= 0.0 {
                    return Err(invalid(
                        "shell reference layer density must be finite and positive",
                    ));
                }
                let width = layer.top_z_m - layer.bottom_z_m;
                let fraction = width / thickness;
                let center = (layer.bottom_z_m + 0.5 * width) / thickness;
                let weight = rho * fraction;
                if !fraction.is_finite() || fraction <= 0.0 || !weight.is_finite() || weight <= 0.0
                {
                    return Err(invalid("reference layer mass fraction is unrepresentable"));
                }
                q[0] += weight;
                q[1] += weight * center;
                q[2] += weight * (center * center + fraction * fraction / 12.0);
                bounds.push([layer.bottom_z_m, layer.top_z_m]);
                densities_kg_m3.push(rho);
            }
            let moments = ReferenceDensityMoments::new(q[0], q[1], q[2])
                .map_err(|e| invalid(e.to_string()))?;
            let native = prepare_reference_mass(&element.geometry, [thickness; 3], [moments; 3])
                .map_err(|e| invalid(e.to_string()))?;
            total_mass += native.reference_mass_kg;
            // Transform local physical axial rotations with the same frame as
            // translations, then apply S^-T M S^-1 for q=(u,L*theta).
            // Reflect the computed upper triangle to preserve exact symmetry.
            for i in 0..18 {
                let bi = 6 * (i / 6) + 3 * ((i % 6) / 3);
                let row = 6 * element.nodes[i / 6] + i % 6;
                for j in i..18 {
                    let bj = 6 * (j / 6) + 3 * ((j % 6) / 3);
                    let column = 6 * element.nodes[j / 6] + j % 6;
                    let mut value = 0.0;
                    for p in 0..3 {
                        for r in 0..3 {
                            let local = native.physical[18 * (bi + p) + bj + r]
                                + native.drilling[18 * (bi + p) + bj + r];
                            value += element.geometry.frame[p][i % 3]
                                * local
                                * element.geometry.frame[r][j % 3];
                        }
                    }
                    if i % 6 >= 3 {
                        value /= self.coordinate_length_m;
                    }
                    if j % 6 >= 3 {
                        value /= self.coordinate_length_m;
                    }
                    if !value.is_finite() {
                        return Err(invalid("scaled shell reference inertia is unrepresentable"));
                    }
                    if value != 0.0 {
                        matrix.push(row, column, value);
                        if i != j {
                            matrix.push(column, row, value);
                        }
                    }
                }
            }
            facets.push(FacetInventory {
                bounds,
                densities_kg_m3,
                thickness_m: thickness,
                moments,
            });
        }
        let mass = ReferenceMass::from_matrix(matrix.assemble(), total_mass, ndofs, cx)?;
        checkpoint(cx)?;
        Ok((facets, mass))
    }

    /// Advance exactly one physical interval from a retained endpoint (or
    /// virgin zero displacement when `initial` is None). `load_path` must be
    /// `[1.0]`: loads and prescribed coordinates are actual endpoint values.
    /// Supports include full consistent-mass coupling. Free bodies need no
    /// artificial restraints. All trial returns use the same old history.
    ///
    /// The reference inventory must remain fixed across accepted steps. Timed
    /// material laws must be constructed for this same `dt_s`; the energetic
    /// wrapper additionally checks duration for its supported native laws.
    /// This first-order rule dissipates numerical kinetic energy. Small-strain
    /// and small-rotation admission remains unchanged; no finite rigid spin or
    /// contact/thermal evolution is implied by adding inertia.
    pub fn advance_backward_euler(
        &self,
        initial: Option<&SmallStrainShellSolution<L::State>>,
        velocity: &[[f64; 6]],
        reference: &PlasticShellReferenceInertia<'_>,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellDynamicStep<L::State>, SmallStrainTetError> {
        let elements = self.prepare(cx)?;
        self.advance_backward_euler_prepared(&elements, initial, velocity, reference, dt_s, cx)
    }

    pub(super) fn advance_backward_euler_prepared(
        &self,
        elements: &[Element],
        initial: Option<&SmallStrainShellSolution<L::State>>,
        velocity: &[[f64; 6]],
        reference: &PlasticShellReferenceInertia<'_>,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellDynamicStep<L::State>, SmallStrainTetError> {
        checkpoint(cx)?;
        if self.load_path != [1.0]
            || velocity.len() != self.mesh.nodes.len()
            || initial.is_some_and(|old| old.coordinates.len() != self.mesh.nodes.len())
        {
            return Err(invalid(
                "shell dynamics requires one factor-one endpoint and complete old motion",
            ));
        }
        reference.admit(self, cx)?;
        let tangent_entries = elements
            .len()
            .checked_mul(324)
            .and_then(|count| count.checked_add(reference.mass.matrix().nnz()))
            .ok_or_else(|| invalid("dynamic shell tangent count overflow"))?;
        if tangent_entries > self.budget.maximum_tangent_entries {
            return Err(invalid(
                "material and inertial shell tangent exceeds allocation budget",
            ));
        }
        let length = self.coordinate_length_m;
        let old_q: Vec<_> = (0..6 * self.mesh.nodes.len())
            .map(|d| {
                initial.map_or(0.0, |old| old.coordinates[d / 6][d % 6])
                    * if d % 6 < 3 { 1.0 } else { length }
            })
            .collect();
        let old_v: Vec<_> = velocity
            .iter()
            .flat_map(|v| {
                v.iter()
                    .enumerate()
                    .map(move |(c, value)| value * if c < 3 { 1.0 } else { length })
            })
            .collect();
        let inertia = BackwardEulerInertia::new(&reference.mass, &old_q, &old_v, dt_s, cx)?;
        let solution = self.solve_prepared(
            elements,
            initial.map(|old| (old.coordinates.as_slice(), old.section_states.as_slice())),
            Some(&inertia),
            cx,
        )?;
        let new_q: Vec<_> = solution
            .coordinates
            .iter()
            .flat_map(|q| {
                q.iter()
                    .enumerate()
                    .map(move |(c, value)| value * if c < 3 { 1.0 } else { length })
            })
            .collect();
        let (scaled_velocity, total) = inertia.motion_flat(&new_q, cx)?;
        let next_velocity: Vec<[f64; 6]> = scaled_velocity
            .chunks_exact(6)
            .map(|v| core::array::from_fn(|c| v[c] / if c < 3 { 1.0 } else { length }))
            .collect();
        if next_velocity.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid("physical shell velocity is unrepresentable"));
        }
        let (physical, drilling_gauge) =
            reference.partition(elements, velocity, &next_velocity, cx)?;
        for (native, split) in [
            (
                total.kinetic_before_j,
                physical.kinetic_before_j + drilling_gauge.kinetic_before_j,
            ),
            (
                total.kinetic_after_j,
                physical.kinetic_after_j + drilling_gauge.kinetic_after_j,
            ),
            (
                total.kinetic_change_j,
                physical.kinetic_change_j + drilling_gauge.kinetic_change_j,
            ),
            (
                total.numerical_loss_j,
                physical.numerical_loss_j + drilling_gauge.numerical_loss_j,
            ),
            (
                total.inertial_work_j,
                physical.inertial_work_j + drilling_gauge.inertial_work_j,
            ),
        ] {
            let scale = total
                .kinetic_before_j
                .max(total.kinetic_after_j)
                .max(total.numerical_loss_j)
                .max(native.abs())
                .max(split.abs());
            if !split.is_finite() || (native - split).abs() > 4096.0 * f64::EPSILON * scale {
                return Err(invalid(
                    "shell physical/gauge partition misses original inertial work",
                ));
            }
        }
        checkpoint(cx)?;
        Ok(PlasticShellDynamicStep {
            solution,
            velocity: next_velocity,
            inertia: ShellInertialEnergy {
                total,
                physical,
                drilling_gauge,
            },
        })
    }
}

impl PlasticShellReferenceInertia<'_> {
    fn partition(
        &self,
        elements: &[Element],
        old: &[[f64; 6]],
        new: &[[f64; 6]],
        cx: &Cx<'_>,
    ) -> Result<(ShellKineticEnergy, ShellKineticEnergy), SmallStrainTetError> {
        let mut accounts = [ShellKineticEnergy::default(); 2];
        for (element, inventory) in elements.iter().zip(&self.facets) {
            checkpoint(cx)?;
            let local = |v: &[[f64; 6]]| -> [f64; 18] {
                core::array::from_fn(|d| {
                    let node = element.nodes[d / 6];
                    let offset = 3 * ((d % 6) / 3);
                    (0..3)
                        .map(|axis| element.geometry.frame[d % 3][axis] * v[node][offset + axis])
                        .sum()
                })
            };
            let before = local(old);
            let after = local(new);
            let (physical, drilling) = reference_kinetic_increment(
                &element.geometry,
                [inventory.thickness_m; 3],
                [inventory.moments; 3],
                &before,
                &after,
            )
            .map_err(|e| invalid(e.to_string()))?;
            for (account, value) in accounts.iter_mut().zip([physical, drilling]) {
                account.kinetic_before_j += value.kinetic_before_j;
                account.kinetic_after_j += value.kinetic_after_j;
                account.kinetic_change_j += value.kinetic_change_j;
                account.numerical_loss_j += value.numerical_loss_j;
                account.inertial_work_j += value.inertial_work_j;
            }
        }
        for account in &accounts {
            let scale = account
                .kinetic_before_j
                .max(account.kinetic_after_j)
                .max(account.numerical_loss_j)
                .max(account.inertial_work_j.abs());
            if [
                account.kinetic_before_j,
                account.kinetic_after_j,
                account.numerical_loss_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
                || !account.kinetic_change_j.is_finite()
                || !account.inertial_work_j.is_finite()
                || (account.inertial_work_j - account.kinetic_change_j - account.numerical_loss_j)
                    .abs()
                    > 4096.0 * f64::EPSILON * scale
            {
                return Err(invalid(
                    "shell kinetic partition is not a finite nonnegative storage/loss account",
                ));
            }
        }
        Ok((accounts[0], accounts[1]))
    }
}
