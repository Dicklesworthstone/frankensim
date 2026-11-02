//! Objective discrete-gradient elasticity with consistent quadratic inertia.
//!
//! Native material energies supply two-state stresses whose work equals their
//! physical storage change and whose material moment is symmetric. Centered
//! reference inertia supplies the physical kinetic change without numerical
//! kinetic loss. Accepted deformation, velocity, energy and time are immutable;
//! a candidate is published only after local, global and accumulated work gates.
//!
//! The true algorithmic tangent is generally NONSYMMETRIC. Native Newton/FGMRES
//! applies it unchanged; a symmetric sparse LDL factor is only a preconditioner.
//! Its inertia is not advertised as a stability certificate for the true system.
//! Reference geometry and signed consistent P2 mass are prepared once and shared.
//!
//! This is a second-order conservative nonlinear time discretization, with
//! explicit constitutive roundoff/work admission. It is not frequency matched:
//! the instrument-render linear splitting and its phase gate remain separate.
//! Loads are held dead forces and supports are absolute endpoint displacements.
//! Contact, follower forces, heat and irreversible material history have their
//! separate owners. Sampled geometric admission does not certify injectivity.

use std::cell::RefCell;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::discrete_gradient::bulk::{DiscreteBulkPoint, DiscreteBulkSettings};
use fs_material::hyper::discrete_gradient::{
    DiscreteGradientPoint, DiscreteGradientSettings, EnergyPart,
};
use fs_solver::nonlinear::{NewtonKrylovState, NewtonReport, NonlinearProblem, SolverRunProgress};
use fs_sparse::{Coo, Csr, DirectOrdering, LdltFactor, LdltOptions, SymbolicLdlt};

use super::element::{MixedPoint, PreparedMixedElement, PreparedMixedMesh};
use super::{MixedLawField, MixedUpError, MixedUpSettings};
use crate::inertia::{MidpointInertia, ReferenceMass};
use crate::linear3::TetAssemblyBudget;

pub use crate::inertia::InertialEnergy;

mod assembly;

/// Explicit material/spatial formulation; no inferred bulk penalty or defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ConservativeFormulation {
    /// Original full elastic energy on P2 displacement. Nearly incompressible
    /// displacement-only locking is not removed by the time discretization.
    FullMaterial,
    /// P2 displacement / continuous P1 pressure enforcing endpoint weak J=1.
    /// Pressure is the interval's algebraic work multiplier, not an endpoint
    /// Cauchy pressure. Native remainder energy is independent of bulk modulus.
    /// Each weak pressure moment is bounded relative to total reference volume.
    ExactIncompressible {
        maximum_relative_volume_residual: f64,
        /// Initial weak volume-rate allowance / total reference volume (1/s).
        /// Admits DAE-compatible input; later velocities satisfy the centered
        /// position update, not an extra endpoint tangent-velocity equation.
        maximum_initial_relative_volume_rate_per_s: f64,
        /// Dimensionless absolute local allowance for Gd:DeltaF = DeltaJ.
        maximum_point_volume_work_defect: f64,
    },
    /// Original finite-compliance bulk energy on retained quadrature volumes.
    /// P1 interval pressure is the original U secant, not an endpoint pressure.
    /// Both geometric and retained volumes must stay in the declared bracket.
    OriginalMaterial {
        bulk: DiscreteBulkSettings,
        maximum_relative_volume_residual: f64,
        maximum_point_volume_work_defect: f64,
    },
}

impl ConservativeFormulation {
    fn volume_controls(self) -> Option<(f64, f64)> {
        match self {
            Self::FullMaterial => None,
            Self::ExactIncompressible {
                maximum_relative_volume_residual,
                maximum_point_volume_work_defect,
                ..
            }
            | Self::OriginalMaterial {
                maximum_relative_volume_residual,
                maximum_point_volume_work_defect,
                ..
            } => Some((
                maximum_relative_volume_residual,
                maximum_point_volume_work_defect,
            )),
        }
    }

    fn bulk_settings(self) -> Option<DiscreteBulkSettings> {
        match self {
            Self::OriginalMaterial { bulk, .. } => Some(bulk),
            _ => None,
        }
    }

    fn admit_geometric_volume(self, volume: f64) -> Result<(), ConservativeTetError> {
        if let Some(settings) = self.bulk_settings() {
            if !volume.is_finite()
                || volume < settings.minimum_volume_ratio
                || volume > settings.maximum_volume_ratio
            {
                return Err(invalid(
                    "geometric volume is outside the original-bulk bracket",
                ));
            }
        }
        Ok(())
    }
}

/// Original reference inputs and explicit preparation/constitutive budgets.
#[derive(Clone, Copy, Debug)]
pub struct ConservativeTetDefinition<'a> {
    pub nodes_m: &'a [[f64; 3]],
    pub elements: &'a [[usize; 10]],
    /// Original laws, uniform, per element or per reference quadrature point.
    pub material: MixedLawField<'a>,
    pub reference_density_kg_m3: &'a [f64],
    pub formulation: ConservativeFormulation,
    pub budget: TetAssemblyBudget,
    pub quadrature_points_per_dir: usize,
    pub maximum_quadrature_points: usize,
    pub maximum_mass_assembly_bytes: usize,
    pub material_settings: DiscreteGradientSettings,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ConservativeTetError {
    Native(MixedUpError),
    InitialVolumeRateGate {
        pressure_node: usize,
        rate_m3_s: f64,
        maximum_abs_m3_s: f64,
    },
    VolumeGate {
        pressure_node: usize,
        residual_m3: f64,
        maximum_abs_m3: f64,
    },
    EnergyGate {
        quantity: &'static str,
        value_j: f64,
        maximum_abs_j: f64,
    },
}

impl From<MixedUpError> for ConservativeTetError {
    fn from(value: MixedUpError) -> Self {
        Self::Native(value)
    }
}
impl core::fmt::Display for ConservativeTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(error) => error.fmt(f),
            Self::InitialVolumeRateGate {
                pressure_node,
                rate_m3_s,
                maximum_abs_m3_s,
            } => write!(
                f,
                "initial conservative volume rate at node {pressure_node}: {rate_m3_s:e} m^3/s exceeds {maximum_abs_m3_s:e} m^3/s"
            ),
            Self::VolumeGate {
                pressure_node,
                residual_m3,
                maximum_abs_m3,
            } => write!(
                f,
                "conservative volume moment at node {pressure_node}: {residual_m3:e} m^3 exceeds {maximum_abs_m3:e} m^3"
            ),
            Self::EnergyGate {
                quantity,
                value_j,
                maximum_abs_j,
            } => write!(
                f,
                "conservative solid {quantity}: {value_j:e} J exceeds {maximum_abs_j:e} J"
            ),
        }
    }
}
impl core::error::Error for ConservativeTetError {}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ConservativeTetError> {
    super::checkpoint(cx).map_err(Into::into)
}
fn invalid(what: impl Into<String>) -> ConservativeTetError {
    super::invalid(what).into()
}
fn inertia_error(error: crate::inertia::InertiaError) -> ConservativeTetError {
    super::inertia_error(error).into()
}
fn material_error(element: usize, point: usize, error: MaterialError) -> ConservativeTetError {
    MixedUpError::Material {
        element,
        point,
        error,
    }
    .into()
}
fn gate(
    quantity: &'static str,
    value_j: f64,
    maximum_abs_j: f64,
) -> Result<(), ConservativeTetError> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        Err(ConservativeTetError::EnergyGate {
            quantity,
            value_j,
            maximum_abs_j,
        })
    } else {
        Ok(())
    }
}

fn deformation(element: &PreparedMixedElement, point: &MixedPoint, u: &[[f64; 3]]) -> [f64; 9] {
    let mut f = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    for a in 1..10 {
        for i in 0..3 {
            for j in 0..3 {
                f[3 * i + j] += (u[element.connectivity[a]][i] - u[element.connectivity[0]][i])
                    * point.gradients[a][j];
            }
        }
    }
    f
}

fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

fn law_at<'a>(
    field: MixedLawField<'a>,
    element: usize,
    point: usize,
    points: usize,
) -> &'a Hyperelastic {
    match field {
        MixedLawField::Uniform(law) => law,
        MixedLawField::PerElement(laws) => &laws[element],
        MixedLawField::PerQuadrature(laws) => &laws[element * points + point],
    }
}

impl<'a> ConservativeTetDefinition<'a> {
    /// Bind actual initial displacement and velocity, without an equilibrium or
    /// zero-strain assumption. Density is original mass per reference volume.
    pub fn prepare(
        &self,
        displacement_m: &[[f64; 3]],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<ConservativeTetSolid<'a>, ConservativeTetError> {
        self.prepare_with_pressure(displacement_m, velocity_mps, &[], cx)
    }

    /// Supply explicit P1 initial pressure in sorted corner-node order. An empty
    /// slice means zero. OriginalMaterial uses the physical ENDPOINT pressure
    /// to initialize each retained volume through U'(Jbar)=p, then admits the
    /// weak geometric-volume constraint. ExactIncompressible uses it only as a
    /// Newton predictor. Initial displacement/velocity are never projected.
    pub fn prepare_with_pressure(
        &self,
        displacement_m: &[[f64; 3]],
        velocity_mps: &[[f64; 3]],
        initial_pressure_pa: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ConservativeTetSolid<'a>, ConservativeTetError> {
        checkpoint(cx)?;
        if self
            .formulation
            .volume_controls()
            .is_some_and(|(relative, local)| {
                [relative, local].iter().any(|v| !v.is_finite() || *v < 0.0)
            })
        {
            return Err(invalid(
                "finite nonnegative mixed-volume allowances required",
            ));
        }
        if displacement_m.len() != self.nodes_m.len()
            || velocity_mps.len() != self.nodes_m.len()
            || displacement_m
                .iter()
                .chain(velocity_mps)
                .flatten()
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "finite full initial displacement and velocity required",
            ));
        }
        let mesh = Arc::new(
            PreparedMixedMesh::new(
                self.nodes_m,
                self.elements,
                self.budget,
                self.quadrature_points_per_dir,
                self.maximum_quadrature_points,
                cx,
            )
            .map_err(MixedUpError::Geometry)?,
        );
        let points = mesh.elements.len() * mesh.points_per_element;
        let count = match self.material {
            MixedLawField::Uniform(_) => None,
            MixedLawField::PerElement(laws) => Some((laws.len(), mesh.elements.len())),
            MixedLawField::PerQuadrature(laws) => Some((laws.len(), points)),
        };
        if count.is_some_and(|(got, expected)| got != expected) {
            return Err(invalid(
                "elastic law field must cover the reference quadrature",
            ));
        }
        let mass = Arc::new(
            ReferenceMass::from_tet10_quadrature(
                mesh.elements.iter().map(|element| {
                    (
                        element.connectivity,
                        element.points.iter().map(|q| (q.shapes, q.weight_m3)),
                    )
                }),
                self.reference_density_kg_m3,
                3 * mesh.nodes_count,
                self.maximum_mass_assembly_bytes,
                cx,
            )
            .map_err(inertia_error)?,
        );
        let np = if self.formulation.volume_controls().is_some() {
            mesh.pressure_nodes.len()
        } else {
            0
        };
        if (!initial_pressure_pa.is_empty() && initial_pressure_pa.len() != np)
            || initial_pressure_pa.iter().any(|p| !p.is_finite())
        {
            return Err(invalid("initial pressure must cover the pressure nodes"));
        }
        let algorithmic_pressure_pa = if initial_pressure_pa.is_empty() {
            vec![0.0; np]
        } else {
            initial_pressure_pa.to_vec()
        };
        let frozen =
            self.freeze_points(&mesh, displacement_m, None, &algorithmic_pressure_pa, cx)?;
        self.admit_volume(&mesh, &frozen.pressure_residual_m3)?;
        self.admit_initial_volume_rate(&mesh, &frozen.points, velocity_mps, cx)?;
        let momentum = mass
            .action(velocity_mps.as_flattened(), cx)
            .map_err(inertia_error)?;
        let kinetic_energy_j = 0.5
            * velocity_mps
                .as_flattened()
                .iter()
                .zip(momentum)
                .fold(0.0, |sum, (v, mv)| v.mul_add(mv, sum));
        let initial_total_energy_j = frozen.stored_energy_j + kinetic_energy_j;
        if !initial_total_energy_j.is_finite()
            || !kinetic_energy_j.is_finite()
            || kinetic_energy_j < 0.0
        {
            return Err(invalid("initial physical energy is unrepresentable"));
        }
        checkpoint(cx)?;
        Ok(ConservativeTetSolid {
            definition: *self,
            mesh,
            mass,
            points: frozen.points,
            bulk_points: frozen.bulk_points,
            bulk_volume_ratios: frozen.bulk_volume_ratios,
            displacement_m: displacement_m.to_vec(),
            velocity_mps: velocity_mps.to_vec(),
            algorithmic_pressure_pa,
            pressure_residual_m3: frozen.pressure_residual_m3,
            stored_energy_j: frozen.stored_energy_j,
            bulk_stored_energy_j: frozen.bulk_stored_energy_j,
            kinetic_energy_j,
            initial_total_energy_j,
            cumulative_work_j: 0.0,
            work_correction_j: 0.0,
            time_s: 0.0,
        })
    }

    fn freeze_points(
        &self,
        mesh: &PreparedMixedMesh,
        u: &[[f64; 3]],
        retained_bulk_volumes: Option<&[f64]>,
        initial_pressure_pa: &[f64],
        cx: &Cx<'_>,
    ) -> Result<FrozenState<'a>, ConservativeTetError> {
        let count = mesh.elements.len() * mesh.points_per_element;
        let mut points = Vec::with_capacity(count);
        let mut energy = 0.0;
        let mixed = self.formulation.volume_controls().is_some();
        let bulk_settings = self.formulation.bulk_settings();
        if retained_bulk_volumes
            .is_some_and(|volumes| bulk_settings.is_none() || volumes.len() != count)
        {
            return Err(invalid(
                "retained bulk volumes must cover the original quadrature",
            ));
        }
        let mut bulk_points = Vec::with_capacity(if bulk_settings.is_some() { count } else { 0 });
        let mut bulk_volume_ratios = Vec::with_capacity(bulk_points.capacity());
        let mut bulk_stored_energy_j = 0.0;
        let mut pressure_residual_m3 = vec![0.0; if mixed { mesh.pressure_nodes.len() } else { 0 }];
        for (e, element) in mesh.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let law = law_at(self.material, e, q, mesh.points_per_element);
                let f = deformation(element, point, u);
                let j = determinant(&f);
                self.formulation.admit_geometric_volume(j)?;
                let frozen = DiscreteGradientPoint::try_new(
                    law,
                    if mixed {
                        EnergyPart::Remainder
                    } else {
                        EnergyPart::Full
                    },
                    f,
                    self.material_settings,
                )
                .map_err(|error| material_error(e, q, error))?;
                energy += point.weight_m3 * frozen.energy_before_pa();
                let bulk_volume = if let Some(settings) = bulk_settings {
                    // A compatible initial endpoint pressure avoids unresolved
                    // volume components that would alternate under a quadratic
                    // secant. Later endpoints retain their solved volume;
                    // replacing it by det(F) would destroy physical bulk work.
                    let bulk = match retained_bulk_volumes {
                        Some(volumes) => DiscreteBulkPoint::try_new(
                            law,
                            volumes[e * mesh.points_per_element + q],
                            settings,
                        ),
                        None => {
                            let pressure = (0..4).fold(0.0, |s, a| {
                                point.pressure_shapes[a]
                                    .mul_add(initial_pressure_pa[element.pressure_indices[a]], s)
                            });
                            DiscreteBulkPoint::from_endpoint_pressure(law, pressure, settings)
                        }
                    }
                    .map_err(|error| material_error(e, q, error))?;
                    let ratio = bulk.accepted_volume_ratio();
                    bulk_stored_energy_j += point.weight_m3 * bulk.energy_before_pa();
                    bulk_points.push(bulk);
                    bulk_volume_ratios.push(ratio);
                    ratio
                } else {
                    1.0
                };
                if mixed {
                    let residual = j - bulk_volume;
                    for a in 0..4 {
                        pressure_residual_m3[element.pressure_indices[a]] +=
                            point.weight_m3 * point.pressure_shapes[a] * residual;
                    }
                }
                points.push(frozen);
            }
        }
        let stored_energy_j = energy + bulk_stored_energy_j;
        if !energy.is_finite() || !bulk_stored_energy_j.is_finite() || !stored_energy_j.is_finite()
        {
            return Err(invalid("physical elastic storage is unrepresentable"));
        }
        Ok(FrozenState {
            points,
            bulk_points,
            bulk_volume_ratios,
            stored_energy_j,
            bulk_stored_energy_j,
            pressure_residual_m3,
        })
    }

    fn admit_volume(
        &self,
        mesh: &PreparedMixedMesh,
        residual: &[f64],
    ) -> Result<(), ConservativeTetError> {
        if let Some((relative, _)) = self.formulation.volume_controls() {
            let maximum_abs_m3 = relative * mesh.total_volume_m3;
            if !maximum_abs_m3.is_finite() {
                return Err(invalid("mixed-volume allowance is unrepresentable"));
            }
            for (a, &residual_m3) in residual.iter().enumerate() {
                if !residual_m3.is_finite() || residual_m3.abs() > maximum_abs_m3 {
                    return Err(ConservativeTetError::VolumeGate {
                        pressure_node: mesh.pressure_nodes[a],
                        residual_m3,
                        maximum_abs_m3,
                    });
                }
            }
        }
        Ok(())
    }

    fn admit_initial_volume_rate(
        &self,
        mesh: &PreparedMixedMesh,
        points: &[DiscreteGradientPoint<'a>],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<(), ConservativeTetError> {
        let ConservativeFormulation::ExactIncompressible {
            maximum_initial_relative_volume_rate_per_s: relative_rate,
            maximum_point_volume_work_defect: local,
            ..
        } = self.formulation
        else {
            return Ok(());
        };
        let maximum_abs_m3_s = relative_rate * mesh.total_volume_m3;
        if !maximum_abs_m3_s.is_finite() || relative_rate < 0.0 {
            return Err(invalid(
                "initial exact-volume rate allowance is unrepresentable",
            ));
        }
        let mut rates = vec![0.0; mesh.pressure_nodes.len()];
        for (e, element) in mesh.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let frozen = &points[e * mesh.points_per_element + q];
                // Reuse the native volume derivative owner at the admitted
                // initial F. No bulk modulus or penalty enters this constraint.
                let volume = frozen
                    .volume_state(frozen.accepted_deformation(), local)
                    .map_err(|error| material_error(e, q, error))?;
                let mut f_rate = [0.0; 9];
                for a in 1..10 {
                    for i in 0..3 {
                        for j in 0..3 {
                            f_rate[3 * i + j] += (velocity_mps[element.connectivity[a]][i]
                                - velocity_mps[element.connectivity[0]][i])
                                * point.gradients[a][j];
                        }
                    }
                }
                let rate = volume
                    .endpoint_cofactor
                    .iter()
                    .zip(f_rate)
                    .fold(0.0, |s, (cofactor, velocity)| cofactor.mul_add(velocity, s));
                for a in 0..4 {
                    rates[element.pressure_indices[a]] +=
                        point.weight_m3 * point.pressure_shapes[a] * rate;
                }
            }
        }
        for (a, rate_m3_s) in rates.into_iter().enumerate() {
            if !rate_m3_s.is_finite() || rate_m3_s.abs() > maximum_abs_m3_s {
                return Err(ConservativeTetError::InitialVolumeRateGate {
                    pressure_node: mesh.pressure_nodes[a],
                    rate_m3_s,
                    maximum_abs_m3_s,
                });
            }
        }
        Ok(())
    }
}

struct FrozenState<'a> {
    points: Vec<DiscreteGradientPoint<'a>>,
    bulk_points: Vec<DiscreteBulkPoint>,
    bulk_volume_ratios: Vec<f64>,
    stored_energy_j: f64,
    bulk_stored_energy_j: f64,
    pressure_residual_m3: Vec<f64>,
}

/// Immutable accepted finite-strain motion and its original reference inventory.
#[derive(Debug)]
pub struct ConservativeTetSolid<'a> {
    definition: ConservativeTetDefinition<'a>,
    mesh: Arc<PreparedMixedMesh>,
    mass: Arc<ReferenceMass>,
    points: Vec<DiscreteGradientPoint<'a>>,
    bulk_points: Vec<DiscreteBulkPoint>,
    bulk_volume_ratios: Vec<f64>,
    displacement_m: Vec<[f64; 3]>,
    velocity_mps: Vec<[f64; 3]>,
    algorithmic_pressure_pa: Vec<f64>,
    pressure_residual_m3: Vec<f64>,
    stored_energy_j: f64,
    bulk_stored_energy_j: f64,
    kinetic_energy_j: f64,
    initial_total_energy_j: f64,
    cumulative_work_j: f64,
    work_correction_j: f64,
    time_s: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ConservativeTetEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Original U at retained mixed quadrature volumes; zero in other modes.
    pub bulk_stored_before_j: f64,
    pub bulk_stored_after_j: f64,
    /// Algorithmic bulk work p*(Jbar_after-Jbar_before).
    pub bulk_work_j: f64,
    pub bulk_constitutive_work_error_j: f64,
    /// Full-material work, or remainder work when pressure owns the bulk.
    pub material_work_j: f64,
    pub pressure_work_j: f64,
    pub dead_load_work_j: f64,
    pub support_work_j: f64,
    pub constitutive_work_error_j: f64,
    /// Actual pressure work minus bulk work and p dotted with the weak-volume change.
    pub constraint_work_error_j: f64,
    pub equilibrium_work_error_j: f64,
    /// External/support work minus physical stored and kinetic changes.
    pub balance_error_j: f64,
    /// Accumulated external/support work minus physical total-energy change.
    pub cumulative_balance_error_j: f64,
}

#[derive(Debug)]
pub struct ConservativeTetStep<'a> {
    pub solid: ConservativeTetSolid<'a>,
    /// None only for fully prescribed FullMaterial motion. Mixed formulations
    /// retain algebraic pressure equations and admit their coupled factorization.
    pub newton: Option<NewtonReport>,
    pub reactions_n: Vec<(usize, f64)>,
    pub energy: ConservativeTetEnergy,
    pub inertia: InertialEnergy,
    pub min_det_f: f64,
    pub pressure_residual_m3: Vec<f64>,
}

struct Evaluation {
    forces_n: Vec<f64>,
    tangent: Option<Csr>,
    stored_j: f64,
    bulk_stored_j: f64,
    bulk_work_j: f64,
    bulk_volume_ratios: Vec<f64>,
    material_work_j: f64,
    pressure_work_j: f64,
    pressure_residual_m3: Vec<f64>,
    min_det_f: f64,
}

impl<'a> ConservativeTetSolid<'a> {
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.displacement_m
    }
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }
    /// Sorted reference corner-node IDs, empty for FullMaterial motion.
    pub fn pressure_nodes(&self) -> &[usize] {
        if self.algorithmic_pressure_pa.is_empty() {
            &[]
        } else {
            &self.mesh.pressure_nodes
        }
    }
    /// Algebraic multipliers for the last accepted interval. At time zero this
    /// returns the supplied initial pressure (zero with `prepare`).
    /// These are not endpoint Cauchy pressures and carry no inertia or storage.
    pub fn algorithmic_pressure_pa(&self) -> &[f64] {
        &self.algorithmic_pressure_pa
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
    /// Original finite-bulk volumes in reference element/quadrature order.
    /// Empty for FullMaterial and ExactIncompressible.
    pub fn bulk_volume_ratios(&self) -> &[f64] {
        &self.bulk_volume_ratios
    }
    pub fn bulk_stored_energy_j(&self) -> f64 {
        self.bulk_stored_energy_j
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.kinetic_energy_j
    }
    pub fn total_energy_j(&self) -> f64 {
        self.stored_energy_j + self.kinetic_energy_j
    }
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        self.definition.reference_density_kg_m3
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Advance one positive interval. Explicit work allowance gates material,
    /// momentum, total and accumulated balance. No accepted field is mutated.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: MixedUpSettings,
        maximum_work_error_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ConservativeTetStep<'a>, ConservativeTetError> {
        checkpoint(cx)?;
        let time_s = self.time_s + dt_s;
        if !maximum_work_error_j.is_finite()
            || maximum_work_error_j < 0.0
            || !time_s.is_finite()
            || time_s <= self.time_s
            || settings.load_steps != 1
            || settings.maximum_newton_iterations == 0
            || settings.maximum_factor_bytes == 0
            || ![settings.length_scale_m, settings.stress_scale_pa]
                .iter()
                .all(|s| s.is_finite() && *s > 0.0)
            || body_force_n_m3.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "one positive finite interval, explicit scales and finite work allowance required",
            ));
        }
        settings
            .newton
            .validate()
            .map_err(MixedUpError::NewtonConfiguration)?;
        // Actual/symmetric-preconditioning COO, sorting, CSR, vectors. The
        // independent native factor budget bounds numeric factor and scratch.
        let np = self.algorithmic_pressure_pa.len();
        let local_dofs = if np == 0 { 30 } else { 34 };
        let estimated = self
            .mesh
            .elements
            .len()
            .checked_mul((local_dofs * local_dofs + 3 * 10 * 10) * 256)
            .and_then(|n| {
                self.mesh
                    .nodes_count
                    .checked_add(np)
                    .and_then(|n| n.checked_mul(1024))
                    .and_then(|v| n.checked_add(v))
            })
            .ok_or_else(|| invalid("conservative assembly size overflow"))?;
        if estimated > settings.maximum_assembly_bytes {
            return Err(invalid(format!(
                "conservative assembly needs at most {estimated} bytes; cap {}",
                settings.maximum_assembly_bytes
            )));
        }
        let inertia = MidpointInertia::new(
            &self.mass,
            self.displacement_m.as_flattened(),
            self.velocity_mps.as_flattened(),
            dt_s,
            cx,
        )
        .map_err(inertia_error)?;
        let ndofs = self.mesh.nodes_count * 3;
        let mut u = self.displacement_m.clone();
        let mut pressure_pa = self.algorithmic_pressure_pa.clone();
        let mut constrained = vec![false; ndofs];
        for &(d, value) in prescribed_m {
            if d >= ndofs || !value.is_finite() || constrained[d] {
                return Err(invalid(
                    "invalid or duplicate prescribed conservative displacement",
                ));
            }
            constrained[d] = true;
            u[d / 3][d % 3] = value;
        }
        let free: Vec<_> = (0..ndofs)
            .filter(|&d| !constrained[d])
            .chain(ndofs..ndofs + np)
            .collect();
        if free.len() > self.definition.budget.maximum_free_dofs {
            return Err(invalid(
                "conservative free displacement/pressure budget exceeded",
            ));
        }
        let mut free_index = vec![None; ndofs + np];
        for (i, &d) in free.iter().enumerate() {
            free_index[d] = Some(i);
        }
        let mut external_n = vec![0.0; ndofs];
        for &(d, value) in nodal_forces_n {
            if d >= ndofs || !value.is_finite() {
                return Err(invalid("invalid conservative nodal force"));
            }
            external_n[d] += value;
        }
        for element in &self.mesh.elements {
            checkpoint(cx)?;
            for point in &element.points {
                for a in 0..10 {
                    for i in 0..3 {
                        external_n[3 * element.connectivity[a] + i] +=
                            point.weight_m3 * point.shapes[a] * body_force_n_m3[i];
                    }
                }
            }
        }
        let work_scale_j = settings.stress_scale_pa * self.mesh.total_volume_m3;
        if !work_scale_j.is_finite()
            || work_scale_j <= 0.0
            || external_n.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "conservative load or work scale is unrepresentable",
            ));
        }
        let mut report = None;
        if !free.is_empty() {
            let endpoint = Endpoint {
                solid: self,
                base_u: &u,
                external_n: &external_n,
                free: &free,
                free_index: &free_index,
                inertia: &inertia,
                settings,
                work_scale_j,
                cx,
                symbolic: RefCell::new(None),
                linearization: RefCell::new(None),
                fatal_error: RefCell::new(None),
                trial_error: RefCell::new(None),
            };
            let initial = free
                .iter()
                .map(|&d| {
                    if d < ndofs {
                        u[d / 3][d % 3] / settings.length_scale_m
                    } else {
                        pressure_pa[d - ndofs] / settings.stress_scale_pa
                    }
                })
                .collect();
            let mut newton =
                NewtonKrylovState::new(&endpoint, initial, settings.newton).map_err(|error| {
                    endpoint
                        .trial_error
                        .borrow_mut()
                        .take()
                        .unwrap_or_else(|| MixedUpError::NewtonConfiguration(error).into())
                })?;
            let target = settings
                .newton
                .absolute_tolerance
                .max(settings.newton.relative_tolerance * newton.residual_norm());
            let result =
                newton.run_controlled(&endpoint, settings.maximum_newton_iterations, &mut || {
                    cx.checkpoint().is_err()
                });
            checkpoint(cx)?;
            if let Some(error) = endpoint.fatal_error.borrow_mut().take() {
                return Err(error);
            }
            if result.progress != SolverRunProgress::Complete || !result.report.converged {
                return Err(MixedUpError::NewtonStalled {
                    load_factor: 1.0,
                    report: result.report,
                }
                .into());
            }
            // A solved initial residual does not establish an admissible mixed
            // factorization. Exact constraints require independent rows; finite
            // bulk retains its original nonzero pressure compliance instead.
            if np > 0 {
                endpoint.ensure_linearization(&newton.x)?;
            }
            let (next_u, next_p) = endpoint.physical_state(&newton.x)?;
            let final_state = self.evaluate(&next_u, &next_p, false, &inertia, cx)?;
            let residual: Vec<_> = free
                .iter()
                .map(|&d| endpoint.scaled_residual(&final_state, d))
                .collect();
            if fs_solver::norm2(&residual) > target {
                return Err(MixedUpError::NewtonStalled {
                    load_factor: 1.0,
                    report: result.report,
                }
                .into());
            }
            report = Some(result.report);
            u = next_u;
            pressure_pa = next_p;
        }
        let evaluation = self.evaluate(&u, &pressure_pa, false, &inertia, cx)?;
        self.definition
            .admit_volume(&self.mesh, &evaluation.pressure_residual_m3)?;
        let reactions_n: Vec<_> = prescribed_m
            .iter()
            .map(|&(d, _)| (d, evaluation.forces_n[d] - external_n[d]))
            .collect();
        let (velocity_mps, inertia_energy) = inertia
            .motion(u.as_flattened(), cx)
            .map_err(inertia_error)?;
        let du: Vec<_> = u
            .as_flattened()
            .iter()
            .zip(self.displacement_m.as_flattened())
            .map(|(a, b)| a - b)
            .collect();
        let dead_load_work_j = external_n
            .iter()
            .zip(&du)
            .fold(0.0, |s, (f, u)| f.mul_add(*u, s));
        let support_work_j = reactions_n
            .iter()
            .fold(0.0, |s, (d, f)| f.mul_add(du[*d], s));
        let work_j = dead_load_work_j + support_work_j;
        let stored_change_j = evaluation.stored_j - self.stored_energy_j;
        let bulk_change_j = evaluation.bulk_stored_j - self.bulk_stored_energy_j;
        let increment = work_j - self.work_correction_j;
        let cumulative_work_j = self.cumulative_work_j + increment;
        let work_correction_j = (cumulative_work_j - self.cumulative_work_j) - increment;
        let algebraic_pressure_work_j = pressure_pa
            .iter()
            .zip(&evaluation.pressure_residual_m3)
            .zip(&self.pressure_residual_m3)
            .fold(0.0, |s, ((p, after), before)| p.mul_add(after - before, s));
        let energy = ConservativeTetEnergy {
            stored_before_j: self.stored_energy_j,
            stored_after_j: evaluation.stored_j,
            bulk_stored_before_j: self.bulk_stored_energy_j,
            bulk_stored_after_j: evaluation.bulk_stored_j,
            bulk_work_j: evaluation.bulk_work_j,
            bulk_constitutive_work_error_j: evaluation.bulk_work_j - bulk_change_j,
            material_work_j: evaluation.material_work_j,
            pressure_work_j: evaluation.pressure_work_j,
            dead_load_work_j,
            support_work_j,
            constitutive_work_error_j: evaluation.material_work_j
                - (stored_change_j - bulk_change_j),
            constraint_work_error_j: evaluation.pressure_work_j
                - evaluation.bulk_work_j
                - algebraic_pressure_work_j,
            equilibrium_work_error_j: work_j
                - evaluation.material_work_j
                - evaluation.pressure_work_j
                - inertia_energy.inertial_work_j,
            balance_error_j: work_j - stored_change_j - inertia_energy.kinetic_change_j,
            cumulative_balance_error_j: cumulative_work_j
                - (evaluation.stored_j + inertia_energy.kinetic_after_j
                    - self.initial_total_energy_j),
        };
        for (quantity, value) in [
            ("constitutive work", energy.constitutive_work_error_j),
            ("original bulk work", energy.bulk_constitutive_work_error_j),
            ("constraint work identity", energy.constraint_work_error_j),
            (
                "pressure storage work",
                energy.pressure_work_j - bulk_change_j,
            ),
            ("momentum work", energy.equilibrium_work_error_j),
            ("kinetic work", inertia_energy.balance_error_j),
            ("total work", energy.balance_error_j),
            ("accumulated total work", energy.cumulative_balance_error_j),
        ] {
            gate(quantity, value, maximum_work_error_j)?;
        }
        let frozen = self.definition.freeze_points(
            &self.mesh,
            &u,
            self.definition
                .formulation
                .bulk_settings()
                .map(|_| evaluation.bulk_volume_ratios.as_slice()),
            &pressure_pa,
            cx,
        )?;
        self.definition
            .admit_volume(&self.mesh, &frozen.pressure_residual_m3)?;
        gate(
            "accepted physical storage",
            frozen.stored_energy_j - evaluation.stored_j,
            maximum_work_error_j,
        )?;
        gate(
            "accepted original bulk storage",
            frozen.bulk_stored_energy_j - evaluation.bulk_stored_j,
            maximum_work_error_j,
        )?;
        checkpoint(cx)?;
        Ok(ConservativeTetStep {
            solid: Self {
                definition: self.definition,
                mesh: Arc::clone(&self.mesh),
                mass: Arc::clone(&self.mass),
                points: frozen.points,
                bulk_points: frozen.bulk_points,
                bulk_volume_ratios: frozen.bulk_volume_ratios,
                displacement_m: u,
                velocity_mps,
                algorithmic_pressure_pa: pressure_pa,
                pressure_residual_m3: frozen.pressure_residual_m3,
                stored_energy_j: frozen.stored_energy_j,
                bulk_stored_energy_j: frozen.bulk_stored_energy_j,
                kinetic_energy_j: inertia_energy.kinetic_after_j,
                initial_total_energy_j: self.initial_total_energy_j,
                cumulative_work_j,
                work_correction_j,
                time_s,
            },
            newton: report,
            reactions_n,
            energy,
            inertia: inertia_energy,
            min_det_f: evaluation.min_det_f,
            pressure_residual_m3: evaluation.pressure_residual_m3,
        })
    }
}

struct Linearization {
    x: Vec<f64>,
    matrix: Csr,
    preconditioner: LdltFactor,
}

struct Endpoint<'p, 'a, 'c, 'scope> {
    solid: &'p ConservativeTetSolid<'a>,
    base_u: &'p [[f64; 3]],
    external_n: &'p [f64],
    free: &'p [usize],
    free_index: &'p [Option<usize>],
    inertia: &'p MidpointInertia<'p>,
    settings: MixedUpSettings,
    work_scale_j: f64,
    cx: &'c Cx<'scope>,
    symbolic: RefCell<Option<SymbolicLdlt>>,
    linearization: RefCell<Option<Linearization>>,
    fatal_error: RefCell<Option<ConservativeTetError>>,
    trial_error: RefCell<Option<ConservativeTetError>>,
}

impl Endpoint<'_, '_, '_, '_> {
    fn variable_scale(&self, dof: usize) -> f64 {
        if dof < 3 * self.solid.mesh.nodes_count {
            self.settings.length_scale_m
        } else {
            self.settings.stress_scale_pa
        }
    }

    fn scaled_residual(&self, state: &Evaluation, dof: usize) -> f64 {
        let nu = 3 * self.solid.mesh.nodes_count;
        let physical = if dof < nu {
            state.forces_n[dof] - self.external_n[dof]
        } else {
            state.pressure_residual_m3[dof - nu]
        };
        physical * self.variable_scale(dof) / self.work_scale_j
    }

    fn physical_state(&self, x: &[f64]) -> Result<(Vec<[f64; 3]>, Vec<f64>), ConservativeTetError> {
        if x.len() != self.free.len() || x.iter().any(|v| !v.is_finite()) {
            return Err(invalid("invalid native conservative Newton state"));
        }
        let mut u = self.base_u.to_vec();
        let nu = 3 * self.solid.mesh.nodes_count;
        let mut p = vec![0.0; self.solid.algorithmic_pressure_pa.len()];
        for (&d, &v) in self.free.iter().zip(x) {
            if d < nu {
                u[d / 3][d % 3] = v * self.settings.length_scale_m;
            } else {
                p[d - nu] = v * self.settings.stress_scale_pa;
            }
        }
        if u.iter().flatten().chain(&p).any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable physical conservative state"));
        }
        Ok((u, p))
    }

    fn ensure_linearization(&self, x: &[f64]) -> Result<(), ConservativeTetError> {
        checkpoint(self.cx)?;
        if self
            .linearization
            .borrow()
            .as_ref()
            .is_some_and(|linear| linear.x == x)
        {
            return Ok(());
        }
        *self.linearization.borrow_mut() = None;
        let (u, p) = self.physical_state(x)?;
        let full = self
            .solid
            .evaluate(&u, &p, true, self.inertia, self.cx)?
            .tangent
            .expect("requested tangent");
        let mut actual = Coo::new(self.free.len(), self.free.len());
        for (row, &d) in self.free.iter().enumerate() {
            if row % 128 == 0 {
                checkpoint(self.cx)?;
            }
            let (columns, values) = full.row(d);
            for (&col, &value) in columns.iter().zip(values) {
                if let Some(j) = self.free_index[col] {
                    let scale =
                        self.variable_scale(d) * self.variable_scale(col) / self.work_scale_j;
                    actual.push(row, j, value * scale);
                }
            }
        }
        let matrix = actual.assemble();
        let mut symmetric = Coo::new(self.free.len(), self.free.len());
        let nu = 3 * self.solid.mesh.nodes_count;
        for row in 0..matrix.nrows() {
            if row % 128 == 0 {
                checkpoint(self.cx)?;
            }
            let (columns, values) = matrix.row(row);
            for (&col, &value) in columns.iter().zip(values) {
                if !value.is_finite() {
                    return Err(invalid("scaled conservative tangent is unrepresentable"));
                }
                if self.free[row] < nu && self.free[col] < nu {
                    symmetric.push(row, col, 0.5 * value);
                    symmetric.push(col, row, 0.5 * value);
                } else if self.free[row] >= nu && self.free[col] < nu {
                    // Precondition with the ENDPOINT constraint derivative in
                    // both crosses, not the average of two distinct derivatives.
                    // A redundant endpoint constraint then remains singular.
                    symmetric.push(row, col, value);
                    symmetric.push(col, row, value);
                } else if self.free[row] >= nu && self.free[col] >= nu {
                    // Original finite bulk has a negative symmetric compliance
                    // block. Keeping it admits driven fully prescribed bodies;
                    // the exact-constraint branch still has zero compliance.
                    symmetric.push(row, col, value);
                }
            }
        }
        let preconditioning_matrix = symmetric.assemble();
        if self.symbolic.borrow().is_none() {
            *self.symbolic.borrow_mut() = Some(
                SymbolicLdlt::analyze(&preconditioning_matrix, DirectOrdering::Amd)
                    .map_err(MixedUpError::Direct)?,
            );
        }
        checkpoint(self.cx)?;
        let preconditioner = self
            .symbolic
            .borrow()
            .as_ref()
            .expect("analyzed pattern")
            .factor(
                &preconditioning_matrix,
                &LdltOptions {
                    max_factor_bytes: Some(self.settings.maximum_factor_bytes),
                    ..LdltOptions::default()
                },
            )
            .map_err(MixedUpError::Direct)?;
        checkpoint(self.cx)?;
        *self.linearization.borrow_mut() = Some(Linearization {
            x: x.to_vec(),
            matrix,
            preconditioner,
        });
        Ok(())
    }

    fn refuse_callback(&self, output: &mut [f64], error: ConservativeTetError) {
        output.fill(f64::NAN);
        *self.fatal_error.borrow_mut() = Some(error);
    }
}

impl NonlinearProblem for Endpoint<'_, '_, '_, '_> {
    fn dimension(&self) -> usize {
        self.free.len()
    }
    fn residual(&self, x: &[f64], residual: &mut [f64]) {
        let value = self
            .physical_state(x)
            .and_then(|(u, p)| self.solid.evaluate(&u, &p, false, self.inertia, self.cx));
        match value {
            Ok(state) => {
                for (r, &d) in residual.iter_mut().zip(self.free) {
                    *r = self.scaled_residual(&state, d);
                }
                *self.trial_error.borrow_mut() = None;
            }
            Err(error) => {
                residual.fill(f64::NAN);
                *self.trial_error.borrow_mut() = Some(error);
            }
        }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        match self.ensure_linearization(x) {
            Ok(()) => self
                .linearization
                .borrow()
                .as_ref()
                .expect("cached true tangent")
                .matrix
                .spmv(direction, output),
            Err(error) => self.refuse_callback(output, error),
        }
    }
    fn preconditioner_apply(
        &self,
        x: &[f64],
        _outer: usize,
        _inner: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        match self.ensure_linearization(x) {
            Ok(()) => output.copy_from_slice(
                &self
                    .linearization
                    .borrow()
                    .as_ref()
                    .expect("cached symmetric preconditioner")
                    .preconditioner
                    .solve(residual),
            ),
            Err(error) => self.refuse_callback(output, error),
        }
    }
}
