//! Causal material memory in the complete original MITC displacement space.
//!
//! The caller supplies the RELAXED native pencil and admitted excess resultant
//! branches. Their original passive factors and quadrature become Maxwell strain
//! rows of fs-phs's existing mechanical owner. No mode is selected, no material
//! time is averaged, and no frequency-sampled plane-stress law is substituted.

use super::relaxation::ShearShellRelaxationBranch;
use super::strain::{PreparedElementStrain, RealElementStrain, ShellStrainPoint, integration_rule};
use super::viscosity::{ShearDissipation, ShearShellViscosity};
use super::{MitcShellModel, ShellMass};
use crate::PlateError;
use fs_exec::Cx;
use fs_math::{dd::Dd, det};
use fs_phs::mechanical::{LinearMechanicalSystem, MAX_MECHANICAL_STATES};
use fs_phs::{PhsError, RelaxationBranch};

#[derive(Debug)]
pub enum MitcTransientError {
    Plate(PlateError),
    Mechanical(PhsError),
    Invalid(&'static str),
    Budget {
        quantity: &'static str,
        requested: u128,
        maximum: u128,
    },
    StrainDomain {
        element: usize,
        observed: f64,
        maximum: f64,
    },
}
impl From<PlateError> for MitcTransientError {
    fn from(value: PlateError) -> Self {
        Self::Plate(value)
    }
}
impl From<PhsError> for MitcTransientError {
    fn from(value: PhsError) -> Self {
        Self::Mechanical(value)
    }
}
impl core::fmt::Display for MitcTransientError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-PLATE-MITC-TRANSIENT: {self:?}")
    }
}
impl std::error::Error for MitcTransientError {}
type Result<T> = core::result::Result<T, MitcTransientError>;

/// Explicit applicability of the equivalent-single-layer resultant source.
/// The norm covers xx, yy, xy, xz and yz with tensor-shear weights. It does not
/// infer an unmeasured normal strain, failure criterion or thermal validity.
#[derive(Clone, Copy, Debug)]
pub struct MitcTransientDomain {
    /// Physical thickness faces measured from the ORIGINAL reference surface.
    pub faces_m: [f64; 2],
    /// Positive reduced-strain Frobenius ceiling, strictly below one.
    pub maximum_reduced_strain: f64,
}

/// Cold source preparation only. Native time propagators have separate budgets.
#[derive(Clone, Copy, Debug)]
pub struct MitcTransientBudget {
    pub maximum_elements: usize,
    pub maximum_coordinates: usize,
    /// Seven quadrature rows per nonzero constitutive factor channel, before
    /// any geometrically exact-zero row is omitted. No near-zero row is dropped.
    pub maximum_memory_channels: usize,
    /// Conservative numerical payload including the native dense source.
    pub maximum_payload_bytes: usize,
    /// Conservative cold scalar-work score, not a timing prediction.
    pub maximum_work: u64,
}

/// Deterministic meaning of one original material-memory coordinate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MitcMemoryChannel {
    pub element: usize,
    pub source_branch: usize,
    pub quadrature_point: usize,
    pub factor_channel: usize,
}

/// An admitted physical snapshot; velocities retain original consistent inertia.
#[derive(Clone, Debug)]
pub struct MitcTransientObservation {
    /// Global translations [m] and physical axial rotations [rad], six per node.
    pub displacement: Vec<f64>,
    pub velocity: Vec<f64>,
    /// Original kinetic, relaxed elastic and Maxwell spring storage [J].
    pub stored_energy_j: f64,
    /// Instantaneous Kelvin + Maxwell material heat rate [W], not a step loss.
    /// Native propagators report integrated heat and numerical balance separately.
    pub material_heat_w: f64,
}

fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| PlateError::Cancelled.into())
}
fn cap(quantity: &'static str, requested: u128, maximum: u128) -> Result<()> {
    if requested > maximum {
        Err(MitcTransientError::Budget {
            quantity,
            requested,
            maximum,
        })
    } else {
        Ok(())
    }
}
fn active(factor: &ShearDissipation) -> [bool; 12] {
    let mut active = [false; 12];
    for component in 0..8 {
        let mut e = [0.0; 8];
        e[component] = 1.0;
        let values = factor.channels([e[0], e[1], e[2]], [e[3], e[4], e[5]], [e[6], e[7]]);
        for (slot, value) in active.iter_mut().zip(values) {
            *slot |= value != 0.0;
        }
    }
    active
}

/// Complete original displacement coordinates and unmerged local material memory.
/// Borrowed geometry/sections remain the original elastic condensation owner.
pub struct MitcTransientSource<'a> {
    model: &'a MitcShellModel,
    mechanics: LinearMechanicalSystem,
    observation: Vec<PreparedElementStrain>,
    domains: Vec<MitcTransientDomain>,
    channels: Vec<MitcMemoryChannel>,
}

impl<'a> MitcTransientSource<'a> {
    /// Bind original consistent M and RELAXED K to independent facet viscosity
    /// and excess spectra in the facet engineering-strain axes. Include source
    /// orientations and shear correction exactly once when authoring resultants.
    /// Existing drilling and connection stiffness remain lossless in K.
    ///
    /// K must be the equilibrium/relaxed source. Supplying instantaneous K and
    /// its excess again would double-count stiffness; no relaxation is inferred
    /// from K. Arbitrary frequency-dependent laminate plane-stress elimination
    /// cannot be converted to these branches by sampling a harmonic response.
    ///
    /// The native pHS source is bounded dense (currently 256 complete states),
    /// with full force ports in original retained support-coordinate order.
    /// No nonzero mapped memory row, weak channel or physical coordinate is truncated.
    /// Exactly zero mapped rows have zero retained history; independent heat from
    /// a pre-existing unobservable history on such a row is outside this source.
    /// Native source factor admission is one bounded atomic operation; all
    /// spatial preparation polls the supplied Cx at facet/point/row boundaries.
    pub fn prepare(
        relaxed_model: &'a MitcShellModel,
        viscosity: &[ShearShellViscosity],
        branches: &[Vec<ShearShellRelaxationBranch>],
        domains: &[MitcTransientDomain],
        budget: MitcTransientBudget,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        let model = relaxed_model;
        let count = model.mesh.tris.len();
        let n = model.pencil.free;
        if viscosity.len() != count
            || branches.len() != count
            || domains.len() != count
            || n == 0
            || model.mass_model != ShellMass::Consistent
        {
            return Err(MitcTransientError::Invalid(
                "complete facet sources and nonempty original consistent-mass coordinates required",
            ));
        }
        cap("elements", count as u128, budget.maximum_elements as u128)?;
        cap("coordinates", n as u128, budget.maximum_coordinates as u128)?;
        let mut memory = 0_u128;
        let mut source_branches = 0_u128;
        for (spectrum, domain) in branches.iter().zip(domains) {
            poll(cx)?;
            if domain.faces_m.iter().any(|z| !z.is_finite())
                || domain.faces_m[0] >= domain.faces_m[1]
                || !domain.maximum_reduced_strain.is_finite()
                || domain.maximum_reduced_strain <= 0.0
                || domain.maximum_reduced_strain >= 1.0
            {
                return Err(MitcTransientError::Invalid(
                    "finite ordered source faces and a reduced-strain limit in (0,1) required",
                ));
            }
            for branch in spectrum {
                poll(cx)?;
                source_branches += 1;
                memory += 7 * active(branch.factor()).into_iter().filter(|v| *v).count() as u128;
            }
        }
        cap(
            "memory channels",
            memory,
            budget.maximum_memory_channels as u128,
        )?;
        let states = 2 * n as u128 + memory;
        cap(
            "native complete states",
            states,
            MAX_MECHANICAL_STATES as u128,
        )?;
        // Covers all simultaneous native source/factor matrices and retained
        // strain observers. Borrowed immutable geometry and stack arrays excluded.
        let bytes = 512 * states * states
            + 128 * (n as u128) * (n as u128)
            + 32 * memory * n as u128
            + count as u128 * (PreparedElementStrain::RETAINED_BYTES + 128) as u128
            + 96 * model.mesh.nodes.len() as u128;
        cap(
            "source payload bytes",
            bytes,
            budget.maximum_payload_bytes as u128,
        )?;
        let work =
            2_000_000 * count as u128 + 131_072 * source_branches + 512 * states * states * states;
        cap("source preparation work", work, budget.maximum_work as u128)?;

        let mut damping = vec![Dd::ZERO; n * n];
        let mut arms = Vec::with_capacity(memory as usize);
        let mut channels = Vec::with_capacity(memory as usize);
        let mut observation = Vec::with_capacity(count);
        let rule = integration_rule();
        if rule.len() != 7 {
            return Err(MitcTransientError::Invalid(
                "original seven-point strain rule required",
            ));
        }
        for (element, tri) in model.mesh.tris.iter().enumerate() {
            poll(cx)?;
            let field = PreparedElementStrain::new(
                &model.mesh,
                element,
                &model.directors,
                &model.sections[element],
                &mut || cx.checkpoint().is_err(),
            )?;
            let maps =
                core::array::from_fn::<_, 18, _>(|i| model.pencil.dof_map[6 * tri[i / 6] + i % 6]);
            let mut basis = Vec::with_capacity(18);
            for column in 0..18 {
                poll(cx)?;
                let mut nodal = [0.0; 18];
                let block = 6 * (column / 6) + 3 * ((column % 6) / 3);
                for axis in 0..3 {
                    nodal[block + axis] = model.frames[tri[column / 6]][column % 3][axis];
                }
                basis.push(field.real_field(&nodal, &mut || cx.checkpoint().is_err())?);
            }
            let material = viscosity[element].dissipation()?;
            let area = model.mesh.facet(element)?.area_m2;
            for (point, &(lambda, weight)) in rule.iter().enumerate() {
                poll(cx)?;
                let mut strain = [ShellStrainPoint {
                    membrane: [0.0; 3],
                    gradient_per_m: [0.0; 3],
                    transverse_shear: [0.0; 2],
                }; 18];
                for (out, field) in strain.iter_mut().zip(&basis) {
                    *out = field.at_barycentric(lambda)?;
                }
                let scale = det::sqrt(area) * det::sqrt(weight);
                let rows = factor_rows(&material, &strain, scale)?;
                for i in 0..18 {
                    poll(cx)?;
                    let Some(a) = maps[i] else {
                        continue;
                    };
                    for j in i..18 {
                        let Some(b) = maps[j] else {
                            continue;
                        };
                        let value = (0..12).fold(Dd::ZERO, |sum, k| {
                            sum + Dd::from_f64(rows[i][k]) * Dd::from_f64(rows[j][k])
                        });
                        damping[n * a + b] = damping[n * a + b] + value;
                        if a != b {
                            damping[n * b + a] = damping[n * b + a] + value;
                        }
                    }
                }
                for (source_branch, branch) in branches[element].iter().enumerate() {
                    poll(cx)?;
                    let rows = factor_rows(branch.factor(), &strain, scale)?;
                    for channel in 0..12 {
                        let mut projection = vec![0.0; n];
                        for column in 0..18 {
                            if let Some(index) = maps[column] {
                                projection[index] += rows[column][channel];
                            }
                        }
                        if projection.iter().any(|v| !v.is_finite()) {
                            return Err(MitcTransientError::Invalid(
                                "material strain row overflow",
                            ));
                        }
                        if projection.iter().all(|v| *v == 0.0) {
                            continue;
                        }
                        // The physical stiffness is already inside this ORIGINAL
                        // square-factor row. Native branch stiffness is unity.
                        arms.push(RelaxationBranch {
                            projection,
                            stiffness: 1.0,
                            relaxation_time_s: branch.relaxation_time_s(),
                        });
                        channels.push(MitcMemoryChannel {
                            element,
                            source_branch,
                            quadrature_point: point,
                            factor_channel: channel,
                        });
                    }
                }
            }
            observation.push(field);
        }
        let damping: Vec<_> = damping.into_iter().map(Dd::to_f64).collect();
        let mut mass = vec![0.0; n * n];
        let mut stiffness = vec![0.0; n * n];
        let mut forces = vec![0.0; n * n];
        for i in 0..n {
            poll(cx)?;
            forces[n * i + i] = 1.0;
            for (source, output) in [
                (&model.pencil.m, &mut mass),
                (&model.pencil.k, &mut stiffness),
            ] {
                let (columns, values) = source.row(i);
                for (&j, &value) in columns.iter().zip(values) {
                    output[n * i + j] = value;
                }
            }
        }
        poll(cx)?;
        let mechanics =
            LinearMechanicalSystem::new(&mass, &damping, &stiffness, &forces, n, n, arms)?;
        poll(cx)?;
        Ok(Self {
            model,
            mechanics,
            observation,
            domains: domains.to_vec(),
            channels,
        })
    }

    pub const fn model(&self) -> &'a MitcShellModel {
        self.model
    }
    /// Original native owner for prepared exact/discrete-gradient propagation
    /// and physical interconnection. Admit every candidate before publication.
    pub const fn mechanics(&self) -> &LinearMechanicalSystem {
        &self.mechanics
    }
    pub fn memory_channels(&self) -> &[MitcMemoryChannel] {
        &self.channels
    }

    /// Initialize q, v and ORIGINAL weighted viscous strains [sqrt(J)].
    /// Zero memory leaves the Maxwell springs loaded by the initial q.
    pub fn state(&self, q: &[f64], v: &[f64], memory: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>> {
        poll(cx)?;
        let state = self.mechanics.state(q, v, memory)?;
        self.admit_state(&state, cx)?;
        Ok(state)
    }
    /// Explicit fully relaxed initialization, without erasing a supplied history.
    pub fn relaxed_state(&self, q: &[f64], v: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>> {
        poll(cx)?;
        let state = self.mechanics.relaxed_state(q, v)?;
        self.admit_state(&state, cx)?;
        Ok(state)
    }
    /// Whole-facet amplitude admission at one actual time. This is not a bound
    /// between time samples; source-domain/time refinement remains necessary.
    pub fn admit_state(&self, state: &[f64], cx: &Cx<'_>) -> Result<()> {
        poll(cx)?;
        if state.len() != self.mechanics.system().state_dim()
            || state.iter().any(|v| !v.is_finite())
        {
            return Err(MitcTransientError::Invalid(
                "complete finite native state required",
            ));
        }
        let q = &state[..self.mechanics.coordinates()];
        for (element, observer) in self.observation.iter().enumerate() {
            poll(cx)?;
            let field = self.strain(element, observer, q, cx)?;
            let domain = self.domains[element];
            let observed =
                field.reduced_elastic_layer_bound(domain.faces_m, [[0.0; 3]; 2], &mut || {
                    cx.checkpoint().is_err()
                })?;
            if observed > domain.maximum_reduced_strain {
                return Err(MitcTransientError::StrainDomain {
                    element,
                    observed,
                    maximum: domain.maximum_reduced_strain,
                });
            }
        }
        poll(cx)
    }
    fn strain(
        &self,
        element: usize,
        observer: &PreparedElementStrain,
        q: &[f64],
        cx: &Cx<'_>,
    ) -> Result<RealElementStrain> {
        let tri = self.model.mesh.tris[element];
        let mut nodal = [0.0; 18];
        for (i, value) in nodal.iter_mut().enumerate() {
            let node = tri[i / 6];
            let block = 3 * ((i % 6) / 3);
            let mut sum = Dd::ZERO;
            for axis in 0..3 {
                if let Some(index) = self.model.pencil.dof_map[6 * node + block + axis] {
                    sum = sum
                        + Dd::from_f64(self.model.frames[node][axis][i % 3])
                            * Dd::from_f64(q[index]);
                }
            }
            *value = sum.to_f64();
        }
        Ok(observer.real_field(&nodal, &mut || cx.checkpoint().is_err())?)
    }

    /// Physical motion and instantaneous ORIGINAL material dissipation only.
    /// Drilling/connection elastic storage remains in H; it creates no heat.
    pub fn observe(&self, state: &[f64], cx: &Cx<'_>) -> Result<MitcTransientObservation> {
        self.admit_state(state, cx)?;
        let velocity = self.mechanics.velocities(state)?;
        let system = self.mechanics.system();
        let stored_energy_j = system.hamiltonian(state);
        let material_heat_w = system.resistance_power_from_effort(&system.effort(state))?;
        if !stored_energy_j.is_finite() || !material_heat_w.is_finite() || material_heat_w < 0.0 {
            return Err(MitcTransientError::Invalid(
                "unrepresentable physical storage or heat",
            ));
        }
        let displacement = self.model.expand(&state[..self.mechanics.coordinates()])?;
        let velocity = self.model.expand(&velocity)?;
        poll(cx)?;
        Ok(MitcTransientObservation {
            displacement,
            velocity,
            stored_energy_j,
            material_heat_w,
        })
    }
}

fn factor_rows(
    factor: &ShearDissipation,
    strain: &[ShellStrainPoint; 18],
    scale: f64,
) -> Result<[[f64; 12]; 18]> {
    let original =
        strain.map(|e| factor.channels(e.membrane, e.gradient_per_m, e.transverse_shear));
    let rows = original.map(|row| row.map(|v| scale * v));
    if !scale.is_finite()
        || scale <= 0.0
        || rows.iter().flatten().any(|v| !v.is_finite())
        || original
            .iter()
            .flatten()
            .zip(rows.iter().flatten())
            .any(|(source, scaled)| *source != 0.0 && *scaled == 0.0)
    {
        return Err(MitcTransientError::Invalid(
            "unrepresentable physical material factor",
        ));
    }
    Ok(rows)
}
