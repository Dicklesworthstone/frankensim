//! Quasistatic normal recontact of explicitly completed tensile fractures.
//!
//! Native P1 equilibrium and the original `fs-contact` Signorini natural map
//! are solved simultaneously in displacement and normal-force coordinates.
//! No inverse of an isolated, possibly fully softened fragment is required.
//! The existing Newton--FGMRES owner handles the indefinite mixed system;
//! original geometry and each selected native tangent are reused by Krylov.
//! All material trials start from the same accepted history.
//!
//! This is supported, small-strain equilibrium on fixed paired traces. Loaded
//! pairs must already touch and remain tangentially matched. Closing an open
//! gap requires a separate resolved endpoint before applying contact pressure.
//! Free rigid modes, impact, finite sliding, collision discovery and branch
//! stability are not supplied by this endpoint. Numerical scales introduce no
//! physical compliance, stored contact energy or artificial fracture loss.
//!
//! The normal law is `p >= 0, g >= 0, p*g = 0`, equivalently
//! `p - max(p - k*g, 0) = 0`; see Kunisch and Stadler (2005),
//! DOI 10.1051/m2an:2005036, and the shared sustained-contact owner. Here
//! `g = n.(u_b-u_a)` and the work-conjugate nodal forces are `(-p*n,+p*n)`.

use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;

use fs_contact::sustained::{ContactLinearization, ContactPoint, NaturalMapScales, linearize};
use fs_exec::Cx;
use fs_material::brittle::BrittleSmallStrain;
use fs_solver::{
    NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NonlinearProblem, SolverRunProgress,
};
use fs_sparse::Csr;

use super::ReleasedBrittleFace;
use crate::brittle::solid::{
    BrittleTetControls, BrittleTetEnergy, BrittleTetError, BrittleTetSolid, BrittleTetStep,
    checkpoint, energy_gate,
};
use crate::plastic3::{
    IncrementRecord, IncrementSettings, PreparedSmallStrainTet, SmallStrainLawField,
    SmallStrainTetProblem, SmallStrainTetSolution,
};

/// SI conditioning, bounded solver work and independent physical gates.
#[derive(Clone, Copy, Debug)]
pub struct BrittleFaceEquilibriumControls {
    pub newton: NewtonKrylovConfig,
    pub maximum_iterations: usize,
    pub maximum_operator_evaluations: usize,
    pub maximum_contact_rows: usize,
    pub maximum_workspace_entries: usize,
    pub displacement_scale_m: f64,
    pub force_scale_n: f64,
    pub force_tolerance_n: f64,
    pub activation_tolerance_m: f64,
    pub tangential_tolerance_m: f64,
    pub complementarity_tolerance_j: f64,
    pub contact_work_tolerance_j: f64,
    pub force_moment_tolerance_n_m: f64,
    pub work: BrittleTetControls,
}

/// One distinct physical normal covector, in first-occurrence order.
/// Force is the total for all its exact positive aliases, not an allocation
/// of pressure to individual adjacent faces. Opposite half-spaces never alias.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleFaceEquilibriumRow {
    pub face_id: u64,
    pub vertex: usize,
    pub nodes: [usize; 2],
    pub normal: [f64; 3],
    pub normal_force_n: f64,
    pub gap_m: f64,
}

/// An original face/vertex witness of an aggregate physical contact row.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleFaceEquilibriumAlias {
    pub face_id: u64,
    pub vertex: usize,
    pub nodes: [usize; 2],
    pub normal: [f64; 3],
    pub aggregate_row: usize,
}

#[derive(Debug)]
pub struct BrittleFaceEquilibriumStep<'a> {
    /// Accepted original material state and native work receipt, including the
    /// contact nodal forces. Source/history/topology ownership is unchanged.
    pub native: BrittleTetStep<'a>,
    /// Applied dead-load/support work with internal contact work removed.
    /// Its equilibrium residual independently passes the original work gate.
    pub external_energy: BrittleTetEnergy,
    pub rows: Vec<BrittleFaceEquilibriumRow>,
    /// Every original non-tied pair, including shared-edge duplicates. These
    /// witnesses do not assign individual forces within an aggregate row.
    pub original_pairs: Vec<BrittleFaceEquilibriumAlias>,
    pub nonlinear: NewtonReport,
    pub maximum_force_residual_n: f64,
    pub normal_contact_work_j: f64,
    pub absolute_contact_work_j: f64,
    pub absolute_contact_couple_n_m: f64,
    /// Native force/tangent assemblies, including independent final replay.
    pub operator_evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BrittleFaceEquilibriumError {
    Invalid(&'static str),
    Solid(BrittleTetError),
    Solver(String),
    NotConverged(NewtonReport),
    Budget(&'static str),
    Geometry {
        face: u64,
        vertex: usize,
        reason: &'static str,
    },
    ForceResidual {
        maximum_n: f64,
        tolerance_n: f64,
    },
}
impl core::fmt::Display for BrittleFaceEquilibriumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "released-face equilibrium refused: {self:?}")
    }
}
impl std::error::Error for BrittleFaceEquilibriumError {}
impl From<BrittleTetError> for BrittleFaceEquilibriumError {
    fn from(error: BrittleTetError) -> Self {
        Self::Solid(error)
    }
}

#[derive(Clone, Copy)]
struct Row {
    face: ReleasedBrittleFace,
    vertex: usize,
    nodes: [usize; 2],
}
impl Row {
    fn separation(self, u: &[f64]) -> [f64; 3] {
        core::array::from_fn(|axis| u[3 * self.nodes[1] + axis] - u[3 * self.nodes[0] + axis])
    }
    fn gap(self, u: &[f64]) -> f64 {
        dot(self.face.normal, self.separation(u))
    }
    fn scatter(self, force: f64, output: &mut [f64]) {
        for axis in 0..3 {
            let value = force * self.face.normal[axis];
            output[3 * self.nodes[0] + axis] -= value;
            output[3 * self.nodes[1] + axis] += value;
        }
    }
    fn geometry(self, reason: &'static str) -> BrittleFaceEquilibriumError {
        BrittleFaceEquilibriumError::Geometry {
            face: self.face.id,
            vertex: self.vertex,
            reason,
        }
    }
}

/// Compare full physical covectors exactly, including held coordinates.
/// Modal signed/support-reduced aliases cannot be used here: opposite signs
/// change the unilateral cone, and different held nodes change reactions.
fn group_rows(
    faces: &[ReleasedBrittleFace],
    seeds: &[f64],
    maximum: usize,
    cx: &Cx<'_>,
) -> Result<(Vec<Row>, Vec<BrittleFaceEquilibriumAlias>, Vec<f64>), BrittleFaceEquilibriumError> {
    use BrittleFaceEquilibriumError as E;
    let mut unique = BTreeMap::new();
    let mut rows: Vec<Row> = Vec::new();
    let mut original = Vec::new();
    let mut summed_seeds = Vec::<f64>::new();
    for &face in faces {
        for vertex in 0..3 {
            checkpoint(cx)?;
            let nodes = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
            if nodes[0] == nodes[1] {
                continue;
            }
            if original.len() >= maximum {
                return Err(E::Budget("released contact rows"));
            }
            let seed = *seeds
                .get(original.len())
                .ok_or(E::Invalid("one seed per original non-tied pair required"))?;
            if !seed.is_finite() || seed < 0.0 {
                return Err(E::Invalid("finite nonnegative contact seeds required"));
            }
            let (ordered, orientation) = if nodes[0] < nodes[1] {
                (nodes, 1.0)
            } else {
                ([nodes[1], nodes[0]], -1.0)
            };
            let normal = face.normal.map(|value| {
                let value = orientation * value;
                if value == 0.0 { 0 } else { value.to_bits() }
            });
            let index = match unique.get(&(ordered, normal)).copied() {
                Some(index) => {
                    let prior: &mut Row = &mut rows[index];
                    prior.face.opening_tolerance_m =
                        prior.face.opening_tolerance_m.min(face.opening_tolerance_m);
                    summed_seeds[index] += seed;
                    if !summed_seeds[index].is_finite() {
                        return Err(E::Invalid("unrepresentable aggregate contact seed"));
                    }
                    index
                }
                None => {
                    let index = rows.len();
                    unique.insert((ordered, normal), index);
                    rows.push(Row {
                        face,
                        vertex,
                        nodes,
                    });
                    summed_seeds.push(seed);
                    index
                }
            };
            original.push(BrittleFaceEquilibriumAlias {
                face_id: face.id,
                vertex,
                nodes,
                normal: face.normal,
                aggregate_row: index,
            });
        }
    }
    if seeds.len() != original.len() {
        return Err(E::Invalid("one seed per original non-tied pair required"));
    }
    Ok((rows, original, summed_seeds))
}

impl<'a> BrittleTetSolid<'a> {
    /// Continue a failed supported specimen under normal contact pressure.
    ///
    /// Prescribed coordinates and dead loads are absolute targets. One seed
    /// per ORIGINAL non-tied released pair selects the numerical active set;
    /// a positive seed activates an initially touching pair even when its
    /// independently softened fragment has no normal stiffness. Seeds are not
    /// retained forces, physical preload or history. Exact positive full-row
    /// aliases combine their seeds and normal force, preserving all original
    /// witnesses and their strictest opening allowance. No per-face pressure
    /// distribution is inferred. Other dependencies or unconstrained rigid
    /// modes may refuse the mixed solve.
    ///
    /// The actual candidate must independently satisfy free force equilibrium,
    /// the original opening tolerances, nonnegative contact pressure, physical
    /// complementarity and pairwise work/moment gates. No candidate position,
    /// force or material state is projected. Every refusal leaves `self` intact.
    pub fn advance_released_face_equilibrium(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        initial_normal_forces_n: &[f64],
        controls: BrittleFaceEquilibriumControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceEquilibriumStep<'a>, BrittleFaceEquilibriumError> {
        use BrittleFaceEquilibriumError as E;
        checkpoint(cx)?;
        validate(controls)?;
        self.require_released_faces(&self.solution.displacement_m, cx)?;
        let dofs = self.solution.displacement_m.len() * 3;
        if prescribed_m.len() > dofs || nodal_forces_n.len() > dofs {
            return Err(E::Invalid("load entries exceed original DOF count"));
        }
        let mut held = vec![false; dofs];
        for &(dof, value) in prescribed_m {
            checkpoint(cx)?;
            if dof >= dofs || !value.is_finite() || held[dof] {
                return Err(E::Invalid("finite unique prescribed coordinates required"));
            }
            held[dof] = true;
        }
        let free: Vec<_> = (0..dofs).filter(|&dof| !held[dof]).collect();
        let mut external = vec![0.0; dofs];
        for &(dof, force) in nodal_forces_n {
            checkpoint(cx)?;
            if dof >= dofs || !force.is_finite() {
                return Err(E::Invalid("finite existing load coordinates required"));
            }
            external[dof] += force;
            if !external[dof].is_finite() {
                return Err(E::Invalid("unrepresentable dead load"));
            }
        }
        let (rows, original_pairs, seeds) = group_rows(
            &self.released_faces,
            initial_normal_forces_n,
            controls.maximum_contact_rows,
            cx,
        )?;
        if free.is_empty() || rows.is_empty() {
            return Err(E::Invalid("free coordinates and released pairs required"));
        }
        let dimension = free
            .len()
            .checked_add(rows.len())
            .ok_or(E::Budget("mixed dimension"))?;
        // Original P1 geometry/CSR is bounded by TetAssemblyBudget. This
        // additional allowance covers mixed vectors, contact cache and Krylov.
        let restart = controls.newton.linear_restart;
        let workspace = restart
            .checked_add(1)
            .and_then(|r| {
                dimension
                    .checked_mul(2)
                    .and_then(|n| n.checked_add(restart))
                    .and_then(|w| r.checked_mul(w))
            })
            .and_then(|k| dimension.checked_mul(24).and_then(|v| k.checked_add(v)))
            .and_then(|k| dofs.checked_mul(8).and_then(|v| k.checked_add(v)))
            .and_then(|k| {
                controls
                    .maximum_iterations
                    .checked_mul(32)
                    .and_then(|v| k.checked_add(v))
            })
            // FGMRES retains residual history across ALL configured restart
            // cycles; construction and Newton receipts can clone that history.
            .and_then(|k| {
                restart
                    .checked_add(1)
                    .and_then(|r| r.checked_mul(controls.newton.max_linear_cycles))
                    .and_then(|h| h.checked_mul(8))
                    .and_then(|h| k.checked_add(h))
            })
            .and_then(|k| {
                original_pairs
                    .len()
                    .checked_mul(24)
                    .and_then(|v| k.checked_add(v))
            })
            .ok_or(E::Budget("mixed workspace overflow"))?;
        if workspace > controls.maximum_workspace_entries {
            return Err(E::Budget("mixed workspace"));
        }
        let laws: Vec<_> = self
            .materials
            .iter()
            .map(|source| source.small_strain_law())
            .collect();
        let native = SmallStrainTetProblem {
            nodes_m: self.geometry.nodes_m(),
            tetrahedra: self.geometry.tetrahedra(),
            laws: SmallStrainLawField::PerElement(&laws),
            prescribed_m,
            nodal_forces_n,
            load_path: &[1.0],
            budget: self.budget,
            settings: IncrementSettings::default(),
        };
        let prepared = native
            .prepare_evaluator(cx)
            .map_err(BrittleTetError::from)?;
        let mut initial: Vec<_> = free
            .iter()
            .map(|&dof| {
                self.solution.displacement_m.as_flattened()[dof] / controls.displacement_scale_m
            })
            .collect();
        initial.extend(seeds.iter().map(|p| p / controls.force_scale_n));
        let problem = Equilibrium {
            solid: self,
            prepared: &prepared,
            prescribed_m,
            external: &external,
            free,
            rows,
            controls,
            cx,
            evaluations: Cell::new(0),
            fault: RefCell::new(None),
            cached: RefCell::new(None),
            scratch: RefCell::new(Scratch {
                force: vec![0.0; dofs],
                displacement_direction: vec![0.0; dofs],
            }),
        };
        let solver = NewtonKrylovState::new(&problem, initial, controls.newton);
        problem.poll()?;
        let mut solver = solver.map_err(|e| E::Solver(e.to_string()))?;
        let run = solver.run_controlled(&problem, controls.maximum_iterations, &mut || {
            problem.poll().is_err()
        });
        problem.poll()?;
        if run.progress != SolverRunProgress::Complete {
            return Err(E::Solver(format!("{:?}", run.progress)));
        }
        if !run.report.converged {
            return Err(E::NotConverged(run.report));
        }
        // Force a fresh original constitutive evaluation after convergence.
        *problem.cached.borrow_mut() = None;
        let mut residual = vec![0.0; dimension];
        problem.action(&solver.x, None, &mut residual)?;
        let cache = problem.cached.borrow();
        let final_state = cache
            .as_ref()
            .expect("fresh residual retains original response");
        let old_u = self.solution.displacement_m.as_flattened();
        let new_u = &final_state.displacement;
        let mut contact_load = vec![0.0; dofs];
        let mut records = Vec::with_capacity(problem.rows.len());
        let mut normal_work = 0.0;
        let mut absolute_work = 0.0;
        let mut absolute_couple = 0.0;
        for (i, &row) in problem.rows.iter().enumerate() {
            checkpoint(cx)?;
            let force = solver.x[problem.free.len() + i] * controls.force_scale_n;
            let old_gap = row.gap(old_u);
            let gap = row.gap(new_u);
            if !force.is_finite() || force < 0.0 {
                return Err(row.geometry("nonnegative finite normal force required"));
            }
            if !gap.is_finite() || gap < -row.face.opening_tolerance_m {
                return Err(
                    row.geometry("released pair overlaps beyond original opening allowance")
                );
            }
            energy_gate(
                "normal complementarity",
                force * gap,
                controls.complementarity_tolerance_j,
            )?;
            if force > controls.force_tolerance_n {
                if old_gap.abs() > controls.activation_tolerance_m
                    || gap.abs() > controls.activation_tolerance_m
                {
                    return Err(row.geometry("loaded contact must touch at both endpoints"));
                }
                for u in [old_u, new_u.as_slice()] {
                    let separation = row.separation(u);
                    let tangent = norm(core::array::from_fn(|axis| {
                        separation[axis] - row.gap(u) * row.face.normal[axis]
                    }));
                    if !tangent.is_finite() || tangent > controls.tangential_tolerance_m {
                        return Err(
                            row.geometry("loaded reference pairs must remain tangentially matched")
                        );
                    }
                }
            }
            let work = force * (gap - old_gap);
            normal_work += work;
            absolute_work += work.abs();
            let initial_couple = norm(cross(row.separation(old_u), row.face.normal)) * force;
            let final_couple = norm(cross(row.separation(new_u), row.face.normal)) * force;
            absolute_couple += initial_couple.max(final_couple);
            row.scatter(force, &mut contact_load);
            records.push(BrittleFaceEquilibriumRow {
                face_id: row.face.id,
                vertex: row.vertex,
                nodes: row.nodes,
                normal: row.face.normal,
                normal_force_n: force,
                gap_m: gap,
            });
        }
        energy_gate(
            "absolute normal contact work",
            absolute_work,
            controls.contact_work_tolerance_j,
        )?;
        energy_gate(
            "absolute normal force couple",
            absolute_couple,
            controls.force_moment_tolerance_n_m,
        )?;
        let total: Vec<_> = external
            .iter()
            .zip(&contact_load)
            .enumerate()
            .map(|(dof, (&f, &c))| (dof, f + c))
            .collect();
        if total.iter().any(|(_, f)| !f.is_finite()) {
            return Err(E::Invalid("unrepresentable final total load"));
        }
        // Recheck in actual newtons, without recovering a physical residual
        // from a possibly underflowed normalized solver equation.
        let mut maximum_force_residual_n = 0.0_f64;
        for &dof in &problem.free {
            checkpoint(cx)?;
            let residual_n = (final_state.internal[dof] - external[dof]) - contact_load[dof];
            if !residual_n.is_finite() {
                return Err(E::Invalid("nonfinite physical free-force residual"));
            }
            maximum_force_residual_n = maximum_force_residual_n.max(residual_n.abs());
        }
        for block in &final_state.contact {
            maximum_force_residual_n = maximum_force_residual_n.max(block.residual_n[0].abs());
        }
        if maximum_force_residual_n > controls.force_tolerance_n {
            return Err(E::ForceResidual {
                maximum_n: maximum_force_residual_n,
                tolerance_n: controls.force_tolerance_n,
            });
        }
        let mut actual_contact_work = 0.0;
        for dof in 0..dofs {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            actual_contact_work =
                contact_load[dof].mul_add(new_u[dof] - old_u[dof], actual_contact_work);
        }
        energy_gate(
            "contact force/jump work transpose",
            actual_contact_work - normal_work,
            controls.contact_work_tolerance_j,
        )?;
        let strains = prepared
            .strain_field(new_u, cx)
            .map_err(BrittleTetError::from)?;
        let mut states = Vec::with_capacity(strains.len());
        let mut stresses = Vec::with_capacity(strains.len());
        let mut evolving_elements = 0;
        for (element, (source, strain)) in self.materials.iter().zip(&strains).enumerate() {
            checkpoint(cx)?;
            let old = &self.solution.element_states[element];
            let trial = source
                .law()
                .trial(strain, old)
                .map_err(|source| BrittleTetError::Material { element, source })?;
            evolving_elements += usize::from(trial.state != *old);
            states.push(trial.state);
            stresses.push(trial.stress_pa);
        }
        let solution = SmallStrainTetSolution {
            displacement_m: new_u.as_chunks::<3>().0.to_vec(),
            element_strain: strains,
            element_stress: stresses,
            element_states: states,
            increments: vec![IncrementRecord {
                load_factor: 1.0,
                residual_history: vec![maximum_force_residual_n],
                reactions_n: prescribed_m
                    .iter()
                    .map(|&(dof, _)| {
                        (final_state.internal[dof] - external[dof]) - contact_load[dof]
                    })
                    .collect(),
                evolving_elements,
                indefinite_fallbacks: 0,
            }],
        };
        let accepted =
            self.accept_equilibrium(solution, prescribed_m, &total, controls.work, cx)?;
        let mut external_energy = accepted.energy;
        external_energy.dead_load_endpoint_work_j -= actual_contact_work;
        external_energy.equilibrium_work_residual_j = external_energy.dead_load_endpoint_work_j
            + external_energy.support_endpoint_work_j
            - external_energy.internal_endpoint_work_j;
        energy_gate(
            "external support/dead-load work without internal contact",
            external_energy.equilibrium_work_residual_j,
            controls.work.equilibrium_work_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(BrittleFaceEquilibriumStep {
            native: accepted,
            external_energy,
            rows: records,
            original_pairs,
            nonlinear: run.report,
            maximum_force_residual_n,
            normal_contact_work_j: normal_work,
            absolute_contact_work_j: absolute_work,
            absolute_contact_couple_n_m: absolute_couple,
            operator_evaluations: problem.evaluations.get(),
        })
    }
}

struct Cached {
    point: Vec<f64>,
    displacement: Vec<f64>,
    internal: Vec<f64>,
    contact: Vec<ContactLinearization>,
    tangent: Option<Csr>,
}

struct Scratch {
    force: Vec<f64>,
    displacement_direction: Vec<f64>,
}

struct Equilibrium<'s, 'p, 'a, 'c, 'cx> {
    solid: &'s BrittleTetSolid<'a>,
    prepared: &'s PreparedSmallStrainTet<'p, 's, BrittleSmallStrain<'a>>,
    prescribed_m: &'s [(usize, f64)],
    external: &'s [f64],
    free: Vec<usize>,
    rows: Vec<Row>,
    controls: BrittleFaceEquilibriumControls,
    cx: &'c Cx<'cx>,
    evaluations: Cell<usize>,
    fault: RefCell<Option<BrittleFaceEquilibriumError>>,
    cached: RefCell<Option<Cached>>,
    scratch: RefCell<Scratch>,
}

impl Equilibrium<'_, '_, '_, '_, '_> {
    fn fail(&self, error: BrittleFaceEquilibriumError) {
        if self.fault.borrow().is_none() {
            *self.fault.borrow_mut() = Some(error);
        }
    }
    fn poll(&self) -> Result<(), BrittleFaceEquilibriumError> {
        if let Some(error) = self.fault.borrow().as_ref() {
            return Err(error.clone());
        }
        checkpoint(self.cx)?;
        Ok(())
    }
    fn charge(&self) -> Result<(), BrittleFaceEquilibriumError> {
        self.poll()?;
        if self.evaluations.get() >= self.controls.maximum_operator_evaluations {
            return Err(BrittleFaceEquilibriumError::Budget(
                "native force/tangent assemblies",
            ));
        }
        self.evaluations.set(self.evaluations.get() + 1);
        Ok(())
    }
    fn ensure(&self, x: &[f64], tangent: bool) -> Result<(), BrittleFaceEquilibriumError> {
        use BrittleFaceEquilibriumError as E;
        self.poll()?;
        if x.len() != self.dimension() || x.iter().any(|v| !v.is_finite()) {
            return Err(E::Invalid("finite complete mixed coordinates required"));
        }
        let matches = self
            .cached
            .borrow()
            .as_ref()
            .is_some_and(|cache| cache.point == x);
        if !matches {
            self.charge()?;
            let mut displacement = self.solid.solution.displacement_m.as_flattened().to_vec();
            for &(dof, value) in self.prescribed_m {
                displacement[dof] = value;
            }
            for (&dof, &value) in self.free.iter().zip(x) {
                displacement[dof] = value * self.controls.displacement_scale_m;
            }
            let internal = self
                .prepared
                .internal_force(&displacement, &self.solid.solution.element_states, self.cx)
                .map_err(BrittleTetError::from)?;
            let k = self.controls.force_scale_n / self.controls.displacement_scale_m;
            let mut contact = Vec::with_capacity(self.rows.len());
            for (i, &row) in self.rows.iter().enumerate() {
                checkpoint(self.cx)?;
                let point = ContactPoint {
                    force_n: [
                        x[self.free.len() + i] * self.controls.force_scale_n,
                        0.0,
                        0.0,
                    ],
                    gap_m: row.gap(&displacement),
                    slip_increment_m: [0.0; 2],
                    friction_coefficient: 0.0,
                };
                contact.push(
                    linearize(
                        point,
                        NaturalMapScales {
                            normal_n_per_m: k,
                            tangent_n_per_m: k,
                        },
                    )
                    .map_err(|_| E::Invalid("unrepresentable native normal contact law"))?,
                );
            }
            *self.cached.borrow_mut() = Some(Cached {
                point: x.to_vec(),
                displacement,
                internal,
                contact,
                tangent: None,
            });
        }
        if tangent
            && self
                .cached
                .borrow()
                .as_ref()
                .is_some_and(|cache| cache.tangent.is_none())
        {
            self.charge()?;
            let mut cache = self.cached.borrow_mut();
            let cache = cache.as_mut().expect("prepared original response");
            cache.tangent = Some(
                self.prepared
                    .tangent_matrix(
                        &cache.displacement,
                        &self.solid.solution.element_states,
                        self.cx,
                    )
                    .map_err(BrittleTetError::from)?,
            );
        }
        Ok(())
    }
    fn action(
        &self,
        x: &[f64],
        direction: Option<&[f64]>,
        output: &mut [f64],
    ) -> Result<(), BrittleFaceEquilibriumError> {
        if output.len() != self.dimension()
            || direction
                .is_some_and(|v| v.len() != self.dimension() || v.iter().any(|x| !x.is_finite()))
        {
            return Err(BrittleFaceEquilibriumError::Invalid(
                "mixed operator dimensions",
            ));
        }
        self.ensure(x, direction.is_some())?;
        let cached = self.cached.borrow();
        let cache = cached.as_ref().expect("successful original response");
        let scale = self.controls.force_scale_n;
        let mut scratch = self.scratch.borrow_mut();
        let Scratch {
            force,
            displacement_direction: du,
        } = &mut *scratch;
        force.fill(0.0);
        du.fill(0.0);
        if let Some(v) = direction {
            for (&dof, &value) in self.free.iter().zip(v) {
                du[dof] = value * self.controls.displacement_scale_m;
            }
            let tangent = cache.tangent.as_ref().expect("prepared native tangent");
            for (index, &dof) in self.free.iter().enumerate() {
                if index % 256 == 0 {
                    checkpoint(self.cx)?;
                }
                let (columns, values) = tangent.row(dof);
                for (&column, &coefficient) in columns.iter().zip(values) {
                    force[dof] = coefficient.mul_add(du[column], force[dof]);
                }
            }
        } else {
            for (dof, value) in force.iter_mut().enumerate() {
                *value = cache.internal[dof] - self.external[dof];
            }
        }
        for (i, &row) in self.rows.iter().enumerate() {
            checkpoint(self.cx)?;
            let block = &cache.contact[i];
            let contact_force = direction.unwrap_or(x)[self.free.len() + i] * scale;
            // Equilibrium is f_int - f_external - B^T p = 0.
            row.scatter(-contact_force, force);
            output[self.free.len() + i] = match direction {
                None => block.residual_n[0] / scale,
                Some(_) => {
                    (block.force_derivative[0][0] * contact_force
                        + block.path_derivative[0][0] * row.gap(du))
                        / scale
                }
            };
        }
        for (value, &dof) in output.iter_mut().zip(&self.free) {
            *value = force[dof] / scale;
        }
        if output.iter().any(|v| !v.is_finite()) {
            return Err(BrittleFaceEquilibriumError::Invalid(
                "unrepresentable mixed residual/tangent",
            ));
        }
        self.poll()
    }
}

impl NonlinearProblem for Equilibrium<'_, '_, '_, '_, '_> {
    fn dimension(&self) -> usize {
        self.free.len() + self.rows.len()
    }
    fn residual(&self, x: &[f64], output: &mut [f64]) {
        if let Err(error) = self.action(x, None, output) {
            self.fail(error);
            output.fill(f64::NAN);
        }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], output: &mut [f64]) {
        if let Err(error) = self.action(x, Some(direction), output) {
            self.fail(error);
            output.fill(f64::NAN);
        }
    }
}

fn validate(c: BrittleFaceEquilibriumControls) -> Result<(), BrittleFaceEquilibriumError> {
    use BrittleFaceEquilibriumError as E;
    c.newton.validate().map_err(|e| E::Solver(e.to_string()))?;
    if c.maximum_iterations == 0
        || c.maximum_operator_evaluations == 0
        || c.maximum_contact_rows == 0
        || c.maximum_workspace_entries == 0
    {
        return Err(E::Budget(
            "positive work, row and workspace allowances required",
        ));
    }
    if [c.displacement_scale_m, c.force_scale_n]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || !(c.force_scale_n / c.displacement_scale_m).is_finite()
        || c.force_scale_n / c.displacement_scale_m == 0.0
    {
        return Err(E::Invalid("representable positive SI scales required"));
    }
    if [
        c.force_tolerance_n,
        c.activation_tolerance_m,
        c.tangential_tolerance_m,
        c.complementarity_tolerance_j,
        c.contact_work_tolerance_j,
        c.force_moment_tolerance_n_m,
        c.work.material_energy_tolerance_j,
        c.work.equilibrium_work_tolerance_j,
        c.work.maximum_endpoint_work_defect_j,
    ]
    .iter()
    .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(E::Invalid(
            "finite nonnegative physical tolerances required",
        ));
    }
    Ok(())
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}
fn norm(a: [f64; 3]) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
