//! The nonlinear enthalpy body on the composition contract.
//!
//! [`EnthalpyBody`] puts one [`EnthalpyBackwardEuler`] owner, with its
//! reference nodal masses, its caloric chart `T(h)` and its problem's
//! conductivity, behind the heat ports of a
//! [`ConductionBody`](super::ConductionBody): each named boundary patch is a
//! mean port or a trace port ([`PatchPort`]), and its heat reaches the owner
//! through the same Neumann rows. Prepared by
//! `MaterialEnthalpyTransport::prepare_uniform`, the chart integrates the
//! card's Cp(T), the masses hold the card's reference density, and the
//! conductivity is the card's k(T): every coefficient of the body comes from
//! one material card.
//!
//! The state is the nodal specific enthalpy [J/kg]; the ports return the END
//! temperature `T(h)`, a mean port its patch's area mean. The stored energy
//! is the owner's inventory `Σ m_v h_v`. The step is the owner's unchanged
//! `advance`. The tangent is the owner's `linearize_step` of that same
//! endpoint, whose public adjoint is a transpose solve: the body forms the
//! port Jacobian one row for each port entry, by one transpose solve each,
//! and applies that one matrix and its transpose, so the two products are
//! exact transposes. The body also offers the step as the owner's Newton
//! residual in its end enthalpy (`StageResidual`), for a monolithic solve:
//! its derivative is the tangent the owner's Newton applies, and `accept`
//! books the receipt `propose` books.
//!
//! The entropy generated inside the body is the conduction term alone, as for
//! a conduction body, over the owner's endpoint operator with its phase
//! conductivity. The owner's refusals stand: a step whose endpoint leaves the
//! chart, or whose Newton solve misses its policy, refuses the step. Robin,
//! Dirichlet and contact rows, volumetric sources and temperature-dependent
//! sheet conduction are refused, as for a conduction body. Latent storage
//! in the chart is the owner's, and is not tested through ports. No contact
//! between the traces of two enthalpy bodies is built here.
//!
//! The state has a snapshot codec, as a conduction body's has. Its identity
//! binds the owner's caloric chart by the card it came from and the chart's
//! own identity, which holds its receipts. A solver with more than its
//! mesh, masses and one chart (nodal charts, surface storage, a phase
//! multiplier on the conductivity from `with_phase_conductivity`, sourced
//! transport) offers no identity of those parts, so a body on one refuses
//! to be written or read, as does a body whose boundary carries sheet
//! conduction.

use fs_blake3::DomainHasher;
use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::component::{
    Class, Component, PortSpec, PortTangent, Proposal, Receipt, StageResidual,
};

use super::{
    Patch, PatchPort, admit_boundary, admit_ports, bind_policy, bind_problem, bind_values,
    boundary_for, check_width, decode_field, encode_field, heat_in_j, identity_of, input_bars,
    port_readings, reading_weights, unsheeted,
};
use crate::transient::enthalpy::adjoint::EnthalpyAdjointError;
use crate::transient::enthalpy::{EnthalpyBackwardEuler, EnthalpyStepConfig, EnthalpyStepSolution};
use crate::{ConductionError, ConductionProblem, LinearConfig};

/// The enthalpy owner behind heat ports on named boundary patches.
///
/// Port `i` has the id of region `ports[i]`, as for a
/// [`ConductionBody`](super::ConductionBody): a mean port's input is the heat
/// rate into the body [W] and its output the patch's area-mean end
/// temperature [K]; a trace port's inputs are the heat loads into the body at
/// the patch's vertices [W] and its outputs their end temperatures [K]. The
/// state is the nodal specific enthalpy [J/kg] in vertex order.
#[derive(Debug)]
pub struct EnthalpyBody<'m, 'c> {
    solver: EnthalpyBackwardEuler<'m, 'c>,
    problem: ConductionProblem<'m>,
    config: EnthalpyStepConfig,
    /// The transpose solve of every row of the port Jacobian.
    adjoint: LinearConfig,
    ports: Vec<PortSpec>,
    patches: Vec<Patch>,
    /// Entries of all ports together, in the inputs and in the outputs.
    width: usize,
    /// The owner's reference nodal masses [kg]: the stored energy of a field
    /// `h` is `Σ_v m_v h_v`.
    masses: Vec<f64>,
}

impl<'m, 'c> EnthalpyBody<'m, 'c> {
    /// Put `solver` and `problem` behind heat ports on the regions `ports`.
    ///
    /// The boundary and source must be those a conduction body admits: every
    /// region a uniform zero-flux Neumann row, no temperature-dependent sheet
    /// conduction, a zero uniform source. `config` is the owner's step policy,
    /// used for every step and tangent; `adjoint` the policy of each transpose
    /// solve of the tangent.
    ///
    /// # Errors
    /// [`ConductionError::Config`] naming the refused port, region or source;
    /// a solver prepared on another mesh.
    pub fn new(
        solver: EnthalpyBackwardEuler<'m, 'c>,
        problem: ConductionProblem<'m>,
        ports: &[PatchPort<'_>],
        config: EnthalpyStepConfig,
        adjoint: LinearConfig,
    ) -> Result<Self, ConductionError> {
        if !std::ptr::eq(solver.mesh(), problem.mesh) {
            return Err(super::invalid(
                "the solver and the problem must use the same mesh".into(),
            ));
        }
        admit_boundary(&problem)?;
        let (ports, patches, width) = admit_ports(&problem, ports)?;
        let masses = solver.reference_nodal_masses_kg().to_vec();
        Ok(Self {
            solver,
            problem,
            config,
            adjoint,
            ports,
            patches,
            width,
            masses,
        })
    }

    /// Area [m²] of port `port`'s patch.
    #[must_use]
    pub fn patch_area_m2(&self, port: usize) -> f64 {
        self.patches[port].area_m2
    }

    /// The mesh vertices of trace port `port`, in the order of its inputs and
    /// outputs; `None` for a mean port.
    #[must_use]
    pub fn trace_vertices(&self, port: usize) -> Option<&[usize]> {
        self.patches[port]
            .trace
            .as_ref()
            .map(|trace| trace.vertices.as_slice())
    }

    /// The row-major port Jacobian `d outputs / d inputs` [K/W] of the
    /// owner's endpoint from `old` under `inputs`: row `i` by one transpose
    /// solve seeded with output `i`'s reading weights.
    fn port_jacobian(
        &self,
        old: &[f64],
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        self.check_width("inputs", inputs.len())?;
        let boundary = boundary_for(&self.patches, &self.problem, inputs)?;
        let problem = ConductionProblem {
            boundary: &boundary,
            ..self.problem
        };
        let step = self
            .solver
            .linearize_step(cx, problem, None, old, dt_s, self.config)?;
        let mut rows = vec![0.0; self.width * self.width];
        let mut unit = vec![0.0; self.width];
        let mut weights = vec![0.0; self.masses.len()];
        for (i, row) in rows.chunks_exact_mut(self.width).enumerate() {
            unit.fill(0.0);
            unit[i] = 1.0;
            weights.fill(0.0);
            reading_weights(&self.patches, &unit, &mut weights);
            let seed = step.temperature_pullback(cx, &weights)?;
            let gradient = step.pullback(cx, &seed, self.adjoint)?;
            let nodal = step.nodal_load_pullback(cx, &gradient.adjoint)?;
            input_bars(&self.patches, &nodal, row);
        }
        Ok(rows)
    }

    fn check_width(&self, vector: &str, length: usize) -> Result<(), ConductionError> {
        check_width(self.width, self.ports.len(), vector, length)
    }

    /// Whether the identity binds all of this body: the solver is fixed by
    /// its mesh, masses and one chart, and the boundary carries no sheet.
    fn has_codec(&self) -> bool {
        self.solver.is_plain() && unsheeted(&self.problem)
    }
}

impl Component for EnthalpyBody<'_, '_> {
    type State = Vec<f64>;
    type Error = EnthalpyAdjointError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-conduction::component::EnthalpyBody",
        state_family: "nodal-specific-enthalpy-field",
        schema_grammar: "the specific enthalpy [J/kg] at every vertex of the body's mesh, in vertex order",
        codec_grammar: "the vertex count, then each specific enthalpy as the bits of an f64",
        codec_version: 1,
    });

    /// Refuses, writing nothing, for a solver with more than its mesh,
    /// masses and one chart, and for a boundary with sheet conduction: the
    /// identity cannot bind them.
    fn encode_state(&self, state: &Vec<f64>, words: &mut Vec<u64>) -> bool {
        if !self.has_codec() {
            return false;
        }
        encode_field(state, words);
        true
    }

    /// Refuses a count that is not the mesh's, a specific enthalpy that is
    /// not finite, and everything for a body `encode_state` refuses.
    /// Whether a value lies on the chart is the owner's to refuse at the
    /// next step.
    fn decode_state(&self, words: &[u64]) -> Option<Vec<f64>> {
        if !self.has_codec() {
            return None;
        }
        decode_field(words, self.masses.len(), |_| true)
    }

    /// One hash, as four words, of what fixes the body's steps: the mesh,
    /// every element's conductivity and the ports, the reference nodal
    /// masses, the owner's caloric chart (the card it came from and the
    /// chart's own identity, which holds its receipts), the step policy and
    /// the policy of the tangent's transpose solves.
    fn identity_words(&self, words: &mut Vec<u64>) {
        let mut hasher = DomainHasher::new("fs-conduction/enthalpy-body/v1");
        bind_problem(&mut hasher, &self.problem, &self.patches);
        bind_values(&mut hasher, &self.masses);
        let chart = self.solver.phase_curve();
        hasher.update(chart.material_card_identity().as_bytes());
        hasher.update(chart.identity().as_bytes());
        bind_policy(&mut hasher, &self.config);
        bind_policy(&mut hasher, &self.adjoint);
        identity_of(&hasher, words);
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    /// Nonlinear in the port inputs through `T(h)` and k(T), and smooth
    /// between the chart's knots.
    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Vec<f64>) -> Result<f64, EnthalpyAdjointError> {
        if state.len() != self.masses.len() {
            return Err(ConductionError::FieldLength {
                field: "body specific enthalpy",
                expected: self.masses.len(),
                found: state.len(),
            }
            .into());
        }
        let energy: f64 = self.masses.iter().zip(state).map(|(m, h)| m * h).sum();
        if energy.is_finite() {
            Ok(energy)
        } else {
            Err(ConductionError::NonFinite {
                field: "body stored energy",
                bits: energy.to_bits(),
            }
            .into())
        }
    }

    fn propose(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Vec<f64>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, EnthalpyAdjointError> {
        self.check_width("inputs", inputs.len())?;
        self.check_width("outputs", outputs.len())?;
        let boundary = boundary_for(&self.patches, &self.problem, inputs)?;
        let problem = ConductionProblem {
            boundary: &boundary,
            ..self.problem
        };
        let step = self
            .solver
            .advance(cx, problem, None, old, dt_s, self.config)?;
        port_readings(&self.patches, &step.temperature, outputs);
        let receipt = self.step_receipt(cx, problem, old, dt_s, inputs, &step)?;
        new.clear();
        new.extend_from_slice(&step.specific_enthalpy_j_kg);
        Ok(Proposal::Step(receipt))
    }
}

impl EnthalpyBody<'_, '_> {
    /// The receipt of the owner's step from `old` to `step`'s endpoint under
    /// `inputs`: `propose`'s, and `accept`'s on the monolithic path.
    fn step_receipt(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        step: &EnthalpyStepSolution,
    ) -> Result<Receipt, EnthalpyAdjointError> {
        let end = &step.specific_enthalpy_j_kg;
        let generated_j_per_k = self
            .solver
            .endpoint_generation_j_k(cx, problem, old, end, dt_s)?;
        let heat_j = heat_in_j(&self.patches, dt_s, inputs);
        let stored_old_j = self.stored_energy_j(old)?;
        let stored_new_j = self.stored_energy_j(&step.specific_enthalpy_j_kg)?;
        // Each stored energy is a sum of `n` positive terms, and the owner's
        // change one of `n` terms no larger: `(n + 16) eps` of both ends
        // covers the rounding of all three, as for a conduction body.
        let units = (self.masses.len() + 16) as f64 * f64::EPSILON;
        // The owner accepts its endpoint when the Newton residual is at most
        // the larger of the absolute tolerance and the relative tolerance
        // times the step's first residual: that target is the receipt's.
        let first = step
            .newton
            .history
            .first()
            .map_or(step.newton.residual_norm, |attempt| attempt.residual_before);
        let newton = self.config.newton;
        let receipt = Receipt {
            stored_change_j: step.stored_energy_change_j,
            dissipated_j: 0.0,
            numerical_loss_j: heat_j - step.stored_energy_change_j,
            entropy_generated_j_per_k: generated_j_per_k,
            energy_tolerance_j: units * (stored_old_j.abs() + stored_new_j.abs() + heat_j.abs()),
            residual: step.newton.residual_norm,
            tolerance: newton
                .absolute_tolerance
                .max(newton.relative_tolerance * first),
        };
        Ok(receipt)
    }
}

/// The owner's Newton residual `m (h - h_old) + dt (A(T) T - b)` [J] in the
/// end enthalpy, with `T = T(h)`, and its tangent. `accept` runs the owner's
/// accounts and joule gate on the kernel's endpoint and books `propose`'s
/// receipt, with a zero residual: the owner solves nothing there.
impl StageResidual for EnthalpyBody<'_, '_> {
    fn unknowns(&self) -> usize {
        self.masses.len()
    }

    fn initial_guess(
        &self,
        old: &Vec<f64>,
        _dt_s: f64,
        _inputs: &[f64],
        z: &mut [f64],
    ) -> Result<(), EnthalpyAdjointError> {
        if old.len() != z.len() || z.len() != self.masses.len() {
            return Err(super::invalid("a stage guess needs one entry per vertex".into()).into());
        }
        z.copy_from_slice(old);
        Ok(())
    }

    fn residual(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        r: &mut [f64],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), EnthalpyAdjointError> {
        self.check_width("inputs", inputs.len())?;
        self.check_width("outputs", outputs.len())?;
        let boundary = boundary_for(&self.patches, &self.problem, inputs)?;
        let problem = ConductionProblem {
            boundary: &boundary,
            ..self.problem
        };
        let endpoint =
            self.solver
                .evaluate_endpoint(cx, problem, None, old, z, dt_s, self.config, None)?;
        if r.len() != endpoint.residual_j.len() {
            return Err(super::invalid("a residual needs one entry per vertex".into()).into());
        }
        r.copy_from_slice(&endpoint.residual_j);
        port_readings(&self.patches, &endpoint.temperature, outputs);
        Ok(())
    }

    fn jvp_z(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        _inputs: &[f64],
        z: &[f64],
        dz: &[f64],
        dr: &mut [f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), EnthalpyAdjointError> {
        self.check_width("tangent outputs", d_outputs.len())?;
        // The port loads enter `b` alone: the derivative holds no input.
        let solver = &self.solver;
        let slope = solver.stage_jacobian_apply(cx, self.problem, old, z, dt_s, false, dz, dr)?;
        let d_temperature: Vec<f64> = slope.iter().zip(dz).map(|(s, d)| s * d).collect();
        port_readings(&self.patches, &d_temperature, d_outputs);
        Ok(())
    }

    fn vjp_z(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        _inputs: &[f64],
        z: &[f64],
        bar_r: &[f64],
        bar_outputs: &[f64],
        bar_z: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), EnthalpyAdjointError> {
        self.check_width("adjoint outputs", bar_outputs.len())?;
        let solver = &self.solver;
        let slope =
            solver.stage_jacobian_apply(cx, self.problem, old, z, dt_s, true, bar_r, bar_z)?;
        let mut weights = vec![0.0; slope.len()];
        reading_weights(&self.patches, bar_outputs, &mut weights);
        for ((bar, s), w) in bar_z.iter_mut().zip(&slope).zip(&weights) {
            *bar += s * w;
        }
        Ok(())
    }

    fn accept(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        z: &[f64],
        new: &mut Vec<f64>,
        cx: &Cx<'_>,
    ) -> Result<Receipt, EnthalpyAdjointError> {
        self.check_width("inputs", inputs.len())?;
        let boundary = boundary_for(&self.patches, &self.problem, inputs)?;
        let problem = ConductionProblem {
            boundary: &boundary,
            ..self.problem
        };
        let step = self
            .solver
            .stage_endpoint(cx, problem, old, z, dt_s, self.config)?;
        let receipt = self.step_receipt(cx, problem, old, dt_s, inputs, &step)?;
        new.clear();
        new.extend_from_slice(&step.specific_enthalpy_j_kg);
        Ok(receipt)
    }
}

impl PortTangent for EnthalpyBody<'_, '_> {
    fn jvp(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), EnthalpyAdjointError> {
        self.check_width("tangent inputs", d_inputs.len())?;
        self.check_width("tangent outputs", d_outputs.len())?;
        let rows = self.port_jacobian(old, dt_s, inputs, cx)?;
        for (output, row) in d_outputs.iter_mut().zip(rows.chunks_exact(self.width)) {
            *output = row.iter().zip(d_inputs).map(|(j, d)| j * d).sum();
        }
        Ok(())
    }

    fn vjp(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), EnthalpyAdjointError> {
        self.check_width("adjoint outputs", bar_outputs.len())?;
        self.check_width("adjoint inputs", bar_inputs.len())?;
        let rows = self.port_jacobian(old, dt_s, inputs, cx)?;
        for (j, bar) in bar_inputs.iter_mut().enumerate() {
            *bar = rows
                .chunks_exact(self.width)
                .zip(bar_outputs)
                .map(|(row, b)| row[j] * b)
                .sum();
        }
        Ok(())
    }
}
