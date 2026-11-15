//! The transient finite-element body on the composition contract.
//!
//! [`ConductionBody`] puts one [`BackwardEuler`] owner, with its mesh, its
//! conductivity and its boundary partition, behind the Component contract.
//! Each named boundary patch of the body is one scalar heat port: the port
//! receives the heat rate `Q` [W] into the body and applies it on its faces as
//! the uniform outward flux `-Q / A`, through the owner's own Neumann rows. It
//! returns the patch's area mean of the END temperature of the step, the
//! temperature at which a uniform flux does its work (`∫ T q dA = Q T̄`) and
//! backward Euler's own reading of it. The step is the owner's unchanged
//! `advance`; the tangent is the owner's `linearize_step` with no Robin port,
//! driven by the same nodal loads.
//!
//! No claim: the entropy generated inside the body (conduction between nodes
//! at different temperatures) is reported as zero, a lower bound. A nonuniform
//! flux over a patch, Robin, Dirichlet and contact traces, volumetric sources,
//! temperature-dependent conductivity and sheet conduction are refused, not
//! wrapped. The state has no snapshot codec.

use fs_exec::Cx;
use fs_ports::StableId;
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};

use crate::transient::backward_euler::{BackwardEuler, StepConfig, StepLinearization};
use crate::{ConductionError, ConductionProblem, ScalarField, ThermalBc, ThermalBoundary};

/// One port's boundary patch: its region, its area and its faces.
#[derive(Debug)]
struct Patch {
    region: usize,
    area_m2: f64,
    /// Vertices and area of every face of the patch, in boundary-slot order.
    faces: Vec<([usize; 3], f64)>,
}

/// A backward-Euler finite-element body whose ports are named boundary patches.
///
/// Port `i` is a [`PortSpec::heat`] with [`PortRole::FlowIn`] and the id of
/// region `ports[i]`: its input is the heat rate into the body [W], its output
/// the patch's area-mean end temperature [K]. The state is the nodal
/// temperature field [K] in vertex order. The stored energy is `1ᵀ C T`, with
/// the owner's prepared capacitance `C`, so its change between two states is
/// the owner's own `stored_energy_change_j` up to rounding.
#[derive(Debug)]
pub struct ConductionBody<'m> {
    solver: BackwardEuler<'m>,
    problem: ConductionProblem<'m>,
    config: StepConfig,
    ports: Vec<PortSpec>,
    patches: Vec<Patch>,
    /// `C 1`: the stored energy of a field `T` is `Σ_v weight_v T_v`.
    weights: Vec<f64>,
}

impl<'m> ConductionBody<'m> {
    /// Put `solver` and `problem` behind heat ports on the regions `ports`.
    ///
    /// Every region of `problem.boundary` must be a uniform zero-flux Neumann
    /// row: a port's region because the port's heat replaces that zero, any
    /// other because heat through it would cross no port. The source must be
    /// `ScalarField::Uniform(0.0)` and the conductivity independent of
    /// temperature. `config` is the owner's step configuration, used for
    /// every step and tangent.
    ///
    /// # Errors
    /// [`ConductionError::Config`] naming the refused port, region, source or
    /// material; a solver prepared on another mesh.
    pub fn new(
        solver: BackwardEuler<'m>,
        problem: ConductionProblem<'m>,
        ports: &[&str],
        config: StepConfig,
    ) -> Result<Self, ConductionError> {
        if !std::ptr::eq(solver.mesh(), problem.mesh) {
            return Err(invalid(
                "the solver and the problem must use the same mesh".into(),
            ));
        }
        let boundary = problem.boundary;
        for (name, condition) in boundary.region_names().iter().zip(boundary.conditions()) {
            if !insulated(condition) {
                return Err(invalid(format!(
                    "region {name} is a {} row that is not a uniform zero flux; heat through it would cross no port",
                    condition.tag()
                )));
            }
        }
        if boundary.has_temperature_dependent_surface_conduction() {
            return Err(invalid(
                "temperature-dependent sheet conduction is not wrapped".into(),
            ));
        }
        if !matches!(problem.source, ScalarField::Uniform(s) if *s == 0.0) {
            return Err(invalid(
                "a volumetric source would add heat through no port".into(),
            ));
        }
        let nonlinear = match problem.element_materials {
            Some(materials) => (0..problem.mesh.element_count()).try_fold(false, |any, e| {
                Ok::<_, ConductionError>(any || materials.model_for(e)?.is_temperature_dependent())
            })?,
            None => problem.material.is_temperature_dependent(),
        };
        if nonlinear {
            return Err(invalid(
                "temperature-dependent conductivity is not wrapped".into(),
            ));
        }
        let mut specs = Vec::with_capacity(ports.len());
        let mut patches = Vec::with_capacity(ports.len());
        for (index, &name) in ports.iter().enumerate() {
            if ports[..index].contains(&name) {
                return Err(invalid(format!("port {name} is named twice")));
            }
            let region = boundary
                .region_names()
                .iter()
                .position(|region| region == name)
                .ok_or_else(|| invalid(format!("port {name} names no boundary region")))?;
            let faces: Vec<([usize; 3], f64)> = problem
                .mesh
                .boundary()
                .iter()
                .enumerate()
                .filter(|&(slot, _)| boundary.region_for(slot) == Some(region))
                .map(|(_, face)| (face.vertices.map(|v| v as usize), face.area))
                .collect();
            let area_m2: f64 = faces.iter().map(|&(_, area)| area).sum();
            if !(area_m2.is_finite() && area_m2 > 0.0) {
                return Err(invalid(format!("port {name} has no boundary area")));
            }
            let id =
                StableId::new(name).map_err(|error| invalid(format!("port {name}: {error:?}")))?;
            let spec = PortSpec::heat(id, PortRole::FlowIn)
                .map_err(|error| invalid(format!("port {name}: {error:?}")))?;
            specs.push(spec);
            patches.push(Patch {
                region,
                area_m2,
                faces,
            });
        }
        if specs.is_empty() {
            return Err(invalid("a body needs at least one port".into()));
        }
        let ones = vec![1.0; problem.mesh.vertex_count()];
        let mut weights = vec![0.0; ones.len()];
        solver.capacitance().spmv(&ones, &mut weights);
        Ok(Self {
            solver,
            problem,
            config,
            ports: specs,
            patches,
            weights,
        })
    }

    /// Area [m²] of port `port`'s patch: the flux of a heat rate `Q` on it is
    /// `Q / patch_area_m2(port)`.
    #[must_use]
    pub fn patch_area_m2(&self, port: usize) -> f64 {
        self.patches[port].area_m2
    }

    /// The owner's boundary with each port's heat rate applied on its patch.
    fn boundary_for(&self, heat_w: &[f64]) -> Result<ThermalBoundary, ConductionError> {
        let rows: Vec<(usize, f64)> = self
            .patches
            .iter()
            .zip(heat_w)
            .map(|(patch, &q)| (patch.region, -q / patch.area_m2))
            .collect();
        self.problem.boundary.with_uniform_neumann_additions(&rows)
    }

    /// The owner's tangent at `old` under the held heat rates.
    fn linearize<'a>(
        &'a self,
        old: &[f64],
        dt_s: f64,
        heat_w: &[f64],
        cx: &Cx<'_>,
    ) -> Result<StepLinearization<'a>, ConductionError> {
        self.check_width("inputs", heat_w.len())?;
        let boundary = self.boundary_for(heat_w)?;
        let problem = ConductionProblem {
            boundary: &boundary,
            ..self.problem
        };
        self.solver
            .linearize_step(cx, problem, None, old, dt_s, self.config, None, &[])
    }

    /// Area mean of a nodal field over every patch, into `means`.
    fn patch_means(&self, field: &[f64], means: &mut [f64]) {
        for (mean, patch) in means.iter_mut().zip(&self.patches) {
            let mut sum = 0.0;
            for &(vertices, area) in &patch.faces {
                for vertex in vertices {
                    sum += (area / patch.area_m2 / 3.0) * field[vertex];
                }
            }
            *mean = sum;
        }
    }

    /// The nodal loads [W] of the heat rates `heat_w` spread uniformly over
    /// their patches: `A_f / 3` of each face's flux on each of its vertices.
    fn nodal_loads(&self, heat_w: &[f64], loads: &mut [f64]) {
        for (patch, &q) in self.patches.iter().zip(heat_w) {
            for &(vertices, area) in &patch.faces {
                for vertex in vertices {
                    loads[vertex] += (area / 3.0) * (q / patch.area_m2);
                }
            }
        }
    }

    fn check_width(&self, vector: &str, length: usize) -> Result<(), ConductionError> {
        if length == self.ports.len() {
            Ok(())
        } else {
            Err(invalid(format!(
                "{vector} have {length} entries for {} ports",
                self.ports.len()
            )))
        }
    }
}

fn insulated(condition: &ThermalBc) -> bool {
    matches!(condition, ThermalBc::Neumann { outward_flux: ScalarField::Uniform(flux) } if *flux == 0.0)
}

fn invalid(what: String) -> ConductionError {
    ConductionError::Config {
        parameter: "conduction body",
        what,
    }
}

impl Component for ConductionBody<'_> {
    type State = Vec<f64>;
    type Error = ConductionError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    /// Affine in the heat rates with slopes that depend on the step alone,
    /// but only to the tolerance of the owner's Krylov solve.
    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Vec<f64>) -> Result<f64, ConductionError> {
        if state.len() != self.weights.len() {
            return Err(ConductionError::FieldLength {
                field: "body temperature",
                expected: self.weights.len(),
                found: state.len(),
            });
        }
        let energy: f64 = self.weights.iter().zip(state).map(|(w, t)| w * t).sum();
        if energy.is_finite() {
            Ok(energy)
        } else {
            Err(ConductionError::NonFinite {
                field: "body stored energy",
                bits: energy.to_bits(),
            })
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
    ) -> Result<Proposal, ConductionError> {
        self.check_width("inputs", inputs.len())?;
        self.check_width("outputs", outputs.len())?;
        let boundary = self.boundary_for(inputs)?;
        let problem = ConductionProblem {
            boundary: &boundary,
            ..self.problem
        };
        let step = self
            .solver
            .advance(cx, problem, None, old, dt_s, self.config)?;
        self.patch_means(&step.temperature, outputs);
        // The heat in, summed as the kernel sums the port energies.
        let mut heat_j = 0.0;
        for &q in inputs {
            heat_j += dt_s * q;
        }
        let stored_old_j = self.stored_energy_j(old)?;
        let stored_new_j = self.stored_energy_j(&step.temperature)?;
        // Each stored energy is a sum of `n` positive terms, and the owner's
        // change one of `n` terms no larger: `(n + 16) eps` of both ends
        // covers the rounding of all three.
        let units = (self.weights.len() + 16) as f64 * f64::EPSILON;
        let receipt = Receipt {
            stored_change_j: step.stored_energy_change_j,
            dissipated_j: 0.0,
            numerical_loss_j: heat_j - step.stored_energy_change_j,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: units * (stored_old_j.abs() + stored_new_j.abs() + heat_j.abs()),
            residual: step.relative_residual,
            tolerance: self.config.linear.tolerance,
        };
        new.clear();
        new.extend_from_slice(&step.temperature);
        Ok(Proposal::Step(receipt))
    }
}

impl PortTangent for ConductionBody<'_> {
    fn jvp(
        &self,
        old: &Vec<f64>,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), ConductionError> {
        self.check_width("tangent inputs", d_inputs.len())?;
        self.check_width("tangent outputs", d_outputs.len())?;
        let step = self.linearize(old, dt_s, inputs, cx)?;
        let mut direction = step.zero_direction();
        self.nodal_loads(d_inputs, &mut direction.nodal_load_w);
        let change = step.apply(cx, &direction)?;
        self.patch_means(&change.temperature_k, d_outputs);
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
    ) -> Result<(), ConductionError> {
        self.check_width("adjoint outputs", bar_outputs.len())?;
        self.check_width("adjoint inputs", bar_inputs.len())?;
        let step = self.linearize(old, dt_s, inputs, cx)?;
        // The patch means are the transpose of the nodal loads up to each
        // patch's area: weights `A_f / (3 A_p)` on its vertices.
        let mut weights = vec![0.0; self.weights.len()];
        for (patch, &bar) in self.patches.iter().zip(bar_outputs) {
            for &(vertices, area) in &patch.faces {
                for vertex in vertices {
                    weights[vertex] += (area / patch.area_m2 / 3.0) * bar;
                }
            }
        }
        let gradient = step.pullback(cx, &weights, &[], &[])?;
        for (bar, patch) in bar_inputs.iter_mut().zip(&self.patches) {
            let mut sum = 0.0;
            for &(vertices, area) in &patch.faces {
                for vertex in vertices {
                    sum += (area / 3.0) * (gradient.nodal_load[vertex] / patch.area_m2);
                }
            }
            *bar = sum;
        }
        Ok(())
    }
}
