//! Exact sparse linearization of the original layered backward-Euler residual.
//!
//! Unknowns are cell energy increments, in the original owner's CellId order.
//! For R = delta_E - dt Q(T(E)), J = I - dt (dQ/dT) diag(1/C).
//! Capacity is the original per-cell C(T); flux tangents differentiate the same
//! integrated conductivity/contact equations as the trial coefficients. There
//! is no derivative of an iterative solver, finite-difference material query,
//! frozen secant-conductance approximation, or second residual/acceptance owner.
//!
//! This operator admits a whole isolated layered-body graph with its original
//! fixed cell sources, linear convection and reservoirs. Larger coupled graphs
//! continue to use `BoundLayeredThermalBody::step_coefficients` and their own
//! complete Jacobian assembly. A tangent is local to one exact candidate; a
//! transpose action alone does not supply a trajectory adjoint or error bound.

use super::BoundLayeredThermalBody;
use crate::thermal::conduction::{LayerConductionConfig, LayerConductionResponse};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation,
    backward_euler::{BackwardEulerStep, BackwardEulerTrial, StepCoefficients},
};
use std::cell::RefCell;

/// One frozen original thermal step and its source-qualified material binding.
/// This supplies residuals and derivatives to a caller's nonlinear solver; it
/// never advances or accepts a physical state on its own.
pub struct LayeredThermalStep<'s, 'm, 'a> {
    step: BackwardEulerStep<'s, 'm>,
    binding: &'s BoundLayeredThermalBody<'a>,
    external_conductances: Vec<f64>,
    config: LayerConductionConfig,
}

impl<'s, 'm, 'a> LayeredThermalStep<'s, 'm, 'a> {
    pub fn new(
        binding: &'s BoundLayeredThermalBody<'a>,
        initial: &'s DepositionSimulation<'m>,
        target_s: f64,
        config: LayerConductionConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        let n = binding.body.cell_count();
        if initial.cells().len() != n || initial.links().len() != n - 1 {
            return Err(bad(
                "implicit layered operator requires exactly its own cells and links",
            ));
        }
        // The existing adapter owns source, geometry, material-object and link
        // admission. Do not maintain a second weaker binding check here.
        binding.step_coefficients(initial, config, &mut checkpoint)?;
        let step = BackwardEulerStep::new(initial, target_s, &mut checkpoint)?;
        let mut external_conductances = Vec::with_capacity(n);
        for &id in step.cell_ids() {
            checkpoint()?;
            let boundary = initial.cells()[&id].boundary;
            // Use the original owner's actual exposed area, including authored
            // contact patches and reservoir coverage. Boundary events within
            // this interval are refused by BackwardEulerStep::new.
            let g = initial.exposed_area_m2(id)? * boundary.htc_w_m2_k
                + boundary.reservoir.map_or(0.0, |r| r.conductance_w_k);
            if !g.is_finite() || g < 0.0 {
                return Err(DepositionError::Arithmetic);
            }
            external_conductances.push(g);
        }
        checkpoint()?;
        Ok(Self {
            step,
            binding,
            external_conductances,
            config,
        })
    }

    pub fn cell_ids(&self) -> &[CellId] {
        self.step.cell_ids()
    }
    pub const fn duration_s(&self) -> f64 {
        self.step.duration_s()
    }
    pub const fn target_s(&self) -> f64 {
        self.step.target_s()
    }

    /// Evaluate the original physical residual with a tangent of THAT candidate.
    pub fn evaluate(
        &self,
        energy_increments_j: &[f64],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, LayeredThermalLinearization>, DepositionError> {
        self.evaluate_inner(energy_increments_j, false, checkpoint)
    }

    /// Same original residual, with the original owner's optional entropy audit.
    pub fn evaluate_with_entropy(
        &self,
        energy_increments_j: &[f64],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, LayeredThermalLinearization>, DepositionError> {
        self.evaluate_inner(energy_increments_j, true, checkpoint)
    }

    fn evaluate_inner(
        &self,
        energy_increments_j: &[f64],
        with_entropy: bool,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, LayeredThermalLinearization>, DepositionError> {
        // The owner invokes the two closures sequentially. Share one bounded
        // checkpoint without weakening cancellation inside material evaluation.
        let checkpoint = RefCell::new(&mut checkpoint);
        let resolve = |candidate: &DepositionSimulation<'m>| {
            let StepCoefficients {
                sources_w,
                conductances_w_k,
                payload: links,
            } = self
                .binding
                .step_coefficients(candidate, self.config, || (*checkpoint.borrow_mut())())?;
            let linearization =
                self.linearize(candidate, links, || (*checkpoint.borrow_mut())())?;
            Ok(StepCoefficients {
                sources_w,
                conductances_w_k,
                payload: linearization,
            })
        };
        let poll = || (*checkpoint.borrow_mut())();
        if with_entropy {
            self.step
                .evaluate_with_entropy(energy_increments_j, resolve, poll)
        } else {
            self.step.evaluate(energy_increments_j, resolve, poll)
        }
    }

    fn linearize(
        &self,
        candidate: &DepositionSimulation<'_>,
        links: Vec<LayerConductionResponse>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<LayeredThermalLinearization, DepositionError> {
        let n = self.cell_ids().len();
        let dt = self.duration_s();
        let mut temperatures = Vec::with_capacity(n);
        let mut capacities = Vec::with_capacity(n);
        let mut diagonal = Vec::with_capacity(n);
        for (i, &id) in self.cell_ids().iter().enumerate() {
            checkpoint()?;
            let t = candidate.temperature_k(id)?;
            let layer = self.binding.body.cells[i].layer;
            let capacity = self.binding.body.storages[layer]
                .heat_capacity_j_k(t)
                .map_err(|_| DepositionError::MaterialDomain)?;
            if !capacity.is_finite() || capacity <= 0.0 {
                return Err(DepositionError::Arithmetic);
            }
            let boundary = nonnegative(dt * (self.external_conductances[i] / capacity))?;
            diagonal.push(positive(1.0 + boundary)?);
            temperatures.push(t);
            capacities.push(capacity);
        }
        let mut lower = Vec::with_capacity(n - 1);
        let mut upper = Vec::with_capacity(n - 1);
        for (i, link) in links.iter().enumerate() {
            checkpoint()?;
            let [dq_dta, dq_dtb] = link.endpoint_derivative_w_k;
            let left = positive(dt * (dq_dta / capacities[i]))?;
            let right = positive(dt * (-dq_dtb / capacities[i + 1]))?;
            diagonal[i] = positive(diagonal[i] + left)?;
            diagonal[i + 1] = positive(diagonal[i + 1] + right)?;
            // R_left += dt Q and R_right -= dt Q. Columns differentiate
            // THEIR cell energies, so unequal capacities make J nonsymmetric.
            lower.push(-left);
            upper.push(-right);
        }
        checkpoint()?;
        Ok(LayeredThermalLinearization {
            diagonal,
            lower,
            upper,
            temperatures,
            capacities,
            links,
        })
    }
}

/// Tridiagonal derivative of one exact original thermal trial in joule inputs.
/// Assembly and retention require O(cells) storage. Each apply/transpose/diagonal
/// inverse action is O(cells), allocates nothing, and publishes no output until
/// all arithmetic and cancellation checks pass. Positive diagonal and finite
/// coefficients do not constitute a conditioning or convergence certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct LayeredThermalLinearization {
    diagonal: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    temperatures: Vec<f64>,
    capacities: Vec<f64>,
    links: Vec<LayerConductionResponse>,
}

#[derive(Clone, Copy)]
enum Action {
    Forward,
    Transpose,
    InverseDiagonal,
}

impl LayeredThermalLinearization {
    pub fn dimension(&self) -> usize {
        self.diagonal.len()
    }
    pub fn diagonal(&self) -> &[f64] {
        &self.diagonal
    }
    pub fn lower_diagonal(&self) -> &[f64] {
        &self.lower
    }
    pub fn upper_diagonal(&self) -> &[f64] {
        &self.upper
    }
    pub fn cell_temperatures_k(&self) -> &[f64] {
        &self.temperatures
    }
    pub fn heat_capacities_j_k(&self) -> &[f64] {
        &self.capacities
    }
    pub fn link_responses(&self) -> &[LayerConductionResponse] {
        &self.links
    }

    pub fn apply_into(
        &self,
        direction: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        self.apply(direction, output, Action::Forward, checkpoint)
    }

    /// The exact transpose, including unequal capacities and nonsymmetric flux
    /// tangents. It never substitutes the forward action for a symmetry guess.
    pub fn apply_transpose_into(
        &self,
        direction: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        self.apply(direction, output, Action::Transpose, checkpoint)
    }

    /// Jacobi inverse action at this frozen candidate, not a system solve.
    pub fn apply_inverse_diagonal_into(
        &self,
        direction: &[f64],
        output: &mut [f64],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        self.apply(direction, output, Action::InverseDiagonal, checkpoint)
    }

    fn apply(
        &self,
        direction: &[f64],
        output: &mut [f64],
        action: Action,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        checkpoint()?;
        if direction.len() != self.dimension() || output.len() != self.dimension() {
            return Err(bad(
                "thermal tangent vector length differs from its original cells",
            ));
        }
        for (i, value) in direction.iter().enumerate() {
            if i % 64 == 0 {
                checkpoint()?;
            }
            if !value.is_finite() {
                return Err(DepositionError::Arithmetic);
            }
        }
        // Validate the bounded action before publishing. Repeating the same
        // immutable row expression avoids a heap scratch buffer. The final
        // copy is bounded by MAX_LAYERED_THERMAL_CELLS (512).
        for i in 0..self.dimension() {
            if i % 64 == 0 {
                checkpoint()?;
            }
            if !self.row(i, direction, action).is_finite() {
                return Err(DepositionError::Arithmetic);
            }
        }
        checkpoint()?;
        for (i, value) in output.iter_mut().enumerate() {
            *value = self.row(i, direction, action);
        }
        Ok(())
    }

    fn row(&self, i: usize, x: &[f64], action: Action) -> f64 {
        let (lower, upper) = match action {
            Action::Forward => (&self.lower, &self.upper),
            Action::Transpose => (&self.upper, &self.lower),
            Action::InverseDiagonal => return x[i] / self.diagonal[i],
        };
        let mut value = self.diagonal[i] * x[i];
        if i > 0 {
            value += lower[i - 1] * x[i - 1];
        }
        if i + 1 < self.dimension() {
            value += upper[i] * x[i + 1];
        }
        value
    }
}

fn nonnegative(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
fn positive(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
fn bad(reason: &'static str) -> DepositionError {
    DepositionError::InvalidInput(reason)
}
