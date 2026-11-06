//! Exact held-current motion through the existing quadratic-flow owner.
//!
//! The original material, consistent rod inertia, reciprocal coupling and
//! charge remain unchanged. Cold preparation caches the full coupled matrix
//! exponential and voltage integral. Warm stepping uses fixed caller-owned
//! buffers, without a Newton solve, an exponential or a material clone.
//! "Exact" means the linear held-input model up to floating-point error; it
//! does not remove the rod's spatial reduction or certify a real-time deadline.

use super::{PiezoelectricBody, PiezoelectricCouplingError, Result, invalid, poll, positive};
use fs_phs::PreparedStepError;
use fs_phs::linear_exact::prepared::PreparedQuadraticFlow;
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget, split_quadratic};

/// Independent electrical and energetic acceptance for an exact-flow window.
#[derive(Clone, Copy, Debug)]
pub struct ExactPiezoelectricControls {
    /// Absolute combined stored-energy minus source-work tolerance [J].
    pub energy_tolerance_j: f64,
    /// Absolute free-charge increment minus integrated current tolerance [C].
    pub charge_tolerance_c: f64,
}

/// Physical observables integrated on the same continuous linear trajectory.
#[derive(Clone, Copy, Debug)]
pub struct ExactPiezoelectricStep {
    /// Fixed preparation interval [s].
    pub duration_s: f64,
    /// Actual energy change, including mechanical and electric storage [J].
    pub energy_change_j: f64,
    /// Held mechanical force times the actual extension increment [J].
    pub mechanical_supply_j: f64,
    /// Held current times the integrated terminal voltage [J].
    pub electrical_supply_j: f64,
    /// Stored-energy change minus the two independently evaluated supplies [J].
    pub energy_residual_j: f64,
    /// Free-charge increment minus current times duration [C].
    pub charge_residual_c: f64,
    /// Time-mean terminal voltage [V], not an endpoint average.
    pub mean_voltage_v: f64,
    /// Time-mean tip velocity [m/s].
    pub mean_velocity_m_s: f64,
    /// Conservative computed full-window energy bound [J]. This is a numerical
    /// norm estimate, not an outward-rounded interval certificate.
    pub interval_energy_bound_j: f64,
}

/// A reusable fixed-clock, source-bound sensor/actuator with current control.
/// Zero current gives an isolated electrode, retaining its original charge.
/// Every accepted step preserves the entire material/mechanical state. A
/// conservative whole-window energy envelope also checks the small-signal
/// strain/field limits, so an excursion between equal endpoints cannot hide.
/// Cancellation or refusal leaves the body unchanged; only scratch changes.
pub struct PreparedPiezoelectricCurrent {
    body: PiezoelectricBody,
    flow: PreparedQuadraticFlow,
    duration_s: f64,
    augmented: [f64; 6],
    balanced: [f64; 3],
    endpoint: [f64; 3],
    integrals: [f64; 1],
}

impl PreparedPiezoelectricCurrent {
    /// Cache the complete three-state flow under explicit cold work/memory
    /// budgets. The original source law supplies exact quadratic metadata.
    pub fn new(
        body: PiezoelectricBody,
        duration_s: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        if !positive(duration_s) {
            return Err(invalid("exact piezoelectric duration must be finite and positive"));
        }
        body.stored_energy_j()?;
        let system = body.rod.port_hamiltonian()?;
        let split = split_quadratic(&system, split_budget)?;
        let capacitance = body.rod.axially_clamped_capacitance_f();
        let voltage = vec![
            -body.rod.coupling_c_m() / capacitance,
            0.0,
            1.0 / capacitance,
        ];
        let flow = PreparedQuadraticFlow::new_controlled(
            &system,
            &split,
            duration_s,
            &[],
            &[voltage],
            dense_budget,
            &mut cancelled,
        )
        .map_err(|error| match error {
            PreparedStepError::Cancelled => PiezoelectricCouplingError::Cancelled,
            PreparedStepError::Solver(error) => PiezoelectricCouplingError::Native(error),
        })?;
        poll(&mut cancelled)?;
        Ok(Self {
            body,
            flow,
            duration_s,
            augmented: [0.0; 6],
            balanced: [0.0; 3],
            endpoint: [0.0; 3],
            integrals: [0.0; 1],
        })
    }

    /// Accepted original body and material; observation does not advance it.
    pub const fn body(&self) -> &PiezoelectricBody {
        &self.body
    }

    /// Prepared interval [s]. Changing it requires new preparation.
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Return the physical state to another compatible boundary owner.
    pub fn into_body(self) -> PiezoelectricBody {
        self.body
    }

    /// Advance a held force/current and publish only after charge, work,
    /// constitutive-domain and cancellation checks pass. No warm allocation is
    /// performed by this owner or its prepared numerical kernel.
    pub fn step(
        &mut self,
        force_n: f64,
        current_a: f64,
        controls: ExactPiezoelectricControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExactPiezoelectricStep> {
        poll(&mut cancelled)?;
        if !force_n.is_finite()
            || !current_a.is_finite()
            || !positive(controls.energy_tolerance_j)
            || !positive(controls.charge_tolerance_c)
        {
            return Err(invalid(
                "finite held piezoelectric sources and positive work/charge tolerances required",
            ));
        }
        let state = self.body.state;
        let inputs = [force_n, current_a];
        let norm = self.flow.interval_coordinate_norm_bound(&state, &inputs)?;
        // H = k_E*x^2/2 + C*V^2/2 + p^2/(2*m). Bounding the original
        // positive storage bounds strain and field simultaneously.
        let rod = &self.body.rod;
        let extension_bound = norm / fs_math::det::sqrt(rod.short_circuit_stiffness_n_m());
        let voltage_bound = norm / fs_math::det::sqrt(rod.axially_clamped_capacitance_f());
        let interval_energy_bound_j = 0.5 * norm * norm;
        if !extension_bound.is_finite()
            || !voltage_bound.is_finite()
            || !interval_energy_bound_j.is_finite()
            || extension_bound / rod.length_m > rod.limits.max_abs_strain
            || voltage_bound / rod.length_m > rod.limits.max_abs_electric_field_v_m
        {
            return Err(invalid(
                "complete exact-flow window exceeds the declared piezoelectric small-signal envelope",
            ));
        }
        self.flow.evaluate_into(
            &state,
            &inputs,
            &mut self.augmented,
            &mut self.balanced,
            &mut self.endpoint,
            &mut self.integrals,
        )?;
        let before = rod.energy(&state)?;
        let after = rod.energy(&self.endpoint)?;
        let extension = self.endpoint[0] - state[0];
        let charge = self.endpoint[2] - state[2];
        let mechanical_supply_j = force_n * extension;
        let electrical_supply_j = current_a * self.integrals[0];
        let energy_change_j = after - before;
        let energy_residual_j = energy_change_j - mechanical_supply_j - electrical_supply_j;
        let charge_residual_c = charge - self.duration_s * current_a;
        let mean_voltage_v = self.integrals[0] / self.duration_s;
        let mean_velocity_m_s = extension / self.duration_s;
        if [
            energy_change_j,
            mechanical_supply_j,
            electrical_supply_j,
            energy_residual_j,
            charge_residual_c,
            mean_voltage_v,
            mean_velocity_m_s,
        ]
        .iter()
        .any(|x| !x.is_finite())
        {
            return Err(invalid("exact piezoelectric source work or charge is unrepresentable"));
        }
        if charge_residual_c.abs() > controls.charge_tolerance_c {
            return Err(invalid("exact piezoelectric endpoint failed free-charge conservation"));
        }
        if energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(PiezoelectricCouplingError::Balance {
                residual_j: energy_residual_j,
                tolerance_j: controls.energy_tolerance_j,
            });
        }
        poll(&mut cancelled)?;
        self.body.state = self.endpoint;
        Ok(ExactPiezoelectricStep {
            duration_s: self.duration_s,
            energy_change_j,
            mechanical_supply_j,
            electrical_supply_j,
            energy_residual_j,
            charge_residual_c,
            mean_voltage_v,
            mean_velocity_m_s,
            interval_energy_bound_j,
        })
    }
}
