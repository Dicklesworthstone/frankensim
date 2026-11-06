//! Exact fixed-voltage constraint elimination on the original rod.
//! Q is reconstructed from its original anchor, never projected or integrated
//! separately. Both mechanical coordinates remain in the native quadratic-flow owner.

use super::{
    DenseExponentialBudget, ExactPiezoelectricControls, ExactPiezoelectricStep, LinearSplitBudget,
    PiezoelectricBody, PiezoelectricCouplingError, PreparedQuadraticFlow, PreparedStepError,
    Result, invalid, poll, positive, split_quadratic,
};
use fs_phs::{PortHamiltonian, QuadraticStorage};

/// Reusable exact motion under the original body's fixed terminal voltage.
/// The complete rod/material and initial [extension, momentum, charge] survive
/// preparation. Only the algebraic voltage constraint is eliminated: the native
/// two-state generator retains stiffness k_E, consistent mass and force
/// F + theta*V. No physical mechanical mode or independent gain is introduced.
///
/// Voltage is read from the original body at preparation and cannot be changed
/// here. Switching to an incompatible ideal clamp requires an actual circuit or
/// current-drive model; it must not reset charge or discard stored energy.
pub struct PreparedPiezoelectricVoltage {
    body: PiezoelectricBody,
    flow: PreparedQuadraticFlow,
    duration_s: f64,
    voltage_v: f64,
    anchor_extension_m: f64,
    anchor_charge_c: f64,
    augmented: [f64; 4],
    balanced: [f64; 2],
    endpoint: [f64; 2],
    integrals: [f64; 0],
}

impl PreparedPiezoelectricVoltage {
    /// Prepare the exact constrained oscillator under explicit cold budgets.
    /// The voltage is inherited from the supplied state, without projection.
    pub fn new(
        body: PiezoelectricBody,
        duration_s: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        if !positive(duration_s) {
            return Err(invalid(
                "exact piezoelectric duration must be finite and positive",
            ));
        }
        body.stored_energy_j()?;
        let voltage_v = body.voltage_v()?;
        if !voltage_v.is_finite() {
            return Err(invalid(
                "original fixed piezoelectric voltage is unrepresentable",
            ));
        }
        let rod = &body.rod;
        let system = PortHamiltonian::new(
            2,
            1,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0; 4],
            vec![0.0, 1.0],
            Box::new(QuadraticStorage::new(
                vec![
                    rod.short_circuit_stiffness_n_m(),
                    0.0,
                    0.0,
                    1.0 / rod.modal_mass_kg(),
                ],
                2,
            )?),
        )?;
        let split = split_quadratic(&system, split_budget)?;
        let flow = PreparedQuadraticFlow::new_controlled(
            &system,
            &split,
            duration_s,
            &[],
            &[],
            dense_budget,
            &mut cancelled,
        )
        .map_err(|error| match error {
            PreparedStepError::Cancelled => PiezoelectricCouplingError::Cancelled,
            PreparedStepError::Solver(error) => PiezoelectricCouplingError::Native(error),
        })?;
        poll(&mut cancelled)?;
        Ok(Self {
            anchor_extension_m: body.extension_m(),
            anchor_charge_c: body.charge_c(),
            body,
            flow,
            duration_s,
            voltage_v,
            augmented: [0.0; 4],
            balanced: [0.0; 2],
            endpoint: [0.0; 2],
            integrals: [],
        })
    }

    /// Accepted complete original material, mechanical state and electrode charge.
    pub const fn body(&self) -> &PiezoelectricBody {
        &self.body
    }

    /// Original terminal voltage, fixed for this prepared owner's lifetime [V].
    pub const fn voltage_v(&self) -> f64 {
        self.voltage_v
    }

    /// Prepared interval [s]. Changing it requires new preparation.
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Return the complete state to another compatible electrical boundary owner.
    pub fn into_body(self) -> PiezoelectricBody {
        self.body
    }

    /// Advance a held mechanical force using fixed warm buffers. Before
    /// publication, independently admit the original rod's full energy, both
    /// source works, induced charge and voltage constraint at both endpoints.
    /// The native mechanical norm bounds the entire interval's strain; the
    /// constant voltage, including its explicit tolerance, bounds electric field.
    /// Cancellation or refusal changes only scratch, never the accepted body.
    pub fn step(
        &mut self,
        force_n: f64,
        voltage_tolerance_v: f64,
        controls: ExactPiezoelectricControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExactPiezoelectricStep> {
        poll(&mut cancelled)?;
        if !force_n.is_finite()
            || !positive(voltage_tolerance_v)
            || !positive(controls.energy_tolerance_j)
            || !positive(controls.charge_tolerance_c)
        {
            return Err(invalid(
                "finite force and positive fixed-voltage, work and charge tolerances required",
            ));
        }
        let before_voltage = self.body.voltage_v()?;
        if !before_voltage.is_finite()
            || (before_voltage - self.voltage_v).abs() > voltage_tolerance_v
        {
            return Err(invalid("accepted body violates its original fixed voltage"));
        }
        let rod = &self.body.rod;
        let state = self.body.state;
        let mechanical = [state[0], state[1]];
        let theta = rod.coupling_c_m();
        let inputs = [force_n + theta * self.voltage_v];
        let norm = self
            .flow
            .interval_coordinate_norm_bound(&mechanical, &inputs)?;
        let extension_bound = norm / fs_math::det::sqrt(rod.short_circuit_stiffness_n_m());
        let voltage_bound = self.voltage_v.abs() + voltage_tolerance_v;
        let interval_energy_bound_j = 0.5 * norm * norm
            + 0.5 * rod.axially_clamped_capacitance_f() * voltage_bound * voltage_bound;
        if !extension_bound.is_finite()
            || !voltage_bound.is_finite()
            || !interval_energy_bound_j.is_finite()
            || extension_bound / rod.length_m > rod.limits.max_abs_strain
            || voltage_bound / rod.length_m > rod.limits.max_abs_electric_field_v_m
        {
            return Err(invalid(
                "complete exact fixed-voltage window exceeds the piezoelectric small-signal envelope",
            ));
        }
        self.flow.evaluate_into(
            &mechanical,
            &inputs,
            &mut self.augmented,
            &mut self.balanced,
            &mut self.endpoint,
            &mut self.integrals,
        )?;
        let next = [
            self.endpoint[0],
            self.endpoint[1],
            self.anchor_charge_c + theta * (self.endpoint[0] - self.anchor_extension_m),
        ];
        let before = rod.energy(&state)?;
        let after = rod.energy(&next)?;
        let after_voltage = rod.length_m * rod.response(next[0], next[2])?.electric_field_v_m;
        if !after_voltage.is_finite()
            || (after_voltage - self.voltage_v).abs() > voltage_tolerance_v
        {
            return Err(invalid(
                "exact endpoint violated the original fixed-voltage constraint",
            ));
        }
        let extension = next[0] - state[0];
        let charge = next[2] - state[2];
        let energy_change_j = after - before;
        let mechanical_supply_j = force_n * extension;
        let electrical_supply_j = self.voltage_v * charge;
        let energy_residual_j = energy_change_j - mechanical_supply_j - electrical_supply_j;
        let charge_residual_c = charge - theta * extension;
        let mean_velocity_m_s = extension / self.duration_s;
        if [
            energy_change_j,
            mechanical_supply_j,
            electrical_supply_j,
            energy_residual_j,
            charge_residual_c,
            mean_velocity_m_s,
        ]
        .iter()
        .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "exact fixed-voltage work or induced charge is unrepresentable",
            ));
        }
        if charge_residual_c.abs() > controls.charge_tolerance_c {
            return Err(invalid(
                "exact fixed-voltage step failed induced-charge conservation",
            ));
        }
        if energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(PiezoelectricCouplingError::Balance {
                residual_j: energy_residual_j,
                tolerance_j: controls.energy_tolerance_j,
            });
        }
        poll(&mut cancelled)?;
        self.body.state = next;
        Ok(ExactPiezoelectricStep {
            duration_s: self.duration_s,
            energy_change_j,
            mechanical_supply_j,
            electrical_supply_j,
            energy_residual_j,
            charge_residual_c,
            mean_voltage_v: self.voltage_v,
            mean_velocity_m_s,
            interval_energy_bound_j,
        })
    }
}
