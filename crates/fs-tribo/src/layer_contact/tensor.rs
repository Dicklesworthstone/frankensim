//! Full two-axis recoverable and viscous shear on the original layer law.
//!
//! Each equilibrium/Maxwell spring retains its own tensor; different arms
//! need not share principal directions. The plastic branch is explicitly the
//! original isotropic Bingham law, not an inferred anisotropic yield surface.
//! The held-coefficient midpoint balance is unchanged. One bounded scalar
//! plastic compliance solves the tensor return without dropping cross terms.
//! Kelvin viscosity acts on the recoverable shear rate in parallel with those
//! springs, before the series plastic slider. It has no stored state or energy.

use super::*;
use crate::tangential::transport::{TangentialHistoryRotation, TangentialTransportConfig};

/// A positive-semidefinite shear stiffness, represented by its original Gram
/// rows. This admits rank-deficient memory without eigenvalue clipping.
#[derive(Clone, Debug, PartialEq)]
pub struct ShearStiffness {
    rows: Vec<[f64; 2]>,
    matrix: [[f64; 2]; 2],
}
impl ShearStiffness {
    /// At most six physical factor rows, in sqrt(N/m). Zero rows are allowed.
    pub fn from_gram_rows(rows: Vec<[f64; 2]>) -> Result<Self, TriboError> {
        if rows.is_empty() || rows.len() > 6 {
            return Err(invalid("layer_tensor_factor_rows"));
        }
        let mut matrix = [[0.0; 2]; 2];
        for row in &rows {
            for value in row {
                finite(*value, "layer_tensor_factor")?;
            }
            for i in 0..2 {
                for j in i..2 {
                    matrix[i][j] = checked_add(
                        matrix[i][j],
                        checked_mul(row[i], row[j], "layer_tensor_stiffness")?,
                        "layer_tensor_stiffness",
                    )?;
                }
            }
        }
        matrix[1][0] = matrix[0][1];
        Ok(Self { rows, matrix })
    }
    /// Strict SPD matrix in N/m. Symmetry and positive pivots are checked;
    /// roundoff-sized negative pivots are never promoted into physical data.
    pub fn positive_definite(matrix: [[f64; 2]; 2]) -> Result<Self, TriboError> {
        for value in matrix.iter().flatten() {
            finite(*value, "layer_tensor_stiffness")?;
        }
        if matrix[0][1] != matrix[1][0] || matrix[0][0] <= 0.0 {
            return Err(invalid("layer_tensor_positive_definite"));
        }
        let first = det::sqrt(matrix[0][0]);
        let cross = matrix[0][1] / first;
        let remainder = matrix[1][1] - cross * cross;
        positive_finite(remainder, "layer_tensor_positive_pivot")?;
        Self::from_gram_rows(vec![[first, cross], [0.0, det::sqrt(remainder)]])
    }
    pub const fn matrix_n_m(&self) -> [[f64; 2]; 2] {
        self.matrix
    }
    pub fn gram_rows(&self) -> &[[f64; 2]] {
        &self.rows
    }
    fn quadratic(&self, vector: [f64; 2]) -> f64 {
        self.rows
            .iter()
            .map(|r| {
                let v = dot(*r, vector);
                v * v
            })
            .sum()
    }
    fn rotated(&self, rotation: TangentialHistoryRotation) -> Result<Self, TriboError> {
        Self::from_gram_rows(
            self.rows
                .iter()
                .map(|&row| rotation.apply(row))
                .collect::<Result<Vec<_>, _>>()?,
        )
    }
}

/// Positive-semidefinite Kelvin shear viscosity in N s/m, retaining the original
/// Gram rows in sqrt(N s/m). This is a rate law, never a stored spring.
#[derive(Clone, Debug, PartialEq)]
pub struct ShearViscosity {
    factor: ShearStiffness,
}
impl ShearViscosity {
    /// At most six physical rows; zero and rank-deficient viscosity are valid.
    pub fn from_gram_rows(rows: Vec<[f64; 2]>) -> Result<Self, TriboError> {
        Ok(Self {
            factor: ShearStiffness::from_gram_rows(rows)?,
        })
    }
    pub const fn matrix_n_s_m(&self) -> [[f64; 2]; 2] {
        self.factor.matrix
    }
    pub fn gram_rows(&self) -> &[[f64; 2]] {
        &self.factor.rows
    }
    fn rotated(&self, rotation: TangentialHistoryRotation) -> Result<Self, TriboError> {
        Ok(Self {
            factor: self.factor.rotated(rotation)?,
        })
    }
}

/// One original Maxwell tensor and relaxation time; branch order is retained.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorLayerRelaxationArm {
    pub stiffness: ShearStiffness,
    pub relaxation_time_s: f64,
}

/// Complete affine derivatives of the retained equilibrium and Maxwell
/// stiffness tensors [N/(m K)]. Signed cross terms and indefinite derivatives
/// are physical possibilities; only the actual stored tensors must be PSD.
/// These numerical inputs do not themselves establish material-source authority.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorLayerThermalSlopes {
    pub equilibrium_n_m_k: [[f64; 2]; 2],
    pub relaxation_n_m_k: Vec<[[f64; 2]; 2]>,
}
#[derive(Clone, Debug, PartialEq)]
pub(super) struct TensorLayerPotential {
    equilibrium: ShearStiffness,
    pub(super) arms: Vec<TensorLayerRelaxationArm>,
    viscosity: Option<ShearViscosity>,
}

/// A joint zero-work re-expression of the constitutive tensors and all history.
/// The source/geometry owner must admit the corresponding physical orientation.
pub struct TensorLayerFrameTransfer {
    pub law: ThirdBodyLayerContact,
    pub state: ThirdBodyLayerState,
    pub stored_energy_residual_j: f64,
}

impl ThirdBodyLayerContact {
    /// Complete tensor potential, with the original scalar plastic flow law.
    /// The positive equilibrium tensor and all PSD arms are explicit inputs.
    /// This constructor carries no measured-material or thermal-law authority.
    #[allow(clippy::too_many_arguments)]
    pub fn new_tensor(
        equilibrium: ShearStiffness,
        arms: Vec<TensorLayerRelaxationArm>,
        yield_force_n: f64,
        plastic_resistance_n_s_m: f64,
        temperature_kelvin: f64,
        maximum_elastic_slip_m: f64,
    ) -> Result<Self, TriboError> {
        ShearStiffness::positive_definite(equilibrium.matrix)?;
        if arms.len() > MAX_LAYER_RELAXATION_ARMS {
            return Err(invalid("layer_relaxation_count"));
        }
        let mut total = equilibrium.matrix;
        for arm in &arms {
            positive_finite(arm.relaxation_time_s, "layer_relaxation_time_s")?;
            total = add(total, arm.stiffness.matrix, 1.0)?;
        }
        let mut law = Self::new(
            equilibrium.matrix[0][0],
            yield_force_n,
            plastic_resistance_n_s_m,
            temperature_kelvin,
            maximum_elastic_slip_m,
        )?;
        law.tensor = Some(TensorLayerPotential {
            equilibrium,
            arms,
            viscosity: None,
        });
        Ok(law)
    }
    /// Add the original Kelvin rate tensor to this tensor material. Its strain
    /// increment is total shear minus plastic shear; applying it to total slip
    /// would change the network topology and double count plastic dissipation.
    pub fn with_kelvin_viscosity(mut self, viscosity: ShearViscosity) -> Result<Self, TriboError> {
        let tensor = self
            .tensor
            .as_mut()
            .ok_or_else(|| invalid("layer_kelvin_requires_tensor_owner"))?;
        if tensor.viscosity.is_some() {
            return Err(invalid("layer_kelvin_already_bound"));
        }
        tensor.viscosity = Some(viscosity);
        Ok(self)
    }
    /// Complete rate tensor, absent for a law without Kelvin viscosity.
    pub fn tensor_kelvin_viscosity(&self) -> Option<&ShearViscosity> {
        self.tensor
            .as_ref()
            .and_then(|value| value.viscosity.as_ref())
    }
    /// Complete equilibrium tensor in the checkpoint's tangent coordinates.
    pub fn stiffness_tensor_n_m(&self) -> [[f64; 2]; 2] {
        self.tensor.as_ref().map_or(
            [[self.stiffness_n_m, 0.0], [0.0, self.stiffness_n_m]],
            |value| value.equilibrium.matrix,
        )
    }
    /// Tensor arms, absent for the separately preserved isotropic constructor.
    pub fn tensor_relaxation(&self) -> Option<&[TensorLayerRelaxationArm]> {
        self.tensor.as_ref().map(|value| value.arms.as_slice())
    }

    /// S=-e^T K' e/2-sum(r_j^T K'_j r_j/2), U=Psi+T*S, using the
    /// SAME retained physical spring strains as the tensor force law. Every
    /// branch derivative is required in source order, including a zero slope.
    /// The coupled owner admits the full affine interval and accounts for both
    /// coefficient work and delta(T*S). Neither Kelvin viscosity nor plastic
    /// kinetics has stored free energy, even when either varies with temperature.
    pub fn tensor_thermodynamics(&self, state: &ThirdBodyLayerState,
        slopes: &TensorLayerThermalSlopes)
        -> Result<ThirdBodyLayerThermodynamics, TriboError>
    {
        let potential = self.tensor_state(state)?;
        if slopes.relaxation_n_m_k.len() != potential.arms.len() {
            return Err(invalid("layer_thermal_memory_derivative_count"));
        }
        let contribution = |matrix: [[f64; 2]; 2], slip: [f64; 2]| {
            if matrix[0][1] != matrix[1][0]
                || matrix.iter().flatten().any(|v| !v.is_finite()) {
                return Err(invalid("layer_tensor_temperature_derivative"));
            }
            let mut entropy = 0.0;
            for i in 0..2 {
                for j in 0..2 {
                    let term = checked_mul(checked_mul(-0.5*matrix[i][j], slip[i],
                        "layer_mechanical_entropy_j_k")?, slip[j], "layer_mechanical_entropy_j_k")?;
                    entropy = checked_add(entropy, term, "layer_mechanical_entropy_j_k")?;
                }
            }
            Ok(entropy)
        };
        let mut entropy_j_k = contribution(slopes.equilibrium_n_m_k, state.elastic_slip_m)?;
        for (matrix, &slip) in slopes.relaxation_n_m_k.iter().zip(&state.relaxation_slip_m) {
            entropy_j_k = checked_add(entropy_j_k, contribution(*matrix, slip)?,
                "layer_mechanical_entropy_j_k")?;
        }
        let free_energy_j = self.tensor_energy(state);
        nonnegative_finite(free_energy_j, "layer_free_energy_j")?;
        let internal_energy_j = checked_add(free_energy_j,
            checked_mul(self.temperature_kelvin, entropy_j_k, "layer_entropy_energy_j")?,
            "layer_internal_energy_j")?;
        Ok(ThirdBodyLayerThermodynamics { free_energy_j, entropy_j_k, internal_energy_j })
    }
    fn tensor_state<'a>(
        &'a self,
        state: &ThirdBodyLayerState,
    ) -> Result<&'a TensorLayerPotential, TriboError> {
        let potential = self
            .tensor
            .as_ref()
            .ok_or_else(|| invalid("layer_tensor_owner"))?;
        if state.owner != *self || state.relaxation_slip_m.len() != potential.arms.len() {
            return Err(invalid("layer_checkpoint_owner"));
        }
        self.admit_elastic_slip(state.elastic_slip_m)?;
        for value in state.plastic_slip_m {
            finite(value, "layer_plastic_slip_m")?;
        }
        for &slip in &state.relaxation_slip_m {
            self.admit_elastic_slip(slip)?;
        }
        Ok(potential)
    }
    pub(super) fn tensor_energy(&self, state: &ThirdBodyLayerState) -> f64 {
        let Ok(potential) = self.tensor_state(state) else {
            return f64::NAN;
        };
        let mut energy = 0.5 * potential.equilibrium.quadratic(state.elastic_slip_m);
        for (arm, &slip) in potential.arms.iter().zip(&state.relaxation_slip_m) {
            energy += 0.5 * arm.stiffness.quadratic(slip);
        }
        energy
    }
    pub(super) fn continue_tensor_state(
        &self,
        previous: &ThirdBodyLayerState,
        coefficient_work: bool,
    ) -> Result<ThirdBodyLayerContinuation, TriboError> {
        let old = previous.owner.tensor_state(previous)?;
        let new = self
            .tensor
            .as_ref()
            .ok_or_else(|| invalid("layer_tensor_owner"))?;
        if old.arms.len() != new.arms.len()
            || self.maximum_elastic_slip_m.to_bits()
                != previous.owner.maximum_elastic_slip_m.to_bits()
        {
            return Err(invalid(
                "layer_continuation_requires_unchanged_elastic_domain",
            ));
        }
        if !coefficient_work
            && (old.equilibrium != new.equilibrium
                || old
                    .arms
                    .iter()
                    .zip(&new.arms)
                    .any(|(a, b)| a.stiffness != b.stiffness))
        {
            return Err(invalid(
                "layer_continuation_requires_unchanged_elastic_potential",
            ));
        }
        let mut state = previous.clone();
        state.owner = self.clone();
        self.tensor_state(&state)?;
        let change = self.tensor_energy(&state) - previous.owner.tensor_energy(previous);
        finite(change, "layer_coefficient_work_j")?;
        Ok(ThirdBodyLayerContinuation {
            state,
            coefficient_free_energy_change_j: change,
        })
    }
    /// Rotate each original tensor AND its conjugate retained strains together.
    /// This is a coordinate transfer, not relative slip or a thermal step.
    pub fn transport_tensor_state(
        &self,
        state: &ThirdBodyLayerState,
        rotation: TangentialHistoryRotation,
        tolerance: TangentialTransportConfig,
    ) -> Result<TensorLayerFrameTransfer, TriboError> {
        tolerance
            .validate()
            .map_err(|_| invalid("layer_tensor_transport_tolerance"))?;
        let potential = self.tensor_state(state)?;
        if rotation.cosine_sine() == [1.0, 0.0] {
            return Ok(TensorLayerFrameTransfer {
                law: self.clone(),
                state: state.clone(),
                stored_energy_residual_j: 0.0,
            });
        }
        let mut law = Self::new_tensor(
            potential.equilibrium.rotated(rotation)?,
            potential
                .arms
                .iter()
                .map(|arm| {
                    Ok(TensorLayerRelaxationArm {
                        stiffness: arm.stiffness.rotated(rotation)?,
                        relaxation_time_s: arm.relaxation_time_s,
                    })
                })
                .collect::<Result<Vec<_>, TriboError>>()?,
            self.yield_force_n,
            self.plastic_resistance_n_s_m,
            self.temperature_kelvin,
            self.maximum_elastic_slip_m,
        )?;
        if let Some(viscosity) = &potential.viscosity {
            law = law.with_kelvin_viscosity(viscosity.rotated(rotation)?)?;
        }
        let transformed = ThirdBodyLayerState {
            elastic_slip_m: rotation.apply(state.elastic_slip_m)?,
            plastic_slip_m: rotation.apply(state.plastic_slip_m)?,
            relaxation_slip_m: state
                .relaxation_slip_m
                .iter()
                .map(|&r| rotation.apply(r))
                .collect::<Result<Vec<_>, _>>()?,
            owner: law.clone(),
        };
        let before = self.tensor_energy(state);
        let after = law.tensor_energy(&transformed);
        let residual = after - before;
        if !residual.is_finite()
            || residual.abs()
                > tolerance.absolute_energy_tolerance_j
                    + tolerance.relative_energy_tolerance * before.max(after)
        {
            return Err(invalid("layer_tensor_transport_changed_storage"));
        }
        Ok(TensorLayerFrameTransfer {
            law,
            state: transformed,
            stored_energy_residual_j: residual,
        })
    }
    pub(super) fn tensor_dissipation(
        &self,
        old: &ThirdBodyLayerState,
        new: &ThirdBodyLayerState,
        dt_s: f64,
    ) -> Result<ThirdBodyLayerDissipation, TriboError> {
        positive_finite(dt_s, "layer_dissipation_dt_s")?;
        let potential = self.tensor_state(old)?;
        self.tensor_state(new)?;
        let distance = det::hypot(
            new.plastic_slip_m[0] - old.plastic_slip_m[0],
            new.plastic_slip_m[1] - old.plastic_slip_m[1],
        );
        let plastic_j = if distance == 0.0 {
            0.0
        } else {
            checked_mul(
                self.yield_force_n + self.plastic_resistance_n_s_m * (distance / dt_s),
                distance,
                "layer_plastic_dissipation_j",
            )?
        };
        let mut relaxation_j = 0.0;
        for ((arm, before), after) in potential
            .arms
            .iter()
            .zip(&old.relaxation_slip_m)
            .zip(&new.relaxation_slip_m)
        {
            let mid = [
                0.5 * before[0] + 0.5 * after[0],
                0.5 * before[1] + 0.5 * after[1],
            ];
            relaxation_j = checked_add(
                relaxation_j,
                checked_mul(
                    dt_s / arm.relaxation_time_s,
                    arm.stiffness.quadratic(mid),
                    "layer_relaxation_heat_j",
                )?,
                "layer_relaxation_heat_j",
            )?;
        }
        nonnegative_finite(plastic_j, "layer_plastic_dissipation_j")?;
        nonnegative_finite(relaxation_j, "layer_relaxation_heat_j")?;
        let delta = [
            new.elastic_slip_m[0] - old.elastic_slip_m[0],
            new.elastic_slip_m[1] - old.elastic_slip_m[1],
        ];
        let kelvin_j = potential
            .viscosity
            .as_ref()
            .map_or(0.0, |c| c.factor.quadratic(delta) / dt_s);
        nonnegative_finite(kelvin_j, "layer_kelvin_heat_j")?;
        Ok(ThirdBodyLayerDissipation {
            plastic_j,
            relaxation_j,
            kelvin_j,
        })
    }

    pub(super) fn advance_tensor(
        &self,
        state: &ThirdBodyLayerState,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<ThirdBodyLayerState>, TriboError> {
        input.solve.validate()?;
        let potential = self.tensor_state(state)?;
        positive_finite(
            input.normal_n,
            "layer_requires_retained_compressive_contact",
        )?;
        positive_finite(input.dt_s, "dt_s")?;
        nonnegative_finite(input.compliance_m_per_n, "compliance_m_per_n")?;
        if input.contact_temperature_k.map(f64::to_bits) != Some(self.temperature_kelvin.to_bits())
        {
            return Err(invalid("layer_temperature_kelvin"));
        }
        if let Some(p) = input.nominal_pressure_pa {
            nonnegative_finite(p, "contact_pressure_pa")?;
        }
        for value in input.slip_increment_m {
            finite(value, "slip_increment_m")?;
        }
        let mut effective = potential.equilibrium.matrix;
        let mut history = apply(effective, state.elastic_slip_m);
        let mut weights = Vec::with_capacity(potential.arms.len());
        for (arm, &slip) in potential.arms.iter().zip(&state.relaxation_slip_m) {
            let ratio = input.dt_s / arm.relaxation_time_s;
            if !ratio.is_finite() || ratio > 2.0 {
                return Err(TriboError::OutsideApplicability {
                    field: "layer_dt_over_relaxation_time",
                    value: ratio,
                    minimum: 0.0,
                    maximum: 2.0,
                });
            }
            let weight = 1.0 / (1.0 + 0.5 * ratio);
            effective = add(effective, arm.stiffness.matrix, weight)?;
            let force = apply(arm.stiffness.matrix, slip);
            for i in 0..2 {
                history[i] = checked_add(history[i], weight * force[i], "layer_history_force_n")?;
            }
            weights.push(weight);
        }
        // F = history + (K_eff/2 + C/dt) Delta_e, with
        // Delta_e = free_slip - (body_compliance + plastic_compliance) F.
        // Keeping 2*C/dt inside the existing half-stiffness return preserves
        // all noncommuting tensor cross terms and the same scalar solve.
        if let Some(c) = &potential.viscosity {
            if c.factor.matrix.iter().flatten().any(|&v| v != 0.0) {
                effective = add(effective, c.factor.matrix, 2.0 / input.dt_s)?;
            }
        }
        let free = apply(effective, input.slip_increment_m);
        let rhs = [history[0] + 0.5 * free[0], history[1] + 0.5 * free[1]];
        for value in rhs {
            finite(value, "layer_trial_force_n")?;
        }
        let flow = input.dt_s / self.plastic_resistance_n_s_m;
        positive_finite(flow, "layer_plastic_step_compliance")?;
        let force_at = |plastic_compliance: f64| -> Result<[f64; 2], TriboError> {
            let load = 0.5 * (input.compliance_m_per_n + plastic_compliance);
            let matrix = add([[1.0, 0.0], [0.0, 1.0]], effective, load)?;
            solve_spd(matrix, rhs)
        };
        let elastic_trial = force_at(0.0)?;
        let yielding = det::hypot(elastic_trial[0], elastic_trial[1]) > self.yield_force_n;
        let mut force = elastic_trial;
        let mut iterations = 0;
        if yielding {
            let (mut low, mut high) = (0.0, flow);
            let mut residual = f64::INFINITY;
            let mut tolerance = 0.0;
            for iteration in 1..=input.solve.max_iterations {
                iterations = iteration;
                let compliance = 0.5 * low + 0.5 * high;
                force = force_at(compliance)?;
                let norm = det::hypot(force[0], force[1]);
                let physical = if norm > self.yield_force_n {
                    flow * (1.0 - self.yield_force_n / norm)
                } else {
                    0.0
                };
                let mismatch = apply(effective, force).map(|v| 0.5 * (compliance - physical) * v);
                residual = mismatch[0].abs().max(mismatch[1].abs());
                tolerance = input.solve.tolerance(norm, norm)?;
                if residual.is_finite() && residual <= tolerance {
                    break;
                }
                if compliance < physical {
                    low = compliance;
                } else {
                    high = compliance;
                }
            }
            if !residual.is_finite() || residual > tolerance {
                return Err(TriboError::TangentialSolve {
                    residual_n: residual,
                    tolerance_n: tolerance,
                    iterations,
                });
            }
        }
        let magnitude = det::hypot(force[0], force[1]);
        finite(magnitude, "layer_return_force_n")?;
        if magnitude > input.solve.maximum_force_n {
            return Err(TriboError::OutsideApplicability {
                field: "maximum_force_n",
                value: magnitude,
                minimum: 0.0,
                maximum: input.solve.maximum_force_n,
            });
        }
        let plastic_compliance = if yielding && magnitude > self.yield_force_n {
            flow * (1.0 - self.yield_force_n / magnitude)
        } else {
            0.0
        };
        let slip = [
            input.slip_increment_m[0] - input.compliance_m_per_n * force[0],
            input.slip_increment_m[1] - input.compliance_m_per_n * force[1],
        ];
        let mut candidate = state.clone();
        for i in 0..2 {
            let plastic = plastic_compliance * force[i];
            candidate.elastic_slip_m[i] += slip[i] - plastic;
            candidate.plastic_slip_m[i] += plastic;
        }
        let delta = [
            candidate.elastic_slip_m[0] - state.elastic_slip_m[0],
            candidate.elastic_slip_m[1] - state.elastic_slip_m[1],
        ];
        let elastic_mid = [
            state.elastic_slip_m[0] + 0.5 * delta[0],
            state.elastic_slip_m[1] + 0.5 * delta[1],
        ];
        let mut midpoint_force = apply(potential.equilibrium.matrix, elastic_mid);
        if let Some(c) = &potential.viscosity {
            let viscous_force = apply(c.factor.matrix, delta);
            for i in 0..2 {
                midpoint_force[i] += viscous_force[i] / input.dt_s;
            }
        }
        for (j, arm) in potential.arms.iter().enumerate() {
            let mid = [
                weights[j] * (state.relaxation_slip_m[j][0] + 0.5 * delta[0]),
                weights[j] * (state.relaxation_slip_m[j][1] + 0.5 * delta[1]),
            ];
            candidate.relaxation_slip_m[j] = [
                2.0 * mid[0] - state.relaxation_slip_m[j][0],
                2.0 * mid[1] - state.relaxation_slip_m[j][1],
            ];
            let branch_force = apply(arm.stiffness.matrix, mid);
            for i in 0..2 {
                midpoint_force[i] += branch_force[i];
            }
        }
        self.tensor_state(&candidate)?;
        let residual = (force[0] - midpoint_force[0])
            .abs()
            .max((force[1] - midpoint_force[1]).abs());
        let tolerance = input.solve.tolerance(magnitude, magnitude)?;
        if !residual.is_finite() || residual > tolerance {
            return Err(TriboError::TangentialSolve {
                residual_n: residual,
                tolerance_n: tolerance,
                iterations,
            });
        }
        let old = self.tensor_energy(state);
        let new = self.tensor_energy(&candidate);
        nonnegative_finite(old, "layer_old_energy_j")?;
        nonnegative_finite(new, "layer_new_energy_j")?;
        let stored_change_j = new - old;
        let loss = self.tensor_dissipation(state, &candidate, input.dt_s)?;
        let dissipated_j = checked_add(
            checked_add(loss.plastic_j, loss.relaxation_j, "layer_dissipated_j")?,
            loss.kelvin_j,
            "layer_dissipated_j",
        )?;
        let closure_residual_j = stored_change_j + dissipated_j - dot(force, slip);
        finite(closure_residual_j, "layer_closure_residual_j")?;
        Ok(TangentialStep {
            traction_n: force,
            slip_increment_m: slip,
            state: candidate,
            regime: if yielding {
                FrictionRegime::Sliding
            } else {
                FrictionRegime::Sticking
            },
            stored_change_j,
            dissipated_j,
            closure_residual_j,
            compatibility_residual_m: [0.0; 2],
            iterations,
        })
    }
}

fn invalid(field: &'static str) -> TriboError {
    TriboError::InvalidInput { field }
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn apply(matrix: [[f64; 2]; 2], vector: [f64; 2]) -> [f64; 2] {
    [dot(matrix[0], vector), dot(matrix[1], vector)]
}
fn add(a: [[f64; 2]; 2], b: [[f64; 2]; 2], scale: f64) -> Result<[[f64; 2]; 2], TriboError> {
    let mut out = a;
    for i in 0..2 {
        for j in 0..2 {
            out[i][j] = checked_add(
                a[i][j],
                checked_mul(scale, b[i][j], "layer_tensor_matrix")?,
                "layer_tensor_matrix",
            )?;
        }
    }
    Ok(out)
}
fn solve_spd(a: [[f64; 2]; 2], b: [f64; 2]) -> Result<[f64; 2], TriboError> {
    let l = det::sqrt(a[0][0]);
    positive_finite(l, "layer_tensor_load_line")?;
    let off = a[1][0] / l;
    let pivot = a[1][1] - off * off;
    positive_finite(pivot, "layer_tensor_load_line")?;
    let second = (b[1] - off * (b[0] / l)) / pivot;
    let first = (b[0] / l - off * second) / l;
    finite(first, "layer_tensor_force")?;
    finite(second, "layer_tensor_force")?;
    Ok([first, second])
}
