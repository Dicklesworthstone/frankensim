//! One held incident wave at an explicitly selected original physical port.
//! The reference impedance is a boundary relation, never material viscosity.
use super::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, finite,
    positive_heat::PositiveHeat, prepared::PreparedQuadraticFlow, refuse, split_quadratic,
};
use crate::{PhsError, PortHamiltonian, storage_energy_change};
use fs_math::det;

/// Independent physical accounts over the complete prepared duration.
#[derive(Clone, Copy, Debug)]
pub struct ScatteringStepRecord {
    /// Original constitutive storage change [J].
    pub delta_h: f64,
    /// Integral of the ORIGINAL physical resistance, excluding reference loading [J].
    pub physical_heat: f64,
    /// Work through all unselected original ports [J].
    pub external_work: f64,
    /// Net incoming boundary work, including outgoing-wave exchange [J].
    pub boundary_work: f64,
    /// Integral of a²/Z for the held incident wave [J].
    pub incident_energy: f64,
    /// Independently integrated (y-a)²/Z, with no energy subtraction [J].
    pub outgoing_energy: f64,
    /// Interval mean outgoing wave b = y-a, in original port-output units.
    pub outgoing_wave: f64,
}

/// Complete endpoint and ORIGINAL mean port outputs, in their original order.
pub struct ScatteringStep {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub account: ScatteringStepRecord,
}

struct Scratch {
    augmented: Vec<f64>,
    balanced: Vec<f64>,
    endpoint: Vec<f64>,
    integrals: Vec<f64>,
    output: Vec<f64>,
    midpoint: Vec<f64>,
    effort: Vec<f64>,
    derivative: Vec<f64>,
    heat_coordinates: Vec<f64>,
}

/// Exact complete quadratic flow with one selected scattering boundary.
///
/// Original input/output are (u,y), with power u*y. The selected input slot
/// instead holds an incident wave a in y's units: u=(2a-y)/Z, b=y-a.
/// Z has the original y/u units (acoustic pressure/volume-flow impedance).
/// Every other input retains its original meaning; every returned output is
/// the ORIGINAL y. All states, storage, material memory and physical R remain
/// represented. Preparation uses the existing exact flow and positive heat
/// owners; warm calls do no exponentials or constitutive model inference.
///
/// The owner stores numerical scratch, not accepted physical history. Its
/// complete candidate publishes only on success. Clock/source bandwidth and
/// nonlinear exciter closure remain the consumer's responsibility.
pub struct PreparedScatteringZoh {
    original: PortHamiltonian,
    robin: PortHamiltonian,
    flow: PreparedQuadraticFlow,
    physical_heat: PositiveHeat,
    reference: PositiveHeat,
    port: usize,
    dt: f64,
    inverse_impedance: f64,
    root_inverse_impedance: f64,
    active_port: bool,
    pressure_state: Vec<f64>,
    pressure_input: Vec<f64>,
    scratch: Scratch,
}

impl PreparedScatteringZoh {
    /// Own the original physical system and one fixed positive-clock cache.
    /// Only an explicitly declared exact quadratic storage is admitted.
    /// Existing complete state/port and dense preparation limits apply.
    pub fn new_owned(
        mut system: PortHamiltonian,
        port: usize,
        impedance: f64,
        duration: f64,
        budget: LinearSplitBudget,
    ) -> Result<Self, PhsError> {
        if port >= system.m {
            return Err(PhsError::BadPortPairing);
        }
        let inverse_impedance = 1.0 / impedance;
        if !impedance.is_finite()
            || impedance <= 0.0
            || !inverse_impedance.is_finite()
            || inverse_impedance <= 0.0
            || !duration.is_finite()
            || duration <= 0.0
        {
            return Err(refuse(
                "scattering requires a representable positive impedance and clock",
            ));
        }
        // Admit the original storage and complete capacity before copying.
        let original_split = split_quadratic(&system, budget)?;
        if !original_split.is_exact_quadratic() {
            return Err(refuse(
                "scattering requires declared identically quadratic storage",
            ));
        }
        let (n, m) = (system.n, system.m);
        let root_inverse_impedance = det::sqrt(inverse_impedance);
        let mut reference_root = vec![0.0; n * n];
        for i in 0..n {
            reference_root[i] = root_inverse_impedance * system.g[i * m + port];
        }
        let active_port = reference_root.iter().any(|value| *value != 0.0);
        let mut reference_r = vec![0.0; n * n];
        let mut robin = system.share_storage();
        for i in 0..n {
            for j in 0..n {
                let value = reference_root[i] * reference_root[j];
                reference_r[i * n + j] = value;
                robin.r[i * n + j] += value;
            }
            robin.g[i * m + port] = 2.0 * inverse_impedance * system.g[i * m + port];
            if system.g[i * m + port] != 0.0
                && (reference_r[i * n + i] == 0.0 || robin.g[i * m + port] == 0.0)
            {
                return Err(refuse("scattering port scale is unrepresentable"));
            }
        }
        if !finite(&reference_r)
            || !finite(&reference_root)
            || !finite(&robin.r)
            || !finite(&robin.g)
        {
            return Err(refuse("scattering boundary structure is nonfinite"));
        }
        // The modified R is original material plus an explicit rank-one passive
        // boundary. The reference contribution is accounted separately below.
        let split = split_quadratic(&robin, budget)?;
        let q = split.quadratic_matrix();
        let outputs: Vec<Vec<f64>> = (0..m)
            .map(|p| {
                (0..n)
                    .map(|j| (0..n).map(|k| system.g[k * m + p] * q[k * n + j]).sum())
                    .collect()
            })
            .collect();
        let flow = PreparedQuadraticFlow::new(
            &robin,
            &split,
            duration,
            &[],
            &outputs,
            DenseExponentialBudget {
                max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                max_dense_bytes: usize::MAX,
                max_multiply_terms: u64::MAX,
            },
        )?;
        let (mut pressure_state, mut pressure_input) = flow.linear_integral_map(port)?;
        for value in pressure_state.iter_mut().chain(&mut pressure_input) {
            *value /= duration;
        }
        if !finite(&pressure_state) || !finite(&pressure_input) {
            return Err(refuse("scattering mean output map is nonfinite"));
        }
        let physical_heat = PositiveHeat::new(&robin, &split, duration, &system.r)?;
        let reference =
            PositiveHeat::from_factor(&robin, &split, duration, &reference_r, reference_root)?;
        Ok(Self {
            original: system,
            robin,
            flow,
            physical_heat,
            reference,
            port,
            dt: duration,
            inverse_impedance,
            root_inverse_impedance,
            active_port,
            pressure_state,
            pressure_input,
            scratch: Scratch {
                augmented: vec![0.0; n + m + 1],
                balanced: vec![0.0; n],
                endpoint: vec![0.0; n],
                integrals: vec![0.0; m],
                output: vec![0.0; m],
                midpoint: vec![0.0; n],
                effort: vec![0.0; n],
                derivative: vec![0.0; n],
                heat_coordinates: vec![0.0; 2 * n],
            },
        })
    }

    /// Original physical system; reference impedance is not part of its R.
    pub fn system(&self) -> &PortHamiltonian {
        &self.original
    }
    pub fn duration(&self) -> f64 {
        self.dt
    }
    pub fn selected_port(&self) -> usize {
        self.port
    }

    /// Current-window relation y_bar = free + beta*a, in original output units.
    /// The selected input value is ignored; all other authored inputs enter
    /// free. This is one O(n+m) contraction of the already prepared exact map.
    pub fn average_relation(&self, state: &[f64], inputs: &[f64]) -> Result<(f64, f64), PhsError> {
        self.admit_layout(state, inputs)?;
        let free = self
            .pressure_state
            .iter()
            .zip(state)
            .map(|(a, x)| a * x)
            .sum::<f64>()
            + self
                .pressure_input
                .iter()
                .zip(inputs)
                .enumerate()
                .filter(|(p, _)| *p != self.port)
                .map(|(_, (a, u))| a * u)
                .sum::<f64>();
        if !free.is_finite() {
            return Err(refuse("scattering free mean output is nonfinite"));
        }
        Ok((free, self.pressure_input[self.port]))
    }

    /// Allocating complete candidate; accepted source state is never modified.
    pub fn step(&mut self, state: &[f64], inputs: &[f64]) -> Result<ScatteringStep, PhsError> {
        let account = self.evaluate(state, inputs)?;
        Ok(ScatteringStep {
            x: self.scratch.endpoint.clone(),
            y: self.scratch.output.clone(),
            account,
        })
    }

    /// Same candidate and physical accounts into complete caller buffers.
    /// Successful warm calls allocate no owner memory; storage callbacks must
    /// obey their existing allocation-free contract. Refusal changes neither output.
    pub fn step_into(
        &mut self,
        state: &[f64],
        inputs: &[f64],
        next: &mut [f64],
        average_outputs: &mut [f64],
    ) -> Result<ScatteringStepRecord, PhsError> {
        if next.len() != self.original.n || average_outputs.len() != self.original.m {
            return Err(refuse(
                "scattering output buffers require complete original layout",
            ));
        }
        let account = self.evaluate(state, inputs)?;
        next.copy_from_slice(&self.scratch.endpoint);
        average_outputs.copy_from_slice(&self.scratch.output);
        Ok(account)
    }

    fn admit_layout(&self, state: &[f64], inputs: &[f64]) -> Result<(), PhsError> {
        if state.len() != self.original.n
            || inputs.len() != self.original.m
            || !finite(state)
            || !finite(inputs)
        {
            return Err(refuse(
                "scattering requires complete finite original state and input",
            ));
        }
        Ok(())
    }

    fn evaluate(
        &mut self,
        state: &[f64],
        inputs: &[f64],
    ) -> Result<ScatteringStepRecord, PhsError> {
        self.admit_layout(state, inputs)?;
        let s = &mut self.scratch;
        self.flow.evaluate_into(
            state,
            inputs,
            &mut s.augmented,
            &mut s.balanced,
            &mut s.endpoint,
            &mut s.integrals,
        )?;
        for (out, integral) in s.output.iter_mut().zip(&s.integrals) {
            *out = integral / self.dt;
        }
        let physical_heat = self.physical_heat.evaluate_into(
            &self.robin,
            state,
            inputs,
            &mut s.effort,
            &mut s.derivative,
            &mut s.heat_coordinates,
        )?;
        let reference = self.reference.evaluate_into(
            &self.robin,
            state,
            inputs,
            &mut s.effort,
            &mut s.derivative,
            &mut s.heat_coordinates,
        )?;
        let a = inputs[self.port];
        let shifted = a * self.root_inverse_impedance;
        let incident_amplitude = shifted * det::sqrt(self.dt);
        let incident_energy = incident_amplitude * incident_amplitude;
        let outgoing_energy = if self.active_port {
            self.reference.evaluate_shifted_into(
                &self.robin,
                state,
                inputs,
                shifted,
                &mut s.effort,
                &mut s.derivative,
                &mut s.heat_coordinates,
            )?
        } else {
            incident_energy
        };
        let selected_work = (2.0 * self.inverse_impedance * a) * s.integrals[self.port];
        let boundary_work = selected_work - reference;
        let external_work: f64 = inputs
            .iter()
            .zip(&s.integrals)
            .enumerate()
            .filter(|(p, _)| *p != self.port)
            .map(|(_, (u, integral))| u * integral)
            .sum();
        let delta_h = storage_energy_change(
            self.original.storage.as_ref(),
            state,
            &s.endpoint,
            &mut s.midpoint,
            &mut s.effort,
        );
        let scale = super::prepared_zoh::energetic_scale_into(&self.original, state, &mut s.effort)
            .max(super::prepared_zoh::energetic_scale_into(
                &self.original,
                &s.endpoint,
                &mut s.effort,
            ))
            .max(delta_h.abs())
            .max(physical_heat.abs())
            .max(external_work.abs())
            .max(selected_work.abs())
            .max(reference.abs());
        let wave_scale = scale.max(incident_energy.abs()).max(outgoing_energy.abs());
        let balance = delta_h + physical_heat - external_work - boundary_work;
        let wave_balance = boundary_work - (incident_energy - outgoing_energy);
        let outgoing_wave = s.output[self.port] - a;
        if !finite(&s.output)
            || !self.original.hamiltonian(&s.endpoint).is_finite()
            || !finite(&[
                physical_heat,
                reference,
                incident_energy,
                outgoing_energy,
                selected_work,
                external_work,
                delta_h,
                scale,
                wave_scale,
                balance,
                wave_balance,
                outgoing_wave,
            ])
            || physical_heat < 0.0
            || outgoing_energy < 0.0
            || balance.abs() > 1e-12 * scale
            || wave_balance.abs() > 1e-12 * wave_scale
        {
            return Err(refuse(
                "scattering original physical work and wave balance refused",
            ));
        }
        Ok(ScatteringStepRecord {
            delta_h,
            physical_heat,
            external_work,
            boundary_work,
            incident_energy,
            outgoing_energy,
            outgoing_wave,
        })
    }
}
