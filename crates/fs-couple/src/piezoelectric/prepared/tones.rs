//! Continuous sinusoidal force/current drives on the original piezoelectric rod.
//!
//! One shared cosine/sine carrier per frequency drives both original physical
//! ports. Carriers are numerical source coordinates with no stored energy.
//! The public native extended quadratic-flow owner propagates the complete
//! original rod and integrates F(t)*v(t), I(t)*V(t), I(t) and V(t). Signed work
//! observations are never relabeled as heat or inferred from an energy deficit.
//!
//! Frequencies and amplitudes are fixed after preparation. DC is an explicit
//! zero-frequency cosine term. Accepted carrier phase and the original body
//! persist together; no phase normalization, energy rescaling or mechanical
//! mode replacement is performed. Exact refers to this linear source model
//! up to floating-point error, not continuum accuracy or a real-time deadline.

use super::{
    ExactPiezoelectricControls, PiezoelectricBody, PiezoelectricCouplingError,
    Result, invalid, poll, positive,
};
use fs_math::det;
use fs_phs::PreparedStepError;
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
    split_quadratic,
};
use fs_phs::linear_exact::prepared::PreparedQuadraticFlow;

/// One fixed carrier, with independent force and current phase offsets:
/// F(t)=F_cos*cos(omega*t)+F_sin*sin(omega*t), and likewise for I(t).
/// Time zero is this prepared owner's initial state. Frequencies must be
/// distinct; combine coefficients for a shared frequency before preparation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PiezoelectricTone {
    pub angular_frequency_rad_s: f64,
    pub force_cosine_n: f64,
    pub force_sine_n: f64,
    pub current_cosine_a: f64,
    pub current_sine_a: f64,
}

/// Physical source integrals and acceptance residuals for one continuous window.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExactPiezoelectricToneStep {
    pub duration_s: f64,
    /// Rounded views of accepted interval count times the fixed duration [s].
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub energy_change_j: f64,
    /// Native signed integral of force times physical velocity [J].
    pub mechanical_supply_j: f64,
    /// Native signed integral of current times terminal voltage [J].
    pub electrical_supply_j: f64,
    /// Native linear integral of the original total current waveform [C].
    pub integrated_current_c: f64,
    pub energy_residual_j: f64,
    pub charge_residual_c: f64,
    pub mean_voltage_v: f64,
    pub mean_velocity_m_s: f64,
    /// Conservative computed full-window physical energy bound. It can refuse
    /// a valid highly oscillatory window and is not an outward-rounded certificate.
    pub interval_energy_bound_j: f64,
}

/// Fixed-clock continuous tones acting on the complete original rod.
/// Warm stepping allocates no owner buffers and performs no exponential,
/// factorization, Newton solve or material clone. Original body, carrier phase,
/// accepted clock and last report publish atomically after all physical checks.
pub struct PreparedPiezoelectricTones {
    body: PiezoelectricBody,
    tones: Vec<PiezoelectricTone>,
    amplitudes: Vec<[f64; 2]>,
    flow: PreparedQuadraticFlow,
    duration_s: f64,
    root_mass: f64,
    root_capacitance: f64,
    root_stiffness: f64,
    accepted_intervals: u64,
    elapsed_time_s: f64,
    // Original [extension, momentum, charge], then [cosine, sine] per carrier.
    // The original body alone owns all physical storage.
    state: Vec<f64>,
    endpoint: Vec<f64>,
    augmented: Vec<f64>,
    balanced: Vec<f64>,
    integrals: [f64; 4],
    last_step: Option<ExactPiezoelectricToneStep>,
}

impl PreparedPiezoelectricTones {
    /// Cache the complete original flow and four independent observations.
    /// Both physical and numerical source coordinates count against max_modes;
    /// the native dense budget also includes every source and work observation.
    /// Empty tones give an isolated, unforced rod. DC sine coefficients must
    /// be zero because the zero-frequency sine does not drive either port.
    pub fn new(
        body: PiezoelectricBody,
        tones: &[PiezoelectricTone],
        duration_s: f64,
        split_budget: LinearSplitBudget,
        dense_budget: DenseExponentialBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        if !positive(duration_s) {
            return Err(invalid("piezoelectric tones require a positive finite duration"));
        }
        let n = tones.len().checked_mul(2).and_then(|extra| extra.checked_add(3))
            .ok_or_else(|| invalid("piezoelectric tone state extent overflows"))?;
        let augmented_dimension = n.checked_add(3)
            .ok_or_else(|| invalid("piezoelectric tone augmented extent overflows"))?;
        if n > split_budget.max_modes
            || augmented_dimension > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || augmented_dimension > dense_budget.max_dimension
            || n.checked_mul(n).and_then(|square| square.checked_mul(3*8))
                .is_none_or(|bytes| bytes > dense_budget.max_dense_bytes)
        {
            return Err(invalid("complete piezoelectric rod/source layout exceeds the explicit preparation budget"));
        }
        body.stored_energy_j()?;
        let root_mass = det::sqrt(body.rod.modal_mass_kg());
        let root_capacitance = det::sqrt(body.rod.axially_clamped_capacitance_f());
        let root_stiffness = det::sqrt(body.rod.short_circuit_stiffness_n_m());
        let mut amplitudes = Vec::with_capacity(tones.len());
        let mut extra_scales = Vec::with_capacity(n-3);
        for (index, tone) in tones.iter().enumerate() {
            poll(&mut cancelled)?;
            let frequency = tone.angular_frequency_rad_s;
            if [frequency, tone.force_cosine_n, tone.force_sine_n,
                tone.current_cosine_a, tone.current_sine_a].iter().any(|v| !v.is_finite())
                || frequency < 0.0
                || tones[..index].iter().any(|other|
                    other.angular_frequency_rad_s == frequency)
                || frequency == 0.0 && (tone.force_sine_n != 0.0 || tone.current_sine_a != 0.0)
                || !(frequency*duration_s).is_finite()
                || frequency != 0.0 && frequency*duration_s == 0.0
            {
                return Err(invalid("piezoelectric tones need distinct finite nonnegative frequencies and representable finite coefficients/phase"));
            }
            let force = tone.force_cosine_n.hypot(tone.force_sine_n);
            let current = tone.current_cosine_a.hypot(tone.current_sine_a);
            let force_rate = force/root_mass;
            let current_rate = current/root_capacitance;
            let rate = force_rate.hypot(current_rate);
            let scale = if force == 0.0 && current == 0.0 {
                1.0
            } else {
                duration_s*rate
            };
            if !force.is_finite() || !current.is_finite() || !positive(scale)
                || force != 0.0 && force_rate == 0.0
                || current != 0.0 && current_rate == 0.0
            {
                return Err(invalid("piezoelectric source energy-coordinate scale is unrepresentable"));
            }
            amplitudes.push([force,current]);
            extra_scales.extend_from_slice(&[scale,scale]);
        }
        let system = body.rod.port_hamiltonian()?;
        let (j,r,g) = system.structure();
        if system.state_dim() != 3 || system.port_dim() != 2
            || r.iter().any(|&value| value != 0.0)
        {
            return Err(invalid("continuous tones require the original lossless three-state two-port rod"));
        }
        let split = split_quadratic(&system,split_budget)?;
        let q = split.quadratic_matrix();
        let mut generator = vec![0.0; n*n];
        for i in 0..3 {
            for k in 0..3 {
                for column in 0..3 {
                    generator[i*n+column] += (j[i*3+k]-r[i*3+k])*q[k*3+column];
                }
            }
        }
        // Physical output rows G^T Q: velocity and terminal voltage.
        let output: [[f64; 3]; 2] = core::array::from_fn(|port|
            core::array::from_fn(|column| (0..3)
                .map(|i| g[i*2+port]*q[i*3+column]).sum()));
        let mut work = [vec![0.0; n*n],vec![0.0; n*n]];
        let mut current = vec![0.0; n];
        let mut voltage = vec![0.0; n];
        voltage[..3].copy_from_slice(&output[1]);
        let mut state = vec![0.0; n];
        state[..3].copy_from_slice(&body.state);
        for (index, tone) in tones.iter().enumerate() {
            poll(&mut cancelled)?;
            let carrier = 3+2*index;
            state[carrier] = 1.0;
            generator[carrier*n+carrier+1] = -tone.angular_frequency_rad_s;
            generator[(carrier+1)*n+carrier] = tone.angular_frequency_rad_s;
            for (offset, coefficients) in [
                [tone.force_cosine_n,tone.current_cosine_a],
                [tone.force_sine_n,tone.current_sine_a],
            ].into_iter().enumerate() {
                let source = carrier+offset;
                current[source] = coefficients[1];
                for i in 0..3 {
                    generator[i*n+source] = g[i*2]*coefficients[0]+g[i*2+1]*coefficients[1];
                    for port in 0..2 {
                        // z^T W z = source_waveform * original physical output.
                        let (a,b) = (coefficients[port],output[port][i]);
                        // Halve the larger factor before multiplying, retaining
                        // representable work across mixed force/current scales.
                        let value = if a.abs() >= b.abs() { (0.5*a)*b } else { a*(0.5*b) };
                        work[port][i*n+source] = value;
                        work[port][source*n+i] = value;
                    }
                }
            }
        }
        let flow = PreparedQuadraticFlow::new_extended_controlled(
            &system,&split,&generator,&extra_scales,duration_s,
            &work,&[current,voltage],dense_budget,&mut cancelled,
        ).map_err(|error| match error {
            PreparedStepError::Cancelled => PiezoelectricCouplingError::Cancelled,
            PreparedStepError::Solver(error) => PiezoelectricCouplingError::Native(error),
        })?;
        poll(&mut cancelled)?;
        Ok(Self {
            body, tones: tones.to_vec(), amplitudes, flow, duration_s,
            root_mass, root_capacitance, root_stiffness,
            accepted_intervals: 0, elapsed_time_s: 0.0,
            state, endpoint: vec![0.0; n], augmented: vec![0.0; augmented_dimension],
            balanced: vec![0.0; n], integrals: [0.0; 4], last_step: None,
        })
    }

    pub const fn body(&self) -> &PiezoelectricBody { &self.body }
    pub const fn duration_s(&self) -> f64 { self.duration_s }
    pub const fn accepted_intervals(&self) -> u64 { self.accepted_intervals }
    pub const fn elapsed_time_s(&self) -> f64 { self.elapsed_time_s }
    /// Original coefficients at this owner's time zero, in caller frequency order.
    pub fn tones(&self) -> &[PiezoelectricTone] { &self.tones }
    /// Accepted [cosine, sine] carrier coordinates in the same order as tones().
    /// Checkpoints must retain these together with the body and accepted clock.
    /// Small numerical phase-norm drift is retained and included in domain bounds.
    pub fn phase_quadratures(&self) -> &[f64] { &self.state[3..] }
    pub const fn last_step(&self) -> Option<&ExactPiezoelectricToneStep> {
        self.last_step.as_ref()
    }
    pub fn into_body(self) -> PiezoelectricBody { self.body }

    /// Advance the continuous program for one fixed interval. Neither source
    /// power is an endpoint product, and integrated current is independently
    /// compared with the original electrode charge increment.
    pub fn step(
        &mut self,
        controls: ExactPiezoelectricControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExactPiezoelectricToneStep> {
        poll(&mut cancelled)?;
        if !positive(controls.energy_tolerance_j) || !positive(controls.charge_tolerance_c) {
            return Err(invalid("continuous piezoelectric sources require positive finite energy/charge tolerances"));
        }
        let next_interval = self.accepted_intervals.checked_add(1)
            .ok_or_else(|| invalid("piezoelectric source clock is exhausted"))?;
        let next_time = next_interval as f64*self.duration_s;
        if !next_time.is_finite() || next_time <= self.elapsed_time_s {
            return Err(invalid("piezoelectric source clock increment is unrepresentable"));
        }
        let before = self.body.stored_energy_j()?;
        let mut force_bound = 0.0;
        let mut current_bound = 0.0;
        for (index, amplitude) in self.amplitudes.iter().enumerate() {
            poll(&mut cancelled)?;
            let carrier = 3+2*index;
            // Evaluate the ACTUAL accepted carrier, including accumulated
            // floating-point norm drift, rather than resetting its amplitude.
            let radius = self.state[carrier].hypot(self.state[carrier+1]);
            force_bound += amplitude[0]*radius;
            current_bound += amplitude[1]*radius;
        }
        let rod = &self.body.rod;
        let force_rate = force_bound/self.root_mass;
        let current_rate = current_bound/self.root_capacitance;
        let source_increment = self.duration_s*force_rate.hypot(current_rate);
        // H = k_E*x^2/2 + m*v^2/2 + C*V^2/2. Original power conjugacy gives
        // d sqrt(2H)/dt <= hypot(|F|/sqrt(m),|I|/sqrt(C)).
        // The sum of each carrier's actual amplitude bounds both waveforms
        // at every interior time, including resonance and DC.
        let norm = (det::sqrt(2.0*before)+source_increment)
            * (1.0+128.0*f64::EPSILON*(self.tones.len()+1) as f64);
        let interval_energy_bound_j = (0.5*norm)*norm;
        let extension_bound = norm/self.root_stiffness;
        let voltage_bound = norm/self.root_capacitance;
        let limits = rod.limits();
        if [force_bound,current_bound,source_increment,norm,interval_energy_bound_j,
            extension_bound,voltage_bound].iter().any(|v| !v.is_finite())
            || force_bound != 0.0 && force_rate == 0.0
            || current_bound != 0.0 && current_rate == 0.0
            || (force_bound != 0.0 || current_bound != 0.0) && source_increment == 0.0
            || extension_bound/rod.length_m() > limits.max_abs_strain
            || voltage_bound/rod.length_m() > limits.max_abs_electric_field_v_m
        {
            return Err(invalid("complete continuous-tone window exceeds the original piezoelectric source envelope"));
        }
        self.flow.evaluate_into(
            &self.state,&[0.0,0.0],&mut self.augmented,&mut self.balanced,
            &mut self.endpoint,&mut self.integrals,
        )?;
        poll(&mut cancelled)?;
        let after = rod.energy(&self.endpoint[..3])?;
        let energy_change_j = after-before;
        let mechanical_supply_j = self.integrals[0];
        let electrical_supply_j = self.integrals[1];
        let integrated_current_c = self.integrals[2];
        let energy_residual_j = energy_change_j-mechanical_supply_j-electrical_supply_j;
        let charge_residual_c = self.endpoint[2]-self.state[2]-integrated_current_c;
        let mean_voltage_v = self.integrals[3]/self.duration_s;
        let mean_velocity_m_s = (self.endpoint[0]-self.state[0])/self.duration_s;
        if self.endpoint.iter().any(|v| !v.is_finite())
            || [energy_change_j,mechanical_supply_j,electrical_supply_j,integrated_current_c,
                energy_residual_j,charge_residual_c,mean_voltage_v,mean_velocity_m_s]
                .iter().any(|v| !v.is_finite())
        {
            return Err(invalid("continuous piezoelectric phase, charge or source work is unrepresentable"));
        }
        if charge_residual_c.abs() > controls.charge_tolerance_c {
            return Err(invalid("continuous piezoelectric source failed integrated free-charge conservation"));
        }
        if energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(PiezoelectricCouplingError::Balance {
                residual_j: energy_residual_j, tolerance_j: controls.energy_tolerance_j,
            });
        }
        let report = ExactPiezoelectricToneStep {
            duration_s: self.duration_s, start_time_s: self.elapsed_time_s,
            end_time_s: next_time, energy_change_j, mechanical_supply_j, electrical_supply_j,
            integrated_current_c, energy_residual_j, charge_residual_c, mean_voltage_v,
            mean_velocity_m_s, interval_energy_bound_j,
        };
        poll(&mut cancelled)?;
        self.body.state.copy_from_slice(&self.endpoint[..3]);
        self.state.copy_from_slice(&self.endpoint);
        self.accepted_intervals = next_interval;
        self.elapsed_time_s = next_time;
        self.last_step = Some(report);
        Ok(report)
    }
}
