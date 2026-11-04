//! The original complete coupled pencil under continuous tone forcing.
use super::*;
use fs_phs::linear_exact::{DenseExponentialBudget, PreparedExactTones};

// Same independent-channel reconstruction allowance as the prepared owner.
const CHANNEL_SUM_RELATIVE_TOLERANCE: f64 = 1e-11;

/// Private because the preparation and its exact source must remain together.
pub(crate) struct PreparedJointTones {
    flow: PreparedExactTones<'static>,
    frequencies: Vec<f64>,
    groups: Vec<Vec<usize>>,
    original_to_group: Vec<usize>,
    rest_tone: Option<usize>,
    rest_force: Vec<f64>,
    link_channels: Vec<Option<usize>>,
    initial: Vec<f64>,
    endpoint: Vec<f64>,
    cosine: Vec<f64>,
    sine: Vec<f64>,
    cosine_end: Vec<f64>,
    sine_end: Vec<f64>,
    mean_velocity: Vec<f64>,
    channel_heat: Vec<f64>,
    endpoint_external: Vec<f64>,
    force_envelope: Vec<f64>,
}

impl PreparedJointTones {
    pub(crate) fn new(
        network: &CoupledModalSystem,
        frequencies: &[f64],
        budget: DenseExponentialBudget,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        Self::prepare(network, frequencies, budget, gate, false)
    }

    pub(crate) fn new_affine(
        network: &CoupledModalSystem,
        frequencies: &[f64],
        budget: DenseExponentialBudget,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        Self::prepare(network, frequencies, budget, gate, true)
    }

    fn prepare(
        network: &CoupledModalSystem,
        frequencies: &[f64],
        budget: DenseExponentialBudget,
        gate: &CancelGate,
        affine: bool,
    ) -> Result<Self, ModalCouplingError> {
        poll(Some(gate))?;
        network.require_time_domain()?;
        let count = network.mode_count();
        let maximum_frequency =
            core::f64::consts::PI * network.config.nyquist_guard_fraction / network.dt;
        if frequencies.is_empty()
            || frequencies.len() > fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || frequencies
                .iter()
                .any(|w| !w.is_finite() || *w < 0.0 || *w > maximum_frequency)
            || network
                .models
                .iter()
                .flat_map(|model| model.modes())
                .any(|mode| mode.pressure_per_modal_velocity.im != 0.0)
        {
            return Err(invalid(
                "coupled tones require bounded resolved nonnegative frequencies and real instantaneous pressure maps",
            ));
        }
        let mut unique = Vec::new();
        let mut groups: Vec<Vec<usize>> = Vec::new();
        let mut original_to_group = Vec::with_capacity(frequencies.len());
        for (index, &frequency) in frequencies.iter().enumerate() {
            let group = if let Some(group) = unique.iter().position(|&w| w == frequency) {
                group
            } else {
                unique.push(frequency);
                groups.push(Vec::new());
                unique.len() - 1
            };
            groups[group].push(index);
            original_to_group.push(group);
        }
        // A spring rest length is a physical constant load in the quadratic
        // source pencil. Its work is removed again when restoring the original
        // spring energy account; it is never counted as authored forcing.
        let rest_tone = if network
            .connections
            .iter()
            .any(|link| link.stiffness_n_m != 0.0 && link.rest_extension_m != 0.0)
        {
            Some(if let Some(index) = unique.iter().position(|&w| w == 0.0) {
                index
            } else {
                unique.push(0.0);
                groups.push(Vec::new());
                unique.len() - 1
            })
        } else {
            None
        };
        let n = count
            .checked_mul(2)
            .and_then(|n| n.checked_add(network.material_memory.len()))
            .and_then(|n| n.checked_add(network.radiation_state.len()))
            .ok_or_else(|| invalid("complete coupled tone state extent"))?;
        let coefficients = count
            .checked_mul(unique.len())
            .ok_or_else(|| invalid("complete coupled tone force extent"))?;
        let force_coordinates = coefficients
            .checked_mul(if affine { 4 } else { 2 })
            .ok_or_else(|| invalid("complete coupled tone quadrature extent"))?;
        let augmented = n
            .checked_add(force_coordinates)
            .and_then(|v| v.checked_add(1))
            .ok_or_else(|| invalid("complete coupled tone generator extent"))?;
        let heat_dimension = n
            .checked_mul(2)
            .and_then(|v| v.checked_add(force_coordinates))
            .ok_or_else(|| invalid("complete coupled tone heat extent"))?;
        let channels = 3usize
            .checked_add(
                network
                    .connections
                    .iter()
                    .filter(|link| link.damping_n_s_m != 0.0)
                    .count(),
            )
            .ok_or_else(|| invalid("complete coupled tone loss channel extent"))?;
        // Admit the shared source assembler and its original channel witnesses
        // BEFORE allocating dense matrices. The core separately admits its full
        // tone/heat preparation within the remaining allowance.
        let square = n
            .checked_mul(n)
            .ok_or_else(|| invalid("coupled tone source square extent"))?;
        let source_bytes = square
            .checked_mul(16 + 2 * channels)
            .and_then(|v| v.checked_mul(core::mem::size_of::<f64>()))
            .and_then(|v| {
                v.checked_add(if affine {
                    coefficients.checked_mul(2 * core::mem::size_of::<f64>())?
                } else {
                    0
                })
            })
            .ok_or_else(|| invalid("coupled tone source byte extent"))?;
        let source_work = (n as u64)
            .checked_pow(3)
            .and_then(|v| v.checked_mul(128 + 4 * channels as u64))
            .and_then(|v| {
                v.checked_add(
                    (square as u64)
                        .checked_mul((network.connections.len() + channels + 16) as u64)?,
                )
            })
            .ok_or_else(|| invalid("coupled tone source work extent"))?;
        if augmented > fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || heat_dimension > fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
            || augmented > budget.max_dimension
            || heat_dimension > budget.max_dimension
            || source_bytes >= budget.max_dense_bytes
            || source_work >= budget.max_multiply_terms
        {
            return Err(invalid(
                "complete coupled tone source exceeds its dimension/byte/work allowance; no coordinate was removed",
            ));
        }
        let modes: Vec<_> = network
            .models
            .iter()
            .flat_map(|model| model.modes().iter().copied())
            .collect();
        let damping_rates: Vec<_> = network.models.iter()
            .flat_map(|model| model.damping_rates_per_s().iter().copied()).collect();
        let poles: Vec<_> = network.radiation.as_ref().map_or_else(Vec::new, |air| {
            air.model
                .poles
                .iter()
                .map(|pole| {
                    let mut coupling = vec![0.0; count];
                    for (&index, &gain) in air.source_modes.iter().zip(&pole.coupling) {
                        coupling[index] = gain;
                    }
                    crate::render::plate::impact::radiation::Pole {
                        omega: pole.omega,
                        zeta: pole.zeta,
                        coupling,
                    }
                })
                .collect()
        });
        // Preserve an anchored dashpot's original coupled-basis direction.
        // Springs/rest lengths continue to use their own original columns.
        let damping_columns: Vec<_> = (0..network.connections.len())
            .map(|index| {
                network
                    .connection_damping_column(index)
                    .map(<[f64]>::to_vec)
            })
            .collect();
        if network
            .connections
            .iter()
            .zip(&damping_columns)
            .any(|(link, column)| {
                link.constant_loss_reference_omega_rad_s.is_some() && column.is_none()
            })
        {
            return Err(invalid(
                "complete original constant-loss damping columns required",
            ));
        }
        let base = ExactBlock::new_with_damping_columns(
            &modes,
            &damping_rates,
            &network.columns,
            &network.connections,
            network.relaxation_branches.clone(),
            &poles,
            Some(&damping_columns),
            gate,
        )?;
        let mut witnesses = Vec::with_capacity(channels);
        witnesses.push(diagonal_witness(&base.component_r, n)?);
        witnesses.push(diagonal_witness(&base.material_r, n)?);
        witnesses.push(diagonal_witness(&base.radiation_r, n)?);
        let mut link_channels = Vec::with_capacity(network.connections.len());
        for row in &base.link_loss_rows {
            poll(Some(gate))?;
            link_channels.push(if let Some(row) = row {
                let index = witnesses.len();
                let mut matrix = vec![0.0; square];
                for i in 0..n {
                    for j in 0..n {
                        matrix[i * n + j] = finite(row[i] * row[j])?;
                    }
                }
                witnesses.push(
                    ResistanceWitness::try_from_source_rows(n, matrix, row.clone())
                        .map_err(ModalCouplingError::Propagation)?,
                );
                Some(index)
            } else {
                None
            });
        }
        let rest_force = base.rest_force;
        let split_budget = LinearSplitBudget {
            max_modes: n,
            max_omega_dt: core::f64::consts::PI,
        };
        let dense_budget = DenseExponentialBudget {
            max_dense_bytes: budget.max_dense_bytes - source_bytes,
            max_multiply_terms: budget.max_multiply_terms - source_work,
            ..budget
        };
        let flow = if affine {
            PreparedExactTones::new_owned_affine_with_resistance_channels_controlled(
                base.system,
                &unique,
                network.dt,
                split_budget,
                dense_budget,
                &witnesses,
                || gate.is_requested(),
            )
        } else {
            PreparedExactTones::new_owned_with_resistance_channels_controlled(
                base.system,
                &unique,
                network.dt,
                split_budget,
                dense_budget,
                &witnesses,
                || gate.is_requested(),
            )
        }
        .map_err(|error| match error {
            fs_phs::PreparedStepError::Cancelled => ModalCouplingError::Cancelled,
            fs_phs::PreparedStepError::Solver(source) => ModalCouplingError::Propagation(source),
        })?;
        poll(Some(gate))?;
        Ok(Self {
            flow,
            frequencies: frequencies.to_vec(),
            groups,
            original_to_group,
            rest_tone,
            rest_force,
            link_channels,
            initial: vec![0.0; n],
            endpoint: vec![0.0; n],
            cosine: vec![0.0; coefficients],
            sine: vec![0.0; coefficients],
            cosine_end: vec![0.0; if affine { coefficients } else { 0 }],
            sine_end: vec![0.0; if affine { coefficients } else { 0 }],
            mean_velocity: vec![0.0; count],
            channel_heat: vec![0.0; channels],
            endpoint_external: vec![0.0; count],
            force_envelope: vec![0.0; count],
        })
    }

    pub(crate) fn angular_frequencies_rad_s(&self) -> &[f64] {
        &self.frequencies
    }

    pub(crate) fn is_affine(&self) -> bool {
        self.flow.is_affine()
    }

    /// A causal observer uses the exact same original source and grouping.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn prepare_observer(
        &self,
        states: usize,
        a: &[f64],
        b: &[f64],
        c: &[f64],
        d: &[f64],
        budget: DenseExponentialBudget,
        gate: &CancelGate,
    ) -> Result<fs_phs::linear_exact::PreparedToneObserver, ModalCouplingError> {
        fs_phs::linear_exact::PreparedToneObserver::new_controlled(
            &self.flow,
            states,
            1,
            a,
            b,
            c,
            d,
            budget,
            || gate.is_requested(),
        )
        .map_err(|error| match error {
            fs_phs::PreparedStepError::Cancelled => ModalCouplingError::Cancelled,
            fs_phs::PreparedStepError::Solver(source) => ModalCouplingError::Propagation(source),
        })
    }

    /// Preview only after stage has packed this candidate's exact forces and
    /// initial state. No physical state, work or loss channel is recomputed.
    pub(crate) fn preview_observer(
        &self,
        observer: &mut fs_phs::linear_exact::PreparedToneObserver,
        initial: &[f64],
        next: &mut [f64],
        output: &mut [f64],
        gate: Option<&CancelGate>,
    ) -> Result<(), ModalCouplingError> {
        let result = if self.flow.is_affine() {
            observer.preview_affine_into_controlled(
                &self.flow,
                &self.initial,
                initial,
                &self.cosine,
                &self.sine,
                &self.cosine_end,
                &self.sine_end,
                next,
                output,
                || gate.is_some_and(CancelGate::is_requested),
            )
        } else {
            observer.preview_into_controlled(
                &self.flow,
                &self.initial,
                initial,
                &self.cosine,
                &self.sine,
                next,
                output,
                || gate.is_some_and(CancelGate::is_requested),
            )
        };
        result.map_err(|error| match error {
            fs_phs::PreparedStepError::Cancelled => ModalCouplingError::Cancelled,
            fs_phs::PreparedStepError::Solver(source) => ModalCouplingError::Propagation(source),
        })
    }

    pub(crate) fn advance_force_phase(
        &self,
        tone: usize,
        cosine: f64,
        sine: f64,
    ) -> Result<(f64, f64), ModalCouplingError> {
        let &group = self
            .original_to_group
            .get(tone)
            .ok_or_else(|| invalid("original coupled tone index"))?;
        self.flow
            .advance_force_phase(group, cosine, sine)
            .map_err(ModalCouplingError::Propagation)
    }

    pub(crate) fn stage(
        &mut self,
        network: &mut CoupledModalSystem,
        cosine_forces: &[f64],
        sine_forces: &[f64],
        gate: Option<&CancelGate>,
    ) -> Result<(), ModalCouplingError> {
        self.stage_inner(network, cosine_forces, sine_forces, None, gate)
    }

    pub(crate) fn stage_affine(
        &mut self,
        network: &mut CoupledModalSystem,
        cosine_start: &[f64],
        sine_start: &[f64],
        cosine_end: &[f64],
        sine_end: &[f64],
        gate: Option<&CancelGate>,
    ) -> Result<(), ModalCouplingError> {
        self.stage_inner(
            network,
            cosine_start,
            sine_start,
            Some((cosine_end, sine_end)),
            gate,
        )
    }

    fn stage_inner(
        &mut self,
        network: &mut CoupledModalSystem,
        cosine_forces: &[f64],
        sine_forces: &[f64],
        end_forces: Option<(&[f64], &[f64])>,
        gate: Option<&CancelGate>,
    ) -> Result<(), ModalCouplingError> {
        poll(gate)?;
        let count = network.mode_count();
        let input_tones = self.frequencies.len();
        let tones = self.groups.len();
        if cosine_forces.len() != count * input_tones
            || sine_forces.len() != cosine_forces.len()
            || cosine_forces
                .iter()
                .chain(sine_forces)
                .any(|v| !v.is_finite())
            || self.flow.is_affine() != end_forces.is_some()
            || end_forces.is_some_and(|(c, s)| {
                c.len() != cosine_forces.len()
                    || s.len() != cosine_forces.len()
                    || c.iter().chain(s).any(|v| !v.is_finite())
            })
        {
            return Err(invalid(
                "complete finite tone endpoints must match the prepared source family",
            ));
        }
        let mut next = 0;
        for model in &network.models {
            for state in model.states() {
                self.initial[next] = state.displacement_m_sqrt_kg;
                self.initial[next + 1] = state.velocity_m_sqrt_kg_per_s;
                next += 2;
            }
        }
        self.initial[next..next + network.material_memory.len()]
            .copy_from_slice(&network.material_memory);
        next += network.material_memory.len();
        self.initial[next..].copy_from_slice(&network.radiation_state);
        for mode in 0..count {
            poll(gate)?;
            let mut endpoint_sum = Sum::default();
            let mut envelope = 0.0;
            for (tone, group) in self.groups.iter().enumerate() {
                let mut cosine = Sum::default();
                let mut sine = Sum::default();
                for &original in group {
                    cosine.add(cosine_forces[mode * input_tones + original])?;
                    sine.add(sine_forces[mode * input_tones + original])?;
                }
                let fc = cosine.value()?;
                let fs = sine.value()?;
                // Coalesce each endpoint before differencing in the core.
                // A small envelope survives large exactly cancelling tones.
                let (fc_end, fs_end) = if let Some((c, s)) = end_forces {
                    let mut cosine = Sum::default();
                    let mut sine = Sum::default();
                    for &original in group {
                        cosine.add(c[mode * input_tones + original])?;
                        sine.add(s[mode * input_tones + original])?;
                    }
                    (cosine.value()?, sine.value()?)
                } else {
                    (fc, fs)
                };
                let (end, _) = self
                    .flow
                    .advance_force_phase(tone, fc_end, fs_end)
                    .map_err(ModalCouplingError::Propagation)?;
                endpoint_sum.add(end)?;
                let total_fc = finite(
                    fc + if self.rest_tone == Some(tone) {
                        self.rest_force[mode]
                    } else {
                        0.0
                    },
                )?;
                self.cosine[mode * tones + tone] = total_fc;
                self.sine[mode * tones + tone] = fs;
                let mut amplitude = if self.flow.angular_frequencies_rad_s()[tone] == 0.0 {
                    total_fc.abs()
                } else {
                    det::hypot(total_fc, fs)
                };
                if end_forces.is_some() {
                    let total_end = finite(
                        fc_end
                            + if self.rest_tone == Some(tone) {
                                self.rest_force[mode]
                            } else {
                                0.0
                            },
                    )?;
                    self.cosine_end[mode * tones + tone] = total_end;
                    self.sine_end[mode * tones + tone] = fs_end;
                    // Norm is convex along each affine coefficient segment.
                    // Endpoint maxima therefore bound the whole force interval.
                    amplitude =
                        amplitude.max(if self.flow.angular_frequencies_rad_s()[tone] == 0.0 {
                            total_end.abs()
                        } else {
                            det::hypot(total_end, fs_end)
                        });
                }
                envelope = finite(envelope + amplitude)?;
            }
            self.endpoint_external[mode] = endpoint_sum.value()?;
            self.force_envelope[mode] = envelope;
        }
        let record = if end_forces.is_some() {
            self.flow.step_affine_observed_into(
                &self.initial,
                &self.cosine,
                &self.sine,
                &self.cosine_end,
                &self.sine_end,
                &mut self.endpoint,
                &mut self.mean_velocity,
                &mut self.channel_heat,
            )
        } else {
            self.flow.step_observed_into(
                &self.initial,
                &self.cosine,
                &self.sine,
                &mut self.endpoint,
                &mut self.mean_velocity,
                &mut self.channel_heat,
            )
        }
        .map_err(ModalCouplingError::Propagation)?;
        poll(gate)?;
        let mut connection_loss = 0.0;
        let mut transfers = Vec::with_capacity(network.connections.len());
        for ((column, law), channel) in network
            .columns
            .iter()
            .zip(&network.connections)
            .zip(&self.link_channels)
        {
            poll(gate)?;
            let heat = channel.map_or(0.0, |i| self.channel_heat[i]);
            connection_loss = finite(connection_loss + heat)?;
            let mut before = -law.rest_extension_m;
            let mut after = before;
            for (mode, &gain) in column.iter().enumerate() {
                before = finite(before + gain * self.initial[2 * mode])?;
                after = finite(after + gain * self.endpoint[2 * mode])?;
            }
            transfers.push(crate::port_adoption::PortTransfer::work(
                finite(-0.5 * law.stiffness_n_m * (after - before) * (after + before) - heat)?,
                false,
            ));
        }
        let losses = [
            self.channel_heat[0],
            connection_loss,
            self.channel_heat[1],
            self.channel_heat[2],
        ];
        let independent_total = finite(losses.iter().sum())?;
        if (independent_total - record.dissipated).abs()
            > CHANNEL_SUM_RELATIVE_TOLERANCE
                * independent_total
                    .abs()
                    .max(record.dissipated.abs())
                    .max(f64::MIN_POSITIVE)
        {
            return Err(invalid(
                "continuous coupled source loss channels disagree with total heat",
            ));
        }
        let mut rest_work = Sum::default();
        for (mode, &rest) in self.rest_force.iter().enumerate() {
            rest_work.add(finite(
                rest * (self.endpoint[2 * mode] - self.initial[2 * mode]),
            )?)?;
        }
        let external_work = finite(record.supplied - rest_work.value()?)?;
        let initial_bound_energy = finite(self.flow.system().hamiltonian(&self.initial))?;
        // The original full-interval force bound receives a conservative bound
        // on the actual time-varying effort, not merely its endpoint sample.
        // The same shared admission owns all endpoint/energy/pressure/port gates.
        stage_candidate(
            network,
            &self.endpoint_external,
            PropagatedCandidate {
                endpoint: self.endpoint.clone(),
                start: self.initial.clone(),
                input: self.force_envelope.clone(),
                losses,
                initial_bound_energy,
                external_work: Some(external_work),
                connection_transfers: Some(transfers),
            },
            gate,
        )
    }
}

fn diagonal_witness(matrix: &[f64], n: usize) -> Result<ResistanceWitness, ModalCouplingError> {
    let matrix = if matrix.is_empty() {
        vec![0.0; n * n]
    } else {
        matrix.to_vec()
    };
    let mut rows = Vec::new();
    for i in 0..n {
        if matrix[i * n + i] == 0.0 {
            continue;
        }
        let start = rows.len();
        rows.resize(start + n, 0.0);
        rows[start + i] = finite(det::sqrt(matrix[i * n + i]))?;
    }
    ResistanceWitness::try_from_source_rows(n, matrix, rows)
        .map_err(ModalCouplingError::Propagation)
}

/// Original-order compensated addition, including exact opposite sources.
#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), ModalCouplingError> {
        let next = finite(self.sum + value)?;
        self.correction = finite(
            self.correction
                + if self.sum.abs() >= value.abs() {
                    (self.sum - next) + value
                } else {
                    (value - next) + self.sum
                },
        )?;
        self.sum = next;
        Ok(())
    }
    fn value(&self) -> Result<f64, ModalCouplingError> {
        finite(self.sum + self.correction)
    }
}
