//! Joint held-wave closure at original physical ports of one complete system.
//! Reference loading is boundary exchange, never constitutive absorption.
use super::{
    DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM, finite,
    positive_heat::PositiveHeat, prepared::PreparedQuadraticFlow, refuse, split_quadratic,
};
use crate::{PhsError, PortHamiltonian, storage_energy_change};
use fs_math::det;

/// Independent complete-system accounts over the prepared interval [J].
#[derive(Clone, Copy, Debug)]
pub struct MultiportScatteringStepRecord {
    /// Change in original constitutive stored energy [J].
    pub delta_h: f64,
    /// Original physical resistance only; excludes every reference impedance.
    pub physical_heat: f64,
    /// Work at unselected original ports.
    pub external_work: f64,
    /// Sum of the selected physical boundary works [J].
    pub boundary_work: f64,
    /// Sum of the independently integrated incident wave energies [J].
    pub incident_energy: f64,
    /// Sum of the independently integrated outgoing wave energies [J].
    pub outgoing_energy: f64,
}

/// Independent accounts for one selected port, in the declared selection order.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ScatteringPortRecord {
    /// Integral of the original physical boundary power u*y [J].
    pub boundary_work: f64,
    /// Integral of a*a/Z [J].
    pub incident_energy: f64,
    /// Independently integrated (y-a)^2/Z [J], never an energy difference.
    pub outgoing_energy: f64,
    /// Interval MEAN b=y-a, in original output units. Its square is not energy.
    pub outgoing_wave: f64,
}

struct ReferencePort {
    inverse: f64,
    root_inverse: f64,
    active: bool,
    energy: PositiveHeat,
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
    accounts: Vec<ScatteringPortRecord>,
}

/// Exact complete quadratic flow with simultaneous scattering boundaries.
///
/// Original power is u*y. Each selected input slot holds an incident wave a
/// in its ORIGINAL output's units, imposing u=(2a-y)/Z and b=y-a. Unselected
/// inputs and every returned mean output retain their original meaning/order.
/// All selected ports load ONE coupled generator; sequential scalar closures
/// would propagate a different physical system.
///
/// Preparation retains all states, storage, physical damping and port coupling.
/// Existing complete-state/port and dense-exponential capacities apply; no
/// modes are discarded to fit them. Each port retains a positive factor of its
/// own outgoing-wave integral, including small representable matched reflection.
/// Cold factor work grows with the selected-port count. Warm calls reuse all
/// owner scratch and never exponentiate, infer materials, or allocate.
///
/// This owner holds numerical scratch, not accepted physical history. Consumer
/// source bandwidth, spatial convergence and nonlinear closure remain explicit
/// responsibilities. Refusal publishes none of the caller's candidate buffers.
pub struct PreparedMultiportScatteringZoh {
    original: PortHamiltonian,
    robin: PortHamiltonian,
    flow: PreparedQuadraticFlow,
    physical_heat: PositiveHeat,
    ports: Vec<(usize, f64)>,
    selected: Vec<bool>,
    references: Vec<ReferencePort>,
    dt: f64,
    mean_state: Vec<f64>,
    mean_input: Vec<f64>,
    direct: Vec<f64>,
    scratch: Scratch,
}

impl PreparedMultiportScatteringZoh {
    /// Own one original system and fixed positive-duration cache.
    /// Selected slots must be distinct, with positive representable impedances.
    pub fn new_owned(
        mut system: PortHamiltonian,
        ports: &[(usize, f64)],
        duration: f64,
        budget: LinearSplitBudget,
    ) -> Result<Self, PhsError> {
        if ports.is_empty() || !duration.is_finite() || duration <= 0.0 {
            return Err(refuse(
                "scattering requires selected ports and a positive clock",
            ));
        }
        for &(port, impedance) in ports {
            if port >= system.m {
                return Err(PhsError::BadPortPairing);
            }
            let inverse = 1.0 / impedance;
            if !impedance.is_finite() || impedance <= 0.0 || !inverse.is_finite() || inverse <= 0.0
            {
                return Err(refuse(
                    "scattering impedance must be positive and representable",
                ));
            }
        }
        // Admit complete dimensions before allocating port-dependent caches.
        let original_split = split_quadratic(&system, budget)?;
        if !original_split.is_exact_quadratic() {
            return Err(refuse(
                "scattering requires declared identically quadratic storage",
            ));
        }
        let (n, m, k) = (system.n, system.m, ports.len());
        let mut selected = vec![false; m];
        for &(port, _) in ports {
            if selected[port] {
                return Err(refuse("scattering selected ports must be distinct"));
            }
            selected[port] = true;
        }
        let mut robin = system.share_storage();
        for &(port, impedance) in ports {
            let inverse = 1.0 / impedance;
            let root_inverse = det::sqrt(inverse);
            for i in 0..n {
                let root_i = root_inverse * system.g[i * m + port];
                for j in 0..n {
                    robin.r[i * n + j] += root_i * (root_inverse * system.g[j * m + port]);
                }
                robin.g[i * m + port] = 2.0 * inverse * system.g[i * m + port];
                if system.g[i * m + port] != 0.0
                    && (root_i * root_i == 0.0 || robin.g[i * m + port] == 0.0)
                {
                    return Err(refuse("scattering port scale is unrepresentable"));
                }
            }
        }
        if !finite(&robin.r) || !finite(&robin.g) {
            return Err(refuse("scattering boundary structure is nonfinite"));
        }
        // Every reference factor below uses this SAME jointly closed generator.
        let split = split_quadratic(&robin, budget)?;
        let q = split.quadratic_matrix();
        let outputs: Vec<Vec<f64>> = (0..m)
            .map(|p| {
                (0..n)
                    .map(|j| (0..n).map(|i| system.g[i * m + p] * q[i * n + j]).sum())
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
        let physical_heat = PositiveHeat::new(&robin, &split, duration, &system.r)?;
        let mut references = Vec::with_capacity(k);
        let mut mean_state = Vec::with_capacity(k * n);
        let mut mean_input = Vec::with_capacity(k * m);
        for &(port, impedance) in ports {
            let inverse = 1.0 / impedance;
            let root_inverse = det::sqrt(inverse);
            // PositiveHeat's first factor row is the original y/sqrt(Z).
            let mut root = vec![0.0; n * n];
            for i in 0..n {
                root[i] = root_inverse * system.g[i * m + port];
            }
            let active = root.iter().any(|value| *value != 0.0);
            let mut resistance = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    resistance[i * n + j] = root[i] * root[j];
                }
            }
            let energy = PositiveHeat::from_factor(&robin, &split, duration, &resistance, root)?;
            references.push(ReferencePort {
                inverse,
                root_inverse,
                active,
                energy,
            });
            let (state_map, input_map) = flow.linear_integral_map(port)?;
            mean_state.extend(state_map.into_iter().map(|value| value / duration));
            mean_input.extend(input_map.into_iter().map(|value| value / duration));
        }
        if !finite(&mean_state) || !finite(&mean_input) {
            return Err(refuse("scattering mean output map is nonfinite"));
        }
        let mut direct = vec![0.0; k * k];
        for row in 0..k {
            for (col, &(port, _)) in ports.iter().enumerate() {
                direct[row * k + col] = mean_input[row * m + port];
            }
        }
        Ok(Self {
            original: system,
            robin,
            flow,
            physical_heat,
            ports: ports.to_vec(),
            selected,
            references,
            dt: duration,
            mean_state,
            mean_input,
            direct,
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
                accounts: vec![ScatteringPortRecord::default(); k],
            },
        })
    }

    /// Original physical system, excluding all reference loading from its R.
    pub fn system(&self) -> &PortHamiltonian {
        &self.original
    }
    /// Fixed prepared interval [s].
    pub fn duration(&self) -> f64 {
        self.dt
    }
    /// Original slot and impedance, in the returned account/map order.
    pub fn ports(&self) -> &[(usize, f64)] {
        &self.ports
    }

    /// Current-window relation y_bar=free+direct*a at selected ORIGINAL outputs.
    /// The row-major direct matrix and free vector follow ports() order.
    /// Selected input slots are ignored in free; other authored inputs remain.
    /// Refusal preserves both caller buffers; this call allocates no heap memory.
    pub fn average_relation_into(
        &self,
        state: &[f64],
        inputs: &[f64],
        free: &mut [f64],
        direct: &mut [f64],
    ) -> Result<(), PhsError> {
        self.admit_layout(state, inputs)?;
        let (n, m, k) = (self.original.n, self.original.m, self.ports.len());
        if free.len() != k || direct.len() != k * k {
            return Err(refuse(
                "scattering relation requires complete selected-port buffers",
            ));
        }
        // The admitted 2*(n+m+1)<=256 bound also bounds k, before this cache exists.
        let mut candidate = [0.0; MAX_LINEAR_EXPONENTIAL_BLOCK_DIM];
        for (row, value) in candidate[..k].iter_mut().enumerate() {
            *value = self.mean_state[row * n..(row + 1) * n]
                .iter()
                .zip(state)
                .map(|(a, x)| a * x)
                .sum::<f64>()
                + self.mean_input[row * m..(row + 1) * m]
                    .iter()
                    .zip(inputs)
                    .enumerate()
                    .filter(|(p, _)| !self.selected[*p])
                    .map(|(_, (a, u))| a * u)
                    .sum::<f64>();
        }
        if !finite(&candidate[..k]) {
            return Err(refuse("scattering free mean output is nonfinite"));
        }
        free.copy_from_slice(&candidate[..k]);
        direct.copy_from_slice(&self.direct);
        Ok(())
    }

    /// Complete original endpoint/mean outputs plus independent selected-port
    /// accounts. Publication is transactional across ALL three output buffers.
    pub fn step_into(
        &mut self,
        state: &[f64],
        inputs: &[f64],
        next: &mut [f64],
        mean_outputs: &mut [f64],
        port_accounts: &mut [ScatteringPortRecord],
    ) -> Result<MultiportScatteringStepRecord, PhsError> {
        if next.len() != self.original.n
            || mean_outputs.len() != self.original.m
            || port_accounts.len() != self.ports.len()
        {
            return Err(refuse(
                "scattering output buffers require complete original layout",
            ));
        }
        let account = self.evaluate(state, inputs)?;
        next.copy_from_slice(&self.scratch.endpoint);
        mean_outputs.copy_from_slice(&self.scratch.output);
        port_accounts.copy_from_slice(&self.scratch.accounts);
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
    ) -> Result<MultiportScatteringStepRecord, PhsError> {
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
        let external_work: f64 = inputs
            .iter()
            .zip(&s.integrals)
            .enumerate()
            .filter(|(p, _)| !self.selected[*p])
            .map(|(_, (u, integral))| u * integral)
            .sum();
        let delta_h = storage_energy_change(
            self.original.storage.as_ref(),
            state,
            &s.endpoint,
            &mut s.midpoint,
            &mut s.effort,
        );
        let mut scale =
            super::prepared_zoh::energetic_scale_into(&self.original, state, &mut s.effort)
                .max(super::prepared_zoh::energetic_scale_into(
                    &self.original,
                    &s.endpoint,
                    &mut s.effort,
                ))
                .max(delta_h.abs())
                .max(physical_heat.abs())
                .max(external_work.abs());
        let mut boundary_work = 0.0;
        let mut incident_energy = 0.0;
        let mut outgoing_energy = 0.0;
        for (slot, (&(port, _), reference)) in self.ports.iter().zip(&self.references).enumerate() {
            let reference_energy = reference.energy.evaluate_into(
                &self.robin,
                state,
                inputs,
                &mut s.effort,
                &mut s.derivative,
                &mut s.heat_coordinates,
            )?;
            let a = inputs[port];
            let shifted = a * reference.root_inverse;
            let incident_amplitude = shifted * det::sqrt(self.dt);
            let incoming = incident_amplitude * incident_amplitude;
            let outgoing = if reference.active {
                reference.energy.evaluate_shifted_into(
                    &self.robin,
                    state,
                    inputs,
                    shifted,
                    &mut s.effort,
                    &mut s.derivative,
                    &mut s.heat_coordinates,
                )?
            } else {
                incoming
            };
            let selected_work = (2.0 * reference.inverse * a) * s.integrals[port];
            let work = selected_work - reference_energy;
            let outgoing_wave = s.output[port] - a;
            scale = scale.max(selected_work.abs()).max(reference_energy.abs());
            let wave_scale = scale.max(incoming.abs()).max(outgoing.abs());
            let wave_balance = work - (incoming - outgoing);
            if !finite(&[
                reference_energy,
                incoming,
                outgoing,
                selected_work,
                work,
                outgoing_wave,
                wave_scale,
                wave_balance,
            ]) || outgoing < 0.0
                || wave_balance.abs() > 1e-12 * wave_scale
            {
                return Err(refuse("scattering independent port wave balance refused"));
            }
            s.accounts[slot] = ScatteringPortRecord {
                boundary_work: work,
                incident_energy: incoming,
                outgoing_energy: outgoing,
                outgoing_wave,
            };
            boundary_work += work;
            incident_energy += incoming;
            outgoing_energy += outgoing;
        }
        scale = scale.max(boundary_work.abs());
        let balance = delta_h + physical_heat - external_work - boundary_work;
        let wave_scale = scale.max(incident_energy.abs()).max(outgoing_energy.abs());
        let wave_balance = boundary_work - (incident_energy - outgoing_energy);
        if !finite(&s.output)
            || !self.original.hamiltonian(&s.endpoint).is_finite()
            || !finite(&[
                physical_heat,
                external_work,
                delta_h,
                boundary_work,
                incident_energy,
                outgoing_energy,
                scale,
                balance,
                wave_scale,
                wave_balance,
            ])
            || physical_heat < 0.0
            || balance.abs() > 1e-12 * scale
            || wave_balance.abs() > 1e-12 * wave_scale
        {
            return Err(refuse(
                "scattering original physical work and wave balance refused",
            ));
        }
        Ok(MultiportScatteringStepRecord {
            delta_h,
            physical_heat,
            external_work,
            boundary_work,
            incident_energy,
            outgoing_energy,
        })
    }
}
