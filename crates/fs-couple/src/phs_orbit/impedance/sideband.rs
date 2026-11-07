//! Loaded original-state probe response through the existing signed-sideband solver.
//!
//! The base and every shifted boundary port belong to one explicitly authored
//! physical source. Samples are not interpolated or realized in time. Real-system
//! frequency symmetry converts exp(-iwt) data at |w| to signed exp(+iwt) ports;
//! this is neither a causality certificate nor a Floquet stability calculation.

use super::{
    C64, ComplexSum, Cx, HbOrbit, OrbitProblem, PeriodicPortImpedance, PhsImpedanceOrbit,
    PortImpedanceError, PreparedPhsImpedance, Result, admit_passive_port_impedance, det, finite,
    invalid, physical_operator, poll,
};
use crate::phs_orbit::sideband::{
    self as physical, ObservationBudget, PhsPortSidebandCarrierDerivativeResponse,
    PhsPortSidebandControls, PhsPortSidebandResponse,
};
use fs_orbit::sideband::{
    PreparedPeriodicLinearization, PreparedSidebandResponse, SidebandResponse, SidebandScales,
    sensitivity::{SidebandCarrierDerivativeControls, SidebandCarrierDerivativeResponse},
};

mod frequency;
use frequency::signed_frequency_entry;

/// Explicit nonnegative-frequency data, in the same physical chart as the base load.
#[derive(Clone, Debug)]
pub struct PortImpedanceFrequencySample {
    /// Nonnegative physical angular frequency, in radians per unit time.
    pub angular_frequency: f64,
    /// Complete row-major exp(-iwt) impedance; DC must be real.
    pub impedance: Vec<C64>,
    /// Authored dZ_-(w)/dw for w >= 0, not a difference between table entries.
    /// A two-sided real-system derivative at DC must be purely imaginary.
    /// None permits responses but refuses derivative queries at this frequency.
    pub angular_derivative: Option<Vec<C64>>,
}

/// Additional table and physical-map resources; generic sideband and observer
/// caps remain separate. Caller/source storage and allocator metadata are excluded.
#[derive(Clone, Copy, Debug)]
pub struct ImpedanceSidebandLimits {
    /// Authored absolute-frequency samples, additionally capped at 4096.
    pub maximum_frequency_samples: usize,
    /// Distinct declared signed carriers, additionally capped at 128.
    pub maximum_carriers: usize,
    /// Pointwise physical admission plus original dense G*Z and G*Z' maps.
    pub maximum_products: u64,
    /// Retained table/map/base payload and physical admission scratch; hard cap 256 MiB.
    pub maximum_workspace_bytes: usize,
}

/// One admitted finite probe band, reusable across its declared signed carriers.
///
/// `carriers` and sample frequencies must be strictly increasing. Every retained
/// |carrier+k*omega| must be supplied here or exactly coincide with an admitted
/// base harmonic. Samples at shared base frequencies must agree with the base
/// matrix entry-for-entry. No constant-load fallback or tolerance-based frequency
/// matching is used. Additional samples can provide derivatives at base harmonics.
pub struct SidebandPortImpedance<'load> {
    base: &'load PeriodicPortImpedance,
    carriers: Vec<f64>,
    harmonics: usize,
    samples: Vec<PortImpedanceFrequencySample>,
    limits: ImpedanceSidebandLimits,
}

impl<'load> SidebandPortImpedance<'load> {
    /// Admit complete source coverage and pointwise passivity before any original
    /// storage/Hessian evaluation. Values and optional slopes must describe the
    /// SAME physical load; admission cannot establish continuity or causality.
    pub fn new(
        base: &'load PeriodicPortImpedance,
        carriers: Vec<f64>,
        harmonics: usize,
        samples: Vec<PortImpedanceFrequencySample>,
        limits: ImpedanceSidebandLimits,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        let p = base.ports;
        if harmonics > 64
            || carriers.is_empty()
            || carriers.len() > limits.maximum_carriers.min(128)
            || samples.len() > limits.maximum_frequency_samples.min(4096)
            || carriers.iter().any(|v| !v.is_finite())
            || carriers.windows(2).any(|w| w[0] >= w[1])
            || samples
                .windows(2)
                .any(|w| w[0].angular_frequency >= w[1].angular_frequency)
            || samples.iter().any(|sample| {
                let w = sample.angular_frequency;
                !w.is_finite()
                    || w < 0.0
                    || sample.impedance.len() != p * p
                    || sample
                        .angular_derivative
                        .as_ref()
                        .is_some_and(|slope| slope.len() != p * p)
            })
        {
            return Err(invalid(
                "bounded sorted carriers/frequencies, complete finite matrices and real DC symmetry required",
            ));
        }
        let table = Self {
            base,
            carriers,
            harmonics,
            samples,
            limits,
        };
        table.resources(0, 0)?;
        for sample in &table.samples {
            poll(cx)?;
            if sample
                .impedance
                .iter()
                .any(|&z| !finite(z) || (sample.angular_frequency == 0.0 && z.im != 0.0))
                || sample.angular_derivative.as_ref().is_some_and(|slope| {
                    slope
                        .iter()
                        .any(|&z| !finite(z) || (sample.angular_frequency == 0.0 && z.re != 0.0))
                })
            {
                return Err(invalid(
                    "finite physical values and two-sided real DC symmetry required",
                ));
            }
            if let Some(h) = base_harmonic(base, sample.angular_frequency) {
                if sample.impedance != base.matrices[h] {
                    return Err(invalid(
                        "probe and base impedance disagree at a shared physical frequency",
                    ));
                }
            }
        }
        // Coverage precedes PSD callbacks and all later physical source work.
        for &carrier in &table.carriers {
            poll(cx)?;
            for k in -(harmonics as isize)..=harmonics as isize {
                let w = (k as f64).mul_add(base.omega, carrier);
                if !w.is_finite() || table.matrix(w.abs()).is_none() {
                    return Err(invalid(
                        "every declared signed sideband needs an explicit physical load value",
                    ));
                }
            }
        }
        for sample in &table.samples {
            poll(cx)?;
            admit_passive_port_impedance(p, &sample.impedance, cx).map_err(|source| {
                if source == fs_modal::harmonic::HarmonicResponseError::Cancelled {
                    PortImpedanceError::Cancelled
                } else {
                    PortImpedanceError::FrequencyAdmission {
                        angular_frequency: sample.angular_frequency,
                        source,
                    }
                }
            })?;
        }
        Ok(table)
    }

    fn sample(&self, frequency: f64) -> Option<usize> {
        self.samples
            .binary_search_by(|sample| sample.angular_frequency.total_cmp(&frequency))
            .ok()
            .or_else(|| {
                (frequency == 0.0
                    && self
                        .samples
                        .first()
                        .is_some_and(|s| s.angular_frequency == 0.0))
                .then_some(0)
            })
    }
    fn matrix(&self, frequency: f64) -> Option<&[C64]> {
        self.sample(frequency)
            .map(|i| self.samples[i].impedance.as_slice())
            .or_else(|| {
                base_harmonic(self.base, frequency).map(|h| self.base.matrices[h].as_slice())
            })
    }
    fn derivative(&self, frequency: f64) -> Option<&[C64]> {
        self.samples
            .get(self.sample(frequency)?)?
            .angular_derivative
            .as_deref()
    }
    fn resources(&self, states: usize, base_width: usize) -> Result<()> {
        let (q, p, d) = (
            self.samples.len() as u128,
            self.base.ports as u128,
            (states + self.base.ports) as u128,
        );
        let slopes = self
            .samples
            .iter()
            .filter(|s| s.angular_derivative.is_some())
            .count() as u128;
        let maps = if states == 0 { 0 } else { q + slopes };
        let work = 8 * q * p * p * p + maps * (states as u128 * p * p + d * d);
        let retained_entries: u128 = self
            .samples
            .iter()
            .map(|sample| {
                sample.impedance.capacity() as u128
                    + sample
                        .angular_derivative
                        .as_ref()
                        .map_or(0, |slope| slope.capacity() as u128)
            })
            .sum();
        let bytes = 16 * retained_entries
            + 8 * (self.samples.capacity() as u128 + self.carriers.capacity() as u128)
            + 16 * maps * d * d
            + 16 * d * base_width as u128
            + 256 * p * p
            + 128 * p;
        if work > self.limits.maximum_products as u128
            || bytes > self.limits.maximum_workspace_bytes.min(256 << 20) as u128
        {
            return Err(PortImpedanceError::Budget(
                "loaded sideband physical table/map resources exceeded",
            ));
        }
        Ok(())
    }
}

fn base_harmonic(base: &PeriodicPortImpedance, frequency: f64) -> Option<usize> {
    let h = (frequency / base.omega).round();
    (frequency.is_finite()
        && frequency >= 0.0
        && h <= base.harmonic_count() as f64
        && h * base.omega == frequency)
        .then_some(h as usize)
}

struct LoadedProblem<'loaded, 'model, 'scope> {
    loaded: &'loaded PreparedPhsImpedance<'model, 'scope>,
    table: &'loaded SidebandPortImpedance<'model>,
    operators: Vec<Vec<C64>>,
    derivatives: Vec<Option<Vec<C64>>>,
}
impl LoadedProblem<'_, '_, '_> {
    fn state_jacobian(&self, time: f64, state: &[f64], out: &mut [f64]) -> bool {
        self.loaded.state_jacobian(time, state, out)
    }
    fn derivative(&self, frequency: f64, out: &mut [C64]) -> bool {
        let Some(index) = self.table.sample(frequency.abs()) else {
            return false;
        };
        let Some(block) = &self.derivatives[index] else {
            return false;
        };
        if out.len() != block.len() || !self.loaded.source_poll() {
            return false;
        }
        for (dst, &value) in out.iter_mut().zip(block) {
            *dst = signed_frequency_entry(value, frequency, true);
        }
        true
    }
}
impl OrbitProblem for LoadedProblem<'_, '_, '_> {
    fn dim(&self) -> usize {
        self.loaded.dim()
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        self.loaded.island(time, state, out);
    }
    fn island_state_jacobian(&self) -> Option<fs_orbit::IslandStateJacobian<Self>> {
        self.loaded
            .original
            .tangent
            .available()
            .then_some(Self::state_jacobian)
    }
    fn port(&self, s: C64) -> Vec<C64> {
        if !finite(s) || s.re != 0.0 || !self.loaded.source_poll() {
            return vec![C64::new(f64::NAN, f64::NAN); self.dim() * self.dim()];
        }
        let block = if let Some(index) = self.table.sample(s.im.abs()) {
            &self.operators[index]
        } else if let Some(h) = base_harmonic(self.loaded.load, s.im.abs()) {
            &self.loaded.operators[h]
        } else {
            return vec![C64::new(f64::NAN, f64::NAN); self.dim() * self.dim()];
        };
        block
            .iter()
            .map(|&z| signed_frequency_entry(z, s.im, false))
            .collect()
    }
    // A sampled port has no inferred time-domain RHS or Floquet realization.
}

/// Complete original physical response; generic residuals still assess the
/// augmented [state; flow] equations, including the boundary constraint.
#[derive(Clone, Debug)]
pub struct PhsImpedanceSidebandResponse {
    /// Original states and independently refined G^T Hxx dx physical output.
    /// This feeds the existing acoustic sideband observer without changing units.
    pub response: PhsPortSidebandResponse,
    /// Retained auxiliary boundary flow, sideband-major at the STATE order.
    pub flow_coefficients: Vec<C64>,
    /// Opposing boundary effort Z(w)*flow, in the same positive-time peak chart.
    pub effort_coefficients: Vec<C64>,
    /// Per-port L2 discrepancy from the independently observed physical flow.
    pub retained_flow_mismatch: Vec<f64>,
}

/// Fixed-base derivative with the same LU and original physical observer.
#[derive(Clone, Debug)]
pub struct PhsImpedanceSidebandCarrierDerivativeResponse {
    /// Original physical response/slope, compatible with the acoustic observer.
    /// Both residuals still check the complete augmented equations.
    pub response: PhsPortSidebandCarrierDerivativeResponse,
    /// Retained auxiliary flow in the state-order positive-time peak chart.
    pub flow_coefficients: Vec<C64>,
    /// Retained opposing boundary effort in the same chart and order.
    pub effort_coefficients: Vec<C64>,
    /// Derivative of the auxiliary flow per signed carrier angular frequency.
    pub flow_angular_derivative: Vec<C64>,
    /// Includes both Z*y' and Z'*y, with signed-frequency chain rules.
    pub effort_angular_derivative: Vec<C64>,
    /// Per-port primal auxiliary-versus-original-observation L2 discrepancy.
    pub retained_flow_mismatch: Vec<f64>,
    /// Corresponding derivative discrepancy in port units per angular frequency.
    pub derivative_flow_mismatch: Vec<f64>,
}

/// Original loaded base plus one pre-admitted finite probe band. Construction
/// does no storage work. Preparation independently reassesses the full base
/// equations and the native full Hessian; public result diagnostics are untrusted.
/// The exclusive loaded-source borrow permits one live physical map cache.
pub struct PhsImpedanceSidebandSource<'loaded, 'model, 'scope> {
    problem: LoadedProblem<'loaded, 'model, 'scope>,
    orbit: &'loaded PhsImpedanceOrbit,
    augmented: HbOrbit,
    controls: PhsPortSidebandControls,
}
impl Drop for PhsImpedanceSidebandSource<'_, '_, '_> {
    fn drop(&mut self) {}
}

impl<'loaded, 'model, 'scope> PhsImpedanceSidebandSource<'loaded, 'model, 'scope> {
    /// Bind the exact prepared base load, input and original state ordering.
    /// Output order must cover the retained flow order; its envelope tolerances
    /// also limit auxiliary-flow disagreement in original physical port units.
    /// Only base state and retained flow coefficients define the operating orbit;
    /// stored effort/evidence fields are neither trusted nor used as a substitute
    /// for the independent complete source-balance check during preparation.
    pub fn new(
        loaded: &'loaded mut PreparedPhsImpedance<'model, 'scope>,
        orbit: &'loaded PhsImpedanceOrbit,
        table: &'loaded SidebandPortImpedance<'model>,
        controls: &PhsPortSidebandControls,
    ) -> Result<Self> {
        poll(loaded.cx)?;
        let system = loaded.original.system;
        physical::admit_source(system, &orbit.state, controls)?;
        let (n, p, d, width) = (
            system.state_dim(),
            system.port_dim(),
            loaded.dim(),
            orbit.state.coeffs[0].len(),
        );
        if !core::ptr::eq(loaded.load, table.base)
            || orbit.state.omega != loaded.load.omega
            || width > loaded.load.matrices.len()
            || table.harmonics != controls.sideband.harmonics
            || controls.output_harmonics < table.harmonics
            || controls.sideband.samples_per_period as u128
                <= 2 * loaded.input.highest_nonzero_harmonic() as u128
            || orbit.flow_coefficients.len() != p
            || orbit.flow_coefficients.iter().any(|row| {
                row.len() != width || row[0].im != 0.0 || row.iter().any(|&z| !finite(z))
            })
        {
            return Err(invalid(
                "same complete loaded base, resolved input and flow/output orders required",
            ));
        }
        let bands = 2 * table.harmonics + 1;
        if d as u128 * bands as u128 > controls.sideband.maximum_unknowns.min(1024) as u128 {
            return Err(PortImpedanceError::Budget(
                "all original states and physical flows must fit the sideband cap",
            ));
        }
        for &carrier in &table.carriers {
            physical::admit_carrier(&orbit.state, carrier, controls)?;
        }
        table.resources(n, width)?;
        // Reserve four boundary spectra, full augmented vector capacities after
        // state compaction, both mismatch vectors and owned output tolerances.
        // Original-state/input/output arrays keep their existing cap.
        let extra = 96 * p as u128 * bands as u128 + 24 * p as u128;
        if extra > controls.max_workspace_bytes as u128 {
            return Err(PortImpedanceError::Budget(
                "loaded boundary spectra exceed the adapter payload cap",
            ));
        }
        let mut controls = controls.clone();
        controls.max_workspace_bytes -= extra as usize;
        physical::admit_storage_with_retained(n, p, &controls, controls.initial_samples, 0)?;
        let mut operators = Vec::with_capacity(table.samples.len());
        let mut derivatives = Vec::with_capacity(table.samples.len());
        for sample in &table.samples {
            operators.push(physical_operator(
                system,
                sample.angular_frequency,
                &sample.impedance,
                false,
                loaded.cx,
            )?);
            derivatives.push(
                sample
                    .angular_derivative
                    .as_ref()
                    .map(|slope| {
                        physical_operator(system, sample.angular_frequency, slope, true, loaded.cx)
                    })
                    .transpose()?,
            );
        }
        let mut coefficients = orbit.state.coeffs.clone();
        coefficients.extend(orbit.flow_coefficients.iter().cloned());
        let augmented = HbOrbit {
            omega: orbit.state.omega,
            coeffs: coefficients,
            residual_trace: Vec::new(),
            residual: orbit.state.residual,
            aft_samples: orbit.state.aft_samples,
            sampling_defect: orbit.state.sampling_defect,
        };
        Ok(Self {
            problem: LoadedProblem {
                loaded,
                table,
                operators,
                derivatives,
            },
            orbit,
            augmented,
            controls,
        })
    }

    /// Reuse the existing independent base/tangent sampling owner. Scales cover
    /// x followed by physical flow; added equation scales have flow units.
    pub fn prepare(
        &mut self,
        scales: &SidebandScales,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhsImpedanceSidebands<'_, 'loaded, 'model, 'scope>> {
        let loaded = self.problem.loaded;
        loaded.cancelled.set(false);
        loaded.port_refused.set(false);
        loaded.original.source_refused.set(false);
        let result = PreparedPeriodicLinearization::new(
            &self.problem,
            &self.augmented,
            self.controls.sideband,
            scales,
            || cancelled() || loaded.cx.checkpoint().is_err(),
        );
        loaded.after_solve()?;
        Ok(PreparedPhsImpedanceSidebands {
            source: self,
            linearization: result?,
        })
    }

    fn forcing(&self, forcing: &[C64], cancelled: &mut impl FnMut() -> bool) -> Result<Vec<C64>> {
        let system = self.problem.loaded.original.system;
        physical::admit_forcing(system, forcing, &self.controls)?;
        let mut rhs =
            physical::input_forcing(system, forcing, self.controls.sideband.harmonics, cancelled)?;
        let (n, d) = (system.state_dim(), self.problem.dim());
        rhs.resize(d * (2 * self.controls.sideband.harmonics + 1), C64::ZERO);
        for band in (0..2 * self.controls.sideband.harmonics + 1).rev() {
            rhs.copy_within(band * n..(band + 1) * n, band * d);
            rhs[band * d + n..(band + 1) * d].fill(C64::ZERO);
        }
        Ok(rhs)
    }
    fn observation_budget(&self, pair: bool, forcing_slope: bool) -> Result<ObservationBudget> {
        let (n, p, d) = (
            self.problem.loaded.original.system.state_dim(),
            self.problem.table.base.ports,
            self.problem.dim(),
        );
        let c = &self.controls;
        let (h, envelopes) = (c.sideband.harmonics, if pair { 2_u128 } else { 1 });
        let bands = (2 * h + 1) as u128;
        let work = (1 + u128::from(forcing_slope)) * physical::input_work(n, p, h)
            + (if pair { 3 } else { 1 }) * 64 * p as u128 * p as u128 * bands
            + 16 * envelopes * d as u128 * bands;
        if 14 * envelopes * c.initial_samples as u128 > c.max_hessian_actions as u128
            || work + envelopes * physical::minimum_observation_work(n, p, c.initial_samples)
                > c.max_observation_work as u128
        {
            return Err(PortImpedanceError::Budget(
                "loaded response/derivative exceed their cumulative physical observation caps",
            ));
        }
        Ok(ObservationBudget {
            actions: 0,
            work,
            retained_bytes: 0,
        })
    }
    fn split(&self, mut response: SidebandResponse) -> (SidebandResponse, Vec<C64>) {
        let n = self.problem.loaded.original.system.state_dim();
        let (p, d) = (self.problem.table.base.ports, self.problem.dim());
        let mut flow = Vec::with_capacity(p * (2 * self.controls.sideband.harmonics + 1));
        let mut write = 0;
        for band in 0..2 * self.controls.sideband.harmonics + 1 {
            flow.extend_from_slice(&response.coefficients[band * d + n..(band + 1) * d]);
            response
                .coefficients
                .copy_within(band * d..band * d + n, write);
            write += n;
        }
        response.coefficients.truncate(write);
        (response, flow)
    }
    fn boundary(
        &self,
        carrier: f64,
        flow: &[C64],
        output: &[C64],
        primal: Option<&[C64]>,
        tolerances: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<C64>, Vec<f64>)> {
        let table = self.problem.table;
        let (p, h, output_h) = (
            table.base.ports,
            table.harmonics,
            self.controls.output_harmonics,
        );
        let mut effort = vec![C64::ZERO; flow.len()];
        let mut mismatch = vec![0.0; p];
        for band in 0..2 * h + 1 {
            physical::poll(cancelled)?;
            let k = band as isize - h as isize;
            let w = (k as f64).mul_add(table.base.omega, carrier);
            let z = table
                .matrix(w.abs())
                .ok_or_else(|| invalid("missing retained physical impedance"))?;
            let slope = if primal.is_some() {
                Some(
                    table
                        .derivative(w.abs())
                        .ok_or_else(|| invalid("complete authored load slopes required"))?,
                )
            } else {
                None
            };
            for a in 0..p {
                let mut sum = ComplexSum::default();
                for b in 0..p {
                    let value = signed_frequency_entry(z[a * p + b].conj(), w, false);
                    sum.product(value, flow[band * p + b]);
                    if let (Some(base), Some(slope)) = (primal, slope) {
                        let derivative = signed_frequency_entry(slope[a * p + b].conj(), w, true);
                        sum.product(derivative, base[band * p + b]);
                    }
                }
                effort[band * p + a] = sum.value()?;
                let observed = output[(band + output_h - h) * p + a];
                mismatch[a] = det::hypot(mismatch[a], (flow[band * p + a] - observed).abs());
            }
        }
        if mismatch
            .iter()
            .zip(tolerances)
            .any(|(&norm, &limit)| !norm.is_finite() || norm > limit)
        {
            return Err(PortImpedanceError::FlowMismatch {
                derivative: primal.is_some(),
                norms: mismatch,
            });
        }
        Ok((effort, mismatch))
    }
}

/// One complete loaded periodic tangent, reusable across the declared carriers.
pub struct PreparedPhsImpedanceSidebands<'source, 'loaded, 'model, 'scope> {
    source: &'source PhsImpedanceSidebandSource<'loaded, 'model, 'scope>,
    linearization: PreparedPeriodicLinearization<'source, LoadedProblem<'loaded, 'model, 'scope>>,
}
impl Drop for PreparedPhsImpedanceSidebands<'_, '_, '_, '_> {
    fn drop(&mut self) {}
}
impl<'loaded, 'model, 'scope> PreparedPhsImpedanceSidebands<'_, 'loaded, 'model, 'scope> {
    /// Independently reassessed full augmented base residual on both source grids.
    pub fn base_residual(&self) -> f64 {
        self.linearization.base_residual
    }
    /// Difference of independently sampled complete normalized native tangents.
    pub fn tangent_sampling_defect(&self) -> f64 {
        self.linearization.tangent_sampling_defect
    }
    /// Factor a declared signed carrier without reevaluating the base/tangent.
    pub fn factor(
        &mut self,
        carrier: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FactoredPhsImpedanceSidebands<'_, 'loaded, 'model, 'scope>> {
        if !self.source.problem.table.carriers.contains(&carrier) {
            return Err(invalid(
                "carrier was not included in the admitted physical probe band",
            ));
        }
        let result = self.linearization.factor(carrier, || {
            cancelled() || self.source.problem.loaded.cx.checkpoint().is_err()
        });
        self.source.problem.loaded.after_solve()?;
        let factor = result?;
        Ok(FactoredPhsImpedanceSidebands {
            source: self.source,
            factor,
            base_residual: self.base_residual(),
            tangent_sampling_defect: self.tangent_sampling_defect(),
        })
    }
}

/// One full loaded sideband operator/LU, reusable for original physical drives.
pub struct FactoredPhsImpedanceSidebands<'prepared, 'loaded, 'model, 'scope> {
    source: &'prepared PhsImpedanceSidebandSource<'loaded, 'model, 'scope>,
    factor: PreparedSidebandResponse<'prepared>,
    base_residual: f64,
    tangent_sampling_defect: f64,
}
impl Drop for FactoredPhsImpedanceSidebands<'_, '_, '_, '_> {
    fn drop(&mut self) {}
}
impl FactoredPhsImpedanceSidebands<'_, '_, '_, '_> {
    /// Signed original probe carrier in radians per unit time.
    pub fn carrier_angular_frequency(&self) -> f64 {
        self.factor.carrier_angular_frequency()
    }
    /// Solve an infinitesimal effort perturbation. Both signed sidebands remain
    /// explicit; no factor-of-two conversion or conjugate folding is performed.
    pub fn solve(
        &mut self,
        forcing: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsImpedanceSidebandResponse> {
        let source = self.source;
        let loaded = source.problem.loaded;
        let mut cancelled = || cancelled() || loaded.cx.checkpoint().is_err();
        physical::poll(&mut cancelled)?;
        let mut budget = source.observation_budget(false, false)?;
        let state = {
            let rhs = source.forcing(forcing, &mut cancelled)?;
            self.factor.solve(&rhs, &mut cancelled)?
        };
        let (state_response, flow_coefficients) = source.split(state);
        let system = loaded.original.system;
        let c = &source.controls;
        let (port_coefficients, observation_trace, hessian_actions) =
            physical::observe_with_budget(
                &loaded.original,
                &source.orbit.state,
                &state_response.coefficients,
                c,
                &c.envelope_tolerances,
                &mut budget,
                &mut cancelled,
            )?;
        let carrier = self.carrier_angular_frequency();
        let (effort_coefficients, retained_flow_mismatch) = source.boundary(
            carrier,
            &flow_coefficients,
            &port_coefficients,
            None,
            &c.envelope_tolerances,
            &mut cancelled,
        )?;
        physical::poll(&mut cancelled)?;
        Ok(PhsImpedanceSidebandResponse {
            response: PhsPortSidebandResponse {
                omega_rad_s: source.orbit.state.omega,
                carrier_angular_frequency: carrier,
                state_harmonics: c.sideband.harmonics,
                output_harmonics: c.output_harmonics,
                state_dimension: system.state_dim(),
                port_count: system.port_dim(),
                state_response,
                port_coefficients,
                base_residual: self.base_residual,
                tangent_sampling_defect: self.tangent_sampling_defect,
                observation_trace,
                hessian_actions,
            },
            flow_coefficients,
            effort_coefficients,
            retained_flow_mismatch,
        })
    }

    /// Differentiate at the FIXED nonlinear base, reusing this LU. Every shifted
    /// frequency needs an explicitly authored load slope, including real DC.
    /// For w<0, Z_+(w)=Z_-(|w|) and Z_+'(w)=-Z_-'(|w|); for w>=0 both are
    /// conjugated. Purely imaginary DC slopes make these two limits agree.
    /// No finite difference, finite-amplitude perturbation or new factor is used.
    pub fn solve_with_carrier_derivative(
        &mut self,
        forcing: &[C64],
        forcing_angular_derivative: Option<&[C64]>,
        controls: SidebandCarrierDerivativeControls,
        derivative_envelope_tolerances: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsImpedanceSidebandCarrierDerivativeResponse> {
        let source = self.source;
        let loaded = source.problem.loaded;
        let system = loaded.original.system;
        let c = &source.controls;
        let (p, h) = (system.port_dim(), c.sideband.harmonics);
        let carrier = self.carrier_angular_frequency();
        let mut cancelled = || cancelled() || loaded.cx.checkpoint().is_err();
        physical::poll(&mut cancelled)?;
        physical::admit_forcing(system, forcing, c)?;
        if let Some(slope) = forcing_angular_derivative {
            physical::admit_forcing(system, slope, c)?;
        }
        physical::admit_port_tolerances(p, derivative_envelope_tolerances)?;
        // Complete derivative coverage is checked before any new input/source work.
        for k in -(h as isize)..=h as isize {
            if source
                .problem
                .table
                .derivative((k as f64).mul_add(source.orbit.state.omega, carrier).abs())
                .is_none()
            {
                return Err(invalid(
                    "every retained signed frequency requires an authored load derivative",
                ));
            }
        }
        let mut budget = source.observation_budget(true, forcing_angular_derivative.is_some())?;
        let pair = {
            let result = self.factor.prepare_carrier_derivative(
                controls,
                |frequency, out| source.problem.derivative(frequency, out),
                &mut cancelled,
            );
            loaded.after_solve()?;
            let workspace = result?;
            let rhs = source.forcing(forcing, &mut cancelled)?;
            let rhs_derivative = match forcing_angular_derivative {
                Some(slope) => source.forcing(slope, &mut cancelled)?,
                None => vec![C64::ZERO; source.problem.dim() * (2 * h + 1)],
            };
            workspace.solve(&rhs, &rhs_derivative, &mut cancelled)?
        };
        let (response, flow_coefficients) = source.split(pair.response);
        let (angular_derivative, flow_angular_derivative) = source.split(pair.angular_derivative);
        let response = physical::observe_original_pair(
            &loaded.original,
            &source.orbit.state,
            carrier,
            SidebandCarrierDerivativeResponse {
                response,
                angular_derivative,
            },
            c,
            derivative_envelope_tolerances,
            self.base_residual,
            self.tangent_sampling_defect,
            &mut budget,
            &mut cancelled,
        )?;
        let (effort_coefficients, retained_flow_mismatch) = source.boundary(
            carrier,
            &flow_coefficients,
            &response.response.port_coefficients,
            None,
            &c.envelope_tolerances,
            &mut cancelled,
        )?;
        let (effort_angular_derivative, derivative_flow_mismatch) = source.boundary(
            carrier,
            &flow_angular_derivative,
            &response.port_angular_derivative,
            Some(&flow_coefficients),
            derivative_envelope_tolerances,
            &mut cancelled,
        )?;
        physical::poll(&mut cancelled)?;
        Ok(PhsImpedanceSidebandCarrierDerivativeResponse {
            response,
            flow_coefficients,
            effort_coefficients,
            flow_angular_derivative,
            effort_angular_derivative,
            retained_flow_mismatch,
            derivative_flow_mismatch,
        })
    }
}
