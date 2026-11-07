//! Small-signal sidebands of the complete original host and rational-load ODE.
//!
//! The source's actual A/B/C/D coordinates remain in the periodic tangent and
//! every solve. No impedance fit, sampled replacement, D inverse or hidden pole
//! elimination is used. The host uses its full original storage Hessian.
//!
//! Responses use peak envelopes Re(exp(i*nu*t) sum_k X_k exp(i*k*omega*t)).
//! Signed frequencies are never folded. The exact explicit-ODE carrier
//! derivative is iI on BOTH host and load coordinates, including at zero.
//!
//! The original positive source band qualifies all retained, observed and
//! finite omitted-equation probe bands. Exact zero uses the supplied rational
//! model's DC; this does not establish physical DC or transient qualification.
//! No source storage, heat, stability or infinite-tail certificate is inferred.
//! Original flow/direct-effort and Fourier interfaces remain f64; the shared
//! checks do not promise arbitrary relative accuracy at coupled cancellation.

use super::{PhsMimoError, PreparedPhsMimo};
use crate::phs_orbit::sideband::{
    self as physical, ObservationBudget, PhsPortSidebandCheck, PhsPortSidebandControls,
    PhsPortSidebandError, PhsPortSidebandResponse,
};
use fs_math::{c64::C64, dd::Dd};
use fs_orbit::{
    HbOrbit,
    sideband::{
        PreparedPeriodicLinearization, PreparedSidebandResponse, SidebandError, SidebandResponse,
        SidebandScales,
        sensitivity::SidebandCarrierDerivativeControls,
    },
};

/// Original causal-owner and shared physical-observer refusals remain distinct.
#[derive(Debug)]
pub enum PhsMimoSidebandError {
    Causal(PhsMimoError),
    Physical(PhsPortSidebandError),
}
impl From<PhsMimoError> for PhsMimoSidebandError {
    fn from(value: PhsMimoError) -> Self {
        Self::Causal(value)
    }
}
impl From<PhsPortSidebandError> for PhsMimoSidebandError {
    fn from(value: PhsPortSidebandError) -> Self {
        Self::Physical(value)
    }
}
impl From<SidebandError> for PhsMimoSidebandError {
    fn from(value: SidebandError) -> Self {
        Self::Physical(value.into())
    }
}
type Result<T> = core::result::Result<T, PhsMimoSidebandError>;

fn invalid(message: &'static str) -> PhsMimoSidebandError {
    PhsMimoSidebandError::Physical(PhsPortSidebandError::Invalid(message))
}
fn budget(message: &'static str) -> PhsMimoSidebandError {
    PhsMimoSidebandError::Physical(PhsPortSidebandError::Budget(message))
}

/// Host-only physical chart plus EVERY original retained source coordinate.
/// Together the two state vectors form the original joint response. The
/// state_response residuals in response certify the COMPLETE joint equations,
/// including source rows, rather than a host-only equation or fitted reduction.
#[derive(Clone, Debug)]
pub struct PhsMimoSidebandResponse {
    pub response: PhsPortSidebandResponse,
    pub source_dimension: usize,
    /// Sideband-major original source coordinates, (k+state_harmonics)*q+j.
    /// Real-pole and real/negative-imaginary pair ordering is unchanged.
    pub source_coefficients: Vec<C64>,
    pub boundary_port_count: usize,
    /// S^T Hxx dx, sideband-major over the reported OUTPUT order.
    pub boundary_flow_coefficients: Vec<C64>,
    /// C dz + D dy, over the same output order. No D or B inverse is needed.
    pub boundary_effort_coefficients: Vec<C64>,
    pub boundary_flow_trace: Vec<PhsPortSidebandCheck>,
    /// The independently sampled D dy trace also covers the total effort:
    /// output order covers the retained C dz polynomial exactly.
    pub boundary_effort_trace: Vec<PhsPortSidebandCheck>,
    /// Unchanged positive source qualification band, in radians/time.
    pub source_band_rad_s: [f64; 2],
    /// Whether a retained, reported or finite omitted-equation probe band is
    /// exactly zero. Such a probe uses model DC, not a physical DC certificate.
    pub probe_includes_model_dc: bool,
}

/// The second complete chart is the analytic derivative with respect to signed
/// carrier angular frequency at FIXED base orbit, source and probe envelopes.
/// Its coefficients/tolerances have their original units per radians/time.
/// A supplied forcing derivative is included without differentiating the base.
#[derive(Clone, Debug)]
pub struct PhsMimoSidebandCarrierDerivativeResponse {
    pub response: PhsMimoSidebandResponse,
    pub angular_derivative: PhsMimoSidebandResponse,
    /// Shared original Hessian actions across both physical observations.
    pub total_hessian_actions: usize,
}

/// One immutable causal source/base and its physical observation limits.
/// The exclusive constructor borrow prevents concurrent source/tangent queries;
/// preparation and factorization preserve that discipline through their lifetimes.
pub struct PhsMimoSidebandSource<'loaded, 'model, 'scope> {
    loaded: &'loaded PreparedPhsMimo<'model, 'scope>,
    orbit: &'loaded HbOrbit,
    host_orbit: HbOrbit,
    controls: PhsPortSidebandControls,
    tolerances: Vec<f64>,
}
impl Drop for PhsMimoSidebandSource<'_, '_, '_> {
    fn drop(&mut self) {}
}

impl<'loaded, 'model, 'scope> PhsMimoSidebandSource<'loaded, 'model, 'scope> {
    /// Bind a COMPLETE original joint orbit. Stored residual evidence is not
    /// trusted: prepare independently re-evaluates the actual source balance.
    /// Separate absolute envelope tolerances are required for boundary flow and
    /// effort; each includes retained sampling change AND observed omitted tail.
    /// Host port tolerances and all cumulative caps come from controls.
    pub fn new(
        loaded: &'loaded mut PreparedPhsMimo<'model, 'scope>,
        orbit: &'loaded HbOrbit,
        controls: &PhsPortSidebandControls,
        boundary_flow_tolerances: &[f64],
        boundary_effort_tolerances: &[f64],
    ) -> Result<Self> {
        loaded.reset()?;
        let system = loaded.original.system;
        let (n, d, p, m) = (
            loaded.host_dimension(), loaded.dimension(),
            loaded.boundary.ports, system.port_dim(),
        );
        let width = orbit.coeffs.first().map_or(0, Vec::len);
        if orbit.omega != loaded.input.omega_rad_s()
            || !orbit.omega.is_finite()
            || orbit.omega <= 0.0
            || orbit.coeffs.len() != d
            || width == 0
            || width > loaded.limits.maximum_harmonics + 1
            || controls.sideband.harmonics > 64
            || controls.output_harmonics < controls.sideband.harmonics
            || controls.output_harmonics > loaded.limits.maximum_harmonics
            || controls.sideband.samples_per_period as u128
                <= 2 * loaded.input.highest_nonzero_harmonic() as u128
            || orbit.coeffs.iter().any(|row| {
                row.len() != width || row[0].im != 0.0
                    || row.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            })
        {
            return Err(invalid(
                "complete original joint orbit, fixed input clock and resolved output/input orders required",
            ));
        }
        physical::admit_port_tolerances(m, &controls.envelope_tolerances)?;
        physical::admit_port_tolerances(p, boundary_flow_tolerances)?;
        physical::admit_port_tolerances(p, boundary_effort_tolerances)?;
        let bands = 2 * controls.sideband.harmonics + 1;
        if d as u128 * bands as u128 > controls.sideband.maximum_unknowns.min(1024) as u128 {
            return Err(budget("all original host and source states must fit the sideband cap"));
        }
        let retained = retained_payload(n, d, m + 2 * p, width, controls, 1);
        physical::admit_storage_with_retained(
            n, m + 2 * p, controls, controls.initial_samples, retained,
        )?;
        let host_orbit = HbOrbit {
            omega: orbit.omega,
            coeffs: orbit.coeffs[..n].to_vec(),
            residual_trace: Vec::new(),
            residual: 0.0,
            aft_samples: 0,
            sampling_defect: 0.0,
        };
        physical::admit_source(system, &host_orbit, controls)?;
        if width > 1 {
            loaded.admit_harmonics(width - 1)?;
        }
        let mut tolerances = Vec::with_capacity(m + 2 * p);
        tolerances.extend_from_slice(&controls.envelope_tolerances);
        tolerances.extend_from_slice(boundary_flow_tolerances);
        tolerances.extend_from_slice(boundary_effort_tolerances);
        Ok(Self {
            loaded, orbit, host_orbit, controls: controls.clone(), tolerances,
        })
    }

    /// Cache one independently admitted full periodic tangent. Scales contain
    /// host coordinates followed by EVERY retained source coordinate.
    pub fn prepare(
        &mut self,
        scales: &SidebandScales,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhsMimoSidebands<'_, 'loaded, 'model, 'scope>> {
        self.loaded.reset()?;
        let result = PreparedPeriodicLinearization::new(
            self.loaded, self.orbit, self.controls.sideband, scales,
            || cancelled() || self.loaded.cx.checkpoint().is_err(),
        );
        self.loaded.after()?;
        Ok(PreparedPhsMimoSidebands {
            source: self,
            linearization: result?,
        })
    }

    fn admit_carrier(&self, carrier: f64) -> Result<bool> {
        physical::admit_carrier(&self.host_orbit, carrier, &self.controls)?;
        // These are the same finite omitted-equation guard bands measured by
        // fs-orbit, in addition to every retained and reported physical band.
        let maximum = (2 * self.controls.sideband.harmonics + 2)
            .max(self.controls.output_harmonics);
        let mut dc = false;
        for k in -(maximum as isize)..=maximum as isize {
            self.loaded.cx.checkpoint().map_err(|_| PhsMimoError::Cancelled)?;
            let frequency = (k as f64).mul_add(self.orbit.omega, carrier);
            if !frequency.is_finite() {
                return Err(invalid("finite original signed probe frequencies required"));
            }
            if frequency == 0.0 {
                dc = true;
            } else {
                self.loaded.source.continuous_operator().admit_frequency(frequency.abs())
                    .map_err(PhsMimoError::from)?;
            }
        }
        Ok(dc)
    }

    fn forcing(&self, forcing: &[C64], cancelled: &mut impl FnMut() -> bool) -> Result<Vec<C64>> {
        let system = self.loaded.original.system;
        physical::admit_forcing(system, forcing, &self.controls)?;
        let h = self.controls.sideband.harmonics;
        let mut rhs = physical::input_forcing(system, forcing, h, cancelled)?;
        let (n, d) = (self.loaded.host_dimension(), self.loaded.dimension());
        rhs.resize(d * (2 * h + 1), C64::ZERO);
        for band in (0..2 * h + 1).rev() {
            rhs.copy_within(band * n..(band + 1) * n, band * d);
            rhs[band * d + n..(band + 1) * d].fill(C64::ZERO);
        }
        Ok(rhs)
    }

    fn observation_budget(&self, pair: bool, forcing_slope: bool) -> Result<ObservationBudget> {
        let c = &self.controls;
        let (n, d, p, m) = (
            self.loaded.host_dimension(), self.loaded.dimension(),
            self.loaded.boundary.ports, self.loaded.original.system.port_dim(),
        );
        let outputs = m + 2 * p;
        let envelopes = if pair { 2_u128 } else { 1 };
        let bands = (2 * c.sideband.harmonics + 1) as u128;
        let work = (1 + u128::from(forcing_slope))
            * physical::input_work(n, m, c.sideband.harmonics)
            + envelopes * (64 * p as u128 * self.loaded.source_dimension() as u128 * bands
                + 32 * d as u128 * bands
                + 32 * outputs as u128 * (2 * c.output_harmonics + 1) as u128);
        let minimum = physical::minimum_observation_work(n, outputs, c.initial_samples)
            + 7 * c.initial_samples as u128 * self.extra_work_per_sample();
        if 14 * envelopes * c.initial_samples as u128 > c.max_hessian_actions as u128
            || work + envelopes * minimum > c.max_observation_work as u128
        {
            return Err(budget(
                "joint response and physical observations exceed their cumulative action/work caps",
            ));
        }
        let retained_bytes = retained_payload(
            n, d, outputs, self.host_orbit.coeffs[0].len(), c, envelopes,
        );
        physical::admit_storage_with_retained(
            n, outputs, c, c.initial_samples, retained_bytes,
        )?;
        Ok(ObservationBudget { actions: 0, work, retained_bytes })
    }

    fn extra_work_per_sample(&self) -> u128 {
        // The shared observer charges the G and S maps. D is a separate p-by-p
        // contraction for both real and imaginary directions at every sample.
        128 * (self.loaded.boundary.ports as u128).pow(2)
    }

    fn split(&self, mut state: SidebandResponse) -> (SidebandResponse, Vec<C64>) {
        let (n, d, q) = (
            self.loaded.host_dimension(), self.loaded.dimension(), self.loaded.source_dimension(),
        );
        let bands = 2 * self.controls.sideband.harmonics + 1;
        let mut source = Vec::with_capacity(q * bands);
        for band in 0..bands {
            source.extend_from_slice(&state.coefficients[band * d + n..(band + 1) * d]);
            state.coefficients.copy_within(band * d..band * d + n, band * n);
        }
        state.coefficients.truncate(n * bands);
        (state, source)
    }

    #[allow(clippy::too_many_arguments)]
    fn observe(
        &self,
        state: SidebandResponse,
        carrier: f64,
        dc: bool,
        tolerances: &[f64],
        base_residual: f64,
        tangent_sampling_defect: f64,
        budget: &mut ObservationBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PhsMimoSidebandResponse> {
        let (state_response, source_coefficients) = self.split(state);
        let (system, boundary, c) = (
            self.loaded.original.system, &self.loaded.boundary, &self.controls,
        );
        let (m, p) = (system.port_dim(), boundary.ports);
        let (mapped, trace, hessian_actions) = physical::observe_mapped_with_budget(
            &self.loaded.original,
            &self.host_orbit,
            &state_response.coefficients,
            c,
            tolerances,
            m + 2 * p,
            self.extra_work_per_sample(),
            &|effort, output| {
                let (ports, rest) = output.split_at_mut(m);
                let (flow, direct) = rest.split_at_mut(p);
                system.output_from_effort_into(effort, ports).is_ok()
                    && boundary.flow_into(effort, flow)
                    && boundary.direct_effort_into(flow, direct)
            },
            budget,
            cancelled,
        )?;
        let output_bands = 2 * c.output_harmonics + 1;
        let mut port_coefficients = Vec::with_capacity(m * output_bands);
        let mut boundary_flow_coefficients = Vec::with_capacity(p * output_bands);
        let mut boundary_effort_coefficients = Vec::with_capacity(p * output_bands);
        for (band, values) in mapped.chunks_exact(m + 2 * p).enumerate() {
            physical::poll(cancelled)?;
            port_coefficients.extend_from_slice(&values[..m]);
            boundary_flow_coefficients.extend_from_slice(&values[m..m + p]);
            let k = band as isize - c.output_harmonics as isize;
            for port in 0..p {
                let direct = values[m + p + port];
                let (mut re, mut im) = (Dd::from_f64(direct.re), Dd::from_f64(direct.im));
                if k.unsigned_abs() <= c.sideband.harmonics {
                    let first = (k + c.sideband.harmonics as isize) as usize * boundary.states;
                    for j in 0..boundary.states {
                        let coefficient = boundary.c[port * boundary.states + j];
                        if coefficient != 0.0 {
                            let z = source_coefficients[first + j];
                            re = re + Dd::from_f64(coefficient) * Dd::from_f64(z.re);
                            im = im + Dd::from_f64(coefficient) * Dd::from_f64(z.im);
                        }
                    }
                }
                let effort = C64::new(re.to_f64(), im.to_f64());
                if !effort.re.is_finite() || !effort.im.is_finite() {
                    return Err(PhsPortSidebandError::SourceRefused.into());
                }
                boundary_effort_coefficients.push(effort);
            }
        }
        let (mut port_trace, mut flow_trace, mut effort_trace) =
            (Vec::new(), Vec::new(), Vec::new());
        for check in trace {
            let part = |start: usize, len: usize| PhsPortSidebandCheck {
                samples: check.samples,
                omitted_envelope_norm: check.omitted_envelope_norm[start..start + len].to_vec(),
                sampling_change_norm: check.sampling_change_norm.as_ref()
                    .map(|values| values[start..start + len].to_vec()),
            };
            port_trace.push(part(0, m));
            flow_trace.push(part(m, p));
            effort_trace.push(part(m + p, p));
        }
        physical::poll(cancelled)?;
        Ok(PhsMimoSidebandResponse {
            response: PhsPortSidebandResponse {
                omega_rad_s: self.orbit.omega,
                carrier_angular_frequency: carrier,
                state_harmonics: c.sideband.harmonics,
                output_harmonics: c.output_harmonics,
                state_dimension: self.loaded.host_dimension(),
                port_count: m,
                state_response,
                port_coefficients,
                base_residual,
                tangent_sampling_defect,
                observation_trace: port_trace,
                hessian_actions,
            },
            source_dimension: self.loaded.source_dimension(),
            source_coefficients,
            boundary_port_count: p,
            boundary_flow_coefficients,
            boundary_effort_coefficients,
            boundary_flow_trace: flow_trace,
            boundary_effort_trace: effort_trace,
            source_band_rad_s: self.loaded.source.limits().band_rad_s,
            probe_includes_model_dc: dc,
        })
    }
}

// Conservative live adapter reserve: host-base copy, owned tolerances, mapped
// forcing, complete joint/source state spectra (including compacted Vec capacity),
// both observations' retained outputs, split traces and real/complex map scratch.
// Generic tangent/LU/derivative storage retains its OWN separate caps.
fn retained_payload(
    n: usize, d: usize, outputs: usize, width: usize,
    c: &PhsPortSidebandControls, envelopes: u128,
) -> u128 {
    let (n, d, outputs, width) = (n as u128, d as u128, outputs as u128, width as u128);
    32 * n * width + 64 * n + 64 * outputs
        + envelopes * (
            64 * d * (2 * c.sideband.harmonics + 1) as u128
            + 64 * outputs * (2 * c.output_harmonics + 1) as u128
            + 1024 * outputs + 16_384
        )
}

/// One full original-state periodic tangent, reusable across admitted carriers.
pub struct PreparedPhsMimoSidebands<'source, 'loaded, 'model, 'scope> {
    source: &'source PhsMimoSidebandSource<'loaded, 'model, 'scope>,
    linearization: PreparedPeriodicLinearization<'source, PreparedPhsMimo<'model, 'scope>>,
}
impl Drop for PreparedPhsMimoSidebands<'_, '_, '_, '_> {
    fn drop(&mut self) {}
}
impl<'loaded, 'model, 'scope> PreparedPhsMimoSidebands<'_, 'loaded, 'model, 'scope> {
    pub fn base_residual(&self) -> f64 {
        self.linearization.base_residual
    }
    pub fn tangent_sampling_defect(&self) -> f64 {
        self.linearization.tangent_sampling_defect
    }
    /// Factor one signed carrier while retaining every original source state.
    pub fn factor(
        &mut self,
        carrier: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FactoredPhsMimoSidebands<'_, 'loaded, 'model, 'scope>> {
        self.source.loaded.reset()?;
        physical::poll(&mut cancelled)?;
        let dc = self.source.admit_carrier(carrier)?;
        let base_residual = self.base_residual();
        let tangent_sampling_defect = self.tangent_sampling_defect();
        let result = self.linearization.factor(
            carrier, || cancelled() || self.source.loaded.cx.checkpoint().is_err(),
        );
        self.source.loaded.after()?;
        Ok(FactoredPhsMimoSidebands {
            source: self.source, factor: result?,
            base_residual, tangent_sampling_defect, dc,
        })
    }
}

/// One complete joint factor, reusable for original physical port probes.
pub struct FactoredPhsMimoSidebands<'prepared, 'loaded, 'model, 'scope> {
    source: &'prepared PhsMimoSidebandSource<'loaded, 'model, 'scope>,
    factor: PreparedSidebandResponse<'prepared>,
    base_residual: f64,
    tangent_sampling_defect: f64,
    dc: bool,
}
impl Drop for FactoredPhsMimoSidebands<'_, '_, '_, '_> {
    fn drop(&mut self) {}
}
impl FactoredPhsMimoSidebands<'_, '_, '_, '_> {
    pub fn carrier_angular_frequency(&self) -> f64 {
        self.factor.carrier_angular_frequency()
    }

    /// Peak complex external effort envelope, sideband-major in original G
    /// port order. Source states are driven only through their actual B/S map.
    pub fn solve(
        &mut self,
        forcing: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsMimoSidebandResponse> {
        let source = self.source;
        source.loaded.reset()?;
        let mut cancelled = || cancelled() || source.loaded.cx.checkpoint().is_err();
        physical::poll(&mut cancelled)?;
        let mut budget = source.observation_budget(false, false)?;
        let state = {
            let rhs = source.forcing(forcing, &mut cancelled)?;
            self.factor.solve(&rhs, &mut cancelled)?
        };
        let result = source.observe(
            state, self.carrier_angular_frequency(), self.dc, &source.tolerances,
            self.base_residual, self.tangent_sampling_defect, &mut budget, &mut cancelled,
        );
        source.loaded.after()?;
        result
    }

    /// Solve the response and its analytic signed carrier derivative on the
    /// SAME full joint LU. No frequency differences, replacement impedance or
    /// extra base/source preparation is performed.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_carrier_derivative(
        &mut self,
        forcing: &[C64],
        forcing_angular_derivative: Option<&[C64]>,
        controls: SidebandCarrierDerivativeControls,
        port_derivative_tolerances: &[f64],
        boundary_flow_derivative_tolerances: &[f64],
        boundary_effort_derivative_tolerances: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhsMimoSidebandCarrierDerivativeResponse> {
        let source = self.source;
        source.loaded.reset()?;
        let system = source.loaded.original.system;
        let c = &source.controls;
        let (m, p, d) = (system.port_dim(), source.loaded.boundary.ports, source.loaded.dimension());
        let mut cancelled = || cancelled() || source.loaded.cx.checkpoint().is_err();
        physical::poll(&mut cancelled)?;
        physical::admit_forcing(system, forcing, c)?;
        if let Some(slope) = forcing_angular_derivative {
            physical::admit_forcing(system, slope, c)?;
        }
        physical::admit_port_tolerances(m, port_derivative_tolerances)?;
        physical::admit_port_tolerances(p, boundary_flow_derivative_tolerances)?;
        physical::admit_port_tolerances(p, boundary_effort_derivative_tolerances)?;
        let mut budget = source.observation_budget(true, forcing_angular_derivative.is_some())?;
        let pair = {
            let workspace = self.factor.prepare_carrier_derivative(
                controls,
                |frequency, out| {
                    if !frequency.is_finite() || out.len() != d * d {
                        return false;
                    }
                    out.fill(C64::ZERO);
                    for j in 0..d {
                        out[j * d + j] = C64::new(0.0, 1.0);
                    }
                    true
                },
                &mut cancelled,
            )?;
            let rhs = source.forcing(forcing, &mut cancelled)?;
            let rhs_derivative = match forcing_angular_derivative {
                Some(slope) => source.forcing(slope, &mut cancelled)?,
                None => vec![C64::ZERO; d * (2 * c.sideband.harmonics + 1)],
            };
            workspace.solve(&rhs, &rhs_derivative, &mut cancelled)?
        };
        let carrier = self.carrier_angular_frequency();
        let response = source.observe(
            pair.response, carrier, self.dc, &source.tolerances,
            self.base_residual, self.tangent_sampling_defect, &mut budget, &mut cancelled,
        )?;
        let mut derivative_tolerances = Vec::with_capacity(m + 2 * p);
        derivative_tolerances.extend_from_slice(port_derivative_tolerances);
        derivative_tolerances.extend_from_slice(boundary_flow_derivative_tolerances);
        derivative_tolerances.extend_from_slice(boundary_effort_derivative_tolerances);
        let derivative = source.observe(
            pair.angular_derivative, carrier, self.dc, &derivative_tolerances,
            self.base_residual, self.tangent_sampling_defect, &mut budget, &mut cancelled,
        ).map_err(|error| match error {
            PhsMimoSidebandError::Physical(PhsPortSidebandError::ObservationUnresolved { last }) =>
                PhsMimoSidebandError::Physical(
                    PhsPortSidebandError::DerivativeObservationUnresolved { last },
                ),
            other => other,
        });
        source.loaded.after()?;
        Ok(PhsMimoSidebandCarrierDerivativeResponse {
            response,
            angular_derivative: derivative?,
            total_hessian_actions: budget.actions,
        })
    }
}
