//! Passive series R-L-C impedance at a characteristic-wave port.
//!
//! This is the bilinear (implicit-midpoint) realization of
//! `p = R U + L dU/dt + p_c`, `C dp_c/dt = U`, not a new integrator.
//! Eliminate the wave port `p = 2 a - Z U` in the SAME solve. The two
//! retained physical states give H = (L U^2 + C p_c^2)/2; only R dissipates.
//! See J. O. Smith, Physical Audio Signal Processing, "Digitizing Elementary
//! Reflectances by Bilinear Transform". A passive circuit is not evidence that
//! these caller-declared coefficients describe a particular physical outlet.

use crate::waveguide::WaveguideError;

/// Positive-real series impedance `R + s L + 1/(s C)` in acoustic SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeriesImpedanceSpec {
    /// Resistance [Pa s/m^3], nonnegative. Zero is lossless, not a default fit.
    pub resistance_pa_s_m3: f64,
    /// Inertance [Pa s^2/m^3], nonnegative. Zero omits inertial storage.
    pub inertance_pa_s2_m3: f64,
    /// Compliance [m^3/Pa], strictly positive when present. None omits the
    /// spring term; it does NOT mean a rigid or closed termination.
    pub compliance_m3_pa: Option<f64>,
}

impl SeriesImpedanceSpec {
    /// Check physical coefficients without choosing a timestep or load.
    ///
    /// # Errors
    /// Negative/nonfinite resistance or inertance, invalid compliance.
    pub fn validate(&self) -> Result<(), WaveguideError> {
        if !self.resistance_pa_s_m3.is_finite() || self.resistance_pa_s_m3 < 0.0
            || !self.inertance_pa_s2_m3.is_finite() || self.inertance_pa_s2_m3 < 0.0
            || self.compliance_m3_pa.is_some_and(|c| !c.is_finite() || c <= 0.0)
        {
            return Err(WaveguideError("series load requires passive finite R/L and positive optional C"));
        }
        Ok(())
    }
}

/// End-of-step physical storage coordinates, not midpoint port observations.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ImpedanceState {
    /// Inertive branch flow [m^3/s]; zero when the inertance is absent.
    pub inertive_flow_m3_s: f64,
    /// Compliance pressure [Pa]; zero when the compliance is absent.
    pub compliance_pressure_pa: f64,
}

/// Immutable candidate. A preview never advances storage or the load clock.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ImpedanceFrame {
    /// End state, accepted only after every coupled participant succeeds.
    pub state: ImpedanceState,
    /// Reflected characteristic pressure wave [Pa].
    pub reflected_pressure_pa: f64,
    /// Midpoint pressure across the complete series load [Pa].
    pub pressure_pa: f64,
    /// Midpoint volume flow INTO the load [m^3/s].
    pub flow_m3_s: f64,
    /// Actual inertance/compliance storage after the step [J].
    pub stored_energy_j: f64,
    /// Difference in actual storage [J].
    pub storage_change_j: f64,
    /// R U_mid^2 dt [J], never total terminal work or a fitted loss.
    pub dissipated_energy_j: f64,
    /// Midpoint pressure times flow times dt [J]; may be negative on release.
    pub supplied_work_j: f64,
}

impl ImpedanceFrame {
    /// Uncorrected local energy-balance residual [J].
    #[must_use]
    pub fn balance_residual_j(&self) -> f64 {
        self.storage_change_j + self.dissipated_energy_j - self.supplied_work_j
    }
}

/// Constant-coefficient, zero-initialized impedance; no per-step allocation.
/// Bilinear frequency warping is explicit: s = (2/dt) (z-1)/(z+1).
#[derive(Debug, Clone, Copy)]
pub struct SeriesImpedance {
    spec: SeriesImpedanceSpec,
    state: ImpedanceState,
    port_impedance: f64,
    dt: f64,
    inverse_sum: f64,
    inertia_weight: f64,
    compliance_coefficient: f64,
    sqrt_inertance: f64,
    sqrt_compliance: f64,
}

impl SeriesImpedance {
    /// Bind an explicit load to a wave impedance and timestep at zero energy.
    ///
    /// # Errors
    /// Invalid coefficients/port/time or unrepresentable discrete coefficients.
    pub fn new(spec: SeriesImpedanceSpec, port_impedance: f64, dt: f64)
        -> Result<Self, WaveguideError>
    {
        spec.validate()?;
        if !port_impedance.is_finite() || port_impedance <= 0.0 || !dt.is_finite() || dt <= 0.0 {
            return Err(WaveguideError("series load requires positive finite port impedance and timestep"));
        }
        let zm = 2.0 * (spec.inertance_pa_s2_m3 / dt);
        let zc = spec.compliance_m3_pa.map_or(0.0, |c| 0.5 * (dt / c));
        // Scale the positive coefficient sum, avoiding overflow of otherwise
        // representable admittance weights. No negative impedance is repaired.
        let scale = port_impedance.max(spec.resistance_pa_s_m3).max(zm).max(zc);
        let sum = port_impedance / scale + spec.resistance_pa_s_m3 / scale + zm / scale + zc / scale;
        let inverse_sum = (1.0 / scale) / sum;
        let inertia_weight = (zm / scale) / sum;
        if ![zm, zc, scale, sum, inverse_sum, inertia_weight].iter().all(|x| x.is_finite())
            || inverse_sum <= 0.0
            || (spec.inertance_pa_s2_m3 > 0.0 && (zm == 0.0 || inertia_weight == 0.0))
            || (spec.compliance_m3_pa.is_some() && zc == 0.0)
        {
            return Err(WaveguideError("series load discrete coefficients are not representable"));
        }
        Ok(Self {
            spec, state: ImpedanceState::default(), port_impedance, dt,
            inverse_sum, inertia_weight, compliance_coefficient: zc,
            sqrt_inertance: fs_math::det::sqrt(spec.inertance_pa_s2_m3),
            sqrt_compliance: fs_math::det::sqrt(spec.compliance_m3_pa.unwrap_or(0.0)),
        })
    }

    /// Immutable coefficient declaration; changing it would require material work.
    #[must_use]
    pub const fn spec(&self) -> &SeriesImpedanceSpec { &self.spec }

    /// Accepted storage coordinates, unchanged by previews and refusals.
    #[must_use]
    pub const fn state(&self) -> ImpedanceState { self.state }

    fn energy_at(&self, state: ImpedanceState) -> f64 {
        let u = self.sqrt_inertance * state.inertive_flow_m3_s;
        let p = self.sqrt_compliance * state.compliance_pressure_pa;
        0.5 * u * u + 0.5 * p * p
    }

    /// Independently evaluated quadratic storage [J], not accumulated port work.
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 { self.energy_at(self.state) }

    /// Solve the reactive boundary and characteristic port simultaneously.
    ///
    /// # Errors
    /// A nonfinite incident wave, state, reflected wave, energy or work.
    pub fn preview_step(&self, incident: f64) -> Result<ImpedanceFrame, WaveguideError> {
        let old = self.state;
        // (Z + R + 2L/dt + dt/(2C)) U_mid = 2a + 2L/dt U_0 - p_c0.
        let flow = 2.0 * (incident * self.inverse_sum)
            + self.inertia_weight * old.inertive_flow_m3_s
            - old.compliance_pressure_pa * self.inverse_sum;
        let reflected = incident - self.port_impedance * flow;
        let pressure = incident + reflected;
        let state = ImpedanceState {
            inertive_flow_m3_s: if self.spec.inertance_pa_s2_m3 > 0.0 {
                2.0 * flow - old.inertive_flow_m3_s
            } else { 0.0 },
            compliance_pressure_pa: if self.spec.compliance_m3_pa.is_some() {
                old.compliance_pressure_pa + 2.0 * (self.compliance_coefficient * flow)
            } else { 0.0 },
        };
        let stored_energy_j = self.energy_at(state);
        let frame = ImpedanceFrame {
            state, reflected_pressure_pa: reflected, pressure_pa: pressure, flow_m3_s: flow,
            stored_energy_j, storage_change_j: stored_energy_j - self.stored_energy_j(),
            dissipated_energy_j: (self.spec.resistance_pa_s_m3 * flow) * (flow * self.dt),
            supplied_work_j: pressure * (flow * self.dt),
        };
        if ![incident, flow, reflected, pressure, state.inertive_flow_m3_s,
            state.compliance_pressure_pa, stored_energy_j, frame.storage_change_j,
            frame.dissipated_energy_j, frame.supplied_work_j, frame.balance_residual_j()]
            .iter().all(|x| x.is_finite())
        {
            return Err(WaveguideError("series load step left the finite set"));
        }
        Ok(frame)
    }

    // Private to this module and its network owner: never accept a caller-forged
    // or stale preview through the public API.
    pub(crate) fn accept_frame(&mut self, frame: ImpedanceFrame) { self.state = frame.state; }

    /// Advance after an immutable preview succeeds. No allocation or partial state.
    ///
    /// # Errors
    /// Same numerical admission as `preview_step`.
    pub fn step(&mut self, incident: f64) -> Result<ImpedanceFrame, WaveguideError> {
        let frame = self.preview_step(incident)?;
        self.accept_frame(frame);
        Ok(frame)
    }
}
