//! Independent bore and exterior-mouth flow coordinates in one joint load.
//! The mouth flow is supplied by its physical owner, not inferred from the bore.

use super::*;
use crate::dg_port::terminal::IndependentMouthsRelation;

/// Mean efforts of the original coupled bore/exterior owners.
/// Coordinates are [flow INTO the bore, flow OUT into the exterior mouth],
/// both in m³/s. Efforts are [bore inlet pressure, exterior mouth pressure].
/// Their work pairing is the sum of the two pressure-volume-flow products.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgIndependentMouthResponse {
    pub history_pressure_pa: [f64; 2],
    /// Row-major 2-by-2 effort response [Pa s/m³].
    pub impedance_pa_s_m3: [[f64; 2]; 2],
}

/// Complete mean effort law after the original terminal/branch elimination.
/// Coordinates are bore inflow followed by each declared mouth's outward flow.
#[derive(Clone, Debug, PartialEq)]
pub struct DgIndependentMouthsResponse {
    pub history_pressure_pa: Vec<f64>,
    /// Flat row-major square matrix with width history_pressure_pa.len().
    pub impedance_pa_s_m3: Vec<f64>,
}

impl DgIndependentMouthsResponse {
    pub(in crate::dg_port) fn singleton(&self)
        -> Result<DgIndependentMouthResponse, AcousticRealizeError>
    {
        if self.history_pressure_pa.len() != 2 || self.impedance_pa_s_m3.len() != 4 {
            return Err(invalid("scalar independent mouth response requires exactly one mouth"));
        }
        Ok(DgIndependentMouthResponse {
            history_pressure_pa: [self.history_pressure_pa[0], self.history_pressure_pa[1]],
            impedance_pa_s_m3: [[self.impedance_pa_s_m3[0], self.impedance_pa_s_m3[1]],
                [self.impedance_pa_s_m3[2], self.impedance_pa_s_m3[3]]],
        })
    }
}

/// One physical mouth's signed exchange with the shared radiation history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgIndependentMouthExchange {
    pub outward_flow_m3_s: f64,
    pub mean_pressure_pa: f64,
    pub load_work_j: f64,
}

/// Complete source work, counting the original bore inlet once and every mouth
/// once. Rational work is signed exchange, not a dissipation certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct DgIndependentMouthsFrame {
    pub mouths: Vec<DgIndependentMouthExchange>,
    pub inlet_pressure_pa: f64,
    pub source_work_j: f64,
    pub pressure_work_residual_j: f64,
}

impl DgIndependentMouthsFrame {
    pub(in crate::dg_port) fn singleton(&self) -> Option<DgIndependentMouthFrame> {
        let [mouth] = self.mouths.as_slice() else { return None; };
        Some(DgIndependentMouthFrame {
            outward_flow_m3_s: mouth.outward_flow_m3_s,
            mean_pressure_pa: mouth.mean_pressure_pa,
            inlet_pressure_pa: self.inlet_pressure_pa,
            load_work_j: mouth.load_work_j,
            source_work_j: self.source_work_j,
            pressure_work_residual_j: self.pressure_work_residual_j,
        })
    }

    pub(in crate::dg_port) fn new(
        inlet_flow: f64, mouth_flows: &[f64], mouth_pressures: &[f64],
        inlet_pressure: f64, inlet_work: f64,
        response: &DgIndependentMouthsResponse, dt: f64,
    ) -> Result<Self, AcousticRealizeError> {
        let width = mouth_flows.len() + 1;
        if mouth_flows.is_empty() || mouth_pressures.len() != mouth_flows.len()
            || response.history_pressure_pa.len() != width
            || response.impedance_pa_s_m3.len() != width * width
        {
            return Err(invalid("independent mouth work requires the complete original coordinate basis"));
        }
        // The one-mouth path keeps the exact original arithmetic and checks.
        if mouth_flows.len() == 1 {
            let frame = DgIndependentMouthFrame::new(inlet_flow, mouth_flows[0],
                mouth_pressures[0], inlet_pressure, inlet_work, &response.singleton()?, dt)?;
            return Ok(Self {
                mouths: vec![DgIndependentMouthExchange {
                    outward_flow_m3_s: frame.outward_flow_m3_s,
                    mean_pressure_pa: frame.mean_pressure_pa,
                    load_work_j: frame.load_work_j,
                }],
                inlet_pressure_pa: frame.inlet_pressure_pa,
                source_work_j: frame.source_work_j,
                pressure_work_residual_j: frame.pressure_work_residual_j,
            });
        }
        let flows: Vec<_> = core::iter::once(inlet_flow).chain(mouth_flows.iter().copied()).collect();
        let pressures: Vec<_> = core::iter::once(inlet_pressure).chain(mouth_pressures.iter().copied()).collect();
        let mut expected = response.history_pressure_pa.clone();
        let mut scales = vec![0.0; width];
        for row in 0..width {
            scales[row] = pressures[row].abs() + expected[row].abs();
            for (column, flow) in flows.iter().enumerate() {
                let term = response.impedance_pa_s_m3[row * width + column] * flow;
                expected[row] += term;
                scales[row] += term.abs();
            }
        }
        let mouths: Vec<_> = mouth_flows.iter().zip(mouth_pressures).map(|(&flow, &pressure)|
            DgIndependentMouthExchange { outward_flow_m3_s: flow, mean_pressure_pa: pressure,
                load_work_j: pressure * (flow * dt) }).collect();
        let source_work_j = mouths.iter().fold(inlet_work, |work, mouth| work + mouth.load_work_j);
        let expected_work = flows.iter().zip(&expected).map(|(q, p)| p * (q * dt)).sum::<f64>();
        let pressure_work_residual_j = source_work_j - expected_work;
        let work_scale = inlet_work.abs() + mouths.iter().map(|mouth| mouth.load_work_j.abs()).sum::<f64>()
            + flows.iter().zip(&expected).map(|(q, p)| (p * (q * dt)).abs()).sum::<f64>()
            + flows.iter().zip(&scales).map(|(q, scale)| scale * q.abs() * dt).sum::<f64>();
        if !flows.iter().chain(&pressures).chain(&expected).chain(&scales).all(|value| value.is_finite())
            || ![source_work_j, expected_work, pressure_work_residual_j, work_scale]
                .into_iter().all(f64::is_finite)
            || pressures.iter().zip(&expected).zip(&scales)
                .any(|((actual, expected), scale)| (actual - expected).abs() > 1e-10 * scale.max(1e-300))
            || pressure_work_residual_j.abs() > 1e-10 * work_scale.max(1e-300)
        {
            return Err(invalid("joint DG/mouth mean efforts disagree with original complete pressure work"));
        }
        Ok(Self { mouths, inlet_pressure_pa: inlet_pressure, source_work_j, pressure_work_residual_j })
    }
}

/// Accepted independent mouth exchange from the original joint radiation owner.
/// The rational work is signed exchange, not a dissipation/storage certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DgIndependentMouthFrame {
    pub outward_flow_m3_s: f64,
    pub mean_pressure_pa: f64,
    pub inlet_pressure_pa: f64,
    pub load_work_j: f64,
    /// Original bore inlet work plus original exterior-mouth load work.
    pub source_work_j: f64,
    pub pressure_work_residual_j: f64,
}

impl DgIndependentMouthFrame {
    pub(in crate::dg_port) fn new(
        inlet_flow: f64, mouth_flow: f64, mouth_pressure: f64,
        inlet_pressure: f64, inlet_work: f64,
        response: &DgIndependentMouthResponse, dt: f64,
    ) -> Result<Self, AcousticRealizeError> {
        let flows = [inlet_flow, mouth_flow];
        let pressures = [inlet_pressure, mouth_pressure];
        let mut expected = response.history_pressure_pa;
        let mut scales = [0.0; 2];
        for row in 0..2 {
            scales[row] = pressures[row].abs() + expected[row].abs();
            for (column, flow) in flows.iter().enumerate() {
                let term = response.impedance_pa_s_m3[row][column] * flow;
                expected[row] += term;
                scales[row] += term.abs();
            }
        }
        let load_work_j = mouth_pressure * (mouth_flow * dt);
        let source_work_j = inlet_work + load_work_j;
        let expected_work = expected[0] * (inlet_flow * dt)
            + expected[1] * (mouth_flow * dt);
        let pressure_work_residual_j = source_work_j - expected_work;
        let work_scale = inlet_work.abs() + load_work_j.abs()
            + (expected[0] * (inlet_flow * dt)).abs()
            + (expected[1] * (mouth_flow * dt)).abs()
            + scales[0] * inlet_flow.abs() * dt
            + scales[1] * mouth_flow.abs() * dt;
        if !flows.into_iter().chain(pressures).chain(expected).chain(scales)
            .chain([load_work_j, source_work_j, expected_work,
                pressure_work_residual_j, work_scale]).all(f64::is_finite)
            || (pressures[0] - expected[0]).abs() > 1e-10 * scales[0].max(1e-300)
            || (pressures[1] - expected[1]).abs() > 1e-10 * scales[1].max(1e-300)
            || pressure_work_residual_j.abs() > 1e-10 * work_scale.max(1e-300)
        {
            return Err(invalid("independent DG/mouth mean efforts disagree with original pressure work"));
        }
        Ok(Self { outward_flow_m3_s: mouth_flow, mean_pressure_pa: mouth_pressure,
            inlet_pressure_pa: inlet_pressure, load_work_j, source_work_j,
            pressure_work_residual_j })
    }
}

impl DgPortLoad {
    /// Attach the computed rectangular mouth and retain all self/mutual paths
    /// in the SAME original rational load as the existing DG outlets.
    ///
    /// This is a cold construction. The supplied rectangle is an explicit
    /// uniform mouth mode in the shared rigid baffle. A resolved fluid adapter
    /// must validate its own geometry, projection and pressure-work error.
    /// This prepares independent-flow use; the legacy scalar API retains its
    /// documented inlet-minus-mouth effort and opposite-flow convention.
    pub fn with_independent_mouth(
        self, mouth: DrivenMouthSurface, mechanics_rate_hz: u32, gate: &CancelGate,
    ) -> Result<Self, MimoRadiationInputError> {
        self.with_driven_mouth(mouth, mechanics_rate_hz, gate)
    }

    /// Attach every declared rectangle to the SAME computed radiation owner.
    /// Their mutual terms and every circular-outlet cross term remain live.
    pub fn with_independent_mouths(
        self, mouths: &[DrivenMouthSurface], mechanics_rate_hz: u32, gate: &CancelGate,
    ) -> Result<Self, MimoRadiationInputError> {
        self.with_driven_mouths(mouths, mechanics_rate_hz, gate)
    }

    /// Two independent mean efforts from the current original bore, branch,
    /// channel and radiation histories. This query advances no owner.
    pub fn independent_mouth_response(&self)
        -> Result<DgIndependentMouthResponse, AcousticRealizeError>
    {
        self.independent_mouths_relation()?.efforts.singleton()
    }

    pub fn independent_mouths_response(&self)
        -> Result<DgIndependentMouthsResponse, AcousticRealizeError>
    {
        Ok(self.independent_mouths_relation()?.efforts)
    }

    pub(in crate::dg_port) fn independent_mouths_relation(&self)
        -> Result<IndependentMouthsRelation, AcousticRealizeError>
    {
        let coupling = self.terminal_coupling.as_ref().filter(|value| value.driven_mouth)
            .ok_or_else(|| invalid("independent mouth flow requires its original joint exterior"))?;
        let terminal = self.terminal.as_ref()
            .ok_or_else(|| invalid("independent mouth flow lost its original rational owner"))?;
        coupling.independent_mouths_relation(self.multiport_response()?, terminal,
            &self.terminal_channels)
    }

    /// Advance the two declared physical flows together, retaining original
    /// bore, branch, channel and rational states on refusal or cancellation.
    pub fn step_with_mouth_flow(
        &mut self, bore_inflow_m3_s: f64, mouth_outflow_m3_s: f64, gate: &CancelGate,
    ) -> Result<crate::dg_port::DgStepOutput, AcousticRealizeError> {
        if gate.is_requested() { return Err(invalid("independent DG/mouth step cancelled")); }
        let trial = self.preview_step_with_mouth_flow(bore_inflow_m3_s, mouth_outflow_m3_s)?;
        if gate.is_requested() { return Err(invalid("independent DG/mouth step cancelled")); }
        Ok(trial.accept())
    }

    pub const fn last_independent_mouth_frame(&self) -> Option<&DgIndependentMouthFrame> {
        self.last_independent_mouth_frame.as_ref()
    }

    /// Advance the complete supplied flow vector on one original physical clock.
    /// Wrong arity, nonfinite input, owner refusal and cancellation publish none
    /// of the duct, branch, channel or shared radiation histories.
    pub fn step_with_mouth_flows(
        &mut self, bore_inflow_m3_s: f64, mouth_outflows_m3_s: &[f64], gate: &CancelGate,
    ) -> Result<crate::dg_port::DgStepOutput, AcousticRealizeError> {
        if gate.is_requested() { return Err(invalid("joint DG/mouth step cancelled")); }
        let trial = self.preview_step_with_mouth_flows(bore_inflow_m3_s, mouth_outflows_m3_s)?;
        if gate.is_requested() { return Err(invalid("joint DG/mouth step cancelled")); }
        Ok(trial.accept())
    }

    pub const fn last_independent_mouths_frame(&self) -> Option<&DgIndependentMouthsFrame> {
        self.last_independent_mouths_frame.as_ref()
    }
}
