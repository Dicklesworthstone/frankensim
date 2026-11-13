//! The original thermal Carreau pipe on the composition contract.
//!
//! Two hydraulic terminals share one fully developed no-slip radial profile.
//! Their pressure difference drives its signed volume flow. The heat port's
//! temperature re-evaluates the same authored material law, and all original
//! viscous heat leaves that port. A connected thermal storage can therefore
//! change the viscosity without a pairing-specific evolution loop.
//!
//! This is a quasisteady, incompressible resistance: it adds no fluid storage,
//! inertia, entrance loss, transition criterion or thermal transport. Uniform
//! temperature and the authored wall-rate ceiling remain explicit assumptions.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_ports::component::{
    Class, Component, EndOfStep, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};

use super::{
    CarreauPipeControls, CarreauPipeProfile, LaminarCarreauPipe, MaterialError, state_error,
};

/// A source-owned nonlinear hydraulic resistance with thermal feedback.
///
/// Inputs, in port order, are pressures `a`, `b` [Pa] and the uniform material
/// temperature `heat` [K]. Outputs are volume flows `Q`, `-Q` [m³/s] into the
/// pipe and the negative original viscous power [W]. The state is `()`; all
/// deposited heat belongs to the connected thermal owner.
///
/// Proposals and analytic port tangents allocate no heap storage on success.
/// The original bounded constitutive solve and positive radial quadrature
/// determine every response; unsupported temperature, shifted shear support,
/// wall rate, or exhausted work refuses without publishing outputs.
#[derive(Clone, Debug)]
pub struct CarreauPipeComponent {
    pipe: LaminarCarreauPipe,
    controls: CarreauPipeControls,
    ports: [PortSpec; 3],
}

fn poll(cx: &Cx<'_>) -> Result<(), MaterialError> {
    cx.checkpoint()
        .map_err(|_| state_error("Carreau pipe component cancelled"))
}

fn values(values: &[f64]) -> Result<[f64; 3], MaterialError> {
    match values {
        [a, b, c] if values.iter().all(|value| value.is_finite()) => Ok([*a, *b, *c]),
        _ => Err(state_error(
            "Carreau pipe component requires three finite port values",
        )),
    }
}

impl CarreauPipeComponent {
    /// Wrap the original admitted pipe and its finite numerical work budgets.
    ///
    /// Each proposal rebinds its uniform temperature to the `heat` input,
    /// retaining the original law, geometry and physical wall-rate ceiling.
    pub fn new(
        pipe: LaminarCarreauPipe,
        controls: CarreauPipeControls,
    ) -> Result<Self, MaterialError> {
        controls.admit()?;
        let hydraulic = |name| {
            PortSpec::power(
                StableId::new(name).expect("static port name"),
                PortKind::FluidPressureFlux,
                PortRole::EffortIn,
                NonZeroUsize::MIN,
            )
            .expect("scalar hydraulic power port")
        };
        Ok(Self {
            pipe,
            controls,
            ports: [
                hydraulic("a"),
                hydraulic("b"),
                PortSpec::heat(
                    StableId::new("heat").expect("static port name"),
                    PortRole::EffortIn,
                )
                .expect("thermal outlet"),
            ],
        })
    }

    /// The unchanged original law, geometry and authored wall-rate envelope.
    #[must_use]
    pub const fn pipe(&self) -> &LaminarCarreauPipe {
        &self.pipe
    }

    /// Original local-solve and quadrature controls.
    #[must_use]
    pub const fn controls(&self) -> CarreauPipeControls {
        self.controls
    }

    fn trial(
        &self,
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(CarreauPipeProfile, Receipt), MaterialError> {
        poll(cx)?;
        let [a, b, temperature] = values(inputs)?;
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(state_error(
                "Carreau pipe component requires a positive finite timestep",
            ));
        }
        let pressure = a - b;
        // Rebinding is the source's own admission of both temperature and
        // the complete physical rate interval, including the zero-shear core.
        let pipe = self.pipe.at_temperature(temperature)?;
        let profile = pipe.profile_at_pressure_poll(pressure, self.controls, || poll(cx))?;
        let flow = profile.flow();
        let work_a = dt_s * a * flow.volumetric_flow_m3_s;
        let work_b = -dt_s * b * flow.volumetric_flow_m3_s;
        let heat = dt_s * flow.dissipated_power_w;
        // The accepted wall-stress solve may reconstruct a pressure a little
        // different from the supplied pressure. Keep its ORIGINAL heat and
        // expose that measured defect in the audit allowance, never as heat
        // invented by the adapter or as a manufactured numerical loss.
        let source_defect =
            dt_s * flow.pressure_residual_pa.abs() * flow.volumetric_flow_m3_s.abs();
        let rounding = 64.0 * f64::EPSILON;
        let allowance = source_defect
            + rounding * work_a.abs()
            + rounding * work_b.abs()
            + rounding * heat
            + f64::MIN_POSITIVE;
        let entropy = heat / temperature;
        let tolerance =
            (self.controls.relative_constitutive_tolerance + 16.0 * f64::EPSILON) * pressure.abs();
        if ![work_a, work_b, heat, allowance, entropy, tolerance]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(state_error(
                "Carreau pipe component energy account is unrepresentable",
            ));
        }
        poll(cx)?;
        Ok((
            profile,
            Receipt {
                stored_change_j: 0.0,
                dissipated_j: heat,
                numerical_loss_j: 0.0,
                entropy_generated_j_per_k: entropy,
                energy_tolerance_j: allowance,
                residual: flow.pressure_residual_pa.abs(),
                tolerance,
            },
        ))
    }

    fn jacobian(
        &self,
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<[[f64; 3]; 3], MaterialError> {
        let (profile, _) = self.trial(dt_s, inputs, cx)?;
        let conductance = profile.pressure_conductance_m3_s_pa()?;
        let thermal = profile.temperature_flow_derivative_m3_s_k()?;
        let flow = profile.flow();
        let heat_pressure = flow.volumetric_flow_m3_s + flow.pressure_drop_pa * conductance;
        let heat_temperature = flow.pressure_drop_pa * thermal;
        let matrix = [
            [conductance, -conductance, thermal],
            [-conductance, conductance, -thermal],
            [-heat_pressure, heat_pressure, -heat_temperature],
        ];
        if matrix.iter().flatten().any(|value| !value.is_finite()) {
            return Err(state_error(
                "Carreau pipe component tangent is unrepresentable",
            ));
        }
        Ok(matrix)
    }
}

impl Component for CarreauPipeComponent {
    type State = ();
    type Error = MaterialError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, MaterialError> {
        Ok(0.0)
    }

    fn end_of_step(&self, _port: usize) -> EndOfStep {
        // This algebraic, held-input response has no incoming flow history
        // from which the midpoint endpoint relation could be reconstructed.
        EndOfStep::Untold
    }

    fn propose(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, MaterialError> {
        if outputs.len() != 3 {
            return Err(state_error(
                "Carreau pipe component requires three output slots",
            ));
        }
        let (profile, receipt) = self.trial(dt_s, inputs, cx)?;
        let flow = profile.flow();
        let next = [
            flow.volumetric_flow_m3_s,
            -flow.volumetric_flow_m3_s,
            -flow.dissipated_power_w,
        ];
        poll(cx)?;
        *new = ();
        outputs.copy_from_slice(&next);
        Ok(Proposal::Step(receipt))
    }
}

impl PortTangent for CarreauPipeComponent {
    fn jvp(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), MaterialError> {
        if d_outputs.len() != 3 {
            return Err(state_error(
                "Carreau pipe component requires three tangent output slots",
            ));
        }
        let direction = values(d_inputs)?;
        let matrix = self.jacobian(dt_s, inputs, cx)?;
        let mut next = [0.0; 3];
        for (value, row) in next.iter_mut().zip(matrix) {
            *value = row.iter().zip(direction).map(|(j, d)| j * d).sum();
        }
        if next.iter().any(|value| !value.is_finite()) {
            return Err(state_error(
                "Carreau pipe component tangent action is unrepresentable",
            ));
        }
        poll(cx)?;
        d_outputs.copy_from_slice(&next);
        Ok(())
    }

    fn vjp(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), MaterialError> {
        if bar_inputs.len() != 3 {
            return Err(state_error(
                "Carreau pipe component requires three tangent input slots",
            ));
        }
        let direction = values(bar_outputs)?;
        let matrix = self.jacobian(dt_s, inputs, cx)?;
        let mut next = [0.0; 3];
        for (column, value) in next.iter_mut().enumerate() {
            *value = matrix
                .iter()
                .zip(direction)
                .map(|(row, d)| row[column] * d)
                .sum();
        }
        if next.iter().any(|value| !value.is_finite()) {
            return Err(state_error(
                "Carreau pipe component adjoint action is unrepresentable",
            ));
        }
        poll(cx)?;
        bar_inputs.copy_from_slice(&next);
        Ok(())
    }
}
