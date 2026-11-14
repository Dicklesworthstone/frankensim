//! Original temperature-dependent layer conduction on the composition contract.
//!
//! The heat ports receive the temperatures chosen by their thermal owners.
//! The unchanged layer law integrates each conductivity chart and solves its
//! contact jump. It returns one heat transfer and its two analytic endpoint
//! derivatives, without freezing conductivity at the initial state. Thermal
//! storage, time integration and phase evolution remain with the connected
//! owners; these source-free, quasisteady layers have no independent state.

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_ports::StableId;
use fs_ports::component::{
    Class, Component, EndOfStep, PortRole, PortSpec, PortTangent, Proposal, Receipt,
};

use super::{
    ConductivityCurve, DepositionError, LayerConductionConfig, LayerConductionResponse,
    LayerGeometry,
};

/// A two-terminal conductor using the original supported material charts.
///
/// Inputs are the temperatures `a`, `b` [K]; outputs are `Q`, `-Q` [W],
/// positive into the conductor. The component owns the two conductivity
/// curves, so it can enter a kernel System without borrowing or leaking
/// material data. Curves cloned from a resolved source retain their exact
/// knots and temperature support; source authority and uncertainty are not
/// inferred from the numbers. A trial queries no database and allocates no
/// heap storage.
///
/// [`LayerGeometry`] supplies the existing planar, cylindrical or spherical
/// resistance coordinates. A heat port's temperature can be an endpoint,
/// midpoint or entropy-conjugate temperature according to its thermal
/// owner's discretization. The adapter does not change that choice or infer
/// phase from temperature, and it is not a transient field solver.
#[derive(Clone, Debug)]
pub struct LayerConductor {
    curves: [ConductivityCurve; 2],
    geometry: LayerGeometry,
    contact_resistance_m2_k_w: f64,
    config: LayerConductionConfig,
    ports: [PortSpec; 2],
}

fn poll(cx: &Cx<'_>) -> Result<(), DepositionError> {
    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
}

fn pair(values: &[f64]) -> Result<[f64; 2], DepositionError> {
    match values {
        [a, b] if a.is_finite() && b.is_finite() => Ok([*a, *b]),
        _ => Err(DepositionError::InvalidInput(
            "layer conductor requires two finite port values",
        )),
    }
}

fn output_width(values: &[f64]) -> Result<(), DepositionError> {
    if values.len() == 2 {
        Ok(())
    } else {
        Err(DepositionError::InvalidInput(
            "layer conductor requires two output slots",
        ))
    }
}

impl LayerConductor {
    /// Own the original curves and geometry, with a finite local solve budget.
    ///
    /// # Errors
    /// A negative/nonfinite contact resistance, zero iteration budget, or
    /// negative/nonfinite heat tolerance. Source support is checked at each
    /// trial's actual bulk and solved interface temperatures.
    pub fn new(
        curves: [ConductivityCurve; 2],
        geometry: LayerGeometry,
        contact_resistance_m2_k_w: f64,
        config: LayerConductionConfig,
    ) -> Result<Self, DepositionError> {
        geometry.conduction([&curves[0], &curves[1]], contact_resistance_m2_k_w)?;
        if config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(DepositionError::InvalidInput(
                "invalid layer-conduction heat tolerance or iteration budget",
            ));
        }
        let heat = |name| {
            StableId::new(name)
                .map_err(|_| DepositionError::InvalidInput("invalid layer port name"))
                .and_then(|id| {
                    PortSpec::heat(id, PortRole::EffortIn)
                        .map_err(|_| DepositionError::InvalidInput("invalid layer heat port"))
                })
        };
        Ok(Self {
            curves,
            geometry,
            contact_resistance_m2_k_w,
            config,
            ports: [heat("a")?, heat("b")?],
        })
    }

    /// The original ordered conductivity charts, including their exact support.
    #[must_use]
    pub fn conductivity_curves(&self) -> &[ConductivityCurve; 2] {
        &self.curves
    }

    /// The physical geometry, before lowering to planar resistance coordinates.
    #[must_use]
    pub const fn geometry(&self) -> LayerGeometry {
        self.geometry
    }

    /// The original scalar solve's work and heat-continuity controls.
    #[must_use]
    pub const fn config(&self) -> LayerConductionConfig {
        self.config
    }

    fn trial(
        &self,
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(LayerConductionResponse, Receipt), DepositionError> {
        poll(cx)?;
        let [a, b] = pair(inputs)?;
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(DepositionError::InvalidInput(
                "layer conductor requires a positive finite timestep",
            ));
        }
        // Geometry metrics are already prepared; this borrowed view copies
        // only references and scalars, not the owned charts or their knots.
        let law = self.geometry.conduction(
            [&self.curves[0], &self.curves[1]],
            self.contact_resistance_m2_k_w,
        )?;
        let response = law.evaluate([a, b], self.config, || poll(cx))?;
        let heat_j = dt_s * response.power_a_to_b_w;
        // Q and Ta-Tb have the same sign. Avoid subtracting nearly equal
        // reciprocals, and avoid forming the product Ta*Tb.
        let entropy = (heat_j / a) * ((a - b) / b);
        let allowance = 64.0 * f64::EPSILON * heat_j.abs() + f64::MIN_POSITIVE;
        // The native audit compares both layer powers. Recover the other
        // one from its retained signed hot-minus-cold residual.
        let hot_power = response.power_a_to_b_w.abs();
        let cold_power = hot_power - response.layer_heat_residual_w;
        let tolerance = self.config.absolute_heat_tolerance_w
            + self.config.relative_heat_tolerance * hot_power.max(cold_power.abs());
        if ![heat_j, entropy, allowance, tolerance]
            .into_iter()
            .all(f64::is_finite)
            || entropy < 0.0
        {
            return Err(DepositionError::Arithmetic);
        }
        Ok((
            response,
            Receipt {
                stored_change_j: 0.0,
                dissipated_j: 0.0,
                numerical_loss_j: 0.0,
                entropy_generated_j_per_k: entropy,
                energy_tolerance_j: allowance,
                residual: response.layer_heat_residual_w.abs(),
                tolerance,
            },
        ))
    }
}

impl Component for LayerConductor {
    type State = ();
    type Error = DepositionError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, DepositionError> {
        Ok(0.0)
    }

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-material::thermal::conduction::component::LayerConductor",
        state_family: "stateless-layer-conduction",
        schema_grammar: "unit state on the original conductivity charts and geometry",
        codec_grammar: "zero u64 words",
        codec_version: 1,
    });

    fn encode_state(&self, _state: &(), _words: &mut Vec<u64>) -> bool {
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<()> {
        words.is_empty().then_some(())
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        for curve in &self.curves {
            words.push(curve.knots().len() as u64);
            for knot in curve.knots() {
                words.extend(knot.map(f64::to_bits));
            }
        }
        words.extend(self.geometry.equivalent_distances_m().map(f64::to_bits));
        words.push(self.geometry.contact_area_m2().to_bits());
        // Equal resistance coordinates do not make different original
        // geometries interchangeable. Retain the geometry kind and inputs.
        match self.geometry.radial_coordinates_m() {
            None => words.push(0),
            Some((radii, contact)) => {
                words.push(if self.geometry.cylinder_length_m().is_some() {
                    1
                } else {
                    2
                });
                words.extend(radii.map(f64::to_bits));
                words.push(contact.to_bits());
                if let Some(length) = self.geometry.cylinder_length_m() {
                    words.push(length.to_bits());
                }
            }
        }
        words.extend([
            self.contact_resistance_m2_k_w.to_bits(),
            self.config.max_iterations as u64,
            self.config.absolute_heat_tolerance_w.to_bits(),
            self.config.relative_heat_tolerance.to_bits(),
        ]);
    }

    fn end_of_step(&self, _port: usize) -> EndOfStep {
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
    ) -> Result<Proposal, DepositionError> {
        output_width(outputs)?;
        let (response, receipt) = self.trial(dt_s, inputs, cx)?;
        poll(cx)?;
        *new = ();
        outputs.copy_from_slice(&[response.power_a_to_b_w, -response.power_a_to_b_w]);
        Ok(Proposal::Step(receipt))
    }
}

impl PortTangent for LayerConductor {
    fn jvp(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), DepositionError> {
        output_width(d_outputs)?;
        let [da, db] = pair(d_inputs)?;
        let (response, _) = self.trial(dt_s, inputs, cx)?;
        let [qa, qb] = response.endpoint_derivative_w_k;
        let tangent = qa * da + qb * db;
        if !tangent.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        poll(cx)?;
        d_outputs.copy_from_slice(&[tangent, -tangent]);
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
    ) -> Result<(), DepositionError> {
        output_width(bar_inputs)?;
        let [ba, bb] = pair(bar_outputs)?;
        let (response, _) = self.trial(dt_s, inputs, cx)?;
        let [qa, qb] = response.endpoint_derivative_w_k;
        // Temperature-dependent k makes qa != -qb in general. The
        // transpose is not another call to the forward action.
        let tangent = [qa * (ba - bb), qb * (ba - bb)];
        if tangent.iter().any(|value| !value.is_finite()) {
            return Err(DepositionError::Arithmetic);
        }
        poll(cx)?;
        bar_inputs.copy_from_slice(&tangent);
        Ok(())
    }
}
