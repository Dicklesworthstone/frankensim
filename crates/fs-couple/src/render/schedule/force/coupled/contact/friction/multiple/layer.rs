//! Source-bound thin material shear layers on the original simultaneous owner.
//!
//! The constitutive owner is fs-tribo's midpoint elastic/viscoplastic simple-shear
//! return map. This adapter supplies the actual held-force body compliance,
//! retains the material claims and publishes history with the original network.
//! Area, thickness, temperature and tangent frame are fixed. An explicitly
//! left-attached film can retain its material patch while the counterface moves.
//! Normal opening, layer transport, wear and pressure-dependent yield are
//! outside this retained-compressive-contact rung; none is replaced by a reset.
//! The explicit thermal adapter continues this same state with affine elastic
//! softening and its reversible heat; ordinary stepping remains isothermal.

/// Closed finite-cell heating and source-backed dissipative thermal feedback.
pub mod thermal;

use super::*;
use fs_material::state_point::ResolvedShearLayer;
use fs_tribo::layer_contact::{ThirdBodyLayerContact, ThirdBodyLayerState, ThirdBodyLayerRelaxationArm};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStep, TangentialStepInput};

/// Share the original source-to-law admission with the two-axis modal owner.
pub(in crate::render::schedule::force::coupled::contact::friction) fn material_layer_law(
    source: &ResolvedShearLayer,
) -> Result<ThirdBodyLayerContact, ModalFrictionError> {
    Ok(MaterialShearLayer::new(source.clone())?.law)
}

/// Independent reversible storage and physical plastic loss of one material
/// layer. The parent port's dissipation is this loss, never all traction work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShearLayerFrame {
    /// Whether the accepted interval produced irreversible plastic slip.
    pub plastic_flow: bool,
    /// Accepted elastic simple-shear displacement [m], in the retained tangent.
    pub elastic_slip_m: f64,
    /// Accepted accumulated plastic slip [m], in the same signed tangent.
    pub plastic_slip_m: f64,
    /// Accepted reversible layer energy [J].
    pub stored_energy_j: f64,
    /// Change of reversible layer energy during this accepted sample [J].
    pub stored_change_j: f64,
    /// Nonnegative yield, plastic-viscosity and Maxwell relaxation loss [J].
    pub dissipated_j: f64,
    /// Actual sub-yield/over-yield Maxwell-arm heat, included in dissipated_j.
    pub relaxation_dissipated_j: f64,
    /// Constitutive storage + loss - constitutive traction work [J].
    pub closure_residual_j: f64,
}

#[derive(Clone, Debug)]
pub(super) struct MaterialShearLayer {
    source: ResolvedShearLayer,
    law: ThirdBodyLayerContact,
    state: ThirdBodyLayerState,
    staged: Option<TangentialStep<ThirdBodyLayerState>>,
    /// Explicit material-side ownership, not permission to move the film itself.
    attached_to_left: bool,
}
impl MaterialShearLayer {
    fn new(source: ResolvedShearLayer) -> Result<Self, ModalFrictionError> {
        let law = ThirdBodyLayerContact::new(
            source.stiffness_n_m(),
            source.yield_force_n(),
            source.plastic_resistance_n_s_m(),
            source.temperature_kelvin(),
            source.maximum_elastic_slip_m(),
        )?.with_relaxation(source.relaxation().iter().map(|&(stiffness_n_m, relaxation_time_s)|
            ThirdBodyLayerRelaxationArm { stiffness_n_m, relaxation_time_s }).collect())?;
        let state = law.initial_state();
        Ok(Self {
            source,
            law,
            state,
            staged: None,
            attached_to_left: false,
        })
    }
    pub(super) fn energy_j(&self) -> f64 {
        self.law.stored_energy_j(&self.state)
    }
    pub(super) fn staged_energy_j(&self) -> f64 {
        self.staged
            .as_ref()
            .map_or(self.energy_j(), |s| self.law.stored_energy_j(&s.state))
    }
    fn trial(
        &self,
        free: f64,
        compliance: f64,
        normal: f64,
        dt: f64,
        config: ModalFrictionConfig,
    ) -> Result<TangentialStep<ThirdBodyLayerState>, ModalFrictionError> {
        Ok(self.law.advance(
            &self.state,
            &TangentialStepInput {
                normal_n: normal,
                slip_increment_m: [free, 0.0],
                dt_s: dt,
                contact_temperature_k: Some(self.source.temperature_kelvin()),
                nominal_pressure_pa: Some(finite(normal / self.source.area_m2())?),
                compliance_m_per_n: compliance,
                solve: TangentialSolveConfig {
                    max_iterations: config.max_iterations,
                    maximum_force_n: config.maximum_force_n,
                    force_absolute_tolerance_n: config.force_absolute_tolerance_n,
                    force_relative_tolerance: config.force_relative_tolerance,
                },
            },
        )?)
    }
    pub(super) fn solve(
        &self,
        free: f64,
        compliance: f64,
        normal: f64,
        dt: f64,
        config: ModalFrictionConfig,
    ) -> Result<(f64, ModalFrictionRegime, usize, Option<FrictionRoots>), ModalFrictionError> {
        // Every coordinate trial starts from the SAME accepted history.
        let step = self.trial(free, compliance, normal, dt, config)?;
        Ok((
            step.traction_n[0],
            ModalFrictionRegime::Shearing {
                plastic_flow: step.regime == fs_tribo::FrictionRegime::Sliding,
            },
            step.iterations,
            None,
        ))
    }
    fn actual(
        &self,
        slip: f64,
        traction: f64,
        normal: f64,
        dt: f64,
        config: ModalFrictionConfig,
        iterations: usize,
    ) -> Result<(TangentialStep<ThirdBodyLayerState>, f64), ModalFrictionError> {
        limit(
            "material-layer slip speed",
            finite(slip / dt)?.abs(),
            config.maximum_slip_speed_m_s,
        )?;
        limit(
            "material-layer traction",
            traction.abs(),
            config.maximum_force_n,
        )?;
        // Re-evaluate the constitutive owner on the ACTUAL shared endpoint,
        // with zero body compliance because the body motion is now known.
        let step = self.trial(slip, 0.0, normal, dt, config)?;
        let residual = finite(traction - step.traction_n[0])?;
        let tolerance = scalar::force_tolerance(config, traction, step.traction_n[0])?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::ContactSolve {
                residual_n: residual,
                tolerance_n: tolerance,
                iterations,
            }
            .into());
        }
        Ok((step, residual))
    }
    pub(super) fn check(
        &self,
        slip: f64,
        traction: f64,
        normal: f64,
        dt: f64,
        config: ModalFrictionConfig,
        iterations: usize,
    ) -> Result<(), ModalFrictionError> {
        self.actual(slip, traction, normal, dt, config, iterations)
            .map(|_| ())
    }
    pub(super) fn stage(
        &mut self,
        slip: f64,
        traction: f64,
        normal: f64,
        dt: f64,
        config: ModalFrictionConfig,
        iterations: usize,
    ) -> Result<(ShearLayerFrame, f64), ModalFrictionError> {
        let (step, residual) = self.actual(slip, traction, normal, dt, config, iterations)?;
        let frame = ShearLayerFrame {
            plastic_flow: step.regime == fs_tribo::FrictionRegime::Sliding,
            elastic_slip_m: step.state.elastic_slip_m()[0],
            plastic_slip_m: step.state.plastic_slip_m()[0],
            stored_energy_j: self.law.stored_energy_j(&step.state),
            stored_change_j: step.stored_change_j,
            dissipated_j: step.dissipated_j,
            relaxation_dissipated_j: self.law.step_dissipation(&self.state,
                &step.state, dt)?.relaxation_j,
            closure_residual_j: step.closure_residual_j,
        };
        self.staged = Some(step);
        Ok((frame, residual))
    }
    pub(super) fn publish(&mut self) {
        if let Some(step) = self.staged.take() {
            self.state = step.state;
        }
    }
}

impl MultiFrictionModalSystem {
    /// Replace the active dry traction mechanism of explicitly selected ports
    /// with their declared isothermal material layers, before playback.
    ///
    /// Each tuple is (original port, retained material + area + thickness +
    /// elastic-strain limit). Other contacts keep their original laws. The old
    /// dry law remains inspectable as construction data but is not evaluated
    /// on a layer-bound port; a dry-interface source cannot coexist there.
    ///
    /// The original full compliance, optional series hair drive, reciprocal
    /// receiver, normal-contact reaction and acoustic observer all participate.
    /// Source claims must cover static storage and the full represented temporal
    /// band through Nyquist. This admits a discretized constitutive model, not
    /// spectral accuracy or an identified rosin material.
    ///
    /// A layer needs positive compression throughout every accepted sample.
    /// Opening refuses without erasing storage; layer fracture/detachment is
    /// not implemented. Dynamic relocation refuses. Bound nonlinear mechanical
    /// storage uses the existing joint discrete-gradient iteration, retaining
    /// this layer's accepted history and reversible energy in every trial.
    pub fn with_material_shear_layers(
        self,
        bindings: Vec<(usize, ResolvedShearLayer)>,
        gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        self.bind_material_shear_layers(bindings, false, gate)
    }

    /// Retain each film on the original LEFT material patch and tangent frame.
    ///
    /// Only the opposing body's velocity/normal maps may migrate. The film's
    /// area, thickness, caloric cell, elastic/plastic and Maxwell history stay
    /// with the same left material. A geometry producer must preserve physical
    /// co-location, parallel contact directions and the finite film footprint.
    /// This supplies no film advection, wear, detachment or area evolution.
    /// Any mapped normal must use the same distinct pair of components.
    pub fn with_material_shear_layers_attached_to_left(
        self,
        bindings: Vec<(usize, ResolvedShearLayer)>,
        gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        self.bind_material_shear_layers(bindings, true, gate)
    }

    fn bind_material_shear_layers(
        mut self, bindings: Vec<(usize, ResolvedShearLayer)>, attached_to_left: bool,
        gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        poll(Some(gate))?;
        if bindings.is_empty() {
            return Ok(self);
        }
        if self.samples_rendered() != 0
            || self.shear_layers.iter().any(Option::is_some)
        {
            return Err(invalid("material shear layers require one complete binding before playback").into());
        }
        for (port, source) in bindings {
            poll(Some(gate))?;
            if port >= self.port_count()
                || self.shear_layers[port].is_some()
                || self.interface(port).is_some()
            {
                return Err(invalid("one material shear layer per actual port, without a competing dry-interface source").into());
            }
            source
                .require_angular_frequency_rad_s(0.0)
                .map_err(|_| invalid("material shear layer must support static elastic storage"))?;
            source.require_angular_frequency_rad_s(core::f64::consts::PI / self.sample_period_s())
                .map_err(|_| invalid("material shear layer must support the represented temporal band through Nyquist"))?;
            if attached_to_left {
                let pair = &self.points[port].port;
                if pair.right.as_ref().is_none_or(|right| right.component == pair.left.component) {
                    return Err(invalid("a retained left film needs a distinct declared counterface").into());
                }
                if let Some(normal) = self.contact_loads[port].and_then(|i| self.obstacle(i)) {
                    if normal.left.component != pair.left.component
                        || Some(normal.right.component) != pair.right.as_ref().map(|r| r.component) {
                        return Err(invalid("retained film normal and tangent must name the same material pair").into());
                    }
                }
            }
            let mut layer = MaterialShearLayer::new(source)?;
            layer.attached_to_left = attached_to_left;
            self.shear_layers[port] = Some(layer);
        }
        poll(Some(gate))?;
        Ok(self)
    }

    /// Original complete source and geometry for this layer-bound port.
    #[must_use]
    pub fn material_shear_layer(&self, port: usize) -> Option<&ResolvedShearLayer> {
        self.shear_layers
            .get(port)
            .and_then(Option::as_ref)
            .map(|layer| &layer.source)
    }

    /// Whether this film explicitly retains the original left material patch.
    #[must_use]
    pub fn material_shear_layer_attached_to_left(&self, port: usize) -> bool {
        self.shear_layers.get(port).and_then(Option::as_ref)
            .is_some_and(|layer| layer.attached_to_left)
    }

    /// Accepted constitutive history only; refused trials never appear here.
    #[must_use]
    pub fn material_shear_state(&self, port: usize) -> Option<&ThirdBodyLayerState> {
        self.shear_layers
            .get(port)
            .and_then(Option::as_ref)
            .map(|layer| &layer.state)
    }
}
