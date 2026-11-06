//! Held kinetic stress with source-resolved molecular power at each thermal state.
//! The original constitutive root and enthalpy timeline remain authoritative.

use super::*;

/// Source and physical stress-power sample at one thermal evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceViscousFeedbackSample {
    pub coordinates: [usize; 3],
    /// Original pinned material response at this temperature and strain rate.
    pub material: ResolvedLiquidState,
    /// Force-corrected sqrt(2 dev(D):dev(D)) [1/s], using this same viscosity.
    pub shear_rate_per_s: f64,
    /// Physical molecular collision coefficient used in that reconstruction.
    pub collision_tau: f64,
    /// Located original liquid mass / reference density [m3].
    pub liquid_volume_m3: f64,
    /// Shear plus tied bulk power `2 mu D:D` [W/m3], excluding numerical loss.
    pub power_w_m3: f64,
}

/// Accepted molecular heat on the original thermal interval, before VOF advection.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceViscousFeedbackReport {
    /// Sum of the actual accepted substep power integrals [J].
    pub heat_j: f64,
    /// Last source samples actually applied by the thermal integrator.
    pub last_applied: Vec<FreeSurfaceViscousFeedbackSample>,
    /// Fresh source-supported endpoint samples; these supply no second heat.
    pub thermal_endpoint: Vec<FreeSurfaceViscousFeedbackSample>,
}

#[derive(Clone)]
struct StressCell {
    id: CellId,
    coordinates: [usize; 3],
    volume_m3: f64,
    /// Guo-corrected deviatoric gamma_lattice * tau, held during thermal splitting.
    numerator: f64,
    /// Independent spherical part of the held stress, with the same tied tau.
    dilation_numerator: f64,
}

#[derive(Clone)]
pub(super) struct ViscousFeedback {
    cells: Vec<StressCell>,
    dx: f64,
    dt: f64,
    maximum_speed: f64,
    // These cache only the immutable source/stress binding of this candidate.
    // Never reuse a response across a different enthalpy bit pattern or step.
    last_applied_h: Vec<u64>,
    thermal_endpoint_h: Vec<u64>,
    pub(super) report: FreeSurfaceViscousFeedbackReport,
}

impl ViscousFeedback {
    fn cached(&self, id: CellId, h: f64) -> Option<&FreeSurfaceViscousFeedbackSample> {
        // The original thermal owner assigns consecutive liquid IDs before
        // any finite-body IDs. Verify the binding instead of trusting an index.
        let index = usize::try_from(id.0).ok()?;
        if self.cells.get(index)?.id != id {
            return None;
        }
        if self.last_applied_h.get(index) == Some(&h.to_bits()) {
            self.report.last_applied.get(index)
        } else if self.thermal_endpoint_h.get(index) == Some(&h.to_bits()) {
            self.report.thermal_endpoint.get(index)
        } else {
            None
        }
    }

    pub(super) fn cached_material(&self, id: CellId, h: f64) -> Option<&ResolvedLiquidState> {
        self.cached(id, h).map(|sample| &sample.material)
    }

    pub(super) fn prepare(
        liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
        ids: &BTreeMap<MassLocation3, CellId>,
        thermal: &DepositionSimulation<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        let mut cells = Vec::with_capacity(ids.len());
        for (&location, &id) in ids {
            poll(gate)?;
            let coordinates = FreeSurface3::transport_coordinates(location).ok_or(
                Error::Invalid("molecular power requires a located liquid cell"),
            )?;
            let [numerator, dilation_numerator] =
                liquid.flow.lattice.transport_strain_rates(location, 1.0)?;
            if [numerator, dilation_numerator].iter().any(|value| !value.is_finite() || *value < 0.0) {
                return Err(Error::Invalid("unrepresentable liquid stress numerator"));
            }
            cells.push(StressCell {
                id,
                coordinates,
                volume_m3: thermal.cells()[&id].geometry.volume_m3,
                numerator,
                dilation_numerator,
            });
        }
        Ok(Self {
            cells,
            dx: liquid.flow.scaling.cell_size_m(),
            dt: liquid.flow.scaling.timestep_s(),
            maximum_speed: liquid.flow.scaling.flow.u_lattice * liquid.flow.scaling.velocity_scale,
            last_applied_h: Vec::new(),
            thermal_endpoint_h: Vec::new(),
            report: FreeSurfaceViscousFeedbackReport {
                heat_j: 0.0,
                last_applied: Vec::new(),
                thermal_endpoint: Vec::new(),
            },
        })
    }
}

impl Conduction<'_, '_> {
    fn sample_viscous_feedback(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<FreeSurfaceViscousFeedbackSample>, Error> {
        let Some(feedback) = &self.viscous_feedback else {
            return Ok(Vec::new());
        };
        let mut samples = Vec::with_capacity(feedback.cells.len());
        for cell in &feedback.cells {
            check()?;
            let h = thermal.cells()[&cell.id].specific_h_j_kg;
            if let Some(sample) = feedback.cached(cell.id, h) {
                samples.push(sample.clone());
                continue;
            }
            let (material, rate, tau) = if let Some(model) = &self.rheology {
                let resolved = model.sample(h, self.rheology_stress[&cell.id], || {
                    check().map_err(Error::from)
                })?;
                (
                    resolved.material,
                    resolved.sample.shear_rate_per_s,
                    resolved.sample.collision_tau,
                )
            } else {
                let (material, transport) =
                    resolve_liquid(self.law, self.caloric, h, self.density, self.buoyancy)?;
                let (scaling, _) = plan_momentum_scaling(
                    transport.kinematic_viscosity_m2_s(),
                    feedback.dx,
                    feedback.dt,
                    feedback.maximum_speed,
                )?;
                let rate = (cell.numerator / scaling.tau) / feedback.dt;
                (material, rate, scaling.tau)
            };
            let dilation = (cell.dilation_numerator / tau) / feedback.dt;
            let mu = material.dynamic_viscosity_pa_s();
            let shear_power = (mu * rate) * rate;
            let bulk_power = (mu * dilation) * dilation;
            let power = shear_power + bulk_power;
            if !rate.is_finite()
                || rate < 0.0
                || (cell.numerator != 0.0 && rate == 0.0)
                || !dilation.is_finite()
                || (cell.dilation_numerator != 0.0 && dilation == 0.0)
                || !power.is_finite()
                || power < 0.0
                || (rate != 0.0 && shear_power == 0.0)
                || (dilation != 0.0 && bulk_power == 0.0)
            {
                return Err(Error::Invalid(
                    "liquid molecular feedback power is unrepresentable",
                ));
            }
            samples.push(FreeSurfaceViscousFeedbackSample {
                coordinates: cell.coordinates,
                material,
                shear_rate_per_s: rate,
                collision_tau: tau,
                liquid_volume_m3: cell.volume_m3,
                power_w_m3: power,
            });
        }
        Ok(samples)
    }

    pub(super) fn add_viscous_feedback(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        powers: &mut BTreeMap<CellId, f64>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        let samples = self.sample_viscous_feedback(thermal, check)?;
        let Some(feedback) = &mut self.viscous_feedback else {
            return Ok(());
        };
        for (cell, sample) in feedback.cells.iter().zip(&samples) {
            check()?;
            let watts = environment::finite(sample.power_w_m3 * sample.liquid_volume_m3)?;
            if sample.power_w_m3 != 0.0 && watts == 0.0 {
                return Err(Error::Invalid("liquid molecular feedback power underflows"));
            }
            if watts != 0.0 {
                let total = powers.entry(cell.id).or_insert(0.0);
                *total = environment::finite(*total + watts)?;
            }
        }
        feedback.last_applied_h = feedback
            .cells
            .iter()
            .map(|cell| thermal.cells()[&cell.id].specific_h_j_kg.to_bits())
            .collect();
        feedback.report.last_applied = samples;
        Ok(())
    }

    pub(super) fn accept_viscous_feedback(
        &mut self,
        duration_s: f64,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        let endpoint = self.sample_viscous_feedback(thermal, check)?;
        let Some(feedback) = &mut self.viscous_feedback else {
            return Ok(());
        };
        for sample in &feedback.report.last_applied {
            check()?;
            let heat = (sample.power_w_m3 * sample.liquid_volume_m3) * duration_s;
            if !heat.is_finite() || heat < 0.0 || (sample.power_w_m3 != 0.0 && heat == 0.0) {
                return Err(Error::Invalid(
                    "liquid molecular feedback heat is unrepresentable",
                ));
            }
            feedback.report.heat_j = environment::finite(feedback.report.heat_j + heat)?;
        }
        feedback.thermal_endpoint_h = feedback
            .cells
            .iter()
            .map(|cell| thermal.cells()[&cell.id].specific_h_j_kg.to_bits())
            .collect();
        feedback.report.thermal_endpoint = endpoint;
        Ok(())
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Conduct heat and refresh molecular viscous power on every thermal substep.
    ///
    /// This opt-in holds the original force-corrected nonequilibrium stress
    /// numerator and, for pressure-dependent rheology, original mechanical
    /// pressure throughout the thermal split. Each actual h resolves the SAME
    /// material claims and constitutive root, so mu(T,gamma), collision tau and
    /// gamma share one source state. Physical strain is not held fixed while
    /// viscosity changes. The original thermal integrator deposits mu*gamma²
    /// once per accepted substep; endpoint diagnostics add no second heat.
    ///
    /// The accepted thermal endpoint then drives one original momentum/VOF
    /// step and conservative mass*h transport. This is first-order splitting,
    /// not simultaneous kinetic/thermal evolution or exact total-energy closure.
    /// The explicit thermal step still requires nonlinear time refinement.
    /// Authored heat must exclude molecular work. Bulk viscosity, pressure
    /// work, phase heat and numerical damping are not inferred. Source,
    /// thermal, native-flow and cancellation refusals publish nothing.
    pub fn step_with_conduction_and_viscous_feedback(
        &mut self,
        faces: &[FreeSurfaceConductionFace],
        heat: &[FreeSurfaceCellHeat],
        controls: FreeSurfaceConductionControls,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceConductionReport, Error> {
        self.advance_with_viscous_feedback(
            faces,
            heat,
            controls,
            &[],
            &mut [],
            &[],
            &[],
            &[],
            &[],
            &[],
            FreeSurfaceInterfaceSolve::default(),
            None,
            None,
            None,
            true,
            gate,
        )
    }
}
