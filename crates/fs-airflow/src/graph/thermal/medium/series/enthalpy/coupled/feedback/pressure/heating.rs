//! Hydraulic dissipation and wall exchange on the same liquid trajectory.
//! Constant-density, equal-area/elevation flow converts the declared pressure
//! loss entirely to fluid heat. Viscosity uses the regional arithmetic bulk
//! mean; convection retains the caller's bulk/film convention. The original
//! caloric exponential integrates both contributions together.

use fs_convection::ThermalDirection;
use fs_exec::Cx;
use fs_material::phase::EquilibriumPhaseState;
use fs_qty::{MassFlowRate, Pressure, Temperature, VolumetricFlowRate};

use super::super::{
    LiquidFeedbackBudget, LiquidFeedbackError, LiquidHeatMarch, LiquidHeatPath, LiquidHeatRegion,
    Section, checkpoint, positive,
};
use super::{
    LiquidHydraulicRegion, LiquidPressureState, LiquidResistanceSection, LiquidResistanceModel, PressureDrivenBudget,
    PressureDrivenError, PressureDrivenSolution, solve_pressure_root,
};
use crate::graph::thermal::medium::FrozenFlowMedium;
use crate::graph::thermal::medium::series::enthalpy::{
    admit_dissipative_density, march_region_with_heat, state_at_temperature,
};

/// A jointly resolved heat/pressure trajectory at a prescribed mass flow.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidDissipativeMarch<'a> {
    /// Wall-only exchanges and source-backed fluid endpoint states.
    pub thermal: LiquidHeatMarch<'a>,
    /// Original pressure-loss evaluations whose heat entered that trajectory.
    pub hydraulics: LiquidPressureState,
    /// Independent actual endpoint `mdot*(h_out-h_in)`, W.
    pub enthalpy_gain_w: f64,
    /// Actual advected energy minus wall and hydraulic heat reports, W.
    /// This is numerical reporting error, never additional physical heat.
    pub balance_residual_w: f64,
}

impl<'a, 'card> LiquidHeatPath<'a, 'card> {
    /// Resolve material properties, hydraulic heat and wall exchange together.
    ///
    /// Choosing this model declares equal elevation, constant source density,
    /// equal flow areas, unchanged bulk kinetic energy and complete local
    /// conversion of the declared loss into fluid heat. It does not model pump
    /// inefficiency, gas throttling or heat retained by a porous skeleton.
    /// Each regional solve starts at the inlet and uses predicted outlets to
    /// find a sign bracket, then bisects. It never probes an unrelated caloric
    /// maximum. All trials still require original source/correlation support.
    ///
    /// # Errors
    /// Original source/Reynolds/correlation refusal, cancellation, invalid
    /// geometry, or failure to close temperature, HTC and hydraulic power.
    pub fn march_with_dissipation(
        &self,
        cx: &Cx<'_>,
        walls: &[f64],
        resistances: &[LiquidResistanceSection<'_, impl LiquidResistanceModel>],
        budget: LiquidFeedbackBudget,
    ) -> Result<LiquidDissipativeMarch<'a>, PressureDrivenError> {
        checkpoint(cx)?;
        budget.admit()?;
        self.admit_resistances(resistances)?;
        admit_dissipative_density(self.source, self.law.reference_state().density_kg_m3())?;
        let area = self.sections[0].geometry.flow_area;
        if self.sections.iter().any(|s| s.geometry.flow_area != area) {
            return Err(PressureDrivenError::Invalid(
                "hydraulic heating requires equal flow areas",
            ));
        }
        if walls.len() != self.sections.len() {
            return Err(PressureDrivenError::Invalid(
                "one wall temperature is required per region",
            ));
        }
        for &wall in walls {
            state_at_temperature(self.source.curve(), wall)?;
        }
        let mut inlet = self.inlet;
        let mut thermal = Vec::with_capacity(walls.len());
        let mut hydraulic = Vec::with_capacity(walls.len());
        let (mut wall_heat, mut pressure, mut hydraulic_heat) = (0.0, 0.0, 0.0);
        for ((section, &wall), resistance) in self.sections.iter().zip(walls).zip(resistances) {
            let (heat, loss) =
                self.solve_heated_region(cx, section, inlet, wall, resistance, budget)?;
            wall_heat += heat.heat_rate_w;
            pressure += loss.loss.pressure_drop_pa;
            hydraulic_heat += loss.loss.dissipated_power_w;
            if !wall_heat.is_finite() || !pressure.is_finite() || !hydraulic_heat.is_finite() {
                return Err(PressureDrivenError::Invalid(
                    "aggregate heated-flow report is unrepresentable",
                ));
            }
            inlet = heat.outlet;
            thermal.push(heat);
            hydraulic.push(loss);
        }
        let gain =
            self.mass_flow * (inlet.specific_enthalpy_j_kg() - self.inlet.specific_enthalpy_j_kg());
        let (partial, a) = fs_math::eft::two_sum(gain, -wall_heat);
        let (rounded, b) = fs_math::eft::two_sum(partial, -hydraulic_heat);
        let residual = rounded + (a + b);
        if !gain.is_finite() || !residual.is_finite() {
            return Err(PressureDrivenError::Invalid(
                "advected energy or heat balance is unrepresentable",
            ));
        }
        self.admit_porous_trajectory(cx, resistances, &thermal, &mut hydraulic)?;
        checkpoint(cx)?;
        Ok(LiquidDissipativeMarch {
            thermal: LiquidHeatMarch {
                source: self.source,
                regions: thermal,
                mass_flow_kg_per_s: self.mass_flow,
                total_heat_rate_w: wall_heat,
            },
            hydraulics: LiquidPressureState {
                regions: hydraulic,
                total_pressure_drop_pa: pressure,
                total_dissipated_power_w: hydraulic_heat,
            },
            enthalpy_gain_w: gain,
            balance_residual_w: residual,
        })
    }

    /// Select mass flow by pressure while retaining all declared loss as heat.
    /// Uses the original bounded pressure root, source laws and local feedback.
    /// The omitted-power allowance must be exactly zero for this opt-in path.
    ///
    /// # Errors
    /// As `march_with_dissipation`, plus original pressure/bracket/work refusals.
    pub fn solve_pressure_driven_with_dissipation(
        &self,
        cx: &Cx<'_>,
        walls: &[f64],
        resistances: &[LiquidResistanceSection<'_, impl LiquidResistanceModel>],
        drive: Pressure,
        local: LiquidFeedbackBudget,
        budget: PressureDrivenBudget,
    ) -> Result<PressureDrivenSolution<LiquidHeatMarch<'a>>, PressureDrivenError> {
        checkpoint(cx)?;
        budget.admit(drive)?;
        local.admit()?;
        self.admit_resistances(resistances)?;
        if budget.maximum_neglected_hydraulic_power_w != 0.0 {
            return Err(PressureDrivenError::Invalid(
                "fully dissipative heat path requires zero omitted-power allowance",
            ));
        }
        solve_pressure_root(cx, drive, budget, |mass| {
            let path = self.at_pressure_trial_mass(cx, mass)?;
            let march = path.march_with_dissipation(cx, walls, resistances, local)?;
            Ok((march.thermal, march.hydraulics))
        })
    }

    fn solve_heated_region(
        &self,
        cx: &Cx<'_>,
        section: &Section,
        inlet: EquilibriumPhaseState,
        wall: f64,
        resistance: &LiquidResistanceSection<'_, impl LiquidResistanceModel>,
        budget: LiquidFeedbackBudget,
    ) -> Result<(LiquidHeatRegion, LiquidHydraulicRegion), PressureDrivenError> {
        let mut guess = inlet.temperature_k();
        let (mut low, mut high): (Option<f64>, Option<f64>) = (None, None);
        for _ in 0..=budget.max_region_iterations {
            let (candidate, hydraulic) =
                self.evaluate_heated_region(cx, section, inlet, wall, guess, resistance)?;
            let residual = candidate.temperature_residual_k;
            let predicted = candidate.outlet.temperature_k();
            if residual.abs() <= budget.temperature_tolerance_k {
                let (mut checked, checked_loss) =
                    self.evaluate_heated_region(cx, section, inlet, wall, predicted, resistance)?;
                let before = candidate.boundary.coefficient().value.value();
                let after = checked.boundary.coefficient().value.value();
                checked.coefficient_residual_relative = (after - before).abs() / before.max(after);
                let p0 = hydraulic.loss.dissipated_power_w;
                let p1 = checked_loss.loss.dissipated_power_w;
                if checked.temperature_residual_k.abs() <= budget.temperature_tolerance_k
                    && checked.coefficient_residual_relative
                        <= budget.coefficient_relative_tolerance
                    && (p1 - p0).abs() / p0.max(p1) <= budget.coefficient_relative_tolerance
                {
                    return Ok((checked, checked_loss));
                }
            }
            if residual >= 0.0 {
                low = Some(guess);
            } else {
                high = Some(guess);
            }
            let next = match (low, high) {
                (Some(a), Some(b)) => 0.5 * a + 0.5 * b,
                _ => predicted,
            };
            if next == guess {
                break;
            }
            guess = next;
        }
        Err(LiquidFeedbackError::NotConverged {
            region: section.name.clone(),
        }
        .into())
    }

    fn evaluate_heated_region(
        &self,
        cx: &Cx<'_>,
        section: &Section,
        inlet: EquilibriumPhaseState,
        wall: f64,
        guess: f64,
        resistance: &LiquidResistanceSection<'_, impl LiquidResistanceModel>,
    ) -> Result<(LiquidHeatRegion, LiquidHydraulicRegion), PressureDrivenError> {
        checkpoint(cx)?;
        let bulk = 0.5 * inlet.temperature_k() + 0.5 * guess;
        let hydraulic_material = self
            .law
            .at_temperature(Temperature::new(bulk))
            .map_err(LiquidFeedbackError::from)?;
        let loss = resistance
            .resistance
            .evaluate_material(&hydraulic_material, bulk, MassFlowRate::new(self.mass_flow))
            .map_err(|source| PressureDrivenError::Resistance {
                region: section.name.clone(),
                source,
            })?;
        let property_temperature_k =
            self.property_temperature
                .evaluate(inlet.temperature_k(), guess, wall);
        let material = if property_temperature_k == bulk {
            hydraulic_material.clone()
        } else {
            self.law
                .at_temperature(Temperature::new(property_temperature_k))
                .map_err(LiquidFeedbackError::from)?
        };
        let volume_flow = self.mass_flow / material.density_kg_m3();
        let channel = FrozenFlowMedium::Liquid(&material)
            .channel(
                &section.name,
                VolumetricFlowRate::new(volume_flow),
                Temperature::new(inlet.temperature_k()),
                section.geometry,
                section.correlation,
                // Pressure heating can make an initially isothermal wall cool
                // the fluid. Direction-sensitive correlations use the same
                // regional bulk estimate as this accepted trajectory.
                if wall >= bulk {
                    ThermalDirection::HeatingFluid
                } else {
                    ThermalDirection::CoolingFluid
                },
            )
            .map_err(LiquidFeedbackError::from)?;
        let segment = &channel.path().segments()[0];
        let outlet = march_region_with_heat(
            cx,
            self.source.curve(),
            inlet,
            self.mass_flow,
            segment,
            wall,
            loss.dissipated_power_w,
        )?;
        let heat_rate_w = self.mass_flow
            * (outlet.specific_enthalpy_j_kg() - inlet.specific_enthalpy_j_kg())
            - loss.dissipated_power_w;
        let reference = wall - heat_rate_w / (segment.htc_w_per_m2_k() * segment.area_m2());
        let velocity_m_s = volume_flow / section.geometry.flow_area.value();
        if !heat_rate_w.is_finite() || !positive(reference) || !positive(velocity_m_s) {
            return Err(PressureDrivenError::Invalid(
                "heated regional heat, reference or velocity is unrepresentable",
            ));
        }
        let boundary = channel
            .robin_boundary(Temperature::new(reference))
            .map_err(LiquidFeedbackError::from)?;
        let correlation = channel.correlation().clone();
        checkpoint(cx)?;
        Ok((
            LiquidHeatRegion {
                region: section.name.clone(),
                material,
                property_temperature_k,
                volume_flow_m3_s: volume_flow,
                velocity_m_s,
                inlet,
                outlet,
                heat_rate_w,
                boundary,
                correlation,
                temperature_residual_k: outlet.temperature_k() - guess,
                coefficient_residual_relative: 0.0,
            },
            LiquidHydraulicRegion {
                region: section.name.clone(),
                material: hydraulic_material,
                loss,
                porous: None,
            },
        ))
    }
}
