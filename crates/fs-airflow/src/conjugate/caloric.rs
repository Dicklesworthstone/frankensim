//! Source-exact Cp(T) in steady, single-phase wall heat exchange.
//!
//! The original material law integrates
//! `integral Cp(T)/(T_wall-T) dT = h*A/mass_flow` on its admitted isobar.
//! Each segment reports `Q = mass_flow * integral Cp(T) dT` and the Robin
//! reference `T_wall-Q/(h*A)`. The same wall/stream transfer therefore has one
//! physical heat rate even when Cp changes strongly along the channel.
//!
//! Mass flow, pressure, composition and h are prescribed. Only actual fluid
//! temperatures require source support; a hotter or colder wall may lie outside
//! it. This does not add pressure work, phase change, h(T), hydraulic feedback,
//! mixing junctions or hydraulic/solid feedback derivatives. Its inlet and wall
//! responses differentiate this same original source-caloric march. The
//! constant-Cp adjoint/transfer APIs accept `AirPath`, never this distinct path.

/// Finite, source-caloric co-current and counterflow heat exchangers.
pub mod exchanger;

mod linearization;
pub use linearization::{CaloricPathLinearization, CaloricPathTangent, CaloricSegmentTangent};

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::thermal::transport::isobaric::{IsobaricHeatCapacityLaw, LinearHeatExchange};

use super::{
    AirMarch, AirPath, AirSegment, AirflowError, ConjugateConfig, ConjugateSolution, IqnIls,
    IqnIlsConfig, SegmentState, SolidRegionState, conjugate_iqn_error, effectiveness,
    finite_coupling, finite_positive, solve_conjugate_driver_with_march,
};

/// One immutable original Cp source and the native ordered heat-exchange path.
/// Cloning shares the source and its prepared caloric index.
#[derive(Clone, Debug, PartialEq)]
pub struct CaloricAirPath {
    geometry: AirPath,
    caloric: Arc<IsobaricHeatCapacityLaw>,
}

/// Analytic inlet response of the original source-caloric march with fixed
/// wall temperatures, mass flow, segment conductances and material source.
/// These derivatives do not include solid, hydraulic or composition feedback.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricInletResponse {
    /// d h_out / d h_in, dimensionless and in [0,1]. Strong heat exchange can
    /// underflow to zero; rounded equality with a wall does not itself zero it.
    pub enthalpy_gain: f64,
    /// d T_out / d T_in = enthalpy_gain * Cp(T_in) / Cp(T_out).
    /// Unlike the enthalpy gain, this can exceed one for a decreasing Cp law.
    pub temperature_gain: f64,
}

impl CaloricAirPath {
    /// Bind a fixed mass flow to the original source-supported Cp(T) isobar.
    /// Region ownership and geometric/flow representability use `AirPath`'s
    /// existing admission. No source uncertainty or extrapolation is inferred.
    pub fn new(
        inlet_temperature_k: f64,
        mass_flow_kg_per_s: f64,
        caloric: Arc<IsobaricHeatCapacityLaw>,
        segments: Vec<AirSegment>,
    ) -> Result<Self, AirflowError> {
        let cp = caloric
            .specific_heat_cp_j_kg_k(inlet_temperature_k)
            .map_err(AirflowError::CaloricSource)?;
        let geometry = AirPath::new(inlet_temperature_k, mass_flow_kg_per_s, cp, segments)?;
        Ok(Self { geometry, caloric })
    }

    /// Original source identity, pressure, support, receipts and uncertainty.
    #[must_use]
    pub fn caloric_law(&self) -> &IsobaricHeatCapacityLaw {
        &self.caloric
    }

    /// Native wall regions, in actual inlet-to-outlet order.
    #[must_use]
    pub fn segments(&self) -> &[AirSegment] {
        self.geometry.segments()
    }

    /// Prescribed inlet absolute temperature [K].
    #[must_use]
    pub const fn inlet_temperature_k(&self) -> f64 {
        self.geometry.inlet_temperature_k()
    }

    /// Prescribed positive mass flow [kg/s].
    #[must_use]
    pub const fn mass_flow_kg_per_s(&self) -> f64 {
        self.geometry.mass_flow_kg_per_s
    }

    /// Integrate the original Cp(T) across every traversed source interval.
    /// Relative enthalpy increments retain small heat transfers independently
    /// of the source's enthalpy datum. A late source refusal returns no march.
    pub fn march(&self, wall_temperatures_k: &[f64]) -> Result<AirMarch, AirflowError> {
        self.march_with_enthalpy_gain(wall_temperatures_k)
            .map(|(march, _)| march)
    }

    /// March once and retain the analytic fixed-wall inlet response.
    /// Each source segment uses d h_out/d h_in = (T_wall-T_out)/(T_wall-T_in)
    /// before rounding absolute temperatures, with the analytic exponential
    /// limit at equilibrium. Original affine Cp integration retains very small
    /// gains even when an outlet rounds to its wall temperature. No numerical
    /// differentiation, frozen Cp or extrapolated wall calorimetry is used.
    pub fn march_with_inlet_response(
        &self,
        wall_temperatures_k: &[f64],
    ) -> Result<(AirMarch, CaloricInletResponse), AirflowError> {
        let (march, enthalpy_gain) = self.march_with_enthalpy_gain(wall_temperatures_k)?;
        let inlet_cp = self
            .caloric
            .specific_heat_cp_j_kg_k(self.inlet_temperature_k())
            .map_err(AirflowError::CaloricSource)?;
        let outlet_cp = self
            .caloric
            .specific_heat_cp_j_kg_k(march.outlet_temperature_k)
            .map_err(AirflowError::CaloricSource)?;
        let temperature_gain = finite_coupling(
            "source-caloric inlet temperature response",
            (enthalpy_gain * inlet_cp) / outlet_cp,
        )?;
        Ok((
            march,
            CaloricInletResponse {
                enthalpy_gain,
                temperature_gain,
            },
        ))
    }

    pub(crate) fn march_with_enthalpy_gain(
        &self,
        wall_temperatures_k: &[f64],
    ) -> Result<(AirMarch, f64), AirflowError> {
        self.march_with_response(wall_temperatures_k, |_, _| Ok(()))
    }

    fn march_with_response(
        &self,
        wall_temperatures_k: &[f64],
        mut retain: impl FnMut(&LinearHeatExchange, f64) -> Result<(), AirflowError>,
    ) -> Result<(AirMarch, f64), AirflowError> {
        if wall_temperatures_k.len() != self.segments().len() {
            return Err(AirflowError::SolidResponseArity {
                expected: self.segments().len(),
                found: wall_temperatures_k.len(),
            });
        }
        let mut inlet = self.inlet_temperature_k();
        let mut states = Vec::with_capacity(self.segments().len());
        let mut total = 0.0;
        let mut correction = 0.0;
        let mut enthalpy_gain = 1.0;
        for (segment, &wall) in self.segments().iter().zip(wall_temperatures_k) {
            finite_positive("caloric wall temperature", wall)?;
            let conductance = finite_positive(
                "air segment conductance",
                segment.area_m2() * segment.htc_w_per_m2_k(),
            )?;
            let exposure = finite_positive(
                "caloric exchange per mass flow",
                conductance / self.mass_flow_kg_per_s(),
            )?;
            let cp = self
                .caloric
                .specific_heat_cp_j_kg_k(inlet)
                .map_err(AirflowError::CaloricSource)?;
            let ntu = finite_positive("caloric inlet NTU", exposure / cp)?;
            let exchange = self
                .caloric
                .linear_heat_exchange(inlet, wall, exposure)
                .map_err(AirflowError::CaloricSource)?;
            let outlet = exchange.temperature_k;
            enthalpy_gain *= exchange.inlet_enthalpy_gain;
            // Retain the source integral on the relative temperature change.
            // Recovering heat from an already rounded absolute outlet amplifies
            // that rounding when a weak exchange forms its Robin reference.
            let heat = exchange.specific_enthalpy_change_j_kg;
            let heat_rate_w =
                finite_coupling("source-caloric heat rate", self.mass_flow_kg_per_s() * heat)?;
            if heat != 0.0 && heat_rate_w == 0.0 {
                return Err(AirflowError::InvalidConjugateInput {
                    field: "unrepresentable source-caloric heat rate",
                    value_bits: heat_rate_w.to_bits(),
                });
            }
            let driving = finite_coupling("source-caloric temperature difference", wall - inlet)?;
            let eps = if driving == 0.0 {
                effectiveness(ntu).0
            } else {
                (outlet - inlet) / driving
            };
            let reference = finite_positive(
                "source-caloric Robin reference",
                exchange.mean_temperature_k,
            )?;
            // This reference is a positive-weight mean of the actual stream
            // temperatures. Reject a resolved inconsistency without clipping heat.
            let rounding = 8.0 * f64::EPSILON * inlet.max(outlet).max(wall);
            if !(0.0..=1.0).contains(&eps)
                || reference < inlet.min(outlet) - rounding
                || reference > inlet.max(outlet) + rounding
            {
                return Err(AirflowError::InvalidConjugateInput {
                    field: "source-caloric exchange resolution",
                    value_bits: reference.to_bits(),
                });
            }
            retain(&exchange, exposure)?;
            states.push(SegmentState {
                region: segment.region().to_owned(),
                inlet_temperature_k: inlet,
                outlet_temperature_k: outlet,
                reference_temperature_k: reference,
                ntu,
                effectiveness: eps,
                heat_rate_w,
            });
            let increment = heat_rate_w - correction;
            let next = finite_coupling("total source-caloric heat rate", total + increment)?;
            correction = (next - total) - increment;
            total = next;
            inlet = outlet;
        }
        Ok((
            AirMarch {
                segments: states,
                outlet_temperature_k: inlet,
                total_heat_rate_w: total,
            },
            enthalpy_gain,
        ))
    }

    /// Solve the original solid/stream fixed point with source Cp(T), retaining
    /// native scalar relaxation, cancellation, region checks and watt gates.
    pub fn solve<F>(
        &self,
        cx: &Cx<'_>,
        config: &ConjugateConfig,
        solid: F,
    ) -> Result<ConjugateSolution, AirflowError>
    where
        F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
    {
        let initial = vec![self.inlet_temperature_k(); self.segments().len()];
        self.solve_from(cx, config, &initial, solid)
    }

    /// Resume from explicit reference temperatures using the same native driver.
    /// Fixed relaxation replays an identically partitioned continuation;
    /// Aitken history is reinitialized as in `solve_conjugate_from`.
    pub fn solve_from<F>(
        &self,
        cx: &Cx<'_>,
        config: &ConjugateConfig,
        initial_references_k: &[f64],
        solid: F,
    ) -> Result<ConjugateSolution, AirflowError>
    where
        F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
    {
        solve_conjugate_driver_with_march(
            cx,
            &self.geometry,
            config,
            initial_references_k,
            None,
            solid,
            |walls| self.march(walls),
        )
    }

    /// Apply native bounded IQN-ILS to the same source-caloric fixed point.
    /// Both temperature and per-region watt convergence remain mandatory.
    pub fn solve_iqn<F>(
        &self,
        cx: &Cx<'_>,
        config: &ConjugateConfig,
        acceleration: IqnIlsConfig,
        solid: F,
    ) -> Result<ConjugateSolution, AirflowError>
    where
        F: FnMut(&Cx<'_>, &[f64]) -> Result<Vec<SolidRegionState>, AirflowError>,
    {
        let initial = vec![self.inlet_temperature_k(); self.segments().len()];
        let accelerator =
            IqnIls::new(self.segments().len(), acceleration).map_err(conjugate_iqn_error)?;
        solve_conjugate_driver_with_march(
            cx,
            &self.geometry,
            config,
            &initial,
            Some(accelerator),
            solid,
            |walls| self.march(walls),
        )
    }
}
