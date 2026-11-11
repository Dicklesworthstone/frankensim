//! Generalized-Newtonian source response on the retained native MAC timeline.
//!
//! The same symmetric-stress quadrature owns the physical rate invariant and
//! molecular heat. Source viscosity is evaluated at accepted-start T/shear;
//! native inertia and the existing thermal split remain unchanged. No BGK tau,
//! cell speed, LES loss, Darcy loss or numerical time loss is a shear rate.

use super::*;

mod quasisteady;

/// Numerical closure of source response against the actual native strain.
/// These controls apply to both stationary flow and an endpoint step from a
/// fixed accepted momentum history; they do not assert constitutive uniqueness.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RheologyControls {
    pub max_iterations: usize,
    /// Each cell must meet this absolute OR the relative physical-rate gate.
    pub absolute_rate_tolerance_per_s: f64,
    pub relative_rate_tolerance: f64,
    /// Re-query the actual rate: viscosity, conductivity, Cp and any active
    /// buoyancy density departure must agree within this relative bound.
    pub relative_transport_tolerance: f64,
    /// Convex relaxation of the next source-rate query, in (0,1].
    pub relaxation: f64,
}

/// Closure of the coefficients actually used by native momentum and heat.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RheologyReport {
    pub iterations: usize,
    pub maximum_rate_mismatch_per_s: f64,
    pub maximum_relative_rate_mismatch: f64,
    pub maximum_relative_transport_mismatch: f64,
}

// The stationary public names denote these same numerical controls/receipts.
pub use RheologyControls as QuasisteadyRheologyControls;
pub use RheologyReport as QuasisteadyRheologyReport;

impl RheologyControls {
    pub(super) fn admit(self) -> Result<(), DepositionError> {
        if self.max_iterations == 0
            || !self.absolute_rate_tolerance_per_s.is_finite()
            || self.absolute_rate_tolerance_per_s < 0.0
            || !self.relative_rate_tolerance.is_finite()
            || !(0.0 < self.relative_rate_tolerance && self.relative_rate_tolerance <= 1.0)
            || !self.relative_transport_tolerance.is_finite()
            || !(0.0 < self.relative_transport_tolerance
                && self.relative_transport_tolerance <= 1.0)
            || !self.relaxation.is_finite()
            || !(0.0 < self.relaxation && self.relaxation <= 1.0)
        {
            return Err(invalid(
                "rheology needs a positive iteration budget, finite rate/transport gates and relaxation in (0,1]",
            ));
        }
        Ok(())
    }
}

fn relative_difference(a: f64, b: f64) -> f64 {
    if a == b {
        0.0
    } else {
        (a - b).abs() / a.abs().max(b.abs())
    }
}

/// Smallest nonzero shear-rate magnitude a source's surface query admits: 2^-200,
/// the lower bound of fs-matdb's surface coordinate domain (`MIN_COORDINATE` in its
/// `surface.rs`). A rate of exactly zero is admitted separately.
const MIN_QUERY_RATE_PER_S: f64 = f64::from_bits((1023 - 200) << 52);

/// The rate a source is queried at for the native rate `rate`. The strain quadrature
/// can leave a nonzero tail far below any physical rate (8.3e-162 1/s in a cell of
/// the buoyancy fixture), which the surface guard refuses because it is under 2^-200.
/// Such a rate is queried as exactly 0, which is a legal query (cold starts and solids
/// are queried there): the law's value at 8e-162 differs from its value at 0 by about
/// (lambda x 8e-162)^a, far below any tolerance. Only the query changes: the stored
/// native rate keeps its value.
fn query_rate(rate: f64) -> f64 {
    if rate.abs() < MIN_QUERY_RATE_PER_S {
        0.0
    } else {
        rate
    }
}

/// Immutable constitutive source on actual caloric and kinematic coordinates.
/// Implementations must retain original source claims and refuse unsupported
/// states. They must not mutate shared external physical history. Exact solids
/// retain their phase-first source semantics and receive zero native rate.
pub trait FlowRheologySource<F>: Sync {
    fn resolve(
        &self,
        voxel: usize,
        cell: &DepositedCell,
        temperature_k: f64,
        shear_rate_per_s: f64,
    ) -> Result<F, MaterialError>;
}
impl<F, Q> FlowRheologySource<F> for Q
where
    Q: Fn(usize, &DepositedCell, f64, f64) -> Result<F, MaterialError> + Sync,
{
    fn resolve(
        &self,
        voxel: usize,
        cell: &DepositedCell,
        temperature_k: f64,
        shear_rate_per_s: f64,
    ) -> Result<F, MaterialError> {
        self(voxel, cell, temperature_k, shear_rate_per_s)
    }
}

#[derive(Clone)]
pub(super) struct RheologyBinding<'a, F> {
    source: &'a dyn FlowRheologySource<F>,
    pub(super) rates: Vec<f64>,
    pub(super) quasisteady: Option<QuasisteadyRheologyControls>,
    report: Option<QuasisteadyRheologyReport>,
}

impl<'a, F: Clone, R, O: DepositionObserver, T: FlowTransportSnapshot>
    MaterialFlowFeedback<'a, F, R, O, T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Bind T/shear-dependent material response to the existing retained flow.
    ///
    /// Native symmetric-strain quadrature yields the cell-volume RMS
    /// `sqrt(2 dev(D):dev(D))` in 1/s, including graded edges and no-slip wall
    /// half cells. It observes the exact retained face velocities; no
    /// cell-speed reconstruction or constitutive relaxation parameter enters.
    /// Source coefficients are frozen over each accepted interval. The SAME
    /// viscosity drives full-stress momentum and its optional molecular heat;
    /// source k and admitted Darcy closure remain on the same material epoch.
    ///
    /// This explicitly time-lagged generalized-Newtonian constitutive law needs
    /// timestep/mesh refinement. It adds no elastic/yield history, nonlinear
    /// implicit solve, phase kinetics, moving solid skeleton or free surface.
    /// Attach retained momentum first; quasisteady flow has no accepted native
    /// starting velocity and refuses. Both pre-advance and actual endpoint T/rate
    /// are source-admitted. Any late failure discards the whole original thermal,
    /// mechanical and wrapped-observer candidate, including finite gas exchange.
    pub fn with_native_rheology(
        mut self,
        thermal: &DepositionSimulation<'_>,
        source: &'a dyn FlowRheologySource<F>,
    ) -> Result<Self, MaterialFlowError<O::Error>> {
        check_gate(self.gate, &mut || Ok(()))?;
        self.check_map(thermal, &mut || Ok(()))?;
        self.check_unsteady_start(thermal, &mut || Ok(()))?;
        if self.unsteady.is_none() || self.rheology.is_some() {
            return Err(
                invalid("native rheology needs one source and retained unsteady momentum").into(),
            );
        }
        self.rheology = Some(RheologyBinding {
            source,
            rates: Vec::new(),
            quasisteady: None,
            report: None,
        });
        self.endpoint = self.sample_materials(thermal, &mut || Ok(()))?;
        admit_density(&self.endpoint, self.density_kg_m3)?;
        self.mushy_drag(&self.endpoint, &mut || Ok(()))?;
        self.check_porous_endpoint(&mut || Ok(()))?;
        self.check_phase_conduction(&self.endpoint, &mut || Ok(()))?;
        Ok(self)
    }

    /// Accepted endpoint rates in native voxel order. Source-start rates are
    /// retained separately on `applied()`, matching its physical coefficients.
    /// Before a quasisteady advance this is the explicitly admitted initial
    /// seed; after it, this is the actual solved rate, not a new thermal-endpoint
    /// equilibrium solve. The next source epoch closes that equilibrium.
    pub fn endpoint_shear_rate_per_s(&self) -> Option<&[f64]> {
        self.rheology
            .as_ref()
            .map(|binding| binding.rates.as_slice())
    }

    pub(super) fn sample_materials(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<FlowMaterialSample<F, T>>, MaterialFlowError<O::Error>> {
        self.sample_materials_at_rates(thermal, None, check)
    }

    /// Explicit rate iterates are numerical source coordinates. They never
    /// replace accepted native velocity or advance its momentum history.
    pub(super) fn sample_materials_at_rates(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        rates: Option<Vec<f64>>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<FlowMaterialSample<F, T>>, MaterialFlowError<O::Error>> {
        let Some(binding) = &self.rheology else {
            return sample(
                thermal,
                self.domain,
                self.input,
                self.gate,
                &mut self.resolve,
                self.exact_phase.then_some(self.mushy).flatten(),
                self.snapshot,
                check,
            );
        };
        let source = binding.source;
        let mut rates = if let Some(rates) = rates {
            rates
        } else if binding.quasisteady.is_some() {
            binding.rates.clone()
        } else {
            let owner = self
                .unsteady
                .as_ref()
                .ok_or(invalid("native rheology lost its momentum owner"))?;
            let active = self.momentum_domain(&self.endpoint, check)?;
            let domain = active.as_ref().unwrap_or(self.domain);
            if domain.fluid_count() == 0 {
                vec![0.0; domain.cell_count()]
            } else {
                let fluid = self.reference.map(|_| self.reference_fluid()).transpose()?;
                owner.shear_rates(domain, fluid.as_ref(), self.config, self.gate)?
            }
        };
        let phase_law = self.exact_phase.then_some(self.mushy).flatten();
        let mut resolve = |voxel: usize, cell: &DepositedCell, temperature: f64| {
            // Native strain still uses the previously admitted topology until
            // the phase handoff below. A newly exact-solid cell can inherit a
            // neighbouring edge contribution even when its own faces are at
            // rest. Its CURRENT caloric phase owns this source query: solids
            // receive zero, without clipping any live momentum or its energy.
            if let Some(law) = phase_law {
                if matches!(
                    law.classify_at_specific_enthalpy(cell.specific_h_j_kg)?,
                    PhaseMobility::ImmobileSolid { .. }
                ) {
                    rates[voxel] = 0.0;
                }
            }
            source.resolve(voxel, cell, temperature, query_rate(rates[voxel]))
        };
        let samples = sample(
            thermal,
            self.domain,
            self.input,
            self.gate,
            &mut resolve,
            phase_law,
            self.snapshot,
            check,
        )?;
        self.rheology.as_mut().expect("binding retained").rates = rates;
        Ok(samples)
    }

    /// One shared admission and mismatch calculation for stationary and
    /// unsteady constitutive iterations. Observed sources remain distinct from
    /// the queried samples that supplied the candidate's momentum and power.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(super) fn observe_rheology(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        epoch: &MaterialFlowEpoch<F, T>,
        rates: Vec<f64>,
        iteration: usize,
        controls: RheologyControls,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(RheologyReport, bool, Vec<FlowMaterialSample<F, T>>), MaterialFlowError<O::Error>>
    {
        let queried = epoch
            .shear_rate_per_s
            .as_ref()
            .ok_or(invalid("bound rheology lacks source rates"))?;
        if rates.len() != queried.len() || rates.len() != self.domain.cell_count() {
            return Err(
                invalid("native and source shear fields must share the original grid").into(),
            );
        }
        let mut report = RheologyReport {
            iterations: iteration,
            maximum_rate_mismatch_per_s: 0.0,
            maximum_relative_rate_mismatch: 0.0,
            maximum_relative_transport_mismatch: 0.0,
        };
        let mut rate_closed = true;
        for (&actual, &query) in rates.iter().zip(queried) {
            check_gate(self.gate, check)?;
            let error = (actual - query).abs();
            let relative = relative_difference(actual, query);
            report.maximum_rate_mismatch_per_s = report.maximum_rate_mismatch_per_s.max(error);
            report.maximum_relative_rate_mismatch =
                report.maximum_relative_rate_mismatch.max(relative);
            rate_closed &= error <= controls.absolute_rate_tolerance_per_s
                || relative <= controls.relative_rate_tolerance;
        }
        // Source-admit the ACTUAL physical rate before any relaxation.
        let observed = self.sample_materials_at_rates(thermal, Some(rates), check)?;
        admit_density(&observed, self.density_kg_m3)?;
        self.check_phase_conduction(&observed, check)?;
        for (applied, actual) in epoch.samples.iter().zip(&observed) {
            check_gate(self.gate, check)?;
            if let Some((gravity, density)) = self.buoyancy {
                if gravity.iter().any(|&g| g != 0.0) {
                    report.maximum_relative_transport_mismatch = report
                        .maximum_relative_transport_mismatch
                        .max(relative_difference(
                            density(&applied.material) - self.density_kg_m3,
                            density(&actual.material) - self.density_kg_m3,
                        ));
                }
            }
            report.maximum_relative_transport_mismatch = report
                .maximum_relative_transport_mismatch
                .max(relative_difference(
                    applied.transport.thermal_conductivity_w_m_k(),
                    actual.transport.thermal_conductivity_w_m_k(),
                ));
            match (
                applied.transport.mobile_transport(),
                actual.transport.mobile_transport(),
            ) {
                (Some(a), Some(b)) => {
                    for (x, y) in [
                        (a.dynamic_viscosity_pa_s(), b.dynamic_viscosity_pa_s()),
                        (a.specific_heat_j_kg_k(), b.specific_heat_j_kg_k()),
                    ] {
                        report.maximum_relative_transport_mismatch = report
                            .maximum_relative_transport_mismatch
                            .max(relative_difference(x, y));
                    }
                }
                (None, None) => {}
                _ => {
                    return Err(invalid(
                        "constitutive iteration changed occupancy at fixed caloric state",
                    )
                    .into());
                }
            }
        }
        let closed = rate_closed
            && report.maximum_relative_transport_mismatch <= controls.relative_transport_tolerance;
        if closed {
            self.rheology.as_mut().expect("binding retained").report = Some(report);
        }
        Ok((report, closed, observed))
    }
}

#[cfg(test)]
mod tests {
    use super::{MIN_QUERY_RATE_PER_S, query_rate};

    #[test]
    fn a_rate_tail_under_the_surface_domain_queries_as_zero_and_a_rate_inside_it_is_unchanged() {
        assert_eq!(MIN_QUERY_RATE_PER_S, 2.0_f64.powi(-200));
        // 8.3e-162 is the tail the buoyancy fixture's strain quadrature leaves (measured).
        for tail in [
            1e-161,
            8.316801694027249e-162,
            2.0_f64.powi(-201),
            f64::MIN_POSITIVE,
        ] {
            assert_eq!(query_rate(tail).to_bits(), 0.0_f64.to_bits(), "{tail}");
        }
        // 2^-200 is the smallest admitted magnitude; 2^-199 and everything above pass unchanged.
        for kept in [2.0_f64.powi(-200), 2.0_f64.powi(-199), 1e-12, 0.5, 300.0] {
            assert_eq!(query_rate(kept).to_bits(), kept.to_bits(), "{kept}");
        }
        assert_eq!(query_rate(0.0).to_bits(), 0.0_f64.to_bits());
    }
}
