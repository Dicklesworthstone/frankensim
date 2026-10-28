//! Material viscosity at the local physical temperature and shear rate.
//!
//! The existing non-equilibrium stress estimators give `gamma_dot * tau` in
//! lattice units. Resolve `tau = 1/2 + 3 nu(T, gamma_dot) dt / dx^2` together
//! with `gamma_dot = numerator / (tau dt)`; the background tau is only an
//! initial iterate. This avoids interpreting a state with the wrong viscosity
//! after step-local transport restores the background relaxation field.
//!
//! This is an explicit generalized-Newtonian, fixed-storage BGK coupling.
//! Local constitutive convergence is not implicit time integration, viscoelastic
//! memory, a stability guarantee or experimental validation. The base entry
//! adds no viscous heat; [`heating`] explicitly thermalizes the same samples.

/// Opt-in conversion of this constitutive stress power into sensible heat.
pub mod heating;

use super::{
    HeatSourceLimits, MaterialSlab, MaterialTransportError, MaterialTransportStep, ThermalLbm,
    ThermalLbm3, TransportSample, advance_with_sampler, poll, sample_temperatures_with,
    temperatures,
};
use crate::material::plan_scaling_from_material;
use crate::thermal::sources::Slab;
use fs_exec::CancelGate;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
};

/// Exact physical coordinates supplied to the constitutive material owner.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RheologyPoint {
    /// Linear fluid-cell index, with x fastest; wall reservoirs are omitted.
    pub cell: usize,
    /// Absolute temperature [K] at the sampled time level.
    pub temperature_k: f64,
    /// `sqrt(2 D:D)` [1/s], where `D = (grad u + grad u^T)/2`.
    pub shear_rate_per_s: f64,
}

/// Bounded constitutive solve with a damped search and bracketed safeguards.
#[derive(Debug, Clone, Copy)]
pub struct RheologySolve {
    /// Maximum material evaluations per fluid cell and time level.
    pub maximum_iterations: usize,
    /// Fraction of an unbracketed tau update, in `(0, 1]`. Once evaluated
    /// residuals bracket a root, safeguarded secant/bisection takes over.
    pub relaxation: f64,
    /// Absolute residual tolerance for dimensionless lattice relaxation time.
    pub absolute_relaxation_tolerance: f64,
    /// Relative tau residual tolerance, in `[0, 1)`.
    pub relative_relaxation_tolerance: f64,
}
impl RheologySolve {
    pub(crate) fn validate(self) -> Result<(), MaterialTransportError> {
        if self.maximum_iterations == 0
            || !self.relaxation.is_finite()
            || self.relaxation <= 0.0
            || self.relaxation > 1.0
            || !self.absolute_relaxation_tolerance.is_finite()
            || self.absolute_relaxation_tolerance < 0.0
            || !self.relative_relaxation_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_relaxation_tolerance)
            || (self.absolute_relaxation_tolerance == 0.0
                && self.relative_relaxation_tolerance == 0.0)
        {
            return Err(MaterialTransportError::Invalid(
                "invalid local rheology solve controls",
            ));
        }
        Ok(())
    }
}

/// The same constitutive root is used by fixed-domain and free-surface flow.
pub(crate) struct LocalRheology<F> {
    pub value: F,
    pub shear_rate_per_s: f64,
    pub iterations: usize,
    pub relaxation_defect: f64,
}

pub(crate) enum LocalRheologyError<E> {
    Evaluation(E),
    Invalid(&'static str),
    Convergence {
        iterations: usize,
        relaxation_defect: f64,
    },
}

/// Resolve the force-corrected stress numerator against the material's physical
/// BGK relaxation. Callers admit controls and units and retain source authority.
/// A bracket is learned only from the existing damped search: no global tau
/// endpoints or extra source-domain assumptions are introduced. Subsequent tau
/// candidates remain strictly between sampled values; the source owner still
/// admits every resulting physical rate.
pub(crate) fn solve_local<F, E>(
    numerator: f64,
    timestep_s: f64,
    initial_tau: f64,
    solve: RheologySolve,
    mut check: impl FnMut() -> Result<(), E>,
    mut evaluate: impl FnMut(f64) -> Result<(F, f64), E>,
) -> Result<LocalRheology<F>, LocalRheologyError<E>> {
    let mut tau = initial_tau;
    let mut relaxation_defect = f64::INFINITY;
    let mut previous: Option<(f64, f64)> = None;
    let mut bracket: Option<[(f64, f64); 2]> = None;
    let mut previous_width = f64::INFINITY;
    for iterations in 1..=solve.maximum_iterations {
        check().map_err(LocalRheologyError::Evaluation)?;
        let rate = (numerator / tau) / timestep_s;
        if !rate.is_finite() || rate < 0.0 {
            return Err(LocalRheologyError::Invalid(
                "physical shear rate is unrepresentable",
            ));
        }
        let (value, next) = evaluate(rate).map_err(LocalRheologyError::Evaluation)?;
        let residual = tau - next;
        relaxation_defect = residual.abs();
        let tolerance = solve.absolute_relaxation_tolerance
            + solve.relative_relaxation_tolerance * next.max(tau);
        if !tolerance.is_finite() {
            return Err(LocalRheologyError::Invalid(
                "rheology tolerance is unrepresentable",
            ));
        }
        if relaxation_defect <= tolerance {
            return Ok(LocalRheology {
                value,
                shear_rate_per_s: rate,
                iterations,
                relaxation_defect,
            });
        }
        if let Some(ends) = &mut bracket {
            // The candidate lies inside the sampled interval. Replace the
            // endpoint with the same residual sign, retaining opposite signs
            // without assuming a particular constitutive slope or continuity.
            let same_as_lower = (residual > 0.0) == (ends[0].1 > 0.0);
            ends[usize::from(!same_as_lower)] = (tau, residual);
        } else if let Some(old) = previous
            && (residual > 0.0) != (old.1 > 0.0)
        {
            bracket = Some(if old.0 < tau {
                [old, (tau, residual)]
            } else {
                [(tau, residual), old]
            });
        }
        previous = Some((tau, residual));
        if let Some([(lower, f_lower), (upper, f_upper)]) = bracket {
            let width = upper - lower;
            // Scale first to avoid overflowing the residual sum. The secant
            // point is a convex combination of evaluated, admissible taus.
            let scale = f_lower.abs().max(f_upper.abs());
            let a = (f_lower / scale).abs();
            let b = (f_upper / scale).abs();
            let secant = lower + (a / (a + b)) * width;
            let candidate = if width <= 0.5 * previous_width && secant > lower && secant < upper {
                secant
            } else {
                // Do not let regula falsi stagnate at a steep curve or kink:
                // the interval must halve at least every two evaluations.
                lower + 0.5 * width
            };
            if candidate <= lower || candidate >= upper {
                // A tiny interval is not evidence of a small constitutive
                // residual (e.g. for a discontinuous source response).
                return Err(LocalRheologyError::Convergence {
                    iterations,
                    relaxation_defect,
                });
            }
            previous_width = width;
            tau = candidate;
        } else {
            tau += solve.relaxation * (next - tau);
        }
    }
    Err(LocalRheologyError::Convergence {
        iterations: solve.maximum_iterations,
        relaxation_defect,
    })
}

/// Original material response and the converged physical coordinates it used.
#[derive(Debug, Clone, PartialEq)]
pub struct RheologyMaterial<F> {
    source: F,
    point: RheologyPoint,
    transport: FluidTransportState,
    iterations: usize,
    relaxation_defect: f64,
}
impl<F> RheologyMaterial<F> {
    /// Original material owner, including its property/source receipts.
    pub const fn source(&self) -> &F {
        &self.source
    }
    /// Actual source query; no nominal or output-tuned rate is substituted.
    pub const fn point(&self) -> RheologyPoint {
        self.point
    }
    /// Number of local material evaluations at this accepted time level.
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    /// Absolute fixed-point residual in tau at the accepted query.
    pub const fn relaxation_defect(&self) -> f64 {
        self.relaxation_defect
    }
}
impl<F> FluidTransport for RheologyMaterial<F> {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        // A concrete source owner is read once for each material query.
        Ok(self.transport)
    }
}

trait StressSlab<const D: usize, const Q: usize>: Slab<D, Q> {
    fn shear_numerator(&self, cell: usize) -> f64;
}
impl StressSlab<2, 5> for ThermalLbm {
    fn shear_numerator(&self, cell: usize) -> f64 {
        let m = self.grid.moments(cell);
        let force = std::array::from_fn(|d| self.grid.g[d].mul_add(m.rho, self.grid.fext[cell][d]));
        crate::core2::force_corrected_shear_rate(&self.grid.f[cell], m.rho, m.u, force, 1.0)
    }
}
impl StressSlab<3, 7> for ThermalLbm3 {
    fn shear_numerator(&self, cell: usize) -> f64 {
        let (rho, u) = self.grid.moments3(cell);
        crate::d3q19::force_corrected_shear_rate3(
            &self.grid.f[cell],
            rho,
            u,
            self.grid.fext[cell],
            1.0,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn sample<S: StressSlab<D, Q>, F: FluidTransport, const D: usize, const Q: usize>(
    mapping: &MaterialSlab<S>,
    state: &S,
    cells: usize,
    plane: usize,
    range: [f64; 2],
    solve: RheologySolve,
    gate: &CancelGate,
    resolve: &mut impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
) -> Result<Vec<TransportSample<RheologyMaterial<F>>>, MaterialTransportError> {
    let temperature = temperatures::<S, D, Q>(mapping, state, range, gate)?;
    // Validate the complete stress field before asking any source to resolve it.
    let mut numerator = vec![0.0; cells];
    for cell in plane..cells - plane {
        poll(gate)?;
        numerator[cell] = state.shear_numerator(cell);
        if !numerator[cell].is_finite() || numerator[cell] < 0.0 {
            return Err(MaterialTransportError::Invalid(
                "unrepresentable non-equilibrium shear stress",
            ));
        }
    }
    let units = mapping.force_units;
    sample_temperatures_with::<S, RheologyMaterial<F>, D, Q>(
        mapping,
        &temperature,
        cells,
        plane,
        range,
        gate,
        &mut |cell, temperature_k, retained| {
            let root = solve_local(
                numerator[cell],
                units.timestep_s,
                state.momentum_relaxation()[cell],
                solve,
                || poll(gate),
                |rate| {
                    let point = RheologyPoint {
                        cell,
                        temperature_k,
                        shear_rate_per_s: rate,
                    };
                    let error = |source| MaterialTransportError::State {
                        cell,
                        temperature_k,
                        source,
                    };
                    let source = resolve(point, retained).map_err(|e| error(e.into()))?;
                    let transport = source.transport_state().map_err(|e| error(e.into()))?;
                    let scale = plan_scaling_from_material(
                        &transport,
                        units.cell_size_m,
                        units.timestep_s,
                        0.0,
                    )
                    .map_err(error)?;
                    Ok(((source, transport), scale.flow().tau))
                },
            )
            .map_err(|error| match error {
                LocalRheologyError::Evaluation(error) => error,
                LocalRheologyError::Invalid(what) => MaterialTransportError::Invalid(what),
                LocalRheologyError::Convergence {
                    iterations,
                    relaxation_defect,
                } => MaterialTransportError::RheologyConvergence {
                    cell,
                    iterations,
                    relaxation_defect,
                },
            })?;
            Ok(RheologyMaterial {
                source: root.value.0,
                transport: root.value.1,
                point: RheologyPoint {
                    cell,
                    temperature_k,
                    shear_rate_per_s: root.shear_rate_per_s,
                },
                iterations: root.iterations,
                relaxation_defect: root.relaxation_defect,
            })
        },
    )
}

#[allow(clippy::too_many_arguments)]
fn advance<S: StressSlab<D, Q>, F: FluidTransport, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>,
    power: &[f64],
    forces: &[[f64; D]],
    depth_m: f64,
    limits: HeatSourceLimits,
    solve: RheologySolve,
    gate: &CancelGate,
    mut resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
) -> Result<MaterialTransportStep<RheologyMaterial<F>>, MaterialTransportError> {
    poll(gate)?;
    solve.validate()?;
    advance_with_sampler::<S, RheologyMaterial<F>, D, Q>(
        slab,
        power,
        forces,
        depth_m,
        limits,
        gate,
        |mapping, state, cells, plane| {
            sample::<S, F, D, Q>(
                mapping,
                state,
                cells,
                plane,
                limits.temperature_range_k,
                solve,
                gate,
                &mut resolve,
            )
        },
    )
}

impl MaterialSlab<ThermalLbm> {
    /// Advance the same 2D heat/flow kernels with viscosity resolved at both
    /// actual temperature [K] and self-consistent shear rate [1/s].
    ///
    /// The callback owns material assignment, fixed pressure/composition and
    /// source policy. It may return `ResolvedLiquidState` from a card whose
    /// viscosity curve uses the supplied rate. Density and heat capacity must
    /// retain their reference values, as in `step_with_material_transport`.
    /// No viscosity floor, cap, material-name switch or rheology fit is added.
    ///
    /// Both initial and candidate fields must converge inside source support.
    /// Coefficients at the initial field drive one explicit step; endpoint
    /// queries check that candidate, rather than retroactively changing it.
    /// A refused source, nonconvergence or cancellation preserves all solver
    /// populations, loads, coefficients and the physical-unit mapping. Local
    /// iteration polls cancellation between source evaluations; callbacks must
    /// themselves return. Viscous heating is not added to the sensible-heat law.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_rheology<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 2]],
        depth_m: f64,
        limits: HeatSourceLimits,
        solve: RheologySolve,
        gate: &CancelGate,
        resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<MaterialTransportStep<RheologyMaterial<F>>, MaterialTransportError> {
        advance::<ThermalLbm, F, 2, 5>(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            solve,
            gate,
            resolve,
        )
    }
}
impl MaterialSlab<ThermalLbm3> {
    /// 3D counterpart of `MaterialSlab::step_with_material_rheology` using
    /// the same constitutive solve, unit conversions, balances and publication.
    pub fn step_with_material_rheology<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        limits: HeatSourceLimits,
        solve: RheologySolve,
        gate: &CancelGate,
        resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<MaterialTransportStep<RheologyMaterial<F>>, MaterialTransportError> {
        let depth_m = self.force_units.cell_size_m;
        advance::<ThermalLbm3, F, 3, 7>(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            solve,
            gate,
            resolve,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::{LocalRheologyError, RheologySolve, solve_local};
    use std::cell::Cell;

    fn controls(relaxation: f64) -> RheologySolve {
        RheologySolve {
            maximum_iterations: 100,
            relaxation,
            absolute_relaxation_tolerance: 1e-13,
            relative_relaxation_tolerance: 1e-12,
        }
    }

    // Continuous positive viscosity [Pa s], with a threefold thickening
    // transition. At rho=1000 kg/m^3 and dx=dt=0.01 SI, tau=0.5+0.3*mu.
    fn thickening_tau(rate: f64) -> f64 {
        let mu = (1.0 + 10.0 * (rate - 1.2)).clamp(1.0, 3.0);
        0.5 + 0.3 * mu
    }

    #[test]
    fn strong_material_thickening_converges_inside_sampled_source_support() {
        for (relaxation, initial_tau) in [(0.5, 0.8), (0.75, 0.8), (1.0, 0.8), (0.75, 1.2)] {
            let polls = Cell::new(0usize);
            let mut calls = 0;
            let mut queried_rates = Vec::new();
            let result = solve_local(
                0.016,
                0.01,
                initial_tau,
                controls(relaxation),
                || {
                    polls.set(polls.get() + 1);
                    Ok(())
                },
                |rate| {
                    // Blind probes at a numerical tau floor would leave this
                    // source's physical support despite its admissible root.
                    if !(1.1..=2.0).contains(&rate) {
                        return Err("outside material support");
                    }
                    calls += 1;
                    queried_rates.push(rate);
                    Ok(((calls, rate), thickening_tau(rate)))
                },
            );
            let root = match result {
                Ok(root) => root,
                Err(LocalRheologyError::Convergence {
                    iterations,
                    relaxation_defect,
                }) => panic!("{iterations} queries, tau defect {relaxation_defect}"),
                Err(_) => panic!("an admissible material query was refused"),
            };
            // gamma*tau=1.6/s, mu=10*gamma-11 in the transition:
            // tau^2+2.8*tau-4.8=0, so tau=1.2 and gamma=4/3 /s.
            assert!((root.shear_rate_per_s - 4.0 / 3.0).abs() < 2e-12);
            assert!((thickening_tau(root.shear_rate_per_s) - 1.2).abs() < 2e-12);
            assert!(root.relaxation_defect < 2e-12);
            assert_eq!(root.iterations, calls);
            assert_eq!(polls.get(), calls);
            assert_eq!(root.value, (calls, *queried_rates.last().unwrap()));
            assert_eq!(root.value.1, root.shear_rate_per_s);
            assert!(calls <= 20, "{calls} material queries");
            if initial_tau == 1.2 {
                assert_eq!(calls, 1, "an already consistent material needs one query");
            }
        }
    }

    #[test]
    fn bracketed_material_solve_retains_budget_refusal_and_cancellation() {
        let mut calls = 0;
        let result = solve_local(
            0.016,
            0.01,
            0.8,
            RheologySolve {
                maximum_iterations: 2,
                ..controls(0.75)
            },
            || Ok::<_, &'static str>(()),
            |rate| {
                calls += 1;
                Ok(((), thickening_tau(rate)))
            },
        );
        assert_eq!(calls, 2);
        assert!(matches!(
            result,
            Err(LocalRheologyError::Convergence { iterations: 2, .. })
        ));

        for cancelled in [false, true] {
            let polls = Cell::new(0usize);
            let mut calls = 0;
            let result = solve_local(
                0.016,
                0.01,
                0.8,
                controls(0.75),
                || {
                    polls.set(polls.get() + 1);
                    if cancelled && polls.get() == 3 {
                        Err("cancelled")
                    } else {
                        Ok(())
                    }
                },
                |rate| {
                    calls += 1;
                    if calls == 3 {
                        Err("source refused interior query")
                    } else {
                        Ok(((), thickening_tau(rate)))
                    }
                },
            );
            let expected = if cancelled {
                "cancelled"
            } else {
                "source refused interior query"
            };
            assert!(matches!(result, Err(LocalRheologyError::Evaluation(e)) if e == expected));
            assert_eq!(polls.get(), 3);
            assert_eq!(calls, if cancelled { 2 } else { 3 });
        }
    }

    #[test]
    fn narrow_material_bracket_cannot_hide_a_finite_constitutive_defect() {
        let mut calls = 0;
        let result = solve_local(
            0.016,
            0.01,
            0.8,
            controls(0.75),
            || Ok::<_, ()>(()),
            |rate| {
                calls += 1;
                // Positive viscosities, but a jump across the requested
                // stress: there is no self-consistent constitutive root.
                Ok(((), if rate > 4.0 / 3.0 { 1.4 } else { 0.8 }))
            },
        );
        assert!(matches!(
            result,
            Err(LocalRheologyError::Convergence {
                iterations,
                relaxation_defect,
            }) if iterations == calls && relaxation_defect > 0.19
        ));
        assert!(calls <= controls(0.75).maximum_iterations);
    }
}
