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

/// Bounded damped fixed-point controls for each local constitutive query.
#[derive(Debug, Clone, Copy)]
pub struct RheologySolve {
    /// Maximum material evaluations per fluid cell and time level.
    pub maximum_iterations: usize,
    /// Fraction of the proposed tau update, in `(0, 1]`.
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
    for iterations in 1..=solve.maximum_iterations {
        check().map_err(LocalRheologyError::Evaluation)?;
        let rate = (numerator / tau) / timestep_s;
        if !rate.is_finite() || rate < 0.0 {
            return Err(LocalRheologyError::Invalid(
                "physical shear rate is unrepresentable",
            ));
        }
        let (value, next) = evaluate(rate).map_err(LocalRheologyError::Evaluation)?;
        relaxation_defect = (next - tau).abs();
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
        tau += solve.relaxation * (next - tau);
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

// Add the half-force transformation before reusing the existing stress
// estimators. Its second moment is +(u F^T + F u^T)/2. This is the raw Guo
// source, WITHOUT the collision prefactor (1 - 1/(2 tau)). See Premnath,
// Pattison and Banerjee, arXiv:0901.2786, equations (5), (7), and (B23):
// https://arxiv.org/pdf/0901.2786 . This restores the physical non-equilibrium
// second moment of the time-explicit transformed populations.
fn force_corrected<const D: usize, const Q: usize>(
    mut f: [f64; Q],
    velocities: [[f64; D]; Q],
    weights: [f64; Q],
    u: [f64; D],
    force: [f64; D],
) -> [f64; Q] {
    for q in 0..Q {
        let eu: f64 = velocities[q].iter().zip(u).map(|(e, u)| e * u).sum();
        let raw: f64 = (0..D)
            .map(|d| (3.0 * (velocities[q][d] - u[d]) + 9.0 * eu * velocities[q][d]) * force[d])
            .sum();
        f[q] += 0.5 * weights[q] * raw;
    }
    f
}

trait StressSlab<const D: usize, const Q: usize>: Slab<D, Q> {
    fn shear_numerator(&self, cell: usize) -> f64;
}
impl StressSlab<2, 5> for ThermalLbm {
    fn shear_numerator(&self, cell: usize) -> f64 {
        let m = self.grid.moments(cell);
        let force = std::array::from_fn(|d| self.grid.g[d].mul_add(m.rho, self.grid.fext[cell][d]));
        let corrected = force_corrected(
            self.grid.f[cell],
            crate::E.map(|e| [f64::from(e.0), f64::from(e.1)]),
            crate::W,
            m.u,
            force,
        );
        crate::core2::shear_rate(
            &corrected,
            &crate::equilibrium(m.rho, m.u[0], m.u[1]),
            m.rho,
            1.0,
        )
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
