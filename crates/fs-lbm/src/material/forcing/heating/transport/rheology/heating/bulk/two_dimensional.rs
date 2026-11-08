//! Physical planar shear/bulk response on the original D2Q9 thermal step.
//!
//! With D_zz=0, the in-plane stress is 2 mu dev_2(D) + (zeta_3+mu/3) theta I.
//! Therefore tau_trace = 1/2 + 3 (zeta_3+mu/3) dt/(rho dx²), while shear keeps
//! tau_s = 1/2 + 3 mu dt/(rho dx²). This is the dimensional reduction of
//! Dubois & Lallemand, Eq. (17), and its dimension-dependent viscous split:
//! <https://www.imo.universite-paris-saclay.fr/~francois.dubois/travaux/reseau/reseau-documents-2020s/POF24-AR-16191-20mars2025.pdf>.
//! Physical heat remains 2 mu dev_3(D):dev_3(D) + zeta_3 theta².

use super::super::super::{
    LocalRheologyError, sample_temperatures_with, solve_local, temperatures,
};
use super::super::{
    DilationRate, FluidTransport, FluidTransportState, HeatSourceLimits, MaterialError,
    MaterialSlab, MaterialTransportError, RheologyMaterial, RheologyPoint, RheologySolve,
    ThermalLbm, TransportSample, ViscousHeating, ViscousRheologyStep,
    advance_with_sampler_sources_and_candidate, candidate_with_prepared_input, poll, source_inputs,
};
use super::BulkViscousMaterial;
use crate::core2::{CollisionModel2, material_strain_components2};
use crate::material::{MaterialLbmError, plan_scaling_from_material};
use fs_exec::CancelGate;

fn trace_tau(
    mapping: &MaterialSlab<ThermalLbm>,
    shear_tau: f64,
    bulk_viscosity_pa_s: f64,
) -> Result<f64, &'static str> {
    let units = mapping.force_units;
    let bulk_lattice = ((bulk_viscosity_pa_s / units.transport.density_kg_m3())
        / units.velocity_scale)
        / units.cell_size_m;
    let bulk_excess = 3.0 * bulk_lattice;
    let tau = 0.5 + ((shear_tau - 0.5) / 3.0 + bulk_excess);
    if !bulk_viscosity_pa_s.is_finite()
        || bulk_viscosity_pa_s < 0.0
        || !bulk_excess.is_finite()
        || (bulk_viscosity_pa_s > 0.0 && 0.5 + bulk_excess <= 0.5)
        || !shear_tau.is_finite()
        || shear_tau <= 0.5
        || !tau.is_finite()
        || tau <= 0.5
    {
        return Err("physical planar bulk viscosity or trace relaxation is unrepresentable");
    }
    Ok(tau)
}

fn norm(a: f64, b: f64) -> f64 {
    let scale = a.max(b);
    if scale == 0.0 {
        0.0
    } else {
        let x = a / scale;
        let y = b / scale;
        scale * (x * x + y * y).sqrt()
    }
}

#[allow(clippy::too_many_arguments)]
fn sample<F: FluidTransport>(
    mapping: &MaterialSlab<ThermalLbm>,
    state: &ThermalLbm,
    cells: usize,
    plane: usize,
    range: [f64; 2],
    solve: RheologySolve,
    gate: &CancelGate,
    resolve: &mut impl FnMut(
        RheologyPoint,
        &FluidTransportState,
    ) -> Result<BulkViscousMaterial<F>, MaterialError>,
) -> Result<Vec<TransportSample<RheologyMaterial<BulkViscousMaterial<F>>>>, MaterialTransportError>
{
    let temperature = temperatures::<ThermalLbm, 2, 5>(mapping, state, range, gate)?;
    let mut stresses = vec![[0.0; 2]; cells];
    for cell in plane..cells - plane {
        poll(gate)?;
        let m = state.grid.moments(cell);
        let force =
            core::array::from_fn(|d| state.grid.g[d].mul_add(m.rho, state.grid.fext[cell][d]));
        // Unit trace tau retains A=gamma_dev2*tau_s and
        // B=sqrt(2/3)*|theta|*tau_trace in the original Guo force frame.
        stresses[cell] =
            material_strain_components2(&state.grid.f[cell], m.rho, m.u, force, Some(1.0));
        if stresses[cell].iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(MaterialTransportError::Invalid(
                "unrepresentable planar shear/bulk stress",
            ));
        }
    }
    let units = mapping.force_units;
    sample_temperatures_with::<ThermalLbm, RheologyMaterial<BulkViscousMaterial<F>>, 2, 5>(
        mapping,
        &temperature,
        cells,
        plane,
        range,
        gate,
        &mut |cell, temperature_k, retained| {
            let [a, b] = stresses[cell];
            // Normalize before forming N or inverse relaxation, avoiding N/g
            // overflow and squared-component underflow at tiny physical stress.
            let scale = a.max(b);
            let (numerator, weights) = if scale == 0.0 {
                (0.0, [1.0, 0.0])
            } else {
                let x = a / scale;
                let y = (b / scale) * core::f64::consts::FRAC_1_SQRT_2;
                let n = norm(x, y);
                (scale * n, [x / n, y / n])
            };
            if !numerator.is_finite() || (scale != 0.0 && numerator == 0.0) {
                return Err(MaterialTransportError::Invalid(
                    "unrepresentable planar strain invariant",
                ));
            }
            let error = |source| MaterialTransportError::State {
                cell,
                temperature_k,
                source,
            };
            let root = solve_local(
                numerator,
                units.timestep_s,
                state.grid.tau[cell],
                solve,
                || poll(gate),
                |rate| {
                    if numerator != 0.0 && rate == 0.0 {
                        return Err(MaterialTransportError::Invalid(
                            "physical planar shear rate is unrepresentable",
                        ));
                    }
                    let point = RheologyPoint {
                        cell,
                        temperature_k,
                        shear_rate_per_s: rate,
                    };
                    let source = resolve(point, retained).map_err(|e| error(e.into()))?;
                    let conversion = plan_scaling_from_material(
                        &source.transport,
                        units.cell_size_m,
                        units.timestep_s,
                        0.0,
                    )
                    .map_err(error)?;
                    let shear_tau = conversion.flow().tau;
                    let bulk_tau = trace_tau(mapping, shear_tau, source.bulk_viscosity_pa_s)
                        .map_err(|what| error(MaterialLbmError::Invalid(what)))?;
                    // gamma*dt = N/lambda. Re-evaluating the same source gives
                    // 1/lambda_next = hypot(w_s/tau_s, w_b/tau_trace).
                    // Trace-only and shear-only flow reduce to their own tau.
                    let equivalent_tau = if numerator == 0.0 {
                        shear_tau
                    } else {
                        1.0 / norm(weights[0] / shear_tau, weights[1] / bulk_tau)
                    };
                    if !equivalent_tau.is_finite() || equivalent_tau <= 0.5 {
                        return Err(MaterialTransportError::Invalid(
                            "unrepresentable equivalent planar relaxation",
                        ));
                    }
                    Ok(((source, bulk_tau), equivalent_tau))
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
            let (source, bulk_tau) = root.value;
            let dilation = (b / bulk_tau) / units.timestep_s;
            if !dilation.is_finite() || (b != 0.0 && dilation == 0.0) {
                return Err(MaterialTransportError::Invalid(
                    "physical planar dilation is unrepresentable",
                ));
            }
            Ok(RheologyMaterial {
                transport: source.transport,
                source,
                point: RheologyPoint {
                    cell,
                    temperature_k,
                    shear_rate_per_s: root.shear_rate_per_s,
                },
                iterations: root.iterations,
                relaxation_defect: root.relaxation_defect,
                dilation_rate_per_s: Some(DilationRate::Independent(dilation)),
            })
        },
    )
}

impl MaterialSlab<ThermalLbm> {
    /// Resolve physical planar mu(T,gamma) and zeta_3(T,gamma) together, and
    /// convert their original stress power into heat on the native transaction.
    ///
    /// gamma=sqrt(2 dev_3(D):dev_3(D)), with D_zz=0. The actual local trace
    /// coefficient is zeta_2=zeta_3+mu/3, so both sourced coefficients enter
    /// reconstruction of planar extension. Let A be the in-plane deviatoric
    /// stress numerator, B the sqrt(2/3) trace numerator, and N=hypot(A,B/sqrt(2)).
    /// The existing bounded scalar solver uses gamma=N/(lambda*dt), where
    /// 1/lambda=hypot((A/N)/tau_s, (B/(sqrt(2)*N))/tau_trace).
    /// Its relaxation_defect measures this dimensionless equivalent lambda.
    /// Zero stress queries gamma=0; zeta_3=0 retains physical Stokes viscosity
    /// through zeta_2=mu/3. Only the actual source-derived shear tau
    /// and trace tau are installed in collision; lambda is never installed.
    ///
    /// The same accepted source snapshot supplies both coefficients and
    /// mu*gamma²+zeta_3*(div u)², without double-counting the planar trace.
    /// Initial coefficients and power are held over one explicit step;
    /// endpoint source queries admit the resulting candidate before publication.
    /// The accepted grid retains the applied local trace field. Repeat this
    /// method for material feedback; no molecular memory or implicit heat solve
    /// is inferred. Density/Cp, low-Mach and source-domain restrictions remain.
    /// All source, convergence, conversion, heat-balance and cancellation
    /// refusals preserve the original populations, loads and trace field.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_local_bulk_rheology_and_viscous_heating<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 2]],
        depth_m: f64,
        limits: HeatSourceLimits,
        solve: RheologySolve,
        gate: &CancelGate,
        mut resolve: impl FnMut(
            RheologyPoint,
            &FluidTransportState,
        ) -> Result<BulkViscousMaterial<F>, MaterialError>,
    ) -> Result<ViscousRheologyStep<BulkViscousMaterial<F>>, MaterialTransportError> {
        poll(gate)?;
        solve.validate()?;
        let (transport, viscous) = advance_with_sampler_sources_and_candidate::<
            ThermalLbm,
            RheologyMaterial<BulkViscousMaterial<F>>,
            ViscousHeating,
            2,
            5,
        >(
            self,
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            gate,
            |mapping, state, cells, plane| {
                sample(
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
            |mapping, original, applied| {
                source_inputs(
                    mapping,
                    power_w_m3,
                    forces_n_m3,
                    depth_m,
                    limits,
                    original,
                    applied,
                    |material| Some(material.source().bulk_viscosity_pa_s),
                    gate,
                )
            },
            |mapping, input, applied| {
                let mut taus = vec![1.0; mapping.lattice.grid.f.len()];
                for value in applied {
                    poll(gate)?;
                    taus[value.cell] = trace_tau(
                        mapping,
                        value.momentum_relaxation,
                        value.material.source().bulk_viscosity_pa_s,
                    )
                    .map_err(MaterialTransportError::Invalid)?;
                }
                candidate_with_prepared_input::<ThermalLbm, _, 2, 5>(
                    mapping,
                    input,
                    applied,
                    gate,
                    move |trial| {
                        trial.grid.collision = CollisionModel2::CentralMomentBulk { bulk_tau: 1.0 };
                        trial.grid.central_moment_bulk_tau = Some(taus);
                    },
                )
            },
        )?;
        Ok(ViscousRheologyStep { transport, viscous })
    }

    /// Constant explicitly authored physical zeta_3 with source-local shear
    /// viscosity. The same planar coupled root and local trace field are used:
    /// even constant zeta_3 gives spatially varying zeta_2 when mu varies.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_rheology_and_bulk_viscous_heating<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 2]],
        depth_m: f64,
        bulk_viscosity_pa_s: f64,
        limits: HeatSourceLimits,
        solve: RheologySolve,
        gate: &CancelGate,
        mut resolve: impl FnMut(RheologyPoint, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<ViscousRheologyStep<BulkViscousMaterial<F>>, MaterialTransportError> {
        poll(gate)?;
        trace_tau(self, 1.0, bulk_viscosity_pa_s).map_err(MaterialTransportError::Invalid)?;
        self.step_with_local_bulk_rheology_and_viscous_heating(
            power_w_m3,
            forces_n_m3,
            depth_m,
            limits,
            solve,
            gate,
            |point, retained| {
                BulkViscousMaterial::try_new(resolve(point, retained)?, bulk_viscosity_pa_s)
            },
        )
    }
}
