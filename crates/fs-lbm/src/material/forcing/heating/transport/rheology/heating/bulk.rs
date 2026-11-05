//! Local physical shear/bulk viscosity from one actual material query.

use super::super::{LocalRheologyError, sample_temperatures_with, solve_local, temperatures};
use super::{
    FluidTransport, FluidTransportState, HeatSourceLimits, MaterialError, MaterialSlab,
    MaterialTransportError, RheologyMaterial, RheologyPoint, RheologySolve, ThermalLbm3,
    TransportSample, ViscousHeating, ViscousRheologyStep,
    advance_with_sampler_sources_and_candidate, candidate_with_prepared_input, physical_bulk_tau,
    poll, source_inputs,
};
use crate::d3q19::{BulkRelaxation3, CollisionModel3, rheology_stress_components3};
use crate::material::{MaterialLbmError, plan_scaling_from_material};
use fs_exec::CancelGate;

/// One original transport response and its explicitly declared physical bulk
/// viscosity. Neither a bulk/shear ratio nor a numerical damping law is inferred.
#[derive(Debug, Clone, PartialEq)]
pub struct BulkViscousMaterial<F> {
    source: F,
    transport: FluidTransportState,
    bulk_viscosity_pa_s: f64,
}
impl<F: FluidTransport> BulkViscousMaterial<F> {
    /// Read the original transport owner once and retain that exact snapshot.
    /// Positive finite bulk viscosity is required by the receiving collision
    /// family; its SI-to-lattice representability is checked at the actual query.
    pub fn try_new(source: F, bulk_viscosity_pa_s: f64) -> Result<Self, MaterialError> {
        if !bulk_viscosity_pa_s.is_finite() || bulk_viscosity_pa_s <= 0.0 {
            return Err(MaterialError::Parameters {
                what: "local physical bulk viscosity must be positive and finite".into(),
            });
        }
        let transport = source.transport_state()?;
        Ok(Self {
            source,
            transport,
            bulk_viscosity_pa_s,
        })
    }
}
impl<F> BulkViscousMaterial<F> {
    /// Original source, including its material query and property receipts.
    pub const fn source(&self) -> &F {
        &self.source
    }
    /// Physical coefficient returned at this exact temperature/strain query.
    pub const fn bulk_viscosity_pa_s(&self) -> f64 {
        self.bulk_viscosity_pa_s
    }
}
impl<F> FluidTransport for BulkViscousMaterial<F> {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.transport)
    }
}

fn norm([a, b]: [f64; 2]) -> f64 {
    let scale = a.max(b);
    if scale == 0.0 {
        return 0.0;
    }
    let a = a / scale;
    let b = b / scale;
    scale * (a * a + b * b).sqrt()
}

// With N=hypot(Adev,Abulk), gamma=N/(tau_eff dt), where
// 1/tau_eff=hypot((Adev/N)/tau_s,(Abulk/N)/tau_b). This is a scalar
// constitutive coordinate, never a substitute collision relaxation time.
// Scale by the smaller actual tau; pure modes keep their original arithmetic.
fn effective_tau(weights: [f64; 2], tau: [f64; 2]) -> f64 {
    if weights[1] == 0.0 {
        return tau[0];
    }
    if weights[0] == 0.0 {
        return tau[1];
    }
    let scale = tau[0].min(tau[1]);
    scale / norm([weights[0] * (scale / tau[0]), weights[1] * (scale / tau[1])])
}

#[allow(clippy::too_many_arguments)]
fn sample<F: FluidTransport>(
    mapping: &MaterialSlab<ThermalLbm3>,
    state: &ThermalLbm3,
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
    let temperature = temperatures::<ThermalLbm3, 3, 7>(mapping, state, range, gate)?;
    // Unit trace relaxation extracts the unrelaxed dilation numerator. Both
    // physical relaxation times will come from the same constitutive query.
    let unit = CollisionModel3::CentralMomentBulk {
        shear_rate: 1.0,
        bulk_rate: 1.0,
        higher_order_rate: 1.0,
    };
    let mut stresses = vec![[0.0; 3]; cells];
    for cell in plane..cells - plane {
        poll(gate)?;
        let (rho, velocity) = state.grid.moments3(cell);
        let [a, b] = rheology_stress_components3(
            &state.grid.f[cell],
            rho,
            velocity,
            state.grid.fext[cell],
            unit,
        );
        stresses[cell] = [a, b, norm([a, b])];
        if stresses[cell]
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(MaterialTransportError::Invalid(
                "unrepresentable shear/bulk stress",
            ));
        }
    }
    let units = mapping.force_units;
    sample_temperatures_with::<ThermalLbm3, RheologyMaterial<BulkViscousMaterial<F>>, 3, 7>(
        mapping,
        &temperature,
        cells,
        plane,
        range,
        gate,
        &mut |cell, temperature_k, retained| {
            let [a, b, magnitude] = stresses[cell];
            let weights = if magnitude == 0.0 {
                [1.0, 0.0]
            } else {
                [a / magnitude, b / magnitude]
            };
            let previous_shear = state.grid.tau[cell];
            let previous_bulk = match state.grid.collision_model_at(cell) {
                CollisionModel3::CentralMomentBulk { bulk_rate, .. } => 1.0 / bulk_rate,
                _ => previous_shear,
            };
            let error = |source| MaterialTransportError::State {
                cell,
                temperature_k,
                source,
            };
            let root = solve_local(
                magnitude,
                units.timestep_s,
                effective_tau(weights, [previous_shear, previous_bulk]),
                solve,
                || poll(gate),
                |rate| {
                    let point = RheologyPoint {
                        cell,
                        temperature_k,
                        shear_rate_per_s: rate,
                    };
                    let source = resolve(point, retained).map_err(|e| error(e.into()))?;
                    let scale = plan_scaling_from_material(
                        &source.transport,
                        units.cell_size_m,
                        units.timestep_s,
                        0.0,
                    )
                    .map_err(error)?;
                    let bulk_tau = physical_bulk_tau(mapping, source.bulk_viscosity_pa_s)
                        .map_err(|what| error(MaterialLbmError::Invalid(what)))?;
                    let tau = [scale.flow().tau, bulk_tau];
                    let effective = effective_tau(weights, tau);
                    if !effective.is_finite() || effective <= 0.0 {
                        return Err(error(MaterialLbmError::Invalid(
                            "unrepresentable coupled constitutive relaxation",
                        )));
                    }
                    Ok(((source, tau), effective))
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
            let (source, tau) = root.value;
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
                // The accepted physical rates drive BOTH collision and power.
                // Their norm agrees with the queried rate to the admitted
                // effective-tau residual, including when shear alone is constant.
                strain_components_per_s: Some([
                    (a / tau[0]) / units.timestep_s,
                    (b / tau[1]) / units.timestep_s,
                ]),
            })
        },
    )
}

impl MaterialSlab<ThermalLbm3> {
    /// Advance with source-local `mu(T,gamma)` AND `zeta(T,gamma)` and their
    /// physical viscous heat, where `gamma=sqrt(2 D:D)` in inverse seconds.
    ///
    /// Each callback returns one original transport response and explicit bulk
    /// viscosity for the same actual query. A single bounded scalar solve
    /// resolves both strain contributions against both material rates. Its
    /// `relaxation_defect` measures the effective constitutive coordinate;
    /// collision retains the separate physical shear and bulk relaxations.
    ///
    /// Initial local coefficients are held over the native thermal step. The
    /// same samples supply `2mu D_dev:D_dev+zeta(div u)^2` once. Endpoint source
    /// responses are checked before publication and retained without advancing
    /// again. The accepted grid retains the APPLIED bulk-rate field; repeating
    /// this call closes explicit temperature/shear/bulk feedback. Higher-order
    /// relaxation is uniform and supplies no heat. External power excludes
    /// these molecular losses.
    ///
    /// All source, conversion, constitutive, heat-balance and cancellation
    /// refusals preserve the complete accepted state and collision field.
    /// This remains low-Mach fixed-density/Cp sensible-heat transport, not
    /// molecular relaxation memory, compressible total-energy conservation or
    /// an implicit thermal solve. Zero/unrepresentable viscosity refuses.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_local_bulk_rheology_and_viscous_heating<F: FluidTransport>(
        &mut self,
        power_w_m3: &[f64],
        forces_n_m3: &[[f64; 3]],
        higher_order_rate: f64,
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
        if (CollisionModel3::CentralMomentBulk {
            shear_rate: 1.0,
            bulk_rate: 1.0,
            higher_order_rate,
        })
        .validate()
        .is_err()
        {
            return Err(MaterialTransportError::Invalid(
                "invalid higher-order central-moment relaxation",
            ));
        }
        let depth_m = self.force_units.cell_size_m;
        let (transport, viscous) = advance_with_sampler_sources_and_candidate::<
            ThermalLbm3,
            RheologyMaterial<BulkViscousMaterial<F>>,
            ViscousHeating,
            3,
            7,
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
                let mut rates = vec![1.0; mapping.lattice.grid.f.len()];
                for value in applied {
                    poll(gate)?;
                    rates[value.cell] = 1.0
                        / physical_bulk_tau(mapping, value.material.source().bulk_viscosity_pa_s)
                            .map_err(MaterialTransportError::Invalid)?;
                }
                candidate_with_prepared_input::<ThermalLbm3, _, 3, 7>(
                    mapping,
                    input,
                    applied,
                    gate,
                    move |trial| {
                        trial.grid.trt_magic = None;
                        trial.grid.central_moment_bulk = Some(BulkRelaxation3::Field {
                            rates,
                            higher_order_rate,
                        });
                    },
                )
            },
        )?;
        Ok(ViscousRheologyStep { transport, viscous })
    }
}
