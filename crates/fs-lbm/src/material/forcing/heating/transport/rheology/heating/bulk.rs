//! Local physical shear/bulk viscosity from one actual material query.

mod two_dimensional;

use super::super::{LocalRheologyError, sample_temperatures_with, solve_local, temperatures};
use super::{
    DilationRate, FluidTransport, FluidTransportState, HeatSourceLimits, MaterialError, MaterialSlab,
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
    /// Nonnegative finite physical bulk viscosity is admitted here. The
    /// receiving collision checks its own representability: planar flow can
    /// retain zeta_3=0 through mu/3, while the 3-D trace collision requires
    /// strictly positive bulk relaxation excess.
    pub fn try_new(source: F, bulk_viscosity_pa_s: f64) -> Result<Self, MaterialError> {
        if !bulk_viscosity_pa_s.is_finite() || bulk_viscosity_pa_s < 0.0 {
            return Err(MaterialError::Parameters {
                what: "local physical bulk viscosity must be nonnegative and finite".into(),
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
    let mut stresses = vec![[0.0; 2]; cells];
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
        stresses[cell] = [a, b];
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
            let [a, b] = stresses[cell];
            let previous_shear = state.grid.tau[cell];
            let error = |source| MaterialTransportError::State {
                cell,
                temperature_k,
                source,
            };
            let root = solve_local(
                a,
                units.timestep_s,
                previous_shear,
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
                    let shear_tau = scale.flow().tau;
                    Ok(((source, bulk_tau), shear_tau))
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
                // Bulk stress still supplies heat, using the independently
                // sourced coefficient at the accepted deviatoric shear query.
                dilation_rate_per_s: Some(DilationRate::Independent(
                    (b / bulk_tau) / units.timestep_s,
                )),
            })
        },
    )
}

impl MaterialSlab<ThermalLbm3> {
    /// Advance with source-local `mu(T,gamma)` AND `zeta(T,gamma)` and their
    /// physical viscous heat, where `gamma=sqrt(2 dev(D):dev(D))` in inverse seconds.
    ///
    /// Each callback returns one original transport response and explicit bulk
    /// viscosity for the same actual shear query. The bounded scalar solve
    /// resolves deviatoric strain against the material's shear relaxation;
    /// `relaxation_defect` measures that shear-tau residual. Dilation does not
    /// alter the query, but its separately relaxed stress supplies bulk heat.
    /// Collision retains the separate physical shear and bulk relaxations.
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
