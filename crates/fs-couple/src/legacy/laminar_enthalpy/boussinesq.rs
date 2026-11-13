//! Implicit source-density buoyancy with original reference caloric inventory.
//!
//! The SAME endpoint material supplies density, viscosity and conductivity.
//! Only `(rho_source - rho0) * gravity` enters momentum; the hydrostatic
//! reference load has already been removed from pressure and traction. Native
//! mass and enthalpy storage keep their original reference inventory. Buoyancy
//! work is mechanical work, not an additional heat source.

use super::*;
use fs_flux::bdm::eval_basis;
pub use fs_lbm::material::transport::boussinesq::{BoussinesqConfig, BoussinesqMaterial};
use fs_math::dd::Dd;

/// The shared ORIGINAL Boussinesq projector retains the frozen thermal source.
/// Its `material().source()` is the original source object, and
/// `material().thermal_state()` is its unprojected coherent thermal view.
pub type LaminarBoussinesqProjection<G> = BoussinesqMaterial<LaminarThermalMaterial<G>>;

/// Joint native thermal/flow endpoint retaining the original bounded projection.
pub type LaminarBoussinesqCandidate<'m, G, T> =
    LaminarEnthalpyCandidate<'m, LaminarBoussinesqProjection<G>, T>;

/// Actual density-buoyancy loads and work, distinct from the thermal account.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminarBuoyancyReport {
    pub config: BoussinesqConfig<2>,
    pub reference_density_kg_m3: f64,
    pub maximum_relative_density_departure: f64,
    pub maximum_relative_specific_heat_departure: f64,
    /// Accepted original cell forces [N/m³], excluding the held external force.
    pub cell_force_density_n_m3: Vec<[f64; 2]>,
    /// Integral of accepted buoyancy force dot original affine BDM velocity.
    /// Includes the declared extrusion depth, hence watts rather than W/m.
    pub endpoint_power_w: f64,
    /// Backward-Euler endpoint power times the actual interval. This is not
    /// deposited heat, gravitational-potential change or a passivity proof.
    pub endpoint_work_j: f64,
}

impl<'s, 'mesh, 'source, F: FluidTransport + ?Sized> LaminarEnthalpy<'s, 'mesh, 'source, F> {
    /// Joint BDF1 natural convection, rheology and original enthalpy advection.
    ///
    /// Every thermal/strain query resolves one coherent original source and
    /// admits its density and Cp departures through `BoussinesqMaterial`.
    /// The original rho0 inventory/inertia and reference transport Cp remain;
    /// viscosity and conductivity are the SAME source's actual local values.
    /// The supplied original caloric law still owns enthalpy and temperature;
    /// no Cp integration or EOS fit is inferred from this projection. The
    /// caller must associate it with the source and choose physically small
    /// density/Cp bounds. Those bounds are not an approximation-error estimate.
    ///
    /// `input.force_density` supplies additional NONBUOYANT force only. Do not
    /// include rho0*g, a separate beta*(T-T0) load, or this density difference
    /// again. Open traction and returned pressure use the reduced variable
    /// `p_reduced = p_physical - rho0*g.(x-x_reference)`; any desired physical
    /// pressure/reaction reconstruction needs an explicit hydrostatic datum.
    /// This owner never invents that datum or modifies supplied boundary loads.
    ///
    /// The native final momentum residual includes both endpoint viscosity
    /// and buoyancy. Sign-split enthalpy transport, physical strain/slip heat,
    /// finite original thermal acceptance and both-history atomicity are the
    /// same as `advance`. The projected thermal view carries no inferred
    /// expansion coefficient; the original source view remains retained.
    /// No buoyancy work, pressure work or numerical loss is added as heat.
    /// This remains fixed-volume, first-order Boussinesq flow, not variable-mass
    /// transport, a compressible EOS or whole-system energy certification.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_boussinesq<'m, G: ThermalFluidTransport, T, E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        input: LaminarEnthalpyInput<'_>,
        config: BoussinesqConfig<2>,
        controls: LaminarEnthalpyControls,
        energy_guess_j: &[f64],
        resolve_material: impl Fn(MaterialFluxRheologyPoint, f64, &F) -> Result<G, MaterialError>,
        exchange: impl Fn(
            &Cx<'_>,
            &DepositionSimulation<'m>,
            &MaterialNavierStokesRheologyStep<
                LaminarThermalMaterial<LaminarBoussinesqProjection<G>>,
            >,
        ) -> Result<LaminarHeatExchange<T>, E>,
    ) -> Result<LaminarBoussinesqCandidate<'m, G, T>, LaminarEnthalpyError<E>> {
        checkpoint(cx)?;
        config.validate().map_err(|error| {
            LaminarEnthalpyError::Momentum(MaterialFluxError::Material(MaterialError::State {
                what: error.to_string(),
            }))
        })?;
        self.advance_with_flow(
            cx,
            thermal,
            target_s,
            input,
            controls,
            energy_guess_j,
            |temperatures, duration_s| {
                self.resolve_boussinesq_flow(
                    cx,
                    input,
                    config,
                    controls,
                    temperatures,
                    duration_s,
                    &resolve_material,
                )
            },
            exchange,
        )
    }
    // Ordinary and driven endpoints retain the SAME original projection path.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(super) fn resolve_boussinesq_flow<G: ThermalFluidTransport>(
        &self,
        cx: &Cx<'_>,
        input: LaminarEnthalpyInput<'_>,
        config: BoussinesqConfig<2>,
        controls: LaminarEnthalpyControls,
        temperatures: &[f64],
        duration_s: f64,
        resolve_material: &impl Fn(MaterialFluxRheologyPoint, f64, &F) -> Result<G, MaterialError>,
    ) -> Result<
        (
            MaterialNavierStokesRheologyStep<
                LaminarThermalMaterial<LaminarBoussinesqProjection<G>>,
            >,
            Option<LaminarBuoyancyReport>,
        ),
        MaterialFluxError,
    > {
        let (flow, forces) = self
            .flow
            .bdf1_navier_stokes_mixed_rheology_with_cell_forces(
                input.force_density,
                input.boundaries,
                input.history.velocity_m_s,
                duration_s,
                controls.momentum,
                |point, reference| {
                    if cx.checkpoint().is_err() {
                        return Err(material_error(
                            "coupled buoyancy/enthalpy operation cancelled",
                        ));
                    }
                    let temperature = temperatures[point.cell];
                    let source = resolve_material(point, temperature, self.flow.material())?;
                    let original = thermal_source(source, temperature)?;
                    let projected = BoussinesqMaterial::try_new(original, reference, config)?;
                    let departure =
                        projected.resolved_transport().density_kg_m3() - reference.density_kg_m3();
                    let force = config.gravity_m_s2.map(|g| departure * g);
                    if force.iter().enumerate().any(|(axis, f)| {
                        !f.is_finite()
                            || (*f == 0.0 && departure != 0.0 && config.gravity_m_s2[axis] != 0.0)
                    }) {
                        return Err(material_error("source-density buoyancy is unrepresentable"));
                    }
                    let state = ThermalFluidTransportState::try_new(
                        *projected.reference_storage_transport(),
                        temperature,
                        None,
                    )?;
                    Ok((
                        LaminarThermalMaterial {
                            source: projected,
                            state,
                        },
                        force,
                    ))
                },
            )?;
        let report = buoyancy_report(
            cx,
            self.system,
            &flow,
            forces,
            self.flow.transport(),
            config,
            self.binding.depth_m(),
            duration_s,
        )?;
        Ok((flow, Some(report)))
    }
}

#[allow(clippy::too_many_arguments)]
fn buoyancy_report<G>(
    cx: &Cx<'_>,
    system: &FluxSystem<'_>,
    flow: &MaterialNavierStokesRheologyStep<LaminarThermalMaterial<LaminarBoussinesqProjection<G>>>,
    forces: Vec<[f64; 2]>,
    reference: &FluidTransportState,
    config: BoussinesqConfig<2>,
    depth_m: f64,
    duration_s: f64,
) -> Result<LaminarBuoyancyReport, MaterialFluxError> {
    if forces.len() != system.mesh.tris.len() || flow.samples().len() != forces.len() {
        return Err(MaterialFluxError::Invalid(
            "incomplete accepted buoyancy field",
        ));
    }
    let mut maximum_relative_density_departure = 0.0_f64;
    let mut maximum_relative_specific_heat_departure = 0.0_f64;
    let mut power = Dd::from_f64(0.0);
    let velocity = flow.step.solution.velocity_dofs();
    for (cell, (sample, force)) in flow.samples().iter().zip(&forces).enumerate() {
        if cx.checkpoint().is_err() {
            return Err(MaterialFluxError::Invalid(
                "buoyancy work evaluation cancelled",
            ));
        }
        let resolved = sample.material().source().resolved_transport();
        maximum_relative_density_departure = maximum_relative_density_departure.max(
            ((resolved.density_kg_m3() - reference.density_kg_m3()) / reference.density_kg_m3())
                .abs(),
        );
        maximum_relative_specific_heat_departure = maximum_relative_specific_heat_departure.max(
            ((resolved.specific_heat_j_kg_k() - reference.specific_heat_j_kg_k())
                / reference.specific_heat_j_kg_k())
            .abs(),
        );
        // Cell force is constant and the original BDM1 velocity is affine:
        // their volume integral is exactly area*depth times the centroid value.
        // Evaluate the ORIGINAL basis directly, with no global point lookup.
        let mut mean = [Dd::from_f64(0.0); 2];
        for i in 0..6 {
            let (edge, _) = system.mesh.tri_edges[cell][i / 2];
            let value = eval_basis(&system.bases[cell], i, system.mesh.centroids[cell]);
            let dof = Dd::from_f64(velocity[2 * edge + i % 2]);
            for axis in 0..2 {
                mean[axis] = mean[axis] + dof * Dd::from_f64(value[axis]);
            }
        }
        let work = Dd::from_f64(force[0]) * mean[0] + Dd::from_f64(force[1]) * mean[1];
        power = power + work * Dd::from_f64(system.mesh.areas[cell]) * Dd::from_f64(depth_m);
    }
    let endpoint_power_w = power.to_f64();
    let endpoint_work_j = (power * Dd::from_f64(duration_s)).to_f64();
    if !endpoint_power_w.is_finite() || !endpoint_work_j.is_finite() {
        return Err(MaterialFluxError::Invalid(
            "buoyancy work is unrepresentable",
        ));
    }
    Ok(LaminarBuoyancyReport {
        config,
        reference_density_kg_m3: reference.density_kg_m3(),
        maximum_relative_density_departure,
        maximum_relative_specific_heat_departure,
        cell_force_density_n_m3: forces,
        endpoint_power_w,
        endpoint_work_j,
    })
}
