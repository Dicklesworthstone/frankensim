//! Source-density-owned external hydraulic loads for resolved Stokes flow.

use super::{MaterialStokesBoundary, MaterialStokesSolution, with_physical_boundaries};
use crate::material::{MaterialFluxError, MaterialFluxSystem};
use crate::ns::boundary::{StokesAssembly, StokesResistancePort, StokesSolution};
use fs_material::fluid_transport::FluidTransport;

/// Physical reservoir and external resistance attached to one integral port.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialStokesResistancePort {
    /// Reservoir pressure increment relative to the original total traction [Pa].
    pub reservoir_pressure_pa: f64,
    /// Finite nonnegative pressure/flow-per-depth resistance [Pa s/m²].
    /// p_port = p_reservoir + R*Q_outward. For a 3-D network resistance
    /// in Pa s/m³, multiply by the modeled out-of-plane depth in metres.
    pub resistance_pa_s_per_m2: f64,
}

/// Physical laminar flow with all external resistance laws satisfied together.
pub struct MaterialStokesResistanceSolution<S> {
    /// Resolved field retaining its original material sources and reactions.
    pub flow: S,
    /// Uniform pressure increments added to original port traction [Pa].
    pub pressure_offsets_pa: Vec<f64>,
    /// Signed outward volume flow per unit out-of-plane depth [m²/s].
    pub outward_volume_rates_m2_s: Vec<f64>,
    /// Each external resistance's R*Q² [W/m], separate from bulk fluid heat.
    pub resistance_dissipation_w_per_m: Vec<f64>,
    /// Largest accepted |p_port - p_reservoir - R*Q_outward| [Pa].
    pub maximum_pressure_residual_pa: f64,
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Resolve steady laminar Stokes flow coupled to external resistances.
    ///
    /// Force density is N/m³, original traction and reservoir increments are
    /// Pa; resistance divides pressure by volume flow per unit depth. Source
    /// dynamic viscosity and the SAME reference density determine the flow,
    /// pressure, mechanical reactions and external R*Q² loss. Port velocity
    /// profiles remain unknown. Signed flow reversal is preserved.
    ///
    /// Every traction port may be controlled; reservoir laws select pressure.
    /// The law is linear hydraulic resistance, without external inertia,
    /// compliance, turbulence or automatic thermal deposition. It does not
    /// update the material's thermodynamic pressure or temperature.
    ///
    /// # Errors
    /// Preserves material, boundary, native equation and resistance-law refusals.
    /// No partial field or source mutation is published.
    pub fn solve_stokes_mixed_resistance_ports(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[MaterialStokesResistancePort],
        absolute_pressure_tolerance_pa: f64,
    ) -> Result<MaterialStokesResistanceSolution<MaterialStokesSolution>, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        let assembly = with_physical_boundaries(boundaries, density, |primitive| {
            self.system.assemble_stokes_mixed(
                self.params,
                &|point| force_density(point).map(|f| f / density),
                primitive,
            )
        })
        .map_err(MaterialFluxError::Boundary)?;
        solve_resistance_assembly(
            assembly,
            density,
            ports,
            loads,
            absolute_pressure_tolerance_pa,
            |flow, _| self.material_stokes_solution(flow),
        )
    }

    /// Advance source-owned unsteady Stokes with endpoint reservoir resistances.
    ///
    /// Uses native variable-step BDF2; `older_history=None` starts with BDF1.
    /// Every endpoint pressure is solved together with the inertial flow and
    /// original source viscosity, rather than imposed from lagged flow. Units,
    /// datum and external-model scope match the steady resistance entry.
    /// External heat is reported separately and no temperature is advanced.
    ///
    /// # Errors
    /// Preserves immutable histories, native interval/ratio admission, material
    /// and boundary refusals, and both continuum and port-law convergence.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_stokes_mixed_variable_step_resistance_ports(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[MaterialStokesResistancePort],
        previous_velocity: &[f64],
        older_history: Option<(&[f64], f64)>,
        timestep_s: f64,
        absolute_pressure_tolerance_pa: f64,
    ) -> Result<MaterialStokesResistanceSolution<MaterialStokesSolution>, MaterialFluxError> {
        let assembly = self.assemble_material_stokes_bdf2(
            force_density,
            boundaries,
            previous_velocity,
            older_history,
            timestep_s,
        )?;
        solve_resistance_assembly(
            assembly,
            self.transport.density_kg_m3(),
            ports,
            loads,
            absolute_pressure_tolerance_pa,
            |flow, _| self.material_stokes_solution(flow),
        )
    }
}

pub(super) fn solve_resistance_assembly<S>(
    assembly: StokesAssembly<'_, '_>,
    density: f64,
    ports: &[&[usize]],
    loads: &[MaterialStokesResistancePort],
    tolerance_pa: f64,
    recover: impl FnOnce(StokesSolution, f64) -> Result<S, MaterialFluxError>,
) -> Result<MaterialStokesResistanceSolution<S>, MaterialFluxError> {
    // Validate physical values before normalization: a negative or positive
    // resistance must not underflow to zero and become a different model.
    let scaled = |value: f64| {
        let result = value / density;
        value.is_finite() && result.is_finite() && (value == 0.0 || result != 0.0)
    };
    if loads.len() != ports.len()
        || tolerance_pa < 0.0
        || !scaled(tolerance_pa)
        || loads.iter().any(|load| {
            load.resistance_pa_s_per_m2 < 0.0
                || !scaled(load.resistance_pa_s_per_m2)
                || !scaled(load.reservoir_pressure_pa)
        })
    {
        return Err(MaterialFluxError::Invalid(
            "invalid or unrepresentable physical hydraulic loads or tolerance",
        ));
    }
    let native: Vec<_> = loads
        .iter()
        .map(|load| StokesResistancePort {
            reservoir_pressure_m2_s2: load.reservoir_pressure_pa / density,
            resistance_s_inv: load.resistance_pa_s_per_m2 / density,
        })
        .collect();
    let prepared = assembly
        .prepare_pressure_response(ports)
        .map_err(MaterialFluxError::Boundary)?;
    let solved = prepared
        .solve_resistance_ports(&native, tolerance_pa / density)
        .map_err(MaterialFluxError::Boundary)?;
    let kinetic_energy = density
        * prepared
            .kinetic_energy(&solved.flow.field.x)
            .map_err(MaterialFluxError::Boundary)?;
    let pressure_offsets_pa: Vec<_> = solved
        .pressure_offsets_m2_s2
        .iter()
        .map(|p| density * p)
        .collect();
    let resistance_dissipation_w_per_m: Vec<_> = solved
        .resistance_dissipation_m4_s3
        .iter()
        .map(|loss| density * loss)
        .collect();
    let maximum_pressure_residual_pa = density * solved.maximum_pressure_residual_m2_s2;
    if pressure_offsets_pa
        .iter()
        .chain(&resistance_dissipation_w_per_m)
        .chain([&kinetic_energy, &maximum_pressure_residual_pa])
        .any(|value| !value.is_finite())
    {
        return Err(MaterialFluxError::Invalid(
            "unrepresentable physical hydraulic result",
        ));
    }
    Ok(MaterialStokesResistanceSolution {
        flow: recover(solved.flow, kinetic_energy)?,
        pressure_offsets_pa,
        outward_volume_rates_m2_s: solved.outward_volume_rates_m2_s,
        resistance_dissipation_w_per_m,
        maximum_pressure_residual_pa,
    })
}
