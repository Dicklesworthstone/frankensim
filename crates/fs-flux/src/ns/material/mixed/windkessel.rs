//! One source-density normalization for fluid inertia and external RCR storage.

use super::{MaterialStokesBoundary, MaterialStokesSolution};
use crate::material::{MaterialFluxError, MaterialFluxSystem};
use crate::ns::boundary::{
    StokesAssembly, StokesSolution, StokesWindkesselBalance, StokesWindkesselHistory,
    StokesWindkesselPort, StokesWindkesselSolution, solve_windkessel_assembly,
};
use fs_material::fluid_transport::FluidTransport;

#[path = "windkessel/nonlinear.rs"]
mod nonlinear;

/// Physical RCR/Windkessel load with optional series inertance, per unit depth.
/// R, C and L are explicit external device properties; they are not inferred
/// from the fluid's transport viscosity or from an invented compressibility.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialStokesWindkesselPort {
    /// Endpoint reservoir pressure increment relative to original traction [Pa].
    pub reservoir_pressure_pa: f64,
    /// Nonnegative proximal resistance [Pa s/m²].
    pub proximal_resistance_pa_s_per_m2: f64,
    /// Strictly positive distal resistance [Pa s/m²].
    pub distal_resistance_pa_s_per_m2: f64,
    /// Strictly positive volume-per-depth compliance [m²/Pa].
    /// Divide a 3-D volume compliance [m³/Pa] by the modeled depth [m].
    pub compliance_m2_per_pa: f64,
    /// Nonnegative series inertance [Pa s²/m²]. Zero gives ordinary RCR.
    /// Multiply 3-D resistance/inertance [Pa s/m³, Pa s²/m³] by depth [m].
    pub inertance_pa_s2_per_m2: f64,
}

/// Native fluid velocity and external compliance at one physical time level.
#[derive(Clone, Copy)]
pub struct MaterialStokesWindkesselHistory<'a> {
    /// Original velocity moments [m/s], excluding pressure coefficients.
    pub velocity: &'a [f64],
    /// Compliance pressures relative to the same reservoirs [Pa], in port order.
    pub compliance_pressures_pa: &'a [f64],
}

/// Physical storage/work account for one accepted external hydraulic load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialStokesWindkesselBalance {
    /// Accepted compliance pressure relative to its reservoir [Pa].
    pub compliance_pressure_pa: f64,
    /// Stored compliance energy C*x²/2 [J/m].
    pub compliance_energy_j_per_m: f64,
    /// Stored series-inertance energy L*Q²/2 [J/m].
    pub inertance_energy_j_per_m: f64,
    /// Proximal resistor loss Rp*Q² [W/m].
    pub proximal_dissipation_w_per_m: f64,
    /// Distal resistor loss x²/Rd [W/m].
    pub distal_dissipation_w_per_m: f64,
    /// Reservoir pressure increment times outward flow [W/m], signed.
    pub reservoir_power_w_per_m: f64,
    /// Added port pressure times outward flow [W/m], entering the load.
    pub port_power_w_per_m: f64,
    /// Actual-interval change in stored external energy [W/m].
    pub stored_energy_change_rate_w_per_m: f64,
    /// Signed numerical BDF energy defect [W/m]. This is NOT physical heat.
    pub bdf_energy_defect_w_per_m: f64,
    /// Original port-power identity residual [W/m].
    pub power_balance_residual_w_per_m: f64,
    /// Signed pressure/inertance-law residual [Pa].
    pub pressure_law_residual_pa: f64,
    /// Signed compliance-continuity residual [m²/s].
    pub compliance_flow_residual_m2_s: f64,
}

/// Source-owned resolved flow and implicitly coupled physical hydraulic storage.
pub struct MaterialStokesWindkesselSolution<S> {
    /// Original material flow, source snapshots, continuum heat and reactions.
    pub flow: S,
    /// Added uniform port pressure increments [Pa].
    pub pressure_offsets_pa: Vec<f64>,
    /// Signed outward volume flow per unit depth [m²/s].
    pub outward_volume_rates_m2_s: Vec<f64>,
    /// External state, storage, physical resistor losses and numerical defect.
    pub balances: Vec<MaterialStokesWindkesselBalance>,
    /// Largest accepted original pressure/inertance-law defect [Pa].
    pub maximum_pressure_residual_pa: f64,
    /// Largest accepted original compliance-continuity defect [m²/s].
    pub maximum_flow_residual_m2_s: f64,
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Advance source-owned Stokes and RCR/series-inertance loads together.
    ///
    /// Fluid momentum, external Q inertia and compliance use the SAME admitted
    /// BDF1/startup or unequal-step BDF2 derivative. Flow histories are integrated
    /// from the supplied native velocity moments, never provided independently.
    /// Pressures and R,L divide by the single source density; C multiplies it.
    /// This preserves physical pressure, source viscosity, inertia and work.
    ///
    /// Reservoir increments are held over each interval. The external parameters,
    /// compliance reference, port ordering and fluid density remain fixed through
    /// supplied histories. Each controlled component needs an uncontrolled
    /// traction datum. Storage/losses remain separate from resolved fluid heat;
    /// no thermal state, thermodynamic pressure or phase model is advanced.
    ///
    /// # Errors
    /// Preserves native time/history, original equation and both dynamic-law
    /// checks; refuses nonpassive/unrepresentable SI inputs or outputs. Neither
    /// source owners nor histories are mutated, even if any port is refused.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_stokes_mixed_variable_step_windkessel_ports(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[MaterialStokesWindkesselPort],
        previous: MaterialStokesWindkesselHistory<'_>,
        older_history: Option<(MaterialStokesWindkesselHistory<'_>, f64)>,
        timestep_s: f64,
        absolute_pressure_tolerance_pa: f64,
        absolute_flow_tolerance_m2_s: f64,
    ) -> Result<MaterialStokesWindkesselSolution<MaterialStokesSolution>, MaterialFluxError> {
        let assembly = self.assemble_material_stokes_bdf2(
            force_density,
            boundaries,
            previous.velocity,
            older_history.map(|(history, dt)| (history.velocity, dt)),
            timestep_s,
        )?;
        solve_material_windkessel_assembly(
            assembly,
            self.transport.density_kg_m3(),
            ports,
            loads,
            previous,
            older_history,
            timestep_s,
            absolute_pressure_tolerance_pa,
            absolute_flow_tolerance_m2_s,
            |flow, _| self.material_stokes_solution(flow),
        )
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve_material_windkessel_assembly<S>(
    assembly: StokesAssembly<'_, '_>,
    density: f64,
    ports: &[&[usize]],
    loads: &[MaterialStokesWindkesselPort],
    previous: MaterialStokesWindkesselHistory<'_>,
    older_history: Option<(MaterialStokesWindkesselHistory<'_>, f64)>,
    timestep_s: f64,
    tolerance_pa: f64,
    tolerance_m2_s: f64,
    recover: impl FnOnce(StokesSolution, f64) -> Result<S, MaterialFluxError>,
) -> Result<MaterialStokesWindkesselSolution<S>, MaterialFluxError> {
    let solved = with_kinematic_windkessel(
        density,
        ports,
        loads,
        previous,
        older_history,
        tolerance_pa,
        tolerance_m2_s,
        |native, previous, older| {
            solve_windkessel_assembly(
                assembly,
                ports,
                native,
                previous,
                older,
                timestep_s,
                tolerance_pa / density,
                tolerance_m2_s,
            )
            .map_err(MaterialFluxError::Boundary)
        },
    )?;
    physical_windkessel(solved, density, recover)
}

#[allow(clippy::too_many_arguments)]
fn with_kinematic_windkessel<T>(
    density: f64,
    ports: &[&[usize]],
    loads: &[MaterialStokesWindkesselPort],
    previous: MaterialStokesWindkesselHistory<'_>,
    older_history: Option<(MaterialStokesWindkesselHistory<'_>, f64)>,
    tolerance_pa: f64,
    tolerance_m2_s: f64,
    solve: impl FnOnce(
        &[StokesWindkesselPort],
        StokesWindkesselHistory<'_>,
        Option<(StokesWindkesselHistory<'_>, f64)>,
    ) -> Result<T, MaterialFluxError>,
) -> Result<T, MaterialFluxError> {
    let normalized = |value: f64| {
        let result = value / density;
        value.is_finite() && result.is_finite() && (value == 0.0 || result != 0.0)
    };
    if !density.is_finite()
        || density <= 0.0
        || loads.len() != ports.len()
        || tolerance_pa < 0.0
        || !normalized(tolerance_pa)
        || !tolerance_m2_s.is_finite()
        || tolerance_m2_s < 0.0
        || loads.iter().any(|load| {
            load.proximal_resistance_pa_s_per_m2 < 0.0
                || load.distal_resistance_pa_s_per_m2 <= 0.0
                || load.inertance_pa_s2_per_m2 < 0.0
                || !normalized(load.reservoir_pressure_pa)
                || !normalized(load.proximal_resistance_pa_s_per_m2)
                || !normalized(load.distal_resistance_pa_s_per_m2)
                || !normalized(load.inertance_pa_s2_per_m2)
                || !(load.compliance_m2_per_pa.is_finite()
                    && load.compliance_m2_per_pa > 0.0
                    && (load.compliance_m2_per_pa * density).is_finite()
                    && load.compliance_m2_per_pa * density > 0.0)
        })
    {
        return Err(MaterialFluxError::Invalid(
            "invalid or unrepresentable physical RCR/inertance parameters",
        ));
    }
    let pressures =
        |history: MaterialStokesWindkesselHistory<'_>| -> Result<Vec<f64>, MaterialFluxError> {
            if history.compliance_pressures_pa.len() != ports.len()
                || history
                    .compliance_pressures_pa
                    .iter()
                    .any(|&p| !normalized(p))
            {
                return Err(MaterialFluxError::Invalid(
                    "complete representable physical compliance histories required",
                ));
            }
            Ok(history
                .compliance_pressures_pa
                .iter()
                .map(|p| p / density)
                .collect())
        };
    let recent_pressure = pressures(previous)?;
    let older_pressure = older_history
        .map(|(history, _)| pressures(history))
        .transpose()?;
    let native: Vec<_> = loads
        .iter()
        .map(|load| StokesWindkesselPort {
            reservoir_pressure_m2_s2: load.reservoir_pressure_pa / density,
            proximal_resistance_s_inv: load.proximal_resistance_pa_s_per_m2 / density,
            distal_resistance_s_inv: load.distal_resistance_pa_s_per_m2 / density,
            compliance_s2: load.compliance_m2_per_pa * density,
            inertance: load.inertance_pa_s2_per_m2 / density,
        })
        .collect();
    let old = older_history
        .zip(older_pressure.as_deref())
        .map(|((history, dt), pressure)| {
            (
                StokesWindkesselHistory {
                    velocity: history.velocity,
                    compliance_pressures_m2_s2: pressure,
                },
                dt,
            )
        });
    solve(
        &native,
        StokesWindkesselHistory {
            velocity: previous.velocity,
            compliance_pressures_m2_s2: &recent_pressure,
        },
        old,
    )
}

fn physical_windkessel<N, S>(
    solved: StokesWindkesselSolution<N>,
    density: f64,
    recover: impl FnOnce(N, f64) -> Result<S, MaterialFluxError>,
) -> Result<MaterialStokesWindkesselSolution<S>, MaterialFluxError> {
    let kinetic_energy = density * solved.kinetic_energy_m4_s2;
    let pressure_offsets_pa: Vec<_> = solved
        .pressure_offsets_m2_s2
        .iter()
        .map(|p| density * p)
        .collect();
    let maximum_pressure_residual_pa = density * solved.maximum_pressure_residual_m2_s2;
    let balances = solved
        .balances
        .into_iter()
        .map(|balance| physical_balance(balance, density))
        .collect::<Result<Vec<_>, _>>()?;
    if pressure_offsets_pa
        .iter()
        .chain([&kinetic_energy, &maximum_pressure_residual_pa])
        .any(|v| !v.is_finite())
    {
        return Err(MaterialFluxError::Invalid(
            "unrepresentable physical hydraulic field or pressure",
        ));
    }
    Ok(MaterialStokesWindkesselSolution {
        flow: recover(solved.flow, kinetic_energy)?,
        pressure_offsets_pa,
        outward_volume_rates_m2_s: solved.outward_volume_rates_m2_s,
        balances,
        maximum_pressure_residual_pa,
        maximum_flow_residual_m2_s: solved.maximum_flow_residual_m2_s,
    })
}

fn physical_balance(
    b: StokesWindkesselBalance,
    density: f64,
) -> Result<MaterialStokesWindkesselBalance, MaterialFluxError> {
    let result = MaterialStokesWindkesselBalance {
        compliance_pressure_pa: density * b.compliance_pressure_m2_s2,
        compliance_energy_j_per_m: density * b.compliance_energy_m4_s2,
        inertance_energy_j_per_m: density * b.inertance_energy_m4_s2,
        proximal_dissipation_w_per_m: density * b.proximal_dissipation_m4_s3,
        distal_dissipation_w_per_m: density * b.distal_dissipation_m4_s3,
        reservoir_power_w_per_m: density * b.reservoir_power_m4_s3,
        port_power_w_per_m: density * b.port_power_m4_s3,
        stored_energy_change_rate_w_per_m: density * b.stored_energy_change_rate_m4_s3,
        bdf_energy_defect_w_per_m: density * b.bdf_energy_defect_m4_s3,
        power_balance_residual_w_per_m: density * b.power_balance_residual_m4_s3,
        pressure_law_residual_pa: density * b.pressure_law_residual_m2_s2,
        compliance_flow_residual_m2_s: b.compliance_flow_residual_m2_s,
    };
    if ![
        result.compliance_pressure_pa,
        result.compliance_energy_j_per_m,
        result.inertance_energy_j_per_m,
        result.proximal_dissipation_w_per_m,
        result.distal_dissipation_w_per_m,
        result.reservoir_power_w_per_m,
        result.port_power_w_per_m,
        result.stored_energy_change_rate_w_per_m,
        result.bdf_energy_defect_w_per_m,
        result.power_balance_residual_w_per_m,
        result.pressure_law_residual_pa,
        result.compliance_flow_residual_m2_s,
    ]
    .iter()
    .all(|v| v.is_finite())
    {
        return Err(MaterialFluxError::Invalid(
            "unrepresentable physical hydraulic storage or power",
        ));
    }
    Ok(result)
}
