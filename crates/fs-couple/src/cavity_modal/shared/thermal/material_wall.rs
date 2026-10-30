//! Conjugate radial gas/solid heat conduction in a compact cylindrical cavity.
//!
//! Cell capacities and cylindrical Fourier resistances come from the original
//! gas and selected wall material. The outside wall temperature is prescribed
//! at the reference gas temperature; endcaps are insulated. Pressure is uniform.
//! For cell temperatures theta, C theta_dot - v p_dot = -L theta, where v contains
//! gas cell volumes and zero for the solid and local molecular reservoirs.
//! Eliminating p = -p0*deltaV/V + p0*v'theta/(T0*V) produces the SAME positive
//! compression realization as the isothermal-wall molecular owner.
//!
//! This is a convergent radial finite-volume discretization, not an effective
//! heat-transfer coefficient. Gas cells and wall cells are independent spatial
//! refinement axes. Every resulting memory coordinate joins the existing exact
//! mechanical propagator; no additional thermal energy ledger is introduced.
//! The pressure-driven caloric equation and linearized equation of state use
//! the standard thermoviscous equations (with ideal-gas alpha=1/T):
//! https://www.comsol.com/blogs/theory-of-thermoviscous-acoustics-thermal-and-viscous-losses/

use super::*;
use fs_material::thermal::transport::{PreparedThermalTransport, ThermalTransportSample};

/// Actual wall geometry and explicit cold-work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CylindricalMaterialWallOptions {
    /// Inner radius [m]; the original cavity volume determines cylinder length.
    pub radius_m: f64,
    /// Positive radial thickness of the solid wall [m].
    pub thickness_m: f64,
    /// Uniform radial gas cells; increase independently to resolve gas diffusion.
    pub gas_cells: usize,
    /// Uniform radial solid cells; increase independently to resolve the wall.
    pub wall_cells: usize,
    /// Declared harmonic-response band [rad/s], within the runtime Nyquist guard.
    /// This is a domain declaration, not a spatial accuracy certificate.
    pub band_rad_s: [f64; 2],
    /// Upper bound on gas temperatures, local reservoirs and wall temperatures
    /// together. The dense cold eigensolve is additionally capped at 128.
    pub maximum_coordinates: usize,
}

/// Material-driven gas/wall compression, retaining both original source owners.
#[derive(Clone, Debug, PartialEq)]
pub struct CylindricalMaterialWallCompression {
    options: CylindricalMaterialWallOptions,
    length_m: f64,
    wall: PreparedThermalTransport,
    wall_sample: ThermalTransportSample,
    molecular_relaxation: Option<RelaxationState>,
    equilibrium_stiffness_n_m: f64,
    instantaneous_stiffness_n_m: f64,
    arms: Vec<ThermalCompressionArm>,
    thermal_coordinates: usize,
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

impl CylindricalMaterialWallCompression {
    /// Freeze the original wall source at the gas reference temperature.
    ///
    /// The selected wall conductivity must describe radial heat flow and its
    /// selected heat capacity must apply to the wall's fixed-geometry caloric
    /// constraint. A scalar source does not establish isotropy. No phase change,
    /// mean heating, wall expansion, contact resistance or exterior convection
    /// is inferred. Rebuild when the reference temperature, gas or wall changes.
    ///
    /// Optional molecular reservoirs couple locally to EACH gas temperature
    /// before pressure elimination; they never become a second compression law.
    /// They retain their original local Debye rates and do not diffuse separately.
    pub fn new(
        source: &HelmholtzModalCavity,
        options: CylindricalMaterialWallOptions,
        wall: PreparedThermalTransport,
        molecular_relaxation: Option<RelaxationState>,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled);
        }
        if !positive(options.radius_m)
            || !positive(options.thickness_m)
            || options.gas_cells == 0
            || options.wall_cells == 0
            || options.maximum_coordinates == 0
            || options.maximum_coordinates > spectrum::MAX_COORDINATES
            || !positive(options.band_rad_s[0])
            || !positive(options.band_rad_s[1])
            || options.band_rad_s[0] >= options.band_rad_s[1]
        {
            return Err(network_invalid(
                "material cavity wall needs positive geometry, cells, band and bounded work",
            ));
        }
        let gas = source.gas();
        if molecular_relaxation
            .as_ref()
            .is_some_and(|law| law.gas() != gas)
        {
            return Err(network_invalid(
                "material wall and molecular owners need the same gas snapshot",
            ));
        }
        let wall_sample = wall.at_temperature(gas.temperature).map_err(|_| {
            network_invalid(
                "original wall transport source cannot represent the cavity temperature",
            )
        })?;
        let processes: Vec<_> = molecular_relaxation
            .as_ref()
            .into_iter()
            .flat_map(|law| law.processes())
            .filter(|p| p.heat_capacity_j_kg_k > 0.0)
            .collect();
        let width = processes
            .len()
            .checked_add(1)
            .ok_or_else(|| network_invalid("material cavity wall coordinate count overflow"))?;
        let gas_coordinates = options
            .gas_cells
            .checked_mul(width)
            .ok_or_else(|| network_invalid("material cavity wall coordinate count overflow"))?;
        let n = gas_coordinates
            .checked_add(options.wall_cells)
            .filter(|&n| n <= options.maximum_coordinates)
            .ok_or_else(|| {
                network_invalid("material cavity wall exceeds its thermal coordinate budget")
            })?;
        let radius = options.radius_m;
        let outer = radius + options.thickness_m;
        let volume = source.geometry().volume_m3;
        let length_m = volume / (core::f64::consts::PI * radius * radius);
        let cp_inf = molecular_relaxation
            .as_ref()
            .map_or(gas.specific_heat_cp, |law| {
                law.frozen_cv_j_kg_k() + gas.specific_gas_constant
            });
        if !positive(length_m)
            || !positive(outer)
            || outer <= radius
            || !positive(cp_inf)
            || !positive(gas.thermal_conductivity)
        {
            return Err(network_invalid(
                "material cavity thermal geometry or gas transport is unrepresentable",
            ));
        }
        let ring = core::f64::consts::PI * length_m;
        let mut capacity = vec![0.0; n];
        let mut forcing = vec![0.0; n];
        let mut conductance = vec![0.0; n * n];
        let mut centres = Vec::with_capacity(options.gas_cells + options.wall_cells);
        let mut radial_indices = Vec::with_capacity(centres.capacity());
        for cell in 0..options.gas_cells {
            let inner_r = radius * (cell as f64 / options.gas_cells as f64);
            let outer_r = radius * ((cell + 1) as f64 / options.gas_cells as f64);
            let cell_volume = ring * (outer_r - inner_r) * (outer_r + inner_r);
            let i = cell * width;
            capacity[i] = gas.density * cp_inf * cell_volume;
            forcing[i] = cell_volume;
            centres.push(det::sqrt(0.5 * (inner_r * inner_r + outer_r * outer_r)));
            radial_indices.push(i);
            for (j, process) in processes.iter().enumerate() {
                let k = i + j + 1;
                capacity[k] = gas.density * process.heat_capacity_j_kg_k * cell_volume;
                add_link(
                    &mut conductance,
                    n,
                    i,
                    k,
                    capacity[k] / process.relaxation_time_s,
                )?;
            }
        }
        for cell in 0..options.wall_cells {
            let inner_r = radius + options.thickness_m * (cell as f64 / options.wall_cells as f64);
            let outer_r =
                radius + options.thickness_m * ((cell + 1) as f64 / options.wall_cells as f64);
            let cell_volume = ring * (outer_r - inner_r) * (outer_r + inner_r);
            let i = gas_coordinates + cell;
            capacity[i] = wall_sample.volumetric_heat_capacity_j_m3_k() * cell_volume;
            centres.push(det::sqrt(0.5 * (inner_r * inner_r + outer_r * outer_r)));
            radial_indices.push(i);
        }
        if capacity.iter().any(|v| !positive(*v)) || centres.iter().any(|v| !positive(*v)) {
            return Err(network_invalid(
                "material cavity cell capacities or radial centres are unrepresentable",
            ));
        }
        let wall_k = wall_sample.conductivity_w_m_k();
        for cell in 1..centres.len() {
            if gate.is_requested() {
                return Err(ModalCouplingError::Cancelled);
            }
            let resistance = if cell == options.gas_cells {
                det::ln(radius / centres[cell - 1]) / gas.thermal_conductivity
                    + det::ln(centres[cell] / radius) / wall_k
            } else {
                let k = if cell < options.gas_cells {
                    gas.thermal_conductivity
                } else {
                    wall_k
                };
                det::ln(centres[cell] / centres[cell - 1]) / k
            };
            add_link(
                &mut conductance,
                n,
                radial_indices[cell - 1],
                radial_indices[cell],
                2.0 * ring / resistance,
            )?;
        }
        let last = n - 1;
        let exterior = 2.0 * ring * wall_k / det::ln(outer / centres[centres.len() - 1]);
        if !positive(exterior) {
            return Err(network_invalid(
                "material cavity exterior radial resistance is unrepresentable",
            ));
        }
        conductance[last * n + last] += exterior;
        let roots: Vec<_> = capacity.iter().map(|&c| det::sqrt(c)).collect();
        for i in 0..n {
            forcing[i] /= roots[i];
            for j in 0..n {
                conductance[i * n + j] /= roots[i] * roots[j];
            }
        }
        let spectrum = spectrum::realize(
            gas.pressure,
            gas.pressure / (gas.temperature * volume),
            conductance,
            &forcing,
            &mut || gate.is_requested(),
        )
        .map_err(|reason| {
            if gate.is_requested() {
                ModalCouplingError::Cancelled
            } else {
                network_invalid(reason)
            }
        })?;
        let area = source.neck_area_m2();
        let scale = area * area / volume;
        let equilibrium_stiffness_n_m = gas.pressure * scale;
        let instantaneous_stiffness_n_m = spectrum.frozen_bulk_modulus_pa * scale;
        if !positive(equilibrium_stiffness_n_m) || !positive(instantaneous_stiffness_n_m) {
            return Err(network_invalid(
                "material cavity compression scaling is unrepresentable",
            ));
        }
        let mut arms = Vec::with_capacity(spectrum.arms.len());
        for arm in spectrum.arms {
            let stiffness = arm.increment_pa * scale;
            if !positive(stiffness) {
                return Err(network_invalid(
                    "material cavity compression arm is unrepresentable",
                ));
            }
            arms.push(ThermalCompressionArm {
                stiffness_n_m: stiffness,
                relaxation_time_s: arm.relaxation_time_s,
                root_stiffness: det::sqrt(stiffness),
            });
        }
        Ok(Self {
            options,
            length_m,
            wall,
            wall_sample,
            molecular_relaxation,
            equilibrium_stiffness_n_m,
            instantaneous_stiffness_n_m,
            arms,
            thermal_coordinates: n,
        })
    }

    /// Original cylindrical geometry, response domain and refinement choices.
    pub const fn options(&self) -> CylindricalMaterialWallOptions {
        self.options
    }
    /// Cylinder length derived from original cavity volume [m].
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    /// Original source path, material identity and query receipts.
    pub fn wall(&self) -> &PreparedThermalTransport {
        &self.wall
    }
    /// Original wall source evaluated at the reference gas temperature.
    pub const fn wall_sample(&self) -> ThermalTransportSample {
        self.wall_sample
    }
    /// Original local molecular source, if enabled.
    pub fn molecular_relaxation(&self) -> Option<&RelaxationState> {
        self.molecular_relaxation.as_ref()
    }
    /// Number of physical thermal coordinates used in the bounded cold solve.
    pub const fn thermal_coordinates(&self) -> usize {
        self.thermal_coordinates
    }
    /// Isothermal static stiffness [N/m], including heat escape through the wall.
    pub const fn equilibrium_stiffness_n_m(&self) -> f64 {
        self.equilibrium_stiffness_n_m
    }
    /// Instantaneous molecularly frozen adiabatic stiffness [N/m].
    pub const fn instantaneous_stiffness_n_m(&self) -> f64 {
        self.instantaneous_stiffness_n_m
    }
    /// Positive gas, wall and molecular compression memories, in runtime order.
    pub fn arms(&self) -> &[ThermalCompressionArm] {
        &self.arms
    }
    /// Complex stiffness [N/m], with the +i omega convention in the declared band.
    pub fn stiffness_at(&self, omega: f64) -> Result<C64, ModalCouplingError> {
        compression_stiffness(
            self.equilibrium_stiffness_n_m,
            &self.arms,
            self.options.band_rad_s,
            omega,
        )
    }
    pub(super) fn observation(
        &self,
        extension: f64,
        memory: &[f64],
    ) -> Result<(f64, f64), ModalCouplingError> {
        compression_observation(
            self.equilibrium_stiffness_n_m,
            &self.arms,
            extension,
            memory,
        )
    }
}

fn add_link(
    matrix: &mut [f64],
    n: usize,
    i: usize,
    j: usize,
    g: f64,
) -> Result<(), ModalCouplingError> {
    if !positive(g) {
        return Err(network_invalid(
            "material cavity radial or molecular conductance is unrepresentable",
        ));
    }
    matrix[i * n + i] += g;
    matrix[j * n + j] += g;
    matrix[i * n + j] -= g;
    matrix[j * n + i] -= g;
    Ok(())
}

impl SharedCavity {
    /// Couple material wall conduction, local gas reservoirs and moving surfaces
    /// in ONE compression solve. Wall/gas swaps change source data, not this code.
    ///
    /// The exterior wall is held at the reference gas temperature and endcaps
    /// are insulated. This remains a small-signal, spatially uniform-pressure
    /// cavity with fixed geometry. The wall source selects radial conductivity
    /// and fixed-geometry heat capacity; coefficients freeze at gas temperature.
    /// `None` initializes the joint memory fully relaxed; zeros describe a sudden
    /// frozen-adiabatic displacement. All arms use the original accepted-state,
    /// exact propagation, energy, dissipation and rollback machinery.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_cylindrical_material_wall(
        source: HelmholtzModalCavity,
        surfaces: Vec<CavitySurface>,
        sample_rate_hz: u32,
        neck_initial: ModalAcousticState,
        neck_budget: ModalAcousticTimeBudget,
        config: ModalCouplingConfig,
        neck_options: Option<ViscousNeckOptions>,
        options: CylindricalMaterialWallOptions,
        wall: PreparedThermalTransport,
        molecular_relaxation: Option<RelaxationState>,
        initial_thermal_memory: Option<Vec<f64>>,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled);
        }
        if options.band_rad_s[1]
            > core::f64::consts::PI * sample_rate_hz as f64 * config.nyquist_guard_fraction
        {
            return Err(network_invalid(
                "material cavity thermal band exceeds the original Nyquist guard",
            ));
        }
        let thermal = CylindricalMaterialWallCompression::new(
            &source,
            options,
            wall,
            molecular_relaxation,
            gate,
        )?;
        let viscous = neck_options
            .map(|options| {
                Self::prepare_viscous_neck(&source, sample_rate_hz, config, options, gate)
            })
            .transpose()?;
        Self::new_inner(
            source,
            surfaces,
            sample_rate_hz,
            neck_initial,
            neck_budget,
            config,
            viscous,
            Some((Compression::Material(thermal), initial_thermal_memory)),
            None,
            gate,
        )
    }
}

