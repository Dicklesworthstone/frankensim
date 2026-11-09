//! Conjugate thermal boundaries of the original annular acoustic section.
//!
//! The original radial P1 mass/stiffness and native tridiagonal solve now retain
//! the thermal boundary nodes. Each wall is isothermal, adiabatic, or has an
//! explicitly supplied areal heat capacity C_A and conductance G to the
//! unchanged background bath. Perfect local gas/wall thermal contact gives
//! `q_out = (G - i omega C_A) T_wall` under exp(-i omega t).
//!
//! The gas obeys `rho cp(s) s T - conductivity Laplacian(T) = s p` with its
//! original, possibly molecular, Cp. Wall storage enters that SAME equation;
//! it is not an additional acoustic loss or an equilibrium thermal shunt.
//! The finite wall is a local lumped-capacity boundary. Its through-thickness
//! uniformity and external conductance require caller justification; neither
//! a material value, contact resistance nor an exterior environment is guessed.
//!
//! Fourier, reservoir and bath-exchange terms below are acoustic availability
//! losses. They do not by themselves assign mean caloric heating to gas/wall.
//! Oscillatory heat flux and wall thermal storage are reported separately.
//! Mechanical walls remain rigid/no-slip; axial wall conduction, radiation,
//! thermoelastic motion and background thermal evolution are not modeled.

/// Source-bound finite-thickness solid conduction in the same gas boundary solve.
pub mod conduction;

use super::{
    AnnularAdmission, AnnularSection, FluidDuctError, RadialSolution, bad, filled, finite, poll,
    positive, radial_mass, relative,
};
use crate::{DuctError, SegmentWave, compose::TwoPort};
use fs_la::eigen_complex::solve_tridiagonal_complex;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::gas::GasState;
use fs_material::relaxation::RelaxationState;
use fs_math::{c64::C64, dd::Dd, det};

/// Independently authored inner/outer thermal boundary, at the original gas's
/// background temperature. Coefficients carry no inferred material authority.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AnnularThermalWall {
    /// Zero temperature perturbation: the original isothermal limit.
    Isothermal,
    /// Zero outward heat flux, retaining no-slip mechanical shear.
    Adiabatic,
    /// One local wall temperature, in perfect thermal contact with the gas.
    /// G transfers heat to the background bath; C_A stores it reversibly.
    Finite {
        /// Conductance to the fixed background bath [W/(m² K)], nonnegative.
        conductance_w_m2_k: f64,
        /// Wall heat capacity per interface area [J/(m² K)], nonnegative.
        areal_heat_capacity_j_m2_k: f64,
    },
}

impl AnnularThermalWall {
    pub(crate) fn coefficients(self, omega: f64) -> Result<C64, DuctError> {
        let value = match self {
            Self::Isothermal | Self::Adiabatic => C64::ZERO,
            Self::Finite {
                conductance_w_m2_k,
                areal_heat_capacity_j_m2_k,
            } => {
                if !conductance_w_m2_k.is_finite()
                    || conductance_w_m2_k < 0.0
                    || !areal_heat_capacity_j_m2_k.is_finite()
                    || areal_heat_capacity_j_m2_k < 0.0
                {
                    return Err(bad(
                        "annular wall conductance and heat capacity must be finite and nonnegative",
                    ));
                }
                C64::new(conductance_w_m2_k, omega * areal_heat_capacity_j_m2_k)
            }
        };
        if !finite(value) {
            return Err(bad("annular wall thermal admittance is unrepresentable"));
        }
        Ok(value)
    }

    pub(crate) fn capacity(self) -> f64 {
        match self {
            Self::Finite {
                areal_heat_capacity_j_m2_k,
                ..
            } => areal_heat_capacity_j_m2_k,
            _ => 0.0,
        }
    }
}

/// Actual P1 gas fields in one original annular cell, with native peak phasors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnnularThermalWallField {
    pub radius_m: f64,
    pub pressure_pa: C64,
    pub axial_velocity_m_s: C64,
    pub radial_velocity_gradient_per_s: C64,
    pub gas_temperature_k: C64,
    pub radial_temperature_gradient_k_m: C64,
    /// Direct irreversible shear work [W/m³].
    pub viscous_dissipation_w_m3: f64,
    /// Fourier acoustic availability loss [W/m³].
    pub conduction_dissipation_w_m3: f64,
    /// Original molecular-reservoir acoustic availability loss [W/m³].
    pub molecular_dissipation_w_m3: f64,
}

/// Wall quantities per physical interface area or axial metre, never conflated
/// with the gas's local caloric source. Side 0 is inner; side 1 is outer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnnularThermalWallObservation {
    pub radius_m: f64,
    pub temperature_k: C64,
    /// Outward from fluid into this wall [W/m²]. Isothermal walls use the
    /// original P1 weak boundary reaction, not a guessed gas-film thickness.
    pub fluid_to_wall_heat_flux_w_m2: C64,
    /// Oscillatory heat delivered to the fixed bath [W/m²].
    pub wall_to_bath_heat_flux_w_m2: C64,
    /// Positive cycle-mean acoustic availability loss to this bath [W/m].
    pub bath_dissipation_w_m: f64,
    /// Cycle-mean quadratic wall thermal storage C_A |T'|²/(4 T0) [J/m].
    pub thermal_storage_j_m: f64,
}

/// Independent positive acoustic losses per axial metre, at one coherent p/Q.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AnnularThermalWallDissipation {
    pub viscous_w_m: f64,
    pub conduction_w_m: f64,
    pub molecular_w_m: f64,
    /// Bath-exchange loss on [inner, outer] walls [W/m].
    pub walls_w_m: [f64; 2],
    pub total_w_m: f64,
}

/// Original geometry/material bound to a finite-temperature thermal solve.
/// Preparation is O(N), subsequent cell/wall/whole-section queries O(1), with
/// no allocation or material queries. The original radial mesh is the spatial
/// accuracy axis; algebraic/power residuals are not discretization certificates.
#[derive(Debug)]
pub struct PreparedAnnularThermalWalls<'geometry, 'material, M: ?Sized> {
    section: &'geometry AnnularSection,
    source: &'material M,
    state: ThermoviscousFluidState,
    omega: f64,
    walls: [AnnularThermalWall; 2],
    wave: SegmentWave,
    viscous: RadialSolution,
    temperature_response: Vec<C64>,
    mean_temperature_response: C64,
    velocity_scale: C64,
    temperature_scale: C64,
    local_loss_roots: [f64; 3],
    total_loss_roots: [f64; 5],
    wall_storage_roots: [f64; 2],
    wall_flux_per_pressure: [C64; 2],
    wall_bath_flux_per_pressure: [C64; 2],
    relative_backward_residual: f64,
    relative_thermal_work_defect: f64,
}

impl AnnularSection {
    /// Original gas with independent inner/outer thermal boundaries.
    pub fn prepare_thermal_wall_fields<'g, 'm>(
        &'g self,
        gas: &'m GasState,
        omega: f64,
        walls: [AnnularThermalWall; 2],
        admission: AnnularAdmission,
    ) -> Result<PreparedAnnularThermalWalls<'g, 'm, GasState>, FluidDuctError> {
        self.prepare_thermal_wall_fields_controlled(gas, omega, walls, admission, || false)
    }

    /// Same native O(N) solve, polling every 64 cells and before/after the
    /// bounded tridiagonal tile (at most 4097 nodes). The explicit payload cap
    /// must cover 320*(elements+1)+4096 bytes, including original geometry,
    /// retained fields and native solve scratch, excluding allocator metadata.
    pub fn prepare_thermal_wall_fields_controlled<'g, 'm>(
        &'g self,
        gas: &'m GasState,
        omega: f64,
        walls: [AnnularThermalWall; 2],
        admission: AnnularAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAnnularThermalWalls<'g, 'm, GasState>, FluidDuctError> {
        self.admit_thermal_walls(omega, walls, admission, &mut cancelled)?;
        let state = gas.thermoviscous_state()?;
        prepare(
            self,
            gas,
            gas,
            state,
            C64::from_re(state.transport().specific_heat_j_kg_k()),
            omega,
            walls,
            admission,
            &mut cancelled,
        )
    }

    /// Original molecular Cp and wall temperature solved together. No separate
    /// equilibrium thermal, Stokes--Kirchhoff or bulk-viscosity term is added.
    pub fn prepare_molecular_thermal_wall_fields<'g, 'm>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        walls: [AnnularThermalWall; 2],
        admission: AnnularAdmission,
    ) -> Result<PreparedAnnularThermalWalls<'g, 'm, RelaxationState>, FluidDuctError> {
        self.prepare_molecular_thermal_wall_fields_controlled(
            material,
            omega,
            walls,
            admission,
            || false,
        )
    }

    /// Molecular preparation under the same cancellation and payload bounds.
    pub fn prepare_molecular_thermal_wall_fields_controlled<'g, 'm>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        walls: [AnnularThermalWall; 2],
        admission: AnnularAdmission,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAnnularThermalWalls<'g, 'm, RelaxationState>, FluidDuctError> {
        self.admit_thermal_walls(omega, walls, admission, &mut cancelled)?;
        let gas = material.gas();
        let state = gas.thermoviscous_state()?;
        let cp = if material.processes().is_empty() {
            C64::from_re(state.transport().specific_heat_j_kg_k())
        } else {
            material.caloric_response(omega)?.cp_j_kg_k
        };
        poll(&mut cancelled, 0)?;
        prepare(
            self,
            material,
            gas,
            state,
            cp,
            omega,
            walls,
            admission,
            &mut cancelled,
        )
    }

    fn admit_thermal_walls(
        &self,
        omega: f64,
        walls: [AnnularThermalWall; 2],
        admission: AnnularAdmission,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DuctError> {
        poll(cancelled, 0)?;
        self.admit(omega, admission)?;
        if 320 * (self.elements() + 1) + 4096 > admission.maximum_payload_bytes {
            return Err(bad(
                "annular finite-wall fields exceed the declared payload budget",
            ));
        }
        for wall in walls {
            wall.coefficients(omega)?;
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn prepare<'g, 'm, M: ?Sized>(
    section: &'g AnnularSection,
    source: &'m M,
    gas: &GasState,
    state: ThermoviscousFluidState,
    cp: C64,
    omega: f64,
    walls: [AnnularThermalWall; 2],
    admission: AnnularAdmission,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedAnnularThermalWalls<'g, 'm, M>, FluidDuctError> {
    let transport = state.transport();
    let k0 = omega / state.acoustic().sound_speed_m_s();
    if !positive(k0) || k0 * section.outer_m > admission.maximum_transverse_phase {
        return Err(bad("annular finite-wall gas exceeds its plane-pressure envelope").into());
    }
    let viscous = section.diffusion_solution(
        omega,
        transport.kinematic_viscosity_m2_s(),
        admission,
        cancelled,
    )?;
    let conductivity = transport.thermal_conductivity_w_m_k();
    let argument = C64::new(cp.im, cp.re)
        .scale((omega * gas.density / conductivity) * section.gap_m * section.gap_m);
    let phase = det::sqrt(argument.abs()) * section.maximum_cell_fraction;
    if !finite(argument)
        || !positive(argument.im)
        || argument.re < 0.0
        || !phase.is_finite()
        || phase > admission.maximum_cell_diffusion_phase
    {
        return Err(bad("annular mesh does not resolve finite-wall thermal diffusion").into());
    }
    let wall_admittance = [walls[0].coefficients(omega)?, walls[1].coefficients(omega)?];
    let mean_radius = section.inner_m + 0.5 * section.gap_m;
    let radii = [section.inner_m, section.outer_m];
    let boundary = core::array::from_fn::<_, 2, _>(|side| {
        wall_admittance[side].scale(section.gap_m / conductivity * (radii[side] / mean_radius))
    });
    if !boundary.into_iter().all(finite) {
        return Err(bad("annular scaled thermal boundary is unrepresentable").into());
    }
    let (nodes, backward) = solve_thermal(section, argument, boundary, walls, cancelled)?;
    let n = section.elements();
    let (mut mean_re, mut mean_im) = (Dd::ZERO, Dd::ZERO);
    for (i, value) in nodes.iter().enumerate() {
        if i % 64 == 0 {
            poll(cancelled, i)?;
        }
        let weight = nodal_load(section, i);
        mean_re = mean_re + Dd::from_f64(weight) * Dd::from_f64(value.re);
        mean_im = mean_im + Dd::from_f64(weight) * Dd::from_f64(value.im);
    }
    let mean = C64::new(mean_re.to_f64(), mean_im.to_f64());
    let temperature_scale = cp.recip() / C64::from_re(gas.density);
    let velocity_scale = viscous.response.mean.conj().recip() / C64::from_re(section.area_m2);
    let local_loss_roots = [
        det::sqrt(transport.dynamic_viscosity_pa_s()) * core::f64::consts::FRAC_1_SQRT_2,
        det::sqrt(conductivity) / det::sqrt(gas.temperature) * core::f64::consts::FRAC_1_SQRT_2,
        det::sqrt(gas.density) * det::sqrt(omega) * det::sqrt(cp.im) / det::sqrt(gas.temperature)
            * core::f64::consts::FRAC_1_SQRT_2,
    ];
    let (mut viscous_form, mut conduction_form, mut molecular_form) =
        (Dd::ZERO, Dd::ZERO, Dd::ZERO);
    for i in 0..n {
        if i % 64 == 0 {
            poll(cancelled, i)?;
        }
        viscous_form = viscous_form
            + Dd::from_f64(section.stiffness_edges[i])
                * Dd::from_f64(viscous.difference(i).norm_sq());
        conduction_form = conduction_form
            + Dd::from_f64(section.stiffness_edges[i])
                * Dd::from_f64((nodes[i + 1] - nodes[i]).norm_sq());
        let x = section.fractions[i];
        let width = section.fractions[i + 1] - x;
        let left = section.inner_m / mean_radius + section.gap_m / mean_radius * x;
        let right = section.inner_m / mean_radius + section.gap_m / mean_radius * (x + width);
        let (diagonal, off) = radial_mass(width, left, right);
        let ratio = off / diagonal[0];
        let form = (nodes[i] + nodes[i + 1].scale(ratio)).norm_sq() * diagonal[0]
            + nodes[i + 1].norm_sq() * (diagonal[1] - off * ratio);
        molecular_form = molecular_form + Dd::from_f64(form);
    }
    let area_root = det::sqrt(section.area_m2);
    let mut total_loss_roots = [
        velocity_scale.scale(local_loss_roots[0]).abs() * area_root / section.gap_m
            * det::sqrt(viscous_form.to_f64()),
        temperature_scale.scale(local_loss_roots[1]).abs() * area_root / section.gap_m
            * det::sqrt(conduction_form.to_f64()),
        temperature_scale.scale(local_loss_roots[2]).abs()
            * area_root
            * det::sqrt(molecular_form.to_f64()),
        0.0,
        0.0,
    ];
    let mut wall_storage_roots = [0.0; 2];
    let mut wall_flux = [C64::ZERO; 2];
    let mut bath_flux = [C64::ZERO; 2];
    for side in 0..2 {
        let node = side * n;
        let t = nodes[node].conj() * temperature_scale;
        total_loss_roots[3 + side] = t.abs()
            * det::sqrt(
                core::f64::consts::PI * radii[side] * wall_admittance[side].re / gas.temperature,
            );
        wall_storage_roots[side] = t.abs()
            * det::sqrt(
                core::f64::consts::PI * radii[side] * walls[side].capacity()
                    / (2.0 * gas.temperature),
            );
        wall_flux[side] = if walls[side] == AnnularThermalWall::Isothermal {
            let edge = if side == 0 { 0 } else { n - 1 };
            let other = if side == 0 { 1 } else { n - 1 };
            let off = C64::from_re(-section.stiffness_edges[edge])
                + argument.scale(section.mass_edges[edge]);
            let reaction = argument.scale(nodal_load(section, node)) - off * nodes[other];
            reaction.conj()
                * temperature_scale
                * C64::from_re(conductivity / section.gap_m * mean_radius / radii[side])
        } else {
            t * wall_admittance[side].conj()
        };
        bath_flux[side] = if walls[side] == AnnularThermalWall::Isothermal {
            wall_flux[side]
        } else {
            t.scale(wall_admittance[side].re)
        };
    }
    let compliance = (C64::ONE - mean * cp.conj().recip().scale(gas.specific_gas_constant))
        .scale(1.0 / (gas.density * gas.specific_gas_constant * gas.temperature));
    let s = C64::new(0.0, omega);
    let series = s.scale(gas.density / section.area_m2) / viscous.response.mean;
    let shunt = s.scale(section.area_m2) * compliance;
    let (wavenumber, characteristic_impedance) =
        crate::helmholtz_wave_from_telegraph(series, shunt);
    let heat = total_loss_roots[1..].iter().map(|r| r * r).sum::<f64>();
    let heat_from_port = 0.5 * shunt.re;
    let work_defect = relative((heat - heat_from_port).abs(), heat + heat_from_port.abs());
    let wave = SegmentWave {
        wavenumber,
        characteristic_impedance,
        specific_impedance: characteristic_impedance.scale(section.area_m2),
        // Circular axes use the physical radius; annuli retain half their gap.
        shear_number: (if section.inner_m == 0.0 {
            section.outer_m
        } else {
            0.5 * section.gap_m
        }) * det::sqrt(omega / transport.kinematic_viscosity_m2_s()),
    };
    if ![
        mean,
        temperature_scale,
        velocity_scale,
        compliance,
        wavenumber,
        characteristic_impedance,
        wave.specific_impedance,
    ]
    .into_iter()
    .all(finite)
        || compliance.re <= 0.0
        || wavenumber.re <= 0.0
        || wavenumber.im < 0.0
        || characteristic_impedance.re <= 0.0
        || k0.max(wavenumber.abs()) * section.outer_m > admission.maximum_transverse_phase
        || !total_loss_roots
            .into_iter()
            .chain(wall_storage_roots)
            .all(|r| r.is_finite() && r >= 0.0)
        || !wall_flux.into_iter().chain(bath_flux).all(finite)
        || backward > admission.maximum_relative_residual
        || work_defect > admission.maximum_relative_residual
    {
        return Err(bad(
            "annular finite-wall field, wave or energy residual exceeds its admission",
        )
        .into());
    }
    poll(cancelled, n)?;
    Ok(PreparedAnnularThermalWalls {
        section,
        source,
        state,
        omega,
        walls,
        wave,
        viscous,
        temperature_response: nodes,
        mean_temperature_response: mean,
        velocity_scale,
        temperature_scale,
        local_loss_roots,
        total_loss_roots,
        wall_storage_roots,
        wall_flux_per_pressure: wall_flux,
        wall_bath_flux_per_pressure: bath_flux,
        relative_backward_residual: backward,
        relative_thermal_work_defect: work_defect,
    })
}

fn nodal_load(section: &AnnularSection, node: usize) -> f64 {
    if node == 0 {
        section.boundary_weights[0]
    } else if node == section.elements() {
        section.boundary_weights[1]
    } else {
        section.load[node - 1]
    }
}

// Same weighted P1 operator as the original Dirichlet owner, now retaining
// boundary nodes. Natural boundary terms add G + i omega C_A, once.
fn solve_thermal(
    section: &AnnularSection,
    argument: C64,
    boundary: [C64; 2],
    walls: [AnnularThermalWall; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<C64>, f64), DuctError> {
    let n = section.elements();
    // The exact constant field satisfies both insulated boundaries for ANY
    // original Cp. Preserve its zero Fourier/bath work instead of solving an
    // increasingly ill-conditioned Neumann Laplacian as frequency vanishes.
    if walls.iter().all(|w| *w != AnnularThermalWall::Isothermal) && boundary == [C64::ZERO; 2] {
        poll(cancelled, n)?;
        return Ok((filled(n + 1, C64::ONE)?, 0.0));
    }
    let scale = 1.0
        / argument
            .abs()
            .max(boundary[0].abs())
            .max(boundary[1].abs())
            .max(1.0);
    let z = argument.scale(scale);
    let mut diagonal = filled(n + 1, C64::ZERO)?;
    let mut off = filled(n, C64::ZERO)?;
    let mut values = filled(n + 1, C64::ZERO)?;
    for i in 0..=n {
        if i % 64 == 0 {
            poll(cancelled, i)?;
        }
        let (stiffness, mass) = if i == 0 {
            (
                section.stiffness_edges[0],
                section.boundary_weights[0] - section.mass_edges[0],
            )
        } else if i == n {
            (
                section.stiffness_edges[n - 1],
                section.boundary_weights[1] - section.mass_edges[n - 1],
            )
        } else {
            (
                section.stiffness_edges[i - 1] + section.stiffness_edges[i],
                section.mass_diagonal[i - 1],
            )
        };
        diagonal[i] = C64::from_re(scale * stiffness) + z.scale(mass);
        values[i] = z.scale(nodal_load(section, i));
        if i < n {
            off[i] =
                C64::from_re(-scale * section.stiffness_edges[i]) + z.scale(section.mass_edges[i]);
        }
    }
    for side in 0..2 {
        let node = side * n;
        if walls[side] == AnnularThermalWall::Isothermal {
            diagonal[node] = C64::ONE;
            off[if side == 0 { 0 } else { n - 1 }] = C64::ZERO;
            values[node] = C64::ZERO;
        } else {
            diagonal[node] = diagonal[node] + boundary[side].scale(scale);
        }
    }
    if diagonal
        .iter()
        .chain(&off)
        .chain(&values)
        .any(|&v| !finite(v))
    {
        return Err(bad("annular finite-wall thermal matrix is unrepresentable"));
    }
    poll(cancelled, n)?;
    solve_tridiagonal_complex(&off, &diagonal, &off, &mut values, 1)
        .map_err(|_| DuctError::Singular)?;
    poll(cancelled, n)?;
    let mut residual = 0.0_f64;
    for i in 0..=n {
        if i % 64 == 0 {
            poll(cancelled, i)?;
        }
        let is_fixed = (i == 0 && walls[0] == AnnularThermalWall::Isothermal)
            || (i == n && walls[1] == AnnularThermalWall::Isothermal);
        let drive = if is_fixed {
            C64::ZERO
        } else {
            z.scale(nodal_load(section, i))
        };
        let mut error = diagonal[i] * values[i] - drive;
        let mut magnitude = diagonal[i].abs() * values[i].abs() + drive.abs();
        if i > 0 {
            error = error + off[i - 1] * values[i - 1];
            magnitude += off[i - 1].abs() * values[i - 1].abs();
        }
        if i < n {
            error = error + off[i] * values[i + 1];
            magnitude += off[i].abs() * values[i + 1].abs();
        }
        residual = residual.max(relative(error.abs(), magnitude));
    }
    Ok((values, residual))
}

impl<'g, 'm, M: ?Sized> PreparedAnnularThermalWalls<'g, 'm, M> {
    pub const fn source(&self) -> &'m M {
        self.source
    }
    pub const fn section(&self) -> &'g AnnularSection {
        self.section
    }
    pub const fn state(&self) -> &ThermoviscousFluidState {
        &self.state
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn walls(&self) -> &[AnnularThermalWall; 2] {
        &self.walls
    }
    pub const fn wave(&self) -> &SegmentWave {
        &self.wave
    }
    /// Mean normalized thermal response, in the radial solve's +i omega convention.
    pub const fn mean_temperature_response(&self) -> C64 {
        self.mean_temperature_response
    }
    pub const fn relative_backward_residual(&self) -> f64 {
        self.relative_backward_residual
    }
    pub const fn relative_thermal_work_defect(&self) -> f64 {
        self.relative_thermal_work_defect
    }
    /// Same native finite-length pressure/flow scattering as the original owner.
    pub fn two_port(&self, length_m: f64, reference: f64) -> Result<TwoPort, DuctError> {
        TwoPort::from_uniform_wave(&self.wave, self.section.area_m2, length_m, reference)
    }

    /// Local fields in one original P1 cell; fraction is in [0,1]. A node's
    /// adjacent-cell derivatives differ at finite mesh resolution.
    pub fn field_in_cell(
        &self,
        cell: usize,
        fraction: f64,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<AnnularThermalWallField, DuctError> {
        amplitudes(pressure_pa, volume_flow_m3_s)?;
        if cell >= self.section.elements()
            || !fraction.is_finite()
            || !(0.0..=1.0).contains(&fraction)
        {
            return Err(bad(
                "annular finite-wall field cell or interpolation fraction is invalid",
            ));
        }
        let x = self.section.fractions[cell];
        let width = self.section.fractions[cell + 1] - x;
        let length = self.section.gap_m * width;
        let interpolate =
            |left: C64, right: C64| left.scale(1.0 - fraction) + right.scale(fraction);
        let u = interpolate(self.viscous.node(cell), self.viscous.node(cell + 1)).conj()
            * self.velocity_scale
            * volume_flow_m3_s;
        let du = self.viscous.difference(cell).conj() * self.velocity_scale / C64::from_re(length)
            * volume_flow_m3_s;
        let t = interpolate(
            self.temperature_response[cell],
            self.temperature_response[cell + 1],
        )
        .conj()
            * self.temperature_scale
            * pressure_pa;
        let dt = (self.temperature_response[cell + 1] - self.temperature_response[cell]).conj()
            * self.temperature_scale
            / C64::from_re(length)
            * pressure_pa;
        let losses = [
            du.scale(self.local_loss_roots[0]).norm_sq(),
            dt.scale(self.local_loss_roots[1]).norm_sq(),
            t.scale(self.local_loss_roots[2]).norm_sq(),
        ];
        if ![u, du, t, dt].into_iter().all(finite) || !losses.into_iter().all(f64::is_finite) {
            return Err(bad("annular finite-wall local field is unrepresentable"));
        }
        Ok(AnnularThermalWallField {
            radius_m: self.section.inner_m + self.section.gap_m * (x + fraction * width),
            pressure_pa,
            axial_velocity_m_s: u,
            radial_velocity_gradient_per_s: du,
            gas_temperature_k: t,
            radial_temperature_gradient_k_m: dt,
            viscous_dissipation_w_m3: losses[0],
            conduction_dissipation_w_m3: losses[1],
            molecular_dissipation_w_m3: losses[2],
        })
    }

    /// Side 0 is inner; side 1 is outer. Flux is outward FROM the fluid at
    /// either side, so an insulated wall has exactly zero flux in this view.
    pub fn wall_observation(
        &self,
        side: usize,
        pressure_pa: C64,
    ) -> Result<AnnularThermalWallObservation, DuctError> {
        amplitudes(pressure_pa, C64::ZERO)?;
        if side > 1 {
            return Err(bad("annular thermal wall side must be zero or one"));
        }
        let temperature = self.temperature_response[side * self.section.elements()].conj()
            * self.temperature_scale
            * pressure_pa;
        let flux = self.wall_flux_per_pressure[side] * pressure_pa;
        let bath_flux = self.wall_bath_flux_per_pressure[side] * pressure_pa;
        let loss = pressure_pa.scale(self.total_loss_roots[3 + side]).norm_sq();
        let storage = pressure_pa.scale(self.wall_storage_roots[side]).norm_sq();
        if ![temperature, flux, bath_flux].into_iter().all(finite)
            || !loss.is_finite()
            || !storage.is_finite()
        {
            return Err(bad("annular wall thermal observation is unrepresentable"));
        }
        Ok(AnnularThermalWallObservation {
            radius_m: if side == 0 {
                self.section.inner_m
            } else {
                self.section.outer_m
            },
            temperature_k: temperature,
            fluid_to_wall_heat_flux_w_m2: flux,
            wall_to_bath_heat_flux_w_m2: bath_flux,
            bath_dissipation_w_m: loss,
            thermal_storage_j_m: storage,
        })
    }

    /// Positive analytic P1 integrals, independent of the acoustic power
    /// difference used to check their consistency during preparation.
    pub fn dissipation_per_axial_length(
        &self,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<AnnularThermalWallDissipation, DuctError> {
        amplitudes(pressure_pa, volume_flow_m3_s)?;
        let values = core::array::from_fn::<_, 5, _>(|i| {
            (if i == 0 {
                volume_flow_m3_s
            } else {
                pressure_pa
            })
            .scale(self.total_loss_roots[i])
            .norm_sq()
        });
        let total = values.into_iter().sum::<f64>();
        if !values.into_iter().all(f64::is_finite) || !total.is_finite() {
            return Err(bad(
                "annular finite-wall integrated acoustic loss is unrepresentable",
            ));
        }
        Ok(AnnularThermalWallDissipation {
            viscous_w_m: values[0],
            conduction_w_m: values[1],
            molecular_w_m: values[2],
            walls_w_m: [values[3], values[4]],
            total_w_m: total,
        })
    }
}

fn amplitudes(p: C64, q: C64) -> Result<(), DuctError> {
    if !finite(p) || !finite(q) {
        Err(bad("annular pressure and volume flow must be finite"))
    } else {
        Ok(())
    }
}
