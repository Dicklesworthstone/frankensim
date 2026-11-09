//! Independent thermal storage/exchange on the four original rectangular walls.
//!
//! The gas solves rho cp(s) s T - k Laplacian(T) = s p on the actual 2-D
//! rectangle. Each physical face is isothermal, insulated, or has authored
//! conductance G and areal capacity C, q_out=(G-i omega C)T. Original molecular
//! Cp enters once. The original Fourier owner still supplies no-slip shear.
//!
//! A declared bilinear mesh resolves thermal corners and unequal walls; native
//! bounded sparse LDLT solves its original Q1 mass/stiffness/boundary forms.
//! The wall is locally uniform through its thickness, in perfect gas contact;
//! G/C require caller authority. No planar solid, empirical data, gas-film
//! coefficient, wall accommodation, or thermoelastic motion is fabricated.

/// Source-bound solid jackets with resolved transverse and corner conduction.
pub mod conduction;

use super::{
    FluidDuctError, RectangularAdmission, RectangularSection, admit_tail, bad,
    field::{DiffusionField, amplitudes},
    finite, poll, positive,
};
use crate::{DuctError, SegmentWave, compose::TwoPort};
use fs_material::{
    fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState},
    gas::GasState,
    relaxation::RelaxationState,
};
use fs_math::{c64::C64, det};
use fs_sparse::{LdltError, direct::DirectLimits};

/// The same explicitly authored local G/C law used by annular/circular walls.
/// Order is [short-coordinate zero, short-coordinate end, long zero, long end].
pub use super::super::annular::thermal_wall::AnnularThermalWall as RectangularThermalWall;

mod physics;
use physics::ThermalField;

/// Explicit Q1 thermal mesh and native sparse work/storage limits.
/// This is a spatial refinement choice, not a continuum-error certificate.
#[derive(Clone, Copy, Debug)]
pub struct RectangularThermalLimits {
    /// Gas elements along [short, long] sides; each >=2. The complete gas or
    /// conjugate gas/solid mesh admits at most 16384 total physical nodes.
    pub cells: [usize; 2],
    /// max(cell widths)*sqrt(|i omega rho Cp/k|), in (0,1]. A conjugate solid
    /// uses its original volumetric capacity and smallest conductivity eigenvalue.
    pub maximum_cell_diffusion_phase: f64,
    /// Original algebraic and physical thermal-work residual, in (0,1).
    pub maximum_relative_residual: f64,
    /// Conservative numerical payload: thermal staging, sparse analysis/factor,
    /// field arrays and original Fourier shear field. Excludes caller sources
    /// and allocator/runtime overhead. Native fill may refuse before reaching
    /// the node cap; no allocation or work limit is raised automatically.
    pub maximum_payload_bytes: usize,
    pub sparse: DirectLimits,
}

#[derive(Debug)]
pub enum RectangularThermalWallError {
    Fluid(FluidDuctError),
    Duct(DuctError),
    Sparse(LdltError),
    Transport(fs_material::thermal::transport::ThermalTransportError),
}
impl From<fs_material::thermal::transport::ThermalTransportError> for RectangularThermalWallError {
    fn from(v: fs_material::thermal::transport::ThermalTransportError) -> Self {
        Self::Transport(v)
    }
}
impl From<FluidDuctError> for RectangularThermalWallError {
    fn from(v: FluidDuctError) -> Self {
        Self::Fluid(v)
    }
}
impl From<DuctError> for RectangularThermalWallError {
    fn from(v: DuctError) -> Self {
        Self::Duct(v)
    }
}
impl From<LdltError> for RectangularThermalWallError {
    fn from(v: LdltError) -> Self {
        Self::Sparse(v)
    }
}
impl From<fs_material::MaterialError> for RectangularThermalWallError {
    fn from(v: fs_material::MaterialError) -> Self {
        Self::Fluid(v.into())
    }
}
impl core::fmt::Display for RectangularThermalWallError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rectangular thermal wall: {self:?}")
    }
}
impl std::error::Error for RectangularThermalWallError {}

/// Original physical gas fields, using peak exp(-i omega t) phasors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectangularThermalField {
    pub coordinates_m: [f64; 2],
    pub pressure_pa: C64,
    pub axial_velocity_m_s: C64,
    pub transverse_velocity_gradient_per_s: [C64; 2],
    pub gas_temperature_k: C64,
    pub transverse_temperature_gradient_k_m: [C64; 2],
    pub viscous_dissipation_w_m3: f64,
    pub conduction_dissipation_w_m3: f64,
    pub molecular_dissipation_w_m3: f64,
}

/// Wall quantities per actual surface area, never a volume-heating assignment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectangularThermalWallObservation {
    pub coordinates_m: [f64; 2],
    pub temperature_k: C64,
    /// Original Robin law for finite/insulated walls. At an isothermal wall
    /// this is the one-sided Q1 gradient flux and converges by mesh refinement;
    /// it is not a redistributed discrete corner reaction.
    pub fluid_to_wall_heat_flux_w_m2: C64,
    /// Flux from the one-sided reconstructed gas gradient, for independent
    /// local refinement. It need not equal the weak Robin law pointwise.
    pub reconstructed_gas_heat_flux_w_m2: C64,
    pub wall_to_bath_heat_flux_w_m2: C64,
    pub bath_dissipation_w_m2: f64,
    pub thermal_storage_j_m2: f64,
}

/// Independent positive original physical losses integrated over the section.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RectangularThermalDissipation {
    pub viscous_w_m: f64,
    pub conduction_w_m: f64,
    pub molecular_w_m: f64,
    pub walls_w_m: [f64; 4],
    pub total_w_m: f64,
}

/// Source-bound prepared gas/thermal field; no material query at observation.
/// Preparation uses original Fourier modes plus bounded sparse factorization.
/// A point query costs O(Fourier terms) for shear and O(1) for Q1 thermal fields.
/// Integrated queries are O(1), and all observations allocate nothing.
#[derive(Debug)]
pub struct PreparedRectangularThermalWalls<'g, 'm, M: ?Sized> {
    section: &'g RectangularSection,
    source: &'m M,
    state: ThermoviscousFluidState,
    omega: f64,
    walls: [RectangularThermalWall; 4],
    wave: SegmentWave,
    viscous: DiffusionField,
    thermal: ThermalField,
    velocity_scale: C64,
    temperature_scale: C64,
    loss_roots: [f64; 3],
    integrated_loss_roots: [f64; 7],
    wall_storage_roots: [f64; 4],
    wall_coefficients: [C64; 4],
    background_temperature_k: f64,
    relative_thermal_work_defect: f64,
    peak_payload_bytes: usize,
}

impl RectangularSection {
    pub fn prepare_thermal_wall_fields<'g, 'm>(
        &'g self,
        gas: &'m GasState,
        omega: f64,
        walls: [RectangularThermalWall; 4],
        admission: RectangularAdmission,
        limits: RectangularThermalLimits,
    ) -> Result<PreparedRectangularThermalWalls<'g, 'm, GasState>, RectangularThermalWallError>
    {
        self.prepare_thermal_wall_fields_controlled(gas, omega, walls, admission, limits, || false)
    }
    pub fn prepare_thermal_wall_fields_controlled<'g, 'm>(
        &'g self,
        gas: &'m GasState,
        omega: f64,
        walls: [RectangularThermalWall; 4],
        admission: RectangularAdmission,
        limits: RectangularThermalLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRectangularThermalWalls<'g, 'm, GasState>, RectangularThermalWallError>
    {
        poll(&mut cancelled, 0)?;
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
            limits,
            &mut cancelled,
        )
    }
    pub fn prepare_molecular_thermal_wall_fields<'g, 'm>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        walls: [RectangularThermalWall; 4],
        admission: RectangularAdmission,
        limits: RectangularThermalLimits,
    ) -> Result<PreparedRectangularThermalWalls<'g, 'm, RelaxationState>, RectangularThermalWallError>
    {
        self.prepare_molecular_thermal_wall_fields_controlled(
            material,
            omega,
            walls,
            admission,
            limits,
            || false,
        )
    }
    pub fn prepare_molecular_thermal_wall_fields_controlled<'g, 'm>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        walls: [RectangularThermalWall; 4],
        admission: RectangularAdmission,
        limits: RectangularThermalLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRectangularThermalWalls<'g, 'm, RelaxationState>, RectangularThermalWallError>
    {
        poll(&mut cancelled, 0)?;
        let gas = material.gas();
        let state = gas.thermoviscous_state()?;
        let cp = if material.processes().is_empty() {
            C64::from_re(state.transport().specific_heat_j_kg_k())
        } else {
            material.caloric_response(omega)?.cp_j_kg_k
        };
        prepare(
            self,
            material,
            gas,
            state,
            cp,
            omega,
            walls,
            admission,
            limits,
            &mut cancelled,
        )
    }
}
#[allow(clippy::too_many_arguments)]
fn prepare<'g, 'm, M: ?Sized>(
    section: &'g RectangularSection,
    source: &'m M,
    gas: &GasState,
    state: ThermoviscousFluidState,
    cp: C64,
    omega: f64,
    walls: [RectangularThermalWall; 4],
    admission: RectangularAdmission,
    limits: RectangularThermalLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedRectangularThermalWalls<'g, 'm, M>, RectangularThermalWallError> {
    if !positive(omega)
        || !positive(admission.maximum_relative_tail)
        || admission.maximum_relative_tail >= 1.0
        || !positive(admission.maximum_transverse_phase)
        || admission.maximum_transverse_phase >= core::f64::consts::PI
        || !finite(cp)
        || cp.re <= 0.0
        || cp.im < 0.0
    {
        return Err(bad(
            "rectangular thermal frequency, material or plane-pressure limits are invalid",
        )
        .into());
    }
    let k0 = omega / state.acoustic().sound_speed_m_s();
    if !positive(k0) || k0 * section.long_m > admission.maximum_transverse_phase {
        return Err(bad("rectangular thermal gas exceeds its plane-pressure envelope").into());
    }
    let transport = state.transport();
    let coefficients = [
        walls[0].coefficients(omega)?,
        walls[1].coefficients(omega)?,
        walls[2].coefficients(omega)?,
        walls[3].coefficients(omega)?,
    ];
    let fourier_bytes = section
        .terms()
        .checked_mul(128)
        .and_then(|v| v.checked_add(4096))
        .ok_or_else(|| bad("rectangular shear payload extent overflow"))?;
    let thermal_bytes = limits
        .maximum_payload_bytes
        .checked_sub(fourier_bytes)
        .ok_or_else(|| bad("rectangular thermal payload cannot retain its original shear field"))?;
    // Admit thermal extents before preparing the original Fourier field.
    let conductivity = transport.thermal_conductivity_w_m_k();
    let z = C64::new(cp.im, cp.re).scale(omega * gas.density / conductivity);
    let thermal = physics::prepare(
        [section.short_m, section.long_m],
        z,
        coefficients.map(|v| v.scale(1.0 / conductivity)),
        walls,
        RectangularThermalLimits {
            maximum_payload_bytes: thermal_bytes,
            ..limits
        },
        cancelled,
    )?;
    let diffusion = section.diffusion(omega, transport.kinematic_viscosity_m2_s(), cancelled)?;
    admit_tail(diffusion, admission.maximum_relative_tail)?;
    let viscous = section.prepare_diffusion_fields(
        C64::new(0.0, omega / transport.kinematic_viscosity_m2_s()),
        cancelled,
    )?;
    let velocity_scale = diffusion.mean.conj().recip() / C64::from_re(section.area_m2);
    let temperature_scale = cp.recip() / C64::from_re(gas.density);
    let loss_roots = [
        det::sqrt(transport.dynamic_viscosity_pa_s() * 0.5),
        det::sqrt(conductivity / gas.temperature * 0.5),
        det::sqrt(gas.density) * det::sqrt(omega) * det::sqrt(cp.im / gas.temperature * 0.5),
    ];
    let mut integrated_loss_roots = [0.0; 7];
    integrated_loss_roots[0] = velocity_scale.scale(loss_roots[0]).abs() * viscous.gradient_norm;
    integrated_loss_roots[1] =
        temperature_scale.scale(loss_roots[1]).abs() * det::sqrt(thermal.gradient_form);
    integrated_loss_roots[2] =
        temperature_scale.scale(loss_roots[2]).abs() * det::sqrt(thermal.field_form);
    let mut wall_storage_roots = [0.0; 4];
    for side in 0..4 {
        integrated_loss_roots[3 + side] = temperature_scale.abs()
            * det::sqrt(coefficients[side].re * thermal.wall_forms[side] / (2.0 * gas.temperature));
        wall_storage_roots[side] = temperature_scale.abs()
            * det::sqrt(
                walls[side].capacity() * thermal.wall_forms[side] / (4.0 * gas.temperature),
            );
    }
    let compliance = (C64::ONE - thermal.mean * cp.conj().recip().scale(gas.specific_gas_constant))
        .scale(1.0 / (gas.density * gas.specific_gas_constant * gas.temperature));
    let s = C64::new(0.0, omega);
    let series = s.scale(gas.density / section.area_m2) / diffusion.mean;
    let shunt = s.scale(section.area_m2) * compliance;
    let (wavenumber, characteristic_impedance) =
        crate::helmholtz_wave_from_telegraph(series, shunt);
    let heat = integrated_loss_roots[1..]
        .iter()
        .map(|r| r * r)
        .sum::<f64>();
    let port_heat = 0.5 * shunt.re;
    let relative_thermal_work_defect = if heat == 0.0 && port_heat == 0.0 {
        0.0
    } else {
        (heat - port_heat).abs() / (heat + port_heat.abs()).max(f64::MIN_POSITIVE)
    };
    let wave = SegmentWave {
        wavenumber,
        characteristic_impedance,
        specific_impedance: characteristic_impedance.scale(section.area_m2),
        shear_number: 0.5
            * section.short_m
            * det::sqrt(omega / transport.kinematic_viscosity_m2_s()),
    };
    if ![
        compliance,
        wave.wavenumber,
        wave.characteristic_impedance,
        wave.specific_impedance,
        velocity_scale,
        temperature_scale,
    ]
    .into_iter()
    .all(finite)
        || compliance.re <= 0.0
        || wave.wavenumber.re <= 0.0
        || wave.wavenumber.im < 0.0
        || wave.characteristic_impedance.re <= 0.0
        || k0.max(wave.wavenumber.abs()) * section.long_m > admission.maximum_transverse_phase
        || !relative_thermal_work_defect.is_finite()
        || relative_thermal_work_defect > limits.maximum_relative_residual
        || integrated_loss_roots
            .into_iter()
            .chain(wall_storage_roots)
            .any(|v| !v.is_finite() || v < 0.0)
    {
        return Err(bad(
            "rectangular thermal field, wave or original work balance exceeds admission",
        )
        .into());
    }
    poll(cancelled, section.terms())?;
    let peak_payload_bytes = thermal.peak_payload_bytes + fourier_bytes;
    Ok(PreparedRectangularThermalWalls {
        section,
        source,
        state,
        omega,
        walls,
        wave,
        viscous,
        thermal,
        velocity_scale,
        temperature_scale,
        loss_roots,
        integrated_loss_roots,
        wall_storage_roots,
        wall_coefficients: coefficients,
        background_temperature_k: gas.temperature,
        relative_thermal_work_defect,
        peak_payload_bytes,
    })
}
impl<'g, 'm, M: ?Sized> PreparedRectangularThermalWalls<'g, 'm, M> {
    pub const fn source(&self) -> &'m M {
        self.source
    }
    pub const fn section(&self) -> &'g RectangularSection {
        self.section
    }
    pub const fn state(&self) -> &ThermoviscousFluidState {
        &self.state
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn walls(&self) -> &[RectangularThermalWall; 4] {
        &self.walls
    }
    pub const fn wave(&self) -> &SegmentWave {
        &self.wave
    }
    /// Mean h under the shared +i omega diffusion convention.
    pub fn mean_temperature_response(&self) -> C64 {
        self.thermal.mean
    }
    pub fn thermal_cells(&self) -> [usize; 2] {
        self.thermal.cells
    }
    pub fn relative_backward_residual(&self) -> f64 {
        self.thermal.relative_residual
    }
    pub const fn relative_thermal_work_defect(&self) -> f64 {
        self.relative_thermal_work_defect
    }
    pub const fn peak_payload_bytes(&self) -> usize {
        self.peak_payload_bytes
    }
    pub fn thermal_factor_nonzeros(&self) -> usize {
        self.thermal.factor_nonzeros
    }
    pub fn two_port(&self, length_m: f64, reference: f64) -> Result<TwoPort, DuctError> {
        TwoPort::from_uniform_wave(&self.wave, self.section.area_m2, length_m, reference)
    }
    pub fn field_at(
        &self,
        coordinates_m: [f64; 2],
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<RectangularThermalField, RectangularThermalWallError> {
        self.field_at_controlled(coordinates_m, pressure_pa, volume_flow_m3_s, || false)
    }
    pub fn field_at_controlled(
        &self,
        coordinates_m: [f64; 2],
        pressure_pa: C64,
        volume_flow_m3_s: C64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RectangularThermalField, RectangularThermalWallError> {
        poll(&mut cancelled, 0)?;
        amplitudes(pressure_pa, volume_flow_m3_s)?;
        let t = self
            .thermal
            .at(coordinates_m)?
            .map(|v| v.conj() * self.temperature_scale * pressure_pa);
        let v = self
            .viscous
            .at(
                coordinates_m[0] / self.section.short_m,
                coordinates_m[1],
                self.section.long_m,
                &mut cancelled,
            )?
            .map(|v| v.conj() * self.velocity_scale * volume_flow_m3_s);
        let viscous = v[1..]
            .iter()
            .map(|v| v.scale(self.loss_roots[0]).norm_sq())
            .sum();
        let conduction = t[1..]
            .iter()
            .map(|v| v.scale(self.loss_roots[1]).norm_sq())
            .sum();
        let molecular = t[0].scale(self.loss_roots[2]).norm_sq();
        if !v.into_iter().chain(t).all(finite)
            || [viscous, conduction, molecular]
                .into_iter()
                .any(|v: f64| !v.is_finite())
        {
            return Err(bad(
                "rectangular thermal physical fields or local work are unrepresentable",
            )
            .into());
        }
        poll(&mut cancelled, self.section.terms())?;
        Ok(RectangularThermalField {
            coordinates_m,
            pressure_pa,
            axial_velocity_m_s: v[0],
            transverse_velocity_gradient_per_s: [v[1], v[2]],
            gas_temperature_k: t[0],
            transverse_temperature_gradient_k_m: [t[1], t[2]],
            viscous_dissipation_w_m3: viscous,
            conduction_dissipation_w_m3: conduction,
            molecular_dissipation_w_m3: molecular,
        })
    }
    pub fn wall_observation(
        &self,
        side: usize,
        fraction: f64,
        pressure_pa: C64,
    ) -> Result<RectangularThermalWallObservation, RectangularThermalWallError> {
        if side >= 4
            || !fraction.is_finite()
            || !(0.0..=1.0).contains(&fraction)
            || !finite(pressure_pa)
        {
            return Err(bad(
                "rectangular thermal wall observation needs an original face and finite amplitude",
            )
            .into());
        }
        let point = match side {
            0 => [0.0, self.section.long_m * fraction],
            1 => [self.section.short_m, self.section.long_m * fraction],
            2 => [self.section.short_m * fraction, 0.0],
            _ => [self.section.short_m * fraction, self.section.long_m],
        };
        let t = self
            .thermal
            .at(point)?
            .map(|v| v.conj() * self.temperature_scale * pressure_pa);
        let axis = usize::from(side >= 2);
        let normal = if side % 2 == 0 { -1.0 } else { 1.0 };
        let reconstructed =
            t[1 + axis].scale(-normal * self.state.transport().thermal_conductivity_w_m_k());
        let flux = if self.walls[side] == RectangularThermalWall::Isothermal {
            reconstructed
        } else {
            t[0] * self.wall_coefficients[side].conj()
        };
        let bath = if self.walls[side] == RectangularThermalWall::Isothermal {
            flux
        } else {
            t[0].scale(self.wall_coefficients[side].re)
        };
        let loss = t[0]
            .scale(det::sqrt(
                self.wall_coefficients[side].re / (2.0 * self.background_temperature_k),
            ))
            .norm_sq();
        let storage = t[0]
            .scale(det::sqrt(
                self.walls[side].capacity() / (4.0 * self.background_temperature_k),
            ))
            .norm_sq();
        if ![t[0], reconstructed, flux, bath].into_iter().all(finite)
            || !loss.is_finite()
            || !storage.is_finite()
        {
            return Err(bad("rectangular thermal wall field or work is unrepresentable").into());
        }
        Ok(RectangularThermalWallObservation {
            coordinates_m: point,
            temperature_k: t[0],
            fluid_to_wall_heat_flux_w_m2: flux,
            reconstructed_gas_heat_flux_w_m2: reconstructed,
            wall_to_bath_heat_flux_w_m2: bath,
            bath_dissipation_w_m2: loss,
            thermal_storage_j_m2: storage,
        })
    }
    pub fn dissipation_per_axial_length(
        &self,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<RectangularThermalDissipation, DuctError> {
        amplitudes(pressure_pa, volume_flow_m3_s)?;
        let values = core::array::from_fn::<_, 7, _>(|i| {
            (if i == 0 {
                volume_flow_m3_s
            } else {
                pressure_pa
            })
            .scale(self.integrated_loss_roots[i])
            .norm_sq()
        });
        let total = values.iter().sum();
        if values.iter().any(|v| !v.is_finite()) || !f64::is_finite(total) {
            return Err(bad(
                "rectangular integrated thermal work is unrepresentable",
            ));
        }
        Ok(RectangularThermalDissipation {
            viscous_w_m: values[0],
            conduction_w_m: values[1],
            molecular_w_m: values[2],
            walls_w_m: [values[3], values[4], values[5], values[6]],
            total_w_m: total,
        })
    }
    pub fn wall_storage_per_axial_length(&self, pressure_pa: C64) -> Result<[f64; 4], DuctError> {
        amplitudes(pressure_pa, C64::ZERO)?;
        let value = self
            .wall_storage_roots
            .map(|r| pressure_pa.scale(r).norm_sq());
        if value.iter().any(|v| !v.is_finite()) {
            return Err(bad("rectangular wall thermal storage is unrepresentable"));
        }
        Ok(value)
    }
}
#[cfg(test)]
mod tests;
