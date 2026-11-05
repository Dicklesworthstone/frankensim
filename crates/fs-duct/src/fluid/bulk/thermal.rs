//! Full one-dimensional Navier--Stokes--Fourier propagation at constant area.
//!
//! Both acoustic and thermal diffusion modes are retained, in exp(-i omega t).
//! The conjugate ports are total compressive traction/volume flow and
//! normalized temperature theta=T'/T0 paired with heat power. Their rightward
//! energy flux is `Re(t conj(U) + theta conj(Q))/2`. Thermal end conditions remain
//! separate from mechanical loads. No lateral wall layers, transverse modes,
//! nonlinear thermodynamics or molecular relaxation are inferred. Explicit
//! original molecular reservoirs enter through
//! [`BulkThermalCylinder::from_relaxation`]. Authored
//! homogeneous regions may share ideal mechanical/thermal contacts through
//! [`BulkThermalCylinder::from_regions`]. Radius supplies area only: this is
//! planar/slip-wall propagation, not a no-slip viscothermal bore model.

use super::super::{bad, wave_from_state};
use super::BulkTractionResponse;
use crate::compose::{ComposeOptions, TwoPort, compose_chain, redheffer_star};
use crate::{DuctError, Segment};
use fs_la::eigen_complex::lu_complex;
use fs_material::fluid_acoustics::thermoviscous::absorption::BulkThermoviscousEquilibrium;
use fs_math::{c64::C64, det};

mod regional;
pub use regional::BulkThermalRegion;
mod molecular;

/// Homogeneous thermal condition at one axial end, independent of its load.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThermalEndCondition {
    /// Fixed reference temperature: T'=0, with heat exchange allowed.
    Isothermal,
    /// Insulated end: normal heat flux=0, with temperature allowed to vary.
    Adiabatic,
}
impl ThermalEndCondition {
    fn reflection(self) -> f64 {
        match self {
            Self::Isothermal => -1.0,
            Self::Adiabatic => 1.0,
        }
    }
}

/// An original material, duct or exponential owner refused; no fallback is selected.
#[derive(Debug, Clone, PartialEq)]
pub enum BulkThermalError {
    Material(fs_material::MaterialError),
    Duct(DuctError),
    Exponential(fs_phs::PhsError),
}
impl From<fs_material::MaterialError> for BulkThermalError {
    fn from(error: fs_material::MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<DuctError> for BulkThermalError {
    fn from(error: DuctError) -> Self {
        Self::Duct(error)
    }
}
impl From<fs_phs::PhsError> for BulkThermalError {
    fn from(error: fs_phs::PhsError) -> Self {
        Self::Exponential(error)
    }
}
impl core::fmt::Display for BulkThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "{error}"),
            Self::Duct(error) => write!(f, "{error}"),
            Self::Exponential(error) => write!(f, "FS-DUCT-BULK-THERMAL: {error}"),
        }
    }
}
impl std::error::Error for BulkThermalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material(error) => Some(error),
            Self::Duct(error) => Some(error),
            Self::Exponential(error) => Some(error),
        }
    }
}

/// Cached complete four-field cylinder and its thermal-boundary reduction.
/// The private mechanical two-port always carries traction, never pressure.
#[derive(Clone, Debug)]
pub struct BulkThermalCylinder {
    network: ThermalNetwork,
    port: TwoPort,
    thermal_drive: [C64; 4],
    ends: [ThermalEndCondition; 2],
}

// A complete operator BEFORE thermal termination. Different regional endpoints
// retain their own constitutive pressure maps and the references used to scale
// those maps; a numerical wave-reference change never changes the material.
#[derive(Clone, Debug)]
struct ThermalNetwork {
    scattering: [C64; 16],
    references: [f64; 2],
    temperature_k: f64,
    endpoints: [PressureEndpoint; 2],
    condition: f64,
}

#[derive(Clone, Copy, Debug)]
struct PressureEndpoint {
    references: [f64; 2],
    coefficients: [C64; 2],
}
impl PressureEndpoint {
    fn pressure(self, traction: C64, theta: C64) -> C64 {
        let rm = det::sqrt(self.references[0]);
        let rt = det::sqrt(self.references[1]);
        (self.coefficients[0] * traction.scale(1.0 / rm)
            + self.coefficients[1] * theta.scale(1.0 / rt))
        .scale(rm)
    }
}

/// Endpoint response to unit inlet volume flow. Heat is positive rightward at
/// BOTH ends. All values are complex small-signal amplitudes, not mean heating.
#[derive(Clone, Copy, Debug)]
pub struct BulkThermalResponse {
    pub mechanical: BulkTractionResponse,
    /// Inlet T' / inlet volume flow [K s/m³].
    pub inlet_temperature_over_inlet_flow_k_s_m3: C64,
    /// Outlet T' / inlet volume flow [K s/m³].
    pub outlet_temperature_over_inlet_flow_k_s_m3: C64,
    /// Inlet rightward heat power / inlet volume flow [Pa].
    pub inlet_heat_power_over_inlet_flow_pa: C64,
    /// Outlet rightward heat power / inlet volume flow [Pa].
    pub outlet_heat_power_over_inlet_flow_pa: C64,
}

impl ThermalNetwork {
    fn new(
        material: &BulkThermoviscousEquilibrium,
        radius_m: f64,
        length_m: f64,
        omega_rad_s: f64,
        maximum_condition: f64,
    ) -> Result<Self, BulkThermalError> {
        Self::with_caloric_response(
            material,
            radius_m,
            length_m,
            omega_rad_s,
            maximum_condition,
            None,
            &mut || false,
        )
    }

    // A molecular gas changes only the heat equation's dynamic Cp. The real
    // equilibrium owns EOS, transport, references and pressure reconstruction.
    fn with_caloric_response(
        material: &BulkThermoviscousEquilibrium,
        radius_m: f64,
        length_m: f64,
        omega_rad_s: f64,
        maximum_condition: f64,
        molecular_cp_j_kg_k: Option<C64>,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Self, BulkThermalError> {
        poll(cancelled, 0)?;
        Segment::Cylinder {
            radius: radius_m,
            length: length_m,
        }
        .validate()?;
        if !maximum_condition.is_finite() || maximum_condition < 1.0 {
            return Err(bad("thermal cylinder needs a finite condition bound >=1").into());
        }
        let fluid = material.absorption().state();
        let acoustic = fluid.acoustic();
        let transport = fluid.transport();
        let mechanical = wave_from_state(&acoustic, radius_m, omega_rad_s)?;
        let area = core::f64::consts::PI * radius_m * radius_m;
        let c = acoustic.sound_speed_m_s();
        let chi = transport.thermal_diffusivity_m2_s();
        let thermal_rate = det::sqrt(omega_rad_s) / det::sqrt(chi);
        let thermal_reference = 1.0
            / area
            / material.temperature_k()
            / transport.thermal_conductivity_w_m_k()
            / thermal_rate;
        let references = [mechanical.characteristic_impedance.re, thermal_reference];
        let nu = (4.0 / 3.0) * transport.kinematic_viscosity_m2_s()
            + material.absorption().bulk_viscosity_pa_s() / acoustic.density_kg_m3();
        let x = omega_rad_s * (nu / c) / c;
        let gamma = fluid.heat_capacity_ratio();
        let d = C64::new(1.0, -gamma * x);
        // Signed coupling is read from the original thermal state, never
        // recovered as sqrt(gamma-1). The equilibrium owner already checks
        // alpha²*c²*T0/Cp = gamma-1 to floating-point accuracy.
        let coupling = material.isobaric_expansion_per_k()
            * c
            * det::sqrt(material.temperature_k() / transport.specific_heat_j_kg_k())
            * det::sqrt(c / (det::sqrt(chi) * det::sqrt(omega_rad_s)));
        if references.iter().any(|r| !r.is_finite() || *r <= 0.0)
            || !thermal_rate.is_finite()
            || thermal_rate <= 0.0
            || !x.is_finite()
            || !d.im.is_finite()
            || !coupling.is_finite()
        {
            return Err(bad("thermal cylinder reference scales are unrepresentable").into());
        }
        // y=[t/sqrt(Rm), theta/sqrt(Rt), U*sqrt(Rm), Q*sqrt(Rt)].
        // Momentum, Fourier conduction, continuity and enthalpy give y'=G y.
        // The symmetric lower-left block uses the original thermodynamic
        // identity to avoid subtracting large viscous terms in the heat row.
        let mut generator = [C64::ZERO; 16];
        generator[2] = C64::new(0.0, omega_rad_s / c);
        generator[7] = C64::from_re(-thermal_rate);
        generator[8] = (generator[2] / d).scale(gamma);
        generator[9] = -(generator[2] / d).scale(coupling);
        generator[12] = generator[9];
        generator[13] = (C64::new(1.0, -x) / d) * C64::new(0.0, thermal_rate);
        if let Some(dynamic_cp) = molecular_cp_j_kg_k {
            let cp = transport.specific_heat_j_kg_k();
            let relative = C64::new(dynamic_cp.re / cp, dynamic_cp.im / cp);
            // Q' gains -s*rho*A*T0*delta_Cp*theta, s=-i*omega. In the
            // original real power coordinates its coefficient is
            // i*sqrt(omega/chi)*delta_Cp/Cp_eq. Evaluate the combined row
            // without subtracting nearly equal equilibrium/frozen Cp values.
            // No EOS term changes; (1-i*x)/d = 1+i*x*(gamma-1)/d.
            generator[13] = C64::new(0.0, thermal_rate)
                * (relative + (C64::new(0.0, x) / d).scale(gamma - 1.0));
        }
        let norm = matrix_norm(&generator, 4);
        if generator.iter().any(|&z| !finite(z)) || !norm.is_finite() || norm <= 0.0 {
            return Err(bad("thermal cylinder generator is unrepresentable").into());
        }
        let mut span = length_m;
        let mut doublings = 0;
        while norm * span > 0.25 {
            if doublings == 48 {
                return Err(bad("thermal cylinder exceeds 48 bounded span doublings").into());
            }
            span *= 0.5;
            doublings += 1;
        }
        let transfer = exponential(&generator, span)?;
        poll(cancelled, 1)?;
        let (mut scattering, mut condition) = transfer_scattering(&transfer, maximum_condition)?;
        for _ in 0..doublings {
            let (next, estimate) = redheffer_star(
                &scattering,
                &scattering,
                2,
                maximum_condition,
                &mut *cancelled,
            )?;
            scattering.copy_from_slice(&next);
            condition = condition.max(estimate);
        }
        let endpoint = PressureEndpoint {
            references,
            // p/sqrt(Rm) = (y0-i*x*coupling*y1)/d. Divide BEFORE
            // multiplying so large viscosity cannot erase a small pressure.
            coefficients: [d.recip(), (C64::new(0.0, -x) / d).scale(coupling)],
        };
        poll(cancelled, doublings + 2)?;
        Ok(Self {
            scattering,
            references,
            temperature_k: material.temperature_k(),
            endpoints: [endpoint; 2],
            condition,
        })
    }
}

impl BulkThermalCylinder {
    /// Prepare from one checked equilibrium, in SI, at positive frequency.
    /// A power-scaled short-span generator is exponentiated once. At most 48
    /// native Redheffer doublings retain both modes without a growing
    /// full-length transfer matrix. `maximum_condition` bounds every local
    /// inverse and feedback cancellation; nonfinite or unresolved work refuses.
    pub fn new(
        material: &BulkThermoviscousEquilibrium,
        radius_m: f64,
        length_m: f64,
        omega_rad_s: f64,
        ends: [ThermalEndCondition; 2],
        maximum_condition: f64,
    ) -> Result<Self, BulkThermalError> {
        let network =
            ThermalNetwork::new(material, radius_m, length_m, omega_rad_s, maximum_condition)?;
        Ok(Self::from_network(network, ends, maximum_condition)?)
    }

    fn from_network(
        mut network: ThermalNetwork,
        ends: [ThermalEndCondition; 2],
        maximum_condition: f64,
    ) -> Result<Self, DuctError> {
        let (port, thermal_drive, terminal_condition) = terminate_thermal(
            &network.scattering,
            network.references[0],
            ends,
            maximum_condition,
        )?;
        network.condition = network.condition.max(terminal_condition);
        Ok(Self {
            network,
            port,
            thermal_drive,
            ends,
        })
    }

    /// Full row-major 4x4 scattering, ordered [left mechanical, left thermal,
    /// right mechanical, right thermal]. At either end e=sqrt(R)(a+b) and
    /// inward j=(a-b)/sqrt(R); the rightward flows reverse at the right end.
    pub const fn scattering(&self) -> [C64; 16] {
        self.network.scattering
    }

    /// Positive real references [traction/volume flow (Pa s/m³), theta/heat
    /// power (1/W)]. Thermal effort is T'/T0, not temperature in kelvin.
    pub const fn references(&self) -> [f64; 2] {
        self.network.references
    }

    /// Mechanical traction scattering after the declared thermal ends close.
    pub const fn traction_scattering(&self) -> [C64; 4] {
        self.port.scattering
    }

    /// Apply a traction/volume-flow load [Pa s/m³]; None closes outlet flow.
    /// Thermal endpoints remain as prepared. No eigenmode pressure conversion
    /// is used: pressure follows the original constitutive equations with the
    /// actual endpoint temperature. Numerical owner refusals propagate.
    pub fn response(
        &self,
        traction_load_pa_s_m3: Option<C64>,
        options: ComposeOptions,
    ) -> Result<BulkThermalResponse, DuctError> {
        if self.network.condition > options.maximum_condition {
            return Err(DuctError::Singular);
        }
        let response = compose_chain(&[self.port], traction_load_pa_s_m3, options)?;
        let rm = det::sqrt(self.network.references[0]);
        let rt = det::sqrt(self.network.references[1]);
        let mechanical_incident = [
            (response.z_in.scale(1.0 / rm) + C64::from_re(rm)).scale(0.5),
            (response.p_out_over_u_in.scale(1.0 / rm) - response.u_out_over_u_in.scale(rm))
                .scale(0.5),
        ];
        let mut theta = [C64::ZERO; 2];
        let mut heat = [C64::ZERO; 2];
        for i in 0..2 {
            let incoming = self.thermal_drive[2 * i] * mechanical_incident[0]
                + self.thermal_drive[2 * i + 1] * mechanical_incident[1];
            let r = self.ends[i].reflection();
            theta[i] = incoming.scale((1.0 + r) * rt);
            heat[i] = incoming.scale((1.0 - r) / rt * if i == 0 { 1.0 } else { -1.0 });
        }
        let pressures = [
            self.network.endpoints[0].pressure(response.z_in, theta[0]),
            self.network.endpoints[1].pressure(response.p_out_over_u_in, theta[1]),
        ];
        let temperature = theta.map(|v| v.scale(self.network.temperature_k));
        if pressures
            .iter()
            .chain(&temperature)
            .chain(&heat)
            .any(|&z| !finite(z))
        {
            return Err(bad("thermal cylinder endpoint response is unrepresentable"));
        }
        Ok(BulkThermalResponse {
            mechanical: BulkTractionResponse {
                inlet_traction_impedance_pa_s_m3: response.z_in,
                inlet_pressure_over_flow_pa_s_m3: pressures[0],
                outlet_traction_over_inlet_flow_pa_s_m3: response.p_out_over_u_in,
                outlet_pressure_over_inlet_flow_pa_s_m3: pressures[1],
                outlet_over_inlet_flow: response.u_out_over_u_in,
                condition_estimate: response.condition_estimate.max(self.network.condition),
            },
            inlet_temperature_over_inlet_flow_k_s_m3: temperature[0],
            outlet_temperature_over_inlet_flow_k_s_m3: temperature[1],
            inlet_heat_power_over_inlet_flow_pa: heat[0],
            outlet_heat_power_over_inlet_flow_pa: heat[1],
        })
    }
}

fn poll(cancelled: &mut dyn FnMut() -> bool, completed_steps: usize) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps })
    } else {
        Ok(())
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn matrix_norm(a: &[C64], n: usize) -> f64 {
    a.chunks_exact(n)
        .map(|row| row.iter().map(|z| z.abs()).sum::<f64>())
        .fold(0.0, f64::max)
}

// Only this SHORT transfer is formed; the full cylinder lives in scattering.
fn exponential(generator: &[C64; 16], span: f64) -> Result<[C64; 16], fs_phs::PhsError> {
    let mut real = [0.0; 64];
    for i in 0..4 {
        for j in 0..4 {
            let z = generator[4 * i + j].scale(span);
            real[8 * i + j] = z.re;
            real[8 * i + j + 4] = -z.im;
            real[8 * (i + 4) + j] = z.im;
            real[8 * (i + 4) + j + 4] = z.re;
        }
    }
    let e = fs_phs::linear_exact::exponential(&real, 8)?;
    Ok(core::array::from_fn(|index| {
        C64::new(
            e[(index / 4) * 8 + index % 4],
            e[(index / 4 + 4) * 8 + index % 4],
        )
    }))
}

// Fixed tiny dimensions; all pivoting/solves are the native complex LU owner.
fn inverse(a: &[C64], n: usize, limit: f64) -> Result<(Vec<C64>, f64), DuctError> {
    let factor = lu_complex(a, n).map_err(|_| DuctError::Singular)?;
    let mut inverse = vec![C64::ZERO; n * n];
    for j in 0..n {
        let mut column = vec![C64::ZERO; n];
        column[j] = C64::ONE;
        factor.solve(&mut column);
        for i in 0..n {
            inverse[i * n + j] = column[i];
        }
    }
    let condition = matrix_norm(a, n) * matrix_norm(&inverse, n);
    if inverse.iter().any(|&z| !finite(z)) || !condition.is_finite() || condition > limit {
        return Err(DuctError::Singular);
    }
    Ok((inverse, condition))
}

fn transfer_scattering(e: &[C64; 16], limit: f64) -> Result<([C64; 16], f64), DuctError> {
    let mut incoming = [C64::ZERO; 16];
    let mut outgoing = incoming;
    for i in 0..2 {
        incoming[4 * i + i] = C64::from_re(0.5);
        incoming[4 * i + i + 2] = C64::from_re(0.5);
        outgoing[4 * i + i] = C64::from_re(0.5);
        outgoing[4 * i + i + 2] = C64::from_re(-0.5);
        for j in 0..4 {
            incoming[4 * (i + 2) + j] = (e[4 * i + j] - e[4 * (i + 2) + j]).scale(0.5);
            outgoing[4 * (i + 2) + j] = (e[4 * i + j] + e[4 * (i + 2) + j]).scale(0.5);
        }
    }
    let (inv, condition) = inverse(&incoming, 4, limit)?;
    let s = core::array::from_fn(|index| {
        (0..4).fold(C64::ZERO, |sum, k| {
            sum + outgoing[(index / 4) * 4 + k] * inv[k * 4 + index % 4]
        })
    });
    if s.iter().any(|&z| !finite(z)) {
        return Err(DuctError::Singular);
    }
    Ok((s, condition))
}

fn terminate_thermal(
    s: &[C64; 16],
    reference: f64,
    ends: [ThermalEndCondition; 2],
    limit: f64,
) -> Result<(TwoPort, [C64; 4], f64), DuctError> {
    let mut feedback = [C64::ZERO; 4];
    let mut drive = feedback;
    for i in 0..2 {
        for j in 0..2 {
            let r = ends[i].reflection();
            feedback[2 * i + j] = C64::from_re(if i == j { 1.0 } else { 0.0 })
                - s[(2 * i + 1) * 4 + 2 * j + 1].scale(r);
            drive[2 * i + j] = s[(2 * i + 1) * 4 + 2 * j].scale(r);
        }
    }
    let (inv, mut condition) = inverse(&feedback, 2, limit)?;
    // Include feedback subtraction as well as matrix conditioning: a tiny
    // scalar/diagonal denominator is not well resolved just because cond=1.
    let scale = 1.0
        + (0..2)
            .map(|i| {
                (0..2)
                    .map(|j| s[(2 * i + 1) * 4 + 2 * j + 1].abs())
                    .sum::<f64>()
            })
            .fold(0.0, f64::max);
    condition = condition.max(scale * matrix_norm(&inv, 2));
    if !condition.is_finite() || condition > limit {
        return Err(DuctError::Singular);
    }
    let heat = core::array::from_fn(|index| {
        (0..2).fold(C64::ZERO, |sum, k| {
            sum + inv[(index / 2) * 2 + k] * drive[k * 2 + index % 2]
        })
    });
    let effective = core::array::from_fn(|index| {
        (0..2).fold(s[(index / 2) * 8 + 2 * (index % 2)], |sum, k| {
            sum + s[(index / 2) * 8 + 2 * k + 1] * heat[k * 2 + index % 2]
        })
    });
    Ok((TwoPort::new(effective, reference)?, heat, condition))
}

#[cfg(test)]
mod tests {
    use super::super::BulkTractionCylinder;
    use super::*;
    use ThermalEndCondition::{Adiabatic, Isothermal};
    use fs_material::fluid_acoustics::{
        FluidAcousticState, thermoviscous::ThermoviscousFluidState,
    };
    use fs_material::fluid_transport::{FluidTransportState, ThermalFluidTransportState};

    fn medium(alpha: f64, shear: f64, bulk: f64) -> BulkThermoviscousEquilibrium {
        let transport = FluidTransportState::try_new(2.0, 8.0, 0.6, shear).unwrap();
        let beta_s = 1.0 / 200.0;
        let beta_t = beta_s + alpha * alpha * 12.5 / 16.0;
        let fluid = ThermoviscousFluidState::try_new(
            FluidAcousticState::try_new(2.0, 10.0).unwrap(),
            transport,
            8.0 * beta_s / beta_t,
        )
        .unwrap();
        BulkThermoviscousEquilibrium::from_states(
            fluid,
            ThermalFluidTransportState::try_new(transport, 12.5, Some(alpha)).unwrap(),
            bulk,
        )
        .unwrap()
    }
    fn near(a: C64, b: C64, relative: f64) {
        assert!(
            (a - b).abs() <= relative * b.abs().max(1e-100),
            "{a:?} != {b:?}"
        );
    }
    fn prepare(
        material: &BulkThermoviscousEquilibrium,
        length: f64,
        omega: f64,
        ends: [ThermalEndCondition; 2],
    ) -> BulkThermalCylinder {
        BulkThermalCylinder::new(material, 0.02, length, omega, ends, 1e12).unwrap()
    }

    #[test]
    fn g1_uncoupled_heat_matches_exact_viscous_acoustic_mode() {
        for (shear, bulk) in [(0.003, 0.002), (3000.0, 2000.0)] {
            let material = medium(0.0, shear, bulk);
            let acoustic =
                BulkTractionCylinder::new(material.absorption(), 0.02, 0.04, 100.0).unwrap();
            let load = Some(C64::new(1.4, 0.3).scale(acoustic.traction_reference_pa_s_m3()));
            let reference = acoustic.response(load, ComposeOptions::default()).unwrap();
            for ends in [
                [Isothermal, Isothermal],
                [Adiabatic, Adiabatic],
                [Isothermal, Adiabatic],
            ] {
                let full = prepare(&material, 0.04, 100.0, ends)
                    .response(load, ComposeOptions::default())
                    .unwrap();
                near(
                    full.mechanical.inlet_traction_impedance_pa_s_m3,
                    reference.inlet_traction_impedance_pa_s_m3,
                    2e-10,
                );
                near(
                    full.mechanical.inlet_pressure_over_flow_pa_s_m3,
                    reference.inlet_pressure_over_flow_pa_s_m3,
                    2e-10,
                );
                near(
                    full.mechanical.outlet_pressure_over_inlet_flow_pa_s_m3,
                    reference.outlet_pressure_over_inlet_flow_pa_s_m3,
                    2e-10,
                );
                near(
                    full.mechanical.outlet_over_inlet_flow,
                    reference.outlet_over_inlet_flow,
                    2e-10,
                );
                assert_eq!(full.inlet_temperature_over_inlet_flow_k_s_m3, C64::ZERO);
                assert_eq!(full.outlet_heat_power_over_inlet_flow_pa, C64::ZERO);
            }
        }
    }

    #[test]
    fn g1_thermal_ends_select_isothermal_or_adiabatic_static_compliance() {
        let material = medium(0.08, 0.003, 0.002);
        let area = core::f64::consts::PI * 0.02 * 0.02;
        let (length, omega) = (0.04, 1e-5);
        for ends in [
            [Isothermal, Isothermal],
            [Isothermal, Adiabatic],
            [Adiabatic, Adiabatic],
        ] {
            let response = prepare(&material, length, omega, ends)
                .response(None, ComposeOptions::default())
                .unwrap();
            let modulus = response.mechanical.inlet_traction_impedance_pa_s_m3
                * C64::new(0.0, -omega * area * length);
            let expected = if ends == [Adiabatic, Adiabatic] {
                200.0
            } else {
                100.0
            };
            near(modulus, C64::from_re(expected), 2e-6);
            assert_eq!(response.mechanical.outlet_over_inlet_flow, C64::ZERO);
        }
    }

    // Independent physical fields [t,theta,U,Q]. These are the original
    // dimensional mass, momentum, heat and stress equations, not the normalized
    // generator or a transfer/scattering rearrangement from production.
    fn physical_rhs(y: [C64; 4], alpha: f64, area: f64, omega: f64) -> ([C64; 4], f64) {
        let eta = 0.006;
        let pressure = (y[0] + C64::new(0.0, -omega * eta * alpha * 12.5) * y[1])
            / C64::new(1.0, -omega * eta * 0.01);
        let ux = C64::new(0.0, omega) * (pressure.scale(0.01) - y[1].scale(alpha * 12.5));
        let theta_x = y[3].scale(-1.0 / (0.6 * area * 12.5));
        let derivative = [
            C64::new(0.0, omega * 2.0 / area) * y[2],
            theta_x,
            ux.scale(area),
            C64::new(0.0, omega * area)
                * (y[1].scale(2.0 * 8.0 * 12.5) - pressure.scale(alpha * 12.5)),
        ];
        let heat = 0.5 * area * (eta * ux.norm_sq() + 0.6 * 12.5 * theta_x.norm_sq());
        (derivative, heat)
    }

    #[test]
    fn g1_full_end_response_closes_original_equations_and_positive_heat() {
        let area = core::f64::consts::PI * 0.02 * 0.02;
        let (length, omega) = (0.04, 100.0);
        for alpha in [-0.08, 0.08] {
            for ends in [
                [Isothermal, Isothermal],
                [Isothermal, Adiabatic],
                [Adiabatic, Adiabatic],
            ] {
                let cylinder = prepare(&medium(alpha, 0.003, 0.002), length, omega, ends);
                let load = Some(C64::new(1.4, 0.3).scale(cylinder.references()[0]));
                let result = cylinder.response(load, ComposeOptions::default()).unwrap();
                let m = result.mechanical;
                let inlet = [
                    m.inlet_traction_impedance_pa_s_m3,
                    result
                        .inlet_temperature_over_inlet_flow_k_s_m3
                        .scale(1.0 / 12.5),
                    C64::ONE,
                    result.inlet_heat_power_over_inlet_flow_pa,
                ];
                let outlet = [
                    m.outlet_traction_over_inlet_flow_pa_s_m3,
                    result
                        .outlet_temperature_over_inlet_flow_k_s_m3
                        .scale(1.0 / 12.5),
                    m.outlet_over_inlet_flow,
                    result.outlet_heat_power_over_inlet_flow_pa,
                ];
                for (i, y) in [inlet, outlet].iter().enumerate() {
                    match ends[i] {
                        Isothermal => assert_eq!(y[1], C64::ZERO),
                        Adiabatic => assert_eq!(y[3], C64::ZERO),
                    }
                    let pressure = (y[0] + C64::new(0.0, -omega * 0.006 * alpha * 12.5) * y[1])
                        / C64::new(1.0, -omega * 0.006 * 0.01);
                    near(
                        [
                            m.inlet_pressure_over_flow_pa_s_m3,
                            m.outlet_pressure_over_inlet_flow_pa_s_m3,
                        ][i],
                        pressure,
                        2e-12,
                    );
                }
                // Independent RK4 integration of the primitive equations and
                // Simpson integration of eta*A|u'|²/2 + kappa*T0*A|theta'|²/2.
                let steps = 2048;
                let dx = length / steps as f64;
                let mut y = inlet;
                let mut heat = 0.0;
                for step in 0..=steps {
                    let (k1, source) = physical_rhs(y, alpha, area, omega);
                    heat += source
                        * if step == 0 || step == steps {
                            1.0
                        } else if step % 2 == 0 {
                            2.0
                        } else {
                            4.0
                        };
                    if step == steps {
                        break;
                    }
                    let advance =
                        |k: [C64; 4], h: f64| core::array::from_fn(|i| y[i] + k[i].scale(h));
                    let k2 = physical_rhs(advance(k1, dx * 0.5), alpha, area, omega).0;
                    let k3 = physical_rhs(advance(k2, dx * 0.5), alpha, area, omega).0;
                    let k4 = physical_rhs(advance(k3, dx), alpha, area, omega).0;
                    y = core::array::from_fn(|i| {
                        y[i] + (k1[i] + k2[i].scale(2.0) + k3[i].scale(2.0) + k4[i]).scale(dx / 6.0)
                    });
                }
                let endpoint_scale = [
                    outlet[0].abs(),
                    inlet[1].abs().max(outlet[1].abs()).max(1.0),
                    outlet[2].abs(),
                    inlet[3].abs().max(outlet[3].abs()).max(1.0),
                ];
                for i in 0..4 {
                    assert!(
                        (y[i] - outlet[i]).abs() < 2e-9 * endpoint_scale[i],
                        "field {i}: {:?} != {:?}",
                        y[i],
                        outlet[i]
                    );
                }
                heat *= dx / 3.0;
                let flux = |z: [C64; 4]| 0.5 * (z[0] * z[2].conj() + z[1] * z[3].conj()).re;
                let loss = flux(inlet) - flux(outlet);
                assert!(
                    heat > 0.0 && (loss / heat - 1.0).abs() < 2e-9,
                    "work {loss}, heat {heat}"
                );
            }
        }
    }

    #[test]
    fn g3_full_scattering_stays_stable_when_thermal_and_acoustic_spans_are_opaque() {
        let material = medium(0.08, 0.003, 0.002);
        for length in [0.04, 10000.0] {
            let full = prepare(&material, length, 1000.0, [Isothermal, Adiabatic]);
            let half = prepare(&material, length * 0.5, 1000.0, [Isothermal, Adiabatic]);
            let (split, _) =
                redheffer_star(&half.scattering(), &half.scattering(), 2, 1e12, || false).unwrap();
            for (a, b) in full.scattering().iter().zip(&split) {
                assert!((*a - *b).abs() < 2e-12);
            }
            for i in 0..4 {
                assert!(
                    (0..4)
                        .map(|j| full.scattering()[4 * j + i].norm_sq())
                        .sum::<f64>()
                        <= 1.0 + 2e-11
                );
                for j in 0..4 {
                    assert!(
                        (full.scattering()[4 * i + j] - full.scattering()[4 * j + i]).abs() < 2e-12
                    );
                }
            }
            let positive = full
                .response(
                    Some(C64::from_re(full.references()[0])),
                    ComposeOptions::default(),
                )
                .unwrap();
            let negative = prepare(
                &medium(-0.08, 0.003, 0.002),
                length,
                1000.0,
                [Isothermal, Adiabatic],
            )
            .response(
                Some(C64::from_re(full.references()[0])),
                ComposeOptions::default(),
            )
            .unwrap();
            near(
                positive.mechanical.inlet_traction_impedance_pa_s_m3,
                negative.mechanical.inlet_traction_impedance_pa_s_m3,
                2e-12,
            );
            near(
                positive.inlet_heat_power_over_inlet_flow_pa,
                -negative.inlet_heat_power_over_inlet_flow_pa,
                2e-12,
            );
            assert!(positive.mechanical.inlet_traction_impedance_pa_s_m3.re > 0.0);
        }
    }

    #[test]
    fn g0_full_thermal_cylinder_refuses_invalid_or_unresolved_inputs() {
        let material = medium(0.08, 0.003, 0.002);
        let ends = [Isothermal, Adiabatic];
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
            assert!(BulkThermalCylinder::new(&material, bad, 0.04, 100.0, ends, 1e12).is_err());
            assert!(BulkThermalCylinder::new(&material, 0.02, bad, 100.0, ends, 1e12).is_err());
            assert!(BulkThermalCylinder::new(&material, 0.02, 0.04, bad, ends, 1e12).is_err());
        }
        for limit in [0.0, 1.0, f64::NAN, f64::INFINITY] {
            assert!(BulkThermalCylinder::new(&material, 0.02, 0.04, 100.0, ends, limit).is_err());
        }
    }
}
