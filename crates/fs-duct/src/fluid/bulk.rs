//! Uniform-cylinder acoustic propagation with total-compressive-traction ports.
//!
//! The original material supplies the full homogeneous NSF acoustic mode.
//! The stable two-port owner propagates traction and volume flow; pressure is
//! reconstructed separately. [`BulkTractionCylinder`] keeps only the two
//! traveling acoustic waves and their mechanical power. [`thermal`] also
//! retains thermal diffusion, physical temperature/heat ports and thermal end
//! conditions, with ideal thermal/mechanical contacts between authored regions.
//! Neither model adds no-slip wall loss, transverse modes or radius changes.

use super::{FluidDuctError, bad, wave_from_state};
use crate::compose::{ComposeOptions, TwoPort, compose_chain};
use crate::{DuctError, Segment, SegmentWave};
use fs_material::fluid_acoustics::thermoviscous::absorption::{
    BulkThermoviscousWave, ClassicalBulkAbsorption,
};
use fs_math::c64::C64;

/// Full acoustic AND thermal diffusion propagation, with physical temperature/
/// heat ports, coherent regional contacts and independent thermal end conditions.
pub mod thermal;

mod absorption;
pub use absorption::BulkTractionAbsorptionResponse;

/// Prepared acoustic-mode cylinder. The underlying pressure-named two-port
/// algebra is private: every exposed mechanical quantity here is traction.
#[derive(Clone, Copy, Debug)]
pub struct BulkTractionCylinder {
    wave: BulkThermoviscousWave,
    port: TwoPort,
    pressure_over_traction: C64,
    material: ClassicalBulkAbsorption,
    area_m2: f64,
    length_m: f64,
    omega_rad_s: f64,
}

/// Mechanical load response for unit inlet volume flow, with pressure
/// observables reconstructed from the same acoustic eigenmode.
#[derive(Clone, Copy, Debug)]
pub struct BulkTractionResponse {
    /// Inlet total compressive traction / inlet volume flow [Pa s/m³].
    pub inlet_traction_impedance_pa_s_m3: C64,
    /// Inlet thermodynamic pressure / inlet volume flow [Pa s/m³].
    pub inlet_pressure_over_flow_pa_s_m3: C64,
    /// Outlet total compressive traction / inlet volume flow [Pa s/m³].
    pub outlet_traction_over_inlet_flow_pa_s_m3: C64,
    /// Outlet thermodynamic pressure / inlet volume flow [Pa s/m³].
    pub outlet_pressure_over_inlet_flow_pa_s_m3: C64,
    /// Outlet / inlet rightward volume flow [dimensionless].
    pub outlet_over_inlet_flow: C64,
    pub condition_estimate: f64,
}

impl BulkTractionCylinder {
    /// Freeze the response at positive frequency [rad/s] and radius/length [m].
    /// The supplied material snapshot retains its explicit shear, bulk and
    /// thermal coefficients and acoustic-branch validity. No material query,
    /// wall approximation or attenuation fit is introduced here.
    pub fn new(
        material: &ClassicalBulkAbsorption,
        radius_m: f64,
        length_m: f64,
        omega_rad_s: f64,
    ) -> Result<Self, FluidDuctError> {
        let segment = Segment::Cylinder {
            radius: radius_m,
            length: length_m,
        };
        segment.validate()?;
        let base = wave_from_state(&material.state().acoustic(), radius_m, omega_rad_s)?;
        let wave = material.thermoviscous_wave(omega_rad_s)?;
        let traction = wave.traction_impedance_pa_s_m();
        let area = core::f64::consts::PI * radius_m * radius_m;
        let port = TwoPort::from_wave(
            &segment,
            &SegmentWave {
                wavenumber: wave.wavenumber_rad_m(),
                characteristic_impedance: traction.scale(1.0 / area),
                specific_impedance: traction,
                shear_number: f64::INFINITY,
            },
            base.characteristic_impedance.re,
        )?;
        let pressure_over_traction = wave.pressure_impedance_pa_s_m() / traction;
        if !pressure_over_traction.re.is_finite() || !pressure_over_traction.im.is_finite() {
            return Err(bad("bulk pressure/traction ratio is unrepresentable").into());
        }
        Ok(Self {
            wave,
            port,
            pressure_over_traction,
            material: *material,
            area_m2: area,
            length_m,
            omega_rad_s,
        })
    }

    pub const fn wave(&self) -> BulkThermoviscousWave {
        self.wave
    }

    /// Row-major [S11,S12,S21,S22] for traction/volume-flow waves, with the
    /// positive real reference below. This is not a pressure-port matrix.
    pub const fn traction_scattering(&self) -> [C64; 4] {
        self.port.scattering
    }

    /// Mechanical wave reference [Pa s/m³], equal to the original rho*c/area.
    pub const fn traction_reference_pa_s_m3(&self) -> f64 {
        self.port.reference_impedance
    }

    /// Apply total compressive traction / volume-flow load [Pa s/m³]. `None`
    /// imposes zero outlet flow within this acoustic subspace; it does not
    /// impose a thermal end condition. Pressure loads require explicit
    /// conversion and are not silently reinterpreted as traction loads.
    pub fn response(
        &self,
        traction_load_pa_s_m3: Option<C64>,
        options: ComposeOptions,
    ) -> Result<BulkTractionResponse, DuctError> {
        let result = compose_chain(&[self.port], traction_load_pa_s_m3, options)?;
        let inlet_pressure = self.pressure_over_traction * result.z_in;
        let outlet_pressure = self.pressure_over_traction * result.p_out_over_u_in;
        if [inlet_pressure, outlet_pressure]
            .iter()
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(bad("bulk pressure reconstruction is unrepresentable"));
        }
        Ok(BulkTractionResponse {
            inlet_traction_impedance_pa_s_m3: result.z_in,
            inlet_pressure_over_flow_pa_s_m3: inlet_pressure,
            outlet_traction_over_inlet_flow_pa_s_m3: result.p_out_over_u_in,
            outlet_pressure_over_inlet_flow_pa_s_m3: outlet_pressure,
            outlet_over_inlet_flow: result.u_out_over_u_in,
            condition_estimate: result.condition_estimate,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::fluid_acoustics::{
        FluidAcousticState, thermoviscous::ThermoviscousFluidState,
    };
    use fs_material::fluid_transport::FluidTransportState;
    use fs_math::det;

    fn material(gamma: f64) -> ClassicalBulkAbsorption {
        let state = ThermoviscousFluidState::try_new(
            FluidAcousticState::try_new(2.0, 10.0).unwrap(),
            FluidTransportState::try_new(2.0, 8.0, 0.6, 0.003).unwrap(),
            8.0 / gamma,
        )
        .unwrap();
        ClassicalBulkAbsorption::from_material(&state, 0.002).unwrap()
    }
    fn near(a: C64, b: C64, tolerance: f64) {
        assert!(
            (a - b).abs() <= tolerance * b.abs().max(1e-100),
            "{a:?} != {b:?}"
        );
    }
    fn exponential(k: C64, length: f64) -> C64 {
        let (s, c) = det::sin_cos(k.re * length);
        C64::new(c, s).scale(det::exp(-k.im * length))
    }

    #[test]
    fn g1_bulk_cylinder_matching_preserves_pressure_and_stable_splitting() {
        let law = material(2.0);
        let radius = 0.02;
        let area = core::f64::consts::PI * radius * radius;
        for length in [0.04, 10000.0] {
            let cylinder = BulkTractionCylinder::new(&law, radius, length, 1000.0).unwrap();
            let half = BulkTractionCylinder::new(&law, radius, length * 0.5, 1000.0).unwrap();
            let wave = cylinder.wave();
            let load = wave.traction_impedance_pa_s_m().scale(1.0 / area);
            let response = cylinder
                .response(Some(load), ComposeOptions::default())
                .unwrap();
            let split = compose_chain(
                &[half.port, half.port],
                Some(load),
                ComposeOptions::default(),
            )
            .unwrap();
            let transmitted = exponential(wave.wavenumber_rad_m(), length);
            near(response.inlet_traction_impedance_pa_s_m3, load, 2e-14);
            near(
                response.inlet_pressure_over_flow_pa_s_m3,
                wave.pressure_impedance_pa_s_m().scale(1.0 / area),
                2e-14,
            );
            near(response.outlet_over_inlet_flow, transmitted, 2e-14);
            near(
                response.outlet_pressure_over_inlet_flow_pa_s_m3,
                response.inlet_pressure_over_flow_pa_s_m3 * transmitted,
                2e-14,
            );
            near(response.inlet_traction_impedance_pa_s_m3, split.z_in, 2e-14);
            near(
                response.outlet_over_inlet_flow,
                split.u_out_over_u_in,
                2e-14,
            );
        }
    }

    #[test]
    fn g1_bulk_cylinder_viscous_heating_equals_mechanical_work_loss() {
        let (radius, length, omega) = (0.02, 0.04, 1000.0);
        let area = core::f64::consts::PI * radius * radius;
        let cylinder = BulkTractionCylinder::new(&material(1.0), radius, length, omega).unwrap();
        let load = C64::new(1.4, 0.3).scale(cylinder.traction_reference_pa_s_m3());
        for boundary in [Some(load), None] {
            let response = cylinder
                .response(boundary, ComposeOptions::default())
                .unwrap();
            let k = cylinder.wave().wavenumber_rad_m();
            let t = response.inlet_traction_impedance_pa_s_m3;
            let velocity = C64::from_re(1.0 / area);
            let a = (velocity + t / cylinder.wave().traction_impedance_pa_s_m()).scale(0.5);
            let b = (velocity - t / cylinder.wave().traction_impedance_pa_s_m()).scale(0.5);
            // Exact integral of eta*A*|du/dx|²/2, including interference of
            // forward/backward waves; independent of scattering or a loss fit.
            let forward = -det::expm1(-2.0 * k.im * length) / (2.0 * k.im);
            let backward = det::expm1(2.0 * k.im * length) / (2.0 * k.im);
            let cross = (exponential(C64::from_re(2.0 * k.re), length) - C64::ONE)
                / C64::new(0.0, 2.0 * k.re);
            let heat = 0.5
                * 0.006
                * area
                * k.norm_sq()
                * (a.norm_sq() * forward + b.norm_sq() * backward
                    - 2.0 * (a * b.conj() * cross).re);
            let outlet_work = (response.outlet_traction_over_inlet_flow_pa_s_m3
                * response.outlet_over_inlet_flow.conj())
            .re;
            let lost_work = 0.5 * (t.re - outlet_work);
            assert!(heat > 0.0 && (lost_work / heat - 1.0).abs() < 2e-12);
        }
    }

    #[test]
    fn g0_bulk_cylinder_rejects_invalid_geometry_frequency_and_load() {
        let law = material(2.0);
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
            assert!(BulkTractionCylinder::new(&law, bad, 0.04, 1000.0).is_err());
            assert!(BulkTractionCylinder::new(&law, 0.02, bad, 1000.0).is_err());
            assert!(BulkTractionCylinder::new(&law, 0.02, 0.04, bad).is_err());
        }
        let cylinder = BulkTractionCylinder::new(&law, 0.02, 0.04, 1000.0).unwrap();
        assert!(
            cylinder
                .response(Some(C64::new(f64::NAN, 0.0)), ComposeOptions::default())
                .is_err()
        );
    }
}
