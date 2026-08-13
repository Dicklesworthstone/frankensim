//! Certified compact radiators from a thin plate.
//!
//! Geometry + orthotropic section go in; `fs-plate` assembles the DKT
//! pencil and `fs-modal` returns inertia-certified eigenpairs. A
//! guitar top, a bulkhead, and a panel are the same object. Radiation
//! reaction is the baffled-piston small-`ka` series fitted by
//! `fs-vfit`, not a named instrument radiator.

use crate::acoustic_realize::AcousticRealizeError;
use fs_material::elastic::OrthotropicElastic;
use fs_material::gas::GasState;
use fs_material::visco::RayleighDamping;
use fs_math::c64::C64;
use fs_math::det;
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateError, PlateMesh, PlateSection, assemble, modes,
};
use fs_scenario::{RadiatingPlate, ThinPlate};
use fs_vfit::FitOptions;
use fs_vfit::discretize::{DigitalFilter, DigitalFilterState, realize_tabulated_impedance};

/// One driven compact radiator harvested from a certified plate mode.
#[derive(Debug, Clone)]
pub struct CompactBody {
    /// Radiating monopole area [m²] (`∫ φ dA`).
    pub area_m2: f64,
    /// Modal mass [kg] (1 for a mass-normalized mode).
    pub mass_kg: f64,
    /// Angular frequency [rad/s].
    pub omega: f64,
    /// Viscous damping ratio.
    pub zeta: f64,
    y: f64,
    v: f64,
    rad: Option<(DigitalFilter, DigitalFilterState)>,
}

impl CompactBody {
    /// From a caller-supplied compact monopole.
    ///
    /// # Errors
    /// Non-physical parameters.
    pub fn from_radiator(spec: RadiatingPlate) -> Result<Self, AcousticRealizeError> {
        if !(spec.area_m2 > 0.0
            && spec.mass_kg > 0.0
            && spec.frequency_hz > 0.0
            && spec.damping_ratio >= 0.0)
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "compact radiator parameters must be physical",
            });
        }
        Ok(Self {
            area_m2: spec.area_m2,
            mass_kg: spec.mass_kg,
            omega: core::f64::consts::TAU * spec.frequency_hz,
            zeta: spec.damping_ratio,
            y: 0.0,
            v: 0.0,
            rad: None,
        })
    }

    /// Advance under a generalized force and return acceleration.
    pub fn drive(&mut self, force_n: f64, dt: f64) -> f64 {
        let f_rad = if let Some((filter, state)) = self.rad.as_mut() {
            match filter.step(state, self.v) {
                Ok(p_face) => -p_face * self.area_m2,
                Err(_) => 0.0,
            }
        } else {
            0.0
        };
        let acc = (force_n + f_rad) / self.mass_kg
            - 2.0 * self.zeta * self.omega * self.v
            - self.omega * self.omega * self.y;
        self.v += dt * acc;
        self.y += dt * self.v;
        acc
    }

    /// Compact monopole pressure at distance `listener_m`.
    #[must_use]
    pub fn radiate(&self, acc: f64, rho: f64, listener_m: f64) -> f64 {
        rho * self.area_m2 * acc / (4.0 * core::f64::consts::PI * listener_m)
    }

    /// Volume velocity of the monopole [m³/s].
    #[must_use]
    pub fn volume_velocity(&self) -> f64 {
        self.area_m2 * self.v
    }

    /// Drive and radiate in one step.
    pub fn drive_and_radiate(&mut self, force_n: f64, dt: f64, rho: f64, listener_m: f64) -> f64 {
        let acc = self.drive(force_n, dt);
        self.radiate(acc, rho, listener_m)
    }

    fn attach_piston_load(&mut self, gas: &GasState, sample_rate_hz: u32) {
        if self.rad.is_some() {
            return;
        }
        let radius = (self.area_m2 / core::f64::consts::PI).sqrt();
        if let Some(filter) = piston_radiation_filter(radius, gas, sample_rate_hz) {
            let state = filter.zero_state();
            self.rad = Some((filter, state));
        }
    }
}

/// A von Karman simply-supported plate as one pHS (isotropic).
pub struct VkBody {
    sys: fs_phs::PortHamiltonian,
    x: Vec<f64>,
    areas: Vec<f64>,
}

impl VkBody {
    /// Build from an isotropic thin plate.
    ///
    /// # Errors
    /// Non-isotropic section or nlmodal admission.
    pub fn from_plate(plate: ThinPlate) -> Result<Self, AcousticRealizeError> {
        if (plate.e1_pa - plate.e2_pa).abs() > 1.0e-6 * plate.e1_pa.abs() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "von Karman plate is isotropic; e1 must equal e2",
            });
        }
        let n = plate.n_modes.max(1).min(3);
        let disp = odd_odd_modes(n);
        let stress = vec![fs_nlmodal::SineMode { m: 2, n: 1 }];
        let model = fs_nlmodal::von_karman_ss_plate(
            &fs_nlmodal::VkPlateParams {
                lx: plate.length_m,
                ly: plate.width_m,
                h: plate.thickness_m,
                young: plate.e1_pa,
                nu: plate.nu12,
                rho: plate.density_kg_m3,
            },
            &disp,
            &stress,
        )
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let zetas = vec![plate.damping_ratio; disp.len()];
        let mut areas = Vec::with_capacity(disp.len());
        let mut drive = Vec::with_capacity(disp.len());
        let norm = (2.0
            / (plate.density_kg_m3 * plate.thickness_m * plate.length_m * plate.width_m))
            .sqrt();
        for md in &disp {
            let mf = md.m as f64;
            let nf = md.n as f64;
            let pi = core::f64::consts::PI;
            let a = if md.m % 2 == 1 && md.n % 2 == 1 {
                norm * (2.0 * plate.length_m / (mf * pi)) * (2.0 * plate.width_m / (nf * pi))
            } else {
                0.0
            };
            areas.push(a);
            let phi = norm * det::sin(mf * pi * 0.25) * det::sin(nf * pi * 0.5);
            drive.push(phi);
        }
        let sys = fs_nlmodal::assemble(model.storage, &zetas, &drive)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let n_state = 2 * disp.len();
        Ok(Self {
            sys,
            x: vec![0.0; n_state],
            areas,
        })
    }

    fn n_modes(&self) -> usize {
        self.areas.len()
    }

    /// Step under a physical force and return radiated pressure + volume velocity.
    pub fn drive_and_radiate(
        &mut self,
        force_n: f64,
        dt: f64,
        rho: f64,
        listener_m: f64,
    ) -> Result<(f64, f64), AcousticRealizeError> {
        let rec = fs_phs::step(&self.sys, &self.x, &[force_n], dt)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let mut p = 0.0;
        let mut u = 0.0;
        for k in 0..self.n_modes() {
            let v0 = self.x[2 * k + 1];
            let v1 = rec.x[2 * k + 1];
            let acc = (v1 - v0) / dt;
            p += rho * self.areas[k] * acc / (4.0 * core::f64::consts::PI * listener_m);
            u += self.areas[k] * v1;
        }
        self.x = rec.x;
        Ok((p, u))
    }

    /// Modal volume velocity of the von Karman bank [m³/s].
    #[must_use]
    pub fn volume_velocity(&self) -> f64 {
        (0..self.n_modes())
            .map(|k| self.areas[k] * self.x[2 * k + 1])
            .sum()
    }
}

fn odd_odd_modes(n: usize) -> Vec<fs_nlmodal::SineMode> {
    let mut out = Vec::new();
    for sum in (2..20).step_by(2) {
        for m in (1..sum).step_by(2) {
            let nn = sum - m;
            if nn % 2 == 1 {
                out.push(fs_nlmodal::SineMode { m, n: nn });
                if out.len() == n {
                    return out;
                }
            }
        }
    }
    out
}

/// Linear compact radiators plus an optional von Karman pHS.
#[derive(Default)]
pub struct PlateBank {
    /// Linear modal monopoles.
    pub linear: Vec<CompactBody>,
    /// Isotropic von Karman plate, if requested.
    pub vk: Option<VkBody>,
}

impl PlateBank {
    /// True when nothing will radiate or load the waveguide.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.linear.is_empty() && self.vk.is_none()
    }

    /// Volume velocity of every radiator.
    #[must_use]
    pub fn volume_velocity(&self) -> f64 {
        let linear: f64 = self.linear.iter().map(CompactBody::volume_velocity).sum();
        linear + self.vk.as_ref().map_or(0.0, VkBody::volume_velocity)
    }

    /// Drive every radiator and sum compact-monopole pressures.
    pub fn drive_and_radiate(
        &mut self,
        force_n: f64,
        dt: f64,
        rho: f64,
        listener_m: f64,
    ) -> Result<f64, AcousticRealizeError> {
        let mut p = 0.0;
        for body in &mut self.linear {
            p += body.drive_and_radiate(force_n, dt, rho, listener_m);
        }
        if let Some(vk) = &mut self.vk {
            p += vk.drive_and_radiate(force_n, dt, rho, listener_m)?.0;
        }
        Ok(p)
    }

    /// Fit a baffled-piston radiation impedance onto every linear
    /// compact radiator. A failed identification leaves that body
    /// unloaded — it does not invent a damper.
    pub fn attach_radiation_loads(&mut self, gas: &GasState, sample_rate_hz: u32) {
        for body in &mut self.linear {
            body.attach_piston_load(gas, sample_rate_hz);
        }
    }
}

/// Certified compact radiators of a thin plate.
///
/// # Errors
/// Section, mesh, or modal-window refusals.
pub fn certified_radiators(plate: ThinPlate) -> Result<Vec<CompactBody>, AcousticRealizeError> {
    if plate.n_modes == 0 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "thin plate needs at least one mode",
        });
    }
    let law = OrthotropicElastic::new(
        [plate.e1_pa, plate.e2_pa, plate.e2_pa],
        [plate.nu12, plate.nu12, plate.nu12],
        [plate.g12_pa, plate.g12_pa, plate.g12_pa],
        1.0,
    )
    .map_err(|_| AcousticRealizeError::InvalidDescription {
        what: "plate elastic constants refused",
    })?;
    let section = PlateSection::orthotropic(&law, plate.thickness_m, plate.density_kg_m3)
        .map_err(map_plate)?;
    let (nx, ny) = (5, 4);
    let mesh = PlateMesh::rectangle(plate.length_m, plate.width_m, nx, ny);
    let boundary = PlateMesh::rectangle_boundary(nx, ny);
    let model = assemble(
        &mesh,
        &section,
        &boundary,
        &[],
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::SimplySupported,
        },
    )
    .map_err(map_plate)?;
    let omega11 = ss_omega11(&section, plate.length_m, plate.width_m);
    let lo = (0.25 * omega11).powi(2);
    let hi = (omega11 * (plate.n_modes as f64 + 2.0) * 4.0).powi(2);
    let report = modes(&model, (lo, hi), &fs_modal::SliceOptions::default()).map_err(map_plate)?;
    if report.modes.is_empty() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "plate modal window returned no certified modes",
        });
    }
    let n_keep = plate.n_modes.min(report.modes.len());
    let nn = mesh.node_count();
    let nodal_area = plate.length_m * plate.width_m / nn as f64;
    let mut out = Vec::with_capacity(n_keep);
    for pair in report.modes.iter().take(n_keep) {
        let omega = pair.lambda.max(0.0).sqrt();
        if !(omega > 0.0) {
            continue;
        }
        let mut monopole = 0.0;
        let mut drive = 0.0;
        let mut drive_w = 0.0;
        let mid_y = 0.5 * plate.width_m;
        for (inode, &(x, y)) in mesh.nodes.iter().enumerate() {
            let Some(reduced) = model.dof_map[3 * inode] else {
                continue;
            };
            let phi = pair.phi.get(reduced).copied().unwrap_or(0.0);
            monopole += phi * nodal_area;
            if x <= plate.length_m / (nx as f64).max(1.0) {
                let w = 1.0 / (1.0 + (y - mid_y).abs());
                drive += phi * w;
                drive_w += w;
            }
        }
        let phi_drive = if drive_w > 0.0 { drive / drive_w } else { 1.0 };
        let mass = if phi_drive.abs() > 1.0e-12 {
            1.0 / phi_drive.abs()
        } else {
            1.0
        };
        out.push(CompactBody {
            area_m2: monopole.abs().max(1.0e-8),
            mass_kg: mass,
            omega,
            zeta: plate.damping_ratio,
            y: 0.0,
            v: 0.0,
            rad: None,
        });
    }
    if out.len() >= 2 {
        let w0 = out[0].omega;
        // Same authored ratio at ω0 and 4ω0: the Rayleigh bowl then
        // splits the higher certified modes. Pinning at the first and
        // last kept frequencies would assign them identical zetas.
        if let Ok(rayleigh) =
            RayleighDamping::from_two_points(w0, plate.damping_ratio, 4.0 * w0, plate.damping_ratio)
        {
            for body in &mut out {
                body.zeta = rayleigh.zeta_at(body.omega);
            }
        }
    }
    if out.is_empty() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "no usable plate radiators after harvesting",
        });
    }
    Ok(out)
}

/// Baffled-piston `Z(ω) = p/v` fitted as a passive discrete filter.
fn piston_radiation_filter(
    radius: f64,
    gas: &GasState,
    sample_rate_hz: u32,
) -> Option<DigitalFilter> {
    if !(radius > 0.0 && sample_rate_hz > 0) {
        return None;
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    let nyquist = core::f64::consts::PI / dt;
    let omega_lo = 40.0_f64.min(0.05 * nyquist);
    let omega_hi = (0.40 * nyquist).max(omega_lo * 4.0);
    if !(omega_hi > omega_lo) {
        return None;
    }
    let n = 24usize;
    // Small-ka baffled-piston series (the fs-bem Rayleigh-integral
    // oracle). fs-bem cannot be a production dep of this crate —
    // couple → bem → solver → feec → couple. The series is the same
    // object: z/(ρc) = (ka)²/2 − i 8ka/(3π) on e^{-iωt}.
    let mut omega = Vec::with_capacity(n);
    let mut z = Vec::with_capacity(n);
    let zc = gas.density * gas.sound_speed;
    for k in 0..n {
        let t = k as f64 / (n as f64 - 1.0);
        let w = omega_lo * det::exp(t * det::ln(omega_hi / omega_lo));
        let ka = w * radius / gas.sound_speed;
        let re = zc * 0.5 * ka * ka;
        let im = -zc * 8.0 * ka / (3.0 * core::f64::consts::PI);
        // Acoustics `e^{-iωt}` → vfit `e^{+iωt}`.
        omega.push(w);
        z.push(C64::new(re, -im));
    }
    if omega.len() < 8 {
        return None;
    }
    let mut opts = FitOptions::new(3);
    opts.fit_e = false;
    opts.iterations = 8;
    realize_tabulated_impedance(&omega, &z, dt, &opts, omega[omega.len() / 2]).ok()
}

fn ss_omega11(section: &PlateSection, a: f64, b: f64) -> f64 {
    let pi = core::f64::consts::PI;
    let d = section.d[0].abs().max(1.0e-18);
    let rho_h = (section.density * section.thickness).max(1.0e-18);
    pi * pi * (1.0 / (a * a) + 1.0 / (b * b)) * (d / rho_h).sqrt()
}

fn map_plate(_err: PlateError) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription {
        what: "thin-plate modal solve refused",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simply_supported_first_mode_is_near_the_analytic_value() {
        let plate = ThinPlate {
            length_m: 0.20,
            width_m: 0.15,
            thickness_m: 0.002,
            density_kg_m3: 7800.0,
            e1_pa: 200e9,
            e2_pa: 200e9,
            nu12: 0.3,
            g12_pa: 200e9 / (2.0 * 1.3),
            damping_ratio: 0.01,
            n_modes: 2,
            geometric_nonlinearity: false,
        };
        let bodies = certified_radiators(plate).expect("modes");
        let section = PlateSection::isotropic(200e9, 0.3, 0.002, 7800.0).expect("section");
        let want = ss_omega11(&section, 0.20, 0.15);
        let got = bodies[0].omega;
        assert!(
            (got - want).abs() / want < 0.20,
            "certified ω={got:.1} vs SS analytic {want:.1}"
        );
        assert!(
            (bodies[0].zeta - bodies[1].zeta).abs() > 1.0e-12,
            "two-point Rayleigh must split modal zetas"
        );
    }
}
