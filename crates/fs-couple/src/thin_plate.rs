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
        // Baffled half-space (same piston as the self-load), not free-space.
        rho * self.area_m2 * acc / (2.0 * core::f64::consts::PI * listener_m)
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

/// A von Karman plate as one pHS.
///
/// Isotropic simply-supported plates use analytic sine modes.
/// Clamped or orthotropic bending uses DKT-sampled displacement
/// with the same sine Airy membrane channel.
pub struct VkBody {
    sys: fs_phs::PortHamiltonian,
    x: Vec<f64>,
    /// Modal monopole areas [m²].
    pub areas: Vec<f64>,
}

impl VkBody {
    /// Build from a thin plate.
    ///
    /// # Errors
    /// Section, modal-window, or nlmodal admission refusals.
    pub fn from_plate(plate: ThinPlate) -> Result<Self, AcousticRealizeError> {
        Self::from_plate_ports(plate, false)
    }

    fn from_plate_ports(plate: ThinPlate, with_area: bool) -> Result<Self, AcousticRealizeError> {
        let isotropic = (plate.e1_pa - plate.e2_pa).abs() <= 1.0e-6 * plate.e1_pa.abs();
        if isotropic && !plate.clamped {
            Self::from_ss_sine(plate, with_area)
        } else {
            Self::from_sampled_fe(plate, with_area)
        }
    }

    fn from_ss_sine(plate: ThinPlate, with_area: bool) -> Result<Self, AcousticRealizeError> {
        let n = plate.n_modes.clamp(1, 3);
        let disp = odd_odd_modes(n);
        // Extra Airy channels above (2,1) trip the nlmodal quadrature
        // certificate on this mesh; more displacement modes still
        // couple through the one certified membrane channel.
        let stress = vec![fs_nlmodal::SineMode { m: 2, n: 1 }];
        let model = fs_nlmodal::von_karman_ss_plate(
            &fs_nlmodal::VkPlateParams {
                lx: plate.length_m,
                ly: plate.width_m,
                h: plate.thickness_m,
                young: plate.e1_pa,
                nu: plate.nu12,
                rho: plate.density_kg_m3,
                pretension_n_m: plate.pretension_n_m,
            },
            &disp,
            &stress,
        )
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let zetas = vk_zetas(
            plate.density_kg_m3,
            plate.thickness_m,
            plate.damping_ratio,
            &model.storage.omegas,
        );
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
        finish_vk(model.storage, &zetas, &drive, areas, with_area)
    }

    #[allow(clippy::too_many_lines)] // one coherent FE sampling stage
    fn from_sampled_fe(plate: ThinPlate, with_area: bool) -> Result<Self, AcousticRealizeError> {
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
        let (nx_c, ny_c) = (8, 8);
        let mesh = PlateMesh::rectangle(plate.length_m, plate.width_m, nx_c, ny_c);
        let boundary = PlateMesh::rectangle_boundary(nx_c, ny_c);
        let model = assemble(
            &mesh,
            &section,
            &boundary,
            &[],
            &AssemblyOptions {
                pretension: plate.pretension_n_m,
                support: if plate.clamped {
                    EdgeSupport::Clamped
                } else {
                    EdgeSupport::SimplySupported
                },
            },
        )
        .map_err(map_plate)?;
        let omega11 = ss_omega11(&section, plate.length_m, plate.width_m);
        let lo = (0.25 * omega11).powi(2);
        let hi = (omega11 * (plate.n_modes as f64 + 2.0) * 4.0).powi(2);
        let report =
            modes(&model, (lo, hi), &fs_modal::SliceOptions::default()).map_err(map_plate)?;
        if report.modes.is_empty() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "plate modal window returned no certified modes",
            });
        }
        let n_keep = plate.n_modes.min(report.modes.len()).min(3);
        let nx = nx_c + 1;
        let ny = ny_c + 1;
        let mut disp = Vec::with_capacity(n_keep);
        for pair in report.modes.iter().take(n_keep) {
            let omega = pair.lambda.max(0.0).sqrt();
            if !(omega > 0.0) {
                continue;
            }
            let mut w = vec![0.0; nx * ny];
            for j in 0..ny {
                for i in 0..nx {
                    let node = j * nx + i;
                    w[j * nx + i] = match model.dof_map.get(3 * node).copied().flatten() {
                        Some(r) => pair.phi.get(r).copied().unwrap_or(0.0),
                        None => 0.0,
                    };
                }
            }
            disp.push(fs_nlmodal::SampledPlateMode { omega, w, nx, ny });
        }
        if disp.is_empty() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "no usable FE displacement samples",
            });
        }
        let vk = fs_nlmodal::von_karman_sampled_plate(
            &fs_nlmodal::VkPlateParams {
                lx: plate.length_m,
                ly: plate.width_m,
                h: plate.thickness_m,
                young: f64::midpoint(plate.e1_pa, plate.e2_pa),
                nu: plate.nu12,
                rho: plate.density_kg_m3,
                pretension_n_m: plate.pretension_n_m,
            },
            &disp,
            &[
                fs_nlmodal::SineMode { m: 1, n: 1 },
                fs_nlmodal::SineMode { m: 2, n: 2 },
            ],
        )
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let dx = plate.length_m / (nx - 1) as f64;
        let dy = plate.width_m / (ny - 1) as f64;
        let area_el = dx * dy;
        let rho_h = plate.density_kg_m3 * plate.thickness_m;
        let mut areas = Vec::with_capacity(disp.len());
        let mut drive = Vec::with_capacity(disp.len());
        let sx = 0.25 * plate.length_m;
        let sy = 0.5 * plate.width_m;
        for mode in &disp {
            let energy: f64 = mode.w.iter().map(|v| v * v).sum::<f64>() * area_el * rho_h;
            if !(energy > 0.0) {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "sampled FE mode has no L2 mass",
                });
            }
            let scale = 1.0 / energy.sqrt();
            areas.push(mode.w.iter().map(|v| v * scale * area_el).sum());
            drive.push(
                bilinear_sample(&mode.w, nx, ny, plate.length_m, plate.width_m, sx, sy) * scale,
            );
        }
        let zetas = vk_zetas(
            plate.density_kg_m3,
            plate.thickness_m,
            plate.damping_ratio,
            &vk.storage.omegas,
        );
        finish_vk(vk.storage, &zetas, &drive, areas, with_area)
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
            p += rho * self.areas[k] * acc / (2.0 * core::f64::consts::PI * listener_m);
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

fn vk_zetas(density: f64, thickness: f64, damping_ratio: f64, omegas: &[f64]) -> Vec<f64> {
    let mut zetas = vec![damping_ratio; omegas.len()];
    if let Ok(te) = thermoelastic_for_density(density, 293.15) {
        for (z, &w) in zetas.iter_mut().zip(omegas) {
            *z += fs_material::visco::loss_factor_to_zeta(te.loss_factor(w, thickness));
        }
    }
    zetas
}

fn finish_vk(
    storage: fs_nlmodal::SosModalStorage,
    zetas: &[f64],
    drive: &[f64],
    areas: Vec<f64>,
    with_area: bool,
) -> Result<VkBody, AcousticRealizeError> {
    let n = areas.len();
    let m = if with_area { 2 } else { 1 };
    let mut g = vec![0.0; (2 * n) * m];
    for k in 0..n {
        g[(2 * k + 1) * m] = drive[k];
        if with_area {
            g[(2 * k + 1) * m + 1] = areas[k];
        }
    }
    let omegas = storage.omegas.clone();
    let sys = fs_nlmodal::assemble_storage(n, &omegas, zetas, m, g, Box::new(storage))
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    Ok(VkBody {
        sys,
        x: vec![0.0; 2 * n],
        areas,
    })
}

/// Von Karman plate as a 1- or 2-port pHS (bridge, optional face).
///
/// # Errors
/// Same as [`VkBody::from_plate`].
pub fn vk_plate_phs(
    plate: ThinPlate,
    with_area: bool,
) -> Result<fs_phs::PortHamiltonian, AcousticRealizeError> {
    Ok(VkBody::from_plate_ports(plate, with_area)?.sys)
}

fn bilinear_sample(w: &[f64], nx: usize, ny: usize, lx: f64, ly: f64, x: f64, y: f64) -> f64 {
    if nx < 2 || ny < 2 || w.len() != nx * ny {
        return 0.0;
    }
    let gx = (x / lx.max(1.0e-30)) * (nx - 1) as f64;
    let gy = (y / ly.max(1.0e-30)) * (ny - 1) as f64;
    let i0 = gx.floor().clamp(0.0, (nx - 2) as f64) as usize;
    let j0 = gy.floor().clamp(0.0, (ny - 2) as f64) as usize;
    let tx = (gx - i0 as f64).clamp(0.0, 1.0);
    let ty = (gy - j0 as f64).clamp(0.0, 1.0);
    let at = |i: usize, j: usize| w[j * nx + i];
    (1.0 - tx) * (1.0 - ty) * at(i0, j0)
        + tx * (1.0 - ty) * at(i0 + 1, j0)
        + (1.0 - tx) * ty * at(i0, j0 + 1)
        + tx * ty * at(i0 + 1, j0 + 1)
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

/// Optional flow-driven Helmholtz volume facing the plate monopoles.
struct PlateCavity {
    sys: fs_phs::PortHamiltonian,
    x: Vec<f64>,
}

/// Linear compact radiators plus an optional von Karman pHS.
#[derive(Default)]
pub struct PlateBank {
    /// Linear modal monopoles.
    pub linear: Vec<CompactBody>,
    /// Von Karman plate, if requested.
    pub vk: Option<VkBody>,
    /// Lumped Helmholtz volume, if the assembly declared one.
    cavity: Option<PlateCavity>,
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
        let u_vol = self.volume_velocity();
        let p_cav = if let Some(cav) = self.cavity.as_mut() {
            let rec = fs_phs::step(&cav.sys, &cav.x, &[u_vol], dt)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
            cav.x = rec.x;
            cav.sys.output(&cav.x)[0]
        } else {
            0.0
        };
        let mut p = 0.0;
        for body in &mut self.linear {
            let f_cav = p_cav * body.area_m2;
            p += body.drive_and_radiate(force_n + f_cav, dt, rho, listener_m);
        }
        if let Some(vk) = &mut self.vk {
            let a_vk: f64 = vk.areas.iter().sum();
            p += vk
                .drive_and_radiate(force_n + p_cav * a_vk, dt, rho, listener_m)?
                .0;
        }
        Ok(p)
    }

    /// Face the bank with a flow-driven Helmholtz volume whose damper
    /// is the compact-mouth radiation resistance at `ω₀`.
    ///
    /// # Errors
    /// Non-physical cavity or pHS admission.
    pub fn attach_cavity(
        &mut self,
        cavity: fs_scenario::HelmholtzCavity,
        gas: &GasState,
    ) -> Result<(), AcousticRealizeError> {
        if !(cavity.volume_m3 > 0.0 && cavity.neck_radius_m > 0.0 && cavity.neck_length_m >= 0.0) {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "Helmholtz cavity geometry must be physical",
            });
        }
        let pi = core::f64::consts::PI;
        let area = pi * cavity.neck_radius_m * cavity.neck_radius_m;
        let l_eff = cavity.neck_length_m + 2.0 * (8.0 / (3.0 * pi)) * cavity.neck_radius_m;
        let omega0 = gas.sound_speed * (area / (cavity.volume_m3 * l_eff)).sqrt();
        let r_rad = fs_phs::compact_radiation_impedance(
            gas.density,
            gas.sound_speed,
            cavity.neck_radius_m,
            omega0,
            fs_phs::MouthFlange::Unflanged,
        )
        .map_or(0.0, |(r, _)| r);
        let sys = fs_phs::helmholtz_resonator_flow(
            cavity.volume_m3,
            cavity.neck_radius_m,
            cavity.neck_length_m,
            gas.density,
            gas.sound_speed,
            r_rad,
        )
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        self.cavity = Some(PlateCavity {
            x: vec![0.0; sys.state_dim()],
            sys,
        });
        Ok(())
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
#[allow(clippy::too_many_lines)] // one coherent certification stage
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
            pretension: plate.pretension_n_m,
            support: if plate.clamped {
                EdgeSupport::Clamped
            } else {
                EdgeSupport::SimplySupported
            },
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
    if let Ok(te) = thermoelastic_for_density(plate.density_kg_m3, 293.15) {
        for body in &mut out {
            body.zeta += fs_material::visco::loss_factor_to_zeta(
                te.loss_factor(body.omega, plate.thickness_m),
            );
        }
    }
    if out.is_empty() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "no usable plate radiators after harvesting",
        });
    }
    Ok(out)
}

fn thermoelastic_for_density(
    density: f64,
    t0: f64,
) -> Result<fs_material::visco::ThermoelasticZener, fs_material::visco::ViscoError> {
    if density > 5_000.0 {
        fs_material::visco::ThermoelasticZener::structural_steel(t0)
    } else {
        fs_material::visco::ThermoelasticZener::aluminum(t0)
    }
}

/// Rayleigh-integral baffled-piston face impedance `Z = p/v` under
/// `e^{-iωt}` (mass-like `Im Z < 0`). This is the same half-space
/// kernel as `fs_bem::helmholtz::baffled_piston_impedance`, written
/// here so couple does not depend on bem (cycle through feec).
fn baffled_piston_impedance(radius: f64, omega: f64, gas: &GasState, rings: usize) -> Option<C64> {
    if !(radius > 0.0 && omega > 0.0 && rings > 0) {
        return None;
    }
    let k = omega / gas.sound_speed;
    if !(k > 0.0 && k.is_finite()) {
        return None;
    }
    let mut cells: Vec<([f64; 2], f64)> = Vec::new();
    for m in 0..rings {
        let r0 = radius * m as f64 / rings as f64;
        let r1 = radius * (m + 1) as f64 / rings as f64;
        let rc = f64::midpoint(r0, r1);
        let sectors = 6 * (m + 1);
        let band_area = core::f64::consts::PI * (r1 * r1 - r0 * r0);
        for sct in 0..sectors {
            let th = core::f64::consts::TAU * (sct as f64 + 0.5) / sectors as f64;
            cells.push((
                [rc * det::cos(th), rc * det::sin(th)],
                band_area / sectors as f64,
            ));
        }
    }
    let omega_rho = omega * gas.density;
    let total_area = core::f64::consts::PI * radius * radius;
    let mut mean_p = C64::new(0.0, 0.0);
    for (i, &(xi, ai)) in cells.iter().enumerate() {
        let mut integral = C64::new(0.0, 0.0);
        for (j, &(yj, aj)) in cells.iter().enumerate() {
            if i == j {
                let ac = (ai / core::f64::consts::PI).sqrt();
                let ka_c = k * ac;
                let self_term = C64::new(
                    det::sin(ka_c) / k.max(1.0e-18),
                    (1.0 - det::cos(ka_c)) / k.max(1.0e-18),
                );
                integral = integral + self_term.scale(2.0 * core::f64::consts::PI);
            } else {
                let dx = xi[0] - yj[0];
                let dy = xi[1] - yj[1];
                let r = det::sqrt(dx * dx + dy * dy);
                let kr = k * r;
                integral = integral + C64::new(det::cos(kr), det::sin(kr)).scale(aj / r);
            }
        }
        let p_i = integral * C64::new(0.0, -omega_rho / (2.0 * core::f64::consts::PI));
        mean_p = mean_p + p_i.scale(ai / total_area);
    }
    Some(mean_p)
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
    // Rayleigh-integral face impedance (same kernel as fs-bem, no
    // bem dep). Small-ka series is the fallback if a ring solve
    // refuses. Acoustics `e^{-iωt}` → vfit `e^{+iωt}` via conjugate.
    let mut omega = Vec::with_capacity(n);
    let mut z = Vec::with_capacity(n);
    let zc = gas.density * gas.sound_speed;
    for k in 0..n {
        let t = k as f64 / (n as f64 - 1.0);
        let w = omega_lo * det::exp(t * det::ln(omega_hi / omega_lo));
        let z_face = baffled_piston_impedance(radius, w, gas, 8).unwrap_or_else(|| {
            let ka = w * radius / gas.sound_speed;
            C64::new(
                zc * 0.5 * ka * ka,
                -zc * 8.0 * ka / (3.0 * core::f64::consts::PI),
            )
        });
        omega.push(w);
        z.push(C64::new(z_face.re, -z_face.im));
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
            pretension_n_m: 0.0,
            clamped: false,
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

    #[test]
    fn clamped_and_pretension_raise_the_certified_frequency() {
        let ss = ThinPlate {
            length_m: 0.20,
            width_m: 0.15,
            thickness_m: 0.002,
            density_kg_m3: 7800.0,
            e1_pa: 200e9,
            e2_pa: 200e9,
            nu12: 0.3,
            g12_pa: 200e9 / (2.0 * 1.3),
            damping_ratio: 0.01,
            n_modes: 1,
            geometric_nonlinearity: false,
            pretension_n_m: 0.0,
            clamped: false,
        };
        let mut clamp = ss;
        clamp.clamped = true;
        let mut taut = ss;
        taut.pretension_n_m = 5.0e4;
        let w_ss = certified_radiators(ss).expect("ss")[0].omega;
        let w_cl = certified_radiators(clamp).expect("clamped")[0].omega;
        let w_t = certified_radiators(taut).expect("taut")[0].omega;
        assert!(w_cl > w_ss, "clamped ω {w_cl} must exceed SS {w_ss}");
        assert!(w_t > w_ss, "pretension ω {w_t} must exceed SS {w_ss}");
    }

    #[test]
    fn rayleigh_piston_is_passive_and_mass_like() {
        let gas = GasState::try_new(
            &fs_material::gas::GasSpec::dry_air_ussa1976(),
            293.15,
            101_325.0,
        )
        .expect("air");
        let z = baffled_piston_impedance(0.05, 2.0e3, &gas, 8).expect("z");
        assert!(z.re > 0.0, "resistance must be positive");
        assert!(
            z.im < 0.0,
            "mass-like reactance is negative under e^{{-iωt}}"
        );
    }

    #[test]
    fn clamped_and_orthotropic_von_karman_use_fe_modes() {
        let mut plate = ThinPlate {
            length_m: 0.20,
            width_m: 0.15,
            thickness_m: 0.002,
            density_kg_m3: 7800.0,
            e1_pa: 200e9,
            e2_pa: 200e9,
            nu12: 0.3,
            g12_pa: 200e9 / (2.0 * 1.3),
            damping_ratio: 0.01,
            n_modes: 1,
            geometric_nonlinearity: true,
            pretension_n_m: 0.0,
            clamped: true,
        };
        let mut clamped = VkBody::from_plate(plate).expect("clamped VK");
        let (p_cl, _) = clamped
            .drive_and_radiate(1.0, 1.0e-5, 1.2, 1.0)
            .expect("drive clamped");
        assert!(p_cl.is_finite() && p_cl.abs() > 0.0);

        plate.clamped = false;
        plate.e2_pa = 0.6 * plate.e1_pa;
        let mut ortho = VkBody::from_plate(plate).expect("orthotropic VK");
        let (p_or, _) = ortho
            .drive_and_radiate(1.0, 1.0e-5, 1.2, 1.0)
            .expect("drive ortho");
        assert!(p_or.is_finite() && p_or.abs() > 0.0);
        assert!(
            (p_cl - p_or).abs() > 1.0e-16,
            "clamped and orthotropic FE banks must not be identical"
        );
    }
}
