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
use fs_material::state_point::IsotropicThermoelasticStatePoint;
use fs_material::visco::{RayleighDamping, ThermoelasticZener};
use fs_math::c64::C64;
use fs_math::det;
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateError, PlateMesh, PlateSection, assemble, modes,
};
use fs_scenario::{IsotropicPlateThermal, RadiatingPlate, ThinPlate};
use fs_vfit::FitOptions;
use fs_vfit::discretize::{DigitalFilter, DigitalFilterState, realize_tabulated_impedance};

/// One driven compact radiator harvested from a certified plate mode.
#[derive(Debug, Clone)]
pub struct CompactBody {
    /// Signed modal monopole area [m²] (`∫ φ dA`). Zero is silent
    /// in this compact monopole approximation, not in all radiation models.
    pub area_m2: f64,
    /// Modal mass [kg] (`φᵀ M φ`), in the same basis as both ports.
    pub mass_kg: f64,
    /// Signed projection of the unit-total-force drive footprint onto φ.
    /// One for a caller-authored lumped radiator.
    pub drive_participation: f64,
    /// Angular frequency [rad/s].
    pub omega: f64,
    /// Viscous damping ratio.
    pub zeta: f64,
    y: f64,
    v: f64,
    // Physical face area, independent of modal normalization and cancellation.
    piston_area_m2: f64,
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
            drive_participation: 1.0,
            omega: core::f64::consts::TAU * spec.frequency_hz,
            zeta: spec.damping_ratio,
            y: 0.0,
            v: 0.0,
            piston_area_m2: spec.area_m2,
            rad: None,
        })
    }

    /// Advance under a generalized force and return acceleration.
    ///
    /// # Errors
    /// Propagates the radiation-filter step refusal instead of silently
    /// demoting the reaction force to zero (which would mask model
    /// instability as benign decay).
    pub fn drive(&mut self, force_n: f64, dt: f64) -> Result<f64, AcousticRealizeError> {
        let f_rad = if let Some((filter, state)) = self.rad.as_mut() {
            // p = Z_face * mean surface velocity, Q = area_modal * q_dot.
            // Generalized reaction is -area_modal * p, preserving work and
            // eigenvector scaling even for negative or nearly cancelling modes.
            match filter.step(state, self.v * self.area_m2 / self.piston_area_m2) {
                Ok(p_face) => -p_face * self.area_m2,
                Err(e) => return Err(AcousticRealizeError::Nonlinear(e.to_string())),
            }
        } else {
            0.0
        };
        let acc = (force_n + f_rad) / self.mass_kg
            - 2.0 * self.zeta * self.omega * self.v
            - self.omega * self.omega * self.y;
        self.v += dt * acc;
        self.y += dt * self.v;
        Ok(acc)
    }

    /// Compact monopole pressure at distance `listener_m`.
    #[must_use]
    pub fn radiate(&self, acc: f64, rho: f64, listener_m: f64) -> f64 {
        // Baffled half-space (same piston as the self-load), not free-space.
        rho * self.area_m2 * acc / (2.0 * core::f64::consts::PI * listener_m)
    }

    /// Observe an acceleration in the pHS mass-normalized coordinate.
    pub(crate) fn radiate_mass_normalized(&self, acc: f64, rho: f64, listener_m: f64) -> f64 {
        self.radiate(acc / self.mass_kg.sqrt(), rho, listener_m)
    }

    /// Volume velocity of the monopole [m³/s].
    #[must_use]
    pub fn volume_velocity(&self) -> f64 {
        self.area_m2 * self.v
    }

    /// Drive and radiate in one step.
    ///
    /// # Errors
    /// Propagates [`CompactBody::drive`] radiation refusals.
    pub fn drive_and_radiate(
        &mut self,
        force_n: f64,
        dt: f64,
        rho: f64,
        listener_m: f64,
    ) -> Result<f64, AcousticRealizeError> {
        let acc = self.drive(force_n, dt)?;
        Ok(self.radiate(acc, rho, listener_m))
    }

    fn attach_piston_load(&mut self, gas: &GasState, sample_rate_hz: u32) {
        if self.rad.is_some() || self.area_m2 == 0.0 {
            return;
        }
        let radius = (self.piston_area_m2 / core::f64::consts::PI).sqrt();
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
        admitted_thermoelastic(&plate)?;
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
        let zetas = vk_zetas(plate, &model.storage.omegas)?;
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
        let zetas = vk_zetas(plate, &vk.storage.omegas)?;
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

fn vk_zetas(plate: ThinPlate, omegas: &[f64]) -> Result<Vec<f64>, AcousticRealizeError> {
    let mut zetas = vec![plate.damping_ratio; omegas.len()];
    if let Some(te) = admitted_thermoelastic(&plate)? {
        for (z, &w) in zetas.iter_mut().zip(omegas) {
            *z += thermoelastic_zeta(te, w, plate.thickness_m)?;
        }
    }
    Ok(zetas)
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
            p += body.drive_and_radiate(
                force_n * body.drive_participation + f_cav,
                dt,
                rho,
                listener_m,
            )?;
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
    let thermoelastic = admitted_thermoelastic(&plate)?;
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
    let projection = PlateProjection::new(&mesh, plate.length_m / nx as f64);
    let mut out = Vec::with_capacity(n_keep);
    for pair in report.modes.iter().take(n_keep) {
        let omega = pair.lambda.max(0.0).sqrt();
        if !(omega > 0.0) {
            continue;
        }
        out.push(projection.mode(&model, &pair.phi, omega, plate.damping_ratio)?);
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
    if let Some(te) = thermoelastic {
        for body in &mut out {
            body.zeta += thermoelastic_zeta(te, body.omega, plate.thickness_m)?;
        }
    }
    if out.is_empty() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "no usable plate radiators after harvesting",
        });
    }
    Ok(out)
}

// P1 surface quadrature on the DKT nodal displacement trace. The production
// rectangle's first cell strip is an explicit uniform traction footprint;
// its integral is one unit of total bridge force. This does not reconstruct
// a higher-order DKT displacement field or a finite-wavelength radiator.
struct PlateProjection {
    nodal_area: Vec<f64>,
    unit_force: Vec<f64>,
    area: f64,
}

impl PlateProjection {
    fn new(mesh: &PlateMesh, drive_strip_width: f64) -> Self {
        let mut nodal_area = vec![0.0; mesh.node_count()];
        let mut unit_force = vec![0.0; mesh.node_count()];
        let mut drive_area = 0.0;
        for &[i, j, k] in &mesh.tris {
            let (a, b, c) = (mesh.nodes[i], mesh.nodes[j], mesh.nodes[k]);
            let area = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0));
            // Strip boundary follows element edges on the generated rectangle.
            let driven = (a.0 + b.0 + c.0) / 3.0 < drive_strip_width;
            for node in [i, j, k] {
                nodal_area[node] += area / 3.0;
                if driven {
                    unit_force[node] += area / 3.0;
                }
            }
            if driven {
                drive_area += area;
            }
        }
        for weight in &mut unit_force {
            *weight /= drive_area;
        }
        let area = nodal_area.iter().sum();
        Self {
            nodal_area,
            unit_force,
            area,
        }
    }

    fn mode(
        &self,
        model: &fs_plate::PlateModel,
        phi: &[f64],
        omega: f64,
        zeta: f64,
    ) -> Result<CompactBody, AcousticRealizeError> {
        let refuse = || AcousticRealizeError::InvalidDescription {
            what: "plate mode needs finite ports and a positive finite mass in the same basis",
        };
        if phi.len() != model.free || phi.iter().any(|v| !v.is_finite()) {
            return Err(refuse());
        }
        let mut m_phi = vec![0.0; phi.len()];
        model.m.spmv(phi, &mut m_phi);
        let mass: f64 = phi.iter().zip(m_phi).map(|(p, m)| p * m).sum();
        let (mut area, mut drive) = (0.0, 0.0);
        for node in 0..self.nodal_area.len() {
            if let Some(r) = model.dof_map[3 * node] {
                area += self.nodal_area[node] * phi[r];
                drive += self.unit_force[node] * phi[r];
            }
        }
        if !(mass > 0.0
            && mass.is_finite()
            && area.is_finite()
            && drive.is_finite()
            && self.area > 0.0
            && self.area.is_finite())
        {
            return Err(refuse());
        }
        Ok(CompactBody {
            area_m2: area,
            mass_kg: mass,
            drive_participation: drive,
            omega,
            zeta,
            y: 0.0,
            v: 0.0,
            piston_area_m2: self.area,
            rad: None,
        })
    }
}

/// Bind one resolved isotropic material state to an isotropic plate description.
/// All elastic and thermal inputs move together; geometry and excitation remain
/// caller-owned. An orthotropic description refuses instead of losing anisotropy.
pub fn with_isotropic_thermoelastic_state(
    mut plate: ThinPlate,
    state: &IsotropicThermoelasticStatePoint,
) -> Result<ThinPlate, AcousticRealizeError> {
    require_isotropic_thermoelastic(&plate)?;
    let law = state.law();
    plate.density_kg_m3 = law.rho;
    plate.e1_pa = law.e;
    plate.e2_pa = law.e;
    plate.nu12 = state.poisson_ratio();
    plate.g12_pa = law.e / (2.0 * (1.0 + plate.nu12));
    plate.thermoelastic = Some(IsotropicPlateThermal {
        temperature_k: law.t0,
        linear_expansion_per_k: law.alpha_t,
        specific_heat_j_kg_k: law.cp,
        conductivity_w_m_k: law.conductivity,
        state_identity: Some(state.resolved().identity()),
    });
    admitted_thermoelastic(&plate)?;
    Ok(plate)
}

fn require_isotropic_thermoelastic(plate: &ThinPlate) -> Result<(), AcousticRealizeError> {
    let g_iso = plate.e1_pa / (2.0 * (1.0 + plate.nu12));
    if !(plate.e1_pa > 0.0
        && plate.e1_pa.is_finite()
        && plate.nu12 > -1.0
        && plate.nu12 < 0.5
        && (plate.e2_pa / plate.e1_pa - 1.0).abs() <= 1.0e-6
        && (plate.g12_pa / g_iso - 1.0).abs() <= 1.0e-6)
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "isotropic thermoelastic loss requires isotropic E, nu and G; anisotropic loss is unavailable",
        });
    }
    Ok(())
}

fn admitted_thermoelastic(
    plate: &ThinPlate,
) -> Result<Option<ThermoelasticZener>, AcousticRealizeError> {
    let Some(thermal) = plate.thermoelastic else {
        return Ok(None);
    };
    require_isotropic_thermoelastic(plate)?;
    if !thermal.linear_expansion_per_k.is_finite()
        || [
            plate.density_kg_m3,
            plate.thickness_m,
            thermal.temperature_k,
            thermal.specific_heat_j_kg_k,
            thermal.conductivity_w_m_k,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "thermoelastic loss requires finite expansion and positive finite rho, thickness, T, cp and conductivity",
        });
    }
    let law = ThermoelasticZener {
        e: plate.e1_pa,
        alpha_t: thermal.linear_expansion_per_k,
        t0: thermal.temperature_k,
        rho: plate.density_kg_m3,
        cp: thermal.specific_heat_j_kg_k,
        conductivity: thermal.conductivity_w_m_k,
    };
    let tau = law.relaxation_time(plate.thickness_m);
    if !(law.relaxation_strength().is_finite() && tau > 0.0 && tau.is_finite()) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "thermoelastic relaxation strength or time is unrepresentable",
        });
    }
    Ok(Some(law))
}

fn thermoelastic_zeta(
    law: ThermoelasticZener,
    omega: f64,
    thickness: f64,
) -> Result<f64, AcousticRealizeError> {
    let eta = law.loss_factor(omega, thickness);
    if !(eta >= 0.0 && eta.is_finite()) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "thermoelastic modal loss is unrepresentable",
        });
    }
    Ok(fs_material::visco::loss_factor_to_zeta(eta))
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

    fn projection_fixture() -> (PlateMesh, fs_plate::PlateModel, PlateProjection) {
        let mesh = PlateMesh::rectangle(2.0, 1.0, 2, 1);
        let section = PlateSection::isotropic(70e9, 0.3, 0.01, 1000.0).unwrap();
        let model = assemble(
            &mesh,
            &section,
            &[],
            &[],
            &AssemblyOptions {
                pretension: 0.0,
                support: EdgeSupport::SimplySupported,
            },
        )
        .unwrap();
        let projection = PlateProjection::new(&mesh, 1.0);
        (mesh, model, projection)
    }

    #[test]
    fn g1_plate_projection_integrates_mass_area_and_unit_force_footprint() {
        let (mesh, model, projection) = projection_fixture();
        let mut phi = vec![0.0; model.free];
        for (node, &(x, _)) in mesh.nodes.iter().enumerate() {
            // Manufactured affine trace w=x, wx=1, wy=0, not an eigenmode.
            phi[model.dof_map[3 * node].unwrap()] = x;
            phi[model.dof_map[3 * node + 1].unwrap()] = 1.0;
        }
        let body = projection.mode(&model, &phi, 100.0, 0.0).unwrap();
        // Translation uses the owner's lumped nodal mass: sum(w_i² A_i)=3.
        // The constant x slope also carries the owner's rotary inertia.
        let expected_mass = 1000.0 * 0.01 * 3.0 + 1000.0 * 0.01_f64.powi(3) / 12.0 * 2.0;
        assert!((body.mass_kg / expected_mass - 1.0).abs() < 1.0e-14);
        assert!((body.area_m2 - 2.0).abs() < 1.0e-14); // integral x over [0,2]×[0,1]
        assert!((body.drive_participation - 0.5).abs() < 1.0e-14); // mean x on [0,1]×[0,1]
        assert!((projection.unit_force.iter().sum::<f64>() - 1.0).abs() < 1.0e-14);
        assert!(
            projection
                .mode(&model, &vec![0.0; model.free], 100.0, 0.0)
                .is_err()
        );
        assert!(
            projection
                .mode(&model, &phi[..phi.len() - 1], 100.0, 0.0)
                .is_err()
        );
    }

    #[test]
    fn g3_plate_projection_retains_monopole_cancellation() {
        let (mesh, model, projection) = projection_fixture();
        let mut phi = vec![0.0; model.free];
        for (node, &(x, _)) in mesh.nodes.iter().enumerate() {
            phi[model.dof_map[3 * node].unwrap()] = x - 1.0;
        }
        let mut body = projection.mode(&model, &phi, 100.0, 0.0).unwrap();
        assert!(
            body.area_m2.abs() < 1.0e-15,
            "opposite lobes cancel, without an area floor"
        );
        assert!((body.drive_participation + 0.5).abs() < 1.0e-14);
        let acc = body.drive(body.drive_participation, 1.0e-4).unwrap();
        assert!(acc.abs() > 0.01, "this mode is mechanically excited");
        assert!(body.radiate(acc, 1.2, 1.0).abs() < 1.0e-15);
    }

    #[test]
    fn g3_plate_projection_sign_and_scale_preserve_pressure_energy_and_self_load() {
        let (mesh, model, projection) = projection_fixture();
        let mut phi = vec![0.0; model.free];
        for (node, &(x, _)) in mesh.nodes.iter().enumerate() {
            phi[model.dof_map[3 * node].unwrap()] = 1.0 - 1.5 * x;
        }
        let gas = GasState::try_new(
            &fs_material::gas::GasSpec::dry_air_ussa1976(),
            293.15,
            101_325.0,
        )
        .unwrap();
        let base = projection.mode(&model, &phi, 100.0, 0.01).unwrap();
        assert!(
            base.area_m2 < 0.0 && base.drive_participation > 0.0,
            "signed transfer falsifier"
        );
        for loaded in [false, true] {
            let mut histories = Vec::new();
            for scale in [1.0, -4.0, 0.125] {
                let scaled: Vec<_> = phi.iter().map(|v| scale * v).collect();
                let mut body = projection.mode(&model, &scaled, 100.0, 0.01).unwrap();
                if loaded {
                    body.attach_piston_load(&gas, 48_000);
                    assert!(body.rad.is_some());
                }
                let mut history = Vec::new();
                for j in 0..256 {
                    let force = (j as f64 * 0.03).cos();
                    let acc = body
                        .drive(force * body.drive_participation, 1.0 / 48_000.0)
                        .unwrap();
                    let energy =
                        0.5 * body.mass_kg * (body.v.powi(2) + (body.omega * body.y).powi(2));
                    history.push([
                        body.radiate(acc, gas.density, 1.0),
                        body.volume_velocity(),
                        energy,
                    ]);
                }
                histories.push(history);
            }
            for other in &histories[1..] {
                for (a, b) in histories[0].iter().zip(other) {
                    for k in 0..3 {
                        assert!(
                            (a[k] - b[k]).abs() <= 1.0e-10 * a[k].abs().max(1.0e-20),
                            "loaded={loaded}, quantity={k}: {a:?} vs {b:?}"
                        );
                    }
                }
            }
        }
    }

    fn thermal_plate() -> ThinPlate {
        ThinPlate {
            length_m: 0.20,
            width_m: 0.15,
            thickness_m: 0.002,
            density_kg_m3: 4999.0,
            e1_pa: 70e9,
            e2_pa: 70e9,
            nu12: 0.3,
            g12_pa: 70e9 / 2.6,
            damping_ratio: 0.0,
            thermoelastic: Some(IsotropicPlateThermal {
                temperature_k: 300.0,
                linear_expansion_per_k: 20e-6,
                specific_heat_j_kg_k: 600.0,
                conductivity_w_m_k: 100.0,
                state_identity: None, // authored numerical fixture, not measured material data
            }),
            n_modes: 1,
            geometric_nonlinearity: false,
            pretension_n_m: 0.0,
            clamped: false,
        }
    }

    fn independent_zener_zeta(plate: ThinPlate, omega: f64) -> f64 {
        let t = plate.thermoelastic.expect("explicit thermal inputs");
        let heat_per_volume = plate.density_kg_m3 * t.specific_heat_j_kg_k;
        let tau = plate.thickness_m.powi(2) * heat_per_volume
            / (core::f64::consts::PI.powi(2) * t.conductivity_w_m_k);
        let delta =
            plate.e1_pa * t.linear_expansion_per_k.powi(2) * t.temperature_k / heat_per_volume;
        0.5 * delta * omega * tau / (1.0 + (omega * tau).powi(2))
    }

    #[test]
    fn g1_thermoelastic_loss_reaches_linear_and_nonlinear_plate_operators() {
        let base = thermal_plate();
        let mut expansion_twin = base;
        expansion_twin
            .thermoelastic
            .as_mut()
            .unwrap()
            .linear_expansion_per_k *= 2.0;
        let mut density_twin = base;
        density_twin.density_kg_m3 = 5001.0;
        let mut hot = base;
        hot.thermoelastic.as_mut().unwrap().temperature_k = 450.0;
        let mut losses = Vec::new();
        for plate in [base, expansion_twin, density_twin, hot] {
            let body = certified_radiators(plate).expect("linear plate").remove(0);
            let expected = independent_zener_zeta(plate, body.omega);
            assert!((body.zeta / expected - 1.0).abs() < 1e-12);
            let vk = VkBody::from_plate(plate).expect("nonlinear plate");
            // Read the assembled pHS, not the helper that manufactured its zeta.
            // The infinitesimal tangent removes the von Karman quartic term.
            let q = 1e-12;
            let omega = (vk.sys.effort(&[q, 0.0])[0] / q).sqrt();
            let (_, r, _) = vk.sys.structure();
            let nonlinear_zeta = r[3] / (2.0 * omega);
            let expected = independent_zener_zeta(plate, omega);
            assert!((nonlinear_zeta / expected - 1.0).abs() < 1e-10);
            let x = [0.0, 1e-4];
            let step = fs_phs::step(&vk.sys, &x, &[0.0], 1e-5).expect("damped step");
            assert!(vk.sys.hamiltonian(&step.x) < vk.sys.hamiltonian(&x));
            losses.push((body.zeta, nonlinear_zeta));
        }
        for column in [0, 1] {
            let loss = |i: usize| {
                if column == 0 {
                    losses[i].0
                } else {
                    losses[i].1
                }
            };
            assert!(
                (loss(1) / loss(0) - 4.0).abs() < 1e-10,
                "equal density must not erase different expansion coefficients"
            );
            assert!(
                (loss(2) / loss(0) - 1.0).abs() < 0.005,
                "crossing 5000 kg/m3 must not switch material laws"
            );
            assert!(
                (loss(3) / loss(0) - 1.5).abs() < 1e-10,
                "specimen temperature must reach both damping operators"
            );
        }
    }

    #[test]
    fn g0_thermoelastic_loss_is_explicit_and_refuses_anisotropy_or_bad_inputs() {
        let mut missing = thermal_plate();
        missing.thermoelastic = None;
        assert_eq!(certified_radiators(missing).unwrap()[0].zeta, 0.0);
        let vk = VkBody::from_plate(missing).expect("unclaimed thermal loss");
        assert!(vk.sys.structure().1.iter().all(|value| *value == 0.0));
        for defect in 0..5 {
            let mut invalid = thermal_plate();
            match defect {
                0 => invalid.e2_pa *= 0.8,
                1 => invalid.g12_pa *= 0.8, // E1 == E2 does not establish isotropy
                2 => invalid.thermoelastic.as_mut().unwrap().conductivity_w_m_k = 0.0,
                3 => invalid.thermoelastic.as_mut().unwrap().temperature_k = f64::NAN,
                _ => invalid.thermoelastic.as_mut().unwrap().specific_heat_j_kg_k = -1.0,
            }
            assert!(matches!(
                certified_radiators(invalid),
                Err(AcousticRealizeError::InvalidDescription { .. })
            ));
            assert!(matches!(
                VkBody::from_plate(invalid),
                Err(AcousticRealizeError::InvalidDescription { .. })
            ));
        }
    }

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
            thermoelastic: None,
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
            thermoelastic: None,
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
            thermoelastic: None,
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
