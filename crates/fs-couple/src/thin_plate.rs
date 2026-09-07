//! Certified compact radiators from a thin plate.
//!
//! Geometry + orthotropic section go in; `fs-plate` assembles the DKT
//! pencil and `fs-modal` returns inertia-certified eigenpairs. A
//! guitar top, a bulkhead, and a panel are the same object. Radiation
//! reaction is the baffled-piston small-`ka` series fitted by
//! `fs-vfit`, not a named instrument radiator.

use crate::acoustic_realize::AcousticRealizeError;
use crate::string_specimen::KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::EvaluationDecision;
use fs_material::gas::GasState;
use fs_material::state_point::{
    DENSITY_PROPERTY, IsotropicThermoelasticStatePoint, ORTHOTROPIC_POISSON_RATIO_PROPERTIES,
    ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES, ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES,
    POISSON_RATIO_PROPERTY, ResolvedMaterialStatePoint, YOUNG_MODULUS_PROPERTY,
};
use fs_material::visco::{RayleighDamping, ThermoelasticZener};
use fs_math::c64::C64;
use fs_math::det;
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateChart, PlateError, PlateMesh, PlateModel, PlateRegion,
    PlateSection, assemble, modes,
};
use fs_qty::{Density, Dims, DynViscosity, Pressure, QuantitySpec};
use fs_scenario::{
    IsotropicPlateBendingViscosity, IsotropicPlateThermal, RadiatingPlate, ThinPlate,
};
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
        plane_stress_section(plate)?;
        PlateBendingLaw::new(&plate)?;
        let isotropic = has_isotropic_elasticity(&plate);
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
        let mut zetas = vec![plate.damping_ratio; disp.len()];
        let loss = PlateBendingLaw::new(&plate)?;
        if loss.is_active() {
            let d = plane_stress_section(plate)?.d[0];
            for ((z, md), &omega) in zetas.iter_mut().zip(&disp).zip(&model.storage.omegas) {
                let k2 = (md.m as f64 * core::f64::consts::PI / plate.length_m).powi(2)
                    + (md.n as f64 * core::f64::consts::PI / plate.width_m).powi(2);
                let fraction = bending_energy_fraction(d * k2, d * k2 + plate.pretension_n_m)?;
                *z += fraction * loss.zeta(omega)?;
            }
        }
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
        let section = plane_stress_section(plate)?;
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
        let bending_loss = PlateBendingLoss::new(plate, &mesh, &section, &boundary)?;
        let nx = nx_c + 1;
        let ny = ny_c + 1;
        let mut disp = Vec::with_capacity(n_keep);
        let mut zetas = Vec::with_capacity(n_keep);
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
            zetas.push(
                plate.damping_ratio
                    + match &bending_loss {
                        Some(loss) => loss.modal_zeta(&model, &pair.phi, omega)?,
                        None => 0.0,
                    },
            );
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
    PlateBendingLaw::new(&plate)?;
    if plate.n_modes == 0 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "thin plate needs at least one mode",
        });
    }
    let section = plane_stress_section(plate)?;
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
    let projection = PlateProjection::new(&mesh, plate.length_m / nx as f64);
    let bending_loss = PlateBendingLoss::new(plate, &mesh, &section, &boundary)?;
    harvest_radiators(
        &model,
        &projection,
        (lo, hi),
        plate.n_modes,
        plate.damping_ratio,
        bending_loss.as_ref(),
    )
}

/// Mechanical and modal controls for an arbitrary, possibly heterogeneous chart.
/// Regional materials do not imply a damping law: the scalar ratio here is
/// authored separately. Thermal transport and nonlinear membranes are absent.
#[derive(Debug, Clone)]
pub struct PlateChartRadiation {
    /// Uniform prestress and boundary support, applied to the chart boundary.
    pub assembly: AssemblyOptions,
    /// Eigenvalue search interval [(rad/s)²], not prescribed modal frequencies.
    /// Only this interval is searched; it need not contain the fundamental.
    pub eigenvalue_window: (f64, f64),
    /// Maximum number of certified in-window modes to retain (at least one).
    pub n_modes: usize,
    /// Authored nonnegative viscous damping ratio at omega0 and 4 omega0.
    pub damping_ratio: f64,
    /// Nodal transverse forces per unit total drive force, in chart node order.
    /// Finite signed weights must sum to one. Forces on supported w DOFs are
    /// reacted by the support; they are not redistributed to free nodes.
    pub unit_force_weights: Vec<f64>,
}

/// Derive radiators from the chart's actual element sections and geometry.
/// Uses the same modal mass, signed area, force projection and piston model as
/// [`certified_radiators`]. This is a perfectly bonded, flat linear plate with
/// a compact baffled monopole observer, not a full acoustic field solution.
///
/// # Errors
/// Invalid chart, force footprint, modal controls or eigenproblem refuse.
pub fn certified_chart_radiators(
    chart: &PlateChart,
    options: &PlateChartRadiation,
) -> Result<Vec<CompactBody>, AcousticRealizeError> {
    let model = chart.assemble(&[], &options.assembly).map_err(map_plate)?;
    let projection = PlateProjection::with_force_weights(&chart.mesh, &options.unit_force_weights)?;
    harvest_radiators(
        &model,
        &projection,
        options.eigenvalue_window,
        options.n_modes,
        options.damping_ratio,
        None,
    )
}

fn harvest_radiators(
    model: &PlateModel,
    projection: &PlateProjection,
    window: (f64, f64),
    n_modes: usize,
    damping_ratio: f64,
    bending_loss: Option<&PlateBendingLoss>,
) -> Result<Vec<CompactBody>, AcousticRealizeError> {
    if n_modes == 0
        || !damping_ratio.is_finite()
        || damping_ratio < 0.0
        || !window.0.is_finite()
        || !window.1.is_finite()
        || window.0 < 0.0
        || window.1 <= window.0
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "plate needs positive modal budget, nonnegative finite damping and ordered nonnegative finite window",
        });
    }
    let report = modes(model, window, &fs_modal::SliceOptions::default()).map_err(map_plate)?;
    if report.modes.is_empty() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "plate modal window returned no certified modes",
        });
    }
    let n_keep = n_modes.min(report.modes.len());
    let mut out = Vec::with_capacity(n_keep);
    let mut bending_zetas = Vec::with_capacity(n_keep);
    for pair in report.modes.iter().take(n_keep) {
        let omega = pair.lambda.max(0.0).sqrt();
        if !(omega > 0.0) {
            continue;
        }
        out.push(projection.mode(model, &pair.phi, omega, damping_ratio)?);
        bending_zetas.push(match bending_loss {
            Some(loss) => loss.modal_zeta(model, &pair.phi, omega)?,
            None => 0.0,
        });
    }
    if out.len() >= 2 {
        let w0 = out[0].omega;
        // Same authored ratio at ω0 and 4ω0: the Rayleigh bowl then
        // splits the higher certified modes. Pinning at the first and
        // last kept frequencies would assign them identical zetas.
        if let Ok(rayleigh) =
            RayleighDamping::from_two_points(w0, damping_ratio, 4.0 * w0, damping_ratio)
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
    for (body, bending_zeta) in out.iter_mut().zip(bending_zetas) {
        body.zeta += bending_zeta;
    }
    Ok(out)
}

// Weak-loss modal projection: only bending stores the strain energy associated
// with the material bending laws. Static geometric stiffness changes
// the mode and its total energy, but must not acquire that loss angle.
// See Fedorov et al., https://arxiv.org/abs/1807.07086, Eq. (1), Sec. II.
struct PlateBendingLaw {
    thermal: Option<ThermoelasticZener>,
    viscous: Option<IsotropicPlateBendingViscosity>,
    viscous_time_s: f64,
    thickness: f64,
}

impl PlateBendingLaw {
    fn new(plate: &ThinPlate) -> Result<Self, AcousticRealizeError> {
        let thermal = admitted_thermoelastic(plate)?;
        let mut viscous_time_s = 0.0;
        if let Some(viscous) = plate.kelvin_voigt_bending {
            let (lo, hi) = viscous.omega_band_rad_s;
            viscous_time_s = viscous.viscosity_pa_s / plate.e1_pa;
            if !has_isotropic_elasticity(plate) || plate.damping_ratio != 0.0 {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "Kelvin-Voigt plate loss requires isotropic elasticity and zero authored damping",
                });
            }
            if !viscous.viscosity_pa_s.is_finite()
                || viscous.viscosity_pa_s < 0.0
                || !viscous_time_s.is_finite()
                || (viscous.viscosity_pa_s > 0.0 && viscous_time_s == 0.0)
                || !lo.is_finite()
                || !hi.is_finite()
                || lo < 0.0
                || hi < lo
            {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "plate viscosity and angular-frequency applicability must be finite and nonnegative with an ordered band",
                });
            }
        }
        Ok(Self {
            thermal,
            viscous: plate.kelvin_voigt_bending,
            viscous_time_s,
            thickness: plate.thickness_m,
        })
    }

    fn is_active(&self) -> bool {
        self.thermal.is_some() || self.viscous.is_some()
    }

    fn zeta(&self, omega: f64) -> Result<f64, AcousticRealizeError> {
        let mut zeta = match self.thermal {
            Some(law) => thermoelastic_zeta(law, omega, self.thickness)?,
            None => 0.0,
        };
        if let Some(viscous) = self.viscous {
            let (lo, hi) = viscous.omega_band_rad_s;
            if !omega.is_finite() || omega < lo || omega > hi {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "retained plate frequency lies outside material bending-loss applicability",
                });
            }
            // Moment = D (curvature + (eta/E) curvature_rate), hence
            // C_b = (eta/E) K_b. The caller projects the bending/total energy
            // ratio, so this does not damp conservative geometric stiffness.
            zeta += 0.5 * self.viscous_time_s * omega;
        }
        if !zeta.is_finite() {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "plate bending loss is unrepresentable",
            });
        }
        Ok(zeta)
    }
}

struct PlateBendingLoss {
    law: PlateBendingLaw,
    bending: Option<PlateModel>,
}

impl PlateBendingLoss {
    fn new(
        plate: ThinPlate,
        mesh: &PlateMesh,
        section: &PlateSection,
        boundary: &[usize],
    ) -> Result<Option<Self>, AcousticRealizeError> {
        let law = PlateBendingLaw::new(&plate)?;
        if !law.is_active() {
            return Ok(None);
        }
        let bending = if plate.pretension_n_m == 0.0 {
            None
        } else {
            Some(
                assemble(
                    mesh,
                    section,
                    boundary,
                    &[],
                    &AssemblyOptions {
                        pretension: 0.0,
                        support: if plate.clamped {
                            EdgeSupport::Clamped
                        } else {
                            EdgeSupport::SimplySupported
                        },
                    },
                )
                .map_err(map_plate)?,
            )
        };
        Ok(Some(Self { law, bending }))
    }

    fn modal_zeta(
        &self,
        total: &PlateModel,
        phi: &[f64],
        omega: f64,
    ) -> Result<f64, AcousticRealizeError> {
        let fraction = if let Some(bending) = &self.bending {
            let mut k_phi = vec![0.0; phi.len()];
            bending.k.spmv(phi, &mut k_phi);
            let bending_energy = phi.iter().zip(&k_phi).map(|(p, k)| p * k).sum();
            total.k.spmv(phi, &mut k_phi);
            let total_energy = phi.iter().zip(&k_phi).map(|(p, k)| p * k).sum();
            bending_energy_fraction(bending_energy, total_energy)?
        } else {
            1.0
        };
        Ok(fraction * self.law.zeta(omega)?)
    }
}

fn bending_energy_fraction(bending: f64, total: f64) -> Result<f64, AcousticRealizeError> {
    let fraction = bending / total;
    if !bending.is_finite()
        || bending < 0.0
        || !total.is_finite()
        || total <= 0.0
        || !fraction.is_finite()
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "bending-loss participation requires finite nonnegative bending and positive total modal energy",
        });
    }
    // Stable compressive prestress may produce a ratio above one; do not clamp it.
    Ok(fraction)
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
    fn with_force_weights(mesh: &PlateMesh, weights: &[f64]) -> Result<Self, AcousticRealizeError> {
        let total: f64 = weights.iter().sum();
        if weights.len() != mesh.node_count()
            || weights.iter().any(|w| !w.is_finite())
            || !total.is_finite()
            || (total - 1.0).abs() > 1.0e-12
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "plate drive needs one finite force weight per node summing to one",
            });
        }
        let mut nodal_area = vec![0.0; mesh.node_count()];
        for &[i, j, k] in &mesh.tris {
            let (a, b, c) = (mesh.nodes[i], mesh.nodes[j], mesh.nodes[k]);
            let area = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0));
            for node in [i, j, k] {
                nodal_area[node] += area / 3.0;
            }
        }
        Ok(Self {
            area: nodal_area.iter().sum(),
            nodal_area,
            unit_force: weights.to_vec(),
        })
    }

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

fn plane_stress_section(plate: ThinPlate) -> Result<PlateSection, AcousticRealizeError> {
    PlateSection::orthotropic_plane_stress_at_angle(
        plate.e1_pa,
        plate.e2_pa,
        plate.nu12,
        plate.g12_pa,
        plate.thickness_m,
        plate.density_kg_m3,
        plate.material_angle_rad,
    )
    .map_err(map_plate)
}

/// Elastic approximation and mapping from card axes to section material axes.
/// `ThinPlate::material_angle_rad` then rotates that section in the rectangle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlateMaterialModel {
    /// Isotropic density, Young's modulus and Poisson ratio; derive G.
    Isotropic,
    /// Card axes 1,2 map to section axes 1,2; axis 3 is normal.
    Orthotropic12,
    /// Exchange card axes 1,2 in the section; apply Poisson reciprocity.
    Orthotropic21,
}

/// Independent thickness prescription for material comparison at fixed x,y size.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PlateThicknessConstraint {
    /// Current thickness [m] is fixed; mass changes with density.
    FixedThickness(f64),
    /// Assigned mass [kg] is fixed; derive thickness as `m/(rho area)`.
    FixedMass(f64),
}

/// A named, perfectly bonded plate region at one resolved material state.
/// Its mass prescription applies to this region's actual triangles, not a
/// bounding rectangle. Material state and authored mechanical choices remain
/// separate; no thermal strain or history is inferred from a state-point query.
#[derive(Clone, Debug, PartialEq)]
pub struct PlateRegionMaterial {
    /// Unique nonempty label and triangle indices, in the supplied mesh.
    pub region: PlateRegion,
    /// Exact resolved properties and their original source receipts.
    pub material: ResolvedMaterialStatePoint,
    /// Constitutive approximation and material-axis mapping.
    pub model: PlateMaterialModel,
    /// Fixed local thickness or fixed total mass of this region.
    pub thickness_constraint: PlateThicknessConstraint,
    /// Counterclockwise material-axis angle in the chart plane [rad].
    pub material_angle_rad: f64,
}

/// One resolved region, exposed read-only by [`ResolvedPlateChart`].
#[derive(Clone, Debug)]
pub struct ResolvedPlateRegion {
    /// Retained assignment, including source data and the original mass policy.
    pub input: PlateRegionMaterial,
    /// Sum of the assigned mid-surface triangle areas [m²].
    pub area_m2: f64,
    /// Mass derived from this region's section and area [kg].
    pub mass_kg: f64,
    /// Rotated section consumed by every triangle in this region.
    pub section: PlateSection,
}

/// Compiled plate geometry and regional sources that cannot drift independently.
/// All access is immutable; material substitution recompiles the same geometry
/// and retained constraints. A cloned numerical chart is an authored copy and
/// does not retain this binding's source authority automatically.
#[derive(Clone, Debug)]
pub struct ResolvedPlateChart {
    chart: PlateChart,
    regions: Vec<ResolvedPlateRegion>,
    region_for_element: Vec<usize>,
    mass_kg: f64,
    identity: ContentHash,
}

impl ResolvedPlateChart {
    /// Numerical chart accepted by the existing plate and acoustic operators.
    #[must_use]
    pub const fn chart(&self) -> &PlateChart {
        &self.chart
    }

    /// Original assignments and their resolved section/area/mass.
    #[must_use]
    pub fn regions(&self) -> &[ResolvedPlateRegion] {
        &self.regions
    }

    /// Exact material source for a triangle, or `None` for an invalid index.
    #[must_use]
    pub fn material_at_element(&self, element: usize) -> Option<&ResolvedMaterialStatePoint> {
        self.region_for_element
            .get(element)
            .map(|&i| &self.regions[i].input.material)
    }

    /// Total mass integrated in mesh element order [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }

    /// Geometry, boundary-node set and resolved per-element material identity.
    /// Region labels/order and repartitioning with identical resolved
    /// per-element inputs do not affect it.
    /// It excludes solver controls and unresolved state evolution; it is not
    /// an identity for a complete acoustic scenario or the authored mass policy.
    #[must_use]
    pub const fn specimen_identity(&self) -> ContentHash {
        self.identity
    }

    /// Replace one named region's material and rebuild geometry-derived inputs.
    /// Other sources, triangle assignments, orientations and mass/thickness
    /// prescriptions are retained. A failed rebind leaves `self` unchanged.
    ///
    /// # Errors
    /// Unknown names or inadmissible replacement data refuse.
    pub fn with_region_material(
        &self,
        name: &str,
        material: &ResolvedMaterialStatePoint,
    ) -> Result<Self, AcousticRealizeError> {
        let i = self
            .regions
            .iter()
            .position(|r| r.input.region.name == name)
            .ok_or_else(|| plate_assignment_error(Some(name), None, "unknown material region"))?;
        let mut inputs: Vec<_> = self.regions.iter().map(|r| r.input.clone()).collect();
        inputs[i].material = material.clone();
        compile_plate_material_chart(
            self.chart.mesh.clone(),
            self.chart.boundary_nodes.clone(),
            inputs,
        )
    }
}

/// Compile complete, disjoint regional material assignments into the existing
/// plate chart. Shared nodes bond regions perfectly. All local sections use the
/// same plane-stress reduction as uniform specimens, with h = m/(rho A_region)
/// for a fixed regional mass. No homogenization, interfacial law, thermal strain,
/// nonlinear membrane, or continuous time evolution is synthesized.
///
/// # Errors
/// Invalid geometry, duplicate/empty region labels, overlaps, gaps, invalid
/// triangle indices or physical material/section inputs refuse. All triangles
/// must be assigned exactly once; there is no fallback material.
pub fn compile_plate_material_chart(
    mesh: PlateMesh,
    mut boundary_nodes: Vec<usize>,
    inputs: Vec<PlateRegionMaterial>,
) -> Result<ResolvedPlateChart, AcousticRealizeError> {
    let mesh = PlateMesh::from_unstructured(mesh.nodes, mesh.tris).map_err(map_plate)?;
    let mut assigned = vec![None; mesh.tris.len()];
    let mut names = std::collections::BTreeSet::new();
    for (r, input) in inputs.iter().enumerate() {
        let name = input.region.name.as_str();
        if name.trim().is_empty() || !names.insert(name) {
            return Err(plate_assignment_error(
                Some(name),
                None,
                "region labels must be nonempty and unique",
            ));
        }
        if input.region.triangle_indices.is_empty() {
            return Err(plate_assignment_error(
                Some(name),
                None,
                "material region is empty",
            ));
        }
        for &element in &input.region.triangle_indices {
            let slot = assigned.get_mut(element).ok_or_else(|| {
                plate_assignment_error(Some(name), Some(element), "triangle index outside mesh")
            })?;
            if slot.replace(r).is_some() {
                return Err(plate_assignment_error(
                    Some(name),
                    Some(element),
                    "triangle assigned more than once",
                ));
            }
        }
    }
    let region_for_element: Vec<usize> = assigned
        .into_iter()
        .enumerate()
        .map(|(e, r)| {
            r.ok_or_else(|| {
                plate_assignment_error(None, Some(e), "triangle has no material assignment")
            })
        })
        .collect::<Result<_, _>>()?;
    let mut areas = vec![0.0; inputs.len()];
    let mut element_areas = Vec::with_capacity(mesh.tris.len());
    // Fixed mesh order makes all regional integrals independent of the authored
    // region order and of the order of triangle indices inside each region.
    for (e, t) in mesh.tris.iter().enumerate() {
        let (a, b, c) = (mesh.nodes[t[0]], mesh.nodes[t[1]], mesh.nodes[t[2]]);
        let area = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0));
        areas[region_for_element[e]] += area;
        element_areas.push(area);
    }
    let mut regions = Vec::with_capacity(inputs.len());
    for (input, area_m2) in inputs.into_iter().zip(areas) {
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(plate_assignment_error(
                Some(&input.region.name),
                None,
                "region area is not finite and positive",
            ));
        }
        let rho = plate_property(&input.material, DENSITY_PROPERTY, Density::DIMS)?;
        let (e1, e2, nu, g) = plate_elasticity(&input.material, input.model)?;
        let h = match input.thickness_constraint {
            PlateThicknessConstraint::FixedThickness(h) => h,
            PlateThicknessConstraint::FixedMass(m) => m / area_m2 / rho,
        };
        let section = PlateSection::orthotropic_plane_stress_at_angle(
            e1,
            e2,
            nu,
            g,
            h,
            rho,
            input.material_angle_rad,
        )
        .map_err(map_plate)?;
        let mass_kg = rho * h * area_m2;
        if !mass_kg.is_finite() || mass_kg <= 0.0 {
            return Err(plate_assignment_error(
                Some(&input.region.name),
                None,
                "region mass overflows or underflows",
            ));
        }
        regions.push(ResolvedPlateRegion {
            input,
            area_m2,
            mass_kg,
            section,
        });
    }
    let sections: Vec<_> = region_for_element
        .iter()
        .map(|&r| regions[r].section)
        .collect();
    let mass_kg: f64 = sections
        .iter()
        .zip(element_areas)
        .map(|(s, a)| s.density * s.thickness * a)
        .sum();
    if !mass_kg.is_finite() || mass_kg <= 0.0 {
        return Err(plate_assignment_error(
            None,
            None,
            "total plate mass overflows or underflows",
        ));
    }
    boundary_nodes.sort_unstable();
    boundary_nodes.dedup();
    let chart = PlateChart::with_boundary_and_regions(
        mesh,
        sections[0],
        boundary_nodes,
        regions.iter().map(|r| r.input.region.clone()).collect(),
    )
    .map_err(map_plate)?
    .with_element_sections(sections)
    .map_err(map_plate)?;
    let mut identity = DomainHasher::new("org.frankensim.fs-couple.regional-plate-specimen.v1");
    identity.update(&(chart.mesh.nodes.len() as u64).to_le_bytes());
    for &(x, y) in &chart.mesh.nodes {
        identity.update(&x.to_bits().to_le_bytes());
        identity.update(&y.to_bits().to_le_bytes());
    }
    identity.update(&(chart.mesh.tris.len() as u64).to_le_bytes());
    for (t, &r) in chart.mesh.tris.iter().zip(&region_for_element) {
        for &node in t {
            identity.update(&(node as u64).to_le_bytes());
        }
        let region = &regions[r];
        identity.update(region.input.material.identity().as_bytes());
        identity.update(&[match region.input.model {
            PlateMaterialModel::Isotropic => 0,
            PlateMaterialModel::Orthotropic12 => 1,
            PlateMaterialModel::Orthotropic21 => 2,
        }]);
        identity.update(&region.section.thickness.to_bits().to_le_bytes());
        identity.update(&region.input.material_angle_rad.to_bits().to_le_bytes());
    }
    identity.update(&(chart.boundary_nodes.len() as u64).to_le_bytes());
    for &node in &chart.boundary_nodes {
        identity.update(&(node as u64).to_le_bytes());
    }
    Ok(ResolvedPlateChart {
        chart,
        regions,
        region_for_element,
        mass_kg,
        identity: identity.finalize(),
    })
}

fn plate_assignment_error(
    region: Option<&str>,
    element: Option<usize>,
    what: &'static str,
) -> AcousticRealizeError {
    AcousticRealizeError::PlateAssignment {
        region: region.map(str::to_owned),
        element,
        what,
    }
}

/// Uniform material-bound plate with the original property-use receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedPlateSpecimen {
    plate: ThinPlate,
    material: ResolvedMaterialStatePoint,
    model: PlateMaterialModel,
    thickness_constraint: PlateThicknessConstraint,
    mass_kg: f64,
    identity: ContentHash,
}

impl ResolvedPlateSpecimen {
    /// Select proportional isotropic Kelvin–Voigt bending at the resolved state.
    /// Requires `kelvin_voigt_bending_viscosity` [Pa s] and constant rho/E/nu
    /// claims. The viscosity claim must declare an omega band [rad/s], intersected
    /// with every coefficient's band. Other state coordinates remain frozen.
    /// The viscous tensor shares the elastic Poisson ratio; separate bulk/shear
    /// relaxation, nonlinear membrane viscosity and thermal evolution are absent.
    /// Original property receipts remain available through [`Self::material`].
    ///
    /// # Errors
    /// Missing/nonconstant coefficients, wrong dimensions, anisotropy, authored
    /// damping, or invalid viscosity/applicability refuse. Actual equilibrium
    /// modal frequencies are checked by realization, including zero viscosity.
    pub fn with_kelvin_voigt_bending_loss(mut self) -> Result<Self, AcousticRealizeError> {
        let refuse = |what| AcousticRealizeError::InvalidDescription { what };
        if self.model != PlateMaterialModel::Isotropic {
            return Err(refuse(
                "material plate bending viscosity requires isotropic elasticity",
            ));
        }
        let property = self
            .material
            .property(KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY)
            .ok_or_else(|| refuse("material plate loss needs kelvin_voigt_bending_viscosity"))?;
        let viscosity_pa_s = plate_property(
            &self.material,
            KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
            DynViscosity::DIMS,
        )?;
        let mut band = property
            .answer()
            .evidence
            .model
            .validity
            .bound("omega")
            .ok_or_else(|| {
                refuse("material plate viscosity needs an explicit omega band in rad/s")
            })?;
        for key in [
            DENSITY_PROPERTY,
            YOUNG_MODULUS_PROPERTY,
            POISSON_RATIO_PROPERTY,
            KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
        ] {
            let answer = self
                .material
                .property(key)
                .ok_or_else(|| refuse("material plate loss needs density, E, nu and viscosity"))?
                .answer();
            if answer.receipt.decision != EvaluationDecision::ConstantWithinValidity {
                return Err(refuse(
                    "Kelvin-Voigt plate coefficients must be validity-wide scalar constants",
                ));
            }
            if let Some((lo, hi)) = answer.evidence.model.validity.bound("omega") {
                band = (band.0.max(lo), band.1.min(hi));
            }
        }
        self.plate.kelvin_voigt_bending = Some(IsotropicPlateBendingViscosity {
            viscosity_pa_s,
            omega_band_rad_s: band,
            material_state_identity: Some(self.material.identity()),
        });
        PlateBendingLaw::new(&self.plate)?;
        Ok(self)
    }

    /// Rotated constitutive section for an element or region of a plate chart.
    /// This projection carries numeric operator inputs; retain this specimen
    /// alongside the chart to retain its material receipts. Its full-rectangle
    /// mass is not the mass of a region: assembly integrates the local rho*h.
    ///
    /// # Errors
    /// Forwards section admission errors (the resolved specimen already passed it).
    pub fn section(&self) -> Result<PlateSection, AcousticRealizeError> {
        plane_stress_section(self.plate)
    }

    /// Description consumed by the existing acoustic assembly and plate operators.
    #[must_use]
    pub const fn plate(&self) -> ThinPlate {
        self.plate
    }

    /// Exact material resolution, without promoting source uncertainty or authority.
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// Authored elastic approximation and uniform material-axis mapping.
    #[must_use]
    pub const fn model(&self) -> PlateMaterialModel {
        self.model
    }

    /// Authored thickness/mass prescription, separate from material authority.
    #[must_use]
    pub const fn thickness_constraint(&self) -> PlateThicknessConstraint {
        self.thickness_constraint
    }

    /// Uniform specimen mass [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }

    /// Material, resolved geometry, axis mapping and angle identity. Excludes
    /// supports, loading, damping and solver controls; not a complete scenario hash.
    #[must_use]
    pub const fn specimen_identity(&self) -> ContentHash {
        self.identity
    }
}

/// Bind uniform plate geometry to the shared requirement-driven material API.
///
/// Isotropic binding needs `density`, `young_modulus`, `poisson_ratio`.
/// Orthotropic binding needs `density`, `young_modulus_{1,2}`,
/// `poisson_ratio_12`, `shear_modulus_12`, in the material's principal frame.
/// No out-of-plane or unrelated optical/thermal properties are required.
/// Values describe one frozen material state; no history or mass is transferred
/// between specimens. Geometry, supports, pretension and authored damping remain
/// separate from material claims. The template's uniform material angle is
/// retained. Only linear orthotropic bending is supported: the existing
/// nonlinear membrane approximation is isotropic.
///
/// # Errors
/// Missing/mismatched quantities, invalid geometry/model or unrepresentable
/// stiffness/mass refuse. A template with thermoelastic loss refuses: use a
/// freshly resolved thermoelastic bundle instead of retaining stale coefficients.
pub fn with_uniform_plate_material_state(
    mut plate: ThinPlate,
    state: &ResolvedMaterialStatePoint,
    model: PlateMaterialModel,
    thickness_constraint: PlateThicknessConstraint,
) -> Result<ResolvedPlateSpecimen, AcousticRealizeError> {
    let refuse = |what| AcousticRealizeError::InvalidDescription { what };
    // A selected viscosity law survives a material swap, its numeric inputs do
    // not. Re-resolve below from this exact state's property receipts.
    let bind_viscosity = plate.kelvin_voigt_bending.take().is_some();
    if plate.thermoelastic.is_some() {
        return Err(refuse(
            "elastic plate rebind requires a fresh thermal-loss binding",
        ));
    }
    if plate.geometric_nonlinearity && model != PlateMaterialModel::Isotropic {
        return Err(refuse(
            "material-bound nonlinear orthotropic membrane response is unavailable",
        ));
    }
    if [plate.length_m, plate.width_m]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || plate.n_modes == 0
        || [plate.damping_ratio, plate.pretension_n_m]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(refuse(
            "plate geometry, modal budget or mechanical inputs are invalid",
        ));
    }
    plate.density_kg_m3 = plate_property(state, DENSITY_PROPERTY, Density::DIMS)?;
    bind_plate_elasticity(&mut plate, state, model)?;
    plate.thickness_m = match thickness_constraint {
        PlateThicknessConstraint::FixedThickness(h) => h,
        PlateThicknessConstraint::FixedMass(m) => {
            m / plate.length_m / plate.width_m / plate.density_kg_m3
        }
    };
    plane_stress_section(plate)?;
    let mass_kg = plate.density_kg_m3 * plate.thickness_m * plate.length_m * plate.width_m;
    if !mass_kg.is_finite() || mass_kg <= 0.0 {
        return Err(refuse("plate mass overflows or underflows"));
    }
    let mut identity = DomainHasher::new("org.frankensim.fs-couple.uniform-plate-specimen.v2");
    identity.update(state.identity().as_bytes());
    identity.update(&[match model {
        PlateMaterialModel::Isotropic => 0,
        PlateMaterialModel::Orthotropic12 => 1,
        PlateMaterialModel::Orthotropic21 => 2,
    }]);
    for value in [
        plate.length_m,
        plate.width_m,
        plate.thickness_m,
        plate.material_angle_rad,
    ] {
        identity.update(&value.to_bits().to_le_bytes());
    }
    let specimen = ResolvedPlateSpecimen {
        plate,
        material: state.clone(),
        model,
        thickness_constraint,
        mass_kg,
        identity: identity.finalize(),
    };
    if bind_viscosity {
        specimen.with_kelvin_voigt_bending_loss()
    } else {
        Ok(specimen)
    }
}

/// Bind uniform isotropic plate mechanics and thermal loss to one resolved state.
///
/// Recomputes thickness/mass, elasticity and all thermal inputs. Existing thermal
/// coefficients are replaced; their presence is a law selection, not reusable
/// material data. If bending viscosity is selected, the same resolved bundle must
/// include its claim, and it too is rebound. The result retains the shared state
/// identity and original receipts for both mechanisms. Temperature is the explicit
/// specimen query coordinate, not inferred from the surrounding gas.
///
/// # Errors
/// Invalid geometry, missing viscosity data or inadmissible loss coefficients
/// refuse without modifying the input. This is a frozen isotropic specimen,
/// not a thermal history, anisotropic loss law or phase-change model.
pub fn with_uniform_isotropic_thermoelastic_material_state(
    mut plate: ThinPlate,
    state: &IsotropicThermoelasticStatePoint,
    thickness_constraint: PlateThicknessConstraint,
) -> Result<ResolvedPlateSpecimen, AcousticRealizeError> {
    let bind_viscosity = plate.kelvin_voigt_bending.take().is_some();
    plate.thermoelastic = None;
    let mut specimen = with_uniform_plate_material_state(
        plate,
        state.resolved(),
        PlateMaterialModel::Isotropic,
        thickness_constraint,
    )?;
    specimen.plate = with_isotropic_thermoelastic_state(specimen.plate, state)?;
    if bind_viscosity {
        specimen.with_kelvin_voigt_bending_loss()
    } else {
        Ok(specimen)
    }
}

fn bind_plate_elasticity(
    plate: &mut ThinPlate,
    state: &ResolvedMaterialStatePoint,
    model: PlateMaterialModel,
) -> Result<(), AcousticRealizeError> {
    let (e1, e2, nu12, g12) = plate_elasticity(state, model)?;
    plate.e1_pa = e1;
    plate.e2_pa = e2;
    plate.nu12 = nu12;
    plate.g12_pa = g12;
    Ok(())
}

fn plate_elasticity(
    state: &ResolvedMaterialStatePoint,
    model: PlateMaterialModel,
) -> Result<(f64, f64, f64, f64), AcousticRealizeError> {
    let (e1, e2, nu12, g12) = match model {
        PlateMaterialModel::Isotropic => {
            let e = plate_property(state, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?;
            let nu = plate_property(state, POISSON_RATIO_PROPERTY, Dims::NONE)?;
            if !(nu > -1.0 && nu < 0.5) {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "isotropic plate Poisson ratio must be in (-1, 0.5)",
                });
            }
            (e, e, nu, e / (2.0 * (1.0 + nu)))
        }
        PlateMaterialModel::Orthotropic12 | PlateMaterialModel::Orthotropic21 => {
            let e1 = plate_property(
                state,
                ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[0],
                Pressure::DIMS,
            )?;
            let e2 = plate_property(
                state,
                ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[1],
                Pressure::DIMS,
            )?;
            let nu = plate_property(state, ORTHOTROPIC_POISSON_RATIO_PROPERTIES[0], Dims::NONE)?;
            let g = plate_property(
                state,
                ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[0],
                Pressure::DIMS,
            )?;
            if model == PlateMaterialModel::Orthotropic21 {
                (e2, e1, nu * e2 / e1, g)
            } else {
                (e1, e2, nu, g)
            }
        }
    };
    Ok((e1, e2, nu12, g12))
}

fn plate_property(
    state: &ResolvedMaterialStatePoint,
    key: &str,
    dims: Dims,
) -> Result<f64, AcousticRealizeError> {
    let p = state
        .property(key)
        .ok_or(AcousticRealizeError::InvalidDescription {
            what: "plate material is missing a required in-plane property",
        })?;
    if p.requirement().quantity() != QuantitySpec::dimensional(dims) || !p.value_si().is_finite() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "plate properties require finite dimension-only SI scalars; semantic aliases cannot be erased",
        });
    }
    Ok(p.value_si())
}

/// Bind one resolved isotropic material state to an isotropic plate description.
/// All elastic and thermal inputs move together; geometry and excitation remain
/// caller-owned. An orthotropic description refuses instead of losing anisotropy.
pub fn with_isotropic_thermoelastic_state(
    mut plate: ThinPlate,
    state: &IsotropicThermoelasticStatePoint,
) -> Result<ThinPlate, AcousticRealizeError> {
    require_isotropic_thermoelastic(&plate)?;
    if plate.kelvin_voigt_bending.is_some() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "thermal plate rebind requires a fresh material bending-loss binding",
        });
    }
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
    if !has_isotropic_elasticity(plate) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "isotropic thermoelastic loss requires isotropic E, nu and G; anisotropic loss is unavailable",
        });
    }
    Ok(())
}

fn has_isotropic_elasticity(plate: &ThinPlate) -> bool {
    let g_iso = plate.e1_pa / (2.0 * (1.0 + plate.nu12));
    plate.e1_pa > 0.0
        && plate.e1_pa.is_finite()
        && plate.nu12 > -1.0
        && plate.nu12 < 0.5
        && (plate.e2_pa / plate.e1_pa - 1.0).abs() <= 1.0e-6
        && (plate.g12_pa / g_iso - 1.0).abs() <= 1.0e-6
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
    // Sine (1,1) Rayleigh quotient sets a search scale, not a certified
    // eigenvalue. D16/D26 cross integrals vanish for this trial function,
    // but they remain in the actual FE pencil. It is exact only for an
    // aligned, simply supported continuum without prestress/rotary inertia.
    let kx = core::f64::consts::PI / a;
    let ky = core::f64::consts::PI / b;
    let d = &section.d;
    ((d[0] * kx.powi(4) + 2.0 * (d[1] + 2.0 * d[8]) * kx.powi(2) * ky.powi(2) + d[4] * ky.powi(4))
        / (section.density * section.thickness))
        .sqrt()
}

fn map_plate(err: PlateError) -> AcousticRealizeError {
    AcousticRealizeError::Plate(err)
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
            material_angle_rad: 0.0,
            damping_ratio: 0.0,
            thermoelastic: Some(IsotropicPlateThermal {
                temperature_k: 300.0,
                linear_expansion_per_k: 20e-6,
                specific_heat_j_kg_k: 600.0,
                conductivity_w_m_k: 100.0,
                state_identity: None, // authored numerical fixture, not measured material data
            }),
            kelvin_voigt_bending: None,
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
    fn g1_thermoelastic_pretension_dilution_and_ringdown_match_sine_solution() {
        let mut plate = thermal_plate();
        plate.pretension_n_m = 50_000.0;
        plate.thermoelastic.as_mut().unwrap().linear_expansion_per_k = 2e-4;
        let k2 = core::f64::consts::PI.powi(2)
            * (plate.length_m.recip().powi(2) + plate.width_m.recip().powi(2));
        let d = plate.e1_pa * plate.thickness_m.powi(3) / (12.0 * (1.0 - plate.nu12.powi(2)));
        let fraction = d * k2 / (d * k2 + plate.pretension_n_m);
        let omega = ((d * k2.powi(2) + plate.pretension_n_m * k2)
            / (plate.density_kg_m3 * plate.thickness_m))
            .sqrt();
        let zeta = fraction * independent_zener_zeta(plate, omega);
        assert!(fraction < 0.5, "fixture must distinguish undiluted loss");
        let body = VkBody::from_plate(plate).unwrap();
        let (_, r, _) = body.sys.structure();
        assert!((r[3] / (2.0 * omega * zeta) - 1.0).abs() < 1e-12);

        // Infinitesimal free vibration removes quartic membrane effects. Check
        // the actual stepped state against the independent damped oscillator,
        // including its phase; a post-hoc amplitude envelope cannot pass this.
        let q0 = 1e-10;
        let duration = 0.01;
        let wd = omega * (1.0 - zeta * zeta).sqrt();
        let envelope = (-zeta * omega * duration).exp();
        let q_exact =
            q0 * envelope * ((wd * duration).cos() + zeta * omega / wd * (wd * duration).sin());
        let v_exact = -q0 * envelope * omega * omega / wd * (wd * duration).sin();
        let mut errors = Vec::new();
        for steps in [1000, 2000] {
            let dt = duration / steps as f64;
            let mut x = vec![q0, 0.0];
            for _ in 0..steps {
                x = fs_phs::step(&body.sys, &x, &[0.0], dt).unwrap().x;
            }
            errors
                .push(((x[0] - q_exact).powi(2) + ((x[1] - v_exact) / omega).powi(2)).sqrt() / q0);
        }
        assert!(
            errors[1] < errors[0] / 3.5,
            "order-two ringdown: {errors:?}"
        );
        assert!(
            errors[1] < 5e-4,
            "absolute phase/state accuracy: {errors:?}"
        );
    }

    #[test]
    fn g1_kelvin_voigt_sine_ringdown_matches_constitutive_equation() {
        for pretension in [0.0, 50_000.0] {
            let mut plate = thermal_plate();
            plate.thermoelastic = None;
            plate.pretension_n_m = pretension;
            plate.kelvin_voigt_bending = Some(IsotropicPlateBendingViscosity {
                viscosity_pa_s: 2e6,
                omega_band_rad_s: (0.0, 1e6),
                material_state_identity: None,
            });
            let k2 = core::f64::consts::PI.powi(2)
                * (plate.length_m.recip().powi(2) + plate.width_m.recip().powi(2));
            let section_factor = plate.thickness_m.powi(3) / (12.0 * (1.0 - plate.nu12.powi(2)));
            let mass_per_area = plate.density_kg_m3 * plate.thickness_m;
            let omega2 =
                (plate.e1_pa * section_factor * k2.powi(2) + pretension * k2) / mass_per_area;
            // Direct plane-stress Kelvin-Voigt equation: rho*h*w_tt +
            // eta*h^3/[12(1-nu^2)] Laplacian^2(w_t) + D Laplacian^2(w) - T Laplacian(w) = 0.
            let c = 2e6 * section_factor * k2.powi(2) / mass_per_area;
            let body = VkBody::from_plate(plate).unwrap();
            assert!((body.sys.structure().1[3] / c - 1.0).abs() < 1e-12);
            let q0 = 1e-10;
            let duration = 0.01;
            let wd = (omega2 - 0.25 * c * c).sqrt();
            let envelope = (-0.5 * c * duration).exp();
            let q_exact =
                q0 * envelope * ((wd * duration).cos() + c / (2.0 * wd) * (wd * duration).sin());
            let v_exact = -q0 * envelope * omega2 / wd * (wd * duration).sin();
            let mut errors = Vec::new();
            for steps in [1000, 2000] {
                let mut x = vec![q0, 0.0];
                let initial_energy = body.sys.hamiltonian(&x);
                for _ in 0..steps {
                    x = fs_phs::step(&body.sys, &x, &[0.0], duration / steps as f64)
                        .unwrap()
                        .x;
                }
                assert!(body.sys.hamiltonian(&x) < initial_energy);
                errors.push(
                    ((x[0] - q_exact).powi(2) + (x[1] - v_exact).powi(2) / omega2).sqrt() / q0,
                );
            }
            assert!(
                errors[1] < errors[0] / 3.5,
                "order-two state/phase convergence: {errors:?}"
            );
            assert!(errors[1] < 5e-4, "analytic ringdown: {errors:?}");
        }
    }

    #[test]
    fn g0_kelvin_voigt_plate_checks_all_retained_frequencies_and_inputs() {
        let mut plate = thermal_plate();
        plate.thermoelastic = None;
        plate.kelvin_voigt_bending = Some(IsotropicPlateBendingViscosity {
            viscosity_pa_s: 0.0, // zero strength still retains physical applicability
            omega_band_rad_s: (0.0, 1e6),
            material_state_identity: None,
        });
        for defect in 0..6 {
            let mut invalid = plate;
            match defect {
                0 => {
                    invalid
                        .kelvin_voigt_bending
                        .as_mut()
                        .unwrap()
                        .viscosity_pa_s = -1.0
                }
                1 => {
                    invalid
                        .kelvin_voigt_bending
                        .as_mut()
                        .unwrap()
                        .viscosity_pa_s = f64::NAN
                }
                2 => {
                    invalid
                        .kelvin_voigt_bending
                        .as_mut()
                        .unwrap()
                        .omega_band_rad_s = (10.0, 1.0)
                }
                3 => invalid.e2_pa *= 0.5,
                4 => invalid.g12_pa *= 0.5,
                _ => invalid.damping_ratio = 0.002,
            }
            assert!(certified_radiators(invalid).is_err());
            assert!(VkBody::from_plate(invalid).is_err());
        }
        for clamped in [false, true] {
            plate.clamped = clamped;
            let mut narrow = plate;
            let linear = certified_radiators(plate).unwrap();
            narrow
                .kelvin_voigt_bending
                .as_mut()
                .unwrap()
                .omega_band_rad_s
                .1 = linear[0].omega * 1.01;
            assert!(certified_radiators(narrow).is_ok());
            narrow.n_modes = 2;
            assert!(
                certified_radiators(narrow).is_err(),
                "higher linear mode must be checked"
            );

            let vk = VkBody::from_plate(plate).unwrap();
            let omega = (vk.sys.effort(&[1e-12, 0.0])[0] / 1e-12).sqrt();
            narrow.n_modes = 1;
            narrow
                .kelvin_voigt_bending
                .as_mut()
                .unwrap()
                .omega_band_rad_s
                .1 = omega * 1.01;
            assert!(VkBody::from_plate(narrow).is_ok());
            narrow.n_modes = 2;
            assert!(
                VkBody::from_plate(narrow).is_err(),
                "higher nonlinear reference mode must be checked"
            );
        }
    }

    #[test]
    fn g1_thermoelastic_fe_dilution_uses_actual_prestressed_mode_energy() {
        for (clamped, sampled_vk) in [(false, false), (true, false), (true, true)] {
            let mut plate = thermal_plate();
            plate.clamped = clamped;
            plate.pretension_n_m = 50_000.0;
            plate.damping_ratio = 0.003; // must not be diluted with thermal loss
            let (nx, ny) = if sampled_vk { (8, 8) } else { (5, 4) };
            let mesh = PlateMesh::rectangle(plate.length_m, plate.width_m, nx, ny);
            let section = plane_stress_section(plate).unwrap();
            let model = assemble(
                &mesh,
                &section,
                &PlateMesh::rectangle_boundary(nx, ny),
                &[],
                &AssemblyOptions {
                    pretension: plate.pretension_n_m,
                    support: if clamped {
                        EdgeSupport::Clamped
                    } else {
                        EdgeSupport::SimplySupported
                    },
                },
            )
            .unwrap();
            let w11 = ss_omega11(&section, plate.length_m, plate.width_m);
            let report = modes(
                &model,
                ((0.25 * w11).powi(2), (12.0 * w11).powi(2)),
                &fs_modal::SliceOptions::default(),
            )
            .unwrap();
            let pair = &report.modes[0];
            let omega = pair.lambda.sqrt();
            // Independent geometric-energy integral of the P1 displacement
            // gradient. Subtract it from total energy; do not call the new
            // bending-matrix projection or infer a fraction from frequency alone.
            for scale in [1.0, -3.0] {
                let phi: Vec<_> = pair.phi.iter().map(|p| scale * p).collect();
                let mut kp = vec![0.0; phi.len()];
                model.k.spmv(&phi, &mut kp);
                let total: f64 = phi.iter().zip(kp).map(|(p, k)| p * k).sum();
                let mut geometric = 0.0;
                for tri in &mesh.tris {
                    let [a, b, c] = tri.map(|i| mesh.nodes[i]);
                    let twice_area = (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0);
                    let [wa, wb, wc] = tri.map(|i| model.dof_map[3 * i].map_or(0.0, |r| phi[r]));
                    let gx = (wa * (b.1 - c.1) + wb * (c.1 - a.1) + wc * (a.1 - b.1)) / twice_area;
                    let gy = (wa * (c.0 - b.0) + wb * (a.0 - c.0) + wc * (b.0 - a.0)) / twice_area;
                    geometric += plate.pretension_n_m * 0.5 * twice_area * (gx * gx + gy * gy);
                }
                assert!(geometric / total > 0.1);
                for (thermal, viscous) in [(true, false), (false, true), (true, true)] {
                    let mut input = plate;
                    let mut undiluted = if thermal {
                        independent_zener_zeta(plate, omega)
                    } else {
                        0.0
                    };
                    if !thermal {
                        input.thermoelastic = None;
                    }
                    if viscous {
                        input.damping_ratio = 0.0;
                        input.kelvin_voigt_bending = Some(IsotropicPlateBendingViscosity {
                            viscosity_pa_s: 2e6,
                            omega_band_rad_s: (0.0, 1e6),
                            material_state_identity: None,
                        });
                        undiluted += 0.5 * 2e6 / plate.e1_pa * omega;
                    }
                    let expected = (1.0 - geometric / total) * undiluted;
                    let actual = if sampled_vk {
                        let body = VkBody::from_plate(input).unwrap();
                        let (_, r, _) = body.sys.structure();
                        r[3] / (2.0 * omega) - input.damping_ratio
                    } else {
                        certified_radiators(input).unwrap()[0].zeta - input.damping_ratio
                    };
                    assert!(
                        (actual / expected - 1.0).abs() < 1e-9,
                        "clamped={clamped}, sampled={sampled_vk}, scale={scale}, thermal={thermal}, viscous={viscous}: {actual} vs {expected}"
                    );
                }
            }
        }
    }

    #[test]
    fn g1_shear_anisotropy_reaches_nonlinear_bending_and_bad_angles_refuse() {
        let mut plate = thermal_plate();
        plate.thermoelastic = None;
        let isotropic = VkBody::from_plate(plate).unwrap();
        let q = 1e-12;
        let iso_omega2 = isotropic.sys.effort(&[q, 0.0])[0] / q;
        plate.g12_pa *= 0.5;
        let anisotropic = VkBody::from_plate(plate).unwrap();
        let aniso_omega2 = anisotropic.sys.effort(&[q, 0.0])[0] / q;
        assert!(
            aniso_omega2 < 0.98 * iso_omega2,
            "shear anisotropy must lower the actual bending tangent even when E1 == E2"
        );
        for angle in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            plate.material_angle_rad = angle;
            assert!(certified_radiators(plate).is_err());
            assert!(VkBody::from_plate(plate).is_err());
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
            material_angle_rad: 0.0,
            damping_ratio: 0.01,
            thermoelastic: None,
            kelvin_voigt_bending: None,
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
            material_angle_rad: 0.0,
            damping_ratio: 0.01,
            thermoelastic: None,
            kelvin_voigt_bending: None,
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
            material_angle_rad: 0.0,
            damping_ratio: 0.01,
            thermoelastic: None,
            kelvin_voigt_bending: None,
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
