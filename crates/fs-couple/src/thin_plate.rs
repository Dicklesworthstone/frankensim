//! Certified compact radiators from a thin plate.
//!
//! Geometry + orthotropic section go in; `fs-plate` assembles the DKT
//! pencil and `fs-modal` returns inertia-certified eigenpairs. A
//! guitar top, a bulkhead, and a panel are the same object. Radiation
//! reaction is the shared exact baffled-piston impedance fitted by
//! `fs-vfit`, not a named instrument radiator.

use crate::acoustic_realize::AcousticRealizeError;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, advance_exact_zoh};
use crate::plate_viscosity::bending_energy_fraction;
use crate::string_specimen::KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY;
use fs_blake3::{ContentHash, DomainHasher};
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

/// Full DKT coordinates from the SAME solve as a material image's ports.
/// Internal consumers use rotations for constitutive projection; a transverse
/// nodal trace alone cannot reconstruct curvature or material memory.
#[derive(Clone, Debug)]
pub(crate) struct RetainedPlateBasis {
    pub(crate) source_identity: ContentHash,
    pub(crate) dof_map: Vec<Option<usize>>,
    pub(crate) mass_normalized_modes: Vec<Vec<f64>>,
}

/// Original DKT in-plane bending strains at one area quadrature point.
/// Homogeneous thickness is integrated analytically: strain is evaluated at
/// z=h/sqrt(12), with the full volume weight. Surface strain is sqrt(3) times
/// this strain. This identity is restricted to linear pure bending.
pub(crate) struct RetainedPlateBendingPoint {
    pub(crate) region: usize,
    pub(crate) columns: Vec<[f64; 3]>,
    pub(crate) volume_m3: f64,
}

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
    // Producer-owned identity: modes harvested from one physical face share
    // one uniform-piston volume-flow load. Equal areas never imply identity.
    piston_face: Option<std::sync::Arc<()>>,
    rad: Option<PistonLoad>,
}

/// The numerical filter and the physical state/clock from which it was built.
/// Keeping these together prevents a retained reaction from silently using a
/// different atmosphere or sample rate than the pressure observer.
#[derive(Debug, Clone)]
struct PistonLoad {
    filter: DigitalFilter,
    state: DigitalFilterState,
    gas: GasState,
    sample_rate_hz: u32,
}

/// One exact-ZOH step of a [`CompactBody`]: its start state and held force.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct HeldStep {
    y: f64,
    v: f64,
    force_n: f64,
}

impl HeldStep {
    /// A body at rest under no force (the quiescent past before a render).
    pub(crate) const REST: Self = Self { y: 0.0, v: 0.0, force_n: 0.0 };
}

impl CompactBody {
    // Shared geometric reducers already use M-orthonormal physical vectors.
    // Unit modal mass preserves their signed volume and force projections;
    // the face area remains independent of modal cancellation/normalization.
    pub(crate) fn from_mass_normalized_projection(
        omega: f64, area: f64, drive: f64, face_area: f64,
    ) -> Result<Self, AcousticRealizeError> {
        if !omega.is_finite() || omega <= 0.0 || !area.is_finite() || !drive.is_finite()
            || !face_area.is_finite() || face_area < 0.0 || (face_area == 0.0 && area != 0.0) {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "geometric radiator requires finite physical mode projections",
            });
        }
        Ok(Self { area_m2:area, mass_kg:1.0, drive_participation:drive,
            omega, zeta:0.0, y:0.0, v:0.0, piston_area_m2:face_area, piston_face:None, rad:None })
    }
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
            piston_face: None,
            rad: None,
        })
    }

    /// Advance under a held generalized force and return endpoint acceleration.
    /// Uses the shared exact-ZOH oscillator transition. When present, the
    /// radiation filter supplies a reaction held over this step; the coupled
    /// filter/structure evolution is still a sampled approximation.
    ///
    /// # Errors
    /// Refuses invalid force, duration, mass, frequency or damping and
    /// non-finite motion. A filter already advanced before a motion refusal
    /// is not rolled back; callers must discard a refused coupled run.
    /// A retained digital radiation load requires its original sample period;
    /// a clock mismatch refuses before either the filter or body advances.
    /// Propagates the radiation-filter step refusal instead of silently
    /// demoting the reaction force to zero (which would mask model
    /// instability as benign decay).
    pub fn drive(&mut self, force_n: f64, dt: f64) -> Result<f64, AcousticRealizeError> {
        self.drive_held(force_n, dt).map(|(acc, _)| acc)
    }

    /// [`CompactBody::drive`], also returning the step's start state and the
    /// total generalized force (reaction included) held over it.
    pub(crate) fn drive_held(
        &mut self,
        force_n: f64,
        dt: f64,
    ) -> Result<(f64, HeldStep), AcousticRealizeError> {
        if !(force_n.is_finite()
            && dt.is_finite()
            && dt > 0.0
            && self.mass_kg.is_finite()
            && self.mass_kg > 0.0
            && self.omega.is_finite()
            && self.omega > 0.0
            && self.zeta.is_finite()
            && self.zeta >= 0.0)
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "compact radiator step requires finite force, positive duration/mass/frequency and nonnegative damping",
            });
        }
        let f_rad = if let Some(load) = self.rad.as_mut() {
            if dt.to_bits() != f64::from(load.sample_rate_hz).recip().to_bits() {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "piston radiation load and mechanical step require the same sample clock",
                });
            }
            // p = Z_face * mean surface velocity, Q = area_modal * q_dot.
            // Generalized reaction is -area_modal * p, preserving work and
            // eigenvector scaling even for negative or nearly cancelling modes.
            match load.filter.step(&mut load.state, self.v * self.area_m2 / self.piston_area_m2) {
                Ok(p_face) => -p_face * self.area_m2,
                Err(e) => return Err(AcousticRealizeError::Nonlinear(e.to_string())),
            }
        } else {
            0.0
        };
        let step = HeldStep { y: self.y, v: self.v, force_n: force_n + f_rad };
        let (y, v, acc) = self.transition(step, dt)?;
        self.v = v;
        self.y = y;
        Ok((acc, step))
    }

    /// Exact acceleration `elapsed_s` into a recorded step of this body, under
    /// the same exact-ZOH law; the body's own state is not touched. At the full
    /// step duration it is the acceleration [`CompactBody::drive`] returned.
    ///
    /// # Errors
    /// Invalid elapsed time, or non-finite motion.
    pub(crate) fn acceleration_within(
        &self,
        step: HeldStep,
        elapsed_s: f64,
    ) -> Result<f64, AcousticRealizeError> {
        if !(elapsed_s.is_finite() && elapsed_s > 0.0) {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "retarded radiator time must lie inside a recorded step",
            });
        }
        self.transition(step, elapsed_s).map(|(_, _, acc)| acc)
    }

    fn transition(&self, step: HeldStep, dt: f64) -> Result<(f64, f64, f64), AcousticRealizeError> {
        let root_mass = self.mass_kg.sqrt();
        let next = advance_exact_zoh(
            ModalAcousticMode {
                angular_frequency_rad_s: self.omega,
                damping_ratio: self.zeta,
                pressure_per_modal_velocity: C64::new(0.0, 0.0),
            },
            ModalAcousticState {
                displacement_m_sqrt_kg: step.y * root_mass,
                velocity_m_sqrt_kg_per_s: step.v * root_mass,
            },
            step.force_n / root_mass,
            dt,
        );
        let y = next.displacement_m_sqrt_kg / root_mass;
        let v = next.velocity_m_sqrt_kg_per_s / root_mass;
        let acc = step.force_n / self.mass_kg
            - 2.0 * self.zeta * self.omega * v
            - self.omega * self.omega * y;
        if !(y.is_finite() && v.is_finite() && acc.is_finite()) {
            return Err(AcousticRealizeError::Nonlinear(
                "compact radiator transition is not finite".to_string(),
            ));
        }
        Ok((y, v, acc))
    }

    /// Compact monopole pressure at distance `listener_m`.
    #[must_use]
    pub fn radiate(&self, acc: f64, rho: f64, listener_m: f64) -> f64 {
        // Baffled half-space (same piston as the self-load), not free-space.
        rho * self.area_m2 * acc / (2.0 * core::f64::consts::PI * listener_m)
    }

    /// Observe an acceleration in the pHS mass-normalized coordinate.
    #[must_use]
    pub fn radiate_mass_normalized(&self, acc: f64, rho: f64, listener_m: f64) -> f64 {
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

    /// Load this body's radiating piston with the baffled-piston radiation of `gas` at
    /// `sample_rate_hz`: the per-body digital filter that the one-port bowed termination
    /// (`Termination::PlateOnePort`) attaches to its body. A body with no radiating
    /// area is left unloaded, and a body that already carries a load keeps it.
    ///
    /// # Errors
    /// An incoherent `gas`; a retained load whose gas or clock differs from these; a
    /// non-positive piston face area or a zero sample clock; the filter's own refusal.
    pub fn attach_piston_load(
        &mut self,
        gas: &GasState,
        sample_rate_hz: u32,
    ) -> Result<(), AcousticRealizeError> {
        gas.validate_coherence().map_err(|e| AcousticRealizeError::Ambient(e.to_string()))?;
        self.validate_radiation_load(gas, sample_rate_hz)?;
        if !self.piston_area_m2.is_finite() || self.piston_area_m2 <= 0.0 || sample_rate_hz == 0 {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "piston radiation loading requires a positive physical face area and sample clock",
            });
        }
        if self.rad.is_some() || self.area_m2 == 0.0 {
            return Ok(());
        }
        let radius = (self.piston_area_m2 / core::f64::consts::PI).sqrt();
        let filter = piston_radiation_filter(radius, gas, sample_rate_hz)?;
        let state = filter.zero_state();
        self.rad = Some(PistonLoad { filter, state, gas: *gas, sample_rate_hz });
        Ok(())
    }

    /// Admit a retained load against the exact environment and clock of its
    /// consumer without clearing filter memory or rebuilding material physics.
    pub(crate) fn validate_radiation_load(
        &self, gas: &GasState, sample_rate_hz: u32,
    ) -> Result<(), AcousticRealizeError> {
        if let Some(load) = &self.rad {
            if load.gas != *gas || load.sample_rate_hz != sample_rate_hz {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "retained piston radiation load differs from the declared gas or sample clock",
                });
            }
        }
        Ok(())
    }
}

mod cavity;
mod bank;
mod piston;
mod prepared;
mod radiation;
pub use bank::{PlateBank, PlateBankImage, PlateBankObservation, PlateBankStep};
pub use piston::{BaffledPistonRadiation, PistonRadiationFit};
pub use cavity::{VkCavityInitial, VkCavityLoad, VkCavityObservation};
pub use prepared::{PreparedVkBody, PreparedVkStep};

/// A von Karman plate as one pHS.
///
/// Isotropic simply-supported plates use analytic sine modes.
/// Clamped or orthotropic bending uses DKT-sampled displacement
/// with the same sine Airy membrane channel.
pub struct VkBody {
    sys: fs_phs::PortHamiltonian,
    linear: fs_phs::linear_exact::LinearSplit,
    x: Vec<f64>,
    /// Signed surface projections in mass-normalized coordinates [m²/sqrt(kg)].
    pub areas: Vec<f64>,
}

impl VkBody {
    /// Realize a source-reduced nonlinear plate with its complete viscous
    /// damping matrix, including cross-mode terms and undamped subspaces.
    ///
    /// All inputs use the SAME mass-normalized coordinates as `storage`:
    /// `modal_damping_per_s` is row-major C in `q'' + C q' + dV/dq = b F`,
    /// `force_participation` is the physical force projection b [1/sqrt(kg)],
    /// and `areas` contains the signed surface integrals [m²/sqrt(kg)].
    /// A source reducer must supply these projections; no damping ratio or
    /// modal area is fitted, inferred, or added to the supplied operator.
    ///
    /// The existing pHS owner checks passive damping and advances the same
    /// nonlinear storage. Radiation remains the compact baffled monopole
    /// observation of this body, with no radiation back-reaction inferred.
    ///
    /// # Errors
    /// Inconsistent/nonfinite projections, non-passive damping, invalid
    /// storage, or the original pHS/time-resolution admission refusal.
    pub fn from_modal_storage_with_damping(
        storage: fs_nlmodal::SosModalStorage,
        modal_damping_per_s: &[f64],
        force_participation: &[f64],
        areas: Vec<f64>,
    ) -> Result<Self, AcousticRealizeError> {
        finish_vk_damping(
            storage, VkDamping::Matrix(modal_damping_per_s),
            force_participation, areas, false,
        )
    }

    /// Add causal Maxwell material memory to a body at rest.
    ///
    /// Each branch uses the original mass-normalized state order
    /// `[q0, p0, q1, p1, ...]` and may project only displacement entries.
    /// The existing pHS owner adds `(sqrt(k) b.q - z)^2 / 2` to the
    /// equilibrium storage and evolves `z` with the declared relaxation time.
    /// Cross-mode strain, frequency-dependent stiffness and dissipation are
    /// retained together; the branch is never replaced by a damping ratio.
    /// For a source projection already scaled by sqrt(stiffness), use k = 1.
    /// The source reducer must supply branches in this body's SAME basis and
    /// within their constitutive applicability domain.
    ///
    /// This is one cold attachment: the original state must be zero and have
    /// no material memories. All new histories start at zero, which is both
    /// relaxed and unrelaxed at zero displacement. Existing equilibrium
    /// stiffness must exclude the added relaxing stiffness. An empty list is
    /// an exact no-op. `maximum_memory_states` bounds the complete branch set
    /// (at most 256); the shared render owner also admits total state capacity
    /// and the complete mechanical/relaxation rate when stepping.
    ///
    /// # Errors
    /// Refuses motion or existing memory, exhausted capacity, momentum-loaded
    /// or malformed branches, nonphysical coefficients, and shared pHS/render
    /// admission failures. No branch or off-diagonal coupling is discarded.
    pub fn with_relaxation_branches(
        self,
        branches: Vec<fs_phs::RelaxationBranch>,
        maximum_memory_states: usize,
    ) -> Result<Self, AcousticRealizeError> {
        if branches.is_empty() {
            return Ok(self);
        }
        let base_dim = 2 * self.n_modes();
        if maximum_memory_states > 256 || branches.len() > maximum_memory_states
            || self.x.len() != base_dim || self.sys.state_dim() != base_dim
            || self.x.iter().any(|x| *x != 0.0)
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "nonlinear material memory needs one zero-state attachment within its declared capacity",
            });
        }
        if branches.iter().any(|branch| branch.projection.len() != base_dim
            || branch.projection.iter().enumerate().any(|(i, b)|
                !b.is_finite() || i % 2 == 1 && *b != 0.0))
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "nonlinear material memory needs finite displacement-only projections on the original modal basis",
            });
        }
        let sys = self.sys.with_relaxation_branches(branches)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let linear = prepare_render_linear(&sys)?;
        let x = vec![0.0; sys.state_dim()];
        Ok(Self { sys, linear, x, areas: self.areas })
    }

    /// Accepted state: mass-normalized `[q0, p0, q1, p1, ...]`, followed by
    /// constitutive and acoustic/cavity histories in their attachment order.
    /// Each load retains its declared coordinates; original modal pairs stay first.
    #[must_use]
    pub fn state(&self) -> &[f64] {
        &self.x
    }

    /// Complete mechanical, constitutive and attached acoustic storage [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.sys.hamiltonian(&self.x)
    }

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
        let (drive, areas) = ss_sine_ports(
            plate.length_m, plate.width_m,
            plate.density_kg_m3 * plate.thickness_m, &disp,
        );
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
        let rec = fs_phs::linear_exact::step_render(&self.sys, &self.linear, &self.x, &[force_n], dt)
            .map_err(render_step_error)?;
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

// The same mass-normalized sine basis as fs-nlmodal's nonlinear storage.
// Point drive is at (Lx/4, Ly/2); the face port integrates the entire surface.
fn ss_sine_ports(
    length_m: f64,
    width_m: f64,
    surface_density_kg_m2: f64,
    modes: &[fs_nlmodal::SineMode],
) -> (Vec<f64>, Vec<f64>) {
    // Integral sin²(m*pi*x/Lx) sin²(n*pi*y/Ly) dA = Lx*Ly/4.
    // Using the one-dimensional factor 2 halves physical force-to-flow gain.
    let norm = (4.0 / (surface_density_kg_m2 * length_m * width_m)).sqrt();
    let mut drive = Vec::with_capacity(modes.len());
    let mut areas = Vec::with_capacity(modes.len());
    for mode in modes {
        let mf = mode.m as f64;
        let nf = mode.n as f64;
        let pi = core::f64::consts::PI;
        areas.push(if mode.m % 2 == 1 && mode.n % 2 == 1 {
            norm * (2.0 * length_m / (mf * pi)) * (2.0 * width_m / (nf * pi))
        } else {
            0.0
        });
        drive.push(norm * det::sin(mf * pi * 0.25) * det::sin(nf * pi * 0.5));
    }
    (drive, areas)
}

fn finish_vk(
    storage: fs_nlmodal::SosModalStorage,
    zetas: &[f64],
    drive: &[f64],
    areas: Vec<f64>,
    with_area: bool,
) -> Result<VkBody, AcousticRealizeError> {
    finish_vk_damping(storage, VkDamping::Ratios(zetas), drive, areas, with_area)
}

enum VkDamping<'a> {
    Ratios(&'a [f64]),
    Matrix(&'a [f64]),
}

fn finish_vk_damping(
    storage: fs_nlmodal::SosModalStorage,
    damping: VkDamping<'_>,
    drive: &[f64],
    areas: Vec<f64>,
    with_area: bool,
) -> Result<VkBody, AcousticRealizeError> {
    storage.validate().map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    let n = storage.omegas.len();
    if drive.len() != n || areas.len() != n
        || drive.iter().chain(&areas).any(|x| !x.is_finite())
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "nonlinear plate needs finite force and surface projections on every storage mode",
        });
    }
    let m = if with_area { 2 } else { 1 };
    let mut g = vec![0.0; (2 * n) * m];
    for k in 0..n {
        g[(2 * k + 1) * m] = drive[k];
        if with_area {
            g[(2 * k + 1) * m + 1] = areas[k];
        }
    }
    let omegas = storage.omegas.clone();
    let sys = match damping {
        VkDamping::Ratios(zetas) =>
            fs_nlmodal::assemble_storage(n, &omegas, zetas, m, g, Box::new(storage)),
        VkDamping::Matrix(matrix) =>
            fs_nlmodal::assemble_storage_with_damping(n, &omegas, matrix, m, g, Box::new(storage)),
    }.map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    let linear = prepare_render_linear(&sys)?;
    Ok(VkBody {
        sys,
        linear,
        x: vec![0.0; 2 * n],
        areas,
    })
}

fn prepare_render_linear(sys: &fs_phs::PortHamiltonian)
    -> Result<fs_phs::linear_exact::LinearSplit, AcousticRealizeError>
{
    fs_phs::linear_exact::split_quadratic(sys, fs_phs::linear_exact::LinearSplitBudget {
        max_modes: sys.state_dim(), max_omega_dt: core::f64::consts::PI,
    }).map_err(render_step_error)
}

fn render_step_error(error: fs_phs::PhsError) -> AcousticRealizeError {
    match error {
        fs_phs::PhsError::TimeResolution(band) => AcousticRealizeError::TimeResolution(band),
        other => AcousticRealizeError::Nonlinear(other.to_string()),
    }
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
        bending_loss
            .as_ref()
            .map(|loss| loss as &dyn ModalBendingLoss),
    )
}

/// Mechanical and modal controls for an arbitrary, possibly heterogeneous chart.
/// Regional materials do not automatically select a damping law. The scalar
/// ratio is authored separately, and must be zero when source-bound bending
/// loss is selected. Thermal transport and nonlinear membranes are absent.
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
    bending_loss: Option<&dyn ModalBendingLoss>,
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
    let face = std::sync::Arc::new(());
    for (body, bending_zeta) in out.iter_mut().zip(bending_zetas) {
        body.zeta += bending_zeta;
        body.piston_face = Some(face.clone());
    }
    Ok(out)
}

// Weak-loss modal projection: only bending stores the strain energy associated
// with the material bending laws. Static geometric stiffness changes
// the mode and its total energy, but must not acquire that loss angle.
// See Fedorov et al., https://arxiv.org/abs/1807.07086, Eq. (1), Sec. II.
pub(crate) struct PlateBendingLaw {
    pub(crate) thermal: Option<ThermoelasticZener>,
    pub(crate) viscous: Option<IsotropicPlateBendingViscosity>,
    pub(crate) viscous_time_s: f64,
    pub(crate) thickness: f64,
}

impl PlateBendingLaw {
    fn new(plate: &ThinPlate) -> Result<Self, AcousticRealizeError> {
        let thermal = admitted_thermoelastic(plate)?;
        let mut viscous_time_s = 0.0;
        if let Some(viscous) = plate.kelvin_voigt_bending {
            if !has_isotropic_elasticity(plate) || plate.damping_ratio != 0.0 {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "Kelvin-Voigt plate loss requires isotropic elasticity and zero authored damping",
                });
            }
            viscous_time_s = material_loss::viscous_time(viscous, plate.e1_pa)?;
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

    pub(crate) fn zeta(&self, omega: f64) -> Result<f64, AcousticRealizeError> {
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
            zeta += crate::plate_viscosity::kelvin_voigt_zeta(self.viscous_time_s, omega);
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
}

// Uniform and regional losses share the same certified mode harvesting.
pub(crate) trait ModalBendingLoss {
    fn modal_zeta(
        &self,
        total: &PlateModel,
        phi: &[f64],
        omega: f64,
    ) -> Result<f64, AcousticRealizeError>;
}

impl ModalBendingLoss for PlateBendingLoss {
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
            piston_face: None,
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
    bending_viscosity: Option<Vec<RegionalBendingViscosity>>,
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

    /// Region index of a triangle in declared region order; the index must be
    /// a valid element of this chart.
    pub(crate) fn region_of_element(&self, element: usize) -> usize {
        self.region_for_element[element]
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
        let replacement = compile_plate_material_chart(
            self.chart.mesh.clone(),
            self.chart.boundary_nodes.clone(),
            inputs,
        )?;
        if self.bending_viscosity.is_some() {
            replacement.with_kelvin_voigt_bending_loss()
        } else {
            Ok(replacement)
        }
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
        bending_viscosity: None,
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
    /// Requires `kelvin_voigt_bending_viscosity` [Pa s] and rho/E/nu claims
    /// constant at fixed positive absolute temperature (scalars or T-only curves).
    /// The viscosity claim must declare a frequency band, intersected
    /// with every coefficient's band. Other state coordinates remain frozen.
    /// One typed frequency axis admits Hz or rad/s; legacy `omega` is rad/s.
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
        self.plate.kelvin_voigt_bending = Some(material_loss::resolve_viscosity(&self.material)?);
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

pub(crate) fn plate_elasticity(
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

pub(crate) fn plate_property(
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

/// Baffled-piston `Z(ω) = p/v` fitted as a passive discrete filter.
fn piston_radiation_filter(
    radius: f64,
    gas: &GasState,
    sample_rate_hz: u32,
) -> Result<DigitalFilter, AcousticRealizeError> {
    if !(radius > 0.0 && radius.is_finite() && sample_rate_hz > 0) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "piston radiation radius and sample rate",
        });
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    let nyquist = core::f64::consts::PI / dt;
    let omega_lo = 40.0_f64.min(0.05 * nyquist);
    let omega_hi = (0.40 * nyquist).max(omega_lo * 4.0);
    if !(omega_hi > omega_lo) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "piston radiation fit band",
        });
    }
    let n = 24usize;
    // Acoustic p/U → face p/v by the geometric area. Acoustics
    // e^{-iωt} → vfit e^{+iωt} by conjugation, without a fallback.
    let mut omega = Vec::with_capacity(n);
    let mut z = Vec::with_capacity(n);
    let area = core::f64::consts::PI * radius * radius;
    for k in 0..n {
        let t = k as f64 / (n as f64 - 1.0);
        let w = omega_lo * det::exp(t * det::ln(omega_hi / omega_lo));
        let (r, x) = fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, radius, w)
            .map_err(|_| AcousticRealizeError::InvalidDescription {
                what: "piston radiation impedance refused",
            })?;
        let z_face = C64::new(r, x).scale(area);
        omega.push(w);
        z.push(C64::new(z_face.re, -z_face.im));
    }
    if omega.len() < 8 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "piston radiation fit samples",
        });
    }
    let mut opts = FitOptions::new(3);
    opts.fit_e = false;
    opts.iterations = 8;
    realize_tabulated_impedance(&omega, &z, dt, &opts, omega[omega.len() / 2])
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
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

    /// G1: the 2-D modal mass integral and virtual work use the same basis.
    #[test]
    fn nonlinear_sine_ports_use_the_storage_mass_normalization() {
        let modes = [
            fs_nlmodal::SineMode { m: 1, n: 1 },
            fs_nlmodal::SineMode { m: 3, n: 1 },
            fs_nlmodal::SineMode { m: 2, n: 1 },
        ];
        let (lx, ly, mu) = (0.4, 0.3, 2.1);
        let (drive, area) = ss_sine_ports(lx, ly, mu, &modes);
        let pi = core::f64::consts::PI;
        let norm = drive[0] / (pi/4.0).sin();
        assert!((mu * lx * ly * norm * norm / 4.0 - 1.0).abs() < 2e-13);
        assert!((area[0] - norm*4.0*lx*ly/(pi*pi)).abs() < 2e-13);
        assert!((area[1] - area[0]/3.0).abs() < 2e-13);
        assert_eq!(area[2], 0.0);
        // At rest, a unit physical point force gives this face acceleration.
        // It contains BOTH projections; a mismatched normalization loses 1/2.
        let expected = 16.0 * (pi/4.0).sin() / (mu*pi*pi);
        assert!((drive[0] * area[0] - expected).abs() < 2e-13);
        // Same specimen dimensions and twice the surface density: all modes
        // keep unit modal mass, while force-to-flow gain scales inversely.
        let (heavy_drive, heavy_area) = ss_sine_ports(lx, ly, 2.0*mu, &modes);
        assert!((heavy_drive[0]*heavy_area[0] / (drive[0]*area[0]) - 0.5).abs() < 2e-13);
    }

    /// Complete baffled-piston face impedance [Pa s/m] under e^{+iwt},
    /// rho c [1 - J1(2ka)/ka + i H1(2ka)/ka], continued to complex ka by the
    /// J1 and Struve H1 power series (no fitted digits).
    fn piston_face_impedance(rho_c: f64, ka: C64) -> C64 {
        let q = -(ka * ka);
        let (mut j1, mut tj) = (C64::ZERO, ka);
        // Gamma(3/2) Gamma(5/2) = 3 pi / 8.
        let (mut h1, mut th) = (C64::ZERO, (ka * ka).scale(8.0 / (3.0 * core::f64::consts::PI)));
        for m in 0..80 {
            let m = f64::from(m);
            j1 = j1 + tj;
            h1 = h1 + th;
            tj = tj * q.scale(1.0 / ((m + 1.0) * (m + 2.0)));
            th = th * q.scale(1.0 / ((m + 1.5) * (m + 2.5)));
        }
        (C64::ONE - j1 / ka + C64::new(0.0, 1.0) * h1 / ka).scale(rho_c)
    }

    /// The fitted filter's face impedance at a complex frequency s, through
    /// its own biquads at z^-1 = e^{-s dt}.
    fn filter_at(filter: &DigitalFilter, s: C64) -> C64 {
        let e = (-s).scale(filter.t_s);
        let w = C64::new(det::exp(e.re) * det::cos(e.im), det::exp(e.re) * det::sin(e.im));
        filter.sections.iter().fold(C64::from_re(filter.direct), |acc, b| {
            acc + (C64::from_re(b.b[0]) + w.scale(b.b[1]) + (w * w).scale(b.b[2]))
                / (C64::ONE + w.scale(b.a[0]) + (w * w).scale(b.a[1]))
        })
    }

    /// Free-decay root s = -sigma + i omega of m s^2 + s A Z(s) + m w0^2 = 0
    /// for a face impedance Z, by Newton from `guess`.
    fn decay_root(mass: f64, omega0: f64, area: f64, guess: C64, z: impl Fn(C64) -> C64) -> C64 {
        let f = |s: C64| s * s.scale(mass) + s * z(s).scale(area) + C64::from_re(mass * omega0 * omega0);
        let mut s = guess;
        for _ in 0..100 {
            let h = s.scale(1e-7);
            let step = f(s) * h / (f(s + h) - f(s));
            s = s - step;
            if step.abs() <= 1e-15 * s.abs() {
                break;
            }
        }
        s
    }

    /// Ring a body down from one impulse step and measure its decay rate
    /// (least squares on the log of each cycle's peak) and frequency (from
    /// interpolated upward zero crossings), over 60 ms.
    fn ring_down(mut body: CompactBody, rate: u32) -> (f64, f64) {
        let dt = 1.0 / f64::from(rate);
        body.drive(1.0, dt).unwrap();
        let v: Vec<f64> = (0..(rate as usize * 6 / 100))
            .map(|_| {
                body.drive(0.0, dt).unwrap();
                body.v
            })
            .collect();
        let mut crossings = Vec::new();
        for i in 1..v.len() {
            if v[i - 1] < 0.0 && v[i] >= 0.0 {
                crossings.push((i - 1) as f64 + v[i - 1] / (v[i - 1] - v[i]));
            }
        }
        let cycles = crossings.len() - 1;
        let omega = core::f64::consts::TAU * cycles as f64 / ((crossings[cycles] - crossings[0]) * dt);
        let peaks: Vec<(f64, f64)> = crossings
            .windows(2)
            .map(|w| {
                let (a, b) = (w[0].ceil() as usize, w[1].floor() as usize);
                let (k, peak) = (a..=b)
                    .map(|k| (k, v[k]))
                    .fold((a, f64::MIN), |m, p| if p.1 > m.1 { p } else { m });
                (k as f64 * dt, peak.ln())
            })
            .collect();
        let n = peaks.len() as f64;
        let (st, sy) = peaks.iter().fold((0.0, 0.0), |s, p| (s.0 + p.0, s.1 + p.1));
        let (stt, sty) = peaks.iter().fold((0.0, 0.0), |s, p| (s.0 + p.0 * p.0, s.1 + p.0 * p.1));
        (-(n * sty - st * sy) / (n * stt - st * st), omega)
    }

    /// Amount of air loading (.362, rulings 6231 and 6237). A radiation-dominated
    /// body (0.05 m^2, 0.05 kg, 400 Hz, no structural damping, ka near 0.9 in
    /// air) is loaded by each gas's piston, whose radius sqrt(A / pi) is the
    /// declared face area's equivalent circle.
    /// - Oracle: the exact complex root of m s^2 + s A Z(s) + k = 0, with the
    ///   closed-form piston impedance continued to complex frequency. The
    ///   first-order real-axis value is printed beside it: it is biased by up
    ///   to 3.3 %, which is why the first declared band (decay 3 %, frequency
    ///   0.3 % against it at 192 kHz) failed at 2 atm (+3.05 %, -0.31 %).
    /// - Fit term: the same root through the fitted filter alone, measured on
    ///   the fit, not on the render. Its R and X errors at the loaded
    ///   frequency are printed, and so is their maximum over the fit band,
    ///   which is uninformative at the band's low edge where R -> 0.
    /// - Lag term: residual minus fit, the held reaction's O(dt) lag. On the
    ///   96/192/384/768 kHz ladder it must shrink at observed order >= 0.8 per
    ///   doubling, for decay and frequency, in every gas.
    /// - Gate at 768 kHz: |residual| <= |fit| + |lag(384) - lag(768)| /
    ///   (2^0.8 - 1), the remaining lag by Richardson at the declared minimum
    ///   order 0.8 (ruling 6253). The first form, |lag(384) - lag(768)|, is
    ///   that bound at order exactly 1, a knife edge for a first-order lag: it
    ///   failed by roundoff on air's frequency (3e-7 absolute) and He-O2's
    ///   decay (0.0014 %). Direction strict at every rate: 2 atm and CO2 damp
    ///   more than air, He-O2 less.
    /// - Planted negative: at 1e-9 of air's density the load vanishes and the
    ///   body rings at its structural zeta omega0.
    #[test]
    fn piston_load_damps_and_mass_loads_the_body_by_the_gas_impedance() {
        use fs_material::gas_species::{GasSpecies, resolve_gas_mixture};
        use fs_qty::semantic::{Composition, CompositionBasis};
        let air = fs_material::gas::GasSpec::dry_air_ussa1976();
        let mixture = |species: &[GasSpecies], fractions: Vec<f64>| {
            resolve_gas_mixture(
                species,
                &Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
                293.15,
                101_325.0,
            )
            .unwrap()
        };
        let gases = [
            ("air", GasState::try_new(&air, 293.15, 101_325.0).unwrap()),
            ("2 atm", GasState::try_new(&air, 293.15, 202_650.0).unwrap()),
            ("CO2", mixture(&[GasSpecies::CarbonDioxide], vec![1.0])),
            ("He-O2", mixture(&[GasSpecies::Helium, GasSpecies::Oxygen], vec![0.8, 0.2])),
        ];
        let rates = [96_000_u32, 192_000, 384_000, 768_000];
        let (area, mass, omega0) = (0.05, 0.05, core::f64::consts::TAU * 400.0);
        let radius = (area / core::f64::consts::PI).sqrt();
        let body = |zeta: f64| {
            CompactBody::from_radiator(RadiatingPlate {
                area_m2: area,
                mass_kg: mass,
                frequency_hz: 400.0,
                damping_ratio: zeta,
            })
            .unwrap()
        };
        // measured[gas][rate] = (sigma, omega)
        let mut measured = vec![Vec::new(); gases.len()];
        for (g, (name, gas)) in gases.iter().enumerate() {
            let rho_c = gas.density * gas.sound_speed;
            let ka = |s: C64| (C64::new(0.0, -1.0) * s).scale(radius / gas.sound_speed);
            // First-order real-axis oracle, kept for the record.
            let mut first = (0.0, omega0);
            for _ in 0..50 {
                let (r, x) = fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, radius, first.1).unwrap();
                let m_eff = mass - x * area * area / first.1;
                let sigma = r * area * area / (2.0 * m_eff);
                first = (sigma, (mass * omega0 * omega0 / m_eff - sigma * sigma).sqrt());
            }
            let exact = decay_root(mass, omega0, area, C64::new(-first.0, first.1), |s| piston_face_impedance(rho_c, ka(s)));
            let (sigma_e, omega_e) = (-exact.re, exact.im);
            let closed_real = piston_face_impedance(rho_c, C64::from_re(omega_e * radius / gas.sound_speed));
            let (mut lag_sigma, mut lag_omega, mut fits) = (Vec::new(), Vec::new(), Vec::new());
            for &rate in &rates {
                let mut loaded = body(0.0);
                loaded.attach_piston_load(gas, rate).unwrap();
                let filter = &loaded.rad.as_ref().unwrap().filter;
                let fitted = decay_root(mass, omega0, area, exact, |s| filter_at(filter, s));
                let at = filter.eval(omega_e);
                // The fit's own R and X errors over its band (40 rad/s to 0.4 Nyquist).
                let (lo, hi) = (40.0_f64, 0.4 * core::f64::consts::PI * f64::from(rate));
                let (mut band_r, mut band_x) = (0.0_f64, 0.0_f64);
                for i in 0..=200 {
                    let w = lo * (hi / lo).powf(f64::from(i) / 200.0);
                    let closed = piston_face_impedance(rho_c, C64::from_re(w * radius / gas.sound_speed));
                    let fit = filter.eval(w);
                    band_r = band_r.max((fit.re - closed.re).abs() / closed.re.abs());
                    band_x = band_x.max((fit.im - closed.im).abs() / closed.im.abs());
                }
                let (sigma_m, omega_m) = ring_down(loaded, rate);
                let residual = (sigma_m / sigma_e - 1.0, omega_m / omega_e - 1.0);
                let fit = (-fitted.re / sigma_e - 1.0, fitted.im / omega_e - 1.0);
                lag_sigma.push(residual.0 - fit.0);
                lag_omega.push(residual.1 - fit.1);
                fits.push((residual, fit));
                measured[g].push((sigma_m, omega_m));
                println!(
                    "{{\"case\":\"piston-air-loading\",\"gas\":\"{name}\",\"rate\":{rate},\"rho\":{},\"c\":{},\"ka\":{},\"sigma\":{{\"measured\":{sigma_m},\"exact_root\":{sigma_e},\"first_order\":{},\"fit_root\":{}}},\"omega\":{{\"measured\":{omega_m},\"exact_root\":{omega_e},\"first_order\":{},\"fit_root\":{}}},\"residual\":[{},{}],\"fit_term\":[{},{}],\"lag\":[{},{}],\"fit_error_at_omega\":[{},{}],\"fit_error_over_band\":[{band_r},{band_x}]}}",
                    gas.density,
                    gas.sound_speed,
                    omega_e * radius / gas.sound_speed,
                    first.0,
                    -fitted.re,
                    first.1,
                    fitted.im,
                    residual.0,
                    residual.1,
                    fit.0,
                    fit.1,
                    residual.0 - fit.0,
                    residual.1 - fit.1,
                    (at.re - closed_real.re) / closed_real.re,
                    (at.im - closed_real.im) / closed_real.im
                );
            }
            for lags in [&lag_sigma, &lag_omega] {
                for pair in lags.windows(2) {
                    let order = (pair[0].abs() / pair[1].abs()).log2();
                    assert!(order >= 0.8, "{name}: lag ladder {lags:?} has order {order}");
                }
            }
            let ((residual, fit), last) = (fits[3], rates.len() - 1);
            let richardson = 2.0_f64.powf(0.8) - 1.0;
            assert!(
                residual.0.abs()
                    <= fit.0.abs() + (lag_sigma[last - 1] - lag_sigma[last]).abs() / richardson,
                "{name}: decay residual {} against fit {} and lag {lag_sigma:?}",
                residual.0,
                fit.0
            );
            assert!(
                residual.1.abs()
                    <= fit.1.abs() + (lag_omega[last - 1] - lag_omega[last]).abs() / richardson,
                "{name}: frequency residual {} against fit {} and lag {lag_omega:?}",
                residual.1,
                fit.1
            );
        }
        for r in 0..rates.len() {
            let sigma = |g: usize| measured[g][r].0;
            assert!(sigma(1) > sigma(0) && sigma(2) > sigma(0) && sigma(3) < sigma(0), "rate {}", rates[r]);
        }
        // Planted negative: no gas to speak of, only the structure's own loss.
        let thin = GasState::try_new(&air, 293.15, 101_325.0*1e-9).unwrap();
        let mut structural = body(0.02);
        structural.attach_piston_load(&thin, 192_000).unwrap();
        let (sigma, _) = ring_down(structural, 192_000);
        assert!((sigma / (0.02 * omega0) - 1.0).abs() <= 1e-4, "structural decay {sigma}");
    }

    #[test]
    fn explicit_no_exterior_surface_retains_mechanics_without_a_fictitious_piston() {
        use super::*;
        let mut body = CompactBody::from_mass_normalized_projection(1000.0, 0.0, 2.0, 0.0).unwrap();
        let acceleration = body.drive(1.0, 1.0 / 48_000.0).unwrap();
        assert!(acceleration != 0.0);
        assert_eq!(body.radiate_mass_normalized(acceleration, 1.2, 1.0), 0.0);
        let gas = GasState::try_new(&fs_material::gas::GasSpec::dry_air_ussa1976(),
            293.15, 101_325.0).unwrap();
        assert!(body.attach_piston_load(&gas, 48_000).is_err());
        assert!(CompactBody::from_mass_normalized_projection(1000.0, 0.1, 2.0, 0.0).is_err());
    }

    #[test]
    fn plate_bank_cavity_island_is_invariant_under_time_partition() {
        use super::*;
        let gas = GasState::try_new(&fs_material::gas::GasSpec::dry_air_ussa1976(),293.15,101_325.0).unwrap();
        let make = || {
            let mut bank = PlateBank::default();
            bank.attach_cavity(fs_scenario::HelmholtzCavity {volume_m3:0.01,
                neck_radius_m:0.01,neck_length_m:0.02,inner_end:fs_scenario::NeckEnd::Flanged,
                outer_end:fs_scenario::NeckEnd::Flanged},&gas).unwrap();
            bank.cavity.as_mut().unwrap().x[0] = 1e-6;
            bank
        };
        for dt in [1e-4,5e-4,1e-3] {
            let mut whole = make(); let mut half = make();
            whole.drive_and_radiate(0.0,dt,gas.density,1.0).unwrap();
            half.drive_and_radiate(0.0,dt/2.0,gas.density,1.0).unwrap();
            half.drive_and_radiate(0.0,dt/2.0,gas.density,1.0).unwrap();
            for (a,b) in whole.cavity.as_ref().unwrap().x.iter().zip(&half.cavity.as_ref().unwrap().x) {
                assert!((a-b).abs()<=1e-12*a.abs().max(b.abs()).max(1e-12),"{a} {b}");
            }
        }
    }

    #[test]
    fn vk_render_has_an_exact_linear_map_and_a_typed_nonlinear_band() {
        use fs_nlmodal::{SosModalStorage, StressChannel};
        use super::*;
        let make = |nonlinear| finish_vk(SosModalStorage {
            omegas: vec![1.0],
            channels: if nonlinear {vec![StressChannel {coefficient:1.0,coupling:vec![1.0]}]}
                else {vec![]},
        }, &[0.0], &[1.0], vec![1.0], false).unwrap();
        for z in [0.1, 1.0, 2.5] {
            let mut body = make(false);
            body.x[0] = 1e-4;
            body.drive_and_radiate(0.0,z,1.2,1.0).unwrap();
            assert!((body.x[0]/1e-4-det::cos(z)).abs()<1e-12);
            assert!((body.x[1]/1e-4+det::sin(z)).abs()<1e-12);
        }
        let cap = fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP;
        for z in [cap/4.0,cap/2.0,cap] {
            let mut body = make(true);
            body.x[0] = 1e-4;
            body.drive_and_radiate(0.0,z,1.2,1.0).unwrap();
            let phase = det::atan2(-body.x[1],body.x[0]);
            let cents = 1200.0*det::ln(phase/z)/core::f64::consts::LN_2;
            assert!(cents.abs()<=1.0,"{cents}");
        }
        let mut body = make(true);
        let before = body.x.clone();
        assert!(matches!(body.drive_and_radiate(0.0,cap+1e-12,1.2,1.0),
            Err(AcousticRealizeError::TimeResolution(_))));
        assert_eq!(body.x,before);
    }

    use super::*;

    fn compact_zoh_fixture(zeta: f64) -> CompactBody {
        CompactBody::from_radiator(RadiatingPlate {
            area_m2: 0.03,
            mass_kg: 0.25,
            frequency_hz: 1000.0 / core::f64::consts::TAU,
            damping_ratio: zeta,
        })
        .unwrap()
    }

    /// G1: independent constant-force solution, including endpoint pressure.
    /// The step exceeds the undamped Euler stability limit omega*dt=2.
    #[test]
    fn compact_zoh_matches_undamped_oscillator_and_endpoint_pressure() {
        let mut body = compact_zoh_fixture(0.0);
        let dt = 0.003;
        let force = 2.0;
        let acceleration = force / body.mass_kg;
        let angle = body.omega * dt;
        let expected_y = acceleration / body.omega.powi(2) * (1.0 - angle.cos());
        let expected_v = acceleration / body.omega * angle.sin();
        let expected_acc = acceleration * angle.cos();
        let pressure = body.drive_and_radiate(force, dt, 1.2, 2.0).unwrap();
        assert!((body.y - expected_y).abs() < 1.0e-15);
        assert!((body.v - expected_v).abs() < 1.0e-12);
        let expected_pressure = 1.2 * 0.03 * expected_acc / (4.0 * core::f64::consts::PI);
        assert!((pressure - expected_pressure).abs() < 1.0e-12);
        // Prior semi-implicit Euler gives y = a*dt^2, not this solution.
        assert!((acceleration * dt * dt - expected_y).abs() > 1.0e-5);
    }

    /// G1/G3: damped analytical step and held-force subdivision invariance.
    /// Subdivision equivalence applies without a sampled radiation filter.
    #[test]
    fn compact_zoh_damped_solution_and_subdivision_agree() {
        let dt = 0.003;
        for zeta in [0.0, 0.2, 1.0, 2.0] {
            let mut whole = compact_zoh_fixture(zeta);
            let mut split = whole.clone();
            let acc = whole.drive(2.0, dt).unwrap();
            let mut split_acc = 0.0;
            for _ in 0..16 {
                split_acc = split.drive(2.0, dt / 16.0).unwrap();
            }
            assert!((whole.y - split.y).abs() < 1.0e-14, "zeta={zeta}");
            assert!((whole.v - split.v).abs() < 1.0e-11, "zeta={zeta}");
            assert!((acc - split_acc).abs() < 1.0e-8, "zeta={zeta}");
            if zeta == 0.2 {
                let wd = whole.omega * (1.0 - zeta * zeta).sqrt();
                let decay = (-zeta * whole.omega * dt).exp();
                let expected_y = 2.0 / (whole.mass_kg * whole.omega.powi(2))
                    * (1.0 - decay * ((wd * dt).cos() + zeta * whole.omega / wd * (wd * dt).sin()));
                let expected_v = 2.0 / (whole.mass_kg * wd) * decay * (wd * dt).sin();
                assert!((whole.y - expected_y).abs() < 1.0e-14);
                assert!((whole.v - expected_v).abs() < 1.0e-11);
            }
        }
    }

    #[test]
    fn compact_zoh_refuses_invalid_step_before_changing_motion() {
        let mut body = compact_zoh_fixture(0.1);
        body.drive(2.0, 0.001).unwrap();
        let before = (body.y.to_bits(), body.v.to_bits());
        for (force, dt) in [
            (f64::NAN, 0.001),
            (1.0, 0.0),
            (1.0, -0.1),
            (1.0, f64::INFINITY),
        ] {
            assert!(body.drive(force, dt).is_err());
            assert_eq!((body.y.to_bits(), body.v.to_bits()), before);
        }
    }

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
                    body.attach_piston_load(&gas, 48_000).expect("piston load");
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
        let (r, x) =
            fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, 0.05, 2.0e3).expect("z");
        assert!(r > 0.0, "resistance must be positive");
        assert!(x < 0.0, "mass-like reactance is negative under e^{{-iωt}}");
    }

    /// G0: an invalid gas used to fall through the Option into the low-ka
    /// formula. It now refuses by name; a physical gas at the same ka still
    /// has a well-defined closed-form load.
    #[test]
    fn thin_plate_refuses_instead_of_low_ka_fallback() {
        let mut gas = GasState::try_new(
            &fs_material::gas::GasSpec::dry_air_ussa1976(),
            293.15,
            101_325.0,
        )
        .unwrap();
        let ka = 1e-3;
        assert!(
            fs_phs::baffled_piston_impedance(
                gas.density,
                gas.sound_speed,
                0.05,
                ka * gas.sound_speed / 0.05
            )
            .is_ok()
        );
        // The former ring path returned None at non-finite wave number,
        // then silently manufactured a low-ka value using the same gas.
        gas.sound_speed = 0.0;
        let result = piston_radiation_filter(0.05, &gas, 48_000);
        assert!(matches!(
            result,
            Err(AcousticRealizeError::InvalidDescription {
                what: "piston radiation impedance refused"
            })
        ));
        assert!(piston_radiation_filter(f64::NAN, &gas, 48_000).is_err());
        assert!(piston_radiation_filter(0.05, &gas, 0).is_err());
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

// Source-bound regional bending loss uses the same DKT and modal owners.
mod material_loss;
pub use material_loss::{
    OrthotropicPlateBendingViscosity, RegionalBendingViscosity, certified_material_chart_radiators,
};
pub(crate) use material_loss::viscous_time;
