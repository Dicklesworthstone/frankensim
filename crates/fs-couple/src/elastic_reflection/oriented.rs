//! Principal mechanical supports rotated into the original elastic surface frame.

use super::{ElasticSupportError, SupportStorage, displacement, support_storage};
use fs_ga::So3;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::elastic::acoustic::interface::ElasticMode;
use fs_material::elastic::acoustic::interface::surface::{
    ElasticSurfaceCondition, ElasticSurfaceImpedanceTensor, ElasticSurfaceResponse,
    IsotropicElasticSurface, PreparedElasticSurface,
};
use fs_math::{c64::C64, dd::Dd, det};
use fs_phs::{WallPin, wall_specific_impedance};

/// Three original mass/spring/damper supports in a shared principal frame.
///
/// `principal_to_surface` maps principal components into Cartesian surface
/// coordinates (x/y tangential, z outward normal). Its columns are the support
/// directions. The original `So3` owns rotation admission; the original wall
/// law owns pin/frequency admission when the support is evaluated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrientedOscillatorSupport {
    principal: [WallPin; 3],
    principal_to_surface: So3,
}

impl OrientedOscillatorSupport {
    /// Retain the authored pins and checked rotation without fitting or rescaling.
    pub const fn new(principal: [WallPin; 3], principal_to_surface: So3) -> Self {
        Self {
            principal,
            principal_to_surface,
        }
    }

    pub const fn principal_supports(&self) -> &[WallPin; 3] {
        &self.principal
    }

    pub const fn principal_to_surface(&self) -> So3 {
        self.principal_to_surface
    }

    /// The original reciprocal support tensor `R diag(Z_i) R^T` [Pa s/m].
    /// This changes reflected polarization when principal support directions
    /// are not aligned with the surface's normal and incident-wave plane.
    pub fn impedance_at_frequency(
        &self,
        omega_rad_s: f64,
    ) -> Result<ElasticSurfaceImpedanceTensor, ElasticSupportError> {
        let mut principal = [C64::ZERO; 3];
        for (pin, value) in self.principal.iter().zip(&mut principal) {
            *value = wall_specific_impedance(pin, omega_rad_s)?;
        }
        let rotation = self.principal_to_surface.matrix().m;
        // Keep the authored dissipative channels. Refactoring the rounded
        // Cartesian resistance could invent a weak positive eigenvalue in an
        // exactly lossless principal direction.
        let gram = core::array::from_fn(|row| {
            core::array::from_fn(|axis| rotation[row * 3 + axis] * det::sqrt(principal[axis].re))
        });
        let mut reactance = [[0.0; 3]; 3];
        for row in 0..3 {
            for col in row..3 {
                let mut imaginary = Dd::ZERO;
                for axis in 0..3 {
                    let weight = Dd::from_f64(rotation[row * 3 + axis])
                        * Dd::from_f64(rotation[col * 3 + axis]);
                    imaginary = imaginary + weight * Dd::from_f64(principal[axis].im);
                }
                reactance[row][col] = imaginary.to_f64();
                reactance[col][row] = imaginary.to_f64();
            }
        }
        Ok(ElasticSurfaceImpedanceTensor::from_resistance_gram(
            gram, reactance,
        )?)
    }
}

/// Cached P/SV/SH reflection from three oriented original mechanical supports.
///
/// This retains the same stationary, linear, locally reacting boundary domain
/// as the two-pin owner. It adds direction-dependent response and cross-mode
/// conversion; it does not resolve support layers or advance a transient clock.
#[derive(Debug)]
pub struct OrientedOscillatorSupportedSurface<'a> {
    state: &'a ElasticAcousticState,
    support: OrientedOscillatorSupport,
    omega_rad_s: f64,
    rotation: [f64; 9],
    reflection: PreparedElasticSurface<'a>,
}

/// Coherent native reflection and storage resolved along the actual supports.
#[derive(Clone, Debug)]
pub struct OrientedSupportedElasticResponse {
    /// The original native field and sole dissipative support-power account.
    pub reflection: ElasticSurfaceResponse,
    /// Peak velocities along the three principal support directions [m/s].
    pub principal_velocity_m_s: [C64; 3],
    /// Peak displacements along those same principal directions [m].
    pub principal_displacement_m: [C64; 3],
    /// Peak displacement in the Cartesian surface frame [m].
    pub support_displacement_m: [C64; 3],
    /// Cycle-mean kinetic storage of each original principal pin [J/m²].
    pub mean_kinetic_energy_j_m2: [f64; 3],
    /// Cycle-mean spring storage of each original principal pin [J/m²].
    pub mean_elastic_energy_j_m2: [f64; 3],
    /// Total kinetic and spring storage over the three supports [J/m²].
    pub mean_stored_energy_j_m2: f64,
}

impl<'a> OrientedOscillatorSupportedSurface<'a> {
    pub fn prepare(
        state: &'a ElasticAcousticState,
        support: OrientedOscillatorSupport,
        omega_rad_s: f64,
        mode: ElasticMode,
        incidence_cosine: f64,
    ) -> Result<Self, ElasticSupportError> {
        Self::prepare_controlled(
            state,
            support,
            omega_rad_s,
            mode,
            incidence_cosine,
            &mut || false,
        )
    }

    /// Admit the original tensor and prepare one bounded native reflection solve.
    /// Cancellation publishes no partially prepared owner.
    pub fn prepare_controlled(
        state: &'a ElasticAcousticState,
        support: OrientedOscillatorSupport,
        omega_rad_s: f64,
        mode: ElasticMode,
        incidence_cosine: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElasticSupportError> {
        if cancelled() {
            return Err(ElasticSupportError::Cancelled);
        }
        let impedance = support.impedance_at_frequency(omega_rad_s)?;
        let reflection = IsotropicElasticSurface::new(
            state,
            ElasticSurfaceCondition::TensorImpedance(impedance),
        )?
        .prepare_incidence_cancellable(mode, incidence_cosine, cancelled)?;
        Ok(Self {
            state,
            support,
            omega_rad_s,
            rotation: support.principal_to_surface.matrix().m,
            reflection,
        })
    }

    pub const fn state(&self) -> &'a ElasticAcousticState {
        self.state
    }

    pub const fn support(&self) -> &OrientedOscillatorSupport {
        &self.support
    }

    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    pub const fn reflection(&self) -> &PreparedElasticSurface<'a> {
        &self.reflection
    }

    /// Reuse the complete native modal/velocity maps and the original pin energies.
    /// This performs fixed-size matrix products only, without a new solve or allocation.
    pub fn drive(
        &self,
        incoming_velocity_m_s: [C64; 3],
    ) -> Result<OrientedSupportedElasticResponse, ElasticSupportError> {
        let reflection = self.reflection.drive(incoming_velocity_m_s)?;
        let motion = support_motion(
            self.support.principal,
            self.rotation,
            reflection.boundary_velocity_m_s,
            self.omega_rad_s,
        )?;
        Ok(OrientedSupportedElasticResponse {
            reflection,
            principal_velocity_m_s: motion.principal_velocity,
            principal_displacement_m: motion.storage.displacement,
            support_displacement_m: motion.surface_displacement,
            mean_kinetic_energy_j_m2: motion.storage.kinetic,
            mean_elastic_energy_j_m2: motion.storage.elastic,
            mean_stored_energy_j_m2: motion.storage.total,
        })
    }
}

pub(super) struct SupportMotion {
    pub principal_velocity: [C64; 3],
    pub surface_displacement: [C64; 3],
    pub storage: SupportStorage,
}

// Direct and coated supports use exactly the same original pin coordinates and
// potentials. The caller supplies the actual SUPPORT velocity, which differs
// from the incident face velocity when a finite coating intervenes.
pub(super) fn support_motion(
    principal: [WallPin; 3],
    rotation: [f64; 9],
    velocity: [C64; 3],
    omega: f64,
) -> Result<SupportMotion, ElasticSupportError> {
    let principal_velocity = core::array::from_fn(|axis| {
        let mut real = Dd::ZERO;
        let mut imaginary = Dd::ZERO;
        for row in 0..3 {
            let weight = Dd::from_f64(rotation[row * 3 + axis]);
            real = real + weight * Dd::from_f64(velocity[row].re);
            imaginary = imaginary + weight * Dd::from_f64(velocity[row].im);
        }
        C64::new(real.to_f64(), imaginary.to_f64())
    });
    let storage = support_storage(principal, principal_velocity, omega)?;
    let mut surface_displacement = [C64::ZERO; 3];
    for (value, velocity) in surface_displacement.iter_mut().zip(velocity) {
        *value = displacement(velocity, omega)?;
    }
    Ok(SupportMotion {
        principal_velocity,
        surface_displacement,
        storage,
    })
}
