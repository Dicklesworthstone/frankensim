//! Coherent and transient elastic reflection from original passive wall mechanics.
//!
//! The solid keeps its source-bound stiffness, density, identity and strain domain.
//! Normal and tangential supports use the same [`WallPin`] and
//! [`wall_specific_impedance`] as native port-Hamiltonian acoustic walls. Their
//! mass, spring and dashpot act in parallel on the surface displacement, with
//! `traction = -Z velocity` and `Z = r - i omega sigma + i K/omega`.
//!
//! These are stationary, locally reacting, linear boundaries. The two-pin
//! owner has isotropic tangential response; the oriented owner retains three
//! independent principal supports and their original rotation. A pin is an
//! explicitly authored reduced support. The coating owner resolves an original
//! finite elastic layer in front of that support, including its full interior
//! fields; the pins themselves do not imply a measured material or finite shell.
//! The normal-incidence transient owner evolves their original displacement and
//! velocity, including causal re-emission and integrated physical work.
//! The viscoelastic coating owner also retains original bulk/shear material heat
//! and converts the pin phasor convention explicitly for native fluid incidence.

mod coating;
pub use coating::{OscillatorBackedElasticCoating, SupportedElasticCoatingResponse};
mod normal;
pub use normal::{
    NormalOscillatorBackedElasticCoating, NormalOscillatorSupportedSurface,
    NormalSupportedElasticCoatingResponse, NormalSupportedElasticResponse,
};
mod oriented;
pub use oriented::{
    OrientedOscillatorSupport, OrientedOscillatorSupportedSurface, OrientedSupportedElasticResponse,
};
mod transient;
pub use transient::{
    PreparedTransientElasticReflection, TransientElasticReflectionError,
    TransientElasticReflectionFrame, TransientElasticReflectionLimits,
    TransientElasticReflectionWork, TransientElasticReflectionWorkspace,
    TransientElasticSupportEnergy, TransientElasticSupportState,
};
mod viscoelastic;
pub use viscoelastic::{
    OscillatorBackedViscoelasticHeating, OscillatorBackedViscoelasticLayers,
    SupportedViscoelasticResponse,
};

use fs_material::MaterialError;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::elastic::acoustic::interface::ElasticMode;
use fs_material::elastic::acoustic::viscoelastic::ViscoelasticSlabError;
use fs_material::elastic::acoustic::interface::surface::{
    ElasticSurfaceCondition, ElasticSurfaceImpedance, ElasticSurfaceResponse,
    IsotropicElasticSurface, PreparedElasticSurface,
};
use fs_math::{c64::C64, dd::Dd, det};
use fs_phs::{PhsError, WallPin, wall_specific_impedance};

/// A native material/support refusal, cancellation, or unrepresentable storage.
#[derive(Debug)]
pub enum ElasticSupportError {
    Material(MaterialError),
    Viscoelastic(ViscoelasticSlabError),
    Support(PhsError),
    Cancelled,
    NumericalRange,
}

impl From<MaterialError> for ElasticSupportError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<ViscoelasticSlabError> for ElasticSupportError {
    fn from(error: ViscoelasticSlabError) -> Self {
        Self::Viscoelastic(error)
    }
}
impl From<PhsError> for ElasticSupportError {
    fn from(error: PhsError) -> Self {
        Self::Support(error)
    }
}
impl core::fmt::Display for ElasticSupportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "elastic support material: {error}"),
            Self::Viscoelastic(error) => write!(f, "elastic support coating: {error}"),
            Self::Support(error) => write!(f, "elastic support mechanics: {error}"),
            Self::Cancelled => f.write_str("elastic support preparation cancelled"),
            Self::NumericalRange => f.write_str("elastic support storage is unrepresentable"),
        }
    }
}
impl core::error::Error for ElasticSupportError {}

/// Prepared reflection retaining both original mechanical supports and frequency.
/// The boundary factorization is reused for every coherent P/SV/SH drive.
#[derive(Debug)]
pub struct OscillatorSupportedSurface<'a> {
    state: &'a ElasticAcousticState,
    normal: WallPin,
    tangential: WallPin,
    omega_rad_s: f64,
    reflection: PreparedElasticSurface<'a>,
}

/// Native reflected waves, support work and cycle-mean mechanical storage.
/// Cartesian arrays follow the surface frame: x/y tangential, z outward normal.
#[derive(Clone, Debug)]
pub struct SupportedElasticResponse {
    /// Original modal amplitudes, physical boundary traces and absorbed W/m².
    /// Its positive support absorption is the sole dissipative power account.
    pub reflection: ElasticSurfaceResponse,
    /// Peak displacement from the same accepted velocity, `x = i v / omega` [m].
    pub support_displacement_m: [C64; 3],
    /// Cycle-mean `sigma |v|²/4` in each local support axis [J/m²].
    pub mean_kinetic_energy_j_m2: [f64; 3],
    /// Cycle-mean `K |x|²/4` in each local support axis [J/m²].
    pub mean_elastic_energy_j_m2: [f64; 3],
    /// Sum of original kinetic and elastic storage over all support axes [J/m²].
    pub mean_stored_energy_j_m2: f64,
}

impl<'a> OscillatorSupportedSurface<'a> {
    /// Prepare the native solid boundary at the original support frequency.
    /// Pin validity and frequency admission belong to `fs-phs`; material and
    /// wave admission belong to `fs-material`. Resonance is allowed, including
    /// a lossless pin's exact zero impedance. No resistance floor is introduced.
    pub fn prepare(
        state: &'a ElasticAcousticState,
        normal: WallPin,
        tangential: WallPin,
        omega_rad_s: f64,
        mode: ElasticMode,
        incidence_cosine: f64,
    ) -> Result<Self, ElasticSupportError> {
        Self::prepare_controlled(
            state,
            normal,
            tangential,
            omega_rad_s,
            mode,
            incidence_cosine,
            &mut || false,
        )
    }

    /// Transactional preparation with the native bounded cancellation points.
    pub fn prepare_controlled(
        state: &'a ElasticAcousticState,
        normal: WallPin,
        tangential: WallPin,
        omega_rad_s: f64,
        mode: ElasticMode,
        incidence_cosine: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ElasticSupportError> {
        if cancelled() {
            return Err(ElasticSupportError::Cancelled);
        }
        let support = ElasticSurfaceImpedance::new(
            wall_specific_impedance(&normal, omega_rad_s)?,
            wall_specific_impedance(&tangential, omega_rad_s)?,
        )?;
        let reflection =
            IsotropicElasticSurface::new(state, ElasticSurfaceCondition::Impedance(support))?
                .prepare_incidence_cancellable(mode, incidence_cosine, cancelled)?;
        Ok(Self {
            state,
            normal,
            tangential,
            omega_rad_s,
            reflection,
        })
    }

    pub const fn state(&self) -> &'a ElasticAcousticState {
        self.state
    }
    pub const fn normal_support(&self) -> &WallPin {
        &self.normal
    }
    pub const fn tangential_support(&self) -> &WallPin {
        &self.tangential
    }
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    pub const fn reflection(&self) -> &PreparedElasticSurface<'a> {
        &self.reflection
    }

    /// Apply peak coherent modal velocities [m/s] and retain mechanical storage.
    /// Reuses the original reflection and work solve, without allocation or a
    /// second dissipative model. This frequency-domain result advances no clock
    /// and does not convert reactive stored energy into heat.
    pub fn drive(
        &self,
        incoming_velocity_m_s: [C64; 3],
    ) -> Result<SupportedElasticResponse, ElasticSupportError> {
        let reflection = self.reflection.drive(incoming_velocity_m_s)?;
        let storage = support_storage(
            [self.tangential, self.tangential, self.normal],
            reflection.boundary_velocity_m_s,
            self.omega_rad_s,
        )?;
        Ok(SupportedElasticResponse {
            reflection,
            support_displacement_m: storage.displacement,
            mean_kinetic_energy_j_m2: storage.kinetic,
            mean_elastic_energy_j_m2: storage.elastic,
            mean_stored_energy_j_m2: storage.total,
        })
    }
}

struct SupportStorage {
    displacement: [C64; 3],
    kinetic: [f64; 3],
    elastic: [f64; 3],
    total: f64,
}

// Both surface owners use the original pin potentials in their physical axes.
fn support_storage(
    pins: [WallPin; 3],
    velocity: [C64; 3],
    omega: f64,
) -> Result<SupportStorage, ElasticSupportError> {
    let mut storage = SupportStorage {
        displacement: [C64::ZERO; 3],
        kinetic: [0.0; 3],
        elastic: [0.0; 3],
        total: 0.0,
    };
    let mut total = Dd::ZERO;
    for axis in 0..3 {
        storage.displacement[axis] = displacement(velocity[axis], omega)?;
        storage.kinetic[axis] = mean_energy(pins[axis].surface_density, velocity[axis].abs())?;
        storage.elastic[axis] = mean_energy(
            pins[axis].stiffness_per_area,
            storage.displacement[axis].abs(),
        )?;
        total = total + Dd::from_f64(storage.kinetic[axis]) + Dd::from_f64(storage.elastic[axis]);
    }
    storage.total = total.to_f64();
    if !storage.total.is_finite() {
        return Err(ElasticSupportError::NumericalRange);
    }
    Ok(storage)
}

fn displacement(velocity: C64, omega: f64) -> Result<C64, ElasticSupportError> {
    let value = C64::new(-velocity.im / omega, velocity.re / omega);
    if !value.re.is_finite()
        || !value.im.is_finite()
        || (velocity != C64::ZERO && value == C64::ZERO)
    {
        return Err(ElasticSupportError::NumericalRange);
    }
    Ok(value)
}

fn mean_energy(coefficient: f64, amplitude: f64) -> Result<f64, ElasticSupportError> {
    let root = det::sqrt(coefficient) * (0.5 * amplitude);
    let energy = root * root;
    if !energy.is_finite() || (amplitude != 0.0 && energy == 0.0) {
        return Err(ElasticSupportError::NumericalRange);
    }
    Ok(energy)
}
