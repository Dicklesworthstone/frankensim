//! Intrinsic solid absorption from the original causal bulk/shear Maxwell laws.
//!
//! A laterally infinite homogeneous isotropic slab at NORMAL incidence has
//! only longitudinal displacement. Its modulus is M*=K*+4G*/3, including
//! both constitutive loss channels. The exp(+i omega t) convention has positive
//! loss moduli and outgoing k with nonpositive imaginary part.
//! See Meza-Fajardo & Lai (2007), GJI 171, 1247–1257, equations (5), (10)–(11):
//! <https://doi.org/10.1111/j.1365-246X.2007.03590.x>.
//!
//! The original scalar Maxwell owners evaluate K* and G*: no fitted attenuation,
//! constant loss factor, static-to-dynamic modulus conversion or weak-loss
//! expansion is introduced. Bounded inward waves retain coherent two-face
//! driving and thick slabs. Positive bulk/shear strain-work integrals remain
//! independent of reflected/transmitted powers, including below roundoff loss.
//! This frozen, small-strain, steady harmonic response does not evolve a
//! thermal state, infer measured source authority, resolve oblique P/S waves,
//! move a porous frame or model a finite aperture.

use crate::MaterialError;
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use crate::visco::GeneralizedMaxwell;
use fs_blake3::ContentHash;
use fs_math::c64::C64;

mod kernel;

/// Coherent bonded solid stacks with the original loss law in every layer.
pub mod layers;

/// Bounded total number of bulk plus shear Maxwell branches at preparation.
pub const MAX_VISCOELASTIC_SLAB_BRANCHES: usize = 64;

/// Explicit original material definition; the scalar laws mean K and G, not E.
/// Density and support must describe this same frozen material state.
#[derive(Clone, Copy, Debug)]
pub struct ViscoelasticSlabSource<'a> {
    /// Original bulk-modulus relaxation law [Pa].
    pub bulk: &'a GeneralizedMaxwell,
    /// Original shear-modulus relaxation law [Pa].
    pub shear: &'a GeneralizedMaxwell,
    /// Original mass per solid volume [kg/m³].
    pub density_kg_m3: f64,
    /// Closed supported angular-frequency interval [rad/s].
    pub frequency_band_rad_s: (f64, f64),
    /// Positive peak absolute longitudinal strain ceiling of the original law.
    pub maximum_absolute_strain: f64,
    /// Nonzero identity of the authored material state, not inferred evidence.
    pub source_identity: ContentHash,
}

/// Typed refusal with no partial driven field or power report.
#[derive(Debug)]
pub enum ViscoelasticSlabError {
    /// Invalid source, geometry, wave, drive or unrepresentable numerical result.
    Invalid(&'static str),
    /// The native exterior fluid source refused its own state.
    Material(MaterialError),
    /// Original material's frequency or conservative peak-strain support failed.
    OutsideSourceSupport,
    /// Authored work allowance or hard branch limit was exceeded.
    Budget,
    /// Cancellation was requested at a bounded work boundary.
    Cancelled,
}
impl core::fmt::Display for ViscoelasticSlabError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "viscoelastic slab refused: {self:?}")
    }
}
impl core::error::Error for ViscoelasticSlabError {}
impl From<MaterialError> for ViscoelasticSlabError {
    fn from(value: MaterialError) -> Self {
        Self::Material(value)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ViscoelasticSlabError> {
    if cancelled() {
        Err(ViscoelasticSlabError::Cancelled)
    } else {
        Ok(())
    }
}

fn source_moduli(
    source: ViscoelasticSlabSource<'_>,
    angular_frequency_rad_s: f64,
    maximum_branches: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[C64; 2], ViscoelasticSlabError> {
    let count = source
        .bulk
        .terms
        .len()
        .checked_add(source.shear.terms.len());
    if maximum_branches > MAX_VISCOELASTIC_SLAB_BRANCHES
        || count.is_none_or(|n| n > maximum_branches)
    {
        return Err(ViscoelasticSlabError::Budget);
    }
    let (lo, hi) = source.frequency_band_rad_s;
    if !lo.is_finite()
        || !hi.is_finite()
        || lo < 0.0
        || hi < lo
        || !source.maximum_absolute_strain.is_finite()
        || source.maximum_absolute_strain <= 0.0
        || source.source_identity == ContentHash([0; 32])
    {
        return Err(ViscoelasticSlabError::Invalid(
            "material identity and support must be explicit",
        ));
    }
    if !angular_frequency_rad_s.is_finite()
        || angular_frequency_rad_s <= 0.0
        || angular_frequency_rad_s < lo
        || angular_frequency_rad_s > hi
    {
        return Err(ViscoelasticSlabError::OutsideSourceSupport);
    }
    for law in [source.bulk, source.shear] {
        if !law.e_inf.is_finite() || law.e_inf <= 0.0 {
            return Err(ViscoelasticSlabError::Invalid(
                "original equilibrium modulus must be positive",
            ));
        }
        for &(modulus, time) in &law.terms {
            poll(cancelled)?;
            if !modulus.is_finite() || modulus < 0.0 || !time.is_finite() || time <= 0.0 {
                return Err(ViscoelasticSlabError::Invalid(
                    "original Maxwell branch is invalid",
                ));
            }
        }
    }
    Ok([source.bulk, source.shear].map(|law| {
        let (storage, loss) = law.modulus(angular_frequency_rad_s);
        C64::new(storage, loss)
    }))
}

/// One source-bound slab at one frequency, between two lossless exterior fluids.
/// Preparation queries each distinct exterior source once. Driving and point
/// observation allocate nothing and perform constant work with no source query.
pub struct ViscoelasticSlab<'a> {
    source: ViscoelasticSlabSource<'a>,
    exterior_sources: [&'a dyn FluidAcoustics; 2],
    exterior_states: [FluidAcousticState; 2],
    wave: kernel::Slab,
}

/// Independent mean powers over the full authored slab area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscoelasticSlabPower {
    /// Incoming cycle-mean watts at the left and right faces, respectively.
    pub incident_power_w: [f64; 2],
    /// Outgoing cycle-mean watts at the left and right faces, respectively.
    pub outgoing_power_w: [f64; 2],
    /// Spatially integrated original bulk relaxation heat [W].
    pub bulk_absorbed_power_w: f64,
    /// Spatially integrated original shear relaxation heat [W].
    pub shear_absorbed_power_w: f64,
    /// Incoming minus outgoing minus independently integrated heat [W].
    pub balance_residual_w: f64,
    /// Conservative all-depth peak-strain bound used for source admission.
    pub maximum_absolute_strain_bound: f64,
}

/// Physical peak phasors and local cycle-mean heat inside the solid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscoelasticSlabField {
    /// Compressive normal traction, -sigma_zz [Pa]; not mean hydrostatic pressure.
    pub compressive_traction_pa: C64,
    /// Normal particle velocity [m/s], positive from left to right.
    pub normal_velocity_m_s: C64,
    /// Longitudinal strain epsilon_zz; both transverse normal strains are zero.
    pub longitudinal_strain: C64,
    /// Original bulk relaxation heat density [W/m³].
    pub bulk_heat_density_w_m3: f64,
    /// Original shear relaxation heat density [W/m³].
    pub shear_heat_density_w_m3: f64,
}

/// Reusable coherent driven field, retaining the same prepared material owner.
pub struct DrivenViscoelasticSlab<'s, 'a> {
    slab: &'s ViscoelasticSlab<'a>,
    drive: kernel::Drive,
}

impl<'a> ViscoelasticSlab<'a> {
    /// Prepare a normal-incidence slab without copying either original law.
    /// The explicit allowance bounds the combined branch count, up to 64.
    /// Exterior sources provide lossless acoustic snapshots; no thermal or
    /// molecular exterior loss is inferred. Zero thickness is an exact interface
    /// identity but still validates source support and branch parameters.
    pub fn try_new(
        source: ViscoelasticSlabSource<'a>,
        exterior_sources: [&'a dyn FluidAcoustics; 2],
        angular_frequency_rad_s: f64,
        thickness_m: f64,
        maximum_branches: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ViscoelasticSlabError> {
        poll(cancelled)?;
        let moduli = source_moduli(source, angular_frequency_rad_s, maximum_branches, cancelled)?;
        let left = exterior_sources[0].acoustic_state()?;
        let right = if core::ptr::eq(exterior_sources[0], exterior_sources[1]) {
            left
        } else {
            exterior_sources[1].acoustic_state()?
        };
        poll(cancelled)?;
        let wave = kernel::Slab::new(
            source.density_kg_m3,
            moduli,
            angular_frequency_rad_s,
            thickness_m,
            [
                left.characteristic_impedance_pa_s_m(),
                right.characteristic_impedance_pa_s_m(),
            ],
        )
        .map_err(ViscoelasticSlabError::Invalid)?;
        poll(cancelled)?;
        Ok(Self {
            source,
            exterior_sources,
            exterior_states: [left, right],
            wave,
        })
    }

    /// Original source references, domain, density and material-state identity.
    pub const fn source(&self) -> ViscoelasticSlabSource<'a> {
        self.source
    }
    /// Original lossless exterior sources, without re-evaluation.
    pub const fn exterior_sources(&self) -> [&'a dyn FluidAcoustics; 2] {
        self.exterior_sources
    }
    /// Coherent snapshots used by the two physical boundary conditions.
    pub const fn exterior_states(&self) -> [FluidAcousticState; 2] {
        self.exterior_states
    }
    /// Prepared original bulk and shear complex moduli, in that order [Pa].
    pub const fn moduli_pa(&self) -> [C64; 2] {
        self.wave.moduli
    }
    /// Passive forward longitudinal wavenumber [rad/m], in exp(+i omega t).
    pub const fn wavenumber_rad_m(&self) -> C64 {
        self.wave.k
    }
    /// Physical solid thickness [m].
    pub const fn thickness_m(&self) -> f64 {
        self.wave.length
    }
    /// Row-major pressure scattering matrix: outgoing [left,right] from incoming
    /// [left,right]. Different exterior impedances require their own power weights.
    pub const fn pressure_scattering(&self) -> [C64; 4] {
        self.wave.scattering
    }

    /// Apply coherent peak incident pressures at the two exterior faces.
    /// Both illuminate the SAME complete positive area. The conservative bound
    /// (|A|+|B|)/|M*| admits all interior peak strains; it may refuse a cancellation-
    /// dominated field that a tighter spatial maximum would admit. No result is
    /// published after source-domain, balance or cancellation refusal.
    pub fn drive<'s>(
        &'s self,
        incident_peak_pressure_pa: [C64; 2],
        area_m2: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DrivenViscoelasticSlab<'s, 'a>, ViscoelasticSlabError> {
        poll(cancelled)?;
        let drive = self
            .wave
            .drive(incident_peak_pressure_pa, area_m2)
            .map_err(ViscoelasticSlabError::Invalid)?;
        if drive.strain_bound > self.source.maximum_absolute_strain {
            return Err(ViscoelasticSlabError::OutsideSourceSupport);
        }
        poll(cancelled)?;
        Ok(DrivenViscoelasticSlab { slab: self, drive })
    }
}

impl DrivenViscoelasticSlab<'_, '_> {
    /// Independent acoustic and irreversible-heat power accounts.
    pub const fn power(&self) -> ViscoelasticSlabPower {
        ViscoelasticSlabPower {
            incident_power_w: self.drive.incident_power,
            outgoing_power_w: self.drive.outgoing_power,
            bulk_absorbed_power_w: self.drive.heat[0],
            shear_absorbed_power_w: self.drive.heat[1],
            balance_residual_w: self.drive.residual,
            maximum_absolute_strain_bound: self.drive.strain_bound,
        }
    }
    /// Outgoing peak pressures at the original left and right reference faces.
    pub const fn outgoing_pressure_pa(&self) -> [C64; 2] {
        self.drive.outgoing
    }
    /// Observe original solid fields and separate intrinsic heat mechanisms.
    /// Depth includes both faces. A zero-thickness identity has no interior and
    /// refuses this query. No spatial quadrature or material query is performed.
    pub fn field(
        &self,
        depth_m: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ViscoelasticSlabField, ViscoelasticSlabError> {
        poll(cancelled)?;
        let field = self
            .slab
            .wave
            .field(self.drive.inward, self.drive.parity, depth_m)
            .map_err(ViscoelasticSlabError::Invalid)?;
        poll(cancelled)?;
        Ok(ViscoelasticSlabField {
            compressive_traction_pa: field.pressure,
            normal_velocity_m_s: field.velocity,
            longitudinal_strain: field.strain,
            bulk_heat_density_w_m3: field.heat_density[0],
            shear_heat_density_w_m3: field.heat_density[1],
        })
    }
}
