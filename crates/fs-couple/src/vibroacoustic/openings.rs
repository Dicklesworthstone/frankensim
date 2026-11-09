//! Joint radiation of the moving structure and actual cavity openings.
//!
//! A rigid-wall pressure basis retains its uniform compression mode. Each
//! opening retains physical end-flow unknowns, not a second Helmholtz mode.
//! The exterior sees [original structural velocities, opening volume flows]
//! through ONE complete sampled impedance. The reciprocal cavity map supplies
//! both opening pressure and outward flow in the acoustic weak equations.

use super::{AcousticMedium, C64, ModalRadiationLoading, VibroError, VibroacousticModel};
use fs_exec::Cx;
use fs_material::gas::GasState;

mod physics;
pub use physics::CavityOpeningResponse;
use physics::OpeningFrequency;
pub(super) use physics::{flow_delta, relative_flow_delta};
mod distributed;
#[cfg(test)]
pub(super) use distributed::DistributedOpeningFrequency;
pub use distributed::{DistributedCavityOpeningRadiation, DistributedOpeningControls};
pub(super) use distributed::{
    DuctRow, Material as DuctMaterial, observe_duct, prepare_cylindrical_rows,
};
#[cfg(test)]
pub(super) use physics::OpeningFrequency as CompactOpeningFrequency;

/// Same coupled solve; the finite duct retains both of its physical end flows.
pub(super) enum CoupledOpeningFrequency<'a> {
    Compact(OpeningFrequency<'a>),
    Distributed(distributed::DistributedOpeningFrequency<'a>),
}
impl CoupledOpeningFrequency<'_> {
    pub(super) fn unknown_count(&self) -> usize {
        match self {
            Self::Compact(law) => law.count(),
            Self::Distributed(law) => law.unknown_count(),
        }
    }
    pub(super) fn add_operator(
        &self,
        omega: f64,
        bulk: f64,
        lambdas: &[f64],
        matrix: &mut [C64],
        keep_s: usize,
        keep_a: usize,
        n: usize,
    ) -> Result<(), VibroError> {
        match self {
            Self::Compact(law) => law.add_operator(omega, bulk, lambdas, matrix, keep_s, keep_a, n),
            Self::Distributed(law) => {
                law.add_operator(omega, bulk, lambdas, matrix, keep_s, keep_a, n)
            }
        }
    }
    pub(super) fn response(
        &self,
        omega: f64,
        displacement: &[C64],
        pressure: &[C64],
        flow: &[C64],
    ) -> Result<CavityOpeningResponse, VibroError> {
        match self {
            Self::Compact(law) => law.response(omega, displacement, pressure, flow),
            Self::Distributed(law) => law.response(omega, displacement, pressure, flow),
        }
    }
}

/// Maximum number of explicitly retained compact or distributed openings.
pub const MAX_CAVITY_OPENINGS: usize = 64;

/// A straight cylindrical neck inside the exterior radiation boundary.
/// Compact and distributed owners impose their own axial propagation domains.
///
/// The exterior matrix contains the OUTER end, including radiation mass and
/// loss, exactly once. The inner entrance either belongs to the resolved cavity
/// pressure field or to an explicitly selected compact fs-phs end correction.
#[derive(Clone, Copy, Debug)]
pub struct CylindricalCavityNeck {
    pub radius_m: f64,
    pub length_m: f64,
    /// None: the retained cavity pressure/flow field owns entrance inertia.
    /// Some: add the native compact correction for this flange. This is
    /// admitted only with ONE uniform cavity mode, which cannot resolve the
    /// entrance field; combining it with higher cavity modes would count
    /// resolved inner-end inertia again.
    pub inner_end: Option<fs_phs::MouthFlange>,
}

#[derive(Clone, Debug)]
struct NeckRow {
    omega: f64,
    impedance: Vec<C64>,
}

/// Complete mutual exterior loading and the original internal neck laws.
///
/// Exterior port order is [structural modal velocities, outward volume flows].
/// A BEM/Rayleigh source row for opening j is its outward velocity per unit
/// volume flow, with integral 1 over that mouth; it is NOT unit velocity.
/// The corresponding effort is pressure [Pa]. Structural source rows retain
/// their original mass normalization. Consequently the matrix includes every
/// wall/mouth and mouth/mouth reaction in their physical port units.
///
/// `pressure_projection[q][j] = integral psi_q chi_j dA` uses that same
/// unit-volume-flow profile on the CAVITY side. It maps opening flows to
/// outward modal flow; its TRANSPOSE maps modal pressure to opening pressure.
/// For uniform flow it is the area mean of psi_q, not its area integral.
/// The caller owns matching the actual source geometry, basis and normals,
/// and keeping the opening distinct from moving-wall/source/boundary patches.
///
/// The native Bessel duct owner supplies the complex series impedance
/// `-i k Zc L` at every exact exterior frequency. An explicitly unresolved
/// inner entrance may add `-i omega rho delta_inner / area`; a resolved cavity
/// must use no added end correction. This retains viscous diffusion
/// in both neck inertia and heat, without a resonance-frozen resistance.
/// One flow describes the whole SHORT neck: its compressibility, thermal
/// compression loss, transverse modes, edge-flow loss, mean flow and nonlinear
/// jet separation are not resolved. The explicit small-phase admission is an
/// approximation domain, not an error certificate.
#[derive(Clone, Debug)]
pub struct CavityOpeningRadiation {
    structural_modes: usize,
    pressure_projection: Vec<Vec<f64>>,
    necks: Vec<CylindricalCavityNeck>,
    gas: GasState,
    maximum_wavenumber_length: f64,
    exterior: ModalRadiationLoading,
    rows: Vec<NeckRow>,
}

fn invalid(what: &'static str) -> VibroError {
    VibroError::BadParameter { what }
}

fn poll(cx: &Cx<'_>) -> Result<(), VibroError> {
    cx.checkpoint().map_err(|_| VibroError::Cancelled)
}

impl CavityOpeningRadiation {
    /// Prepare neck response on the COMPLETE original joint exterior grid.
    ///
    /// Supply one already admitted `ModalRadiationLoading` with
    /// `structural_modes + necks.len()` ports; `from_bem_bake` retains the
    /// producer's complete grid and declaration. The additional neck storage
    /// is linear in that existing table's frequencies and opening count.
    /// No source direction or mutual entry is discarded.
    ///
    /// At every frequency and neck,
    /// `|k| max(radius, length + delta_inner) <= maximum_wavenumber_length <= 0.25`.
    /// A tighter caller limit is allowed. Longer ducts require their distributed
    /// pressure/flow owner; this constructor never extrapolates a compact neck.
    ///
    /// # Errors
    /// Shape, source-gas, geometry or compact-domain refusal; an unsupported
    /// native duct/end response; allocation failure or cancellation. Admission
    /// is complete before an immutable loading is returned.
    pub fn from_cylindrical_necks(
        structural_modes: usize,
        pressure_projection: Vec<Vec<f64>>,
        necks: &[CylindricalCavityNeck],
        gas: &GasState,
        exterior: ModalRadiationLoading,
        maximum_wavenumber_length: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, VibroError> {
        poll(cx)?;
        let count = necks.len();
        if structural_modes == 0
            || count == 0
            || count > MAX_CAVITY_OPENINGS
            || structural_modes.checked_add(count) != Some(exterior.mode_count())
            || pressure_projection.is_empty()
            || pressure_projection.iter().any(|row| row.len() != count)
        {
            return Err(VibroError::ShapeMismatch {
                what: "joint loading needs original structure, cavity map and 1..=64 opening ports",
            });
        }
        if pressure_projection.iter().flatten().any(|v| !v.is_finite())
            || !maximum_wavenumber_length.is_finite()
            || maximum_wavenumber_length <= 0.0
            || maximum_wavenumber_length > 0.25
        {
            return Err(invalid(
                "finite opening map and compact phase limit in (0, 0.25] required",
            ));
        }
        AcousticMedium::try_from_gas(gas)?;
        if pressure_projection.len() != 1 && necks.iter().any(|neck| neck.inner_end.is_some()) {
            return Err(invalid(
                "compact inner-end correction is exclusive to one uniform cavity mode",
            ));
        }
        let mut inner_lengths = Vec::with_capacity(count);
        for neck in necks {
            poll(cx)?;
            if !neck.radius_m.is_finite()
                || neck.radius_m <= 0.0
                || !neck.length_m.is_finite()
                || neck.length_m < 0.0
            {
                return Err(invalid(
                    "cavity neck needs positive radius and nonnegative finite length",
                ));
            }
            let inner = match neck.inner_end {
                Some(flange) => fs_phs::neck_end_correction(neck.radius_m, flange)
                    .map_err(|_| invalid("native inner neck end correction refused"))?,
                None => 0.0,
            };
            inner_lengths.push(inner);
        }
        let mut rows = Vec::new();
        rows.try_reserve_exact(exterior.angular_frequencies_rad_s().len())
            .map_err(|_| invalid("cavity neck frequency table allocation failed"))?;
        for omega in exterior.angular_frequencies_rad_s() {
            poll(cx)?;
            let mut impedance = Vec::with_capacity(count);
            for (neck, &inner) in necks.iter().zip(&inner_lengths) {
                poll(cx)?;
                let wave =
                    fs_duct::segment_wave(gas, neck.radius_m, omega, fs_duct::LossModel::Bessel)
                        .map_err(|_| invalid("native Bessel neck response refused"))?;
                let phase = wave.wavenumber.abs() * neck.radius_m.max(neck.length_m + inner);
                if !phase.is_finite() || phase > maximum_wavenumber_length {
                    return Err(invalid(
                        "cavity neck exceeds its declared compact wavenumber-length domain",
                    ));
                }
                let area = core::f64::consts::PI * neck.radius_m * neck.radius_m;
                let z = (C64::new(0.0, -1.0) * wave.wavenumber * wave.characteristic_impedance)
                    .scale(neck.length_m)
                    + C64::new(0.0, -omega * (gas.density * inner / area));
                if !z.re.is_finite() || !z.im.is_finite() || z.re < 0.0 || z.im > 0.0 {
                    return Err(invalid("native neck impedance is not finite and passive"));
                }
                impedance.push(z);
            }
            rows.push(NeckRow { omega, impedance });
        }
        poll(cx)?;
        Ok(Self {
            structural_modes,
            pressure_projection,
            necks: necks.to_vec(),
            gas: *gas,
            maximum_wavenumber_length,
            exterior,
            rows,
        })
    }

    pub fn opening_count(&self) -> usize {
        self.necks.len()
    }
    pub fn cavity_mode_count(&self) -> usize {
        self.pressure_projection.len()
    }
    pub const fn structural_mode_count(&self) -> usize {
        self.structural_modes
    }
    pub fn gas(&self) -> &GasState {
        &self.gas
    }
    pub fn necks(&self) -> &[CylindricalCavityNeck] {
        &self.necks
    }
    pub fn pressure_projection(&self) -> &[Vec<f64>] {
        &self.pressure_projection
    }
    pub const fn maximum_wavenumber_length(&self) -> f64 {
        self.maximum_wavenumber_length
    }
    pub fn exterior_loading(&self) -> &ModalRadiationLoading {
        &self.exterior
    }

    /// Original internal-only impedance [Pa s/m^3], in opening order.
    /// An absent frequency refuses; no interpolation or frozen neck law.
    pub fn neck_impedances(&self, omega: f64) -> Result<&[C64], VibroError> {
        Ok(self.at_frequency(omega)?.neck)
    }

    pub(super) fn at_frequency(&self, omega: f64) -> Result<OpeningFrequency<'_>, VibroError> {
        let exterior = self.exterior.at_omega(omega)?;
        let row = self
            .rows
            .binary_search_by(|row| row.omega.total_cmp(&omega))
            .map_err(|_| invalid("opening neck has no original frequency sample"))?;
        Ok(OpeningFrequency {
            structural_modes: self.structural_modes,
            pressure_projection: &self.pressure_projection,
            exterior,
            neck: &self.rows[row].impedance,
        })
    }
}

impl VibroacousticModel {
    /// Add actual radiating openings to a RIGID-WALL cavity pressure basis.
    ///
    /// The first acoustic mode must be the zero-frequency uniform compression
    /// mode. In particular, do not use `helmholtz_resonator_mode`: that already
    /// includes a neck resonance, which would be counted again. The ordinary
    /// rectangular cavity constructor supplies the required sealed basis.
    ///
    /// This joint exterior replaces a structural-only exterior attachment;
    /// both cannot coexist. Existing structural material, cavity molecular loss
    /// and disjoint locally reacting walls remain in the same solve. All
    /// openings stay present during modal truncation, with the same original
    /// pressure map and full exterior stride. The unloaded eigenfrequency
    /// diagnostic remains the sealed-cavity reference; use the FRF for vented
    /// resonances.
    pub fn with_cavity_openings(
        mut self,
        openings: CavityOpeningRadiation,
    ) -> Result<Self, VibroError> {
        if openings.structural_modes != self.omegas_s.len()
            || openings.cavity_mode_count() != self.omegas_a.len()
        {
            return Err(VibroError::ShapeMismatch {
                what: "cavity openings must retain the original structural and pressure bases",
            });
        }
        if self.z_modal.is_some() || self.openings.is_some() || self.distributed_openings.is_some()
        {
            return Err(invalid(
                "one complete joint exterior replaces separate structural radiation",
            ));
        }
        if self.omegas_a[0] != 0.0 {
            return Err(invalid(
                "radiating openings require the rigid-wall cavity's uniform zero mode",
            ));
        }
        if self
            .internal_ducts
            .as_ref()
            .is_some_and(|ducts| ducts.gas() != openings.gas())
        {
            return Err(invalid(
                "internal and exterior opening ducts must retain the same original gas",
            ));
        }
        let medium = AcousticMedium::try_from_gas(&openings.gas)?;
        if medium.rho0 != self.medium_a.rho0
            || medium.c0 != self.medium_a.c0
            || self
                .cavity_relaxation
                .as_ref()
                .is_some_and(|law| law.gas() != &openings.gas)
        {
            return Err(invalid(
                "internal cavity necks and pressure basis must retain the same gas",
            ));
        }
        self.openings = Some(openings);
        Ok(self)
    }

    pub fn cavity_openings(&self) -> Option<&CavityOpeningRadiation> {
        self.openings.as_ref()
    }
}

#[cfg(test)]
mod tests;
