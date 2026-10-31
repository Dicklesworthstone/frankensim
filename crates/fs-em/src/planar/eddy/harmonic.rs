//! Bias-dependent winding impedance from resolved magnetoquasistatic fields.
//!
//! Peak phasors use `exp(-i omega t)`, matching `fs_phs::harmonic`. The original
//! smooth DC operating point supplies its complete magnetic Hessian K_t; the
//! original floating-conductor projection supplies C_sigma. Each unit winding
//! solves `(K_t - i omega C_sigma) a = f`. No scalar permeability replacement,
//! fitted core loss or dense whole-field factorization is introduced.
//!
//! The same signed P1 winding loads recover linkage and terminal impedance
//! `Z = R_copper - i omega (F^T a + L_end)`. Stranded winding cells and resolved
//! eddy cells must be disjoint. Winding resistance comes from original material
//! sources and represented volume. Explicit [`SolidWindingRegion`] bindings
//! additionally resolve skin/proximity redistribution in driven solid sections,
//! using their prescribed integrated current and the same conductivity operator.
//! Stranded profiles stay fixed; their internal redistribution is not inferred.
//!
//! Geometry, temperature, boundary trace, remanence and DC bias stay fixed.
//! These are infinitesimal perturbations, with incremental stored energy and
//! cycle-mean AC heat; DC copper heat is additional. Hysteresis, finite-amplitude
//! saturation cycles, displacement current and a causal time-domain reduction
//! of the frequency-dependent terminal matrix remain outside this owner.

/// Native circuit-network composition of the field-derived terminal matrix.
pub mod circuit;
mod response;
/// Resolved solid winding currents, skin/proximity redistribution and heat.
pub mod solid;
mod solve;
pub use solid::{HarmonicSolidCell, SolidWindingCurrent, SolidWindingRegion};

use super::{Conduction, EddyConductor, EddyError, bind, load, poll, positive, regions};
use crate::{
    planar::{
        MagneticError, PlanarMagnetostatics,
        nonlinear::{IsotropicMagneticMaterial, NonlinearMagneticError},
    },
    winding::{
        incremental::PlanarIncrementalWindingBasis,
        ohmic::{OhmicWindingError, WindingConductor, bind_sources},
    },
};
use fs_material::field_storage::{admit_fixed_magnetic_refresh, linear_magnetic_permeability};
use fs_math::c64::C64;
use fs_solver::SolveReport;
use fs_sparse::Csr;

/// Numerical work and physical admission limits, held across a frequency sweep.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicWindingControls {
    pub relative_tolerance: f64,
    /// Actual weak-equation target is max(this [A m], relative * norm(rhs)).
    pub absolute_residual_a_m: f64,
    pub restart: usize,
    pub max_cycles: usize,
    /// Existing FGMRES workspace, including restart history and Hessenberg data.
    pub max_krylov_scalars: usize,
    /// Conservative additional field/source/response numerical payload. Existing
    /// magnetic assembly limits and the separate Krylov budget still apply.
    pub max_field_scalars: usize,
    /// Complex reciprocity and source/field power-pairing allowance [ohm].
    pub impedance_tolerance_ohm: f64,
    /// Cycle-average real and reactive power-accounting allowance [W or var].
    pub power_tolerance_w: f64,
    pub current_balance_tolerance_a: f64,
    /// Caller-admitted maximum cell diameter / smallest local skin depth,
    /// using BOTH differential and transverse permeability at the DC bias.
    /// This is a spatial admission limit, not an error certificate.
    pub max_cell_skin_ratio: f64,
}

#[derive(Debug)]
pub enum HarmonicWindingError {
    Invalid(&'static str),
    Magnetic(NonlinearMagneticError),
    Eddy(EddyError),
    Winding(OhmicWindingError),
    Nonconvergence(SolveReport),
    Balance(&'static str),
    Cancelled,
}
impl core::fmt::Display for HarmonicWindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "harmonic winding field: {self:?}")
    }
}
impl std::error::Error for HarmonicWindingError {}
impl From<EddyError> for HarmonicWindingError {
    fn from(error: EddyError) -> Self {
        match error {
            EddyError::Cancelled | EddyError::Magnetic(MagneticError::Cancelled) => Self::Cancelled,
            other => Self::Eddy(other),
        }
    }
}
impl From<MagneticError> for HarmonicWindingError {
    fn from(error: MagneticError) -> Self {
        Self::from(EddyError::Magnetic(error))
    }
}
impl From<NonlinearMagneticError> for HarmonicWindingError {
    fn from(error: NonlinearMagneticError) -> Self {
        match error {
            NonlinearMagneticError::Magnetic(MagneticError::Cancelled) => Self::Cancelled,
            other => Self::Magnetic(other),
        }
    }
}
impl From<OhmicWindingError> for HarmonicWindingError {
    fn from(error: OhmicWindingError) -> Self {
        match error {
            OhmicWindingError::Cancelled => Self::Cancelled,
            other => Self::Winding(other),
        }
    }
}
type Result<T> = core::result::Result<T, HarmonicWindingError>;

fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(HarmonicWindingError::Invalid(
            "unrepresentable harmonic arithmetic",
        ))
    }
}
fn phasor(value: C64) -> Result<C64> {
    finite(value.re)?;
    finite(value.im)?;
    Ok(value)
}

/// Prepared sparse differential field and original electrical sources. Multiple
/// frequencies reuse source admission, the DC Hessian and signed winding loads.
pub struct PlanarHarmonicWindings<'b, 'm> {
    basis: &'b PlanarIncrementalWindingBasis<'m>,
    tangent: Csr,
    free: Vec<usize>,
    loads: Vec<Vec<f64>>,
    conduction: Conduction,
    eddy_conductors: Vec<EddyConductor>,
    winding_conductors: Vec<WindingConductor>,
    conductor_windings: Vec<usize>,
    resistance_ohm: Vec<f64>,
    solid: solid::SolidWindings,
    skin_ratio_per_sqrt_rad_s: Vec<f64>,
    controls: HarmonicWindingControls,
    restart: usize,
}

/// One retained derivative field per terminal; all prescribed perturbation
/// potentials are zero even when the original DC boundary trace is nonzero.
#[derive(Clone, Debug)]
pub struct HarmonicWindingUnitResponse {
    potential_per_current_h_m: Vec<C64>,
    report: SolveReport,
}
impl HarmonicWindingUnitResponse {
    pub fn potential_per_current_h_m(&self) -> &[C64] {
        &self.potential_per_current_h_m
    }
    pub fn report(&self) -> &SolveReport {
        &self.report
    }
}

/// Frequency-specific terminal matrix and resolved unit fields. Repeated drives
/// reuse these fields; only the at-most-16-winding terminal matrix is dense.
pub struct HarmonicWindingFrequency<'p, 'b, 'm> {
    prepared: &'p PlanarHarmonicWindings<'b, 'm>,
    omega_rad_s: f64,
    impedance_ohm: Vec<C64>,
    effective_inductance_h: Vec<f64>,
    unit_responses: Vec<HarmonicWindingUnitResponse>,
    maximum_reciprocity_defect_ohm: f64,
    maximum_power_pairing_defect_ohm: f64,
    maximum_cell_skin_ratio: f64,
}
impl HarmonicWindingFrequency<'_, '_, '_> {
    pub fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    /// Row-major Z in original winding order; inductive imaginary parts are
    /// NEGATIVE under the shared exp(-i omega t) convention.
    pub fn impedance_ohm(&self) -> &[C64] {
        &self.impedance_ohm
    }
    /// Row-major -Im(Z)/omega, including explicit end-storage diagonals.
    pub fn effective_inductance_h(&self) -> &[f64] {
        &self.effective_inductance_h
    }
    pub fn unit_responses(&self) -> &[HarmonicWindingUnitResponse] {
        &self.unit_responses
    }
    pub fn maximum_reciprocity_defect_ohm(&self) -> f64 {
        self.maximum_reciprocity_defect_ohm
    }
    pub fn maximum_power_pairing_defect_ohm(&self) -> f64 {
        self.maximum_power_pairing_defect_ohm
    }
    pub fn maximum_cell_skin_ratio(&self) -> f64 {
        self.maximum_cell_skin_ratio
    }
}

#[derive(Clone, Debug)]
pub struct HarmonicEddyCell {
    pub source: EddyConductor,
    pub volume_m3: f64,
    /// P1 peak values in original local vertex order, including the floating
    /// scalar-potential gradient. J = sigma E at this same material state.
    pub electric_field_v_m: [C64; 3],
    pub current_density_a_m2: [C64; 3],
    pub mean_power_w: f64,
}
#[derive(Clone, Debug)]
pub struct HarmonicWindingCell {
    pub source: WindingConductor,
    pub winding: usize,
    pub volume_m3: f64,
    pub current_density_a_m2: C64,
    pub electric_field_v_m: C64,
    pub mean_power_w: f64,
}

/// Actual phase-sensitive superposition and independently reconstructed spatial
/// loss. All phasors are peaks; energies/powers are means over one AC cycle.
#[derive(Clone, Debug)]
pub struct HarmonicWindingResponse {
    pub omega_rad_s: f64,
    pub currents_a: Vec<C64>,
    pub terminal_voltages_v: Vec<C64>,
    pub potential_wb_m: Vec<C64>,
    pub flux_density_t: Vec<[C64; 2]>,
    pub eddy_cells: Vec<HarmonicEddyCell>,
    pub winding_cells: Vec<HarmonicWindingCell>,
    /// Actual total current and heat in driven solid sections, deposited once.
    pub solid_cells: Vec<HarmonicSolidCell>,
    /// Signed current targets and measured integrals in driven solid sections.
    pub solid_region_currents: Vec<SolidWindingCurrent>,
    /// Unbound floating regions only; driven sections have separate receipts.
    pub region_currents_a: Vec<(u32, C64)>,
    /// Local quadratic perturbation storage, including declared end leakage.
    pub mean_magnetic_energy_j: f64,
    /// Unbound floating conductor cells only.
    pub mean_eddy_power_w: f64,
    /// Fixed-profile stranded conductor cells only.
    pub mean_winding_power_w: f64,
    /// Total heat in driven solid cells, including all local cross terms.
    pub mean_solid_power_w: f64,
    pub mean_supplied_power_w: f64,
    /// Im(sum(conj(I) V))/2; negative for this convention's inductive load.
    pub reactive_power_var: f64,
    /// Spatial Joule power minus terminal real power.
    pub power_residual_w: f64,
    /// Terminal reactive power + 2 omega * mean magnetic perturbation storage.
    pub reactive_power_residual_var: f64,
}
impl HarmonicWindingResponse {
    /// Physical cycle-mean heat across disjoint floating, stranded and solid
    /// source cells. No DC or circulating-current subtotal is deposited twice.
    pub fn mean_dissipated_power_w(&self) -> f64 {
        self.mean_eddy_power_w + self.mean_winding_power_w + self.mean_solid_power_w
    }
}

impl<'b, 'm> PlanarHarmonicWindings<'b, 'm> {
    /// Bind every cell's explicit role. `non_eddy_triangles` includes stranded
    /// winding cells: their intrinsic resistance is supplied separately by
    /// `winding_conductors`, exactly once. No end-turn resistance is inferred.
    pub fn new(
        basis: &'b PlanarIncrementalWindingBasis<'m>,
        eddy_conductors: Vec<EddyConductor>,
        non_eddy_triangles: &[usize],
        winding_conductors: Vec<WindingConductor>,
        controls: HarmonicWindingControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::new_with_solid_regions(
            basis,
            eddy_conductors,
            non_eddy_triangles,
            winding_conductors,
            vec![],
            controls,
            cancelled,
        )
    }

    /// Add explicitly driven solid conductor sections. Their existing DC
    /// profiles must equal n*sigma/integral(sigma dA); this constructor never
    /// substitutes a different winding load or operating field. Unbound
    /// conducting sections remain floating. Stranded and solid pieces may
    /// belong to one winding, but may never occupy the same source cell.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_solid_regions(
        basis: &'b PlanarIncrementalWindingBasis<'m>,
        eddy_conductors: Vec<EddyConductor>,
        non_eddy_triangles: &[usize],
        winding_conductors: Vec<WindingConductor>,
        solid_regions: Vec<SolidWindingRegion>,
        controls: HarmonicWindingControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let problem = basis.operating_field().problem();
        let mesh = problem.mesh();
        let n = mesh.vertices().len();
        let cells = mesh.faces().len();
        let m = basis.winding_count();
        let histories = basis
            .unit_responses()
            .iter()
            .try_fold(0usize, |sum, response| {
                sum.checked_add(
                    response
                        .report()
                        .history
                        .len()
                        .max(controls.max_cycles)
                        .checked_add(2)?,
                )
            });
        // Includes all retained unit fields, loads, a returned full response,
        // sparse tangent payload (at most nine entries per cell), heat receipts
        // retained solve histories and reused application buffers. The source owner's assembly budget
        // separately bounds constitutive assembly scratch.
        let payload = n
            .checked_mul(m)
            .and_then(|v| v.checked_mul(3))
            .and_then(|v| n.checked_mul(32).and_then(|w| v.checked_add(w)))
            .and_then(|v| cells.checked_mul(120).and_then(|w| v.checked_add(w)))
            .and_then(|v| cells.checked_mul(m).and_then(|w| v.checked_add(w)))
            .and_then(|v| {
                m.checked_mul(m)
                    .and_then(|w| w.checked_mul(8))
                    .and_then(|w| v.checked_add(w))
            })
            .and_then(|v| histories.and_then(|h| v.checked_add(h)));
        if payload.is_none_or(|v| v > controls.max_field_scalars)
            || !positive(controls.relative_tolerance)
            || controls.relative_tolerance >= 1.0
            || !positive(controls.absolute_residual_a_m)
            || !positive(controls.impedance_tolerance_ohm)
            || !positive(controls.power_tolerance_w)
            || !positive(controls.current_balance_tolerance_a)
            || !positive(controls.max_cell_skin_ratio)
            || controls.restart == 0
            || controls.max_cycles == 0
        {
            return Err(HarmonicWindingError::Invalid(
                "invalid harmonic controls or field payload budget",
            ));
        }
        let prescribed = PlanarMagnetostatics::admit_mesh_boundary(
            mesh,
            basis.boundary_potential_wb_m(),
            cancelled,
        )?;
        let free: Vec<_> = (0..n).filter(|&v| prescribed[v].is_none()).collect();
        let dimensions = free
            .len()
            .checked_mul(2)
            .ok_or(HarmonicWindingError::Invalid("harmonic dimension overflow"))?;
        let restart = controls.restart.min(dimensions.max(1));
        let workspace = restart
            .checked_mul(2)
            .and_then(|v| v.checked_add(12))
            .and_then(|v| v.checked_mul(dimensions))
            .and_then(|v| {
                restart
                    .checked_add(8)
                    .and_then(|r| r.checked_mul(restart))
                    .and_then(|h| v.checked_add(h))
            })
            .and_then(|v| {
                controls
                    .max_cycles
                    .checked_mul(4)
                    .and_then(|h| v.checked_add(h))
            });
        if workspace.is_none_or(|v| v > controls.max_krylov_scalars)
            || restart.checked_mul(controls.max_cycles).is_none()
        {
            return Err(HarmonicWindingError::Invalid(
                "harmonic Krylov budget exceeded",
            ));
        }
        let (eddy_conductors, conduction) = bind(problem, eddy_conductors, None, cancelled)?;
        let components = regions(problem, &eddy_conductors, non_eddy_triangles, cancelled)?;
        let solid = solid::SolidWindings::new(
            basis,
            &conduction,
            &eddy_conductors,
            solid_regions,
            controls,
            cancelled,
        )?;
        let mut is_eddy = vec![false; cells];
        for source in &eddy_conductors {
            is_eddy[source.triangle] = true;
        }
        if winding_conductors
            .iter()
            .any(|s| s.triangle < cells && is_eddy[s.triangle])
        {
            return Err(HarmonicWindingError::Invalid(
                "stranded winding and resolved eddy conductor cells must be disjoint",
            ));
        }
        // Only stranded pieces pass through the original fixed-profile copper
        // owner. Solid pieces are already source-bound by the conduction owner.
        let mut stranded_profiles = Vec::new();
        let mut stranded_windings = Vec::new();
        for (winding, profile) in basis.profiles_per_m2().iter().enumerate() {
            poll(cancelled)?;
            if profile
                .iter()
                .enumerate()
                .any(|(cell, &j)| j != 0.0 && !is_eddy[cell])
            {
                stranded_windings.push(winding);
                stranded_profiles.push(
                    profile
                        .iter()
                        .enumerate()
                        .map(|(cell, &j)| if is_eddy[cell] { 0.0 } else { j })
                        .collect::<Vec<_>>(),
                );
            }
        }
        let (winding_conductors, mut conductor_windings, stranded_resistance) = bind_sources(
            mesh,
            &stranded_profiles,
            winding_conductors,
            None,
            |source| {
                let material = source.conductor.material();
                match &problem.materials()[source.triangle] {
                    IsotropicMagneticMaterial::Linear(original) => {
                        admit_fixed_magnetic_refresh(original, material)?
                    }
                    IsotropicMagneticMaterial::Reversible(law) => {
                        law.admit_fixed_material_refresh(material)?
                    }
                    IsotropicMagneticMaterial::Recoil(law) => {
                        law.material().admit_fixed_material_refresh(material)?
                    }
                }
                Ok(())
            },
            cancelled,
        )?;
        for winding in &mut conductor_windings {
            *winding = stranded_windings[*winding];
        }
        let mut resistance_ohm = solid.resistance_ohm.clone();
        for (&winding, &resistance) in stranded_windings.iter().zip(&stranded_resistance) {
            resistance_ohm[winding] = finite(resistance_ohm[winding] + resistance)?;
        }
        if resistance_ohm.iter().any(|&r| !positive(r)) {
            return Err(HarmonicWindingError::Invalid(
                "every winding requires represented positive DC resistance",
            ));
        }
        let mut loads = Vec::with_capacity(m);
        for profile in basis.profiles_per_m2() {
            loads.push(load(
                problem,
                profile,
                &components,
                controls.current_balance_tolerance_a,
                cancelled,
            )?);
        }
        let mut skin_ratio_per_sqrt_rad_s = Vec::with_capacity(conduction.cells.len());
        for (source, cell) in eddy_conductors.iter().zip(&conduction.cells) {
            poll(cancelled)?;
            let permeability = match &problem.materials()[source.triangle] {
                IsotropicMagneticMaterial::Linear(material) => {
                    linear_magnetic_permeability(material).map_err(EddyError::Material)?
                }
                IsotropicMagneticMaterial::Recoil(law) => law.material().permeability_h_m(),
                IsotropicMagneticMaterial::Reversible(law) => {
                    let b = basis.operating_field().magnetic_flux_density_t()[source.triangle];
                    let response = law
                        .evaluate(b[0].hypot(b[1]))
                        .map_err(EddyError::Material)?;
                    response
                        .secant_reluctivity_m_h
                        .min(response.differential_reluctivity_m_h)
                        .recip()
                }
            };
            let points = cell.vertices.map(|v| mesh.vertices()[v]);
            let diameter = (0..3)
                .map(|k| {
                    let p = points[k];
                    let q = points[(k + 1) % 3];
                    (p[0] - q[0]).hypot(p[1] - q[1])
                })
                .fold(0.0_f64, f64::max);
            let coefficient = diameter * (0.5 * permeability * cell.conductivity_s_m).sqrt();
            if !positive(coefficient) {
                return Err(HarmonicWindingError::Invalid(
                    "unrepresentable differential skin depth",
                ));
            }
            skin_ratio_per_sqrt_rad_s.push(coefficient);
        }
        let (_, tangent) =
            problem.magnetic_linearization(basis.operating_field().potential_wb_m(), cancelled)?;
        poll(cancelled)?;
        Ok(Self {
            basis,
            tangent,
            free,
            loads,
            conduction,
            eddy_conductors,
            winding_conductors,
            conductor_windings,
            resistance_ohm,
            solid,
            skin_ratio_per_sqrt_rad_s,
            controls,
            restart,
        })
    }
    pub fn basis(&self) -> &'b PlanarIncrementalWindingBasis<'m> {
        self.basis
    }
    /// All spatially resolved conductive cells, including driven solid sections.
    pub fn eddy_conductors(&self) -> &[EddyConductor] {
        &self.eddy_conductors
    }
    pub fn winding_conductors(&self) -> &[WindingConductor] {
        &self.winding_conductors
    }
    /// Source-integrated stranded plus solid DC resistance. Frequency-dependent
    /// circulating-current loss is resolved by the field, not added here.
    pub fn resistance_ohm(&self) -> &[f64] {
        &self.resistance_ohm
    }
    pub fn solid_regions(&self) -> &[SolidWindingRegion] {
        &self.solid.bindings
    }
}
