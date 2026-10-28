//! Regional material viscosity with ALL retained modal cross terms.
//!
//! The original DKT operator projects C = Phi^T C_FE Phi on the SAME
//! mass-normalized modes used by force and pressure ports. Cholesky columns
//! express this tensor as passive dashpot channels in the existing coupled
//! modal solver. No diagonal damping approximation, modal refit, independent
//! eigensolve, added stiffness or alternative time integrator is introduced.
//! This is fixed-state, homogeneous-through-thickness, linear plate bending.
//! Coupled dashpot time discretization still needs refinement; material source
//! validity is neither a temporal accuracy certificate nor measured timbre.
//!
//! This module is also the ONE owner of the diagonal (energy-quotient) image of
//! the same plate viscosity: [`DiagonalPlateViscosity`] gives one decay ratio per
//! mode for scalar regional loss, and [`kelvin_voigt_zeta`] the uniform plate's
//! viscous term. Both keep the legacy arithmetic, element order and refusals.
//!
//! [`hysteretic`] is the same owner's frequency-independent (constant loss
//! factor) image for sourced wood loss: the plate stiffness split into named
//! shares, each share's modal strain energy projected once on the complete
//! source slice, and the diagonal loss a board reduction consumes, with the
//! neglected off-diagonal coupling measured.

use crate::acoustic_realize::AcousticRealizeError;
use crate::thin_plate::{CompactBody, ModalBendingLoss, PlateBendingLaw, PlateMaterialModel,
    RegionalBendingViscosity, ResolvedPlateChart, ResolvedPlateRegion, plate_property, viscous_time};
use crate::render::schedule::force::coupled::{ModalAttachment, ModalConnection};
use fs_exec::CancelGate;
use fs_material::plane_stress_viscosity::{MaterialPlane, PlaneStressViscosity};
use fs_material::state_point::YOUNG_MODULUS_PROPERTY;
use fs_blake3::ContentHash;
use fs_plate::PlateModel;
use fs_qty::Pressure;

/// Source-bound hereditary bending on the same original plate basis.
pub mod relaxation;
/// Constant-loss-factor (hysteretic) material loss by modal strain-energy shares.
pub mod hysteretic;

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn owner(error: impl core::fmt::Display) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(error.to_string())
}
fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() { Err(AcousticRealizeError::Cancelled) } else { Ok(()) }
}

/// One region's reduced plane-stress viscosity law, its applicability band
/// (viscous AND elastic/density coefficients) and its rotated DKT bending
/// tensor [N m s] in mesh coordinates. The single resolution used by the full
/// tensor image and by the diagonal image; Orthotropic21 exchanges 11/22 and
/// 16/26 like the elastic chart. No missing or zero coefficient is inferred.
pub(crate) fn resolve_region_viscosity(region: &ResolvedPlateRegion, properties: [&str; 6])
    -> Result<(PlaneStressViscosity, (f64, f64), [f64; 9]), AcousticRealizeError>
{
    let (plane, names) = if region.input.model == PlateMaterialModel::Orthotropic21 {
        (MaterialPlane::Axes21, [properties[1],properties[0],properties[2],properties[3],properties[5],properties[4]])
    } else { (MaterialPlane::Axes12, properties) };
    let law = PlaneStressViscosity::from_material(region.input.material.clone(),plane,names).map_err(owner)?;
    let mut band = law.omega_band_rad_s();
    let elastic: &[&str] = match region.input.model {
        PlateMaterialModel::Isotropic => &["density","young_modulus","poisson_ratio"],
        _ => &["density","young_modulus_1","young_modulus_2","poisson_ratio_12","shear_modulus_12"],
    };
    for key in elastic {
        let p = region.input.material.property(key).ok_or_else(|| invalid("viscous plate is missing its original elastic coefficient"))?;
        if !p.is_constant_at_fixed_temperature() {
            return Err(invalid("viscous plate needs fixed-state elastic and inertial coefficients"));
        }
        if let Some((lo,hi)) = p.angular_frequency_band_rad_s().map_err(owner)? {
            band = (band.0.max(lo),band.1.min(hi));
        }
    }
    if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0 {
        return Err(invalid("plate elastic and viscous applicability bands do not intersect"));
    }
    let bending = law.bending_matrix_n_m_s(region.section.thickness,region.input.material_angle_rad).map_err(owner)?;
    Ok((law, band, bending))
}

/// Viscous term of the uniform Kelvin-Voigt plate law at one frequency:
/// moment = D (curvature + tau curvature_rate), so zeta = tau omega / 2 before
/// the caller's bending/total energy fraction. Legacy arithmetic, moved here.
pub(crate) fn kelvin_voigt_zeta(viscous_time_s: f64, omega: f64) -> f64 {
    0.5 * viscous_time_s * omega
}

/// Bending (or viscously weighted bending) over total modal strain energy.
/// Stable compressive prestress may produce a ratio above one; it is not clamped.
pub(crate) fn bending_energy_fraction(bending: f64, total: f64) -> Result<f64, AcousticRealizeError> {
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

struct DiagonalElement {
    dofs: [Option<usize>; 9],
    stiffness: [f64; 81],
    region: usize,
}
enum DiagonalLaw {
    // Proportional isotropic bending: element matrix is the ELASTIC bending
    // tensor, weighted by the law's zeta(omega) = tau omega / 2.
    Isotropic(PlateBendingLaw),
    // Full rotated viscous tensor: element matrix is the VISCOUS bending tensor,
    // weighted by omega / 2, so the weighted energy is omega q^T C_e q / 2.
    Orthotropic((f64, f64)),
}

/// Diagonal weak-loss image of the selected regional bending viscosity on the
/// existing modes: zeta_i = sum_e w_r(omega) q_e^T K_e q_e / (phi^T K phi), with
/// K the original (prestressed) model. For an isotropic region K_e is the elastic
/// DKT bending matrix and w = tau omega / 2 (legacy, bit for bit); for a region
/// with a full reduced plane-stress viscosity K_e is the DKT matrix of its
/// rotated viscous tensor and w = omega / 2, which equals C_ii / (2 omega_i) of
/// the full projection for an exact eigenpair. Off-diagonal modal damping is
/// NOT represented; use [`PlateViscosity`] for the full tensor.
pub(crate) struct DiagonalPlateViscosity {
    laws: Vec<DiagonalLaw>,
    elements: Vec<DiagonalElement>,
}
impl DiagonalPlateViscosity {
    pub(crate) fn new(specimen: &ResolvedPlateChart, viscosities: &[RegionalBendingViscosity],
        model: &PlateModel) -> Result<Self, AcousticRealizeError>
    {
        let mut laws = Vec::with_capacity(viscosities.len());
        for (region, law) in specimen.regions().iter().zip(viscosities) {
            match law {
                RegionalBendingViscosity::Isotropic(viscous) => {
                    let viscous = *viscous;
                    let young = plate_property(&region.input.material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?;
                    laws.push(DiagonalLaw::Isotropic(PlateBendingLaw { thermal: None, viscous: Some(viscous),
                        viscous_time_s: viscous_time(viscous, young)?, thickness: region.section.thickness }));
                }
                RegionalBendingViscosity::Orthotropic(tensor) => {
                    laws.push(DiagonalLaw::Orthotropic(tensor.omega_band_rad_s()));
                }
            }
        }
        let tensors: Vec<&[f64; 9]> = viscosities.iter().enumerate().map(|(region, law)| match law {
            RegionalBendingViscosity::Isotropic(_) => &specimen.regions()[region].section.d,
            RegionalBendingViscosity::Orthotropic(tensor) => tensor.bending_matrix_n_m_s(),
        }).collect();
        Self::assemble(specimen, laws, &tensors, model)
    }

    /// The diagonal image of a full tensor selection, on the same modes. Every
    /// region takes the `Orthotropic` branch of [`Self::new`]: its rotated
    /// viscous DKT matrix with `w = omega / 2`. A region with the same tensor
    /// and band therefore gives the same bits through either constructor.
    pub(crate) fn from_tensor(specimen: &ResolvedPlateChart, viscosity: &PlateViscosity,
        model: &PlateModel) -> Result<Self, AcousticRealizeError>
    {
        if viscosity.chart_identity != specimen.specimen_identity() {
            return Err(invalid("diagonal tensor image needs its exact chart"));
        }
        let laws = viscosity.regions.iter().map(|r| DiagonalLaw::Orthotropic(r.band)).collect();
        let tensors: Vec<&[f64; 9]> = viscosity.regions.iter().map(|r| &r.bending).collect();
        Self::assemble(specimen, laws, &tensors, model)
    }

    // One element loop for both constructors.
    fn assemble(specimen: &ResolvedPlateChart, laws: Vec<DiagonalLaw>, tensors: &[&[f64; 9]],
        model: &PlateModel) -> Result<Self, AcousticRealizeError>
    {
        let mesh = &specimen.chart().mesh;
        let mut elements = Vec::with_capacity(mesh.tris.len());
        for (element, tri) in mesh.tris.iter().enumerate() {
            let region = specimen.region_of_element(element);
            let x = tri.map(|node| mesh.nodes[node].0);
            let y = tri.map(|node| mesh.nodes[node].1);
            // Original DKT element, same section, geometry and boundary map as
            // the modal pencil. No geometric stiffness enters viscous storage.
            let tensor = tensors[region];
            let stiffness = if tensor.iter().all(|v| *v == 0.0) { [0.0; 81] } else {
                fs_plate::dkt_stiffness(&x, &y, tensor, element).map_err(AcousticRealizeError::Plate)?.0
            };
            let dofs = std::array::from_fn(|i| model.dof_map[3 * tri[i / 3] + i % 3]);
            elements.push(DiagonalElement { dofs, stiffness, region });
        }
        Ok(Self { laws, elements })
    }
}
impl ModalBendingLoss for DiagonalPlateViscosity {
    fn modal_zeta(&self, total: &PlateModel, phi: &[f64], omega: f64)
        -> Result<f64, AcousticRealizeError>
    {
        // Check ALL source bands, including silent/zero-viscosity regions.
        let zetas: Vec<f64> = self.laws.iter().map(|law| match law {
            DiagonalLaw::Isotropic(law) => law.zeta(omega),
            DiagonalLaw::Orthotropic((lo, hi)) => {
                if !omega.is_finite() || omega < *lo || omega > *hi {
                    return Err(AcousticRealizeError::InvalidDescription {
                        what: "retained plate frequency lies outside material bending-loss applicability",
                    });
                }
                Ok(0.5 * omega)
            }
        }).collect::<Result<_, _>>()?;
        let mut k_phi = vec![0.0; phi.len()];
        total.k.spmv(phi, &mut k_phi);
        let total_energy = phi.iter().zip(k_phi).map(|(p, k)| p * k).sum();
        let mut weighted_bending = 0.0;
        // Fixed element order, not region-declaration order or hash order.
        for element in &self.elements {
            let q: [f64; 9] = element.dofs.map(|dof| dof.map_or(0.0, |i| phi[i]));
            let mut energy = 0.0;
            for r in 0..9 {
                for c in 0..9 {
                    energy += q[r] * element.stiffness[9 * r + c] * q[c];
                }
            }
            weighted_bending += zetas[element.region] * energy;
        }
        bending_energy_fraction(weighted_bending, total_energy)
    }
}

/// Cold projection limits and an explicit coordinate scale for dashpot channels.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateViscosityLimits {
    /// Maximum retained plate modes; bounded above by the 64 connection slots.
    pub max_modes: usize,
    /// Bound on element projections, retained vectors and cubic factor work.
    /// Counted conservatively as elements*81*modes^2 + DOFs*modes + 4*modes^3.
    pub max_projection_terms: usize,
    /// Positive resistance scale [N s/m]. Columns scale inversely with its
    /// square root so C is unchanged; connection force ceilings use this scale.
    pub reference_resistance_n_s_m: f64,
}

/// Original per-region viscosity, including source records and elastic support.
#[derive(Clone, Debug)]
pub struct ViscousPlateRegion {
    name: String,
    law: PlaneStressViscosity,
    band: (f64, f64),
    bending: [f64; 9],
}
impl ViscousPlateRegion {
    /// Original material region label.
    #[must_use] pub fn name(&self) -> &str { &self.name }
    /// Shared material owner and its exact source/axis/property records.
    #[must_use] pub const fn law(&self) -> &PlaneStressViscosity { &self.law }
    /// Viscosity AND elastic/density applicability [rad/s].
    #[must_use] pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.band }
    /// Original DKT engineering-curvature tensor [N m s] in mesh coordinates.
    #[must_use] pub const fn bending_matrix_n_m_s(&self) -> &[f64; 9] { &self.bending }
}

/// Complete source-bound loss selection for an existing immutable plate chart.
/// Its chart identity binds all original source states and element geometry.
#[derive(Clone, Debug)]
pub struct PlateViscosity {
    chart_identity: ContentHash,
    properties: [String; 6],
    regions: Vec<ViscousPlateRegion>,
    element_regions: Vec<usize>,
}
impl PlateViscosity {
    /// Select six explicit reduced plane-stress coefficients from EVERY region's
    /// original material bundle, in material 11,22,66,12,16,26 order. No missing
    /// or zero term is inferred. Orthotropic21 exchanges 11/22 and 16/26 using
    /// the same ordered-plane convention as the elastic chart. The material
    /// law rotates to each region's actual orientation and integrates its h^3.
    ///
    /// Existing scalar bending loss must be absent. Scalars or temperature-only
    /// curves at the original fixed state are allowed; frequency curves are not
    /// broadband constitutive laws. All elastic coefficients also limit the band.
    pub fn from_chart(chart: &ResolvedPlateChart, properties: [&str; 6], gate: &CancelGate)
        -> Result<Self, AcousticRealizeError>
    {
        poll(gate)?;
        if chart.bending_viscosity().is_some() || chart.regions().len() > 1024
            || chart.chart().mesh.tris.len() > 100_000 {
            return Err(invalid("tensor viscosity requires no scalar loss and bounded original regions/elements"));
        }
        let mut regions = Vec::with_capacity(chart.regions().len());
        let mut element_regions = vec![usize::MAX; chart.chart().mesh.tris.len()];
        for (i, region) in chart.regions().iter().enumerate() {
            poll(gate)?;
            let (law, band, bending) = resolve_region_viscosity(region, properties)?;
            for &element in &region.input.region.triangle_indices {
                poll(gate)?;
                let slot = element_regions.get_mut(element).ok_or_else(|| invalid("viscosity region names an absent element"))?;
                if *slot != usize::MAX { return Err(invalid("viscosity repeats an element")); }
                *slot = i;
            }
            regions.push(ViscousPlateRegion { name:region.input.region.name.clone(),law,band,bending });
        }
        if regions.is_empty() || element_regions.contains(&usize::MAX) {
            return Err(invalid("viscosity requires complete nonempty regional coverage"));
        }
        Ok(Self { chart_identity:chart.specimen_identity(), properties:properties.map(str::to_string), regions, element_regions })
    }
    /// Re-resolve the SAME law selection after a cold material/geometry rebind.
    /// Missing replacement coefficients refuse; no prior material loss is copied.
    pub fn rebind(&self, chart: &ResolvedPlateChart, gate: &CancelGate) -> Result<Self, AcousticRealizeError> {
        Self::from_chart(chart,self.properties.each_ref().map(String::as_str),gate)
    }
    /// Original data and orientation-resolved tensors in declared region order.
    #[must_use] pub fn regions(&self) -> &[ViscousPlateRegion] { &self.regions }
    /// Exact chart/source identity; not a scientific validity certificate.
    #[must_use] pub const fn chart_identity(&self) -> ContentHash { self.chart_identity }
    /// True only for explicitly sourced zero tensors in every region.
    #[must_use] pub fn is_zero(&self) -> bool { self.regions.iter().all(|r|r.bending.iter().all(|v|*v == 0.0)) }

    pub(crate) fn preflight(&self, chart:&ResolvedPlateChart, modes:usize,
        limits:PlateViscosityLimits, gate:&CancelGate) -> Result<(),AcousticRealizeError>
    {
        poll(gate)?;
        if chart.specimen_identity() != self.chart_identity || chart.bending_viscosity().is_some()
            || !(1..=64).contains(&limits.max_modes) || modes == 0 || modes > limits.max_modes
            || !limits.reference_resistance_n_s_m.is_finite() || limits.reference_resistance_n_s_m <= 0.0 {
            return Err(invalid("tensor plate reduction needs its exact chart, bounded modes and a positive resistance scale"));
        }
        let work = chart.chart().mesh.tris.len().checked_mul(81)
            .and_then(|w|w.checked_mul(modes)).and_then(|w|w.checked_mul(modes))
            .and_then(|w|chart.chart().mesh.nodes.len().checked_mul(3*modes).and_then(|v|w.checked_add(v)))
            .and_then(|w|w.checked_add(4*modes*modes*modes))
            .ok_or_else(||invalid("tensor plate projection work overflows"))?;
        if work > limits.max_projection_terms { return Err(invalid("tensor plate projection exceeds its declared work budget")); }
        Ok(())
    }

    // Called only from the original eigensolve, with its actual normalized vectors.
    pub(crate) fn project(&self, chart:&ResolvedPlateChart, model:&PlateModel,
        phi:&[Vec<f64>], omega:&[f64], limits:PlateViscosityLimits, gate:&CancelGate)
        -> Result<ModalPlateViscosity,AcousticRealizeError>
    {
        let n = self.admit(chart,phi,omega,limits,gate)?;
        let mut matrix = vec![0.0; n*n];
        let mesh = &chart.chart().mesh;
        for (e,tri) in mesh.tris.iter().enumerate() {
            poll(gate)?;
            let d = self.regions[self.element_regions[e]].bending;
            if d.iter().all(|v|*v == 0.0) { continue; }
            let (c,_) = fs_plate::dkt_stiffness(&tri.map(|v|mesh.nodes[v].0),
                &tri.map(|v|mesh.nodes[v].1),&d,e).map_err(owner)?;
            let dofs:[Option<usize>;9] = std::array::from_fn(|i|model.dof_map[3*tri[i/3]+i%3]);
            let local:Vec<[f64;9]> = phi.iter().map(|p|dofs.map(|v|v.map_or(0.0,|i|p[i]))).collect();
            for i in 0..n {
                poll(gate)?;
                for j in i..n {
                    let mut value = 0.0;
                    for a in 0..9 { for b in 0..9 { value += local[i][a]*c[9*a+b]*local[j][b]; } }
                    matrix[n*i+j] += value; matrix[n*j+i] = matrix[n*i+j];
                }
            }
        }
        ModalPlateViscosity::factor(matrix,n,limits.reference_resistance_n_s_m,gate)
    }

    // The actual retained modes every projection admits.
    fn admit(&self, chart:&ResolvedPlateChart, phi:&[Vec<f64>], omega:&[f64],
        limits:PlateViscosityLimits, gate:&CancelGate) -> Result<usize,AcousticRealizeError>
    {
        let n = phi.len(); self.preflight(chart,n,limits,gate)?;
        if omega.len() != n || phi.iter().any(|v|v.iter().any(|x|!x.is_finite())) {
            return Err(invalid("viscous projection needs aligned finite modal coordinates"));
        }
        for region in &self.regions {
            poll(gate)?;
            for &w in omega {
                if !w.is_finite() || w <= 0.0 || w < region.band.0 || w > region.band.1 {
                    return Err(invalid("actual retained plate mode is outside original elastic/viscous source support"));
                }
            }
        }
        Ok(n)
    }

    /// Windowed sibling of `project` on the SAME modes (IRP M0.2-W5 scalable
    /// path). Entries are formed only for `|w_i - w_k| < window` [rad/s]; the
    /// rest of each row is bounded, never formed.
    ///
    /// `x_j = C_FE phi_j` comes from the same element loop, with a bound on its
    /// own rounding. `fs_modal::residual_far_row_bounds` then forms the window
    /// entries and bounds every other entry of the row by Bessel on the
    /// lumped-mass residual `r_j = M^-1 x_j - sum_{k in near+j} phi_k C_kj`,
    /// charging the measured defect of `Phi^T M Phi = I`.
    pub(crate) fn project_windowed(&self, chart:&ResolvedPlateChart, model:&PlateModel,
        phi:&[Vec<f64>], omega:&[f64], window:f64, limits:PlateViscosityLimits, gate:&CancelGate)
        -> Result<WindowedPlateViscosity,AcousticRealizeError>
    {
        let n = self.admit(chart,phi,omega,limits,gate)?;
        let free = model.free;
        if phi.iter().any(|v|v.len() != free) { return Err(invalid("windowed viscosity needs reduced-DOF modal vectors")); }
        let mut mass = vec![0.0; free];
        for (d,slot) in mass.iter_mut().enumerate() {
            let (cols,vals) = model.m.row(d);
            for (&c,&v) in cols.iter().zip(vals) {
                if c == d { *slot = v; }
                else if v != 0.0 { return Err(invalid("windowed viscosity bound needs the lumped diagonal plate mass")); }
            }
            if !slot.is_finite() || *slot <= 0.0 { return Err(invalid("windowed viscosity bound needs a positive plate mass")); }
        }
        // x_j and |C_FE||phi_j|, with the per-DOF term count, for its rounding.
        let mesh = &chart.chart().mesh;
        let (mut x, mut size) = (vec![vec![0.0; free]; n], vec![vec![0.0; free]; n]);
        let mut terms = vec![0usize; free];
        for (e,tri) in mesh.tris.iter().enumerate() {
            poll(gate)?;
            let d = self.regions[self.element_regions[e]].bending;
            if d.iter().all(|v|*v == 0.0) { continue; }
            let (c,_) = fs_plate::dkt_stiffness(&tri.map(|v|mesh.nodes[v].0),
                &tri.map(|v|mesh.nodes[v].1),&d,e).map_err(owner)?;
            let dofs:[Option<usize>;9] = std::array::from_fn(|i|model.dof_map[3*tri[i/3]+i%3]);
            for &row in dofs.iter().flatten() { terms[row] += 9; }
            for j in 0..n {
                let local = dofs.map(|v|v.map_or(0.0,|i|phi[j][i]));
                for (a,row) in dofs.iter().enumerate() {
                    let Some(row) = *row else { continue };
                    for b in 0..9 { let t = c[9*a+b]*local[b]; x[j][row] += t; size[j][row] += t.abs(); }
                }
            }
        }
        // |dx_d| <= gamma_K (|C_FE||phi_j|)_d with K terms summed into any DOF.
        let k_eps = terms.iter().copied().max().unwrap_or(0) as f64*f64::EPSILON;
        let applied_error:Vec<f64> = size.iter().map(|s| {
            let norm = s.iter().zip(&mass).map(|(v,m)|v*v/m).sum::<f64>().sqrt();
            1.01*k_eps/(1.0-k_eps)*norm
        }).collect();
        poll(gate)?;
        let bounds = fs_modal::residual_far_row_bounds(omega,window,&mass,phi,&x,&applied_error).map_err(owner)?;
        poll(gate)?;
        Ok(WindowedPlateViscosity { omega:omega.to_vec(), window, near:bounds.near, far:bounds.far_row_norm,
            column:bounds.column_norm, defect:bounds.basis_defect })
    }
}

/// Windowed material damping on the retained plate modes (IRP M0.2-W5). Only
/// entries with `|w_i - w_k| < window` exist; Bessel on the mass-norm residual
/// bounds the rest of each row. The dense modal matrix is never formed.
#[derive(Clone,Debug)]
pub struct WindowedPlateViscosity {
    omega:Vec<f64>, window:f64, near:Vec<Vec<(usize,f64)>>, far:Vec<f64>, column:Vec<f64>, defect:f64,
}
impl WindowedPlateViscosity {
    /// Undamped frequencies of the retained modes [rad/s].
    #[must_use] pub fn omega_rad_s(&self)->&[f64] {&self.omega}
    /// Coupling window [rad/s].
    #[must_use] pub const fn window_rad_s(&self)->f64 {self.window}
    /// `C_ij` [1/s] when it was formed: the diagonal or inside the window.
    #[must_use] pub fn entry(&self,i:usize,j:usize)->Option<f64> {
        let row = self.near.get(i)?;
        row.binary_search_by_key(&j,|&(k,_)|k).ok().map(|at|row[at].1)
    }
    /// Bound on the norm of each row's unformed off-diagonal entries [1/s].
    #[must_use] pub fn far_row_norm_bound(&self)->&[f64] {&self.far}
    /// Bound on each column norm `||C e_j||` over every plate mode [1/s].
    #[must_use] pub fn column_norm_bound(&self)->&[f64] {&self.column}
    /// Bound on `||Phi^T M Phi - I||_2` of the retained modes, charged to both bounds.
    #[must_use] pub const fn basis_defect(&self)->f64 {self.defect}
    /// Certify from the formed entries and the residual far bound.
    ///
    /// # Errors
    /// Any `fs_modal::certify_modal_damping_windowed_far_rows` refusal.
    pub fn certify(&self,budget:&fs_modal::ModalDampingBudget)
        -> Result<fs_modal::CertifiedModalDamping,fs_modal::ModalError>
    {
        let near = |i:usize,j:usize| self.entry(i,j).unwrap_or(f64::NAN);
        fs_modal::certify_modal_damping_windowed_far_rows(&self.omega,&near,&self.far,self.window,budget)
    }
}

/// Full mass-normalized C [1/s] and its equivalent passive modal channels.
/// The tensor is not reduced to one decay rate per mode.
#[derive(Clone,Debug)]
pub struct ModalPlateViscosity {
    matrix:Vec<f64>, columns:Vec<Vec<f64>>, resistance:f64, modes:usize, diagonal_zeta:Option<Vec<f64>>,
}
impl ModalPlateViscosity {
    pub(crate) fn factor(matrix:Vec<f64>,n:usize,resistance:f64,gate:&CancelGate)->Result<Self,AcousticRealizeError> {
        if matrix.iter().any(|v|!v.is_finite()) { return Err(invalid("projected material damping is not finite")); }
        if matrix.iter().all(|v|*v == 0.0) {
            return Ok(Self {matrix,columns:Vec::new(),resistance,modes:n,diagonal_zeta:None});
        }
        let scale:Vec<_> = (0..n).map(|i|matrix[n*i+i].sqrt()).collect();
        if scale.iter().any(|s|!s.is_finite() || *s <= 0.0) {
            return Err(invalid("nonzero projected viscosity requires resolved positive-definite modal damping"));
        }
        let normalized:Vec<_> = (0..n*n).map(|k|matrix[k]/scale[k/n]/scale[k%n]).collect();
        let root = fs_la::factor::cholesky(&normalized,n).map_err(owner)?;
        let r = resistance.sqrt(); let mut columns = Vec::with_capacity(n);
        for j in 0..n {
            poll(gate)?;
            let column:Vec<_> = (0..n).map(|i|scale[i]*root.l(i,j)/r).collect();
            if column.iter().any(|v|!v.is_finite()) || column[j] == 0.0 {
                return Err(invalid("material damping channel scale is unrepresentable"));
            }
            columns.push(column);
        }
        // Check the realized C, not only its factorization; never drop cross terms.
        for i in 0..n { for j in 0..n {
            poll(gate)?;
            let actual:f64 = columns.iter().map(|b|(b[i]*r)*(b[j]*r)).sum();
            let allowance = 1e-10*(matrix[n*i+i].sqrt()*matrix[n*j+j].sqrt());
            if !actual.is_finite() || !allowance.is_finite() || (actual-matrix[n*i+j]).abs() > allowance {
                return Err(invalid("material damping factor does not reconstruct its original full tensor"));
            }
        } }
        Ok(Self {matrix,columns,resistance,modes:n,diagonal_zeta:None})
    }
    /// Attach the one diagonal evaluator's zeta per mode (IRP M0.2-W3/W5),
    /// computed by the caller on the SAME eigenpairs.
    pub(crate) fn with_diagonal_zeta(mut self, zeta:Vec<f64>)->Self { self.diagonal_zeta = Some(zeta); self }
    /// `DiagonalPlateViscosity`'s zeta on these modes, when the projecting
    /// eigensolve supplied it; certified Diagonal blocks report exactly these.
    #[must_use] pub fn diagonal_zeta(&self)->Option<&[f64]> { self.diagonal_zeta.as_deref() }
    /// Dense symmetric C [1/s], row-major in the original mass-normalized basis.
    #[must_use] pub fn matrix_per_s(&self)->&[f64] { &self.matrix }
    /// Mode count, including modes with explicit zero dissipation.
    #[must_use] pub const fn mode_count(&self)->usize {self.modes}
    /// Number of passive channels. Zero tensors allocate no artificial damping.
    #[must_use] pub fn channel_count(&self)->usize {self.columns.len()}
    /// Borrow the admitted passive columns without rebuilding connection objects.
    /// Their reference resistance is the one supplied to this factor's constructor.
    pub(crate) fn channel_shapes(&self)->&[Vec<f64>] { &self.columns }
    /// Reference resistance of these passive channel shapes [N s/m].
    pub(crate) fn reference_resistance_n_s_m(&self)->f64 { self.resistance }
    /// Realize C through the original bilateral dashpot law on one component.
    /// Its coordinate origin is zero; no added spring, mass, or external forcing.
    /// The same component on both sides keeps the complete coordinate order.
    #[must_use] pub fn connections(&self,component:usize)->Vec<ModalConnection> {
        self.columns.iter().map(|b|ModalConnection {
            left:ModalAttachment {component,shapes:b.clone()},
            right:ModalAttachment {component,shapes:vec![0.0;self.modes]},
            stiffness_n_m:0.0,damping_n_s_m:self.resistance,rest_extension_m:0.0,
        }).collect()
    }

    /// Factor a symmetric positive-semidefinite (PSD) modal damping matrix $C_{\text{modal}}$ (`n × n`,
    /// row-major [1/s]) into passive dashpot channels at reference resistance `resistance` [$\text{N}\cdot\text{s/m}$].
    ///
    /// Unlike bulk bending viscosity (where every retained elastic mode has strictly positive diagonal
    /// damping), localized interface/glue damping may be rank-deficient or vanish on rigid-body modes.
    /// Complete pivoted Cholesky extracts one passive channel per active rank direction and verifies
    /// that $\sum_k (\sqrt{R}\,\mathbf{b}_k)(\sqrt{R}\,\mathbf{b}_k)^T$ reconstructs $C_{\text{modal}}$ within $10^{-10}$.
    pub fn from_psd_modal_matrix(
        matrix: Vec<f64>,
        n: usize,
        resistance: f64,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        if n == 0
            || matrix.len() != n * n
            || !resistance.is_finite()
            || resistance <= 0.0
            || matrix.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "PSD modal damping requires finite n*n matrix and positive resistance scale",
            ));
        }
        let mut max_diag = 0.0_f64;
        for i in 0..n {
            let d_ii = matrix[i * n + i];
            if d_ii < -1e-12 {
                return Err(invalid("PSD modal damping has a negative diagonal entry"));
            }
            max_diag = max_diag.max(d_ii.max(0.0));
        }
        if max_diag == 0.0 {
            return Ok(Self {
                matrix: vec![0.0; n * n],
                columns: Vec::new(),
                resistance,
                modes: n,
                diagonal_zeta: None,
            });
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if (matrix[i * n + j] - matrix[j * n + i]).abs() > 1e-10 * max_diag {
                    return Err(invalid("PSD modal damping matrix is not symmetric"));
                }
            }
        }

        // Complete diagonal-pivoted Cholesky on the symmetric PSD matrix.
        let mut rem = matrix.clone();
        for i in 0..n {
            if rem[i * n + i] < 0.0 && rem[i * n + i] >= -1e-12 * max_diag {
                rem[i * n + i] = 0.0;
            }
        }
        let r_sqrt = resistance.sqrt();
        // Pivots below this factor's own PSD, symmetry and reconstruction
        // tolerance (1e-10 of the largest diagonal) are roundoff in the modal
        // projection, not damping channels. Rigid-mode rows of a glue box
        // measured 4e-13 and -7e-14 relative, and at 1e-13 one became a
        // spurious seventh channel.
        let pivot_tol = 1e-10 * max_diag;
        let mut columns: Vec<Vec<f64>> = Vec::new();
        let mut used = vec![false; n];

        for _step in 0..n {
            poll(gate)?;
            let mut best_p = None;
            let mut best_diag = pivot_tol;
            for i in 0..n {
                if !used[i] {
                    let d_i = rem[i * n + i];
                    if d_i < -1e-10 * max_diag {
                        return Err(invalid(
                            "modal damping matrix is not positive semidefinite (negative Cholesky pivot)",
                        ));
                    }
                    if d_i > best_diag {
                        best_diag = d_i;
                        best_p = Some(i);
                    }
                }
            }
            let Some(p) = best_p else {
                break;
            };
            used[p] = true;
            let root_p = best_diag.sqrt();
            let mut col = vec![0.0; n];
            col[p] = root_p / r_sqrt;
            for i in 0..n {
                if !used[i] {
                    col[i] = (rem[i * n + p] / root_p) / r_sqrt;
                }
            }
            for i in 0..n {
                if !used[i] {
                    let l_ip = col[i] * r_sqrt;
                    for j in i..n {
                        if !used[j] {
                            let l_jp = col[j] * r_sqrt;
                            rem[i * n + j] -= l_ip * l_jp;
                            rem[j * n + i] = rem[i * n + j];
                        }
                    }
                }
            }
            rem[p * n + p] = 0.0;
            columns.push(col);
        }

        // Verify remaining diagonal entries are non-negative within tolerance and full matrix reconstructs.
        for i in 0..n {
            if !used[i] && rem[i * n + i] < -1e-10 * max_diag {
                return Err(invalid(
                    "modal damping matrix is not positive semidefinite (negative residual diagonal)",
                ));
            }
        }
        let allowance = 1e-10 * max_diag;
        for i in 0..n {
            for j in 0..n {
                let actual: f64 = columns
                    .iter()
                    .map(|b| (b[i] * r_sqrt) * (b[j] * r_sqrt))
                    .sum();
                if !actual.is_finite() || (actual - matrix[i * n + j]).abs() > allowance {
                    return Err(invalid(
                        "PSD modal damping factorization does not reconstruct the full modal tensor",
                    ));
                }
            }
        }
        Ok(Self {
            matrix,
            columns,
            resistance,
            modes: n,
            diagonal_zeta: None,
        })
    }

    /// Project the Kelvin–Voigt glue damping operator $C_{\text{glue}}$ of a [`fs_plate::shell::junction::TiedShellModel`]
    /// onto mass-normalized reduced mode shapes `modes` through this single modal viscous owner (`M0.2-W3` / `M1.2-W2b`).
    pub fn from_tied_shell_glue(
        tied: &fs_plate::shell::junction::TiedShellModel,
        modes: &[Vec<f64>],
        resistance: f64,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let modal_c = tied.project_glue_damping(modes).map_err(owner)?;
        Self::from_psd_modal_matrix(modal_c, modes.len(), resistance, gate)
    }
}

/// Canonical `fs-matdb` claim key for area-specific normal glue stiffness [$\text{N/m}^3$].
pub const GLUE_NORMAL_STIFFNESS_PROPERTY: &str = "glue-normal-stiffness";
/// Canonical `fs-matdb` claim key for area-specific seam-tangent shear glue stiffness [$\text{N/m}^3$].
pub const GLUE_SHEAR_STIFFNESS_1_PROPERTY: &str = "glue-shear-stiffness-1";
/// Canonical `fs-matdb` claim key for area-specific transverse shear glue stiffness [$\text{N/m}^3$].
pub const GLUE_SHEAR_STIFFNESS_2_PROPERTY: &str = "glue-shear-stiffness-2";
/// Canonical `fs-matdb` claim key for area-specific normal glue viscosity [$\text{N}\cdot\text{s/m}^3$].
pub const GLUE_NORMAL_VISCOSITY_PROPERTY: &str = "glue-normal-viscosity";
/// Canonical `fs-matdb` claim key for area-specific seam-tangent shear glue viscosity [$\text{N}\cdot\text{s/m}^3$].
pub const GLUE_SHEAR_VISCOSITY_1_PROPERTY: &str = "glue-shear-viscosity-1";
/// Canonical `fs-matdb` claim key for area-specific transverse shear glue viscosity [$\text{N}\cdot\text{s/m}^3$].
pub const GLUE_SHEAR_VISCOSITY_2_PROPERTY: &str = "glue-shear-viscosity-2";

/// Evidence-bearing resolution of an ordered [`fs_matdb::InterfaceSystemCard`] into an
/// [`fs_plate::shell::junction::GlueLine`].
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedGlueLine {
    /// Physical Kelvin–Voigt glue line parameters consumed by `fs-plate::shell::junction`.
    pub glue_line: fs_plate::shell::junction::GlueLine,
    /// Canonical content hash of the source [`fs_matdb::InterfaceSystemCard`].
    pub card_identity: ContentHash,
    /// Ordered surface A descriptor.
    pub surface_a: fs_matdb::SurfaceSpec,
    /// Ordered surface B descriptor.
    pub surface_b: fs_matdb::SurfaceSpec,
    /// Third-body adhesive identifier from `SystemContext.third_body`.
    pub adhesive_third_body: String,
    /// Named history state from `SystemContext.history`.
    pub history: String,
    /// Usage receipts for `[normal_n_m3, shear_n_m3[0], shear_n_m3[1]]`.
    pub stiffness_receipts: [fs_matdb::PropertyUsageReceipt; 3],
    /// Source uncertainty models for `[normal_n_m3, shear_n_m3[0], shear_n_m3[1]]`.
    pub stiffness_uncertainties: [fs_matdb::UncertaintyModel; 3],
    /// Optional usage receipts for `[viscosity_n_s_m3[0..3]]`.
    pub viscosity_receipts: Option<[fs_matdb::PropertyUsageReceipt; 3]>,
    /// Optional source uncertainty models for `[viscosity_n_s_m3[0..3]]`.
    pub viscosity_uncertainties: Option<[fs_matdb::UncertaintyModel; 3]>,
}

impl ResolvedGlueLine {
    /// Resolve an ordered [`fs_matdb::InterfaceSystemCard`] into a [`fs_plate::shell::junction::GlueLine`]
    /// with full query receipts, unit checks (`N/m^3` and `N*s/m^3` via `fs-qty`), and third-body/history validation.
    ///
    /// - `expected_surfaces`: ordered `(surface_a_material, surface_b_material)` pair; a swapped pair refuses.
    /// - `require_viscosity`: if `true`, resolves the three Kelvin–Voigt viscosity claims; if `false`, resolves
    ///   stiffness-only (`viscosity_n_s_m3 = None`).
    /// - `width_m`: physical glue-line width [m] (must be finite and $> 0$).
    pub fn from_card(
        card: &fs_matdb::InterfaceSystemCard,
        expected_surfaces: (&fs_matdb::MaterialStateId, &fs_matdb::MaterialStateId),
        point: &fs_matdb::QueryPoint,
        policy: fs_matdb::SelectionPolicy,
        width_m: f64,
        require_viscosity: bool,
    ) -> Result<Self, AcousticRealizeError> {
        if &card.surface_a().material != expected_surfaces.0
            || &card.surface_b().material != expected_surfaces.1
        {
            return Err(AcousticRealizeError::Nonlinear(format!(
                "ordered glue card surface mismatch: expected ({:?}, {:?}), found ({:?}, {:?})",
                expected_surfaces.0,
                expected_surfaces.1,
                card.surface_a().material,
                card.surface_b().material
            )));
        }
        let Some(adhesive) = card.third_body() else {
            return Err(AcousticRealizeError::Nonlinear(
                "glue interface card requires SystemContext.third_body naming the adhesive".to_string(),
            ));
        };
        if adhesive.trim().is_empty() || card.history().trim().is_empty() {
            return Err(AcousticRealizeError::Nonlinear(
                "glue interface card requires non-empty adhesive third_body and named history".to_string(),
            ));
        }
        if !width_m.is_finite() || width_m <= 0.0 {
            return Err(AcousticRealizeError::Nonlinear(
                "glue line width_m must be finite and strictly positive".to_string(),
            ));
        }

        let k_dims = crate::shear_shell::connections::stiffness_per_area_dims();
        let c_dims = crate::shear_shell::connections::viscosity_per_area_dims();

        let query_claim = |prop_name: &str,
                           expected_dims: fs_qty::Dims,
                           allow_zero: bool|
         -> Result<
            (
                f64,
                fs_matdb::PropertyUsageReceipt,
                fs_matdb::UncertaintyModel,
            ),
            AcousticRealizeError,
        > {
            let ans = card.claims().query(prop_name, point, policy).map_err(|e| {
                AcousticRealizeError::Nonlinear(format!(
                    "missing or invalid glue property {prop_name:?}: {e}"
                ))
            })?;
            let sample = &ans.evidence.value;
            if sample.dims != expected_dims {
                return Err(AcousticRealizeError::Nonlinear(format!(
                    "unit mismatch for glue property {prop_name:?}: expected dims {:?}, found {:?}",
                    expected_dims.0, sample.dims.0
                )));
            }
            let val = sample.value;
            let ok = if allow_zero {
                val.is_finite() && val >= 0.0
            } else {
                val.is_finite() && val > 0.0
            };
            if !ok {
                return Err(AcousticRealizeError::Nonlinear(format!(
                    "glue property {prop_name:?} value {val} is out of admissible range"
                )));
            }
            Ok((val, ans.receipt, sample.uncertainty.clone()))
        };

        let (kn, rk_n, uk_n) = query_claim(GLUE_NORMAL_STIFFNESS_PROPERTY, k_dims, false)?;
        let (ks1, rk_s1, uk_s1) = query_claim(GLUE_SHEAR_STIFFNESS_1_PROPERTY, k_dims, false)?;
        let (ks2, rk_s2, uk_s2) = query_claim(GLUE_SHEAR_STIFFNESS_2_PROPERTY, k_dims, false)?;

        let (visc_vals, visc_receipts, visc_uncert) = if require_viscosity {
            let (vn, rv_n, uv_n) = query_claim(GLUE_NORMAL_VISCOSITY_PROPERTY, c_dims, true)?;
            let (vs1, rv_s1, uv_s1) = query_claim(GLUE_SHEAR_VISCOSITY_1_PROPERTY, c_dims, true)?;
            let (vs2, rv_s2, uv_s2) = query_claim(GLUE_SHEAR_VISCOSITY_2_PROPERTY, c_dims, true)?;
            (
                Some([vn, vs1, vs2]),
                Some([rv_n, rv_s1, rv_s2]),
                Some([uv_n, uv_s1, uv_s2]),
            )
        } else {
            (None, None, None)
        };

        let glue_line = fs_plate::shell::junction::GlueLine {
            normal_n_m3: kn,
            shear_n_m3: [ks1, ks2],
            viscosity_n_s_m3: visc_vals,
            width_m,
        };
        glue_line.validate().map_err(owner)?;

        Ok(Self {
            glue_line,
            card_identity: card.content_hash(),
            surface_a: card.surface_a().clone(),
            surface_b: card.surface_b().clone(),
            adhesive_third_body: adhesive.to_string(),
            history: card.history().to_string(),
            stiffness_receipts: [rk_n, rk_s1, rk_s2],
            stiffness_uncertainties: [uk_n, uk_s1, uk_s2],
            viscosity_receipts: visc_receipts,
            viscosity_uncertainties: visc_uncert,
        })
    }
}

#[cfg(test)]
mod glue_card_tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
        PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
        SystemContext, UncertaintyModel,
    };

    fn fixture_surfaces() -> (MaterialStateId, MaterialStateId) {
        (
            MaterialStateId {
                chemistry: "Picea-abies-spruce".to_string(),
                phase: "cellular-orthotropic-wood".to_string(),
                process: "quarter-sawn-air-dried-12pct-mc".to_string(),
                revision: 1,
            },
            MaterialStateId {
                chemistry: "Acer-pseudoplatanus-maple".to_string(),
                phase: "cellular-orthotropic-wood".to_string(),
                process: "quarter-sawn-rib-bent-12pct-mc".to_string(),
                revision: 1,
            },
        )
    }

    fn make_claim(name: &str, value: f64, dims: fs_qty::Dims) -> PropertyClaim {
        PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("temperature", 270.0, 330.0),
            uncertainty: UncertaintyModel::RelativeHalfWidth {
                fraction: 0.1,
                confidence: 0.95,
            },
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "authored-fixture:hide-glue-spruce-maple-v1".to_string(),
                license: "CC0-1.0".to_string(),
                artifact: None,
            },
        }
    }

    fn assemble_fixture_card(omit_key: Option<&str>, bad_dims_key: Option<&str>) -> InterfaceSystemCard {
        let (mat_a, mat_b) = fixture_surfaces();
        let k_dims = crate::shear_shell::connections::stiffness_per_area_dims();
        let c_dims = crate::shear_shell::connections::viscosity_per_area_dims();
        let wrong_dims = fs_qty::Pressure::DIMS; // Pa = N/m^2 instead of N/m^3

        let mut claims = ClaimSet::new();
        let entries = [
            (GLUE_NORMAL_STIFFNESS_PROPERTY, 2.0e11, k_dims),
            (GLUE_SHEAR_STIFFNESS_1_PROPERTY, 1.0e11, k_dims),
            (GLUE_SHEAR_STIFFNESS_2_PROPERTY, 1.0e11, k_dims),
            (GLUE_NORMAL_VISCOSITY_PROPERTY, 4.0e5, c_dims),
            (GLUE_SHEAR_VISCOSITY_1_PROPERTY, 2.0e5, c_dims),
            (GLUE_SHEAR_VISCOSITY_2_PROPERTY, 2.0e5, c_dims),
        ];
        for (key, val, dims) in entries {
            if omit_key == Some(key) {
                continue;
            }
            let d = if bad_dims_key == Some(key) {
                wrong_dims
            } else {
                dims
            };
            claims.insert_claim(make_claim(key, val, d)).unwrap();
        }
        InterfaceSystemCard::assemble(
            SurfaceSpec {
                material: mat_a,
                texture_frame: "planed-toothed-0.1mm-L".to_string(),
            },
            SurfaceSpec {
                material: mat_b,
                texture_frame: "scraped-rib-edge-L".to_string(),
            },
            SystemContext {
                medium: "cured-bondline-0.08mm".to_string(),
                third_body: Some("authored-fixture:hot-hide-glue-192g-bloom".to_string()),
                environment: "air-293K-50pctRH".to_string(),
                history: "cured-72h-293K-virgin".to_string(),
            },
            claims,
            Vec::new(),
        )
        .unwrap()
    }

    #[test]
    fn g0_glue_card_missing_property_refuses_by_name() {
        let (mat_a, mat_b) = fixture_surfaces();
        let point = QueryPoint::new().with("temperature", 293.15).unwrap();
        for missing in [
            GLUE_NORMAL_STIFFNESS_PROPERTY,
            GLUE_SHEAR_STIFFNESS_1_PROPERTY,
            GLUE_SHEAR_STIFFNESS_2_PROPERTY,
            GLUE_NORMAL_VISCOSITY_PROPERTY,
            GLUE_SHEAR_VISCOSITY_1_PROPERTY,
            GLUE_SHEAR_VISCOSITY_2_PROPERTY,
        ] {
            let card = assemble_fixture_card(Some(missing), None);
            let err = ResolvedGlueLine::from_card(
                &card,
                (&mat_a, &mat_b),
                &point,
                SelectionPolicy::SingleClaimOnly,
                0.003,
                true,
            )
            .unwrap_err();
            let msg = err.to_string();
            assert!(
                msg.contains(missing),
                "expected refusal message to name missing property {missing:?}, got: {msg}"
            );
        }
    }

    #[test]
    fn g0_glue_card_units_checked() {
        let (mat_a, mat_b) = fixture_surfaces();
        let point = QueryPoint::new().with("temperature", 293.15).unwrap();
        for bad_key in [
            GLUE_NORMAL_STIFFNESS_PROPERTY,
            GLUE_SHEAR_STIFFNESS_1_PROPERTY,
            GLUE_NORMAL_VISCOSITY_PROPERTY,
        ] {
            let card = assemble_fixture_card(None, Some(bad_key));
            let err = ResolvedGlueLine::from_card(
                &card,
                (&mat_a, &mat_b),
                &point,
                SelectionPolicy::SingleClaimOnly,
                0.003,
                true,
            )
            .unwrap_err();
            let msg = err.to_string();
            assert!(
                msg.contains("unit mismatch") && msg.contains(bad_key),
                "expected unit mismatch naming {bad_key:?}, got: {msg}"
            );
        }
    }
}

/// One eigensolve's physical radiators and full damping, with original sources.
#[derive(Clone,Debug)]
pub struct ViscousPlateImage {
    pub(crate) source:PlateViscosity,
    pub(crate) radiators:Vec<CompactBody>,
    pub(crate) damping:ModalPlateViscosity,
}
impl ViscousPlateImage {
    /// Material assignments that produced this exact damping and geometry image.
    #[must_use] pub const fn source(&self)->&PlateViscosity {&self.source}
    /// Modes/masses/signed force and radiation projections from the SAME solve.
    /// These have zero scalar material loss; use damping connections exactly once.
    #[must_use] pub fn radiators(&self)->&[CompactBody] {&self.radiators}
    /// Full nonproportional damping projected on these exact modes.
    #[must_use] pub const fn damping(&self)->&ModalPlateViscosity {&self.damping}
    /// Certified modal damping (IRP M0.2-W5) of this image. It is certified on
    /// the full matrix, and every Diagonal block reports the one diagonal
    /// owner's zeta (`DiagonalPlateViscosity`), with its radius grown to keep
    /// the exact pole enclosed.
    ///
    /// # Errors
    /// Any `fs_modal::certify_modal_damping_owned_diagonal` refusal, or a
    /// `DampingInput` when this image carries no diagonal zeta.
    pub fn certified_modal_damping(&self, budget:&fs_modal::ModalDampingBudget)
        -> Result<fs_modal::CertifiedModalDamping,fs_modal::ModalError>
    {
        let zeta = self.damping.diagonal_zeta().ok_or(fs_modal::ModalError::DampingInput {
            field:"plate image without its diagonal zeta", index:0 })?;
        let omega:Vec<f64> = self.radiators.iter().map(|r|r.omega).collect();
        fs_modal::certify_modal_damping_owned_diagonal(&omega,self.damping.matrix_per_s(),zeta,budget)
    }
}
