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

use crate::acoustic_realize::AcousticRealizeError;
use crate::thin_plate::{CompactBody, PlateMaterialModel, ResolvedPlateChart};
use crate::render::schedule::force::coupled::{ModalAttachment, ModalConnection};
use fs_exec::CancelGate;
use fs_material::plane_stress_viscosity::{MaterialPlane, PlaneStressViscosity};
use fs_blake3::ContentHash;
use fs_plate::PlateModel;

/// Source-bound hereditary bending on the same original plate basis.
pub mod relaxation;

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn owner(error: impl core::fmt::Display) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(error.to_string())
}
fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() { Err(AcousticRealizeError::Cancelled) } else { Ok(()) }
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
    matrix:Vec<f64>, columns:Vec<Vec<f64>>, resistance:f64, modes:usize,
}
impl ModalPlateViscosity {
    pub(crate) fn factor(matrix:Vec<f64>,n:usize,resistance:f64,gate:&CancelGate)->Result<Self,AcousticRealizeError> {
        if matrix.iter().any(|v|!v.is_finite()) { return Err(invalid("projected material damping is not finite")); }
        if matrix.iter().all(|v|*v == 0.0) {
            return Ok(Self {matrix,columns:Vec::new(),resistance,modes:n});
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
        Ok(Self {matrix,columns,resistance,modes:n})
    }
    /// Dense symmetric C [1/s], row-major in the original mass-normalized basis.
    #[must_use] pub fn matrix_per_s(&self)->&[f64] { &self.matrix }
    /// Mode count, including modes with explicit zero dissipation.
    #[must_use] pub const fn mode_count(&self)->usize {self.modes}
    /// Number of passive channels. Zero tensors allocate no artificial damping.
    #[must_use] pub fn channel_count(&self)->usize {self.columns.len()}
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
}
