//! Objective material loss on the original total-Lagrangian membrane.
//!
//! This is an explicit reference-volume Kelvin--Voigt extension: second-Piola
//! viscous stress is conjugate to engineering Green-strain rate. The existing
//! reduced material tensor supplies its coefficients, in the original
//! warp/weft axes. This objective extension does not turn a small-strain source
//! into a calibrated finite-strain law. Source validity and the membrane's
//! original geometric/strain envelope both remain necessary.
//!
//! Reference area times reference thickness integrates power, with the original
//! elastic owner's natural-volume pullback for shifted references. Thermal
//! references require the explicit same-claim endpoint constructor.
//! Neither current area, inertia, geometric prestress nor an invented decay rate
//! is damped.
//! In particular, a flat prestretched sheet has ZERO linear transverse loss
//! from this mechanism; curved geometry or in-plane strain rate enables loss.

use super::vibration::MembraneModes;
use super::{Element, MembraneError, MembranePatch, dot, finite_vectors, poll};
use fs_material::plane_stress_viscosity::{
    MaterialPlane, PlaneStressViscosity, factor_passive_resultant,
};
use fs_material::sheet::{MembraneLaw, THICKNESS};

/// Original membrane/source refusal or an incompatible material binding.
#[derive(Debug)]
pub enum MembraneViscosityError {
    /// Original geometry, material, work-limit or cancellation refusal.
    Membrane(MembraneError),
    /// Original dissipative constitutive refusal with element attribution.
    Material {
        /// Original triangle index.
        element: usize,
        /// Original material-law error.
        source: fs_material::MaterialError,
    },
    /// Original pinned-claim or endpoint material query refusal.
    Source {
        /// Original triangle index.
        element: usize,
        /// Original material-state error, without source substitution.
        source: fs_material::state_point::MaterialStatePointError,
    },
    /// Source identity, axes, band or finite-value admission failed.
    Invalid(&'static str),
}

impl From<MembraneError> for MembraneViscosityError {
    fn from(error: MembraneError) -> Self {
        Self::Membrane(error)
    }
}

impl core::fmt::Display for MembraneViscosityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "membrane viscosity: {self:?}")
    }
}

impl std::error::Error for MembraneViscosityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Membrane(error) => Some(error),
            Self::Material { source, .. } => Some(source),
            Self::Source { source, .. } => Some(source),
            Self::Invalid(_) => None,
        }
    }
}

type Result<T> = core::result::Result<T, MembraneViscosityError>;

fn invalid(what: &'static str) -> MembraneViscosityError {
    MembraneViscosityError::Invalid(what)
}

/// Physical full-space damping and independently integrated irreversible power.
#[derive(Clone, Debug)]
pub struct MembraneViscousResponse {
    /// The resisting term C(u)v [N] on the LEFT of the momentum equation.
    /// Actual viscous nodal force is its negative; v dot C(u)v equals power.
    pub resisting_force_n: Vec<[f64; 3]>,
    /// Nonnegative irreversible power [W], without numerical integration error.
    pub power_w: f64,
    /// Same power by original reference triangle [W], with no nodal heat split.
    pub element_power_w: Vec<f64>,
}

/// Frozen original sources and their physical strain-rate geometry owner.
#[derive(Clone, Debug)]
pub struct MembraneViscosity<'a> {
    patch: &'a MembranePatch,
    laws: Vec<PlaneStressViscosity>,
    reference_laws: Option<Vec<PlaneStressViscosity>>,
    volumes_m3: Vec<f64>,
    factors: Vec<[[f64; 6]; 3]>,
    natural_references: Vec<Option<NaturalReferencePullback>>,
    resolved_band_rad_s: (f64, f64),
}

// The elastic owner uses Fe = F G^-1 and W_mesh = det(G) W_natural.
// Engineering Green rates therefore transform by D = diag(g1^-2,
// g2^-2, (g1*g2)^-1), and their conjugate stresses by det(G) D.
// None preserves the original arithmetic for an unchanged reference.
#[derive(Clone, Copy, Debug)]
struct NaturalReferencePullback {
    rate: [f64; 3],
    stress: [f64; 3],
    volume_ratio: f64,
}

impl NaturalReferencePullback {
    fn from_stretches(stretches: [f64; 3]) -> Result<Option<Self>> {
        if stretches == [1.0; 3] {
            return Ok(None);
        }
        let inverse = stretches.map(f64::recip);
        let rate = [
            inverse[0] * inverse[0],
            inverse[1] * inverse[1],
            inverse[0] * inverse[1],
        ];
        let volume_ratio = stretches[0] * stretches[1] * stretches[2];
        let stress = rate.map(|scale| volume_ratio * scale);
        if rate.iter().chain(&stress).chain([&volume_ratio])
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(invalid("unrepresentable natural-reference viscosity pullback"));
        }
        Ok(Some(Self { rate, stress, volume_ratio }))
    }
}

// A resolved nonzero direction must not silently disappear during a change
// of material reference. Exact source/geometric null directions stay zero.
fn pullback_scale(value: f64, scale: f64) -> Result<f64> {
    let transformed = value * scale;
    if !transformed.is_finite() || (value != 0.0 && transformed == 0.0) {
        return Err(invalid("unrepresentable natural-reference viscous value"));
    }
    Ok(transformed)
}

/// Complete modal material resistance, with source rows proving its PSD form.
/// No modes or off-diagonal couplings are discarded.
#[derive(Clone, Debug)]
pub struct MembraneModalViscosity {
    matrix_per_s: Vec<f64>,
    dissipation_rows: Vec<f64>,
    mode_count: usize,
    row_count: usize,
}

impl MembraneModalViscosity {
    /// Row-major full C for mass-normalized coordinates, in inverse seconds.
    pub fn matrix_per_s(&self) -> &[f64] {
        &self.matrix_per_s
    }

    /// Row-major S, row_count by mode_count, with C = S transpose S.
    /// Each original facet contributes six material-factor rows in fixed order.
    /// Their squared action on modal velocity is physical power [W].
    pub fn dissipation_rows(&self) -> &[f64] {
        &self.dissipation_rows
    }

    /// Number of retained source rows, including exact zero rows.
    pub const fn row_count(&self) -> usize {
        self.row_count
    }

    /// Complete original mode count.
    pub const fn mode_count(&self) -> usize {
        self.mode_count
    }
}

// Difference-based gradients annihilate rigid translation algebraically.
fn rate(element: &Element, f: [[f64; 3]; 2], velocity: [[f64; 3]; 3]) -> [f64; 3] {
    let mut fdot = [[0.0; 3]; 2];
    for a in 1..3 {
        for axis in 0..3 {
            let difference = velocity[a][axis] - velocity[0][axis];
            for (alpha, column) in fdot.iter_mut().enumerate() {
                column[axis] += difference * element.gradients[a][alpha];
            }
        }
    }
    [
        dot(f[0], fdot[0]),
        dot(f[1], fdot[1]),
        dot(f[0], fdot[1]) + dot(f[1], fdot[0]),
    ]
}

impl<'a> MembraneViscosity<'a> {
    /// Bind one source tensor per original triangle, without inventing losses.
    ///
    /// Resolve sheet and viscosity properties together in ONE complete material
    /// state bundle, then derive both owners from that bundle. Exact equality
    /// includes the query, source selections and fixed environment. Axes12
    /// means the sheet's ordered warp/weft axes; the element already rotates
    /// its gradients into this frame, so viscosity is not rotated a second time.
    ///
    /// This constructor explicitly selects reference Green-rate Kelvin--Voigt
    /// damping. A shifted nonthermal hyperelastic reference uses Fe = F G^-1
    /// and det(G) times natural-volume power, exactly as the selected elastic
    /// energy. The original source stress and passive factors are pulled back
    /// to mesh coordinates without changing the source tensor or its nulls.
    /// Thermal references require [`Self::new_with_thermal_endpoints`], which
    /// refreshes only the pinned viscosity claims at the elastic endpoint.
    /// Reference thickness comes from the original sheet, independently of mass.
    ///
    /// The finite ordered resolved band must lie inside EVERY participating
    /// viscosity, areal-density, thickness and active membrane-law source band.
    /// Unused properties retained in the complete bundle do not constrain it.
    /// A request outside their intersection refuses; it is never silently clipped.
    /// This is a declared applicability interval, not a signal bandwidth proof.
    /// Transient consumers must separately admit their forcing and sampling.
    pub fn new(
        patch: &'a MembranePatch,
        laws: Vec<PlaneStressViscosity>,
        resolved_band_rad_s: (f64, f64),
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::assemble(patch, laws, None, resolved_band_rad_s, cancelled)
    }

    // Endpoint laws can enter only through the private thermal source binder;
    // their complete reference laws retain the original exact sheet binding.
    fn assemble(
        patch: &'a MembranePatch,
        laws: Vec<PlaneStressViscosity>,
        reference_laws: Option<Vec<PlaneStressViscosity>>,
        resolved_band_rad_s: (f64, f64),
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let (lo, hi) = resolved_band_rad_s;
        if !lo.is_finite() || !hi.is_finite() || lo < 0.0 || hi < lo {
            return Err(invalid("finite ordered membrane viscosity band required"));
        }
        if laws.len() != patch.elements.len()
            || reference_laws.as_ref().is_some_and(|sources| sources.len() != laws.len())
        {
            return Err(invalid("one source viscosity per original triangle required"));
        }
        let mut volumes_m3 = Vec::with_capacity(laws.len());
        let mut factors = Vec::with_capacity(laws.len());
        let mut natural_references = Vec::with_capacity(laws.len());
        let mut source_band = patch.mechanical_frequency_band_rad_s(&mut cancelled)?
            .unwrap_or((0.0_f64, f64::INFINITY));
        for (index, ((element, sheet), law)) in
            patch.elements.iter().zip(&patch.sheets).zip(&laws).enumerate()
        {
            poll(&mut cancelled)?;
            let reference_law = reference_laws.as_ref().map_or(law, |sources| &sources[index]);
            if reference_law.source() != sheet.material()
                || reference_law.plane() != MaterialPlane::Axes12
                || law.plane() != MaterialPlane::Axes12
            {
                return Err(invalid("viscosity needs the identical sheet source and warp/weft axes"));
            }
            let natural_reference = if let MembraneLaw::Hyperelastic(elastic) = &element.law {
                if elastic.mean_thermal_expansion().is_some() && reference_laws.is_none() {
                    return Err(invalid("thermal membrane references need endpoint viscosity sources"));
                }
                NaturalReferencePullback::from_stretches(elastic.reference_stretches())?
            } else {
                None
            };
            for omega in [lo, hi] {
                law.validate_frequency(omega)
                    .map_err(|source| MembraneViscosityError::Material { element: index, source })?;
            }
            let (loss_lo, loss_hi) = law.omega_band_rad_s();
            source_band = (source_band.0.max(loss_lo), source_band.1.min(loss_hi));
            // Reference-volume dissipation also consumes thickness, including
            // for membrane elastic laws whose energy is already per area.
            let thickness = sheet.material().property(THICKNESS)
                .ok_or_else(|| invalid("missing original membrane thickness source"))?;
            if let Some((low, high)) = thickness.angular_frequency_band_rad_s()
                .map_err(|source| MembraneViscosityError::Source { element: index, source })?
            {
                source_band = (source_band.0.max(low), source_band.1.min(high));
            }
            if source_band.1 < source_band.0 || lo < source_band.0 || hi > source_band.1 {
                return Err(invalid("declared membrane band exceeds combined mechanical and viscosity source support"));
            }
            let volume = element.area * sheet.thickness_m();
            if !volume.is_finite() || volume <= 0.0 {
                return Err(invalid("unrepresentable reference membrane volume"));
            }
            let mut factor = factor_passive_resultant(law.matrix_pa_s())
                .map_err(|source| MembraneViscosityError::Material { element: index, source })?;
            if let Some(reference) = natural_reference {
                let root_volume = reference.volume_ratio.sqrt();
                for (row, scale) in factor.iter_mut().zip(reference.rate) {
                    let scale = root_volume * scale;
                    for value in row {
                        *value = pullback_scale(*value, scale)?;
                    }
                }
            }
            volumes_m3.push(volume);
            factors.push(factor);
            natural_references.push(natural_reference);
        }
        poll(&mut cancelled)?;
        Ok(Self { patch, laws, reference_laws, volumes_m3, factors, natural_references, resolved_band_rad_s })
    }

    /// Exact original geometry, material and reference-inertia owner.
    pub const fn patch(&self) -> &'a MembranePatch {
        self.patch
    }

    /// Actual immutable loss laws used by force, heat and modal projection.
    /// Thermal facets use the selected elastic endpoint context; coefficients
    /// remain frozen there. Laws stay in their original material axes.
    pub fn laws(&self) -> &[PlaneStressViscosity] {
        &self.laws
    }

    /// Exact original sheet-bundle laws, before any thermal endpoint query.
    /// Their mass/thickness context and every original usage receipt remain
    /// unchanged. For [`Self::new`], these are also the actual [`Self::laws`].
    pub fn reference_laws(&self) -> &[PlaneStressViscosity] {
        self.reference_laws.as_deref().unwrap_or(&self.laws)
    }

    /// Declared resolved interval inside all participating mechanical/loss bands [rad/s].
    pub const fn resolved_band_rad_s(&self) -> (f64, f64) {
        self.resolved_band_rad_s
    }

    /// Evaluate objective damping and its physical power at the original state.
    /// Outputs allocate O(nodes + elements); failure publishes no partial result.
    /// The declared frequency band is NOT inferred from one instantaneous state.
    pub fn evaluate(
        &self,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MembraneViscousResponse> {
        poll(&mut cancelled)?;
        let mut resisting_force_n = vec![[0.0; 3]; self.patch.mesh.nodes.len()];
        let mut element_power_w = vec![0.0; self.laws.len()];
        let power_w = self.evaluate_into(
            displacement_m, velocity_m_s, &mut resisting_force_n,
            &mut element_power_w, &mut cancelled,
        )?;
        Ok(MembraneViscousResponse { resisting_force_n, power_w, element_power_w })
    }

    // Only private trial arrays may be partially changed on refusal. Fixed-q
    // time/linear-system consumers can reuse these arrays without allocations.
    pub(super) fn evaluate_into(
        &self,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        force_n: &mut [[f64; 3]],
        element_power_w: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64> {
        poll(cancelled)?;
        let nodes = self.patch.mesh.nodes.len();
        if !finite_vectors(displacement_m, nodes) || !finite_vectors(velocity_m_s, nodes)
            || force_n.len() != nodes || element_power_w.len() != self.laws.len()
        {
            return Err(invalid("complete finite membrane state and work arrays required"));
        }
        force_n.fill([0.0; 3]);
        let mut power_w = 0.0;
        for (index, element) in self.patch.elements.iter().enumerate() {
            poll(cancelled)?;
            let kinematics = self.patch.kinematics(element, displacement_m);
            let f = kinematics.gradient;
            self.patch.response(index, kinematics)?;
            let local_velocity = element.nodes.map(|node| velocity_m_s[node]);
            let mut strain_rate = rate(element, f, local_velocity);
            if let Some(reference) = self.natural_references[index] {
                for (value, scale) in strain_rate.iter_mut().zip(reference.rate) {
                    *value = pullback_scale(*value, scale)?;
                }
            }
            let (mut stress, mut power) = self.laws[index].stress_and_power(strain_rate)
                .map_err(|source| MembraneViscosityError::Material { element: index, source })?;
            if let Some(reference) = self.natural_references[index] {
                for (value, scale) in stress.iter_mut().zip(reference.stress) {
                    *value = pullback_scale(*value, scale)?;
                }
                power = pullback_scale(power, reference.volume_ratio)?;
            }
            let [s11, s22, s12] = stress;
            let volume = self.volumes_m3[index];
            let heat = volume * power;
            if !heat.is_finite() || (power != 0.0 && heat == 0.0) {
                return Err(invalid("unrepresentable membrane dissipation"));
            }
            element_power_w[index] = heat;
            power_w += heat;
            for (a, &node) in element.nodes.iter().enumerate() {
                let [gx, gy] = element.gradients[a];
                for axis in 0..3 {
                    force_n[node][axis] += volume * (
                        (s11 * gx + s12 * gy) * f[0][axis]
                        + (s12 * gx + s22 * gy) * f[1][axis]
                    );
                }
            }
        }
        if !power_w.is_finite() || !finite_vectors(force_n, nodes) {
            return Err(invalid("unrepresentable membrane damping force or power"));
        }
        poll(cancelled)?;
        Ok(power_w)
    }

    /// Six physical source rows at one admitted held geometry, in local-node order.
    /// Their squared action on nodal velocity is the facet's dissipated power.
    /// The caller admits the full state once; local checks keep a complete cache
    /// linear in the number of facets instead of rescanning every node per facet.
    pub(super) fn element_dissipation_rows(
        &self,
        index: usize,
        displacement_m: &[[f64; 3]],
    ) -> Result<[[[f64; 3]; 3]; 6]> {
        let element = self.patch.elements.get(index)
            .ok_or_else(|| invalid("original membrane triangle index required"))?;
        if displacement_m.len() != self.patch.mesh.nodes.len()
            || element.nodes.iter().any(|&node|
                displacement_m[node].iter().any(|value| !value.is_finite()))
        {
            return Err(invalid("complete locally finite membrane displacement required"));
        }
        let kinematics = self.patch.kinematics(element, displacement_m);
        let f = kinematics.gradient;
        self.patch.response(index, kinematics)?;
        let scale = self.volumes_m3[index].sqrt();
        let factor = &self.factors[index];
        let mut rows = [[[0.0; 3]; 3]; 6];
        for (k, row) in rows.iter_mut().enumerate() {
            for a in 1..3 {
                let [gx, gy] = element.gradients[a];
                let along_first = factor[0][k] * gx + factor[2][k] * gy;
                let along_second = factor[2][k] * gx + factor[1][k] * gy;
                for axis in 0..3 {
                    let value = along_first * f[0][axis] + along_second * f[1][axis];
                    let scaled = scale * value;
                    if !scaled.is_finite() || (value != 0.0 && scaled == 0.0) {
                        return Err(invalid("unrepresentable membrane dissipation row"));
                    }
                    row[a][axis] = scaled;
                }
            }
            // Balance the original triangle algebraically, including roundoff
            // in its stored gradients; no net translational resistance is added.
            for axis in 0..3 {
                row[0][axis] = -(row[1][axis] + row[2][axis]);
                if !row[0][axis].is_finite() {
                    return Err(invalid("unrepresentable membrane dissipation row"));
                }
            }
        }
        Ok(rows)
    }

    /// Project full material resistance through the ORIGINAL complete basis.
    ///
    /// Geometry is retained by the modal owner; foreign patches refuse.
    /// Every retained eigenfrequency enclosure must lie in the declared source band.
    /// This does not qualify arbitrary transients, unresolved loading or tails.
    ///
    /// The source Gram rows preserve directional and geometric nullspaces and
    /// retain every cross-mode term. They reuse the material owner's passive
    /// factorization tolerance; no new positive floor or modal fit is introduced.
    /// `max_terms` bounds projection contractions/storage before allocation.
    /// Constitutive state admission retains the original law's own iteration
    /// controls; hyperelastic local solves are outside this term count.
    pub fn project_modes(
        &self,
        modes: &MembraneModes<'_>,
        max_terms: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MembraneModalViscosity> {
        poll(&mut cancelled)?;
        if !core::ptr::eq(self.patch, modes.patch()) {
            return Err(invalid("viscosity and modes must share the original membrane"));
        }
        let n = modes.modes().len();
        let full = self.patch.mesh.nodes.len().checked_mul(3);
        let rows = self.laws.len().checked_mul(6);
        let square = n.checked_mul(n);
        let entries = rows.and_then(|rows| rows.checked_mul(n));
        let terms = full.and_then(|full| {
            self.laws.len().checked_mul(n.checked_mul(96)?)
                .and_then(|v| v.checked_add(rows?.checked_mul(square?)?))
                .and_then(|v| v.checked_add(entries?))
                .and_then(|v| v.checked_add(full))
                .and_then(|v| v.checked_add(square?))
        });
        if terms.is_none_or(|work| work > max_terms) {
            return Err(MembraneError::Budget("membrane viscosity projection").into());
        }
        let (full, rows, square, entries) = (
            full.unwrap(), rows.unwrap(), square.unwrap(), entries.unwrap(),
        );
        let (lo, hi) = self.resolved_band_rad_s;
        for mode in modes.modes() {
            poll(&mut cancelled)?;
            for lambda in [mode.interval.0, mode.lambda, mode.interval.1] {
                let omega = lambda.sqrt();
                if !omega.is_finite() || omega < lo || omega > hi {
                    return Err(invalid("membrane frequency enclosure exceeds the declared material band"));
                }
            }
        }
        let mut index = vec![usize::MAX; full];
        for (position, &dof) in modes.free_dofs().iter().enumerate() {
            index[dof] = position;
        }
        let mut dissipation_rows = vec![0.0; entries];
        let displacement = modes.linearization_displacement_m();
        for (element_index, element) in self.patch.elements.iter().enumerate() {
            poll(&mut cancelled)?;
            let kinematics = self.patch.kinematics(element, displacement);
            let f = kinematics.gradient;
            self.patch.response(element_index, kinematics)?;
            let scale = self.volumes_m3[element_index].sqrt();
            for (column, mode) in modes.modes().iter().enumerate() {
                let local = element.nodes.map(|node| core::array::from_fn(|axis| {
                    let dof = index[3 * node + axis];
                    if dof == usize::MAX { 0.0 } else { mode.phi[dof] }
                }));
                let strain_rate = rate(element, f, local);
                if strain_rate.iter().any(|value| !value.is_finite()) {
                    return Err(invalid("unrepresentable modal membrane strain rate"));
                }
                for k in 0..6 {
                    let value: f64 = (0..3).map(|a|
                        self.factors[element_index][a][k] * strain_rate[a]).sum();
                    let scaled = scale * value;
                    if !scaled.is_finite() || (value != 0.0 && scaled == 0.0) {
                        return Err(invalid("unrepresentable modal membrane dissipation row"));
                    }
                    dissipation_rows[(6 * element_index + k) * n + column] = scaled;
                }
            }
        }
        let mut matrix_per_s = vec![0.0; square];
        for i in 0..n {
            for j in i..n {
                let mut value = 0.0;
                for row in 0..rows {
                    if row % 256 == 0 { poll(&mut cancelled)?; }
                    value += dissipation_rows[row * n + i] * dissipation_rows[row * n + j];
                }
                if !value.is_finite() {
                    return Err(invalid("unrepresentable modal membrane resistance"));
                }
                matrix_per_s[i * n + j] = value;
                matrix_per_s[j * n + i] = value;
            }
        }
        poll(&mut cancelled)?;
        Ok(MembraneModalViscosity { matrix_per_s, dissipation_rows, mode_count: n, row_count: rows })
    }
}

// Explicit pinned-claim endpoint construction; no thermal evolution or rebind.
mod thermal;

/// Exact coupled evolution of the complete source-derived modal resistance.
pub mod modal;
