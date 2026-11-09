//! Complete harmonic laminate response in the original MITC shell trial space.
//!
//! The original elastic director-bubble map is a FIXED Galerkin approximation.
//! Every frequency uses its same eight strain rows, physical mass, support axes
//! and separately retained numerical drilling/physical connection stiffness.
//! Actual full-tensor layers supply all membrane/bending/transverse-shear cross
//! terms. Neither the bubble coordinates nor the source shear profile is
//! re-equilibrated at frequency. No frequency-dependent mass, modal fit,
//! prestress, thermal strain or temperature shift is inferred.
//!
//! Sources are immutable material snapshots at their authored reference state;
//! callers must resolve any changed temperature/state before constructing them.
//! Native source frequency bands and full-cycle strain domains remain binding.

use super::strain::{PreparedElementStrain, ShellStrainPoint, integration_rule};
use super::{MitcShellModel, ShearShellSection, ShellMass, element};
use crate::PlateError;
use fs_exec::Cx;
use fs_material::laminate::LayerLaw;
use fs_material::laminate::harmonic::{
    HarmonicLaminate, HarmonicLaminateBudget, HarmonicLaminateError, PreparedHarmonicLaminate,
};
use fs_material::visco::tensor::NormalStressHarmonic;
use fs_math::{c64::C64, dd::Dd, det};
use fs_modal::harmonic::{
    HarmonicResponse, HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_sparse::{Coo, Csr};
use std::collections::BTreeMap;
use std::sync::Arc;

mod supports;
pub use supports::PreparedSupportedHarmonicMitc;

/// Native spatial, source-material, budget or shared harmonic refusal.
#[derive(Debug)]
pub enum MitcHarmonicError {
    Plate(PlateError),
    /// The named original facet failed source/fixed-trial or inertia binding.
    SourceMismatch {
        element: usize,
        what: &'static str,
    },
    /// The original material refusal, including its original layer index.
    Material {
        element: usize,
        source: HarmonicLaminateError,
    },
    Budget {
        quantity: &'static str,
        requested: u128,
        maximum: u128,
    },
    Harmonic(HarmonicResponseError),
}
impl core::fmt::Display for MitcHarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-PLATE-MITC-HARMONIC: {self:?}")
    }
}
impl std::error::Error for MitcHarmonicError {}
impl From<PlateError> for MitcHarmonicError {
    fn from(value: PlateError) -> Self {
        Self::Plate(value)
    }
}
impl From<HarmonicResponseError> for MitcHarmonicError {
    fn from(value: HarmonicResponseError) -> Self {
        Self::Harmonic(value)
    }
}
type Result<T> = core::result::Result<T, MitcHarmonicError>;

/// Cold preparation and per-frequency constitutive-work allowances.
/// The existing model and original Arc-owned material sources are borrowed.
/// Dense factor/solve storage is separately admitted by HarmonicResponseControls.
#[derive(Clone, Copy, Debug)]
pub struct MitcHarmonicBudget {
    pub maximum_elements: usize,
    /// Sum of original layer visits over facets, including repeated sections.
    pub maximum_layers: usize,
    /// Additional cached strain rows and source indices, or prepared materials.
    /// Logical numerical payload; allocator metadata and fixed stacks excluded.
    pub maximum_payload_bytes: usize,
    /// Conservative work score per preparation, excluding the shared modal LU.
    pub maximum_work: u64,
}

#[derive(Clone, Copy)]
struct Point {
    rows: [[f64; 18]; 8],
    measure: f64,
}
struct Facet {
    observation: PreparedElementStrain,
    points: Vec<Point>,
    coordinates: [Option<usize>; 18],
    material: usize,
    sin: f64,
    cos: f64,
}
struct Material {
    source: Arc<HarmonicLaminate>,
    first_element: usize,
}

/// Original native shell, complete actual layers and immutable strain maps.
///
/// Supply one source per original facet. Repeated Arc identities share their
/// expensive per-frequency material preparation, while each original facet
/// retains its own directors, orientation and support-coordinate map.
/// The reference shear-profile section must reproduce the original elastic
/// A/B/D/H exactly. Actual full-tensor stiffness may differ, but its original
/// I0/I1/I2 must match every retained reference-inertia corner of that facet.
/// The mass must be the native consistent choice, including its original gauge.
pub struct HarmonicMitcShell<'a> {
    model: &'a MitcShellModel,
    materials: Vec<Material>,
    facets: Vec<Facet>,
    budget: MitcHarmonicBudget,
}

fn checkpoint(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| PlateError::Cancelled.into())
}
fn invalid(what: &'static str) -> MitcHarmonicError {
    HarmonicResponseError::Invalid(what).into()
}
fn cap(quantity: &'static str, requested: u128, maximum: u128) -> Result<()> {
    if requested > maximum {
        Err(MitcHarmonicError::Budget {
            quantity,
            requested,
            maximum,
        })
    } else {
        Ok(())
    }
}
fn mismatch(element: usize, what: &'static str) -> MitcHarmonicError {
    MitcHarmonicError::SourceMismatch { element, what }
}
fn material_budget(budget: MitcHarmonicBudget) -> HarmonicLaminateBudget {
    HarmonicLaminateBudget {
        maximum_layers: budget.maximum_layers,
        maximum_payload_bytes: budget.maximum_payload_bytes,
        maximum_work: budget.maximum_work,
    }
}

// Exactly the work-conjugate engineering rotation used by ShearShellSection.
// A common facet-to-section rotation acts on membrane, gradient AND shear.
fn section_strain(point: ShellStrainPoint, sin: f64, cos: f64) -> [f64; 8] {
    let rotate = |[x, y, xy]: [f64; 3]| {
        [
            cos * cos * x + sin * sin * y + cos * sin * xy,
            sin * sin * x + cos * cos * y - cos * sin * xy,
            -2.0 * cos * sin * x + 2.0 * cos * sin * y + (cos * cos - sin * sin) * xy,
        ]
    };
    let membrane = rotate(point.membrane);
    let gradient = rotate(point.gradient_per_m);
    let [xz, yz] = point.transverse_shear;
    [
        membrane[0],
        membrane[1],
        membrane[2],
        gradient[0],
        gradient[1],
        gradient[2],
        cos * xz + sin * yz,
        -sin * xz + cos * yz,
    ]
}

fn admit_binding(model: &MitcShellModel, element: usize, source: &HarmonicLaminate) -> Result<()> {
    let actual = &model.sections[element];
    let expected = ShearShellSection::from_laminate(
        source.shear_trial().section(),
        actual.material_angle_rad,
    )?;
    if actual.coupled_matrices() != expected.coupled_matrices() {
        return Err(mismatch(
            element,
            "original elastic section, reference surface and shear trial must agree",
        ));
    }
    let expected_inertia = source.inertia_moments();
    let corners = model
        .reference_inertia_moments(element)
        .ok_or_else(|| mismatch(element, "original reference inertia is missing"))?;
    if corners.iter().any(|corner| {
        [
            corner.areal_mass_kg_m2(),
            corner.mass_first_moment_kg_m(),
            corner.rotary_inertia_kg(),
        ] != expected_inertia
    }) {
        return Err(mismatch(
            element,
            "actual laminate density moments must match every original reference-inertia corner",
        ));
    }
    Ok(())
}

impl<'a> HarmonicMitcShell<'a> {
    /// Bind original native geometry/supports/mass to complete actual sources.
    /// No global mesh clone, material replacement or factorization occurs.
    pub fn try_new(
        model: &'a MitcShellModel,
        sources: Vec<Arc<HarmonicLaminate>>,
        budget: MitcHarmonicBudget,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        checkpoint(cx)?;
        let count = model.mesh.tris.len();
        if sources.len() != count {
            return Err(PlateError::SectionCount {
                expected: count,
                actual: sources.len(),
            }
            .into());
        }
        if model.mass_model != ShellMass::Consistent {
            return Err(invalid(
                "harmonic laminate MITC requires the original consistent mass",
            ));
        }
        cap(
            "MITC harmonic facets",
            count as u128,
            budget.maximum_elements as u128,
        )?;
        let mut layers = 0_u128;
        for source in &sources {
            checkpoint(cx)?;
            layers += source.section().layers().len() as u128;
        }
        cap(
            "MITC harmonic layer visits",
            layers,
            budget.maximum_layers as u128,
        )?;
        // Seven native integration points plus the six original Bernstein maps.
        // The extra 128 bytes/facet covers temporary identity-map/index entries.
        let bytes = count as u128
            * (core::mem::size_of::<Facet>()
                + 7 * core::mem::size_of::<Point>()
                + core::mem::size_of::<Material>()
                + 128) as u128;
        cap(
            "MITC harmonic cached payload bytes",
            bytes,
            budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "MITC harmonic preparation work",
            2_000_000 * count as u128,
            budget.maximum_work as u128,
        )?;
        let rule = integration_rule();
        if rule.len() != 7 {
            return Err(invalid(
                "original MITC degree-five integration rule required",
            ));
        }
        let mut identities = BTreeMap::new();
        let mut materials = Vec::new();
        let mut facets = Vec::with_capacity(count);
        for (element, source) in sources.into_iter().enumerate() {
            checkpoint(cx)?;
            admit_binding(model, element, &source)?;
            // Address identity only avoids repeated preparation of the SAME
            // immutable Arc; source output order remains original facet order.
            let identity = Arc::as_ptr(&source) as usize;
            let material = if let Some(&index) = identities.get(&identity) {
                index
            } else {
                let index = materials.len();
                identities.insert(identity, index);
                materials.push(Material {
                    source,
                    first_element: element,
                });
                index
            };
            let observation = PreparedElementStrain::new(
                &model.mesh,
                element,
                &model.directors,
                &model.sections[element],
                &mut || cx.checkpoint().is_err(),
            )?;
            let geometry = model.mesh.facet(element)?;
            let tri = model.mesh.tris[element];
            let angle = model.sections[element].material_angle_rad;
            let (sin, cos) = (det::sin(angle), det::cos(angle));
            let mut points: Vec<_> = rule
                .iter()
                .map(|(_, weight)| Point {
                    rows: [[0.0; 18]; 8],
                    measure: geometry.area_m2 * weight,
                })
                .collect();
            for column in 0..18 {
                checkpoint(cx)?;
                let mut basis = [0.0; 18];
                basis[column] = 1.0;
                let field = observation.real_field(&basis, &mut || cx.checkpoint().is_err())?;
                for (point, &(lambda, _)) in points.iter_mut().zip(&rule) {
                    let strain = section_strain(field.at_barycentric(lambda)?, sin, cos);
                    for component in 0..8 {
                        point.rows[component][column] = strain[component];
                    }
                }
            }
            // Original physical GLOBAL nodal rows -> original support axes.
            for point in &mut points {
                checkpoint(cx)?;
                let global = point.rows;
                for component in 0..8 {
                    for column in 0..18 {
                        let node = tri[column / 6];
                        let block = 6 * (column / 6) + 3 * ((column % 6) / 3);
                        let mut sum = Dd::ZERO;
                        for axis in 0..3 {
                            sum = sum
                                + Dd::from_f64(model.frames[node][column % 3][axis])
                                    * Dd::from_f64(global[component][block + axis]);
                        }
                        point.rows[component][column] = sum.to_f64();
                    }
                }
                if !point.measure.is_finite()
                    || point.measure <= 0.0
                    || point.rows.iter().flatten().any(|v| !v.is_finite())
                {
                    return Err(invalid(
                        "native MITC harmonic strain map is unrepresentable",
                    ));
                }
            }
            let coordinates =
                core::array::from_fn(|i| model.pencil.dof_map[6 * tri[i / 6] + i % 6]);
            facets.push(Facet {
                observation,
                points,
                coordinates,
                material,
                sin,
                cos,
            });
        }
        checkpoint(cx)?;
        Ok(Self {
            model,
            materials,
            facets,
            budget,
        })
    }

    pub const fn model(&self) -> &'a MitcShellModel {
        self.model
    }

    /// Complete original source associated with one original facet.
    pub fn laminate(&self, element: usize) -> Option<&HarmonicLaminate> {
        self.facets
            .get(element)
            .map(|facet| self.materials[facet.material].source.as_ref())
    }

    fn frequency_payload(
        &self,
        omega: f64,
        controls: HarmonicResponseControls,
        factor: bool,
        cx: &Cx<'_>,
    ) -> Result<usize> {
        checkpoint(cx)?;
        if !omega.is_finite() || omega <= 0.0 {
            return Err(invalid(
                "MITC harmonic frequency must be finite and positive",
            ));
        }
        let n = self.model.pencil.free;
        cap(
            "MITC harmonic coordinates",
            n as u128,
            controls.maximum_coordinates as u128,
        )?;
        if !controls.maximum_relative_residual.is_finite()
            || controls.maximum_relative_residual <= 0.0
            || controls.maximum_relative_residual >= 1.0
        {
            return Err(invalid(
                "MITC harmonic residual allowance must lie in (0,1)",
            ));
        }
        let mut bytes = self.materials.len() as u128
            * core::mem::size_of::<PreparedHarmonicLaminate<'_>>() as u128;
        let mut work = 2_000_000 * self.facets.len() as u128;
        for material in &self.materials {
            checkpoint(cx)?;
            for layer in material.source.section().layers() {
                bytes += core::mem::size_of::<NormalStressHarmonic<'_>>() as u128;
                work += 500_000;
                if let LayerLaw::Prony(law) = &layer.law {
                    work += 4096 * law.branches().len() as u128;
                }
            }
        }
        cap(
            "MITC harmonic constitutive payload bytes",
            bytes,
            self.budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "MITC harmonic frequency work",
            work,
            self.budget.maximum_work as u128,
        )?;
        // COO capacity growth, stable-sort workspace and both finished CSR
        // arrays are bounded before staging. Borrowed native matrices excluded.
        let triplets = 648 * self.facets.len() as u128 + self.model.fixed_stiffness.nnz() as u128;
        let sparse_peak = 128 * triplets + 32 * (n as u128 + 1) + bytes;
        cap(
            "MITC harmonic operator payload bytes",
            sparse_peak,
            controls.maximum_payload_bytes as u128,
        )?;
        if factor {
            // Three source dense arrays coexist with the shared factor owner.
            let dense = n
                .checked_mul(n)
                .and_then(|v| v.checked_mul(168))
                .and_then(|v| n.checked_mul(128).and_then(|w| v.checked_add(w)))
                .ok_or_else(|| invalid("MITC harmonic dense payload extent overflows"))?;
            cap(
                "MITC harmonic factor payload bytes",
                dense as u128 + bytes,
                controls.maximum_payload_bytes as u128,
            )?;
        }
        usize::try_from(bytes).map_err(|_| invalid("MITC harmonic material payload overflows"))
    }

    /// Full sparse storage and viscosity at one original material frequency.
    /// Bare resonances do not require a factor here. A coupled owner must call
    /// admit_displacement before publishing its solved physical response.
    pub fn prepare_harmonic_operator(
        &self,
        omega: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedHarmonicMitcOperator<'_, 'a>> {
        self.frequency_payload(omega, controls, false, cx)?;
        let mut materials = Vec::with_capacity(self.materials.len());
        for material in &self.materials {
            checkpoint(cx)?;
            materials.push(
                material
                    .source
                    .prepare_harmonic(omega, material_budget(self.budget), || {
                        cx.checkpoint().is_err()
                    })
                    .map_err(|source| MitcHarmonicError::Material {
                        element: material.first_element,
                        source,
                    })?,
            );
        }
        let n = self.model.pencil.free;
        let mut storage = Coo::new(n, n);
        let mut damping = Coo::new(n, n);
        for row in 0..n {
            checkpoint(cx)?;
            let (columns, values) = self.model.fixed_stiffness.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                if value != 0.0 {
                    storage.push(row, column, value);
                }
            }
        }
        for facet in &self.facets {
            checkpoint(cx)?;
            let material = &materials[facet.material];
            let mut stiffness = [[Dd::ZERO; 18]; 18];
            let mut loss = [[Dd::ZERO; 18]; 18];
            for point in &facet.points {
                element::full_constitutive_work(
                    &point.rows,
                    material.storage(),
                    point.measure,
                    &mut stiffness,
                    &mut || cx.checkpoint().is_err(),
                )?;
                element::full_constitutive_work(
                    &point.rows,
                    material.loss(),
                    point.measure,
                    &mut loss,
                    &mut || cx.checkpoint().is_err(),
                )?;
            }
            for i in 0..18 {
                checkpoint(cx)?;
                let Some(row) = facet.coordinates[i] else {
                    continue;
                };
                for j in i..18 {
                    let Some(column) = facet.coordinates[j] else {
                        continue;
                    };
                    let k = stiffness[i][j].to_f64();
                    let original_loss = loss[i][j].to_f64();
                    let c = original_loss / omega;
                    if !c.is_finite() || (original_loss != 0.0 && c == 0.0) {
                        return Err(invalid(
                            "MITC source loss divided by frequency is unrepresentable",
                        ));
                    }
                    for (matrix, value) in [(&mut storage, k), (&mut damping, c)] {
                        if value != 0.0 {
                            matrix.push(row, column, value);
                            if row != column {
                                matrix.push(column, row, value);
                            }
                        }
                    }
                }
            }
        }
        checkpoint(cx)?;
        let stiffness = storage.assemble();
        let damping = damping.assemble();
        for row in 0..n {
            checkpoint(cx)?;
            if stiffness
                .row(row)
                .1
                .iter()
                .chain(damping.row(row).1)
                .any(|v| !v.is_finite())
            {
                return Err(invalid("assembled full MITC harmonic operator overflows"));
            }
        }
        Ok(PreparedHarmonicMitcOperator {
            source: self,
            omega,
            materials,
            stiffness,
            damping,
        })
    }

    /// Factor the complete original physical pencil K'(omega)-omega²M-iK''.
    /// Repeated drives reuse the same native fs-modal factor and material maps.
    pub fn prepare_harmonic_response(
        &self,
        omega: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedHarmonicMitc<'_, 'a>> {
        self.frequency_payload(omega, controls, true, cx)?;
        self.prepare_harmonic_operator(omega, controls, cx)?
            .factor(controls, cx)
    }
}

/// A complete original free-coordinate material pencil, without a bare solve.
pub struct PreparedHarmonicMitcOperator<'a, 'model> {
    source: &'a HarmonicMitcShell<'model>,
    omega: f64,
    materials: Vec<PreparedHarmonicLaminate<'a>>,
    stiffness: Csr,
    damping: Csr,
}
impl<'a, 'model> PreparedHarmonicMitcOperator<'a, 'model> {
    pub const fn source(&self) -> &'a HarmonicMitcShell<'model> {
        self.source
    }
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }
    /// The original native consistent mass, including its original gauge.
    pub fn mass(&self) -> &Csr {
        &self.source.model.pencil.m
    }
    /// Full actual material storage plus original lossless drilling/connections.
    pub const fn stiffness(&self) -> &Csr {
        &self.stiffness
    }
    /// Complete material loss/omega, including signed cross-coupling entries.
    pub const fn damping(&self) -> &Csr {
        &self.damping
    }

    /// Certify the entire native facet polynomial, every original layer's
    /// thickness and every common harmonic phase. Six in-plane Bernstein
    /// controls feed the material owner's own thickness/cycle admission.
    /// This sufficient convex bound may conservatively refuse valid motion.
    pub fn admit_displacement(&self, displacement: &[C64], cx: &Cx<'_>) -> Result<()> {
        admit_displacement(self.source, &self.materials, displacement, cx)
    }

    /// Consume these source operators into the existing dense harmonic owner.
    /// No independent eigensolver, damping fit or changed mass is introduced.
    pub fn factor(
        self,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedHarmonicMitc<'a, 'model>> {
        let material_bytes = self
            .source
            .frequency_payload(self.omega, controls, true, cx)?;
        let n = self.source.model.pencil.free;
        let source_bytes = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(24))
            .and_then(|v| v.checked_add(material_bytes))
            .ok_or_else(|| invalid("MITC harmonic dense source extent overflows"))?;
        let mass = dense(&self.source.model.pencil.m, n, cx)?;
        let storage = dense(&self.stiffness, n, cx)?;
        let viscosity = dense(&self.damping, n, cx)?;
        let Self {
            source,
            omega,
            materials,
            stiffness,
            damping,
        } = self;
        drop(stiffness);
        drop(damping);
        let response = PreparedHarmonicResponse::try_new(
            &mass,
            &viscosity,
            &storage,
            n,
            omega,
            HarmonicResponseControls {
                maximum_payload_bytes: controls.maximum_payload_bytes - source_bytes,
                ..controls
            },
            cx,
        )?;
        checkpoint(cx)?;
        Ok(PreparedHarmonicMitc {
            source,
            omega,
            materials,
            response,
        })
    }
}

/// One native factor and the exact material/domain owners used to assemble it.
pub struct PreparedHarmonicMitc<'a, 'model> {
    source: &'a HarmonicMitcShell<'model>,
    omega: f64,
    materials: Vec<PreparedHarmonicLaminate<'a>>,
    response: PreparedHarmonicResponse,
}
impl<'a, 'model> PreparedHarmonicMitc<'a, 'model> {
    pub const fn source(&self) -> &'a HarmonicMitcShell<'model> {
        self.source
    }
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }

    /// Solve original free/support-coordinate force and moment phasors with
    /// exp(-i omega t). Full source-strain admission precedes publication.
    /// Shared average_viscous_power_w reports actual material loss; numerical
    /// drilling and fixed physical springs introduce no heat. This is a fixed
    /// material-state, linear forced response, not thermal/prestressed motion.
    pub fn solve(&self, force: &[C64], cx: &Cx<'_>) -> Result<HarmonicResponse> {
        let response = self.response.solve(force, cx)?;
        admit_displacement(self.source, &self.materials, &response.displacement, cx)?;
        Ok(response)
    }
}

fn admit_displacement(
    source: &HarmonicMitcShell<'_>,
    materials: &[PreparedHarmonicLaminate<'_>],
    displacement: &[C64],
    cx: &Cx<'_>,
) -> Result<()> {
    checkpoint(cx)?;
    if displacement.len() != source.model.pencil.free
        || displacement
            .iter()
            .any(|v| !v.re.is_finite() || !v.im.is_finite())
    {
        return Err(invalid(
            "complete finite original MITC displacement phasors required",
        ));
    }
    for (element, facet) in source.facets.iter().enumerate() {
        checkpoint(cx)?;
        let tri = source.model.mesh.tris[element];
        let mut real = [0.0; 18];
        let mut imaginary = [0.0; 18];
        for i in 0..18 {
            let node = tri[i / 6];
            let block = 3 * ((i % 6) / 3);
            let mut value = C64::ZERO;
            for axis in 0..3 {
                if let Some(index) = source.model.pencil.dof_map[6 * node + block + axis] {
                    value =
                        value + displacement[index].scale(source.model.frames[node][axis][i % 3]);
                }
            }
            real[i] = value.re;
            imaginary[i] = value.im;
        }
        let field = facet
            .observation
            .harmonic_field(&real, &imaginary, &mut || cx.checkpoint().is_err())?;
        for control in field.bernstein_controls() {
            let re = section_strain(control.real, facet.sin, facet.cos);
            let im = section_strain(control.imaginary, facet.sin, facet.cos);
            let generalized = core::array::from_fn(|i| C64::new(re[i], im[i]));
            materials[facet.material]
                .admit_generalized_strain(&generalized, || cx.checkpoint().is_err())
                .map_err(|source| MitcHarmonicError::Material { element, source })?;
        }
    }
    checkpoint(cx)
}

fn dense(matrix: &Csr, n: usize, cx: &Cx<'_>) -> Result<Vec<f64>> {
    let mut output = vec![0.0; n * n];
    for row in 0..n {
        checkpoint(cx)?;
        let (columns, values) = matrix.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            output[n * row + column] = value;
        }
    }
    Ok(output)
}
