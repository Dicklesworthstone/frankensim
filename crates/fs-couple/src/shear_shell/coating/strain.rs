//! Source-bound local film-strain envelopes on an unchanged physical mode basis.
//!
//! The original MITC owner condenses each coated element once. Each query first
//! sums the complex modal displacement coherently, then bounds the recovered
//! tensor over every coated triangle and both exact film faces. These are
//! Bernstein convex-hull bounds on the modeled field, not quadrature samples,
//! continuum discretization certificates or material-failure predictions.

use super::ResolvedShellCoating;
use crate::acoustic_realize::AcousticRealizeError;
use crate::shear_shell::{
    MaterialShellImage, ResolvedShearShell, assembly::MaterialShellAssemblyImage, invalid, owner,
    poll,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_plate::shell::mitc::{
    NodalDirectors,
    strain::{HarmonicElementStrain, PreparedElementStrain, RealElementStrain, TensorStrainBound},
};
use std::mem::size_of;

/// Explicit preparation, storage and complete planned query allowance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ShellCoatingStrainBudget {
    /// Prepared storage plus one returned report and this adapter's element-local
    /// query arrays. Original borrowed source/mode storage, allocator bookkeeping,
    /// the original element owner's fixed-size stack frames and caller-retained
    /// older reports are excluded. No source film, element or retained mode is
    /// omitted to fit.
    pub maximum_bytes: usize,
    /// Preparation plus `evaluations` complete queries, in bounded scalar work.
    pub maximum_work: usize,
    /// Positive number of successful evaluations admitted before preparation.
    pub evaluations: usize,
}

/// Admitted resource counts for this exact source, mode basis and query count.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ShellCoatingStrainResources {
    pub retained_bytes: usize,
    pub evaluation_bytes: usize,
    pub preparation_work: usize,
    pub per_evaluation_work: usize,
    pub total_work: usize,
    pub evaluations: usize,
    pub coated_elements: usize,
    pub film_element_pairs: usize,
}

/// Original film identity; component and film order are never sorted or merged.
#[derive(Clone, Copy, Debug)]
pub struct ShellCoatingFilmSource<'a> {
    /// `None` for a single shell, otherwise the original assembly component.
    pub component_index: Option<usize>,
    pub component_name: Option<&'a str>,
    pub film_index: usize,
    /// Exact source state, region, side, authored strain limit and film faces.
    pub source: &'a ResolvedShellCoating,
}

/// Whole-film amplitude envelope for supplied peak complex body motion, or
/// absolute spatial envelope for an explicitly observed real snapshot.
#[derive(Clone, Debug)]
pub struct ShellCoatingFilmStrain<'a> {
    pub component_index: Option<usize>,
    pub component_name: Option<&'a str>,
    pub film_index: usize,
    pub source: &'a ResolvedShellCoating,
    pub covered_elements: usize,
    /// Upper bound on the peak full-tensor Frobenius norm over the complete film.
    /// Plane-stress epsilon_zz uses the film's authored real Poisson ratio.
    /// Engineering shear becomes tensor shear; no energy shear correction is
    /// applied to strain. This is the original equivalent-single-layer MITC
    /// assumed shear, not a recovered 3D interlaminar stress/strain field.
    pub frobenius_upper_bound: f64,
    /// Component amplitude envelopes [xx, yy, zz, xy, yz, zx] in facet axes.
    /// These separately maximized values need not occur at a common phase/site.
    pub component_amplitudes: [f64; 6],
    /// Original triangle owning the ruling Bernstein control bound. A control
    /// coefficient is not generally a value at a physical maximum location.
    pub ruling_element: usize,
    pub ruling_control_point: usize,
    /// 0/1 selects the exact bottom/top face from `source.faces_m()`.
    pub ruling_z_face: usize,
    /// Bound divided by the source's authored small-strain envelope. Above one
    /// means this sufficient bound cannot establish envelope compliance; it
    /// does not establish actual overstrain, damage or failure.
    pub limit_ratio: f64,
    /// Sufficient nonnegative common multiplier of the complete complex motion.
    /// `None` means the response is exactly zero, with no finite strain cap.
    pub sufficient_amplitude_scale: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct ShellCoatingStrainReport<'a> {
    pub films: Vec<ShellCoatingFilmStrain<'a>>,
    /// Every whole-film bound lies within its authored small-strain envelope.
    /// No unmodeled preload, cure strain or other static strain is included.
    pub within_authored_envelope: bool,
    /// Minimum sufficient scale over all films, preserving every modal phase.
    /// `None` means there are no films or every observed film strain is zero.
    pub sufficient_common_amplitude_scale: Option<f64>,
}

/// Numeric film envelope, independent of a report borrow. A prepared transient
/// observer reuses this storage at every accepted endpoint.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ShellCoatingFilmStrainValue {
    pub covered_elements: usize,
    pub frobenius_upper_bound: f64,
    pub component_amplitudes: [f64; 6],
    pub ruling_element: usize,
    pub ruling_control_point: usize,
    pub ruling_z_face: usize,
    pub limit_ratio: f64,
    pub sufficient_amplitude_scale: Option<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
/// Common applicability result for one complete instantaneous or harmonic query.
pub struct ShellCoatingStrainSummary {
    pub within_authored_envelope: bool,
    pub sufficient_common_amplitude_scale: Option<f64>,
}

struct Element {
    index: usize,
    region: usize,
    nodes: [usize; 3],
    field: PreparedElementStrain,
}

struct Part {
    // These immutable source allocations survive moves of their owning bank.
    // Only identity is compared; addresses are never dereferenced or exposed.
    source_storage: [usize; 4],
    film_regions: Vec<usize>,
    elements: Vec<Element>,
}

/// Prepared local observation tied to one original shell/assembly source and
/// its retained modes. No foreign directors, sections or material limits can
/// be substituted. Query errors and cancellation leave its allowance intact.
pub struct PreparedShellCoatingStrain<'a> {
    source: CoatingStrainSource<'a>,
    plan: OwnedShellCoatingStrain,
}

/// Internal detached plan, bound to the immutable original image allocations.
/// Keeping source borrows out of this type permits the original performance to
/// own it without cloning its shell, modes or constitutive histories.
pub(crate) struct OwnedShellCoatingStrain {
    parts: Vec<Part>,
    values: Vec<ShellCoatingFilmStrainValue>,
    resources: ShellCoatingStrainResources,
    modes: usize,
    films: usize,
    assembly: bool,
    snapshot: bool,
    completed_evaluations: usize,
}

#[derive(Clone, Copy)]
pub(crate) enum CoatingStrainSource<'a> {
    Shell(&'a MaterialShellImage),
    Assembly(&'a MaterialShellAssemblyImage),
}

struct SourcePart<'a> {
    index: Option<usize>,
    name: Option<&'a str>,
    shell: &'a ResolvedShearShell,
    directors: &'a NodalDirectors,
    modes: &'a [Vec<f64>],
}

impl SourcePart<'_> {
    fn storage_identity(&self) -> [usize; 4] {
        [
            self.modes.as_ptr() as usize,
            self.modes.len(),
            self.shell.mesh().nodes.as_ptr() as usize,
            self.shell.coatings().as_ptr() as usize,
        ]
    }
}

impl<'a> CoatingStrainSource<'a> {
    fn len(&self) -> usize {
        match self {
            Self::Shell(_) => 1,
            Self::Assembly(image) => image.source().components.len(),
        }
    }

    fn part(&self, i: usize) -> SourcePart<'a> {
        match self {
            Self::Shell(image) => SourcePart {
                index: None,
                name: None,
                shell: image.source(),
                directors: image.directors(),
                modes: image.physical_modes(),
            },
            Self::Assembly(image) => {
                let part = &image.source().components[i];
                SourcePart {
                    index: Some(i),
                    name: Some(&part.name),
                    shell: &part.specimen,
                    directors: &part.directors,
                    modes: &image.physical_modes()[i],
                }
            }
        }
    }

    pub(crate) fn film_sources(self) -> impl Iterator<Item = ShellCoatingFilmSource<'a>> {
        (0..self.len()).flat_map(move |i| {
            let part = self.part(i);
            part.shell
                .coatings()
                .iter()
                .enumerate()
                .map(move |(film_index, source)| ShellCoatingFilmSource {
                    component_index: part.index,
                    component_name: part.name,
                    film_index,
                    source,
                })
        })
    }

    pub(crate) fn prepare_snapshot(
        self,
        budget: ShellCoatingStrainBudget,
        gate: &CancelGate,
    ) -> Result<OwnedShellCoatingStrain, AcousticRealizeError> {
        prepare(self, budget, true, gate)
    }
}

impl MaterialShellImage {
    /// Prepare complete local film observation with the original modal basis,
    /// elastic sections and assembly directors. Valid for frozen or hereditary
    /// films; this observes supplied motion and does not solve a constitutive law.
    pub fn prepare_coating_strain(
        &self,
        budget: ShellCoatingStrainBudget,
        gate: &CancelGate,
    ) -> Result<PreparedShellCoatingStrain<'_>, AcousticRealizeError> {
        let source = CoatingStrainSource::Shell(self);
        Ok(PreparedShellCoatingStrain {
            source,
            plan: prepare(source, budget, false, gate)?,
        })
    }
}

impl MaterialShellAssemblyImage {
    /// Observe every component film in the SAME assembled modal coordinates.
    pub fn prepare_coating_strain(
        &self,
        budget: ShellCoatingStrainBudget,
        gate: &CancelGate,
    ) -> Result<PreparedShellCoatingStrain<'_>, AcousticRealizeError> {
        let source = CoatingStrainSource::Assembly(self);
        Ok(PreparedShellCoatingStrain {
            source,
            plan: prepare(source, budget, false, gate)?,
        })
    }
}

fn add_product(total: &mut usize, a: usize, b: usize) -> Result<(), AcousticRealizeError> {
    *total = a
        .checked_mul(b)
        .and_then(|n| total.checked_add(n))
        .ok_or_else(|| invalid("complete shell coating strain resources overflow"))?;
    Ok(())
}

fn film_region(
    shell: &ResolvedShearShell,
    film: &ResolvedShellCoating,
) -> Result<usize, AcousticRealizeError> {
    shell
        .regions()
        .iter()
        .position(|region| region.input().region.region.name == film.input().region)
        .ok_or_else(|| invalid("coating strain requires every original source film region"))
}

fn strain_owner(error: fs_plate::PlateError) -> AcousticRealizeError {
    match error {
        fs_plate::PlateError::Cancelled => AcousticRealizeError::Cancelled,
        error => owner(error),
    }
}

fn prepare(
    source: CoatingStrainSource<'_>,
    budget: ShellCoatingStrainBudget,
    snapshot: bool,
    gate: &CancelGate,
) -> Result<OwnedShellCoatingStrain, AcousticRealizeError> {
    poll(gate)?;
    if budget.evaluations == 0 || source.len() == 0 {
        return Err(invalid(
            "coating strain requires a source and positive planned evaluations",
        ));
    }
    let modes = source.part(0).modes.len();
    let mut films = 0usize;
    let mut resources = ShellCoatingStrainResources {
        retained_bytes: if snapshot {
            size_of::<OwnedShellCoatingStrain>()
        } else {
            size_of::<PreparedShellCoatingStrain<'_>>()
        },
        evaluation_bytes: size_of::<TensorStrainBound>()
            + if snapshot {
                RealElementStrain::RETAINED_BYTES + 18 * size_of::<f64>()
            } else {
                size_of::<ShellCoatingStrainReport<'_>>()
                    + HarmonicElementStrain::RETAINED_BYTES
                    + 36 * size_of::<f64>()
            },
        preparation_work: 1,
        per_evaluation_work: if snapshot { modes } else { 2 * modes },
        total_work: 0,
        evaluations: budget.evaluations,
        coated_elements: 0,
        film_element_pairs: 0,
    };
    // Count exact coverage without allocating an element field or result vector.
    // Regions partition the original mesh; several films may share one region.
    for i in 0..source.len() {
        poll(gate)?;
        let part = source.part(i);
        let shell = part.shell;
        if modes == 0
            || part.modes.len() != modes
            || part
                .modes
                .iter()
                .any(|m| m.len() != 6 * shell.mesh().nodes.len())
        {
            return Err(invalid(
                "coating strain requires the complete original physical mode basis",
            ));
        }
        let count = shell.coatings().len();
        add_product(&mut films, count, 1)?;
        add_product(&mut resources.retained_bytes, 1, size_of::<Part>())?;
        add_product(&mut resources.retained_bytes, count, size_of::<usize>())?;
        add_product(&mut resources.preparation_work, modes, 1)?;
        add_product(
            &mut resources.preparation_work,
            count,
            4 * shell.regions().len(),
        )?;
        add_product(
            &mut resources.preparation_work,
            shell.mesh().tris.len(),
            count + 1,
        )?;
        add_product(&mut resources.per_evaluation_work, count, 64)?;
        add_product(&mut resources.per_evaluation_work, 1, 32)?;
        for film in shell.coatings() {
            let region = film_region(shell, film)?;
            let pairs = shell.regions()[region]
                .input()
                .region
                .region
                .triangle_indices
                .len();
            add_product(&mut resources.film_element_pairs, pairs, 1)?;
            add_product(
                &mut resources.per_evaluation_work,
                pairs,
                if snapshot {
                    PreparedElementStrain::REAL_LAYER_BOUND_WORK + 32
                } else {
                    PreparedElementStrain::LAYER_BOUND_WORK + 32
                },
            )?;
        }
        for region in shell.regions() {
            if !shell
                .coatings()
                .iter()
                .any(|film| film.input().region == region.input().region.region.name)
            {
                continue;
            }
            let elements = region.input().region.region.triangle_indices.len();
            add_product(&mut resources.coated_elements, elements, 1)?;
            add_product(
                &mut resources.retained_bytes,
                elements,
                size_of::<Element>(),
            )?;
            add_product(
                &mut resources.preparation_work,
                elements,
                PreparedElementStrain::PREPARATION_WORK,
            )?;
            add_product(
                &mut resources.per_evaluation_work,
                elements,
                if snapshot {
                    36 * modes + PreparedElementStrain::REAL_FIELD_WORK + count
                } else {
                    72 * modes + PreparedElementStrain::HARMONIC_FIELD_WORK + count
                },
            )?;
        }
    }
    add_product(
        &mut resources.retained_bytes,
        films,
        size_of::<ShellCoatingFilmStrainValue>(),
    )?;
    if !snapshot {
        add_product(
            &mut resources.evaluation_bytes,
            films,
            size_of::<ShellCoatingFilmStrain<'_>>(),
        )?;
    }
    resources.total_work = resources.preparation_work;
    add_product(
        &mut resources.total_work,
        budget.evaluations,
        resources.per_evaluation_work,
    )?;
    let total_bytes = resources
        .retained_bytes
        .checked_add(resources.evaluation_bytes);
    if total_bytes.is_none_or(|bytes| bytes > budget.maximum_bytes)
        || resources.total_work > budget.maximum_work
    {
        return Err(invalid(
            "complete local coating strain preparation and planned queries exceed their declared allowance",
        ));
    }
    let mut parts = Vec::with_capacity(source.len());
    for i in 0..source.len() {
        poll(gate)?;
        let part = source.part(i);
        let mut film_regions = Vec::with_capacity(part.shell.coatings().len());
        for film in part.shell.coatings() {
            film_regions.push(film_region(part.shell, film)?);
        }
        let n = part
            .shell
            .regions()
            .iter()
            .enumerate()
            .filter(|(r, _)| film_regions.contains(r))
            .map(|(_, region)| region.input().region.region.triangle_indices.len())
            .sum();
        let mut elements = Vec::with_capacity(n);
        for (index, nodes) in part.shell.mesh().tris.iter().copied().enumerate() {
            if index % 64 == 0 {
                poll(gate)?;
            }
            let region = part.shell.element_regions[index];
            if film_regions.contains(&region) {
                let field = PreparedElementStrain::new(
                    part.shell.mesh(),
                    index,
                    part.directors,
                    &part.shell.sections()[index],
                    &mut || gate.is_requested(),
                )
                .map_err(strain_owner)?;
                elements.push(Element {
                    index,
                    region,
                    nodes,
                    field,
                });
            }
        }
        parts.push(Part {
            source_storage: part.storage_identity(),
            film_regions,
            elements,
        });
    }
    poll(gate)?;
    Ok(OwnedShellCoatingStrain {
        parts,
        values: vec![ShellCoatingFilmStrainValue::default(); films],
        resources,
        modes,
        films,
        assembly: matches!(source, CoatingStrainSource::Assembly(_)),
        snapshot,
        completed_evaluations: 0,
    })
}

impl<'a> PreparedShellCoatingStrain<'a> {
    pub const fn resources(&self) -> &ShellCoatingStrainResources {
        &self.plan.resources
    }
    pub const fn film_count(&self) -> usize {
        self.plan.films
    }
    pub const fn completed_evaluations(&self) -> usize {
        self.plan.completed_evaluations
    }

    pub fn film_sources(&self) -> impl Iterator<Item = ShellCoatingFilmSource<'a>> + '_ {
        self.source.film_sources()
    }

    /// Observe one supplied peak phasor in the original mass-normalized body
    /// coordinates [m sqrt(kg)]. The coherent sum precedes strain norms; modal
    /// amplitudes or separate-force bounds must never be added as if in phase.
    /// Successful reports consume one admitted query. The physical image and
    /// accepted constitutive histories are never modified.
    pub fn evaluate(
        &mut self,
        displacement_m_sqrt_kg: &[C64],
        gate: &CancelGate,
    ) -> Result<ShellCoatingStrainReport<'a>, AcousticRealizeError> {
        let summary =
            self.plan
                .evaluate(self.source, Motion::Harmonic(displacement_m_sqrt_kg), gate)?;
        Ok(self.report(summary))
    }

    /// Observe an actual real body displacement at one instant. Static offsets
    /// already present in these coordinates participate coherently. Removed
    /// prestress/cure strain and motion between observed instants are excluded.
    /// The original harmonic preparation allowance conservatively admits these
    /// cheaper real queries; a successful query consumes one planned evaluation.
    pub fn evaluate_snapshot(
        &mut self,
        displacement_m_sqrt_kg: &[f64],
        gate: &CancelGate,
    ) -> Result<ShellCoatingStrainReport<'a>, AcousticRealizeError> {
        let summary = self
            .plan
            .evaluate_snapshot(self.source, displacement_m_sqrt_kg, gate)?;
        Ok(self.report(summary))
    }

    fn report(&self, summary: ShellCoatingStrainSummary) -> ShellCoatingStrainReport<'a> {
        let films = self
            .source
            .film_sources()
            .zip(&self.plan.values)
            .map(|(source, value)| ShellCoatingFilmStrain {
                component_index: source.component_index,
                component_name: source.component_name,
                film_index: source.film_index,
                source: source.source,
                covered_elements: value.covered_elements,
                frobenius_upper_bound: value.frobenius_upper_bound,
                component_amplitudes: value.component_amplitudes,
                ruling_element: value.ruling_element,
                ruling_control_point: value.ruling_control_point,
                ruling_z_face: value.ruling_z_face,
                limit_ratio: value.limit_ratio,
                sufficient_amplitude_scale: value.sufficient_amplitude_scale,
            })
            .collect();
        ShellCoatingStrainReport {
            films,
            within_authored_envelope: summary.within_authored_envelope,
            sufficient_common_amplitude_scale: summary.sufficient_common_amplitude_scale,
        }
    }
}

#[derive(Clone, Copy)]
enum Motion<'a> {
    Harmonic(&'a [C64]),
    Snapshot(&'a [f64]),
}

impl Motion<'_> {
    fn valid(&self, modes: usize) -> bool {
        match self {
            Self::Harmonic(q) => {
                q.len() == modes && q.iter().all(|z| z.re.is_finite() && z.im.is_finite())
            }
            Self::Snapshot(q) => q.len() == modes && q.iter().all(|x| x.is_finite()),
        }
    }
}

impl OwnedShellCoatingStrain {
    pub(crate) const fn resources(&self) -> &ShellCoatingStrainResources {
        &self.resources
    }
    pub(crate) const fn film_count(&self) -> usize {
        self.films
    }
    pub(crate) const fn mode_count(&self) -> usize {
        self.modes
    }
    pub(crate) const fn completed_evaluations(&self) -> usize {
        self.completed_evaluations
    }
    /// Scratch values from the most recent query; publish them only after that
    /// query succeeds. Failure changes no physical state or successful count.
    pub(crate) fn values(&self) -> &[ShellCoatingFilmStrainValue] {
        &self.values
    }

    pub(crate) fn evaluate_snapshot(
        &mut self,
        source: CoatingStrainSource<'_>,
        displacement_m_sqrt_kg: &[f64],
        gate: &CancelGate,
    ) -> Result<ShellCoatingStrainSummary, AcousticRealizeError> {
        self.evaluate(source, Motion::Snapshot(displacement_m_sqrt_kg), gate)
    }

    fn evaluate(
        &mut self,
        source: CoatingStrainSource<'_>,
        motion: Motion<'_>,
        gate: &CancelGate,
    ) -> Result<ShellCoatingStrainSummary, AcousticRealizeError> {
        poll(gate)?;
        if self.completed_evaluations == self.resources.evaluations {
            return Err(invalid(
                "coating strain observation exhausted its planned evaluations",
            ));
        }
        if !motion.valid(self.modes) || (self.snapshot && matches!(motion, Motion::Harmonic(_))) {
            return Err(invalid(
                "coating strain requires admitted finite displacement on every original body mode",
            ));
        }
        if source.len() != self.parts.len()
            || self.assembly != matches!(source, CoatingStrainSource::Assembly(_))
            || self
                .parts
                .iter()
                .enumerate()
                .any(|(i, part)| source.part(i).storage_identity() != part.source_storage)
        {
            return Err(invalid(
                "coating strain observation requires its unchanged original source image",
            ));
        }
        self.values.fill(ShellCoatingFilmStrainValue::default());
        let mut offset = 0;
        for (index, part) in self.parts.iter().enumerate() {
            poll(gate)?;
            let source_part = source.part(index);
            let films = source_part.shell.coatings();
            let values = &mut self.values[offset..offset + films.len()];
            offset += films.len();
            for element in &part.elements {
                poll(gate)?;
                match motion {
                    Motion::Snapshot(q) => {
                        let mut real = [0.0; 18];
                        for (mode, amplitude) in source_part.modes.iter().zip(q) {
                            for (local_node, &node) in element.nodes.iter().enumerate() {
                                for dof in 0..6 {
                                    real[6 * local_node + dof] += mode[6 * node + dof] * amplitude;
                                }
                            }
                        }
                        let field = element
                            .field
                            .real_field(&real, &mut || gate.is_requested())
                            .map_err(strain_owner)?;
                        observe_element(
                            element,
                            &part.film_regions,
                            films,
                            values,
                            |faces, nu| {
                                field.isotropic_layer_bound(faces, nu, &mut || gate.is_requested())
                            },
                        )?;
                    }
                    Motion::Harmonic(q) => {
                        let mut real = [0.0; 18];
                        let mut imaginary = [0.0; 18];
                        for (mode, amplitude) in source_part.modes.iter().zip(q) {
                            for (local_node, &node) in element.nodes.iter().enumerate() {
                                for dof in 0..6 {
                                    let shape = mode[6 * node + dof];
                                    real[6 * local_node + dof] += shape * amplitude.re;
                                    imaginary[6 * local_node + dof] += shape * amplitude.im;
                                }
                            }
                        }
                        let field = element
                            .field
                            .harmonic_field(&real, &imaginary, &mut || gate.is_requested())
                            .map_err(strain_owner)?;
                        observe_element(
                            element,
                            &part.film_regions,
                            films,
                            values,
                            |faces, nu| {
                                field.isotropic_layer_bound(faces, nu, &mut || gate.is_requested())
                            },
                        )?;
                    }
                }
            }
        }
        let mut sufficient_common_amplitude_scale: Option<f64> = None;
        let mut within_authored_envelope = true;
        for (film, source) in self.values.iter_mut().zip(source.film_sources()) {
            let limit = source.source.input().material.strain_limit();
            film.limit_ratio = film.frobenius_upper_bound / limit;
            if film.covered_elements == 0 || !film.limit_ratio.is_finite() {
                return Err(invalid(
                    "coating strain requires complete finite source-film bounds",
                ));
            }
            within_authored_envelope &= film.limit_ratio <= 1.0;
            if film.frobenius_upper_bound > 0.0 {
                // Downward rounding avoids claiming the next representable
                // scale above the computed limit; positive overflow is capped.
                let scale = (limit / film.frobenius_upper_bound)
                    .min(f64::MAX)
                    .next_down()
                    .max(0.0);
                film.sufficient_amplitude_scale = Some(scale);
                sufficient_common_amplitude_scale =
                    Some(sufficient_common_amplitude_scale.map_or(scale, |s| s.min(scale)));
            }
        }
        poll(gate)?;
        self.completed_evaluations += 1;
        Ok(ShellCoatingStrainSummary {
            within_authored_envelope,
            sufficient_common_amplitude_scale,
        })
    }
}

fn observe_element(
    element: &Element,
    film_regions: &[usize],
    films: &[ResolvedShellCoating],
    values: &mut [ShellCoatingFilmStrainValue],
    mut layer: impl FnMut([f64; 2], f64) -> Result<TensorStrainBound, fs_plate::PlateError>,
) -> Result<(), AcousticRealizeError> {
    for ((film, &region), value) in films.iter().zip(film_regions).zip(values) {
        if region != element.region {
            continue;
        }
        let faces = film.faces_m();
        let bound = layer([faces.0, faces.1], film.input().material.poisson_ratio())
            .map_err(strain_owner)?;
        if value.covered_elements == 0 || bound.frobenius_upper_bound > value.frobenius_upper_bound
        {
            value.frobenius_upper_bound = bound.frobenius_upper_bound;
            value.ruling_element = element.index;
            value.ruling_control_point = bound.ruling_control_point;
            value.ruling_z_face = bound.ruling_z_face;
        }
        for (stored, candidate) in value
            .component_amplitudes
            .iter_mut()
            .zip(bound.component_amplitudes)
        {
            *stored = stored.max(candidate);
        }
        value.covered_elements += 1;
    }
    Ok(())
}
