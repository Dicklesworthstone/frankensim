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
    strain::{HarmonicElementStrain, PreparedElementStrain, TensorStrainBound},
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

/// Whole-film amplitude envelope for the supplied peak complex body motion.
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

struct Element {
    index: usize,
    region: usize,
    nodes: [usize; 3],
    field: PreparedElementStrain,
}

struct Part<'a> {
    component_index: Option<usize>,
    component_name: Option<&'a str>,
    source: &'a ResolvedShearShell,
    modes: &'a [Vec<f64>],
    film_regions: Vec<usize>,
    elements: Vec<Element>,
}

/// Prepared local observation tied to one original shell/assembly source and
/// its retained modes. No foreign directors, sections or material limits can
/// be substituted. Query errors and cancellation leave its allowance intact.
pub struct PreparedShellCoatingStrain<'a> {
    parts: Vec<Part<'a>>,
    resources: ShellCoatingStrainResources,
    modes: usize,
    films: usize,
    completed_evaluations: usize,
}

enum Source<'a> {
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

impl<'a> Source<'a> {
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
        prepare(Source::Shell(self), budget, gate)
    }
}

impl MaterialShellAssemblyImage {
    /// Observe every component film in the SAME assembled modal coordinates.
    pub fn prepare_coating_strain(
        &self,
        budget: ShellCoatingStrainBudget,
        gate: &CancelGate,
    ) -> Result<PreparedShellCoatingStrain<'_>, AcousticRealizeError> {
        prepare(Source::Assembly(self), budget, gate)
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

fn prepare<'a>(
    source: Source<'a>,
    budget: ShellCoatingStrainBudget,
    gate: &CancelGate,
) -> Result<PreparedShellCoatingStrain<'a>, AcousticRealizeError> {
    poll(gate)?;
    if budget.evaluations == 0 || source.len() == 0 {
        return Err(invalid(
            "coating strain requires a source and positive planned evaluations",
        ));
    }
    let modes = source.part(0).modes.len();
    let mut films = 0usize;
    let mut resources = ShellCoatingStrainResources {
        retained_bytes: size_of::<PreparedShellCoatingStrain<'a>>(),
        evaluation_bytes: size_of::<ShellCoatingStrainReport<'a>>()
            + HarmonicElementStrain::RETAINED_BYTES
            + size_of::<TensorStrainBound>()
            + 36 * size_of::<f64>(),
        preparation_work: 1,
        per_evaluation_work: 2 * modes,
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
        add_product(&mut resources.retained_bytes, 1, size_of::<Part<'a>>())?;
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
                PreparedElementStrain::LAYER_BOUND_WORK + 32,
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
                72 * modes + PreparedElementStrain::HARMONIC_FIELD_WORK + count,
            )?;
        }
    }
    add_product(
        &mut resources.evaluation_bytes,
        films,
        size_of::<ShellCoatingFilmStrain<'a>>(),
    )?;
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
            component_index: part.index,
            component_name: part.name,
            source: part.shell,
            modes: part.modes,
            film_regions,
            elements,
        });
    }
    poll(gate)?;
    Ok(PreparedShellCoatingStrain {
        parts,
        resources,
        modes,
        films,
        completed_evaluations: 0,
    })
}

impl<'a> PreparedShellCoatingStrain<'a> {
    pub const fn resources(&self) -> &ShellCoatingStrainResources {
        &self.resources
    }
    pub const fn film_count(&self) -> usize {
        self.films
    }
    pub const fn completed_evaluations(&self) -> usize {
        self.completed_evaluations
    }

    pub fn film_sources(&self) -> impl Iterator<Item = ShellCoatingFilmSource<'a>> + '_ {
        self.parts.iter().flat_map(|part| {
            part.source
                .coatings()
                .iter()
                .enumerate()
                .map(move |(film_index, source)| ShellCoatingFilmSource {
                    component_index: part.component_index,
                    component_name: part.component_name,
                    film_index,
                    source,
                })
        })
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
        poll(gate)?;
        let q = displacement_m_sqrt_kg;
        if self.completed_evaluations == self.resources.evaluations {
            return Err(invalid(
                "coating strain observation exhausted its planned evaluations",
            ));
        }
        if q.len() != self.modes || q.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
            return Err(invalid(
                "coating strain requires one finite complex displacement per original body mode",
            ));
        }
        let mut films = Vec::with_capacity(self.films);
        for part in &self.parts {
            poll(gate)?;
            let offset = films.len();
            for (film_index, source) in part.source.coatings().iter().enumerate() {
                films.push(ShellCoatingFilmStrain {
                    component_index: part.component_index,
                    component_name: part.component_name,
                    film_index,
                    source,
                    covered_elements: 0,
                    frobenius_upper_bound: 0.0,
                    component_amplitudes: [0.0; 6],
                    ruling_element: 0,
                    ruling_control_point: 0,
                    ruling_z_face: 0,
                    limit_ratio: 0.0,
                    sufficient_amplitude_scale: None,
                });
            }
            for element in &part.elements {
                poll(gate)?;
                let mut real = [0.0; 18];
                let mut imag = [0.0; 18];
                for (mode, amplitude) in part.modes.iter().zip(q) {
                    for (local_node, &node) in element.nodes.iter().enumerate() {
                        for dof in 0..6 {
                            let shape = mode[6 * node + dof];
                            real[6 * local_node + dof] += shape * amplitude.re;
                            imag[6 * local_node + dof] += shape * amplitude.im;
                        }
                    }
                }
                let field = element
                    .field
                    .harmonic_field(&real, &imag, &mut || gate.is_requested())
                    .map_err(strain_owner)?;
                for (film_index, &region) in part.film_regions.iter().enumerate() {
                    if region != element.region {
                        continue;
                    }
                    let report = &mut films[offset + film_index];
                    let faces = report.source.faces_m();
                    let bound = field
                        .isotropic_layer_bound(
                            [faces.0, faces.1],
                            report.source.input().material.poisson_ratio(),
                            &mut || gate.is_requested(),
                        )
                        .map_err(strain_owner)?;
                    if report.covered_elements == 0
                        || bound.frobenius_upper_bound > report.frobenius_upper_bound
                    {
                        report.frobenius_upper_bound = bound.frobenius_upper_bound;
                        report.ruling_element = element.index;
                        report.ruling_control_point = bound.ruling_control_point;
                        report.ruling_z_face = bound.ruling_z_face;
                    }
                    for (stored, candidate) in report
                        .component_amplitudes
                        .iter_mut()
                        .zip(bound.component_amplitudes)
                    {
                        *stored = stored.max(candidate);
                    }
                    report.covered_elements += 1;
                }
            }
        }
        let mut sufficient_common_amplitude_scale: Option<f64> = None;
        let mut within_authored_envelope = true;
        for film in &mut films {
            let limit = film.source.input().material.strain_limit();
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
        Ok(ShellCoatingStrainReport {
            films,
            within_authored_envelope,
            sufficient_common_amplitude_scale,
        })
    }
}
