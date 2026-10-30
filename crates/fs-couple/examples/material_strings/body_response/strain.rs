//! Optional observation of the original local coating field. All mechanics,
//! condensation, material closure and spatial bounds belong to the shell owner.
use super::*;
use fs_couple::shear_shell::coating::strain::{
    PreparedShellCoatingStrain, ShellCoatingStrainBudget, ShellCoatingStrainReport,
};

/// The same numeric bound is either a unit-force transfer or the strain of
/// the complete authored physical load. Its sufficient scale follows that input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Units {
    PerNewton,
    Physical,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Options {
    pub maximum_bytes: usize,
    pub maximum_work: usize,
}

impl Options {
    pub(super) fn parse(values: &[OsString]) -> Result<Self> {
        if values.len() != 2 {
            return Err("--body-response-strain needs MAX_BYTES MAX_WORK".into());
        }
        let integer = |value: &OsString| -> Result<usize> {
            Ok(value.to_str().ok_or("strain budgets must be UTF-8 integers")?.parse()?)
        };
        let maximum_bytes = integer(&values[0])?;
        let maximum_work = integer(&values[1])?;
        if maximum_bytes == 0 || maximum_work == 0 {
            return Err("coating strain requires positive explicit byte and work allowances".into());
        }
        Ok(Self { maximum_bytes, maximum_work })
    }

    pub(super) fn prepare<'a>(
        self,
        shell: Option<&'a MaterialShellImage>,
        assembly: Option<&'a MaterialShellAssemblyImage>,
        evaluations: usize,
        gate: &fs_exec::CancelGate,
    ) -> Result<PreparedShellCoatingStrain<'a>> {
        let budget = ShellCoatingStrainBudget {
            maximum_bytes: self.maximum_bytes,
            maximum_work: self.maximum_work,
            evaluations,
        };
        Ok(match (shell, assembly) {
            (Some(image), None) => image.prepare_coating_strain(budget, gate)?,
            (None, Some(image)) => image.prepare_coating_strain(budget, gate)?,
            _ => return Err("coating strain requires one original shell or assembly image".into()),
        })
    }
}

pub(super) fn additional_columns(films: usize) -> Result<usize> {
    films.checked_mul(12).and_then(|n| n.checked_add(3))
        .ok_or_else(|| "coating strain report shape overflow".into())
}

pub(super) fn metadata(text: &mut String, observer: &PreparedShellCoatingStrain<'_>, units: Units) -> Result<()> {
    writeln!(text, "# coating strain: peak full-tensor Frobenius bound on the original retained MITC harmonic increment; complete coated triangles and both film faces")?;
    writeln!(text, "# strain model: original elastic bubble condensation; plane-stress zz from film Poisson ratio; tensor shear=engineering shear/2; no shear-energy correction applied to strain")?;
    writeln!(text, "# strain scope: authored constitutive envelope, not a measured failure threshold; excludes unmodeled static/preload/cure strain, transients, missing modes, continuum and 3D interlaminar recovery; floating-point bounds retain the original owner's roundoff scope")?;
    writeln!(text, "# strain bound above limit means compliance is unresolved by this sufficient bound; the ruling Bernstein control need not be a physical maximum location")?;
    match units {
        Units::PerNewton => writeln!(text, "# sufficient force: independent 1 N drive in each row; scales the complete coherent response at that frequency; inf means no finite cap from these film bounds, not unlimited instrument strength")?,
        Units::Physical => writeln!(text, "# sufficient amplitude scale: dimensionless common multiplier of every authored complex bridge force, preserving relative phase; bounds and limit ratios describe the actual summed load; inf means no finite cap from these film bounds, not unlimited instrument strength")?,
    }
    writeln!(text, "# coating strain resources: {:?}", observer.resources())?;
    for (index, film) in observer.film_sources().enumerate() {
        let source = film.source;
        let input = source.input();
        let (bottom, top) = source.faces_m();
        writeln!(text,
            "# strain film {index}: component_index={:?} component_name={:?} source_film={} region={} source={} bottom_m={bottom} top_m={top} poisson_ratio={} authored_frobenius_limit={}",
            film.component_index, film.component_name, film.film_index, input.region,
            input.material.source().identity(), input.material.poisson_ratio(), input.material.strain_limit())?;
    }
    Ok(())
}

pub(super) fn header(text: &mut String, films: usize, units: Units) -> Result<()> {
    let physical = units == Units::Physical;
    if physical {
        write!(text, "\tcoating_strain_max_limit_ratio\tcoating_strain_sufficient_common_amplitude_scale\tcoating_strain_bound_within_envelope")?;
    } else {
        write!(text, "\tcoating_strain_max_limit_ratio_at_1N\tcoating_strain_sufficient_peak_force_N\tcoating_strain_bound_within_envelope_at_1N")?;
    }
    let suffix = if physical { "" } else { "_per_N" };
    for film in 0..films {
        if physical {
            write!(text, "\tfilm{film}_strain_frobenius_bound\tfilm{film}_strain_limit_ratio\tfilm{film}_strain_sufficient_amplitude_scale")?;
        } else {
            write!(text, "\tfilm{film}_strain_frobenius_bound_per_N\tfilm{film}_strain_limit_ratio_at_1N\tfilm{film}_strain_sufficient_peak_force_N")?;
        }
        for component in ["xx", "yy", "zz", "xy", "yz", "zx"] {
            write!(text, "\tfilm{film}_strain_{component}_amplitude_bound{suffix}")?;
        }
        write!(text, "\tfilm{film}_strain_ruling_element\tfilm{film}_strain_ruling_control\tfilm{film}_strain_ruling_z_face")?;
    }
    Ok(())
}

fn scale(text: &mut String, value: Option<f64>) -> Result<()> {
    match value {
        Some(value) => write!(text, "\t{value:.17e}")?,
        None => write!(text, "\tinf")?,
    }
    Ok(())
}

pub(super) fn row(text: &mut String, report: &ShellCoatingStrainReport<'_>) -> Result<()> {
    let ratio = report.films.iter().map(|film| film.limit_ratio).fold(0.0_f64, f64::max);
    write!(text, "\t{ratio:.17e}")?;
    scale(text, report.sufficient_common_amplitude_scale)?;
    write!(text, "\t{}", u8::from(report.within_authored_envelope))?;
    for film in &report.films {
        write!(text, "\t{:.17e}\t{:.17e}", film.frobenius_upper_bound, film.limit_ratio)?;
        scale(text, film.sufficient_amplitude_scale)?;
        for bound in film.component_amplitudes {
            write!(text, "\t{bound:.17e}")?;
        }
        write!(text, "\t{}\t{}\t{}", film.ruling_element, film.ruling_control_point, film.ruling_z_face)?;
    }
    Ok(())
}
