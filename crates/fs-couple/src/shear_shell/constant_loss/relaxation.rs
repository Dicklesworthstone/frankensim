//! Causal constant-loss substrate memory in the existing MITC shell basis.
//!
//! The source's viscosity anchor defines loss stiffness H = omega_ref V.
//! It does not identify the elastic stiffness's reference. The caller supplies
//! that independent interpretation explicitly. A calibrated storage tensor is
//! reduced to its positive equilibrium tensor before the original eigensolve;
//! the same five Maxwell branches then supply both storage and dissipation.
//! This is a finite-band constitutive approximation with an explicitly selected
//! passive continuation, not additional measured material evidence.

use super::{BRANCHES, ConstantLossSpectrum};
use crate::acoustic_realize::AcousticRealizeError;
use crate::shear_shell::{
    PlateMaterialModel, ResolvedShearShell, ResolvedShellRegion,
    compile_material_shear_shell, invalid, owner, plate_elasticity, plate_property, poll,
};
use fs_exec::CancelGate;
use fs_material::sheet::SheetMatrix;
use fs_plate::shell::mitc::{
    MitcBudget, NodalDirectors, ShearShellSection,
    viscosity::{ShearShellViscosity, project_viscosity},
};
use fs_qty::{Density, Pressure};

/// Independent interpretation of the supplied elastic section, selected by
/// the caller rather than inferred from the viscosity encoding frequency.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SubstrateStorageReference {
    /// The supplied elastic tensor is the fully relaxed (zero-frequency) tensor.
    Relaxed,
    /// The supplied elastic tensor is storage at this angular frequency [rad/s].
    /// Subtract the Maxwell storage at this frequency before assembly. The
    /// original tensor is recovered at the reference without double stiffness.
    AtAngularFrequency(f64),
}

/// Selecting this policy explicitly supplies the causal Prony continuation
/// outside the original finite source band. Source support is not widened.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SubstrateRelaxationOptions {
    pub storage_reference: SubstrateStorageReference,
    /// Complete substrate coordinates after exact equal-rate aggregation.
    /// For assemblies each selected component's ceiling applies to the total.
    pub max_memory_coordinates: usize,
}

impl SubstrateRelaxationOptions {
    pub(crate) fn validate(self) -> Result<(), AcousticRealizeError> {
        if !(1..=8192).contains(&self.max_memory_coordinates)
            || matches!(self.storage_reference,
                SubstrateStorageReference::AtAngularFrequency(w) if !w.is_finite() || w <= 0.0)
        {
            return Err(invalid("substrate memory needs an explicit positive storage reference and 1..=8192 coordinates"));
        }
        Ok(())
    }

    pub(crate) fn projection_terms(self, modes: usize) -> Option<usize> {
        // The single positive loss form is factored once, then its rows are
        // scaled by sqrt(g_j). Every branch shares the same five source rates.
        let square = modes.checked_mul(modes)?;
        square.checked_mul(modes)?.checked_mul(4)?
            .checked_add(square.checked_mul(BRANCHES + 2)?)
    }
}

/// Complete positive substrate spectrum on one unchanged equilibrium basis.
#[derive(Clone, Debug)]
pub struct ShellSubstrateRelaxation {
    loss_stiffness: Vec<f64>,
    rows: Vec<Vec<f64>>,
    times: Vec<f64>,
    modes: usize,
    band: (f64, f64),
    spectrum: ConstantLossSpectrum,
}
impl ShellSubstrateRelaxation {
    pub const fn mode_count(&self) -> usize { self.modes }
    pub fn coupling_rows_per_s(&self) -> &[Vec<f64>] { &self.rows }
    pub fn relaxation_times_s(&self) -> &[f64] { &self.times }
    /// Original finite material/fit support, distinct from the continuation.
    pub const fn source_band_rad_s(&self) -> (f64, f64) { self.band }
    /// H [1/s²] in the original mass-normalized basis. The complex stiffness
    /// increment is (storage_shape(omega) - i shape(omega)) H for exp(-i omega t).
    pub fn loss_stiffness_per_s2(&self) -> &[f64] { &self.loss_stiffness }
    /// Storage and positive loss multipliers of the retained loss-stiffness form.
    pub fn modulus_shapes(&self, omega_rad_s: f64) -> (f64, f64) {
        (self.spectrum.storage_shape(omega_rad_s), self.spectrum.shape(omega_rad_s))
    }

    pub(crate) fn from_projected(
        parts: Vec<ProjectedSubstrate>, omega: &[f64], max_work: usize, gate: &CancelGate,
    ) -> Result<Option<Self>, AcousticRealizeError> {
        poll(gate)?;
        if parts.is_empty() { return Ok(None); }
        let n = omega.len();
        if n == 0 || n > crate::shear_shell::MAX_SHELL_MODES {
            return Err(invalid("substrate memory needs a complete bounded body basis"));
        }
        let mut form = vec![0.0; n*n];
        let mut band = (0.0_f64, f64::MAX);
        let mut maximum = usize::MAX;
        for part in parts {
            poll(gate)?;
            part.options.validate()?;
            if part.form.len() != n*n
                || part.options.projection_terms(n).is_none_or(|work| work > max_work) {
                return Err(invalid("substrate memory exceeds its original projection allowance"));
            }
            maximum = maximum.min(part.options.max_memory_coordinates);
            band.0 = band.0.max(part.band.0);
            band.1 = band.1.min(part.band.1);
            for (out, value) in form.iter_mut().zip(part.form) { *out += value; }
        }
        if band.1 < band.0 || omega.iter().any(|w|
            !w.is_finite() || *w <= 0.0 || *w < band.0 || *w > band.1)
            || form.iter().any(|value| !value.is_finite()) {
            return Err(invalid("substrate memory has unsupported frequencies or unrepresentable stiffness"));
        }
        let factors = crate::plate_viscosity::relaxation::factor(&form, n, gate)?;
        let count = factors.len().checked_mul(BRANCHES)
            .ok_or_else(|| invalid("substrate memory coordinate extent overflow"))?;
        if count > maximum {
            return Err(invalid("complete substrate memory exceeds its declared coordinate ceiling"));
        }
        let mut rows = Vec::with_capacity(count);
        let mut times = Vec::with_capacity(count);
        let spectrum = ConstantLossSpectrum::fit();
        for (rate, weight) in spectrum.branches() {
            poll(gate)?;
            if !rate.is_finite() || rate <= 0.0 || !weight.is_finite() || weight <= 0.0 {
                return Err(invalid("substrate Maxwell spectrum is not positive and representable"));
            }
            let scale = fs_math::det::sqrt(weight);
            for factor in &factors {
                let row: Vec<_> = factor.iter().map(|value| scale * value).collect();
                if row.iter().any(|value| !value.is_finite())
                    || row.iter().all(|value| *value == 0.0) {
                    return Err(invalid("substrate material row is unrepresentable"));
                }
                rows.push(row);
                times.push(1.0 / rate);
            }
        }
        poll(gate)?;
        Ok(Some(Self { loss_stiffness: form, rows, times, modes: n, band, spectrum }))
    }
}

pub(crate) struct ProjectedSubstrate {
    form: Vec<f64>,
    band: (f64, f64),
    options: SubstrateRelaxationOptions,
}

impl ResolvedShearShell {
    pub const fn substrate_relaxation_options(&self) -> Option<SubstrateRelaxationOptions> {
        self.substrate_relaxation
    }

    /// Replace ONLY anchored substrate dashpots by their complete causal
    /// spectrum. Unanchored Kelvin--Voigt regions, films and interfaces keep
    /// their independent laws. The source records remain unchanged.
    ///
    /// An anchored stiffness is converted to equilibrium before the original
    /// shell assembly. Both the bare substrate and the full coated section
    /// must remain positive definite; a coating cannot hide an active core.
    /// The reference is an authored constitutive interpretation, not a new
    /// claim about how an unannotated material measurement was made.
    ///
    /// Body-bank constructors attach these explicitly selected histories to
    /// the existing mechanical owner, initially relaxed at its displacement.
    pub fn with_substrate_relaxation(
        &self, options: SubstrateRelaxationOptions, gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        options.validate()?;
        let substrate = compile_material_shear_shell(self.mesh.clone(),
            self.regions.iter().map(|region| region.input.clone()).collect(), self.budget, gate)?;
        let coated = crate::shear_shell::coating::apply_coatings(substrate,
            self.coatings.iter().map(|film| film.input().clone()).collect(),
            self.coating_relaxation, gate)?;
        let mut rebuilt = apply_selection(coated, Some(options), gate)?.with_interfaces(
            self.interfaces.iter().map(|link| link.input().clone()).collect(),
            self.interface_budget, gate)?;
        rebuilt.pressure_surfaces = self.pressure_surfaces.clone();
        poll(gate)?;
        Ok(rebuilt)
    }

    pub(crate) fn project_substrate_form(
        &self, directors: &NodalDirectors, physical_modes: &[Vec<f64>], omega: &[f64],
        budget: MitcBudget, max_projection_terms: usize, gate: &CancelGate,
    ) -> Result<Option<ProjectedSubstrate>, AcousticRealizeError> {
        poll(gate)?;
        let Some(options) = self.substrate_relaxation else { return Ok(None); };
        let n = physical_modes.len();
        if n == 0 || n > crate::shear_shell::MAX_SHELL_MODES || omega.len() != n
            || self.loss_projection_terms(n).is_none_or(|work| work > max_projection_terms) {
            return Err(invalid("substrate projection needs complete bounded modes and work"));
        }
        let mut form = vec![0.0; n*n];
        let mut band = (0.0_f64, f64::MAX);
        for region in self.regions.iter().filter(|region| region.loss_reference_omega.is_some()) {
            band.0 = band.0.max(region.omega_band.0);
            band.1 = band.1.min(region.omega_band.1);
        }
        if omega.iter().any(|w| !w.is_finite() || *w <= 0.0 || *w < band.0 || *w > band.1) {
            return Err(invalid("substrate mode exceeds its original material/fit support"));
        }
        for (reference, viscosity) in self.loss_groups()? {
            let Some(reference) = reference else { continue; };
            poll(gate)?;
            let part = project_viscosity(&self.mesh, directors, &self.sections, &viscosity,
                physical_modes, budget, max_projection_terms, || gate.is_requested()).map_err(owner)?;
            for (out, value) in form.iter_mut().zip(part) { *out += reference * value; }
        }
        if form.iter().any(|value| !value.is_finite()) {
            return Err(invalid("substrate loss stiffness projection overflow"));
        }
        poll(gate)?;
        Ok(Some(ProjectedSubstrate { form, band, options }))
    }

    pub fn project_substrate_relaxation(
        &self, directors: &NodalDirectors, physical_modes: &[Vec<f64>], omega: &[f64],
        budget: MitcBudget, max_projection_terms: usize, gate: &CancelGate,
    ) -> Result<Option<ShellSubstrateRelaxation>, AcousticRealizeError> {
        let forms = self.project_substrate_form(directors, physical_modes, omega,
            budget, max_projection_terms, gate)?.into_iter().collect();
        ShellSubstrateRelaxation::from_projected(forms, omega, max_projection_terms, gate)
    }
}

/// Apply after films are assembled, so their unchanged reference laws, mass
/// moments and B coupling survive. Every rebind starts from original sources.
pub(crate) fn apply_selection(
    mut shell: ResolvedShearShell, options: Option<SubstrateRelaxationOptions>, gate: &CancelGate,
) -> Result<ResolvedShearShell, AcousticRealizeError> {
    poll(gate)?;
    let Some(options) = options else { return Ok(shell); };
    options.validate()?;
    if shell.substrate_relaxation.is_some()
        || !shell.regions.iter().any(|region| region.loss_reference_omega.is_some()) {
        return Err(invalid("substrate memory needs a fresh shell with an explicit source loss anchor"));
    }
    let spectrum = ConstantLossSpectrum::fit();
    for region in &mut shell.regions {
        poll(gate)?;
        let Some(reference) = region.loss_reference_omega else { continue; };
        let storage = match options.storage_reference {
            SubstrateStorageReference::Relaxed => 0.0,
            SubstrateStorageReference::AtAngularFrequency(omega) => {
                if omega < region.omega_band.0 || omega > region.omega_band.1 {
                    return Err(invalid("substrate storage reference lies outside its original material support"));
                }
                reference * spectrum.storage_shape(omega)
            }
        };
        if !storage.is_finite() {
            return Err(invalid("substrate equilibrium storage subtraction overflow"));
        }
        if storage != 0.0 {
            // Validate the substrate by itself, before the coating could mask
            // negative equilibrium storage in any physical strain direction.
            let bare = original_substrate_section(region)?;
            subtract_storage(&bare, &region.viscosity, storage)?;
            region.section = subtract_storage(&region.section, &region.viscosity, storage)?;
        }
    }
    for (element, &region_index) in shell.element_regions.iter().enumerate() {
        poll(gate)?;
        let region = &shell.regions[region_index];
        let angle = fs_plate::shell::survey::material_axis_angle(
            &shell.mesh.facet(element).map_err(owner)?, region.input.axis1_global).map_err(owner)?
            + region.input.region.material_angle_rad;
        shell.sections[element] = region.section.with_material_angle_rad(angle).map_err(owner)?;
    }
    shell.substrate_relaxation = Some(options);
    poll(gate)?;
    Ok(shell)
}

fn original_substrate_section(region: &ResolvedShellRegion)
    -> Result<ShearShellSection, AcousticRealizeError> {
    let input = &region.input.region;
    let source = &input.material;
    let rho = plate_property(source, "density", Density::DIMS)?;
    let (e1, e2, nu, g12) = plate_elasticity(source, input.model)?;
    let (g13, g23) = if input.model.is_isotropic() { (g12, g12) } else {
        let a = plate_property(source, "shear_modulus_13", Pressure::DIMS)?;
        let b = plate_property(source, "shear_modulus_23", Pressure::DIMS)?;
        if input.model == PlateMaterialModel::Orthotropic21 { (b, a) } else { (a, b) }
    };
    ShearShellSection::from_orthotropic_solid(e1, e2, nu, g12, g13, g23,
        region.thickness_m, rho, 0.0, region.input.shear_correction).map_err(owner)
}

fn subtract_storage(
    section: &ShearShellSection, viscosity: &ShearShellViscosity, scale: f64,
) -> Result<ShearShellSection, AcousticRealizeError> {
    let (a, d) = section.principal_matrices();
    let (va, vb, vd, vs) = viscosity.coupled_matrices();
    let source_shear = section.principal_transverse_shear_n_m();
    let subtract = |elastic: &[f64; 9], loss: &[f64; 9]| {
        core::array::from_fn(|i| elastic[i] - scale * loss[i])
    };
    let shear = [source_shear[0] - scale * vs[0][0],
        source_shear[1] - scale * vs[0][1], source_shear[2] - scale * vs[1][1]];
    ShearShellSection::new_coupled(
        SheetMatrix::from_coefficients(subtract(a.coefficients(), &va)).map_err(owner)?,
        subtract(section.coupling_pa_m2(), &vb),
        SheetMatrix::from_coefficients(subtract(d.coefficients(), &vd)).map_err(owner)?,
        shear, 1.0, section.areal_mass_kg_m2(), section.mass_first_moment_kg_m(),
        section.rotary_inertia_kg(), 0.0,
    ).map_err(owner)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn storage_reference_is_recovered_without_changing_mass_or_double_stiffness() {
        let source = ShearShellSection::from_orthotropic_solid(
            12e9, 0.8e9, 0.22, 0.7e9, 0.6e9, 0.08e9,
            0.003, 450.0, 0.0, 5.0/6.0).unwrap();
        let (a, d) = source.principal_matrices();
        let shear = source.principal_transverse_shear_n_m();
        let reference = core::f64::consts::TAU * 200.0;
        let anchor = core::f64::consts::TAU * 1000.0;
        let loss = 0.02 / reference;
        let viscosity = ShearShellViscosity::new(
            a.coefficients().map(|v| v*loss),
            d.coefficients().map(|v| v*loss), shear.map(|v| v*loss)).unwrap();
        let scale = reference * ConstantLossSpectrum::fit().storage_shape(anchor);
        let equilibrium = subtract_storage(&source, &viscosity, scale).unwrap();
        assert_eq!(equilibrium.areal_mass_kg_m2(), source.areal_mass_kg_m2());
        assert_eq!(equilibrium.rotary_inertia_kg(), source.rotary_inertia_kg());
        for (eq, original) in equilibrium.principal_matrices().0.coefficients().iter()
            .zip(a.coefficients()) {
            let restored = eq + original * loss * scale;
            assert!((restored - original).abs() <= 4.0*f64::EPSILON*original.abs());
        }
        assert!(subtract_storage(&source, &viscosity, 2.0/loss).is_err());
    }

    #[test]
    fn projected_cross_terms_become_the_same_positive_causal_spectrum() {
        let gate = CancelGate::new_clock_free();
        let form = vec![4.0, -6.0, -6.0, 9.0];
        let options = SubstrateRelaxationOptions {
            storage_reference: SubstrateStorageReference::Relaxed,
            max_memory_coordinates: BRANCHES,
        };
        let make = |options| ProjectedSubstrate {
            form: form.clone(), band: (100.0, 10_000.0), options,
        };
        let image = ShellSubstrateRelaxation::from_projected(
            vec![make(options)], &[200.0, 2000.0], 100_000, &gate).unwrap().unwrap();
        assert_eq!(image.coupling_rows_per_s().len(), BRANCHES);
        let spectrum = ConstantLossSpectrum::fit();
        for omega in [200.0_f64, 1000.0, 8000.0] {
            let mut storage = vec![0.0; 4];
            let mut loss = vec![0.0; 4];
            for (row, &tau) in image.rows.iter().zip(&image.times) {
                let x = omega*tau;
                for i in 0..2 { for j in 0..2 {
                    storage[2*i+j] += row[i]*row[j]*super::super::debye_storage(x);
                    loss[2*i+j] += row[i]*row[j]*super::super::debye(x);
                }}
                // Constant-strain free relaxation: spring storage lost is
                // exactly the integrated positive dashpot dissipation.
                let strain = row[0] - 2.0*row[1];
                let dt = 0.4*tau;
                let remaining = strain*fs_math::det::exp(-dt/tau);
                let released = 0.5*(strain*strain - remaining*remaining);
                let heat = 0.5*strain*strain*(1.0-fs_math::det::exp(-2.0*dt/tau));
                assert!(heat >= 0.0);
                assert!((released-heat).abs() <= 2e-15*strain*strain);
            }
            for i in 0..4 {
                assert!((storage[i] - spectrum.storage_shape(omega)*form[i]).abs() < 1e-11);
                assert!((loss[i] - spectrum.shape(omega)*form[i]).abs() < 1e-11);
            }
        }
        let mut too_small = options;
        too_small.max_memory_coordinates = BRANCHES-1;
        assert!(ShellSubstrateRelaxation::from_projected(
            vec![make(too_small)], &[200.0, 2000.0], 100_000, &gate).is_err());
        assert!(ShellSubstrateRelaxation::from_projected(
            vec![make(options)], &[50.0, 2000.0], 100_000, &gate).is_err());
    }
}
