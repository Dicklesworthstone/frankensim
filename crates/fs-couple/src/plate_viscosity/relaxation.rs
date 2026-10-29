//! Regional Maxwell bending in the original DKT coordinates.
//!
//! A region's selected excess Young modulus scales its ORIGINAL rotated
//! equilibrium bending tensor. Different regions retain different spectra.
//! This is a fixed-state proportional tensor within each region, not a general
//! three-dimensional tensor Prony reduction. Prestress is never relaxed.
use super::{AcousticRealizeError, ResolvedPlateChart, invalid, owner, poll};
use crate::thin_plate::{PlateMaterialModel, RetainedPlateBasis};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::visco::GeneralizedMaxwell;
use fs_math::det;
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

/// Complete source tensors with hereditary normal-traction elimination.
pub mod tensor;

/// Explicit initial material history, independent of equilibrium displacement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InitialPlateMemory {
    /// Each arm initially has zero elastic extension at the original q.
    Relaxed,
    /// Every viscous strain initially vanishes; displaced material stores energy.
    Unrelaxed,
}

/// Cold projection and complete constitutive-state budgets. Null strain
/// directions carry no memory; every resolved positive channel is retained.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PlateRelaxationLimits {
    pub max_modes: usize,
    pub max_projection_terms: usize,
    pub max_memory_states: usize,
}

/// Immutable original chart selection and all source-resolved regional spectra.
#[derive(Clone, Debug)]
pub struct PlateBendingRelaxation {
    chart_identity: ContentHash,
    properties: Vec<(String, String)>,
    spectra: Vec<GeneralizedMaxwell>,
    bands: Vec<(f64, f64)>,
}
impl PlateBendingRelaxation {
    /// Select modulus [Pa] / time [s] pairs on EVERY original material region.
    /// The equilibrium modulus remains in the existing elastic pencil. No
    /// density, damping, geometry or material coefficient is inferred by name.
    /// Temperature-only curves are frozen at the original resolved state;
    /// frequency-sampled coefficients cannot masquerade as a hereditary law.
    pub fn from_chart(chart: &ResolvedPlateChart, properties: &[(String, String)],
        gate: &CancelGate) -> Result<Self, AcousticRealizeError>
    {
        poll(gate)?;
        if chart.bending_viscosity().is_some() || properties.is_empty() || properties.len() > 8
            || chart.regions().is_empty() || chart.regions().len() > 1024 {
            return Err(invalid("regional memory needs original regions, 1..=8 source pairs and no selected scalar material damping"));
        }
        let mut pairs = std::collections::BTreeSet::new();
        let mut moduli = std::collections::BTreeSet::new();
        let mut times = std::collections::BTreeSet::new();
        for (modulus, time) in properties {
            if modulus.is_empty() || time.is_empty() || modulus == time
                || !pairs.insert((modulus.as_str(), time.as_str()))
                || times.contains(modulus.as_str()) || moduli.contains(time.as_str()) {
                return Err(invalid("regional memory needs nonempty distinct pairs with consistent physical quantities"));
            }
            moduli.insert(modulus.as_str()); times.insert(time.as_str());
        }
        let mut spectra = Vec::with_capacity(chart.regions().len());
        let mut bands = Vec::with_capacity(chart.regions().len());
        for region in chart.regions() {
            poll(gate)?;
            let source = &region.input.material;
            let mut band = (0.0, f64::MAX);
            let elastic: &[(&str, Dims)] = match region.input.model {
                PlateMaterialModel::Isotropic => &[("density", Density::DIMS),
                    ("young_modulus", Pressure::DIMS), ("poisson_ratio", Dims::NONE)],
                _ => &[("density", Density::DIMS), ("young_modulus_1", Pressure::DIMS),
                    ("young_modulus_2", Pressure::DIMS), ("poisson_ratio_12", Dims::NONE),
                    ("shear_modulus_12", Pressure::DIMS)],
            };
            for &(key, dims) in elastic { coefficient(source, key, dims, &mut band)?; }
            let equilibrium_key = if region.input.model == PlateMaterialModel::Isotropic {
                "young_modulus"
            } else { "young_modulus_1" };
            let equilibrium = coefficient(source, equilibrium_key, Pressure::DIMS, &mut band)?;
            let terms = properties.iter().map(|(modulus, time)| Ok((
                coefficient(source, modulus, Pressure::DIMS, &mut band)?,
                coefficient(source, time, Dims([0, 0, 1, 0, 0, 0]), &mut band)?,
            ))).collect::<Result<Vec<_>, AcousticRealizeError>>()?;
            let law = GeneralizedMaxwell::new(equilibrium, terms).map_err(owner)?;
            if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0
                || law.terms.iter().any(|&(e, _)| e > 0.0 &&
                    (!(e/equilibrium).is_finite() || e/equilibrium == 0.0)) {
                return Err(invalid("regional material spectrum or intersected source band is unrepresentable"));
            }
            spectra.push(law);
            bands.push(band);
        }
        Ok(Self { chart_identity: chart.specimen_identity(), properties: properties.to_vec(), spectra, bands })
    }
    #[must_use] pub const fn chart_identity(&self) -> ContentHash { self.chart_identity }
    #[must_use] pub fn material_spectra(&self) -> &[GeneralizedMaxwell] { &self.spectra }
    #[must_use] pub fn properties(&self) -> &[(String, String)] { &self.properties }
    #[must_use] pub fn omega_bands_rad_s(&self) -> &[(f64, f64)] { &self.bands }
    /// A cold material replacement re-resolves every coefficient and band.
    pub fn rebind(&self, chart: &ResolvedPlateChart, gate: &CancelGate) -> Result<Self, AcousticRealizeError> {
        Self::from_chart(chart, &self.properties, gate)
    }
    pub(crate) fn validate_band(&self, low: f64, high: f64, dt: f64) -> Result<(), AcousticRealizeError> {
        if !low.is_finite() || !high.is_finite() || low < 0.0 || high < low
            || !dt.is_finite() || dt <= 0.0 || self.bands.iter().any(|&(a,b)| low < a || high > b)
            || self.spectra.iter().flat_map(|s| &s.terms).any(|&(e,t)| e > 0.0 && dt/t > 2.0) {
            return Err(invalid("complete coupled frequency band or clock is outside the original material relaxation support"));
        }
        Ok(())
    }
    pub(crate) fn project(&self, chart: &ResolvedPlateChart, basis: &RetainedPlateBasis,
        limits: PlateRelaxationLimits, gate: &CancelGate) -> Result<ModalPlateRelaxation, AcousticRealizeError>
    {
        poll(gate)?;
        let n = basis.mass_normalized_modes.len();
        if chart.specimen_identity() != self.chart_identity || n == 0 || n > limits.max_modes
            || !(1..=64).contains(&limits.max_modes) || limits.max_memory_states > 256 {
            return Err(invalid("regional memory projection needs the exact original chart and bounded complete modes/states"));
        }
        let forms = basis.regional_bending_forms(chart, &vec![1.0; n], limits.max_projection_terms, gate)?;
        if forms.len() != self.spectra.len() { return Err(invalid("regional material source layout changed")); }
        let mut rows = Vec::new();
        let mut times = Vec::new();
        for (form, law) in forms.iter().zip(&self.spectra) {
            poll(gate)?;
            let positive = law.terms.iter().filter(|(e,_)| *e > 0.0).count();
            if positive == 0 { continue; }
            let factors = factor(form, n, gate)?;
            if factors.len().checked_mul(positive).and_then(|a| rows.len().checked_add(a))
                .is_none_or(|a| a > limits.max_memory_states) {
                return Err(invalid("complete regional bending history exceeds its explicit state budget"));
            }
            for &(e,tau) in &law.terms {
                if e == 0.0 { continue; }
                let scale = det::sqrt(e/law.e_inf);
                for row in &factors {
                    let row: Vec<_> = row.iter().map(|c| scale*c).collect();
                    if row.iter().any(|v| !v.is_finite()) || row.iter().all(|v| *v == 0.0) {
                        return Err(invalid("regional material projection lost a positive strain channel"));
                    }
                    rows.push(row); times.push(tau);
                }
            }
        }
        poll(gate)?;
        Ok(ModalPlateRelaxation { source: self.clone(), rows, times, modes: n, limits })
    }
}

fn coefficient(source: &ResolvedMaterialStatePoint, key: &str, dims: Dims,
    band: &mut (f64,f64)) -> Result<f64, AcousticRealizeError>
{
    let p = source.property(key).ok_or_else(|| invalid("regional memory lacks a selected source coefficient"))?;
    if p.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !p.is_constant_at_fixed_temperature() || !p.value_si().is_finite() {
        return Err(invalid("regional memory needs dimensioned constant coefficients at the original material state"));
    }
    if let Some((lo,hi)) = p.angular_frequency_band_rad_s().map_err(owner)? {
        band.0 = band.0.max(lo); band.1 = band.1.min(hi);
    }
    Ok(p.value_si())
}

/// Positive strain rows c in ORIGINAL mass-normalized body coordinates.
/// Every channel stores `(c^T q-z)^2/2` [J] and dissipates that mismatch²/tau.
#[derive(Clone, Debug)]
pub struct ModalPlateRelaxation {
    source: PlateBendingRelaxation,
    rows: Vec<Vec<f64>>,
    times: Vec<f64>,
    modes: usize,
    limits: PlateRelaxationLimits,
}
impl ModalPlateRelaxation {
    #[must_use] pub const fn source(&self) -> &PlateBendingRelaxation { &self.source }
    #[must_use] pub fn coupling_rows_per_s(&self) -> &[Vec<f64>] { &self.rows }
    #[must_use] pub fn relaxation_times_s(&self) -> &[f64] { &self.times }
    #[must_use] pub const fn mode_count(&self) -> usize { self.modes }
    #[must_use] pub const fn limits(&self) -> PlateRelaxationLimits { self.limits }
}

/// Shared positive factor of a complete regional bending form. Localized
/// regions can have nonzero rows while leaving combinations of retained modes
/// unstrained. Those null directions carry no memory. Both reed and body
/// consumers retain the SAME complete form, including its cross terms, within
/// the existing dimensionless reconstruction tolerance; no stiffness is added.
pub(crate) fn factor(form: &[f64], n: usize, gate: &CancelGate)
    -> Result<Vec<Vec<f64>>, AcousticRealizeError>
{
    poll(gate)?;
    if n == 0 || n > crate::shear_shell::MAX_SHELL_MODES || form.len() != n*n || form.iter().any(|v| !v.is_finite()) {
        return Err(invalid("regional bending factor needs a finite complete bounded square form"));
    }
    let active: Vec<_> = (0..n).filter(|&i| form[n*i..n*(i+1)].iter().any(|v| *v != 0.0)).collect();
    let size = active.len();
    if size == 0 { return Ok(Vec::new()); }
    let scales: Vec<_> = active.iter().map(|&i| det::sqrt(form[n*i+i])).collect();
    if scales.iter().any(|s| !s.is_finite() || *s <= 0.0) {
        return Err(invalid("regional bending memory needs a resolved positive strain subspace"));
    }
    let normalized: Vec<_> = (0..size*size).map(|i|
        form[n*active[i/size]+active[i%size]]/scales[i/size]/scales[i%size]).collect();
    let factors = match fs_la::factor::cholesky(&normalized,size) {
        // Preserve successful positive-definite projections and their channel
        // order exactly. A material patch need not be positive definite on the
        // complete instrument's modes, however: use the SAME pivoted PSD owner
        // as localized interface loss if the unpivoted factor refuses.
        Ok(root) => (0..size).map(|channel|
            (0..size).map(|i| root.l(i,channel)).collect::<Vec<_>>()).collect(),
        Err(_) => super::ModalPlateViscosity::from_psd_modal_matrix(
            normalized, size, 1.0, gate,
        )?.columns,
    };
    let mut rows = Vec::with_capacity(factors.len());
    for factor in factors {
        poll(gate)?;
        let mut row = vec![0.0;n];
        for (i,&mode) in active.iter().enumerate() { row[mode] = scales[i]*factor[i]; }
        if row.iter().any(|v| !v.is_finite()) || row.iter().all(|v| *v == 0.0) {
            return Err(invalid("regional material strain factor is unrepresentable"));
        }
        rows.push(row);
    }
    for i in 0..n { for j in 0..n {
        poll(gate)?;
        let actual: f64 = rows.iter().map(|row| row[i]*row[j]).sum();
        let scale = det::sqrt(form[n*i+i])*det::sqrt(form[n*j+j]);
        if !actual.is_finite() || !scale.is_finite() || (actual-form[n*i+j]).abs() > 1e-10*scale {
            return Err(invalid("regional material channels do not reconstruct the original bending form"));
        }
    }}
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g0_localized_strain_retains_its_nullspace_and_complete_scaled_work() {
        let gate = CancelGate::new_clock_free();
        // One physical strain channel on three coordinates, including a
        // completely unstrained coordinate and disparate coordinate scales.
        let form = [4e12, -6e6, 0.0, -6e6, 9.0, 0.0, 0.0, 0.0, 0.0];
        let rows = factor(&form, 3, &gate).unwrap();
        assert_eq!(rows.len(), 1);
        for q in [[3e-6_f64, 2.0, 7.0], [1e-6, -1.0, 0.0]] {
            let expected = 0.5 * (2e6 * q[0] - 3.0 * q[1]).powi(2);
            let actual: f64 = rows.iter().map(|row| {
                0.5 * row.iter().zip(q).map(|(c, x)| c * x).sum::<f64>().powi(2)
            }).sum();
            assert!((actual - expected).abs() <= 1e-12 * expected.max(1.0));
        }
        assert!(rows.iter().all(|row| row[2] == 0.0));
        assert!(factor(&[1.0, 2.0, 2.0, 1.0], 2, &gate).is_err());
        assert!(factor(&[0.0, 1.0, 1.0, 0.0], 2, &gate).is_err());
    }

    #[test]
    fn g0_positive_definite_regional_channels_keep_original_arithmetic() {
        let gate = CancelGate::new_clock_free();
        let form = [4.0, 1.0, 1.0, 9.0];
        let scale = [det::sqrt(4.0), det::sqrt(9.0)];
        let normalized: Vec<_> = (0..4).map(|i| form[i] / scale[i / 2] / scale[i % 2]).collect();
        let original = fs_la::factor::cholesky(&normalized, 2).unwrap();
        let expected: Vec<Vec<_>> = (0..2).map(|j|
            (0..2).map(|i| scale[i] * original.l(i, j)).collect()).collect();
        assert_eq!(factor(&form, 2, &gate).unwrap(), expected);
    }
}
