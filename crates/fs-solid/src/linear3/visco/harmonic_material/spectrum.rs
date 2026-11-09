//! Coherent finite-interval material energy on the original solid quadrature.

use super::{ViscoHarmonicMaterialBudget, cap};
use crate::linear3::visco::{ViscoElementError, ViscoElements, ViscoHarmonicError, checkpoint};
use fs_exec::Cx;
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_material::visco::{
    MAX_HARMONIC_MATERIAL_LINES, TensorPronyHarmonicSpectrum, TensorPronyHarmonicWindowPoint,
    TensorPronySpectrumWindow, ViscoError,
};
use fs_math::{c64::C64, dd::Dd};

/// All authored strain tones at one original constitutive point. The material
/// account includes their interference over the actual observation interval.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscoHarmonicMaterialSpectrumPoint {
    pub element: usize,
    pub quadrature_point: usize,
    /// Positive original reference-volume weight [m3].
    pub volume_m3: f64,
    /// Peak exp(-i omega t) Mandel strains, in the original input tone order.
    pub strain_mandel: Vec<[C64; 6]>,
    /// Original compatible strains when the constitutive field uses B-bar.
    /// The complete compatible waveform independently satisfies the source limit.
    pub compatible_strain_mandel: Option<Vec<[C64; 6]>>,
    /// Physical storage, work and Kelvin/Maxwell heat densities, including
    /// equal-frequency interference and all sum/difference-frequency terms.
    pub material: TensorPronyHarmonicWindowPoint,
}

/// Complete original-quadrature material energy for a coherent line spectrum.
/// Different frequencies need not share a common period. Every phasor uses the
/// same clock; unknown phases are not interpreted as independent random phases.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscoHarmonicMaterialSpectrumField {
    pub window: HarmonicWindow,
    /// The unchanged positive angular frequencies [rad/s], in input order.
    pub angular_frequencies_rad_s: Vec<f64>,
    pub points: Vec<ViscoHarmonicMaterialSpectrumPoint>,
    pub start_stored_energy_j: f64,
    pub end_stored_energy_j: f64,
    /// Sum of direct pointwise increments; rounded endpoint totals are never
    /// subtracted to recover a small physical storage change.
    pub stored_energy_change_j: f64,
    pub kelvin_heat_j: f64,
    pub maxwell_heat_j: f64,
    pub material_heat_j: f64,
    pub stress_work_j: f64,
    /// Independently computed stress work minus storage change minus heat.
    pub energy_balance_residual_j: f64,
}

/// Immutable original material spectra and free/support maps. This observes
/// supplied free-coordinate displacement phasors from the existing mechanical
/// solvers; it does not refit the source or replace its hereditary constitutive law.
pub struct ViscoHarmonicMaterialSpectrumObserver<'a> {
    source: &'a ViscoElements,
    angular_frequencies_rad_s: Vec<f64>,
    materials: Vec<TensorPronyHarmonicSpectrum<'a>>,
    local_maps: Vec<Vec<Option<usize>>>,
    points: usize,
    local_extent: usize,
}

impl ViscoElements {
    /// Prepare coherent multitone material observation on every original point.
    ///
    /// Each frequency is admitted by each original material law. Fixed supports
    /// remain zero, matching the harmonic mechanical owners. Coherent repeated
    /// frequencies combine inside the material owner before waveform admission;
    /// distinct-frequency peak envelopes conservatively admit the entire motion.
    ///
    /// The existing point/payload/work budget covers preparation and one complete
    /// observation. Costs use the AUTHORED tone count F, including duplicates:
    /// work <= 262144*F^2*sum(points*(1+branches)) +
    /// 8192*F*sum(1+branches). No frequency, point or branch is dropped to fit.
    pub fn prepare_harmonic_material_spectrum_observer(
        &self,
        angular_frequencies_rad_s: &[f64],
        budget: ViscoHarmonicMaterialBudget,
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicMaterialSpectrumObserver<'_>, ViscoHarmonicError> {
        checkpoint(cx)?;
        if angular_frequencies_rad_s.is_empty()
            || angular_frequencies_rad_s.len() > MAX_HARMONIC_MATERIAL_LINES
            || angular_frequencies_rad_s
                .iter()
                .any(|omega| !omega.is_finite() || *omega <= 0.0)
        {
            return Err(ViscoElementError::Material(ViscoError::Parameters {
                what: "coherent material observation requires 1..=64 positive finite frequencies",
            })
            .into());
        }
        let tones = angular_frequencies_rad_s.len() as u128;
        let (mut points, mut cache_bytes, mut field_bytes, mut work) =
            (0_u128, 0_u128, 0_u128, 0_u128);
        let mut local_extent = 0;
        for element in &self.elements {
            checkpoint(cx)?;
            let count = element.points.len() as u128;
            let branches = element.law.branches().len() as u128;
            points += count;
            local_extent = local_extent.max(element.dofs.len());
            // Conservative upper bounds include spectrum vector headers,
            // original source weights, tone grouping and temporal pair caches.
            cache_bytes += core::mem::size_of::<TensorPronyHarmonicSpectrum<'_>>() as u128
                + 128 * tones * (1 + branches)
                + core::mem::size_of::<Vec<Option<usize>>>() as u128
                + element.dofs.len() as u128 * core::mem::size_of::<Option<usize>>() as u128;
            for point in &element.points {
                let strain_fields = if point.compatible_gradients.is_some() {
                    2
                } else {
                    1
                };
                field_bytes += core::mem::size_of::<ViscoHarmonicMaterialSpectrumPoint>() as u128
                    + strain_fields * tones * core::mem::size_of::<[C64; 6]>() as u128;
            }
            work += 262144 * tones * tones * count * (1 + branches) + 8192 * tones * (1 + branches);
        }
        // One temporal cache is live at a time, while all original material
        // spectra, the returned field and one local nodal gather remain live.
        let bytes = cache_bytes
            + field_bytes
            + 2 * tones * core::mem::size_of::<f64>() as u128
            + local_extent as u128 * tones * core::mem::size_of::<C64>() as u128
            + core::mem::size_of::<TensorPronySpectrumWindow<'_, '_>>() as u128
            + 256 * tones
            + 64 * tones * tones;
        cap(
            "harmonic material spectrum quadrature points",
            points,
            budget.maximum_points as u128,
        )?;
        cap(
            "harmonic material spectrum numerical payload bytes",
            bytes,
            budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "harmonic material spectrum observation work",
            work,
            budget.maximum_work as u128,
        )?;
        let mut materials = Vec::with_capacity(self.elements.len());
        let mut local_maps = Vec::with_capacity(self.elements.len());
        for element in &self.elements {
            checkpoint(cx)?;
            materials.push(
                element
                    .law
                    .prepare_harmonic_spectrum(angular_frequencies_rad_s)
                    .map_err(ViscoElementError::Material)?,
            );
            local_maps.push(
                element
                    .dofs
                    .iter()
                    .map(|dof| self.equilibrium.free_dofs.binary_search(dof).ok())
                    .collect(),
            );
        }
        checkpoint(cx)?;
        Ok(ViscoHarmonicMaterialSpectrumObserver {
            source: self,
            angular_frequencies_rad_s: angular_frequencies_rad_s.to_vec(),
            materials,
            local_maps,
            points: points as usize,
            local_extent,
        })
    }
}

impl<'a> ViscoHarmonicMaterialSpectrumObserver<'a> {
    /// Exact immutable geometry, source laws and support map being observed.
    pub const fn source(&self) -> &'a ViscoElements {
        self.source
    }

    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        &self.angular_frequencies_rad_s
    }

    /// Integrate physical material work/heat over an unweighted time interval.
    /// `displacements[tone][free_coordinate]` uses peak exp(-i omega t), with
    /// every tone sharing the source clock and original free-coordinate map.
    ///
    /// One temporal cache per element and one nodal buffer serve all points.
    /// This is the settled linear particular response, without startup memory,
    /// changing temperature, DC preload, or a new mechanical equilibrium claim.
    /// Any source, budget or cancellation refusal publishes no partial field.
    pub fn observe_window(
        &self,
        displacements: &[&[C64]],
        window: HarmonicWindow,
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicMaterialSpectrumField, ViscoHarmonicError> {
        checkpoint(cx)?;
        if window.profile() != HarmonicWindowProfile::Uniform {
            return Err(ViscoElementError::Material(ViscoError::Parameters {
                what: "physical harmonic material heat requires a uniform time interval",
            })
            .into());
        }
        let tones = self.angular_frequencies_rad_s.len();
        if displacements.len() != tones
            || displacements.iter().any(|values| {
                values.len() != self.source.equilibrium.free_dofs.len()
                    || values
                        .iter()
                        .any(|z| !z.re.is_finite() || !z.im.is_finite())
            })
        {
            return Err(ViscoElementError::InvalidHistory {
                what: "coherent material observation needs one finite free-coordinate vector per tone",
            }
            .into());
        }
        let mut points = Vec::with_capacity(self.points);
        let mut local = Vec::with_capacity(tones * self.local_extent);
        let mut totals = [Dd::ZERO; 7];
        for (element, ((source, material), map)) in self
            .source
            .elements
            .iter()
            .zip(&self.materials)
            .zip(&self.local_maps)
            .enumerate()
        {
            checkpoint(cx)?;
            let prepared = material
                .prepare_window(window)
                .map_err(ViscoElementError::Material)?;
            local.clear();
            for displacement in displacements {
                checkpoint(cx)?;
                local.extend(
                    map.iter()
                        .map(|index| index.map_or(C64::ZERO, |i| displacement[i])),
                );
            }
            for (quadrature_point, point) in source.points.iter().enumerate() {
                checkpoint(cx)?;
                let mut strain = Vec::with_capacity(tones);
                let mut compatible = point
                    .compatible_gradients
                    .as_ref()
                    .map(|_| Vec::with_capacity(tones));
                for tone in 0..tones {
                    let nodal = &local[tone * map.len()..(tone + 1) * map.len()];
                    let re = point.strain_from_local(|i| nodal[i].re);
                    let im = point.strain_from_local(|i| nodal[i].im);
                    strain.push(core::array::from_fn(|i| C64::new(re[i], im[i])));
                    if let Some(compatible) = &mut compatible {
                        let re = point.compatible_strain_from_local(re, |i| nodal[i].re);
                        let im = point.compatible_strain_from_local(im, |i| nodal[i].im);
                        compatible.push(core::array::from_fn(|i| C64::new(re[i], im[i])));
                    }
                }
                if let Some(compatible) = &compatible {
                    material
                        .validate_strains(compatible)
                        .map_err(ViscoElementError::Material)?;
                }
                let response = prepared
                    .evaluate_with_checkpoint::<ViscoElementError>(&strain, || checkpoint(cx))?;
                let volume = Dd::from_f64(point.volume_m3);
                for (total, value) in totals.iter_mut().zip([
                    response.start_stored_energy_j_m3,
                    response.end_stored_energy_j_m3,
                    response.stored_energy_change_j_m3,
                    response.kelvin_heat_j_m3,
                    response.maxwell_heat_j_m3,
                    response.material_heat_j_m3,
                    response.stress_work_j_m3,
                ]) {
                    *total = *total + volume * Dd::from_f64(value);
                }
                points.push(ViscoHarmonicMaterialSpectrumPoint {
                    element,
                    quadrature_point,
                    volume_m3: point.volume_m3,
                    strain_mandel: strain,
                    compatible_strain_mandel: compatible,
                    material: response,
                });
            }
        }
        let balance = (totals[6] - totals[2] - totals[5]).to_f64();
        let [start, end, change, kelvin, maxwell, heat, work] = totals.map(Dd::to_f64);
        if [start, end, change, kelvin, maxwell, heat, work, balance]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(ViscoElementError::NonFiniteResult.into());
        }
        checkpoint(cx)?;
        Ok(ViscoHarmonicMaterialSpectrumField {
            window,
            angular_frequencies_rad_s: self.angular_frequencies_rad_s.clone(),
            points,
            start_stored_energy_j: start,
            end_stored_energy_j: end,
            stored_energy_change_j: change,
            kelvin_heat_j: kelvin,
            maxwell_heat_j: maxwell,
            material_heat_j: heat,
            stress_work_j: work,
            energy_balance_residual_j: balance,
        })
    }
}
