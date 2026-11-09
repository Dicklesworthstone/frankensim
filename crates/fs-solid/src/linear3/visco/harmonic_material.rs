//! Full original-quadrature stress and physical material energy in harmonic motion.

use super::{ViscoElementError, ViscoElements, ViscoHarmonicError, checkpoint};
use fs_exec::Cx;
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_material::visco::{
    TensorPronyHarmonic, TensorPronyHarmonicPoint, TensorPronyHarmonicWindowPoint, ViscoError,
};
use fs_math::{c64::C64, dd::Dd};
use fs_modal::harmonic::HarmonicResponseError;

mod spectrum;
pub use spectrum::{
    ViscoHarmonicMaterialSpectrumField, ViscoHarmonicMaterialSpectrumObserver,
    ViscoHarmonicMaterialSpectrumPoint,
};
mod causal;
pub use causal::{
    ViscoHarmonicMaterialCausalField, ViscoHarmonicMaterialCausalObserver,
    ViscoHarmonicMaterialCausalStep, ViscoHarmonicMaterialHistory,
};

/// Cold cache plus returned field admission. The original geometry/laws and
/// displacement input are borrowed, not included in this additional payload.
#[derive(Clone, Copy, Debug)]
pub struct ViscoHarmonicMaterialBudget {
    /// Original constitutive quadrature points; none may be subsampled.
    pub maximum_points: usize,
    /// Conservative heap numerical payload for the cache, local scratch and
    /// complete returned field, excluding fixed stack and allocator metadata.
    pub maximum_payload_bytes: usize,
    /// Conservative preparation plus one complete observation work allowance.
    /// Cycle means require 65536*sum(points*(1+branches)) +
    /// 4096*sum(1+branches). A finite-window observation checks its own larger
    /// allowance: 262144*sum(points*(1+branches)) + 8192*sum(1+branches).
    /// An observer admitted only for cycle means remains usable for those means.
    pub maximum_work: u64,
}

/// One original element/quadrature point, in the original global Mandel frame.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscoHarmonicMaterialPoint {
    pub element: usize,
    pub quadrature_point: usize,
    /// Positive original constitutive reference-volume weight [m3].
    pub volume_m3: f64,
    /// Actual constitutive strain used by original K'/K'' assembly.
    pub strain_mandel: [C64; 6],
    /// Original compatible strain when the constitutive field uses B-bar.
    /// Both fields independently satisfy the source's complete-cycle limit.
    pub compatible_strain_mandel: Option<[C64; 6]>,
    pub material: TensorPronyHarmonicPoint,
}

/// Original physical field and independently integrated material totals.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscoHarmonicMaterialField {
    pub points: Vec<ViscoHarmonicMaterialPoint>,
    pub mean_stored_energy_j: f64,
    pub mean_material_heat_w: f64,
    pub mean_stress_power_w: f64,
    pub power_balance_residual_w: f64,
}

/// Physical finite-interval material response at one original quadrature point.
/// All strain phasors use the same peak exp(-i omega t) source clock.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscoHarmonicMaterialWindowPoint {
    pub element: usize,
    pub quadrature_point: usize,
    pub volume_m3: f64,
    /// The strain used by the original constitutive assembly.
    pub strain_mandel: [C64; 6],
    /// Original compatible strain, retained and independently admitted for B-bar.
    pub compatible_strain_mandel: Option<[C64; 6]>,
    pub material: TensorPronyHarmonicWindowPoint,
}

/// Independently integrated physical storage, work and heat over one interval.
/// These are total energies [J], not cycle means or weighted detector exposures.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscoHarmonicMaterialWindowField {
    pub window: HarmonicWindow,
    pub points: Vec<ViscoHarmonicMaterialWindowPoint>,
    pub start_stored_energy_j: f64,
    pub end_stored_energy_j: f64,
    /// Sum of direct pointwise increments, retaining changes below the precision
    /// of subtracting the separately rounded total endpoint energies.
    pub stored_energy_change_j: f64,
    pub kelvin_heat_j: f64,
    pub maxwell_heat_j: f64,
    pub material_heat_j: f64,
    pub stress_work_j: f64,
    /// Signed stress work minus storage change minus physical material heat.
    pub energy_balance_residual_j: f64,
}

/// Immutable per-element material spectra and original free/support maps.
/// This observes a supplied displacement phasor; it does not certify that the
/// displacement solves an external force problem. It can consume either the
/// condensed or complete-memory harmonic solver without another mechanical solve.
pub struct ViscoHarmonicMaterialObserver<'a> {
    source: &'a ViscoElements,
    omega: f64,
    materials: Vec<TensorPronyHarmonic<'a>>,
    local_maps: Vec<Vec<Option<usize>>>,
    points: usize,
    local_extent: usize,
    budget: ViscoHarmonicMaterialBudget,
    window_payload_bytes: u128,
    window_work: u128,
}

fn cap(what: &'static str, requested: u128, maximum: u128) -> Result<(), ViscoHarmonicError> {
    if requested > maximum {
        return Err(HarmonicResponseError::Budget {
            quantity: what,
            requested: usize::try_from(requested).unwrap_or(usize::MAX),
            cap: usize::try_from(maximum).unwrap_or(usize::MAX),
        }
        .into());
    }
    Ok(())
}

impl ViscoElements {
    /// Prepare complete stress/heat observation from the same original strain
    /// operators, material laws, supports and reference quadrature as the solid.
    /// Frequency and material domain admission remain source-owned. Preparation
    /// caches free/support lookup and spectra once for repeated excitation cases.
    pub fn prepare_harmonic_material_observer(
        &self,
        omega: f64,
        budget: ViscoHarmonicMaterialBudget,
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicMaterialObserver<'_>, ViscoHarmonicError> {
        checkpoint(cx)?;
        let (mut points, mut bytes, mut work) = (0_u128, 0_u128, 0_u128);
        let mut window_work = 0_u128;
        let mut local_extent = 0;
        for element in &self.elements {
            checkpoint(cx)?;
            let branches = element.law.branches().len() as u128;
            let count = element.points.len() as u128;
            points += count;
            local_extent = local_extent.max(element.dofs.len());
            bytes += core::mem::size_of::<TensorPronyHarmonic<'_>>() as u128
                + 128 * branches // all cached binary spectrum weights, including padding
                + core::mem::size_of::<Vec<Option<usize>>>() as u128
                + element.dofs.len() as u128 * core::mem::size_of::<Option<usize>>() as u128;
            work += 65536 * count * (1 + branches) + 4096 * (1 + branches);
            // Two endpoint stores, direct storage increments, physical heat,
            // and independent stress work use the same original source law.
            window_work += 262144 * count * (1 + branches) + 8192 * (1 + branches);
        }
        let window_payload_bytes = bytes
            + points * core::mem::size_of::<ViscoHarmonicMaterialWindowPoint>() as u128
            + local_extent as u128 * core::mem::size_of::<C64>() as u128;
        bytes += points * core::mem::size_of::<ViscoHarmonicMaterialPoint>() as u128
            + local_extent as u128 * core::mem::size_of::<C64>() as u128;
        cap(
            "harmonic material quadrature points",
            points,
            budget.maximum_points as u128,
        )?;
        cap(
            "harmonic material numerical payload bytes",
            bytes,
            budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "harmonic material observation work",
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
                    .prepare_harmonic_material(omega)
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
        Ok(ViscoHarmonicMaterialObserver {
            source: self,
            omega,
            materials,
            local_maps,
            points: points as usize,
            local_extent,
            budget,
            window_payload_bytes,
            window_work,
        })
    }
}

impl<'a> ViscoHarmonicMaterialObserver<'a> {
    /// Immutable original source; no law or geometry is rebound during observation.
    pub const fn source(&self) -> &'a ViscoElements {
        self.source
    }
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }

    /// Observe all original points from peak displacement exp(-i omega t).
    /// Homogeneous prescribed supports are zero, as in the harmonic solvers.
    /// Per-point constitutive evaluation allocates nothing; the one returned
    /// field and one reusable local nodal buffer retain the prepared budget.
    /// Cancellation/domain refusal returns no partial field or material history.
    pub fn observe(
        &self,
        displacement: &[C64],
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicMaterialField, ViscoHarmonicError> {
        checkpoint(cx)?;
        if displacement.len() != self.source.equilibrium.free_dofs.len()
            || displacement
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(ViscoElementError::InvalidHistory {
                what: "harmonic observation needs one finite phasor per original free coordinate",
            }
            .into());
        }
        let mut points = Vec::with_capacity(self.points);
        let mut local = Vec::with_capacity(self.local_extent);
        let (mut energy, mut heat, mut work) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
        for (element, ((source, material), map)) in self
            .source
            .elements
            .iter()
            .zip(&self.materials)
            .zip(&self.local_maps)
            .enumerate()
        {
            checkpoint(cx)?;
            local.clear();
            local.extend(
                map.iter()
                    .map(|index| index.map_or(C64::ZERO, |i| displacement[i])),
            );
            for (quadrature_point, point) in source.points.iter().enumerate() {
                checkpoint(cx)?;
                let re = point.strain_from_local(|i| local[i].re);
                let im = point.strain_from_local(|i| local[i].im);
                let strain = core::array::from_fn(|i| C64::new(re[i], im[i]));
                let compatible = if point.compatible_gradients.is_some() {
                    let re = point.compatible_strain_from_local(re, |i| local[i].re);
                    let im = point.compatible_strain_from_local(im, |i| local[i].im);
                    let strain = core::array::from_fn(|i| C64::new(re[i], im[i]));
                    source
                        .law
                        .validate_harmonic_mandel_strain(&strain)
                        .map_err(ViscoElementError::Material)?;
                    Some(strain)
                } else {
                    None
                };
                let response = material
                    .evaluate(&strain)
                    .map_err(ViscoElementError::Material)?;
                let volume = Dd::from_f64(point.volume_m3);
                energy = energy + volume * Dd::from_f64(response.mean_stored_energy_j_m3);
                heat = heat + volume * Dd::from_f64(response.mean_material_heat_w_m3);
                work = work + volume * Dd::from_f64(response.mean_stress_power_w_m3);
                points.push(ViscoHarmonicMaterialPoint {
                    element,
                    quadrature_point,
                    volume_m3: point.volume_m3,
                    strain_mandel: strain,
                    compatible_strain_mandel: compatible,
                    material: response,
                });
            }
        }
        let balance = (work - heat).to_f64();
        let (energy, heat, work) = (energy.to_f64(), heat.to_f64(), work.to_f64());
        if [energy, heat, work, balance].iter().any(|v| !v.is_finite()) {
            return Err(ViscoElementError::NonFiniteResult.into());
        }
        checkpoint(cx)?;
        Ok(ViscoHarmonicMaterialField {
            points,
            mean_stored_energy_j: energy,
            mean_material_heat_w: heat,
            mean_stress_power_w: work,
            power_balance_residual_w: balance,
        })
    }

    /// Observe unweighted physical heat, work and endpoint storage over `window`.
    /// Every original quadrature point and compatible B-bar strain is retained.
    /// One allocation-free material window is prepared per element and reused
    /// for all its points; the result field and one local nodal buffer use the
    /// originally supplied payload cap. The larger window work allowance is
    /// checked before allocation, without changing ordinary cycle-mean admission.
    ///
    /// Only uniform intervals describe physical heat and endpoint conservation.
    /// This observes the periodic particular displacement, with source material
    /// conditions fixed. It neither advances startup histories nor temperature.
    /// Refusal or cancellation returns no partial field and changes no source.
    pub fn observe_window(
        &self,
        displacement: &[C64],
        window: HarmonicWindow,
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicMaterialWindowField, ViscoHarmonicError> {
        checkpoint(cx)?;
        cap(
            "harmonic material window numerical payload bytes",
            self.window_payload_bytes,
            self.budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "harmonic material window observation work",
            self.window_work,
            self.budget.maximum_work as u128,
        )?;
        if window.profile() != HarmonicWindowProfile::Uniform {
            return Err(ViscoElementError::Material(ViscoError::Parameters {
                what: "physical harmonic material heat requires a uniform time interval",
            })
            .into());
        }
        if displacement.len() != self.source.equilibrium.free_dofs.len()
            || displacement
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(ViscoElementError::InvalidHistory {
                what: "harmonic observation needs one finite phasor per original free coordinate",
            }
            .into());
        }
        let mut points = Vec::with_capacity(self.points);
        let mut local = Vec::with_capacity(self.local_extent);
        // Independent total storage, direct increment, physical heat and work.
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
            local.extend(
                map.iter()
                    .map(|index| index.map_or(C64::ZERO, |i| displacement[i])),
            );
            for (quadrature_point, point) in source.points.iter().enumerate() {
                checkpoint(cx)?;
                let re = point.strain_from_local(|i| local[i].re);
                let im = point.strain_from_local(|i| local[i].im);
                let strain = core::array::from_fn(|i| C64::new(re[i], im[i]));
                let compatible = if point.compatible_gradients.is_some() {
                    let re = point.compatible_strain_from_local(re, |i| local[i].re);
                    let im = point.compatible_strain_from_local(im, |i| local[i].im);
                    let strain = core::array::from_fn(|i| C64::new(re[i], im[i]));
                    source
                        .law
                        .validate_harmonic_mandel_strain(&strain)
                        .map_err(ViscoElementError::Material)?;
                    Some(strain)
                } else {
                    None
                };
                let response = prepared
                    .evaluate(&strain)
                    .map_err(ViscoElementError::Material)?;
                let values = [
                    response.start_stored_energy_j_m3,
                    response.end_stored_energy_j_m3,
                    response.stored_energy_change_j_m3,
                    response.kelvin_heat_j_m3,
                    response.maxwell_heat_j_m3,
                    response.material_heat_j_m3,
                    response.stress_work_j_m3,
                ];
                let volume = Dd::from_f64(point.volume_m3);
                for (total, value) in totals.iter_mut().zip(values) {
                    *total = *total + volume * Dd::from_f64(value);
                }
                points.push(ViscoHarmonicMaterialWindowPoint {
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
        Ok(ViscoHarmonicMaterialWindowField {
            window,
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
