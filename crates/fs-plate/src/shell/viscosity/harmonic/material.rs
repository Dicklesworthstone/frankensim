//! Original layer stress, spring energy and physical heat on the CST/DKT field.

use super::{
    C64, Cx, Dd, HarmonicResponseError, PreparedTensorPronyShell, TensorPronyShell,
    TensorShellHarmonicError, checkpoint, compensated_products, invalid, strain_rows,
};
use fs_material::visco::{TensorPronyHarmonic, TensorPronyHarmonicPoint};

mod causal;
pub use causal::{
    TensorShellMaterialCausalBudget, TensorShellMaterialCausalField,
    TensorShellMaterialCausalHistory, TensorShellMaterialCausalObserver,
    TensorShellMaterialCausalPoint, TensorShellMaterialCausalStep,
    TensorShellViscoelasticCausalField, TensorShellViscoelasticCausalObserver,
    TensorShellViscoelasticCausalPoint, TensorShellViscoelasticCausalStep,
    TensorShellViscoelasticEnergyDensity,
};
mod spectrum;
pub use spectrum::{
    TensorShellMaterialSpectrumField, TensorShellMaterialSpectrumObserver,
    TensorShellMaterialSpectrumPoint,
};

/// Additional material cache and one complete returned field. Original source
/// geometry, the prepared mechanical factor and fixed-size stack are excluded.
#[derive(Clone, Copy, Debug)]
pub struct TensorShellMaterialBudget {
    /// Six positive quadrature points per original physical layer and facet.
    pub maximum_points: usize,
    /// Additional heap numerical payload, excluding allocator metadata.
    pub maximum_payload_bytes: usize,
    /// Single-frequency mean-field preparation plus one observation score:
    /// 4096*sum(1+branches) + 6*65536*sum(1+branches).
    /// Coherent spectra use the authored-line bound documented by
    /// `TensorPronyShell::prepare_material_spectrum_observer`.
    pub maximum_work: u64,
}

/// One original triangle/layer integration point in that layer's material axes.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorShellMaterialPoint {
    pub element: usize,
    pub layer: usize,
    pub triangle_quadrature_point: usize,
    pub thickness_quadrature_point: usize,
    /// Original triangle barycentric coordinates, ordered 01, 12, 20.
    pub triangle_barycentric: [f64; 3],
    /// Signed thickness coordinate from the declared reference surface [m].
    pub thickness_position_m: f64,
    /// Positive original reference volume of this quadrature point [m3].
    pub volume_m3: f64,
    /// Full recovered Mandel strain in ORIGINAL layer material axes.
    pub strain_mandel: [C64; 6],
    /// Full stress, spring storage and independent physical heat at this point.
    pub material: TensorPronyHarmonicPoint,
}

/// Complete original layer field and independently integrated physical totals.
/// Numerical drilling, kinetic energy and boundary work are not material energy.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorShellMaterialField {
    pub points: Vec<TensorShellMaterialPoint>,
    pub mean_equilibrium_energy_j: f64,
    pub mean_maxwell_energy_j: f64,
    pub mean_stored_energy_j: f64,
    pub mean_kelvin_heat_w: f64,
    pub mean_maxwell_heat_w: f64,
    pub mean_material_heat_w: f64,
    pub mean_stress_power_w: f64,
    pub power_balance_residual_w: f64,
}

/// Original geometry and support contraction, shared by repeated harmonic solves
/// and material observation. All CST/DKT rows come from the existing owner.
pub(super) struct PreparedFacet {
    frame: [[f64; 3]; 3],
    coordinates: [Option<usize>; 18],
    rows: [[[f64; 18]; 6]; 3],
    area_m2: f64,
}
impl PreparedFacet {
    pub(super) fn new(source: &TensorPronyShell, element: usize)
        -> Result<Self, TensorShellHarmonicError>
    {
        let g = source.mesh.facet(element)?;
        let tri = source.mesh.tris[element];
        Ok(Self {
            frame: g.frame,
            coordinates: core::array::from_fn(|i| source.model.dof_map[6 * tri[i / 6] + i % 6]),
            rows: strain_rows(&g, element)?,
            area_m2: g.area_m2,
        })
    }

    pub(super) fn strains(&self, displacement: &[C64]) -> [[C64; 6]; 3] {
        self.strains_by(|r| displacement[r])
    }

    pub(super) fn strains_by(&self, mut displacement: impl FnMut(usize) -> C64) -> [[C64; 6]; 3] {
        let local: [C64; 18] = core::array::from_fn(|i| {
            let base = 6 * (i / 6) + 3 * ((i % 6) / 3);
            let terms: [(f64, C64); 3] = core::array::from_fn(|p| (
                self.frame[i % 3][p],
                self.coordinates[base + p].map_or(C64::ZERO, |r| displacement(r)),
            ));
            C64::new(
                compensated_products(terms.iter().map(|(a, b)| (*a, b.re))),
                compensated_products(terms.iter().map(|(a, b)| (*a, b.im))),
            )
        });
        self.rows.map(|point| core::array::from_fn(|r| C64::new(
            compensated_products((0..18).map(|i| (point[r][i], local[i].re))),
            compensated_products((0..18).map(|i| (point[r][i], local[i].im))),
        )))
    }
}

/// Borrows one prepared mechanical frequency and all its original layer maps.
/// Preparation adds only original constitutive spectra. Repeated observations
/// reuse geometry, supports, plane-stress recovery and material factors.
pub struct TensorShellMaterialObserver<'prepared, 'source> {
    prepared: &'prepared PreparedTensorPronyShell<'source>,
    materials: Vec<Vec<TensorPronyHarmonic<'source>>>,
    points: usize,
}

fn cap(quantity: &'static str, requested: u128, maximum: u128)
    -> Result<(), TensorShellHarmonicError>
{
    if requested > maximum {
        return Err(HarmonicResponseError::Budget {
            quantity,
            requested: usize::try_from(requested).unwrap_or(usize::MAX),
            cap: usize::try_from(maximum).unwrap_or(usize::MAX),
        }.into());
    }
    Ok(())
}

impl<'source> PreparedTensorPronyShell<'source> {
    /// Prepare complete original layer fields at this SAME mechanical frequency.
    /// No second mechanical solve/factorization, geometry reconstruction or
    /// material inverse is performed by point observation. The original source
    /// law, orientation, reference surface and frequency domain stay fixed.
    pub fn prepare_material_observer<'prepared>(
        &'prepared self,
        budget: TensorShellMaterialBudget,
        cx: &Cx<'_>,
    ) -> Result<TensorShellMaterialObserver<'prepared, 'source>, TensorShellHarmonicError> {
        checkpoint(cx)?;
        let source: &'source TensorPronyShell = self.source;
        let (mut points, mut bytes, mut work) = (0_u128, 0_u128, 0_u128);
        for section in &source.sections {
            checkpoint(cx)?;
            bytes += core::mem::size_of::<Vec<TensorPronyHarmonic<'source>>>() as u128;
            for layer in 0..section.layer_count() {
                let branches = section.layer_law(layer)?.branches().len() as u128;
                points += 6;
                bytes += core::mem::size_of::<TensorPronyHarmonic<'source>>() as u128
                    + 128 * branches;
                work += (4096 + 6 * 65536) * (1 + branches);
            }
        }
        bytes += points * core::mem::size_of::<TensorShellMaterialPoint>() as u128;
        cap("shell material field points", points, budget.maximum_points as u128)?;
        cap("shell material field payload bytes", bytes, budget.maximum_payload_bytes as u128)?;
        cap("shell material field work", work, budget.maximum_work as u128)?;
        let mut materials = Vec::with_capacity(source.sections.len());
        for (element, section) in source.sections.iter().enumerate() {
            checkpoint(cx)?;
            let mut layers = Vec::with_capacity(section.layer_count());
            for layer in 0..section.layer_count() {
                checkpoint(cx)?;
                layers.push(section.layer_law(layer)?
                    .prepare_harmonic_material(self.omega_rad_s)
                    .map_err(|error| section.material_error(element, layer, error))?);
            }
            materials.push(layers);
        }
        checkpoint(cx)?;
        Ok(TensorShellMaterialObserver { prepared: self, materials, points: points as usize })
    }
}

impl<'prepared, 'source> TensorShellMaterialObserver<'prepared, 'source> {
    pub const fn source(&self) -> &'source TensorPronyShell {
        self.prepared.source
    }
    pub const fn omega_rad_s(&self) -> f64 {
        self.prepared.omega_rad_s
    }

    /// Observe peak physical displacement/rotation phasors under exp(-i omega t).
    /// The original midpoint triangle rule and two positive Gauss points through
    /// each layer integrate this affine CST/DKT field's quadratic energy exactly
    /// in real arithmetic. Every layer face at every facet vertex independently
    /// admits the complete strain cycle BEFORE its integration points are emitted.
    ///
    /// This observes the supplied field without certifying its force equation.
    /// No startup history, temperature evolution, transverse-traction closure or
    /// spatial/model error estimate is inferred. Refusal returns no partial field.
    pub fn observe(&self, displacement: &[C64], cx: &Cx<'_>)
        -> Result<TensorShellMaterialField, TensorShellHarmonicError>
    {
        checkpoint(cx)?;
        if displacement.len() != self.prepared.source.model.free
            || displacement.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
        {
            return Err(invalid("one finite phasor per original free shell coordinate required"));
        }
        const BARYCENTRIC: [[f64; 3]; 3] =
            [[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]];
        const GAUSS: f64 = 0.5773502691896258;
        let mut points = Vec::with_capacity(self.points);
        // Original equilibrium/Maxwell energy, Kelvin/Maxwell heat, stress work.
        let mut sums = [Dd::ZERO; 5];
        for element in 0..self.materials.len() {
            checkpoint(cx)?;
            let facet = &self.prepared.facets[element];
            let section = &self.prepared.source.sections[element];
            let recovery = &self.prepared.materials[element];
            let midpoint = facet.strains(displacement);
            recovery.admit_facet(section, element, &midpoint, cx)?;
            for (layer, material) in self.materials[element].iter().enumerate() {
                checkpoint(cx)?;
                let [bottom, top] = recovery.layer_faces(layer)
                    .ok_or_else(|| invalid("original shell material layer is missing"))?;
                let half = 0.5 * (top - bottom);
                let center = bottom + half;
                let volume = (facet.area_m2 / 3.0) * half;
                if !half.is_finite() || half <= 0.0 || !volume.is_finite() || volume <= 0.0 {
                    return Err(invalid("original shell layer quadrature volume is unrepresentable"));
                }
                for triangle_quadrature_point in 0..3 {
                    for (thickness_quadrature_point, sign) in [-1.0, 1.0].into_iter().enumerate() {
                        checkpoint(cx)?;
                        let z = center + sign * GAUSS * half;
                        let strain = recovery.recover_strain(
                            layer, &midpoint[triangle_quadrature_point], z,
                        ).map_err(|error| section.material_error(element, layer, error))?;
                        let response = material.evaluate(&strain)
                            .map_err(|error| section.material_error(element, layer, error))?;
                        let channels = [
                            response.mean_equilibrium_energy_j_m3,
                            response.mean_maxwell_energy_j_m3,
                            response.mean_kelvin_heat_w_m3,
                            response.mean_maxwell_heat_w_m3,
                            response.mean_stress_power_w_m3,
                        ];
                        for (sum, value) in sums.iter_mut().zip(channels) {
                            *sum = *sum + Dd::from_f64(volume) * Dd::from_f64(value);
                        }
                        points.push(TensorShellMaterialPoint {
                            element, layer, triangle_quadrature_point, thickness_quadrature_point,
                            triangle_barycentric: BARYCENTRIC[triangle_quadrature_point],
                            thickness_position_m: z,
                            volume_m3: volume,
                            strain_mandel: strain,
                            material: response,
                        });
                    }
                }
            }
        }
        let energy = (sums[0] + sums[1]).to_f64();
        let heat = (sums[2] + sums[3]).to_f64();
        let residual = (sums[4] - sums[2] - sums[3]).to_f64();
        let values = sums.map(|s| s.to_f64());
        if values.iter().chain([energy, heat, residual].iter()).any(|v| !v.is_finite())
            || values[..4].iter().any(|v| *v < 0.0)
        {
            return Err(invalid("original shell material totals are unrepresentable"));
        }
        checkpoint(cx)?;
        Ok(TensorShellMaterialField {
            points,
            mean_equilibrium_energy_j: values[0],
            mean_maxwell_energy_j: values[1],
            mean_stored_energy_j: energy,
            mean_kelvin_heat_w: values[2],
            mean_maxwell_heat_w: values[3],
            mean_material_heat_w: heat,
            mean_stress_power_w: values[4],
            power_balance_residual_w: residual,
        })
    }
}
