//! Finite-interval coherent material heat and work on original shell layers.

use super::super::{
    C64, Cx, Dd, MaterialPlane, TensorPronyShell, TensorShellHarmonicError, checkpoint, invalid,
    section::principal_strain,
};
use super::{PreparedFacet, TensorShellMaterialBudget, cap};
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_material::visco::{
    MAX_HARMONIC_MATERIAL_LINES, PlaneStressHarmonic, PlaneStressHarmonicSpectrum,
    PlaneStressSpectrumWindow, TensorPronyHarmonicWindowPoint, ViscoError,
};

/// One original layer volume with all coherent engineering strain lines.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorShellMaterialSpectrumPoint {
    pub element: usize,
    pub layer: usize,
    pub triangle_quadrature_point: usize,
    pub thickness_quadrature_point: usize,
    pub triangle_barycentric: [f64; 3],
    /// Signed position relative to the original section reference surface [m].
    pub thickness_position_m: f64,
    pub volume_m3: f64,
    /// [epsilon_x, epsilon_y, gamma_xy] in the original layer's ordered axes,
    /// in authored frequency order. The material owner recovers all six complex
    /// Mandel strains with its frequency-specific transverse constraint.
    pub in_plane_strain: Vec<[C64; 3]>,
    /// Relaxed static engineering strain in the same original layer axes.
    /// The original zero-frequency constraint recovers its transverse strains.
    pub in_plane_bias: [f64; 3],
    /// Full original constitutive storage, stress work and Kelvin/Maxwell heat
    /// densities [J/m3], including every temporal interference term.
    pub material: TensorPronyHarmonicWindowPoint,
}

/// Actual interval quantities for the complete original shell material field.
/// Drilling stabilization, kinetic energy and boundary work are separate owners.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorShellMaterialSpectrumField {
    pub window: HarmonicWindow,
    pub angular_frequencies_rad_s: Vec<f64>,
    pub points: Vec<TensorShellMaterialSpectrumPoint>,
    pub start_stored_energy_j: f64,
    pub end_stored_energy_j: f64,
    /// Original direct pointwise increments; never rounded endpoint subtraction.
    pub stored_energy_change_j: f64,
    pub kelvin_heat_j: f64,
    pub maxwell_heat_j: f64,
    pub material_heat_j: f64,
    pub stress_work_j: f64,
    /// Independently integrated stress work minus storage change minus heat.
    pub energy_balance_residual_j: f64,
}

struct LayerSpectrum<'a> {
    material: PlaneStressHarmonicSpectrum<'a>,
    faces_m: [f64; 2],
    to_principal: [[f64; 3]; 3],
}

enum EvaluationError {
    Material(ViscoError),
    Cancelled,
}
impl From<ViscoError> for EvaluationError {
    fn from(value: ViscoError) -> Self {
        Self::Material(value)
    }
}

/// Immutable CST/DKT geometry, support maps and every original layer spectrum.
/// Supplied displacement banks may come from the bare or acoustically loaded
/// owners; preparation neither assembles nor factors a mechanical pencil.
pub struct TensorShellMaterialSpectrumObserver<'a> {
    source: &'a TensorPronyShell,
    angular_frequencies_rad_s: Vec<f64>,
    facets: Vec<PreparedFacet>,
    layers: Vec<Vec<LayerSpectrum<'a>>>,
    points: usize,
}

impl TensorPronyShell {
    /// Prepare coherent observation of all original homogeneous/laminate layers.
    ///
    /// Each authored frequency uses the complete complex plane-stress reduction,
    /// retaining anisotropy, Kelvin viscosity and all original Maxwell branches.
    /// Equal frequencies combine coherently before transverse recovery and source
    /// admission. Frequencies need no common period; an empty bank is zero AC.
    /// Original supports remain stationary, exactly as in the source shell.
    ///
    /// The existing material budget covers retained geometry/spectra, one full
    /// output and temporary temporal/strain buffers. For F authored lines, L
    /// physical layers and E facets, conservative work is bounded by
    /// sum_layers((6*262144*F² + 65536*(F+1))*(1+branches)) + (32768+32768*F)*E.
    /// No layer, point or branch is discarded to satisfy a budget.
    pub fn prepare_material_spectrum_observer(
        &self,
        angular_frequencies_rad_s: &[f64],
        budget: TensorShellMaterialBudget,
        cx: &Cx<'_>,
    ) -> Result<TensorShellMaterialSpectrumObserver<'_>, TensorShellHarmonicError> {
        checkpoint(cx)?;
        if angular_frequencies_rad_s.len() > MAX_HARMONIC_MATERIAL_LINES
            || angular_frequencies_rad_s
                .iter()
                .any(|w| !w.is_finite() || *w <= 0.0)
        {
            return Err(invalid(
                "shell material spectrum requires at most 64 positive finite frequencies",
            ));
        }
        let tones = angular_frequencies_rad_s.len() as u128;
        let mut points = 0_u128;
        let mut bytes = self.sections.len() as u128
            * (core::mem::size_of::<PreparedFacet>()
                + core::mem::size_of::<Vec<LayerSpectrum<'_>>>()) as u128;
        let mut work = (32768 + 32768 * tones) * self.sections.len() as u128;
        for section in &self.sections {
            checkpoint(cx)?;
            for layer in 0..section.layer_count() {
                checkpoint(cx)?;
                let branches = section.layer_law(layer)?.branches().len() as u128;
                points += 6;
                bytes += core::mem::size_of::<LayerSpectrum<'_>>() as u128
                    + 128 * tones * (1 + branches)
                    + tones * core::mem::size_of::<PlaneStressHarmonic>() as u128;
                work += (6 * 262144 * tones * tones + 65536 * (tones + 1)) * (1 + branches);
            }
        }
        bytes += points
            * (core::mem::size_of::<TensorShellMaterialSpectrumPoint>() as u128
                + tones * core::mem::size_of::<[C64; 3]>() as u128)
            + 2 * tones * core::mem::size_of::<f64>() as u128
            + tones
                * (core::mem::size_of::<[[C64; 6]; 3]>() + core::mem::size_of::<[C64; 3]>())
                    as u128
            + core::mem::size_of::<PlaneStressSpectrumWindow<'_, '_>>() as u128
            + 256 * tones
            + 64 * tones * tones;
        cap(
            "shell material spectrum points",
            points,
            budget.maximum_points as u128,
        )?;
        cap(
            "shell material spectrum payload bytes",
            bytes,
            budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "shell material spectrum work",
            work,
            budget.maximum_work as u128,
        )?;
        let mut facets = Vec::with_capacity(self.sections.len());
        let mut layers = Vec::with_capacity(self.sections.len());
        for (element, section) in self.sections.iter().enumerate() {
            checkpoint(cx)?;
            facets.push(PreparedFacet::new(self, element)?);
            let mut spectra = Vec::with_capacity(section.layer_count());
            for layer in 0..section.layer_count() {
                checkpoint(cx)?;
                let (faces_m, to_principal) = section.layer_geometry(element, layer)?;
                let material = section
                    .layer_law(layer)?
                    .prepare_plane_stress_harmonic_spectrum(
                        MaterialPlane::Axes12,
                        angular_frequencies_rad_s,
                    )
                    .map_err(|e| section.material_error(element, layer, e))?;
                spectra.push(LayerSpectrum {
                    material,
                    faces_m,
                    to_principal,
                });
            }
            layers.push(spectra);
        }
        checkpoint(cx)?;
        Ok(TensorShellMaterialSpectrumObserver {
            source: self,
            angular_frequencies_rad_s: angular_frequencies_rad_s.to_vec(),
            facets,
            layers,
            points: points as usize,
        })
    }
}

impl<'a> TensorShellMaterialSpectrumObserver<'a> {
    pub const fn source(&self) -> &'a TensorPronyShell {
        self.source
    }
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        &self.angular_frequencies_rad_s
    }

    /// Observe one uniform physical interval on the common exp(-i omega t) clock.
    ///
    /// Each input contains all original free translations/rotations for its tone.
    /// At every original facet vertex and layer face, the recovered full 3D
    /// waveform must satisfy the source strain ceiling before point publication.
    /// The conservative distinct-frequency envelope bounds the complete affine
    /// CST/DKT field; it is neither a time-sampled nor an exact multitone maximum.
    ///
    /// Three triangle midpoint weights and two positive Gauss points per layer
    /// integrate this affine spatial field's quadratic material quantities exactly
    /// in real arithmetic. Frequency-specific transverse relaxation is retained.
    /// No kinetic/drilling/boundary energy, startup memory, thermal evolution,
    /// force-equation validation or spatial-convergence certificate is inferred.
    pub fn observe_window(
        &self,
        displacements: &[&[C64]],
        window: HarmonicWindow,
        cx: &Cx<'_>,
    ) -> Result<TensorShellMaterialSpectrumField, TensorShellHarmonicError> {
        self.observe(displacements, None, window, cx)
    }

    /// Observe the same coherent motion about an explicitly relaxed static field.
    ///
    /// `bias` contains one real displacement/rotation per original free coordinate.
    /// Each layer recovers its static transverse strain using the full equilibrium
    /// tensor at zero frequency, separately from every dispersive AC constraint.
    /// Equilibrium DC/AC cross terms enter physical storage and stress work;
    /// relaxed bias adds no Kelvin/Maxwell heat. The complete offset waveform must
    /// satisfy the original material strain ceiling at every vertex/layer face.
    /// Nonzero material bias requires original DC source coverage. This does not establish
    /// static force equilibrium, geometric prestress, thermal strain or startup memory.
    pub fn observe_window_with_bias(
        &self,
        displacements: &[&[C64]],
        bias: &[f64],
        window: HarmonicWindow,
        cx: &Cx<'_>,
    ) -> Result<TensorShellMaterialSpectrumField, TensorShellHarmonicError> {
        self.observe(displacements, Some(bias), window, cx)
    }

    fn observe(
        &self,
        displacements: &[&[C64]],
        bias: Option<&[f64]>,
        window: HarmonicWindow,
        cx: &Cx<'_>,
    ) -> Result<TensorShellMaterialSpectrumField, TensorShellHarmonicError> {
        checkpoint(cx)?;
        let tones = self.angular_frequencies_rad_s.len();
        if window.profile() != HarmonicWindowProfile::Uniform
            || displacements.len() != tones
            || bias.is_some_and(|q| {
                q.len() != self.source.model.free || q.iter().any(|v| !v.is_finite())
            })
            || displacements.iter().any(|q| {
                q.len() != self.source.model.free
                    || q.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            })
        {
            return Err(invalid(
                "physical shell spectrum observation needs a uniform interval and complete finite displacement vectors",
            ));
        }
        const BARYCENTRIC: [[f64; 3]; 3] = [[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]];
        const GAUSS: f64 = 0.5773502691896258;
        let mut points = Vec::with_capacity(self.points);
        let mut midpoint = Vec::with_capacity(tones);
        let mut in_plane = Vec::with_capacity(tones);
        let mut sums = [Dd::ZERO; 7];
        for (element, (facet, layers)) in self.facets.iter().zip(&self.layers).enumerate() {
            checkpoint(cx)?;
            let section = &self.source.sections[element];
            let bias_midpoint = bias.map_or([[C64::ZERO; 6]; 3], |q| {
                facet.strains_by(|i| C64::from_re(q[i]))
            });
            midpoint.clear();
            for q in displacements {
                checkpoint(cx)?;
                midpoint.push(facet.strains(q));
            }
            for (layer, spectrum) in layers.iter().enumerate() {
                checkpoint(cx)?;
                for vertex in 0..3 {
                    for z in spectrum.faces_m {
                        checkpoint(cx)?;
                        in_plane.clear();
                        for tone in &midpoint {
                            let strain = core::array::from_fn(|r| {
                                tone[vertex][r] + tone[(vertex + 2) % 3][r]
                                    - tone[(vertex + 1) % 3][r]
                            });
                            in_plane.push(principal_strain(&strain, z, &spectrum.to_principal));
                        }
                        let strain = core::array::from_fn(|r| {
                            bias_midpoint[vertex][r] + bias_midpoint[(vertex + 2) % 3][r]
                                - bias_midpoint[(vertex + 1) % 3][r]
                        });
                        let bias =
                            principal_strain(&strain, z, &spectrum.to_principal).map(|v| v.re);
                        spectrum
                            .material
                            .validate_strains_with_bias(&in_plane, &bias)
                            .map_err(|e| section.material_error(element, layer, e))?;
                    }
                }
                let [bottom, top] = spectrum.faces_m;
                let half = 0.5 * (top - bottom);
                let center = bottom + half;
                let volume = (facet.area_m2 / 3.0) * half;
                if !half.is_finite() || half <= 0.0 || !volume.is_finite() || volume <= 0.0 {
                    return Err(invalid(
                        "original shell layer quadrature volume is unrepresentable",
                    ));
                }
                let temporal = spectrum
                    .material
                    .prepare_window(window)
                    .map_err(|e| section.material_error(element, layer, e))?;
                for triangle_quadrature_point in 0..3 {
                    for (thickness_quadrature_point, sign) in [-1.0, 1.0].into_iter().enumerate() {
                        checkpoint(cx)?;
                        let z = center + sign * GAUSS * half;
                        in_plane.clear();
                        for tone in &midpoint {
                            in_plane.push(principal_strain(
                                &tone[triangle_quadrature_point],
                                z,
                                &spectrum.to_principal,
                            ));
                        }
                        // Preserve original element/layer diagnostics while the
                        // native material checks cancellation between branches.
                        let bias = principal_strain(
                            &bias_midpoint[triangle_quadrature_point],
                            z,
                            &spectrum.to_principal,
                        )
                        .map(|v| v.re);
                        let response = temporal
                            .evaluate_with_bias_and_checkpoint::<EvaluationError>(
                                &in_plane,
                                &bias,
                                || cx.checkpoint().map_err(|_| EvaluationError::Cancelled),
                            )
                            .map_err(|e| match e {
                                EvaluationError::Material(error) => {
                                    section.material_error(element, layer, error)
                                }
                                EvaluationError::Cancelled => {
                                    super::super::PlateError::Cancelled.into()
                                }
                            })?;
                        for (sum, value) in sums.iter_mut().zip([
                            response.start_stored_energy_j_m3,
                            response.end_stored_energy_j_m3,
                            response.stored_energy_change_j_m3,
                            response.kelvin_heat_j_m3,
                            response.maxwell_heat_j_m3,
                            response.material_heat_j_m3,
                            response.stress_work_j_m3,
                        ]) {
                            *sum = *sum + Dd::from_f64(volume) * Dd::from_f64(value);
                        }
                        points.push(TensorShellMaterialSpectrumPoint {
                            element,
                            layer,
                            triangle_quadrature_point,
                            thickness_quadrature_point,
                            triangle_barycentric: BARYCENTRIC[triangle_quadrature_point],
                            thickness_position_m: z,
                            volume_m3: volume,
                            in_plane_strain: in_plane.clone(),
                            in_plane_bias: bias,
                            material: response,
                        });
                    }
                }
            }
        }
        let balance = (sums[6] - sums[2] - sums[5]).to_f64();
        let [start, end, change, kelvin, maxwell, heat, work] = sums.map(Dd::to_f64);
        if [start, end, change, kelvin, maxwell, heat, work, balance]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "original shell material spectrum totals are not representable",
            ));
        }
        checkpoint(cx)?;
        Ok(TensorShellMaterialSpectrumField {
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
