//! Shared gas composition and refractivity into the original spectral renderer.
//!
//! The renderer consumes the material owner's prepared mass coefficient and
//! depolarization; no renderer-side gas law, species preset or spectral table
//! approximation is introduced. Scattering and absorption remain distinct:
//! this elastic molecular model emits no heat and supplies no absorption.

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_material::gas_scattering::{GasScatteringError, RAYLEIGH_REFERENCE_NM, RayleighGasState};

use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};

use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
};
use crate::volumes::transport::{
    AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum, VolumeTransportError,
};

/// Original gas/optical state and its exact visible-band render coefficients.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialGasScatteringBinding {
    state: RayleighGasState,
    scattering: MassScatteringSpectrum,
    phase: PhaseFunction,
}

impl MaterialGasScatteringBinding {
    /// Prepare exact inverse-fourth spectral scattering and the analytic
    /// depolarized phase from the admitted material model.
    pub fn try_new(state: RayleighGasState) -> Result<Self, VolumeTransportError> {
        let scattering = MassScatteringSpectrum::rayleigh(
            RAYLEIGH_REFERENCE_NM,
            state.mass_scattering_reference_m2_kg(),
        )?;
        let phase = PhaseFunction::rayleigh_depolarized(state.depolarization_ratio())?;
        Ok(Self {
            state,
            scattering,
            phase,
        })
    }

    /// The exact original state, source references and gas composition identity.
    #[must_use]
    pub const fn state(&self) -> &RayleighGasState {
        &self.state
    }

    /// Coefficients for a caller-owned density field of this frozen composition.
    /// Density is kg/m³; no dimensionless LBM scaling or species inference occurs.
    /// A different composition requires a fresh material evaluation and binding.
    #[must_use]
    pub const fn scattering_spectrum(&self) -> MassScatteringSpectrum {
        self.scattering
    }

    /// Exact scalar molecular phase with the original depolarization convention.
    #[must_use]
    pub const fn phase(&self) -> &PhaseFunction {
        &self.phase
    }

    /// Bind a smooth density reconstruction of this SAME frozen gas composition
    /// to native curved-ray propagation. `density_field` supplies density in
    /// kg/m³ and its actual spatial gradient in kg/m⁴; a cell-average grid is
    /// not automatically a differentiable field. The caller identifies the
    /// accepted snapshot and reconstruction through `snapshot_identity`.
    ///
    /// Molecular refractivity comes from the original optical references,
    /// not a fitted renderer coefficient. This adapter supplies phase/group
    /// index and zero true absorption. It does not turn elastic scattering
    /// into absorption or integrate scattered radiance along the bent path.
    pub fn bind_refractive_field<'a, F: GasDensityGradientField + ?Sized>(
        &self,
        density_field: &'a F,
        snapshot_identity: ContentHash,
    ) -> Result<MaterialGasRefractiveField<'a, F>, GasRefractiveFieldError> {
        if snapshot_identity == ContentHash([0; 32]) {
            return Err(GasRefractiveFieldError::MissingSnapshotIdentity);
        }
        let mut hash = DomainHasher::new("org.frankensim.render.material-gas-refractive-field.v1");
        hash.update(self.state.content_identity().as_bytes());
        hash.update(snapshot_identity.as_bytes());
        Ok(MaterialGasRefractiveField {
            state: self.state,
            density_field,
            identity: hash.finalize(),
        })
    }

    /// Bind an actual homogeneous gas snapshot, borrowing its original density
    /// buffer. Every cell must equal the accepted material density; a later
    /// temperature/pressure state must be reevaluated before rebinding.
    ///
    /// This is explicitly the elastic, nonabsorbing molecular component. It
    /// does not claim spectral gas absorption, thermal emission, refraction,
    /// polarization or a spatially varying composition. For spatially varying
    /// mass density at fixed composition, use the exposed spectrum and phase
    /// with `ScatteringGrid`'s independently admitted physical field.
    pub fn bind_homogeneous_grid<'a>(
        &self,
        cx: &Cx<'_>,
        grid: VolumeGrid<'a>,
        budget: AbsorbingGridBudget,
    ) -> Result<MaterialGasScatteringGrid<'a>, GasGridBindingError> {
        cx.checkpoint()
            .map_err(|_| GasGridBindingError::Transport(VolumeTransportError::Cancelled))?;
        if grid.data.len() > budget.max_grid_cells {
            return Err(GasGridBindingError::Transport(VolumeTransportError::Budget));
        }
        for (chunk_index, chunk) in grid.data.chunks(1024).enumerate() {
            cx.checkpoint()
                .map_err(|_| GasGridBindingError::Transport(VolumeTransportError::Cancelled))?;
            if let Some(index) = chunk
                .iter()
                .position(|density| *density != self.state.gas().density)
            {
                return Err(GasGridBindingError::DensityMismatch {
                    index: chunk_index * 1024 + index,
                });
            }
        }
        let absorption = AbsorbingGrid::try_new(
            cx,
            grid,
            MassAbsorptionSpectrum::constant(0.0).map_err(GasGridBindingError::Transport)?,
            None,
            budget,
        )
        .map_err(GasGridBindingError::Transport)?;
        let medium = ScatteringGrid::try_new(cx, absorption, self.scattering, self.phase.clone())
            .map_err(GasGridBindingError::Transport)?;
        let mut hash = DomainHasher::new("org.frankensim.render.material-gas-scattering-grid.v1");
        hash.update(self.state.content_identity().as_bytes());
        hash.update(medium.content_identity().as_bytes());
        Ok(MaterialGasScatteringGrid {
            state: self.state,
            medium,
            identity: hash.finalize(),
        })
    }
}

/// One value and gradient from a caller-owned smooth physical gas snapshot.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasDensityGradientSample {
    /// Actual nonnegative mass density [kg/m³].
    pub density_kg_m3: f64,
    /// Spatial derivative of that same reconstruction [kg/m⁴].
    pub gradient_kg_m4: Vec3,
}

/// Explicit reconstruction seam; no finite differences, smoothing, LBM unit
/// conversion or outside-domain vacuum are inferred by the renderer.
pub trait GasDensityGradientField {
    /// Evaluate the immutable snapshot, refusing points outside its support.
    fn sample_density(
        &self,
        position_m: Point3,
    ) -> Result<GasDensityGradientSample, GradientIndexError>;
}

impl<F> GasDensityGradientField for F
where
    F: Fn(Point3) -> Result<GasDensityGradientSample, GradientIndexError>,
{
    fn sample_density(
        &self,
        position_m: Point3,
    ) -> Result<GasDensityGradientSample, GradientIndexError> {
        self(position_m)
    }
}

/// Native GRIN input tied to the original gas sources and one physical snapshot.
/// The borrowed provider is evaluated without allocations at each ray stage.
pub struct MaterialGasRefractiveField<'a, F: GasDensityGradientField + ?Sized> {
    state: RayleighGasState,
    density_field: &'a F,
    identity: ContentHash,
}

impl<F: GasDensityGradientField + ?Sized> MaterialGasRefractiveField<'_, F> {
    /// Original accepted gas and the exact source-bound molecular optics.
    #[must_use]
    pub const fn state(&self) -> &RayleighGasState {
        &self.state
    }

    /// Original gas/source identity plus caller-declared snapshot/reconstruction.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Evaluate with the original material refusal available to diagnostics.
    /// The `GradientIndexField` entry point maps a material-domain refusal to
    /// `OutsideFieldDomain` so the native integrator never extrapolates it.
    pub fn sample_optics(
        &self,
        position_m: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GasRefractiveFieldError> {
        if ![position_m.x, position_m.y, position_m.z]
            .iter()
            .all(|x| x.is_finite())
        {
            return Err(GasRefractiveFieldError::Field(
                GradientIndexError::InvalidSample,
            ));
        }
        let density = self
            .density_field
            .sample_density(position_m)
            .map_err(GasRefractiveFieldError::Field)?;
        if ![
            density.gradient_kg_m4.x,
            density.gradient_kg_m4.y,
            density.gradient_kg_m4.z,
        ]
        .iter()
        .all(|x| x.is_finite())
        {
            return Err(GasRefractiveFieldError::Field(
                GradientIndexError::InvalidSample,
            ));
        }
        let index = self
            .state
            .refractive_index_at_density(density.density_kg_m3, wavelength_nm)
            .map_err(GasRefractiveFieldError::Gas)?;
        GradientIndexSample::try_new(
            index.phase_index,
            density.gradient_kg_m4.scale(index.derivative_density_m3_kg),
            0.0,
            Some(index.phase_index),
        )
        .map_err(GasRefractiveFieldError::Field)
    }
}

impl<F: GasDensityGradientField + ?Sized> GradientIndexField for MaterialGasRefractiveField<'_, F> {
    fn sample(
        &self,
        position: Point3,
        wavelength_nm: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        self.sample_optics(position, wavelength_nm)
            .map_err(|error| match error {
                GasRefractiveFieldError::Field(error) => error,
                GasRefractiveFieldError::Gas(GasScatteringError::RefractivityDomain) => {
                    GradientIndexError::OutsideFieldDomain
                }
                GasRefractiveFieldError::Gas(_)
                | GasRefractiveFieldError::MissingSnapshotIdentity => {
                    GradientIndexError::InvalidSample
                }
            })
    }
}

/// A missing source/snapshot or an original field/material refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum GasRefractiveFieldError {
    /// The caller did not identify the physical density reconstruction.
    MissingSnapshotIdentity,
    /// Original molecular optics refusal, including the dilute-use boundary.
    Gas(GasScatteringError),
    /// Original provider or native GRIN admission refusal.
    Field(GradientIndexError),
}

impl core::fmt::Display for GasRefractiveFieldError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material gas refractive field refused: {self:?}")
    }
}
impl core::error::Error for GasRefractiveFieldError {}

/// A borrowed rendering grid retaining its original gas and optical source state.
pub struct MaterialGasScatteringGrid<'a> {
    state: RayleighGasState,
    medium: ScatteringGrid<'a>,
    identity: ContentHash,
}

impl<'a> MaterialGasScatteringGrid<'a> {
    /// Direct input to `render_scattering` and physical-grid photometry.
    #[must_use]
    pub const fn medium(&self) -> &ScatteringGrid<'a> {
        &self.medium
    }
    /// Original material state; retained independently of the numerical field hash.
    #[must_use]
    pub const fn state(&self) -> &RayleighGasState {
        &self.state
    }
    /// Complete source/state binding and numerical rendering snapshot identity.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }
}

/// Atomic refusal before a mismatched gas density enters light transport.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GasGridBindingError {
    /// The borrowed cell does not carry the accepted homogeneous density.
    DensityMismatch {
        /// Index in the original borrowed simulation buffer.
        index: usize,
    },
    /// Original grid geometry, optical transport, budget or cancellation refusal.
    Transport(VolumeTransportError),
}

impl core::fmt::Display for GasGridBindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material gas scattering grid refused: {self:?}")
    }
}
impl core::error::Error for GasGridBindingError {}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_geom::{Point3, Vec3};
    use fs_material::{
        gas::GasSpec,
        gas_scattering::{MolecularRayleigh, MolecularRayleighMixture},
    };
    use fs_qty::semantic::{Composition, CompositionBasis};

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            f(
                &Cx::new(
                    &gate,
                    arena,
                    StreamKey {
                        seed: 9,
                        kernel_id: 1,
                        tile: 0,
                        iteration: 0,
                    },
                    Budget::INFINITE,
                    ExecMode::Deterministic,
                ),
                &gate,
            )
        })
    }

    fn binding(pressure: f64, source: u8) -> MaterialGasScatteringBinding {
        let optical =
            MolecularRayleigh::try_new(0.0003, 2.5e25, 550.0, 0.03, ContentHash([source; 32]))
                .unwrap();
        let model = MolecularRayleighMixture::try_new(
            &[GasSpec::dry_air_ussa1976()],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
            &[optical],
            None,
        )
        .unwrap();
        MaterialGasScatteringBinding::try_new(model.evaluate(300.0, pressure).unwrap()).unwrap()
    }

    fn budget() -> AbsorbingGridBudget {
        AbsorbingGridBudget {
            max_grid_cells: 2,
            max_steps_per_segment: 8,
        }
    }

    #[test]
    fn g1_material_gas_reaches_exact_native_spectral_transport() {
        with_cx(|cx, _| {
            for pressure in [50_000.0, 100_000.0, 200_000.0] {
                let binding = binding(pressure, 1);
                let density = [binding.state().gas().density; 2];
                let grid = binding
                    .bind_homogeneous_grid(
                        cx,
                        VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0, 1.0, 1000.0]),
                        budget(),
                    )
                    .unwrap();
                let wavelengths = [380.0, 417.3, 550.0, 709.7, 780.0];
                let result = grid
                    .medium()
                    .transmittance(
                        cx,
                        Point3::new(0.5, 0.5, 0.0),
                        Vec3::new(0.0, 0.0, 1.0),
                        2000.0,
                        &wavelengths,
                    )
                    .unwrap();
                for (got, lambda) in result.into_iter().zip(wavelengths) {
                    let expected =
                        (-2000.0 * binding.state().scattering_per_m(lambda).unwrap()).exp();
                    assert!(
                        (got - expected).abs() < 2e-13,
                        "{lambda}: {got} != {expected}"
                    );
                }
                for cosine in [-1.0, -0.6, 0.0, 0.4, 1.0] {
                    let phase = binding.phase().eval_cos(cosine).unwrap();
                    assert!((phase - binding.state().phase_per_sr(cosine).unwrap()).abs() < 2e-15);
                }
                assert!(!grid.medium().has_thermal_emission());
                assert_eq!(grid.state(), binding.state());
            }
        });
    }

    #[test]
    fn g0_gas_grid_refuses_density_mismatch_and_retains_optical_provenance() {
        with_cx(|cx, _| {
            let a = binding(100_000.0, 1);
            let b = binding(100_000.0, 2);
            let density = [a.state().gas().density];
            let make = || VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0; 3]);
            let ga = a.bind_homogeneous_grid(cx, make(), budget()).unwrap();
            let gb = b.bind_homogeneous_grid(cx, make(), budget()).unwrap();
            assert_eq!(
                ga.medium().content_identity(),
                gb.medium().content_identity()
            );
            assert_ne!(ga.content_identity(), gb.content_identity());
            let bad = [a.state().gas().density * 2.0];
            assert!(matches!(
                a.bind_homogeneous_grid(
                    cx,
                    VolumeGrid::new([1; 3], &bad, [0.0; 3], [1.0; 3]),
                    budget()
                ),
                Err(GasGridBindingError::DensityMismatch { index: 0 })
            ));
        });
    }

    #[test]
    fn g1_original_gas_sources_supply_refractive_index_and_density_gradient() {
        let other_source = binding(100_000.0, 2);
        let binding = binding(100_000.0, 1);
        let rho = binding.state().gas().density;
        let reconstruction = |p: Point3| {
            Ok(GasDensityGradientSample {
                density_kg_m3: rho * (1.0 + 0.1 * p.y),
                gradient_kg_m4: Vec3::new(0.0, 0.1 * rho, 0.0),
            })
        };
        let field = binding
            .bind_refractive_field(&reconstruction, ContentHash([4; 32]))
            .unwrap();
        let p = Point3::new(0.5, 0.0, 0.5);
        for wavelength in [380.0, 550.0, 780.0] {
            let sample = field.sample(p, wavelength).unwrap();
            let original = binding.state().refractive_index(wavelength).unwrap();
            assert_eq!(sample.phase_index(), original.phase_index);
            assert_eq!(sample.group_index(), Some(original.phase_index));
            assert_eq!(sample.absorption_per_m(), 0.0);
            assert_eq!(
                sample.phase_gradient_per_m(),
                Vec3::new(0.0, 0.1 * rho * original.derivative_density_m3_kg, 0.0)
            );
            let h = 0.01;
            let a = field.sample(Point3::new(p.x, -h, p.z), wavelength).unwrap();
            let b = field.sample(Point3::new(p.x, h, p.z), wavelength).unwrap();
            let finite_difference = (b.phase_index() - a.phase_index()) / (2.0 * h);
            assert!((finite_difference / sample.phase_gradient_per_m().y - 1.0).abs() < 1.0e-8);
        }
        assert_eq!(field.state(), binding.state());
        let updated = binding
            .bind_refractive_field(&reconstruction, ContentHash([5; 32]))
            .unwrap();
        assert_ne!(field.content_identity(), updated.content_identity());
        let changed = other_source
            .bind_refractive_field(&reconstruction, ContentHash([4; 32]))
            .unwrap();
        assert_ne!(field.content_identity(), changed.content_identity());
    }

    #[test]
    fn g0_gas_refractive_field_preserves_domain_and_invalid_gradient_refusals() {
        let binding = binding(100_000.0, 1);
        let outside = |_p: Point3| Err(GradientIndexError::OutsideFieldDomain);
        let field = binding
            .bind_refractive_field(&outside, ContentHash([4; 32]))
            .unwrap();
        assert_eq!(
            field.sample(Point3::new(0.0, 0.0, 0.0), 550.0),
            Err(GradientIndexError::OutsideFieldDomain)
        );
        assert!(matches!(
            binding.bind_refractive_field(&outside, ContentHash([0; 32])),
            Err(GasRefractiveFieldError::MissingSnapshotIdentity)
        ));
        for (density, gradient, expected) in [
            (
                -1.0,
                Vec3::new(0.0, 0.0, 0.0),
                GradientIndexError::OutsideFieldDomain,
            ),
            (
                1.0,
                Vec3::new(f64::NAN, 0.0, 0.0),
                GradientIndexError::InvalidSample,
            ),
            (
                1.0e20,
                Vec3::new(0.0, 0.0, 0.0),
                GradientIndexError::OutsideFieldDomain,
            ),
        ] {
            let invalid = |_p: Point3| {
                Ok(GasDensityGradientSample {
                    density_kg_m3: density,
                    gradient_kg_m4: gradient,
                })
            };
            let field = binding
                .bind_refractive_field(&invalid, ContentHash([4; 32]))
                .unwrap();
            assert_eq!(
                field.sample(Point3::new(0.0, 0.0, 0.0), 550.0),
                Err(expected)
            );
        }
    }

    #[test]
    fn g1_gas_density_gradient_bends_native_ray_and_retains_stratified_invariant() {
        use crate::charts::Ray;
        use crate::gradient_index::{GradientIndexOptions, propagate};

        let binding = binding(100_000.0, 1);
        let rho = binding.state().gas().density;
        let reconstruction = |p: Point3| {
            Ok(GasDensityGradientSample {
                density_kg_m3: rho * (1.0 + 0.2 * p.y),
                gradient_kg_m4: Vec3::new(0.0, 0.2 * rho, 0.0),
            })
        };
        let field = binding
            .bind_refractive_field(&reconstruction, ContentHash([4; 32]))
            .unwrap();
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
        };
        let options = GradientIndexOptions {
            initial_step_m: 0.05,
            max_step_m: 0.1,
            ..Default::default()
        };
        let result = propagate(&field, ray, Vec3::new(0.0, 0.0, 1.0), 550.0, 2.0, options).unwrap();
        assert!(result.ray.origin.y > 1.0e-5 && result.ray.dir.y > 1.0e-5);
        assert_eq!(result.ray.origin.z, 0.0);
        // Translational symmetry in x conserves the x optical momentum.
        let invariant = result.final_sample.phase_index() * result.ray.dir.x;
        assert!((invariant - result.initial_sample.phase_index()).abs() < 1.0e-10);
        assert_eq!(
            result.group_optical_path_length_m,
            Some(result.optical_path_length_m)
        );
        assert_eq!(result.power_transmittance(), 1.0);
        assert!(result.optical_path_length_m > 2.0);
    }
}
