//! Shared gas composition and refractivity into the original spectral renderer.
//!
//! The renderer consumes the material owner's prepared mass coefficient and
//! depolarization; no renderer-side gas law, species preset or spectral table
//! approximation is introduced. Scattering and absorption remain distinct:
//! this elastic molecular model emits no heat and supplies no absorption.

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::gas_scattering::{RAYLEIGH_REFERENCE_NM, RayleighGasState};

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
}
