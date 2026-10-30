//! Source-preserving fluid admission to the existing time-domain DG owner.
//!
//! The same DG masses, pressure/flow ports, wall coordinates and exact
//! exponential operate on a fluid's actual rho and adiabatic sound speed.
//! Generic fluid admission creates no gas EOS or mutable second material.

use super::*;
use crate::fluid::FluidDuctError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};

mod update;

mod sealed {
    pub trait Medium {}
}

/// Read-only coefficients of an original admitted DG medium.
///
/// Implementations are the existing gas state and the source-preserving fluid
/// wrapper produced by [`DgDuct::from_material`]. This bound keeps arbitrary
/// mutable coefficient providers outside the assembled physical owner.
pub trait DgMedium: sealed::Medium + Clone {
    /// Equilibrium density [kg/m³].
    fn density_kg_m3(&self) -> f64;
    /// Adiabatic small-signal sound speed [m/s].
    fn sound_speed_m_s(&self) -> f64;
    /// Original admission, before constructing any physical operator.
    #[doc(hidden)]
    fn admit(&self) -> Result<(), DuctError>;
}

impl sealed::Medium for GasState {}
impl DgMedium for GasState {
    fn density_kg_m3(&self) -> f64 {
        self.density
    }
    fn sound_speed_m_s(&self) -> f64 {
        self.sound_speed
    }
    fn admit(&self) -> Result<(), DuctError> {
        crate::validate_gas(self)
    }
}

/// Immutable original fluid owner and its admitted local acoustic response.
/// Source receipts, EOS/phase policy and query coordinates stay on `source`.
#[derive(Clone, Debug)]
pub struct DgFluidMaterial<M> {
    source: M,
    acoustic: FluidAcousticState,
    // Frozen caloric/transport tuple that produced the actual diffusion arms.
    // A live refresh must not reevaluate a changed source to infer old storage.
    thermoviscous: Option<ThermoviscousFluidState>,
}

impl<M> DgFluidMaterial<M> {
    /// Original complete material evaluation, not a synthesized gas surrogate.
    pub const fn source(&self) -> &M {
        &self.source
    }
    /// Coefficients actually used in the acoustic mass and characteristic maps.
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
}
impl<M: Clone> sealed::Medium for DgFluidMaterial<M> {}
impl<M: Clone> DgMedium for DgFluidMaterial<M> {
    fn density_kg_m3(&self) -> f64 {
        self.acoustic.density_kg_m3()
    }
    fn sound_speed_m_s(&self) -> f64 {
        self.acoustic.sound_speed_m_s()
    }
    fn admit(&self) -> Result<(), DuctError> {
        Ok(())
    }
}

impl sealed::Medium for FluidAcousticState {}
impl DgMedium for FluidAcousticState {
    fn density_kg_m3(&self) -> f64 {
        (*self).density_kg_m3()
    }
    fn sound_speed_m_s(&self) -> f64 {
        (*self).sound_speed_m_s()
    }
    fn admit(&self) -> Result<(), DuctError> {
        Ok(())
    }
}

impl DgDuct {
    /// Bind one original uniform fluid evaluation to the native DG operator.
    ///
    /// Cylinders, cones and profiles use the SAME geometric quadrature, mass
    /// transforms, boundaries and time propagator as gas ducts. All original
    /// source properties and receipts are retained. A later source refresh
    /// must be explicitly rebound; changing this frozen operator is not a heat,
    /// species or mean-flow evolution. The live `rebind_uniform_material` and
    /// `rebind_uniform_viscothermal_material` methods preserve physical state
    /// and account the work of an explicitly prescribed background change.
    ///
    /// The resulting owner supports ordinary native state initialization,
    /// propagation, exact pressure/flow work, side ports, continuous harmonic
    /// response, compliant walls and optional thermoviscous memory. Gas-only
    /// molecular and gas-field update methods remain on the gas specialization.
    pub fn from_material<M: FluidAcoustics + Clone>(
        duct: &Duct,
        material: &M,
        mesh: DgMesh,
        inlet: DgBoundary,
        outlet: DgBoundary,
    ) -> Result<DgDuct<DgFluidMaterial<M>>, FluidDuctError> {
        let acoustic = material.acoustic_state()?;
        let medium = DgFluidMaterial {
            source: material.clone(),
            acoustic,
            thermoviscous: None,
        };
        Ok(DgDuct::assemble(
            duct,
            &medium,
            &medium,
            mesh,
            inlet,
            outlet,
            None,
            &mut |_| Ok::<_, DuctError>(acoustic),
        )?)
    }
}

impl<M: Clone> DgDuct<DgFluidMaterial<M>> {
    /// Complete original source retained by this uniform material duct.
    pub const fn material(&self) -> &M {
        self.gas.source()
    }
    /// Actual frozen acoustic coefficients used by the native operator.
    pub const fn acoustic_material(&self) -> FluidAcousticState {
        self.gas.acoustic()
    }
    /// Exact frozen tuple used by the current diffusion arms, if attached.
    pub const fn thermoviscous_material(&self) -> Option<ThermoviscousFluidState> {
        self.gas.thermoviscous
    }
}

impl<M: ThermoviscousFluid + Clone> DgDuct<DgFluidMaterial<M>> {
    /// Attach the existing positive Bessel-derived wall diffusion realization
    /// using this SAME original fluid's viscosity, conductivity, Cp and Cp/Cv.
    ///
    /// Uniform circular geometry is required. Every pole remains in the native
    /// dense exponential and energy state; no separate time integrator, fitted
    /// gas identity or duplicated inviscid storage enters the system. The
    /// original sampled band/error admission and complete-state DOF cap apply.
    /// Cp == Cv creates no thermal arms, while retaining viscous diffusion.
    pub fn with_viscothermal(
        mut self,
        poles: usize,
        band_hz: [f64; 2],
        maximum_relative_error: f64,
    ) -> Result<Self, FluidDuctError> {
        let state = self.material().thermoviscous_state()?;
        if state.acoustic() != self.gas.acoustic {
            return Err(DuctError::BadParameter {
                what: "DG thermoviscous and acoustic views must be the same original material state",
            }.into());
        }
        let radius = self.uniform_radius_m.ok_or(DuctError::BadParameter {
            what: "DG thermoviscous fluid memory requires identical circular cylinders",
        })?;
        let transport = state.transport();
        let model = fs_phs::ZkDiffusive::new(
            radius,
            state.acoustic().density_kg_m3(),
            transport.dynamic_viscosity_pa_s(),
            state.heat_capacity_ratio(),
            transport.prandtl(),
            poles,
            band_hz.map(|f| core::f64::consts::TAU * f),
            maximum_relative_error,
        )
        .map_err(|_| DuctError::BadParameter {
            what: "DG fluid viscothermal Bessel pole data or declared band is unresolved",
        })?;
        self.gas.thermoviscous = Some(state);
        self.outlet_gas.thermoviscous = Some(state);
        Ok(self.with_uniform_viscothermal_model(model)?)
    }
}
