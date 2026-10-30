//! Source tables in temperature and physical shear rate drive the actual VOF
//! collision and molecular heat. The fixed-domain rheology owner supplies the
//! same local constitutive solve; no second viscosity law or stress estimator
//! is introduced. This is generalized Newtonian response, without elastic
//! memory, yield-surface history, phase mobility or capillary energy. A sourced
//! constant sigma can accompany the same complete thermal/rheology bundle; the
//! original native capillary operator retains its own mechanical accuracy.

use super::*;
use crate::material::free_surface::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceInitialCell3, MaterialFreeSurfaceScaling,
};
use crate::material::transport::rheology::{LocalRheologyError, solve_local};
use fs_exec::Cx;
use fs_material::liquid::surface::ResolvedLiquidSurfaceState;

pub use crate::material::transport::rheology::RheologySolve;

/// Accepted constitutive query at one actual wet collision cell. The original
/// material tuple and its source receipts are available through `material_at`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceRheologySample {
    /// Source query `sqrt(2 D:D)` [1/s], not a lattice or cyclic frequency.
    pub shear_rate_per_s: f64,
    /// `1/2 + 3 (mu/rho0) dt/dx^2` from this SAME source state.
    pub collision_tau: f64,
    pub iterations: usize,
    /// Absolute constitutive fixed-point defect in dimensionless tau.
    pub relaxation_defect: f64,
}

#[derive(Clone, Debug)]
pub(in crate::material::free_surface) struct RheologyBinding {
    rate_axis: String,
    solve: RheologySolve,
}

#[derive(Clone, Copy)]
pub(super) struct RheologyStress {
    /// Force-corrected lattice `gamma_dot * tau`, before the thermal split.
    numerator: f64,
    initial_tau: f64,
}

#[derive(Clone)]
pub(super) struct RheologyModel<'m> {
    binding: RheologyBinding,
    law: &'m LiquidTemperatureLaw<'m>,
    caloric: &'m EquilibriumEnthalpyPhaseCurve,
    density: f64,
    buoyancy: Option<BoussinesqLiquid>,
    dx: f64,
    dt: f64,
    max_speed: f64,
}

pub(super) struct RheologicalState {
    pub material: ResolvedLiquidState,
    pub sample: FreeSurfaceRheologySample,
}

impl RheologyBinding {
    pub(in crate::material::free_surface) fn try_new(
        source: &ResolvedLiquidSurfaceState,
        law: &LiquidTemperatureLaw<'_>,
        solve: RheologySolve,
    ) -> Result<Self, Error> {
        solve
            .validate()
            .map_err(|_| Error::Invalid("invalid liquid rheology controls"))?;
        let (temperature_axis, rate_axis) = source.rheology_axes().ok_or(Error::Invalid(
            "moving-liquid rheology needs an explicitly evolving source bundle",
        ))?;
        if temperature_axis != law.temperature_axis() {
            return Err(Error::Invalid(
                "rheology and caloric temperature axes must agree",
            ));
        }
        // Thermal conductivity and Cp may belong to a broader law bundle than
        // the momentum ingress. They cannot introduce an unevolved source axis.
        ResolvedLiquidSurfaceState::try_for_shear_rate_evolution_with_capillarity(
            law.reference_state().material().clone(),
            law.viscosity_basis(),
            source.capillarity(),
            temperature_axis,
            rate_axis,
        )?;
        if ["density", "specific_heat_capacity"].iter().any(|name| {
            law.reference_state()
                .material()
                .property(name)
                .is_some_and(|property| {
                    property
                        .source_axes()
                        .iter()
                        .any(|(axis, _)| axis == rate_axis)
                })
        }) {
            return Err(Error::Invalid(
                "shear-dependent density or caloric storage needs an additional thermodynamic state owner",
            ));
        }
        Ok(Self {
            rate_axis: rate_axis.to_owned(),
            solve,
        })
    }
}

impl RheologyModel<'_> {
    pub(super) fn sample(
        &self,
        h: f64,
        stress: RheologyStress,
        check: impl FnMut() -> Result<(), Error>,
    ) -> Result<RheologicalState, Error> {
        let temperature = buoyancy::liquid_temperature(self.caloric, h, self.density)?;
        let root = solve_local(
            stress.numerator,
            self.dt,
            stress.initial_tau,
            self.binding.solve,
            check,
            |rate| {
                let point = self
                    .law
                    .reference_point()
                    .clone()
                    .with(self.law.temperature_axis(), temperature)
                    .and_then(|point| point.with(&self.binding.rate_axis, rate))
                    .map_err(|_| Error::Invalid("invalid liquid temperature/shear coordinates"))?;
                let material = self
                    .law
                    .reference_state()
                    .at_query_point(self.law.card(), &point)?;
                let transport = buoyancy::liquid_transport(&material, self.density, self.buoyancy)?;
                let (scaling, _) = plan_momentum_scaling(
                    transport.kinematic_viscosity_m2_s(),
                    self.dx,
                    self.dt,
                    self.max_speed,
                )?;
                Ok(((material, scaling.tau), scaling.tau))
            },
        )
        .map_err(|error| match error {
            LocalRheologyError::Evaluation(error) => error,
            LocalRheologyError::Invalid(what) => Error::Invalid(what),
            LocalRheologyError::Convergence {
                iterations,
                relaxation_defect,
            } => Error::RheologyConvergence {
                iterations,
                relaxation_defect,
            },
        })?;
        Ok(RheologicalState {
            material: root.value.0,
            sample: FreeSurfaceRheologySample {
                shear_rate_per_s: root.shear_rate_per_s,
                collision_tau: root.value.1,
                iterations: root.iterations,
                relaxation_defect: root.relaxation_defect,
            },
        })
    }
}

impl<'s> MaterialFreeSurfaceScaling<'s> {
    /// Initialize moving liquid from one source-backed T/shear constitutive law.
    ///
    /// Every occupied/ghost cell immediately solves physical shear against its
    /// own molecular viscosity. The SAME coefficient drives native collision
    /// and optional molecular heat. Later thermal substeps refresh all transport
    /// properties against the held initial stress, then VOF evolves once and
    /// admits the actual candidate stress before publication. Source-domain,
    /// convergence and cancellation failures return no initialized candidate.
    ///
    /// The ordinary frozen constructors refuse this source. Caloric density is
    /// fixed and all cells are pure liquid. Capillarity is neglected or uses
    /// source-declared constant sigma in the same complete thermal bundle.
    /// The numerical conversion carry transports enthalpy but has no assigned
    /// strain, viscosity or molecular heat until it returns to a resolved cell.
    #[allow(clippy::too_many_arguments)]
    pub fn free_surface3_from_cells_with_rheology<'m>(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        controls: FreeSurfaceEnthalpyControls,
        solve: RheologySolve,
        initial_h: impl Fn([usize; 3]) -> f64,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceEnthalpy3<'s, 'm>, Error> {
        self.free_surface3_from_cells_with_rheology_model(
            domain, initial, law, caloric, controls, solve, None, initial_h, cx, gate,
        )
    }

    /// Initialize one moving-liquid owner with original buoyancy AND rheology.
    ///
    /// The same pinned T/rate bundle supplies raw rho, Cp, mu and k. The existing
    /// Boussinesq owner bounds rho/Cp against the original reference while rho0
    /// owns transported mass, caloric storage and collision inertia. The native
    /// constitutive solve uses force-corrected stress and the original viscosity;
    /// attaching gravity does not freeze the rate or introduce another law.
    ///
    /// Initial force rebinding preserves physical momentum. Later heat, VOF
    /// transport, conduction and molecular heating all continue through the
    /// existing combined owner. Caloric states must remain pure liquid; no
    /// phase-volume change or solid-skeleton dynamics is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn free_surface3_from_cells_with_boussinesq_rheology<'m>(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        controls: FreeSurfaceEnthalpyControls,
        solve: RheologySolve,
        buoyancy: crate::material::transport::boussinesq::BoussinesqConfig<3>,
        initial_h: impl Fn([usize; 3]) -> f64,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceEnthalpy3<'s, 'm>, Error> {
        self.free_surface3_from_cells_with_rheology_model(
            domain, initial, law, caloric, controls, solve, Some(buoyancy), initial_h, cx, gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn free_surface3_from_cells_with_rheology_model<'m>(
        self,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        controls: FreeSurfaceEnthalpyControls,
        solve: RheologySolve,
        buoyancy: Option<crate::material::transport::boussinesq::BoussinesqConfig<3>>,
        initial_h: impl Fn([usize; 3]) -> f64,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceEnthalpy3<'s, 'm>, Error> {
        poll(gate)?;
        let binding = RheologyBinding::try_new(self.material(), law, solve)?;
        let flow = self.free_surface3_from_cells_inner(domain, initial, cx)?;
        let liquid = flow.with_bound_thermal_transport(
            law, caloric, controls, initial_h, buoyancy, Some(binding), gate,
        )?;
        cx.checkpoint().map_err(|_| FreeSurfaceError3::Cancelled)?;
        Ok(liquid)
    }
}

impl<'s> MaterialFreeSurface3<'s> {
    // A conservative transfer initializes native populations through the same
    // owner before attaching its original buoyancy and rheology models.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::material::free_surface) fn with_bound_thermal_transport<'m>(
        self,
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        controls: FreeSurfaceEnthalpyControls,
        initial_h: impl Fn([usize; 3]) -> f64,
        buoyancy: Option<crate::material::transport::boussinesq::BoussinesqConfig<3>>,
        binding: Option<RheologyBinding>,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceEnthalpy3<'s, 'm>, Error> {
        self.with_enthalpy_model(law, caloric, controls, initial_h, buoyancy, binding, None, gate)
    }
}

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Current rate, physical collision coefficient and local convergence check.
    #[must_use]
    pub fn rheology_at(&self, coordinates: [usize; 3]) -> Option<&FreeSurfaceRheologySample> {
        self.location_at(coordinates)
            .and_then(|key| self.rheology_samples.get(&key))
    }

    pub(super) fn rheology_model(&self) -> Option<RheologyModel<'m>> {
        self.rheology.as_ref().map(|binding| RheologyModel {
            binding: binding.clone(),
            law: self.law,
            caloric: self.caloric,
            density: self.flow.scaling.material().density_kg_m3(),
            buoyancy: self.buoyancy,
            dx: self.flow.scaling.cell_size_m(),
            dt: self.flow.scaling.timestep_s(),
            max_speed: self.flow.scaling.flow.u_lattice * self.flow.scaling.velocity_scale,
        })
    }

    pub(super) fn rheology_stress(
        &self,
        lattice: &FreeSurface3,
        location: MassLocation3,
    ) -> Result<RheologyStress, Error> {
        Ok(RheologyStress {
            numerator: lattice.transport_shear_rate(location, 1.0)?,
            initial_tau: self
                .rheology_samples
                .get(&location)
                .map_or(self.flow.scaling.flow.tau, |sample| sample.collision_tau),
        })
    }
}
