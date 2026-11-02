//! Isothermal source rheology on the actual sparse VOF trajectory.
//!
//! The same local constitutive solve as thermal flow closes physical strain
//! rate against source viscosity. Native collision, streaming, mass exchange,
//! conversion and uniform capillary reconstruction remain the numerical owners.
//! No caloric property, viscosity floor or fitted material-name law is invented.

use std::collections::BTreeMap;

use fs_exec::{CancelGate, Cx};
use fs_matdb::MaterialCard;
use fs_material::liquid::surface::{
    LiquidCapillarity, ResolvedLiquidSurfaceState, SURFACE_TENSION_DIMS, SURFACE_TENSION_PROPERTY,
};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{MaterialStatePointError, ResolvedMaterialStatePoint};
use fs_qty::QuantitySpec;

use super::{
    MaterialFreeSurface3, MaterialFreeSurfaceDomain3, MaterialFreeSurfaceInitialCell3,
    MaterialFreeSurfaceSample3, capillary_scaling, plan_free_surface_from_material,
};
use crate::d3q19::freesurface3::pressure::GasPressureField3;
use crate::d3q19::freesurface3::transport::MassLocation3;
use crate::d3q19::freesurface3::{FreeSurface3, FreeSurfaceError3};
use crate::material::transport::rheology::{LocalRheologyError, solve_local};
use crate::material::{MaterialLbmError, plan_momentum_scaling};

pub use crate::material::transport::rheology::RheologySolve;

/// Source, native flow or bounded constitutive solve refused this candidate.
#[derive(Clone, Debug, PartialEq)]
pub enum FreeSurfaceRheologyError {
    Material(MaterialStatePointError),
    Flow(MaterialLbmError),
    Convergence {
        coordinates: [usize; 3],
        iterations: usize,
        relaxation_defect: f64,
    },
    Invalid(&'static str),
}
impl core::fmt::Display for FreeSurfaceRheologyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(e) => write!(f, "isothermal free-surface source: {e}"),
            Self::Flow(e) => write!(f, "isothermal free-surface flow: {e}"),
            Self::Convergence {
                coordinates,
                iterations,
                relaxation_defect,
            } => write!(
                f,
                "free-surface rheology at {coordinates:?} did not converge in {iterations} queries (tau defect {relaxation_defect})",
            ),
            Self::Invalid(e) => write!(f, "isothermal free-surface rheology: {e}"),
        }
    }
}
impl std::error::Error for FreeSurfaceRheologyError {}
impl From<MaterialStatePointError> for FreeSurfaceRheologyError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Material(e)
    }
}
impl From<MaterialLbmError> for FreeSurfaceRheologyError {
    fn from(e: MaterialLbmError) -> Self {
        Self::Flow(e)
    }
}
impl From<FreeSurfaceError3> for FreeSurfaceRheologyError {
    fn from(e: FreeSurfaceError3) -> Self {
        Self::Flow(e.into())
    }
}
type Error = FreeSurfaceRheologyError;

/// One immutable material assignment at fixed temperature and external context.
///
/// The default law evolves only nonnegative physical shear rate. The explicit
/// [`Self::try_new_with_pressure`] law additionally resolves viscosity at the
/// actual absolute mechanical pressure. Density and surface tension cannot
/// depend on either evolving coordinate. Temperature, composition, environment
/// and source selections stay fixed. This is an isothermal approximation, not
/// a finite heat reservoir or a compressible liquid equation of state.
pub struct IsothermalRheologyLaw<'m> {
    card: &'m MaterialCard,
    reference: ResolvedLiquidSurfaceState,
    temperature_axis: String,
    rate_axis: String,
    pressure_axis: Option<String>,
    surface_tension_n_m: Option<f64>,
}

impl<'m> IsothermalRheologyLaw<'m> {
    /// Adopt the same complete material bundle already used by other physics.
    /// Thermal conductivity and heat capacity are not required. The viscosity
    /// basis and the decision to include capillarity are explicit scenario data.
    /// Other dynamic response/applicability axes refuse through the shared
    /// material admission. Every later query pins the SAME selected claims.
    pub fn try_new(
        card: &'m MaterialCard,
        material: ResolvedMaterialStatePoint,
        viscosity_basis: LiquidViscosityBasis,
        capillarity: LiquidCapillarity,
        temperature_axis: &str,
        shear_rate_axis: &str,
    ) -> Result<Self, Error> {
        Self::try_new_inner(
            card,
            material,
            viscosity_basis,
            capillarity,
            temperature_axis,
            shear_rate_axis,
            None,
        )
    }

    /// Evolve the same local shear law at solved absolute liquid pressure [Pa].
    ///
    /// This opt-in uses `p_ref + cs²*(rho_lattice-1)*rho_source*(dx/dt)²`
    /// solely as the mechanical pressure coordinate of the viscosity source.
    /// Physical density stays fixed; pressure-dependent density or surface
    /// tension refuses rather than implying compressibility or unsupported
    /// surface-stress evolution. No analytic pressure-viscosity law or floor
    /// substitutes for the original material data.
    ///
    /// The flow must be constructed with an explicit pressure reference equal
    /// to this bundle's initial absolute pressure. Every accepted endpoint
    /// re-queries the SAME source selections at its own pressure and shear rate.
    /// Therefore a physical common shift of pressure/reference can change
    /// viscosity; it is not a gauge symmetry of a pressure-dependent material.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_pressure(
        card: &'m MaterialCard,
        material: ResolvedMaterialStatePoint,
        viscosity_basis: LiquidViscosityBasis,
        capillarity: LiquidCapillarity,
        temperature_axis: &str,
        shear_rate_axis: &str,
        pressure_axis: &str,
    ) -> Result<Self, Error> {
        Self::try_new_inner(
            card,
            material,
            viscosity_basis,
            capillarity,
            temperature_axis,
            shear_rate_axis,
            Some(pressure_axis),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn try_new_inner(
        card: &'m MaterialCard,
        material: ResolvedMaterialStatePoint,
        viscosity_basis: LiquidViscosityBasis,
        capillarity: LiquidCapillarity,
        temperature_axis: &str,
        shear_rate_axis: &str,
        pressure_axis: Option<&str>,
    ) -> Result<Self, Error> {
        if material.card_identity() != card.content_hash() {
            return Err(Error::Invalid(
                "rheology requires the original material card",
            ));
        }
        let reference = match pressure_axis {
            Some(axis) => ResolvedLiquidSurfaceState::try_for_pressure_shear_rate_evolution(
                material,
                viscosity_basis,
                LiquidCapillarity::Neglected,
                temperature_axis,
                shear_rate_axis,
                axis,
            )?,
            None => ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
                material,
                viscosity_basis,
                temperature_axis,
                shear_rate_axis,
            )?,
        };
        let material = reference.material();
        // The reference-density momentum/VOF owner has no compressible storage
        // law, and its capillary boundary has no variable-sigma traction law.
        for name in ["density", SURFACE_TENSION_PROPERTY] {
            if name == SURFACE_TENSION_PROPERTY && capillarity == LiquidCapillarity::Neglected {
                continue;
            }
            if material.property(name).is_some_and(|property| {
                property
                    .source_axes()
                    .iter()
                    .any(|(axis, _)| axis == shear_rate_axis || pressure_axis == Some(axis.as_str()))
            }) {
                return Err(Error::Invalid(
                    "isothermal VOF needs density and surface tension independent of evolving shear and pressure",
                ));
            }
        }
        let surface_tension_n_m = surface_tension(material, capillarity)?;
        Ok(Self {
            card,
            reference,
            temperature_axis: temperature_axis.to_owned(),
            rate_axis: shear_rate_axis.to_owned(),
            pressure_axis: pressure_axis.map(str::to_owned),
            surface_tension_n_m,
        })
    }

    /// Complete original source bundle, including any sourced surface tension.
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        self.reference.material()
    }

    /// Fixed sourced capillary coefficient [N/m], or explicitly neglected.
    #[must_use]
    pub const fn surface_tension_n_m(&self) -> Option<f64> {
        self.surface_tension_n_m
    }

    /// Named source coordinate evolved from solved absolute pressure, if opted in.
    #[must_use]
    pub fn pressure_axis(&self) -> Option<&str> {
        self.pressure_axis.as_deref()
    }

    fn resolve(
        &self,
        rate: f64,
        pressure: Option<f64>,
    ) -> Result<ResolvedLiquidSurfaceState, Error> {
        let mut point = self
            .material()
            .query()
            .with(&self.rate_axis, rate)
            .map_err(|_| Error::Invalid("unrepresentable physical shear coordinate"))?;
        if let Some(axis) = &self.pressure_axis {
            let pressure = pressure.ok_or(Error::Invalid("missing solved absolute pressure"))?;
            if !pressure.is_finite() || pressure <= 0.0 {
                return Err(Error::Invalid(
                    "solved absolute pressure must be finite and positive",
                ));
            }
            point = point
                .with(axis, pressure)
                .map_err(|_| Error::Invalid("unrepresentable absolute pressure coordinate"))?;
        }
        let material = self.material().at_query_point(self.card, &point)?;
        let resolved = match &self.pressure_axis {
            Some(axis) => ResolvedLiquidSurfaceState::try_for_pressure_shear_rate_evolution(
                material,
                self.reference.viscosity_basis(),
                LiquidCapillarity::Neglected,
                &self.temperature_axis,
                &self.rate_axis,
                axis,
            )?,
            None => ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
                material,
                self.reference.viscosity_basis(),
                &self.temperature_axis,
                &self.rate_axis,
            )?,
        };
        let capillarity = if self.surface_tension_n_m.is_some() {
            LiquidCapillarity::Sourced
        } else {
            LiquidCapillarity::Neglected
        };
        if resolved.density_kg_m3().to_bits() != self.reference.density_kg_m3().to_bits()
            || surface_tension(resolved.material(), capillarity)? != self.surface_tension_n_m
        {
            return Err(Error::Invalid(
                "fixed density or capillarity changed during rheology",
            ));
        }
        Ok(resolved)
    }
}

/// Actual physical constitutive coordinates at an accepted wet collision cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsothermalRheologySample {
    /// `sqrt(2 D:D)` [1/s], inferred from the force-corrected native stress.
    pub shear_rate_per_s: f64,
    /// Molecular `1/2 + 3 nu dt/dx^2`, used by the next collision.
    pub collision_tau: f64,
    pub iterations: usize,
    /// Absolute local constitutive defect in dimensionless tau.
    pub relaxation_defect: f64,
}

struct CellMaterial {
    material: ResolvedLiquidSurfaceState,
    sample: IsothermalRheologySample,
}

/// The existing sparse VOF with source viscosity at each actual wet-cell rate.
///
/// The immutable law, original physical units, native mass ledger and accepted
/// source samples travel together. Source, speed, convergence or cancellation
/// refusals preserve populations, fill, conversions, mass carry and receipts.
/// This is explicit generalized-Newtonian BGK evolution. It does not evolve
/// heat, viscoelastic memory, yield history, phase change or surface energy.
pub struct MaterialFreeSurfaceRheology3<'s, 'm> {
    flow: MaterialFreeSurface3<'s>,
    law: &'s IsothermalRheologyLaw<'m>,
    solve: RheologySolve,
    samples: BTreeMap<MassLocation3, CellMaterial>,
}

impl<'s, 'm> MaterialFreeSurfaceRheology3<'s, 'm> {
    /// Initialize from supplied fill and physical velocity using the native
    /// constructor. Source coefficients are resolved against every initial wet
    /// cell's actual stress, including zero-mass interface collision ghosts.
    /// Domain allocation is bounded by `domain.max_cells`; local constitutive
    /// work is bounded by `solve.maximum_iterations` at each occupied cell.
    #[allow(clippy::too_many_arguments)]
    pub fn from_cells(
        law: &'s IsothermalRheologyLaw<'m>,
        cell_size_m: f64,
        timestep_s: f64,
        max_speed_m_s: f64,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        solve: RheologySolve,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        Self::from_cells_inner(
            law,
            cell_size_m,
            timestep_s,
            max_speed_m_s,
            domain,
            initial,
            solve,
            None,
            cx,
            gate,
        )
    }

    /// Initialize rheology with a fixed absolute-pressure datum [Pa].
    ///
    /// Uses the original material scaling's pressure-reference owner. Numerical
    /// density one represents `reference_pressure_pa`; initial populations,
    /// velocities and mass are otherwise identical to [`Self::from_cells`].
    /// This does not initialize a hydrostatic or capillary pressure field.
    /// The caller owns initial pressure/load compatibility. For the default law
    /// the source pressure context stays fixed. An explicitly pressure-dependent
    /// law requires this datum to equal its initial source pressure and samples
    /// the actual local pressure on every accepted state. Numerical pressure is
    /// not a liquid equation of state. The datum cannot be rebound later.
    ///
    /// # Errors
    /// The ordinary ingress/source gates apply. Nonpositive or nonfinite
    /// reference pressure refuses before returning a candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn from_cells_with_pressure_reference(
        law: &'s IsothermalRheologyLaw<'m>,
        cell_size_m: f64,
        timestep_s: f64,
        max_speed_m_s: f64,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        solve: RheologySolve,
        reference_pressure_pa: f64,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        Self::from_cells_inner(
            law,
            cell_size_m,
            timestep_s,
            max_speed_m_s,
            domain,
            initial,
            solve,
            Some(reference_pressure_pa),
            cx,
            gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn from_cells_inner(
        law: &'s IsothermalRheologyLaw<'m>,
        cell_size_m: f64,
        timestep_s: f64,
        max_speed_m_s: f64,
        domain: MaterialFreeSurfaceDomain3<'_>,
        initial: &[MaterialFreeSurfaceInitialCell3],
        solve: RheologySolve,
        reference_pressure_pa: Option<f64>,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        poll(gate)?;
        solve
            .validate()
            .map_err(|_| Error::Invalid("invalid liquid rheology controls"))?;
        if let Some(axis) = law.pressure_axis() {
            let initial_pressure = law
                .material()
                .query_point()
                .iter()
                .find(|(name, _)| name == axis)
                .map(|(_, pressure)| *pressure);
            if reference_pressure_pa != initial_pressure {
                return Err(Error::Invalid(
                    "pressure-dependent rheology requires a fixed datum equal to the initial source pressure",
                ));
            }
        }
        let mut scaling = plan_free_surface_from_material(
            &law.reference,
            cell_size_m,
            timestep_s,
            max_speed_m_s,
        )?;
        scaling.surface_tension_lattice = capillary_scaling(
            law.surface_tension_n_m,
            scaling.pressure_unit_pa,
            cell_size_m,
        )?;
        if let Some(reference) = reference_pressure_pa {
            scaling = scaling.with_pressure_reference(reference)?;
        }
        let flow = scaling.free_surface3_from_cells_inner(domain, initial, cx)?;
        let mut candidate = Self {
            flow,
            law,
            solve,
            samples: BTreeMap::new(),
        };
        candidate.samples = candidate.sample_lattice(&candidate.flow.lattice, gate)?;
        poll(gate)?;
        cx.checkpoint().map_err(|_| FreeSurfaceError3::Cancelled)?;
        Ok(candidate)
    }

    /// Actual native collision/mass owner, exposed read-only.
    #[must_use]
    pub const fn lattice(&self) -> &FreeSurface3 {
        &self.flow.lattice
    }

    /// Source assignment, declared evolving coordinates and original claim receipts.
    #[must_use]
    pub const fn law(&self) -> &'s IsothermalRheologyLaw<'m> {
        self.law
    }

    /// Fixed absolute pressure represented by numerical density one [Pa].
    /// None retains the original pressure-perturbation-only convention.
    #[must_use]
    pub const fn reference_pressure_pa(&self) -> Option<f64> {
        self.flow.scaling.reference_pressure_pa()
    }

    /// Current complete source bundle at a wet cell; gas/carry has no rheology.
    #[must_use]
    pub fn material_at(&self, coordinates: [usize; 3]) -> Option<&ResolvedMaterialStatePoint> {
        self.samples
            .get(&FreeSurface3::transport_location(coordinates)?)
            .map(|cell| cell.material.material())
    }

    /// Accepted rate, next collision coefficient and constitutive convergence.
    #[must_use]
    pub fn rheology_at(&self, coordinates: [usize; 3]) -> Option<&IsothermalRheologySample> {
        self.samples
            .get(&FreeSurface3::transport_location(coordinates)?)
            .map(|cell| &cell.sample)
    }

    /// Original conservative VOF ledger, including unregioned carry, in kg.
    pub fn ledger_mass_kg(&self) -> Result<f64, Error> {
        Ok(self.flow.ledger_mass_kg()?)
    }

    /// Wet-cell physical density, Guo-corrected velocity and pressure perturbation.
    pub fn sample_at_si(
        &self,
        x: usize,
        y: usize,
        z: usize,
    ) -> Result<Option<MaterialFreeSurfaceSample3>, Error> {
        Ok(self.flow.sample_at_si(x, y, z)?)
    }

    /// The same necessary native capillary positivity bound as frozen flow.
    #[must_use]
    pub fn capillary_density_excursion_bound(&self) -> f64 {
        self.flow.scaling.capillary_density_excursion_bound()
    }

    /// Advance one physical timestep with current source viscosity, then admit
    /// the new native stress field before publishing either flow or sources.
    /// Endpoint queries validate the explicit step, not an implicit time solve.
    pub fn step_with_gate(&mut self, gate: &CancelGate) -> Result<(), Error> {
        self.step_inner(gate, None)
    }

    /// Advance source-resolved viscosity under prescribed gas pressure [Pa].
    ///
    /// Requires the construction-time fixed pressure reference. The callback
    /// supplies positive absolute gas pressure at exposed liquid cell centers
    /// [m], frozen over this physical timestep. It shares the original physical
    /// coordinate conversion and combined gas/capillary density-excursion bound
    /// in `(0,1)` with [`MaterialFreeSurface3::step_with_gas_pressure`].
    ///
    /// Current source-resolved local relaxation and the pressure reconstruction
    /// act in ONE original native collision/transport step. Its actual endpoint
    /// stress supplies the next bounded constitutive solve before either flow
    /// or material samples publish. No unused thermal mass routes are allocated.
    /// Zero/reference loading retains the ordinary rheology trajectory, and a
    /// common shift of reference and gas pressure leaves that trajectory intact
    /// only when the law keeps its original fixed pressure context.
    ///
    /// For the opt-in pressure-dependent law, actual endpoint liquid pressure
    /// also sets the source viscosity coordinate. The external gas pressure is
    /// never substituted for interior liquid pressure. This is an explicit
    /// pressure/shear coefficient update, not an implicit time solve. Gas
    /// evolution, physical compressibility, pressure-volume work and heat remain
    /// outside this isothermal owner.
    ///
    /// # Errors
    /// Invalid pressure or density bounds, source-domain/convergence failure,
    /// endpoint speed refusal or cancellation leaves populations, mass, carry,
    /// conversion counters and all accepted material/rate samples unchanged.
    pub fn step_with_gas_pressure(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn([f64; 3]) -> f64,
    ) -> Result<(), Error> {
        poll(gate)?;
        let reference = self
            .flow
            .scaling
            .reference_pressure_pa
            .ok_or(Error::Invalid(
                "rheology gas loading requires an explicit fixed pressure reference",
            ))?;
        let field = super::pressure::lattice_pressure_field(
            reference,
            self.flow.scaling.cell_size_m,
            self.flow.scaling.pressure_unit_pa,
            pressure_pa,
        );
        self.step_inner(
            gate,
            Some((maximum_density_excursion, GasPressureField3::Cell(&field))),
        )
    }

    /// Advance source-resolved viscosity with distinct gas pressure per link.
    ///
    /// The positive absolute pressure [Pa] is sampled on the actual native gas
    /// neighbors with the coordinate meaning of [`super::MaterialGasPressureLink3`].
    /// This shares [`Self::step_with_gas_pressure`]'s one collision, endpoint
    /// constitutive solve and atomic publication. A pressure-dependent source
    /// still uses the actual liquid pressure, never a neighboring gas load.
    /// Gas dynamics and pressure-volume heat/work are outside this owner.
    pub fn step_with_gas_link_pressure(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn(super::MaterialGasPressureLink3) -> f64,
    ) -> Result<(), Error> {
        poll(gate)?;
        let reference = self
            .flow
            .scaling
            .reference_pressure_pa
            .ok_or(Error::Invalid(
                "rheology gas loading requires an explicit fixed pressure reference",
            ))?;
        let field = super::pressure::lattice_link_pressure_field(
            reference,
            self.flow.scaling.cell_size_m,
            self.flow.scaling.pressure_unit_pa,
            pressure_pa,
        );
        self.step_inner(
            gate,
            Some((maximum_density_excursion, GasPressureField3::Link(&field))),
        )
    }

    fn step_inner(
        &mut self,
        gate: &CancelGate,
        gas_pressure: Option<(f64, GasPressureField3<'_>)>,
    ) -> Result<(), Error> {
        poll(gate)?;
        let mut tau = BTreeMap::new();
        for (&location, cell) in &self.samples {
            poll(gate)?;
            if let MassLocation3::Cell(key, lane) = location {
                tau.insert((key, lane), cell.sample.collision_tau);
            }
        }
        let mut lattice = self.flow.lattice.clone();
        lattice.step_with_material_relaxation_and_pressure(gate, &tau, gas_pressure)?;
        let samples = self.sample_lattice(&lattice, gate)?;
        poll(gate)?;
        self.flow.lattice = lattice;
        self.samples = samples;
        Ok(())
    }

    fn sample_lattice(
        &self,
        lattice: &FreeSurface3,
        gate: &CancelGate,
    ) -> Result<BTreeMap<MassLocation3, CellMaterial>, Error> {
        let mut samples = BTreeMap::new();
        let scaling = self.flow.scaling;
        for (location, _) in lattice.mass_snapshot() {
            poll(gate)?;
            let Some(coordinates) = FreeSurface3::transport_coordinates(location) else {
                continue;
            };
            let [x, y, z] = coordinates;
            let native = lattice
                .sample_at(x, y, z)
                .ok_or(Error::Invalid("missing wet-cell physical state"))?;
            let pressure = if self.law.pressure_axis().is_some() {
                Some(scaling.absolute_pressure_si(native.density)?)
            } else {
                None
            };
            let speed = native.velocity[0]
                .hypot(native.velocity[1])
                .hypot(native.velocity[2]);
            if !speed.is_finite() || speed > scaling.flow.u_lattice {
                return Err(Error::Invalid(
                    "rheology flow exceeds admitted scenario speed",
                ));
            }
            let initial_tau = self
                .samples
                .get(&location)
                .map_or(scaling.flow.tau, |cell| cell.sample.collision_tau);
            let root = solve_local(
                lattice.transport_shear_rate(location, 1.0)?,
                scaling.timestep_s,
                initial_tau,
                self.solve,
                || poll(gate),
                |rate| {
                    let material = self.law.resolve(rate, pressure)?;
                    let (plan, _) = plan_momentum_scaling(
                        material.kinematic_viscosity_m2_s(),
                        scaling.cell_size_m,
                        scaling.timestep_s,
                        0.0,
                    )?;
                    Ok(((material, plan.tau), plan.tau))
                },
            )
            .map_err(|error| match error {
                LocalRheologyError::Evaluation(e) => e,
                LocalRheologyError::Invalid(what) => Error::Invalid(what),
                LocalRheologyError::Convergence {
                    iterations,
                    relaxation_defect,
                } => Error::Convergence {
                    coordinates,
                    iterations,
                    relaxation_defect,
                },
            })?;
            samples.insert(
                location,
                CellMaterial {
                    material: root.value.0,
                    sample: IsothermalRheologySample {
                        shear_rate_per_s: root.shear_rate_per_s,
                        collision_tau: root.value.1,
                        iterations: root.iterations,
                        relaxation_defect: root.relaxation_defect,
                    },
                },
            );
        }
        Ok(samples)
    }
}

fn surface_tension(
    material: &ResolvedMaterialStatePoint,
    capillarity: LiquidCapillarity,
) -> Result<Option<f64>, Error> {
    if capillarity == LiquidCapillarity::Neglected {
        return Ok(None);
    }
    let property = material
        .property(SURFACE_TENSION_PROPERTY)
        .ok_or(Error::Invalid(
            "sourced capillarity requires original surface tension",
        ))?;
    if property.requirement().quantity() != QuantitySpec::dimensional(SURFACE_TENSION_DIMS)
        || !property.value_si().is_finite()
        || property.value_si() <= 0.0
    {
        return Err(Error::Invalid(
            "surface tension must be a positive scalar in N/m",
        ));
    }
    Ok(Some(property.value_si()))
}

fn poll(gate: &CancelGate) -> Result<(), Error> {
    if gate.is_requested() {
        Err(FreeSurfaceError3::Cancelled.into())
    } else {
        Ok(())
    }
}
