//! Material-dependent developed pipe-wall losses in the existing wall transaction.

use super::{NetworkWallOperator, WallContact, WallDrag, WallExchangeError, WallThermal, poll};
use crate::multicomponent::dissipation::TransportLaw;
use crate::multicomponent::dissipation::evolve::fractions;
use crate::multicomponent::wall_exchange::{
    SourceCoupledDragControl, WallControl, WallUpdate, exchange_cell_with_source_drag,
};
use crate::multicomponent::{Cell, MAX_SPECIES, Mesh1d, MixtureGas, allocated, positive};
use fs_qty::{Length, Pressure, Temperature, Time};

mod annular;
mod churchill;
mod rectangular;
mod turbulent;
pub use annular::AnnularLaminarWall;
pub use rectangular::RectangularLaminarWall;
pub use turbulent::{CircularChurchillWall, CircularGnielinskiWall, CircularTurbulentWall};

/// Explicit circular, developed-laminar approximation for one complete branch.
///
/// The diameter must agree with the original mesh's constant circular area and
/// integrated volumes. The full wetted area of a cell is therefore `4 V / D`;
/// callers cannot accidentally supply the entire pipe area to every cell.
///
/// The Nusselt number is an authored correlation input, not inferred from a gas
/// label. For example, a fully developed constant-wall-temperature approximation
/// may declare 3.66 with its own admitted regime. This type does not establish
/// entrance length, negligible axial conduction, fully developed flow, or a
/// transient boundary-layer model. Zero Reynolds is admitted only if explicitly
/// included in the caller's approximation interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CircularLaminarWall {
    /// Positive diameter of the actual circular branch [m].
    pub diameter: Length,
    /// Nonnegative declared `h D / k`; zero explicitly disables wall heat flow.
    pub nusselt_number: f64,
    /// Inclusive supported Reynolds interval, contained in `[0, 2300]`.
    pub reynolds_number: [f64; 2],
    /// Positive upper Mach bound, at most 0.3 for this low-Mach wall closure.
    pub maximum_mach_number: f64,
}

impl CircularLaminarWall {
    fn admit_mesh(self, mesh: &Mesh1d) -> Result<(), WallExchangeError> {
        let [minimum, maximum] = self.reynolds_number;
        if !positive(self.diameter.value())
            || !self.nusselt_number.is_finite()
            || self.nusselt_number < 0.0
            || !minimum.is_finite()
            || minimum < 0.0
            || !positive(maximum)
            || minimum > maximum
            || maximum > 2300.0
            || !positive(self.maximum_mach_number)
            || self.maximum_mach_number > 0.3
        {
            return Err(WallExchangeError::Invalid(
                "circular wall geometry or applicability",
            ));
        }
        admit_circular_mesh(self.diameter, mesh)
    }
}

/// An explicitly selected developed-flow model on one circular branch.
///
/// Laminar and turbulent branches can share one original transport owner.
/// Regimes are never selected automatically. Churchill explicitly selects a
/// continuous engineering relation; the other models retain their own domains.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CircularWallModel {
    /// Circular Poiseuille drag and an authored Nusselt number.
    Laminar(CircularLaminarWall),
    /// Colebrook-White drag and an authored Nusselt number.
    Turbulent(CircularTurbulentWall),
    /// Smooth-pipe drag and heat transfer from actual Reynolds and Prandtl numbers.
    Gnielinski(CircularGnielinskiWall),
    /// Explicit all-regime pipe friction and an independently authored Nusselt number.
    Churchill(CircularChurchillWall),
}

impl CircularWallModel {
    fn admit_mesh(self, mesh: &Mesh1d) -> Result<(), WallExchangeError> {
        match self {
            Self::Laminar(model) => model.admit_mesh(mesh),
            Self::Turbulent(model) => {
                model.validate()?;
                admit_circular_mesh(model.diameter, mesh)
            }
            Self::Churchill(model) => {
                model.validate()?;
                admit_circular_mesh(model.diameter, mesh)
            }
            Self::Gnielinski(model) => {
                model.validate()?;
                admit_circular_mesh(model.diameter, mesh)?;
                let positions = mesh.positions();
                let length = positions[positions.len() - 1] - positions[0];
                admit(
                    "pipe length / diameter",
                    length / model.diameter.value(),
                    10.0,
                    1e6,
                )
            }
        }
    }

    fn parameters(self) -> (Length, [f64; 2], f64) {
        match self {
            Self::Laminar(model) => (
                model.diameter,
                model.reynolds_number,
                model.maximum_mach_number,
            ),
            Self::Turbulent(model) => (
                model.diameter,
                model.reynolds_number,
                model.maximum_mach_number,
            ),
            Self::Gnielinski(model) => (
                model.diameter,
                model.reynolds_number,
                model.maximum_mach_number,
            ),
            Self::Churchill(model) => (
                model.diameter,
                model.reynolds_number,
                model.maximum_mach_number,
            ),
        }
    }
}

fn admit_circular_mesh(diameter: Length, mesh: &Mesh1d) -> Result<(), WallExchangeError> {
    let area = core::f64::consts::FRAC_PI_4 * diameter.value() * diameter.value();
    let Some(mesh_area) = mesh.uniform_cross_section() else {
        return Err(WallExchangeError::Invalid(
            "developed circular wall requires constant area",
        ));
    };
    if !positive(area) || !same_geometry(area, mesh_area) {
        return Err(WallExchangeError::Invalid(
            "wall diameter differs from original branch area",
        ));
    }
    for (faces, &volume) in mesh.positions().windows(2).zip(mesh.volumes()) {
        if !same_geometry(volume, mesh_area * (faces[1] - faces[0])) {
            return Err(WallExchangeError::Invalid(
                "circular wall volume differs from original geometry",
            ));
        }
    }
    Ok(())
}

/// Explicit geometry and developed-flow law for one original gas branch.
///
/// Circular, rectangular and concentric annular passages share one source-owned
/// transport binding. Each branch retains its own shape, heat relation and applicability.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DevelopedWallModel {
    /// Circular Poiseuille, Colebrook, Gnielinski or Churchill wall model.
    Circular(CircularWallModel),
    /// Full no-slip rectangular Poiseuille shape and an authored Nusselt number.
    Rectangular(RectangularLaminarWall),
    /// Concentric no-slip annular Poiseuille flow and an authored Nusselt number.
    Annular(AnnularLaminarWall),
}

impl DevelopedWallModel {
    fn admit_mesh(self, mesh: &Mesh1d) -> Result<(), WallExchangeError> {
        match self {
            Self::Circular(model) => model.admit_mesh(mesh),
            Self::Rectangular(model) => model.admit_mesh(mesh),
            Self::Annular(model) => model.admit_mesh(mesh),
        }
    }

    fn parameters(self) -> (Length, [f64; 2], f64) {
        match self {
            Self::Circular(model) => model.parameters(),
            Self::Rectangular(model) => (
                model.hydraulic_diameter(),
                model.reynolds_number(),
                model.maximum_mach_number(),
            ),
            Self::Annular(model) => (
                model.hydraulic_diameter(),
                model.reynolds_number(),
                model.maximum_mach_number(),
            ),
        }
    }
}

pub(super) struct DevelopedWallTransport {
    law: Box<dyn TransportLaw + Send + Sync>,
    branches: Vec<Option<DevelopedWallModel>>,
    // Preserve existing typed slice inspectors. Views are prepared only for
    // homogeneous bindings; mixed networks retain only their complete models.
    circular: Option<Vec<Option<CircularWallModel>>>,
    rectangular: Option<Vec<Option<RectangularLaminarWall>>>,
    annular: Option<Vec<Option<AnnularLaminarWall>>>,
    binary: Vec<f64>,
    source_coupled_drag: Option<SourceCoupledDragControl>,
}

impl NetworkWallOperator<'_> {
    /// Bind original transport data to selected developed circular branches.
    ///
    /// `Some` replaces that branch's authored drag and conductance coefficients
    /// with `32 mu/(rho D²)` and `Nu k (4 V/D)/D`. Wall temperatures, finite heat
    /// capacities, and adiabatic choices retain their existing ownership. `None`
    /// leaves the authored branch unchanged; plenums retain their explicit
    /// contacts because they have no developed axial flow.
    ///
    /// The same original `TransportLaw` used by molecular flow supplies current
    /// viscosity and conductivity. In particular, `MaterialTransport` retains
    /// source cards and mixture rules without copying gas-property scalars here.
    /// Full original species/caloric identity must match, not only species names.
    ///
    /// Every source half resolves the actual conservative mixture's T, p and
    /// mole fractions. By default coefficients are frozen over that half
    /// (first-order local constitutive lag); [`Self::with_source_coupled_drag`]
    /// instead resolves resistance throughout mechanical relaxation. Both
    /// endpoints must satisfy the original law and this explicit applicability
    /// domain. Ordinary `advance_for`, reacting, dissipating and source-composed
    /// network steps all use this binding.
    ///
    /// Bind before the first accepted wall step. This consumes the law but no
    /// gas/wall state; it does not silently change an existing accepted report.
    pub fn with_circular_laminar_transport<L: TransportLaw + Send + Sync + 'static>(
        self,
        law: L,
        branches: Vec<Option<CircularLaminarWall>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        let mut models = allocated(branches.len(), None)?;
        for (model, branch) in models.iter_mut().zip(branches) {
            poll(&mut cancelled)?;
            *model = branch.map(CircularWallModel::Laminar);
        }
        self.with_circular_transport(law, models, cancelled)
    }

    /// Resolve laminar or turbulent wall friction from each branch's actual gas.
    ///
    /// The original transport law supplies current viscosity and conductivity;
    /// the selected model supplies geometry, roughness or the smooth-pipe heat
    /// relation, and an explicit applicability domain. The same wall transaction retains heat,
    /// moving-wall work, finite-support recoil and whole-network rollback.
    /// Coefficients are frozen per source interval and checked at both endpoints.
    pub fn with_circular_transport<L: TransportLaw + Send + Sync + 'static>(
        self,
        law: L,
        branches: Vec<Option<CircularWallModel>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        let mut models = allocated(branches.len(), None)?;
        for (model, branch) in models.iter_mut().zip(branches) {
            poll(&mut cancelled)?;
            *model = branch.map(DevelopedWallModel::Circular);
        }
        self.with_developed_transport(law, models, cancelled)
    }

    /// Bind actual gas transport to selected developed rectangular branches.
    ///
    /// The prepared no-slip shape factor supplies `12 mu/(rho h² C)` for
    /// smaller side h, not circular drag at an equivalent diameter. Complete
    /// cell conductance is `Nu k (4 V/Dh)/Dh`; the authored rectangular Nusselt
    /// number retains its own thermal-wall applicability. `None` leaves that
    /// branch's authored contact unchanged; plenums also remain unchanged.
    ///
    /// Source identity, current composition, endpoint domain checks, moving-wall
    /// work, finite supports, heat storage and rollback use the same transaction
    /// as [`Self::with_circular_transport`]. Coefficients are frozen per source
    /// interval, so this is a quasi-steady developed-profile approximation,
    /// not a resolved unsteady transverse viscous boundary layer.
    pub fn with_rectangular_laminar_transport<L: TransportLaw + Send + Sync + 'static>(
        self,
        law: L,
        branches: Vec<Option<RectangularLaminarWall>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        let mut models = allocated(branches.len(), None)?;
        for (model, branch) in models.iter_mut().zip(branches) {
            poll(&mut cancelled)?;
            *model = branch.map(DevelopedWallModel::Rectangular);
        }
        self.with_developed_transport(law, models, cancelled)
    }

    /// Bind current gas transport to selected concentric annular passages.
    ///
    /// Prepared geometry retains exact logarithmic Poiseuille resistance, the
    /// narrow-gap limit and both wetted surfaces. Nu is an explicitly authored
    /// annular relation. Both cylinders share the contact's temperature and
    /// axial velocity; distinct surface states require a different model.
    /// Current material, support work, heat and rollback retain the original
    /// [`Self::with_developed_transport`] owner and source-interval lag.
    pub fn with_annular_laminar_transport<L: TransportLaw + Send + Sync + 'static>(
        self,
        law: L,
        branches: Vec<Option<AnnularLaminarWall>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        let mut models = allocated(branches.len(), None)?;
        for (model, branch) in models.iter_mut().zip(branches) {
            poll(&mut cancelled)?;
            *model = branch.map(DevelopedWallModel::Annular);
        }
        self.with_developed_transport(law, models, cancelled)
    }

    /// Bind one original transport source to mixed circular, rectangular and annular ducts.
    ///
    /// Each selected branch retains its own geometry and admissible flow regime.
    /// `None` preserves authored contact coefficients. The complete network uses
    /// one conservative wall transaction; a refusal on any branch rolls back all
    /// branches and finite walls. Bind before the first accepted wall step.
    pub fn with_developed_transport<L: TransportLaw + Send + Sync + 'static>(
        mut self,
        law: L,
        branches: Vec<Option<DevelopedWallModel>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        if self.transport.is_some() || self.last_report.is_some() || self.trial_active {
            return Err(WallExchangeError::Invalid(
                "wall transport must bind before evolution",
            ));
        }
        if law.source() != self.gas.source()
            || branches.len() != self.branches.len()
            || !branches.iter().any(Option::is_some)
        {
            return Err(WallExchangeError::Invalid(
                "wall transport original source or branch axis",
            ));
        }
        let (mut only_circular, mut only_rectangular, mut only_annular) = (true, true, true);
        for (model, shape) in branches.iter().zip(&self.branches) {
            poll(&mut cancelled)?;
            if let Some(model) = model {
                model.admit_mesh(&shape.mesh)?;
                match model {
                    DevelopedWallModel::Circular(_) => {
                        only_rectangular = false;
                        only_annular = false;
                    }
                    DevelopedWallModel::Rectangular(_) => {
                        only_circular = false;
                        only_annular = false;
                    }
                    DevelopedWallModel::Annular(_) => {
                        only_circular = false;
                        only_rectangular = false;
                    }
                }
            }
        }
        let mut circular = if only_circular {
            Some(allocated(branches.len(), None)?)
        } else {
            None
        };
        let mut rectangular = if only_rectangular {
            Some(allocated(branches.len(), None)?)
        } else {
            None
        };
        let mut annular = if only_annular {
            Some(allocated(branches.len(), None)?)
        } else {
            None
        };
        for (index, model) in branches.iter().enumerate() {
            poll(&mut cancelled)?;
            match (model, &mut circular, &mut rectangular, &mut annular) {
                (Some(DevelopedWallModel::Circular(model)), Some(view), _, _) => {
                    view[index] = Some(*model)
                }
                (Some(DevelopedWallModel::Rectangular(model)), _, Some(view), _) => {
                    view[index] = Some(*model)
                }
                (Some(DevelopedWallModel::Annular(model)), _, _, Some(view)) => {
                    view[index] = Some(*model)
                }
                _ => {}
            }
        }
        let n = self.gas.species_count();
        let binary = allocated(n * n, 0.0)?;
        poll(&mut cancelled)?;
        self.transport = Some(DevelopedWallTransport {
            law: Box::new(law),
            branches,
            circular,
            rectangular,
            annular,
            binary,
            source_coupled_drag: None,
        });
        Ok(self)
    }

    /// Resolve original transport properties throughout each selected wall's
    /// mechanical relaxation instead of freezing its initial drag coefficients.
    ///
    /// Bind after a developed-transport law and before evolution. Every
    /// mechanical trial reads the original mixture's actual T, p and composition,
    /// including slip heating and prescribed-wall work, and re-admits the
    /// declared Reynolds/Mach/source domain. The source controls apply separately
    /// to each selected cell in each source interval; initial and final transport
    /// admission remain additional bounded queries. The native source kernel
    /// exposes its refinement receipt; the network's `maximum_iterations` still
    /// counts only the original thermal solve.
    ///
    /// The resulting mechanical endpoint refreshes the local gas film before
    /// the existing caloric/conduction heat solve. Heat still follows mechanics
    /// with the original first-order composition and a held gas film; shared
    /// thermal ports retain their own authored conductances. No simultaneous
    /// thermoviscous heat integration or unsteady profile memory is inferred.
    ///
    /// Selected walls must be stationary or have a prescribed constant axial
    /// velocity. Finite independent/shared supports on selected branches refuse
    /// in either binding order; their coupled recoil equations are not replaced
    /// by a prescribed-wall solve. Supports on unselected branches are retained.
    pub fn with_source_coupled_drag(
        mut self,
        control: SourceCoupledDragControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        control.validate()?;
        if self.last_report.is_some() || self.trial_active {
            return Err(WallExchangeError::Invalid(
                "source-coupled wall drag must bind before evolution",
            ));
        }
        let transport = self.transport.as_mut().ok_or(WallExchangeError::Invalid(
            "source-coupled wall drag requires original developed transport",
        ))?;
        if transport.source_coupled_drag.is_some()
            || !transport.branches.iter().any(Option::is_some)
        {
            return Err(WallExchangeError::Invalid(
                "source-coupled wall drag binds once to selected branches",
            ));
        }
        for (branch, (model, shape)) in transport.branches.iter().zip(&self.branches).enumerate() {
            poll(&mut cancelled)?;
            if model.is_some()
                && (self.supports.as_ref().is_some_and(|supports| {
                    supports.any_selected(shape.start, shape.mesh.cell_count())
                }) || self
                    .shared_supports
                    .as_ref()
                    .is_some_and(|supports| supports.has_branch(branch)))
            {
                return Err(WallExchangeError::Invalid(
                    "source-coupled wall drag requires prescribed wall motion",
                ));
            }
        }
        poll(&mut cancelled)?;
        transport.source_coupled_drag = Some(control);
        Ok(self)
    }

    /// Per-cell mechanical source accuracy/work controls, when explicitly bound.
    pub fn source_coupled_drag_control(&self) -> Option<SourceCoupledDragControl> {
        self.transport.as_ref()?.source_coupled_drag
    }

    /// Original immutable constitutive owner, including its species source.
    pub fn transport_law(&self) -> Option<&dyn TransportLaw> {
        self.transport
            .as_ref()
            .map(|transport| transport.law.as_ref() as &dyn TransportLaw)
    }

    /// Original circular models for a homogeneous circular binding, otherwise None.
    pub fn circular_walls(&self) -> Option<&[Option<CircularWallModel>]> {
        self.transport.as_ref()?.circular.as_deref()
    }

    /// Original rectangular models for a homogeneous rectangular binding, otherwise None.
    pub fn rectangular_walls(&self) -> Option<&[Option<RectangularLaminarWall>]> {
        self.transport.as_ref()?.rectangular.as_deref()
    }

    /// Original annular models for a homogeneous annular binding, otherwise None.
    pub fn annular_walls(&self) -> Option<&[Option<AnnularLaminarWall>]> {
        self.transport.as_ref()?.annular.as_deref()
    }

    /// Complete original geometry and flow model in branch order, including mixtures.
    pub fn developed_walls(&self) -> Option<&[Option<DevelopedWallModel>]> {
        self.transport
            .as_ref()
            .map(|transport| transport.branches.as_slice())
    }
}

impl DevelopedWallTransport {
    pub(super) fn couples_drag(&self, branch: usize) -> bool {
        self.source_coupled_drag.is_some() && self.branches[branch].is_some()
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_drag(
        &mut self,
        gas: &MixtureGas,
        branch: usize,
        state: Cell<'_>,
        volume: f64,
        duration: Time,
        contact: WallContact,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<WallUpdate>, WallExchangeError> {
        let Some(source_control) = self.source_coupled_drag else {
            return Ok(None);
        };
        if self.branches[branch].is_none() {
            return Ok(None);
        }
        let mut update = exchange_cell_with_source_drag(
            gas,
            state,
            volume,
            duration,
            contact,
            controls,
            source_control,
            |candidate, checkpoint| {
                self.resolve(gas, branch, candidate, volume, contact, &mut || {
                    checkpoint()
                })
                .map(|resolved| resolved.drag)
            },
            &mut *cancelled,
        )?
        .gas;
        // The mechanical source owns no thermal storage. Refresh only the
        // actual endpoint film/drag before the original thermal owner takes over.
        update.contact = self.resolve(
            gas,
            branch,
            Cell {
                momentum: update.momentum,
                energy: update.energy,
                ..state
            },
            volume,
            contact,
            cancelled,
        )?;
        Ok(Some(update))
    }

    pub(super) fn resolve(
        &mut self,
        gas: &MixtureGas,
        branch: usize,
        state: Cell<'_>,
        volume: f64,
        mut contact: WallContact,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<WallContact, WallExchangeError> {
        let Some(model) = self.branches[branch] else {
            return Ok(contact);
        };
        poll(cancelled)?;
        if self.law.source() != gas.source() {
            return Err(WallExchangeError::Invalid(
                "wall transport changed its original source",
            ));
        }
        let primitive = gas.primitive(state, &mut *cancelled)?;
        let mut mole_fractions = [0.0; MAX_SPECIES];
        let x = &mut mole_fractions[..gas.species_count()];
        fractions(gas, state.partial_densities, x)?;
        let coefficients = self.law.evaluate(
            Temperature::new(primitive.temperature),
            Pressure::new(primitive.pressure),
            x,
            &mut self.binary,
            cancelled,
        )?;
        let viscosity = coefficients.dynamic_viscosity.value();
        let conductivity = coefficients.thermal_conductivity.value();
        if !positive(viscosity) || !conductivity.is_finite() || conductivity < 0.0 {
            return Err(WallExchangeError::Invalid(
                "wall law viscosity or conductivity",
            ));
        }
        let (diameter, reynolds_number, maximum_mach_number) = model.parameters();
        let motion = contact.drag.moving_velocity();
        let slip = primitive.velocity - motion.unwrap_or(0.0);
        let reynolds = primitive.density * slip.abs() * diameter.value() / viscosity;
        let mach = slip.abs() / primitive.sound_speed;
        let [minimum, maximum] = reynolds_number;
        admit("Reynolds number", reynolds, minimum, maximum)?;
        admit("Mach number", mach, 0.0, maximum_mach_number)?;
        let (friction_factor, nusselt_number) = match model {
            DevelopedWallModel::Circular(CircularWallModel::Laminar(laminar)) => {
                (None, laminar.nusselt_number)
            }
            DevelopedWallModel::Circular(CircularWallModel::Churchill(churchill)) => {
                (None, churchill.nusselt_number)
            }
            DevelopedWallModel::Circular(CircularWallModel::Turbulent(turbulent)) => (
                Some(turbulent.darcy_factor(reynolds)?),
                turbulent.nusselt_number,
            ),
            DevelopedWallModel::Circular(CircularWallModel::Gnielinski(thermal)) => {
                // Cp, not Cv, defines Pr. Use the original NASA mixture at the
                // current composition and temperature, preserving formation-
                // reference independence and avoiding a constant-gamma proxy.
                let (_, cv_density) =
                    gas.caloric(state.partial_densities, primitive.temperature)?;
                let cp =
                    (cv_density + primitive.pressure / primitive.temperature) / primitive.density;
                let prandtl = viscosity * cp / conductivity;
                let (factor, nu) = thermal.coefficients(reynolds, prandtl)?;
                (Some(factor), nu)
            }
            DevelopedWallModel::Rectangular(rectangular) => (None, rectangular.nusselt_number()),
            DevelopedWallModel::Annular(annular) => (None, annular.nusselt_number()),
        };
        contact.drag = match model {
            DevelopedWallModel::Circular(CircularWallModel::Churchill(churchill)) => {
                let WallDrag::Linear { rate_per_s } = WallDrag::circular_laminar(
                    coefficients.dynamic_viscosity,
                    primitive.density,
                    diameter,
                )?
                else {
                    unreachable!("Poiseuille resistance is linear")
                };
                let excess = churchill.excess(reynolds)?;
                // Retain the exact laminar term. The nonnegative remainder
                // satisfies a*u + b*u*|u| = f_D(Re)*u*|u|/(2D) at this state,
                // while remaining finite and exactly linear at rest.
                let quadratic_per_m = if reynolds == 0.0 {
                    0.0
                } else {
                    32.0 * (excess / reynolds) / diameter.value()
                };
                if !quadratic_per_m.is_finite() {
                    return Err(WallExchangeError::Invalid(
                        "unrepresentable Churchill resistance",
                    ));
                }
                WallDrag::LinearQuadratic {
                    linear_rate_per_s: rate_per_s,
                    quadratic_per_m,
                    wall_velocity_m_s: motion.unwrap_or(0.0),
                }
            }
            DevelopedWallModel::Circular(CircularWallModel::Laminar(_)) => {
                let drag = WallDrag::circular_laminar(
                    coefficients.dynamic_viscosity,
                    primitive.density,
                    diameter,
                )?;
                match (drag, motion) {
                    (WallDrag::Linear { rate_per_s }, Some(wall_velocity_m_s)) => {
                        WallDrag::MovingLinear {
                            rate_per_s,
                            wall_velocity_m_s,
                        }
                    }
                    _ => drag,
                }
            }
            DevelopedWallModel::Circular(
                CircularWallModel::Turbulent(_) | CircularWallModel::Gnielinski(_),
            ) => {
                let friction_factor =
                    friction_factor.expect("turbulent branch evaluated its factor");
                match motion {
                    Some(wall_velocity_m_s) => WallDrag::MovingDarcy {
                        friction_factor,
                        hydraulic_diameter_m: diameter.value(),
                        wall_velocity_m_s,
                    },
                    None => WallDrag::Darcy {
                        friction_factor,
                        hydraulic_diameter_m: diameter.value(),
                    },
                }
            }
            DevelopedWallModel::Rectangular(rectangular) => {
                let rate_per_s = rectangular.linear_rate(viscosity, primitive.density)?;
                match motion {
                    Some(wall_velocity_m_s) => WallDrag::MovingLinear {
                        rate_per_s,
                        wall_velocity_m_s,
                    },
                    None => WallDrag::Linear { rate_per_s },
                }
            }
            DevelopedWallModel::Annular(annular) => {
                let rate_per_s = annular.linear_rate(viscosity, primitive.density)?;
                match motion {
                    Some(wall_velocity_m_s) => WallDrag::MovingLinear {
                        rate_per_s,
                        wall_velocity_m_s,
                    },
                    None => WallDrag::Linear { rate_per_s },
                }
            }
        };
        // The original source solver owns wall storage and energy exchange;
        // this binding supplies only its instantaneous constitutive coefficients.
        let conductance =
            nusselt_number * conductivity * (4.0 * volume / diameter.value()) / diameter.value();
        if !conductance.is_finite() || conductance < 0.0 {
            return Err(WallExchangeError::Invalid(
                "wall conductance is unrepresentable",
            ));
        }
        match &mut contact.thermal {
            WallThermal::Adiabatic => {}
            WallThermal::Reservoir {
                conductance_w_per_k,
                ..
            }
            | WallThermal::Lumped {
                conductance_w_per_k,
                ..
            }
            | WallThermal::Caloric {
                conductance_w_per_k,
                ..
            } => {
                *conductance_w_per_k = conductance;
            }
        }
        poll(cancelled)?;
        Ok(contact)
    }
}

fn same_geometry(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}

fn admit(
    quantity: &'static str,
    value: f64,
    minimum: f64,
    maximum: f64,
) -> Result<(), WallExchangeError> {
    if value.is_finite() && (minimum..=maximum).contains(&value) {
        Ok(())
    } else {
        Err(WallExchangeError::TransportDomain {
            quantity,
            value,
            minimum,
            maximum,
        })
    }
}
