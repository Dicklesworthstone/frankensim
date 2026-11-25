//! The planar nonlinear winding owner as one Component of the composition
//! contract (fs-ports).
//!
//! [`PlanarWinding`] steps its owner's circuit ([`WindingOwner`]; by default
//! [`NonlinearWindingCircuit`]) once per proposal, from the accepted state,
//! with conductors fixed at construction
//! ([`PlanarWinding::new`]) or given by laws at the temperatures its heat
//! ports receive ([`PlanarWinding::at_port_temperatures`]). Ports, in order:
//! `source<p>` for each original external input `p` of the native graph
//! (`EffortIn` for a voltage source, `FlowIn` for a current source),
//! returning the native step-average conjugate `StepRecord::y[p]`; then
//! `heat<k>` for each declared [`HeatReceiver`], receiving temperature and
//! returning minus the mean heat rate of the stranded cells, core cells and
//! native graph resistors mapped to it. Every dissipating element maps to
//! exactly one receiver. The owner is backward Euler with held inputs, so
//! every port answers [`EndOfStep::Untold`].
//!
//! The receipt takes the owner's stored change, its numerical (backward
//! Euler) dissipation and the heat that leaves through the heat ports. A
//! step with mechanical work refuses, and so does one with magnetic
//! parameter work outside the magnetic lane below: the receipt has no term
//! for either. The temperatures the heat ports receive enter
//! the generated entropy and, with laws, the conductors and graph resistors
//! of the step: a [`ConductorLaw`] per receiver for its triangles, and a
//! [`ResistorLaw`] per graph resistor through the owner's
//! `with_graph_resistors`, or none, keeping the circuit's own resistors. A
//! law with no answer refuses the step; the Component never clamps a
//! temperature.
//!
//! In the magnetic lane ([`PlanarWinding::magnetic_at_port_temperatures`]),
//! each triangle with a [`MagneticLaw`] also takes its law's linear material
//! at the temperature its receiver's port gives, at the end of the step as
//! the conductors do, and the state carries those temperatures
//! ([`PlanarWindingState::at_temperatures`]). A port's heat is then its Joule
//! and field heat less `T dS` of its law triangles, and may be negative; the
//! stored energy adds `T s`; the numerical loss adds each law triangle's
//! `-(P + (T1 - T0) s0)`, second order in the step's change of temperature,
//! with `P` the owner's parameter work. Only the saturating owner has this
//! lane, and only for linear materials. The Component declares its stored
//! entropy: the sum of `s = -dpsi/dT` over the law triangles, measured from
//! the state.
//!
//! With the `magnetic-hysteresis` feature, the owner can also be the
//! play-memory circuit (`HystereticWindingCircuit`). Its field heats every
//! magnetic triangle, conducting or not, by pinning and viscous switching, so
//! every one of its triangles maps to a receiver. Its magnetic fits are
//! constant on their source-admitted temperature support, and a step whose
//! port gives a triangle a temperature outside the support between the fit's
//! own temperature and the port's refuses ([`PlanarWindingError::Law`]).
//!
//! The Component borrows the owner's problem for `'m` and leaks nothing. A
//! system of the kernel holds `'static` members, so a caller that adds one
//! to a system gives it a problem that lives for `'static`. The state has no
//! codec: a system holding it cannot be written to a durable snapshot.

use core::any::Any;
use core::fmt;
use core::marker::PhantomData;
use core::ops::Range;
use std::collections::BTreeMap;
use std::rc::Rc;
use std::sync::Arc;

use fs_exec::Cx;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::conductor::{ResolvedConductor, ResolvedConductorState};
use fs_material::field_storage::admit_magnetic_temperature_refresh;
use fs_material::field_storage::thermal::magnetic::ThermalMagneticMaterial;
use fs_phs::circuit::CircuitGraph;
use fs_ports::CoupleError;
use fs_ports::component::{Class, Component, ContractViolation, EndOfStep, PortSpec, Proposal};

pub(crate) mod shared;
pub mod spatial;

use shared::{Account, Element, ElementRefusal, LawCell, Refusal, Reversible, same_bits};

#[cfg(feature = "magnetic-hysteresis")]
use super::hysteretic::{HystereticWindingCircuit, HystereticWindingStep};
use super::nonlinear::{
    NonlinearWindingCircuit, NonlinearWindingControls, NonlinearWindingError, NonlinearWindingStep,
    resistors::GraphResistorMaterial,
};
use super::ohmic::WindingConductor;
use crate::planar::eddy::EddyConductor;
use crate::planar::eddy::nonlinear::{NonlinearEddyControls, NonlinearEddyMagneticCell};
use crate::planar::nonlinear::IsotropicMagneticMaterial;

/// One thermal receiver: the dissipating elements whose heat leaves through
/// one heat port.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct HeatReceiver {
    /// Mesh triangles of stranded winding or core conductors.
    pub triangles: Vec<usize>,
    /// Branch indices of native graph resistors.
    pub resistors: Vec<usize>,
}

/// The conductors of one receiver's triangles at the temperature its heat
/// port received.
pub trait ConductorLaw: fmt::Debug + Send + Sync {
    /// The source-bound state of `triangle`'s conductor at `kelvin`. The
    /// owner refuses a state whose card, basis or claims are not those of the
    /// circuit's own conductor there.
    ///
    /// # Errors
    /// Why the law has no answer there, for example a temperature outside the
    /// validity of the card behind it. The step then refuses with
    /// [`PlanarWindingError::Law`]; the temperature is never clamped.
    fn conductor(&self, triangle: usize, kelvin: f64) -> Result<ResolvedConductorState, String>;

    /// Words that identify the law. They join the Component's identity.
    fn identity_words(&self, words: &mut Vec<u64>);
}

/// One native graph resistor at the temperature its heat port received.
pub trait ResistorLaw: fmt::Debug + Send + Sync {
    /// The source-bound resistor at `kelvin`, a static Ohmic scalar as the
    /// owner's `with_graph_resistors` requires.
    ///
    /// # Errors
    /// As [`ConductorLaw::conductor`].
    fn conductor(&self, kelvin: f64) -> Result<ResolvedConductor, String>;

    /// Words that identify the law. They join the Component's identity.
    fn identity_words(&self, words: &mut Vec<u64>);
}

/// One mesh triangle whose permeability follows the temperature its
/// receiver's heat port gives.
#[derive(Clone, Debug, PartialEq)]
pub struct MagneticLaw {
    /// The mesh triangle.
    pub triangle: usize,
    /// fs-material's linear magnetic law, affine in temperature. Its source
    /// is the field's own material in the triangle.
    pub law: ThermalMagneticMaterial,
}

/// The laws of [`PlanarWinding::magnetic_at_port_temperatures`].
#[derive(Clone, Debug, PartialEq)]
pub struct MagneticLaws {
    /// One per triangle that follows its temperature; every other triangle
    /// keeps the circuit's own material.
    pub triangles: Vec<MagneticLaw>,
    /// The bounds of each step's refreshed field assembly.
    pub assembly: WeightedAssemblyLimits,
}

/// A step's refreshed field materials, one per mesh triangle, and the bounds
/// of their assembly.
type Refresh = (Vec<IsotropicMagneticMaterial>, WeightedAssemblyLimits);

/// The element whose law had no answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LawElement {
    /// A stranded or core conductor's mesh triangle.
    Triangle(usize),
    /// A native graph resistor's branch.
    Resistor(usize),
    /// A play-memory field triangle, whose magnetic fit does not hold there.
    Field(usize),
    /// A [`MagneticLaw`]'s triangle: the law has no answer at the
    /// temperature, or its free energy there is not the field's.
    Permeability(usize),
}

/// Refusal of the winding owner or of its Component boundary.
#[derive(Debug)]
pub enum PlanarWindingError {
    /// The owner's refusal, unchanged.
    Owner(NonlinearWindingError),
    /// Invalid or repeated port metadata.
    Port(ContractViolation),
    /// An invalid port identifier.
    Identifier(CoupleError),
    /// A shape, a value or a construction that is refused.
    Invalid {
        /// What was refused.
        what: &'static str,
    },
    /// The step did mechanical work, or magnetic parameter work outside the
    /// magnetic lane or in a triangle without a law, which the receipt has
    /// no term for.
    UnaccountedWork {
        /// Signed magnetic parameter work \[J\].
        parameter_j: f64,
        /// Signed mechanical work \[J\].
        mechanical_j: f64,
    },
    /// Supplied work, stored change, numerical dissipation and the heat that
    /// leaves do not close within the rounding floor and the declared
    /// tolerance.
    Balance {
        /// Supplied work minus stored change minus heat out minus the
        /// owner's numerical dissipation \[J\].
        defect_j: f64,
        /// The allowance \[J\].
        tolerance_j: f64,
    },
    /// A law had no answer at the temperature its port received, a
    /// play-memory triangle's magnetic fit does not hold there, or a
    /// magnetic law's free energy is not the field's.
    Law {
        /// The heat receiver whose port gave the temperature.
        receiver: usize,
        /// The element asked for.
        element: LawElement,
        /// The temperature \[K\].
        kelvin: f64,
        /// The law's reason.
        message: String,
    },
    /// Cancellation was requested.
    Cancelled,
}

impl fmt::Display for PlanarWindingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Owner(error) => write!(f, "winding owner refused: {error}"),
            Self::Port(violation) => write!(f, "winding port refused: {violation}"),
            Self::Identifier(error) => write!(f, "winding port identifier refused: {error:?}"),
            Self::Invalid { what } => write!(f, "winding component refused: {what}"),
            Self::UnaccountedWork {
                parameter_j,
                mechanical_j,
            } => write!(
                f,
                "winding step did parameter work {parameter_j} J and mechanical work \
                 {mechanical_j} J, which the receipt has no term for"
            ),
            Self::Balance {
                defect_j,
                tolerance_j,
            } => write!(
                f,
                "winding energy account does not close: defect {defect_j} J beyond {tolerance_j} J"
            ),
            Self::Law {
                receiver,
                element,
                kelvin,
                message,
            } => write!(
                f,
                "winding law of {element:?} (receiver {receiver}) has no answer at {kelvin} K: \
                 {message}"
            ),
            Self::Cancelled => f.write_str("winding component cancelled"),
        }
    }
}

impl std::error::Error for PlanarWindingError {}

impl From<NonlinearWindingError> for PlanarWindingError {
    fn from(error: NonlinearWindingError) -> Self {
        Self::Owner(error)
    }
}

const fn invalid(what: &'static str) -> PlanarWindingError {
    PlanarWindingError::Invalid { what }
}

/// A refusal of the shared core, in this Component's words: its mesh
/// elements are triangles.
fn refused(refusal: Refusal) -> PlanarWindingError {
    match refusal {
        Refusal::Invalid(what) => invalid(what),
        Refusal::Element(element) => invalid(match element {
            ElementRefusal::NoConductor => "a receiver names a triangle that carries no conductor",
            ElementRefusal::HeatUnmapped => {
                "the owner reported heat in a triangle no receiver maps"
            }
            ElementRefusal::LawRecordMissing => "the owner's magnetic records omit a law triangle",
            ElementRefusal::NoLawElement => "the state's field has no such law triangle",
        }),
        Refusal::Port(violation) => PlanarWindingError::Port(violation),
        Refusal::Identifier(error) => PlanarWindingError::Identifier(error),
        Refusal::UnaccountedWork {
            parameter_j,
            mechanical_j,
        } => PlanarWindingError::UnaccountedWork {
            parameter_j,
            mechanical_j,
        },
        Refusal::Balance {
            defect_j,
            tolerance_j,
        } => PlanarWindingError::Balance {
            defect_j,
            tolerance_j,
        },
        Refusal::Law {
            receiver,
            element,
            kelvin,
            message,
        } => PlanarWindingError::Law {
            receiver,
            element: match element {
                Element::Resistor(branch) => LawElement::Resistor(branch),
                Element::Permeability(triangle) => LawElement::Permeability(triangle),
            },
            kelvin,
            message,
        },
    }
}

/// Each receiver's triangles and graph resistor branches.
fn maps(receivers: &[HeatReceiver]) -> Vec<(&[usize], &[usize])> {
    receivers
        .iter()
        .map(|receiver| (receiver.triangles.as_slice(), receiver.resistors.as_slice()))
        .collect()
}

/// The winding owner whose accepted circuit a [`PlanarWinding`] steps:
/// [`NonlinearWindingCircuit`], the saturating field, or, with the
/// `magnetic-hysteresis` feature, `HystereticWindingCircuit`, the play-memory
/// field. Sealed: the owners are fs-em's.
pub trait WindingOwner<'m>: owner::Sealed<'m> {}

impl<'m> WindingOwner<'m> for NonlinearWindingCircuit<'m> {}

#[cfg(feature = "magnetic-hysteresis")]
impl<'m> WindingOwner<'m> for HystereticWindingCircuit<'m> {}

mod owner {
    #[cfg(feature = "magnetic-hysteresis")]
    use super::HystereticWindingCircuit;
    #[cfg(feature = "magnetic-hysteresis")]
    use super::shared::{fitted_temperature, same_vectors};
    use super::{
        CircuitGraph, EddyConductor, GraphResistorMaterial, IsotropicMagneticMaterial,
        NonlinearEddyMagneticCell, NonlinearWindingCircuit, NonlinearWindingControls,
        NonlinearWindingError, NonlinearWindingStep, Range, Refresh, WindingConductor, fmt,
    };
    #[cfg(feature = "magnetic-hysteresis")]
    use crate::planar::eddy::hysteretic::{HystereticEddyControls, HystereticEddyStep};
    use crate::planar::eddy::nonlinear::{NonlinearEddyControls, NonlinearEddyStep};

    /// What the Component asks of its owner: one step from an accepted
    /// circuit, the rebinding of its graph resistors, and the accepted
    /// coordinates that its stored energy and state equality read.
    pub trait Sealed<'m>: Clone + Sized {
        /// The field's controls inside [`NonlinearWindingControls`].
        type FieldControls: Copy + fmt::Debug;
        /// The field's record inside [`NonlinearWindingStep`]. A proposed
        /// state keeps the step as `dyn Any`, so it borrows nothing.
        type FieldStep: 'static;
        /// Words that tell this owner's Component apart; none for the
        /// saturating owner, whose words came first.
        const LANE_WORDS: &'static [u64];

        /// The owner's `step_to_candidate` to `target_s`, in parts; with
        /// `materials`, its `step_with_materials_candidate`.
        #[allow(clippy::too_many_arguments)]
        fn step(
            &self,
            target_s: f64,
            materials: Option<Refresh>,
            core: Vec<EddyConductor>,
            windings: Vec<WindingConductor>,
            inputs: &[f64],
            controls: NonlinearWindingControls<Self::FieldControls>,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<(Self, NonlinearWindingStep<Self::FieldStep>), NonlinearWindingError>;

        /// The owner's `with_graph_resistors`.
        fn rebind_graph_resistors(
            &self,
            sources: Vec<GraphResistorMaterial>,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<Self, NonlinearWindingError>;

        /// The mesh triangle and heat \[J\] of each field cell of a step, in
        /// the owner's order.
        fn field_heat(step: &Self::FieldStep) -> Vec<(usize, f64)>;

        /// The mesh triangles whose field cells may dissipate beyond the
        /// conductors: none when only conductors heat.
        fn field_triangles(&self) -> Range<usize>;

        /// Why the field's law in `triangle` does not hold at `kelvin`.
        fn admit_field_temperature(&self, triangle: usize, kelvin: f64) -> Result<(), String>;

        /// Whether two states hold the same field memory beyond the
        /// coordinates below.
        fn same_memory(&self, other: &Self) -> bool;

        /// Circuit, field and end-mode storage, summed as the owner sums
        /// `stored_before_j`.
        fn stored_energy_j(&self) -> Result<f64, NonlinearWindingError>;

        /// The field's materials, one per mesh triangle, when the owner
        /// refreshes them with temperature; none for the play-memory owner.
        fn magnetic_materials(&self) -> Option<&[IsotropicMagneticMaterial]>;

        /// Each mesh triangle's flux density \[T\] and volume \[m³\],
        /// recovered from the accepted field's potential.
        fn flux_densities(&self) -> Result<Vec<([f64; 2], f64)>, NonlinearWindingError>;

        /// A step's magnetic cell records, in mesh order; none for the
        /// play-memory owner.
        fn magnetic_cells(step: &Self::FieldStep) -> &[NonlinearEddyMagneticCell];

        /// Whether two states' problems hold the same mesh, thickness and
        /// materials. A step that refreshes the materials builds its own
        /// problem.
        fn same_problem(&self, other: &Self) -> bool;

        fn graph(&self) -> &CircuitGraph;
        fn time_s(&self) -> f64;
        fn state(&self) -> &[f64];
        fn currents_a(&self) -> &[f64];
        fn flux_linkages_wb(&self) -> &[f64];
        fn resistance_ohm(&self) -> &[f64];
        fn conductors(&self) -> &[WindingConductor];
        /// The field's core conductors.
        fn core_conductors(&self) -> &[EddyConductor];
        /// The address of the field's problem, which equal states share.
        fn problem(&self) -> *const ();
        fn potential_wb_m(&self) -> &[f64];
        fn field_energy_j(&self) -> f64;
        fn field_time_s(&self) -> f64;
    }

    impl<'m> Sealed<'m> for NonlinearWindingCircuit<'m> {
        type FieldControls = NonlinearEddyControls;
        type FieldStep = NonlinearEddyStep;
        const LANE_WORDS: &'static [u64] = &[];

        fn step(
            &self,
            target_s: f64,
            materials: Option<Refresh>,
            core: Vec<EddyConductor>,
            windings: Vec<WindingConductor>,
            inputs: &[f64],
            controls: NonlinearWindingControls,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<(Self, NonlinearWindingStep), NonlinearWindingError> {
            let candidate = match materials {
                Some((materials, assembly)) => self.step_with_materials_candidate(
                    target_s, materials, core, windings, inputs, assembly, controls, cancelled,
                )?,
                None => {
                    self.step_to_candidate(target_s, core, windings, inputs, controls, cancelled)?
                }
            };
            Ok(candidate.into_parts())
        }

        fn rebind_graph_resistors(
            &self,
            sources: Vec<GraphResistorMaterial>,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<Self, NonlinearWindingError> {
            self.with_graph_resistors(sources, cancelled)
        }

        fn field_heat(step: &NonlinearEddyStep) -> Vec<(usize, f64)> {
            step.cells
                .iter()
                .map(|cell| (cell.source.triangle, cell.heat_j))
                .collect()
        }

        fn field_triangles(&self) -> Range<usize> {
            0..0
        }

        fn admit_field_temperature(&self, _: usize, _: f64) -> Result<(), String> {
            Ok(())
        }

        fn same_memory(&self, _: &Self) -> bool {
            true
        }

        fn stored_energy_j(&self) -> Result<f64, NonlinearWindingError> {
            let end = self
                .end_inductance()
                .stored_energy_j(self.currents_a())
                .map_err(NonlinearWindingError::from)?;
            Ok(self.dae().system.hamiltonian(self.state()) + self.field().energy_j() + end)
        }

        fn magnetic_materials(&self) -> Option<&[IsotropicMagneticMaterial]> {
            Some(self.field().problem().materials())
        }

        fn flux_densities(&self) -> Result<Vec<([f64; 2], f64)>, NonlinearWindingError> {
            let field = self.field().field(&mut || false)?;
            let mesh = self.field().problem().mesh();
            field
                .magnetic_flux_density_t()
                .iter()
                .enumerate()
                .map(|(triangle, &b)| {
                    mesh.face_measure(triangle)
                        .map(|volume| (b, volume))
                        .ok_or(NonlinearWindingError::Invalid("missing magnetic volume"))
                })
                .collect()
        }

        fn magnetic_cells(step: &NonlinearEddyStep) -> &[NonlinearEddyMagneticCell] {
            &step.magnetic_cells
        }

        fn same_problem(&self, other: &Self) -> bool {
            let (a, b) = (self.field().problem(), other.field().problem());
            let (x, y) = (a.mesh(), b.mesh());
            core::ptr::eq(a, b)
                || (a.materials() == b.materials()
                    && a.thickness_m().to_bits() == b.thickness_m().to_bits()
                    && x.lineage() == y.lineage()
                    && x.metric() == y.metric()
                    && x.faces() == y.faces()
                    && x.vertex_keys() == y.vertex_keys()
                    && x.vertices().len() == y.vertices().len()
                    && x.vertices()
                        .iter()
                        .zip(y.vertices())
                        .all(|(p, q)| super::same_bits(p, q)))
        }

        fn graph(&self) -> &CircuitGraph {
            Self::graph(self)
        }

        fn time_s(&self) -> f64 {
            Self::time_s(self)
        }

        fn state(&self) -> &[f64] {
            Self::state(self)
        }

        fn currents_a(&self) -> &[f64] {
            Self::currents_a(self)
        }

        fn flux_linkages_wb(&self) -> &[f64] {
            Self::flux_linkages_wb(self)
        }

        fn resistance_ohm(&self) -> &[f64] {
            Self::resistance_ohm(self)
        }

        fn conductors(&self) -> &[WindingConductor] {
            Self::conductors(self)
        }

        fn core_conductors(&self) -> &[EddyConductor] {
            self.field().conductors()
        }

        fn problem(&self) -> *const () {
            core::ptr::from_ref(self.field().problem()).cast()
        }

        fn potential_wb_m(&self) -> &[f64] {
            self.field().potential_wb_m()
        }

        fn field_energy_j(&self) -> f64 {
            self.field().energy_j()
        }

        fn field_time_s(&self) -> f64 {
            self.field().time_s()
        }
    }

    /// The play-memory owner: its field heats each eddy cell and, by pinning
    /// and viscous switching, every magnetic triangle, conducting or not.
    #[cfg(feature = "magnetic-hysteresis")]
    impl<'m> Sealed<'m> for HystereticWindingCircuit<'m> {
        type FieldControls = HystereticEddyControls;
        type FieldStep = HystereticEddyStep;
        const LANE_WORDS: &'static [u64] = &[u64::from_le_bytes(*b"vecplay1")];

        fn step(
            &self,
            target_s: f64,
            materials: Option<Refresh>,
            core: Vec<EddyConductor>,
            windings: Vec<WindingConductor>,
            inputs: &[f64],
            controls: NonlinearWindingControls<HystereticEddyControls>,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<(Self, NonlinearWindingStep<HystereticEddyStep>), NonlinearWindingError>
        {
            if materials.is_some() {
                return Err(NonlinearWindingError::Invalid(
                    "the play-memory field has no temperature refresh",
                ));
            }
            Ok(self
                .step_to_candidate(target_s, core, windings, inputs, controls, cancelled)?
                .into_parts())
        }

        fn rebind_graph_resistors(
            &self,
            sources: Vec<GraphResistorMaterial>,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<Self, NonlinearWindingError> {
            self.with_graph_resistors(sources, cancelled)
        }

        /// Eddy cells, then each magnetic cell's pinning plus viscous heat
        /// (`irreversible_heat_j`), as the thermal module deposits them.
        fn field_heat(step: &HystereticEddyStep) -> Vec<(usize, f64)> {
            step.cells
                .iter()
                .map(|cell| (cell.source.triangle, cell.heat_j))
                .chain(
                    step.magnetic
                        .cells
                        .iter()
                        .map(|cell| (cell.triangle, cell.irreversible_heat_j)),
                )
                .collect()
        }

        fn field_triangles(&self) -> Range<usize> {
            0..self.field().problem().mesh().faces().len()
        }

        /// The fit's coefficients are constant, so they must hold over the
        /// whole interval between the temperature the fit was resolved at
        /// and `kelvin`. At the fit's own temperature they are the fit's.
        fn admit_field_temperature(&self, triangle: usize, kelvin: f64) -> Result<(), String> {
            let law = self
                .field()
                .problem()
                .materials()
                .get(triangle)
                .ok_or_else(|| "the field has no magnetic fit in this triangle".to_owned())?;
            let (axis, fitted) = fitted_temperature(law.material()).ok_or_else(|| {
                "the magnetic fit was resolved at no single absolute temperature".to_owned()
            })?;
            if kelvin.to_bits() == fitted.to_bits() {
                return Ok(());
            }
            law.admit_temperature_interval(axis, [kelvin.min(fitted), kelvin.max(fitted)])
                .map_err(|error| error.to_string())
        }

        /// Every cell's field, flux density, branch polarizations and
        /// stored energy density, the curl of the accepted field and its
        /// stored energy, bit for bit.
        fn same_memory(&self, other: &Self) -> bool {
            let (a, b) = (
                self.field().magnetic_state(),
                other.field().magnetic_state(),
            );
            a.stored_energy_j().to_bits() == b.stored_energy_j().to_bits()
                && same_vectors(a.magnetic_flux_density_t(), b.magnetic_flux_density_t())
                && a.cell_states().len() == b.cell_states().len()
                && a.cell_states().iter().zip(b.cell_states()).all(|(x, y)| {
                    x.stored_energy_j_m3().to_bits() == y.stored_energy_j_m3().to_bits()
                        && same_vectors(
                            &[x.field_a_m(), x.flux_density_t()],
                            &[y.field_a_m(), y.flux_density_t()],
                        )
                        && same_vectors(x.polarizations_t(), y.polarizations_t())
                })
        }

        fn stored_energy_j(&self) -> Result<f64, NonlinearWindingError> {
            let end = self
                .end_inductance()
                .stored_energy_j(self.currents_a())
                .map_err(NonlinearWindingError::from)?;
            Ok(self.dae().system.hamiltonian(self.state()) + self.field().energy_j() + end)
        }

        fn magnetic_materials(&self) -> Option<&[IsotropicMagneticMaterial]> {
            None
        }

        fn flux_densities(&self) -> Result<Vec<([f64; 2], f64)>, NonlinearWindingError> {
            Err(NonlinearWindingError::Invalid(
                "the play-memory field has no temperature refresh",
            ))
        }

        fn magnetic_cells(_: &HystereticEddyStep) -> &[NonlinearEddyMagneticCell] {
            &[]
        }

        /// Never asked: no state of the play-memory owner carries
        /// temperatures.
        fn same_problem(&self, other: &Self) -> bool {
            core::ptr::eq(self.problem(), other.problem())
        }

        fn graph(&self) -> &CircuitGraph {
            Self::graph(self)
        }

        fn time_s(&self) -> f64 {
            Self::time_s(self)
        }

        fn state(&self) -> &[f64] {
            Self::state(self)
        }

        fn currents_a(&self) -> &[f64] {
            Self::currents_a(self)
        }

        fn flux_linkages_wb(&self) -> &[f64] {
            Self::flux_linkages_wb(self)
        }

        fn resistance_ohm(&self) -> &[f64] {
            Self::resistance_ohm(self)
        }

        fn conductors(&self) -> &[WindingConductor] {
            Self::conductors(self)
        }

        fn core_conductors(&self) -> &[EddyConductor] {
            self.field().conductors()
        }

        fn problem(&self) -> *const () {
            core::ptr::from_ref(self.field().problem()).cast()
        }

        fn potential_wb_m(&self) -> &[f64] {
            self.field().potential_wb_m()
        }

        fn field_energy_j(&self) -> f64 {
            self.field().energy_j()
        }

        fn field_time_s(&self) -> f64 {
            self.field().time_s()
        }
    }
}

/// The accepted state: the owner's circuit, which carries the field history;
/// in the magnetic lane, the temperatures \[K\] its materials were evaluated
/// at, one per heat receiver; and, for a state a step proposed, the owner's
/// step that produced it. The step is held as `dyn Any`, not through the
/// owner's associated type, so the state stays covariant in `'m`; a clone
/// shares it.
#[derive(Clone)]
pub struct PlanarWindingState<'m, C = NonlinearWindingCircuit<'m>>(
    C,
    Option<Vec<f64>>,
    PhantomData<&'m ()>,
    Option<Rc<dyn Any>>,
);

impl<'m, C: WindingOwner<'m>> PlanarWindingState<'m, C> {
    /// The state of an accepted circuit.
    #[must_use]
    pub const fn new(circuit: C) -> Self {
        Self(circuit, None, PhantomData, None)
    }

    /// The magnetic lane's state of an accepted circuit whose field
    /// materials were evaluated at `temperatures` \[K\], one per heat
    /// receiver. A step refuses a state whose temperatures do not fit its
    /// Component, and one whose law triangles' free energy at them is not
    /// the field's.
    #[must_use]
    pub const fn at_temperatures(circuit: C, temperatures: Vec<f64>) -> Self {
        Self(circuit, Some(temperatures), PhantomData, None)
    }

    /// The temperatures \[K\] the field materials were evaluated at, in the
    /// magnetic lane.
    #[must_use]
    pub fn temperatures(&self) -> Option<&[f64]> {
        self.1.as_deref()
    }

    /// The accepted circuit.
    #[must_use]
    pub const fn circuit(&self) -> &C {
        &self.0
    }

    /// The accepted circuit, by value.
    #[must_use]
    pub fn into_circuit(self) -> C {
        self.0
    }
}

impl PlanarWindingState<'_> {
    /// The saturating owner's step that produced this state: its cells, graph
    /// resistors, field record, energies and solve report, as the owner's own
    /// `step_to_candidate` gave them. None for a constructed state.
    #[must_use]
    pub fn step(&self) -> Option<&NonlinearWindingStep> {
        self.3.as_deref()?.downcast_ref()
    }
}

#[cfg(feature = "magnetic-hysteresis")]
impl<'m> PlanarWindingState<'m, HystereticWindingCircuit<'m>> {
    /// The play-memory owner's step that produced this state, as the owner's
    /// own `step_to_candidate` gave it. None for a constructed state.
    #[must_use]
    pub fn step(&self) -> Option<&HystereticWindingStep> {
        self.3.as_deref()?.downcast_ref()
    }
}

/// Equal states share the owner's problem (the same allocation) and graph,
/// and hold the same bits in every accepted coordinate: clock, circuit
/// state, currents, linkages, resistances, conductors and field potential,
/// energy and clock, and the owner's field memory (for the play-memory
/// owner, every cell's vector history). The graph is compared by value: a
/// step with resistor laws builds a new one. In the magnetic lane, each step
/// builds its own problem too, so there the problems are compared by mesh,
/// thickness and materials, and the temperatures by their bits; a state
/// with temperatures never equals one without. The owner's step a state
/// keeps is not compared: a proposed state equals the constructed one.
impl<'m, C: WindingOwner<'m>> PartialEq for PlanarWindingState<'m, C> {
    fn eq(&self, other: &Self) -> bool {
        let (a, b) = (&self.0, &other.0);
        let problem = match (&self.1, &other.1) {
            (None, None) => core::ptr::eq(a.problem(), b.problem()),
            (Some(x), Some(y)) => same_bits(x, y) && a.same_problem(b),
            _ => false,
        };
        shared::same_graph(a.graph(), b.graph(), false)
            && problem
            && a.time_s().to_bits() == b.time_s().to_bits()
            && same_bits(a.state(), b.state())
            && same_bits(a.currents_a(), b.currents_a())
            && same_bits(a.flux_linkages_wb(), b.flux_linkages_wb())
            && same_bits(a.resistance_ohm(), b.resistance_ohm())
            && a.conductors() == b.conductors()
            && a.core_conductors() == b.core_conductors()
            && same_bits(a.potential_wb_m(), b.potential_wb_m())
            && a.field_energy_j().to_bits() == b.field_energy_j().to_bits()
            && a.field_time_s().to_bits() == b.field_time_s().to_bits()
            && a.same_memory(b)
    }
}

impl<'m, C: WindingOwner<'m>> fmt::Debug for PlanarWindingState<'m, C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PlanarWindingState")
            .field("time_s", &self.0.time_s())
            .field("currents_a", &self.0.currents_a())
            .field("state", &self.0.state())
            .field("temperatures_k", &self.1)
            .finish_non_exhaustive()
    }
}

/// A step's circuit, when resistor laws rebuilt it, and its core and
/// stranded conductors.
type Supplied<C> = (Option<C>, Vec<EddyConductor>, Vec<WindingConductor>);

/// Where a step's conductors and graph resistances come from.
enum Supply {
    /// Given at construction and used at every step.
    Fixed {
        core: Vec<EddyConductor>,
        windings: Vec<WindingConductor>,
    },
    /// Laws asked at the temperatures the heat ports receive.
    Laws {
        /// Every conductor triangle, with its core region or `None` when
        /// stranded, as the construction's circuit holds them.
        cells: BTreeMap<usize, Option<u32>>,
        /// One per heat receiver.
        conductors: Vec<Arc<dyn ConductorLaw>>,
        /// One per native graph resistor, by branch, or none: the
        /// circuit's own resistors stay.
        resistors: BTreeMap<usize, Arc<dyn ResistorLaw>>,
    },
}

/// The planar winding owner behind the composition contract's source and
/// heat ports. See the module documentation. `C` is the owner's circuit and
/// `F` its field's controls; both default to the saturating owner. Neither
/// is reached through an associated type in a field, so the Component stays
/// covariant in `'m`, as the saturating Component was.
pub struct PlanarWinding<'m, C = NonlinearWindingCircuit<'m>, F = NonlinearEddyControls> {
    graph: CircuitGraph,
    sources: usize,
    supply: Supply,
    receivers: Vec<HeatReceiver>,
    triangle_port: BTreeMap<usize, usize>,
    resistor_port: BTreeMap<usize, usize>,
    /// The owner's field triangles that may dissipate beyond the conductors.
    field: Range<usize>,
    controls: NonlinearWindingControls<F>,
    energy_tolerance_j: f64,
    ports: Vec<PortSpec>,
    /// The magnetic lane's laws, by triangle, and its assembly bounds.
    magnetic: Option<(
        BTreeMap<usize, ThermalMagneticMaterial>,
        WeightedAssemblyLimits,
    )>,
    problem: PhantomData<&'m ()>,
    owner: PhantomData<fn() -> C>,
}

impl<C, F> fmt::Debug for PlanarWinding<'_, C, F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PlanarWinding")
            .field("sources", &self.sources)
            .field("receivers", &self.receivers)
            .field("energy_tolerance_j", &self.energy_tolerance_j)
            .finish_non_exhaustive()
    }
}

impl<'m, C, F> PlanarWinding<'m, C, F>
where
    C: WindingOwner<'m, FieldControls = F>,
    F: Copy + fmt::Debug,
{
    /// The Component of `circuit`'s topology, stepped with these conductors
    /// and controls, its heat leaving through `receivers`.
    ///
    /// # Errors
    /// The energy tolerance is not positive and finite; an external port
    /// drives no source or two; a receiver is empty; a stranded or core
    /// triangle, a triangle of a field that heats beyond its conductors, or a
    /// graph resistor is mapped to no receiver or to two; a receiver names an
    /// element that does not dissipate.
    pub fn new(
        circuit: &C,
        core_conductors: Vec<EddyConductor>,
        winding_conductors: Vec<WindingConductor>,
        receivers: Vec<HeatReceiver>,
        controls: NonlinearWindingControls<F>,
        energy_tolerance_j: f64,
    ) -> Result<Self, PlanarWindingError> {
        let triangles: Vec<usize> = core_conductors
            .iter()
            .map(|c| c.triangle)
            .chain(winding_conductors.iter().map(|c| c.triangle))
            .chain(circuit.field_triangles())
            .collect();
        let supply = Supply::Fixed {
            core: core_conductors,
            windings: winding_conductors,
        };
        Self::build(
            circuit,
            &triangles,
            supply,
            receivers,
            controls,
            energy_tolerance_j,
        )
    }

    /// The Component of `circuit`'s topology whose conductors and graph
    /// resistors follow the temperatures its heat ports receive: receiver
    /// `k`'s law gives the conductor of each of its triangles at port `k`'s
    /// temperature, and each graph resistor's law gives it at the
    /// temperature of the port its heat leaves through. With no resistor
    /// laws, the graph resistors keep `circuit`'s values. The core and
    /// stranded triangles, and the core regions, are `circuit`'s own.
    ///
    /// # Errors
    /// As [`Self::new`]; the number of conductor laws is not the number of
    /// receivers; some graph resistors have laws and others do not, or a law
    /// names a branch that is not a resistor.
    pub fn at_port_temperatures(
        circuit: &C,
        conductor_laws: Vec<Arc<dyn ConductorLaw>>,
        resistor_laws: BTreeMap<usize, Arc<dyn ResistorLaw>>,
        receivers: Vec<HeatReceiver>,
        controls: NonlinearWindingControls<F>,
        energy_tolerance_j: f64,
    ) -> Result<Self, PlanarWindingError> {
        let (supply, triangles) =
            Self::laws(circuit, conductor_laws, resistor_laws, receivers.len())?;
        Self::build(
            circuit,
            &triangles,
            supply,
            receivers,
            controls,
            energy_tolerance_j,
        )
    }

    /// [`Self::at_port_temperatures`] in the magnetic lane (module
    /// documentation): the permeability of each triangle with a
    /// [`MagneticLaw`] also follows the temperature its receiver's port
    /// gives. Each law triangle maps to a receiver, conducting or not. The
    /// Component steps states made by [`PlanarWindingState::at_temperatures`].
    ///
    /// # Errors
    /// As [`Self::at_port_temperatures`]; the owner refreshes no field
    /// materials with temperature (the play-memory owner); there is no law;
    /// two laws name one triangle; a law names a triangle whose field
    /// material is not linear, or whose material is not the law's source
    /// (its card, claims and other coordinates); a law triangle maps to no
    /// receiver.
    pub fn magnetic_at_port_temperatures(
        circuit: &C,
        conductor_laws: Vec<Arc<dyn ConductorLaw>>,
        resistor_laws: BTreeMap<usize, Arc<dyn ResistorLaw>>,
        magnetic: MagneticLaws,
        receivers: Vec<HeatReceiver>,
        controls: NonlinearWindingControls<F>,
        energy_tolerance_j: f64,
    ) -> Result<Self, PlanarWindingError> {
        let materials = circuit
            .magnetic_materials()
            .ok_or_else(|| invalid("this owner refreshes no field materials with temperature"))?;
        let mut laws = BTreeMap::new();
        for MagneticLaw { triangle, law } in magnetic.triangles {
            let Some(IsotropicMagneticMaterial::Linear(material)) = materials.get(triangle) else {
                return Err(invalid(
                    "a magnetic law names a triangle whose field material is not linear",
                ));
            };
            if admit_magnetic_temperature_refresh(material, law.material()).is_err() {
                return Err(invalid(
                    "a magnetic law's source is not its triangle's field material",
                ));
            }
            if laws.insert(triangle, law).is_some() {
                return Err(invalid("two magnetic laws name one triangle"));
            }
        }
        if laws.is_empty() {
            return Err(invalid("the magnetic lane has no law"));
        }
        let (supply, mut triangles) =
            Self::laws(circuit, conductor_laws, resistor_laws, receivers.len())?;
        triangles.extend(laws.keys().copied());
        let mut component = Self::build(
            circuit,
            &triangles,
            supply,
            receivers,
            controls,
            energy_tolerance_j,
        )?;
        component.magnetic = Some((laws, magnetic.assembly));
        Ok(component)
    }

    /// The supply of [`Self::at_port_temperatures`], and the triangles its
    /// receivers must map.
    fn laws(
        circuit: &C,
        conductor_laws: Vec<Arc<dyn ConductorLaw>>,
        resistor_laws: BTreeMap<usize, Arc<dyn ResistorLaw>>,
        receivers: usize,
    ) -> Result<(Supply, Vec<usize>), PlanarWindingError> {
        if conductor_laws.len() != receivers {
            return Err(invalid("one conductor law per heat receiver required"));
        }
        let cells: BTreeMap<usize, Option<u32>> = circuit
            .core_conductors()
            .iter()
            .map(|c| (c.triangle, Some(c.region)))
            .chain(circuit.conductors().iter().map(|c| (c.triangle, None)))
            .collect();
        let triangles: Vec<usize> = cells
            .keys()
            .copied()
            .chain(circuit.field_triangles())
            .collect();
        let supply = Supply::Laws {
            cells,
            conductors: conductor_laws,
            resistors: resistor_laws,
        };
        Ok((supply, triangles))
    }

    fn build(
        circuit: &C,
        triangles: &[usize],
        supply: Supply,
        receivers: Vec<HeatReceiver>,
        controls: NonlinearWindingControls<F>,
        energy_tolerance_j: f64,
    ) -> Result<Self, PlanarWindingError> {
        let laws = match &supply {
            Supply::Laws { resistors, .. } => Some(resistors),
            Supply::Fixed { .. } => None,
        };
        let layout = shared::layout(
            circuit.graph(),
            triangles,
            &maps(&receivers),
            laws,
            energy_tolerance_j,
        )
        .map_err(refused)?;
        Ok(Self {
            graph: circuit.graph().clone(),
            sources: layout.sources,
            supply,
            receivers,
            triangle_port: layout.element_port,
            resistor_port: layout.resistor_port,
            field: circuit.field_triangles(),
            controls,
            energy_tolerance_j,
            ports: layout.ports,
            magnetic: None,
            problem: PhantomData,
            owner: PhantomData,
        })
    }

    /// The owner's clock at the end of a step of `dt_s` from `old`, when the
    /// step, its ports and the state's topology are this Component's.
    fn admit(
        &self,
        old: &PlanarWindingState<'m, C>,
        dt_s: f64,
        inputs: &[f64],
        outputs: usize,
    ) -> Result<f64, PlanarWindingError> {
        shared::admit_ports(dt_s, inputs, outputs, self.sources, self.receivers.len())
            .map_err(refused)?;
        let laws = matches!(&self.supply, Supply::Laws { resistors, .. } if !resistors.is_empty());
        if !shared::same_graph(old.0.graph(), &self.graph, laws)
            || old.0.field_triangles() != self.field
        {
            return Err(invalid(
                "the state's circuit is not this Component's topology",
            ));
        }
        let lane = self.magnetic.as_ref().map(|_| self.receivers.len());
        shared::admit_temperatures(old.1.as_deref(), lane).map_err(refused)?;
        shared::target_time(old.0.time_s(), dt_s).map_err(refused)
    }

    /// The circuit to step, when resistor laws rebuild it, and the core and
    /// stranded conductors of a step whose heat ports received `kelvin`.
    fn supply_at(
        &self,
        old: &C,
        kelvin: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Supplied<C>, PlanarWindingError> {
        let (cells, conductors, resistors) = match &self.supply {
            Supply::Fixed { core, windings } => return Ok((None, core.clone(), windings.clone())),
            Supply::Laws {
                cells,
                conductors,
                resistors,
            } => (cells, conductors, resistors),
        };
        let mut core = Vec::new();
        let mut windings = Vec::new();
        for (&triangle, &region) in cells {
            let receiver = self.triangle_port[&triangle];
            let kelvin = kelvin[receiver];
            let conductor =
                conductors[receiver]
                    .conductor(triangle, kelvin)
                    .map_err(|message| PlanarWindingError::Law {
                        receiver,
                        element: LawElement::Triangle(triangle),
                        kelvin,
                        message,
                    })?;
            match region {
                Some(region) => core.push(EddyConductor {
                    triangle,
                    region,
                    conductor,
                }),
                None => windings.push(WindingConductor {
                    triangle,
                    conductor,
                }),
            }
        }
        if resistors.is_empty() {
            return Ok((None, core, windings));
        }
        let sources =
            shared::graph_resistors(resistors, &self.resistor_port, kelvin).map_err(refused)?;
        let driven = old.rebind_graph_resistors(sources, cancelled)?;
        Ok((Some(driven), core, windings))
    }

    /// The heat of each receiver over a step [J], from the owner's records
    /// of its stranded cells, field cells (`field`) and graph resistors, in
    /// that order.
    fn heat_by_receiver(
        &self,
        step: &NonlinearWindingStep<C::FieldStep>,
        field: &[(usize, f64)],
    ) -> Result<Vec<f64>, PlanarWindingError> {
        shared::heat_by_receiver(
            self.receivers.len(),
            &self.triangle_port,
            &self.resistor_port,
            step.cells
                .iter()
                .map(|cell| (cell.source.triangle, cell.heat_j))
                .chain(field.iter().copied()),
            step.graph_resistors
                .iter()
                .map(|resistor| (resistor.branch, resistor.heat_j)),
        )
        .map_err(refused)
    }

    /// The declared heat receivers, in heat-port order.
    #[must_use]
    pub fn receivers(&self) -> &[HeatReceiver] {
        &self.receivers
    }

    /// In the magnetic lane, the field materials of a step whose heat ports
    /// received `kelvin`: each law triangle's law there, and every other
    /// triangle's material as `old` holds it.
    fn materials_at(&self, old: &C, kelvin: &[f64]) -> Result<Option<Refresh>, PlanarWindingError> {
        let Some((laws, assembly)) = &self.magnetic else {
            return Ok(None);
        };
        let mut materials = old
            .magnetic_materials()
            .ok_or_else(|| invalid("this owner refreshes no field materials with temperature"))?
            .to_vec();
        for (&triangle, law) in laws {
            let receiver = self.triangle_port[&triangle];
            let kelvin = kelvin[receiver];
            let at = law
                .at_temperature(kelvin)
                .map_err(|error| PlanarWindingError::Law {
                    receiver,
                    element: LawElement::Permeability(triangle),
                    kelvin,
                    message: error.to_string(),
                })?;
            *materials
                .get_mut(triangle)
                .ok_or_else(|| invalid("the state's field has no such law triangle"))? =
                IsotropicMagneticMaterial::Linear(at.material().clone());
        }
        Ok(Some((materials, *assembly)))
    }

    /// The magnetic lane's terms of a step from `old`'s temperatures to the
    /// ports' (`kelvin`), read from the owner's magnetic cell records. Each
    /// law's free energy must be the field's at both ends of the step, within
    /// the energy tolerance; a triangle without a law must do no parameter
    /// work.
    fn reversible(
        &self,
        old: &PlanarWindingState<'m, C>,
        cells: &[NonlinearEddyMagneticCell],
        kelvin: &[f64],
    ) -> Result<Option<Reversible>, PlanarWindingError> {
        let Some((laws, _)) = &self.magnetic else {
            return Ok(None);
        };
        let before = old
            .1
            .as_deref()
            .ok_or_else(|| invalid("a state of the magnetic lane carries its temperatures"))?;
        let cells: Vec<LawCell> = cells
            .iter()
            .map(|cell| LawCell {
                element: cell.triangle,
                flux_density_t: cell.flux_density_t.map(|b| b[0].hypot(b[1])),
                free_energy_j: cell.free_energy_j,
                volume_m3: cell.volume_m3,
                parameter_work_j: cell.parameter_work_j,
            })
            .collect();
        shared::reversible(
            laws,
            &self.triangle_port,
            self.receivers.len(),
            before,
            kelvin,
            &cells,
            self.energy_tolerance_j,
        )
        .map(Some)
        .map_err(refused)
    }

    /// In the magnetic lane, a state's `(sum of T s [J], sum of s [J/K])`
    /// over the law triangles, measured from the state alone: each law at
    /// its receiver's temperature, at the flux density the field's potential
    /// gives. Both sums are compensated (`two_sum`).
    fn stored_reversible(
        &self,
        state: &PlanarWindingState<'m, C>,
    ) -> Result<Option<(f64, f64)>, PlanarWindingError> {
        let Some((laws, _)) = &self.magnetic else {
            return Ok(None);
        };
        let temperatures = state
            .1
            .as_deref()
            .filter(|t| t.len() == self.receivers.len())
            .ok_or_else(|| invalid("a state of the magnetic lane carries its temperatures"))?;
        let cells: Vec<(f64, f64)> = state
            .0
            .flux_densities()?
            .into_iter()
            .map(|(b, volume_m3)| (b[0].hypot(b[1]), volume_m3))
            .collect();
        shared::stored_reversible(laws, &self.triangle_port, temperatures, &cells)
            .map(Some)
            .map_err(refused)
    }
}

/// Circuit, field and end-mode storage of an accepted state, summed as the
/// owner sums `stored_before_j`.
fn energy<'m, C: WindingOwner<'m>>(
    state: &PlanarWindingState<'m, C>,
) -> Result<f64, PlanarWindingError> {
    let total = state.0.stored_energy_j()?;
    if !total.is_finite() {
        return Err(invalid("the winding's stored energy is not finite"));
    }
    Ok(total)
}

impl<'m, C, F> Component for PlanarWinding<'m, C, F>
where
    C: WindingOwner<'m, FieldControls = F>,
    F: Copy + fmt::Debug,
{
    type State = PlanarWindingState<'m, C>;
    type Error = PlanarWindingError;

    /// The graph; the fixed conductors' triangles, regions, resistivity bits
    /// and card identities, or the laws' cells and words; the receiver map;
    /// the controls by their `Debug` text (which prints each `f64` so that it
    /// parses back to the same bits); the energy tolerance; then the
    /// owner's lane words (none for the saturating owner, `b"vecplay1"` read
    /// little-endian for the play-memory owner); then, in the magnetic lane
    /// only, `b"magtemp1"`, each law's triangle and words (its card identity,
    /// temperature axis, range and slope, and its source's permeability), and
    /// the assembly bounds. Each list is preceded by its length.
    fn identity_words(&self, words: &mut Vec<u64>) {
        shared::graph_words(&self.graph, words);
        match &self.supply {
            Supply::Fixed { core, windings } => {
                words.extend([0, core.len() as u64]);
                for core in core {
                    words.extend([core.triangle as u64, u64::from(core.region)]);
                    shared::conductor_words(&core.conductor, words);
                }
                words.push(windings.len() as u64);
                for winding in windings {
                    words.push(winding.triangle as u64);
                    shared::conductor_words(&winding.conductor, words);
                }
            }
            Supply::Laws {
                cells,
                conductors,
                resistors,
            } => {
                words.extend([1, cells.len() as u64]);
                for (&triangle, &region) in cells {
                    words.extend([triangle as u64, region.map_or(u64::MAX, u64::from)]);
                }
                words.push(conductors.len() as u64);
                for law in conductors {
                    shared::law_words(words, |words| law.identity_words(words));
                }
                shared::resistor_law_words(resistors, words);
            }
        }
        shared::receiver_words(&maps(&self.receivers), words);
        shared::debug_words(&self.controls, words);
        words.push(self.energy_tolerance_j.to_bits());
        words.extend_from_slice(C::LANE_WORDS);
        if let Some((laws, assembly)) = &self.magnetic {
            shared::lane_words(laws, words);
            words.extend([
                assembly.max_cells as u64,
                assembly.max_dofs as u64,
                assembly.max_triplets as u64,
            ]);
        }
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn end_of_step(&self, _port: usize) -> EndOfStep {
        // Backward Euler with held inputs: no flow ends its step by the
        // midpoint relation, and none is told.
        EndOfStep::Untold
    }

    /// With the magnetic lane's `T s` of each law triangle, measured from
    /// the state.
    fn stored_energy_j(
        &self,
        state: &PlanarWindingState<'m, C>,
    ) -> Result<f64, PlanarWindingError> {
        let total = energy(state)?;
        match self.stored_reversible(state)? {
            None => Ok(total),
            Some((reversible, _)) => Ok(total + reversible),
        }
    }

    /// In the magnetic lane, the sum of `s = -dpsi/dT` over the law
    /// triangles, measured from the state's own field and temperatures;
    /// otherwise none.
    fn stored_entropy_j_per_k(
        &self,
        state: &PlanarWindingState<'m, C>,
    ) -> Result<Option<f64>, PlanarWindingError> {
        Ok(self.stored_reversible(state)?.map(|(_, entropy)| entropy))
    }

    fn propose(
        &self,
        old: &PlanarWindingState<'m, C>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut PlanarWindingState<'m, C>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, PlanarWindingError> {
        cx.checkpoint().map_err(|_| PlanarWindingError::Cancelled)?;
        let target_s = self.admit(old, dt_s, inputs, outputs.len())?;
        let sources = self.sources;
        for triangle in self.field.clone() {
            cx.checkpoint().map_err(|_| PlanarWindingError::Cancelled)?;
            let receiver = self.triangle_port[&triangle];
            let kelvin = inputs[sources + receiver];
            old.0
                .admit_field_temperature(triangle, kelvin)
                .map_err(|message| PlanarWindingError::Law {
                    receiver,
                    element: LawElement::Field(triangle),
                    kelvin,
                    message,
                })?;
        }
        let (driven, core, windings) =
            self.supply_at(&old.0, &inputs[sources..], &mut || cx.checkpoint().is_err())?;
        let materials = self.materials_at(&old.0, &inputs[sources..])?;
        let lane = materials.is_some();
        let (circuit, step) = driven.as_ref().unwrap_or(&old.0).step(
            target_s,
            materials,
            core,
            windings,
            &inputs[..sources],
            self.controls,
            &mut || cx.checkpoint().is_err(),
        )?;
        let kept = Rc::new(step);
        let step = &*kept;
        shared::admit_work(step.parameter_work_j, step.mechanical_work_j, lane).map_err(refused)?;
        let field = C::field_heat(&step.field);
        let heat = self.heat_by_receiver(step, &field)?;
        let kelvin = &inputs[sources..];
        let reversible = self.reversible(old, C::magnetic_cells(&step.field), kelvin)?;
        let account = Account {
            stored_before_j: step.stored_before_j,
            stored_after_j: step.stored_after_j,
            stored_change_j: step.stored_change_j,
            supplied_work_j: step.supplied_work_j,
            numerical_dissipation_j: step.numerical_dissipation_j,
            elements: step.cells.len() + field.len() + step.graph_resistors.len(),
            coordinates: [old.0.state().len(), old.0.potential_wb_m().len()],
        };
        let laws = self.magnetic.as_ref().map_or(0, |(laws, _)| laws.len());
        let closed = shared::close(
            &account,
            &heat,
            kelvin,
            reversible.as_ref().map(|terms| (terms, laws)),
            sources,
            self.energy_tolerance_j,
        )
        .map_err(refused)?;
        cx.checkpoint().map_err(|_| PlanarWindingError::Cancelled)?;
        let receipt = shared::publish(
            &closed,
            &step.report,
            &self.controls.newton,
            &step.native.y,
            sources,
            dt_s,
            outputs,
        );
        let mut state = if lane {
            PlanarWindingState::at_temperatures(circuit, kelvin.to_vec())
        } else {
            PlanarWindingState::new(circuit)
        };
        state.3 = Some(kept);
        *new = state;
        Ok(Proposal::Step(receipt))
    }
}
