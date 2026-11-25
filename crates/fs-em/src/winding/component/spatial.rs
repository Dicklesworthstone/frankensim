//! fs-em's 3D winding circuits as one Component of the composition contract
//! (fs-ports).
//!
//! [`SpatialWinding`] steps its owner's circuit ([`SpatialWindingOwner`]; by
//! default [`SpatialNonlinearWindingCircuit`], the saturating field) once per
//! proposal, from the accepted state, with conductors fixed at construction
//! ([`SpatialWinding::new`]) or given by laws at the temperatures its heat
//! ports receive ([`SpatialWinding::at_port_temperatures`]). Ports, in order:
//! `source<p>` for each original external input `p` of the native graph
//! (`EffortIn` for a voltage source, `FlowIn` for a current source),
//! returning the native step-average conjugate `StepRecord::y[p]`; then
//! `heat<k>` for each declared [`SpatialHeatReceiver`], receiving temperature
//! and returning minus the mean heat rate of the stranded cells, core cells
//! and native graph resistors mapped to it. A mesh cell is a tetrahedron.
//! Every dissipating element maps to exactly one receiver. The owner is
//! backward Euler with held inputs, so every port answers
//! [`EndOfStep::Untold`].
//!
//! The mesh stands still: each step gives the owner a zero velocity at every
//! vertex, and the Component has no mechanical port. The receipt takes the
//! owner's stored change, its numerical (backward Euler) dissipation and the
//! heat that leaves through the heat ports. A step with mechanical work
//! refuses, and so does one with magnetic parameter work outside the
//! magnetic lane below: the receipt has no term for either. The
//! temperatures the heat ports receive enter the generated entropy and, with
//! laws, the conductors and graph resistors of the step: a
//! [`SpatialConductorLaw`] per receiver for its core cells and strands, and a
//! [`ResistorLaw`] per graph resistor through the owner's
//! `with_graph_resistors`, or none, keeping the circuit's own resistors. A
//! law with no answer refuses the step; the Component never clamps a
//! temperature.
//!
//! The owner reads its field energy only by recomputing the field, so the
//! state carries its stored energy: circuit, field and end-mode storage,
//! summed as the owner sums `stored_before_j`, when the state is built
//! ([`SpatialWindingState::new`]), and the owner's `stored_after_j` for a
//! state a step proposed. The owner sums a step's stored change separately,
//! so the kernel's endpoint check compares two sums of the same terms.
//!
//! In the magnetic lane ([`SpatialWinding::magnetic_at_port_temperatures`]),
//! each cell with a [`SpatialMagneticLaw`] also takes its law's linear
//! material at the temperature its receiver's port gives, at the end of the
//! step as the conductors do, and the state carries those temperatures
//! ([`SpatialWindingState::at_temperatures`]). These are `PlanarWinding`'s
//! lane terms, from the core the two share, with a cell's flux density the
//! length of its vector. A port's heat is its Joule heat less `T dS` of its
//! law cells, and may be negative. The stored energy adds `T s`, measured
//! from the state's own field and temperatures. The numerical loss adds each
//! law cell's `-(P + (T1 - T0) s0)`, second order in the step's change of
//! temperature, with `P` the owner's parameter work. The Component declares
//! its stored entropy: the sum of `s = -dpsi/dT` over the law cells,
//! measured from the state.
//!
//! The Component borrows the owner's problem for `'m` and leaks nothing. A
//! system of the kernel holds `'static` members, so a caller that adds one
//! to a system gives it a problem that lives for `'static`. The state has no
//! codec: a system holding it cannot be written to a durable snapshot.

use core::any::Any;
use core::fmt;
use core::marker::PhantomData;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::conductor::ResolvedConductorState;
use fs_material::conductor::spatial::{SpatialConductivity, SpatialConductorMaterial};
use fs_material::field_storage::admit_magnetic_temperature_refresh;
use fs_material::field_storage::thermal::magnetic::ThermalMagneticMaterial;
use fs_phs::circuit::CircuitGraph;
use fs_ports::CoupleError;
use fs_ports::component::{Class, Component, ContractViolation, EndOfStep, PortSpec, Proposal};

use super::ResistorLaw;
use super::shared::{
    self, Account, Element, ElementRefusal, LawCell, Refusal, Reversible, same_bits,
};
use crate::spatial::nonlinear::SpatialNonlinearMaterial;
use crate::spatial::vector_potential::induction::nonlinear::winding::{
    SpatialNonlinearWindingCircuit, SpatialNonlinearWindingError, SpatialNonlinearWindingStep,
    SpatialWindingConductor,
};
use crate::spatial::vector_potential::induction::nonlinear::{
    SpatialNonlinearInductionControls, SpatialNonlinearInductionMagneticCell,
};
use crate::winding::nonlinear::{NonlinearWindingControls, resistors::GraphResistorMaterial};

/// One thermal receiver: the dissipating elements whose heat leaves through
/// one heat port.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SpatialHeatReceiver {
    /// Mesh cells (tetrahedra) of stranded windings or core conductors.
    pub cells: Vec<usize>,
    /// Branch indices of native graph resistors.
    pub resistors: Vec<usize>,
}

/// The conductors of one receiver's cells at the temperature its heat port
/// received.
pub trait SpatialConductorLaw<C: SpatialConductorMaterial>: fmt::Debug + Send + Sync {
    /// The source-bound core conductor of tetrahedron `cell` at `kelvin`. The
    /// owner refuses one whose card, basis, keys or axes are not those of the
    /// field's own conductor there.
    ///
    /// # Errors
    /// Why the law has no answer there, for example a temperature outside the
    /// validity of the card behind it. The step then refuses with
    /// [`SpatialWindingError::Law`]; the temperature is never clamped.
    fn core(&self, cell: usize, kelvin: f64) -> Result<C, String>;

    /// Winding `winding`'s stranded conductor in tetrahedron `cell` at
    /// `kelvin`; several insulated bundles may share a cell. The owner
    /// refuses one whose card or claims are not those of the circuit's own
    /// strand there.
    ///
    /// # Errors
    /// As [`Self::core`].
    fn strand(
        &self,
        winding: usize,
        cell: usize,
        kelvin: f64,
    ) -> Result<ResolvedConductorState, String>;

    /// Words that identify the law. They join the Component's identity.
    fn identity_words(&self, words: &mut Vec<u64>);
}

/// One mesh cell whose permeability follows the temperature its receiver's
/// heat port gives.
#[derive(Clone, Debug, PartialEq)]
pub struct SpatialMagneticLaw {
    /// The tetrahedron.
    pub cell: usize,
    /// fs-material's linear magnetic law, affine in temperature. Its source
    /// is the field's own linear material in the cell.
    pub law: ThermalMagneticMaterial,
}

/// The element whose law had no answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpatialLawElement {
    /// A core conductor's tetrahedron.
    Core(usize),
    /// A winding's strand in a tetrahedron.
    Strand {
        /// The winding.
        winding: usize,
        /// The tetrahedron.
        cell: usize,
    },
    /// A native graph resistor's branch.
    Resistor(usize),
    /// A [`SpatialMagneticLaw`]'s tetrahedron: the law has no answer at the
    /// temperature, or its free energy there is not the field's.
    Permeability(usize),
}

/// Refusal of the 3D winding owner or of its Component boundary.
#[derive(Debug)]
pub enum SpatialWindingError {
    /// The owner's refusal, unchanged (boxed: it is large).
    Owner(Box<SpatialNonlinearWindingError>),
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
    /// magnetic lane or in a cell without a law, which the receipt has no
    /// term for.
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
    /// A law had no answer at the temperature its port received, or a
    /// magnetic law's free energy is not the field's.
    Law {
        /// The heat receiver whose port gave the temperature.
        receiver: usize,
        /// The element asked for.
        element: SpatialLawElement,
        /// The temperature \[K\].
        kelvin: f64,
        /// The law's reason.
        message: String,
    },
    /// Cancellation was requested.
    Cancelled,
}

impl fmt::Display for SpatialWindingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Owner(error) => write!(f, "3D winding owner refused: {error}"),
            Self::Port(violation) => write!(f, "3D winding port refused: {violation}"),
            Self::Identifier(error) => {
                write!(f, "3D winding port identifier refused: {error:?}")
            }
            Self::Invalid { what } => write!(f, "3D winding component refused: {what}"),
            Self::UnaccountedWork {
                parameter_j,
                mechanical_j,
            } => write!(
                f,
                "3D winding step did parameter work {parameter_j} J and mechanical work \
                 {mechanical_j} J, which the receipt has no term for"
            ),
            Self::Balance {
                defect_j,
                tolerance_j,
            } => write!(
                f,
                "3D winding energy account does not close: defect {defect_j} J beyond \
                 {tolerance_j} J"
            ),
            Self::Law {
                receiver,
                element,
                kelvin,
                message,
            } => write!(
                f,
                "3D winding law of {element:?} (receiver {receiver}) has no answer at \
                 {kelvin} K: {message}"
            ),
            Self::Cancelled => f.write_str("3D winding component cancelled"),
        }
    }
}

impl std::error::Error for SpatialWindingError {}

impl From<SpatialNonlinearWindingError> for SpatialWindingError {
    fn from(error: SpatialNonlinearWindingError) -> Self {
        Self::Owner(Box::new(error))
    }
}

const fn invalid(what: &'static str) -> SpatialWindingError {
    SpatialWindingError::Invalid { what }
}

/// A refusal of the shared core, in this Component's words: its mesh
/// elements are cells (tetrahedra).
fn refused(refusal: Refusal) -> SpatialWindingError {
    match refusal {
        Refusal::Invalid(what) => invalid(what),
        Refusal::Element(element) => invalid(match element {
            ElementRefusal::NoConductor => "a receiver names a cell that carries no conductor",
            ElementRefusal::HeatUnmapped => "the owner reported heat in a cell no receiver maps",
            ElementRefusal::LawRecordMissing => "the owner's magnetic records omit a law cell",
            ElementRefusal::NoLawElement => "the state's field has no such law cell",
        }),
        Refusal::Port(violation) => SpatialWindingError::Port(violation),
        Refusal::Identifier(error) => SpatialWindingError::Identifier(error),
        Refusal::UnaccountedWork {
            parameter_j,
            mechanical_j,
        } => SpatialWindingError::UnaccountedWork {
            parameter_j,
            mechanical_j,
        },
        Refusal::Balance {
            defect_j,
            tolerance_j,
        } => SpatialWindingError::Balance {
            defect_j,
            tolerance_j,
        },
        Refusal::Law {
            receiver,
            element,
            kelvin,
            message,
        } => SpatialWindingError::Law {
            receiver,
            element: match element {
                Element::Resistor(branch) => SpatialLawElement::Resistor(branch),
                Element::Permeability(cell) => SpatialLawElement::Permeability(cell),
            },
            kelvin,
            message,
        },
    }
}

/// Each receiver's cells and graph resistor branches.
fn maps(receivers: &[SpatialHeatReceiver]) -> Vec<(&[usize], &[usize])> {
    receivers
        .iter()
        .map(|receiver| (receiver.cells.as_slice(), receiver.resistors.as_slice()))
        .collect()
}

/// The cells that carry a core conductor or a strand, each once.
fn elements(core: impl Iterator<Item = usize>, strands: impl Iterator<Item = usize>) -> Vec<usize> {
    core.chain(strands)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

/// The length of a flux density vector \[T\].
fn magnitude([x, y, z]: [f64; 3]) -> f64 {
    x.hypot(y).hypot(z)
}

/// The winding owner whose accepted circuit a [`SpatialWinding`] steps:
/// [`SpatialNonlinearWindingCircuit`], the saturating field. Sealed: the
/// owners are fs-em's.
pub trait SpatialWindingOwner<'m, C: SpatialConductorMaterial>: owner::Sealed<'m, C> {}

impl<'m, C: SpatialConductorMaterial> SpatialWindingOwner<'m, C>
    for SpatialNonlinearWindingCircuit<'m, C>
{
}

mod owner {
    use super::{
        CircuitGraph, GraphResistorMaterial, NonlinearWindingControls, SpatialConductorMaterial,
        SpatialNonlinearInductionControls, SpatialNonlinearInductionMagneticCell,
        SpatialNonlinearMaterial, SpatialNonlinearWindingCircuit, SpatialNonlinearWindingError,
        SpatialNonlinearWindingStep, SpatialWindingConductor, fmt, magnitude, shared,
    };
    use crate::spatial::vector_potential::induction::nonlinear::SpatialNonlinearInductionStep;

    /// What the Component asks of its owner: one step from an accepted
    /// circuit at standstill, the rebinding of its graph resistors, the
    /// stored energy of an accepted circuit, the magnetic lane's materials
    /// and records, and the accepted coordinates that state equality reads.
    pub trait Sealed<'m, C: SpatialConductorMaterial>: Clone + Sized {
        /// The field's controls inside [`NonlinearWindingControls`].
        type FieldControls: Copy + fmt::Debug;
        /// The field's record inside [`SpatialNonlinearWindingStep`]. A
        /// proposed state keeps the step as `dyn Any`, so it borrows nothing.
        type FieldStep: 'static;

        /// The owner's `step_with_conductors_candidate` to `target_s` with
        /// these conductors and every vertex at rest, in parts; with
        /// `materials`, its `step_with_materials_candidate`.
        #[allow(clippy::too_many_arguments)]
        fn step(
            &self,
            target_s: f64,
            materials: Option<Vec<SpatialNonlinearMaterial>>,
            core: Vec<Option<C>>,
            windings: Vec<SpatialWindingConductor>,
            inputs: &[f64],
            controls: NonlinearWindingControls<Self::FieldControls>,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<
            (Self, SpatialNonlinearWindingStep<Self::FieldStep>),
            SpatialNonlinearWindingError,
        >;

        /// The owner's `with_graph_resistors`.
        fn rebind_graph_resistors(
            &self,
            sources: Vec<GraphResistorMaterial>,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<Self, SpatialNonlinearWindingError>;

        /// The mesh cell and heat \[J\] of each field cell of a step, in the
        /// owner's order.
        fn field_heat(step: &Self::FieldStep) -> Vec<(usize, f64)>;

        /// Field, circuit and end-mode storage, summed in the order the
        /// owner sums `stored_before_j`. The field's energy is recomputed
        /// from its potential.
        fn stored_energy_j(&self) -> Result<f64, SpatialNonlinearWindingError>;

        /// The field's materials, one per mesh cell.
        fn magnetic_materials(&self) -> &[SpatialNonlinearMaterial];

        /// Each mesh cell's flux density magnitude \[T\] and volume \[m³\],
        /// recovered from the accepted field's potential.
        fn flux_densities(&self) -> Result<Vec<(f64, f64)>, SpatialNonlinearWindingError>;

        /// A step's magnetic cell records, in mesh order.
        fn magnetic_cells(step: &Self::FieldStep) -> &[SpatialNonlinearInductionMagneticCell];

        /// Whether two states' problems hold the same mesh (positions and
        /// tetrahedra) and materials. A step that refreshes the materials
        /// builds its own problem.
        fn same_problem(&self, other: &Self) -> bool;

        fn graph(&self) -> &CircuitGraph;
        fn time_s(&self) -> f64;
        fn state(&self) -> &[f64];
        fn currents_a(&self) -> &[f64];
        fn flux_linkages_wb(&self) -> &[f64];
        fn resistance_ohm(&self) -> &[f64];
        fn conductors(&self) -> &[SpatialWindingConductor];
        /// The field's core conductors, one per mesh cell.
        fn core_conductors(&self) -> &[Option<C>];
        /// The address of the field's problem, which equal states share.
        fn problem(&self) -> *const ();
        fn potential_wb(&self) -> &[f64];
        fn field_time_s(&self) -> f64;
    }

    impl<'m, C: SpatialConductorMaterial> Sealed<'m, C> for SpatialNonlinearWindingCircuit<'m, C> {
        type FieldControls = SpatialNonlinearInductionControls;
        type FieldStep = SpatialNonlinearInductionStep;

        fn step(
            &self,
            target_s: f64,
            materials: Option<Vec<SpatialNonlinearMaterial>>,
            core: Vec<Option<C>>,
            windings: Vec<SpatialWindingConductor>,
            inputs: &[f64],
            controls: NonlinearWindingControls<SpatialNonlinearInductionControls>,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<(Self, SpatialNonlinearWindingStep), SpatialNonlinearWindingError> {
            let velocity = vec![[0.0; 3]; self.field().problem().mesh().vertex_count];
            let candidate = match materials {
                Some(materials) => self.step_with_materials_candidate(
                    target_s, &velocity, materials, core, windings, inputs, controls, cancelled,
                )?,
                None => self.step_with_conductors_candidate(
                    target_s, &velocity, core, windings, inputs, controls, cancelled,
                )?,
            };
            Ok(candidate.into_parts())
        }

        fn rebind_graph_resistors(
            &self,
            sources: Vec<GraphResistorMaterial>,
            cancelled: &mut impl FnMut() -> bool,
        ) -> Result<Self, SpatialNonlinearWindingError> {
            self.with_graph_resistors(sources, cancelled)
        }

        fn field_heat(step: &SpatialNonlinearInductionStep) -> Vec<(usize, f64)> {
            step.cells
                .iter()
                .map(|cell| (cell.cell, cell.heat_j))
                .collect()
        }

        fn stored_energy_j(&self) -> Result<f64, SpatialNonlinearWindingError> {
            let field = self.field().field(&mut || false)?.stored_energy_j();
            let end = self.end_inductance().stored_energy_j(self.currents_a())?;
            Ok(field + self.dae().system.hamiltonian(self.state()) + end)
        }

        fn magnetic_materials(&self) -> &[SpatialNonlinearMaterial] {
            self.field().problem().materials()
        }

        fn flux_densities(&self) -> Result<Vec<(f64, f64)>, SpatialNonlinearWindingError> {
            let field = self.field().field(&mut || false)?;
            let problem = self.field().problem();
            field
                .flux_density_t()
                .iter()
                .enumerate()
                .map(|(cell, &b)| {
                    problem
                        .cell_volume_m3(cell)
                        .map(|volume| (magnitude(b), volume))
                        .ok_or(SpatialNonlinearWindingError::Invalid(
                            "missing magnetic volume",
                        ))
                })
                .collect()
        }

        fn magnetic_cells(
            step: &SpatialNonlinearInductionStep,
        ) -> &[SpatialNonlinearInductionMagneticCell] {
            &step.magnetic_cells
        }

        fn same_problem(&self, other: &Self) -> bool {
            let (a, b) = (self.field().problem(), other.field().problem());
            core::ptr::eq(a, b)
                || (a.materials() == b.materials()
                    && a.mesh().tets == b.mesh().tets
                    && shared::same_vectors(a.positions_m(), b.positions_m()))
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

        fn conductors(&self) -> &[SpatialWindingConductor] {
            Self::conductors(self)
        }

        fn core_conductors(&self) -> &[Option<C>] {
            self.field().conductors()
        }

        fn problem(&self) -> *const () {
            core::ptr::from_ref(self.field().problem()).cast()
        }

        fn potential_wb(&self) -> &[f64] {
            self.field().potential_wb()
        }

        fn field_time_s(&self) -> f64 {
            self.field().time_s()
        }
    }
}

/// The accepted state: the owner's circuit, which carries the field history;
/// its stored energy \[J\]; in the magnetic lane, the temperatures its field
/// materials were evaluated at; and, for a state a step proposed, the owner's
/// step that produced it. The step is held as `dyn Any`, not through the
/// owner's associated type, so the state stays covariant in `'m`; a clone
/// shares it.
#[derive(Clone)]
pub struct SpatialWindingState<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
    O = SpatialNonlinearWindingCircuit<'m, C>,
> {
    circuit: O,
    stored_energy_j: f64,
    temperatures: Option<Vec<f64>>,
    step: Option<Rc<dyn Any>>,
    problem: PhantomData<&'m ()>,
    conductor: PhantomData<fn() -> C>,
}

impl<'m, C: SpatialConductorMaterial, O: SpatialWindingOwner<'m, C>> SpatialWindingState<'m, C, O> {
    /// The state of an accepted circuit, with its field, circuit and
    /// end-mode storage summed as the owner sums `stored_before_j`.
    ///
    /// # Errors
    /// The owner cannot recompute its field or end storage, or the sum is not
    /// finite.
    pub fn new(circuit: O) -> Result<Self, SpatialWindingError> {
        let stored_energy_j = circuit.stored_energy_j()?;
        if !stored_energy_j.is_finite() {
            return Err(invalid("the winding's stored energy is not finite"));
        }
        Ok(Self {
            circuit,
            stored_energy_j,
            temperatures: None,
            step: None,
            problem: PhantomData,
            conductor: PhantomData,
        })
    }

    /// The magnetic lane's state of an accepted circuit whose field
    /// materials were evaluated at `temperatures` \[K\], one per heat
    /// receiver, its storage summed as [`Self::new`] sums it. A step refuses
    /// a state whose temperatures do not fit its Component, and one whose
    /// law cells' free energy at them is not the field's.
    ///
    /// # Errors
    /// As [`Self::new`].
    pub fn at_temperatures(
        circuit: O,
        temperatures: Vec<f64>,
    ) -> Result<Self, SpatialWindingError> {
        let mut state = Self::new(circuit)?;
        state.temperatures = Some(temperatures);
        Ok(state)
    }

    /// The temperatures \[K\] the field materials were evaluated at, in the
    /// magnetic lane.
    #[must_use]
    pub fn temperatures(&self) -> Option<&[f64]> {
        self.temperatures.as_deref()
    }

    /// The accepted circuit.
    #[must_use]
    pub const fn circuit(&self) -> &O {
        &self.circuit
    }

    /// The accepted circuit, by value.
    #[must_use]
    pub fn into_circuit(self) -> O {
        self.circuit
    }

    /// The stored energy \[J\]: summed at construction, or the owner's
    /// `stored_after_j` of the step that proposed this state.
    #[must_use]
    pub const fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
}

impl<C: SpatialConductorMaterial> SpatialWindingState<'_, C> {
    /// The saturating owner's step that produced this state: its cells, graph
    /// resistors, field record, energies and solve report, as the owner's own
    /// `step_with_conductors_candidate` gave them. None for a constructed
    /// state.
    #[must_use]
    pub fn step(&self) -> Option<&SpatialNonlinearWindingStep> {
        self.step.as_deref()?.downcast_ref()
    }
}

/// Equal states share the owner's problem (the same allocation) and graph,
/// and hold the same bits in every accepted coordinate: clock, circuit
/// state, currents, linkages, resistances, strands, core conductors and
/// field potential and clock. The graph is compared by value: a step with
/// resistor laws builds a new one. In the magnetic lane, each step builds
/// its own problem too, so there the problems are compared by mesh and
/// materials, and the temperatures by their bits; a state with temperatures
/// never equals one without. Neither the stored energy a state carries nor
/// the owner's step it keeps is compared: a proposed state equals the
/// constructed one.
impl<'m, C: SpatialConductorMaterial, O: SpatialWindingOwner<'m, C>> PartialEq
    for SpatialWindingState<'m, C, O>
{
    fn eq(&self, other: &Self) -> bool {
        let (a, b) = (&self.circuit, &other.circuit);
        let problem = match (&self.temperatures, &other.temperatures) {
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
            && same_bits(a.potential_wb(), b.potential_wb())
            && a.field_time_s().to_bits() == b.field_time_s().to_bits()
    }
}

impl<'m, C: SpatialConductorMaterial, O: SpatialWindingOwner<'m, C>> fmt::Debug
    for SpatialWindingState<'m, C, O>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SpatialWindingState")
            .field("time_s", &self.circuit.time_s())
            .field("currents_a", &self.circuit.currents_a())
            .field("state", &self.circuit.state())
            .field("stored_energy_j", &self.stored_energy_j)
            .field("temperatures_k", &self.temperatures)
            .finish_non_exhaustive()
    }
}

/// A strand of the construction's circuit, whose conductor laws give.
struct Strand {
    winding: usize,
    cell: usize,
    fill_fraction: f64,
}

/// A step's circuit, when resistor laws rebuilt it, and its core conductors
/// and strands.
type Supplied<C, O> = (Option<O>, Vec<Option<C>>, Vec<SpatialWindingConductor>);

/// Where a step's conductors and graph resistances come from.
enum Supply<C: SpatialConductorMaterial> {
    /// Given at construction and used at every step.
    Fixed {
        core: Vec<Option<C>>,
        windings: Vec<SpatialWindingConductor>,
    },
    /// Laws asked at the temperatures the heat ports receive.
    Laws {
        /// The number of mesh cells.
        cells: usize,
        /// The cells with a core conductor, as the construction's circuit
        /// holds them.
        core: Vec<usize>,
        /// The construction's circuit's strands, in its order.
        strands: Vec<Strand>,
        /// One per heat receiver.
        conductors: Vec<Arc<dyn SpatialConductorLaw<C>>>,
        /// One per native graph resistor, by branch, or none: the circuit's
        /// own resistors stay.
        resistors: BTreeMap<usize, Arc<dyn ResistorLaw>>,
    },
}

/// The 3D winding owner behind the composition contract's source and heat
/// ports. See the module documentation. `C` is the core conductors' type,
/// `O` the owner's circuit and `F` its field's controls; `O` and `F` default
/// to the saturating owner.
pub struct SpatialWinding<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
    O = SpatialNonlinearWindingCircuit<'m, C>,
    F = SpatialNonlinearInductionControls,
> {
    graph: CircuitGraph,
    sources: usize,
    supply: Supply<C>,
    receivers: Vec<SpatialHeatReceiver>,
    cell_port: BTreeMap<usize, usize>,
    resistor_port: BTreeMap<usize, usize>,
    controls: NonlinearWindingControls<F>,
    energy_tolerance_j: f64,
    ports: Vec<PortSpec>,
    /// The magnetic lane's laws, by cell.
    magnetic: Option<BTreeMap<usize, ThermalMagneticMaterial>>,
    problem: PhantomData<&'m ()>,
    owner: PhantomData<fn() -> O>,
}

impl<C: SpatialConductorMaterial, O, F> fmt::Debug for SpatialWinding<'_, C, O, F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SpatialWinding")
            .field("sources", &self.sources)
            .field("receivers", &self.receivers)
            .field("energy_tolerance_j", &self.energy_tolerance_j)
            .finish_non_exhaustive()
    }
}

impl<'m, C, O, F> SpatialWinding<'m, C, O, F>
where
    C: SpatialConductorMaterial,
    O: SpatialWindingOwner<'m, C, FieldControls = F>,
    F: Copy + fmt::Debug,
{
    /// The Component of `circuit`'s topology, stepped with these core
    /// conductors (one per mesh cell) and strands and these controls, its
    /// heat leaving through `receivers`.
    ///
    /// # Errors
    /// The energy tolerance is not positive and finite; an external port
    /// drives no source or two; a receiver is empty; a core or stranded cell
    /// or a graph resistor is mapped to no receiver or to two; a receiver
    /// names an element that does not dissipate.
    pub fn new(
        circuit: &O,
        core_conductors: Vec<Option<C>>,
        winding_conductors: Vec<SpatialWindingConductor>,
        receivers: Vec<SpatialHeatReceiver>,
        controls: NonlinearWindingControls<F>,
        energy_tolerance_j: f64,
    ) -> Result<Self, SpatialWindingError> {
        let elements = elements(
            core_conductors
                .iter()
                .enumerate()
                .filter_map(|(cell, conductor)| conductor.as_ref().map(|_| cell)),
            winding_conductors.iter().map(|strand| strand.cell),
        );
        let supply = Supply::Fixed {
            core: core_conductors,
            windings: winding_conductors,
        };
        Self::build(
            circuit,
            &elements,
            supply,
            receivers,
            controls,
            energy_tolerance_j,
        )
    }

    /// The Component of `circuit`'s topology whose conductors and graph
    /// resistors follow the temperatures its heat ports receive: receiver
    /// `k`'s law gives the core conductor and the strands of each of its
    /// cells at port `k`'s temperature, and each graph resistor's law gives
    /// it at the temperature of the port its heat leaves through. With no
    /// resistor laws, the graph resistors keep `circuit`'s values. The core
    /// cells, and each strand's winding, cell and fill, are `circuit`'s own.
    ///
    /// # Errors
    /// As [`Self::new`]; the number of conductor laws is not the number of
    /// receivers; some graph resistors have laws and others do not, or a law
    /// names a branch that is not a resistor.
    pub fn at_port_temperatures(
        circuit: &O,
        conductor_laws: Vec<Arc<dyn SpatialConductorLaw<C>>>,
        resistor_laws: BTreeMap<usize, Arc<dyn ResistorLaw>>,
        receivers: Vec<SpatialHeatReceiver>,
        controls: NonlinearWindingControls<F>,
        energy_tolerance_j: f64,
    ) -> Result<Self, SpatialWindingError> {
        let (supply, elements) =
            Self::laws(circuit, conductor_laws, resistor_laws, receivers.len())?;
        Self::build(
            circuit,
            &elements,
            supply,
            receivers,
            controls,
            energy_tolerance_j,
        )
    }

    /// [`Self::at_port_temperatures`] in the magnetic lane (module
    /// documentation): the permeability of each cell with a
    /// [`SpatialMagneticLaw`] also follows the temperature its receiver's
    /// port gives. Each law cell maps to a receiver, conducting or not. The
    /// Component steps states made by [`SpatialWindingState::at_temperatures`].
    ///
    /// # Errors
    /// As [`Self::at_port_temperatures`]; there is no law; two laws name one
    /// cell; a law names a cell whose field material is not linear, or whose
    /// material is not the law's source (its card, claims and other
    /// coordinates); a law cell maps to no receiver.
    #[allow(clippy::too_many_arguments)]
    pub fn magnetic_at_port_temperatures(
        circuit: &O,
        conductor_laws: Vec<Arc<dyn SpatialConductorLaw<C>>>,
        resistor_laws: BTreeMap<usize, Arc<dyn ResistorLaw>>,
        magnetic: Vec<SpatialMagneticLaw>,
        receivers: Vec<SpatialHeatReceiver>,
        controls: NonlinearWindingControls<F>,
        energy_tolerance_j: f64,
    ) -> Result<Self, SpatialWindingError> {
        let materials = circuit.magnetic_materials();
        let mut laws = BTreeMap::new();
        for SpatialMagneticLaw { cell, law } in magnetic {
            let Some(SpatialNonlinearMaterial::Linear(material)) = materials.get(cell) else {
                return Err(invalid(
                    "a magnetic law names a cell whose field material is not linear",
                ));
            };
            if admit_magnetic_temperature_refresh(material, law.material()).is_err() {
                return Err(invalid(
                    "a magnetic law's source is not its cell's field material",
                ));
            }
            if laws.insert(cell, law).is_some() {
                return Err(invalid("two magnetic laws name one cell"));
            }
        }
        if laws.is_empty() {
            return Err(invalid("the magnetic lane has no law"));
        }
        let (supply, mut elements) =
            Self::laws(circuit, conductor_laws, resistor_laws, receivers.len())?;
        elements.extend(laws.keys().copied());
        let mut component = Self::build(
            circuit,
            &elements,
            supply,
            receivers,
            controls,
            energy_tolerance_j,
        )?;
        component.magnetic = Some(laws);
        Ok(component)
    }

    /// The supply of [`Self::at_port_temperatures`], and the cells its
    /// receivers must map.
    fn laws(
        circuit: &O,
        conductor_laws: Vec<Arc<dyn SpatialConductorLaw<C>>>,
        resistor_laws: BTreeMap<usize, Arc<dyn ResistorLaw>>,
        receivers: usize,
    ) -> Result<(Supply<C>, Vec<usize>), SpatialWindingError> {
        if conductor_laws.len() != receivers {
            return Err(invalid("one conductor law per heat receiver required"));
        }
        let core: Vec<usize> = circuit
            .core_conductors()
            .iter()
            .enumerate()
            .filter_map(|(cell, conductor)| conductor.as_ref().map(|_| cell))
            .collect();
        let strands: Vec<Strand> = circuit
            .conductors()
            .iter()
            .map(|strand| Strand {
                winding: strand.winding,
                cell: strand.cell,
                fill_fraction: strand.fill_fraction,
            })
            .collect();
        let elements = elements(
            core.iter().copied(),
            strands.iter().map(|strand| strand.cell),
        );
        let supply = Supply::Laws {
            cells: circuit.core_conductors().len(),
            core,
            strands,
            conductors: conductor_laws,
            resistors: resistor_laws,
        };
        Ok((supply, elements))
    }

    fn build(
        circuit: &O,
        elements: &[usize],
        supply: Supply<C>,
        receivers: Vec<SpatialHeatReceiver>,
        controls: NonlinearWindingControls<F>,
        energy_tolerance_j: f64,
    ) -> Result<Self, SpatialWindingError> {
        let laws = match &supply {
            Supply::Laws { resistors, .. } => Some(resistors),
            Supply::Fixed { .. } => None,
        };
        let layout = shared::layout(
            circuit.graph(),
            elements,
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
            cell_port: layout.element_port,
            resistor_port: layout.resistor_port,
            controls,
            energy_tolerance_j,
            ports: layout.ports,
            magnetic: None,
            problem: PhantomData,
            owner: PhantomData,
        })
    }

    /// The declared heat receivers, in heat-port order.
    #[must_use]
    pub fn receivers(&self) -> &[SpatialHeatReceiver] {
        &self.receivers
    }

    /// The owner's clock at the end of a step of `dt_s` from `old`, when the
    /// step, its ports and the state's topology are this Component's.
    fn admit(
        &self,
        old: &SpatialWindingState<'m, C, O>,
        dt_s: f64,
        inputs: &[f64],
        outputs: usize,
    ) -> Result<f64, SpatialWindingError> {
        shared::admit_ports(dt_s, inputs, outputs, self.sources, self.receivers.len())
            .map_err(refused)?;
        let laws = matches!(&self.supply, Supply::Laws { resistors, .. } if !resistors.is_empty());
        if !shared::same_graph(old.circuit.graph(), &self.graph, laws) {
            return Err(invalid(
                "the state's circuit is not this Component's topology",
            ));
        }
        let lane = self.magnetic.as_ref().map(|_| self.receivers.len());
        shared::admit_temperatures(old.temperatures.as_deref(), lane).map_err(refused)?;
        shared::target_time(old.circuit.time_s(), dt_s).map_err(refused)
    }

    /// The circuit to step, when resistor laws rebuild it, and the core
    /// conductors and strands of a step whose heat ports received `kelvin`.
    fn supply_at(
        &self,
        old: &O,
        kelvin: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Supplied<C, O>, SpatialWindingError> {
        let (cells, core_cells, strands, conductors, resistors) = match &self.supply {
            Supply::Fixed { core, windings } => return Ok((None, core.clone(), windings.clone())),
            Supply::Laws {
                cells,
                core,
                strands,
                conductors,
                resistors,
            } => (cells, core, strands, conductors, resistors),
        };
        let law = |cell: usize| {
            let receiver = self.cell_port[&cell];
            (receiver, kelvin[receiver], &conductors[receiver])
        };
        let mut core = vec![None; *cells];
        for &cell in core_cells {
            let (receiver, kelvin, law) = law(cell);
            let conductor = law
                .core(cell, kelvin)
                .map_err(|message| SpatialWindingError::Law {
                    receiver,
                    element: SpatialLawElement::Core(cell),
                    kelvin,
                    message,
                })?;
            core[cell] = Some(conductor);
        }
        let mut windings = Vec::with_capacity(strands.len());
        for strand in strands {
            let (receiver, kelvin, law) = law(strand.cell);
            let conductor = law
                .strand(strand.winding, strand.cell, kelvin)
                .map_err(|message| SpatialWindingError::Law {
                    receiver,
                    element: SpatialLawElement::Strand {
                        winding: strand.winding,
                        cell: strand.cell,
                    },
                    kelvin,
                    message,
                })?;
            windings.push(SpatialWindingConductor {
                winding: strand.winding,
                cell: strand.cell,
                fill_fraction: strand.fill_fraction,
                conductor,
            });
        }
        if resistors.is_empty() {
            return Ok((None, core, windings));
        }
        let sources =
            shared::graph_resistors(resistors, &self.resistor_port, kelvin).map_err(refused)?;
        let driven = old.rebind_graph_resistors(sources, cancelled)?;
        Ok((Some(driven), core, windings))
    }

    /// In the magnetic lane, the field materials of a step whose heat ports
    /// received `kelvin`: each law cell's law there, and every other cell's
    /// material as `old` holds it.
    fn materials_at(
        &self,
        old: &O,
        kelvin: &[f64],
    ) -> Result<Option<Vec<SpatialNonlinearMaterial>>, SpatialWindingError> {
        let Some(laws) = &self.magnetic else {
            return Ok(None);
        };
        let mut materials = old.magnetic_materials().to_vec();
        for (&cell, law) in laws {
            let receiver = self.cell_port[&cell];
            let kelvin = kelvin[receiver];
            let at = law
                .at_temperature(kelvin)
                .map_err(|error| SpatialWindingError::Law {
                    receiver,
                    element: SpatialLawElement::Permeability(cell),
                    kelvin,
                    message: error.to_string(),
                })?;
            *materials
                .get_mut(cell)
                .ok_or_else(|| invalid("the state's field has no such law cell"))? =
                SpatialNonlinearMaterial::Linear(at.material().clone());
        }
        Ok(Some(materials))
    }

    /// The magnetic lane's terms of a step from `old`'s temperatures to the
    /// ports' (`kelvin`), read from the owner's magnetic cell records. Each
    /// law's free energy must be the field's at both ends of the step, within
    /// the energy tolerance; a cell without a law must do no parameter work.
    fn reversible(
        &self,
        old: &SpatialWindingState<'m, C, O>,
        cells: &[SpatialNonlinearInductionMagneticCell],
        kelvin: &[f64],
    ) -> Result<Option<Reversible>, SpatialWindingError> {
        let Some(laws) = &self.magnetic else {
            return Ok(None);
        };
        let before = old
            .temperatures
            .as_deref()
            .ok_or_else(|| invalid("a state of the magnetic lane carries its temperatures"))?;
        let cells: Vec<LawCell> = cells
            .iter()
            .map(|cell| LawCell {
                element: cell.cell,
                flux_density_t: cell.flux_density_t.map(magnitude),
                free_energy_j: cell.free_energy_j,
                volume_m3: cell.volume_m3,
                parameter_work_j: cell.parameter_work_j,
            })
            .collect();
        shared::reversible(
            laws,
            &self.cell_port,
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
    /// over the law cells, measured from the state alone: each law at its
    /// receiver's temperature, at the flux density the field's potential
    /// gives. Both sums are compensated (`two_sum`).
    fn stored_reversible(
        &self,
        state: &SpatialWindingState<'m, C, O>,
    ) -> Result<Option<(f64, f64)>, SpatialWindingError> {
        let Some(laws) = &self.magnetic else {
            return Ok(None);
        };
        let temperatures = state
            .temperatures
            .as_deref()
            .filter(|t| t.len() == self.receivers.len())
            .ok_or_else(|| invalid("a state of the magnetic lane carries its temperatures"))?;
        let cells = state.circuit.flux_densities()?;
        shared::stored_reversible(laws, &self.cell_port, temperatures, &cells)
            .map(Some)
            .map_err(refused)
    }
}

/// A core conductor's conductivity (scalar or tensor), principal axes,
/// source basis and keys, and card identity.
fn spatial_conductor_words<C: SpatialConductorMaterial>(conductor: &C, words: &mut Vec<u64>) {
    match conductor.conductivity() {
        SpatialConductivity::Scalar(sigma) => words.extend([0, sigma.to_bits()]),
        SpatialConductivity::Tensor(sigma) => {
            words.push(1);
            words.extend(sigma.iter().flatten().map(|value| value.to_bits()));
        }
    }
    match conductor.principal_axes() {
        None => words.push(0),
        Some(axes) => {
            words.push(1);
            words.extend(axes.iter().flatten().map(|value| value.to_bits()));
        }
    }
    shared::debug_words(&conductor.source_basis(), words);
    shared::debug_words(&conductor.source_keys().collect::<Vec<_>>(), words);
    for chunk in conductor.material().card_identity().0.chunks(8) {
        let mut bytes = [0_u8; 8];
        bytes.copy_from_slice(chunk);
        words.push(u64::from_le_bytes(bytes));
    }
}

impl<'m, C, O, F> Component for SpatialWinding<'m, C, O, F>
where
    C: SpatialConductorMaterial,
    O: SpatialWindingOwner<'m, C, FieldControls = F>,
    F: Copy + fmt::Debug,
{
    type State = SpatialWindingState<'m, C, O>;
    type Error = SpatialWindingError;

    /// `b"spatial1"` read little-endian; the graph; the fixed core
    /// conductors (the number of cells, then each conducting cell and its
    /// conductor's conductivity, axes, source basis and keys, and card
    /// identity) and strands (winding, cell, fill bits, resistivity bits and
    /// card identity), or the laws' core cells, strands and words; the
    /// receiver map; the controls by their `Debug` text (which prints each
    /// `f64` so that it parses back to the same bits); the energy tolerance;
    /// then, in the magnetic lane only, `b"magtemp1"` and each law's cell and
    /// words (its card identity, temperature axis, range and slope, and its
    /// source's permeability). Each list is preceded by its length.
    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(u64::from_le_bytes(*b"spatial1"));
        shared::graph_words(&self.graph, words);
        match &self.supply {
            Supply::Fixed { core, windings } => {
                let conducting = core.iter().filter(|conductor| conductor.is_some()).count();
                words.extend([0, core.len() as u64, conducting as u64]);
                for (cell, conductor) in core.iter().enumerate() {
                    if let Some(conductor) = conductor {
                        words.push(cell as u64);
                        spatial_conductor_words(conductor, words);
                    }
                }
                words.push(windings.len() as u64);
                for strand in windings {
                    words.extend([
                        strand.winding as u64,
                        strand.cell as u64,
                        strand.fill_fraction.to_bits(),
                    ]);
                    shared::conductor_words(&strand.conductor, words);
                }
            }
            Supply::Laws {
                cells,
                core,
                strands,
                conductors,
                resistors,
            } => {
                words.extend([1, *cells as u64, core.len() as u64]);
                words.extend(core.iter().map(|&cell| cell as u64));
                words.push(strands.len() as u64);
                for strand in strands {
                    words.extend([
                        strand.winding as u64,
                        strand.cell as u64,
                        strand.fill_fraction.to_bits(),
                    ]);
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
        if let Some(laws) = &self.magnetic {
            shared::lane_words(laws, words);
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

    /// The energy the state carries: summed at construction, or the owner's
    /// `stored_after_j` of the step that proposed it; in the magnetic lane,
    /// with `T s` of each law cell, measured from the state.
    fn stored_energy_j(
        &self,
        state: &SpatialWindingState<'m, C, O>,
    ) -> Result<f64, SpatialWindingError> {
        match self.stored_reversible(state)? {
            None => Ok(state.stored_energy_j),
            Some((reversible, _)) => Ok(state.stored_energy_j + reversible),
        }
    }

    /// In the magnetic lane, the sum of `s = -dpsi/dT` over the law cells,
    /// measured from the state's own field and temperatures; otherwise none.
    fn stored_entropy_j_per_k(
        &self,
        state: &SpatialWindingState<'m, C, O>,
    ) -> Result<Option<f64>, SpatialWindingError> {
        Ok(self.stored_reversible(state)?.map(|(_, entropy)| entropy))
    }

    fn propose(
        &self,
        old: &SpatialWindingState<'m, C, O>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut SpatialWindingState<'m, C, O>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, SpatialWindingError> {
        cx.checkpoint()
            .map_err(|_| SpatialWindingError::Cancelled)?;
        let target_s = self.admit(old, dt_s, inputs, outputs.len())?;
        let sources = self.sources;
        let kelvin = &inputs[sources..];
        let (driven, core, windings) =
            self.supply_at(&old.circuit, kelvin, &mut || cx.checkpoint().is_err())?;
        let materials = self.materials_at(&old.circuit, kelvin)?;
        let lane = materials.is_some();
        let (circuit, step) = driven.as_ref().unwrap_or(&old.circuit).step(
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
        let reversible = self.reversible(old, O::magnetic_cells(&step.field), kelvin)?;
        let field = O::field_heat(&step.field);
        let heat = shared::heat_by_receiver(
            self.receivers.len(),
            &self.cell_port,
            &self.resistor_port,
            step.cells
                .iter()
                .map(|cell| (cell.cell, cell.heat_j))
                .chain(field.iter().copied()),
            step.graph_resistors
                .iter()
                .map(|resistor| (resistor.branch, resistor.heat_j)),
        )
        .map_err(refused)?;
        let account = Account {
            stored_before_j: step.stored_before_j,
            stored_after_j: step.stored_after_j,
            stored_change_j: step.stored_change_j,
            supplied_work_j: step.supplied_work_j,
            numerical_dissipation_j: step.numerical_dissipation_j,
            elements: step.cells.len() + field.len() + step.graph_resistors.len(),
            coordinates: [old.circuit.state().len(), old.circuit.potential_wb().len()],
        };
        let laws = self.magnetic.as_ref().map_or(0, BTreeMap::len);
        let closed = shared::close(
            &account,
            &heat,
            kelvin,
            reversible.as_ref().map(|terms| (terms, laws)),
            sources,
            self.energy_tolerance_j,
        )
        .map_err(refused)?;
        let stored_energy_j = step.stored_after_j;
        if !stored_energy_j.is_finite() {
            return Err(invalid("the winding's stored energy is not finite"));
        }
        cx.checkpoint()
            .map_err(|_| SpatialWindingError::Cancelled)?;
        let receipt = shared::publish(
            &closed,
            &step.report,
            &self.controls.newton,
            &step.native.y,
            sources,
            dt_s,
            outputs,
        );
        *new = SpatialWindingState {
            circuit,
            stored_energy_j,
            temperatures: lane.then(|| kelvin.to_vec()),
            step: Some(kept),
            problem: PhantomData,
            conductor: PhantomData,
        };
        Ok(Proposal::Step(receipt))
    }
}
