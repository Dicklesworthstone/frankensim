//! Ordered, source-owned contact resistance in the existing implicit wall graph.

use super::*;
use crate::conjugate::FaceContactResistance;
use fs_conduction::enthalpy::conduction::implicit::contact::{
    ImplicitContactLaw, ImplicitContactState,
};

/// One resistance in series with the two actual voxel half cells.
///
/// A source retains its original caloric objects, selected interface response,
/// applicability domain and fixed external conditions, including pressure.
/// Its cell order is the source's A/B order, independent of voxel numbering.
/// Both temperature/phase tangents enter the SAME native implicit residual.
/// Unlisted faces retain perfect contact. Each selected face must touch a
/// participating finite solid; a fixed resistance and a source cannot both
/// own that face. This law has no separate interfacial thermal storage.
#[derive(Clone, Copy)]
pub enum ThermalVolumeContact<'a> {
    Resistance(FaceContactResistance),
    Source {
        cells: [usize; 2],
        law: &'a dyn ImplicitContactLaw,
    },
}

impl From<FaceContactResistance> for ThermalVolumeContact<'_> {
    fn from(value: FaceContactResistance) -> Self {
        Self::Resistance(value)
    }
}

/// Only the view is reordered. The underlying source/card and its fixed
/// conditions remain untouched, as do its analytic enthalpy derivatives.
struct OrderedSource<'a> {
    law: &'a dyn ImplicitContactLaw,
    reversed: bool,
}

impl ImplicitContactLaw for OrderedSource<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        let mut laws = self.law.caloric_laws();
        if self.reversed {
            laws.swap(0, 1);
        }
        laws
    }

    fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        let mut bounds = self.law.enthalpy_bounds_j_kg();
        if self.reversed {
            bounds.swap(0, 1);
        }
        bounds
    }

    fn state_at_specific_enthalpies(
        &self,
        mut h: [f64; 2],
    ) -> Result<ImplicitContactState, DepositionError> {
        if self.reversed {
            h.swap(0, 1);
        }
        let mut state = self.law.state_at_specific_enthalpies(h)?;
        if self.reversed {
            state.resistance_enthalpy_derivatives_m2_kg_k_w_j.swap(0, 1);
        }
        Ok(state)
    }
}

enum PreparedContact<'a> {
    Resistance(f64),
    Source(OrderedSource<'a>),
}

pub(super) struct PreparedContacts<'a> {
    faces: BTreeMap<(usize, usize), PreparedContact<'a>>,
    source_count: usize,
}

impl<'a> PreparedContacts<'a> {
    pub(super) fn new(
        domain: &VoxelDomain,
        contacts: &[ThermalVolumeContact<'a>],
        samples: &[Option<PhaseConductivitySample>],
        wall_at: &[Option<usize>],
        max_faces: usize,
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        if contacts.len() > max_faces {
            return Err(refuse("finite-solid contact work budget exhausted"));
        }
        let mut faces = BTreeMap::new();
        let mut source_count = 0;
        for (index, contact) in contacts.iter().enumerate() {
            poll_every(gate, index)?;
            let cells = match contact {
                ThermalVolumeContact::Resistance(fixed) => fixed.cells,
                ThermalVolumeContact::Source { cells, .. } => *cells,
            };
            let [a, b] = cells;
            if a >= domain.cell_count()
                || b >= domain.cell_count()
                || a == b
                || samples.get(a).is_none_or(Option::is_none)
                || samples.get(b).is_none_or(Option::is_none)
                || (wall_at.get(a).is_none_or(Option::is_none)
                    && wall_at.get(b).is_none_or(Option::is_none))
                || !(0..6).any(|side| domain.neighbor(a, side) == Some(b))
            {
                return Err(refuse(
                    "finite-solid contacts need adjacent active faces touching a participating wall",
                ));
            }
            let value = match contact {
                ThermalVolumeContact::Resistance(fixed) => {
                    if !fixed.resistance_m2_k_w.is_finite() || fixed.resistance_m2_k_w < 0.0 {
                        return Err(refuse("contact resistance must be finite and nonnegative"));
                    }
                    PreparedContact::Resistance(fixed.resistance_m2_k_w)
                }
                ThermalVolumeContact::Source { law, .. } => {
                    source_count += 1;
                    PreparedContact::Source(OrderedSource {
                        law: *law,
                        reversed: a > b,
                    })
                }
            };
            if faces.insert((a.min(b), a.max(b)), value).is_some() {
                return Err(refuse(
                    "one thermal face cannot own duplicate contact resistances",
                ));
            }
        }
        poll(gate)?;
        Ok(Self {
            faces,
            source_count,
        })
    }

    /// Source contacts leave the geometric resistance at zero; the native
    /// graph adds their state-dependent resistance exactly once.
    pub(super) fn fixed_resistance(&self, a: usize, b: usize) -> f64 {
        match self.faces.get(&(a.min(b), a.max(b))) {
            Some(PreparedContact::Resistance(value)) => *value,
            _ => 0.0,
        }
    }

    pub(super) fn has_sources(&self) -> bool {
        self.source_count != 0
    }

    /// The bridge creates native links in increasing voxel order. Native
    /// attachment then checks the exact ordered caloric identities and domains.
    pub(super) fn bindings(
        &self,
        initial: &DepositionSimulation<'_>,
        gate: &CancelGate,
    ) -> Result<Vec<(u64, &dyn ImplicitContactLaw)>, ChtError> {
        poll(gate)?;
        if self.source_count == 0 {
            return Ok(Vec::new());
        }
        let mut bindings = Vec::with_capacity(self.source_count);
        for (index, (&id, link)) in initial.links().iter().enumerate() {
            poll_every(gate, index)?;
            let a =
                usize::try_from(link.a.0).map_err(|_| refuse("invalid contact cell identity"))?;
            let b =
                usize::try_from(link.b.0).map_err(|_| refuse("invalid contact cell identity"))?;
            if let Some(PreparedContact::Source(law)) = self.faces.get(&(a.min(b), a.max(b))) {
                if a >= b {
                    return Err(refuse(
                        "implicit voxel links must retain increasing cell order",
                    ));
                }
                bindings.push((id, law as &dyn ImplicitContactLaw));
            }
        }
        if bindings.len() != self.source_count {
            return Err(refuse(
                "every source contact must own an active native thermal link",
            ));
        }
        poll(gate)?;
        Ok(bindings)
    }

    /// Reconstruct physical heat from the SAME endpoint law. Call before any
    /// zero-conductivity shortcut, so insulating bulk cannot hide a source refusal.
    pub(super) fn resistance_at(
        &self,
        a: usize,
        b: usize,
        mut h: [f64; 2],
        gate: &CancelGate,
    ) -> Result<f64, ChtError> {
        match self.faces.get(&(a.min(b), a.max(b))) {
            None => Ok(0.0),
            Some(PreparedContact::Resistance(value)) => Ok(*value),
            Some(PreparedContact::Source(law)) => {
                poll(gate)?;
                if a > b {
                    h.swap(0, 1);
                }
                if h.iter()
                    .zip(law.enthalpy_bounds_j_kg())
                    .any(|(&value, [lo, hi])| !value.is_finite() || !(lo..=hi).contains(&value))
                {
                    return Err(refuse("contact endpoint left its original enthalpy domain"));
                }
                let state = law
                    .state_at_specific_enthalpies(h)
                    .map_err(|error| native_error(gate, error))?;
                state.admit().map_err(|error| native_error(gate, error))?;
                poll(gate)?;
                Ok(state.resistance_m2_k_w)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Sensible(f64);
    impl ThermalLaw for Sensible {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            Ok(300.0 + self.0 * h)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            self.0
        }
    }

    struct TwoSided<'a>([&'a dyn ThermalLaw; 2]);
    impl ImplicitContactLaw for TwoSided<'_> {
        fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
            self.0
        }
        fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
            [[-2.0, 3.0], [4.0, 9.0]]
        }
        fn state_at_specific_enthalpies(
            &self,
            h: [f64; 2],
        ) -> Result<ImplicitContactState, DepositionError> {
            Ok(ImplicitContactState {
                resistance_m2_k_w: 10.0 + 2.0 * h[0] + 3.0 * h[1],
                resistance_enthalpy_derivatives_m2_kg_k_w_j: [2.0, 3.0],
            })
        }
    }

    #[test]
    fn reversed_source_keeps_original_objects_domains_and_both_derivatives() {
        let a = Sensible(2.0);
        let b = Sensible(3.0);
        let source = TwoSided([&a, &b]);
        let ordered = OrderedSource {
            law: &source,
            reversed: true,
        };
        let curves = ordered.caloric_laws();
        assert!(std::ptr::addr_eq(curves[0], &b as &dyn ThermalLaw));
        assert!(std::ptr::addr_eq(curves[1], &a as &dyn ThermalLaw));
        assert_eq!(ordered.enthalpy_bounds_j_kg(), [[4.0, 9.0], [-2.0, 3.0]]);
        let state = ordered.state_at_specific_enthalpies([6.0, 1.0]).unwrap();
        assert_eq!(state.resistance_m2_k_w, 30.0);
        assert_eq!(
            state.resistance_enthalpy_derivatives_m2_kg_k_w_j,
            [3.0, 2.0]
        );
        let contacts = PreparedContacts {
            faces: BTreeMap::from([((0, 1), PreparedContact::Source(ordered))]),
            source_count: 1,
        };
        let gate = CancelGate::new_clock_free();
        assert_eq!(contacts.fixed_resistance(0, 1), 0.0);
        assert_eq!(
            contacts.resistance_at(0, 1, [6.0, 1.0], &gate).unwrap(),
            30.0
        );
        assert_eq!(
            contacts.resistance_at(1, 0, [1.0, 6.0], &gate).unwrap(),
            30.0
        );
        // The endpoint audit enforces the native domain even when a custom
        // constitutive implementation itself omits a range check.
        assert!(contacts.resistance_at(0, 1, [3.0, 1.0], &gate).is_err());
        assert!(contacts.resistance_at(1, 0, [4.0, 6.0], &gate).is_err());
    }
}
