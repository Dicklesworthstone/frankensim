//! Original fixed source contexts over an independently meshed interface.
//! The A trace owns the map; overlap integration keeps that ownership even
//! when a B face crosses several source patches. Only temperatures evolve.

use super::*;

/// Original source selection for one explicitly named A-side face.
/// The complete query retains pressure, finish and every fixed source axis.
/// Later pointwise temperature queries pin this exact selected claim.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactFaceReference {
    /// Explicit original boundary-face slot on the map-owning A trace.
    pub side_a: usize,
    /// Full admitted reference query, including the fixed local context.
    pub reference_point: QueryPoint,
    /// Selected positive source answer whose claim stays pinned on this face.
    pub resistance: InterfaceResistance,
}

impl<'a> TemperatureContact<'a> {
    /// Assign one original reference query to every A face of a planar
    /// nonmatching contact. The opposite tessellation does not resample the
    /// fixed context map. Each overlap point updates its own two trace
    /// temperatures while retaining the owning A face's pressure and other
    /// coordinates, card identity, and individually selected source claim.
    ///
    /// Entries may arrive in any order and must cover A exactly once. Distinct
    /// interface cards still require distinct physical interfaces. Contexts
    /// are held through a thermal solve; this does not evolve mechanical
    /// pressure, contact geometry, finish or history. Rebinding changed context
    /// remains the responsibility of its physical owner.
    ///
    /// The original pair-test allowance bounds preparation before any local
    /// source queries. Cancellation is checked between source admissions;
    /// common-refinement geometry is prepared once by `with_nonmatching`.
    #[allow(clippy::too_many_arguments)]
    pub fn new_nonmatching_mapped(
        cx: &Cx<'_>,
        name: impl Into<String>,
        mut side_a: Vec<usize>,
        side_b: Vec<usize>,
        card: &'a InterfaceSystemCard,
        mut face_points: Vec<(usize, QueryPoint)>,
        coordinates: TemperatureContactCoordinates,
        selection: SelectionPolicy,
        options: NonmatchingOptions,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        let name = name.into();
        let count = side_a
            .len()
            .checked_add(side_b.len())
            .ok_or_else(|| invalid(&name, "mapped contact face count overflows"))?;
        let tests = count
            .checked_mul(count.saturating_sub(1))
            .map(|n| n / 2)
            .ok_or_else(|| invalid(&name, "mapped contact pair allowance overflows"))?;
        if tests > options.max_pair_tests {
            return Err(invalid(
                &name,
                "mapped contact exceeds the pair allowance before source queries",
            ));
        }
        if side_a.is_empty() || side_b.is_empty() || face_points.len() != side_a.len() {
            return Err(invalid(
                &name,
                "source contexts must name every A face exactly once and no other face",
            ));
        }
        side_a.sort_unstable();
        face_points.sort_by_key(|entry| entry.0);
        if side_a.windows(2).any(|pair| pair[0] == pair[1])
            || face_points
                .iter()
                .zip(&side_a)
                .any(|((face, _), expected)| face != expected)
        {
            return Err(invalid(
                &name,
                "source contexts must name every A face exactly once and no other face",
            ));
        }
        let mut references = Vec::new();
        references
            .try_reserve_exact(face_points.len())
            .map_err(|_| invalid(&name, "mapped contact source allocation refused"))?;
        let mut resistances = Vec::new();
        resistances
            .try_reserve_exact(face_points.len())
            .map_err(|_| invalid(&name, "mapped contact reference allocation refused"))?;
        for (index, (side_a, reference_point)) in face_points.into_iter().enumerate() {
            poll(cx, index)?;
            Self::admit_reference(&name, &reference_point, &coordinates)?;
            let resistance =
                InterfaceResistance::from_card(&name, card, &reference_point, selection)?;
            resistances.push((side_a, resistance.clone()));
            references.push(ContactFaceReference {
                side_a,
                reference_point,
                resistance,
            });
        }
        let surface = NonmatchingSurface::new_mapped(name, side_a, side_b, resistances, options)?;
        let first = &references[0];
        let result = Self {
            card,
            reference_point: first.reference_point.clone(),
            resistance: first.resistance.clone(),
            coordinates,
            reference: ContactSurface::Nonmatching(surface),
            face_references: references,
        };
        poll(cx, result.face_references.len())?;
        Ok(result)
    }

    /// Canonical A-face source contexts. Empty for a uniform-context law.
    /// These are original queries, not averages of later quadrature states.
    #[must_use]
    pub fn reference_face_states(&self) -> &[ContactFaceReference] {
        &self.face_references
    }

    pub(super) fn source_context(
        &self,
        side_a: Option<usize>,
    ) -> Result<(&QueryPoint, &InterfaceResistance), ConductionError> {
        if self.face_references.is_empty() {
            return Ok((&self.reference_point, &self.resistance));
        }
        let side_a = side_a
            .ok_or_else(|| invalid(self.name(), "mapped contact query lost its original A face"))?;
        let index = self
            .face_references
            .binary_search_by_key(&side_a, |entry| entry.side_a)
            .map_err(|_| invalid(self.name(), "mapped contact query names an unbound A face"))?;
        let source = &self.face_references[index];
        Ok((&source.reference_point, &source.resistance))
    }
}
