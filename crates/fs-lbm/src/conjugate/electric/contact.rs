//! Resistive, zero-thickness electrical interfaces in the existing voxel solver.

use std::collections::{BTreeMap, BTreeSet};

use super::{
    ChtError, ElectricSetup, ElectricSolution, Node, Voxel, VoxelDomain, admit,
    finite_positive, poll, poll_every, solve_electric_field, solve_electric_impl,
    solve_electric_field_with_contact_field,
};
use fs_exec::CancelGate;

/// One explicitly located ohmic contact between adjacent voxel cells.
///
/// `resistance_ohm_m2 / face_area` is the additional face resistance. It is
/// distinct from either bulk resistivity (ohm m), and is added once in series.
/// The ordered cells define the heat-partition convention, not current direction.
/// Reversing current must not reverse where resistive heat is deposited.
///
/// The interface has no heat capacity. Its generated heat is assigned to the
/// two existing thermal cells by an explicit partition; no equal split or
/// effusivity model is inferred. Mechanical contact area, pressure dependence,
/// dielectric breakdown, tunnelling and contact dynamics are caller-owned.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElectricalContact {
    /// Two face-adjacent cells, in the caller's surface order.
    pub cells: [usize; 2],
    /// Area-specific contact resistance, ohm m²; finite and non-negative.
    pub resistance_ohm_m2: f64,
    /// Fraction of interface-generated heat assigned to `cells[0]`, in `[0,1]`.
    /// The remaining fraction goes to `cells[1]`. Bulk heat remains local.
    pub heat_fraction_to_first: f64,
}

impl ElectricSetup {
    /// Solve bulk current and contact heating with the shared spatial kernel.
    ///
    /// The bulk callback has the same semantics as [`Self::solve_with_resistivity`].
    /// Contacts are resolved at the supplied temperature state by the caller;
    /// they are not cached across temperature updates. Omitted faces are perfect
    /// electrical contacts, exactly as in the original entry point. A contact
    /// may join two bulk conductor cells or a conductor and an ideal electrode;
    /// its heat may then be deposited in the electrode's thermal cell as well.
    /// It cannot turn an insulating voxel into a conductor. Directly touching
    /// distinct ideal electrodes remain refused by the ordinary topology owner.
    ///
    /// The returned `joule_w` includes both bulk and contact heat exactly once.
    /// Use it directly as thermal power; do not add another `I²R_contact` source.
    /// No material or interface source authority is inferred from raw numbers.
    ///
    /// # Errors
    /// Non-adjacent, duplicated or nonconducting contact faces, invalid contact
    /// data, any ordinary material/electrical refusal, or cancellation. Both
    /// contacts and the setup remain unchanged on failure.
    pub fn solve_with_contacts(
        &self,
        domain: &VoxelDomain,
        temperature: &[f64],
        contacts: &[ElectricalContact],
        mut resolve: impl FnMut(usize, u16, f64) -> Result<f64, ChtError>,
        gate: &CancelGate,
    ) -> Result<ElectricSolution, ChtError> {
        solve_electric_impl(
            domain, self, Some(temperature),
            |cell, material, point, _| {
                resolve(cell, material, point.expect("explicit temperature field"))
            },
            contacts, gate,
        )
    }

    /// Material-assigned fluid/solid conduction with explicit resistive faces.
    ///
    /// This is the contact-capable entry point for [`Self::solve_with_cell_resistivity`].
    /// `Some(rho)` declares conduction; `None` declares insulation independently
    /// of mechanical occupancy. Missing source data must return an error. The
    /// solid preset table is ignored. Electrode cells remain ideal electrical
    /// equipotentials, but may receive contact heat in the thermal problem.
    /// The complete non-electrode temperature field is checked before callbacks.
    ///
    /// # Errors
    /// The refusals of [`Self::solve_with_contacts`] and the per-cell material
    /// path. A contact involving an explicitly insulated cell is refused.
    pub fn solve_with_cell_contacts(
        &self,
        domain: &VoxelDomain,
        temperature: &[f64],
        contacts: &[ElectricalContact],
        resolve: impl FnMut(usize, Voxel, f64) -> Result<Option<f64>, ChtError>,
        gate: &CancelGate,
    ) -> Result<ElectricSolution, ChtError> {
        let (nodes, rho) = resolve_cells(self, domain, temperature, resolve, gate)?;
        solve_electric_field(domain, self, nodes, &rho, contacts, gate)
    }

    /// Solve with explicit open electrical faces as well as resistive contacts.
    ///
    /// An opening removes ONLY that face from electrical reachability, current
    /// balance and Joule deposition. It does not remove either bulk material,
    /// change its heat capacity, or open a thermal/mechanical contact. A face
    /// in neither list remains perfectly connected. A face in both lists or
    /// repeated in either orientation refuses; no ordering can silently win.
    /// Zero resistance means a closed contact, never an opening. Infinity is
    /// not used as a large-resistance approximation.
    ///
    /// Bulk laws are still queried for floating conductors. Their unreferenced
    /// potentials are NaN and stationary current/heat are zero. A disconnected
    /// nonzero ideal current drive refuses; a disconnected zero-current drive
    /// has undefined voltage and zero delivered power. Capacitive transients,
    /// arc/breakdown and motional isolated-current gauges are outside this API.
    /// Distinct touching ideal electrodes retain the original admission refusal.
    ///
    /// # Errors
    /// Existing material/solver refusals, invalid/conflicting face topology,
    /// unrepresentable values, or cancellation. No setup or body is mutated.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_cell_contact_topology(
        &self, domain: &VoxelDomain, temperature: &[f64],
        contacts: &[ElectricalContact], open_faces: &[[usize; 2]],
        resolve: impl FnMut(usize, Voxel, f64) -> Result<Option<f64>, ChtError>,
        gate: &CancelGate,
    ) -> Result<ElectricSolution, ChtError> {
        let (nodes, rho) = resolve_cells(self, domain, temperature, resolve, gate)?;
        let field = ContactField::new(domain, &nodes, contacts, gate)?
            .with_open_faces(domain, open_faces, gate)?;
        solve_electric_field_with_contact_field(domain, self, nodes, &rho, field, None, gate)
    }
}

// One material/topology admission for stationary and moving conductors.
pub(super) fn resolve_cells(
    setup: &ElectricSetup,
    domain: &VoxelDomain,
    temperature: &[f64],
    resolve: impl FnMut(usize, Voxel, f64) -> Result<Option<f64>, ChtError>,
    gate: &CancelGate,
) -> Result<(Vec<Node>, Vec<f64>), ChtError> {
    poll(gate)?;
    let nodes = admit(domain, setup, false)?;
    resolve_admitted_cells(domain, temperature, nodes, resolve, gate)
}

// Numerical gauges never masquerade as electrode cells: every isolated bulk
// cell reaches this same temperature and material evaluation as a wired cell.
pub(super) fn resolve_admitted_cells(
    domain: &VoxelDomain,
    temperature: &[f64],
    mut nodes: Vec<Node>,
    mut resolve: impl FnMut(usize, Voxel, f64) -> Result<Option<f64>, ChtError>,
    gate: &CancelGate,
) -> Result<(Vec<Node>, Vec<f64>), ChtError> {
    poll(gate)?;
    let cells = domain.cell_count();
    if temperature.len() != cells {
        return Err(ChtError::InvalidInput {
            field: "electric.temperature",
            reason: "needs one temperature per cell".into(),
        });
    }
    for c in 0..cells {
        poll_every(gate, c)?;
        if !matches!(nodes[c], Node::Electrode(_)) {
            finite_positive("electric.temperature", temperature[c])?;
        }
    }
    let mut rho = vec![f64::NAN; cells];
    for c in 0..cells {
        poll_every(gate, c)?;
        if matches!(nodes[c], Node::Electrode(_)) { continue; }
        if let Some(value) = resolve(c, domain.voxel_at(c), temperature[c])? {
            finite_positive("electric.resistivity_ohm_m", value)?;
            rho[c] = value;
            nodes[c] = Node::Conductor;
        }
    }
    Ok((nodes, rho))
}

/// Sparse canonical face data. Sorting the key does not discard surface order:
/// the stored heat fraction is transformed to the lower-index cell's role.
pub(super) struct ContactField {
    resistive: BTreeMap<(usize, usize), (f64, f64)>,
    open: BTreeSet<(usize, usize)>,
}

impl ContactField {
    pub(super) fn new(
        domain: &VoxelDomain,
        nodes: &[Node],
        contacts: &[ElectricalContact],
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        let mut field = BTreeMap::new();
        for (index, contact) in contacts.iter().enumerate() {
            if index.is_multiple_of(256) { poll(gate)?; }
            let [a, b] = contact.cells;
            let refuse = |reason: &str| ChtError::InvalidInput {
                field: "electric.contacts",
                reason: format!("contact {index} ({a}, {b}): {reason}"),
            };
            if a >= domain.cell_count() || b >= domain.cell_count()
                || !(0..6).any(|face| domain.neighbor(a, face) == Some(b))
            {
                return Err(refuse("needs two in-domain face-adjacent cells"));
            }
            if !matches!((nodes[a], nodes[b]),
                (Node::Conductor, Node::Conductor | Node::Electrode(_))
                | (Node::Electrode(_), Node::Conductor))
            {
                return Err(refuse("requires a conductor-conductor or conductor-electrode face"));
            }
            let r = contact.resistance_ohm_m2;
            let f = contact.heat_fraction_to_first;
            if !r.is_finite() || r < 0.0 || !f.is_finite() || !(0.0..=1.0).contains(&f) {
                return Err(refuse("requires finite non-negative resistance and a heat fraction in [0,1]"));
            }
            let (key, first_fraction) = if a < b { ((a, b), f) } else { ((b, a), 1.0 - f) };
            if field.insert(key, (r, first_fraction)).is_some() {
                return Err(refuse("the physical face is declared more than once"));
            }
        }
        Ok(Self { resistive: field, open: BTreeSet::new() })
    }

    pub(super) fn with_open_faces(mut self, domain: &VoxelDomain, open_faces: &[[usize; 2]],
        gate: &CancelGate) -> Result<Self, ChtError> {
        if open_faces.len().saturating_add(self.resistive.len())
            > domain.cell_count().saturating_mul(3) {
            return Err(ChtError::InvalidInput {
                field: "electric.open_faces", reason: "face count exceeds the Cartesian domain".into(),
            });
        }
        for (index, &[a, b]) in open_faces.iter().enumerate() {
            if index.is_multiple_of(256) { poll(gate)?; }
            let key = (a.min(b), a.max(b));
            if a >= domain.cell_count() || b >= domain.cell_count()
                || !(0..6).any(|f| domain.neighbor(a, f) == Some(b))
                || self.resistive.contains_key(&key) || !self.open.insert(key) {
                return Err(ChtError::InvalidInput {
                    field: "electric.open_faces",
                    reason: format!("opening {index} ({a}, {b}) needs a unique adjacent face without a resistive assignment"),
                });
            }
        }
        poll(gate)?;
        Ok(self)
    }

    pub(super) fn is_open(&self, a: usize, b: usize) -> bool {
        self.open.contains(&(a.min(b), a.max(b)))
    }

    /// Resistance and heat fraction into `a`, independent of face traversal.
    pub(super) fn at(&self, a: usize, b: usize) -> (f64, f64) {
        let Some(&(r, f)) = self.resistive.get(&(a.min(b), a.max(b))) else { return (0.0, 0.0); };
        (r, if a < b { f } else { 1.0 - f })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conjugate::{Conductor, Electrode, ElectrodeDrive};

    fn bar(current: bool, voltage: f64) -> (VoxelDomain, ElectricSetup) {
        let domain = VoxelDomain::from_fn(4, 1, 1, 1.0, |p| {
            Voxel::Solid(u16::from(p[0] > 2.0))
        }).unwrap();
        let setup = ElectricSetup {
            conductors: vec![Some(Conductor::new(2.0)), Some(Conductor::new(4.0))],
            electrodes: vec![
                Electrode { lo: [0,0,0], hi: [1,1,1], drive: ElectrodeDrive::Voltage(0.0) },
                Electrode { lo: [3,0,0], hi: [4,1,1], drive: if current {
                    ElectrodeDrive::Current(1.0)
                } else { ElectrodeDrive::Voltage(voltage) } },
            ],
        };
        (domain, setup)
    }

    fn solve(domain: &VoxelDomain, setup: &ElectricSetup, contacts: &[ElectricalContact]) -> ElectricSolution {
        setup.solve_with_contacts(domain, &[300.0; 4], contacts,
            |_, material, _| Ok(if material == 0 { 2.0 } else { 4.0 }), &CancelGate::new()).unwrap()
    }

    #[test]
    fn g1_contact_changes_current_and_deposits_only_its_declared_heat_share() {
        for (current, voltage) in [(false, 9.0), (false, -9.0), (true, 9.0)] {
            let (domain, setup) = bar(current, voltage);
            let result = solve(&domain, &setup, &[ElectricalContact {
                cells: [1,2], resistance_ohm_m2: 3.0, heat_fraction_to_first: 0.25,
            }]);
            // R_bulk = 2+4, R_contact = 3, |I| = 1. Heat partition must not
            // change when voltage reverses or when the drive is current-owned.
            assert!((result.electrodes[1].current_a.abs() - 1.0).abs() < 1e-10);
            assert!((result.electrodes[1].voltage_v.abs() - 9.0).abs() < 1e-10);
            assert!((result.joule_w[1] - 2.75).abs() < 1e-10);
            assert!((result.joule_w[2] - 6.25).abs() < 1e-10);
            assert!((result.total_joule_w - 9.0).abs() < 1e-10);
            assert!((result.delivered_w - result.total_joule_w).abs() < 1e-10);
        }
    }

    #[test]
    fn g1_electrode_contact_heat_reaches_both_thermal_cells() {
        let (domain, setup) = bar(false, 9.0);
        let result = solve(&domain, &setup, &[ElectricalContact {
            cells: [0,1], resistance_ohm_m2: 3.0, heat_fraction_to_first: 0.25,
        }]);
        for (heat, expected) in result.joule_w.iter().zip([0.75, 4.25, 4.0, 0.0]) {
            assert!((heat - expected).abs() < 1e-10);
        }
        assert!((result.delivered_w - result.total_joule_w).abs() < 1e-10);
    }

    #[test]
    fn g5_zero_contact_and_reversed_surface_order_preserve_the_same_physics() {
        let (domain, setup) = bar(false, 9.0);
        let original = setup.solve_with_resistivity(&domain, &[300.0; 4],
            |_, m, _| Ok(if m == 0 {2.0} else {4.0}), &CancelGate::new()).unwrap();
        assert_eq!(original, solve(&domain, &setup, &[]));
        assert_eq!(original, solve(&domain, &setup, &[ElectricalContact {
            cells: [1,2], resistance_ohm_m2: 0.0, heat_fraction_to_first: 0.25,
        }]));
        let forward = solve(&domain, &setup, &[ElectricalContact {
            cells: [1,2], resistance_ohm_m2: 3.0, heat_fraction_to_first: 0.25,
        }]);
        let reverse = solve(&domain, &setup, &[ElectricalContact {
            cells: [2,1], resistance_ohm_m2: 3.0, heat_fraction_to_first: 0.75,
        }]);
        assert_eq!(forward, reverse);
    }

    #[test]
    fn g3_material_assigned_fluid_uses_the_same_contact_kernel() {
        let (solid, mut setup) = bar(false, 9.0);
        let contact = [ElectricalContact {cells: [1,2], resistance_ohm_m2: 3.0, heat_fraction_to_first: 0.25}];
        let reference = solve(&solid, &setup, &contact);
        let fluid = VoxelDomain::new(4,1,1,1.0).unwrap();
        setup.conductors.clear();
        let field = setup.solve_with_cell_contacts(&fluid, &[300.0; 4], &contact,
            |cell, occupancy, _| {
                assert_eq!(occupancy, Voxel::Fluid);
                Ok(Some(if cell == 1 {2.0} else {4.0}))
            }, &CancelGate::new()).unwrap();
        assert_eq!(field, reference);
        assert!((0..4).all(|cell| fluid.is_fluid(cell)));
    }

    #[test]
    fn g0_invalid_contact_faces_values_and_cancellation_refuse() {
        let (domain, setup) = bar(false, 9.0);
        let good = ElectricalContact {cells: [1,2], resistance_ohm_m2: 3.0, heat_fraction_to_first: 0.5};
        for contacts in [
            vec![ElectricalContact {cells: [0,2], ..good}],
            vec![ElectricalContact {cells: [1,4], ..good}],
            vec![good, ElectricalContact {cells: [2,1], ..good}],
            vec![ElectricalContact {resistance_ohm_m2: -1.0, ..good}],
            vec![ElectricalContact {resistance_ohm_m2: f64::INFINITY, ..good}],
            vec![ElectricalContact {heat_fraction_to_first: f64::NAN, ..good}],
            vec![ElectricalContact {heat_fraction_to_first: 1.1, ..good}],
        ] {
            assert!(matches!(setup.solve_with_contacts(&domain, &[300.0; 4], &contacts,
                |_, _, _| Ok(2.0), &CancelGate::new()),
                Err(ChtError::InvalidInput {field: "electric.contacts", ..})));
        }
        let gate = CancelGate::new(); gate.request();
        assert_eq!(setup.solve_with_contacts(&domain, &[300.0; 4], &[good],
            |_, _, _| panic!("cancel before material evaluation"), &gate), Err(ChtError::Cancelled));
        let mut insulating = domain.clone(); insulating.set(2,0,0,Voxel::Fluid);
        assert!(matches!(setup.solve_with_contacts(&insulating, &[300.0; 4], &[good],
            |_,_,_| Ok(2.0), &CancelGate::new()), Err(ChtError::InvalidInput {field: "electric.contacts", ..})));
    }
}
