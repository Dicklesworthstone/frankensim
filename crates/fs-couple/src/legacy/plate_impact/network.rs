//! Simultaneous normal contact networks of original spatial plate trajectories.
//!
//! An edge carries opposite impulses at one shared world-x/y point. Its signed
//! incidence assembles the complete contact mobility from each original plate's
//! cached consistent-mass response, including contacts on different interfaces.

use super::multiple::{
    PlateImpactSite, PlateRigidImpactSet, PlateRigidImpactSetBudget, PlateRigidImpactSetControls,
};
use super::plates::{PlateImpactSurface, difference, requests};
use super::{PlateRigidImpactError, PlateTransferAudit, check, dd, poll, vector};
use fs_contact::multicontact::{MultiImpactResult, NormalContactImpact, solve_normal_impacts};
use fs_math::dd::Dd;
use fs_mbd::Vec3;
use fs_plate_transient::PlateLoad;
use fs_time::galpha::SecondOrderState;
use fs_time::galpha::impulse::SecondOrderImpulse;

/// Positive impulse acts along world +z on `positive_plate` and -z on
/// `negative_plate`. Both endpoints are located at `site.point_xy` on their
/// original meshes; each surface declares its own reference face height.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateNetworkContact {
    pub positive_plate: usize,
    pub negative_plate: usize,
    pub site: PlateImpactSite,
}

/// Admission before owner storage, point rows or the shared matrix allocate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PlateImpactNetworkBudget {
    /// At most 128 original owners (two endpoints per admitted contact).
    pub max_plates: usize,
    pub max_contacts: usize,
    /// Sum of original full nodal dimensions, including unused owners.
    pub max_total_full_dofs: usize,
    /// Includes the referenced owners' admitted workspace ceilings, raw nodal
    /// requests, graph rows, staged updates and receipts. Borrowed trajectories
    /// and the contact owner's separately bounded solver scratch are excluded.
    pub max_workspace_components: usize,
}

/// One original plate's actual momentum change and complete support reactions.
/// A redundant law supplies one checked impulse representative; these values
/// do not imply a unique local pressure or thermal partition.
#[derive(Debug)]
pub struct PlateNetworkOwnerBalance {
    pub support_impulses: Vec<(usize, f64)>,
    pub plate_linear_impulse_world: Vec3,
    pub plate_angular_impulse_world: Vec3,
    pub support_linear_impulse_world: Vec3,
    pub support_angular_impulse_world: Vec3,
    pub kinetic_energy_change_j: f64,
}

/// Collective balance from actual native state changes, accumulated before
/// rounding any plate's contribution. Internal edge impulses cancel in world
/// transverse Pz and bending Lx/Ly; support impulses remain external.
#[derive(Debug)]
pub struct PlateNetworkImpulseBalance {
    /// Original owner order; absent when the complete nodal request is zero.
    pub plates: Vec<Option<PlateNetworkOwnerBalance>>,
    pub linear_balance_residual: Vec3,
    pub angular_balance_residual: Vec3,
    pub kinetic_energy_change_j: f64,
    pub relative_impulse_work_j: f64,
    pub work_energy_residual_j: f64,
    pub energy_residual_j: f64,
}

/// Returned after every affected continuing trajectory has been published.
#[derive(Debug)]
pub struct PlateImpactNetworkReport {
    pub contact: MultiImpactResult,
    pub gaps_m: Vec<f64>,
    pub relative_normal_speeds_before_m_s: Vec<f64>,
    pub relative_normal_speeds_after_m_s: Vec<f64>,
    /// Original owner order. Zero nodal requests never restart an owner.
    pub plates: Vec<Option<SecondOrderImpulse>>,
    pub balance: Option<PlateNetworkImpulseBalance>,
}

struct Incident {
    edge: usize,
    sign: f64,
}

struct Owner<'a, L: ?Sized> {
    prepared: PlateRigidImpactSet<'a, L>,
    incidents: Vec<Incident>,
}

struct Edge {
    contact: PlateNetworkContact,
    positive_row: usize,
    negative_row: usize,
    reference_gap: Dd,
}

/// Fixed-normal event graph on independently meshed, supported and damped
/// plates. Use `L = dyn PlateLoad` for heterogeneous original load owners.
/// Preparation caches each original point response; events reuse the complete
/// signed matrix and perform at most one original impulse solve per plate.
pub struct PlateImpactNetwork<'a, L: ?Sized> {
    owners: Vec<Option<Owner<'a, L>>>,
    edges: Vec<Edge>,
    delassus: Vec<f64>,
    controls: PlateRigidImpactSetControls,
}

impl<'a, L: PlateLoad + ?Sized> PlateImpactNetwork<'a, L> {
    /// Each original owner keeps its mass, support map, damping, loads and
    /// restart path. Budgets follow `surfaces` order. Unreferenced owners need
    /// no point cache and their state is not admitted, copied or restarted.
    pub fn prepare(
        surfaces: &[PlateImpactSurface<'a, L>],
        contacts: &[PlateNetworkContact],
        controls: PlateRigidImpactSetControls,
        owner_budgets: &[PlateRigidImpactSetBudget],
        budget: PlateImpactNetworkBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PlateRigidImpactError> {
        poll(cancelled)?;
        let p = surfaces.len();
        let k = contacts.len();
        if p < 2 || p > 128 || p > budget.max_plates {
            return Err(PlateRigidImpactError::Budget("plate-network owner count"));
        }
        if k == 0 || k > 64 || k > budget.max_contacts || k > controls.contact.solver.max_variables
        {
            return Err(PlateRigidImpactError::Budget("plate-network contact count"));
        }
        if owner_budgets.len() != p {
            return Err(PlateRigidImpactError::Invalid(
                "one budget per original plate",
            ));
        }
        for contact in contacts {
            poll(cancelled)?;
            if contact.positive_plate >= p
                || contact.negative_plate >= p
                || contact.positive_plate == contact.negative_plate
            {
                return Err(PlateRigidImpactError::Invalid(
                    "distinct original plate endpoints",
                ));
            }
            let a = surfaces[contact.positive_plate].reference_height_m;
            let b = surfaces[contact.negative_plate].reference_height_m;
            if !a.is_finite() || !b.is_finite() || !(dd(a) - dd(b)).to_f64().is_finite() {
                return Err(PlateRigidImpactError::Invalid(
                    "finite plate face separation",
                ));
            }
        }
        // No owner-sized arrays are needed to admit the whole graph. At most
        // 64 edges are scanned per owner before any graph allocation occurs.
        let mut full_dofs = 0usize;
        let mut owner_workspace = 0usize;
        for (i, surface) in surfaces.iter().enumerate() {
            poll(cancelled)?;
            full_dofs = full_dofs
                .checked_add(surface.dynamics.model().dof_map.len())
                .ok_or(PlateRigidImpactError::Budget("plate-network DOF overflow"))?;
            let count = contacts
                .iter()
                .filter(|c| c.positive_plate == i || c.negative_plate == i)
                .count();
            if count > owner_budgets[i].max_contacts {
                return Err(PlateRigidImpactError::Budget(
                    "plate-network owner contacts",
                ));
            }
            if count
                .checked_mul(count)
                .and_then(|n| n.checked_mul(surface.dynamics.model().free))
                .is_none_or(|n| n > owner_budgets[i].max_mobility_products)
            {
                return Err(PlateRigidImpactError::Budget(
                    "plate-network owner mobility products",
                ));
            }
            if count != 0 {
                owner_workspace = owner_workspace
                    .checked_add(owner_budgets[i].plate.max_workspace_components)
                    .ok_or(PlateRigidImpactError::Budget(
                        "plate-network workspace overflow",
                    ))?;
            }
        }
        if full_dofs > budget.max_total_full_dofs {
            return Err(PlateRigidImpactError::Budget(
                "plate-network total full DOFs",
            ));
        }
        let workspace = (|| {
            owner_workspace
                .checked_add(full_dofs.checked_mul(4)?)?
                .checked_add(k.checked_mul(k)?.checked_mul(2)?)?
                .checked_add(k.checked_mul(48)?)?
                .checked_add(p.checked_mul(40)?)
        })()
        .ok_or(PlateRigidImpactError::Budget(
            "plate-network workspace overflow",
        ))?;
        if workspace > budget.max_workspace_components {
            return Err(PlateRigidImpactError::Budget(
                "plate-network live workspace",
            ));
        }
        let mut edges: Vec<_> = contacts
            .iter()
            .map(|&contact| Edge {
                reference_gap: dd(surfaces[contact.positive_plate].reference_height_m)
                    - dd(surfaces[contact.negative_plate].reference_height_m),
                contact,
                positive_row: 0,
                negative_row: 0,
            })
            .collect();
        let mut owners = Vec::with_capacity(p);
        let mut matrix = vec![Dd::ZERO; k * k];
        for (i, surface) in surfaces.iter().enumerate() {
            poll(cancelled)?;
            let mut incidents = Vec::new();
            let mut sites = Vec::new();
            for (j, edge) in edges.iter_mut().enumerate() {
                let c = edge.contact;
                let sign = if c.positive_plate == i {
                    edge.positive_row = sites.len();
                    1.0
                } else if c.negative_plate == i {
                    edge.negative_row = sites.len();
                    -1.0
                } else {
                    continue;
                };
                incidents.push(Incident { edge: j, sign });
                sites.push(c.site);
            }
            if sites.is_empty() {
                owners.push(None);
                continue;
            }
            let prepared = PlateRigidImpactSet::prepare(
                surface.dynamics,
                surface.mesh,
                surface.full_mass,
                &sites,
                controls,
                owner_budgets[i],
                cancelled,
            )?;
            let count = sites.len();
            for (a, source) in incidents.iter().enumerate() {
                poll(cancelled)?;
                for (b, target) in incidents.iter().enumerate() {
                    let slot = k * source.edge + target.edge;
                    matrix[slot] = matrix[slot]
                        + dd(source.sign * target.sign)
                            * dd(prepared.plate_cross_mobility()[count * a + b]);
                }
            }
            owners.push(Some(Owner {
                prepared,
                incidents,
            }));
        }
        let delassus: Vec<_> = matrix.into_iter().map(Dd::to_f64).collect();
        if delassus.iter().any(|v| !v.is_finite()) {
            return Err(PlateRigidImpactError::Invalid(
                "representable plate-network mobility",
            ));
        }
        poll(cancelled)?;
        Ok(Self {
            owners,
            edges,
            delassus,
            controls,
        })
    }

    /// Cached full row-major B M^-1 B^T across all original owners [kg^-1].
    #[must_use]
    pub fn contact_mobility(&self) -> &[f64] {
        &self.delassus
    }

    /// Apply one global admissible contact law at a shared native physical time.
    /// All referenced owners must share that time, including inactive edges.
    /// Unreferenced owners and exactly zero complete nodal requests retain all
    /// bits. Every other owner publishes only native velocity and acceleration,
    /// after all actual endpoint, momentum, support and kinetic-work checks.
    ///
    /// This event API does not discover collisions, advance persistent contact,
    /// rotate normals, supply friction, or turn mechanical loss into heat.
    pub fn apply(
        &self,
        states: &mut [SecondOrderState],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlateImpactNetworkReport, PlateRigidImpactError> {
        poll(cancelled)?;
        if states.len() != self.owners.len() {
            return Err(PlateRigidImpactError::Invalid(
                "one state per original plate",
            ));
        }
        let mut time = None;
        for (i, owner) in self.owners.iter().enumerate() {
            poll(cancelled)?;
            if let Some(owner) = owner {
                owner.prepared.first.admit_state(&states[i])?;
                if time.is_some_and(|t| t != states[i].t) {
                    return Err(PlateRigidImpactError::Invalid(
                        "plate-network shared physical time",
                    ));
                }
                time = Some(states[i].t);
            }
        }
        let k = self.edges.len();
        let mut contacts = Vec::with_capacity(k);
        for edge in &self.edges {
            poll(cancelled)?;
            let c = edge.contact;
            let positive = self
                .owner(c.positive_plate)
                .prepared
                .row(edge.positive_row)
                .1;
            let negative = self
                .owner(c.negative_plate)
                .prepared
                .row(edge.negative_row)
                .1;
            contacts.push(NormalContactImpact {
                gap: (edge.reference_gap
                    + difference(
                        positive,
                        &states[c.positive_plate].q,
                        negative,
                        &states[c.negative_plate].q,
                        cancelled,
                    )?)
                .to_f64(),
                activation_tolerance: c.site.activation_tolerance_m,
                relative_normal_speed: difference(
                    positive,
                    &states[c.positive_plate].v,
                    negative,
                    &states[c.negative_plate].v,
                    cancelled,
                )?
                .to_f64(),
                restitution: c.site.restitution,
            });
        }
        let contact = solve_normal_impacts(
            &contacts,
            &self.delassus,
            self.controls.contact,
            &mut *cancelled,
        )
        .map_err(PlateRigidImpactError::MultipleContact)?;
        let gaps = contacts.iter().map(|c| c.gap).collect();
        let before: Vec<_> = contacts.iter().map(|c| c.relative_normal_speed).collect();
        let mut receipts: Vec<Option<SecondOrderImpulse>> =
            std::iter::repeat_with(|| None).take(states.len()).collect();
        if contact.normal_impulses.iter().all(|&j| j == 0.0) {
            poll(cancelled)?;
            return Ok(PlateImpactNetworkReport {
                contact,
                gaps_m: gaps,
                relative_normal_speeds_after_m_s: before.clone(),
                relative_normal_speeds_before_m_s: before,
                plates: receipts,
                balance: None,
            });
        }
        let mut trials: Vec<Option<SecondOrderState>> =
            std::iter::repeat_with(|| None).take(states.len()).collect();
        let mut transfers = Vec::with_capacity(states.len());
        for (i, owner) in self.owners.iter().enumerate() {
            poll(cancelled)?;
            let Some(owner) = owner else {
                transfers.push(None);
                continue;
            };
            let impulses: Vec<_> = owner
                .incidents
                .iter()
                .map(|row| row.sign * contact.normal_impulses[row.edge])
                .collect();
            if impulses.iter().all(|&j| j == 0.0) {
                transfers.push(None);
                continue;
            }
            let (full_request, reduced_request) =
                requests(&owner.prepared, &impulses, 1.0, cancelled)?;
            if full_request.iter().all(|&j| j == Dd::ZERO) {
                transfers.push(None);
                continue;
            }
            let (trial, receipt) = owner.prepared.first.stage_plate_impulse(
                &states[i],
                &reduced_request,
                cancelled,
            )?;
            let transfer = owner.prepared.first.account_plate_transfer(
                &states[i],
                &trial,
                &receipt,
                |r| full_request[r],
                cancelled,
            )?;
            trials[i] = Some(trial);
            receipts[i] = Some(receipt);
            transfers.push(Some(transfer));
        }
        let mut after = Vec::with_capacity(k);
        let mut work = Dd::ZERO;
        for (i, edge) in self.edges.iter().enumerate() {
            poll(cancelled)?;
            let c = edge.contact;
            let positive = trials[c.positive_plate]
                .as_ref()
                .unwrap_or(&states[c.positive_plate]);
            let negative = trials[c.negative_plate]
                .as_ref()
                .unwrap_or(&states[c.negative_plate]);
            let speed = difference(
                self.owner(c.positive_plate)
                    .prepared
                    .row(edge.positive_row)
                    .1,
                &positive.v,
                self.owner(c.negative_plate)
                    .prepared
                    .row(edge.negative_row)
                    .1,
                &negative.v,
                cancelled,
            )?
            .to_f64();
            check(
                "plate-network endpoint speed",
                speed - contact.relative_normal_speeds_after[i],
                self.controls.speed_tolerance_m_s
                    + self.controls.impulse.relative_impulse_tolerance
                        * before[i].abs().max(speed.abs()),
            )?;
            work = work + dd(contact.normal_impulses[i]) * (dd(0.5 * before[i]) + dd(0.5 * speed));
            after.push(speed);
        }
        check(
            "original plate-network contact work",
            (work - dd(contact.relative_work)).to_f64(),
            self.controls.energy_tolerance_j,
        )?;
        let balance = self.account_network(transfers, &contact, work)?;
        poll(cancelled)?;
        for (state, trial) in states.iter_mut().zip(trials) {
            if let Some(trial) = trial {
                state.v = trial.v;
                state.a = trial.a;
            }
        }
        Ok(PlateImpactNetworkReport {
            contact,
            gaps_m: gaps,
            relative_normal_speeds_before_m_s: before,
            relative_normal_speeds_after_m_s: after,
            plates: receipts,
            balance: Some(balance),
        })
    }

    fn owner(&self, index: usize) -> &Owner<'a, L> {
        self.owners[index]
            .as_ref()
            .expect("prepared edge references its original plate")
    }

    fn account_network(
        &self,
        transfers: Vec<Option<PlateTransferAudit>>,
        contact: &MultiImpactResult,
        work: Dd,
    ) -> Result<PlateNetworkImpulseBalance, PlateRigidImpactError> {
        let mut linear = [Dd::ZERO; 3];
        let mut angular = [Dd::ZERO; 3];
        let mut energy = Dd::ZERO;
        for transfer in transfers.iter().flatten() {
            for c in 0..3 {
                linear[c] = linear[c] + transfer.linear[c] - transfer.support_linear[c];
                angular[c] = angular[c] + transfer.angular[c] - transfer.support_angular[c];
            }
            energy = energy + transfer.kinetic_energy;
        }
        let scale = contact
            .normal_impulses
            .iter()
            .fold(Dd::ZERO, |sum, &j| sum + dd(j.abs()))
            .to_f64();
        for c in 0..3 {
            check(
                "plate-network linear balance",
                linear[c].to_f64(),
                self.controls.linear_impulse_tolerance_ns
                    + self.controls.impulse.relative_impulse_tolerance * scale,
            )?;
            check(
                "plate-network angular balance",
                angular[c].to_f64(),
                self.controls.angular_impulse_tolerance_nms,
            )?;
        }
        let work_residual = (energy - work).to_f64();
        let energy_residual = (energy + dd(contact.dissipated_energy)).to_f64();
        check(
            "plate-network impact work",
            work_residual,
            self.controls.energy_tolerance_j,
        )?;
        check(
            "plate-network impact energy",
            energy_residual,
            self.controls.energy_tolerance_j,
        )?;
        Ok(PlateNetworkImpulseBalance {
            plates: transfers
                .into_iter()
                .map(|transfer| {
                    transfer.map(|t| PlateNetworkOwnerBalance {
                        support_impulses: t.support_impulses,
                        plate_linear_impulse_world: vector(t.linear.map(Dd::to_f64)),
                        plate_angular_impulse_world: vector(t.angular.map(Dd::to_f64)),
                        support_linear_impulse_world: vector(t.support_linear.map(Dd::to_f64)),
                        support_angular_impulse_world: vector(t.support_angular.map(Dd::to_f64)),
                        kinetic_energy_change_j: t.kinetic_energy.to_f64(),
                    })
                })
                .collect(),
            linear_balance_residual: vector(linear.map(Dd::to_f64)),
            angular_balance_residual: vector(angular.map(Dd::to_f64)),
            kinetic_energy_change_j: energy.to_f64(),
            relative_impulse_work_j: work.to_f64(),
            work_energy_residual_j: work_residual,
            energy_residual_j: energy_residual,
        })
    }
}
