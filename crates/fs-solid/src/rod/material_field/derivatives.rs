//! Source-bound local energy derivatives and the actual rod/link interaction graph.
use super::{SpatialMaterialRod, SpatialRodError};
use crate::rod::Rod;
use crate::rod::derivatives::RodEnergy;
use fs_ga::Vec3;

pub(super) struct SpatialEnergy<'a> {
    owner: &'a SpatialMaterialRod,
    neighbors: Vec<Vec<usize>>,
    links: Vec<Vec<(usize, usize)>>,
}

impl<'a> SpatialEnergy<'a> {
    pub(super) fn new(owner: &'a SpatialMaterialRod) -> Self {
        let count = owner.rod.positions.len();
        let mut neighbors: Vec<Vec<usize>> = (0..count).map(|node|
            (node.saturating_sub(1)..(node + 2).min(count)).collect()).collect();
        let mut links = vec![Vec::new(); count];
        for (index, link) in owner.tension_links.iter().enumerate() {
            let attachments = link.attachments();
            let (a, b) = (attachments[0].node, attachments[1].node);
            neighbors[a].push(b);
            neighbors[b].push(a);
            links[a].push((index, 0));
            links[b].push((index, 1));
        }
        for row in &mut neighbors { row.sort_unstable(); row.dedup(); }
        Self { owner, neighbors, links }
    }
}

impl RodEnergy<SpatialRodError> for SpatialEnergy<'_> {
    fn value(&self, rod: &Rod) -> Result<f64, SpatialRodError> {
        self.owner.equilibrium_energy_for(rod)
    }

    fn derivative(
        &self, base: &Rod, plus: &Rod, minus: &Rod,
        node: usize, component: usize, step: f64,
    ) -> Result<f64, SpatialRodError> {
        let mut result = 0.0;
        for segment in node.saturating_sub(1)..(node + 1).min(self.owner.sections.len()) {
            result += (self.owner.segment_energy_for(plus, segment)?
                - self.owner.segment_energy_for(minus, segment)?) / (2.0 * step);
        }
        // These are already native analytical gradients. Subtract the force
        // ON the rod to obtain the gradient of stored/potential energy.
        for &(index, side) in &self.links[node] {
            let state = self.owner.tension_links[index].state_for(base)?;
            result -= if component < 3 { state.force_n[side][component] }
                else { state.moment_body_nm[side][component - 3] };
        }
        if component < 3 {
            if let Some(forces) = self.owner.gravity_nodal_forces_n() {
                result -= forces[node][component];
            }
        } else if let Some(acceleration) = self.owner.gravity_m_s2() {
            for &(index, side) in &self.links[node] {
                let link = &self.owner.tension_links[index];
                let force = acceleration.map(|g| 0.5 * link.reference_mass_kg() * g);
                let body = base.frames[node].inverse()
                    .rotate(Vec3::new(force[0], force[1], force[2]))
                    .map_err(crate::rod::group_error)?;
                let offset = link.attachments()[side].offset_body_m;
                let moment = Vec3::new(offset[0], offset[1], offset[2]).cross(body);
                result -= [moment.x, moment.y, moment.z][component - 3];
            }
        }
        if !result.is_finite() {
            return Err(SpatialRodError::Invalid("local rod/link/gravity derivative is unrepresentable"));
        }
        Ok(result)
    }

    fn coupled(&self, row_node: usize, column_node: usize) -> bool {
        self.neighbors[row_node].binary_search(&column_node).is_ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::supports::RodSupport;
    use super::super::tension::{MaterialTensionLink, TensionLinkAttachment};
    use crate::rod::{TipLoad, derivatives::perturb};
    use fs_evidence::ValidityDomain;
    use fs_exec::CancelGate;
    use fs_ga::So3;
    use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
    use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection};
    use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point};
    use fs_qty::{Density, Pressure};

    fn section() -> ConcentricShearSection {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, value, dims) in [("density", 1000.0, Density::DIMS),
            ("young_modulus", 1e6, Pressure::DIMS), ("shear_modulus", 4e5, Pressure::DIMS)] {
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(), uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
                provenance: Provenance { source: "authored local rod derivative fixture".into(),
                    license: "test-only".into(), artifact: None },
            }).unwrap();
            requirements.push(ScalarPropertyRequirement::try_new(name, dims,
                ScalarAdmissibility::StrictlyPositive).unwrap());
        }
        let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(),
            phase: "solid".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
        let material = resolve_material_state_point(&card, &QueryPoint::new(), &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
            outer_radius_m: 0.01, material,
        }]).unwrap().with_shear_modulus("shear_modulus", 0.8).unwrap()
    }

    fn linked(offset: [f64; 3], natural: f64) -> SpatialMaterialRod {
        let section = section();
        let rod = SpatialMaterialRod::try_straight(1.0, vec![section.clone(); 4]).unwrap();
        let link = MaterialTensionLink::try_new(section.elastic().clone(), [
            TensionLinkAttachment { node: 1, offset_body_m: offset },
            TensionLinkAttachment { node: 4, offset_body_m: offset },
        ], natural, 1, 0.2).unwrap();
        rod.with_tension_links(vec![link], &CancelGate::new()).unwrap()
    }

    #[test]
    fn nonlocal_link_block_matches_analytical_axial_stiffness() {
        let rod = linked([0.0; 3], 0.7);
        let model = SpatialEnergy::new(&rod);
        assert!(model.coupled(1, 4));
        assert!(model.coupled(4, 1));
        assert!(!model.coupled(2, 4));
        let mut plus = rod.rod.clone();
        let mut minus = rod.rod.clone();
        perturb(&mut plus, 4, 0, 1e-4).unwrap();
        perturb(&mut minus, 4, 0, -1e-4).unwrap();
        let rp = plus.residual_on(&TipLoad::default(), 0.0, &model, &|| Ok(()),
            0.25, &rod.free_dofs, Some(18)).unwrap();
        let rm = minus.residual_on(&TipLoad::default(), 0.0, &model, &|| Ok(()),
            0.25, &rod.free_dofs, Some(18)).unwrap();
        let expected = -rod.sections[0].elastic().axial_stiffness_n() / 0.7;
        assert!(((rp[0] - rm[0]) / 2e-4 - expected).abs() < 1e-7);
    }

    #[test]
    fn nonlocal_tension_link_changes_actual_equilibrium_compliance() {
        let gate = CancelGate::new();
        let mut rod = linked([0.0; 3], 0.75);
        let supports: Vec<_> = (1..=4).map(|node| RodSupport {
            node, position_m: [None, Some(0.0), Some(0.0)], frame: Some(So3::identity()),
        }).collect();
        rod = rod.with_supports(&supports, &gate).unwrap();
        let ea = rod.sections[0].elastic().axial_stiffness_n();
        let load = TipLoad { force: [0.2, 0.0, 0.0], moment: [0.0; 3] };
        rod.solve_static(&load, 2, 1e-8, &gate).unwrap();
        let u1 = load.force[0] * 0.25 / ea;
        let span = load.force[0] / (ea / 0.75 + ea / 0.75);
        for node in 1..=4 {
            let expected = node as f64 * 0.25 + u1 + (node - 1) as f64 * span / 3.0;
            assert!((rod.positions_m()[node][0] - expected).abs() < 1e-9);
        }
        let before = rod.positions_m().to_vec();
        gate.request();
        assert!(rod.solve_static(&load, 1, 1e-8, &gate).is_err());
        assert_eq!(rod.positions_m(), before.as_slice());
    }

    #[test]
    fn analytical_link_and_rotating_weight_gradients_match_complete_potential() {
        let gate = CancelGate::new();
        let rod = linked([0.0, 0.02, 0.01], 0.7)
            .with_gravity([0.0, -9.0, 1.0], &gate).unwrap();
        let mut base = rod.rod.clone();
        perturb(&mut base, 4, 5, 0.1).unwrap();
        base.positions[4][1] += 0.01;
        let model = SpatialEnergy::new(&rod);
        let local = base.residual(&TipLoad::default(), 0.0, &model, &|| Ok(()), 0.25).unwrap();
        let full = base.residual(&TipLoad::default(), 0.0,
            &|r: &Rod| rod.equilibrium_energy_for(r), &|| Ok(()), 0.25).unwrap();
        for (a, b) in local.iter().zip(full) {
            assert!((a - b).abs() < 2e-6 * a.abs().max(1.0), "{a:e} != {b:e}");
        }
    }

    #[test]
    fn slack_links_keep_their_possible_nonlocal_tangent_blocks() {
        let rod = linked([0.0; 3], 0.8);
        assert_eq!(rod.tension_link_states().unwrap()[0].tension_n, 0.0);
        assert!(SpatialEnergy::new(&rod).coupled(1, 4));
    }
}
