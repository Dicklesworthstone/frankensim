//! Source-bound local energy derivatives and the actual rod/link interaction graph.
use super::{SpatialMaterialRod, SpatialRodError};
use crate::rod::Rod;
use crate::rod::derivatives::RodEnergy;
use fs_ga::Vec3;
use crate::rod::derivatives::analytic::Linearization;

mod analytic;

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

    fn linearize(&self, rod: &Rod, free: &[usize], tangent: bool,
        checkpoint: &impl Fn() -> Result<(), SpatialRodError>)
        -> Result<Option<Linearization>, SpatialRodError> {
        Ok(Some(self.assemble_analytic(rod, free, tangent, checkpoint)?))
    }

    // Retained only as an independent local-stencil oracle. Native equilibrium
    // and support reactions use the complete analytical assembly above.
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

    fn section() -> ConcentricShearSection { section_scaled(1.0) }

    fn section_scaled(scale: f64) -> ConcentricShearSection {
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, value, dims) in [("density", 1000.0, Density::DIMS),
            ("young_modulus", 1e6*scale, Pressure::DIMS), ("shear_modulus", 4e5*scale, Pressure::DIMS)] {
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

    #[test]
    fn root_reaction_keeps_small_weight_beside_remote_stored_energy() {
        let gate = CancelGate::new();
        let mut rod = SpatialMaterialRod::try_straight(1.0, vec![section(); 4]).unwrap()
            .with_gravity([0.0, -1e-6, 0.0], &gate).unwrap();
        rod.rod.positions[4][0] = 1e8;
        assert!(rod.energy_j().unwrap() > 1e17);
        let before = rod.positions_m().to_vec();
        let force = rod.gravity_nodal_forces_n().unwrap()[0][1];
        let reaction = rod.root_support_reaction(&gate).unwrap();
        assert_eq!(reaction.force_n[1], Some(-force));
        assert_eq!(rod.positions_m(), before.as_slice());
        gate.request();
        assert!(rod.root_support_reaction(&gate).is_err());
        assert_eq!(rod.positions_m(), before.as_slice());
    }

    #[test]
    fn root_offset_link_and_weight_reactions_match_direct_force_and_moment() {
        let gate = CancelGate::new();
        let section = section();
        let ea = section.elastic().axial_stiffness_n();
        let offset = [0.0, 0.02, 0.01];
        let link = MaterialTensionLink::try_new(section.elastic().clone(), [
            TensionLinkAttachment { node: 0, offset_body_m: offset },
            TensionLinkAttachment { node: 4, offset_body_m: offset },
        ], 0.9, 1, 0.2).unwrap();
        let link_mass = link.reference_mass_kg();
        let rod = SpatialMaterialRod::try_straight(1.0, vec![section; 4]).unwrap()
            .with_tension_links(vec![link], &gate).unwrap()
            .with_gravity([0.0, -9.0, 1.0], &gate).unwrap();
        let tension = ea * ((1.0 - 0.9) / 0.9);
        let weight = rod.gravity_nodal_forces_n().unwrap()[0];
        let expected_force = [-tension, -weight[1], -weight[2]];
        let expected_moment = Vec3::new(offset[0], offset[1], offset[2])
            .cross(Vec3::new(-tension, 4.5 * link_mass, -0.5 * link_mass));
        let reaction = rod.root_support_reaction(&gate).unwrap();
        for c in 0..3 {
            assert!((reaction.force_n[c].unwrap() - expected_force[c]).abs() < 1e-9);
            assert!((reaction.moment_body_nm.unwrap()[c]
                - [expected_moment.x, expected_moment.y, expected_moment.z][c]).abs() < 1e-9);
        }
    }

    #[test]
    fn tension_link_span_and_forces_ignore_large_common_world_translation() {
        let rod = linked([0.0, 0.02, 0.01], 0.7);
        let mut base = rod.rod.clone();
        perturb(&mut base, 4, 5, 0.1).unwrap();
        let link = &rod.tension_links[0];
        let reference = link.state_for(&base).unwrap();
        let mut translated = base.clone();
        for position in &mut translated.positions {
            for value in position { *value += 2.0_f64.powi(40); }
        }
        let moved = link.state_for(&translated).unwrap();
        assert_eq!(moved.length_m, reference.length_m);
        assert_eq!(moved.stored_energy_j, reference.stored_energy_j);
        assert_eq!(moved.force_n, reference.force_n);
        assert_eq!(moved.moment_body_nm, reference.moment_body_nm);
        let rounded = (0..3).map(|c|
            (moved.endpoints_m[1][c] - moved.endpoints_m[0][c]).powi(2)).sum::<f64>().sqrt();
        assert!((rounded - reference.length_m).abs() > 1e-6);
    }

    #[test]
    fn g3_complete_spatial_tangent_covers_heterogeneity_nonlocal_links_and_rotating_weight() {
        use fs_ga::So3Tangent;
        let gate = CancelGate::new();
        let sections = [0.3, 2.0, 0.7, 1.4].map(section_scaled).to_vec();
        let link = MaterialTensionLink::try_new(section().elastic().clone(), [
            TensionLinkAttachment { node: 1, offset_body_m: [0.01, 0.02, -0.01] },
            TensionLinkAttachment { node: 4, offset_body_m: [-0.01, 0.03, 0.02] },
        ], 0.68, 2, 0.2).unwrap();
        let owner = SpatialMaterialRod::try_straight(1.0, sections).unwrap()
            .with_tension_links(vec![link], &gate).unwrap()
            .with_gravity([1.0, -9.0, 0.3], &gate).unwrap();
        let mut base = owner.rod.clone();
        for node in 1..=4 {
            base.positions[node][1] += 0.001 * (node * node) as f64;
            base.frames[node] = base.frames[node].body_plus(So3Tangent::new(
                Vec3::new(0.02*node as f64,-0.01*node as f64,0.03*node as f64))).unwrap();
        }
        let model = SpatialEnergy::new(&owner);
        let exact = model.linearize(&base, &owner.free_dofs, true, &|| Ok(())).unwrap().unwrap();
        let matrix = exact.tangent.as_ref().unwrap();
        for col in 0..24 {
            let h = 2e-6;
            let mut plus = base.clone(); let mut minus = base.clone();
            perturb(&mut plus, col/6+1, col%6, h).unwrap();
            perturb(&mut minus, col/6+1, col%6, -h).unwrap();
            let ep = model.value(&plus).unwrap(); let em = model.value(&minus).unwrap();
            let expected = exact.gradient[col/6+1][col%6];
            assert!(((ep-em)/(2.0*h)-expected).abs() < 3e-6*expected.abs().max(1.0));
            let upper = model.linearize(&plus, &[], false, &|| Ok(())).unwrap().unwrap();
            let lower = model.linearize(&minus, &[], false, &|| Ok(())).unwrap().unwrap();
            for row in 0..24 {
                let probe = (upper.gradient[row/6+1][row%6]-lower.gradient[row/6+1][row%6])/(2.0*h);
                let value = matrix[24*row+col];
                assert!((probe-value).abs() < 3e-6*value.abs().max(1.0), "({row},{col}): {probe} != {value}");
            }
        }
        let selected = [23,7,0,18];
        let projected = model.linearize(&base,&selected,true,&|| Ok(())).unwrap().unwrap();
        assert_eq!(projected.gradient, exact.gradient);
        for (i,&row) in selected.iter().enumerate() {
            for (j,&col) in selected.iter().enumerate() {
                assert_eq!(projected.tangent.as_ref().unwrap()[4*i+j],matrix[24*row+col]);
            }
        }
    }

    #[test]
    fn g3_intrinsic_reference_and_current_material_share_the_analytical_pullback() {
        use super::super::reference::RodReferenceGeometry;
        use fs_ga::So3Tangent;
        let positions = vec![[0.0,0.0,0.0],[0.24,0.03,0.0],[0.46,0.12,0.0]];
        let frames = [0.0,0.2,0.4].map(|a|
            So3::exp(So3Tangent::new(Vec3::new(0.0,0.0,a))).unwrap()).to_vec();
        let reference = RodReferenceGeometry::try_new(0.5,positions,frames).unwrap();
        let owner = SpatialMaterialRod::try_from_reference(reference,
            vec![section_scaled(0.5),section_scaled(2.0)]).unwrap();
        let model = SpatialEnergy::new(&owner);
        let rest = model.linearize(&owner.rod,&owner.free_dofs,true,&|| Ok(())).unwrap().unwrap();
        assert!(rest.gradient.iter().flatten().all(|v| *v == 0.0));
        let mut base = owner.rod.clone();
        base.positions[2][2] += 0.01;
        perturb(&mut base,2,3,0.03).unwrap();
        let exact = model.linearize(&base,&owner.free_dofs,true,&|| Ok(())).unwrap().unwrap();
        for col in 0..12 {
            let h = 1e-6;
            let mut plus=base.clone();let mut minus=base.clone();
            perturb(&mut plus,col/6+1,col%6,h).unwrap();
            perturb(&mut minus,col/6+1,col%6,-h).unwrap();
            let gp=model.linearize(&plus,&[],false,&|| Ok(())).unwrap().unwrap();
            let gm=model.linearize(&minus,&[],false,&|| Ok(())).unwrap().unwrap();
            for row in 0..12 {
                let probe=(gp.gradient[row/6+1][row%6]-gm.gradient[row/6+1][row%6])/(2.0*h);
                let value=exact.tangent.as_ref().unwrap()[12*row+col];
                assert!((probe-value).abs()<3e-6*value.abs().max(1.0));
            }
        }
    }

    #[test]
    fn g0_taut_slack_selection_is_explicit_and_refusals_never_publish_partial_assembly() {
        use std::cell::Cell;
        let owner=linked([0.0;3],0.75);
        let model=SpatialEnergy::new(&owner);
        let at_kink=model.linearize(&owner.rod,&owner.free_dofs,true,&|| Ok(())).unwrap().unwrap();
        // The two endpoints have no rod element in common. Exactly slack uses
        // the zero one-sided link derivative, not a central half-stiffness.
        assert_eq!(at_kink.tangent.as_ref().unwrap()[18],0.0);
        let mut stretched=owner.rod.clone();stretched.positions[4][0]+=0.001;
        let taut=model.linearize(&stretched,&owner.free_dofs,true,&|| Ok(())).unwrap().unwrap();
        let expected=-owner.sections[0].elastic().axial_stiffness_n()/0.75;
        assert!((taut.tangent.as_ref().unwrap()[18]-expected).abs()<1e-10);
        let calls=Cell::new(0);
        let cancel=|| {
            calls.set(calls.get()+1);
            if calls.get()>=3 {Err(SpatialRodError::Cancelled)} else {Ok(())}
        };
        assert!(matches!(model.linearize(&stretched,&owner.free_dofs,true,&cancel),Err(SpatialRodError::Cancelled)));
        assert_eq!(owner.positions_m()[4],[1.0,0.0,0.0]);
        stretched.positions[4][0]=2.0;
        assert!(model.linearize(&stretched,&[],false,&|| Ok(())).is_err());
        assert_eq!(owner.positions_m()[4],[1.0,0.0,0.0]);
    }

}
