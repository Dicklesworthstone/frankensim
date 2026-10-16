//! Synthetic constitutive fixtures: implementation checks, not measured timbre.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};

fn material(young: f64, eta: Option<f64>, upper: f64) -> ResolvedMaterialStatePoint {
    let mut values = vec![("density", Density::DIMS, 7800.0),
        ("young_modulus", Pressure::DIMS, young), ("poisson_ratio", Dims::NONE, 0.3)];
    if let Some(eta) = eta { values.push((KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, DynViscosity::DIMS, eta)); }
    let mut claims = ClaimSet::new(); let mut requirements = Vec::new();
    for (key, dims, value) in values {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 310.0).with("omega", 0.0, upper),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: Provenance { source: "synthetic regional bending law".into(), license: "CC0-1.0".into(), artifact: None },
            observations: vec![],
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "fixture".into(), phase: "solid".into(),
        process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
    resolve_material_state_point(&card, &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn specimen(left: ResolvedMaterialStatePoint, right: ResolvedMaterialStatePoint) -> ResolvedPlateChart {
    let mesh = PlateMesh::rectangle(0.4, 0.3, 4, 4);
    let inputs = [("left", left), ("right", right)].into_iter().enumerate().map(|(r, (name, material))| {
        let triangle_indices = mesh.tris.iter().enumerate().filter_map(|(e, tri)| {
            let x = tri.iter().map(|&i| mesh.nodes[i].0).sum::<f64>() / 3.0;
            ((x >= 0.2) == (r == 1)).then_some(e)
        }).collect();
        PlateRegionMaterial { region: PlateRegion { name: name.into(), triangle_indices }, material,
            model: PlateMaterialModel::Isotropic,
            thickness_constraint: PlateThicknessConstraint::FixedThickness(0.001), material_angle_rad: 0.0 }
    }).collect();
    compile_plate_material_chart(mesh, PlateMesh::rectangle_boundary(4, 4), inputs).unwrap()
}
fn options() -> PlateChartRadiation {
    let mut weights = vec![0.0; 25]; weights[12] = 1.0;
    PlateChartRadiation { assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::SimplySupported },
        eigenvalue_window: (1.0, 1e8), n_modes: 3, damping_ratio: 0.0, unit_force_weights: weights }
}
fn near(a: f64, b: f64) { assert!((a-b).abs() <= 2e-10 * a.abs().max(b.abs()).max(1e-14), "{a:e} != {b:e}"); }

#[test]
fn constant_curvature_uses_regional_bending_energy_and_does_not_damp_prestress() {
    let s = specimen(material(200e9, Some(2e6), 100_000.0), material(100e9, Some(5e6), 100_000.0));
    // Independent patch field on an unconstrained chart; not a manufactured
    // eigenvector computed by the modal/loss assembly being tested.
    let s = compile_plate_material_chart(s.chart.mesh.clone(), vec![],
        s.regions().iter().map(|r| r.input.clone()).collect()).unwrap()
        .with_kelvin_voigt_bending_loss().unwrap();
    let omega = 100.0;
    let mut reference = 0.0;
    for r in s.regions() {
        let eta = plate_property(&r.input.material, KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, DynViscosity::DIMS).unwrap();
        let e = plate_property(&r.input.material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS).unwrap();
        // kappa_xx=1, other curvatures zero: integral kappa^T D kappa.
        let d11 = e * 0.001_f64.powi(3) / (12.0 * (1.0-0.3_f64.powi(2)));
        reference += 0.5 * omega * (eta/e) * r.area_m2 * d11;
    }
    let mut ratios = Vec::new();
    for tension in [0.0, 1000.0] {
        let mut o = options(); o.assembly.pretension = tension;
        let m = s.chart().assemble(&[], &o.assembly).unwrap();
        let mut phi = vec![0.0; m.free];
        for (node, &(x, _)) in s.chart.mesh.nodes.iter().enumerate() {
            for (axis, value) in [0.5*x*x, x, 0.0].into_iter().enumerate() {
                phi[m.dof_map[3*node+axis].unwrap()] = value;
            }
        }
        let loss = RegionalBendingLoss::new(&s, s.bending_viscosity().unwrap(), &m).unwrap();
        let ratio = loss.modal_zeta(&m, &phi, omega).unwrap();
        let mut kphi = vec![0.0; phi.len()]; m.k.spmv(&phi, &mut kphi);
        let storage: f64 = phi.iter().zip(kphi).map(|(p, k)| p*k).sum();
        near(ratio*storage, reference);
        let scaled: Vec<_> = phi.iter().map(|v| -3.0*v).collect();
        near(loss.modal_zeta(&m, &scaled, omega).unwrap(), ratio);
        ratios.push(ratio);
    }
    assert!(ratios[1] < ratios[0], "conservative tension acquired material loss");
}

#[test]
fn common_viscous_time_recovers_the_uniform_law_and_region_order_is_irrelevant() {
    let s = specimen(material(200e9, Some(2e6), 100_000.0), material(100e9, Some(1e6), 100_000.0));
    let selected = s.clone().with_kelvin_voigt_bending_loss().unwrap();
    let mut inputs: Vec<_> = s.regions().iter().map(|r| r.input.clone()).collect(); inputs.reverse();
    for input in &mut inputs { input.region.triangle_indices.reverse(); }
    let reversed = compile_plate_material_chart(s.chart.mesh.clone(), s.chart.boundary_nodes.clone(), inputs)
        .unwrap().with_kelvin_voigt_bending_loss().unwrap();
    assert_eq!(selected.specimen_identity(), reversed.specimen_identity());
    let a = certified_material_chart_radiators(&selected, &options()).unwrap();
    let b = certified_material_chart_radiators(&reversed, &options()).unwrap();
    for (a, b) in a.iter().zip(&b) {
        near(a.zeta, 0.5e-5 * a.omega);
        assert_eq!(a.zeta.to_bits(), b.zeta.to_bits());
        assert_eq!(a.omega.to_bits(), b.omega.to_bits());
        assert_eq!(a.drive_participation.to_bits(), b.drive_participation.to_bits());
    }
}

#[test]
fn viscosity_only_substitution_changes_modal_decay_but_not_mass_pitch_or_ports() {
    let s = specimen(material(200e9, Some(2e6), 100_000.0), material(100e9, Some(1e6), 100_000.0))
        .with_kelvin_voigt_bending_loss().unwrap();
    let identity = s.specimen_identity();
    let other = s.with_region_material("right", &material(100e9, Some(4e6), 100_000.0)).unwrap();
    assert_eq!(s.specimen_identity(), identity);
    assert_eq!(other.bending_viscosity().unwrap()[1].viscosity_pa_s, 4e6);
    assert_eq!(other.regions()[0].input.material, s.regions()[0].input.material);
    let a = certified_material_chart_radiators(&s, &options()).unwrap();
    let b = certified_material_chart_radiators(&other, &options()).unwrap();
    for (a, b) in a.iter().zip(&b) {
        assert_eq!(a.omega, b.omega); assert_eq!(a.mass_kg, b.mass_kg);
        assert_eq!(a.area_m2, b.area_m2); assert_eq!(a.drive_participation, b.drive_participation);
        assert!(b.zeta > a.zeta);
    }
    assert!(s.with_region_material("right", &material(100e9, None, 100_000.0)).is_err());
    assert_eq!(s.specimen_identity(), identity);
}

#[test]
fn zero_loss_keeps_source_domains_and_authored_loss_cannot_be_double_counted() {
    let source = material(200e9, Some(0.0), 10.0);
    let selected = specimen(source.clone(), source).with_kelvin_voigt_bending_loss().unwrap();
    assert!(certified_material_chart_radiators(&selected, &options()).is_err());
    let selected = specimen(material(200e9, Some(0.0), 100_000.0),
        material(200e9, Some(0.0), 100_000.0)).with_kelvin_voigt_bending_loss().unwrap();
    assert!(certified_material_chart_radiators(&selected, &options()).unwrap().iter().all(|r| r.zeta == 0.0));
    let mut o = options(); o.damping_ratio = 0.001;
    assert!(certified_material_chart_radiators(&selected, &o).is_err());
}

#[test]
fn unselected_loss_delegates_legacy_chart_bits_and_missing_material_data_refuse() {
    let s = specimen(material(200e9, None, 100_000.0), material(100e9, None, 100_000.0));
    let mut o = options(); o.damping_ratio = 0.01;
    let a = certified_chart_radiators(s.chart(), &o).unwrap();
    let b = certified_material_chart_radiators(&s, &o).unwrap();
    for (a,b) in a.iter().zip(&b) {
        assert_eq!(a.omega, b.omega); assert_eq!(a.zeta, b.zeta);
        assert_eq!(a.mass_kg, b.mass_kg); assert_eq!(a.area_m2, b.area_m2);
        assert_eq!(a.drive_participation, b.drive_participation);
    }
    assert!(s.with_kelvin_voigt_bending_loss().is_err());
    let s = specimen(material(200e9, Some(-1.0), 100_000.0), material(100e9, Some(1e6), 100_000.0));
    assert!(s.with_kelvin_voigt_bending_loss().is_err());
}

#[path = "tests/coupled.rs"]
mod coupled;
