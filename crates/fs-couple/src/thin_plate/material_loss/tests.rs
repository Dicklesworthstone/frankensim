//! Synthetic constitutive fixtures: implementation checks, not measured timbre.
use super::*;
use crate::plate_viscosity::{PlateViscosity, PlateViscosityLimits};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
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
        let loss = DiagonalPlateViscosity::new(&s, s.bending_viscosity().unwrap(), &m).unwrap();
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
    assert_eq!(other.bending_viscosity().unwrap()[1].isotropic().unwrap().viscosity_pa_s, 4e6);
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

/// Verbatim pre-move copy of the scalar regional evaluator (origin 5278e835c:
/// `thin_plate/material_loss.rs` RegionalBendingLoss, the viscous branch of
/// `thin_plate.rs` PlateBendingLaw::zeta with no thermal law, and
/// bending_energy_fraction). It is the bitwise oracle for the one-owner routing
/// through `plate_viscosity::DiagonalPlateViscosity`; production never calls it.
mod legacy {
    use super::*;
    pub(super) fn modal_zeta(specimen: &ResolvedPlateChart, viscosities: &[IsotropicPlateBendingViscosity],
        model: &PlateModel, phi: &[f64], omega: f64) -> Result<f64, AcousticRealizeError>
    {
        let mut laws = Vec::with_capacity(viscosities.len());
        for (region, &viscous) in specimen.regions.iter().zip(viscosities) {
            let young = plate_property(&region.input.material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?;
            laws.push((viscous, viscous_time(viscous, young)?));
        }
        let mesh = &specimen.chart.mesh;
        let mut elements = Vec::with_capacity(mesh.tris.len());
        for (element, tri) in mesh.tris.iter().enumerate() {
            let region = specimen.region_for_element[element];
            let x = tri.map(|node| mesh.nodes[node].0);
            let y = tri.map(|node| mesh.nodes[node].1);
            let (stiffness, _) = fs_plate::dkt_stiffness(&x, &y,
                &specimen.regions[region].section.d, element).map_err(map_plate)?;
            let dofs: [Option<usize>; 9] = std::array::from_fn(|i| model.dof_map[3 * tri[i / 3] + i % 3]);
            elements.push((dofs, stiffness, region));
        }
        let zetas: Vec<f64> = laws.iter().map(|&(viscous, viscous_time_s)| {
            let mut zeta = 0.0;
            let (lo, hi) = viscous.omega_band_rad_s;
            if !omega.is_finite() || omega < lo || omega > hi {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "retained plate frequency lies outside material bending-loss applicability",
                });
            }
            zeta += 0.5 * viscous_time_s * omega;
            if !zeta.is_finite() {
                return Err(AcousticRealizeError::InvalidDescription { what: "plate bending loss is unrepresentable" });
            }
            Ok(zeta)
        }).collect::<Result<_, _>>()?;
        let mut k_phi = vec![0.0; phi.len()];
        model.k.spmv(phi, &mut k_phi);
        let total_energy: f64 = phi.iter().zip(k_phi).map(|(p, k)| p * k).sum();
        let mut weighted_bending = 0.0;
        for (dofs, stiffness, region) in &elements {
            let q: [f64; 9] = dofs.map(|dof| dof.map_or(0.0, |i| phi[i]));
            let mut energy = 0.0;
            for r in 0..9 {
                for c in 0..9 {
                    energy += q[r] * stiffness[9 * r + c] * q[c];
                }
            }
            weighted_bending += zetas[*region] * energy;
        }
        let fraction = weighted_bending / total_energy;
        if !weighted_bending.is_finite() || weighted_bending < 0.0 || !total_energy.is_finite()
            || total_energy <= 0.0 || !fraction.is_finite() {
            return Err(AcousticRealizeError::InvalidDescription { what: "legacy oracle refused" });
        }
        Ok(fraction)
    }
}

#[test]
fn isotropic_regional_loss_routed_through_plate_viscosity_is_bitwise_legacy() {
    // Uniform (one material in both regions) and two-region isotropic specimens,
    // unstressed and prestressed: zeta bits and radiator bits equal the pre-move code.
    let pairs = [
        (material(200e9, Some(2e6), 100_000.0), material(200e9, Some(2e6), 100_000.0)),
        (material(200e9, Some(2e6), 100_000.0), material(100e9, Some(5e6), 100_000.0)),
    ];
    let mut compared = 0;
    for (left, right) in pairs {
        for tension in [0.0, 1000.0] {
            let s = specimen(left.clone(), right.clone()).with_kelvin_voigt_bending_loss().unwrap();
            let mut o = options(); o.assembly.pretension = tension;
            let laws: Vec<_> = s.bending_viscosity().unwrap().iter().map(|l| *l.isotropic().unwrap()).collect();
            let model = s.chart().assemble(&[], &o.assembly).unwrap();
            let report = modes(&model, o.eigenvalue_window, &fs_modal::SliceOptions::default()).unwrap();
            let loss = DiagonalPlateViscosity::new(&s, s.bending_viscosity().unwrap(), &model).unwrap();
            let radiators = certified_material_chart_radiators(&s, &o).unwrap();
            let mut k = 0;
            for pair in report.modes.iter().take(o.n_modes) {
                let omega = pair.lambda.max(0.0).sqrt();
                if !(omega > 0.0) { continue; }
                let old = legacy::modal_zeta(&s, &laws, &model, &pair.phi, omega).unwrap();
                assert!(old > 0.0);
                assert_eq!(loss.modal_zeta(&model, &pair.phi, omega).unwrap().to_bits(), old.to_bits());
                assert_eq!(radiators[k].zeta.to_bits(), old.to_bits(), "radiator {k} at tension {tension}");
                k += 1;
                compared += 1;
            }
            assert_eq!(k, radiators.len());
        }
    }
    assert_eq!(compared, 12, "two specimens x two tensions x three modes");
    // The uniform plate's viscous term goes through the same owner, bit for bit.
    let law = IsotropicPlateBendingViscosity { viscosity_pa_s: 2e6, omega_band_rad_s: (0.0, 100_000.0),
        material_state_identity: None };
    let tau = viscous_time(law, 200e9).unwrap();
    let uniform = PlateBendingLaw { thermal: None, viscous: Some(law), viscous_time_s: tau, thickness: 0.001 };
    for omega in [1.0, 377.0, 6_283.185_307_179_586, 99_999.0] {
        let mut old = 0.0;
        old += 0.5 * tau * omega;
        assert_eq!(uniform.zeta(omega).unwrap().to_bits(), old.to_bits());
    }
}

// Synthetic orthotropic body with explicit reduced plane-stress viscosities; not measured wood.
fn orthotropic(scale: f64, omit: Option<usize>) -> ResolvedMaterialStatePoint {
    let mut values = vec![("density", Density::DIMS, 600.0),
        ("young_modulus_1", Pressure::DIMS, 8e9), ("young_modulus_2", Pressure::DIMS, 2e9),
        ("poisson_ratio_12", Dims::NONE, 0.2), ("shear_modulus_12", Pressure::DIMS, 0.6e9)];
    for (i, (key, value)) in PLANE_STRESS_VISCOSITY_PROPERTIES.into_iter()
        .zip([3e7, 6e6, 4e6, 1e6, 4e5, -3e5]).enumerate() {
        if omit != Some(i) { values.push((key, DynViscosity::DIMS, value * scale)); }
    }
    let mut claims = ClaimSet::new(); let mut requirements = Vec::new();
    for (key, dims, value) in values {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 310.0).with("omega", 0.0, 100_000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: Provenance { source: "synthetic orthotropic bending viscosity".into(), license: "CC0-1.0".into(), artifact: None },
            observations: vec![],
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "fixture".into(), phase: "solid".into(),
        process: "orthotropic test".into(), revision: 0 }, claims, vec![]).unwrap();
    resolve_material_state_point(&card, &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn orthotropic_specimen(left: ResolvedMaterialStatePoint, right: ResolvedMaterialStatePoint) -> ResolvedPlateChart {
    let mesh = PlateMesh::rectangle(0.4, 0.3, 4, 4);
    let inputs = [("left", left, 0.17), ("right", right, -0.23)].into_iter().enumerate()
        .map(|(r, (name, material, angle))| {
            let triangle_indices = mesh.tris.iter().enumerate().filter_map(|(e, tri)| {
                let x = tri.iter().map(|&i| mesh.nodes[i].0).sum::<f64>() / 3.0;
                ((x >= 0.2) == (r == 1)).then_some(e)
            }).collect();
            PlateRegionMaterial { region: PlateRegion { name: name.into(), triangle_indices }, material,
                model: PlateMaterialModel::Orthotropic12,
                thickness_constraint: PlateThicknessConstraint::FixedThickness(0.002), material_angle_rad: angle }
        }).collect();
    compile_plate_material_chart(mesh, PlateMesh::rectangle_boundary(4, 4), inputs).unwrap()
}

#[test]
fn orthotropic_region_now_admitted_and_matches_full_matrix_diagonal() {
    let plain = orthotropic_specimen(orthotropic(1.0, None), orthotropic(7.0, None));
    // The former isotropic-only refusal now resolves through the tensor owner.
    let selected = plain.clone().with_kelvin_voigt_bending_loss().unwrap();
    let laws = selected.bending_viscosity().unwrap();
    assert!(laws.iter().all(|l| matches!(l, RegionalBendingViscosity::Orthotropic(_)) && l.isotropic().is_none()));
    let o = options();
    let scalar = certified_material_chart_radiators(&selected, &o).unwrap();
    // Independent full projection on the unselected chart (counted once).
    let gate = CancelGate::new_clock_free();
    let tensor = PlateViscosity::from_chart(&plain, PLANE_STRESS_VISCOSITY_PROPERTIES, &gate).unwrap();
    let limits = PlateViscosityLimits { max_modes: 8, max_projection_terms: 1_000_000, reference_resistance_n_s_m: 1.0 };
    let image = plain.certified_viscous_radiators(&o, &tensor, limits, &gate).unwrap();
    let c = image.damping().matrix_per_s();
    let n = image.radiators().len();
    assert_eq!(scalar.len(), n);
    // zeta_i = omega q^T C q / (2 q^T K q) equals C_ii / (2 omega_i) exactly for an
    // exact eigenpair; the shared solve's residual is measured here, not assumed.
    let model = plain.chart().assemble(&[], &o.assembly).unwrap();
    let report = modes(&model, o.eigenvalue_window, &fs_modal::SliceOptions::default()).unwrap();
    for (i, pair) in report.modes.iter().take(n).enumerate() {
        let (mut kp, mut mp) = (vec![0.0; pair.phi.len()], vec![0.0; pair.phi.len()]);
        model.k.spmv(&pair.phi, &mut kp); model.m.spmv(&pair.phi, &mut mp);
        let stiffness: f64 = pair.phi.iter().zip(&kp).map(|(a, b)| a * b).sum();
        let mass: f64 = pair.phi.iter().zip(&mp).map(|(a, b)| a * b).sum();
        let residual = (stiffness - pair.lambda * mass).abs() / (pair.lambda * mass);
        let reference = c[n * i + i] / (2.0 * image.radiators()[i].omega);
        assert_eq!(scalar[i].omega.to_bits(), image.radiators()[i].omega.to_bits());
        let relative = (scalar[i].zeta - reference).abs() / reference;
        assert!(reference > 0.0 && relative <= 1e-13 + 2.0 * residual,
            "mode {i}: relative {relative:e}, eigen residual {residual:e}");
    }
    // A chart with selected scalar loss cannot also select the full tensor.
    assert!(PlateViscosity::from_chart(&selected, PLANE_STRESS_VISCOSITY_PROPERTIES, &gate).is_err());
    // No coefficient is inferred: a missing 16 term refuses.
    assert!(orthotropic_specimen(orthotropic(1.0, Some(4)), orthotropic(7.0, None))
        .with_kelvin_voigt_bending_loss().is_err());
}

#[path = "tests/coupled.rs"]
mod coupled;
