//! Quantitative, source-bound varnish/wood section study.
//!
//! Run `violin_varnish PACK_DIRECTORY [WOOD_THICKNESS_M FILM_THICKNESS_M FILM_NU]`.
//! The directory must contain the seven compiled packs listed in
//! `docs/VIOLIN_VARNISH_PHYSICS.md`. Defaults are explicitly authored geometry:
//! 3 mm wood, one 20 micrometre outer film, and isotropic film nu = 0.30.
//!
//! This is an aligned, flat, perfectly bonded section calculation. Frequency
//! shifts hold the spatial mode shape and sampled material spectrum fixed;
//! decay estimates contain the source wood and film material losses only.
//! There is no whole-violin eigenproblem, radiation, string/contact model,
//! historical recipe identification, or environmental extrapolation here.

#![deny(unsafe_code)]

use fs_la::eigen_complex::lu_complex;
use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, QueryPoint};
use fs_material::{
    coating::ResolvedCoating,
    laminate::{LaminaIntegrals, LaminaLayer, LaminateSection, LayerLaw, ReferenceSurface},
    plane_stress_viscosity::{
        MaterialPlane, PLANE_STRESS_VISCOSITY_PROPERTIES, PlaneStressViscosity,
    },
    state_point::{
        MaterialPropertySelection, OrthotropicPlaneElastic, ResolvedMaterialStatePoint,
        ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_math::{c64::C64, det};
use fs_qty::{
    Density, Dims, DynViscosity, Pressure, QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};
use std::{error::Error, f64::consts::TAU, fmt::Write as _, io::Read, path::Path};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
type Matrix = [C64; 9];
const WOOD: &str = "violin-spruce-norway-viala2020-shell";
const FILMS: [&str; 6] = [
    "varnish-oil-copal-lammlein2021-s1",
    "varnish-oil-copal-lammlein2021-s2",
    "varnish-spirit-shellac-lammlein2021-s1",
    "varnish-spirit-shellac-lammlein2021-s2",
    "varnish-mineral-ground-lammlein2021-s1",
    "varnish-mineral-ground-lammlein2021-s2",
];
const FREQUENCIES_HZ: [f64; 5] = [20., 200., 1000., 5000., 20_000.];
const STRAIN_LIMIT: f64 = 0.0003;
const USAGE: &str = "violin_varnish PACK_DIRECTORY [WOOD_THICKNESS_M FILM_THICKNESS_M FILM_NU]";

fn invalid(message: impl Into<String>) -> Box<dyn Error> {
    std::io::Error::new(std::io::ErrorKind::InvalidInput, message.into()).into()
}

fn load_pack(
    directory: &Path,
    slug: &str,
    report: &mut String,
) -> Result<NormalizedMaterialCardPack> {
    let mut bytes = Vec::new();
    std::fs::File::open(directory.join(format!("{slug}.fsmatpk")))?
        .take(16 * 1024 * 1024 + 1)
        .read_to_end(&mut bytes)?;
    if bytes.len() > 16 * 1024 * 1024 {
        return Err(invalid("a source pack exceeds 16 MiB"));
    }
    let pack = NormalizedPack::from_bytes(&bytes)?;
    if pack.pack_id() != slug {
        return Err(invalid(format!("pack identity does not match {slug}")));
    }
    writeln!(report, "pack\t{slug}\t{}", pack.content_hash())?;
    Ok(NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug.into(),
            phase: "source-identified solid; authored section interpretation".into(),
            process: "original specimen and conditioning retained in source claims".into(),
            revision: 0,
        },
        pack,
    )?)
}

fn record_source(
    pack: &NormalizedMaterialCardPack,
    state: &ResolvedMaterialStatePoint,
    label: &str,
    report: &mut String,
) -> Result<()> {
    writeln!(
        report,
        "card\t{label}\t{}\t{}",
        state.card_identity(),
        state.identity()
    )?;
    for (axis, value) in state.query_point() {
        writeln!(report, "source-axis\t{label}\t{axis}\t{value:.17e}")?;
    }
    for property in state.properties() {
        let answer = property.answer();
        pack.card().claims().verify_receipt(&answer.receipt)?;
        writeln!(
            report,
            "source-property\t{label}\t{}\t{:.17e}\t{}\t{}\t{:?}",
            property.requirement().name(),
            property.value_si(),
            answer.receipt.selected.0,
            answer.receipt.content_hash(),
            answer.evidence.value.uncertainty
        )?;
    }
    Ok(())
}

fn value(source: &ResolvedMaterialStatePoint, name: &str) -> Result<f64> {
    source
        .property(name)
        .map(|p| p.value_si())
        .ok_or_else(|| invalid(format!("missing resolved property {name}")))
}

struct WoodReference {
    layer: LaminaLayer,
    loss_pa: [f64; 9],
    omega_band: (f64, f64),
}

fn resolve_wood(
    pack: &NormalizedMaterialCardPack,
    thickness: f64,
    report: &mut String,
) -> Result<WoodReference> {
    let mut requirements = Vec::new();
    for (name, dims) in [
        ("density", Density::DIMS),
        ("young_modulus_1", Pressure::DIMS),
        ("young_modulus_2", Pressure::DIMS),
        ("poisson_ratio_12", Dims::NONE),
        ("shear_modulus_12", Pressure::DIMS),
        ("shear_modulus_13", Pressure::DIMS),
        ("shear_modulus_23", Pressure::DIMS),
        ("constant_loss_reference_omega", Dims([0, 0, -1, 0, 0, 0])),
    ] {
        requirements.push(ScalarPropertyRequirement::try_new(
            name,
            dims,
            ScalarAdmissibility::StrictlyPositive,
        )?);
    }
    for name in PLANE_STRESS_VISCOSITY_PROPERTIES {
        requirements.push(ScalarPropertyRequirement::try_new(
            name,
            DynViscosity::DIMS,
            ScalarAdmissibility::Finite,
        )?);
    }
    // The existing wood pack uses legacy dimensional axes T [K], omega [rad/s].
    let point = QueryPoint::new()
        .with("T", 293.15)?
        .with("omega", TAU * 1000.)?;
    let source = resolve_material_state_point(
        pack.card(),
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )?;
    record_source(pack, &source, WOOD, report)?;
    let mut band = (0.0_f64, f64::INFINITY);
    for p in source.properties() {
        if !p.is_constant_at_fixed_temperature() {
            return Err(invalid(
                "the section reference requires fixed-state wood coefficients",
            ));
        }
        let (lo, hi) = p.angular_frequency_band_rad_s()?.ok_or_else(|| {
            invalid("every wood coefficient needs its original frequency support")
        })?;
        band = (band.0.max(lo), band.1.min(hi));
    }
    if !(band.0 > 0. && band.1 > band.0 && band.1.is_finite()) {
        return Err(invalid("wood source frequency supports do not intersect"));
    }
    let viscosity = PlaneStressViscosity::from_material(
        source.clone(),
        MaterialPlane::Axes12,
        PLANE_STRESS_VISCOSITY_PROPERTIES,
    )?;
    let reference_omega = value(&source, "constant_loss_reference_omega")?;
    viscosity.validate_frequency(reference_omega)?;
    // This source pack explicitly authors a constant-loss closure. Multiplying
    // V by the actual omega instead would silently turn it into Kelvin--Voigt.
    let loss_pa = viscosity.matrix_pa_s().map(|v| v * reference_omega);
    let plane = OrthotropicPlaneElastic {
        plane: MaterialPlane::Axes12,
        e: [
            value(&source, "young_modulus_1")?,
            value(&source, "young_modulus_2")?,
        ],
        nu: value(&source, "poisson_ratio_12")?,
        g: value(&source, "shear_modulus_12")?,
        strain_limit: STRAIN_LIMIT,
    };
    writeln!(
        report,
        "wood-frequency-support-hz\t{:.17e}\t{:.17e}",
        band.0 / TAU,
        band.1 / TAU
    )?;
    writeln!(
        report,
        "wood-loss-closure\tconstant-loss\tomega-reference-rad-s\t{reference_omega:.17e}"
    )?;
    Ok(WoodReference {
        layer: LaminaLayer {
            thickness_m: thickness,
            orientation_rad: 0.,
            density_kg_m3: value(&source, "density")?,
            law: LayerLaw::PlaneStress {
                stiffness_pa: plane.stiffness(),
                transverse_shear_pa: [
                    value(&source, "shear_modulus_13")?,
                    value(&source, "shear_modulus_23")?,
                ],
            },
            material: Some(source),
        },
        loss_pa,
        omega_band: band,
    })
}

fn resolve_film(
    pack: &NormalizedMaterialCardPack,
    nu: f64,
    label: &str,
    report: &mut String,
) -> Result<ResolvedCoating> {
    let semantic = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            293.15,
        )?
        .with_quantity(
            "frequency",
            semantic(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
            1000.,
        )?
        .with_quantity(
            "conditioning_relative_humidity",
            QuantitySpec::dimensional(Dims::NONE),
            0.35,
        )?
        .with_quantity(
            "source_fitted_response",
            QuantitySpec::dimensional(Dims::NONE),
            1.,
        )?
        .with_quantity(
            "source_density_sample_matched",
            QuantitySpec::dimensional(Dims::NONE),
            0.,
        )?;
    let source = resolve_material_state_point(
        pack.card(),
        &point,
        &ResolvedCoating::requirements()?,
        MaterialPropertySelection::SingleClaimOnly,
    )?;
    record_source(pack, &source, label, report)?;
    Ok(ResolvedCoating::from_material(source, nu, STRAIN_LIMIT)?)
}

fn complex_abd(
    section: &LaminateSection,
    losses_pa: &[[f64; 9]],
) -> Result<(Matrix, Matrix, Matrix)> {
    let (ai, bi, di) = section.integrate_layer_matrices(losses_pa)?;
    let combine = |real: &[f64; 9], imaginary: [f64; 9]| {
        core::array::from_fn(|k| C64::new(real[k], imaginary[k]))
    };
    Ok((
        combine(section.a_pa_m(), ai),
        combine(section.b_pa_m2(), bi),
        combine(section.d_pa_m3(), di),
    ))
}

// Membrane forces are free to equilibrate: epsilon0 = -A^-1 B kappa.
// This is not the bending law for a membrane-constrained or curved violin.
fn free_membrane_bending(a: &Matrix, b: &Matrix, d: &Matrix) -> Result<Matrix> {
    let factor = lu_complex(a, 3)
        .map_err(|error| invalid(format!("complex membrane factorization failed: {error}")))?;
    let mut effective = *d;
    for j in 0..3 {
        let mut column = core::array::from_fn::<_, 3, _>(|k| b[3 * k + j]);
        factor.solve(&mut column);
        for i in 0..3 {
            for k in 0..3 {
                effective[3 * i + j] = effective[3 * i + j] - b[3 * i + k] * column[k];
            }
        }
    }
    if effective
        .iter()
        .any(|z| !z.re.is_finite() || !z.im.is_finite())
    {
        return Err(invalid("nonrepresentable free-membrane bending response"));
    }
    Ok(effective)
}

fn matrix_record(report: &mut String, label: &str, hz: f64, name: &str, m: &Matrix) -> Result<()> {
    // Row-major engineering coordinates [xx, yy, xy]. Each cell is Re then Im.
    write!(report, "matrix\t{label}\t{hz:.17e}\t{name}")?;
    for z in m {
        write!(report, "\t{:.17e}\t{:.17e}", z.re, z.im)?;
    }
    writeln!(report)?;
    Ok(())
}

fn section_records(
    report: &mut String,
    label: &str,
    hz: f64,
    section: &LaminateSection,
    abd: &(Matrix, Matrix, Matrix),
    effective: &Matrix,
) -> Result<()> {
    writeln!(
        report,
        "inertia\t{label}\t{hz:.17e}\t{:.17e}\t{:.17e}\t{:.17e}",
        section.i0_kg_m2(),
        section.i1_kg_m(),
        section.i2_kg()
    )?;
    for (name, matrix) in [
        ("A-Pa-m", &abd.0),
        ("B-Pa-m2", &abd.1),
        ("D-Pa-m3", &abd.2),
        ("D-free-Pa-m3", effective),
    ] {
        matrix_record(report, label, hz, name, matrix)?;
    }
    Ok(())
}

fn form(q: &[f64; 9], left: &[C64; 3], right: &[C64; 3]) -> f64 {
    let mut sum = C64::ZERO;
    for i in 0..3 {
        for j in 0..3 {
            sum = sum + (left[i].conj() * right[j]).scale(q[3 * i + j]);
        }
    }
    sum.re
}

// Integral of e(z)^H Q e(z), e(z)=epsilon0+z*kappa. With complex amplitudes
// this is four times average stored energy for Q', and loss/cycle divided by
// pi for Q''. The common amplitude normalization cancels in every ratio.
fn layer_weight(
    q: &[f64; 9],
    geometry: &LaminaIntegrals,
    epsilon: &[C64; 3],
    curvature: &[C64; 3],
) -> f64 {
    let (z0, z1) = (geometry.z_bottom_m, geometry.z_top_m);
    let m0 = z1 - z0;
    let m1 = m0 * (z0 + z1) / 2.;
    let m2 = m0 * (z0 * z0 + z0 * z1 + z1 * z1) / 3.;
    m0 * form(q, epsilon, epsilon)
        + 2. * m1 * form(q, epsilon, curvature)
        + m2 * form(q, curvature, curvature)
}

fn directional_records(
    report: &mut String,
    label: &str,
    nominal_hz: f64,
    omega: f64,
    bare: &LaminateSection,
    bare_effective: &Matrix,
    section: &LaminateSection,
    losses_pa: &[[f64; 9]],
    abd: &(Matrix, Matrix, Matrix),
    effective: &Matrix,
) -> Result<()> {
    let membrane = lu_complex(&abd.0, 3)
        .map_err(|error| invalid(format!("complex membrane factorization failed: {error}")))?;
    for (direction, index) in [("L", 0), ("R", 1)] {
        let mut curvature = [C64::ZERO; 3];
        curvature[index] = C64::ONE;
        let mut epsilon = core::array::from_fn::<_, 3, _>(|i| -abd.1[3 * i + index]);
        membrane.solve(&mut epsilon);
        let mut stored = [0.; 2];
        let mut lost = [0.; 2];
        for (i, geometry) in section.layer_integrals().iter().enumerate() {
            stored[i] = layer_weight(&geometry.q_bar_pa, geometry, &epsilon, &curvature);
            lost[i] = layer_weight(&losses_pa[i], geometry, &epsilon, &curvature);
        }
        let storage_sum: f64 = stored.iter().sum();
        let loss_sum: f64 = lost.iter().sum();
        let stiffness = effective[4 * index];
        if !(storage_sum > 0. && loss_sum > 0.)
            || stored
                .iter()
                .chain(&lost)
                .any(|v| !v.is_finite() || *v < 0.)
            || (storage_sum - stiffness.re).abs() > 1e-10 * stiffness.re
            || (loss_sum - stiffness.im).abs() > 1e-10 * stiffness.im
        {
            return Err(invalid(
                "layer energy does not reconstruct the condensed passive section",
            ));
        }
        let eta_bare = bare_effective[4 * index].im / bare_effective[4 * index].re;
        let eta_wood = lost[0] / storage_sum;
        let eta_film = lost[1] / storage_sum;
        let eta = eta_wood + eta_film;
        let mass_ratio = section.i0_kg_m2() / bare.i0_kg_m2();
        let stiffness_ratio = stiffness.re / bare_effective[4 * index].re;
        let frequency_ratio = det::sqrt(stiffness_ratio / mass_ratio);
        let first_order = 0.5 * ((stiffness_ratio - 1.) - (mass_ratio - 1.));
        let tau_bare = 2. / (omega * eta_bare);
        let tau_coated = 2. / (omega * eta);
        writeln!(
            report,
            "direction\t{label}\t{nominal_hz:.17e}\t{:.17e}\t{direction}\t{mass_ratio:.17e}\t{stiffness_ratio:.17e}\t{frequency_ratio:.17e}\t{first_order:.17e}\t{eta_bare:.17e}\t{eta_wood:.17e}\t{eta_film:.17e}\t{eta:.17e}\t{tau_bare:.17e}\t{tau_coated:.17e}\t{:.17e}\t{:.17e}\t{:.17e}",
            omega / TAU,
            det::ln(1000.) * tau_bare,
            det::ln(1000.) * tau_coated,
            stored[1] / storage_sum
        )?;
    }
    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).take(5).collect();
    if args.as_slice() == ["--help"] {
        println!("{USAGE}");
        return Ok(());
    }
    if !matches!(args.len(), 1 | 4) {
        return Err(invalid(USAGE));
    }
    let (wood_h, film_h, nu) = if args.len() == 4 {
        (
            args[1].parse::<f64>()?,
            args[2].parse::<f64>()?,
            args[3].parse::<f64>()?,
        )
    } else {
        (0.003, 20e-6, 0.30)
    };
    if !wood_h.is_finite()
        || wood_h <= 0.
        || !film_h.is_finite()
        || film_h < 0.
        || !nu.is_finite()
        || nu <= -1.
        || nu >= 0.5
    {
        return Err(invalid(
            "need wood thickness > 0, film thickness >= 0, -1 < film nu < 0.5; all finite",
        ));
    }
    let directory = Path::new(&args[0]);
    let mut report = String::from("frankensim.violin-varnish-section.v1\n");
    writeln!(
        report,
        "scope\taligned flat section; perfect bond; free membrane force; fixed spatial curvature and translational modal inertia; material loss only"
    )?;
    writeln!(
        report,
        "assumptions\twood-thickness-m\t{wood_h:.17e}\tfilm-thickness-m\t{film_h:.17e}\tfilm-nu\t{nu:.17e}\tstrain-limit\t{STRAIN_LIMIT:.17e}"
    )?;
    writeln!(
        report,
        "source-caution\tmodern fitted film specimens, group-mean unmatched density; wood nu_LR is a Sitka surrogate; separately conditioned wood and film are not one validated coated specimen"
    )?;
    writeln!(
        report,
        "reference-surface\toriginal wood midsurface; external coat is above the wood; no penetration zone added"
    )?;
    writeln!(
        report,
        "matrix-columns\trow-major [xx,yy,xy], each cell Re then Im; inertia columns I0-kg-m-2 I1-kg-m-1 I2-kg"
    )?;
    writeln!(
        report,
        "direction-columns\tpack\tnominal-hz\tevaluation-hz\taxis\tmass-ratio\tbending-ratio\tfixed-spectrum-frequency-ratio\tfirst-order-delta-f-over-f\teta-bare\teta-wood-contribution\teta-film-contribution\teta-total\ttau-amplitude-bare-s\ttau-amplitude-coated-s\tT60-bare-s\tT60-coated-s\tfilm-storage-fraction"
    )?;
    let wood_pack = load_pack(directory, WOOD, &mut report)?;
    let wood = resolve_wood(&wood_pack, wood_h, &mut report)?;
    let reference = ReferenceSurface::Offset(wood_h / 2.);
    let bare = LaminateSection::try_new(vec![wood.layer.clone()], reference)?;
    let bare_abd = complex_abd(&bare, &[wood.loss_pa])?;
    let bare_effective = free_membrane_bending(&bare_abd.0, &bare_abd.1, &bare_abd.2)?;
    section_records(
        &mut report,
        "bare-wood",
        1000.,
        &bare,
        &bare_abd,
        &bare_effective,
    )?;
    for label in FILMS {
        let pack = load_pack(directory, label, &mut report)?;
        let film = resolve_film(&pack, nu, label, &mut report)?;
        let (_, lowered) = film.lower_to_tensor_prony(16, 0.001, 0.001)?;
        writeln!(
            report,
            "prony-grid-certificate\t{label}\t{}\t{}\t{:.17e}\t{:.17e}\t{:.17e}",
            lowered.model.terms.len(),
            lowered.verification_points,
            lowered.sup_rel_err,
            lowered.sup_rel_storage_err,
            lowered.sup_rel_loss_err
        )?;
        for nominal_hz in FREQUENCIES_HZ {
            let (ep, epp) = film.law().modulus(TAU * nominal_hz)?;
            writeln!(
                report,
                "film-spectrum\t{label}\t{nominal_hz:.17e}\t{ep:.17e}\t{epp:.17e}\t{:.17e}",
                epp / ep
            )?;
            let requested = TAU * nominal_hz;
            let omega = if requested >= wood.omega_band.0 && requested <= wood.omega_band.1 {
                requested
            } else if nominal_hz == 5000.
                && requested > wood.omega_band.1
                && (requested - wood.omega_band.1) / requested < 1e-8
            {
                // The existing pack prints its 5000 Hz edge as 31415.9265 rad/s.
                // Evaluate that actual source edge and expose both coordinates.
                writeln!(
                    report,
                    "source-edge-coordinate\t{label}\t{nominal_hz:.17e}\t{:.17e}",
                    wood.omega_band.1 / TAU
                )?;
                wood.omega_band.1
            } else {
                writeln!(
                    report,
                    "section-unavailable\t{label}\t{nominal_hz:.17e}\twood-source-frequency-out-of-band"
                )?;
                continue;
            };
            let mut layers = vec![wood.layer.clone()];
            let mut losses = vec![wood.loss_pa];
            if film_h > 0. {
                let slice = film.at_frequency(film_h, omega)?;
                losses.push(slice.plane_viscosity_pa_s.map(|v| omega * v));
                layers.push(slice.layer);
            }
            let section = LaminateSection::try_new(layers, reference)?;
            let abd = complex_abd(&section, &losses)?;
            let effective = free_membrane_bending(&abd.0, &abd.1, &abd.2)?;
            section_records(&mut report, label, omega / TAU, &section, &abd, &effective)?;
            directional_records(
                &mut report,
                label,
                nominal_hz,
                omega,
                &bare,
                &bare_effective,
                &section,
                &losses,
                &abd,
                &effective,
            )?;
        }
    }
    // No successful report is emitted until every source and section is admitted.
    print!("{report}");
    Ok(())
}
