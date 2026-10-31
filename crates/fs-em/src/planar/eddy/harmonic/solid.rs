//! Current-constrained solid winding sections on the same projected field owner.
//!
//! In one invariant-z connected section, S=int(sigma dA), prescribed integral
//! J dA=n I, and the shared exp(-i omega t) convention give
//! `E = n I/S + i omega (A - mean_sigma(A))`. Eliminating the axial potential
//! gradient leaves the existing `(K-i omega C) A=F I` operator. The reciprocal
//! source profile is `J_unit=n sigma/S`; terminal contribution is
//! `depth*n^2/S * I - i omega int(J_unit A dV)`.
//!
//! Actual total E/J, including local drive/circulating-current cross terms,
//! produce spatial heat once. Only its REGION integral separates into DC
//! resistance plus nonnegative circulating-current loss. Adding those two
//! densities cell by cell would deposit heat in the wrong places.

use super::*;
use std::cell::RefCell;

/// One complete edge-connected conducting section carries `signed_turns * I`.
/// Return sections must be named explicitly with their own sign; the original
/// winding source remains balanced on every magnetic mesh component.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolidWindingRegion {
    pub region: u32,
    pub winding: usize,
    pub signed_turns: f64,
}

/// Original source cell with its ACTUAL total P1 current and irreversible heat.
/// Neither a uniform DC-current heat nor an eddy-only heat is added elsewhere.
#[derive(Clone, Debug)]
pub struct HarmonicSolidCell {
    pub source: EddyConductor,
    pub binding: SolidWindingRegion,
    pub volume_m3: f64,
    pub electric_field_v_m: [C64; 3],
    pub current_density_a_m2: [C64; 3],
    pub mean_power_w: f64,
}

#[derive(Clone, Debug)]
pub struct SolidWindingCurrent {
    pub binding: SolidWindingRegion,
    pub target_current_a: C64,
    pub current_a: C64,
    pub residual_a: C64,
}

pub(super) struct SolidWindings {
    pub(super) bindings: Vec<SolidWindingRegion>,
    /// Parallel to the original conductivity projection's sorted region list.
    region_bindings: Vec<Option<usize>>,
    electric_per_current_v_m_a: Vec<f64>,
    pub(super) resistance_ohm: Vec<f64>,
}

impl SolidWindings {
    pub(super) fn new(
        basis: &PlanarIncrementalWindingBasis<'_>,
        conduction: &Conduction,
        sources: &[EddyConductor],
        mut bindings: Vec<SolidWindingRegion>,
        controls: HarmonicWindingControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let count = basis.winding_count();
        if bindings.len() > conduction.regions.len() {
            return Err(HarmonicWindingError::Invalid(
                "too many driven solid winding regions",
            ));
        }
        bindings.sort_by_key(|binding| binding.region);
        let depth = basis.operating_field().problem().thickness_m();
        let mut region_bindings = vec![None; conduction.regions.len()];
        let mut electric_per_current_v_m_a = vec![0.0; conduction.regions.len()];
        let mut resistance_ohm = vec![0.0; count];
        for (index, binding) in bindings.iter().enumerate() {
            poll(cancelled)?;
            if binding.winding >= count
                || !binding.signed_turns.is_finite()
                || binding.signed_turns == 0.0
                || index > 0 && bindings[index - 1].region == binding.region
            {
                return Err(HarmonicWindingError::Invalid(
                    "invalid or duplicate solid winding region binding",
                ));
            }
            let group = conduction
                .regions
                .binary_search(&binding.region)
                .map_err(|_| {
                    HarmonicWindingError::Invalid(
                        "solid winding region has no admitted conductor section",
                    )
                })?;
            let drive = finite((binding.signed_turns * depth) / conduction.weights[group])?;
            let resistance = finite((binding.signed_turns * depth) * drive)?;
            if drive == 0.0 || !positive(resistance) {
                return Err(HarmonicWindingError::Invalid(
                    "solid winding DC drive or resistance is unrepresentable",
                ));
            }
            region_bindings[group] = Some(index);
            electric_per_current_v_m_a[group] = drive;
            resistance_ohm[binding.winding] = finite(resistance_ohm[binding.winding] + resistance)?;
        }
        // Check the original source rather than silently replacing it. The L1
        // volume-weighted discrepancy bounds every free weak-load component;
        // multiplying by max(1, |I_DC|) also checks the actual retained bias.
        let mut profile_defect_a_m = vec![0.0; count];
        for (source, cell) in sources.iter().zip(&conduction.cells) {
            poll(cancelled)?;
            if let Some(index) = region_bindings[cell.group] {
                let binding = bindings[index];
                let expected =
                    finite(cell.conductivity_s_m * electric_per_current_v_m_a[cell.group])?;
                let original = basis.profiles_per_m2()[binding.winding][source.triangle];
                if expected == 0.0 || original == 0.0 || expected.signum() != original.signum() {
                    return Err(HarmonicWindingError::Invalid(
                        "solid winding profile has the wrong current direction or support",
                    ));
                }
                let defect = finite(
                    (original - expected).abs()
                        * cell.volume_m3
                        * basis.operating_currents_a()[binding.winding].abs().max(1.0),
                )?;
                profile_defect_a_m[binding.winding] =
                    finite(profile_defect_a_m[binding.winding] + defect)?;
                if basis
                    .profiles_per_m2()
                    .iter()
                    .enumerate()
                    .any(|(winding, profile)| {
                        winding != binding.winding && profile[source.triangle] != 0.0
                    })
                {
                    return Err(HarmonicWindingError::Invalid(
                        "one connected solid section cannot carry multiple independent winding coordinates",
                    ));
                }
            } else if basis
                .profiles_per_m2()
                .iter()
                .any(|profile| profile[source.triangle] != 0.0)
            {
                return Err(HarmonicWindingError::Invalid(
                    "active conductor profile requires an explicit solid current binding",
                ));
            }
        }
        let allowance = controls
            .absolute_residual_a_m
            .min(controls.current_balance_tolerance_a * depth);
        if profile_defect_a_m.iter().any(|&defect| defect > allowance) {
            return Err(HarmonicWindingError::Invalid(
                "existing DC winding profile does not match the source-normalized solid current",
            ));
        }
        Ok(Self {
            bindings,
            region_bindings,
            electric_per_current_v_m_a,
            resistance_ohm,
        })
    }
}

pub(super) struct ConductorRecovery {
    pub(super) eddy_cells: Vec<HarmonicEddyCell>,
    pub(super) solid_cells: Vec<HarmonicSolidCell>,
    pub(super) region_currents_a: Vec<(u32, C64)>,
    pub(super) solid_region_currents: Vec<SolidWindingCurrent>,
    pub(super) mean_eddy_power_w: f64,
    pub(super) mean_solid_power_w: f64,
}

pub(super) fn recover(
    model: &PlanarHarmonicWindings<'_, '_>,
    real: &[f64],
    imaginary: &[f64],
    currents: &[C64],
    omega: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ConductorRecovery> {
    let count = model.conduction.regions.len();
    let mut drives = vec![C64::ZERO; count];
    let mut targets = vec![C64::ZERO; count];
    for group in 0..count {
        poll(cancelled)?;
        if let Some(index) = model.solid.region_bindings[group] {
            let binding = model.solid.bindings[index];
            targets[group] = phasor(currents[binding.winding].scale(binding.signed_turns))?;
            drives[group] = phasor(
                currents[binding.winding].scale(model.solid.electric_per_current_v_m_a[group]),
            )?;
        }
    }
    // exp(-i omega t): E=drive+i omega P(A). The shared integrator computes
    // drive-P(increment)/dt. Imaginary A therefore gives real E; real A with
    // minus-imaginary drive gives minus-imaginary E, with identical Joule heat.
    let (electric_real, minus_electric_imaginary) = {
        let callback = RefCell::new(&mut *cancelled);
        let check = || (*callback.borrow_mut())();
        let depth = model.basis.operating_field().problem().thickness_m();
        let real_drives: Vec<_> = drives.iter().map(|v| v.re).collect();
        let minus_imaginary_drives: Vec<_> = drives.iter().map(|v| -v.im).collect();
        let real_targets: Vec<_> = targets.iter().map(|v| v.re).collect();
        let minus_imaginary_targets: Vec<_> = targets.iter().map(|v| -v.im).collect();
        let real_heat = model
            .conduction
            .heat_with_region_drives(
                imaginary,
                omega.recip(),
                depth,
                model.controls.current_balance_tolerance_a,
                &real_drives,
                &real_targets,
                &check,
            )
            .map_err(EddyError::from)?;
        let minus_imaginary_heat = model
            .conduction
            .heat_with_region_drives(
                real,
                omega.recip(),
                depth,
                model.controls.current_balance_tolerance_a,
                &minus_imaginary_drives,
                &minus_imaginary_targets,
                &check,
            )
            .map_err(EddyError::from)?;
        (real_heat, minus_imaginary_heat)
    };
    let mut recovered = ConductorRecovery {
        eddy_cells: vec![],
        solid_cells: vec![],
        region_currents_a: vec![],
        solid_region_currents: vec![],
        mean_eddy_power_w: 0.0,
        mean_solid_power_w: 0.0,
    };
    for (index, (source, cell)) in model
        .eddy_conductors
        .iter()
        .zip(&model.conduction.cells)
        .enumerate()
    {
        poll(cancelled)?;
        let re = &electric_real.cells[index];
        let minus_im = &minus_electric_imaginary.cells[index];
        let mean_power_w = finite(0.5 * re.power_w + 0.5 * minus_im.power_w)?;
        let electric_field_v_m = std::array::from_fn(|k| {
            C64::new(re.electric_field_v_m[k], -minus_im.electric_field_v_m[k])
        });
        let current_density_a_m2 = std::array::from_fn(|k| {
            C64::new(
                re.current_density_a_m2[k],
                -minus_im.current_density_a_m2[k],
            )
        });
        if let Some(binding) = model.solid.region_bindings[cell.group] {
            recovered.mean_solid_power_w = finite(recovered.mean_solid_power_w + mean_power_w)?;
            recovered.solid_cells.push(HarmonicSolidCell {
                source: source.clone(),
                binding: model.solid.bindings[binding],
                volume_m3: re.volume_m3,
                electric_field_v_m,
                current_density_a_m2,
                mean_power_w,
            });
        } else {
            recovered.mean_eddy_power_w = finite(recovered.mean_eddy_power_w + mean_power_w)?;
            recovered.eddy_cells.push(HarmonicEddyCell {
                source: source.clone(),
                volume_m3: re.volume_m3,
                electric_field_v_m,
                current_density_a_m2,
                mean_power_w,
            });
        }
    }
    for (group, (&(region, re), &(_, minus_im))) in electric_real
        .region_currents_a
        .iter()
        .zip(&minus_electric_imaginary.region_currents_a)
        .enumerate()
    {
        poll(cancelled)?;
        let current = phasor(C64::new(re, -minus_im))?;
        let residual = phasor(current - targets[group])?;
        if residual.abs() > model.controls.current_balance_tolerance_a {
            return Err(HarmonicWindingError::Balance(
                "conductor peak current disagrees with its original winding target",
            ));
        }
        if let Some(binding) = model.solid.region_bindings[group] {
            recovered.solid_region_currents.push(SolidWindingCurrent {
                binding: model.solid.bindings[binding],
                target_current_a: targets[group],
                current_a: current,
                residual_a: residual,
            });
        } else {
            recovered.region_currents_a.push((region, current));
        }
    }
    Ok(recovered)
}
