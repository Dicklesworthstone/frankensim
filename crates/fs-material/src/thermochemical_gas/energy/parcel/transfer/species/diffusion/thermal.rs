//! Source-conjugate Maxwell--Stefan, Soret and reciprocal heat transport.
//!
//! The mass-average frame has sum(M_i z_i)=0. With beta=1/T, the original
//! NASA calorics define h_i*=Delta(h_i/T-s_i^0)/Delta(beta). The chemical
//! effort after subtracting h_i* Delta(beta) is R Delta(ln p_i). Species
//! heats of transport H_i [J/kg] add -M_i H_i Delta(beta); the reciprocal
//! reduced heat is -k grad(T)+sum(H_i M_i z_i). Centering H in the mass
//! frame removes its physically irrelevant common offset BEFORE reduction.
//!
//! The existing inventory owner already transports donor partial enthalpy.
//! Therefore its additional heat is Fourier + reciprocal heat +
//! sum(z_i (h_i*-h_i,donor)). Original NASA enthalpy/entropy INCREMENTS form
//! this correction without subtracting arbitrary formation-energy datums.
//! Omitting the correction would break the finite-distance entropy pairing.
//!
//! All coefficients are explicit material claims, including conductivity at
//! ZERO diffusive flux. No heat of transport is inferred from mass or labels.
//! Claims are pinned once, then refreshed at both traces and the central face.
//! This bounded ideal-gas closure owns only temperature/pressure dependence;
//! unowned composition, field, stress or history dependence is refused.
//! Positive common species support is required for logarithmic forces. Exact
//! common zeros remain absent; one-sided zeros refuse, without trace floors.
//!
//! Continuum convention: Van-Brunt, Farrell and Monroe,
//! https://arxiv.org/abs/2109.05082, Eq. 2.8 and Section 3. The divided caloric
//! pairing matches fs-gas's source-conjugate molecular face. A local positive
//! entropy rate does not replace the network's actual finite-step entropy gate.

/// All faces and other material paths use one conservative snapshot update.
pub mod network;

use super::super::entropy::log_ratio;
use super::{
    BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial, ClosedGasParcel,
    GasDiffusionGeometry, GasEnergyError, MAX_DIFFUSION_SPECIES, SpeciesDiffusionFlux, SpeciesId,
    catalog, diffusion_error, fractions, invalid, poll, positive, source_error,
};
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
    resolve_material_state_point,
};
use crate::thermochemical_gas::{ThermochemicalGasModel, energy::updated_point};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_math::dd::Dd;
use fs_qty::{
    AmountConcentration, Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::transport::maxwell_stefan::{
    MaxwellStefanSolveV1, MaxwellStefanV1, MaxwellStefanWorkspaceV1,
};
use fs_thermochem::{Nasa9StandardStateModelV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

/// Mass-specific heat of transport, NOT a species enthalpy or heat capacity.
pub const HEAT_OF_TRANSPORT_PROPERTY: &str = "heat_of_transport";
pub const HEAT_OF_TRANSPORT_DIMS: Dims = Dims([2, 0, -2, 0, 0, 0]);

/// Explicit source binding. The reference point pins one claim; it is not the
/// physical face temperature. Canonical typed temperature/pressure are required.
#[derive(Clone, Debug)]
pub struct GasThermodiffusionMaterial<'a> {
    pub material: &'a MaterialCard,
    pub context: &'a QueryPoint,
    pub selection: MaterialPropertySelection,
}

/// One signed heat-of-transport source for every original species.
#[derive(Clone, Debug)]
pub struct SpeciesHeatOfTransportMaterial<'a> {
    pub species: SpeciesId,
    pub source: GasThermodiffusionMaterial<'a>,
}

/// Explicit face approximation limits and original-row numerical admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasThermodiffusionSolve {
    pub maximum_temperature_difference_k: f64,
    pub maximum_pressure_difference_pa: f64,
    pub diffusion: MaxwellStefanSolveV1,
}
impl GasThermodiffusionSolve {
    fn validate(self) -> Result<(), GasEnergyError> {
        if ![
            self.maximum_temperature_difference_k,
            self.maximum_pressure_difference_pa,
        ]
        .iter()
        .all(|x| x.is_finite() && *x >= 0.0)
            || ![
                self.diffusion.relative_tolerance,
                self.diffusion.absolute_gradient_tolerance_per_m,
                self.diffusion.absolute_mass_flux_tolerance_kg_m2_s,
            ]
            .iter()
            .all(|x| positive(*x))
        {
            return Err(invalid("invalid thermodiffusion support or solve controls"));
        }
        Ok(())
    }
}

/// Same immutable card and selected claim at three actual physical points.
#[derive(Clone, Debug, PartialEq)]
pub struct GasThermodiffusionSample {
    pub value_si: f64,
    pub face: ResolvedMaterialStatePoint,
    pub endpoints: [ResolvedMaterialStatePoint; 2],
}
#[derive(Clone, Debug, PartialEq)]
pub struct GasThermodiffusionPairSample {
    pub species: [SpeciesId; 2],
    pub sample: GasThermodiffusionSample,
}
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesHeatOfTransportSample {
    pub species: SpeciesId,
    /// Original source value, before the physically irrelevant gauge is removed.
    pub sample: GasThermodiffusionSample,
    pub centered_heat_j_kg: f64,
}

/// Signed left-to-right rates, source receipts and original-equation checks.
#[derive(Clone, Debug, PartialEq)]
pub struct GasThermodiffusionFlux {
    pub species: Vec<SpeciesDiffusionFlux>,
    pub pairs: Vec<GasThermodiffusionPairSample>,
    pub heats: Vec<SpeciesHeatOfTransportSample>,
    pub conductivity: GasThermodiffusionSample,
    pub face_temperature_k: f64,
    pub face_pressure_pa: f64,
    pub fourier_heat_flux_w_m2: f64,
    pub reciprocal_heat_flux_w_m2: f64,
    /// Converts existing donor partial-enthalpy carriage to divided calorics.
    pub enthalpy_correction_flux_w_m2: f64,
    /// Apply this ONCE, in addition to the existing species enthalpy transport.
    pub additional_heat_rate_left_to_right_w: f64,
    pub maximum_gradient_residual_per_m: f64,
    pub mass_frame_residual_kg_m2_s: f64,
    pub drag_entropy_production_w_m3_k: Option<f64>,
    /// Paired semidiscrete entropy rate; not a finite-step acceptance certificate.
    pub entropy_production_w_k: Option<f64>,
}

#[derive(Clone, Debug)]
struct PinnedScalar<'a> {
    card: &'a MaterialCard,
    property: &'static str,
    reference: ResolvedMaterialStatePoint,
}
impl<'a> PinnedScalar<'a> {
    fn new(
        source: GasThermodiffusionMaterial<'a>,
        property: &'static str,
        dims: Dims,
        admissibility: ScalarAdmissibility,
    ) -> Result<Self, GasEnergyError> {
        for (axis, kind) in [
            ("temperature", QuantityKind::AbsoluteTemperature),
            ("pressure", QuantityKind::Pressure),
        ] {
            let expected = QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
            if source.context.axis_quantities().get(axis) != Some(&expected) {
                return Err(invalid(
                    "thermodiffusion sources require canonical typed T/p axes",
                ));
            }
        }
        let requirement = ScalarPropertyRequirement::try_new(property, dims, admissibility)
            .map_err(source_error)?;
        let reference = resolve_material_state_point(
            source.material,
            source.context,
            &[requirement],
            source.selection,
        )
        .map_err(source_error)?;
        let value = reference
            .property(property)
            .ok_or(invalid("missing transport property"))?;
        if !value.has_only_fixed_external_conditions(Some("temperature"), false)
            || !value.source_axes().iter().all(|(axis, dims)| {
                axis == "temperature" && *dims == Temperature::DIMS
                    || axis == "pressure" && *dims == Pressure::DIMS
            })
            || !value
                .answer()
                .evidence
                .model
                .validity
                .bounds()
                .keys()
                .all(|axis| axis == "temperature" || axis == "pressure")
        {
            return Err(invalid(
                "thermodiffusion source depends on an unowned coordinate",
            ));
        }
        Ok(Self {
            card: source.material,
            property,
            reference,
        })
    }

    fn at(&self, t: f64, p: f64) -> Result<ResolvedMaterialStatePoint, GasEnergyError> {
        let point = updated_point(&self.reference.query(), t, p)?;
        self.reference
            .at_query_point(self.card, &point)
            .map_err(source_error)
    }

    fn sample(
        &self,
        temperatures: [f64; 2],
        pressures: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<GasThermodiffusionSample, GasEnergyError> {
        poll(cancelled)?;
        let left = self.at(temperatures[0], pressures[0])?;
        poll(cancelled)?;
        let right = self.at(temperatures[1], pressures[1])?;
        poll(cancelled)?;
        let face = self.at(midpoint(temperatures), midpoint(pressures))?;
        let value_si = face
            .property(self.property)
            .ok_or(invalid("missing refreshed transport property"))?
            .value_si();
        Ok(GasThermodiffusionSample {
            value_si,
            face,
            endpoints: [left, right],
        })
    }
}

/// Prepared original species catalog, friction owner and pinned material laws.
/// Preparation and face evaluation do not alter any parcel or source claim.
#[derive(Clone, Debug)]
pub struct GasThermodiffusionModel<'a> {
    calorics: Vec<Nasa9StandardStateModelV1>,
    diffusion: MaxwellStefanV1,
    pairs: Vec<(usize, usize, PinnedScalar<'a>)>,
    heats: Vec<PinnedScalar<'a>>,
    conductivity: PinnedScalar<'a>,
}

/// Reusable matrix scratch; no evolving composition or physical state.
#[derive(Clone, Debug)]
pub struct GasThermodiffusionWorkspace {
    binary: Vec<f64>,
    diffusion: MaxwellStefanWorkspaceV1,
}
impl GasThermodiffusionWorkspace {
    pub fn workspace_bytes(&self) -> usize {
        self.binary.capacity() * core::mem::size_of::<f64>() + self.diffusion.workspace_bytes()
    }
}

impl<'a> GasThermodiffusionModel<'a> {
    /// Pin a complete full-catalog pair matrix and one H_i claim per species.
    /// `conductivity` MUST represent conductivity at zero diffusive species flux,
    /// the Schur complement paired with the authored heats of transport.
    pub fn new(
        reference: &ThermochemicalGasModel,
        pairs: &[BinaryDiffusionMaterial<'a>],
        conductivity: GasThermodiffusionMaterial<'a>,
        heats: &[SpeciesHeatOfTransportMaterial<'a>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, GasEnergyError> {
        poll(&mut cancelled)?;
        let mut calorics = reference.thermodynamics().components().to_vec();
        calorics.sort_by(|a, b| a.species().cmp(b.species()));
        let n = calorics.len();
        if !(1..=MAX_DIFFUSION_SPECIES).contains(&n)
            || pairs.len() != n * (n - 1) / 2
            || heats.len() != n
        {
            return Err(invalid(
                "thermodiffusion requires bounded complete species/pair sources",
            ));
        }
        let diffusion = MaxwellStefanV1::new(
            calorics
                .iter()
                .map(|s| (s.species().clone(), s.molar_mass()))
                .collect(),
        )
        .map_err(diffusion_error)?;
        let mut selected_pairs = Vec::with_capacity(pairs.len());
        for i in 0..n {
            for j in i + 1..n {
                poll(&mut cancelled)?;
                let mut candidates = pairs.iter().filter(|pair| {
                    pair.species
                        == [
                            diffusion.species()[i].clone(),
                            diffusion.species()[j].clone(),
                        ]
                        || pair.species
                            == [
                                diffusion.species()[j].clone(),
                                diffusion.species()[i].clone(),
                            ]
                });
                let pair = candidates
                    .next()
                    .ok_or(invalid("missing thermodiffusion pair"))?;
                if candidates.next().is_some() {
                    return Err(invalid("duplicate thermodiffusion pair"));
                }
                let source = GasThermodiffusionMaterial {
                    material: pair.material,
                    context: pair.context,
                    selection: pair.selection.clone(),
                };
                selected_pairs.push((
                    i,
                    j,
                    PinnedScalar::new(
                        source,
                        BINARY_DIFFUSIVITY_PROPERTY,
                        BINARY_DIFFUSIVITY_DIMS,
                        ScalarAdmissibility::StrictlyPositive,
                    )?,
                ));
            }
        }
        let mut selected_heats = Vec::with_capacity(n);
        for id in diffusion.species() {
            poll(&mut cancelled)?;
            let mut candidates = heats.iter().filter(|heat| &heat.species == id);
            let heat = candidates
                .next()
                .ok_or(invalid("missing heat-of-transport source"))?;
            if candidates.next().is_some() {
                return Err(invalid("duplicate heat-of-transport source"));
            }
            selected_heats.push(PinnedScalar::new(
                heat.source.clone(),
                HEAT_OF_TRANSPORT_PROPERTY,
                HEAT_OF_TRANSPORT_DIMS,
                ScalarAdmissibility::Finite,
            )?);
        }
        let conductivity = PinnedScalar::new(
            conductivity,
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            ScalarAdmissibility::NonNegative,
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            calorics,
            diffusion,
            pairs: selected_pairs,
            heats: selected_heats,
            conductivity,
        })
    }

    pub fn species(&self) -> &[SpeciesId] {
        self.diffusion.species()
    }

    pub fn workspace(&self) -> Result<GasThermodiffusionWorkspace, GasEnergyError> {
        let n = self.species().len();
        Ok(GasThermodiffusionWorkspace {
            binary: vec![0.0; n * n],
            diffusion: self.diffusion.workspace().map_err(diffusion_error)?,
        })
    }

    fn composition(
        &self,
        left: &ClosedGasParcel,
        right: &ClosedGasParcel,
        controls: GasThermodiffusionSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<
        (
            [f64; 2],
            [f64; 2],
            [f64; MAX_DIFFUSION_SPECIES],
            [f64; MAX_DIFFUSION_SPECIES],
            [f64; MAX_DIFFUSION_SPECIES],
        ),
        GasEnergyError,
    > {
        poll(cancelled)?;
        controls.validate()?;
        let temperatures = [
            left.state().state().temperature,
            right.state().state().temperature,
        ];
        let pressures = [
            left.state().state().pressure,
            right.state().state().pressure,
        ];
        if !temperatures.iter().chain(&pressures).all(|x| positive(*x))
            || (temperatures[1] - temperatures[0]).abs() > controls.maximum_temperature_difference_k
            || (pressures[1] - pressures[0]).abs() > controls.maximum_pressure_difference_pa
        {
            return Err(invalid(
                "thermodiffusion endpoint T/p exceeds authored face support",
            ));
        }
        // The existing join also checks pure transport, water and source-pack identities.
        let joined = catalog(left, right)?;
        for source in joined.thermodynamics().components() {
            poll(cancelled)?;
            let index = self
                .species()
                .binary_search(source.species())
                .map_err(|_| invalid("thermodiffusion endpoint adds an unbound species"))?;
            if source != &self.calorics[index] {
                return Err(invalid(
                    "thermodiffusion endpoint replaced its original caloric source",
                ));
            }
        }
        let left_x = fractions(left);
        let right_x = fractions(right);
        let mut xl = [0.0; MAX_DIFFUSION_SPECIES];
        let mut xr = xl;
        let mut x = xl;
        for (i, id) in self.species().iter().enumerate() {
            xl[i] = left_x.get(id).copied().unwrap_or(0.0);
            xr[i] = right_x.get(id).copied().unwrap_or(0.0);
            if (xl[i] == 0.0) != (xr[i] == 0.0) {
                return Err(invalid(
                    "thermodiffusion requires common positive species support",
                ));
            }
            x[i] = midpoint([xl[i], xr[i]]);
            if xl[i] > 0.0 && (!positive(xr[i]) || !positive(x[i])) {
                return Err(invalid(
                    "thermodiffusion trace composition is unrepresentable",
                ));
            }
        }
        // Only the largest component absorbs normalization roundoff, never a trace.
        let n = self.species().len();
        let pivot = largest(&x[..n]);
        let mut rest = Dd::ZERO;
        for (i, value) in x[..n].iter().enumerate() {
            if i != pivot {
                rest = rest + Dd::from_f64(*value);
            }
        }
        let closed = checked(Dd::from_f64(1.0) - rest)?;
        if !positive(closed) || (closed - x[pivot]).abs() > 256.0 * f64::EPSILON {
            return Err(invalid("thermodiffusion face composition cannot be closed"));
        }
        x[pivot] = closed;
        Ok((temperatures, pressures, xl, xr, x))
    }

    /// Evaluate a physical face without allocating another dense solve matrix.
    /// Positive area and distance are required; omit closed connections entirely.
    #[allow(clippy::too_many_arguments)]
    pub fn face_flux(
        &self,
        left: &ClosedGasParcel,
        right: &ClosedGasParcel,
        geometry: GasDiffusionGeometry,
        controls: GasThermodiffusionSolve,
        work: &mut GasThermodiffusionWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasThermodiffusionFlux, GasEnergyError> {
        let n = self.species().len();
        if !positive(geometry.area_m2)
            || !positive(geometry.distance_m)
            || work.binary.len() != n * n
        {
            return Err(invalid("thermodiffusion geometry or workspace dimensions"));
        }
        let (temperatures, pressures, xl, xr, x) =
            self.composition(left, right, controls, &mut cancelled)?;
        let t = midpoint(temperatures);
        let p = midpoint(pressures);
        let delta_t = temperatures[1] - temperatures[0];
        let delta_beta = reciprocal_difference(temperatures)?;
        let pivot = largest(&x[..n]);
        let masses = self.diffusion.molar_masses();
        let mut mean_mass = Dd::ZERO;
        for i in 0..n {
            mean_mass = mean_mass + Dd::from_f64(x[i]) * Dd::from_f64(masses[i].value());
        }
        let mean_mass = checked(mean_mass)?;
        if !positive(mean_mass) {
            return Err(invalid("thermodiffusion mean molar mass"));
        }

        work.binary.fill(0.0);
        let mut pairs = Vec::with_capacity(self.pairs.len());
        for (i, j, source) in &self.pairs {
            let sample = source.sample(temperatures, pressures, &mut cancelled)?;
            work.binary[i * n + j] = sample.value_si;
            work.binary[j * n + i] = sample.value_si;
            pairs.push(GasThermodiffusionPairSample {
                species: [self.species()[*i].clone(), self.species()[*j].clone()],
                sample,
            });
        }
        let conductivity = self
            .conductivity
            .sample(temperatures, pressures, &mut cancelled)?;
        let mut heats = Vec::with_capacity(n);
        for (id, source) in self.species().iter().zip(&self.heats) {
            heats.push(SpeciesHeatOfTransportSample {
                species: id.clone(),
                sample: source.sample(temperatures, pressures, &mut cancelled)?,
                centered_heat_j_kg: 0.0,
            });
        }
        let reference_heat = Dd::from_f64(heats[pivot].sample.value_si);
        let mut mean_heat = Dd::ZERO;
        for i in 0..n {
            mean_heat = mean_heat
                + Dd::from_f64(x[i]) * Dd::from_f64(masses[i].value()) / Dd::from_f64(mean_mass)
                    * (Dd::from_f64(heats[i].sample.value_si) - reference_heat);
        }
        for heat in &mut heats {
            heat.centered_heat_j_kg =
                checked(Dd::from_f64(heat.sample.value_si) - reference_heat - mean_heat)?;
        }

        let mut divided_corrections = [Dd::ZERO; MAX_DIFFUSION_SPECIES];
        let mut efforts = [0.0; MAX_DIFFUSION_SPECIES];
        for (i, model) in self.calorics.iter().enumerate() {
            poll(&mut cancelled)?;
            if x[i] == 0.0 {
                continue;
            }
            let relative = if temperatures[0] == temperatures[1] {
                Dd::ZERO
            } else {
                Dd::from_f64(
                    model
                        .molar_enthalpy_change(
                            Temperature::new(t),
                            Temperature::new(temperatures[1]),
                        )
                        .map_err(source_error)?
                        .value(),
                ) / Dd::from_f64(temperatures[1])
                    - Dd::from_f64(
                        model
                            .molar_enthalpy_change(
                                Temperature::new(t),
                                Temperature::new(temperatures[0]),
                            )
                            .map_err(source_error)?
                            .value(),
                    ) / Dd::from_f64(temperatures[0])
                    - Dd::from_f64(
                        model
                            .molar_entropy_change(
                                Temperature::new(temperatures[0]),
                                Temperature::new(temperatures[1]),
                            )
                            .map_err(source_error)?
                            .value(),
                    )
            };
            let correction = if delta_beta == 0.0 {
                Dd::ZERO
            } else {
                relative / Dd::from_f64(delta_beta)
            };
            checked(correction)?;
            divided_corrections[i] = correction;
            let remainder = relative - correction * Dd::from_f64(delta_beta);
            efforts[i] = checked(
                remainder
                    + Dd::from_f64(R)
                        * (Dd::from_f64(log_ratio(xr[i], xl[i]))
                            + Dd::from_f64(log_ratio(pressures[1], pressures[0]))),
            )?;
        }
        let reference_effort = efforts[pivot];
        for i in 0..n {
            efforts[i] = if i == pivot || x[i] == 0.0 {
                0.0
            } else {
                checked(
                    Dd::from_f64(efforts[i])
                        - Dd::from_f64(masses[i].value()) / Dd::from_f64(masses[pivot].value())
                            * Dd::from_f64(reference_effort),
                )?
            };
        }
        let mut mean_effort = Dd::ZERO;
        for i in 0..n {
            mean_effort = mean_effort + Dd::from_f64(x[i]) * Dd::from_f64(efforts[i]);
        }
        mean_effort = mean_effort / Dd::from_f64(mean_mass);
        let mut forces = [0.0; MAX_DIFFUSION_SPECIES];
        for i in 0..n {
            if x[i] == 0.0 {
                continue;
            }
            let effort = Dd::from_f64(efforts[i])
                - Dd::from_f64(masses[i].value()) * mean_effort
                - Dd::from_f64(masses[i].value())
                    * Dd::from_f64(heats[i].centered_heat_j_kg)
                    * Dd::from_f64(delta_beta);
            forces[i] = checked(
                Dd::from_f64(x[i]) * effort / Dd::from_f64(R) / Dd::from_f64(geometry.distance_m),
            )?;
        }
        let evaluated = self
            .diffusion
            .evaluate_driving_forces(
                &x[..n],
                &forces[..n],
                AmountConcentration::new((p / R) / t),
                &work.binary,
                controls.diffusion,
                &mut work.diffusion,
                &mut cancelled,
            )
            .map_err(diffusion_error)?;
        let fourier = Dd::ZERO
            - Dd::from_f64(conductivity.value_si) * Dd::from_f64(delta_t)
                / Dd::from_f64(geometry.distance_m);
        let mut reciprocal = Dd::ZERO;
        let mut correction = Dd::ZERO;
        let mut species = Vec::with_capacity(n);
        for i in 0..n {
            poll(&mut cancelled)?;
            let mass_flux = evaluated.mass_fluxes()[i];
            let molar_flux = evaluated.molar_fluxes()[i];
            reciprocal =
                reciprocal + Dd::from_f64(heats[i].centered_heat_j_kg) * Dd::from_f64(mass_flux);
            if mass_flux != 0.0 {
                let donor_t = temperatures[usize::from(mass_flux < 0.0)];
                let donor_correction = Dd::from_f64(
                    self.calorics[i]
                        .molar_enthalpy_change(Temperature::new(donor_t), Temperature::new(t))
                        .map_err(source_error)?
                        .value(),
                ) + divided_corrections[i];
                correction = correction + Dd::from_f64(molar_flux) * donor_correction;
            }
            let mass_rate = checked(Dd::from_f64(mass_flux) * Dd::from_f64(geometry.area_m2))?;
            if mass_flux != 0.0 && mass_rate == 0.0 {
                return Err(invalid(
                    "thermodiffusion area would discard a positive trace flux",
                ));
            }
            species.push(SpeciesDiffusionFlux {
                species: self.species()[i].clone(),
                molar_flux_mol_m2_s: molar_flux,
                mass_rate_left_to_right_kg_s: mass_rate,
            });
        }
        let drag = evaluated.drag_entropy_production_w_m3_k();
        let entropy = drag.and_then(|drag| {
            checked(
                (Dd::from_f64(drag) * Dd::from_f64(geometry.distance_m)
                    - fourier * Dd::from_f64(delta_t)
                        / Dd::from_f64(temperatures[0])
                        / Dd::from_f64(temperatures[1]))
                    * Dd::from_f64(geometry.area_m2),
            )
            .ok()
        });
        let report = GasThermodiffusionFlux {
            species,
            pairs,
            heats,
            conductivity,
            face_temperature_k: t,
            face_pressure_pa: p,
            fourier_heat_flux_w_m2: checked(fourier)?,
            reciprocal_heat_flux_w_m2: checked(reciprocal)?,
            enthalpy_correction_flux_w_m2: checked(correction)?,
            additional_heat_rate_left_to_right_w: checked(
                (fourier + reciprocal + correction) * Dd::from_f64(geometry.area_m2),
            )?,
            maximum_gradient_residual_per_m: evaluated.maximum_gradient_residual_per_m(),
            mass_frame_residual_kg_m2_s: evaluated.mass_frame_residual_kg_m2_s(),
            drag_entropy_production_w_m3_k: drag,
            entropy_production_w_k: entropy,
        };
        poll(&mut cancelled)?;
        Ok(report)
    }

    /// Admit the actual combined endpoint against the same source claims and
    /// original caloric catalog. This does not impose separate species ordering:
    /// Soret and multicomponent friction can drive individual species uphill.
    pub fn admit_endpoints(
        &self,
        left: &ClosedGasParcel,
        right: &ClosedGasParcel,
        controls: GasThermodiffusionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), GasEnergyError> {
        let (temperatures, pressures, ..) =
            self.composition(left, right, controls, &mut cancelled)?;
        for (_, _, source) in &self.pairs {
            source.sample(temperatures, pressures, &mut cancelled)?;
        }
        for source in &self.heats {
            source.sample(temperatures, pressures, &mut cancelled)?;
        }
        self.conductivity
            .sample(temperatures, pressures, &mut cancelled)?;
        poll(&mut cancelled)
    }
}

fn midpoint(values: [f64; 2]) -> f64 {
    if values[0] == values[1] {
        values[0]
    } else {
        0.5 * values[0] + 0.5 * values[1]
    }
}
fn largest(values: &[f64]) -> usize {
    (0..values.len())
        .max_by(|&a, &b| values[a].total_cmp(&values[b]))
        .expect("nonempty catalog")
}
fn reciprocal_difference(t: [f64; 2]) -> Result<f64, GasEnergyError> {
    if t[0] == t[1] {
        return Ok(0.0);
    }
    let (lo, hi, sign) = if t[1] > t[0] {
        (t[0], t[1], -1.0)
    } else {
        (t[1], t[0], 1.0)
    };
    let magnitude = ((hi - lo) / hi) / lo;
    if !positive(magnitude) {
        return Err(invalid("unrepresentable reciprocal temperature difference"));
    }
    Ok(sign * magnitude)
}
fn checked(value: Dd) -> Result<f64, GasEnergyError> {
    let result = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && result.is_finite() {
        Ok(result)
    } else {
        Err(invalid(
            "unrepresentable source-conjugate thermodiffusion arithmetic",
        ))
    }
}

#[cfg(test)]
mod tests;
