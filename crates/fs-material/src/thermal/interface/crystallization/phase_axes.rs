//! Original kinetic fractions as live coordinates of a contact resistance law.

use super::{
    CrystallizationBulkConduction, CrystallizationInterfaceContact, CrystallizationInterfaceError,
    CrystallizationPropertyError, CrystallizationSourcePair, IsobaricCrystallizationParcel,
    MaterialPropertySelection, QueryPoint, ThermalInterfacePoint, ThermalMaterialError, invalid,
};
use fs_qty::{Dims, QuantitySpec};

/// Optional crystalline coordinates in the original interface's A/B order.
///
/// All names must identify distinct, existing dimensional scalar axes in the
/// supplied context. Mass coordinates use the absolute crystalline MASS fraction,
/// never X/Xmax or one minus the liquid fraction. Volume coordinates use the
/// same original phase densities as the corresponding kinetic parcel.
///
/// Applying a uniform bulk inventory to an interface is an explicit constitutive
/// assumption: surface enrichment, crystal texture, morphology and interfacial
/// layer history need their own source/model. No resistance mixing rule follows
/// from these coordinates alone.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CrystallizationInterfacePhaseAxes<'a> {
    pub crystalline_mass_fraction: [Option<&'a str>; 2],
    pub crystalline_volume_fraction: [Option<&'a str>; 2],
}

impl<'a> CrystallizationInterfaceContact<'a> {
    /// Prepare contact resistance at both actual temperatures AND crystal states.
    ///
    /// Fractions are refreshed from every candidate, including partial remelting,
    /// before querying the pinned original resistance claim. Failure of its joint
    /// temperature/fraction support aborts the entire original heat transaction.
    /// Frozen input fractions cannot stand in for the solved endpoint inventory.
    ///
    /// Optional bulk conduction supplies the existing series-resistance closure,
    /// with its own live original conductivity sources. Geometry, contact pressure
    /// and separation retain the same fixed conventions as the other constructors.
    /// Either original chart or exact-Cp pairs are retained without resampling;
    /// volume coordinates require their original density sources.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_crystallinity(
        point: ThermalInterfacePoint<'a>,
        pairs: [impl Into<CrystallizationSourcePair<'a>>; 2],
        parcels: [&IsobaricCrystallizationParcel<'_, '_>; 2],
        area_m2: f64,
        fixed_pressure_axis: Option<&str>,
        fixed_separation_axis: Option<&str>,
        phase_axes: CrystallizationInterfacePhaseAxes<'_>,
        selection: MaterialPropertySelection,
        bulk: Option<CrystallizationBulkConduction<'a>>,
    ) -> Result<Self, CrystallizationInterfaceError> {
        let pairs = pairs.map(Into::into);
        let coordinates = PhaseCoordinates {
            mass: phase_axes
                .crystalline_mass_fraction
                .map(|axis| axis.map(str::to_owned)),
            volume: phase_axes
                .crystalline_volume_fraction
                .map(|axis| axis.map(str::to_owned)),
            pairs,
        };
        let mut contact = Self::prepare(
            point,
            pairs,
            parcels,
            area_m2,
            fixed_pressure_axis,
            fixed_separation_axis,
            selection,
            Some(coordinates),
        )?;
        if bulk
            .as_ref()
            .is_some_and(|bulk| contact.bindings != bulk.bindings)
        {
            return Err(invalid(
                "bulk conduction and interface require the same ordered original laws and pressures",
            ));
        }
        contact.bulk = bulk;
        contact.evaluate(parcels[0], parcels[1])?;
        Ok(contact)
    }
}

#[derive(Clone, Debug)]
pub(super) struct PhaseCoordinates<'a> {
    mass: [Option<String>; 2],
    volume: [Option<String>; 2],
    pairs: [CrystallizationSourcePair<'a>; 2],
}

impl PhaseCoordinates<'_> {
    pub(super) fn names(&self) -> impl Iterator<Item = &str> {
        self.mass
            .iter()
            .chain(&self.volume)
            .filter_map(Option::as_deref)
    }

    pub(super) fn admit(&self, context: &QueryPoint) -> Result<(), CrystallizationInterfaceError> {
        for axis in self.names() {
            if context.axis_quantities().get(axis) != Some(&QuantitySpec::dimensional(Dims::NONE)) {
                return Err(invalid(
                    "contact crystallinity requires typed dimensional scalar fraction axes",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn update(
        &self,
        mut query: QueryPoint,
        parcels: [&IsobaricCrystallizationParcel<'_, '_>; 2],
    ) -> Result<QueryPoint, CrystallizationInterfaceError> {
        for (side, parcel) in parcels.into_iter().enumerate() {
            let state = parcel.state();
            let coordinate =
                |error: fs_matdb::MatDbError| ThermalMaterialError::Coordinate(error.to_string());
            if let Some(axis) = &self.mass[side] {
                query = query
                    .with(axis, state.crystalline_mass_fraction())
                    .map_err(coordinate)?;
            }
            if let Some(axis) = &self.volume[side] {
                let volume = self.pairs[side]
                    .volume_state(state.temperature_k(), state.crystalline_mass_fraction())
                    .map_err(CrystallizationPropertyError::from)?;
                query = query
                    .with(axis, volume.crystalline_volume_fraction)
                    .map_err(coordinate)?;
            }
        }
        Ok(query)
    }
}
