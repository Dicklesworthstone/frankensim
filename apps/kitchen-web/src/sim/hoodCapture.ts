/**
 * Three gas broilers under a wall canopy hood that stays where it is: the least air its fan must pull to catch all
 * their smoke, as the broilers stand further back under it and as the gap between their backs and the wall is closed.
 *
 * The rates are the laboratory's, not a plume model's. Three 3-foot gas underfired broilers, cooking, stood under a
 * 10-foot wall canopy hood; the exhaust was turned down until the hot air first escaped past the hood's edges, watched
 * by schlieren imaging: the least exhaust that caught it all is the setup's capture and containment rate. The
 * broilers were pushed back 0, 6 and 12 inches from the hood's front edge, the gap behind them shrinking from 16.5 to
 * 4.5 inches, with that gap open and with it closed by a sheet from the broilers' backs to the wall. All six rates
 * are the Food Service Technology Center's rear seal investigation (Swierczyna, "Applying state of the art commercial
 * kitchen ventilation technologies", 2018); the two open ones nearest the edge are also the ASHRAE Handbook's
 * (HVAC Applications, Kitchen Ventilation, Table 6, from Swierczyna et al. 2006: 790 and 680 L/s per metre of hood).
 *
 * What the fan's air costs: every cubic foot it pulls out is replaced by outdoor air that has to be warmed (or cooled)
 * to the kitchen's. Its sensible heat is the standard-air rule, 1.08 Btu/h per cfm per °F (0.075 lb/ft³ × 0.24
 * Btu/lb·°F × 60 min/h).
 */

/** The hood: 10 feet along the wall, its lower edge 78 inches above the floor. */
export const HOOD_FT = 10;
export const EDGE_IN = 78;
/** The building code's least exhaust for an unlisted wall canopy over heavy-duty appliances such as broilers
 *  (International Mechanical Code, as tabled in the California Energy Commission's Design Guide 1), a foot of hood. */
export const CODE_CFM_FT = 400;
export const CODE_CFM = CODE_CFM_FT * HOOD_FT;

/** Where the broilers stand: their front this far behind the hood's front edge (in), the hood's edge this far in
 *  front of their cooking surface (in), and the gap left between their backs and the wall (in). */
export const SETS = [
  { backIn: 0, cookIn: 6, gapIn: 16.5 },
  { backIn: 6, cookIn: 12, gapIn: 10.5 },
  { backIn: 12, cookIn: 18, gapIn: 4.5 },
] as const;
export type SetIndex = 0 | 1 | 2;

/** The measured capture and containment rates (cfm), by where the broilers stand: the gap open, then sealed. */
export const CC_CFM: Readonly<Record<"open" | "sealed", readonly [number, number, number]>> = {
  open: [5100, 4400, 4000],
  sealed: [3400, 2900, 2800],
};

export interface Setup {
  at: SetIndex;
  sealed: boolean;
}
export const SETUPS: readonly Setup[] = [0, 1, 2].flatMap((at) => [false, true].map((sealed) => ({ at: at as SetIndex, sealed })));

/** The least exhaust (cfm) that catches all the smoke of a setup. */
export const needCfm = (s: Setup): number => CC_CFM[s.sealed ? "sealed" : "open"][s.at];
/** Per foot of hood. */
export const perFoot = (cfm: number): number => cfm / HOOD_FT;
/** Whether a fan pulling `fan` cfm catches it all: at the measured rate or above. */
export const catches = (fan: number, s: Setup): boolean => fan >= needCfm(s);
/** The setups a fan catches, as a count: what moving the broilers and closing the gap buys. */
export const caughtCount = (fan: number): number => SETUPS.filter((s) => catches(fan, s)).length;

/** Litres a second in a cubic foot a minute. */
export const LPS_PER_CFM = 0.47194745;
/** Btu/h to warm `cfm` of outdoor air at `outF` to the kitchen's `inF` (°F): 1.08 · cfm · ΔT, none when it is warmer. */
export const SENSIBLE = 1.08;
export const warmBtuh = (cfm: number, outF: number, inF = 70): number => SENSIBLE * cfm * Math.max(0, inF - outF);
/** The morning the figure warms the air on. */
export const COLD_OUT_F = 40;
