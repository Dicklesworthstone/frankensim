/**
 * CHAPTER 1 — "THE HARDWARE": the showroom turntable (#equipment).
 *
 * Every archetype the kitchens are drawn from (the 67 of the equipment catalog, the 11 more the exemplar kitchens
 * place, and the stove chapter's hotplate) is built here by the procedural model library (models/all.ts) at close-up
 * detail, to its archetype's extents, and shown on a studio turntable: a dark softbox environment for the stainless
 * to reflect, a warm key light with soft shadows, a platter that turns with the model.
 *
 * The controls work the models themselves:
 *   Animate joints   cycles each model's real joints and live quantities (models/*.ts document them): doors swing,
 *                    baskets lift, lids open, fans and blades spin, gas valves and their knobs turn, foods sear.
 *   Thermal preview  paints the nominal set point of each archetype's active zones (equipment_dimensions.tsv), named
 *                    by its manifest's thermal slots, with an illustrative falloff onto the surfaces round them.
 *   Wireframe X-ray  a see-through shell with the model's feature edges; moving joints in cyan, the body in amber.
 *   SI cage          the archetype's W × D × H envelope with dimension lines in millimetres, and its manifest's
 *                    interaction sockets (grasps, stance, pan seats).
 *   Auto-spin        turns the platter; it pauses while you orbit or look close.
 * The spec card is read off the built model: triangles, joints, sockets, extents.
 *
 * Archetype data below is generated from assets/kitchens: equipment/catalog_index.json (ids, extents, MatDB cards),
 * equipment/*.scene3d.v1.json and exemplars/*\/*.scene3d.v1.json (sockets, thermal slots, extents of the exemplar-only
 * archetypes) and specs/equipment_dimensions.tsv (names, utilities, active zones). Frames: sockets and zones are in
 * the fixture frame of models/shop.ts (x across from the left, y up, z from the back at -d to the front at 0).
 */

import * as THREE from "three";
import { LineMaterial } from "three/examples/jsm/lines/LineMaterial.js";
import { LineSegments2 } from "three/examples/jsm/lines/LineSegments2.js";
import { LineSegmentsGeometry } from "three/examples/jsm/lines/LineSegmentsGeometry.js";
import { FLAGSHIP_EQUIPMENT, categorizeCatalogEntry, type CatalogEntry, type EquipmentTabId } from "../sim/exemplars";
import { createScrollSafeOrbit, glOf, observeInView, shaderGate } from "./kit";
import { buildFixture, type Model, type Piece } from "./models/all";
import "./models/electricHob";
import { palette } from "./models/materials";
import { makeRenderer, studioEnvironment } from "./models/stage";

/* ------------------------------------------------------------------------------------------------------------ */
/* Archetypes                                                                                                     */
/* ------------------------------------------------------------------------------------------------------------ */

/** An active zone: a box (min xyz, max xyz) or a sphere (centre xyz, radius), at a nominal set point (°C). */
type Zone = [kind: "b" | "s", a: number, b: number, c: number, d: number, e: number, f: number, tempC: number, key: string];
/** An interaction socket from the manifest: id, kind, position in the fixture frame. */
type Socket = [id: string, kind: string, x: number, y: number, z: number];

interface Arch {
  id: string;
  /** Short label for chips. */
  label: string;
  /** Display name from equipment_dimensions.tsv. */
  name: string;
  /** Width, depth, height in metres. */
  ext: [number, number, number];
  /** MatDB surface card from the catalog ("" for archetypes outside it). */
  card: string;
  /** Utility ports from equipment_dimensions.tsv. */
  util: string;
  cls: string;
  fam: string;
  /** "catalog", the exemplar kitchen that places it, or "stove". */
  src: string;
  sockets: Socket[];
  /** Thermal slots: [slot id, bound channel]. */
  slots: [string, string][];
  zones: Zone[];
}

const ARCH: Arch[] = [
  {"id":"NSF-RANGE-GAS-6B","label":"6-Burner Gas Range","name":"36 in 6-Burner Open-Top Gas Range with Standard Oven","ext":[0.914,0.864,1.473],"card":"surfaces/cookware-cast-iron-seasoned","util":"gas:NG:0.75in:61500W;exhaust:TypeI","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.457,0.704,-0.01],["grip_L","left_hand_grasp",0.2559,0.704,-0.01],["grip_R","right_hand_grasp",0.6581,0.704,-0.01],["stance_front","stance_anchor",0.457,0,0.45],["pan_seat_front_center","pan_seat",0.457,0.914,-0.24]],"slots":[["burner_ring_0","range.burner_heat_w"],["cast_iron_grates","range.grate_temp_k"]],"zones":[["s",0.18,0.914,-0.22,0.11,0,0,310,"burners"],["s",0.457,0.914,-0.22,0.11,0,0,310,"burners"],["s",0.734,0.914,-0.22,0.11,0,0,310,"burners"],["s",0.18,0.914,-0.58,0.11,0,0,310,"burners"],["s",0.457,0.914,-0.58,0.11,0,0,310,"burners"],["s",0.734,0.914,-0.58,0.11,0,0,310,"burners"],["b",0.12,0.22,-0.72,0.79,0.58,-0.12,175,"oven_cavity"]]},
  {"id":"NSF-SALAMANDER-36","label":"36\" Salamander Broiler","name":"36 in Overhead Infrared Salamander Broiler","ext":[0.914,0.457,0.445],"card":"surfaces/cookware-ss304","util":"gas:NG:0.5in:10300W;exhaust:TypeI","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.457,0.165,0],["pan_seat_grid","pan_seat",0.457,0.175,-0.2285]],"slots":[["ir_emitter_deck","salamander.ir_temp_k"]],"zones":[["b",0.06,0.08,-0.38,0.85,0.3,-0.06,300,"broiler_grid"]]},
  {"id":"NSF-GRIDDLE-24","label":"24\" Flat-Top Griddle","name":"24 in Heavy 1 in Steel Plate Flat-Top Griddle","ext":[0.61,0.813,0.381],"card":"surfaces/cookware-carbon-steel","util":"gas:NG:0.75in:17600W;exhaust:TypeI","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","left_hand_grasp",0.11,0.135,-0.02],["stance_front","stance_anchor",0.305,0,0.45],["pan_seat_platen","pan_seat",0.305,0.3014,-0.42]],"slots":[["platen_surface","griddle.platen_temp_k"]],"zones":[["b",0.025,0.28,-0.69,0.585,0.305,-0.08,205,"platen"]]},
  {"id":"NSF-GRIDDLE-36","label":"36\" Flat-Top Griddle","name":"36 in Heavy 1 in Steel Plate Flat-Top Griddle","ext":[0.914,0.813,0.381],"card":"surfaces/cookware-carbon-steel","util":"gas:NG:0.75in:26400W;exhaust:TypeI","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","left_hand_grasp",0.11,0.135,-0.02],["stance_front","stance_anchor",0.457,0,0.45],["pan_seat_platen","pan_seat",0.457,0.3014,-0.42]],"slots":[["platen_surface","griddle.platen_temp_k"]],"zones":[["b",0.025,0.28,-0.69,0.889,0.305,-0.08,205,"platen"]]},
  {"id":"NSF-GRIDDLE-48","label":"48\" Flat-Top Griddle","name":"48 in Heavy 1 in Steel Plate Flat-Top Griddle","ext":[1.219,0.813,0.381],"card":"surfaces/cookware-carbon-steel","util":"gas:NG:0.75in:35200W;exhaust:TypeI","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","left_hand_grasp",0.11,0.135,-0.02],["stance_front","stance_anchor",0.6095,0,0.45],["pan_seat_platen","pan_seat",0.6095,0.3014,-0.42]],"slots":[["platen_surface","griddle.platen_temp_k"]],"zones":[["b",0.025,0.28,-0.69,1.194,0.305,-0.08,205,"platen"]]},
  {"id":"NSF-CHARBROILER-36","label":"36\" Radiant Charbroiler","name":"36 in Radiant-Bar Cast-Iron Grate Charbroiler","ext":[0.914,0.813,0.381],"card":"surfaces/cookware-cast-iron-seasoned","util":"gas:NG:0.75in:35200W;exhaust:TypeI","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["stance_front","stance_anchor",0.457,0,0.45],["pan_seat_grates","pan_seat",0.457,0.357,-0.4065]],"slots":[["radiant_deck","charbroiler.radiant_temp_k"]],"zones":[["b",0.03,0.28,-0.67,0.884,0.33,-0.06,340,"grates"]]},
  {"id":"NSF-FRYER-TWIN-50LB","label":"50-lb Twin-Basket Fryer","name":"50 lb Tube-Fired Twin-Basket Commercial Deep Fryer","ext":[0.394,0.762,1.143],"card":"kinetics/frying-oil-degradation-kinetics","util":"gas:NG:0.75in:35200W;exhaust:TypeI","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["grip_L","left_hand_grasp",0.111,0.93,-0.01],["grip_R","right_hand_grasp",0.281,0.93,-0.01],["handle_pull","right_hand_grasp",0.354,0.36,-0.005],["stance_front","stance_anchor",0.197,0,0.45],["pour_oil","pour_target",0.197,0.84,-0.35]],"slots":[["oil_surface","fryer.oil_temp_k"]],"zones":[["b",0.035,0.56,-0.52,0.359,0.89,-0.12,180,"oil_vat"]]},
  {"id":"NSF-INDUCTION-HOB-2Z","label":"2-Zone Induction Hob","name":"Countertop Dual-Zone Commercial Induction Range","ext":[0.4,0.71,0.145],"card":"surfaces/cookware-ss304","util":"elec:208V:3ph:7000W","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["stance_front","stance_anchor",0.2,0,0.45],["pan_seat_coil_0","pan_seat",0.2,0.145,-0.21],["pan_seat_coil_1","pan_seat",0.2,0.145,-0.51]],"slots":[["induction_glass","induction.coil_power_w"]],"zones":[["s",0.2,0.145,-0.2,0.135,0,0,245,"coils"],["s",0.2,0.145,-0.51,0.135,0,0,245,"coils"]]},
  {"id":"NSF-INDUCTION-SUITE-4Z","label":"4-Zone Induction Suite","name":"36 in Four-Zone Commercial Induction Floor Suite","ext":[0.914,0.813,0.914],"card":"surfaces/cookware-ss304","util":"elec:208V:3ph:14000W","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["stance_front","stance_anchor",0.457,0,0.45],["pan_seat_coil_0","pan_seat",0.24,0.914,-0.24],["pan_seat_coil_1","pan_seat",0.674,0.914,-0.24],["pan_seat_coil_2","pan_seat",0.24,0.914,-0.573],["pan_seat_coil_3","pan_seat",0.674,0.914,-0.573]],"slots":[["induction_glass","induction.coil_power_w"]],"zones":[["s",0.24,0.914,-0.24,0.14,0,0,245,"coils"],["s",0.67,0.914,-0.24,0.14,0,0,245,"coils"],["s",0.24,0.914,-0.57,0.14,0,0,245,"coils"],["s",0.67,0.914,-0.57,0.14,0,0,245,"coils"]]},
  {"id":"NSF-COMBI-6PAN","label":"6-Pan Combi Oven","name":"6-Pan Countertop Electric Combi-Steamer Oven","ext":[0.847,0.776,0.782],"card":"surfaces/doorway-gosney-olama-infiltration","util":"elec:208V:3ph:10800W;water:0.75in;drain:50mm;exhaust:TypeI","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.627,0.391,-0.005],["grip_R","right_hand_grasp",0.627,0.391,-0.005],["stance_front","stance_anchor",0.4235,0,0.55],["pan_seat_mid_rail","pan_seat",0.3285,0.391,-0.42]],"slots":[["cavity_glass","combi.cavity_temp_k_rails_6"]],"zones":[["b",0.11,0.12,-0.64,0.68,0.66,-0.12,235,"cavity"]]},
  {"id":"NSF-COMBI-10PAN","label":"10-Pan Combi Oven","name":"10-Pan Full-Size Combi-Steamer Oven on Stand","ext":[0.847,0.776,1.79],"card":"surfaces/doorway-gosney-olama-infiltration","util":"elec:208V:3ph:18900W;water:0.75in;drain:50mm;exhaust:TypeI","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.627,1.265,-0.005],["grip_R","right_hand_grasp",0.627,1.265,-0.005],["stance_front","stance_anchor",0.4235,0,0.55],["pan_seat_mid_rail","pan_seat",0.3285,1.265,-0.42]],"slots":[["cavity_glass","combi.cavity_temp_k_rails_10"]],"zones":[["b",0.11,0.88,-0.64,0.68,1.62,-0.12,235,"cavity"]]},
  {"id":"NSF-COMBI-20PAN-ROLLIN","label":"20-Pan Roll-In Combi","name":"20-Pan Roll-In Combi-Steamer Oven with Rack Dock","ext":[0.879,0.842,1.872],"card":"surfaces/doorway-gosney-olama-infiltration","util":"elec:480V:3ph:37200W;water:0.75in;drain:50mm;exhaust:TypeI","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.659,0.926,-0.005],["grip_R","right_hand_grasp",0.659,0.926,-0.005],["stance_front","stance_anchor",0.4395,0,0.55],["pan_seat_mid_rail","pan_seat",0.3445,0.926,-0.42],["rack_dock_rollin","rack_dock",0.3445,0.04,-0.42]],"slots":[["cavity_glass","combi.cavity_temp_k_rails_12"]],"zones":[["b",0.11,0.14,-0.7,0.72,1.68,-0.12,230,"cavity"]]},
  {"id":"NSF-TILT-SKILLET-40GAL","label":"40-Gal Tilt Braising Pan","name":"40 Gallon Floor-Model Tilting Braising Pan","ext":[1.219,0.991,1.067],"card":"surfaces/cookware-ss304","util":"elec:208V:3ph:18000W;water:0.5in;exhaust:TypeI","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.6095,0.944,-0.08],["grip_R","right_hand_grasp",1.149,0.72,-0.01],["stance_front","stance_anchor",0.6095,0,0.5],["pour_target","pour_target",0.6095,0.854,-0.5]],"slots":[["braising_bed","skillet.bed_temp_k"]],"zones":[["b",0.1,0.68,-0.76,1.11,0.914,-0.14,185,"pan_bed"]]},
  {"id":"NSF-KETTLE-40GAL","label":"40-Gal Steam Kettle","name":"40 Gallon Floor Steam-Jacketed Kettle with Draw-Off Valve","ext":[0.914,0.965,1.067],"card":"surfaces/cookware-ss304","util":"elec:208V:3ph:18000W;water:0.5in;drain:50mm;exhaust:TypeII","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.457,1,-0.08],["grip_R","right_hand_grasp",0.457,0.34,-0.02],["stance_front","stance_anchor",0.457,0,0.45],["pour_target","pour_target",0.457,0.89,-0.4825]],"slots":[["kettle_liquid","kettle.liquid_temp_k"]],"zones":[["s",0.457,0.72,-0.482,0.33,0,0,105,"bowl_center"]]},
  {"id":"NSF-MICROWAVE-COMM","label":"Commercial Microwave","name":"Heavy-Duty Commercial Countertop Microwave Oven","ext":[0.42,0.56,0.34],"card":"surfaces/cookware-ss304","util":"elec:208V:1ph:3200W","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.285,0.17,0],["pan_seat_cavity","pan_seat",0.189,0.05,-0.28]],"slots":[["cavity_window","microwave.power_w"]],"zones":[["b",0.04,0.06,-0.44,0.34,0.26,-0.08,75,"cavity"]]},
  {"id":"NSF-TOASTER-CONVEYOR","label":"Conveyor Toaster","name":"Commercial Countertop Conveyor Toaster","ext":[0.368,0.457,0.387],"card":"surfaces/cookware-ss304","util":"elec:208V:1ph:2700W","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["pan_seat_feed","pan_seat",0.184,0.18,-0.04],["stance_front","stance_anchor",0.184,0,0.45]],"slots":[["quartz_bar","toaster.element_temp_k"]],"zones":[["b",0.05,0.14,-0.39,0.318,0.22,-0.06,260,"belt"]]},
  {"id":"NSF-AIRFRYER-COMM","label":"Commercial Air Fryer","name":"Commercial Countertop High-Velocity Air Fryer","ext":[0.457,0.559,0.432],"card":"surfaces/cookware-ss304","util":"elec:208V:1ph:3400W","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.2285,0.15,-0.01],["grip_R","right_hand_grasp",0.2285,0.15,-0.01]],"slots":[["airfryer_shell","airfryer.cavity_temp_k"]],"zones":[["b",0.07,0.08,-0.42,0.387,0.28,-0.06,200,"basket_cavity"]]},
  {"id":"NSF-SOUSVIDE-BATH","label":"Sous-Vide Bath","name":"28 L Polycarbonate Sous-Vide Bath with Immersion Circulator","ext":[0.33,0.53,0.36],"card":"surfaces/cookware-ss304","util":"elec:120V:1ph:1100W","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["grip_L","left_hand_grasp",0,0.205,-0.265],["grip_R","right_hand_grasp",0.33,0.205,-0.265],["pour_target","pour_target",0.165,0.18,-0.2385]],"slots":[["bath_water","sousvide.water_temp_k"]],"zones":[["b",0.015,0.02,-0.515,0.315,0.19,-0.015,57,"bath_water"]]},
  {"id":"NSF-HEATLAMP-STRIP-48","label":"48\" Strip Heat Lamp","name":"48 in Commercial Sheathed-Element Infrared Strip Warmer","ext":[1.219,0.152,0.0636],"card":"surfaces/cookware-ss304","util":"elec:120V:1ph:1100W","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["pan_seat_pass","pan_seat",0.6095,-0.35,-0.076]],"slots":[["element_emitter","heatlamp.emitter_temp_k"]],"zones":[["b",0.04,0.01,-0.092,1.179,0.03,-0.06,480,"sheathed_element"]]},
  {"id":"NSF-HOLDING-CABINET","label":"Heated Holding Cabinet","name":"Full-Height Insulated Heated Holding and Proofing Cabinet","ext":[0.66,0.813,1.778],"card":"surfaces/doorway-gosney-olama-infiltration","util":"elec:120V:1ph:2000W","cls":"appliance","fam":"hot_line","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.62,0.9715,-0.005],["stance_front","stance_anchor",0.33,0,0.5]],"slots":[["cabinet_shell","holding.cavity_temp_k"]],"zones":[["b",0.07,0.18,-0.73,0.59,1.62,-0.08,77,"cavity"]]},
  {"id":"NSF-WALKIN-COOLER","label":"8×10 ft Walk-In Cooler","name":"8x10 ft Modular Insulated Walk-In Cooler with Strip Curtain","ext":[2.438,3.048,2.438],"card":"surfaces/lumped-cooling-coil-walkin","util":"elec:208V:1ph:2200W;drain:25mm","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["handle_pull","right_hand_grasp",1.584,1.08,-0.005],["grip_R","right_hand_grasp",1.584,1.08,-0.005],["stance_front","stance_anchor",1.217,0,0.6],["rack_dock_interior","rack_dock",1.219,0.102,-1.524]],"slots":[["evaporator_coil","walkin.coil_frost_mm"],["door_infiltration","walkin.infiltration_w"]],"zones":[["b",0.102,0.102,-2.946,2.336,2.336,-0.102,3.3,"interior"]]},
  {"id":"NSF-WALKIN-FREEZER","label":"8×8 ft Walk-In Freezer","name":"8x8 ft Modular Insulated Walk-In Freezer with Strip Curtain","ext":[2.438,2.438,2.438],"card":"surfaces/lumped-cooling-coil-walkin","util":"elec:208V:1ph:3400W;drain:25mm","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["handle_pull","right_hand_grasp",1.584,1.08,-0.005],["grip_R","right_hand_grasp",1.584,1.08,-0.005],["stance_front","stance_anchor",1.217,0,0.6],["rack_dock_interior","rack_dock",1.219,0.102,-1.219]],"slots":[["evaporator_coil","walkin.coil_frost_mm"],["door_infiltration","walkin.infiltration_w"]],"zones":[["b",0.102,0.102,-2.336,2.336,2.336,-0.102,-18,"interior"]]},
  {"id":"NSF-BLAST-CHILLER-15PAN","label":"15-Pan Blast Chiller","name":"15-Pan Commercial Roll-In and Reach-In Blast Chiller","ext":[0.84,0.88,1.85],"card":"surfaces/lumped-cooling-coil-blast-chiller","util":"elec:208V:3ph:5200W;drain:25mm","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["stance_front","stance_anchor",0.42,0,0.5],["handle_pull","right_hand_grasp",0.8,0.911,-0.005],["grip_R","right_hand_grasp",0.8,0.911,-0.005],["rack_dock_blast","rack_dock",0.42,0.212,-0.45]],"slots":[["thermostat_display","cooler.cavity_temp_k"]],"zones":[["b",0.09,0.22,-0.76,0.75,1.58,-0.1,-4,"cavity"]]},
  {"id":"NSF-REACHIN-1D","label":"1-Door Reach-In","name":"27 in Single-Door Stainless Reach-In Refrigerator","ext":[0.686,0.813,2.057],"card":"surfaces/doorway-gosney-olama-infiltration","util":"elec:120V:1ph:580W","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["stance_front","stance_anchor",0.343,0,0.5],["handle_pull","right_hand_grasp",0.646,1.0145,-0.005],["grip_R","right_hand_grasp",0.646,1.0145,-0.005]],"slots":[["thermostat_display","cooler.cavity_temp_k"]],"zones":[["b",0.07,0.22,-0.73,0.616,1.72,-0.08,3,"cavity"]]},
  {"id":"NSF-REACHIN-2D","label":"2-Door Reach-In","name":"54 in Two-Door Stainless Reach-In Refrigerator","ext":[1.372,0.813,2.057],"card":"surfaces/doorway-gosney-olama-infiltration","util":"elec:120V:1ph:680W","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["stance_front","stance_anchor",0.686,0,0.5],["handle_pull","right_hand_grasp",0.726,1.0145,-0.005],["grip_L","left_hand_grasp",0.646,1.0145,-0.005],["grip_R","right_hand_grasp",0.726,1.0145,-0.005]],"slots":[["thermostat_display","cooler.cavity_temp_k"]],"zones":[["b",0.07,0.22,-0.73,1.3,1.72,-0.08,3,"cavity"]]},
  {"id":"NSF-UNDERCOUNTER-2D","label":"48\" Undercounter Cooler","name":"48 in Two-Door Undercounter Worktop Refrigerator","ext":[1.219,0.762,0.914],"card":"surfaces/doorway-gosney-olama-infiltration","util":"elec:120V:1ph:420W","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["stance_front","stance_anchor",0.6095,0,0.5],["handle_pull","right_hand_grasp",0.6495,0.507,-0.005],["grip_L","left_hand_grasp",0.5695,0.507,-0.005],["grip_R","right_hand_grasp",0.6495,0.507,-0.005]],"slots":[["thermostat_display","cooler.cavity_temp_k"]],"zones":[["b",0.06,0.16,-0.68,1.16,0.82,-0.08,3,"cavity"]]},
  {"id":"NSF-MAKELINE-COLD-72","label":"72\" Cold Makeline","name":"72 in Refrigerated Cold Makeline Prep Table with Flip Lid","ext":[1.829,0.864,1.118],"card":"surfaces/doorway-gosney-olama-infiltration","util":"elec:120V:1ph:740W","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.9145,1.014,-0.36],["claw_zone_L","left_hand_grasp",0.823,0.936,-0.2],["stance_front","stance_anchor",0.9145,0,0.45],["pan_seat_rail","pan_seat",0.9145,0.934,-0.55]],"slots":[["cold_rail","makeline.rail_temp_k"]],"zones":[["b",0.06,0.88,-0.8,1.769,1.02,-0.34,2,"cold_rail"]]},
  {"id":"NSF-STEAMTABLE-HOT-60","label":"60\" Hot Steam Table","name":"60 in Four-Well Electric Hot Food Steamtable","ext":[1.524,0.762,0.914],"card":"surfaces/cookware-ss304","util":"elec:208V:1ph:3000W","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["stance_front","stance_anchor",0.762,0,0.45],["pan_seat_well_0","pan_seat",0.2,0.914,-0.48]],"slots":[["steam_wells","steamtable.well_temp_k"]],"zones":[["b",0.08,0.7,-0.7,0.4,0.914,-0.26,85,"wells"],["b",0.43,0.7,-0.7,0.75,0.914,-0.26,85,"wells"],["b",0.78,0.7,-0.7,1.1,0.914,-0.26,85,"wells"],["b",1.13,0.7,-0.7,1.45,0.914,-0.26,85,"wells"]]},
  {"id":"NSF-ICEMACHINE-AIRCOOLED","label":"Cube Ice Machine","name":"30 in Air-Cooled Commercial Cuber Ice Machine on Bin","ext":[0.762,0.864,1.778],"card":"surfaces/cookware-ss304","util":"elec:208V:1ph:1450W;water:0.5in;drain:25mm","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.381,0.82,-0.01],["stance_front","stance_anchor",0.381,0,0.45]],"slots":[["condenser_louvers","icemachine.waste_heat_w"]],"zones":[["b",0.06,0.18,-0.78,0.702,1.05,-0.08,0,"ice_bin"],["b",0.74,1.15,-0.7,0.762,1.68,-0.16,45,"condenser_louver"]]},
  {"id":"NSF-HOOD-TYPE1-CANOPY","label":"Type-I Canopy Hood (8 ft)","name":"96x48 in Type I Exhaust Canopy with MUA Plenum","ext":[2.438,1.219,0.61],"card":"surfaces/cookware-ss304","util":"exhaust:TypeI:2800CFM;fire_suppression:UL300","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[],"slots":[["hood_plume_capture","hood.exhaust_cfm"]],"zones":[["b",0.15,0.1,-1.15,2.288,0.52,-0.85,68,"baffle_bank"]]},
  {"id":"NSF-HOOD-TYPE2-CONDENSATE","label":"Type-II Condensate Hood","name":"60x42 in Type II Stainless Condensate Canopy Hood","ext":[1.524,1.067,0.457],"card":"surfaces/cookware-ss304","util":"exhaust:TypeII:1100CFM;drain:25mm","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[],"slots":[["hood_plume_capture","hood.exhaust_cfm"]],"zones":[["b",0.0,0.0,-1.067,1.524,0.05,-0.0,54,"perimeter_gutter"]]},
  {"id":"NSF-DISHWASHER-DOOR","label":"Door-Type Dishwasher","name":"High-Temp Hood Door-Type Commercial Dishwasher","ext":[0.64,0.64,1.524],"card":"surfaces/cookware-ss304","util":"elec:208V:3ph:12500W;water_hot:0.75in;drain:50mm;exhaust:TypeII","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.32,0.959,-0.01],["grip_L","left_hand_grasp",0.16,0.959,-0.01],["grip_R","right_hand_grasp",0.48,0.959,-0.01],["rack_dock_wash","rack_dock",0.32,0.864,-0.32]],"slots":[["wash_chamber","dishwasher.rinse_temp_k"]],"zones":[["b",0.06,0.864,-0.58,0.58,1.38,-0.06,82,"wash_chamber"]]},
  {"id":"NSF-DISHWASHER-CONVEYOR","label":"Conveyor Dishwasher","name":"44 in Rack-Conveyor High-Temp Commercial Dishwasher","ext":[1.118,0.787,1.727],"card":"surfaces/cookware-ss304","util":"elec:480V:3ph:26000W;water_hot:0.75in;drain:50mm;exhaust:TypeII","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["handle_pull","right_hand_grasp",0.76,1.124,-0.01],["rack_dock_entry","rack_dock",0.04,0.864,-0.3935]],"slots":[["conveyor_tunnel","conveyor.wash_temp_k"]],"zones":[["b",0.0,0.864,-0.66,1.118,1.42,-0.12,77,"tunnel"]]},
  {"id":"NSF-BOOSTER-HEATER-36KW","label":"36 kW Booster Heater","name":"36 kW Commercial Electric Sanitizing Rinse Booster Heater","ext":[0.457,0.508,0.457],"card":"surfaces/cookware-ss304","util":"elec:208V:3ph:36000W;water_hot:0.75in","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[],"slots":[["booster_tank","booster.outlet_temp_k"]],"zones":[["b",0.04,0.08,-0.468,0.417,0.42,-0.04,82,"tank"]]},
  {"id":"NSF-SINK-3COMP","label":"3-Compartment Sink","name":"72 in Three-Compartment Commercial Sink with Pre-Rinse Sprayer","ext":[1.829,0.762,1.727],"card":"surfaces/cookware-ss304","util":"water_hot:0.5in;water_cold:0.5in;drain:50mm","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["grip_R","right_hand_grasp",0.9145,1.407,-0.4],["stance_front","stance_anchor",0.9145,0,0.45]],"slots":[["sink_basins","sink.wash_temp_k"]],"zones":[["b",0.36,0.56,-0.64,0.72,0.914,-0.12,43,"basins"],["b",0.74,0.56,-0.64,1.1,0.914,-0.12,43,"basins"],["b",1.12,0.56,-0.64,1.48,0.914,-0.12,43,"basins"]]},
  {"id":"NSF-PREPTABLE-48","label":"48\" Worktable","name":"48x30 in 16-Gauge SS304 Prep Table with Undershelf","ext":[1.219,0.762,0.914],"card":"surfaces/cookware-ss304","util":"none","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["stance_front","stance_anchor",0.6095,0,0.45],["pan_seat_center","pan_seat",0.6095,0.914,-0.381],["grip_L","left_hand_grasp",0.4266,0.914,-0.2],["grip_R","right_hand_grasp",0.7923,0.914,-0.2]],"slots":[["worktop","preptable.surface_temp_k"]],"zones":[]},
  {"id":"NSF-PREPTABLE-72","label":"72\" Worktable","name":"72x30 in 16-Gauge SS304 Prep Table with Backsplash","ext":[1.829,0.762,1.016],"card":"surfaces/cookware-ss304","util":"none","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["stance_front","stance_anchor",0.9145,0,0.45],["pan_seat_center","pan_seat",0.9145,0.914,-0.381],["grip_L","left_hand_grasp",0.6401,0.914,-0.2],["grip_R","right_hand_grasp",1.1888,0.914,-0.2]],"slots":[["worktop","preptable.surface_temp_k"]],"zones":[]},
  {"id":"NSF-PREPTABLE-96","label":"96\" Worktable","name":"96x30 in 16-Gauge SS304 Prep Table with Undershelf","ext":[2.438,0.762,0.914],"card":"surfaces/cookware-ss304","util":"none","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["stance_front","stance_anchor",1.219,0,0.45],["pan_seat_center","pan_seat",1.219,0.914,-0.381],["grip_L","left_hand_grasp",0.8533,0.914,-0.2],["grip_R","right_hand_grasp",1.5847,0.914,-0.2]],"slots":[["worktop","preptable.surface_temp_k"]],"zones":[]},
  {"id":"NSF-SPEEDRACK-20PAN","label":"20-Pan Speed Rack","name":"20.5x26 in 20-Tier Mobile Aluminum Speed Rack on Swivel Casters","ext":[0.521,0.66,1.753],"card":"surfaces/cookware-aluminum-anodized","util":"none","cls":"appliance","fam":"cold_and_support","src":"catalog","sockets":[["push_bar_L","left_hand_grasp",0.012,1.05,-0.012],["push_bar_R","right_hand_grasp",0.509,1.05,-0.012],["grip_L","left_hand_grasp",0.012,1.05,-0.648],["grip_R","right_hand_grasp",0.509,1.05,-0.648],["stance_push","stance_anchor",0.2605,0,0.45],["rack_dock_base","rack_dock",0.2605,0,-0.33]],"slots":[["rack_frame","speedrack.frame_temp_k"]],"zones":[]},
  {"id":"NSF-SKILLET-8IN","label":"8\" Cast-Iron Skillet","name":"8 in Seasoned Cast-Iron Commercial Skillet","ext":[0.203,0.335,0.048],"card":"surfaces/cookware-cast-iron-seasoned","util":"none","cls":"cookware","fam":"cookware_and_food","src":"catalog","sockets":[["grip_R","right_hand_grasp",0.1015,0.039,-0.05],["grip_L","left_hand_grasp",0.1015,0.039,-0.11],["pan_seat_interior","pan_seat",0.1015,0.008,-0.2335]],"slots":[["cooking_surface","skillet.surface_temp_k"]],"zones":[["s",0.1015,0.004,-0.1015,0.095,0,0,215,"cooking_disk"]]},
  {"id":"NSF-SKILLET-10IN","label":"10\" Tri-Ply Skillet","name":"10 in Tri-Ply Magnetic Stainless Commercial Skillet","ext":[0.254,0.42,0.052],"card":"surfaces/cookware-ss304","util":"none","cls":"cookware","fam":"cookware_and_food","src":"catalog","sockets":[["grip_R","right_hand_grasp",0.127,0.043,-0.05],["grip_L","left_hand_grasp",0.127,0.043,-0.11],["pan_seat_interior","pan_seat",0.127,0.008,-0.293]],"slots":[["cooking_surface","skillet.surface_temp_k"]],"zones":[["s",0.127,0.004,-0.127,0.118,0,0,215,"cooking_disk"]]},
  {"id":"NSF-SKILLET-12IN","label":"12\" Carbon-Steel Skillet","name":"12 in Heavy Carbon-Steel Commercial Skillet","ext":[0.305,0.475,0.056],"card":"surfaces/cookware-carbon-steel","util":"none","cls":"cookware","fam":"cookware_and_food","src":"catalog","sockets":[["grip_R","right_hand_grasp",0.1525,0.047,-0.05],["grip_L","left_hand_grasp",0.1525,0.047,-0.465],["pan_seat_interior","pan_seat",0.1525,0.008,-0.3225]],"slots":[["cooking_surface","skillet.surface_temp_k"]],"zones":[["s",0.1525,0.004,-0.1525,0.142,0,0,215,"cooking_disk"]]},
  {"id":"NSF-SKILLET-14IN","label":"14\" Nonstick Sauté","name":"14 in PTFE Nonstick Anodized Aluminum Saute Skillet","ext":[0.356,0.545,0.062],"card":"surfaces/cookware-ptfe-nonstick","util":"none","cls":"cookware","fam":"cookware_and_food","src":"catalog","sockets":[["grip_R","right_hand_grasp",0.178,0.053,-0.05],["grip_L","left_hand_grasp",0.178,0.053,-0.535],["pan_seat_interior","pan_seat",0.178,0.008,-0.367]],"slots":[["cooking_surface","skillet.surface_temp_k"]],"zones":[["s",0.178,0.004,-0.178,0.165,0,0,215,"cooking_disk"]]},
  {"id":"NSF-STOCKPOT-8QT","label":"8-Qt Stockpot","name":"8 Quart Stainless Stockpot with Lid, Liquid and Foam Levels","ext":[0.24,0.3,0.215],"card":"surfaces/cookware-ss304","util":"none","cls":"cookware","fam":"cookware_and_food","src":"catalog","sockets":[["grip_L","left_hand_grasp",0.12,0.133,-0.01],["grip_R","right_hand_grasp",0.12,0.133,-0.29],["handle_pull","right_hand_grasp",0.12,0.205,-0.15],["pour_target","pour_target",0.12,0.1396,-0.15]],"slots":[["liquid_level","stockpot.liquid_temp_k"],["bikerman_foam","stockpot.foam_height_m"]],"zones":[["s",0.12,0.004,-0.15,0.116,0,0,100,"interior_center"]]},
  {"id":"NSF-STOCKPOT-20QT","label":"20-Qt Stockpot","name":"20 Quart Heavy Stainless Stockpot with Lid, Liquid and Foam Levels","ext":[0.31,0.38,0.295],"card":"surfaces/cookware-ss304","util":"none","cls":"cookware","fam":"cookware_and_food","src":"catalog","sockets":[["grip_L","left_hand_grasp",0.155,0.213,-0.01],["grip_R","right_hand_grasp",0.155,0.213,-0.37],["handle_pull","right_hand_grasp",0.155,0.285,-0.19],["pour_target","pour_target",0.155,0.202,-0.19]],"slots":[["liquid_level","stockpot.liquid_temp_k"],["bikerman_foam","stockpot.foam_height_m"]],"zones":[["s",0.155,0.004,-0.19,0.15,0,0,100,"interior_center"]]},
  {"id":"NSF-STOCKPOT-40QT","label":"40-Qt Stockpot","name":"40 Quart Commercial Stockpot with Lid, Liquid and Foam Levels","ext":[0.38,0.46,0.375],"card":"surfaces/cookware-ss304","util":"none","cls":"cookware","fam":"cookware_and_food","src":"catalog","sockets":[["grip_L","left_hand_grasp",0.19,0.293,-0.01],["grip_R","right_hand_grasp",0.19,0.293,-0.45],["handle_pull","right_hand_grasp",0.19,0.365,-0.23],["pour_target","pour_target",0.19,0.2644,-0.23]],"slots":[["liquid_level","stockpot.liquid_temp_k"],["bikerman_foam","stockpot.foam_height_m"]],"zones":[["s",0.19,0.005,-0.23,0.184,0,0,100,"interior_center"]]},
  {"id":"NSF-HOTEL-PAN-FULL-2IN","label":"Full Hotel Pan, 2\"","name":"Full-Size GN 1/1 2 in 22-Gauge SS304 Hotel Pan","ext":[0.53,0.325,0.065],"card":"surfaces/cookware-ss304","util":"none","cls":"cookware","fam":"cookware_and_food","src":"catalog","sockets":[["grip_L","left_hand_grasp",0.008,0.063,-0.1625],["grip_R","right_hand_grasp",0.522,0.063,-0.1625],["pour_target","pour_target",0.265,0.0325,-0.1625]],"slots":[["pan_wall","hotel_pan.temp_k"]],"zones":[["b",0.015,0.002,-0.31,0.515,0.063,-0.015,74,"cavity"]]},
  {"id":"NSF-HOTEL-PAN-FULL-4IN","label":"Full Hotel Pan, 4\"","name":"Full-Size GN 1/1 4 in 22-Gauge SS304 Hotel Pan","ext":[0.53,0.325,0.1],"card":"surfaces/cookware-ss304","util":"none","cls":"cookware","fam":"cookware_and_food","src":"catalog","sockets":[["grip_L","left_hand_grasp",0.008,0.098,-0.1625],["grip_R","right_hand_grasp",0.522,0.098,-0.1625],["pour_target","pour_target",0.265,0.05,-0.1625]],"slots":[["pan_wall","hotel_pan.temp_k"]],"zones":[["b",0.015,0.002,-0.31,0.515,0.098,-0.015,74,"cavity"]]},
  {"id":"NSF-HOTEL-PAN-FULL-6IN","label":"Full Hotel Pan, 6\"","name":"Full-Size GN 1/1 6 in 22-Gauge SS304 Hotel Pan","ext":[0.53,0.325,0.15],"card":"surfaces/cookware-ss304","util":"none","cls":"cookware","fam":"cookware_and_food","src":"catalog","sockets":[["grip_L","left_hand_grasp",0.008,0.148,-0.1625],["grip_R","right_hand_grasp",0.522,0.148,-0.1625],["pour_target","pour_target",0.265,0.075,-0.1625]],"slots":[["pan_wall","hotel_pan.temp_k"]],"zones":[["b",0.015,0.002,-0.31,0.515,0.148,-0.015,74,"cavity"]]},
  {"id":"NSF-HOTEL-PAN-HALF-4IN","label":"Half Hotel Pan, 4\"","name":"Half-Size GN 1/2 4 in 22-Gauge SS304 Hotel Pan","ext":[0.325,0.265,0.1],"card":"surfaces/cookware-ss304","util":"none","cls":"cookware","fam":"cookware_and_food","src":"catalog","sockets":[["grip_L","left_hand_grasp",0.008,0.098,-0.1325],["grip_R","right_hand_grasp",0.317,0.098,-0.1325],["pour_target","pour_target",0.1625,0.05,-0.1325]],"slots":[["pan_wall","hotel_pan.temp_k"]],"zones":[["b",0.012,0.002,-0.253,0.313,0.098,-0.012,74,"cavity"]]},
  {"id":"NSF-CARRIER-INSULATED","label":"Insulated Food Carrier","name":"Front-Load Insulated Thermal Food-Pan Carrier with Vent Cap","ext":[0.46,0.635,0.625],"card":"surfaces/doorway-gosney-olama-infiltration","util":"none","cls":"cookware","fam":"cookware_and_food","src":"catalog","sockets":[["grip_L","left_hand_grasp",0.01,0.495,-0.33],["grip_R","right_hand_grasp",0.45,0.495,-0.33],["handle_pull","right_hand_grasp",0.42,0.3125,-0.005]],"slots":[["carrier_vent","carrier.cavity_rh"]],"zones":[["b",0.05,0.05,-0.575,0.41,0.565,-0.06,76,"cavity"]]},
  {"id":"NSF-CUTTING-BOARD-HDPE-24X18","label":"HDPE Cutting Board","name":"24x18x0.75 in High-Density Polyethylene Cutting Board","ext":[0.61,0.457,0.019],"card":"surfaces/atkins-blade-food-cutting","util":"none","cls":"utensil","fam":"cookware_and_food","src":"catalog","sockets":[["claw_zone_L","left_hand_grasp",0.2562,0.019,-0.2285],["pan_seat_board","pan_seat",0.305,0.019,-0.2285]],"slots":[],"zones":[]},
  {"id":"NSF-CHEF-KNIFE-8IN","label":"8\" Chef Knife","name":"8 in Forged High-Carbon Stainless Chef Knife","ext":[0.048,0.335,0.022],"card":"surfaces/atkins-blade-food-cutting","util":"none","cls":"utensil","fam":"cookware_and_food","src":"catalog","sockets":[["pinch_grip_R","right_hand_grasp",0.025,0.011,-0.2298],["grip_R","right_hand_grasp",0.025,0.011,-0.2928],["claw_zone_L","left_hand_grasp",0.024,0.011,-0.0762]],"slots":[["blade_edge","knife.cutting_force_n"]],"zones":[]},
  {"id":"NSF-CHEF-KNIFE-10IN","label":"10\" Chef Knife","name":"10 in Forged High-Carbon Stainless Chef Knife","ext":[0.052,0.385,0.024],"card":"surfaces/atkins-blade-food-cutting","util":"none","cls":"utensil","fam":"cookware_and_food","src":"catalog","sockets":[["pinch_grip_R","right_hand_grasp",0.027,0.012,-0.2622],["grip_R","right_hand_grasp",0.027,0.012,-0.3252],["claw_zone_L","left_hand_grasp",0.026,0.012,-0.0876]],"slots":[["blade_edge","knife.cutting_force_n"]],"zones":[]},
  {"id":"NSF-SLICER-ROTARY","label":"Rotary Slicer","name":"12 in Commercial Gravity-Feed Rotary Food Slicer","ext":[0.56,0.66,0.48],"card":"surfaces/atkins-blade-food-cutting","util":"elec:120V:1ph:420W","cls":"appliance","fam":"cookware_and_food","src":"catalog","sockets":[["grip_R","right_hand_grasp",0.48,0.36,-0.1],["handle_pull","right_hand_grasp",0.48,0.36,-0.1],["grip_L","left_hand_grasp",0.12,0.18,-0.02]],"slots":[["slicer_blade","slicer.cutting_power_w"]],"zones":[]},
  {"id":"NSF-GRATER-BOX","label":"Box Grater","name":"9 in Four-Sided Stainless Commercial Box Grater","ext":[0.115,0.095,0.235],"card":"surfaces/atkins-blade-food-cutting","util":"none","cls":"utensil","fam":"cookware_and_food","src":"catalog","sockets":[["grip_L","left_hand_grasp",0.0575,0.227,-0.0475],["grip_R","right_hand_grasp",0.0575,0.1175,0.02]],"slots":[],"zones":[]},
  {"id":"NSF-PROCESSOR-SBLADE","label":"S-Blade Food Processor","name":"Commercial Continuous-Feed S-Blade Food Processor","ext":[0.24,0.31,0.41],"card":"surfaces/atkins-blade-food-cutting","util":"elec:120V:1ph:1100W","cls":"appliance","fam":"cookware_and_food","src":"catalog","sockets":[["grip_R","right_hand_grasp",0.12,0.39,-0.125],["grip_L","left_hand_grasp",0.12,0.1,-0.02]],"slots":[["s_blade","processor.shaft_power_w"]],"zones":[]},
  {"id":"NSF-TONGS-12IN","label":"12\" Tongs","name":"12 in Heavy-Duty Scalloped Stainless Utility Tongs","ext":[0.055,0.305,0.038],"card":"surfaces/cookware-ss304","util":"none","cls":"utensil","fam":"cookware_and_food","src":"catalog","sockets":[["grip_R","right_hand_grasp",0.0275,0.019,-0.1983],["grip_L","left_hand_grasp",0.0275,0.019,-0.1983]],"slots":[],"zones":[]},
  {"id":"NSF-SPATULA-SLOTTED","label":"Slotted Fish Turner","name":"13 in Slotted Stainless Fish and Griddle Turner Spatula","ext":[0.082,0.33,0.028],"card":"surfaces/cookware-ss304","util":"none","cls":"utensil","fam":"cookware_and_food","src":"catalog","sockets":[["grip_R","right_hand_grasp",0.041,0.019,-0.264]],"slots":[],"zones":[]},
  {"id":"NSF-WHISK-12IN","label":"12\" Balloon Whisk","name":"12 in Heavy Commercial French Wire Balloon Whisk","ext":[0.072,0.305,0.072],"card":"surfaces/cookware-ss304","util":"none","cls":"utensil","fam":"cookware_and_food","src":"catalog","sockets":[["grip_R","right_hand_grasp",0.036,0.036,-0.2379]],"slots":[],"zones":[]},
  {"id":"NSF-LADLE-8OZ","label":"8-oz Ladle","name":"8 oz One-Piece Stainless Commercial Portion Ladle","ext":[0.102,0.355,0.072],"card":"surfaces/cookware-ss304","util":"none","cls":"utensil","fam":"cookware_and_food","src":"catalog","sockets":[["grip_R","right_hand_grasp",0.051,0.062,-0.2769],["pour_target","pour_target",0.051,0.048,-0.051]],"slots":[],"zones":[]},
  {"id":"NSF-DEGLAZE-BOTTLE","label":"Deglazing Wine Bottle","name":"750 mL Culinary Deglaze Wine Bottle and Pouring Vessel","ext":[0.076,0.076,0.295],"card":"surfaces/cookware-copper","util":"none","cls":"utensil","fam":"cookware_and_food","src":"catalog","sockets":[["grip_R","right_hand_grasp",0.038,0.14,0],["pour_target","pour_target",0.038,0.295,-0.038]],"slots":[],"zones":[]},
  {"id":"ARCH-FOOD-STEAK-RIBEYE","label":"Ribeye Steak","name":"Boneless Ribeye Steak Portion that Shrinks as it Sears, with its Fond Layer","ext":[0.165,0.115,0.032],"card":"kinetics/maillard-browning","util":"none","cls":"food","fam":"cookware_and_food","src":"catalog","sockets":[],"slots":[["maillard_crust","steak.maillard_index"],["fond_residue","steak.fond_mass_g"]],"zones":[["b",0.0,0.004,-0.115,0.165,0.032,-0.0,62,"crust_top"],["b",0.0,0.0,-0.115,0.165,0.004,-0.0,150,"fond_layer"]]},
  {"id":"ARCH-FOOD-FRY-LOAD","label":"Basket of Fries","name":"3/8 in Cut French Fry Basket Portion with its Oil Sheen","ext":[0.145,0.285,0.078],"card":"kinetics/crust-glass-transition-gordon-taylor","util":"none","cls":"food","fam":"cookware_and_food","src":"catalog","sockets":[],"slots":[["fry_crust_oil","fries.crust_tg_k"]],"zones":[["b",0.0,0.0,-0.285,0.145,0.078,-0.0,175,"bed"]]},
  {"id":"ARCH-FOOD-EGG-FRIED","label":"Fried Egg","name":"Large Grade-A Fried Egg with Albumen and Yolk Zones","ext":[0.11,0.105,0.022],"card":"kinetics/egg-albumen-two-stage","util":"none","cls":"food","fam":"cookware_and_food","src":"catalog","sockets":[],"slots":[["egg_albumen","egg.albumen_gel_frac"],["egg_yolk","egg.yolk_coag_frac"]],"zones":[["s",0.055,0.011,-0.052,0.024,0,0,68,"yolk_center"]]},
  {"id":"ARCH-FOOD-SOURDOUGH-BOULE","label":"Sourdough Boule","name":"Artisan Sourdough Boule that Rises up to 85% in Oven Spring","ext":[0.18,0.18,0.088],"card":"kinetics/yeast-rosso-fermentation","util":"none","cls":"food","fam":"cookware_and_food","src":"catalog","sockets":[],"slots":[["sourdough_crust","boule.crust_browning"]],"zones":[["s",0.09,0.044,-0.09,0.044,0,0,96,"dome_center"]]},
  {"id":"ARCH-FOOD-BRAISED-PAN-PORTION","label":"Braised Pan Portion","name":"Bulk GN 1/1 Hotel-Pan Braised and Roasted Portion Bed","ext":[0.49,0.29,0.072],"card":"kinetics/meat-collagen-gelatin-hydrolysis","util":"none","cls":"food","fam":"cookware_and_food","src":"catalog","sockets":[],"slots":[["braised_portions","braise.core_temp_k"]],"zones":[["b",0.0,0.0,-0.29,0.49,0.072,-0.0,85,"portion_bed"]]},
  {"id":"NSF-WALKIN-COOLER-14X16","label":"14×16 ft Walk-In Cooler","name":"14x16 ft Commissary Walk-In Cooler (4 C) with Twin Evaporators and 42 in Door","ext":[4.267,4.877,2.438],"card":"","util":"elec:208V:1ph:4200W;drain:25mm","cls":"","fam":"","src":"commissary_catering_3200","sockets":[["handle_pull","right_hand_grasp",2.577,1.08,-0.005],["grip_R","right_hand_grasp",2.577,1.08,-0.005],["stance_front","stance_anchor",2.1335,0,0.6],["rack_dock_interior","rack_dock",2.1335,0.102,-2.4385]],"slots":[["evaporator_coil","walkin-cooler-14x16.evaporator_coil"],["door_infiltration","walkin-cooler-14x16.door_infiltration"]],"zones":[["b",0.102,0.102,-4.775,4.165,2.336,-0.102,3.3,"interior"]]},
  {"id":"NSF-WALKIN-FREEZER-10X14","label":"10×14 ft Walk-In Freezer","name":"10x14 ft Commissary Walk-In Freezer (-18 C) with Heated 42 in Speed-Rack Door","ext":[3.048,4.267,2.438],"card":"","util":"elec:208V:1ph:3800W;drain:25mm","cls":"","fam":"","src":"commissary_catering_3200","sockets":[["handle_pull","right_hand_grasp",1.968,1.08,-0.005],["grip_R","right_hand_grasp",1.968,1.08,-0.005],["stance_front","stance_anchor",1.5245,0,0.6],["rack_dock_interior","rack_dock",1.524,0.102,-2.1335]],"slots":[["evaporator_coil","walkin-freezer-10x14.evaporator_coil"],["door_infiltration","walkin-freezer-10x14.door_infiltration"]],"zones":[["b",0.102,0.102,-4.165,2.946,2.336,-0.102,-18,"interior"]]},
  {"id":"NSF-WIRE-SHELVING-48IN","label":"48\" Wire Shelving","name":"48x24 in 4-Tier Epoxy Wire Dry-Storage Shelving","ext":[1.219,0.61,1.829],"card":"","util":"none","cls":"","fam":"","src":"commissary_catering_3200","sockets":[["stance_front","stance_anchor",0.6095,0,0.45],["pan_seat_mid","pan_seat",0.6095,0.728,-0.305]],"slots":[],"zones":[]},
  {"id":"NSF-VEG-PREP-SINK-48IN","label":"48\" Veg Prep Sink","name":"48 in Single-Bowl Commercial Vegetable Prep Sink with Drainboard","ext":[1.219,0.762,1.118],"card":"","util":"water_hot:0.5in;water_cold:0.5in;drain:50mm","cls":"","fam":"","src":"commissary_catering_3200","sockets":[["grip_R","right_hand_grasp",0.32,1.054,-0.34],["stance_front","stance_anchor",0.6095,0,0.45]],"slots":[],"zones":[]},
  {"id":"NSF-HOOD-TYPE1-CANOPY-144IN","label":"Type-I Canopy Hood (12 ft)","name":"144x48 in 12 ft Commercial Type I Exhaust Canopy with MUA Plenum","ext":[3.658,1.219,0.61],"card":"","util":"exhaust:TypeI:3600CFM;fire_suppression:UL300","cls":"","fam":"","src":"commissary_catering_3200","sockets":[],"slots":[["hood_plume_capture","hood-kettle-12ft.hood_plume_capture"]],"zones":[["b",0.15,0.1,-1.15,3.508,0.52,-0.85,68,"baffle_bank"]]},
  {"id":"NSF-HOOD-TYPE1-CANOPY-288IN","label":"Type-I Canopy Hood (24 ft)","name":"288x48 in 24 ft Commissary Type I Exhaust Canopy with MUA Plenum","ext":[7.315,1.219,0.61],"card":"","util":"exhaust:TypeI:9800CFM;fire_suppression:UL300","cls":"","fam":"","src":"commissary_catering_3200","sockets":[],"slots":[["hood_plume_capture","hood-hotline-24ft.hood_plume_capture"]],"zones":[["b",0.15,0.1,-1.15,7.165,0.52,-0.85,68,"baffle_bank"]]},
  {"id":"NSF-PORTION-SCALE-DIGI","label":"Digital Portion Scale","name":"Commercial Stainless Digital Portion and Staging Scale (50 kg)","ext":[0.33,0.356,0.114],"card":"","util":"elec:120V:1ph:20W","cls":"","fam":"","src":"commissary_catering_3200","sockets":[["pan_seat_platter","pan_seat",0.165,0.114,-0.178]],"slots":[],"zones":[]},
  {"id":"NSF-HANDWASH-SINK-WALL-17IN","label":"Wall Hand Sink","name":"17 in Wall-Mounted SS304 Handwashing Sink with Splash Guards","ext":[0.432,0.381,0.457],"card":"","util":"water_hot:0.5in;water_cold:0.5in;drain:38mm","cls":"","fam":"","src":"commissary_catering_3200","sockets":[["stance_front","stance_anchor",0.216,-0.8,0.45],["grip_R","right_hand_grasp",0.216,0.26,-0.326]],"slots":[["handwash_basin","handwash-sink-prep.handwash_basin"]],"zones":[["b",0.06,0.14,-0.32,0.372,0.3,-0.06,38,"basin"]]},
  {"id":"NSF-HOOD-TYPE1-ISLAND-144IN","label":"Type-I Island Hood (12 ft)","name":"144x84 in 12 ft Double-Island Type I Exhaust Canopy with MUA Plenum","ext":[3.658,2.134,0.61],"card":"","util":"exhaust:TypeI:5200CFM;fire_suppression:UL300","cls":"","fam":"","src":"fine_dining_950","sockets":[],"slots":[["hood_plume_capture","hood-island-1.hood_plume_capture"]],"zones":[["b",0.15,0.1,-1.38,3.508,0.52,-0.75,68,"baffle_bank"]]},
  {"id":"NSF-HOOD-TYPE1-CANOPY-168IN","label":"Type-I Canopy Hood (14 ft)","name":"168x48 in 14 ft Commercial Type I Exhaust Canopy with MUA Plenum","ext":[4.267,1.219,0.61],"card":"","util":"exhaust:TypeI:4200CFM;fire_suppression:UL300","cls":"","fam":"","src":"freestanding_qsr_1450","sockets":[],"slots":[["hood_plume_capture","hood-type1-west.hood_plume_capture"]],"zones":[["b",0.15,0.1,-1.15,4.117,0.52,-0.85,68,"baffle_bank"]]},
  {"id":"NSF-WALKIN-FREEZER-8X10","label":"8×10 ft Walk-In Freezer","name":"8x10 ft Modular Insulated Walk-In Freezer (-18 C) with Strip Curtain","ext":[2.438,3.048,2.438],"card":"","util":"elec:208V:1ph:3400W;drain:25mm","cls":"","fam":"","src":"freestanding_qsr_1450","sockets":[["handle_pull","right_hand_grasp",1.584,1.08,-0.005],["grip_R","right_hand_grasp",1.584,1.08,-0.005],["stance_front","stance_anchor",1.217,0,0.6],["rack_dock_interior","rack_dock",1.219,0.102,-1.524]],"slots":[["evaporator_coil","walkin-freezer-1.evaporator_coil"],["door_infiltration","walkin-freezer-1.door_infiltration"]],"zones":[["b",0.102,0.102,-2.946,2.336,2.336,-0.102,-18,"interior"]]},
  {"id":"NSF-HOTPLATE-ELEC-2B","label":"2-Element Electric Hotplate","name":"Countertop Two-Element Open-Coil Electric Hotplate","ext":[0.305,0.61,0.125],"card":"","util":"elec","cls":"appliance","fam":"hot_line","src":"stove","sockets":[],"slots":[],"zones":[["s",0.1525,0.125,-0.19,0.1,0,0,500,"element"],["s",0.1525,0.125,-0.45,0.1,0,0,500,"element"]]},
];

const BY_ID = new Map(ARCH.map((a) => [a.id, a]));
/** The piece the turntable opens on. */
const FIRST_ID = "NSF-RANGE-GAS-6B";
/** The piece mountEquipmentTurntable builds as it mounts, for main.ts to build ahead (as show() builds it). */
export const EQUIPMENT_PIECES: readonly Piece[] = (() => {
  const ext = BY_ID.get(FIRST_ID)?.ext;
  return ext ? [[FIRST_ID, { w: ext[0], d: ext[1], h: ext[2] }, 2, true]] : [];
})();
const FLAGSHIP_IDS = new Set(FLAGSHIP_EQUIPMENT.map((f) => f.id));

function tabOf(a: Arch): EquipmentTabId {
  if (/SHELVING/.test(a.id)) return "cold_chain";
  if (/PORTION-SCALE/.test(a.id)) return "cookware_food";
  const entry = { archetype_id: a.id, asset_class: a.cls, family: a.fam } as unknown as CatalogEntry;
  return categorizeCatalogEntry(entry);
}

/** "gas:NG:0.75in:61500W;exhaust:TypeI" → "Gas 61.5 kW · Type I hood". */
function utilities(u: string): string {
  if (u === "none") return "No utilities: a passive part";
  if (!u) return "";
  const out: string[] = [];
  const kw = (s: string) => {
    const m = /([\d.]+)W$/.exec(s);
    if (!m) return "";
    const w = Number(m[1]);
    return w >= 1000 ? `${(w / 1000).toFixed(w % 1000 === 0 ? 0 : 1)} kW` : `${w} W`;
  };
  for (const part of u.split(";")) {
    const f = part.split(":");
    const head = f[0] ?? "";
    if (head === "gas") out.push(`Gas ${kw(part)}`);
    else if (head === "elec") out.push(f[1] ? `Electric ${kw(part)} · ${f[1]} ${f[2] === "3ph" ? "3-phase" : "1-phase"}` : "Electric");
    else if (head === "exhaust") out.push(f[2] ? `Type ${f[1]?.replace("Type", "")} exhaust ${Number(f[2].replace("CFM", "")).toLocaleString()} CFM` : `Type ${f[1]?.replace("Type", "")} hood`);
    else if (head === "fire_suppression") out.push(`${f[1]?.replace("UL", "UL ")} suppression`);
    else if (head === "water" || head === "water_hot" || head === "water_cold") out.push(head === "water" ? "Water" : head === "water_hot" ? "Hot water" : "Cold water");
    else if (head === "drain") out.push(`${(f[1] ?? "").replace(/(\d)mm$/, "$1 mm")} drain`);
  }
  return [...new Set(out)].join(" · ");
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Thermal ramp: cold blues through room-temperature indigo to white heat                                          */
/* ------------------------------------------------------------------------------------------------------------ */

const AMBIENT_C = 22;
const RAMP: [number, number][] = [
  [-25, 0xe4f8ff],
  [-12, 0x6fd6ff],
  [0, 0x2f80ff],
  [10, 0x2b3fae],
  [22, 0x1b1736],
  [45, 0x5a1f8c],
  [80, 0xa51a96],
  [130, 0xe2264b],
  [190, 0xf3661b],
  [260, 0xf8b31d],
  [340, 0xfff0a6],
  [460, 0xffffff],
];
const N_ZONES = 10;

const FIELD_GLSL = /* glsl */ `
  uniform vec4 uZA[${N_ZONES}];
  uniform vec4 uZB[${N_ZONES}];
  uniform float uZT[${N_ZONES}];
  uniform int uZN;
  uniform float uAmb;
  uniform float uLambda;
  uniform vec3 uRampC[${RAMP.length}];
  uniform float uRampT[${RAMP.length}];
  float zoneDist(vec4 a, vec4 b, vec3 p) {
    if (a.w < 0.5) {
      vec3 q = max(max(a.xyz - p, p - b.xyz), 0.0);
      return length(q);
    }
    return max(0.0, length(p - a.xyz) - b.x);
  }
  float temperature(vec3 p) {
    float best = 0.0;
    for (int i = 0; i < ${N_ZONES}; i++) {
      if (i >= uZN) break;
      float d = zoneDist(uZA[i], uZB[i], p);
      float dT = (uZT[i] - uAmb) * exp(-d / uLambda);
      if (abs(dT) > abs(best)) best = dT;
    }
    return uAmb + best;
  }
  vec3 ramp(float t) {
    vec3 c = uRampC[0];
    for (int i = 1; i < ${RAMP.length}; i++) {
      float k = clamp((t - uRampT[i - 1]) / max(1e-3, uRampT[i] - uRampT[i - 1]), 0.0, 1.0);
      c = mix(c, uRampC[i], k);
    }
    return c;
  }
`;

const INSPECT_VERT = /* glsl */ `
  uniform mat4 uRootInv;
  varying vec3 vLocal;
  varying vec3 vNormalV;
  varying vec3 vViewPos;
  void main() {
    vec4 wp = modelMatrix * vec4(position, 1.0);
    vLocal = (uRootInv * wp).xyz;
    vec4 mv = viewMatrix * wp;
    vViewPos = mv.xyz;
    vNormalV = normalize(normalMatrix * normal);
    gl_Position = projectionMatrix * mv;
  }
`;

/** Thermal: a solid false-colour surface lit just enough to read its form. X-ray: an additive fresnel shell. */
const INSPECT_FRAG = /* glsl */ `
  ${FIELD_GLSL}
  uniform float uGhost;
  uniform float uThermal;
  uniform vec3 uTint;
  varying vec3 vLocal;
  varying vec3 vNormalV;
  varying vec3 vViewPos;
  void main() {
    vec3 n = normalize(vNormalV);
    if (!gl_FrontFacing) n = -n;
    vec3 v = normalize(-vViewPos);
    float ndv = clamp(dot(n, v), 0.0, 1.0);
    vec3 c = uThermal > 0.5 ? ramp(temperature(vLocal)) : uTint;
    if (uGhost < 0.5) {
      vec3 L = normalize(vec3(0.35, 0.8, 0.45));
      float diff = 0.5 + 0.5 * clamp(dot(n, L), 0.0, 1.0);
      float rim = pow(1.0 - ndv, 3.0) * 0.45;
      gl_FragColor = vec4(c * (diff + rim) + vec3(0.012), 1.0);
    } else {
      // Additive layers sum in display space, so the shell is written there directly rather than encoded: a
      // faint face-on veil, brighter where a surface turns away at its silhouette.
      float f = pow(1.0 - ndv, 2.5);
      gl_FragColor = vec4(pow(c, vec3(1.0 / 2.2)) * (0.014 + 0.24 * f), 1.0);
      return;
    }
    #include <colorspace_fragment>
  }
`;

function rampUniforms(): { c: THREE.Color[]; t: number[] } {
  return { c: RAMP.map(([, hex]) => new THREE.Color(hex)), t: RAMP.map(([t]) => t) };
}

function rampCss(): string {
  return `linear-gradient(90deg, ${RAMP.map(([, hex], i) => `#${hex.toString(16).padStart(6, "0")} ${((i / (RAMP.length - 1)) * 100).toFixed(1)}%`).join(", ")})`;
}
/** Where a temperature falls on the legend bar (stops are evenly spaced), 0..1. */
function rampPos(tC: number): number {
  if (tC <= RAMP[0]![0]) return 0;
  for (let i = 1; i < RAMP.length; i++) {
    const [t0] = RAMP[i - 1]!;
    const [t1] = RAMP[i]!;
    if (tC <= t1) return (i - 1 + (tC - t0) / (t1 - t0)) / (RAMP.length - 1);
  }
  return 1;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Studio                                                                                                         */
/* ------------------------------------------------------------------------------------------------------------ */

/**
 * A dark studio for the stainless to reflect: a big neutral softbox overhead, a warm strip to the right, a cool one
 * to the left, a warm band behind and a soft card in front, the rest near black, so that steel reads as steel.
 */
function showroomEnvironment(r: THREE.WebGLRenderer): THREE.Texture {
  const env = new THREE.Scene();
  env.background = new THREE.Color(0x0b0908);
  const geos: THREE.BufferGeometry[] = [];
  const mats: THREE.Material[] = [];
  const panel = (w: number, h: number, color: number, k: number, at: [number, number, number]) => {
    const g = new THREE.PlaneGeometry(w, h);
    const m = new THREE.MeshBasicMaterial({ color: new THREE.Color(color).multiplyScalar(k), side: THREE.DoubleSide });
    geos.push(g);
    mats.push(m);
    const mesh = new THREE.Mesh(g, m);
    mesh.position.set(at[0], at[1], at[2]);
    mesh.lookAt(0, 0, 0);
    env.add(mesh);
  };
  panel(8, 5, 0xfff3e6, 2.3, [0.6, 7, 1.8]);
  panel(1.8, 8, 0xffcf9c, 2.4, [7.5, 2, 3]);
  panel(1.4, 7, 0xa9c6ff, 1.3, [-7.5, 2.4, 0.5]);
  panel(11, 1.3, 0xffb072, 1.0, [0, 1.4, -8.5]);
  panel(5, 2.6, 0xfff6ec, 0.55, [-1, 1.6, 9]);
  panel(18, 18, 0x231a13, 1, [0, -7, 0]);
  const pm = new THREE.PMREMGenerator(glOf(r));
  const t = pm.fromScene(env, 0.025).texture;
  pm.dispose();
  for (const g of geos) g.dispose();
  for (const m of mats) m.dispose();
  return t;
}

function canvasTexture(size: number, draw: (x: CanvasRenderingContext2D, s: number) => void, srgb = true): THREE.CanvasTexture {
  const c = document.createElement("canvas");
  c.width = c.height = size;
  const x = c.getContext("2d")!;
  draw(x, size);
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = srgb ? THREE.SRGBColorSpace : THREE.NoColorSpace;
  t.anisotropy = 8;
  return t;
}

/** The platter's top: dark satin with fine turned rings, a brass-coloured index ring and degree ticks at its rim. */
function platterTopTexture(): THREE.CanvasTexture {
  return canvasTexture(2048, (x, S) => {
    const c = S / 2;
    const g = x.createRadialGradient(c, c, 0, c, c, c);
    g.addColorStop(0, "#1d1915");
    g.addColorStop(0.75, "#16130f");
    g.addColorStop(1, "#110e0b");
    x.fillStyle = g;
    x.fillRect(0, 0, S, S);
    // Turned rings, as a lathe leaves them.
    for (let r = 6; r < c * 0.93; r += 2.2) {
      x.strokeStyle = `rgba(255,236,210,${(0.012 + 0.018 * Math.abs(Math.sin(r * 0.37))).toFixed(3)})`;
      x.lineWidth = 1;
      x.beginPath();
      x.arc(c, c, r, 0, Math.PI * 2);
      x.stroke();
    }
    // Index ring and ticks every 5°, longer every 30°.
    x.strokeStyle = "rgba(214,160,86,0.85)";
    x.lineWidth = 3;
    x.beginPath();
    x.arc(c, c, c * 0.95, 0, Math.PI * 2);
    x.stroke();
    for (let i = 0; i < 72; i++) {
      const a = (i / 72) * Math.PI * 2;
      const long = i % 6 === 0;
      const r0 = c * (long ? 0.905 : 0.925);
      x.strokeStyle = long ? "rgba(240,190,110,0.9)" : "rgba(214,160,86,0.55)";
      x.lineWidth = long ? 3 : 2;
      x.beginPath();
      x.moveTo(c + Math.cos(a) * r0, c + Math.sin(a) * r0);
      x.lineTo(c + Math.cos(a) * c * 0.948, c + Math.sin(a) * c * 0.948);
      x.stroke();
    }
  });
}

/** A soft dark pad, a rounded rectangle blurred, for contact shading under a model's footprint. */
function contactTexture(): THREE.CanvasTexture {
  return canvasTexture(256, (x, S) => {
    x.clearRect(0, 0, S, S);
    x.filter = "blur(18px)";
    x.fillStyle = "rgba(0,0,0,0.85)";
    const m = S * 0.2;
    x.beginPath();
    x.roundRect(m, m, S - 2 * m, S - 2 * m, S * 0.08);
    x.fill();
  });
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Disposal: per-build resources go, shared ones (the palette, family caches) stay                                */
/* ------------------------------------------------------------------------------------------------------------ */

const disposedOnce = new WeakSet<object>();
const sharedRes = new WeakSet<object>();

function texturesOf(m: THREE.Material): THREE.Texture[] {
  const out: THREE.Texture[] = [];
  for (const v of Object.values(m)) if (v && (v as THREE.Texture).isTexture) out.push(v as THREE.Texture);
  const u = (m as THREE.ShaderMaterial).uniforms;
  if (u) for (const k of Object.keys(u)) {
    const v = u[k]?.value as THREE.Texture | undefined;
    if (v && v.isTexture) out.push(v);
  }
  return out;
}

function resourcesOf(root: THREE.Object3D): { geos: Set<THREE.BufferGeometry>; mats: Set<THREE.Material>; texs: Set<THREE.Texture>; inst: THREE.InstancedMesh[] } {
  const geos = new Set<THREE.BufferGeometry>();
  const mats = new Set<THREE.Material>();
  const texs = new Set<THREE.Texture>();
  const inst: THREE.InstancedMesh[] = [];
  root.traverse((o) => {
    const any = o as THREE.Mesh;
    if ((o as THREE.InstancedMesh).isInstancedMesh) inst.push(o as THREE.InstancedMesh);
    if (any.geometry && !(o as THREE.Sprite).isSprite) geos.add(any.geometry);
    const list = any.material ? (Array.isArray(any.material) ? any.material : [any.material]) : [];
    for (const m of list) {
      mats.add(m);
      for (const t of texturesOf(m)) texs.add(t);
    }
  });
  return { geos, mats, texs, inst };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Choreography: what each model does while it is shown                                                           */
/* ------------------------------------------------------------------------------------------------------------ */

const smooth = (a: number, b: number, x: number) => {
  const k = Math.max(0, Math.min(1, (x - a) / (b - a)));
  return k * k * (3 - 2 * k);
};
/** A door's cycle: shut, swing open, hold, swing shut, hold. Period `p` seconds, phase offset `o` (of a period). */
const cycle = (t: number, p = 9, o = 0) => {
  const u = (((t / p + o) % 1) + 1) % 1;
  return smooth(0.06, 0.3, u) * (1 - smooth(0.62, 0.86, u));
};
/** A smooth triangle 0..1..0 over `p` seconds. */
const tri = (t: number, p = 10, o = 0) => 0.5 - 0.5 * Math.cos(((t / p + o) % 1) * Math.PI * 2);

interface Rig {
  readonly t: number;
  readonly on: boolean;
  readonly model: Model;
  set(name: string, v: number, eps?: number): void;
  open(v: number): void;
  flows(f: (burner: number) => number): void;
}

type Act = (r: Rig) => void;

/** Gas valves: a slow wave across the burners (their knobs turn with them), or all steady. */
const gasWave = (r: Rig, lo = 0.4) => r.flows((i) => (r.on ? lo + (1 - lo) * tri(r.t, 7, -i * 0.13) : 0.82));
const opens = (p = 9, o = 0): Act => (r) => r.open(r.on ? cycle(r.t, p, o) : 0);

/** Pieces that hang overhead and work on their underside. */
const OVERHEAD = /^NSF-(HOOD|HEATLAMP)/;

const CHOREO: [RegExp, Act][] = [
  [/^NSF-RANGE-GAS/, (r) => (gasWave(r), opens(11, 0.55)(r))],
  [/^NSF-GRIDDLE/, (r) => (gasWave(r, 0.5), opens(10, 0.4)(r))],
  [/^NSF-CHARBROILER/, (r) => (gasWave(r, 0.5), r.set("heat", 0.92), opens(10, 0.4)(r))],
  [
    /^NSF-FRYER/,
    (r) => {
      gasWave(r, 0.6);
      r.set("lift_l", r.on ? cycle(r.t, 8, 0) : 0);
      r.set("lift_r", r.on ? cycle(r.t, 8, 0.5) : 0);
      r.open(r.on ? cycle(r.t, 16, 0.75) * 0.9 : 0);
    },
  ],
  [
    /^NSF-SALAMANDER/,
    (r) => {
      r.set("heat", 0.95);
      r.open(r.on ? cycle(r.t, 9, 0.3) : 0);
      r.set("rack", r.on ? Math.floor((r.t / 2.25) % 4) / 3 : 0.34, 0.01);
    },
  ],
  [/^NSF-INDUCTION-HOB/, (r) => [0, 1].forEach((i) => r.set(`zone${i}`, r.on ? 0.15 + 0.85 * tri(r.t, 6, i * 0.5) : 0.7, 0.01))],
  [/^NSF-INDUCTION-SUITE/, (r) => ([0, 1, 2, 3].forEach((i) => r.set(`zone${i}`, r.on ? 0.15 + 0.85 * tri(r.t, 7, i * 0.25) : 0.7, 0.01)), opens(10, 0.6)(r))],
  [
    /^NSF-COMBI/,
    (r) => {
      const u = (((r.t / 12) % 1) + 1) % 1;
      const door = r.on ? smooth(0.05, 0.25, u) * (1 - smooth(0.7, 0.9, u)) : 0;
      r.open(door);
      r.set("fan", r.on ? 1 - door : 0.6);
      r.set("heat", 0.72);
      r.set("steam", r.on ? 0.1 + 0.12 * tri(r.t, 5) * (1 - door) : 0.1, 0.01);
      r.set("dial", r.on ? tri(r.t, 12) : 0.45, 0.005);
      // Only the floor model has a roll-in trolley; it rolls out while the door stands open.
      if (/ROLLIN/.test(r.model.id)) r.set("trolley", r.on ? smooth(0.3, 0.45, u) * (1 - smooth(0.55, 0.68, u)) : 0);
    },
  ],
  [
    /^NSF-TILT-SKILLET/,
    (r) => {
      const u = (((r.t / 12) % 1) + 1) % 1;
      r.set("lid", r.on ? smooth(0.04, 0.18, u) * (1 - smooth(0.82, 0.95, u)) : 0);
      r.set("tilt", r.on ? smooth(0.24, 0.42, u) * (1 - smooth(0.6, 0.76, u)) : 0);
      r.set("heat", 0.8);
    },
  ],
  [/^NSF-KETTLE/, (r) => (r.set("level", 0.72), r.set("heat", 0.8), opens(9)(r))],
  [
    /^NSF-MICROWAVE/,
    (r) => {
      const d = r.on ? cycle(r.t, 8) : 0;
      r.open(d);
      r.set("power", d < 0.02 ? 1 : 0, 0.5);
    },
  ],
  // The belt runs; the crumb tray (the one part that opens) draws out now and then.
  [/^NSF-TOASTER/, (r) => (r.set("heat", 0.9), r.set("speed", r.on ? 0.65 : 0.3), opens(12, 0.5)(r))],
  [
    /^NSF-AIRFRYER/,
    (r) => {
      const d = r.on ? cycle(r.t, 9) : 0;
      r.open(d);
      r.set("fan", d < 0.02 ? 1 : 0);
      r.set("heat", 0.85);
    },
  ],
  [/^NSF-SOUSVIDE/, (r) => (r.set("heat", 1), r.set("flow", r.on ? 0.85 : 0.3))],
  [/^NSF-HEATLAMP/, (r) => r.set("heat", 0.9)],
  [/^NSF-HOLDING-CABINET/, (r) => (opens(10)(r), r.set("heat", 0.7), r.set("humidity", r.on ? 0.3 + 0.5 * tri(r.t, 8) : 0.5, 0.01))],
  [/^NSF-STEAMTABLE/, (r) => (r.set("steam", r.on ? 0.85 : 0.5), r.set("heat", r.on ? 0.55 + 0.4 * tri(r.t, 8) : 0.8, 0.01))],
  [/^NSF-BOOSTER/, (r) => (r.set("heat", 1), r.set("temp", r.on ? 0.62 + 0.12 * tri(r.t, 7) : 0.68, 0.003), r.set("pressure", r.on ? 0.42 + 0.06 * tri(r.t, 5) : 0.45, 0.003))],
  [/^NSF-WALKIN/, (r) => opens(12, 0.5)(r)],
  [/^NSF-ICEMACHINE/, (r) => (opens(9)(r), r.set("ice", r.on ? 0.45 + 0.4 * tri(r.t, 16) : 0.7, 0.02))],
  [/^NSF-(REACHIN|UNDERCOUNTER|BLAST|MAKELINE)/, (r) => opens(9)(r)],
  [/^NSF-SPEEDRACK/, (r) => (r.set("pans", r.on ? Math.round(20 * tri(r.t, 14)) : 14, 0.5), r.set("swivel", r.on ? r.t * 0.7 : 0, 0.01))],
  [/^NSF-WIRE-SHELVING/, (r) => r.set("stock", 1)],
  [/^NSF-HOOD/, (r) => r.set("lights", 1)],
  [/^NSF-DISHWASHER-DOOR/, (r) => (r.set("wash_temp", 66), r.set("rinse_temp", 82), opens(10)(r))],
  [/^NSF-DISHWASHER-CONVEYOR/, (r) => opens(10)(r)],
  [
    /^NSF-(SINK-3COMP|VEG-PREP-SINK)/,
    (r) => {
      const dr = r.on ? cycle(r.t, 10, 0.5) : 0;
      r.set("drain", dr);
      r.set("water", 0.78 * (1 - dr), 0.01);
    },
  ],
  [/^NSF-HANDWASH/, (r) => r.set("flow", r.on && cycle(r.t, 7) > 0.5 ? 1 : 0, 0.5)],
  [/^NSF-CARRIER/, (r) => (r.set("pans", 3, 0.5), opens(11)(r), r.set("vent", r.on ? cycle(r.t, 11, 0.5) : 0))],
  [/^NSF-PROCESSOR/, (r) => opens(7)(r)],
  [/^NSF-SLICER/, (r) => r.set("thickness", r.on ? 0.25 + 0.5 * tri(r.t, 9) : 0.4, 0.01)],
  [/^NSF-PORTION-SCALE/, (r) => r.set("weight", r.on ? Math.round(2350 * cycle(r.t, 8)) / 1000 : 0, 0.0005)],
  [/^NSF-TONGS/, (r) => r.open(r.on ? cycle(r.t, 3.2) : 0)],
  [/^NSF-STOCKPOT/, (r) => (r.set("level", 0.72), r.set("boil", r.on ? 0.85 : 0.2), opens(9)(r))],
  [/^NSF-HOTEL-PAN/, (r) => (r.set("fill", 0.62), r.set("food", r.on ? Math.floor(r.t / 4) % 4 : 0, 0.5))],
  [
    /^ARCH-FOOD-STEAK/,
    (r) => {
      const b = r.on ? tri(r.t, 12) : 0.72;
      r.set("browning", b, 0.02);
      r.set("shrink", 0.22 * b, 0.004);
    },
  ],
  [/^ARCH-FOOD-SOURDOUGH/, (r) => (r.set("spring", r.on ? 1 + 0.85 * tri(r.t, 12) : 1.85, 0.01), r.set("crust", r.on ? smooth(0.2, 1, tri(r.t, 12)) : 0.8, 0.02))],
  [/^ARCH-FOOD-EGG/, (r) => r.set("set", r.on ? 0.15 + 0.85 * tri(r.t, 10) : 0.9, 0.02)],
  [/^ARCH-FOOD-FRY/, (r) => r.set("browning", r.on ? 0.15 + 0.8 * tri(r.t, 10) : 0.75, 0.02)],
  [/^NSF-DEGLAZE/, (r) => r.set("level", r.on ? 0.15 + 0.7 * tri(r.t, 10) : 0.7, 0.01)],
  [
    /^NSF-HOTPLATE/,
    (r) =>
      [0, 1].forEach((i) => {
        const k = r.on ? 0.25 + 0.75 * tri(r.t, 9, i * 0.5) : 0.85;
        r.set(`knob${i}`, k, 0.01);
        r.set(`heat${i}`, k, 0.01);
      }),
  ],
];

/** The live names in CHOREO that carry a part out past the body (the rest turn knobs, tint, fill or reshape). */
const REACHING = new Set(["trolley", "tilt", "lid", "lift_l", "lift_r"]);

/**
 * How far the model reaches from the turntable's axis (cx, cz) over its choreography: the act is run through 24 s of
 * its cycles (its doors, drawers and lids, and the REACHING names; nothing else, since a food's browning or a bottle's
 * level rebuilds geometry) and the meshes measured at each pose; `into` grows to hold them all. The pose is left as the
 * last sample's; the first operate() after puts it back.
 */
function sweptReach(model: Model, meshes: THREE.Mesh[], act: Act | null, cx: number, cz: number, into: THREE.Box3): number {
  const h = model.hooks;
  if (!h.open && !h.set) return 0;
  const tmp = new THREE.Box3();
  let t = 0;
  const rig: Rig = {
    get t() {
      return t;
    },
    on: true,
    model,
    set: (n, v) => {
      if (REACHING.has(n)) h.set?.(n, v);
    },
    open: (v) => h.open?.(v),
    flows: () => {},
  };
  const shownIn = (m: THREE.Object3D) => {
    for (let o: THREE.Object3D | null = m; o && o !== model.group; o = o.parent) if (!o.visible) return false;
    return true;
  };
  let r = 0;
  for (; t < 24; t += 0.5) {
    if (act) act(rig);
    else rig.open(cycle(t, 9));
    model.group.updateMatrixWorld(true);
    for (const m of meshes) {
      if (!shownIn(m)) continue;
      const g = m.geometry;
      if (!g.boundingBox) g.computeBoundingBox();
      tmp.copy(g.boundingBox!).applyMatrix4(m.matrixWorld);
      if (tmp.isEmpty()) continue;
      into.union(tmp);
      r = Math.max(r, Math.hypot(Math.max(Math.abs(tmp.min.x - cx), Math.abs(tmp.max.x - cx)), Math.max(Math.abs(tmp.min.z - cz), Math.abs(tmp.max.z - cz))));
    }
  }
  return r;
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Card text                                                                                                      */
/* ------------------------------------------------------------------------------------------------------------ */

const FINISH: Record<string, string> = {
  ss: "AISI 304 No. 4 stainless",
  ssDull: "galvanised steel",
  ssBright: "polished stainless",
  chrome: "chrome wire",
  alu: "anodised aluminium",
  castIron: "seasoned cast iron",
  griddle: "seasoned carbon steel",
  enamel: "black enamel",
  charcoal: "powder coat",
  panel: "insulated panel skin",
  hdpe: "HDPE",
  glass: "tempered glass",
};

/** FrankenSim's material and physics cards, in plain words. The card id shows on hover. */
const CARD_NAME: Record<string, string> = {
  "surfaces/cookware-ss304": "304 stainless steel",
  "surfaces/cookware-carbon-steel": "Carbon steel",
  "surfaces/cookware-cast-iron-seasoned": "Seasoned cast iron",
  "surfaces/cookware-ptfe-nonstick": "PTFE non-stick",
  "surfaces/cookware-copper": "Copper",
  "surfaces/cookware-aluminum-anodized": "Anodised aluminium",
  "surfaces/doorway-gosney-olama-infiltration": "Warm air through the doorway (Gosney–Olama)",
  "surfaces/atkins-blade-food-cutting": "Blade cutting food (Atkins)",
  "surfaces/lumped-cooling-coil-walkin": "Walk-in cooling coil",
  "surfaces/lumped-cooling-coil-blast-chiller": "Blast-chiller cooling coil",
  "kinetics/maillard-browning": "Maillard browning",
  "kinetics/meat-collagen-gelatin-hydrolysis": "Collagen turning to gelatin",
  "kinetics/egg-albumen-two-stage": "Egg white setting in two stages",
  "kinetics/frying-oil-degradation-kinetics": "Frying oil breaking down",
  "kinetics/yeast-rosso-fermentation": "Yeast fermentation (Rosso)",
  "kinetics/crust-glass-transition-gordon-taylor": "Crust softening (Gordon–Taylor)",
};

/** "knob_0, knob_1, oven_door" → "oven_door, knob ×2". */
function jointSummary(names: string[]): string {
  const groups = new Map<string, number>();
  for (const n of names) {
    const k = n.replace(/_?\d+(_\d+)?$/, "");
    groups.set(k, (groups.get(k) ?? 0) + 1);
  }
  return [...groups].map(([k, c]) => (c > 1 ? `${k} ×${c}` : k)).join(", ");
}

const mm = (m: number) => String(Math.round(m * 1000));

/* ------------------------------------------------------------------------------------------------------------ */
/* Overlay styles (the chapter's own elements inside the stage)                                                   */
/* ------------------------------------------------------------------------------------------------------------ */

const OVERLAY_CSS = `
.eqx-overlay{position:absolute;inset:0;pointer-events:none;overflow:hidden;font-family:"JetBrains Mono Variable",ui-monospace,SFMono-Regular,Menlo,monospace;contain:strict}
.eqx-tag{position:absolute;left:0;top:0;will-change:transform;white-space:nowrap}
.eqx-dim{padding:3px 8px;border-radius:7px;background:rgba(14,10,8,.82);border:1px solid rgba(245,165,36,.6);color:#fde7b0;font-size:11px;font-weight:600;letter-spacing:.02em;line-height:1.2;box-shadow:0 2px 10px rgba(0,0,0,.45)}
.eqx-dim i{font-style:normal;color:#f59e0b;margin-right:5px;font-weight:700}
.eqx-sock{display:flex;align-items:center;gap:5px;font-size:9.5px;color:#a5f3fc;letter-spacing:.02em;text-shadow:0 1px 3px #000,0 0 6px #000}
.eqx-sock{transform-origin:0 50%}.eqx-sock.dot span{display:none}
.eqx-sock::before{content:"";width:9px;height:9px;border-radius:50%;background:#22d3ee;box-shadow:0 0 0 2px rgba(8,6,5,.85),0 0 10px rgba(34,211,238,.9);flex:none}
.eqx-panel{position:absolute;left:12px;top:12px;max-width:min(300px,calc(100% - 24px));padding:9px 11px 10px;border-radius:11px;background:rgba(12,9,7,.8);border:1px solid rgba(255,255,255,.09);color:#e7e2dc;font-size:10.5px;line-height:1.45;backdrop-filter:blur(6px);-webkit-backdrop-filter:blur(6px);box-shadow:0 6px 24px rgba(0,0,0,.4);opacity:0;transform:translateY(-4px);transition:opacity .25s ease,transform .25s ease}
.eqx-panel.on{opacity:1;transform:none}
.eqx-panel h4{margin:0 0 6px;font-size:10px;letter-spacing:.14em;text-transform:uppercase;color:#fbbf24;font-weight:700}
.eqx-bar{position:relative;height:9px;border-radius:5px;margin:2px 0 15px}
.eqx-bar span{position:absolute;top:12px;transform:translateX(-50%);font-size:9px;color:#b8afa6}
.eqx-bar b{position:absolute;top:-3px;width:2px;height:15px;margin-left:-1px;background:#fff;box-shadow:0 0 0 1px rgba(0,0,0,.6);border-radius:1px}
.eqx-zone{display:flex;justify-content:space-between;gap:10px;color:#f3eee8}
.eqx-zone em{font-style:normal;color:#9c938a}
.eqx-fine{margin-top:5px;color:#8f877f;font-size:9.5px}
.eqx-key{display:inline-flex;align-items:center;gap:6px;margin-right:12px}
.eqx-key::before{content:"";width:10px;height:10px;border-radius:3px;background:var(--k)}
.eqx-hint{position:absolute;right:12px;top:12px;padding:5px 9px;border-radius:999px;background:rgba(12,9,7,.62);border:1px solid rgba(255,255,255,.08);color:#b9b0a6;font-size:10px;letter-spacing:.03em;transition:opacity .6s ease}
.eqx-hint.off{opacity:0}
@media (max-width:600px){.eqx-panel{font-size:9.5px;padding:7px 9px;left:8px;top:8px;max-width:min(232px,calc(100% - 16px))}.eqx-panel .eqx-fine:not(.eqx-honest){display:none}.eqx-panel .eqx-bar span:nth-child(even){display:none}.eqx-hint{display:none}.eqx-dim{font-size:10px;padding:2px 6px}.eqx-sock{font-size:8.5px}}
`;

function injectCss(): void {
  if (document.getElementById("eqx-css")) return;
  const s = document.createElement("style");
  s.id = "eqx-css";
  s.textContent = OVERLAY_CSS;
  document.head.appendChild(s);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The chapter                                                                                                    */
/* ------------------------------------------------------------------------------------------------------------ */

interface Shown {
  arch: Arch;
  model: Model;
  /** The fixture frame on the turntable: the model and its cage live in it. */
  frame: THREE.Group;
  meshes: THREE.Mesh[];
  /** Original materials, to restore after thermal or X-ray. */
  orig: Map<THREE.Mesh, THREE.Material | THREE.Material[]>;
  jointMeshes: Set<THREE.Mesh>;
  edges: THREE.LineSegments[] | null;
  hidden: THREE.Object3D[];
  memo: Map<string, number>;
  act: Act | null;
  /** Bounding sphere about the turntable axis: centre height and radius. */
  cy: number;
  R: number;
  /** Platter radius. */
  rp: number;
  /** The model swept round the axis, its moving parts at their full travel too: radius, lowest and highest points. */
  rh: number;
  bot: number;
  top: number;
  /** Camera elevation for this model (radians); below zero for a piece hung overhead, seen from beneath. */
  elev: number;
  /** Hung overhead (hoods, heat lamps): no platter under it, and the camera looks up at its working face. */
  over: boolean;
  born: number;
}

export function mountEquipmentTurntable(root: HTMLElement): void {
  const stageEl = root.querySelector<HTMLElement>("[data-eq-stage]");
  const canvas = root.querySelector<HTMLCanvasElement>("[data-eq-canvas]");
  if (!stageEl || !canvas) return;
  const maybeRenderer = makeRenderer(canvas, { alpha: true, shadows: true });
  if (!maybeRenderer) {
    const note = document.createElement("p");
    note.className = "webgl-fallback";
    note.textContent = "The 3D showroom needs WebGL, which this browser has turned off.";
    stageEl.appendChild(note);
    return;
  }
  const renderer: THREE.WebGLRenderer = maybeRenderer;
  renderer.setClearColor(0x000000, 0);
  renderer.toneMappingExposure = 1.05;
  injectCss();

  const reduceMotion = window.matchMedia?.("(prefers-reduced-motion: reduce)").matches ?? false;

  /* --- Scene, studio, platter ------------------------------------------------------------------------------- */
  const scene = new THREE.Scene();
  const envMode = new URLSearchParams(location.search).get("eqenv");
  scene.environment = envMode === "dark" ? showroomEnvironment(renderer) : studioEnvironment(renderer);
  scene.environmentIntensity = envMode === "dark" ? 1.0 : 0.85;
  const camera = new THREE.PerspectiveCamera(30, 16 / 10, 0.01, 100);
  camera.position.set(0, 1.2, 3.4);
  const gate = shaderGate(renderer, camera);

  const hemi = new THREE.HemisphereLight(0xfff1e2, 0x2a2018, 0.35);
  const key = new THREE.DirectionalLight(0xfff0de, 2.1);
  key.castShadow = true;
  key.shadow.mapSize.set(2048, 2048);
  key.shadow.bias = -0.00025;
  const rim = new THREE.DirectionalLight(0xa9c8ff, 0.9);
  const fill = new THREE.DirectionalLight(0xffe2c4, 0.35);
  // Lights the underside of a piece hung overhead, the way the line below it would.
  const under = new THREE.DirectionalLight(0xfff3e6, 0);
  scene.add(hemi, key, key.target, rim, fill, under, under.target);

  const turntable = new THREE.Group();
  scene.add(turntable);

  const p = palette();
  const platter = new THREE.Group();
  const topTex = platterTopTexture();
  const platterTopMat = new THREE.MeshStandardMaterial({ color: 0xffffff, map: topTex, metalness: 0.5, roughness: 0.46 });
  const platterTop = new THREE.Mesh(new THREE.CircleGeometry(0.965, 128).rotateX(-Math.PI / 2), platterTopMat);
  platterTop.receiveShadow = true;
  const rimGeo = new THREE.LatheGeometry(
    // Bottom centre up to the top edge: LatheGeometry faces outward only when the profile climbs.
    [
      [0.0, -0.04],
      [0.99, -0.04],
      [1.0, -0.034],
      [1.0, -0.014],
      [0.985, -0.004],
      [0.965, 0],
    ].map(([r, y]) => new THREE.Vector2(r, y)),
    128,
  );
  // A satin turned steel edge.
  const rimMat = new THREE.MeshStandardMaterial({ color: 0xd9dde2, metalness: 1, roughness: 0.26 });
  const platterRim = new THREE.Mesh(rimGeo, rimMat);
  platterRim.receiveShadow = true;
  platterRim.castShadow = true;
  const baseGeo = new THREE.CylinderGeometry(0.9, 0.93, 0.03, 96);
  baseGeo.translate(0, -0.055, 0);
  const platterBase = new THREE.Mesh(baseGeo, p.enamel);
  platterBase.receiveShadow = true;
  platter.add(platterTop, platterRim, platterBase);
  turntable.add(platter);

  const aoTex = contactTexture();
  const aoMat = new THREE.MeshBasicMaterial({ map: aoTex, transparent: true, opacity: 0.7, depthWrite: false, polygonOffset: true, polygonOffsetFactor: -2, polygonOffsetUnits: -4, toneMapped: false });
  const aoPad = new THREE.Mesh(new THREE.PlaneGeometry(1, 1).rotateX(-Math.PI / 2), aoMat);
  aoPad.renderOrder = 1;
  turntable.add(aoPad);

  // The studio floor: catches the platter's shadow, fades out to the stage.
  const floorMat = new THREE.ShadowMaterial({ color: 0x000000, opacity: 0.38, depthWrite: false });
  const floor = new THREE.Mesh(new THREE.CircleGeometry(1, 64).rotateX(-Math.PI / 2), floorMat);
  floor.receiveShadow = true;
  scene.add(floor);
  const poolTex = canvasTexture(512, (x, S) => {
    const g = x.createRadialGradient(S / 2, S / 2, 0, S / 2, S / 2, S / 2);
    g.addColorStop(0, "rgba(255,214,170,0.13)");
    g.addColorStop(0.45, "rgba(255,200,150,0.05)");
    g.addColorStop(1, "rgba(255,200,150,0)");
    x.fillStyle = g;
    x.fillRect(0, 0, S, S);
  });
  const pool = new THREE.Mesh(new THREE.PlaneGeometry(1, 1).rotateX(-Math.PI / 2), new THREE.MeshBasicMaterial({ map: poolTex, transparent: true, depthWrite: false, toneMapped: false, blending: THREE.AdditiveBlending }));
  pool.renderOrder = -1;
  scene.add(pool);

  /* --- Inspection materials (thermal and X-ray) ------------------------------------------------------------- */
  const ramp = rampUniforms();
  const makeInspect = (ghost: boolean, side: THREE.Side) =>
    new THREE.ShaderMaterial({
      vertexShader: INSPECT_VERT,
      fragmentShader: INSPECT_FRAG,
      side,
      transparent: ghost,
      depthWrite: !ghost,
      blending: ghost ? THREE.AdditiveBlending : THREE.NormalBlending,
      toneMapped: false,
      uniforms: {
        uRootInv: { value: new THREE.Matrix4() },
        uZA: { value: Array.from({ length: N_ZONES }, () => new THREE.Vector4()) },
        uZB: { value: Array.from({ length: N_ZONES }, () => new THREE.Vector4()) },
        uZT: { value: new Array(N_ZONES).fill(AMBIENT_C) },
        uZN: { value: 0 },
        uAmb: { value: AMBIENT_C },
        uLambda: { value: 0.05 },
        uRampC: { value: ramp.c },
        uRampT: { value: ramp.t },
        uGhost: { value: ghost ? 1 : 0 },
        uThermal: { value: 1 },
        uTint: { value: new THREE.Color(0xf59e0b) },
      },
    });
  const thermalMats = { front: makeInspect(false, THREE.FrontSide), double: makeInspect(false, THREE.DoubleSide) };
  const ghostBody = makeInspect(true, THREE.DoubleSide);
  const ghostJoint = makeInspect(true, THREE.DoubleSide);
  (ghostBody.uniforms.uTint!.value as THREE.Color).set(0xf59e0b);
  (ghostJoint.uniforms.uTint!.value as THREE.Color).set(0x22d3ee);
  const inspectMats = [thermalMats.front, thermalMats.double, ghostBody, ghostJoint];
  const edgeBody = new THREE.LineBasicMaterial({ color: 0xffb347, transparent: true, opacity: 0.42, depthWrite: false, blending: THREE.AdditiveBlending, toneMapped: false });
  const edgeJoint = new THREE.LineBasicMaterial({ color: 0x5eead4, transparent: true, opacity: 0.7, depthWrite: false, blending: THREE.AdditiveBlending, toneMapped: false });

  /* --- Overlay ------------------------------------------------------------------------------------------------ */
  const overlay = document.createElement("div");
  overlay.className = "eqx-overlay";
  overlay.setAttribute("aria-hidden", "true");
  canvas.insertAdjacentElement("afterend", overlay);
  const panel = document.createElement("div");
  panel.className = "eqx-panel";
  overlay.append(panel);
  const hint = stageEl.querySelector<HTMLElement>(".stage-hint") ?? document.createElement("p");
  hint.style.transition = "opacity .6s ease";
  const dimTags = ["W", "D", "H"].map((k) => {
    const el = document.createElement("div");
    el.className = "eqx-tag eqx-dim";
    el.innerHTML = `<i>${k}</i><span></span>`;
    overlay.appendChild(el);
    return el;
  });
  let sockTags: { el: HTMLElement; at: THREE.Vector3; w: number }[] = [];

  /* --- DOM hooks ---------------------------------------------------------------------------------------------- */
  const tabBtns = Array.from(root.querySelectorAll<HTMLButtonElement>("[data-eq-tab]"));
  const pillGrid = root.querySelector<HTMLElement>("[data-eq-pills]");
  const searchInput = root.querySelector<HTMLInputElement>("[data-eq-search]");
  const select = root.querySelector<HTMLSelectElement>("[data-eq-select]");
  const btn = {
    joints: root.querySelector<HTMLButtonElement>("[data-eq-joints]"),
    thermal: root.querySelector<HTMLButtonElement>("[data-eq-thermal]"),
    wire: root.querySelector<HTMLButtonElement>("[data-eq-wire]"),
    cage: root.querySelector<HTMLButtonElement>("[data-eq-cage]"),
    spin: root.querySelector<HTMLButtonElement>("[data-eq-spin]"),
  };
  const meta = {
    name: root.querySelector<HTMLElement>("[data-eq-name]"),
    id: root.querySelector<HTMLElement>("[data-eq-id]"),
    dims: root.querySelector<HTMLElement>("[data-eq-dims]"),
    tris: root.querySelector<HTMLElement>("[data-eq-tris]"),
    card: root.querySelector<HTMLElement>("[data-eq-card]"),
    sockets: root.querySelector<HTMLElement>("[data-eq-sockets]"),
    note: root.querySelector<HTMLElement>("[data-eq-note]"),
  };
  // Defaults come from the markup's aria-pressed, where it says.
  const pressed = (b: HTMLButtonElement | null, dflt: boolean) => (b?.getAttribute("aria-pressed") ?? String(dflt)) === "true";
  const state = {
    animate: pressed(btn.joints, true) && !reduceMotion,
    thermal: pressed(btn.thermal, false),
    wire: pressed(btn.wire, false),
    cage: pressed(btn.cage, true),
    spin: pressed(btn.spin, true) && !reduceMotion,
  };
  const syncButtons = () => {
    btn.joints?.setAttribute("aria-pressed", String(state.animate));
    btn.thermal?.setAttribute("aria-pressed", String(state.thermal));
    btn.wire?.setAttribute("aria-pressed", String(state.wire));
    btn.cage?.setAttribute("aria-pressed", String(state.cage));
    btn.spin?.setAttribute("aria-pressed", String(state.spin));
  };
  syncButtons();

  /* --- Shared resources that stay for the page's life ------------------------------------------------------- */
  const keep = new WeakSet<object>();
  for (const m of Object.values(p) as THREE.Material[]) {
    keep.add(m);
    for (const t of texturesOf(m)) keep.add(t);
  }

  /* --- Cage --------------------------------------------------------------------------------------------------- */
  const cage = new THREE.Group();
  const cageMats: LineMaterial[] = [];
  const lineMat = (color: number, width: number, opacity: number, dashed = false, dash = 0.02) => {
    const m = new LineMaterial({ color, linewidth: width, transparent: true, opacity, dashed, dashSize: dash, gapSize: dash * 0.7, worldUnits: false, toneMapped: false });
    m.resolution.set(canvas.clientWidth || 800, canvas.clientHeight || 500);
    cageMats.push(m);
    return m;
  };
  /** One set of dimension-line variants, shown one at a time, picked to face the camera. */
  interface DimSet {
    variants: { obj: THREE.Group; label: THREE.Vector3; side: THREE.Vector3 }[];
    tag: HTMLElement;
    /** "facing": the variant on the side toward the camera; "edge": the one standing furthest right on screen. */
    pick: "facing" | "edge";
  }
  let dims: DimSet[] = [];
  let stanceRings: LineSegments2[] = [];

  function segs(pts: number[]): LineSegmentsGeometry {
    const g = new LineSegmentsGeometry();
    g.setPositions(pts);
    return g;
  }

  function clearCage(): void {
    for (const c of [...cage.children]) {
      cage.remove(c);
      c.traverse((o) => (o as LineSegments2).geometry?.dispose());
    }
    for (const m of cageMats) m.dispose();
    cageMats.length = 0;
    dims = [];
    stanceRings = [];
    for (const s of sockTags) s.el.remove();
    sockTags = [];
  }

  /** The cage in the fixture frame; (cx, cz) is where the turntable's axis passes, rp the platter's radius. */
  function buildCage(a: Arch, R: number, cx: number, cz: number, rp: number): void {
    clearCage();
    const [w, d, h] = a.ext;
    const s = Math.max(w, d, h);
    const off = Math.max(0.035, 0.09 * Math.max(w, d) ** 0.85) * Math.min(1, 0.4 + s);
    const tk = off * 0.32;
    const y = Math.max(0.0006, R * 0.002);
    // The envelope: dashed edges.
    const P = [
      [0, 0, -d],
      [w, 0, -d],
      [w, 0, 0],
      [0, 0, 0],
      [0, h, -d],
      [w, h, -d],
      [w, h, 0],
      [0, h, 0],
    ];
    const E = [
      [0, 1], [1, 2], [2, 3], [3, 0],
      [4, 5], [5, 6], [6, 7], [7, 4],
      [0, 4], [1, 5], [2, 6], [3, 7],
    ];
    const boxPts: number[] = [];
    for (const [i, j] of E) {
      const p0 = P[i!]!;
      const p1 = P[j!]!;
      boxPts.push(p0[0]!, Math.max(y, p0[1]!), p0[2]!, p1[0]!, Math.max(y, p1[1]!), p1[2]!);
    }
    const env = new LineSegments2(segs(boxPts), lineMat(0xf5b544, 1.25, 0.55, true, s * 0.018));
    env.computeLineDistances();
    cage.add(env);

    const dimMat = lineMat(0xfbbf24, 1.6, 0.95);
    const extMat = lineMat(0xfbbf24, 1.1, 0.55);
    const mk = (main: number[], ext: number[], label: THREE.Vector3, side: THREE.Vector3) => {
      const g = new THREE.Group();
      g.add(new LineSegments2(segs(main), dimMat), new LineSegments2(segs(ext), extMat));
      cage.add(g);
      return { obj: g, label, side };
    };
    const slash = (x: number, yy: number, z: number, ax: "x" | "y" | "z") => {
      const t = tk * 0.55;
      if (ax === "y") return [x - t, yy - t, z + t, x + t, yy + t, z - t];
      return [x - t, yy, z + t, x + t, yy, z - t];
    };
    // Width: along x, in front (z = off) or behind (z = -d - off).
    const wSet: DimSet = { variants: [], tag: dimTags[0]!, pick: "facing" };
    for (const sz of [1, -1]) {
      const z = sz > 0 ? off : -d - off;
      const z0 = sz > 0 ? 0 : -d;
      wSet.variants.push(
        mk([0, y, z, w, y, z, ...slash(0, y, z, "x"), ...slash(w, y, z, "x")], [0, y, z0 + sz * tk * 0.3, 0, y, z + sz * tk, w, y, z0 + sz * tk * 0.3, w, y, z + sz * tk], new THREE.Vector3(w / 2, y, z + sz * tk * 1.8), new THREE.Vector3(0, 0, sz)),
      );
    }
    // Depth: along z, to the right (x = w + off) or left (x = -off).
    const dSet: DimSet = { variants: [], tag: dimTags[1]!, pick: "facing" };
    for (const sx of [1, -1]) {
      const x = sx > 0 ? w + off : -off;
      const x0 = sx > 0 ? w : 0;
      dSet.variants.push(
        mk([x, y, -d, x, y, 0, ...slash(x, y, -d, "z"), ...slash(x, y, 0, "z")], [x0 + sx * tk * 0.3, y, -d, x + sx * tk, y, -d, x0 + sx * tk * 0.3, y, 0, x + sx * tk, y, 0], new THREE.Vector3(x + sx * tk * 1.8, y, -d / 2), new THREE.Vector3(sx, 0, 0)),
      );
    }
    // Height: a vertical line standing off one of the four upright edges, diagonally out from the box; the one at
    // the silhouette's right is shown.
    const hSet: DimSet = { variants: [], tag: dimTags[2]!, pick: "edge" };
    for (const [cx, cz, sx, sz] of [
      [w, 0, 1, 1],
      [0, 0, -1, 1],
      [w, -d, 1, -1],
      [0, -d, -1, -1],
    ] as const) {
      const k = off * 0.72;
      const x = cx + sx * k;
      const z = cz + sz * k;
      const ex = (yy: number) => [cx + sx * tk * 0.25, yy, cz + sz * tk * 0.25, x + sx * tk * 0.7, yy, z + sz * tk * 0.7];
      hSet.variants.push(mk([x, 0, z, x, h, z, ...slash(x, 0, z, "y"), ...slash(x, h, z, "y")], [...ex(h), ...ex(y)], new THREE.Vector3(x + sx * tk * 1.2, h / 2, z + sz * tk * 1.2), new THREE.Vector3(cx, 0, cz)));
    }
    dims = [wSet, dSet, hSet];
    dimTags[0]!.querySelector("span")!.textContent = `${mm(w)} mm`;
    dimTags[1]!.querySelector("span")!.textContent = `${mm(d)} mm`;
    dimTags[2]!.querySelector("span")!.textContent = `${mm(h)} mm`;

    // Sockets: tags at their manifest positions; a stance anchor also gets a ring on the floor where the cook stands.
    const merged = new Map<string, { at: THREE.Vector3; names: string[] }>();
    for (const [sid, kind, sx, sy, sz] of a.sockets) {
      // A stance spot out beyond the platter (where the cook stands) would sit outside the frame: the card lists it.
      if (Math.hypot(sx - cx, sz - cz) > rp * 1.04) continue;
      const k = `${sx.toFixed(3)},${sy.toFixed(3)},${sz.toFixed(3)}`;
      const e = merged.get(k) ?? { at: new THREE.Vector3(sx, Math.max(sy, y), sz), names: [] };
      e.names.push(sid);
      merged.set(k, e);
      if (kind === "stance_anchor") {
        const pts: number[] = [];
        const rr = 0.16;
        for (let i = 0; i < 48; i++) {
          const a0 = (i / 48) * Math.PI * 2;
          const a1 = ((i + 1) / 48) * Math.PI * 2;
          pts.push(sx + Math.cos(a0) * rr, y, sz + Math.sin(a0) * rr, sx + Math.cos(a1) * rr, y, sz + Math.sin(a1) * rr);
        }
        const ring = new LineSegments2(segs(pts), lineMat(0x22d3ee, 1.4, 0.7, true, 0.03));
        ring.computeLineDistances();
        cage.add(ring);
        stanceRings.push(ring);
      }
    }
    for (const { at, names } of merged.values()) {
      const el = document.createElement("div");
      el.className = "eqx-tag eqx-sock";
      const t = document.createElement("span");
      t.textContent = names.join(" / ").replace(/_/g, " ");
      el.appendChild(t);
      overlay.appendChild(el);
      sockTags.push({ el, at, w: 0 });
    }
  }

  /* --- State ---------------------------------------------------------------------------------------------- */
  let shown: Shown | null = null;
  /** Dev only: hold the model as it stands (no choreography), for page checks. */
  let debugHold = false;
  let activeTab: EquipmentTabId = "hot_line";
  let activeId = FIRST_ID;
  let inView = false;
  let dragging = false;
  let lastInteract = -1e9;
  let clock = 0;
  let fitDist = 3;
  const home = new THREE.Vector3();
  let freeH = 1;
  let viewShiftPx = 0;
  let introT = 1;
  let focusAnim: { from: THREE.Vector3; to: THREE.Vector3; d0: number; d1: number; t: number } | null = null;

  const controls = createScrollSafeOrbit(camera, canvas, [0, 0.5, 0], { minDistance: 0.05, maxDistance: 30, maxPolarAngle: Math.PI * 0.485 });
  // OrbitControls writes touch-action: none onto the canvas, which on a phone traps the page's scroll under the
  // stage. Hand vertical swipes back to the page (as the stylesheet intends); sideways drags still turn the model
  // and two-finger pinches still zoom it.
  canvas.style.touchAction = "pan-y";
  controls.zoomToCursor = true;
  controls.enablePan = true;
  controls.screenSpacePanning = true;
  controls.zoomSpeed = 1.1;
  controls.rotateSpeed = 0.8;
  controls.addEventListener("start", () => {
    dragging = true;
    lastInteract = clock;
    focusAnim = null;
    hint.style.opacity = "0";
  });
  controls.addEventListener("end", () => {
    dragging = false;
    lastInteract = clock;
  });

  /* --- Fit ------------------------------------------------------------------------------------------------ */
  /** Measure the space the control chips leave free at the stage's foot, and shift the view to centre in it. */
  function measureFree(): void {
    const sr = stageEl!.getBoundingClientRect();
    let top = sr.bottom;
    for (const b of Object.values(btn)) {
      if (!b) continue;
      const r = b.getBoundingClientRect();
      if (r.height > 0 && r.top > sr.top + sr.height * 0.45 && r.bottom <= sr.bottom + 2) top = Math.min(top, r.top);
    }
    const H = Math.max(1, sr.height);
    const free = Math.max(H * 0.55, top - sr.top - 6);
    freeH = free / H;
    viewShiftPx = (H - free) / 2;
  }

  /** Resizing waits on the GPU: behind a shader compile it is held until the compile is done. */
  function resize(): void {
    gate.defer("resize", resizeNow);
  }
  function resizeNow(): void {
    const w = Math.max(1, Math.floor(stageEl!.clientWidth));
    const h = Math.max(1, Math.floor(stageEl!.clientHeight));
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    measureFree();
    if (viewShiftPx > 0.5) camera.setViewOffset(w, h, 0, viewShiftPx, w, h);
    else camera.clearViewOffset();
    camera.updateProjectionMatrix();
    for (const m of cageMats) m.resolution.set(w, h);
    if (shown) {
      const f = fitView(shown);
      fitDist = f.dist;
      home.set(0, f.y, 0);
    }
  }

  /**
   * The distance and pivot height at which all the turntable can show fits the view from the model's elevation:
   * the model swept a full turn round the axis (a cylinder through its farthest corner, its doors, drawers and lids
   * at full travel) and the platter under it.
   * Sizes differ by three orders of magnitude, from a paring knife to a walk-in; each fills the frame alike.
   */
  const ringA = Array.from({ length: 36 }, (_, i) => (i / 36) * Math.PI * 2);
  function fitView(s: Shown): { dist: number; y: number } {
    const half = THREE.MathUtils.degToRad(camera.fov / 2);
    const tanV = Math.tan(half) * freeH * 0.93;
    const tanH = Math.tan(half) * camera.aspect * 0.93;
    const ce = Math.cos(s.elev);
    const se = Math.sin(s.elev);
    const pl = s.rp * 0.042;
    // A hung piece has no platter round it to pad the frame, so leave its dimension lines room instead.
    const rr = s.over ? s.rh * 1.12 : s.rh;
    const rings: [number, number][] = [
      [rr, s.top],
      [rr, s.over ? s.bot : Math.max(s.bot, 0)],
    ];
    if (!s.over) rings.push([s.rp, 0], [s.rp * 0.985, -pl]);
    let y = (s.top + s.bot) / 2;
    let dist = 1;
    for (let it = 0; it < 4; it++) {
      dist = 0;
      for (const [r, py] of rings) for (const a of ringA) {
        const qx = r * Math.cos(a);
        const qy = py - y;
        const qz = r * Math.sin(a);
        const up = qy * ce - qz * se;
        const back = qy * se + qz * ce;
        dist = Math.max(dist, back + Math.abs(qx) / tanH, back + Math.abs(up) / tanV);
      }
      // Centre the vertical extent: slide the pivot up the axis by the extent's offset.
      let lo = Infinity;
      let hi = -Infinity;
      for (const [r, py] of rings) for (const a of ringA) {
        const qy = py - y;
        const qz = r * Math.sin(a);
        const t = (qy * ce - qz * se) / (dist - (qy * se + qz * ce));
        lo = Math.min(lo, t);
        hi = Math.max(hi, t);
      }
      y += (((hi + lo) / 2) * dist) / Math.max(0.3, ce);
    }
    return { dist, y };
  }

  function frame(s: Shown, animate: boolean): void {
    const f = fitView(s);
    fitDist = f.dist;
    const elev = s.elev;
    home.set(0, f.y, 0);
    controls.target.copy(home);
    const dist = fitDist * (animate && !reduceMotion ? 1.18 : 1);
    camera.position.set(0, f.y + Math.sin(elev) * dist, Math.cos(elev) * dist);
    controls.minDistance = Math.max(0.012, s.R * 0.03);
    controls.maxDistance = fitDist * 3;
    // A floor piece is seen from above its platter; a hung one may be walked round underneath as well.
    controls.minPolarAngle = s.over ? Math.PI * 0.2 : 0;
    controls.maxPolarAngle = s.over ? Math.PI * 0.8 : Math.PI * 0.485;
    controls.update();
    introT = animate && !reduceMotion ? 0 : 1;
  }

  /* --- Build and show --------------------------------------------------------------------------------------- */
  function lightFor(s: Shown): void {
    // The shadow frustum holds the platter and anything the moving parts swing out past it.
    const R = Math.max(s.rp, s.rh) * 1.1;
    const c = new THREE.Vector3(0, s.cy * 0.5, 0);
    key.position.copy(c).add(new THREE.Vector3(R * 0.8, R * 2.8, R * 1.0));
    key.target.position.copy(c);
    const sc = key.shadow.camera;
    sc.left = sc.bottom = -R * 1.2;
    sc.right = sc.top = R * 1.2;
    sc.near = R * 0.5;
    sc.far = R * 6;
    sc.updateProjectionMatrix();
    key.shadow.normalBias = Math.max(0.0008, 0.004 * R);
    rim.position.set(-R * 1.8, R * 1.2, -R * 1.6);
    fill.position.set(-R * 1.4, R * 0.6, R * 1.6);
    under.position.set(R * 0.6, -R * 2.4, R * 1.4);
    under.target.position.copy(c);
    under.intensity = s.over ? 0.55 : 0;
    // Seen from below, metal mirrors the studio's bright ceiling and washes out; turn the studio over for a hung
    // piece, so that its faces read as stainless and its underside catches the softboxes.
    scene.environmentRotation.set(s.over ? Math.PI : 0, 0, 0);
    for (const o of [platter, aoPad, floor, pool]) o.visible = !s.over;
    floor.scale.setScalar(s.rp * 2.4);
    floor.position.y = -0.042 * s.rp - 0.0005;
    pool.scale.setScalar(s.rp * 4.2);
    pool.position.y = floor.position.y - 0.0005;
  }

  function disposeShown(s: Shown): void {
    restoreMaterials(s);
    if (s.edges) for (const e of s.edges) {
      e.parent?.remove(e);
      e.geometry.dispose();
    }
    turntable.remove(s.frame);
    const res = resourcesOf(s.model.group);
    for (const g of res.geos) {
      if (keep.has(g) || sharedRes.has(g)) continue;
      g.dispose();
      disposedOnce.add(g);
    }
    for (const m of res.mats) {
      if (keep.has(m) || sharedRes.has(m) || inspectMats.includes(m as THREE.ShaderMaterial)) continue;
      m.dispose();
      disposedOnce.add(m);
    }
    for (const t of res.texs) {
      if (keep.has(t) || sharedRes.has(t)) continue;
      t.dispose();
      disposedOnce.add(t);
    }
    for (const im of res.inst) im.dispose();
  }

  function show(id: string): void {
    const a = BY_ID.get(id);
    if (!a) return;
    activeId = id;
    const [w, d, h] = a.ext;
    const model = buildFixture(id, { w, d, h }, 2, true);
    // Anything this build shares with one already let go is a family cache: keep it from now on.
    const res = resourcesOf(model.group);
    for (const set of [res.geos, res.mats, res.texs] as Set<object>[]) for (const r of set) if (disposedOnce.has(r)) sharedRes.add(r);

    if (shown) disposeShown(shown);

    // Measure what was built (flames and sprites aside), and centre its footprint on the turntable's axis.
    model.group.updateMatrixWorld(true);
    const box = new THREE.Box3();
    const meshes: THREE.Mesh[] = [];
    model.group.traverse((o) => {
      const m = o as THREE.Mesh;
      if (!m.isMesh || (o as THREE.InstancedMesh).isInstancedMesh || !m.geometry) return;
      meshes.push(m);
      m.geometry.computeBoundingBox();
      box.union(m.geometry.boundingBox!.clone().applyMatrix4(m.matrixWorld));
    });
    if (box.isEmpty()) box.set(new THREE.Vector3(0, 0, -d), new THREE.Vector3(w, h, 0));
    // The archetype's envelope counts too, so that the cage always fits.
    box.union(new THREE.Box3(new THREE.Vector3(0, 0, -d), new THREE.Vector3(w, h, 0)));
    const cx = (box.min.x + box.max.x) / 2;
    const cz = (box.min.z + box.max.z) / 2;
    const rh = 0.5 * Math.hypot(box.max.x - box.min.x, box.max.z - box.min.z);
    const hh = (box.max.y - Math.min(0, box.min.y)) / 2;
    const R = Math.hypot(rh, hh);
    const lift = Math.max(0.00002, R * 0.00015);
    // The moving parts at their full travel (a door swung open, a trolley rolled out, a pan tipped) must stay in the
    // frame as the platter turns: run the choreography through its cycles and take the farthest reach about the axis.
    const act = CHOREO.find(([re]) => re.test(id))?.[1] ?? null;
    const swept = box.clone();
    const reach = sweptReach(model, meshes, act, cx, cz, swept);

    const frameG = new THREE.Group();
    frameG.position.set(-cx, lift, -cz);
    frameG.add(model.group, cage);
    turntable.add(frameG);
    // A three-quarter view from the front right; a long, narrow piece (a knife, a ladle) starts side-on instead,
    // its handle toward the right, so that it fills the frame.
    const bw = box.max.x - box.min.x;
    const bd = box.max.z - box.min.z;
    turntable.rotation.y = bd > 1.6 * bw ? Math.PI / 2 - 0.32 : bw > 3.2 * bd && bw > 2 ? -0.42 : -0.62;

    const jointSet = new Set<THREE.Object3D>(model.joints.values());
    const jointMeshes = new Set<THREE.Mesh>();
    for (const m of meshes) {
      for (let o: THREE.Object3D | null = m.parent; o && o !== model.group; o = o.parent) if (jointSet.has(o)) {
        jointMeshes.add(m);
        break;
      }
    }
    const rp = Math.max(rh * 1.14, Math.max(w, d) * 0.6) + Math.min(0.08, R * 0.04);
    const over = OVERHEAD.test(id);
    platter.scale.setScalar(rp);
    aoPad.scale.set((box.max.x - box.min.x) * 1.22, 1, (box.max.z - box.min.z) * 1.22);
    aoPad.position.set(0, 0, 0);

    const s: Shown = {
      arch: a,
      model,
      frame: frameG,
      meshes,
      orig: new Map(),
      jointMeshes,
      edges: null,
      hidden: [],
      memo: new Map(),
      act,
      cy: hh,
      R,
      rp,
      rh: Math.max(rh, reach),
      bot: Math.min(0, swept.min.y),
      top: swept.max.y,
      // A hung piece is seen as from the line below it: a hood from a cook's eye, a heat lamp from the pass.
      elev: over ? THREE.MathUtils.degToRad(/HEATLAMP/.test(id) ? -30 : -22) : THREE.MathUtils.degToRad(THREE.MathUtils.lerp(32, 14, smooth(0.25, 1.6, (box.max.y - Math.min(0, box.min.y)) / (2 * rh)))),
      over,
      born: clock,
    };
    shown = s;
    buildCage(a, R, cx, cz, rp);
    cage.visible = state.cage;
    lightFor(s);
    frame(s, true);
    applySurface();
    operate(0);
    updateCard(s);
    overlay.dataset.ready = id;
    canvas!.setAttribute("aria-label", `3D turntable: ${a.label}, ${mm(w)} × ${mm(d)} × ${mm(h)} mm. Drag to turn, scroll or pinch to zoom.`);
  }

  /* --- Card ------------------------------------------------------------------------------------------------- */
  function updateCard(s: Shown): void {
    const a = s.arch;
    const [w, d, h] = a.ext;
    if (meta.name) meta.name.textContent = a.label;
    if (meta.id) meta.id.textContent = a.id;
    if (meta.dims) meta.dims.textContent = `${mm(w)} × ${mm(d)} × ${mm(h)} mm`;
    if (meta.tris) {
      meta.tris.textContent = `${Math.round(s.model.tris).toLocaleString("en-US")} tris`;
      meta.tris.title = `${s.meshes.length} merged meshes, built at close-up detail`;
    }
    if (meta.card) {
      meta.card.title = a.card ? `FrankenSim card ${a.card}` : "";
      if (a.card) meta.card.textContent = CARD_NAME[a.card] ?? a.card;
      else {
        // Outside the catalog: name the finish the model is mostly made of.
        const area = new Map<string, number>();
        const rev = new Map<THREE.Material, string>(Object.entries(p).map(([k, v]) => [v as THREE.Material, k]));
        for (const m of s.meshes) {
          const k = rev.get(m.material as THREE.Material);
          if (k && FINISH[k]) area.set(k, (area.get(k) ?? 0) + (m.geometry.index?.count ?? 0));
        }
        const top = [...area].sort((x, y) => y[1] - x[1])[0]?.[0];
        meta.card.textContent = top ? FINISH[top]! : "—";
      }
    }
    const jn = [...s.model.joints.keys()];
    if (meta.sockets) {
      meta.sockets.textContent = `${a.sockets.length} socket${a.sockets.length === 1 ? "" : "s"} · ${jn.length} joint${jn.length === 1 ? "" : "s"}`;
      meta.sockets.title = [jn.length ? `Joints: ${jointSummary(jn)}` : "No moving joints", a.sockets.length ? `Sockets: ${a.sockets.map((x) => x[0]).join(", ")}` : "No sockets"].join("\n");
    }
    if (meta.note) {
      const flag = FLAGSHIP_EQUIPMENT.find((f) => f.id === a.id);
      // A steak or a loaf has no utilities to speak of; a tool or a table says it has none.
      const util = /^ARCH-FOOD/.test(a.id) ? "" : utilities(a.util);
      meta.note.textContent = [a.name || flag?.label || a.label, util].filter(Boolean).join(" · ");
      meta.note.title = jn.length ? `Moving parts: ${jointSummary(jn)}` : "";
    }
  }

  /* --- Surface modes ---------------------------------------------------------------------------------------- */
  function restoreMaterials(s: Shown): void {
    for (const [m, mat] of s.orig) m.material = mat;
    s.orig.clear();
    for (const o of s.hidden) o.visible = true;
    s.hidden = [];
    if (s.edges) for (const e of s.edges) e.visible = false;
  }

  function zonesUniforms(s: Shown): void {
    const z = s.arch.zones.slice(0, N_ZONES);
    const lam = THREE.MathUtils.clamp(0.045 * Math.max(...s.arch.ext), 0.006, 0.09);
    for (const mat of inspectMats) {
      const u = mat.uniforms;
      const za = u.uZA!.value as THREE.Vector4[];
      const zb = u.uZB!.value as THREE.Vector4[];
      const zt = u.uZT!.value as number[];
      z.forEach(([k, a0, b0, c0, d0, e0, f0, t], i) => {
        if (k === "b") {
          za[i]!.set(a0, b0, c0, 0);
          zb[i]!.set(d0, e0, f0, 0);
        } else {
          za[i]!.set(a0, b0, c0, 1);
          zb[i]!.set(d0, 0, 0, 0);
        }
        zt[i] = t;
      });
      u.uZN!.value = z.length;
      u.uLambda!.value = lam;
      u.uThermal!.value = mat === ghostBody || mat === ghostJoint ? (state.thermal ? 1 : 0) : 1;
    }
  }

  /** Dress the shown piece for the current view, then compile any shaders that brings before the next draw. */
  function applySurface(): void {
    applySurfaceNow();
    void gate.refresh(scene);
  }
  function applySurfaceNow(): void {
    const s = shown;
    if (!s) return;
    restoreMaterials(s);
    renderPanel();
    if (!state.thermal && !state.wire) return;
    zonesUniforms(s);
    // Flames, steam and other sprites have no surface to paint.
    s.model.group.traverse((o) => {
      if (((o as THREE.InstancedMesh).isInstancedMesh || (o as THREE.Sprite).isSprite || (o as THREE.Points).isPoints) && o.visible) {
        o.visible = false;
        s.hidden.push(o);
      }
    });
    for (const m of s.meshes) {
      s.orig.set(m, m.material);
      const src = (Array.isArray(m.material) ? m.material[0] : m.material) as THREE.Material;
      if (state.wire) m.material = s.jointMeshes.has(m) ? ghostJoint : ghostBody;
      else m.material = src.side === THREE.DoubleSide || src.transparent ? thermalMats.double : thermalMats.front;
    }
    if (state.wire) {
      if (!s.edges) {
        s.edges = [];
        for (const m of s.meshes) {
          const e = new THREE.LineSegments(new THREE.EdgesGeometry(m.geometry, 32), s.jointMeshes.has(m) ? edgeJoint : edgeBody);
          e.raycast = () => undefined;
          m.add(e);
          s.edges.push(e);
        }
      }
      for (const e of s.edges) e.visible = true;
    }
  }

  function renderPanel(): void {
    const s = shown;
    if (!s || (!state.thermal && !state.wire)) {
      panel.classList.remove("on");
      return;
    }
    const a = s.arch;
    const parts: string[] = [];
    if (state.thermal) {
      const ticks = [-18, 0, 22, 100, 200, 340].map((t) => `<span style="left:${(rampPos(t) * 100).toFixed(1)}%">${t}</span>`).join("");
      const marks = [...new Set(a.zones.map((z) => z[7]))].map((t) => `<b style="left:${(rampPos(t) * 100).toFixed(1)}%"></b>`).join("");
      const seen = new Set<string>();
      const rows = a.zones
        .filter((z) => (seen.has(z[8]) ? false : (seen.add(z[8]), true)))
        .map((z) => `<div class="eqx-zone"><span>${z[8].replace(/_/g, " ")}</span><span>${z[7]} °C</span></div>`)
        .join("");
      const slots = a.slots.length ? `<div class="eqx-fine">Thermal slots: ${a.slots.map(([sid, ch]) => `${sid.replace(/_/g, " ")} <em>(${ch})</em>`).join(", ")}</div>` : "";
      parts.push(
        `<h4>Thermal preview · °C</h4><div class="eqx-bar" style="background:${rampCss()}">${ticks}${marks}</div>` +
          (rows || `<div class="eqx-zone"><span>No heated or chilled zone</span><span>${AMBIENT_C} °C</span></div>`) +
          slots +
          `<div class="eqx-fine eqx-honest">Nominal set points on the archetype's active zones, with an illustrative falloff; not a solved field.</div>`,
      );
    }
    if (state.wire) {
      const j = s.model.joints.size;
      parts.push(
        `<h4${state.thermal ? ' style="margin-top:9px"' : ""}>Wireframe X-ray</h4><span class="eqx-key" style="--k:#f59e0b">body</span><span class="eqx-key" style="--k:#22d3ee">moving joints (${j})</span><div class="eqx-fine">${Math.round(s.model.tris).toLocaleString("en-US")} triangles in ${s.meshes.length} merged meshes; feature edges over 32°.</div>`,
      );
    }
    panel.innerHTML = parts.join("");
    panel.classList.add("on");
  }

  /* --- Operation ------------------------------------------------------------------------------------------- */
  function operate(dt: number): void {
    const s = shown;
    if (!s) return;
    const t = clock - s.born;
    const h = s.model.hooks;
    if (debugHold) {
      h.tick?.(clock, dt, 0);
      return;
    }
    const rig: Rig = {
      t,
      on: state.animate,
      model: s.model,
      set: (name, v, eps = 0.003) => {
        const last = s.memo.get(name);
        if (last !== undefined && Math.abs(last - v) < eps) return;
        s.memo.set(name, v);
        h.set?.(name, v);
      },
      open: (v) => {
        const last = s.memo.get("@open");
        if (last !== undefined && Math.abs(last - v) < 0.002) return;
        s.memo.set("@open", v);
        h.open?.(v);
      },
      flows: (f) => {
        const n = h.burners ?? 0;
        if (!h.flames) return;
        for (let i = 0; i < n; i++) {
          const v = f(i);
          const k = `@flow${i}`;
          const last = s.memo.get(k);
          if (last !== undefined && Math.abs(last - v) < 0.006) continue;
          s.memo.set(k, v);
          h.flames.setFlow(i, v);
        }
        // Pilots and any burner past the valves (an oven's) burn low and steady.
        if (!s.memo.has("@pilot")) {
          s.memo.set("@pilot", 1);
          for (let i = n; i < h.flames.burners; i++) h.flames.setFlow(i, 0.18);
        }
      },
    };
    if (s.act) s.act(rig);
    else rig.open(state.animate ? cycle(t, 9) : 0);
    if (h.glows && !s.memo.has("heat") && s.act === null) rig.set("heat", 0.85);
    h.tick?.(clock, dt, state.animate ? 1 : 0);
  }

  /* --- Close looks: zoom onto surfaces, double-click to focus ------------------------------------------------ */
  const ray = new THREE.Raycaster();
  const ndc = new THREE.Vector2();
  function pick(clientX: number, clientY: number): THREE.Intersection | null {
    const s = shown;
    if (!s) return null;
    const r = canvas!.getBoundingClientRect();
    ndc.set(((clientX - r.left) / r.width) * 2 - 1, -((clientY - r.top) / r.height) * 2 + 1);
    ray.setFromCamera(ndc, camera);
    return ray.intersectObjects(visibleMeshes(s), false)[0] ?? null;
  }
  function visibleMeshes(s: Shown): THREE.Mesh[] {
    return s.meshes.filter((m) => {
      for (let o: THREE.Object3D | null = m; o; o = o.parent) if (!o.visible) return false;
      return true;
    });
  }
  /** Keep the orbit's pivot on the surface in front of the camera, so that zooming in stops at it rather than
   *  passing through. The pivot moves along the line of sight, so the view does not jump. */
  const fwd = new THREE.Vector3();
  function pivotToSurface(): void {
    const s = shown;
    if (!s) return;
    fwd.copy(controls.target).sub(camera.position);
    const dist = fwd.length();
    if (dist < 1e-6) return;
    fwd.divideScalar(dist);
    ray.set(camera.position, fwd);
    ray.far = dist;
    const hit = ray.intersectObjects(visibleMeshes(s), false)[0];
    ray.far = Infinity;
    if (hit && hit.distance < dist - controls.minDistance * 0.5) controls.target.copy(hit.point);
  }
  /**
   * Zooming to the cursor carries the camera along the cursor's ray, past the pivot's surface if another lies across
   * that ray: a walk-in's wall, a door. Stop it short of the first surface on its way from `from`, as far off it as
   * the orbit lets the camera come to its pivot.
   */
  const seg = new THREE.Vector3();
  function stopAtSurface(from: THREE.Vector3): void {
    const s = shown;
    if (!s) return;
    // Only on the way in: backing out through a thin double-sided panel behind the camera is harmless.
    if (camera.position.distanceToSquared(controls.target) >= from.distanceToSquared(controls.target)) return;
    seg.copy(camera.position).sub(from);
    const len = seg.length();
    if (len < 1e-7) return;
    seg.divideScalar(len);
    // The near plane follows the camera's distance to its pivot, so at this standoff the surface is never cut.
    const margin = controls.minDistance;
    ray.set(from, seg);
    ray.far = len + margin;
    const hit = ray.intersectObjects(visibleMeshes(s), false)[0];
    ray.far = Infinity;
    if (hit) camera.position.copy(from).addScaledVector(seg, Math.max(0, hit.distance - margin));
  }
  let pivotQueued = false;
  // Where the camera stood before this frame's zoom: OrbitControls moves it inside its own wheel and pinch handlers,
  // so this is taken in the capture phase, ahead of them.
  const camPrev = new THREE.Vector3();
  let camPrevSet = false;
  const notePrev = () => {
    if (camPrevSet) return;
    camPrev.copy(camera.position);
    camPrevSet = true;
  };
  stageEl.addEventListener("wheel", notePrev, { capture: true, passive: true });
  stageEl.addEventListener(
    "pointermove",
    () => {
      if (pinchPointers >= 2) notePrev();
    },
    { capture: true, passive: true },
  );
  canvas.addEventListener(
    "wheel",
    () => {
      if (!controls.enableZoom) return;
      lastInteract = clock;
      focusAnim = null;
      hint.style.opacity = "0";
      pivotQueued = true;
    },
    { passive: true },
  );
  let pinchPointers = 0;
  canvas.addEventListener("pointerdown", (e) => {
    if (e.pointerType === "touch") pinchPointers++;
  });
  const pointerGone = (e: PointerEvent) => {
    if (e.pointerType === "touch") pinchPointers = Math.max(0, pinchPointers - 1);
  };
  canvas.addEventListener("pointerup", pointerGone);
  canvas.addEventListener("pointercancel", pointerGone);
  canvas.addEventListener("pointermove", () => {
    if (pinchPointers >= 2) pivotQueued = true;
  });

  function focusAt(clientX: number, clientY: number): void {
    const hit = pick(clientX, clientY);
    lastInteract = clock;
    hint.style.opacity = "0";
    if (!hit) {
      // Double-click on empty studio: back to the whole model.
      focusAnim = { from: controls.target.clone(), to: home.clone(), d0: camera.position.distanceTo(controls.target), d1: fitDist, t: 0 };
      return;
    }
    const d0 = camera.position.distanceTo(controls.target);
    focusAnim = { from: controls.target.clone(), to: hit.point.clone(), d0, d1: Math.max(controls.minDistance * 2.5, Math.min(d0, hit.distance) * 0.34), t: 0 };
  }
  canvas.addEventListener("dblclick", (e) => focusAt(e.clientX, e.clientY));
  let lastTap = { t: 0, x: 0, y: 0 };
  canvas.addEventListener("pointerup", (e) => {
    if (e.pointerType !== "touch") return;
    const now = performance.now();
    if (now - lastTap.t < 320 && Math.hypot(e.clientX - lastTap.x, e.clientY - lastTap.y) < 30) {
      focusAt(e.clientX, e.clientY);
      lastTap = { t: 0, x: 0, y: 0 };
    } else lastTap = { t: now, x: e.clientX, y: e.clientY };
  });

  /* --- Chips, search, select ------------------------------------------------------------------------------ */
  function tabItems(tab: EquipmentTabId): Arch[] {
    const items = ARCH.filter((a) => tabOf(a) === tab);
    return [...items.filter((a) => FLAGSHIP_IDS.has(a.id)), ...items.filter((a) => !FLAGSHIP_IDS.has(a.id))];
  }
  function renderPills(query = ""): void {
    if (!pillGrid) return;
    const q = query.trim().toLowerCase();
    pillGrid.innerHTML = "";
    const items = q
      ? ARCH.filter((a) => [a.id, a.label, a.name, a.card, CARD_NAME[a.card] ?? "", a.util].some((f) => f.toLowerCase().includes(q))).slice(0, 16)
      : tabItems(activeTab);
    if (items.length === 0) {
      const empty = document.createElement("span");
      empty.className = "micro";
      empty.textContent = `Nothing in the ${ARCH.length} models matches “${query.trim()}”`;
      pillGrid.appendChild(empty);
      return;
    }
    for (const a of items) {
      const b = document.createElement("button");
      b.type = "button";
      b.className = "chip-btn";
      b.dataset.assetId = a.id;
      b.textContent = a.label;
      b.title = `${a.name || a.label} (${a.id})`;
      b.setAttribute("aria-pressed", String(a.id === activeId));
      b.addEventListener("click", () => choose(a.id));
      pillGrid.appendChild(b);
    }
  }
  function fillSelect(): void {
    if (!select) return;
    select.innerHTML = "";
    const names: Record<EquipmentTabId, string> = { hot_line: "Hot line", ovens_fryers: "Ovens & fryers", cold_chain: "Cold chain", ventilation: "Ventilation, sanitation & prep", cookware_food: "Cookware, tools & food" };
    for (const tab of Object.keys(names) as EquipmentTabId[]) {
      const g = document.createElement("optgroup");
      g.label = names[tab];
      for (const a of tabItems(tab)) {
        const o = document.createElement("option");
        o.value = a.id;
        o.textContent = `${a.label} · ${a.id}`;
        g.appendChild(o);
      }
      select.appendChild(g);
    }
    select.value = activeId;
  }
  function syncChoice(): void {
    if (pillGrid) for (const b of Array.from(pillGrid.querySelectorAll<HTMLButtonElement>("button"))) b.setAttribute("aria-pressed", String(b.dataset.assetId === activeId));
    if (select && select.value !== activeId) select.value = activeId;
    for (const tb of tabBtns) tb.setAttribute("aria-pressed", String(tb.dataset.eqTab === activeTab));
  }
  function choose(id: string): void {
    const a = BY_ID.get(id);
    if (!a) return;
    activeId = id;
    const tab = tabOf(a);
    const tabChanged = tab !== activeTab;
    activeTab = tab;
    if (tabChanged && !(searchInput && searchInput.value.trim())) renderPills("");
    syncChoice();
    // Let the pressed chip paint before the build.
    requestAnimationFrame(() => {
      if (activeId === id) show(id);
    });
  }

  for (const tb of tabBtns) {
    tb.addEventListener("click", () => {
      activeTab = (tb.dataset.eqTab as EquipmentTabId) ?? "hot_line";
      if (searchInput) searchInput.value = "";
      renderPills("");
      const first = tabItems(activeTab)[0];
      if (first) choose(first.id);
      else syncChoice();
    });
  }
  searchInput?.addEventListener("input", () => renderPills(searchInput.value));
  searchInput?.addEventListener("keydown", (e) => {
    if (e.key !== "Enter") return;
    const first = pillGrid?.querySelector<HTMLButtonElement>("button[data-asset-id]");
    if (first?.dataset.assetId) choose(first.dataset.assetId);
  });
  select?.addEventListener("change", () => {
    if (select.value) choose(select.value);
  });

  btn.joints?.addEventListener("click", () => {
    state.animate = !state.animate;
    syncButtons();
  });
  btn.thermal?.addEventListener("click", () => {
    state.thermal = !state.thermal;
    syncButtons();
    applySurface();
  });
  btn.wire?.addEventListener("click", () => {
    state.wire = !state.wire;
    syncButtons();
    applySurface();
  });
  btn.cage?.addEventListener("click", () => {
    state.cage = !state.cage;
    syncButtons();
    cage.visible = state.cage;
  });
  btn.spin?.addEventListener("click", () => {
    state.spin = !state.spin;
    syncButtons();
    lastInteract = -1e9;
  });

  /* --- Per-frame overlay ------------------------------------------------------------------------------------ */
  const vA = new THREE.Vector3();
  const vB = new THREE.Vector3();
  const camLocal = new THREE.Vector3();
  const placed: [number, number, number][] = [];
  /** Put a tag at a point's place on screen: centred on it, or (a socket, of label width `lw`) with its dot on it. */
  function place(el: HTMLElement, world: THREE.Vector3, w: number, h: number, dot = false, lw = 0): boolean {
    vA.copy(world).project(camera);
    if (vA.z > 1 || vA.z < -1 || Math.abs(vA.x) > 1.3 || Math.abs(vA.y) > 1.3) {
      el.style.visibility = "hidden";
      return false;
    }
    const x = (vA.x * 0.5 + 0.5) * w;
    const y = (-vA.y * 0.5 + 0.5) * h;
    el.style.visibility = "visible";
    el.style.transform = `translate3d(${x.toFixed(1)}px, ${y.toFixed(1)}px, 0) ${dot ? "translate(-4.5px, -50%)" : "translate(-50%, -50%)"}`;
    if (dot) {
      // A label that would sit on another keeps only its dot (its name is in the card's list).
      const crowded = placed.some(([px, py, pw]) => Math.abs(py - y) < 15 && x - 6 < px + pw && px - 6 < x + lw);
      el.classList.toggle("dot", crowded);
      placed.push([x, y, crowded ? 12 : lw]);
    }
    return true;
  }
  /**
   * A dimension tag: kept wholly inside the stage, and stepped clear of the tags already placed, so that
   * W, D and H never sit on one another or run off an edge.
   */
  const dimPlaced: [number, number, number, number][] = [];
  function placeDim(el: HTMLElement, world: THREE.Vector3, w: number, h: number): void {
    vA.copy(world).project(camera);
    if (vA.z > 1 || vA.z < -1 || Math.abs(vA.x) > 1.6 || Math.abs(vA.y) > 1.6) {
      el.style.visibility = "hidden";
      return;
    }
    const tw = el.offsetWidth || 90;
    const th = el.offsetHeight || 22;
    const clampX = (x: number) => Math.min(w - tw / 2 - 6, Math.max(tw / 2 + 6, x));
    const clampY = (y: number) => Math.min(h - th / 2 - 6, Math.max(th / 2 + 6, y));
    let x = clampX((vA.x * 0.5 + 0.5) * w);
    let y = clampY((-vA.y * 0.5 + 0.5) * h);
    for (let tries = 0; tries < 3; tries++) {
      const hit = dimPlaced.find(([px, py, pw, ph]) => Math.abs(px - x) < (pw + tw) / 2 + 4 && Math.abs(py - y) < (ph + th) / 2 + 3);
      if (!hit) break;
      const down = y >= hit[1] ? 1 : -1;
      const ny = hit[1] + down * ((hit[3] + th) / 2 + 4);
      y = clampY(ny);
      if (y !== ny) y = clampY(hit[1] - down * ((hit[3] + th) / 2 + 4));
      x = clampX(x);
    }
    dimPlaced.push([x, y, tw, th]);
    el.style.visibility = "visible";
    el.style.transform = `translate3d(${x.toFixed(1)}px, ${y.toFixed(1)}px, 0) translate(-50%, -50%)`;
  }
  function updateOverlay(): void {
    const s = shown;
    const show3 = !!s && state.cage;
    for (const t of dimTags) t.style.display = show3 ? "" : "none";
    for (const t of sockTags) t.el.style.display = show3 ? "" : "none";
    if (!s || !show3) return;
    const w = canvas!.clientWidth;
    const h = canvas!.clientHeight;
    s.frame.updateMatrixWorld(true);
    camLocal.copy(camera.position);
    s.frame.worldToLocal(camLocal);
    const [aw, ad] = s.arch.ext;
    const relX = camLocal.x - aw / 2;
    const relZ = camLocal.z + ad / 2;
    dimPlaced.length = 0;
    for (const set of dims) {
      let best = set.variants[0]!;
      let score = -Infinity;
      for (const v of set.variants) {
        const sc = set.pick === "facing" ? v.side.x * relX + v.side.z * relZ : s.frame.localToWorld(vA.copy(v.side)).project(camera).x;
        if (sc > score) {
          score = sc;
          best = v;
        }
      }
      for (const v of set.variants) v.obj.visible = v === best;
      placeDim(set.tag, s.frame.localToWorld(vB.copy(best.label)), w, h);
    }
    placed.length = 0;
    for (const t of sockTags) {
      if (!t.w && t.el.offsetWidth) t.w = t.el.offsetWidth;
      place(t.el, s.frame.localToWorld(vB.copy(t.at)), w, h, true, t.w || 120);
    }
  }

  /* --- Loop --------------------------------------------------------------------------------------------- */
  resize();
  const ro = new ResizeObserver(() => resize());
  ro.observe(stageEl);
  const firstBtn = Object.values(btn).find(Boolean);
  if (firstBtn?.parentElement) ro.observe(firstBtn.parentElement);
  observeInView(stageEl, (v) => {
    inView = v;
  });

  fillSelect();
  renderPills("");
  syncChoice();
  show(activeId);
  if (import.meta.env.DEV) {
    // A handle for the page checks (memory over a long browse, lighting A/B); absent from production builds.
    (window as unknown as Record<string, unknown>).__eq = {
      renderer,
      scene,
      camera,
      controls,
      key,
      turntable,
      platter,
      floor,
      pool,
      aoPad,
      choose,
      show,
      shown: () => shown,
      fit: () => ({ dist: fitDist, home: home.toArray() }),
      hold: (v: boolean) => (debugHold = v),
    };
  }

  const tmpM = new THREE.Matrix4();
  let memAt = 0;
  let last = performance.now();
  function tick(now: number): void {
    requestAnimationFrame(tick);
    if (!inView || document.hidden) {
      last = now;
      return;
    }
    const dt = Math.min(0.05, Math.max(0, (now - last) / 1000));
    last = now;
    clock += dt;
    const s = shown;
    if (!s) return;

    // Intro: ease in from a little further out.
    if (introT < 1) {
      const k0 = introT;
      introT = Math.min(1, introT + dt / 0.9);
      const e = (x: number) => 1 - Math.pow(1 - x, 3);
      const f0 = 1.18 - 0.18 * e(k0);
      const f1 = 1.18 - 0.18 * e(introT);
      camera.position.sub(controls.target).multiplyScalar(f1 / f0).add(controls.target);
    }
    // Focus: glide the pivot to the picked point and close in.
    if (focusAnim) {
      focusAnim.t = Math.min(1, focusAnim.t + dt / 0.65);
      const e = focusAnim.t * focusAnim.t * (3 - 2 * focusAnim.t);
      const dir = vB.copy(camera.position).sub(controls.target).normalize();
      controls.target.lerpVectors(focusAnim.from, focusAnim.to, e);
      const dist = THREE.MathUtils.lerp(focusAnim.d0, focusAnim.d1, e);
      camera.position.copy(controls.target).addScaledVector(dir, dist);
      if (focusAnim.t >= 1) focusAnim = null;
      lastInteract = clock;
    }
    const camDist = camera.position.distanceTo(controls.target);
    // Zoomed back out past the framing: drift the pivot home again.
    if (!dragging && !focusAnim && camDist > fitDist * 0.92 && controls.target.distanceToSquared(home) > 1e-10) {
      controls.target.lerp(home, 1 - Math.exp(-dt * 2.5));
    }
    const close = camDist < fitDist * 0.8 || controls.target.distanceTo(home) > s.R * 0.25;
    if (state.spin && !dragging && !focusAnim && !close && clock - lastInteract > 2.5) turntable.rotation.y += dt * 0.24;

    operate(dt);
    controls.update();
    if (pivotQueued) {
      pivotQueued = false;
      if (camPrevSet) stopAtSurface(camPrev);
      pivotToSurface();
    }
    camPrevSet = false;

    // Clip planes that follow the camera, tight enough to keep depth precise close up.
    const toCentre = camera.position.distanceTo(vA.set(0, s.cy, 0));
    const reach = Math.max(s.R, s.rp, s.rh);
    const near = THREE.MathUtils.clamp(Math.min(camDist * 0.2, toCentre - reach * 1.05), Math.max(0.0008, controls.minDistance * 0.08), 50);
    const far = toCentre + reach * 4 + 1;
    if (Math.abs(camera.near - near) > near * 0.05 || Math.abs(camera.far - far) > far * 0.05) {
      camera.near = near;
      camera.far = far;
      camera.updateProjectionMatrix();
    }

    if (state.thermal || state.wire) {
      s.frame.updateMatrixWorld(true);
      tmpM.copy(s.frame.matrixWorld).invert();
      for (const m of inspectMats) (m.uniforms.uRootInv!.value as THREE.Matrix4).copy(tmpM);
    }
    // While a new piece's shaders compile the canvas keeps the last frame, and its labels stay with it.
    if (!gate.ready()) return;
    updateOverlay();
    renderer.render(scene, camera);
    if (clock - memAt > 1) {
      memAt = clock;
      overlay.dataset.mem = `${renderer.info.memory.geometries} geometries, ${renderer.info.memory.textures} textures, ${renderer.info.programs?.length ?? 0} programs`;
    }
  }
  requestAnimationFrame(tick);
}
