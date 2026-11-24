/** Every model family, registered. Import this (not the family files) to build fixtures. */
import "./hotline";
import "./cooking";
import "./ovens";
import "./cold";
import "./support";
import "./smallwares";

export { buildFixture, hasBuilder, type Ext, type Hooks, type Model } from "./index";
