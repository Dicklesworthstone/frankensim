import XCTest

final class FrankenSimAppearanceUITests: XCTestCase {
    override func setUpWithError() throws {
        continueAfterFailure = false
    }

    func testAppearanceTogglePersistsLightModeAcrossLaunches() throws {
        let app = XCUIApplication()
        app.launch()

        let toggle = app.buttons["appearance-toggle"]
        XCTAssertTrue(toggle.waitForExistence(timeout: 12))
        XCTAssertTrue(
            ["Switch to light mode", "Switch to dark mode"].contains(toggle.label),
            "Appearance control exposed an unexpected state: \(toggle.label)"
        )

        if toggle.label == "Switch to dark mode" {
            toggle.tap()
            XCTAssertEqual(toggle.label, "Switch to light mode")
        }

        toggle.tap()
        XCTAssertEqual(toggle.label, "Switch to dark mode")
        keepScreenshot(of: app, named: "Remembered light appearance")

        app.terminate()
        app.launch()

        let relaunchedToggle = app.buttons["appearance-toggle"]
        XCTAssertTrue(relaunchedToggle.waitForExistence(timeout: 12))
        XCTAssertEqual(relaunchedToggle.label, "Switch to dark mode")
    }

    func testPhoneStudioExposesAnExplicitRunAndTheFullCatalog() throws {
        let app = XCUIApplication()
        app.launchEnvironment["FSIM_INITIAL_EXPERIMENT"] = "13"
        app.launch()

        let chooser = app.buttons["compact-kernel-chooser"]
        XCTAssertTrue(chooser.waitForExistence(timeout: 12))
        XCTAssertTrue(app.otherElements["Simulation canvas ready"].exists)
        XCTAssertFalse(app.otherElements["Rendered native simulation result"].exists)

        let run = app.buttons["compact-run-button"]
        XCTAssertTrue(run.waitForExistence(timeout: 5))
        XCTAssertTrue(run.isHittable)
        XCTAssertEqual(run.label, "Run Lorenz attractor")
        XCTAssertEqual(run.value as? String, "Ready")
        run.tap()
        XCTAssertTrue(app.otherElements["Rendered native simulation result"].waitForExistence(timeout: 12))
        XCTAssertEqual(run.label, "Run Lorenz attractor")
        XCTAssertEqual(run.value as? String, "Completed")
        XCTAssertTrue(run.isHittable)
        keepScreenshot(of: app, named: "Phone studio after explicit run")

        chooser.tap()
        XCTAssertTrue(app.otherElements["compact-simulation-catalog"].waitForExistence(timeout: 5))
        XCTAssertTrue(app.textFields["catalog-search-field"].exists)
        XCTAssertEqual(app.staticTexts["catalog-result-count"].label, "44 of 44 on-device kernels")
        XCTAssertTrue(app.buttons["catalog-scope-all"].exists)
        XCTAssertTrue(app.buttons["catalog-scope-Foundations"].exists)
        XCTAssertTrue(app.buttons["catalog-scope-Frontier"].exists)
        XCTAssertTrue(app.buttons["catalog-scope-Deep Kernel"].exists)
        XCTAssertTrue(app.buttons["catalog-scope-Campaigns"].exists)
        XCTAssertTrue(app.buttons["catalog-scope-Flagships"].exists)
        XCTAssertTrue(app.buttons["catalog-theory-atlas"].exists)
        XCTAssertTrue(app.staticTexts["Heat diffusion"].exists)
        keepScreenshot(of: app, named: "Phone simulation catalog")

        app.buttons["catalog-theory-atlas"].tap()
        XCTAssertTrue(app.navigationBars["Theory Atlas"].waitForExistence(timeout: 5))
    }

    func testPhoneCatalogSearchAndTierScopesRevealTheFullInventory() throws {
        let app = XCUIApplication()
        app.launchEnvironment["FSIM_SHOW_CATALOG"] = "1"
        app.launch()

        let catalog = app.otherElements["compact-simulation-catalog"]
        XCTAssertTrue(catalog.waitForExistence(timeout: 12))

        let search = app.textFields["catalog-search-field"]
        XCTAssertTrue(search.exists)
        search.tap()
        search.typeText("reed")
        XCTAssertTrue(app.staticTexts["Reed bore"].waitForExistence(timeout: 5))
        XCTAssertEqual(app.staticTexts["catalog-result-count"].label, "1 of 44 on-device kernels")

        app.buttons["catalog-clear-search"].tap()
        let scopes = app.scrollViews["catalog-tier-scopes"]
        let flagships = app.buttons["catalog-scope-Flagships"]
        if !flagships.isHittable { scopes.swipeLeft() }
        XCTAssertTrue(flagships.isHittable)
        flagships.tap()
        XCTAssertEqual(app.staticTexts["catalog-result-count"].label, "3 of 44 on-device kernels")
        XCTAssertTrue(app.staticTexts["Ornithoid aircraft"].exists)
        XCTAssertTrue(app.staticTexts["Laminar vessel"].exists)
        XCTAssertTrue(app.staticTexts["Seismic frame"].exists)
        XCTAssertFalse(app.staticTexts["Heat diffusion"].exists)

        let all = app.buttons["catalog-scope-all"]
        if !all.isHittable { scopes.swipeRight() }
        XCTAssertTrue(all.isHittable)
        all.tap()
        XCTAssertEqual(app.staticTexts["catalog-result-count"].label, "44 of 44 on-device kernels")
    }

    func testPhoneAnimatedResultExposesPlaybackAndFrameInspection() throws {
        let app = XCUIApplication()
        app.launchEnvironment["FSIM_INITIAL_EXPERIMENT"] = "0"
        app.launchEnvironment["FSIM_INITIAL_QUALITY"] = "0.12"
        app.launch()

        let run = app.buttons["compact-run-button"]
        XCTAssertTrue(run.waitForExistence(timeout: 12))
        run.tap()

        let playback = app.otherElements["result-playback-controls"]
        XCTAssertTrue(playback.waitForExistence(timeout: 12))
        XCTAssertTrue(app.sliders["result-playback-slider"].exists)
        XCTAssertTrue(app.buttons["result-playback-restart"].exists)
        XCTAssertTrue(app.buttons["result-playback-speed"].exists)

        let toggle = app.buttons["result-playback-toggle"]
        XCTAssertTrue(toggle.exists)
        if toggle.isEnabled {
            XCTAssertEqual(toggle.label, "Pause result playback")
            toggle.tap()
            XCTAssertEqual(toggle.label, "Play result playback")
        }

        keepScreenshot(of: app, named: "Phone result playback controls")
    }

    private func keepScreenshot(of app: XCUIApplication, named name: String) {
        let attachment = XCTAttachment(screenshot: app.screenshot())
        attachment.name = name
        attachment.lifetime = .keepAlways
        add(attachment)
    }
}
