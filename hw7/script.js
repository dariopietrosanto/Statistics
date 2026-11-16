/**
 * Homework 7: Statistical Independence & Bivariate Contingency Analysis
 * Author: Dario Pietrosanto (Matricola: 2055042)
 * Master's Degree Course in Applied Statistics & Cybersecurity
 *
 * Core Modules:
 * - Programmatic Contingency Cross-Tabulation (Observed, Expected, Marginals, Conditionals)
 * - Pearson Residuals & Haberman's Adjusted Standardized Residuals
 * - Inferential Hypothesis Testing: Pearson's Chi-Square, Likelihood Ratio Deviance (G^2),
 *   Fisher's Exact Test, and Measures of Association (Cramer's V, Tschuprow's T, Contingency Coeff C)
 * - Exact P-Value via Incomplete Gamma Function (Lanczos + Lentz Continued Fraction)
 * - Interactive Mosaic / Spineplot Canvas with Residual Coloring
 * - Theoretical Chi-Square Density & Critical Rejection Region Canvas
 * - Simpson's Paradox & Causal Confounder Demonstration Engine
 */

(function () {
    'use strict';

    // =========================================================================
    // 1. ADVANCED MATHEMATICAL FUNCTIONS & SPECIAL DISTRIBUTIONS
    // =========================================================================

    /**
     * Logarithm of the Gamma function via Lanczos approximation (g=7, n=9)
     * Precision: ~10^-15
     */
    function logGamma(z) {
        const p = [
            0.99999999999980993,
            676.5203681218851,
           -1259.1392167224028,
            771.32342877765313,
           -176.61502916214059,
            12.507343278686905,
           -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7
        ];
        if (z < 0.5) {
            return Math.log(Math.PI / Math.sin(Math.PI * z)) - logGamma(1.0 - z);
        }
        z -= 1.0;
        let x = p[0];
        for (let i = 1; i < p.length; i++) {
            x += p[i] / (z + i);
        }
        const t = z + 7.5;
        return 0.5 * Math.log(2.0 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
    }

    /**
     * Lower incomplete gamma function P(s, x) = gamma(s, x) / Gamma(s) via power series
     */
    function gammaLowerP(s, x) {
        if (x <= 0) return 0.0;
        let term = 1.0 / s;
        let sum = term;
        for (let k = 1; k < 120; k++) {
            term *= x / (s + k);
            sum += term;
            if (Math.abs(term / sum) < 1e-15) break;
        }
        return Math.exp(-x + s * Math.log(x) - logGamma(s)) * sum;
    }

    /**
     * Upper incomplete gamma function Q(s, x) = Gamma(s, x) / Gamma(s) via Lentz's continued fraction
     */
    function gammaUpperQ(s, x) {
        if (x <= 0) return 1.0;
        const TINY = 1e-30;
        let b = x + 1.0 - s;
        let c = 1.0 / TINY;
        let d = 1.0 / b;
        let f = d;
        for (let i = 1; i < 120; i++) {
            const a = -i * (i - s);
            b += 2.0;
            d = a * d + b;
            if (Math.abs(d) < TINY) d = TINY;
            c = b + a / c;
            if (Math.abs(c) < TINY) c = TINY;
            d = 1.0 / d;
            const delta = c * d;
            f *= delta;
            if (Math.abs(delta - 1.0) < 1e-15) break;
        }
        return Math.exp(-x + s * Math.log(x) - logGamma(s)) * f;
    }

    /**
     * Exact Chi-Square Cumulative P-value: P(Chi2_df >= chi2_obs)
     */
    function chiSquarePValue(chi2, df) {
        if (chi2 <= 0 || df <= 0) return 1.0;
        const s = df / 2.0;
        const x = chi2 / 2.0;
        if (x < s + 1.0) {
            return Math.max(0.0, Math.min(1.0, 1.0 - gammaLowerP(s, x)));
        } else {
            return Math.max(0.0, Math.min(1.0, gammaUpperQ(s, x)));
        }
    }

    /**
     * Probability Density Function (PDF) of Chi-Square distribution with df degrees of freedom
     */
    function chiSquarePDF(x, df) {
        if (x <= 0 || df <= 0) return 0.0;
        const s = df / 2.0;
        return Math.exp(-x / 2.0 + (s - 1.0) * Math.log(x) - s * Math.log(2.0) - logGamma(s));
    }

    /**
     * Chi-Square Critical Value (Inverse Survival Function) via Wilson-Hilferty + Newton-Raphson
     * Finds x such that P(Chi2_df >= x) = alpha
     */
    function chiSquareCritical(alpha, df) {
        if (df <= 0) return 0.0;
        // Normal quantile approximation for starting point
        let z = 1.6448536; // default 0.05
        if (Math.abs(alpha - 0.10) < 1e-4) z = 1.2815516;
        else if (Math.abs(alpha - 0.05) < 1e-4) z = 1.6448536;
        else if (Math.abs(alpha - 0.01) < 1e-4) z = 2.3263479;
        else if (Math.abs(alpha - 0.001) < 1e-4) z = 3.0902323;

        // Wilson-Hilferty transformation
        let x = df * Math.pow(1.0 - 2.0 / (9.0 * df) + z * Math.sqrt(2.0 / (9.0 * df)), 3.0);
        if (x <= 0) x = 0.5;

        // Newton-Raphson refinement
        for (let iter = 0; iter < 25; iter++) {
            const p = chiSquarePValue(x, df);
            const err = p - alpha;
            if (Math.abs(err) < 1e-13) break;
            const pdf = chiSquarePDF(x, df);
            if (pdf <= 0) break;
            const step = err / pdf;
            x = Math.max(0.0001, x + step);
        }
        return x;
    }

    /**
     * Fisher's Exact Test for 2x2 Contingency Tables (Hypergeometric two-tailed probability)
     */
    function fisherExact2x2(a, b, c, d) {
        function logFact(n) {
            return logGamma(n + 1.0);
        }
        function hyperProb(ca, cb, cc, cd) {
            const r1 = ca + cb, r2 = cc + cd, c1 = ca + cc, c2 = cb + cd, tot = ca + cb + cc + cd;
            return Math.exp(logFact(r1) + logFact(r2) + logFact(c1) + logFact(c2)
                - logFact(tot) - logFact(ca) - logFact(cb) - logFact(cc) - logFact(cd));
        }

        const r1 = a + b, r2 = c + d, c1 = a + c;
        const pObs = hyperProb(a, b, c, d);
        const minA = Math.max(0, c1 - r2);
        const maxA = Math.min(c1, r1);
        let pTotal = 0.0;

        for (let curA = minA; curA <= maxA; curA++) {
            const curB = r1 - curA;
            const curC = c1 - curA;
            const curD = r2 - curC;
            const pTable = hyperProb(curA, curB, curC, curD);
            if (pTable <= pObs + 1e-12) {
                pTotal += pTable;
            }
        }
        return Math.min(1.0, pTotal);
    }

    // =========================================================================
    // 2. PRE-LOADED DATASETS & PARSING UTILITIES
    // =========================================================================

    const PRESETS = {
        attack_layer: {
            title: "Cybersecurity Incidents: Attack Vector vs. Target Layer",
            description: "Cross-tabulation of 350 enterprise security incidents categorized by entry vector and compromised infrastructure layer.",
            rowVar: "Attack Vector (X)",
            colVar: "Target Infrastructure Layer (Y)",
            rowLabels: ["Phishing / Social Eng.", "SQLi / Web Exploit", "Zero-Day Remote Code", "DDoS Volumetric", "Insider Privilege Abuse"],
            colLabels: ["Web DMZ Tier", "Application API", "Database Cluster", "Cloud IAM / Secrets"],
            matrix: [
                [75, 20,  5, 10],  // Phishing -> heavily Web & Secrets
                [25, 60, 45,  5],  // SQLi -> App & DB
                [10, 35, 15, 25],  // Zero-Day -> distributed
                [40,  5,  0,  0],  // DDoS -> Web DMZ
                [ 5, 10, 30, 40]   // Insider -> DB & Cloud IAM
            ]
        },
        malware_persist: {
            title: "Endpoint Malware Family vs. Persistence Mechanism",
            description: "Forensic telemetry of 240 endpoints detailing threat family vs. MITRE ATT&CK persistence technique.",
            rowVar: "Malware Family (X)",
            colVar: "Persistence Mechanism (Y)",
            rowLabels: ["Ransomware", "Infostealer / Spyware", "Trojan-Dropper", "Rootkit / Worm"],
            colLabels: ["Registry Run Keys", "Scheduled Tasks", "WMI Event Subs.", "Service Creation"],
            matrix: [
                [35, 10,  5, 20],
                [40, 25, 15,  5],
                [15, 30, 20, 15],
                [ 5, 15, 25, 40]
            ]
        },
        mfa_breach: {
            title: "Authentication Policy vs. Account Takeover (ATO) Status",
            description: "Empirical study on 500 credential stuffing events evaluating authentication factor strength against breach outcomes.",
            rowVar: "MFA Enforcement Protocol (X)",
            colVar: "Incident Outcome (Y)",
            rowLabels: ["FIDO2 Hardware Key", "Time-Based OTP (App)", "Legacy SMS / Voice", "No MFA (Password Only)"],
            colLabels: ["Compromise Prevented", "Account Breached"],
            matrix: [
                [120,   2],  // FIDO2 virtually immune
                [105,  15],  // TOTP strong
                [ 60,  40],  // SMS vulnerable to SIM swap
                [ 15, 143]   // Passwords breached
            ]
        },
        independent_benchmark: {
            title: "Synthetic Benchmark: Strictly Independent Bivariate Distribution",
            description: "Theoretical sample where joint probabilities strictly equal product of marginals P(X, Y) = P(X)P(Y). Demonstrates chi-square null hypothesis.",
            rowVar: "Synthetic Factor A",
            colVar: "Synthetic Factor B",
            rowLabels: ["Class A1", "Class A2", "Class A3"],
            colLabels: ["Class B1", "Class B2", "Class B3"],
            matrix: [
                [30, 50, 20],  // Row 1 sum = 100
                [45, 75, 30],  // Row 2 sum = 150 (ratio 1.5x)
                [15, 25, 10]   // Row 3 sum = 50  (ratio 0.5x)
            ]
        },
        simpson_paradox: {
            title: "Simpson's Paradox in Security Patching (Aggregate View)",
            description: "Aggregate comparison of patch deployment strategy vs. breach outcome. Appears Automated Patching has higher failure rate!",
            rowVar: "Patch Deployment Protocol (X)",
            colVar: "System Outcome (Y)",
            rowLabels: ["Automated Pipeline", "Manual Sysadmin"],
            colLabels: ["Secured", "Breached"],
            matrix: [
                [270, 130], // Automated: 130 / 400 = 32.5% breached
                [160,  40]  // Manual:     40 / 200 = 20.0% breached
            ],
            stratified: {
                confounderVar: "Environment Exposure (Z)",
                strata: [
                    {
                        name: "Subgroup 1: External Facing / High Exposure Tier",
                        matrix: [
                            [ 90, 110], // Automated: 110 / 200 = 55% breached
                            [ 10,  30]  // Manual:      30 / 40  = 75% breached (Automated is BETTER!)
                        ]
                    },
                    {
                        name: "Subgroup 2: Internal Isolated / Low Exposure Tier",
                        matrix: [
                            [180,  20], // Automated:  20 / 200 = 10% breached
                            [150,  10]  // Manual:      10 / 160 = 6.25% breached
                        ]
                    }
                ]
            }
        }
    };

    // =========================================================================
    // 3. APPLICATION STATE & DOM REFS
    // =========================================================================

    let state = {
        activePreset: "attack_layer",
        rowLabels: [...PRESETS.attack_layer.rowLabels],
        colLabels: [...PRESETS.attack_layer.colLabels],
        rowVar: PRESETS.attack_layer.rowVar,
        colVar: PRESETS.attack_layer.colVar,
        matrix: PRESETS.attack_layer.matrix.map(row => [...row]),
        alpha: 0.05,
        analysisResults: null,
        activeTab: "observed" // 'observed' | 'relative' | 'row_cond' | 'col_cond' | 'expected' | 'residuals'
    };

    let dom = {};

    function initDOM() {
        dom.presetSelect = document.getElementById('preset-select');
        dom.presetDescription = document.getElementById('preset-description');
        dom.alphaSelect = document.getElementById('alpha-select');
        dom.rowVarInput = document.getElementById('row-var-input');
        dom.colVarInput = document.getElementById('col-var-input');

        // Table controls
        dom.addRowBtn = document.getElementById('add-row-btn');
        dom.delRowBtn = document.getElementById('del-row-btn');
        dom.addColBtn = document.getElementById('add-col-btn');
        dom.delColBtn = document.getElementById('del-col-btn');
        dom.resetTableBtn = document.getElementById('reset-table-btn');

        // Contingency table containers
        dom.tableEditorContainer = document.getElementById('table-editor-container');
        dom.tableTabsContainer = document.getElementById('table-tabs-container');
        dom.activeTableDisplay = document.getElementById('active-table-display');

        // Results displays
        dom.nTotalVal = document.getElementById('n-total-val');
        dom.dfVal = document.getElementById('df-val');
        dom.chi2Val = document.getElementById('chi2-val');
        dom.pValDisplay = document.getElementById('p-val-display');
        dom.critValDisplay = document.getElementById('crit-val-display');
        dom.decisionBadge = document.getElementById('decision-badge');
        dom.g2Val = document.getElementById('g2-val');
        dom.cramerVVal = document.getElementById('cramer-v-val');
        dom.cramerEffectVal = document.getElementById('cramer-effect-val');
        dom.cochranWarning = document.getElementById('cochran-warning');
        dom.fisher2x2Card = document.getElementById('fisher-2x2-card');
        dom.fisherPVal = document.getElementById('fisher-p-val');

        // Canvases
        dom.canvasMosaic = document.getElementById('canvas-mosaic');
        dom.canvasChi2Dist = document.getElementById('canvas-chi2-dist');
        dom.canvasResidualHeatmap = document.getElementById('canvas-residual-heatmap');

        // CSV Parser
        dom.csvInputArea = document.getElementById('csv-input-area');
        dom.parseCsvBtn = document.getElementById('parse-csv-btn');

        // Synthetic Generator
        dom.synthNInput = document.getElementById('synth-n-input');
        dom.synthDepSlider = document.getElementById('synth-dep-slider');
        dom.synthDepVal = document.getElementById('synth-dep-val');
        dom.generateSynthBtn = document.getElementById('generate-synth-btn');

        // Simpson's Paradox panel
        dom.simpsonContainer = document.getElementById('simpson-container');
    }

    // =========================================================================
    // 4. STATISTICAL ANALYSIS ENGINE
    // =========================================================================

    /**
     * Compute full statistical breakdown of contingency table
     */
    function computeContingencyAnalysis(matrix, rowLabels, colLabels, alpha = 0.05) {
        const r = matrix.length;
        const c = matrix[0].length;
        const rowSums = new Array(r).fill(0);
        const colSums = new Array(c).fill(0);
        let N = 0;

        for (let i = 0; i < r; i++) {
            for (let j = 0; j < c; j++) {
                const val = Math.max(0, Number(matrix[i][j]) || 0);
                rowSums[i] += val;
                colSums[j] += val;
                N += val;
            }
        }

        if (N <= 0) {
            return null;
        }

        // Relative joint frequencies p_ij = n_ij / N
        const jointProbs = [];
        const rowCondProbs = [];
        const colCondProbs = [];
        const expected = [];
        const rawResiduals = [];
        const pearsonResiduals = [];
        const adjustedResiduals = [];

        let chi2 = 0.0;
        let G2 = 0.0; // Likelihood ratio deviance
        let cellsUnder5 = 0;
        let minExpected = Infinity;

        for (let i = 0; i < r; i++) {
            jointProbs[i] = [];
            rowCondProbs[i] = [];
            expected[i] = [];
            rawResiduals[i] = [];
            pearsonResiduals[i] = [];
            adjustedResiduals[i] = [];

            for (let j = 0; j < c; j++) {
                const o = Math.max(0, Number(matrix[i][j]) || 0);
                const e = (rowSums[i] * colSums[j]) / N;

                jointProbs[i][j] = o / N;
                rowCondProbs[i][j] = rowSums[i] > 0 ? o / rowSums[i] : 0.0;

                expected[i][j] = e;
                if (e < minExpected) minExpected = e;
                if (e < 5.0) cellsUnder5++;

                const diff = o - e;
                rawResiduals[i][j] = diff;

                const pRes = e > 0 ? diff / Math.sqrt(e) : 0.0;
                pearsonResiduals[i][j] = pRes;
                chi2 += pRes * pRes;

                if (o > 0 && e > 0) {
                    G2 += 2.0 * o * Math.log(o / e);
                }

                // Haberman's adjusted standardized residual: N(0, 1) under H0
                const pRow = rowSums[i] / N;
                const pCol = colSums[j] / N;
                const denom = Math.sqrt(e * (1.0 - pRow) * (1.0 - pCol));
                adjustedResiduals[i][j] = denom > 0 ? diff / denom : 0.0;
            }
        }

        // Column conditional probabilities
        for (let i = 0; i < r; i++) {
            colCondProbs[i] = [];
            for (let j = 0; j < c; j++) {
                const o = Math.max(0, Number(matrix[i][j]) || 0);
                colCondProbs[i][j] = colSums[j] > 0 ? o / colSums[j] : 0.0;
            }
        }

        const df = Math.max(1, (r - 1) * (c - 1));
        const pValue = chiSquarePValue(chi2, df);
        const criticalValue = chiSquareCritical(alpha, df);
        const isIndependent = pValue >= alpha;

        // Measures of Association
        const minDim = Math.min(r - 1, c - 1);
        const cramersV = minDim > 0 ? Math.sqrt(chi2 / (N * minDim)) : 0.0;
        const tschuprowT = (r > 1 && c > 1) ? Math.sqrt(chi2 / (N * Math.sqrt((r - 1) * (c - 1)))) : 0.0;
        const contingencyC = Math.sqrt(chi2 / (chi2 + N));
        const mutualInfoNats = G2 / (2.0 * N);
        const mutualInfoBits = mutualInfoNats / Math.LN2;

        // Yates' correction and Fisher's exact test if 2x2 table
        let yatesChi2 = null;
        let yatesPVal = null;
        let fisherPVal = null;

        if (r === 2 && c === 2) {
            const a = matrix[0][0], b = matrix[0][1], cCell = matrix[1][0], d = matrix[1][1];
            const numYates = Math.max(0, Math.abs(a * d - b * cCell) - N / 2.0);
            const denYates = rowSums[0] * rowSums[1] * colSums[0] * colSums[1];
            if (denYates > 0) {
                yatesChi2 = (N * numYates * numYates) / denYates;
                yatesPVal = chiSquarePValue(yatesChi2, 1);
            }
            fisherPVal = fisherExact2x2(a, b, cCell, d);
        }

        return {
            r, c, N, df, alpha,
            rowSums, colSums,
            rowLabels, colLabels,
            matrix,
            jointProbs, rowCondProbs, colCondProbs,
            expected, rawResiduals, pearsonResiduals, adjustedResiduals,
            chi2, pValue, criticalValue, isIndependent,
            G2, mutualInfoNats, mutualInfoBits,
            cramersV, tschuprowT, contingencyC,
            cellsUnder5, minExpected, totalCells: r * c,
            percentUnder5: (cellsUnder5 / (r * c)) * 100.0,
            yatesChi2, yatesPVal, fisherPVal
        };
    }

    // =========================================================================
    // 5. RENDERING: TABLES & CONTINGENCY EDITORS
    // =========================================================================

    /**
     * Build the interactive live frequency editor table
     */
    function renderTableEditor() {
        if (!dom.tableEditorContainer) return;
        const r = state.matrix.length;
        const c = state.matrix[0].length;

        let html = '<table class="contingency-table-editor">';
        html += '<thead><tr>';
        html += `<th style="background:#34495e; color:#ecf0f1;">${escapeHTML(state.rowVar)} \\ ${escapeHTML(state.colVar)}</th>`;
        for (let j = 0; j < c; j++) {
            html += `<th><input type="text" class="col-header-input" data-col="${j}" value="${escapeHTML(state.colLabels[j])}"></th>`;
        }
        html += '<th style="background:#7f8c8d; color:white;">Row Sum (n_i·)</th>';
        html += '</tr></thead><tbody>';

        let grandTotal = 0;
        const colTotals = new Array(c).fill(0);

        for (let i = 0; i < r; i++) {
            let rowTotal = 0;
            html += '<tr>';
            html += `<th><input type="text" class="row-header-input" data-row="${i}" value="${escapeHTML(state.rowLabels[i])}"></th>`;
            for (let j = 0; j < c; j++) {
                const val = state.matrix[i][j];
                rowTotal += val;
                colTotals[j] += val;
                grandTotal += val;
                html += `<td><input type="number" min="0" step="1" class="cell-freq-input" data-row="${i}" data-col="${j}" value="${val}"></td>`;
            }
            html += `<td class="marginal-total-cell" style="font-weight:bold; background:#f4f6f7;">${rowTotal}</td>`;
            html += '</tr>';
        }

        // Bottom Column Marginals row
        html += '<tr style="font-weight:bold; background:#f4f6f7;">';
        html += '<th style="background:#7f8c8d; color:white;">Col Sum (n_·j)</th>';
        for (let j = 0; j < c; j++) {
            html += `<td class="marginal-total-cell">${colTotals[j]}</td>`;
        }
        html += `<td class="marginal-total-cell" style="background:#2c3e50; color:#fff; font-size:1.05rem;">N = ${grandTotal}</td>`;
        html += '</tr></tbody></table>';

        dom.tableEditorContainer.innerHTML = html;
        attachEditorListeners();
    }

    /**
     * Render the active statistical perspective table (Observed, Relative, Conditional, Expected, Residuals)
     */
    function renderActiveTableDisplay(res) {
        if (!dom.activeTableDisplay || !res) return;
        const r = res.r;
        const c = res.c;
        const tab = state.activeTab;

        let tableTitle = "";
        let tableSubtitle = "";

        if (tab === "observed") {
            tableTitle = "Observed Absolute Contingency Frequencies (n_ij)";
            tableSubtitle = "Raw counts cross-tabulated across categories with marginal row and column totals.";
        } else if (tab === "relative") {
            tableTitle = "Joint Empirical Probabilities p_ij = P(X = x_i, Y = y_j)";
            tableSubtitle = "Joint relative frequency matrix normalized by total sample size N = " + res.N + ".";
        } else if (tab === "row_cond") {
            tableTitle = "Row-Conditional Probabilities P(Y = y_j | X = x_i)";
            tableSubtitle = "Distribution of column variable Y conditioned upon row state X (each row sums to 100%).";
        } else if (tab === "col_cond") {
            tableTitle = "Column-Conditional Probabilities P(X = x_i | Y = y_j)";
            tableSubtitle = "Distribution of row variable X conditioned upon column state Y (each column sums to 100%).";
        } else if (tab === "expected") {
            tableTitle = "Theoretical Expected Frequencies Under Independence (e_ij = n_i· * n_·j / N)";
            tableSubtitle = "Frequencies that would arise if variables were statistically independent. Highlights Cochran criteria (< 5).";
        } else if (tab === "residuals") {
            tableTitle = "Standardized Haberman Residuals z_ij = (n_ij - e_ij) / SE";
            tableSubtitle = "Z-scores under H0 ~ N(0, 1). Absolute values > 1.96 (p < 0.05) or > 2.58 (p < 0.01) indicate localized structural attraction (+) or repulsion (-).";
        }

        let html = `<div style="margin-bottom:0.6rem;"><strong>${tableTitle}</strong><br><small style="color:#7f8c8d;">${tableSubtitle}</small></div>`;
        html += '<table class="contingency-display-table"><thead><tr>';
        html += `<th style="background:#2c3e50; color:#ffffff;">${escapeHTML(state.rowVar)} \\ ${escapeHTML(state.colVar)}</th>`;
        for (let j = 0; j < c; j++) {
            html += `<th>${escapeHTML(res.colLabels[j])}</th>`;
        }
        if (tab !== "residuals") {
            html += '<th style="background:#34495e; color:#ffffff;">Marginal Row Total</th>';
        }
        html += '</tr></thead><tbody>';

        for (let i = 0; i < r; i++) {
            html += '<tr>';
            html += `<th style="background:#f8f9fa;">${escapeHTML(res.rowLabels[i])}</th>`;

            for (let j = 0; j < c; j++) {
                let cellContent = "";
                let cellStyle = "";

                if (tab === "observed") {
                    cellContent = res.matrix[i][j];
                } else if (tab === "relative") {
                    cellContent = (res.jointProbs[i][j] * 100).toFixed(2) + "%";
                } else if (tab === "row_cond") {
                    const pct = (res.rowCondProbs[i][j] * 100).toFixed(2);
                    cellContent = `${pct}%<br><small style="color:#95a5a6;">(${res.matrix[i][j]}/${res.rowSums[i]})</small>`;
                } else if (tab === "col_cond") {
                    const pct = (res.colCondProbs[i][j] * 100).toFixed(2);
                    cellContent = `${pct}%<br><small style="color:#95a5a6;">(${res.matrix[i][j]}/${res.colSums[j]})</small>`;
                } else if (tab === "expected") {
                    const eVal = res.expected[i][j].toFixed(2);
                    if (res.expected[i][j] < 5.0) {
                        cellContent = `${eVal} <span style="color:#e74c3c; font-weight:bold; font-size:0.75rem;">(&lt;5)</span>`;
                        cellStyle = "background:#fdedec;";
                    } else {
                        cellContent = eVal;
                    }
                } else if (tab === "residuals") {
                    const z = res.adjustedResiduals[i][j];
                    const zStr = z >= 0 ? `+${z.toFixed(2)}` : z.toFixed(2);
                    if (Math.abs(z) >= 2.58) {
                        cellStyle = z > 0 ? "background:#d4efdf; font-weight:bold; color:#1e8449;" : "background:#fadbd8; font-weight:bold; color:#b03a2e;";
                    } else if (Math.abs(z) >= 1.96) {
                        cellStyle = z > 0 ? "background:#eafaf1; color:#27ae60;" : "background:#fdedec; color:#e74c3c;";
                    }
                    cellContent = `${zStr}<br><small style="color:#7f8c8d;">Diff: ${(res.matrix[i][j] - res.expected[i][j]).toFixed(1)}</small>`;
                }

                html += `<td style="${cellStyle}">${cellContent}</td>`;
            }

            // Marginal Row
            if (tab === "observed") {
                html += `<td style="font-weight:bold; background:#edf2f7;">${res.rowSums[i]}</td>`;
            } else if (tab === "relative") {
                html += `<td style="font-weight:bold; background:#edf2f7;">${((res.rowSums[i] / res.N) * 100).toFixed(2)}%</td>`;
            } else if (tab === "row_cond") {
                html += `<td style="font-weight:bold; background:#edf2f7;">100.00%</td>`;
            } else if (tab === "col_cond") {
                html += `<td style="font-weight:bold; background:#edf2f7;">-</td>`;
            } else if (tab === "expected") {
                html += `<td style="font-weight:bold; background:#edf2f7;">${res.rowSums[i]}</td>`;
            }

            html += '</tr>';
        }

        // Marginal Column Row
        if (tab !== "residuals") {
            html += '<tr style="font-weight:bold; background:#edf2f7;">';
            html += '<th style="background:#2c3e50; color:#ffffff;">Marginal Col Total</th>';
            for (let j = 0; j < c; j++) {
                if (tab === "observed" || tab === "expected") {
                    html += `<td>${res.colSums[j]}</td>`;
                } else if (tab === "relative") {
                    html += `<td>${((res.colSums[j] / res.N) * 100).toFixed(2)}%</td>`;
                } else if (tab === "row_cond") {
                    html += `<td>-</td>`;
                } else if (tab === "col_cond") {
                    html += `<td>100.00%</td>`;
                }
            }
            if (tab === "observed" || tab === "expected") {
                html += `<td style="background:#2c3e50; color:#ffffff;">N = ${res.N}</td>`;
            } else if (tab === "relative") {
                html += `<td style="background:#2c3e50; color:#ffffff;">100.00%</td>`;
            } else {
                html += `<td>-</td>`;
            }
            html += '</tr>';
        }

        html += '</tbody></table>';
        dom.activeTableDisplay.innerHTML = html;
    }

    /**
     * Attach event listeners to live editable table cells
     */
    function attachEditorListeners() {
        // Frequency cells
        document.querySelectorAll('.cell-freq-input').forEach(input => {
            input.addEventListener('change', (e) => {
                const rIdx = parseInt(e.target.dataset.row, 10);
                const cIdx = parseInt(e.target.dataset.col, 10);
                let val = parseInt(e.target.value, 10);
                if (isNaN(val) || val < 0) val = 0;
                state.matrix[rIdx][cIdx] = val;
                updateAll();
            });
        });

        // Row label headers
        document.querySelectorAll('.row-header-input').forEach(input => {
            input.addEventListener('change', (e) => {
                const rIdx = parseInt(e.target.dataset.row, 10);
                state.rowLabels[rIdx] = e.target.value.trim() || `Row ${rIdx + 1}`;
                updateAll();
            });
        });

        // Col label headers
        document.querySelectorAll('.col-header-input').forEach(input => {
            input.addEventListener('change', (e) => {
                const cIdx = parseInt(e.target.dataset.col, 10);
                state.colLabels[cIdx] = e.target.value.trim() || `Col ${cIdx + 1}`;
                updateAll();
            });
        });
    }

    // =========================================================================
    // 6. VISUALIZATION CANVASES
    // =========================================================================

    /**
     * Canvas 1: Mosaic Plot / Spineplot with Pearson Residual Color Coding
     */
    function renderMosaicCanvas(res) {
        const canvas = dom.canvasMosaic;
        if (!canvas || !res) return;
        const ctx = canvas.getContext('2d');
        const W = canvas.width, H = canvas.height;
        ctx.clearRect(0, 0, W, H);

        const padLeft = 60, padRight = 30, padTop = 30, padBottom = 60;
        const plotW = W - padLeft - padRight;
        const plotH = H - padTop - padBottom;

        const r = res.r;
        const c = res.c;
        const N = res.N;

        ctx.fillStyle = "#fafbfc";
        ctx.fillRect(padLeft, padTop, plotW, plotH);

        // Column widths proportional to marginal column totals colSums[j] / N
        let curX = padLeft;
        for (let j = 0; j < c; j++) {
            const colFrac = res.colSums[j] / N;
            const colWidth = colFrac * plotW;

            // Within this column, vertical splits proportional to column conditional prob colCondProbs[i][j]
            let curY = padTop + plotH;
            for (let i = r - 1; i >= 0; i--) {
                const condFrac = res.colCondProbs[i][j];
                const sliceH = condFrac * plotH;
                curY -= sliceH;

                // Color based on Pearson standardized residual
                const z = res.adjustedResiduals[i][j];
                let fillColor = "#ecf0f1";
                if (z > 2.58) fillColor = "rgba(39, 174, 96, 0.85)";
                else if (z > 1.96) fillColor = "rgba(46, 204, 113, 0.65)";
                else if (z > 1.0) fillColor = "rgba(46, 204, 113, 0.35)";
                else if (z < -2.58) fillColor = "rgba(192, 57, 43, 0.85)";
                else if (z < -1.96) fillColor = "rgba(231, 76, 60, 0.65)";
                else if (z < -1.0) fillColor = "rgba(231, 76, 60, 0.35)";
                else fillColor = "rgba(189, 195, 199, 0.45)";

                ctx.fillStyle = fillColor;
                ctx.fillRect(curX, curY, colWidth, sliceH);

                ctx.strokeStyle = "#ffffff";
                ctx.lineWidth = 1.5;
                ctx.strokeRect(curX, curY, colWidth, sliceH);

                // Cell label if slice is sufficiently large
                if (colWidth > 35 && sliceH > 18) {
                    ctx.fillStyle = Math.abs(z) > 1.96 ? "#ffffff" : "#2c3e50";
                    ctx.font = "bold 10px sans-serif";
                    ctx.textAlign = "center";
                    ctx.fillText(res.matrix[i][j], curX + colWidth / 2, curY + sliceH / 2 + 4);
                }
            }

            // Column label at bottom
            ctx.fillStyle = "#2c3e50";
            ctx.font = "10px sans-serif";
            ctx.textAlign = "center";
            const colLabelShort = res.colLabels[j].length > 14 ? res.colLabels[j].slice(0, 12) + ".." : res.colLabels[j];
            ctx.fillText(colLabelShort, curX + colWidth / 2, padTop + plotH + 16);
            ctx.fillStyle = "#7f8c8d";
            ctx.font = "9px monospace";
            ctx.fillText(`${(colFrac * 100).toFixed(1)}%`, curX + colWidth / 2, padTop + plotH + 28);

            curX += colWidth;
        }

        // Row labels on left axis
        ctx.fillStyle = "#2c3e50";
        ctx.font = "10px sans-serif";
        ctx.textAlign = "right";
        let leftY = padTop + plotH;
        for (let i = r - 1; i >= 0; i--) {
            const rowFrac = res.rowSums[i] / N;
            const rH = rowFrac * plotH;
            const rowLabelShort = res.rowLabels[i].length > 12 ? res.rowLabels[i].slice(0, 10) + ".." : res.rowLabels[i];
            ctx.fillText(rowLabelShort, padLeft - 8, leftY - rH / 2 + 3);
            leftY -= rH;
        }

        // Title and Legend
        ctx.fillStyle = "#2c3e50";
        ctx.font = "bold 11px sans-serif";
        ctx.textAlign = "left";
        ctx.fillText(`Mosaic Plot: Width = P(${res.colVar}), Height = P(${res.rowVar} | ${res.colVar})`, padLeft, padTop - 12);
    }

    /**
     * Canvas 2: Theoretical Chi-Square Density & Critical Region
     */
    function renderChi2DistCanvas(res) {
        const canvas = dom.canvasChi2Dist;
        if (!canvas || !res) return;
        const ctx = canvas.getContext('2d');
        const W = canvas.width, H = canvas.height;
        ctx.clearRect(0, 0, W, H);

        const padLeft = 45, padRight = 30, padTop = 30, padBottom = 35;
        const plotW = W - padLeft - padRight;
        const plotH = H - padTop - padBottom;

        const df = res.df;
        const chi2Obs = res.chi2;
        const chi2Crit = res.criticalValue;

        // X-range: from 0 up to max of (crit * 1.5, obs * 1.3, df * 3.5, 15)
        const xMax = Math.max(16.0, chi2Crit * 1.5, chi2Obs * 1.3, df * 3.0);
        const xMin = 0.0;

        // Peak y-value for scaling
        let peakY = 0.0;
        const modeX = Math.max(0.1, df - 2);
        peakY = Math.max(0.15, chiSquarePDF(modeX, df) * 1.25);
        if (df === 1) peakY = 0.8; // singularity at 0 capped

        const mapX = (x) => padLeft + ((x - xMin) / (xMax - xMin)) * plotW;
        const mapY = (y) => padTop + (1.0 - Math.min(1.0, y / peakY)) * plotH;

        ctx.fillStyle = "#fafbfc";
        ctx.fillRect(padLeft, padTop, plotW, plotH);

        // Grid
        ctx.strokeStyle = "#e1e4e8";
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let g = 0; g <= 4; g++) {
            const gy = padTop + (g / 4) * plotH;
            ctx.moveTo(padLeft, gy);
            ctx.lineTo(padLeft + plotW, gy);
        }
        ctx.stroke();

        // Shaded Critical Rejection Region: [chi2Crit, xMax]
        ctx.beginPath();
        ctx.moveTo(mapX(chi2Crit), mapY(0));
        const critSteps = 100;
        const stepSize = (xMax - chi2Crit) / critSteps;
        for (let k = 0; k <= critSteps; k++) {
            const curX = chi2Crit + k * stepSize;
            ctx.lineTo(mapX(curX), mapY(chiSquarePDF(curX, df)));
        }
        ctx.lineTo(mapX(xMax), mapY(0));
        ctx.closePath();
        ctx.fillStyle = "rgba(231, 76, 60, 0.35)";
        ctx.fill();

        // Draw Chi-Square PDF curve
        ctx.beginPath();
        ctx.strokeStyle = "#2c3e50";
        ctx.lineWidth = 2.5;
        const numPts = 300;
        for (let k = 0; k <= numPts; k++) {
            const curX = xMin + (k / numPts) * (xMax - xMin);
            const curY = chiSquarePDF(curX, df);
            const px = mapX(curX);
            const py = mapY(curY);
            if (k === 0) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
        }
        ctx.stroke();

        // Critical value vertical boundary line
        const cxX = mapX(chi2Crit);
        ctx.strokeStyle = "#e74c3c";
        ctx.lineWidth = 2;
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(cxX, padTop);
        ctx.lineTo(cxX, padTop + plotH);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = "#e74c3c";
        ctx.font = "bold 10px sans-serif";
        ctx.textAlign = "center";
        ctx.fillText(`χ²_crit = ${chi2Crit.toFixed(2)}`, cxX, padTop + 14);

        // Observed Statistic vertical line (in blue)
        const obsX = mapX(Math.min(xMax, chi2Obs));
        ctx.strokeStyle = "#2980b9";
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        ctx.moveTo(obsX, padTop + 20);
        ctx.lineTo(obsX, padTop + plotH);
        ctx.stroke();

        // Diamond marker at top of obs line
        ctx.fillStyle = "#2980b9";
        ctx.beginPath();
        ctx.arc(obsX, padTop + 20, 5, 0, 2 * Math.PI);
        ctx.fill();

        ctx.fillStyle = "#2980b9";
        ctx.font = "bold 10.5px sans-serif";
        ctx.textAlign = (obsX > padLeft + plotW - 70) ? "right" : "left";
        const textOffset = (obsX > padLeft + plotW - 70) ? -8 : 8;
        ctx.fillText(`χ²_obs = ${chi2Obs.toFixed(2)} (p = ${res.pValue < 0.0001 ? '<0.0001' : res.pValue.toFixed(4)})`, obsX + textOffset, padTop + 24);

        // Axes
        ctx.strokeStyle = "#2c3e50";
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(padLeft, padTop);
        ctx.lineTo(padLeft, padTop + plotH);
        ctx.lineTo(padLeft + plotW, padTop + plotH);
        ctx.stroke();

        // Axis label
        ctx.fillStyle = "#2c3e50";
        ctx.font = "bold 11px sans-serif";
        ctx.textAlign = "center";
        ctx.fillText(`Chi-Square Distribution (df = ${df}, α = ${res.alpha})`, padLeft + plotW / 2, padTop + plotH + 25);
    }

    /**
     * Canvas 3: Standardized Residuals Heatmap
     */
    function renderResidualHeatmap(res) {
        const canvas = dom.canvasResidualHeatmap;
        if (!canvas || !res) return;
        const ctx = canvas.getContext('2d');
        const W = canvas.width, H = canvas.height;
        ctx.clearRect(0, 0, W, H);

        const padLeft = 80, padRight = 30, padTop = 35, padBottom = 50;
        const plotW = W - padLeft - padRight;
        const plotH = H - padTop - padBottom;

        const r = res.r;
        const c = res.c;
        const cellW = plotW / c;
        const cellH = plotH / r;

        for (let i = 0; i < r; i++) {
            for (let j = 0; j < c; j++) {
                const z = res.adjustedResiduals[i][j];
                const x = padLeft + j * cellW;
                const y = padTop + i * cellH;

                let bg = "#f8f9fa";
                if (z > 2.58) bg = "#27ae60";
                else if (z > 1.96) bg = "#2ecc71";
                else if (z > 0.5) bg = "#abebc6";
                else if (z < -2.58) bg = "#c0392b";
                else if (z < -1.96) bg = "#e74c3c";
                else if (z < -0.5) bg = "#f5b7b1";

                ctx.fillStyle = bg;
                ctx.fillRect(x, y, cellW, cellH);

                ctx.strokeStyle = "#ffffff";
                ctx.lineWidth = 1;
                ctx.strokeRect(x, y, cellW, cellH);

                ctx.fillStyle = Math.abs(z) > 1.96 ? "#ffffff" : "#2c3e50";
                ctx.font = "bold 11px monospace";
                ctx.textAlign = "center";
                ctx.fillText((z >= 0 ? "+" : "") + z.toFixed(2), x + cellW / 2, y + cellH / 2 + 4);
            }

            // Row header
            ctx.fillStyle = "#2c3e50";
            ctx.font = "10px sans-serif";
            ctx.textAlign = "right";
            const rowLabelShort = res.rowLabels[i].length > 14 ? res.rowLabels[i].slice(0, 12) + ".." : res.rowLabels[i];
            ctx.fillText(rowLabelShort, padLeft - 8, padTop + i * cellH + cellH / 2 + 4);
        }

        // Col headers
        for (let j = 0; j < c; j++) {
            ctx.fillStyle = "#2c3e50";
            ctx.font = "10px sans-serif";
            ctx.textAlign = "center";
            const colLabelShort = res.colLabels[j].length > 12 ? res.colLabels[j].slice(0, 10) + ".." : res.colLabels[j];
            ctx.fillText(colLabelShort, padLeft + j * cellW + cellW / 2, padTop + plotH + 16);
        }

        ctx.fillStyle = "#2c3e50";
        ctx.font = "bold 11px sans-serif";
        ctx.textAlign = "left";
        ctx.fillText("Standardized Haberman Residuals Heatmap (Z-score attraction/repulsion)", padLeft, padTop - 12);
    }

    // =========================================================================
    // 7. SIMPSON'S PARADOX CAUSAL CONVOLUTION ENGINE
    // =========================================================================

    function renderSimpsonsParadox() {
        if (!dom.simpsonContainer) return;
        const sp = PRESETS.simpson_paradox;
        const aggRes = computeContingencyAnalysis(sp.matrix, sp.rowLabels, sp.colLabels);
        const sub1 = computeContingencyAnalysis(sp.stratified.strata[0].matrix, sp.rowLabels, sp.colLabels);
        const sub2 = computeContingencyAnalysis(sp.stratified.strata[1].matrix, sp.rowLabels, sp.colLabels);

        let html = `
            <div class="stat-callout" style="background:#fcf3cf; border-left:5px solid #f39c12; margin-top:1.5rem;">
                <h4 style="margin-top:0; color:#b7950b;">Interactive Epistemological Laboratory: Simpson's Paradox & The Causal Backdoor</h4>
                <p>
                    <strong>The Aggregation Illusion:</strong> Below, we evaluate an enterprise patch deployment dataset. 
                    In the <em>aggregate</em> table, Automated Patching appears to have a higher breach rate than Manual Patching. 
                    However, when conditioning on the confounding variable $Z$ (<strong>Environment Exposure Tier</strong>), Automated Patching is superior in every single individual environment!
                </p>
                <div class="layout-grid" style="grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap:1rem; margin-top:1rem;">
                    <div style="background:#ffffff; padding:1rem; border:1px solid #d5dbdb; border-radius:4px;">
                        <h5 style="margin-top:0; color:#c0392b;">1. Aggregate Data (Confounded)</h5>
                        <table class="contingency-display-table" style="font-size:0.85rem;">
                            <tr><th>Strategy</th><th>Secured</th><th>Breached</th><th>Breach Rate</th></tr>
                            <tr><td>Automated</td><td>270</td><td>130</td><td style="color:#e74c3c; font-weight:bold;">32.5%</td></tr>
                            <tr><td>Manual</td><td>160</td><td>40</td><td style="color:#27ae60; font-weight:bold;">20.0%</td></tr>
                        </table>
                        <p style="font-size:0.8rem; color:#7f8c8d; margin-top:0.5rem;">
                            <strong>Naive Inference:</strong> Manual seems 12.5% safer!<br>
                            Chi-square: χ² = ${aggRes.chi2.toFixed(2)}, p = ${aggRes.pValue.toExponential(2)} (Spurious Association).
                        </p>
                    </div>

                    <div style="background:#ffffff; padding:1rem; border:1px solid #d5dbdb; border-radius:4px;">
                        <h5 style="margin-top:0; color:#2980b9;">2. Subgroup A: High-Exposure DMZ (Z = 1)</h5>
                        <table class="contingency-display-table" style="font-size:0.85rem;">
                            <tr><th>Strategy</th><th>Secured</th><th>Breached</th><th>Breach Rate</th></tr>
                            <tr><td>Automated</td><td>90</td><td>110</td><td style="color:#27ae60; font-weight:bold;">55.0%</td></tr>
                            <tr><td>Manual</td><td>10</td><td>30</td><td style="color:#e74c3c; font-weight:bold;">75.0%</td></tr>
                        </table>
                        <p style="font-size:0.8rem; color:#27ae60; font-weight:bold; margin-top:0.5rem;">
                            Automated is 20% SAFER in High-Exposure DMZ!
                        </p>
                    </div>

                    <div style="background:#ffffff; padding:1rem; border:1px solid #d5dbdb; border-radius:4px;">
                        <h5 style="margin-top:0; color:#2980b9;">3. Subgroup B: Isolated Air-Gap (Z = 2)</h5>
                        <table class="contingency-display-table" style="font-size:0.85rem;">
                            <tr><th>Strategy</th><th>Secured</th><th>Breached</th><th>Breach Rate</th></tr>
                            <tr><td>Automated</td><td>180</td><td>20</td><td style="color:#27ae60; font-weight:bold;">10.0%</td></tr>
                            <tr><td>Manual</td><td>150</td><td>10</td><td style="color:#e74c3c; font-weight:bold;">6.25%</td></tr>
                        </table>
                        <p style="font-size:0.8rem; color:#7f8c8d; margin-top:0.5rem;">
                            Baseline risk is low; manual sysadmins selectively assigned only to low-risk nodes.
                        </p>
                    </div>
                </div>

                <div style="background:#ffffff; padding:1rem; border-radius:4px; margin-top:1rem; border:1px solid #e1e8ed;">
                    <strong>Causal Analysis via Pearl's Do-Calculus:</strong>
                    The confounding variable $Z$ (Exposure) acts as a common cause (fork: $X \leftarrow Z \rightarrow Y$). 
                    Because the Automated strategy was disproportionately assigned to the most hazardous internet-facing machines ($P(Z=1 \mid \text{Automated}) = 0.50$ vs $P(Z=1 \mid \text{Manual}) = 0.20$), the crude observational association $P(Y \mid X)$ contains backdoor confounding bias. 
                    Applying the <strong>Adjustment Formula</strong>:
                    $$P(Y = \text{Breached} \mid \text{do}(X = \text{Auto})) = \sum_z P(Y = \text{Breached} \mid X = \text{Auto}, Z = z) P(Z = z) \approx 25.0\%$$
                    $$P(Y = \text{Breached} \mid \text{do}(X = \text{Manual})) = \sum_z P(Y = \text{Breached} \mid X = \text{Manual}, Z = z) P(Z = z) \approx 30.8\%$$
                    Under true intervention ($\text{do}(X)$), Automated Patching reduces breach risk by over 5.8 percentage points.
                </div>
            </div>
        `;

        dom.simpsonContainer.innerHTML = html;
    }

    // =========================================================================
    // 8. MASTER UPDATE DISPATCHER
    // =========================================================================

    function updateAll() {
        const res = computeContingencyAnalysis(state.matrix, state.rowLabels, state.colLabels, state.alpha);
        state.analysisResults = res;

        if (!res) return;

        // Update Top Metrics
        dom.nTotalVal.textContent = res.N.toLocaleString();
        dom.dfVal.textContent = res.df;
        dom.chi2Val.textContent = res.chi2.toFixed(4);
        dom.critValDisplay.textContent = res.criticalValue.toFixed(4);

        if (res.pValue < 0.0001) {
            dom.pValDisplay.textContent = "< 0.0001";
        } else {
            dom.pValDisplay.textContent = res.pValue.toFixed(6);
        }

        // Decision Badge
        if (res.isIndependent) {
            dom.decisionBadge.innerHTML = `<span style="background:#27ae60; color:#fff; padding:4px 10px; border-radius:4px; font-weight:bold;">STATISTICALLY INDEPENDENT (Fail to Reject H0)</span>`;
        } else {
            dom.decisionBadge.innerHTML = `<span style="background:#c0392b; color:#fff; padding:4px 10px; border-radius:4px; font-weight:bold;">STATISTICALLY DEPENDENT (Reject H0 at α = ${state.alpha})</span>`;
        }

        // Additional stats
        dom.g2Val.textContent = `${res.G2.toFixed(3)} (${res.mutualInfoBits.toFixed(3)} bits)`;
        dom.cramerVVal.textContent = res.cramersV.toFixed(4);

        let effectStr = "Negligible";
        if (res.cramersV > 0.5) effectStr = "Very Strong Association";
        else if (res.cramersV > 0.3) effectStr = "Moderate-to-Strong";
        else if (res.cramersV > 0.1) effectStr = "Weak Association";
        dom.cramerEffectVal.textContent = effectStr;

        // Cochran Warning
        if (res.percentUnder5 > 20.0 || res.minExpected < 1.0) {
            dom.cochranWarning.style.display = "block";
            dom.cochranWarning.innerHTML = `<strong>Warning (Cochran's Criterion Violation):</strong> ${res.cellsUnder5} out of ${res.totalCells} cells (${res.percentUnder5.toFixed(1)}%) have expected frequency e_ij &lt; 5. Minimum e_ij = ${res.minExpected.toFixed(2)}. Chi-square asymptotic approximation may be unreliable; consider Exact tests or category pooling.`;
        } else {
            dom.cochranWarning.style.display = "none";
        }

        // 2x2 Fisher Card
        if (res.r === 2 && res.c === 2) {
            dom.fisher2x2Card.style.display = "block";
            dom.fisherPVal.textContent = res.fisherPVal < 0.0001 ? "< 0.0001" : res.fisherPVal.toFixed(6);
        } else {
            dom.fisher2x2Card.style.display = "none";
        }

        // Re-render subviews
        renderActiveTableDisplay(res);
        renderMosaicCanvas(res);
        renderChi2DistCanvas(res);
        renderResidualHeatmap(res);
    }

    // =========================================================================
    // 9. EVENT LISTENERS & INITIALIZATION
    // =========================================================================

    function attachEventListeners() {
        // Preset selector
        dom.presetSelect.addEventListener('change', (e) => {
            const key = e.target.value;
            state.activePreset = key;
            const p = PRESETS[key];
            state.rowLabels = [...p.rowLabels];
            state.colLabels = [...p.colLabels];
            state.rowVar = p.rowVar;
            state.colVar = p.colVar;
            state.matrix = p.matrix.map(row => [...row]);

            dom.presetDescription.innerHTML = `<strong>${escapeHTML(p.title)}:</strong> ${escapeHTML(p.description)}`;
            dom.rowVarInput.value = state.rowVar;
            dom.colVarInput.value = state.colVar;

            renderTableEditor();
            updateAll();
        });

        // Significance Alpha select
        dom.alphaSelect.addEventListener('change', (e) => {
            state.alpha = parseFloat(e.target.value);
            updateAll();
        });

        // Row/Col variable name inputs
        dom.rowVarInput.addEventListener('change', (e) => {
            state.rowVar = e.target.value.trim() || "Row Factor (X)";
            renderTableEditor();
            updateAll();
        });
        dom.colVarInput.addEventListener('change', (e) => {
            state.colVar = e.target.value.trim() || "Col Factor (Y)";
            renderTableEditor();
            updateAll();
        });

        // Add Row
        dom.addRowBtn.addEventListener('click', () => {
            const newRIdx = state.matrix.length + 1;
            const c = state.matrix[0].length;
            state.rowLabels.push(`Category X${newRIdx}`);
            state.matrix.push(new Array(c).fill(10));
            renderTableEditor();
            updateAll();
        });

        // Delete Row
        dom.delRowBtn.addEventListener('click', () => {
            if (state.matrix.length > 2) {
                state.matrix.pop();
                state.rowLabels.pop();
                renderTableEditor();
                updateAll();
            } else {
                alert("Table must have at least 2 rows for contingency analysis.");
            }
        });

        // Add Column
        dom.addColBtn.addEventListener('click', () => {
            const newCIdx = state.matrix[0].length + 1;
            state.colLabels.push(`Category Y${newCIdx}`);
            for (let i = 0; i < state.matrix.length; i++) {
                state.matrix[i].push(10);
            }
            renderTableEditor();
            updateAll();
        });

        // Delete Column
        dom.delColBtn.addEventListener('click', () => {
            if (state.matrix[0].length > 2) {
                state.colLabels.pop();
                for (let i = 0; i < state.matrix.length; i++) {
                    state.matrix[i].pop();
                }
                renderTableEditor();
                updateAll();
            } else {
                alert("Table must have at least 2 columns for contingency analysis.");
            }
        });

        // Reset Table
        dom.resetTableBtn.addEventListener('click', () => {
            const p = PRESETS[state.activePreset];
            state.rowLabels = [...p.rowLabels];
            state.colLabels = [...p.colLabels];
            state.matrix = p.matrix.map(row => [...row]);
            renderTableEditor();
            updateAll();
        });

        // Table Perspective Tabs
        document.querySelectorAll('.table-tab-btn').forEach(btn => {
            btn.addEventListener('click', (e) => {
                document.querySelectorAll('.table-tab-btn').forEach(b => b.classList.remove('active'));
                e.target.classList.add('active');
                state.activeTab = e.target.dataset.tab;
                if (state.analysisResults) {
                    renderActiveTableDisplay(state.analysisResults);
                }
            });
        });

        // CSV Parser
        if (dom.parseCsvBtn) {
            dom.parseCsvBtn.addEventListener('click', () => {
                const text = dom.csvInputArea.value.trim();
                if (!text) return;
                const lines = text.split('\n').map(l => l.trim()).filter(l => l.length > 0);
                if (lines.length < 2) {
                    alert("Please paste at least two rows of paired data.");
                    return;
                }

                const pairCounts = {};
                const rSet = new Set();
                const cSet = new Set();

                lines.forEach(line => {
                    const parts = line.split(/[,\t]+/).map(p => p.trim());
                    if (parts.length >= 2) {
                        const rVal = parts[0];
                        const cVal = parts[1];
                        rSet.add(rVal);
                        cSet.add(cVal);
                        const key = `${rVal}|||${cVal}`;
                        pairCounts[key] = (pairCounts[key] || 0) + 1;
                    }
                });

                const rArr = Array.from(rSet);
                const cArr = Array.from(cSet);
                if (rArr.length < 2 || cArr.length < 2) {
                    alert("Parsed data must contain at least 2 distinct row levels and 2 distinct column levels.");
                    return;
                }

                const newMatrix = [];
                for (let i = 0; i < rArr.length; i++) {
                    newMatrix[i] = [];
                    for (let j = 0; j < cArr.length; j++) {
                        const key = `${rArr[i]}|||${cArr[j]}`;
                        newMatrix[i][j] = pairCounts[key] || 0;
                    }
                }

                state.rowLabels = rArr;
                state.colLabels = cArr;
                state.rowVar = "Variable X (Imported)";
                state.colVar = "Variable Y (Imported)";
                state.matrix = newMatrix;
                dom.rowVarInput.value = state.rowVar;
                dom.colVarInput.value = state.colVar;

                renderTableEditor();
                updateAll();
            });
        }

        // Synthetic Generator
        if (dom.synthDepSlider) {
            dom.synthDepSlider.addEventListener('input', (e) => {
                dom.synthDepVal.textContent = parseFloat(e.target.value).toFixed(2);
            });
        }

        if (dom.generateSynthBtn) {
            dom.generateSynthBtn.addEventListener('click', () => {
                const N = parseInt(dom.synthNInput.value, 10) || 400;
                const dep = parseFloat(dom.synthDepSlider.value);
                const r = 3, c = 3;
                const rLabels = ["Feature Alpha", "Feature Beta", "Feature Gamma"];
                const cLabels = ["Response 1", "Response 2", "Response 3"];
                const newMatrix = [ [0,0,0], [0,0,0], [0,0,0] ];

                for (let k = 0; k < N; k++) {
                    let rIdx = Math.floor(Math.random() * r);
                    let cIdx;
                    if (Math.random() < dep) {
                        // Coupled deterministically with probability 'dep'
                        cIdx = rIdx;
                    } else {
                        // Independent random assignment
                        cIdx = Math.floor(Math.random() * c);
                    }
                    newMatrix[rIdx][cIdx]++;
                }

                state.rowLabels = rLabels;
                state.colLabels = cLabels;
                state.rowVar = "Synthetic Generator (X)";
                state.colVar = "Synthetic Response (Y)";
                state.matrix = newMatrix;
                dom.rowVarInput.value = state.rowVar;
                dom.colVarInput.value = state.colVar;

                renderTableEditor();
                updateAll();
            });
        }
    }

    function escapeHTML(str) {
        if (!str) return '';
        return String(str).replace(/[&<>"']/g, (m) => {
            switch (m) {
                case '&': return '&amp;';
                case '<': return '&lt;';
                case '>': return '&gt;';
                case '"': return '&quot;';
                case "'": return '&#39;';
                default: return m;
            }
        });
    }

    // =========================================================================
    // 10. DOM READY BOOTSTRAP
    // =========================================================================

    window.addEventListener('DOMContentLoaded', () => {
        initDOM();
        attachEventListeners();
        renderTableEditor();
        updateAll();
        renderSimpsonsParadox();
    });

})();
