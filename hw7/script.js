/**
 * Homework 7: Statistical Independence & Bivariate Contingency Analysis
 * Author: Dario Pietrosanto (Matricola: 2055042)
 * Master's Degree Course in Applied Statistics & Cybersecurity
 *
 * Comprehensive Modules:
 * - Programmatic Contingency Cross-Tabulation (Observed, Expected, Marginals, Conditionals)
 * - Pearson Residuals, Haberman's Adjusted Standardized Residuals, and Cell Chi-Square Contributions
 * - Inferential Hypothesis Testing: Pearson's Chi-Square, Likelihood Ratio Deviance (G^2),
 *   Fisher's Exact Test, and Measures of Association (Cramer's V, Tschuprow's T, Contingency Coeff C)
 * - Exact P-Value via Incomplete Gamma Function (Lanczos + Lentz Continued Fraction)
 * - Non-Parametric Monte Carlo Permutation Test (B = 1,000+ Shuffles)
 * - Dynamic APA-Style Statistical Interpretation & Residual Diagnostics Generator
 * - Interactive Mosaic / Spineplot Canvas with Residual Coloring
 * - Theoretical Chi-Square Density & Critical Rejection Region Canvas
 * - Conditional vs. Marginal Distribution Comparison Canvas (Visual Independence Proof)
 * - Causal DAG Studio: Fork (Confounder), Collider (Berkson's Fallacy), Chain, and Mediator
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
        let z = 1.6448536;
        if (Math.abs(alpha - 0.10) < 1e-4) z = 1.2815516;
        else if (Math.abs(alpha - 0.05) < 1e-4) z = 1.6448536;
        else if (Math.abs(alpha - 0.01) < 1e-4) z = 2.3263479;
        else if (Math.abs(alpha - 0.001) < 1e-4) z = 3.0902323;

        let x = df * Math.pow(1.0 - 2.0 / (9.0 * df) + z * Math.sqrt(2.0 / (9.0 * df)), 3.0);
        if (x <= 0) x = 0.5;

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

    /**
     * Monte Carlo Permutation Test for Contingency Table Independence
     * Computes non-parametric empirical p-value by conditioning on marginal totals
     */
    function monteCarloPermutationTest(matrix, rSums, cSums, N, chi2Obs, B = 1000) {
        const r = matrix.length;
        const c = matrix[0].length;
        
        const colLabelsArr = new Int32Array(N);
        let idx = 0;
        for (let j = 0; j < c; j++) {
            for (let k = 0; k < cSums[j]; k++) {
                colLabelsArr[idx++] = j;
            }
        }

        const expected = [];
        for (let i = 0; i < r; i++) {
            expected[i] = new Float64Array(c);
            for (let j = 0; j < c; j++) {
                expected[i][j] = (rSums[i] * cSums[j]) / N;
            }
        }

        let extremeCount = 0;
        const permMatrix = [];
        for (let i = 0; i < r; i++) permMatrix[i] = new Int32Array(c);

        for (let b = 0; b < B; b++) {
            // Fisher-Yates shuffle
            for (let k = N - 1; k > 0; k--) {
                const j = Math.floor(Math.random() * (k + 1));
                const tmp = colLabelsArr[k];
                colLabelsArr[k] = colLabelsArr[j];
                colLabelsArr[j] = tmp;
            }

            for (let i = 0; i < r; i++) {
                for (let j = 0; j < c; j++) permMatrix[i][j] = 0;
            }

            let pIdx = 0;
            for (let i = 0; i < r; i++) {
                const rowCount = rSums[i];
                for (let k = 0; k < rowCount; k++) {
                    const colCat = colLabelsArr[pIdx++];
                    permMatrix[i][colCat]++;
                }
            }

            let curChi2 = 0;
            for (let i = 0; i < r; i++) {
                for (let j = 0; j < c; j++) {
                    const diff = permMatrix[i][j] - expected[i][j];
                    curChi2 += (diff * diff) / expected[i][j];
                }
            }

            if (curChi2 >= chi2Obs - 1e-9) {
                extremeCount++;
            }
        }

        return (1.0 + extremeCount) / (1.0 + B);
    }

    // =========================================================================
    // 2. PRE-LOADED DATASETS & PARSING UTILITIES
    // =========================================================================

    const PRESETS = {
        attack_layer: {
            title: "Cybersecurity Incidents: Attack Vector vs. Target Layer",
            description: "Cross-tabulation of 460 enterprise security incidents categorized by entry vector and compromised infrastructure layer.",
            rowVar: "Attack Vector (X)",
            colVar: "Target Infrastructure Layer (Y)",
            rowLabels: ["Phishing / Social Eng.", "SQLi / Web Exploit", "Zero-Day Remote Code", "DDoS Volumetric", "Insider Privilege Abuse"],
            colLabels: ["Web DMZ Tier", "Application API", "Database Cluster", "Cloud IAM / Secrets"],
            matrix: [
                [75, 20,  5, 10],  // Phishing -> heavily Web & Secrets (sum: 110)
                [25, 60, 45,  5],  // SQLi -> App & DB (sum: 135)
                [10, 35, 15, 25],  // Zero-Day -> distributed (sum: 85)
                [40,  5,  0,  0],  // DDoS -> Web DMZ (sum: 45)
                [ 5, 10, 30, 40]   // Insider -> DB & Cloud IAM (sum: 85)
            ]
        },
        malware_persist: {
            title: "Endpoint Threat Forensics: Family vs. Persistence Technique",
            description: "Forensic telemetry of 320 endpoints detailing threat family vs. MITRE ATT&CK persistence technique.",
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
            title: "Identity Protection: Authentication Policy vs. Compromise Outcome",
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
            title: "Theoretical Reference: Strictly Independent Bivariate Model",
            description: "Synthetic sample where joint probabilities strictly equal product of marginals P(X, Y) = P(X)P(Y). Demonstrates chi-square null hypothesis.",
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
            title: "Simpson's Paradox in Security Remediation (Confounded Aggregate)",
            description: "Aggregate comparison of patch deployment strategy vs. breach outcome. Appears Automated Patching has higher failure rate!",
            rowVar: "Patch Deployment Protocol (X)",
            colVar: "System Outcome (Y)",
            rowLabels: ["Automated Pipeline", "Manual Sysadmin"],
            colLabels: ["Secured", "Breached"],
            matrix: [
                [260, 140], // Automated: 140 / 400 = 35.0% breached
                [136,  64]  // Manual:      64 / 200 = 32.0% breached
            ],
            stratified: {
                confounderVar: "Environment Exposure (Z)",
                strata: [
                    {
                        name: "Subgroup 1: External Facing / High Exposure Tier",
                        matrix: [
                            [ 80, 120], // Automated: 120 / 200 = 60.0% breached
                            [ 10,  40]  // Manual:      40 / 50  = 80.0% breached (Automated is 20% BETTER!)
                        ]
                    },
                    {
                        name: "Subgroup 2: Internal Isolated / Low Exposure Tier",
                        matrix: [
                            [180,  20], // Automated:  20 / 200 = 10.0% breached
                            [126,  24]  // Manual:      24 / 150 = 16.0% breached (Automated is 6% BETTER!)
                        ]
                    }
                ]
            }
        }
    };

    // =========================================================================
    // 3. APPLICATION STATE & DOM REFERENCES
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
        activeTab: "observed", // 'observed' | 'relative' | 'row_cond' | 'col_cond' | 'expected' | 'residuals' | 'chi2_contrib'
        activeDagModel: "fork" // 'fork' | 'collider' | 'chain' | 'mediator'
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

        // Containers
        dom.tableEditorContainer = document.getElementById('table-editor-container');
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

        // Interpretation Narrative Card
        dom.narrativeInterpretationCard = document.getElementById('narrative-interpretation-card');

        // Monte Carlo Permutation Test
        dom.runPermutationBtn = document.getElementById('run-permutation-btn');
        dom.permResultDisplay = document.getElementById('perm-result-display');

        // Canvases
        dom.canvasMosaic = document.getElementById('canvas-mosaic');
        dom.canvasChi2Dist = document.getElementById('canvas-chi2-dist');
        dom.canvasResidualHeatmap = document.getElementById('canvas-residual-heatmap');
        dom.canvasCondMarginal = document.getElementById('canvas-cond-marginal');
        dom.canvasCausalDag = document.getElementById('canvas-causal-dag');
        dom.dagSelect = document.getElementById('dag-select');
        dom.dagDescription = document.getElementById('dag-description');

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

        if (N <= 0) return null;

        const jointProbs = [];
        const rowCondProbs = [];
        const colCondProbs = [];
        const expected = [];
        const rawResiduals = [];
        const pearsonResiduals = [];
        const adjustedResiduals = [];
        const chi2Contributions = [];

        let chi2 = 0.0;
        let G2 = 0.0;
        let cellsUnder5 = 0;
        let minExpected = Infinity;

        for (let i = 0; i < r; i++) {
            jointProbs[i] = [];
            rowCondProbs[i] = [];
            expected[i] = [];
            rawResiduals[i] = [];
            pearsonResiduals[i] = [];
            adjustedResiduals[i] = [];
            chi2Contributions[i] = [];

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
                const cellChi2 = pRes * pRes;
                chi2Contributions[i][j] = cellChi2;
                chi2 += cellChi2;

                if (o > 0 && e > 0) {
                    G2 += 2.0 * o * Math.log(o / e);
                }

                const pRow = rowSums[i] / N;
                const pCol = colSums[j] / N;
                const denom = Math.sqrt(e * (1.0 - pRow) * (1.0 - pCol));
                adjustedResiduals[i][j] = denom > 0 ? diff / denom : 0.0;
            }
        }

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

        const minDim = Math.min(r - 1, c - 1);
        const cramersV = minDim > 0 ? Math.sqrt(chi2 / (N * minDim)) : 0.0;
        const tschuprowT = (r > 1 && c > 1) ? Math.sqrt(chi2 / (N * Math.sqrt((r - 1) * (c - 1)))) : 0.0;
        const contingencyC = Math.sqrt(chi2 / (chi2 + N));
        const mutualInfoNats = G2 / (2.0 * N);
        const mutualInfoBits = mutualInfoNats / Math.LN2;

        let yatesChi2 = null, yatesPVal = null, fisherPVal = null;
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
            expected, rawResiduals, pearsonResiduals, adjustedResiduals, chi2Contributions,
            chi2, pValue, criticalValue, isIndependent,
            G2, mutualInfoNats, mutualInfoBits,
            cramersV, tschuprowT, contingencyC,
            cellsUnder5, minExpected, totalCells: r * c,
            percentUnder5: (cellsUnder5 / (r * c)) * 100.0,
            yatesChi2, yatesPVal, fisherPVal
        };
    }

    /**
     * Generate dynamic APA-style Statistical Interpretation Narrative
     */
    function generateInterpretationReport(res) {
        if (!res) return "";
        const isSig = res.pValue < res.alpha;
        const pFormatted = res.pValue < 0.0001 ? "p < .0001" : `p = ${res.pValue.toFixed(4)}`;
        const effectDescriptor = res.cramersV > 0.5 ? "a very strong" : (res.cramersV > 0.3 ? "a moderate-to-strong" : (res.cramersV > 0.1 ? "a small-to-moderate" : "a negligible"));

        // Identify most extreme residuals
        let maxPos = { z: -Infinity, r: 0, c: 0 };
        let maxNeg = { z: Infinity, r: 0, c: 0 };
        for (let i = 0; i < res.r; i++) {
            for (let j = 0; j < res.c; j++) {
                const z = res.adjustedResiduals[i][j];
                if (z > maxPos.z) maxPos = { z, r: i, c: j };
                if (z < maxNeg.z) maxNeg = { z, r: i, c: j };
            }
        }

        let narrative = `
            <div style="font-size:0.95rem; line-height:1.6; color:#2c3e50; text-align:left; word-break:normal; overflow-wrap:break-word;">
                <p style="text-align:left; margin-bottom:0.8rem;">
                    <strong>Formal APA-Style Synthesis:</strong><br>
                    A Pearson's chi-square test of independence was evaluated to examine the statistical relationship between 
                    <strong>${escapeHTML(state.rowVar)}</strong> and <strong>${escapeHTML(state.colVar)}</strong> across $N = ${res.N.toLocaleString()}$ total observations. 
                    The association was found to be <strong>${isSig ? 'statistically significant' : 'not statistically significant'}</strong>, 
                    $\\chi^2(${res.df}, N = ${res.N}) = ${res.chi2.toFixed(2)}$, $${pFormatted}$, Cramér's $V = ${res.cramersV.toFixed(3)}$ (representing ${effectDescriptor} association).
                </p>
                <p style="text-align:left; margin-bottom:0.8rem;">
                    <strong>Substantive Interpretation:</strong><br>
                    ${isSig 
                        ? `Because the test statistic exceeds the critical threshold $(\\chi^2_{\\text{obs}} = ${res.chi2.toFixed(2)} > \\chi^2_{\\text{crit}} = ${res.criticalValue.toFixed(2)})$, we reject the null hypothesis of independence ($H_0$). Knowledge of an entity's category in ${escapeHTML(state.rowVar)} substantially alters the conditional probability distribution across ${escapeHTML(state.colVar)}.`
                        : `Because the test statistic does not exceed the critical threshold $(\\chi^2_{\\text{obs}} = ${res.chi2.toFixed(2)} \\le \\chi^2_{\\text{crit}} = ${res.criticalValue.toFixed(2)})$, we fail to reject the null hypothesis ($H_0$). The empirical deviations between observed and expected cell counts are completely compatible with random sampling variability.`
                    }
                </p>
                ${isSig ? `
                <div style="background:#eafaf1; border-left:4px solid #27ae60; padding:0.8rem; border-radius:4px; margin-top:0.6rem; text-align:left; word-break:normal; overflow-wrap:break-word;">
                    <strong>Key Localized Residual Drivers:</strong><br>
                    Inspection of Haberman's adjusted standardized residuals identifies the strongest positive dependency at cell 
                    <strong>[${escapeHTML(res.rowLabels[maxPos.r])} × ${escapeHTML(res.colLabels[maxPos.c])}]</strong> with $z = +${maxPos.z.toFixed(2)}$ ($n = ${res.matrix[maxPos.r][maxPos.c]}$ vs. $e = ${res.expected[maxPos.r][maxPos.c].toFixed(1)}$), 
                    and the strongest negative deficit at cell 
                    <strong>[${escapeHTML(res.rowLabels[maxNeg.r])} × ${escapeHTML(res.colLabels[maxNeg.c])}]</strong> with $z = ${maxNeg.z.toFixed(2)}$ ($n = ${res.matrix[maxNeg.r][maxNeg.c]}$ vs. $e = ${res.expected[maxNeg.r][maxNeg.c].toFixed(1)}$).
                </div>
                ` : ''}
            </div>
        `;
        return narrative;
    }

    // =========================================================================
    // 5. RENDERING: TABLES & CONTINGENCY EDITORS
    // =========================================================================

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

    function renderActiveTableDisplay(res) {
        if (!dom.activeTableDisplay || !res) return;
        const r = res.r;
        const c = res.c;
        const tab = state.activeTab;

        let tableTitle = "", tableSubtitle = "";

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
            tableSubtitle = "Frequencies that would arise under statistical independence. Highlights Cochran criteria (< 5).";
        } else if (tab === "residuals") {
            tableTitle = "Standardized Haberman Residuals z_ij = (n_ij - e_ij) / SE";
            tableSubtitle = "Z-scores under H0 ~ N(0, 1). Absolute values > 1.96 (p < 0.05) or > 2.58 (p < 0.01) indicate localized structural attraction (+) or repulsion (-).";
        } else if (tab === "chi2_contrib") {
            tableTitle = "Cell Contributions to Total Chi-Square Statistic ((n_ij - e_ij)² / e_ij)";
            tableSubtitle = "Component breakdown showing the exact additive contribution of each cell to χ² = " + res.chi2.toFixed(2) + ".";
        }

        let html = `<div style="margin-bottom:0.6rem;"><strong>${tableTitle}</strong><br><small style="color:#7f8c8d;">${tableSubtitle}</small></div>`;
        html += '<table class="contingency-display-table"><thead><tr>';
        html += `<th style="background:#2c3e50; color:#ffffff;">${escapeHTML(state.rowVar)} \\ ${escapeHTML(state.colVar)}</th>`;
        for (let j = 0; j < c; j++) {
            html += `<th>${escapeHTML(res.colLabels[j])}</th>`;
        }
        if (tab !== "residuals" && tab !== "chi2_contrib") {
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
                } else if (tab === "chi2_contrib") {
                    const cVal = res.chi2Contributions[i][j];
                    const pctTotal = res.chi2 > 0 ? ((cVal / res.chi2) * 100).toFixed(1) : 0;
                    if (cVal > 5.0) cellStyle = "background:#fcf3cf; font-weight:bold; color:#b7950b;";
                    cellContent = `${cVal.toFixed(2)}<br><small style="color:#7f8c8d;">(${pctTotal}%)</small>`;
                }

                html += `<td style="${cellStyle}">${cellContent}</td>`;
            }

            if (tab === "observed" || tab === "expected") {
                html += `<td style="font-weight:bold; background:#edf2f7;">${res.rowSums[i]}</td>`;
            } else if (tab === "relative") {
                html += `<td style="font-weight:bold; background:#edf2f7;">${((res.rowSums[i] / res.N) * 100).toFixed(2)}%</td>`;
            } else if (tab === "row_cond") {
                html += `<td style="font-weight:bold; background:#edf2f7;">100.00%</td>`;
            } else if (tab === "col_cond") {
                html += `<td style="font-weight:bold; background:#edf2f7;">-</td>`;
            }

            html += '</tr>';
        }

        if (tab !== "residuals" && tab !== "chi2_contrib") {
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

    function attachEditorListeners() {
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

        document.querySelectorAll('.row-header-input').forEach(input => {
            input.addEventListener('change', (e) => {
                const rIdx = parseInt(e.target.dataset.row, 10);
                state.rowLabels[rIdx] = e.target.value.trim() || `Row ${rIdx + 1}`;
                updateAll();
            });
        });

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
     * Canvas 1: Mosaic Plot
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

        const r = res.r, c = res.c, N = res.N;

        ctx.fillStyle = "#fafbfc";
        ctx.fillRect(padLeft, padTop, plotW, plotH);

        let curX = padLeft;
        for (let j = 0; j < c; j++) {
            const colFrac = res.colSums[j] / N;
            const colWidth = colFrac * plotW;

            let curY = padTop + plotH;
            for (let i = r - 1; i >= 0; i--) {
                const condFrac = res.colCondProbs[i][j];
                const sliceH = condFrac * plotH;
                curY -= sliceH;

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

                if (colWidth > 35 && sliceH > 18) {
                    ctx.fillStyle = Math.abs(z) > 1.96 ? "#ffffff" : "#2c3e50";
                    ctx.font = "bold 10px sans-serif";
                    ctx.textAlign = "center";
                    ctx.fillText(res.matrix[i][j], curX + colWidth / 2, curY + sliceH / 2 + 4);
                }
            }

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

        ctx.fillStyle = "#2c3e50";
        ctx.font = "bold 11px sans-serif";
        ctx.textAlign = "left";
        ctx.fillText(`Mosaic Plot: Width = P(${res.colVar}), Height = P(${res.rowVar} | ${res.colVar})`, padLeft, padTop - 12);
    }

    /**
     * Canvas 2: Chi-Square Density & Critical Region
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

        const xMax = Math.max(16.0, chi2Crit * 1.5, chi2Obs * 1.3, df * 3.0);
        const xMin = 0.0;

        let peakY = 0.0;
        const modeX = Math.max(0.1, df - 2);
        peakY = Math.max(0.15, chiSquarePDF(modeX, df) * 1.25);
        if (df === 1) peakY = 0.8;

        const mapX = (x) => padLeft + ((x - xMin) / (xMax - xMin)) * plotW;
        const mapY = (y) => padTop + (1.0 - Math.min(1.0, y / peakY)) * plotH;

        ctx.fillStyle = "#fafbfc";
        ctx.fillRect(padLeft, padTop, plotW, plotH);

        ctx.strokeStyle = "#e1e4e8";
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let g = 0; g <= 4; g++) {
            const gy = padTop + (g / 4) * plotH;
            ctx.moveTo(padLeft, gy);
            ctx.lineTo(padLeft + plotW, gy);
        }
        ctx.stroke();

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

        const obsX = mapX(Math.min(xMax, chi2Obs));
        ctx.strokeStyle = "#2980b9";
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        ctx.moveTo(obsX, padTop + 20);
        ctx.lineTo(obsX, padTop + plotH);
        ctx.stroke();

        ctx.fillStyle = "#2980b9";
        ctx.beginPath();
        ctx.arc(obsX, padTop + 20, 5, 0, 2 * Math.PI);
        ctx.fill();

        ctx.fillStyle = "#2980b9";
        ctx.font = "bold 10.5px sans-serif";
        ctx.textAlign = (obsX > padLeft + plotW - 70) ? "right" : "left";
        const textOffset = (obsX > padLeft + plotW - 70) ? -8 : 8;
        ctx.fillText(`χ²_obs = ${chi2Obs.toFixed(2)} (p = ${res.pValue < 0.0001 ? '<0.0001' : res.pValue.toFixed(4)})`, obsX + textOffset, padTop + 24);

        ctx.strokeStyle = "#2c3e50";
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(padLeft, padTop);
        ctx.lineTo(padLeft, padTop + plotH);
        ctx.lineTo(padLeft + plotW, padTop + plotH);
        ctx.stroke();

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

        const r = res.r, c = res.c;
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

            ctx.fillStyle = "#2c3e50";
            ctx.font = "10px sans-serif";
            ctx.textAlign = "right";
            const rowLabelShort = res.rowLabels[i].length > 14 ? res.rowLabels[i].slice(0, 12) + ".." : res.rowLabels[i];
            ctx.fillText(rowLabelShort, padLeft - 8, padTop + i * cellH + cellH / 2 + 4);
        }

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

    /**
     * Canvas 4: Conditional vs. Marginal Distribution (Visual Independence Proof)
     * Demonstrates that if X and Y are independent, all conditional profiles P(Y | X = x_i)
     * align identically with the marginal profile P(Y).
     */
    function renderCondMarginalCanvas(res) {
        const canvas = dom.canvasCondMarginal;
        if (!canvas || !res) return;
        const ctx = canvas.getContext('2d');
        const W = canvas.width, H = canvas.height;
        ctx.clearRect(0, 0, W, H);

        const padLeft = 50, padRight = 30, padTop = 30, padBottom = 45;
        const plotW = W - padLeft - padRight;
        const plotH = H - padTop - padBottom;

        const c = res.c;
        const r = res.r;

        ctx.fillStyle = "#fafbfc";
        ctx.fillRect(padLeft, padTop, plotW, plotH);

        // Group width for each column category Y_j
        const groupW = plotW / c;
        const totalBarsPerGroup = r + 1; // 1 marginal + r conditionals
        const barW = Math.max(4, Math.min(22, (groupW * 0.75) / totalBarsPerGroup));

        // Colors for bars: Black for marginal, distinct hues for conditional rows
        const rowColors = [
            "#3498db", "#e67e22", "#9b59b6", "#1abc9c", "#e74c3c", "#34495e"
        ];

        // Draw grid
        ctx.strokeStyle = "#e1e4e8";
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let g = 0; g <= 4; g++) {
            const gy = padTop + (g / 4) * plotH;
            ctx.moveTo(padLeft, gy);
            ctx.lineTo(padLeft + plotW, gy);
            const pct = (1.0 - g / 4) * 100;
            ctx.fillStyle = "#7f8c8d";
            ctx.font = "9px monospace";
            ctx.textAlign = "right";
            ctx.fillText(`${pct.toFixed(0)}%`, padLeft - 6, gy + 3);
        }
        ctx.stroke();

        for (let j = 0; j < c; j++) {
            const groupCenterX = padLeft + j * groupW + groupW / 2;
            const startX = groupCenterX - (totalBarsPerGroup * barW) / 2;

            // 1. Marginal baseline bar P(Y = y_j)
            const margP = res.colSums[j] / res.N;
            const mH = margP * plotH;
            const mX = startX;
            const mY = padTop + plotH - mH;

            ctx.fillStyle = "#2c3e50";
            ctx.fillRect(mX, mY, barW, mH);
            ctx.strokeStyle = "#1a252f";
            ctx.lineWidth = 1;
            ctx.strokeRect(mX, mY, barW, mH);

            // 2. Conditional bars P(Y = y_j | X = x_i)
            for (let i = 0; i < r; i++) {
                const condP = res.rowCondProbs[i][j];
                const cH = condP * plotH;
                const cX = startX + (i + 1) * barW;
                const cY = padTop + plotH - cH;

                ctx.fillStyle = rowColors[i % rowColors.length];
                ctx.fillRect(cX, cY, barW, cH);
                ctx.strokeStyle = "#ffffff";
                ctx.lineWidth = 0.5;
                ctx.strokeRect(cX, cY, barW, cH);
            }

            // Category label
            ctx.fillStyle = "#2c3e50";
            ctx.font = "10px sans-serif";
            ctx.textAlign = "center";
            const colLabelShort = res.colLabels[j].length > 12 ? res.colLabels[j].slice(0, 10) + ".." : res.colLabels[j];
            ctx.fillText(colLabelShort, groupCenterX, padTop + plotH + 16);
        }

        // Axes
        ctx.strokeStyle = "#2c3e50";
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(padLeft, padTop);
        ctx.lineTo(padLeft, padTop + plotH);
        ctx.lineTo(padLeft + plotW, padTop + plotH);
        ctx.stroke();

        ctx.fillStyle = "#2c3e50";
        ctx.font = "bold 11px sans-serif";
        ctx.textAlign = "left";
        ctx.fillText("Conditional Profiles P(Y | X) vs. Marginal Baseline P(Y) [Black Bar]", padLeft, padTop - 12);
    }

    /**
     * Canvas 5: Interactive Causal DAG Studio
     * Visualizes Fork (Confounder), Collider (Berkson's Paradox), Chain, and Mediator
     */
    function renderCausalDagCanvas() {
        const canvas = dom.canvasCausalDag;
        if (!canvas) return;
        const ctx = canvas.getContext('2d');
        const W = canvas.width, H = canvas.height;
        ctx.clearRect(0, 0, W, H);

        ctx.fillStyle = "#fafbfc";
        ctx.fillRect(0, 0, W, H);

        const model = state.activeDagModel;

        function drawNode(x, y, label, subtext, color = "#2980b9") {
            ctx.fillStyle = color;
            ctx.beginPath();
            ctx.arc(x, y, 24, 0, 2 * Math.PI);
            ctx.fill();

            ctx.strokeStyle = "#2c3e50";
            ctx.lineWidth = 2;
            ctx.stroke();

            ctx.fillStyle = "#ffffff";
            ctx.font = "bold 13px sans-serif";
            ctx.textAlign = "center";
            ctx.fillText(label, x, y + 5);

            if (subtext) {
                ctx.fillStyle = "#2c3e50";
                ctx.font = "10px sans-serif";
                ctx.fillText(subtext, x, y + 38);
            }
        }

        function drawArrow(fromX, fromY, toX, toY, label = "") {
            const headlen = 10;
            const angle = Math.atan2(toY - fromY, toX - fromX);

            // Offset from node circles (radius 24)
            const r = 24;
            const startX = fromX + r * Math.cos(angle);
            const startY = fromY + r * Math.sin(angle);
            const endX = toX - r * Math.cos(angle);
            const endY = toY - r * Math.sin(angle);

            ctx.strokeStyle = "#34495e";
            ctx.lineWidth = 2.5;
            ctx.beginPath();
            ctx.moveTo(startX, startY);
            ctx.lineTo(endX, endY);
            ctx.stroke();

            // Arrow head
            ctx.fillStyle = "#34495e";
            ctx.beginPath();
            ctx.moveTo(endX, endY);
            ctx.lineTo(endX - headlen * Math.cos(angle - Math.PI / 6), endY - headlen * Math.sin(angle - Math.PI / 6));
            ctx.lineTo(endX - headlen * Math.cos(angle + Math.PI / 6), endY - headlen * Math.sin(angle + Math.PI / 6));
            ctx.closePath();
            ctx.fill();

            if (label) {
                ctx.fillStyle = "#7f8c8d";
                ctx.font = "10px sans-serif";
                ctx.textAlign = "center";
                ctx.fillText(label, (fromX + toX) / 2, (fromY + toY) / 2 - 8);
            }
        }

        if (model === "fork") {
            // Confounder Fork: X <- Z -> Y
            drawNode(W / 2, 60, "Z", "Common Cause (Confounder)", "#e67e22");
            drawNode(W / 4, 180, "X", "Exposure / Strategy", "#2980b9");
            drawNode((3 * W) / 4, 180, "Y", "Outcome / Breach", "#27ae60");

            drawArrow(W / 2, 60, W / 4, 180, "Causes X");
            drawArrow(W / 2, 60, (3 * W) / 4, 180, "Causes Y");

            // Dashed spurious association arc
            ctx.setLineDash([5, 5]);
            ctx.strokeStyle = "#e74c3c";
            ctx.lineWidth = 2;
            ctx.beginPath();
            ctx.moveTo(W / 4 + 24, 180);
            ctx.lineTo((3 * W) / 4 - 24, 180);
            ctx.stroke();
            ctx.setLineDash([]);

            ctx.fillStyle = "#c0392b";
            ctx.font = "bold 11px sans-serif";
            ctx.textAlign = "center";
            ctx.fillText("Spurious Association (Backdoor Path X ← Z → Y)", W / 2, 215);

            dom.dagDescription.innerHTML = "<strong>Confounder Fork ($X \\leftarrow Z \\rightarrow Y$):</strong> $Z$ acts as a common cause. Observational correlation between $X$ and $Y$ is biased by the backdoor path. Condition on $Z$ or execute an intervention $\\text{do}(X)$ to measure the true causal effect.";
        } else if (model === "collider") {
            // Collider: X -> Z <- Y
            drawNode(W / 4, 60, "X", "Skill / Feature", "#2980b9");
            drawNode((3 * W) / 4, 60, "Y", "Effort / Feature", "#27ae60");
            drawNode(W / 2, 180, "Z", "Collider (e.g. Success)", "#8e44ad");

            drawArrow(W / 4, 60, W / 2, 180, "");
            drawArrow((3 * W) / 4, 60, W / 2, 180, "");

            ctx.setLineDash([5, 5]);
            ctx.strokeStyle = "#e74c3c";
            ctx.lineWidth = 2;
            ctx.beginPath();
            ctx.moveTo(W / 4 + 24, 60);
            ctx.lineTo((3 * W) / 4 - 24, 60);
            ctx.stroke();
            ctx.setLineDash([]);

            ctx.fillStyle = "#c0392b";
            ctx.font = "bold 11px sans-serif";
            ctx.textAlign = "center";
            ctx.fillText("Spurious Dependency Induced Conditioning on Collider Z!", W / 2, 35);

            dom.dagDescription.innerHTML = "<strong>Collider / Inverted Fork ($X \\rightarrow Z \\leftarrow Y$):</strong> $X$ and $Y$ are marginally independent ($X \\perp Y$). Conditioning on collider $Z$ opens the path, creating an artificial negative correlation (Berkson's Fallacy / Selection Bias)!";
        } else if (model === "chain") {
            // Direct Causal Chain: X -> Y
            drawNode(W / 3, 120, "X", "Direct Cause", "#2980b9");
            drawNode((2 * W) / 3, 120, "Y", "Direct Effect", "#27ae60");
            drawArrow(W / 3, 120, (2 * W) / 3, 120, "Direct Causal Link");

            dom.dagDescription.innerHTML = "<strong>Direct Causal Chain ($X \\rightarrow Y$):</strong> $X$ has an unconfounded direct influence on $Y$. Intervening via $\\text{do}(X)$ directly changes the probability distribution of $Y$.";
        } else if (model === "mediator") {
            // Mediator: X -> M -> Y
            drawNode(W / 4, 120, "X", "Exogenous Input", "#2980b9");
            drawNode(W / 2, 120, "M", "Mediator", "#8e44ad");
            drawNode((3 * W) / 4, 120, "Y", "Final Endpoint", "#27ae60");

            drawArrow(W / 4, 120, W / 2, 120, "");
            drawArrow(W / 2, 120, (3 * W) / 4, 120, "");

            dom.dagDescription.innerHTML = "<strong>Mediator Chain ($X \\rightarrow M \\rightarrow Y$):</strong> $X$ causes $Y$ entirely through intermediate mechanism $M$. Conditioning on $M$ d-separates $X$ and $Y$, making them conditionally independent ($X \\perp Y \\mid M$).";
        }

        if (window.renderMathInElement) {
            window.renderMathInElement(dom.dagDescription, {
                delimiters: [{ left: '$', right: '$', display: false }]
            });
        }
    }

    // =========================================================================
    // 7. SIMPSON'S PARADOX CAUSAL CONVOLUTION ENGINE
    // =========================================================================

    function renderSimpsonsParadox() {
        if (!dom.simpsonContainer) return;
        const sp = PRESETS.simpson_paradox;
        const aggRes = computeContingencyAnalysis(sp.matrix, sp.rowLabels, sp.colLabels);

        let html = `
            <div class="stat-callout" style="background:#fcf3cf; border-left:5px solid #f39c12; margin-top:1.5rem;">
                <h4 style="margin-top:0; color:#b7950b;">Interactive Epistemological Laboratory: Simpson's Paradox & The Causal Backdoor</h4>
                <p>
                    <strong>The Aggregation Illusion:</strong> Below, we evaluate an enterprise patch deployment dataset. 
                    In the <em>aggregate</em> table, Automated Patching appears to have a higher breach rate than Manual Patching. 
                    However, when conditioning on the confounding variable $Z$ (<strong>Environment Exposure Tier</strong>), Automated Patching is superior in every single individual environment!
                </p>
                <div class="layout-grid" style="grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap:1rem; margin-top:1rem;">
                    <div style="background:#ffffff; padding:1rem; border:1px solid #d5dbdb; border-radius:4px; box-sizing:border-box;">
                        <h5 style="margin-top:0; color:#c0392b;">1. Aggregate Data (Confounded)</h5>
                        <div style="overflow-x:auto; -webkit-overflow-scrolling:touch;">
                            <table class="subgroup-table">
                                <tr><th>Strategy</th><th>Secured</th><th>Breached</th><th>Breach Rate</th></tr>
                                <tr><td>Automated</td><td>260</td><td>140</td><td style="color:#e74c3c; font-weight:bold;">35.0%</td></tr>
                                <tr><td>Manual</td><td>136</td><td>64</td><td style="color:#27ae60; font-weight:bold;">32.0%</td></tr>
                            </table>
                        </div>
                        <p style="font-size:0.8rem; color:#7f8c8d; margin-top:0.5rem;">
                            <strong>Naive Inference:</strong> Manual seems 3.0% safer in crude aggregate!<br>
                            Chi-square: χ² = ${aggRes.chi2.toFixed(2)}, p = ${aggRes.pValue.toExponential(2)} (Spurious Association).
                        </p>
                    </div>

                    <div style="background:#ffffff; padding:1rem; border:1px solid #d5dbdb; border-radius:4px; box-sizing:border-box;">
                        <h5 style="margin-top:0; color:#2980b9;">2. Subgroup A: High-Exposure DMZ (Z = 1)</h5>
                        <div style="overflow-x:auto; -webkit-overflow-scrolling:touch;">
                            <table class="subgroup-table">
                                <tr><th>Strategy</th><th>Secured</th><th>Breached</th><th>Breach Rate</th></tr>
                                <tr><td>Automated</td><td>80</td><td>120</td><td style="color:#27ae60; font-weight:bold;">60.0%</td></tr>
                                <tr><td>Manual</td><td>10</td><td>40</td><td style="color:#e74c3c; font-weight:bold;">80.0%</td></tr>
                            </table>
                        </div>
                        <p style="font-size:0.8rem; color:#27ae60; font-weight:bold; margin-top:0.5rem;">
                            Automated is 20.0% SAFER in High-Exposure DMZ!
                        </p>
                    </div>

                    <div style="background:#ffffff; padding:1rem; border:1px solid #d5dbdb; border-radius:4px; box-sizing:border-box;">
                        <h5 style="margin-top:0; color:#2980b9;">3. Subgroup B: Isolated Air-Gap (Z = 2)</h5>
                        <div style="overflow-x:auto; -webkit-overflow-scrolling:touch;">
                            <table class="subgroup-table">
                                <tr><th>Strategy</th><th>Secured</th><th>Breached</th><th>Breach Rate</th></tr>
                                <tr><td>Automated</td><td>180</td><td>20</td><td style="color:#27ae60; font-weight:bold;">10.0%</td></tr>
                                <tr><td>Manual</td><td>126</td><td>24</td><td style="color:#e74c3c; font-weight:bold;">16.0%</td></tr>
                            </table>
                        </div>
                        <p style="font-size:0.8rem; color:#27ae60; font-weight:bold; margin-top:0.5rem;">
                            Automated is 6.0% SAFER in Isolated Air-Gap!
                        </p>
                    </div>
                </div>

                <div style="background:#ffffff; padding:1.2rem; border-radius:4px; margin-top:1rem; border:1px solid #e1e8ed; overflow-x:auto; box-sizing:border-box;">
                    <strong>Causal Analysis via Pearl's Do-Calculus:</strong><br>
                    The confounding variable $Z$ (Exposure Tier) acts as a common cause (fork: $X \\leftarrow Z \\rightarrow Y$). 
                    Because the Automated strategy was disproportionately assigned to the most hazardous internet-facing machines ($P(Z=1 \\mid \\text{Automated}) = 200/400 = 50.0\\%$ vs $P(Z=1 \\mid \\text{Manual}) = 50/200 = 25.0\\%$), the crude observational association $P(Y \\mid X)$ contains severe backdoor confounding bias. 
                    Applying the <strong>Backdoor Adjustment Formula</strong> ($P(Z=1) = 250/600$, $P(Z=2) = 350/600$):
                    $$\\begin{aligned}
                    P(Y = \\text{Breached} \\mid \\text{do}(X = \\text{Auto})) &= \\sum_z P(Y = \\text{Breached} \\mid X = \\text{Auto}, Z = z) P(Z = z) \\\\
                    &= 0.60 \\left(\\frac{250}{600}\\right) + 0.10 \\left(\\frac{350}{600}\\right) \\approx \\mathbf{30.83\\%}
                    \\end{aligned}$$
                    $$\\begin{aligned}
                    P(Y = \\text{Breached} \\mid \\text{do}(X = \\text{Manual})) &= \\sum_z P(Y = \\text{Breached} \\mid X = \\text{Manual}, Z = z) P(Z = z) \\\\
                    &= 0.80 \\left(\\frac{250}{600}\\right) + 0.16 \\left(\\frac{350}{600}\\right) \\approx \\mathbf{42.67\\%}
                    \\end{aligned}$$
                    Under true intervention ($\\text{do}(X)$), Automated Patching reduces the population breach risk by <strong>11.84 percentage points</strong>!
                </div>
        `;

        dom.simpsonContainer.innerHTML = html;
        if (window.renderMathInElement) {
            window.renderMathInElement(dom.simpsonContainer, {
                delimiters: [{ left: '$$', right: '$$', display: true }, { left: '$', right: '$', display: false }]
            });
        }
    }

    // =========================================================================
    // 8. MASTER UPDATE DISPATCHER
    // =========================================================================

    function updateAll() {
        const res = computeContingencyAnalysis(state.matrix, state.rowLabels, state.colLabels, state.alpha);
        state.analysisResults = res;

        if (!res) return;

        dom.nTotalVal.textContent = res.N.toLocaleString();
        dom.dfVal.textContent = res.df;
        dom.chi2Val.textContent = res.chi2.toFixed(4);
        dom.critValDisplay.textContent = res.criticalValue.toFixed(4);

        if (res.pValue < 0.0001) {
            dom.pValDisplay.textContent = "< 0.0001";
        } else {
            dom.pValDisplay.textContent = res.pValue.toFixed(6);
        }

        if (res.isIndependent) {
            dom.decisionBadge.innerHTML = `<span style="background:#27ae60; color:#fff; padding:6px 14px; border-radius:4px; font-weight:bold; display:inline-block;">STATISTICALLY INDEPENDENT (Fail to Reject H0 at α = ${state.alpha})</span>`;
        } else {
            dom.decisionBadge.innerHTML = `<span style="background:#c0392b; color:#fff; padding:6px 14px; border-radius:4px; font-weight:bold; display:inline-block;">STATISTICALLY DEPENDENT (Reject H0 at α = ${state.alpha})</span>`;
        }

        dom.g2Val.textContent = `${res.G2.toFixed(3)} (${res.mutualInfoBits.toFixed(3)} bits)`;
        dom.cramerVVal.textContent = res.cramersV.toFixed(4);

        let effectStr = "Negligible";
        if (res.cramersV > 0.5) effectStr = "Very Strong Association";
        else if (res.cramersV > 0.3) effectStr = "Moderate-to-Strong";
        else if (res.cramersV > 0.1) effectStr = "Weak Association";
        dom.cramerEffectVal.textContent = effectStr;

        if (res.percentUnder5 > 20.0 || res.minExpected < 1.0) {
            dom.cochranWarning.style.display = "block";
            dom.cochranWarning.innerHTML = `<strong>Warning (Cochran's Criterion Violation):</strong> ${res.cellsUnder5} out of ${res.totalCells} cells (${res.percentUnder5.toFixed(1)}%) have expected frequency e_ij &lt; 5. Minimum e_ij = ${res.minExpected.toFixed(2)}. Chi-square asymptotic approximation may be unreliable; consider Exact tests, Monte Carlo simulation, or category pooling.`;
        } else {
            dom.cochranWarning.style.display = "none";
        }

        if (res.r === 2 && res.c === 2) {
            dom.fisher2x2Card.style.display = "block";
            dom.fisherPVal.textContent = res.fisherPVal < 0.0001 ? "< 0.0001" : res.fisherPVal.toFixed(6);
        } else {
            dom.fisher2x2Card.style.display = "none";
        }

        // Narrative Report
        if (dom.narrativeInterpretationCard) {
            dom.narrativeInterpretationCard.innerHTML = generateInterpretationReport(res);
            if (window.renderMathInElement) {
                window.renderMathInElement(dom.narrativeInterpretationCard, {
                    delimiters: [{ left: '$', right: '$', display: false }]
                });
            }
        }

        // Reset Monte Carlo result
        if (dom.permResultDisplay) {
            dom.permResultDisplay.innerHTML = `<span style="color:#7f8c8d; font-size:0.9rem;">Click button below to execute 1,000 non-parametric shuffles.</span>`;
        }

        renderActiveTableDisplay(res);
        renderMosaicCanvas(res);
        renderChi2DistCanvas(res);
        renderResidualHeatmap(res);
        renderCondMarginalCanvas(res);
        renderCausalDagCanvas();
    }

    // =========================================================================
    // 9. EVENT LISTENERS & INITIALIZATION
    // =========================================================================

    function attachEventListeners() {
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

        dom.alphaSelect.addEventListener('change', (e) => {
            state.alpha = parseFloat(e.target.value);
            updateAll();
        });

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

        dom.addRowBtn.addEventListener('click', () => {
            const newRIdx = state.matrix.length + 1;
            const c = state.matrix[0].length;
            state.rowLabels.push(`Category X${newRIdx}`);
            state.matrix.push(new Array(c).fill(10));
            renderTableEditor();
            updateAll();
        });

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

        dom.addColBtn.addEventListener('click', () => {
            const newCIdx = state.matrix[0].length + 1;
            state.colLabels.push(`Category Y${newCIdx}`);
            for (let i = 0; i < state.matrix.length; i++) {
                state.matrix[i].push(10);
            }
            renderTableEditor();
            updateAll();
        });

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

        dom.resetTableBtn.addEventListener('click', () => {
            const p = PRESETS[state.activePreset];
            state.rowLabels = [...p.rowLabels];
            state.colLabels = [...p.colLabels];
            state.matrix = p.matrix.map(row => [...row]);
            renderTableEditor();
            updateAll();
        });

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

        // Monte Carlo Permutation Test Button
        if (dom.runPermutationBtn) {
            dom.runPermutationBtn.addEventListener('click', () => {
                if (!state.analysisResults) return;
                dom.permResultDisplay.innerHTML = `<span style="color:#2980b9;">Executing 1,000 Monte Carlo shuffles...</span>`;
                setTimeout(() => {
                    const t0 = performance.now();
                    const pMC = monteCarloPermutationTest(
                        state.matrix,
                        state.analysisResults.rowSums,
                        state.analysisResults.colSums,
                        state.analysisResults.N,
                        state.analysisResults.chi2,
                        1000
                    );
                    const elapsed = (performance.now() - t0).toFixed(1);
                    dom.permResultDisplay.innerHTML = `
                        <strong>Empirical Permutation P-Value (B = 1,000):</strong> 
                        <span style="font-family:monospace; color:${pMC < 0.05 ? '#c0392b' : '#27ae60'}; font-size:1.1rem; font-weight:bold;">${pMC.toFixed(4)}</span> 
                        <small style="color:#7f8c8d;">(${elapsed} ms)</small><br>
                        <small style="color:#2c3e50;">Confirms ${pMC < 0.05 ? 'strong rejection of independence' : 'consistency with independence'} without relying on large-sample normal approximations.</small>
                    `;
                }, 10);
            });
        }

        // Causal DAG Model Selector
        if (dom.dagSelect) {
            dom.dagSelect.addEventListener('change', (e) => {
                state.activeDagModel = e.target.value;
                renderCausalDagCanvas();
            });
        }

        // CSV / Triplet Parser
        if (dom.parseCsvBtn) {
            dom.parseCsvBtn.addEventListener('click', () => {
                const text = dom.csvInputArea.value.trim();
                if (!text) return;
                const lines = text.split('\n').map(l => l.trim()).filter(l => l.length > 0);
                if (lines.length < 2) {
                    alert("Please paste at least two rows of observations.");
                    return;
                }

                const pairCounts = {};
                const rSet = new Set();
                const cSet = new Set();

                lines.forEach(line => {
                    const parts = line.split(/[,\t]+/).map(p => p.trim());
                    if (parts.length >= 3 && !isNaN(Number(parts[2]))) {
                        // Weighted triplet format: X, Y, Count
                        const rVal = parts[0], cVal = parts[1], count = Math.max(0, parseInt(parts[2], 10) || 0);
                        rSet.add(rVal); cSet.add(cVal);
                        const key = `${rVal}|||${cVal}`;
                        pairCounts[key] = (pairCounts[key] || 0) + count;
                    } else if (parts.length >= 2) {
                        // Raw pair format: X, Y
                        const rVal = parts[0], cVal = parts[1];
                        rSet.add(rVal); cSet.add(cVal);
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
                const rLabels = ["Threat Level Low", "Threat Level Med", "Threat Level High"];
                const cLabels = ["Response Standard", "Response Elevated", "Response Critical"];
                const newMatrix = [ [0,0,0], [0,0,0], [0,0,0] ];

                for (let k = 0; k < N; k++) {
                    let rIdx = Math.floor(Math.random() * r);
                    let cIdx;
                    if (Math.random() < dep) {
                        cIdx = rIdx;
                    } else {
                        cIdx = Math.floor(Math.random() * c);
                    }
                    newMatrix[rIdx][cIdx]++;
                }

                state.rowLabels = rLabels;
                state.colLabels = cLabels;
                state.rowVar = "Simulated Threat Category (X)";
                state.colVar = "SOC Response Posture (Y)";
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

    window.addEventListener('DOMContentLoaded', () => {
        initDOM();
        attachEventListeners();
        renderTableEditor();
        updateAll();
        renderSimpsonsParadox();
    });

})();
